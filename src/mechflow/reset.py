"""'mechflow reset': put failed jobs (or, when asked, finished ones) back so that they can run again.

Nothing is deleted: the files of each attempt move into failed_attempts/<date>/ inside the job's folder.

    Step 1, OptTS, Forwards, Backwards   restart from the last geometry of the attempt
    SCAN                                 a new part (SCAN_2, SCAN_3 ...) from the last converged point, as in the
                                         procedure (a SCAN with no converged point runs again as it was)
    IRC, single points (SP)              run again from the same structure
    --from-scratch                       the stage starts again from the beginning: 'mechflow next' writes it again
                                         from the previous stage's result
    --imaginary                          optimisations that finished with imaginary frequencies they should not have
                                         (a minimum: any; a TS: more than one) restart from their last geometry pushed
                                         along those modes (read from orca.hess), so that they can leave the saddle

The inputs of a restarted job are written again from the current project.yaml, so more time or memory set there
applies to the next attempt. Jobs that are queued or running in Slurm are never touched.
"""
import glob
import math
import os
import shutil
from collections import namedtuple

from . import chem
from .orca import META, job_ids, job_state, parse_output, read_hess, read_meta, slurm_problem, write_job
from .workflow import (ACCEPTED, ATTEMPTS, IRC, OPTTS, SCAN, VALIDATE, assembled, complex_parts, restart_scan,
                       scan_folders, species_sp_dir, step_progress, step_sp_dir)

RESTART_XYZ = "last_geometry_restart.xyz"     # the input of an optimisation restarted from its last geometry
AMPLITUDE = 0.3          # Angstrom: --imaginary moves the atom that moves most in the imaginary mode(s) this far
STAGES = ("Step1", "SP", "SCAN", "OptTS", "IRC", "Fwd", "Bwd")
STAGE_NAMES = {"step1": ("Step1",), "sp": ("SP",), "scan": ("SCAN",), "optts": ("OptTS",), "irc": ("IRC",),
               "fwd": ("Fwd",), "forwards": ("Fwd",), "bwd": ("Bwd",), "backwards": ("Bwd",), "4b": ("Fwd", "Bwd")}
LATER = {"SCAN": ("OptTS", "IRC", "Fwd", "Bwd", "SP"), "OptTS": ("IRC", "Fwd", "Bwd", "SP"), "IRC": ("Fwd", "Bwd")}

Job = namedtuple("Job", "item kind stage folder")      # kind: species or step


def parse_stages(text):
    """'OptTS,IRC' -> {'OptTS', 'IRC'} (ValueError for an unknown stage)."""
    out = set()
    for word in str(text).replace(" ", ",").split(","):
        if word.strip():
            if word.strip().lower() not in STAGE_NAMES:
                raise ValueError("unknown stage '%s' (stages: %s)" % (word.strip(), ", ".join(STAGES)))
            out.update(STAGE_NAMES[word.strip().lower()])
    return out


def species_jobs(project, sp):
    if assembled(project, sp):            # a complex built from its molecules has no jobs of its own
        return []
    return [Job(sp.id, "species", "Step1", project.species_job_dir(sp.id)),
            Job(sp.id, "species", "SP", species_sp_dir(project, sp.id))]


def step_jobs(project, st):
    sd = project.step_dir(st.id)
    return [Job(st.id, "step", "SCAN", os.path.join(sd, SCAN)),
            Job(st.id, "step", "OptTS", os.path.join(sd, OPTTS)),
            Job(st.id, "step", "IRC", os.path.join(sd, IRC)),
            Job(st.id, "step", "Fwd", os.path.join(sd, VALIDATE, "Forwards")),
            Job(st.id, "step", "Bwd", os.path.join(sd, VALIDATE, "Backwards")),
            Job(st.id, "step", "SP", step_sp_dir(project, st.id))]


def select_jobs(project, species, steps, stages=None):
    jobs = [j for sp in species for j in species_jobs(project, sp)] + [j for st in steps for j in step_jobs(project, st)]
    return [j for j in jobs if not stages or j.stage in stages]


def downstream(project, job):
    """The later jobs that use the result of a job, so must start again when it is redone."""
    if job.kind == "species":
        if job.stage != "Step1":
            return []
        out = [j for j in species_jobs(project, project.species[job.item]) if j.stage == "SP"]
        for st in project.steps.values():
            start = project.species.get(st.start)
            if st.start == job.item or (start is not None and assembled(project, start) and
                                        job.item in _part_ids(project, start)):
                out += step_jobs(project, st)
        return out
    return [j for j in step_jobs(project, project.steps[job.item]) if j.stage in LATER.get(job.stage, ())]


def _part_ids(project, sp):
    try:
        return [pid for pid, _ in complex_parts(project, sp)]
    except ValueError:
        return []


# ============================================================================== imaginary frequencies

def follow_modes(job, folder=None):
    """The imaginary modes --imaginary pushes a finished optimisation along: every one of a minimum (Step 1, Forwards,
    Backwards), every one but the reaction's (the largest) of a TS with more than one. Returns (hess, [mode numbers])
    - hess None when there is no readable orca.hess."""
    folder = folder or job.folder
    if job.stage not in ("Step1", "Fwd", "Bwd", "OptTS"):
        return None, []
    r = parse_output(os.path.join(folder, "output.out"))
    if not r["terminated_normally"] or not r["imag"]:
        return None, []
    hess = read_hess(os.path.join(folder, "orca.hess"))
    if hess is None:
        return None, []
    imag = sorted((f, m) for m, f in enumerate(hess["freqs"]) if f < -1e-6)      # most negative first
    if job.stage == "OptTS":
        imag = imag[1:]                       # keep the reaction's mode
    return hess, [m for _, m in imag]


def displaced(folder, hess, modes, amplitude=AMPLITUDE):
    """The last geometry of an optimisation (orca.xyz) pushed along normal modes: each mode scaled so that its atom that
    moves most moves 1, summed, and the sum scaled so that its atom that moves most moves `amplitude` Angstrom.
    Returns (Molecule, text)."""
    last = chem.read_xyz(os.path.join(folder, "orca.xyz"))
    if last.symbols != hess["symbols"]:
        raise ValueError("orca.hess and orca.xyz of %s have different atoms" % folder)
    # orca.hess may hold the structure moved (or turned): carry the modes into the frame of orca.xyz
    rot, _ = chem.superpose(hess["coords"], last.coords)
    n = len(last)
    total = [[0.0, 0.0, 0.0] for _ in range(n)]
    for m in modes:
        v = [tuple(hess["modes"][m][3 * a + k] for k in range(3)) for a in range(n)]
        v = chem.transform(v, rot, (0.0, 0.0, 0.0))
        big = max(math.sqrt(sum(x * x for x in d)) for d in v) or 1.0
        for a in range(n):
            for k in range(3):
                total[a][k] += v[a][k] / big
    big = max(math.sqrt(sum(x * x for x in d)) for d in total) or 1.0
    coords = [tuple(c[k] + amplitude * total[a][k] / big for k in range(3)) for a, c in enumerate(last.coords)]
    freqs = ", ".join("%.1f" % hess["freqs"][m] for m in modes)
    text = "pushed %.2f A along the imaginary mode%s (%s cm-1)" % (abs(amplitude), "s" if len(modes) > 1 else "", freqs)
    return chem.Molecule(last.symbols, coords), text


def run_folder(job):
    """The folder in which the job runs: for a SCAN, its last part (SCAN_n) once it has been restarted."""
    return scan_folders(os.path.dirname(job.folder))[-1] if job.stage == "SCAN" else job.folder


def state(job, queue=None):
    return job_state(run_folder(job), queue)


def why_failed(folder):
    marker = slurm_problem(folder)
    if marker:
        m = marker.upper()
        return ("time limit" if "TIME LIMIT" in m else "node failure" if "NODE" in m else
                "cancelled" if "CANCELLED" in m else "out of memory")
    r = parse_output(os.path.join(folder, "output.out"))
    if not r["exists"]:
        return "no ORCA output"
    return "ORCA error" if r["error"] else "ORCA stopped"


def _last_geometry(folder, meta):
    """orca.xyz of an optimisation that got somewhere (the same atoms as its input), else None."""
    if meta.get("kind") not in ("opt", "optts"):
        return None
    last = os.path.join(folder, "orca.xyz")
    try:
        a, b = chem.read_xyz(last), chem.read_xyz(os.path.join(folder, meta["xyz"]))
    except (OSError, chem.XYZError, KeyError):
        return None
    return last if (len(a) and a.symbols == b.symbols) else None


def _scan_points_in(folder):
    return sorted(glob.glob(os.path.join(folder, "orca.[0-9][0-9][0-9].xyz")))


# ============================================================================== which jobs

def plan(project, jobs, queue, redo_done=False, from_scratch=False, downstream_too=False, force=False,
         named=False, stages=None, imaginary=False):
    """What 'reset' would do. Returns (actions, notes): actions = [(job, how, why)] in the order of the jobs, with
    how = 'restart', 'scratch' or 'displace' (imaginary=True: only optimisations that finished with imaginary
    frequencies they should not have, pushed along those modes); notes = lines about the jobs that are left alone."""
    known = queue is not None and queue.ids() is not None
    acts, notes, unknown, cache = {}, [], [], {}

    def S(job):
        if job not in cache:
            cache[job] = state(job, queue)
        return cache[job]

    def blocked(job, s, loud):
        """None if the job may be reset, else a note ('' = leave it quietly). A job that squeue lists is never
        reset, whatever its output says; one whose state is unknown (no squeue) only with --force."""
        if known:
            active = job_ids(run_folder(job)) & queue.ids()
            if not active:
                return None
            return ("%s %s is %s in Slurm (job %s): not reset (scancel it first)" % (
                job.item, job.stage, "queued" if s == "submitted" else "running", ", ".join(sorted(active)))
                    ) if loud else ""
        if s not in ("submitted", "running") or force:
            return None
        if not loud:
            unknown.append(job)
            return ""
        return ("%s %s was submitted, and squeue could not be asked whether it is still in Slurm: not reset (add "
                "--force once you know it is not running)" % (job.item, job.stage))

    for job in jobs:
        s = S(job)
        if s in ("absent", "ready"):
            continue
        if imaginary:
            if s != "finished" or job.stage not in ("Step1", "Fwd", "Bwd", "OptTS"):
                continue
            r = parse_output(os.path.join(job.folder, "output.out"))
            extra = r["imag"] if job.stage != "OptTS" else (r["imag"] if len(r["imag"]) > 1 else [])
            if not extra:
                continue
            hess, modes = follow_modes(job)
            freqs = ", ".join("%.1f" % f for f in r["imag"])
            if not modes:
                notes.append("%s %s has imaginary frequencies (%s) but no readable orca.hess, so it cannot be pushed "
                             "along them: 'mechflow reset %s --stage %s --redo-done' restarts it from its last geometry"
                             % (job.item, job.stage, freqs, job.item, job.stage))
                continue
            acts[job] = ["displace", "%d imaginary frequency(ies): %s" % (len(r["imag"]), freqs)]
            continue
        note = blocked(job, s, named)
        if note is not None:
            if note:
                notes.append(note)
            continue
        if s == "failed":
            why = "failed: " + why_failed(run_folder(job))
            if job.stage == "SCAN" and not from_scratch and not (stages and "SCAN" in stages):
                if step_progress(project, project.steps[job.item], queue)[0] != "SCAN":
                    notes.append("%s SCAN stopped (%s) after passing a maximum, so it is not reset: 'mechflow next' "
                                 "goes on to the OptTS (reset it with --stage SCAN to run it further anyway)" % (
                                     job.item, why_failed(run_folder(job))))
                    continue
        elif s == "finished":
            if not redo_done:
                continue
            why = "finished"
        elif known:                                  # submitted / running, but no longer in the queue
            why = "%s, not in Slurm" % ("queued" if s == "submitted" else "running")
        else:                                        # submitted / running, state unknown, --force
            why = "%s (no squeue; --force)" % s
        how = "scratch" if (from_scratch or (job.stage == "SCAN" and s == "finished")) else "restart"
        acts[job] = [how, why]

    for job in list(acts):                           # the later stages that used the results of these jobs
        later = [d for d in downstream(project, job) if S(d) != "absent"]
        if not downstream_too:
            rest = [d for d in later if d not in acts]
            if rest:
                notes.append("%s %s: %d later stage(s) used its result (%s); add --downstream to reset them too" % (
                    job.item, job.stage, len(rest), ", ".join("%s %s" % (d.item, d.stage) for d in rest)))
            continue
        for d in later:
            if d in acts:
                acts[d][0] = "scratch"               # its start changes, so it starts again from the beginning
                continue
            note = blocked(d, S(d), True)
            if note is not None:
                if note:
                    notes.append(note + ", although it uses the result of %s %s" % (job.item, job.stage))
                continue
            acts[d] = ["scratch", "uses the result of %s %s" % (job.item, job.stage)]
    if unknown:
        notes.append("%d submitted job(s) were left alone because squeue could not be asked whether they are still "
                     "in Slurm; add --force to reset them too (once you know they are not running)" % len(unknown))
    order = {j: k for k, j in enumerate(
        [j for sp in project.species.values() for j in species_jobs(project, sp)] +
        [j for st in project.steps.values() for j in step_jobs(project, st)])}
    return [(j, how, why) for j, (how, why) in sorted(acts.items(), key=lambda kv: order.get(kv[0], 0))], notes


def describe(job, how, amplitude=AMPLITUDE):
    """What resetting the job will do (for --dry-run)."""
    if how == "scratch":
        return "from scratch ('mechflow next' writes it again)"
    if how == "displace":
        hess, modes = follow_modes(job)
        return "restarts from its last geometry %s" % displaced(job.folder, hess, modes, amplitude)[1]
    if job.stage == "SCAN":
        return ("a new SCAN part from its last converged point" if _scan_points_in(run_folder(job))
                else "runs again as it was (no converged SCAN point)")
    meta = read_meta(job.folder) or {}
    if _last_geometry(job.folder, meta):
        return "restarts from the last geometry of that attempt"
    return "runs again from the same structure"


# ============================================================================== resetting

def archive(folder, stamp, keep=(), copy=()):
    """Move everything in a job folder except `keep` into failed_attempts/<stamp>, and copy `copy` there too.
    Returns that folder (None when there was nothing to move)."""
    if not os.path.isdir(folder):
        return None
    names = [f for f in sorted(os.listdir(folder)) if f != ATTEMPTS and f not in keep]
    if not names:
        return None
    dest = os.path.join(folder, ATTEMPTS, stamp)
    k = 2
    while os.path.exists(dest):
        dest = os.path.join(folder, ATTEMPTS, "%s_%d" % (stamp, k))
        k += 1
    os.makedirs(dest)
    for f in copy:
        if os.path.isfile(os.path.join(folder, f)):
            shutil.copy2(os.path.join(folder, f), os.path.join(dest, f))
    for f in names:
        shutil.move(os.path.join(folder, f), os.path.join(dest, f))
    return dest


def _drop_4b_acceptance(project, job, dest, stamp):
    """'mechflow accept S5' was for the old Step 4B result: keep it with the attempt."""
    if job.stage not in ("Fwd", "Bwd"):
        return
    acc = os.path.join(project.step_dir(job.item), VALIDATE, ACCEPTED)
    if os.path.isfile(acc):
        if dest is None:
            dest = os.path.join(job.folder, ATTEMPTS, stamp)
            os.makedirs(dest, exist_ok=True)
        shutil.move(acc, os.path.join(dest, "Step4B" + ACCEPTED))


def _restart(project, job, folder, stamp, amplitude=None):
    """Keep the input, move the attempt's files away, restart an optimisation from its last geometry (amplitude: pushed
    along its imaginary modes), and write orca.inp / submit.sl again from the current settings."""
    meta = read_meta(folder)
    if not meta:
        raise ValueError("no %s in %s, so it cannot be restarted: use --from-scratch" % (META, folder))
    xyz = meta["xyz"]
    last = _last_geometry(folder, meta)
    geometry = pushed = None
    if amplitude is not None:
        hess, modes = follow_modes(job, folder)
        if not modes or not last:
            raise ValueError("no imaginary mode to follow (orca.hess, orca.xyz)")
        mol, pushed = displaced(folder, hess, modes, amplitude)
        geometry = "%d\nlast geometry of the attempt, %s\n%s" % (len(mol), pushed, "".join(
            "%-2s %16.10f %16.10f %16.10f\n" % ((s,) + tuple(c)) for s, c in zip(mol.symbols, mol.coords)))
    elif last:
        with open(last) as fh:
            geometry = fh.read()
    keep = {"orca.inp", "submit.sl", META, xyz, meta.get("first_xyz"), "TS_guess_from_SCAN.txt"}
    dest = archive(folder, stamp, keep={k for k in keep if k}, copy=("orca.inp", "submit.sl", xyz))
    _drop_4b_acceptance(project, job, dest, stamp)
    info = {k: v for k, v in meta.items() if k not in ("kind", "xyz", "charge", "mult", "jobname", "scan", "overrides")}
    info["restarts"] = info.get("restarts", 0) + 1
    if pushed:
        info["pushed"] = info.get("pushed", 0) + 1
    if geometry is not None:
        with open(os.path.join(folder, RESTART_XYZ), "w") as fh:
            fh.write(geometry)
        info.setdefault("first_xyz", xyz)
        xyz = RESTART_XYZ
    write_job(folder, project.settings, meta["kind"], xyz, meta["charge"], meta["mult"], meta["jobname"],
              tuple(meta["scan"]) if meta.get("scan") else None, meta.get("overrides"), info)
    if pushed:
        return dest, "restarts from its last geometry %s (%s)" % (pushed, RESTART_XYZ)
    return dest, ("restarts from the last geometry of that attempt (%s)" % RESTART_XYZ if geometry is not None
                  else "runs again from the same structure")


def apply(project, job, how, stamp, amplitude=AMPLITUDE):
    """Reset one job. Returns (folder holding the attempt's files or None, what was done)."""
    if how == "displace":
        return _restart(project, job, job.folder, stamp, amplitude=amplitude)
    if how == "scratch":
        dest = archive(job.folder, stamp)
        _drop_4b_acceptance(project, job, dest, stamp)
        return dest, "from scratch: 'mechflow next' writes it again"
    if job.stage == "SCAN":
        last = run_folder(job)
        if _scan_points_in(last):
            return None, restart_scan(project, project.steps[job.item])
        return _restart(project, job, last, stamp)
    return _restart(project, job, job.folder, stamp)
