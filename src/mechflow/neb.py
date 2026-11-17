"""Step 2 by NEB-CI (the procedure's Step 2Bb) and a chemist in the loop (approvals).

With neb.enabled, Step 2 of a step is an NEB-CI between its two ends instead of a SCAN:

    Step2Bb_Ends/Left, Right   optimisation of an end of more than one molecule, or of an end file whose comment line
                               says 'optimise first' (an end of one molecule is otherwise used as it is: the optimised
                               structure of that molecule); if its bonds change, the step waits
    Step2Bb_NEB                the NEB-CI between the two ends; its climbing image is the TS guess for the OptTS

The ends are xyz files given in the project (neb.ends, e.g. neb/S5a_left.xyz and neb/S5a_right.xyz): the two sides
of the step's reaction, with the same atoms in the same order.

With approval.enabled, Step 2 (SCAN or NEB) waits until the start is approved and the OptTS waits until the TS guess
is approved: 'mechflow approve <steps> --start' or '--guess' writes approved_start.txt / approved_guess.txt in the
step's folder (or list the steps under approval.start / approval.guess in project.yaml). An approved NEB start
holds only while the end files stay as they were when it was approved.
"""
import getpass
import hashlib
import os
import re
import shutil
import time

from . import chem
from .orca import job_state, parse_output, read_text, slurm_problem, write_job

NEB_ENDS = "Step2Bb_Ends"
NEB = "Step2Bb_NEB"
SIDES = (("left", "Left"), ("right", "Right"))
LEFT_END, RIGHT_END = "left_end.xyz", "right_end.xyz"
CI_FILE = "orca_NEB-CI_converged.xyz"
MEP_FILE = "orca_MEP_trj.xyz"
ALLXYZ = "orca_MEP.allxyz"
RESTART_ALLXYZ = "restart_path.allxyz"
TS_GUESS = "ts_guess.xyz"
APPROVALS = ("start", "guess")
ACCEPTED = ".accepted"
KCAL_KJ = 4.184


# ------------------------------------------------------------------ settings
def enabled(project):
    return bool((project.settings.get("neb") or {}).get("enabled"))


def images(project, st=None):
    """NImages: the step's own neb_images, else neb.images (ORCA's default 8)."""
    own = (st.overrides or {}).get("neb_images") if st is not None else None
    return int(own or (project.settings.get("neb") or {}).get("images") or 8)


def end_file(project, st, side):
    pattern = str((project.settings.get("neb") or {}).get("ends") or "neb/{step}_{side}.xyz")
    return project.path(pattern.format(step=st.id, side=side))


def ts_file(project, st):
    """The step's TS guess for its NEB (neb.ts, e.g. "neb/{step}_ts.xyz": the same atoms in the same order as the two
    ends): ORCA interpolates the starting path from the left end to it and on to the right end. None without the
    setting or the file."""
    pattern = (project.settings.get("neb") or {}).get("ts")
    if not pattern:
        return None
    f = project.path(str(pattern).format(step=st.id))
    return f if os.path.isfile(f) else None


def approval_on(project):
    return bool((project.settings.get("approval") or {}).get("enabled"))


# ------------------------------------------------------------------ approvals
def approval_file(project, st, what):
    return os.path.join(project.step_dir(st.id), "approved_%s.txt" % what)


def ends_digest(project, st):
    """A checksum of the two end files (an approved NEB start holds only while they stay the same)."""
    h = hashlib.sha1()
    for side, _ in SIDES:
        try:
            with open(end_file(project, st, side), "rb") as fh:
                h.update(fh.read())
        except OSError:
            h.update(b"(missing)")
    ts = ts_file(project, st)
    if ts:                                    # (steps without a TS guess keep the checksum they had before neb.ts)
        with open(ts, "rb") as fh:
            h.update(b"(ts)" + fh.read())
    return h.hexdigest()


def approval_state(project, st, what):
    """(approved?, note). Always approved when approval.enabled is false."""
    if not approval_on(project):
        return True, ""
    listed = [str(x) for x in ((project.settings.get("approval") or {}).get(what) or [])]
    if st.id in listed:
        return True, "approval.%s in project.yaml" % what
    f = approval_file(project, st, what)
    if not os.path.isfile(f):
        return False, ""
    text = read_text(f)
    if what == "start" and enabled(project):
        m = re.search(r"^ends: ([0-9a-f]+)", text, re.M)
        if m and m.group(1) != ends_digest(project, st):
            return False, "the end files changed after you approved them"
    return True, "approved (%s)" % os.path.basename(f)


def is_approved(project, st, what):
    return approval_state(project, st, what)[0]


def approve(project, st, what, note=""):
    os.makedirs(project.step_dir(st.id), exist_ok=True)
    lines = ["%s of %s approved by %s on %s" % ("the start" if what == "start" else "the TS guess", st.id,
                                                  getpass.getuser(), time.strftime("%Y-%m-%d %H:%M"))]
    if note:
        lines.append("note: %s" % note)
    if what == "start" and enabled(project):
        lines.append("ends: %s" % ends_digest(project, st))
    with open(approval_file(project, st, what), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def unapprove(project, st, what):
    f = approval_file(project, st, what)
    if os.path.isfile(f):
        os.remove(f)
        return True
    return False


# ------------------------------------------------------------------ the NEB output
def parse_neb(path):
    """What 'mechflow' needs from an NEB output: converged?, the path summary (image, distance, E in Eh, dE in
    kcal/mol, climbing image?), the climbing image."""
    txt = read_text(path)
    r = {"exists": bool(txt), "converged": "THE NEB OPTIMIZATION HAS CONVERGED" in txt,
         "terminated_normally": "ORCA TERMINATED NORMALLY" in txt, "path": [], "ci": None}
    k = txt.rfind("PATH SUMMARY\n")
    if k < 0:
        k = txt.rfind("PATH SUMMARY")
    if k >= 0:
        for ln in txt[k:].splitlines()[1:200]:
            m = re.match(r"\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\S+)\s+(\S+)(.*)", ln)
            if m:
                ci = "CI" in m.group(7)
                r["path"].append({"image": int(m.group(1)), "dist": float(m.group(2)), "E": float(m.group(3)),
                                  "dE_kcal": float(m.group(4)), "ci": ci})
                if ci:
                    r["ci"] = int(m.group(1))
            elif r["path"] and not ln.strip():
                break
    return r


def neb_summary(project, st):
    """A one-line account of a finished NEB (for 'mechflow status')."""
    r = parse_neb(os.path.join(project.step_dir(st.id), NEB, "output.out"))
    if not r["path"]:
        return "NEB converged"
    top = max(r["path"], key=lambda p: p["dE_kcal"])
    ci = [p for p in r["path"] if p["ci"]]
    pick = ci[0] if ci else top
    return "NEB converged: %s image %d of %d, %.1f kJ/mol above the left end (%s)" % (
        "climbing" if ci else "highest", pick["image"], len(r["path"]) - 1, pick["dE_kcal"] * KCAL_KJ, MEP_FILE)


# ------------------------------------------------------------------ the ends
def _molecules(mol, tol):
    return len(chem.fragments(len(mol), chem.bonds(mol, tol)))


def wants_optimisation(mol):
    """An end file whose comment line asks for it ('optimise first') is optimised before the NEB even when it is one
    molecule: e.g. a conformer set up for the reaction rather than the molecule's own Step 1 structure."""
    return "optimise first" in (mol.comment or "").lower()


def end_folder(project, st, label):
    return os.path.join(project.step_dir(st.id), NEB_ENDS, label)


def end_state(project, st, side, label, queue=None):
    """One end of the step: {file, mol, opt (needs an optimisation), state, ok, xyz (the structure the NEB uses),
    note}."""
    from .workflow import bond_change_text
    out = {"file": end_file(project, st, side), "mol": None, "opt": False, "state": None, "ok": False, "xyz": None,
           "note": ""}
    if not os.path.isfile(out["file"]):
        out["note"] = "no NEB end file %s" % os.path.relpath(out["file"], project.root)
        return out
    try:
        out["mol"] = chem.read_xyz(out["file"])
    except (OSError, chem.XYZError) as exc:
        out["note"] = "cannot read %s: %s" % (os.path.relpath(out["file"], project.root), exc)
        return out
    if _molecules(out["mol"], project.tolerance) <= 1 and not wants_optimisation(out["mol"]):
        out.update(ok=True, xyz=out["file"])
        return out
    out["opt"] = True
    folder = end_folder(project, st, label)
    out["state"] = s = job_state(folder, queue)
    if s != "finished":
        return out
    r = parse_output(os.path.join(folder, "output.out"))
    xyz = os.path.join(folder, "orca.xyz")
    if not r["opt_converged"] or not os.path.isfile(xyz):
        out["state"] = "failed"
        out["note"] = "the %s end's optimisation did not converge (%s/%s/output.out)" % (side, NEB_ENDS, label)
        return out
    change = bond_change_text(project, out["mol"], chem.read_xyz(xyz))
    if change and not os.path.isfile(os.path.join(project.step_dir(st.id), NEB_ENDS, ACCEPTED)):
        out["note"] = "the %s end changed during its optimisation (%s): look at %s/%s/orca.xyz" % (
            side, change, NEB_ENDS, label)
        out["changed"] = True
        return out
    out.update(ok=True, xyz=xyz)
    return out


def _check_ends(project, st, ends, queue):
    """A problem with the two end files (or None): the same atoms in the same order, each the molecules of its side."""
    from .workflow import _side_fingerprints, project_expression, species_result
    left, right = ends["left"]["mol"], ends["right"]["mol"]
    if left.symbols != right.symbols:
        return "CHECK", "the two NEB ends do not have the same atoms in the same order (%d and %d atoms)" % (
            len(left), len(right))
    ts = ts_file(project, st)
    if ts and chem.read_xyz(ts).symbols != left.symbols:
        return "CHECK", "the TS guess %s does not have the atoms of the NEB ends in the same order" % (
            os.path.relpath(ts, project.root))
    waiting = [s for s in sorted(set(st.left) | set(st.right))
               if not species_result(project, project.species[s], queue)["ok"]]
    if waiting:
        return "waiting", "for Step 1 of %s (the NEB ends are made of the optimised molecules)" % ", ".join(waiting)
    for side, counter in (("left", st.left), ("right", st.right)):
        want = _side_fingerprints(project, counter, queue)
        got = sorted(chem.fragment_fingerprints(ends[side]["mol"], project.tolerance))
        if got != want:
            return "CHECK", "the %s NEB end (%s) is not %s" % (
                side, " + ".join(f for _, f in got), project_expression(counter))
    return None


def progress(project, st, queue=None):
    """(stage, detail, action) of a step's Step 2 by NEB, up to the TS guess ('prepare_optts' once it is there)."""
    nd = os.path.join(project.step_dir(st.id), NEB)
    sn = job_state(nd, queue)
    if sn == "absent":
        ends = {side: end_state(project, st, side, label, queue) for side, label in SIDES}
        missing = [e["note"] for e in ends.values() if e["mol"] is None]
        if missing:
            return "waiting", "; ".join(missing) + " (neb.ends in project.yaml)", None
        problem = _check_ends(project, st, ends, queue)
        if problem:
            return problem[0], problem[1], None
        ok, why = approval_state(project, st, "start")
        if not ok:
            files = " and ".join(os.path.relpath(ends[s]["file"], project.root) for s, _ in SIDES)
            ts = ts_file(project, st)
            return "approve", "%sNEB ends ready (%s%s): look at them, then 'mechflow approve %s --start'" % (
                (why + "; ") if why else "", files, (", TS guess " + os.path.relpath(ts, project.root)) if ts else "",
                st.id), None
        to_opt = [s for s, _ in SIDES if ends[s]["opt"]]
        if any(ends[s]["state"] == "absent" for s in to_opt):
            why = lambda s_: "marked 'optimise first'" if wants_optimisation(ends[s_]["mol"]) else "more than one molecule"  # noqa: E731
            return "NEB ends", "optimised first: %s" % ", ".join("%s (%s)" % (s_, why(s_)) for s_ in to_opt), "prepare_neb_ends"
        active = [s for s in to_opt if ends[s]["state"] in ("ready", "submitted", "running")]
        if active:
            return "NEB ends", ", ".join("%s %s" % (s, ends[s]["state"]) for s in to_opt), None
        for s in to_opt:
            if ends[s].get("changed"):
                return "CHECK", ends[s]["note"] + "; if it is right: 'mechflow accept %s'" % st.id, None
            if not ends[s]["ok"]:
                return "NEB ends", (ends[s]["note"] or "the %s end failed" % s) + \
                    "; 'mechflow reset %s --stage Ends' runs it on" % st.id, None
        how = ", ".join("%s %s" % (s, "optimised" if ends[s]["opt"] else "as given") for s, _ in SIDES)
        return "NEB", "ends ready (%s)" % how, "prepare_neb"
    if sn in ("ready", "submitted", "running"):
        return "NEB", sn, None
    r = parse_neb(os.path.join(nd, "output.out"))
    if not r["converged"]:
        return "NEB", "did not converge (%s) - see %s/output.out; 'mechflow reset %s --stage NEB' runs it on from its " \
                      "last path" % (slurm_problem(nd) or "ORCA", NEB, st.id), None
    if not os.path.isfile(os.path.join(nd, CI_FILE)):
        return "NEB", "converged without a climbing image (no maximum along the path: barrierless?) - look at " \
                      "%s/%s" % (NEB, MEP_FILE), None
    return "OptTS", neb_summary(project, st), "prepare_optts"


# ------------------------------------------------------------------ writing the jobs
def prepare(project, st, action, image=None):
    from .workflow import OPTTS, _write_step_info, jobname
    sd = project.step_dir(st.id)
    if action == "prepare_neb_ends":
        written = []
        for side, label in SIDES:
            e = end_state(project, st, side, label)
            if not e["opt"] or e["state"] != "absent":
                continue
            folder = end_folder(project, st, label)
            os.makedirs(folder, exist_ok=True)
            name = "%s_end.xyz" % side
            chem.write_xyz(os.path.join(folder, name), e["mol"], "%s %s end (from %s)  charge=%d mult=%d" % (
                st.id, side, os.path.relpath(e["file"], project.root), st.charge, st.mult))
            write_job(folder, project.settings, "nebend", name, st.charge, st.mult, jobname(project, st.id, "End" + label),
                      overrides=st.overrides, info={"item": st.id, "stage": "Ends", "side": side})
            written.append(side)
        _write_step_info(project, st)
        return "optimisation of the %s NEB end%s (more than one molecule)" % (" and ".join(written),
                                                                             "s" if len(written) > 1 else "")
    if action == "prepare_neb":
        folder = os.path.join(sd, NEB)
        os.makedirs(folder, exist_ok=True)
        lines = []
        for (side, label), name in zip(SIDES, (LEFT_END, RIGHT_END)):
            e = end_state(project, st, side, label)
            if not e["ok"]:
                raise ValueError("the %s NEB end is not ready: %s" % (side, e["note"] or e["state"]))
            mol = chem.read_xyz(e["xyz"])
            chem.write_xyz(os.path.join(folder, name), mol, "%s %s end  charge=%d mult=%d" % (st.id, side, st.charge, st.mult))
            lines.append("%s: %s" % (name, os.path.relpath(e["xyz"], project.root)))
        info = {"item": st.id, "stage": "NEB", "neb_end": RIGHT_END}
        ts = ts_file(project, st)
        if ts:                                # the starting path goes through the TS guess (NEB_TS_XYZFILE)
            chem.write_xyz(os.path.join(folder, TS_GUESS), chem.read_xyz(ts), "%s TS guess  charge=%d mult=%d" % (
                st.id, st.charge, st.mult))
            lines.append("%s: %s (the starting path goes through it)" % (TS_GUESS, os.path.relpath(ts, project.root)))
            info["neb_ts"] = TS_GUESS
        with open(os.path.join(folder, "NEB_ends.txt"), "w") as fh:
            fh.write("The NEB-CI of %s runs from %s (left side of the reaction) to %s (right side):\n%s\n" % (
                st.id, LEFT_END, RIGHT_END, "\n".join("  " + x for x in lines)))
        write_job(folder, project.settings, "neb", LEFT_END, st.charge, st.mult, jobname(project, st.id, "NEB"),
                  overrides=st.overrides, info=info)
        _write_step_info(project, st)
        return "NEB-CI with %d images between the two ends%s" % (images(project, st), ", through the TS guess" if ts else "")
    if action == "prepare_optts":
        nd = os.path.join(sd, NEB)
        folder = os.path.join(sd, OPTTS)
        if image is not None:
            frames = chem.read_xyz_frames(os.path.join(nd, MEP_FILE))
            if not 0 <= image < len(frames):
                raise ValueError("the NEB path has images 0-%d" % (len(frames) - 1))
            mol, what = frames[image], "image %d of the NEB path (%s) - chosen with --image" % (image, MEP_FILE)
        else:
            mol, what = chem.read_xyz(os.path.join(nd, CI_FILE)), "the climbing image of the NEB (%s)" % CI_FILE
        os.makedirs(folder, exist_ok=True)
        name = "NEB_TS_guess.xyz"
        chem.write_xyz(os.path.join(folder, name), mol, "%s TS guess from the NEB  charge=%d mult=%d" % (
            st.id, st.charge, st.mult))
        with open(os.path.join(folder, "TS_guess_from_NEB.txt"), "w") as fh:
            fh.write("TS guess: %s\n%s\n" % (what, neb_summary(project, st)))
        write_job(folder, project.settings, "optts", name, st.charge, st.mult, jobname(project, st.id, "OptTS"),
                  overrides=st.overrides, info={"item": st.id, "stage": "OptTS"})
        return "OptTS from %s" % what
    raise ValueError("unknown action %s" % action)


def restart_from_path(folder):
    """Before an NEB runs again: its last path (orca_MEP.allxyz) as the start of the next attempt, if there is one.
    Returns the file name to give as Restart_ALLXYZFile, or None."""
    src = os.path.join(folder, ALLXYZ)
    if os.path.isfile(src):
        shutil.copy(src, os.path.join(folder, RESTART_ALLXYZ))
        return RESTART_ALLXYZ
    return None
