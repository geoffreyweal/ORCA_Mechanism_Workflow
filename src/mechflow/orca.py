"""Writing ORCA inputs and Slurm scripts, reading ORCA outputs, and the state of a job folder.

The inputs have the same layout as the ORCA Mechanism Procedure templates; the output
markers were checked against ORCA 6.1 outputs from the procedure's examples.
"""
import glob
import json
import os
import re
import subprocess

from .defaults import deep_merge

META = ".mechflow.json"          # job parameters, so an input can be rewritten later ('next --refresh')
SUBMITTED = ".submitted"         # job id written when mechflow submits a job

# slurm.print_binding: each task prints its node, the node's NUMA domains and the cores it was given (the
# print_binding_stats step of NeSI's ORCA job scripts). A job is not stopped if this step fails.
PRINT_BINDING = (
    r'''srun --job-name=print_binding_stats bash -c "echo -e \"Task #\${SLURM_PROCID} is running on node '''
    r'''\$(hostname).\n\$(hostname) has the following NUMA configuration:\n\$(lscpu | grep -i --color=none numa)'''
    r'''\nTask #\${SLURM_PROCID} has \$(nproc) CPUs, their core IDs are \$(taskset -c -p \$\$ | awk '{print \$NF}')'''
    r'''\n===========================================\"" || echo "(the cores of each task could not be printed)"''')

# ----------------------------------------------------------------------------- settings for one job


def job_settings(settings, kind, overrides=None):
    """resources + slurm settings for a job of this kind (stage_overrides and per-item overrides applied)."""
    res = dict(settings["resources"])
    slurm = dict(settings["slurm"])
    for extra in ((settings.get("stage_overrides") or {}).get(kind) or {}, overrides or {}):
        for k, v in extra.items():
            if k in res:
                res[k] = v
            else:
                slurm[k] = v
    return res, slurm


def _block(text):
    return (text or "").strip("\n")


def render_input(settings, kind, xyzfile, charge, mult, scan=None, overrides=None):
    """Text of orca.inp for a job of kind 'opt' | 'scan' | 'optts' | 'irc' | 'sp'.

    'sp' is the Step 5 single point: the method and keywords of settings['single_point'] instead of orca.method
    and orca.keywords, with the same SCF and solvent blocks.
    scan = (i, j, first, last, npoints) with 0-based atom numbers, as ORCA expects."""
    o = settings["orca"]
    res, _ = job_settings(settings, kind, overrides)
    if kind == "sp":
        sp = settings.get("single_point") or {}
        out = [x.rstrip() for x in (sp.get("method") or "", sp.get("keywords") or "") if str(x).strip()]
    else:
        out = [o["method"].rstrip(), o["keywords"][kind].rstrip()]
    for blk in (o.get("scf_block"), o.get("solvent_block")):
        if _block(blk):
            out.append(_block(blk))
    out.append("%%PAL NPROCS %s END" % res["nprocs"])
    out.append("%%maxcore %s" % res["maxcore"])
    if kind == "scan":
        i, j, d0, d1, n = scan
        out.append("%GEOM")
        if o.get("scan_maxiter"):
            out.append("    MaxIter %s   # optimisation cycles allowed per SCAN point" % o["scan_maxiter"])
        out.append("    SCAN B %d %d = %.3f, %.3f, %d END" % (i, j, d0, d1, n))
        out.append("END")
    extra = (settings.get("single_point") or {}).get("block") if kind == "sp" else (o.get("blocks") or {}).get(kind)
    if _block(extra):
        out.append(_block(extra))
    out.append("* xyzfile %d %d %s" % (charge, mult, xyzfile))
    return "\n".join(out) + "\n\n\n"


def render_submit(settings, kind, jobname, overrides=None, local=False):
    res, s = job_settings(settings, kind, overrides)
    lines = ["#!/bin/bash -e",
             "#SBATCH --job-name=%s" % jobname[:64],
             "#SBATCH --ntasks=%s" % res["nprocs"],
             "#SBATCH --mem-per-cpu=%s" % s["mem_per_cpu"],
             "#SBATCH --partition=%s" % s["partition"],
             "#SBATCH --time=%s     # Walltime" % s["time"],
             "#SBATCH --output=slurm-%j.out      # %x and %j are replaced by job name and ID",
             "#SBATCH --error=slurm-%j.err",
             "#SBATCH --nodes=1 # OpenMPI can have problems with ORCA over multiple nodes sometimes, depending on your system."]
    lines += [x.strip() for x in (s.get("extra_sbatch") or []) if str(x).strip()]
    lines += [""]
    if s.get("module_purge"):
        lines += ["# Start from a clean environment", "module purge", ""]
    lines += ["# Load ORCA", "module load %s" % s["orca_module"]]
    lines += [str(x) for x in (s.get("extra_lines") or [])]
    if s.get("print_binding"):
        lines += ["", "# The cores each task was given, and the node's NUMA domains", PRINT_BINDING]
    lines += ["",
              "# ORCA under MPI requires that it be called via its full absolute path",
              "orca_exe=$(which orca)",
              "",
              "# Don't use \"srun\" as ORCA does that itself when launching its MPI process."]
    chain = settings.get("chain") or {}
    if not chain.get("enabled"):
        lines += ["${orca_exe} orca.inp > output.out", ""]
        return "\n".join(lines)
    lines += ["orca_status=0",
              "${orca_exe} orca.inp > output.out || orca_status=$?",
              "",
              "# Chained jobs (chain.enabled in project.yaml): if ORCA finished, write and submit every stage that",
              "# can start now. sbatch runs from inside this job, so the job's own SLURM_/SBATCH_ variables go first.",
              "if [ ${orca_status} -eq 0 ] && grep -q \"ORCA TERMINATED NORMALLY\" output.out; then",
              "    (",
              "        for v in $(compgen -e | grep -E '^(SLURM|SBATCH)_'); do unset \"$v\"; done"]
    lines += ["        " + str(x) for x in (chain.get("setup") or [])]
    lines += ["        mechflow chain",
              "    ) || echo \"mechflow chain did not finish: run 'mechflow next' and 'mechflow submit' in the project\"",
              "fi",
              "exit ${orca_status}",
              ""]
    return "\n".join(lines)


def write_job(folder, settings, kind, xyzfile, charge, mult, jobname, scan=None, overrides=None, info=None):
    """Write orca.inp + submit.sl (+ the job parameters in .mechflow.json)."""
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, "orca.inp"), "w") as fh:
        fh.write(render_input(settings, kind, xyzfile, charge, mult, scan, overrides))
    with open(os.path.join(folder, "submit.sl"), "w") as fh:
        fh.write(render_submit(settings, kind, jobname, overrides))
    meta = {"kind": kind, "xyz": xyzfile, "charge": charge, "mult": mult, "jobname": jobname,
            "scan": list(scan) if scan else None, "overrides": overrides or {}}
    meta.update(info or {})
    with open(os.path.join(folder, META), "w") as fh:
        json.dump(meta, fh, indent=1)


def read_meta(folder):
    try:
        with open(os.path.join(folder, META)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def rewrite_job(folder, settings):
    """Rewrite orca.inp / submit.sl of a prepared job from the current settings (True if done)."""
    m = read_meta(folder)
    if not m:
        return False
    write_job(folder, settings, m["kind"], m["xyz"], m["charge"], m["mult"], m["jobname"],
              tuple(m["scan"]) if m.get("scan") else None, m.get("overrides"),
              {k: v for k, v in m.items() if k not in ("kind", "xyz", "charge", "mult", "jobname", "scan", "overrides")})
    return True


# ----------------------------------------------------------------------------- reading outputs

def read_text(path):
    try:
        with open(path, errors="replace") as fh:
            return fh.read()
    except OSError:
        return ""


_FREQ = re.compile(r"\s*(\d+):\s+(-?\d+\.\d+)\s+cm\*\*-1(.*)")


def parse_output(path):
    """What the procedure checks in an ORCA output file."""
    txt = read_text(path)
    r = {"exists": bool(txt)}
    r["terminated_normally"] = "ORCA TERMINATED NORMALLY" in txt
    r["error"] = bool(re.search(r"ORCA finished by error termination|aborting the run|ERROR !!!", txt))
    r["opt_converged"] = "THE OPTIMIZATION HAS CONVERGED" in txt
    r["irc_converged"] = txt.count("THE IRC HAS CONVERGED")
    g = re.findall(r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)", txt)
    r["gibbs"] = float(g[-1]) if g else None
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", txt)
    r["e_el"] = float(e[-1]) if e else None
    # thermal correction of a frequency job: G minus the electronic energy (for G = E(single point) + this)
    c = re.findall(r"G-E\(el\)\s+\.\.\.\s+(-?\d+\.\d+)", txt)
    r["g_minus_eel"] = float(c[-1]) if c else None
    freqs, imag = [], []
    idx = txt.rfind("VIBRATIONAL FREQUENCIES")
    if idx >= 0:
        for ln in txt[idx:].splitlines()[1:2000]:
            m = _FREQ.match(ln)
            if m:
                f = float(m.group(2))
                freqs.append(f)
                if "imaginary" in m.group(3) or f < -1.0:
                    imag.append(f)
            elif freqs and not ln.strip():
                break
    r["freqs"] = freqs
    r["imag"] = imag
    r["scan_steps_started"] = len(re.findall(r"RELAXED SURFACE SCAN STEP\s+\d+", txt))
    return r


SLURM_MARKERS = ("DUE TO TIME LIMIT", "CANCELLED", "oom-kill", "OUT_OF_MEMORY", "Out Of Memory", "out-of-memory",
                 "NODE_FAIL", "DUE TO NODE FAILURE")


def slurm_problem(folder):
    for f in sorted(glob.glob(os.path.join(folder, "slurm-*.err")) + glob.glob(os.path.join(folder, "slurm-*.out"))):
        t = read_text(f)
        for marker in SLURM_MARKERS:
            if marker in t:
                return marker
    return None


# ----------------------------------------------------------------------------- job state

class Queue:
    """Job ids in the Slurm queue (asked once, with a single squeue call)."""

    def __init__(self, enabled=True):
        self.enabled = enabled
        self._ids = "unset"

    def ids(self):
        if self._ids == "unset":
            self._ids = None
            if self.enabled:
                try:
                    res = subprocess.run(["squeue", "-h", "-u", os.environ.get("USER", ""), "-o", "%i"],
                                         capture_output=True, text=True, timeout=60)
                    if res.returncode == 0:
                        self._ids = set(res.stdout.split())
                except (OSError, subprocess.SubprocessError):
                    self._ids = None
        return self._ids


def job_ids(folder):
    ids = set()
    sub = os.path.join(folder, SUBMITTED)
    if os.path.isfile(sub):
        first = read_text(sub).split()[:1]
        if first and first[0].isdigit():
            ids.add(first[0])
    for f in glob.glob(os.path.join(folder, "slurm-*.out")):
        m = re.search(r"slurm-(\d+)\.out$", f)
        if m:
            ids.add(m.group(1))
    return ids


def job_state(folder, queue=None):
    """absent | ready | submitted | running | finished | failed"""
    if not os.path.isfile(os.path.join(folder, "orca.inp")):
        return "absent"
    out = os.path.join(folder, "output.out")
    ids = job_ids(folder)
    act = queue.ids() if (queue is not None and ids) else None
    gone = act is not None and bool(ids) and not (ids & act)
    if not os.path.isfile(out):
        if not os.path.exists(os.path.join(folder, SUBMITTED)) and not ids:
            return "ready"
        return "failed" if gone else "submitted"
    r = parse_output(out)
    if r["terminated_normally"]:
        return "finished"
    if r["error"] or slurm_problem(folder) or gone:
        return "failed"
    return "running"
