"""The ORCA Mechanism Procedure, stage by stage, for every species and step of a project.

    Step 1    optimise every species                       calcs/Species/<ID>/Step1_Geo_Opt
    Step 2Ba  relaxed SCAN of one distance                 calcs/Steps/<ID>/Step2Ba_SCAN
    Step 3    OptTS from the highest SCAN point            calcs/Steps/<ID>/Step3_Opt_TS
    Step 4A   IRC from the transition state                calcs/Steps/<ID>/Step4A_IRC
    Step 4B   optimise both IRC end points                 calcs/Steps/<ID>/Step4B_Validate/{Forwards,Backwards}
    SP        single points (optional: single_point.enabled), e.g. DLPNO-CCSD(T), on each optimised species of one
              molecule and each transition state   calcs/Species/<ID>/SP_Single_Point, calcs/Steps/<ID>/SP_Single_Point
    Step 5    energies                                     mechflow collect
"""
import glob
import math
import os
import re
import shutil

from . import chem
from .orca import (job_state, parse_output, read_meta, read_text, rewrite_job, slurm_problem, write_job)

HARTREE_KJ = 2625.4996
HARTREE_KCAL = 627.5095

STEP1 = "Step1_Geo_Opt"
SCAN = "Step2Ba_SCAN"
OPTTS = "Step3_Opt_TS"
IRC = "Step4A_IRC"
VALIDATE = "Step4B_Validate"
SP = "SP_Single_Point"
ACCEPTED = ".accepted"       # written by 'mechflow accept'
ATTEMPTS = "failed_attempts"   # 'mechflow reset' keeps the files of earlier attempts of a job here


def jobname(project, item, stage):
    return "%s_%s_%s" % (project.name, item, stage)


def is_accepted(folder):
    return os.path.isfile(os.path.join(folder, ACCEPTED))


def _same_geometry(a, b):
    try:
        ma, mb = chem.read_xyz(a), chem.read_xyz(b)
    except (OSError, chem.XYZError):
        return False
    if ma.symbols != mb.symbols:
        return False
    return all(abs(p - q) < 1e-6 for ca, cb in zip(ma.coords, mb.coords) for p, q in zip(ca, cb))


def bond_change_text(project, before, after):
    """'lost C1-O4; new O4-H22' between two structures with the same atom order."""
    b0 = chem.bonds(before, project.tolerance)
    b1 = chem.bonds(after, project.tolerance)
    lab = lambda m, ij: "%s-%s" % (project.atom_label(m, ij[0]), project.atom_label(m, ij[1]))
    parts = []
    if b0 - b1:
        parts.append("lost " + ", ".join(lab(before, ij) for ij in sorted(b0 - b1)))
    if b1 - b0:
        parts.append("new " + ", ".join(lab(after, ij) for ij in sorted(b1 - b0)))
    return "; ".join(parts)


# ============================================================================== species (Step 1)

def species_setup(project, sp, queue=None):
    """Write the Step 1 job of a species when it is missing (or refresh it before it is submitted).

    Returns (message or None, is_warning)."""
    if assembled(project, sp):
        return None, False             # no Step 1: its SCANs start from its molecules' optimised structures
    d = project.species_job_dir(sp.id)
    start = os.path.join(d, "start_structure.xyz")
    state = job_state(d, queue)
    info = {"item": sp.id, "stage": "Step1"}
    if state == "absent":
        os.makedirs(d, exist_ok=True)
        shutil.copy(sp.path, start)
        write_job(d, project.settings, "opt", "start_structure.xyz", sp.charge, sp.mult,
                  jobname(project, sp.id, "Step1"), overrides=sp.overrides, info=info)
        return "Step 1 input written", False
    meta = read_meta(d) or {}
    changed = not _same_geometry(sp.path, start)
    q_changed = (meta.get("charge"), meta.get("mult")) != (sp.charge, sp.mult)
    if state == "ready" and (changed or q_changed):
        shutil.copy(sp.path, start)
        write_job(d, project.settings, "opt", "start_structure.xyz", sp.charge, sp.mult,
                  jobname(project, sp.id, "Step1"), overrides=sp.overrides, info=info)
        return "Step 1 input updated (the species xyz, charge or multiplicity changed)", False
    if state != "ready" and (changed or q_changed):
        return ("species/%s.xyz (or its charge / multiplicity) changed after its Step 1 job was submitted: "
                "'mechflow reset %s --redo-done --from-scratch --downstream' starts that species again from it"
                % (sp.id, sp.id)), True
    return None, False


def species_result(project, sp, queue=None):
    """State of a species' Step 1 optimisation and whether its result can be used."""
    d = project.species_job_dir(sp.id)
    state = job_state(d, queue)
    res = {"state": state, "ok": False, "gibbs": None, "e_el": None, "g_minus_eel": None, "note": "", "xyz": None,
           "imag": []}
    if state == "finished":
        r = parse_output(os.path.join(d, "output.out"))
        res["gibbs"], res["e_el"], res["imag"] = r["gibbs"], r["e_el"], r["imag"]
        res["g_minus_eel"] = r["g_minus_eel"]
        xyz = os.path.join(d, "orca.xyz")
        if not r["opt_converged"]:
            res["note"] = "optimisation did not converge"
        elif r["imag"]:
            res["note"] = "%d imaginary frequency(ies): %s (%s)" % (
                len(r["imag"]), ", ".join("%.1f" % f for f in r["imag"]), imaginary_hint(project, sp.id, d))
        elif r["gibbs"] is None:
            res["note"] = "no Gibbs free energy in output.out"
        elif not os.path.isfile(xyz):
            res["note"] = "orca.xyz missing"
        else:
            res["xyz"] = xyz
            res["ok"] = True
            try:
                change = bond_change_text(project, chem.read_xyz(os.path.join(d, "start_structure.xyz")),
                                          chem.read_xyz(xyz))
            except (OSError, chem.XYZError):
                change = ""
            if change and not is_accepted(d):
                res["ok"] = False
                res["note"] = ("bonds changed during the optimisation (%s): check it, then fix species/%s.xyz "
                               "or run 'mechflow accept %s'" % (change, sp.id, sp.id))
            elif change:
                res["note"] = "bonds changed during the optimisation (%s) - accepted" % change
    elif state == "failed":
        res["note"] = slurm_problem(d) or "ORCA stopped with an error - see output.out"
    return res


# ============================================================================== SP: single points

def sp_enabled(project):
    return bool((project.settings.get("single_point") or {}).get("enabled"))


def sp_auto(project):
    """single_point.auto: 'mechflow next' writes the single points (and 'submit' submits them) by itself. Otherwise only
    'mechflow sp' runs them, when you say so."""
    return sp_enabled(project) and bool((project.settings.get("single_point") or {}).get("auto"))


def species_sp_dir(project, sid):
    return project.path("calcs", "Species", sid, SP)


def step_sp_dir(project, sid):
    return os.path.join(project.step_dir(sid), SP)


def ts_result(project, st):
    """The Step 3 OptTS of a step when it is a usable TS (converged, exactly one imaginary frequency), else None."""
    t3 = os.path.join(project.step_dir(st.id), OPTTS)
    r = parse_output(os.path.join(t3, "output.out"))
    if r["terminated_normally"] and r["opt_converged"] and len(r["imag"]) == 1 and os.path.isfile(os.path.join(t3, "orca.xyz")):
        r["xyz"] = os.path.join(t3, "orca.xyz")
        return r
    return None


def molecules(project, sp):
    """How many molecules a species is (its pieces that are not bonded to each other)."""
    if sp.mol is None:
        return 1
    return len(chem.fragments(len(sp.mol), chem.bonds(sp.mol, project.tolerance)))


def is_state(project, sp):
    """Whether a mechanism has the species in one of its states (other than as the start of a TS's SCAN)."""
    return any(sp.id in s.full for m in project.mechanisms.values() for s in m.states if s.kind != "ts")


def sp_needed(project, sp):
    """Whether a species gets a single point. A species of more than one molecule (a complex, such as the start of a
    SCAN) needs none: a state's energy is the sum of its molecules' energies, and a TS's state counts the TS in place of
    its start. The exceptions: a mechanism has the complex itself as a state, or single_point.complexes is true. A
    complex built from its molecules (complexes.optimise: false) has no optimised structure, so never gets one."""
    if assembled(project, sp):
        return False
    if (project.settings.get("single_point") or {}).get("complexes"):
        return True
    if molecules(project, sp) <= 1:
        return True
    return is_state(project, sp)


# ============================================================================== complexes built from their molecules

def imaginary_settings(project):
    """(auto, max_pushes, amplitude) of the imaginary section"""
    s = project.settings.get("imaginary") or {}
    return bool(s.get("auto")), int(s.get("max_pushes", 2)), float(s.get("amplitude", 0.3))


def pushes(folder):
    """how many times a job was restarted pushed along its imaginary modes"""
    return int((read_meta(folder) or {}).get("pushed", 0))


def imaginary_hint(project, item, folder):
    """what to do about an optimisation with imaginary frequencies it should not have"""
    auto, most, _ = imaginary_settings(project)
    done = pushes(folder)
    if auto and done < most:
        return "'mechflow next' redoes it pushed along the mode (push %d of %d)" % (done + 1, most)
    if done:
        return "pushed %d time(s) already; 'mechflow reset %s --imaginary' pushes it again" % (done, item)
    return "'mechflow reset %s --imaginary' redoes it, pushed along the mode" % item


def optimise_complexes(project):
    """complexes.optimise (default true): give complexes a Step 1 optimisation of their own."""
    return bool((project.settings.get("complexes") or {}).get("optimise", True))


def assembled(project, sp):
    """A complex that gets no Step 1 of its own (complexes.optimise: false): a species of more than one molecule that
    no mechanism has as a state (e.g. the start of a SCAN, such as a substrate hydrogen-bonded to an acid). Its SCANs
    start from its molecules' own optimised structures, put where the complex's xyz file has them (assemble)."""
    if optimise_complexes(project) or sp.mol is None:
        return False
    cache = project.__dict__.setdefault("_assembled", {})
    if sp.id not in cache:
        cache[sp.id] = molecules(project, sp) > 1 and not is_state(project, sp)
    return cache[sp.id]


def _atom_range(idx, base):
    runs, start = [], None
    for k, a in enumerate(idx):
        if start is None:
            start = a
        if k + 1 == len(idx) or idx[k + 1] != a + 1:
            runs.append("%d" % (start + base) if start == a else "%d-%d" % (start + base, a + base))
            start = None
    return ",".join(runs)


def complex_parts(project, sp):
    """[(species ID, atom numbers in the complex)] for each molecule of a complex: the species of one molecule with the
    same atoms and bonds, listed in the same order in its own xyz file. If several species fit (conformers or
    stereoisomers), the one closest in shape is taken. ValueError when a molecule fits none."""
    cache = project.__dict__.setdefault("_complex_parts", {})
    if sp.id in cache:
        if isinstance(cache[sp.id], Exception):
            raise cache[sp.id]
        return cache[sp.id]
    try:
        out = _complex_parts(project, sp)
    except ValueError as exc:
        cache[sp.id] = exc
        raise
    cache[sp.id] = out
    return out


def _complex_parts(project, sp):
    singles = project.__dict__.get("_single_molecules")
    if singles is None:
        singles = {}
        for other in project.species.values():
            if other.mol is not None and molecules(project, other) == 1:
                singles.setdefault(chem.molecule_fingerprint(other.mol, project.tolerance), []).append(other)
        project.__dict__["_single_molecules"] = singles
    bonds = chem.bonds(sp.mol, project.tolerance)
    out = []
    for frag in chem.fragments(len(sp.mol), bonds):
        idx = sorted(frag)
        sub = chem.Molecule([sp.mol.symbols[a] for a in idx], [sp.mol.coords[a] for a in idx])
        what = "the molecule %s of %s (atoms %s)" % (chem.formula_string(chem.formula_of(sub.symbols)), sp.id,
                                                     _atom_range(idx, project.index_base))
        cands = [c for c in singles.get(chem.molecule_fingerprint(sub, project.tolerance), []) if c.id != sp.id]
        if not cands:
            raise ValueError("%s is not a species of the project, so it cannot be optimised on its own: add it as a "
                             "species (an xyz file of its own), or set complexes.optimise: true" % what)
        pos = {a: k for k, a in enumerate(idx)}
        own = set((pos[i], pos[j]) for i, j in bonds if i in pos and j in pos)
        fits = [c for c in cands if c.mol.symbols == sub.symbols and chem.bonds(c.mol, project.tolerance) == own]
        if not fits:
            raise ValueError("%s is %s, but its atoms are listed in another order: list them as in species/%s.xyz" % (
                what, " or ".join(c.id for c in cands), cands[0].id))
        if len(fits) > 1:
            def shape(c):
                rot, t = chem.superpose(c.mol.coords, sub.coords)
                return chem.rmsd(chem.transform(c.mol.coords, rot, t), sub.coords)
            fits.sort(key=shape)
        out.append((fits[0].id, idx))
    return out


CONTACT_SCALE = 3.0      # Angstrom: atoms this much further from the other molecules count e (2.7) times less in the fit


def assemble(project, sp, queue=None):
    """The structure of a complex built from its molecules: each molecule's optimised structure (Step 1) is fitted onto
    the place the complex's xyz file has that molecule, weighted towards its atoms near the other molecules, so that the
    contacts drawn in the complex (e.g. the hydrogen bond a SCAN starts from) stay as they were.

    Returns (Molecule, notes) - or (None, why) while a molecule has no usable optimised structure yet."""
    try:
        parts = complex_parts(project, sp)
    except ValueError as exc:
        return None, "cannot be built from its molecules: %s" % exc
    waiting = []
    for pid, _ in parts:
        r = species_result(project, project.species[pid], queue)
        if not r["ok"]:
            waiting.append("%s (%s)" % (pid, r["note"] or r["state"]))
    if waiting:
        return None, "waiting for Step 1 of %s" % ", ".join(waiting)
    coords = [None] * len(sp.mol)
    fits = []
    for pid, idx in parts:
        opt = chem.read_xyz(species_result(project, project.species[pid], queue)["xyz"])
        mine = set(idx)
        others = [c for a, c in enumerate(sp.mol.coords) if a not in mine]
        target = [sp.mol.coords[a] for a in idx]
        w = [math.exp(-min(math.dist(c, o) for o in others) / CONTACT_SCALE) + 1e-3 for c in target] if others else None
        rot, t = chem.superpose(opt.coords, target, w)
        placed = chem.transform(opt.coords, rot, t)
        for k, a in enumerate(idx):
            coords[a] = placed[k]
        fits.append("%s (moved %.2f A on average from where species/%s.xyz has it)" % (pid, chem.rmsd(placed, target), sp.id))
    mol = chem.Molecule(sp.mol.symbols, coords, "%s built from the optimised %s" % (sp.id, " + ".join(p for p, _ in parts)))
    part_of = {}
    for k, (_, idx) in enumerate(parts):
        for a in idx:
            part_of[a] = k
    clash = [(i, j) for i, j in chem.bonds(mol, project.tolerance) if part_of[i] != part_of[j]]
    if clash:
        i, j = clash[0]
        return None, ("cannot be built from its molecules: their optimised structures put %s and %s %.2f A apart (as "
                      "if bonded); move the molecules further apart in species/%s.xyz" % (
                          project.atom_label(mol, i), project.atom_label(mol, j), chem.distance(mol.coords, i, j), sp.id))
    i, j = min(((a, b) for a in range(len(mol)) for b in range(a + 1, len(mol)) if part_of[a] != part_of[b]),
               key=lambda ab: chem.distance(mol.coords, *ab))
    notes = fits + ["closest contact between molecules: %s...%s %.2f A (%.2f A in species/%s.xyz)" % (
        project.atom_label(mol, i), project.atom_label(mol, j), chem.distance(mol.coords, i, j),
        chem.distance(sp.mol.coords, i, j), sp.id)]
    return mol, notes


def start_structure(project, sid, queue=None):
    """The structure a SCAN starts from: the optimised species, or a complex built from its optimised molecules.
    Returns {"ok", "note", "mol", "xyz" (a file, or None for a built complex), "built" (notes, or None)}."""
    sp = project.species[sid]
    if assembled(project, sp):
        mol, notes = assemble(project, sp, queue)
        if mol is None:
            return {"ok": False, "note": notes, "mol": None, "xyz": None, "built": None}
        return {"ok": True, "note": "", "mol": mol, "xyz": None, "built": notes}
    r = species_result(project, sp, queue)
    return {"ok": r["ok"], "note": r["note"] or r["state"], "mol": chem.read_xyz(r["xyz"]) if r["ok"] else None,
            "xyz": r["xyz"] if r["ok"] else None, "built": None}


def sp_setup(project, item, kind="species", queue=None):
    """Write the single point (SP) of a species (from its Step 1 structure) or of a step's TS (from its Step 3
    structure) once that structure is ready. Returns a message or None."""
    if not sp_enabled(project):
        return None
    if kind == "species":
        sp = item
        if not sp_needed(project, sp):
            return None
        d, r = species_sp_dir(project, sp.id), species_result(project, sp, queue)
        xyz = r["xyz"] if r["ok"] else None
        charge, mult, name, over, label = sp.charge, sp.mult, sp.id, sp.overrides, sp.id
    else:
        st = item
        d, r = step_sp_dir(project, st.id), ts_result(project, st)
        xyz = r["xyz"] if r else None
        charge, mult, name, over, label = st.charge, st.mult, st.id, st.overrides, st.ts
    if xyz is None or job_state(d, queue) != "absent":
        return None
    os.makedirs(d, exist_ok=True)
    shutil.copy(xyz, os.path.join(d, "opt_structure.xyz"))
    write_job(d, project.settings, "sp", "opt_structure.xyz", charge, mult, jobname(project, name, "SP"),
              overrides=over, info={"item": name, "stage": "SP"})
    return "single point (SP) of %s written" % label


def sp_plan(project, species, steps, queue=None, all_ts=False):
    """What 'mechflow sp' would do now: the single points it would write [(kind, item)] - a species once its Step 1 is
    usable (if it needs one: sp_needed), a TS once its step is DONE (with all_ts: once the OptTS has one imaginary
    frequency) - the prepared single points it would submit (folders), and how many structures are not ready yet."""
    write, ready, waiting = [], [], 0
    for sp in species:
        if not sp_needed(project, sp):
            continue
        d = species_sp_dir(project, sp.id)
        s = job_state(d, queue)
        if s == "absent":
            if species_result(project, sp, queue)["ok"]:
                write.append(("species", sp))
            else:
                waiting += 1
        elif s == "ready":
            ready.append(d)
    for st in steps:
        d = step_sp_dir(project, st.id)
        s = job_state(d, queue)
        if s == "absent":
            if ts_result(project, st) is not None and (all_ts or step_progress(project, st, queue)[0] == "DONE"):
                write.append(("step", st))
            else:
                waiting += 1
        elif s == "ready":
            ready.append(d)
    return write, ready, waiting


def sp_energy(folder, queue=None):
    """(state, E in Eh or None) of a single point (SP)"""
    state = job_state(folder, queue)
    if state != "finished":
        return state, None
    return state, parse_output(os.path.join(folder, "output.out"))["e_el"]


def species_geometry(project, sid, queue=None):
    """Optimised structure of a species if available, else its xyz file."""
    sp = project.species[sid]
    r = species_result(project, sp, queue)
    if r["xyz"]:
        return chem.read_xyz(r["xyz"])
    return sp.mol


# ============================================================================== SCAN analysis

def scan_folders(stepdir):
    base = os.path.join(stepdir, SCAN)

    def num(p):
        m = re.search(r"SCAN_(\d+)$", p)
        return int(m.group(1)) if m else 0
    subs = sorted((p for p in glob.glob(os.path.join(base, "SCAN_*")) if os.path.isdir(p)), key=num)
    return subs if subs else [base]


def scan_points(stepdir, i, j):
    """Every converged SCAN point (restarts included, in order) -> list of dicts."""
    pts = []
    for k, folder in enumerate(scan_folders(stepdir)):
        for img in sorted(glob.glob(os.path.join(folder, "orca.[0-9][0-9][0-9].xyz"))):
            try:
                mol = chem.read_xyz(img)
            except (OSError, chem.XYZError):
                continue
            e = chem.energy_from_comment(mol.comment)
            if e is None:
                continue
            pts.append({"folder": folder, "file": img, "image": int(re.findall(r"orca\.(\d+)\.xyz$", img)[0]),
                        "part": k + 1, "d": chem.distance(mol.coords, i, j), "E": e})
    clean = []
    for p in pts:          # a restarted SCAN repeats its first point
        if clean and abs(p["d"] - clean[-1]["d"]) < 1e-3 and p["part"] != clean[-1]["part"]:
            clean[-1] = p
        else:
            clean.append(p)
    return clean


def scan_maximum(pts):
    """Index of the highest interior local maximum, or None."""
    best = None
    for k in range(1, len(pts) - 1):
        if pts[k]["E"] >= pts[k - 1]["E"] and pts[k]["E"] >= pts[k + 1]["E"]:
            if best is None or pts[k]["E"] > pts[best]["E"]:
                best = k
    return best


def scan_profile_text(pts, mark=None, width=46):
    if not pts:
        return "   (no converged SCAN points yet)"
    emin = min(p["E"] for p in pts)
    span = max(max(p["E"] for p in pts) - emin, 1e-12)
    out = ["   %5s %4s %8s %10s  relative energy" % ("image", "part", "d / A", "kJ/mol")]
    for k, p in enumerate(pts):
        bar = "#" * int(round(width * (p["E"] - emin) / span))
        out.append("   %5d %4d %8.3f %10.1f  %s%s" % (p["image"], p["part"], p["d"], (p["E"] - emin) * HARTREE_KJ, bar,
                                                   "  <== TS guess" if k == mark else ""))
    return "\n".join(out)


# ============================================================================== steps

def step_progress(project, st, queue=None):
    """(stage, detail, action) for a step. action is what 'next' would prepare, or None."""
    sd = project.step_dir(st.id)
    scan_dirs = scan_folders(sd)
    s2 = job_state(scan_dirs[-1], queue)
    if s2 == "absent":
        sres = start_structure(project, st.start, queue)
        if sres["built"] is not None:
            return "SCAN", "start structure %s built from its optimised molecules" % st.start, "prepare_scan"
        if sres["ok"]:
            return "SCAN", "start structure %s optimised" % st.start, "prepare_scan"
        if assembled(project, project.species[st.start]):
            return "waiting", "for %s: %s" % (st.start, sres["note"]), None
        return "waiting", "for Step 1 of %s (%s)" % (st.start, sres["note"]), None
    if s2 in ("ready", "submitted", "running"):
        return "SCAN", s2, None
    i, j = st.atoms
    pts = scan_points(sd, i, j)
    t3 = os.path.join(sd, OPTTS)
    s3 = job_state(t3, queue)
    if s3 == "absent":
        if s2 == "failed" and not pts:
            return "SCAN", "failed (%s)" % (slurm_problem(scan_dirs[-1]) or "see output.out"), None
        k = scan_maximum(pts)
        if k is None:
            note = "no energy maximum in %d SCAN point(s)" % len(pts)
            note += ("; the SCAN stopped early: 'mechflow reset %s' runs it on from its last point" % st.id) \
                if s2 == "failed" else \
                " - barrierless here, or scan further / another distance"
            return "SCAN", note, None
        note = "maximum at image %d (d = %.3f A)" % (pts[k]["image"], pts[k]["d"])
        if s2 == "failed":
            note += " [the SCAN stopped early but passed a maximum]"
        return "OptTS", note, "prepare_optts"
    if s3 in ("ready", "submitted", "running"):
        return "OptTS", s3, None
    r3 = parse_output(os.path.join(t3, "output.out"))
    if s3 == "failed" or not r3["opt_converged"]:
        return "OptTS", "did not converge / failed (%s) - see %s/output.out" % (
            slurm_problem(t3) or "ORCA", OPTTS), None
    if len(r3["imag"]) != 1:
        return "OptTS", "%d imaginary frequencies (%s): a TS has exactly 1%s" % (
            len(r3["imag"]), ", ".join("%.0f" % f for f in r3["imag"]),
            " (%s)" % imaginary_hint(project, st.id, t3).replace("the mode", "the extra mode")
            if len(r3["imag"]) > 1 else ""), None
    t4a = os.path.join(sd, IRC)
    s4a = job_state(t4a, queue)
    if s4a == "absent":
        return "IRC", "TS found (%.0f cm-1)" % r3["imag"][0], "prepare_irc"
    if s4a in ("ready", "submitted", "running"):
        return "IRC", s4a, None
    if not (os.path.isfile(os.path.join(t4a, "orca_IRC_F.xyz")) and os.path.isfile(os.path.join(t4a, "orca_IRC_B.xyz"))):
        return "IRC", "failed: no IRC end points written (%s)" % (slurm_problem(t4a) or "see output.out"), None
    fw, bw = os.path.join(sd, VALIDATE, "Forwards"), os.path.join(sd, VALIDATE, "Backwards")
    sf, sb = job_state(fw, queue), job_state(bw, queue)
    if sf == "absent" or sb == "absent":
        r4 = parse_output(os.path.join(t4a, "output.out"))
        note = "IRC converged in both directions" if r4["irc_converged"] >= 2 else \
            "IRC finished (%d of 2 directions converged)" % r4["irc_converged"]
        if sf != sb:
            note += "; %s to write again" % ("Forwards" if sf == "absent" else "Backwards")
        return "Fwd/Bwd", note, "prepare_4b"
    if sf != "finished" or sb != "finished":
        return "Fwd/Bwd", "Forwards %s, Backwards %s" % (sf, sb), None
    ok, msg = validate_4b(project, st, queue)
    if not ok and is_accepted(os.path.join(sd, VALIDATE)):
        return "DONE", "accepted by you: " + msg, None
    return ("DONE" if ok else "CHECK"), msg, None


def _side_fingerprints(project, counter, queue):
    fps = []
    for sid, n in counter.items():
        mol = species_geometry(project, sid, queue)
        fps += chem.fragment_fingerprints(mol, project.tolerance) * n
    return sorted(fps)


def validate_4b(project, st, queue=None):
    """Step 4B checks: both end points converged with no imaginary frequency, the scanned pair is bonded
    at one end only, and the two end points are the two sides of the step's reaction."""
    sd = project.step_dir(st.id)
    i, j = st.atoms
    msgs, ok, ends = [], True, {}
    for side in ("Forwards", "Backwards"):
        d = os.path.join(sd, VALIDATE, side)
        r = parse_output(os.path.join(d, "output.out"))
        if not r["opt_converged"]:
            ok = False
            msgs.append("%s did not converge" % side)
        if r["imag"]:
            ok = False
            msgs.append("%s has %d imaginary frequency(ies) (%s)" % (side, len(r["imag"]), imaginary_hint(project, st.id, d)))
        try:
            ends[side] = chem.read_xyz(os.path.join(d, "orca.xyz"))
        except (OSError, chem.XYZError):
            ok = False
            msgs.append("%s/orca.xyz missing" % side)
    if len(ends) == 2:
        bf = chem.is_bonded(ends["Forwards"], i, j, project.tolerance)
        bb = chem.is_bonded(ends["Backwards"], i, j, project.tolerance)
        df = chem.distance(ends["Forwards"].coords, i, j)
        db = chem.distance(ends["Backwards"].coords, i, j)
        pair = "%s-%s" % (project.atom_label(ends["Forwards"], i), project.atom_label(ends["Forwards"], j))
        if bf == bb:
            ok = False
            msgs.append("%s is %s at BOTH ends (%.2f / %.2f A)" % (pair, "bonded" if bf else "not bonded", df, db))
        else:
            msgs.append("%s bonded %s (%.2f A), not %s (%.2f A)" % (
                pair, "Forwards" if bf else "Backwards", df if bf else db, "Backwards" if bf else "Forwards",
                db if bf else df))
        left = _side_fingerprints(project, st.left, queue)
        right = _side_fingerprints(project, st.right, queue)
        found = {}
        for side, mol in ends.items():
            fp = chem.fragment_fingerprints(mol, project.tolerance)
            found[side] = "left" if fp == left else ("right" if fp == right else None)
        if sorted(v or "" for v in found.values()) == ["left", "right"]:
            fl = "Forwards" if found["Forwards"] == "left" else "Backwards"
            msgs.append("connects %s (%s) and %s (%s)" % (
                project_expression(st.left), fl, project_expression(st.right),
                "Backwards" if fl == "Forwards" else "Forwards"))
        else:
            ok = False
            for side in ("Backwards", "Forwards"):
                if found[side] is None:
                    frags = " + ".join(f for _, f in chem.fragment_fingerprints(ends[side], project.tolerance))
                    msgs.append("the %s end (%s) is neither %s nor %s" % (
                        side, frags, project_expression(st.left), project_expression(st.right)))
            if found["Forwards"] and found["Forwards"] == found["Backwards"]:
                msgs.append("both ends are %s" % project_expression(st.left if found["Forwards"] == "left" else st.right))
    return ok, "; ".join(msgs)


def project_expression(counter):
    from .project import expression
    return expression(counter)


def prepare(project, st, action, image=None):
    """Write the input of the next stage of a step. Returns a one-line message."""
    sd = project.step_dir(st.id)
    i, j = st.atoms
    if action == "prepare_scan":
        sres = start_structure(project, st.start)
        if not sres["ok"]:
            raise ValueError("the start structure %s is not ready: %s" % (st.start, sres["note"]))
        folder = os.path.join(sd, SCAN)
        os.makedirs(folder, exist_ok=True)
        mol = sres["mol"]
        if sres["xyz"]:
            shutil.copy(sres["xyz"], os.path.join(folder, "start_opt.xyz"))
        else:
            chem.write_xyz(os.path.join(folder, "start_opt.xyz"), mol, "%s  charge=%d mult=%d" % (
                mol.comment, st.charge, st.mult))
            with open(os.path.join(folder, "start_built_from_molecules.txt"), "w") as fh:
                fh.write("start_opt.xyz is %s built from the optimised structures (Step 1) of its molecules, each fitted "
                         "onto where species/%s.xyz has it (complexes.optimise: false):\n%s\n" % (
                             st.start, st.start, "\n".join("  " + n for n in sres["built"])))
        d0 = round(chem.distance(mol.coords, i, j), 3)
        n = max(int(round(abs(st.to - d0) / st.scan_step)) + 1, 3)
        d1 = d0 + math.copysign(st.scan_step * (n - 1), st.to - d0)
        write_job(folder, project.settings, "scan", "start_opt.xyz", st.charge, st.mult, jobname(project, st.id, "SCAN"),
                  scan=(i, j, d0, d1, n), overrides=st.overrides, info={"item": st.id, "stage": "SCAN"})
        _write_step_info(project, st)
        return "SCAN %s-%s: %.3f -> %.3f A in %d points" % (project.atom_label(mol, i), project.atom_label(mol, j), d0, d1, n)
    if action == "prepare_optts":
        pts = scan_points(sd, i, j)
        if image is not None:
            sel = [p for p in pts if p["image"] == image]
            if not sel:
                raise ValueError("SCAN image %d not found (see 'mechflow scan-profile %s')" % (image, st.id))
            p = sel[-1]
        else:
            k = scan_maximum(pts)
            if k is None:
                raise ValueError("no energy maximum in the SCAN of %s: choose a point with --image" % st.id)
            p = pts[k]
        folder = os.path.join(sd, OPTTS)
        os.makedirs(folder, exist_ok=True)
        name = os.path.basename(p["file"])
        shutil.copy(p["file"], os.path.join(folder, name))
        with open(os.path.join(folder, "TS_guess_from_SCAN.txt"), "w") as fh:
            fh.write("TS guess: %s (SCAN part %d, image %d, d = %.3f A, E = %.8f Eh)%s\n" % (
                os.path.relpath(p["file"], sd), p["part"], p["image"], p["d"], p["E"],
                " - chosen with --image" if image is not None else " - highest SCAN maximum"))
        write_job(folder, project.settings, "optts", name, st.charge, st.mult, jobname(project, st.id, "OptTS"),
                  overrides=st.overrides, info={"item": st.id, "stage": "OptTS"})
        return "OptTS from SCAN image %d (d = %.3f A)" % (p["image"], p["d"])
    if action == "prepare_irc":
        folder = os.path.join(sd, IRC)
        os.makedirs(folder, exist_ok=True)
        shutil.copy(os.path.join(sd, OPTTS, "orca.xyz"), os.path.join(folder, "TS.xyz"))
        write_job(folder, project.settings, "irc", "TS.xyz", st.charge, st.mult, jobname(project, st.id, "IRC"),
                  overrides=st.overrides, info={"item": st.id, "stage": "IRC"})
        return "IRC from the optimised %s" % st.ts
    if action == "prepare_4b":
        written = []
        for side, f, short in (("Forwards", "orca_IRC_F.xyz", "Fwd"), ("Backwards", "orca_IRC_B.xyz", "Bwd")):
            folder = os.path.join(sd, VALIDATE, side)
            if os.path.isfile(os.path.join(folder, "orca.inp")):
                continue                    # this side is there already (only the other one was reset)
            written.append(side)
            os.makedirs(folder, exist_ok=True)
            shutil.copy(os.path.join(sd, IRC, f), os.path.join(folder, f))
            write_job(folder, project.settings, "opt", f, st.charge, st.mult, jobname(project, st.id, short),
                      overrides=st.overrides, info={"item": st.id, "stage": short})
        if len(written) == 2:
            return "Forwards and Backwards optimisations from the IRC end points"
        return "%s optimisation from its IRC end point" % " ".join(written)
    raise ValueError("unknown action %s" % action)


def _write_step_info(project, st):
    sd = project.step_dir(st.id)
    os.makedirs(sd, exist_ok=True)
    sp = project.species[st.start]
    i, j = st.atoms
    lines = ["%s  (transition state %s)" % (st.id, st.ts), st.name, "",
             "Reaction:      %s" % st.reaction_text(),
             "Used by:       %s" % (", ".join(st.used_by) or "-"),
             ("SCAN starts from %s built from the optimised structures of its molecules (%s), placed as in species/%s.xyz"
              % (st.start, " + ".join(p for p, _ in complex_parts(project, sp)), st.start)
              if assembled(project, sp) else
              "SCAN starts from the optimised structure of %s (calcs/Species/%s/%s/orca.xyz)" % (st.start, st.start, STEP1)),
             "SCAN distance: %s-%s (atoms numbered from %d, as in project.yaml; ORCA's input counts from 0), "
             "%s to %.2f A in %.3f A steps" % (
                 project.atom_label(sp.mol, i), project.atom_label(sp.mol, j), project.index_base,
                 "stretched" if st.mode == "stretch" else "brought together", st.to, st.scan_step),
             "Charge %d, multiplicity %d" % (st.charge, st.mult), "",
             "Run 'mechflow next' after each stage: it checks the stage and writes the next input."]
    with open(os.path.join(sd, "step_info.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def retire(folder):
    """Rename a stage folder to <folder>_old1, _old2 ... so that it can be redone."""
    if os.path.isdir(folder):
        k = 1
        while os.path.exists("%s_old%d" % (folder, k)):
            k += 1
        os.rename(folder, "%s_old%d" % (folder, k))
        return "%s_old%d" % (folder, k)
    return None


def restart_scan(project, st):
    """The procedure's SCAN restart: files move to SCAN_1, SCAN_2 starts from the last converged point."""
    base = os.path.join(project.step_dir(st.id), SCAN)
    if not os.path.isdir(base):
        raise ValueError("%s has no SCAN yet" % st.id)
    subs = [p for p in glob.glob(os.path.join(base, "SCAN_*")) if os.path.isdir(p)]
    if not subs:
        first = os.path.join(base, "SCAN_1")
        os.makedirs(first)
        for f in os.listdir(base):
            if f not in ("SCAN_1", ATTEMPTS):
                shutil.move(os.path.join(base, f), os.path.join(first, f))
        subs = [first]
    last = sorted(subs, key=lambda p: int(re.findall(r"SCAN_(\d+)$", p)[0]))[-1]
    k = int(re.findall(r"SCAN_(\d+)$", last)[0])
    imgs = sorted(glob.glob(os.path.join(last, "orca.[0-9][0-9][0-9].xyz")))
    if not imgs:
        raise ValueError("no converged SCAN point in %s: 'mechflow reset %s' runs it again as it was "
                         "(--from-scratch starts the SCAN again)" % (last, st.id))
    m = re.search(r"SCAN B\s+\d+\s+\d+\s*=\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*(\d+)", read_text(os.path.join(last, "orca.inp")))
    if not m:
        raise ValueError("cannot read the SCAN line of %s/orca.inp" % last)
    d_first, d_final, n_old = float(m.group(1)), float(m.group(2)), int(m.group(3))
    step = abs(d_final - d_first) / max(n_old - 1, 1)
    new = os.path.join(base, "SCAN_%d" % (k + 1))
    os.makedirs(new)
    shutil.copy(imgs[-1], os.path.join(new, "last_converged_SCAN_image_restart.xyz"))
    mol = chem.read_xyz(imgs[-1])
    i, j = st.atoms
    d0 = round(chem.distance(mol.coords, i, j), 3)
    n = max(int(round(abs(d_final - d0) / step)) + 1, 2)
    write_job(new, project.settings, "scan", "last_converged_SCAN_image_restart.xyz", st.charge, st.mult,
              jobname(project, st.id, "SCAN%d" % (k + 1)), scan=(i, j, d0, d_final, n), overrides=st.overrides,
              info={"item": st.id, "stage": "SCAN%d" % (k + 1)})
    return "SCAN_%d from %s: %.3f -> %.3f A in %d points (%.3f A apart)" % (
        k + 1, os.path.basename(imgs[-1]), d0, d_final, n, step)


def job_folders(project, species, steps, with_sp=True):
    """Every job folder of the chosen species and steps, in a sensible order (with_sp=False: not the single points)."""
    for sp in species:
        if not assembled(project, sp):           # (a built complex's Step 1, written earlier, is not submitted)
            yield project.species_job_dir(sp.id)
        if with_sp and sp_needed(project, sp):   # (a complex's single point, written by an older mechflow, is not submitted)
            yield species_sp_dir(project, sp.id)
    for st in steps:
        sd = project.step_dir(st.id)
        for f in scan_folders(sd):
            yield f
        yield os.path.join(sd, OPTTS)
        yield os.path.join(sd, IRC)
        yield os.path.join(sd, VALIDATE, "Forwards")
        yield os.path.join(sd, VALIDATE, "Backwards")
        if with_sp:
            yield os.path.join(sd, SP)


def refresh(project, folders, queue=None):
    """Rewrite the inputs of prepared, not yet submitted jobs from the current settings."""
    n = 0
    for f in folders:
        if job_state(f, queue) == "ready" and rewrite_job(f, project.settings):
            n += 1
    return n
