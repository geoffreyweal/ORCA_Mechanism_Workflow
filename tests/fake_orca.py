"""Fabricate ORCA outputs so the workflow can be tested without ORCA (a test helper only).

finish_ready_jobs(root) writes output.out (and the xyz files ORCA would write) for every job
folder under calcs/ that mechflow has prepared but that has not run yet.
"""
import glob
import hashlib
import json
import math
import os

from mechflow import chem
from mechflow.orca import META, job_state
from mechflow.project import Project


def out_text(G, imag=None, opt=True, irc=0, term=True, scan_steps=0):
    t = "                         Program Version 6.1.0-f.0  -   RELEASE   -\n"
    t += "".join("         *               RELAXED SURFACE SCAN STEP %3d               *\n" % (k + 1)
                 for k in range(scan_steps))
    if opt:
        t += "                    ***        THE OPTIMIZATION HAS CONVERGED     ***\n"
    t += "THE IRC HAS CONVERGED\n" * irc
    t += "FINAL SINGLE POINT ENERGY     %.12f\n" % (G - 0.1)
    t += ("-----------------------\nVIBRATIONAL FREQUENCIES\n-----------------------\n\n"
          "Scaling factor for frequencies =  1.000000000  (already applied!)\n\n")
    for k in range(6):
        t += "   %2d:       0.00 cm**-1\n" % k
    if imag is not None:
        t += "    6:    %.2f cm**-1  ***imaginary mode***\n" % imag
    for k in range(7, 12):
        t += "   %2d:     %6.2f cm**-1\n" % (k, 30.0 + 10 * k)
    t += "\nFinal Gibbs free energy         ...  %.8f Eh\n" % G
    if term:
        t += "\n                             ****ORCA TERMINATED NORMALLY****\n"
    return t


def pseudo_gibbs(mol, charge, salt=""):
    """A made-up but reproducible energy: about -0.5 Eh per electron, plus a small structure-dependent part."""
    fp = chem.molecule_fingerprint(mol) + salt
    jitter = int(hashlib.sha1(fp.encode()).hexdigest()[:6], 16) / 16 ** 6 * 0.02 - 0.01
    return -0.5 * chem.electron_count(mol.symbols, 0) - 0.15 * charge + jitter


def _w(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def _write_mol(path, mol, comment):
    chem.write_xyz(path, mol, comment)


def _separate(mol, i, j, gap=2.5):
    """Pull apart the two pieces that are left when the i-j bond is broken (if there are two)."""
    b = chem.bonds(mol) - {(min(i, j), max(i, j))}
    comps = chem.fragments(len(mol), b)
    comp_j = [c for c in comps if j in c][0]
    if i in comp_j:
        return None
    v = [mol.coords[j][k] - mol.coords[i][k] for k in range(3)]
    L = math.sqrt(sum(x * x for x in v)) or 1.0
    shift = [x / L * gap for x in v]
    coords = [tuple(c[k] + (shift[k] if a in comp_j else 0.0) for k in range(3)) for a, c in enumerate(mol.coords)]
    return chem.Molecule(mol.symbols, coords)


def finish_ready_jobs(root, die_after=None, ts_height=0.03):
    """Fake every prepared job. die_after=N makes SCANs stop (time limit) after N points."""
    project = Project.load(root)
    n = 0
    for meta_path in sorted(glob.glob(os.path.join(root, "calcs", "**", META), recursive=True)):
        folder = os.path.dirname(meta_path)
        if job_state(folder) != "ready":
            continue
        meta = json.load(open(meta_path))
        kind, item = meta["kind"], meta.get("item")
        mol = chem.read_xyz(os.path.join(folder, meta["xyz"]))
        q = meta["charge"]
        n += 1
        if kind == "opt":
            _write_mol(os.path.join(folder, "orca.xyz"), mol, "Coordinates from ORCA-job orca")
            _w(os.path.join(folder, "output.out"), out_text(pseudo_gibbs(mol, q)))
        elif kind == "scan":
            i, j, d0, d1, npts = meta["scan"]
            limit = npts if die_after is None else min(die_after, npts)
            e0 = pseudo_gibbs(mol, q)
            for k in range(limit):
                target = d0 + (d1 - d0) * k / (npts - 1)
                v = [mol.coords[j][c] - mol.coords[i][c] for c in range(3)]
                L = math.sqrt(sum(x * x for x in v))
                new = list(mol.coords)
                new[j] = tuple(mol.coords[i][c] + v[c] / L * target for c in range(3))
                img = chem.Molecule(mol.symbols, new)
                e = e0 + ts_height * math.sin(math.pi * k / (npts - 1))
                _write_mol(os.path.join(folder, "orca.%03d.xyz" % (k + 1)), img,
                           "Coordinates from ORCA-job orca E %.12f" % e)
            if limit == npts:
                _w(os.path.join(folder, "output.out"), out_text(e0, scan_steps=npts))
            else:
                _w(os.path.join(folder, "output.out"), "RELAXED SURFACE SCAN STEP %3d\n" % limit)
                _w(os.path.join(folder, "slurm-123.err"),
                   "slurmstepd: error: *** JOB 123 ON node CANCELLED AT 2026-01-01 DUE TO TIME LIMIT ***\n")
        elif kind == "optts":
            st = project.steps[item]
            _write_mol(os.path.join(folder, "orca.xyz"), mol, "Coordinates from ORCA-job orca")
            g = sum(pseudo_gibbs(project.species[s].mol, project.species[s].charge) * c for s, c in st.left.items())
            _w(os.path.join(folder, "output.out"), out_text(g + ts_height, imag=-350.0))
        elif kind == "irc":
            st = project.steps[item]
            sd = project.step_dir(item)
            start = chem.read_xyz(os.path.join(sd, "Step2Ba_SCAN", "start_opt.xyz")) if os.path.isfile(
                os.path.join(sd, "Step2Ba_SCAN", "start_opt.xyz")) else chem.read_xyz(
                glob.glob(os.path.join(sd, "Step2Ba_SCAN", "SCAN_1", "start_opt.xyz"))[0])
            other_side = st.right if (list(st.left) == [st.start] and sum(st.left.values()) == 1) else st.left
            other = None
            if sum(other_side.values()) == 1:
                cand = project.species[list(other_side)[0]].mol
                if cand.symbols == start.symbols:
                    other = cand
            if other is None:
                i, j = st.atoms
                other = _separate(start, i, j) if st.mode == "stretch" else None
            if other is None:
                imgs = sorted(glob.glob(os.path.join(sd, "Step2Ba_SCAN", "**", "orca.[0-9][0-9][0-9].xyz"),
                                        recursive=True))
                other = chem.read_xyz(imgs[-1])
            _write_mol(os.path.join(folder, "orca_IRC_B.xyz"), start, "IRC backwards end")
            _write_mol(os.path.join(folder, "orca_IRC_F.xyz"), other, "IRC forwards end")
            _w(os.path.join(folder, "output.out"), out_text(-1.0, opt=False, irc=2))
    return n
