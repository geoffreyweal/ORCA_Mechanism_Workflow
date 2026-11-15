"""Complexes built from their optimised molecules (complexes.optimise: false), and 'mechflow reset --imaginary'."""
import json
import math
import os

import fake_orca
import pytest
from mechflow import chem
from mechflow.cli import main
from mechflow.orca import META, job_state
from mechflow.project import Project
from mechflow.workflow import assemble, assembled, complex_parts, start_structure, step_progress


def run(root, *args):
    return main(["-p", root] + list(args))


def write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


# HCN hydrogen-bonded to a water molecule; the same water, turned and moved, is the species W
HCN = [("C", (0.0, 0.0, 0.0)), ("N", (1.153, 0.0, 0.0)), ("H", (-1.0498, 0.1851, 0.0))]
WATER_IN_C9 = [("O", (-2.921, 0.515, 0.0)), ("H", (-3.881, 0.515, 0.0)), ("H", (-2.681, 1.445, 0.0))]


def turned(points, angle=1.1, shift=(4.0, -3.0, 2.5)):
    """the points turned about the axis (1, 2, 2)/3 by `angle` and moved by `shift`"""
    ux, uy, uz = 1 / 3.0, 2 / 3.0, 2 / 3.0
    c, s = math.cos(angle), math.sin(angle)
    rot = [[c + ux * ux * (1 - c), ux * uy * (1 - c) - uz * s, ux * uz * (1 - c) + uy * s],
           [uy * ux * (1 - c) + uz * s, c + uy * uy * (1 - c), uy * uz * (1 - c) - ux * s],
           [uz * ux * (1 - c) - uy * s, uz * uy * (1 - c) + ux * s, c + uz * uz * (1 - c)]]
    return [(e, tuple(sum(rot[i][j] * p[j] for j in range(3)) + shift[i] for i in range(3))) for e, p in points]


def xyz(atoms, comment):
    return "%d\n%s\n%s" % (len(atoms), comment, "".join("%s %.6f %.6f %.6f\n" % ((e,) + tuple(p)) for e, p in atoms))


def add_complex(root, optimise=False, c9=None):
    """W (water), the complex C9 (HCN...OH2), CN- and H3O+, and S2: the proton of HCN moves to the water."""
    sp = os.path.join(root, "species")
    write(os.path.join(sp, "W.xyz"), xyz(turned(WATER_IN_C9), "W water  charge=0 mult=1"))
    write(os.path.join(sp, "C9.xyz"), xyz(c9 or (HCN + WATER_IN_C9), "C9 HCN...OH2  charge=0 mult=1"))
    write(os.path.join(sp, "P5.xyz"), xyz(HCN[:2], "P5 cyanide  charge=-1 mult=1"))
    write(os.path.join(sp, "P6.xyz"), xyz([("O", (0, 0, 0)), ("H", (0.97, 0, 0)), ("H", (-0.32, 0.92, 0)),
                                           ("H", (-0.32, -0.46, 0.8))], "P6 hydronium  charge=1 mult=1"))
    y = os.path.join(root, "project.yaml")
    text = open(y).read().replace("  P1: hydrogen isocyanide\n", "  P1: hydrogen isocyanide\n  W: water\n"
                                  "  C9: HCN hydrogen-bonded to water\n  P5: cyanide\n  P6: hydronium\n")
    text = text.replace("mechanisms:\n", "  S2:\n    name: proton from HCN to water\n    reaction: R1 + W -> P5 + P6\n"
                        "    scan: {from: C9, atoms: [2, 3], to: 1.00}\n\nmechanisms:\n")
    write(y, text + "\ncomplexes:\n  optimise: %s\n" % ("true" if optimise else "false"))


def test_complex_built_from_its_optimised_molecules(hcn, capsys):
    add_complex(hcn)
    run(hcn, "validate")
    out = capsys.readouterr().out
    assert "1 complex(es) are built from their optimised molecules (complexes.optimise: false): C9 = R1 + W" in out
    p = Project.load(hcn)
    assert assembled(p, p.species["C9"]) and not assembled(p, p.species["W"])
    assert complex_parts(p, p.species["C9"]) == [("R1", [0, 1, 2]), ("W", [3, 4, 5])]
    run(hcn, "next", "--no-squeue")
    out = capsys.readouterr().out
    assert "C9 " not in out and not os.path.exists(os.path.join(hcn, "calcs/Species/C9"))       # no Step 1 for C9
    assert step_progress(p, p.steps["S2"])[1].startswith("for C9: waiting for Step 1 of R1 (ready), W (ready)")
    run(hcn, "status", "--no-squeue")
    out = capsys.readouterr().out
    assert "C9      waiting    built from R1 + W (no Step 1 of its own); waiting for Step 1 of R1, W" in out
    assert "(1 complex(es) are built from their optimised molecules" in out
    fake_orca.finish_ready_jobs(hcn)                              # Step 1 of R1, W ...
    p = Project.load(hcn)
    mol, notes = assemble(p, p.species["C9"])
    assert chem.rmsd(mol.coords, p.species["C9"].mol.coords) < 1e-6     # W turned back to where C9 has it
    assert "closest contact between molecules: H2...O3 1.90 A (1.90 A in species/C9.xyz)" in notes[-1]
    run(hcn, "next", "--no-squeue")
    out = capsys.readouterr().out
    assert "S2      SCAN H2-O3: 1.900 -> 1.000 A in 46 points" in out
    scan = os.path.join(hcn, "calcs/Steps/S2/Step2Ba_SCAN")
    start = chem.read_xyz(os.path.join(scan, "start_opt.xyz"))
    assert chem.rmsd(start.coords, p.species["C9"].mol.coords) < 1e-5
    assert "built from the optimised structures (Step 1) of its molecules" in open(
        os.path.join(scan, "start_built_from_molecules.txt")).read()
    assert "S2 built from the optimised structures of its molecules" not in open(
        os.path.join(hcn, "calcs/Steps/S2/step_info.txt")).read()
    assert "SCAN starts from C9 built from the optimised structures of its molecules (R1 + W)" in open(
        os.path.join(hcn, "calcs/Steps/S2/step_info.txt")).read()
    run(hcn, "status", "--no-squeue")
    out = capsys.readouterr().out
    assert "C9      built      built from the optimised R1 + W when its SCAN starts" in out
    # W redone: the SCAN that started from C9 used its structure
    run(hcn, "reset", "W", "--redo-done", "--stage", "Step1", "--no-squeue", "--dry-run")
    assert "W Step1: 1 later stage(s) used its result (S2 SCAN); add --downstream" in capsys.readouterr().out


def test_a_complex_that_is_a_state_is_still_optimised(hcn, capsys):
    """A complex that a mechanism has as a state (here the step starts from it) needs its own energy, so its own Step 1."""
    add_complex(hcn)
    y = os.path.join(hcn, "project.yaml")
    text = open(y).read().replace("reaction: R1 + W -> P5 + P6", "reaction: C9 -> P5 + P6")
    write(y, text.replace("mechanisms:\n", "mechanisms:\n  M2:\n    name: through the complex\n"
                                       "    reference: R1 + W\n    path: [R1 + W, S2, P5 + P6]\n", 1))
    p = Project.load(hcn)
    assert not p.errors, p.errors
    assert not assembled(p, p.species["C9"])
    run(hcn, "next", "--no-squeue")
    assert os.path.isfile(os.path.join(hcn, "calcs/Species/C9/Step1_Geo_Opt/orca.inp"))


def test_complex_problems_are_reported(hcn, capsys):
    add_complex(hcn, c9=HCN + [WATER_IN_C9[1], WATER_IN_C9[0], WATER_IN_C9[2]])    # water listed H, O, H
    try:
        run(hcn, "validate")
    except SystemExit:
        pass
    out = capsys.readouterr().out
    assert "ERROR:   complex C9: the molecule H2O of C9 (atoms 3-5) is W, but its atoms are listed in another order" in out
    p = Project.load(hcn)
    assert "cannot be built from its molecules" in start_structure(p, "C9")["note"]
    nh3 = [("N", (-3.0, 0.5, 0.0)), ("H", (-3.5, 1.3, 0.3)), ("H", (-3.5, -0.3, 0.3)), ("H", (-3.2, 0.5, -1.0))]
    add_complex(hcn, c9=HCN + nh3)
    y = os.path.join(hcn, "project.yaml")
    write(y, open(y).read().replace("reaction: R1 + W -> P5 + P6", "reaction: R1 + W -> P5 + P6").replace(
        "  C9: HCN hydrogen-bonded to water\n", "  C9: HCN hydrogen-bonded to ammonia\n", 1))
    p = Project.load(hcn)
    try:
        complex_parts(p, p.species["C9"])
        assert False, "a molecule that is no species"
    except ValueError as exc:
        assert "is not a species of the project" in str(exc)


def test_old_step1_job_of_a_built_complex(hcn, fake_sbatch, capsys):
    """A complex's Step 1, written and submitted before complexes.optimise was set to false, is left alone, and status
    says it can be cancelled."""
    add_complex(hcn, optimise=True)
    run(hcn, "next", "--no-squeue")
    run(hcn, "submit", "--no-squeue")
    d = os.path.join(hcn, "calcs/Species/C9/Step1_Geo_Opt")
    assert job_state(d) == "submitted"
    y = os.path.join(hcn, "project.yaml")
    write(y, open(y).read().replace("optimise: true", "optimise: false"))
    capsys.readouterr()
    run(hcn, "status", "--no-squeue")
    out = capsys.readouterr().out
    jid = open(os.path.join(d, ".submitted")).read().strip()
    assert "its own Step 1 job (written before) is not needed: scancel %s" % jid in out
    assert "1 of them still have a Step 1 job in Slurm, written before, that is not needed: scancel %s" % jid in out
    run(hcn, "reset", "--no-squeue", "--force", "--dry-run")
    assert "C9" not in capsys.readouterr().out                            # reset leaves it alone too


# ----------------------------------------------------------------------------- reset --imaginary

def hess_text(mol, freqs, modes, shift=(1.0, 0.0, 0.0)):
    """An ORCA .hess file: frequencies, normal modes (columns, in blocks of 5) and the structure in bohr (moved by
    `shift` bohr, as ORCA may write it)."""
    n3 = 3 * len(mol)
    t = "\n$orca_hessian_file\n\n$vibrational_frequencies\n%d\n" % len(freqs)
    t += "".join("%5d %16.6f\n" % (k, f) for k, f in enumerate(freqs))
    t += "\n$normal_modes\n%d %d\n" % (n3, len(freqs))
    for c0 in range(0, len(freqs), 5):
        cols = list(range(c0, min(c0 + 5, len(freqs))))
        t += "      " + "".join("%19d" % c for c in cols) + "\n"
        for r in range(n3):
            t += "%5d   " % r + "".join("%19.10E" % modes[c][r] for c in cols) + "\n"
    t += "\n$atoms\n%d\n" % len(mol)
    for e, c in zip(mol.symbols, mol.coords):
        t += " %-2s %10.5f %18.12f %18.12f %18.12f\n" % ((e, 1.0) + tuple(x / 0.529177210903 + s for x, s in zip(c, shift)))
    return t + "\n$end\n"


def give_imaginary(folder, freqs, modes):
    """Make a finished optimisation (fake ORCA) end with these frequencies (imaginary ones negative) and modes."""
    mol = chem.read_xyz(os.path.join(folder, "orca.xyz"))
    out = open(os.path.join(folder, "output.out")).read()
    imag = [f for f in freqs if f < 0]
    out = "".join(ln for ln in out.splitlines(True) if "imaginary mode" not in ln)
    lines = "".join("    %d:    %.2f cm**-1  ***imaginary mode***\n" % (6 + k, f) for k, f in enumerate(imag))
    out = out.replace("   10:", lines + "   10:")
    write(os.path.join(folder, "output.out"), out)
    write(os.path.join(folder, "orca.hess"), hess_text(mol, freqs, modes))
    return mol


def test_reset_imaginary_pushes_a_minimum_along_its_mode(hcn, capsys):
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    d = Project.load(hcn).species_job_dir("R1")
    mode = [0.0] * 9
    mode[7] = 0.5                                       # the H (atom 2) moves along y
    other = [0.1] * 9
    mol = give_imaginary(d, [0.0] * 6 + [-35.0, 700.0, 2100.0], [other] * 6 + [mode, other, other])
    capsys.readouterr()
    run(hcn, "status", "--no-squeue")
    assert "R1      CHECK      1 imaginary frequency(ies): -35.0 ('mechflow reset R1 --imaginary' redoes it" in \
        capsys.readouterr().out
    run(hcn, "reset", "--imaginary", "--no-squeue", "--dry-run")
    out = capsys.readouterr().out
    assert ("R1      Step1  1 imaginary frequency(ies): -35.0 -> restarts from its last geometry pushed 0.30 A along "
            "the imaginary mode (-35.0 cm-1)") in out
    assert job_state(d) == "finished"
    run(hcn, "reset", "--imaginary", "--no-squeue")
    out = capsys.readouterr().out
    assert "Reset 1 job(s)." in out and job_state(d) == "ready"
    new = chem.read_xyz(os.path.join(d, "last_geometry_restart.xyz"))
    moved = [[b - a for a, b in zip(p, q)] for p, q in zip(mol.coords, new.coords)]
    assert all(abs(x) < 1e-6 for x in moved[0] + moved[1]) and abs(moved[2][1] - 0.3) < 1e-6
    meta = json.load(open(os.path.join(d, META)))
    assert (meta["xyz"], meta["restarts"], meta["pushed"]) == ("last_geometry_restart.xyz", 1, 1)
    [old] = [os.path.join(d, "failed_attempts", f) for f in os.listdir(os.path.join(d, "failed_attempts"))]
    assert os.path.isfile(os.path.join(old, "orca.hess")) and os.path.isfile(os.path.join(old, "output.out"))
    run(hcn, "reset", "--imaginary", "--no-squeue")
    assert "Nothing to reset (no finished optimisation has imaginary frequencies" in capsys.readouterr().out
    fake_orca.finish_ready_jobs(hcn)                    # the new attempt, from the pushed structure
    assert "imaginary" not in open(os.path.join(d, "output.out")).read()
    with pytest.raises(SystemExit):
        run(hcn, "reset", "--imaginary", "--no-squeue", "--amplitude", "0")
    assert "--amplitude must not be 0" in capsys.readouterr().err


def test_reset_imaginary_keeps_the_reaction_mode_of_a_ts(hcn, capsys):
    run(hcn, "next", "--no-squeue")
    for _ in range(2):                                  # Step 1, SCAN; then the OptTS is written
        fake_orca.finish_ready_jobs(hcn)
        run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    p = Project.load(hcn)
    t3 = os.path.join(p.step_dir("S1"), "Step3_Opt_TS")
    reaction, methyl = [0.0] * 9, [0.0] * 9
    reaction[6] = 0.6                                   # H moves along x: the reaction
    methyl[8] = 0.4                                     # H moves along z: a second, unwanted imaginary mode
    mol = give_imaginary(t3, [0.0] * 6 + [-350.0, -40.0, 2100.0], [[0.1] * 9] * 6 + [reaction, methyl, [0.1] * 9])
    capsys.readouterr()
    run(hcn, "status", "--no-squeue")
    assert "2 imaginary frequencies (-350, -40): a TS has exactly 1 ('mechflow reset S1 --imaginary'" in \
        capsys.readouterr().out
    run(hcn, "reset", "S1", "--imaginary", "--no-squeue")
    out = capsys.readouterr().out
    assert "OptTS  2 imaginary frequency(ies): -350.0, -40.0 -> restarts from its last geometry pushed 0.30 A along " \
           "the imaginary mode (-40.0 cm-1)" in out
    new = chem.read_xyz(os.path.join(t3, "last_geometry_restart.xyz"))
    moved = [b - a for p_, q in zip(mol.coords, new.coords) for a, b in zip(p_, q)]
    assert abs(moved[8] - 0.3) < 1e-6 and all(abs(x) < 1e-6 for k, x in enumerate(moved) if k != 8)
    assert "last_geometry_restart.xyz" in open(os.path.join(t3, "orca.inp")).read()


def test_next_pushes_imaginary_modes_by_itself(hcn, capsys):
    """imaginary.auto: 'mechflow next' (and so the chain) restarts a minimum with an imaginary frequency pushed along the
    mode, at most max_pushes times; then it waits for you."""
    with open(os.path.join(hcn, "project.yaml"), "a") as fh:
        fh.write("\nimaginary:\n  auto: true\n  max_pushes: 2\n")
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    d = Project.load(hcn).species_job_dir("R1")
    mode = [0.0] * 9
    mode[7] = 0.5
    for k in (1, 2):
        give_imaginary(d, [0.0] * 6 + [-35.0, 700.0, 2100.0], [[0.1] * 9] * 6 + [mode, [0.1] * 9, [0.1] * 9])
        capsys.readouterr()
        run(hcn, "status", "--no-squeue")
        assert "'mechflow next' redoes it pushed along the mode (push %d of 2)" % k in capsys.readouterr().out
        run(hcn, "next", "--no-squeue")
        out = capsys.readouterr().out
        assert ("R1      Step1 1 imaginary frequency(ies): -35.0: restarts from its last geometry pushed 0.30 A along "
                "the imaginary mode (-35.0 cm-1) (last_geometry_restart.xyz) (push %d of 2)" % k) in out
        assert "1 restart(s) pushed along imaginary modes" in out and job_state(d) == "ready"
        assert json.load(open(os.path.join(d, META)))["pushed"] == k
        fake_orca.finish_ready_jobs(hcn)
    give_imaginary(d, [0.0] * 6 + [-35.0, 700.0, 2100.0], [[0.1] * 9] * 6 + [mode, [0.1] * 9, [0.1] * 9])
    capsys.readouterr()
    run(hcn, "next", "--no-squeue")
    assert "push 3" not in capsys.readouterr().out and job_state(d) == "finished"      # twice is enough: it waits
    run(hcn, "status", "--no-squeue")
    assert "pushed 2 time(s) already; 'mechflow reset R1 --imaginary' pushes it again" in capsys.readouterr().out


def test_calculate_skip_steps(hcn, fake_sbatch, capsys):
    """calculate.skip_steps: next, submit and reset leave those steps alone unless named; status lists them as skipped."""
    with open(os.path.join(hcn, "project.yaml"), "a") as fh:
        fh.write("\ncalculate:\n  skip_steps: [S1]\n")
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    capsys.readouterr()
    run(hcn, "next", "--no-squeue")
    assert "S1" not in capsys.readouterr().out
    assert not os.path.exists(os.path.join(hcn, "calcs/Steps/S1"))
    run(hcn, "status", "--no-squeue")
    out = capsys.readouterr().out
    assert "S1      skipped  calculate.skip_steps in project.yaml" in out and "skipped 1" in out
    run(hcn, "next", "--no-squeue", "--step", "S1")          # named: written all the same
    assert os.path.isfile(os.path.join(hcn, "calcs/Steps/S1/Step2Ba_SCAN/orca.inp"))


def test_bimolecular_irc_ends_are_not_optimised(hcn, capsys):
    """complexes.optimise: false: an IRC end of several molecules gets no Step 4B optimisation; the step is checked from
    the IRC end point as it is. An end of one molecule (S1, HCN -> HNC) is optimised as before."""
    add_complex(hcn)
    run(hcn, "next", "--no-squeue")
    for _ in range(3):                                   # Step 1, SCAN, OptTS: then the IRCs are written
        fake_orca.finish_ready_jobs(hcn)
        run(hcn, "next", "--no-squeue")
    p = Project.load(hcn)
    irc = os.path.join(p.step_dir("S2"), "Step4A_IRC")
    assert os.path.isfile(os.path.join(irc, "orca.inp"))
    c9 = p.species["C9"].mol                             # HCN...OH2: C N H O H H
    X = [list(c) for c in c9.coords]
    o, h = X[3], X[2]
    u = [h[k] - o[k] for k in range(3)]
    L = math.sqrt(sum(x * x for x in u))
    fwd = [[x + 2.5 if k == 0 else x for k, x in enumerate(X[0])], [x + 2.5 if k == 0 else x for k, x in enumerate(X[1])],
           [o[k] + u[k] / L * 0.98 for k in range(3)], X[3], X[4], X[5]]           # CN- ... H3O+
    chem.write_xyz(os.path.join(irc, "orca_IRC_F.xyz"), chem.Molecule(c9.symbols, [tuple(c) for c in fwd]), "IRC F")
    chem.write_xyz(os.path.join(irc, "orca_IRC_B.xyz"), c9, "IRC B")
    write(os.path.join(irc, "output.out"), fake_orca.out_text(-1.0, opt=False, irc=2))
    fake_orca.finish_ready_jobs(hcn)                     # the IRC of S1
    capsys.readouterr()
    run(hcn, "next", "--no-squeue")
    out = capsys.readouterr().out
    assert "S1      Forwards and Backwards optimisations from the IRC end points" in out
    assert "S2" not in out
    assert not os.path.exists(os.path.join(p.step_dir("S2"), "Step4B_Validate"))
    stage, detail, _ = step_progress(Project.load(hcn), p.steps["S2"])
    assert stage == "DONE", detail
    assert "Forwards: the IRC end point as it is (several molecules, not optimised)" in detail
    assert "connects R1 + W (Backwards) and P5 + P6 (Forwards)" in detail


def test_species_keywords(hcn, capsys):
    """A species' own keywords line (e.g. a finer grid for one molecule) is used for its Step 1; a reset after it was
    added rewrites the input with it; other species and single points keep the project's."""
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    d = Project.load(hcn).species_job_dir("R1")
    mode = [0.0] * 9
    mode[7] = 0.5
    give_imaginary(d, [0.0] * 6 + [-7.0, 700.0, 2100.0], [[0.1] * 9] * 6 + [mode, [0.1] * 9, [0.1] * 9])
    y = os.path.join(hcn, "project.yaml")
    write(y, open(y).read().replace("  R1: hydrogen cyanide\n",
                                    '  R1: {name: hydrogen cyanide, keywords: "!OPT FREQ VeryTightOPT TightSCF defgrid3"}\n'))
    p = Project.load(hcn)
    assert not p.errors and p.species["R1"].overrides == {"keywords": "!OPT FREQ VeryTightOPT TightSCF defgrid3"}
    run(hcn, "reset", "R1", "--imaginary", "--no-squeue")
    inp = open(os.path.join(d, "orca.inp")).read()
    assert "!OPT FREQ VeryTightOPT TightSCF defgrid3" in inp and "defgrid2" not in inp
    assert "defgrid2" in open(os.path.join(p.species_job_dir("P1"), "orca.inp")).read()
    write(y, open(y).read().replace('keywords: "!OPT', 'keywords: "OPT'))
    assert any("keywords must be an ORCA '!' line" in e for e in Project.load(hcn).errors)


def test_placement_moves_overlapping_molecules_apart():
    """An optimised molecule fitted on top of the other one (its optimised shape differs from the drawn one) is moved
    and turned as a whole until the drawn contact is as drawn and nothing else is too close; its shape is unchanged."""
    from mechflow import placement
    drawn = [p for _, p in HCN + WATER_IN_C9]
    symbols = [e for e, _ in HCN + WATER_IN_C9]
    water = [(x + 4.3, y - 0.2, z + 0.3) for _, (x, y, z) in WATER_IN_C9]                # its O on the N of HCN
    start = drawn[:3] + water
    bodies = [[0, 1, 2], [3, 4, 5]]
    bonds = chem.bonds(chem.Molecule(symbols, drawn))
    assert any(i < 3 <= j for i, j in chem.bonds(chem.Molecule(symbols, start)))           # overlapping at first
    coords, contacts, left = placement.place(symbols, drawn, start, bodies, bonds, seed="C9")
    assert [(i, j) for i, j, _ in contacts] == [(2, 3)]                                    # the H...O drawn at 1.90 A
    assert abs(math.dist(coords[2], coords[3]) - 1.90) < 0.02
    polar = placement.polar_hydrogens(symbols, bonds)
    assert all(math.dist(coords[i], coords[j]) > 0.98 * placement.minimum(symbols, polar, i, j)
               for i in range(3) for j in range(3, 6) if (i, j) != (2, 3))
    assert coords[:3] == [list(p) for p in drawn[:3]]                                      # the larger/first stays
    for a in range(3, 6):                                                                  # the water is rigid
        for b in range(3, 6):
            assert abs(math.dist(coords[a], coords[b]) - math.dist(water[a - 3], water[b - 3])) < 1e-9
    again = placement.place(symbols, drawn, start, bodies, bonds, seed="C9")[0]
    assert again == coords                                                                 # deterministic


def test_a_complex_with_atoms_too_close_is_not_built(hcn, capsys):
    """placement.too_close finds atoms of different molecules closer than TOO_CLOSE x their minimum (here a water put
    onto the HCN), which makes assemble refuse the complex; the drawn contacts are left out of the check."""
    from mechflow import placement
    symbols = [e for e, _ in HCN + WATER_IN_C9]
    coords = [p for _, p in HCN] + [(x + 0.9, y - 0.1, z) for _, (x, y, z) in WATER_IN_C9]   # H...O 1.0 A, H...H close
    part_of = {0: 0, 1: 0, 2: 0, 3: 1, 4: 1, 5: 1}
    polar = placement.polar_hydrogens(symbols, chem.bonds(chem.Molecule(symbols, coords)))
    worst = placement.too_close(symbols, coords, part_of, polar, contact=[])
    assert worst is not None and worst[2] < placement.TOO_CLOSE * worst[3]
    fine = [p for _, p in HCN + WATER_IN_C9]                                                  # as drawn: H...O 1.90 A
    assert placement.too_close(symbols, fine, part_of, polar, contact=[(2, 3, 1.90)]) is None
