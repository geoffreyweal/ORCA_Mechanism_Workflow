"""Step 2 by NEB-CI (neb.enabled) and a chemist in the loop (approval.enabled)."""
import os

import fake_orca
from mechflow import chem
from mechflow.cli import main
from mechflow.neb import parse_neb
from mechflow.orca import job_state
from mechflow.project import Project
from mechflow.workflow import step_progress
from test_complexes import HCN, WATER_IN_C9, add_complex, write, xyz


def run(root, *args):
    return main(["-p", root] + list(args))


def settings(root, neb=True, approval=True, images=4):
    y = os.path.join(root, "project.yaml")
    text = open(y).read()
    if neb:
        text += "\nneb:\n  enabled: true\n  images: %d\n" % images
    if approval:
        text += "\napproval:\n  enabled: true\n"
    write(y, text)


def hcn_ends(root):
    """the two ends of S1 (HCN -> HNC): R1 and P1 as they are, same atoms in the same order"""
    os.makedirs(os.path.join(root, "neb"), exist_ok=True)
    for side, sid in (("left", "R1"), ("right", "P1")):
        mol = chem.read_xyz(os.path.join(root, "species", sid + ".xyz"))
        chem.write_xyz(os.path.join(root, "neb", "S1_%s.xyz" % side), mol, "S1 %s end" % side)


def progress(root, sid):
    p = Project.load(root)
    return step_progress(p, p.steps[sid])


def test_neb_route_waits_for_the_chemist(hcn, fake_sbatch, capsys):
    settings(hcn)
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)                                   # Step 1 of R1 and P1
    stage, detail, action = progress(hcn, "S1")
    assert stage == "waiting" and "no NEB end file neb/S1_left.xyz" in detail
    hcn_ends(hcn)
    stage, detail, action = progress(hcn, "S1")
    assert stage == "approve" and "--start" in detail and action is None
    run(hcn, "next", "--no-squeue")
    assert job_state(os.path.join(hcn, "calcs/Steps/S1/Step2Bb_NEB")) == "absent"          # nothing without approval
    run(hcn, "status", "--no-squeue")
    assert "1 step(s) wait for you" in capsys.readouterr().out
    run(hcn, "approve", "S1", "--start", "--note", "looked fine", "--no-squeue")
    assert "approved the start" in capsys.readouterr().out
    assert "looked fine" in open(os.path.join(hcn, "calcs/Steps/S1/approved_start.txt")).read()
    stage, detail, action = progress(hcn, "S1")
    assert (stage, action) == ("NEB", "prepare_neb") and "left as given, right as given" in detail
    run(hcn, "next", "--no-squeue")
    inp = open(os.path.join(hcn, "calcs/Steps/S1/Step2Bb_NEB/orca.inp")).read()
    assert "!NEB-CI" in inp and '%NEB\n    NEB_END_XYZFILE "right_end.xyz"\n    NImages 4' in inp and "Interpolation IDPP" in inp
    assert "* xyzfile 0 1 left_end.xyz" in inp
    run(hcn, "submit", "--no-squeue")
    assert "Step2Bb_NEB" in capsys.readouterr().out
    fake_orca.finish_ready_jobs(hcn)                                   # the NEB
    r = parse_neb(os.path.join(hcn, "calcs/Steps/S1/Step2Bb_NEB/output.out"))
    assert r["converged"] and r["ci"] == 2 and len(r["path"]) == 6
    stage, detail, action = progress(hcn, "S1")
    assert stage == "approve" and "--guess" in detail and "climbing image 2 of 5" in detail
    run(hcn, "next", "--no-squeue")
    assert job_state(os.path.join(hcn, "calcs/Steps/S1/Step3_Opt_TS")) == "absent"         # the chemist decides first
    run(hcn, "approve", "S1", "--guess", "--no-squeue")
    run(hcn, "next", "--no-squeue")
    t3 = os.path.join(hcn, "calcs/Steps/S1/Step3_Opt_TS")
    assert job_state(t3) == "ready" and "climbing image" in open(os.path.join(t3, "TS_guess_from_NEB.txt")).read()
    ci = chem.read_xyz(os.path.join(hcn, "calcs/Steps/S1/Step2Bb_NEB/orca_NEB-CI_converged.xyz"))
    assert chem.rmsd(chem.read_xyz(os.path.join(t3, "NEB_TS_guess.xyz")).coords, ci.coords) < 1e-6
    # then on by itself: OptTS -> IRC -> Step 4B
    fake_orca.finish_ready_jobs(hcn)
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    assert progress(hcn, "S1")[0] == "DONE"


def test_an_approved_start_holds_only_for_the_ends_it_was_given_for(hcn, capsys):
    settings(hcn)
    hcn_ends(hcn)
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    run(hcn, "approve", "S1", "--start", "--no-squeue")
    assert progress(hcn, "S1")[0] == "NEB"
    with open(os.path.join(hcn, "neb", "S1_right.xyz"), "a") as fh:
        fh.write("\n")                                                  # the end file changed
    stage, detail, _ = progress(hcn, "S1")
    assert stage == "approve" and "the end files changed after you approved them" in detail
    run(hcn, "approve", "S1", "--start", "--undo", "--no-squeue")
    assert not os.path.exists(os.path.join(hcn, "calcs/Steps/S1/approved_start.txt"))


def test_approvals_listed_in_project_yaml_and_none_needed_without_approval(hcn):
    settings(hcn, approval=False)
    hcn_ends(hcn)
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    assert progress(hcn, "S1")[:3:2] == ("NEB", "prepare_neb")         # approval.enabled false: straight on
    y = os.path.join(hcn, "project.yaml")
    write(y, open(y).read() + "\napproval:\n  enabled: true\n  start: [S1]\n")
    assert progress(hcn, "S1")[0] == "NEB"


def test_neb_ends_of_two_molecules_are_optimised_first(hcn, capsys):
    add_complex(hcn)                                                    # S2: R1 + W -> P5 + P6 (HCN + water -> CN- + H3O+)
    settings(hcn, approval=False)
    hcn_ends(hcn)
    left = HCN + WATER_IN_C9
    right = [HCN[0], HCN[1], ("H", (-2.05, 0.36, 0.0))] + WATER_IN_C9   # the proton on the water: CN- ... H3O+
    write(os.path.join(hcn, "neb", "S2_left.xyz"), xyz(left, "S2 left"))
    write(os.path.join(hcn, "neb", "S2_right.xyz"), xyz(right, "S2 right"))
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)                                   # Step 1
    stage, detail, action = progress(hcn, "S2")
    assert (stage, action) == ("NEB ends", "prepare_neb_ends") and \
        detail == "optimised first: left (more than one molecule), right (more than one molecule)"
    run(hcn, "next", "--no-squeue")
    ends = os.path.join(hcn, "calcs/Steps/S2/Step2Bb_Ends")
    assert "!OPT TightOPT" in open(os.path.join(ends, "Left", "orca.inp")).read()
    assert job_state(os.path.join(ends, "Right")) == "ready"
    fake_orca.finish_ready_jobs(hcn)                                   # the two end optimisations
    # the right end fell back (the proton went back to the carbon): the step waits
    moved = [HCN[0], HCN[1], ("H", (-1.06, 0.19, 0.0))] + WATER_IN_C9
    write(os.path.join(ends, "Right", "orca.xyz"), xyz(moved, "Coordinates from ORCA-job orca"))
    stage, detail, action = progress(hcn, "S2")
    assert stage == "CHECK" and "the right end changed during its optimisation" in detail and "mechflow accept S2" in detail
    run(hcn, "accept", "S2")
    assert "the optimised NEB ends of S2" in capsys.readouterr().out
    stage, detail, action = progress(hcn, "S2")
    assert (stage, action) == ("NEB", "prepare_neb") and "left optimised, right optimised" in detail
    run(hcn, "next", "--no-squeue")
    nd = os.path.join(hcn, "calcs/Steps/S2/Step2Bb_NEB")
    assert chem.rmsd(chem.read_xyz(os.path.join(nd, "left_end.xyz")).coords,
                     chem.read_xyz(os.path.join(ends, "Left", "orca.xyz")).coords) < 1e-6


def test_neb_ends_must_be_the_two_sides(hcn):
    settings(hcn, approval=False)
    hcn_ends(hcn)
    os.replace(os.path.join(hcn, "neb", "S1_left.xyz"), os.path.join(hcn, "neb", "S1_right.xyz"))
    os.makedirs(os.path.join(hcn, "neb"), exist_ok=True)
    write(os.path.join(hcn, "neb", "S1_left.xyz"), xyz(HCN[:2], "only CN"))
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    stage, detail, _ = progress(hcn, "S1")
    assert stage == "CHECK" and "the same atoms in the same order" in detail


def test_an_neb_that_did_not_converge_runs_on_from_its_last_path(hcn, capsys):
    settings(hcn)
    hcn_ends(hcn)
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    run(hcn, "approve", "S1", "--start", "--no-squeue")
    run(hcn, "approve", "S1", "--guess", "--no-squeue")                 # in advance (taken back by the reset)
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn, neb_converged=False)
    stage, detail, _ = progress(hcn, "S1")
    assert stage == "NEB" and "did not converge" in detail and "--stage NEB" in detail
    run(hcn, "reset", "S1", "--stage", "NEB", "--no-squeue", "--force")
    out = capsys.readouterr().out
    assert "runs on from its last path" in out
    nd = os.path.join(hcn, "calcs/Steps/S1/Step2Bb_NEB")
    assert 'Restart_ALLXYZFile "restart_path.allxyz"' in open(os.path.join(nd, "orca.inp")).read()
    assert os.path.isfile(os.path.join(nd, "restart_path.allxyz"))
    assert not os.path.exists(os.path.join(hcn, "calcs/Steps/S1/approved_guess.txt"))


def test_scan_route_with_a_chemist_in_the_loop(hcn, capsys):
    settings(hcn, neb=False)
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    stage, detail, action = progress(hcn, "S1")
    assert stage == "approve" and "start structure R1 optimised" in detail and "--start" in detail
    run(hcn, "approve", "S1", "--start", "--no-squeue")
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)                                   # the SCAN
    stage, detail, action = progress(hcn, "S1")
    assert stage == "approve" and "maximum at image" in detail and "--guess" in detail
    run(hcn, "approve", "S1", "--guess", "--no-squeue")
    assert progress(hcn, "S1")[2] == "prepare_optts"


def test_scan_folders_are_left_alone_on_the_neb_route(hcn, fake_sbatch, capsys):
    """SCANs written (and cancelled) before neb.enabled was switched on are neither reset nor submitted."""
    settings(hcn, neb=False, approval=False)
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    run(hcn, "next", "--no-squeue")                                    # the SCAN of S1 written ...
    scan = os.path.join(hcn, "calcs/Steps/S1/Step2Ba_SCAN")
    write(os.path.join(scan, "slurm-9.err"), "slurmstepd: error: *** JOB 9 ON n CANCELLED AT 2026-01-01 ***\n")
    write(os.path.join(scan, ".submitted"), "9\n")
    write(os.path.join(scan, "output.out"), "cancelled\n")
    settings(hcn, approval=False)                                      # ... then the NEB route switched on
    hcn_ends(hcn)
    capsys.readouterr()
    run(hcn, "reset", "--dry-run", "--no-squeue", "--force")
    assert "SCAN" not in capsys.readouterr().out
    run(hcn, "next", "--no-squeue")
    run(hcn, "submit", "--no-squeue")
    out = capsys.readouterr().out
    assert "Step2Bb_NEB" in out and "Step2Ba_SCAN" not in out


def test_an_end_of_one_molecule_marked_optimise_first_is_optimised(hcn, capsys):
    """'optimise first' in an end file's comment line: an end of one molecule (a conformer set up for the reaction) is
    optimised before the NEB too."""
    settings(hcn, approval=False)
    hcn_ends(hcn)
    mol = chem.read_xyz(os.path.join(hcn, "neb", "S1_left.xyz"))
    chem.write_xyz(os.path.join(hcn, "neb", "S1_left.xyz"), mol, "S1 left end, a conformer set up for it: optimise first")
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    stage, detail, action = progress(hcn, "S1")
    assert (stage, action) == ("NEB ends", "prepare_neb_ends") and detail == "optimised first: left (marked 'optimise first')"
    run(hcn, "next", "--no-squeue")
    ends = os.path.join(hcn, "calcs/Steps/S1/Step2Bb_Ends")
    assert os.path.isfile(os.path.join(ends, "Left", "orca.inp")) and not os.path.exists(os.path.join(ends, "Right"))
    fake_orca.finish_ready_jobs(hcn)
    assert progress(hcn, "S1")[1] == "ends ready (left optimised, right as given)"


def test_a_step_can_have_its_own_number_of_images(hcn, fake_sbatch, capsys):
    settings(hcn, approval=False, images=4)
    y = os.path.join(hcn, "project.yaml")
    text = open(y).read()
    assert "    scan:" in text
    write(y, text.replace("    scan:", "    neb_images: 100\n    scan:", 1))
    p = Project.load(hcn)
    assert not p.errors and p.steps["S1"].overrides == {"neb_images": 100}
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)                                   # Step 1 of R1 and P1
    hcn_ends(hcn)
    run(hcn, "next", "--no-squeue")
    folder = os.path.join(hcn, "calcs/Steps/S1/Step2Bb_NEB")
    inp = open(os.path.join(folder, "orca.inp")).read()
    assert "    NImages 100   # images between the two ends" in inp
    assert "neb_images" not in open(os.path.join(folder, "submit.sl")).read()
    capsys.readouterr()
    run(hcn, "validate")
    assert "4 images (own: S1 100)" in capsys.readouterr().out
    write(y, open(y).read().replace("neb_images: 100", "neb_images: 0"))
    assert any("neb_images must be a whole number" in e for e in Project.load(hcn).errors)


def test_an_neb_can_start_through_a_ts_guess(hcn, fake_sbatch, capsys):
    settings(hcn, approval=True, images=4)
    y = os.path.join(hcn, "project.yaml")
    write(y, open(y).read().replace("  images: 4\n", '  images: 4\n  ts: "neb/{step}_ts.xyz"\n', 1))
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)                                   # Step 1 of R1 and P1
    hcn_ends(hcn)
    run(hcn, "approve", "S1", "--start", "--no-squeue")
    digest = open(os.path.join(hcn, "calcs/Steps/S1/approved_start.txt")).read()
    mol = chem.read_xyz(os.path.join(hcn, "neb", "S1_left.xyz"))
    chem.write_xyz(os.path.join(hcn, "neb", "S1_ts.xyz"), mol, "S1 TS guess")
    stage, detail, _ = progress(hcn, "S1")
    assert stage == "approve" and "TS guess neb/S1_ts.xyz" in detail      # a new file: the approval no longer holds
    assert digest
    run(hcn, "approve", "S1", "--start", "--no-squeue")
    run(hcn, "next", "--no-squeue")
    folder = os.path.join(hcn, "calcs/Steps/S1/Step2Bb_NEB")
    inp = open(os.path.join(folder, "orca.inp")).read()
    assert 'NEB_TS_XYZFILE "ts_guess.xyz"' in inp and os.path.isfile(os.path.join(folder, "ts_guess.xyz"))
    assert "ts_guess.xyz" in open(os.path.join(folder, "NEB_ends.txt")).read()
    capsys.readouterr()
    run(hcn, "validate")
    assert "through a TS guess (neb.ts) for 1 step(s): S1" in capsys.readouterr().out
    bad = chem.Molecule(list(reversed(mol.symbols)), list(reversed(mol.coords)))
    chem.write_xyz(os.path.join(hcn, "neb", "S1_ts.xyz"), bad, "S1 TS guess, atoms in another order")
    os.rename(folder, folder + "_away")
    stage, detail, _ = progress(hcn, "S1")
    assert stage == "CHECK" and "does not have the atoms of the NEB ends" in detail
