"""'mechflow reset': failed jobs set up to run again, the files of each attempt kept in failed_attempts/."""
import glob
import json
import os

import fake_orca
import pytest
from mechflow import chem
from mechflow.cli import main
from mechflow.orca import META, job_state
from mechflow.project import Project
from mechflow.workflow import step_progress

TIME_LIMIT = "slurmstepd: error: *** JOB 555 ON node CANCELLED AT 2026-01-01 DUE TO TIME LIMIT ***\n"


def run(root, *args):
    return main(["-p", root] + list(args))


def stage(root):
    p = Project.load(root)
    return step_progress(p, p.steps["S1"])


def attempts(folder):
    return sorted(glob.glob(os.path.join(folder, "failed_attempts", "*")))


def fail_optimisation(folder, moved=0.05):
    """An optimisation stopped by the time limit after moving its atoms (orca.xyz is its last geometry)."""
    meta = json.load(open(os.path.join(folder, META)))
    mol = chem.read_xyz(os.path.join(folder, meta["xyz"]))
    coords = [(x + moved, y, z) for x, y, z in mol.coords]
    chem.write_xyz(os.path.join(folder, "orca.xyz"), chem.Molecule(mol.symbols, coords), "Coordinates from ORCA-job orca")
    open(os.path.join(folder, "output.out"), "w").write("Program Version 6.1.0\nGEOMETRY OPTIMIZATION CYCLE   7\n")
    open(os.path.join(folder, "slurm-555.err"), "w").write(TIME_LIMIT)
    open(os.path.join(folder, ".submitted"), "w").write("555\n")
    return coords


def full_run(root):
    run(root, "next", "--no-squeue")
    for _ in range(6):
        fake_orca.finish_ready_jobs(root)
        run(root, "next", "--no-squeue")
    assert stage(root)[0] == "DONE"


def test_reset_restarts_an_optimisation_from_its_last_geometry(hcn, capsys):
    run(hcn, "next", "--no-squeue")
    d = Project.load(hcn).species_job_dir("R1")
    start = open(os.path.join(d, "start_structure.xyz")).read()
    coords = fail_optimisation(d)
    assert job_state(d) == "failed"
    capsys.readouterr()
    run(hcn, "status", "--no-squeue")
    assert "1 failed job(s): 'mechflow reset'" in capsys.readouterr().out
    run(hcn, "reset", "--no-squeue", "--dry-run")
    out = capsys.readouterr().out
    assert "R1      Step1  failed: time limit" in out and "restarts from the last geometry" in out
    assert "1 job(s) would be reset" in out
    assert job_state(d) == "failed" and not attempts(d)                  # a dry run changes nothing
    run(hcn, "reset", "--no-squeue")
    out = capsys.readouterr().out
    assert "Reset 1 job(s)." in out and "Next: 'mechflow submit'." in out
    assert job_state(d) == "ready"
    [old] = attempts(d)
    for f in ("output.out", "slurm-555.err", ".submitted", "orca.xyz", "orca.inp", "submit.sl", "start_structure.xyz"):
        assert os.path.isfile(os.path.join(old, f)), f
    restart = chem.read_xyz(os.path.join(d, "last_geometry_restart.xyz"))
    assert all(abs(a - b) < 1e-6 for c1, c2 in zip(restart.coords, coords) for a, b in zip(c1, c2))
    assert "* xyzfile 0 1 last_geometry_restart.xyz" in open(os.path.join(d, "orca.inp")).read()
    assert open(os.path.join(d, "start_structure.xyz")).read() == start      # still the species' own structure
    meta = json.load(open(os.path.join(d, META)))
    assert (meta["xyz"], meta["first_xyz"], meta["restarts"]) == ("last_geometry_restart.xyz", "start_structure.xyz", 1)
    run(hcn, "next", "--no-squeue")                  # 'next' keeps the restart (the species file has not changed)
    assert "updated" not in capsys.readouterr().out
    assert "last_geometry_restart.xyz" in open(os.path.join(d, "orca.inp")).read()
    fake_orca.finish_ready_jobs(hcn)
    assert stage(hcn)[2] == "prepare_scan"           # the restarted Step 1 is used as usual


def test_reset_runs_a_stopped_scan_on(hcn, capsys):
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn, die_after=12)          # the SCAN hits the time limit before its maximum
    assert "'mechflow reset S1'" in stage(hcn)[1]
    capsys.readouterr()
    run(hcn, "reset", "--no-squeue")
    out = capsys.readouterr().out
    assert "S1      SCAN   failed: time limit" in out and "SCAN_2 from orca.012.xyz" in out
    scan = os.path.join(hcn, "calcs/Steps/S1/Step2Ba_SCAN")
    assert job_state(os.path.join(scan, "SCAN_2")) == "ready"
    assert os.path.isfile(os.path.join(scan, "SCAN_1", "slurm-123.err"))       # the stopped part stays as it was
    fake_orca.finish_ready_jobs(hcn)
    assert stage(hcn)[2] == "prepare_optts"
    run(hcn, "reset", "--no-squeue")
    assert "Nothing to reset" in capsys.readouterr().out                      # SCAN_1 is only an earlier part


def test_reset_from_scratch_with_the_later_stages(hcn, capsys):
    full_run(hcn)
    capsys.readouterr()
    run(hcn, "reset", "--no-squeue")
    assert "Nothing to reset" in capsys.readouterr().out                      # nothing failed
    with pytest.raises(SystemExit):
        run(hcn, "reset", "--no-squeue", "--redo-done")                       # finished jobs only when named
    capsys.readouterr()
    run(hcn, "reset", "R1", "--redo-done", "--dry-run", "--no-squeue")
    out = capsys.readouterr().out
    assert "5 later stage(s) used its result (S1 SCAN, S1 OptTS, S1 IRC, S1 Fwd, S1 Bwd)" in out
    run(hcn, "reset", "R1", "--redo-done", "--from-scratch", "--downstream", "--no-squeue")
    out = capsys.readouterr().out
    assert "Reset 6 job(s)." in out and "Next: 'mechflow next'" in out
    assert "S1      OptTS  uses the result of R1 Step1" in out
    p = Project.load(hcn)
    assert job_state(p.species_job_dir("R1")) == "absent" and attempts(p.species_job_dir("R1"))
    assert job_state(p.species_job_dir("P1")) == "finished"                  # not used by S1's SCAN: left alone
    for sub in ("Step2Ba_SCAN", "Step3_Opt_TS", "Step4A_IRC", "Step4B_Validate/Forwards", "Step4B_Validate/Backwards"):
        folder = os.path.join(p.step_dir("S1"), sub)
        assert os.listdir(folder) == ["failed_attempts"], sub
    assert stage(hcn)[0] == "waiting"
    full_run(hcn)                                    # and the whole step runs again


def test_reset_one_side_of_step_4b(hcn, capsys):
    full_run(hcn)
    run(hcn, "reset", "S1", "--stage", "Fwd", "--redo-done", "--from-scratch", "--no-squeue")
    st, detail, action = stage(hcn)
    assert st == "Fwd/Bwd" and action == "prepare_4b" and "Forwards to write again" in detail
    capsys.readouterr()
    run(hcn, "next", "--no-squeue")
    assert "Forwards optimisation from its IRC end point" in capsys.readouterr().out
    sd = Project.load(hcn).step_dir("S1")
    assert job_state(os.path.join(sd, "Step4B_Validate/Forwards")) == "ready"
    assert job_state(os.path.join(sd, "Step4B_Validate/Backwards")) == "finished"
    fake_orca.finish_ready_jobs(hcn)
    assert stage(hcn)[0] == "DONE"


def test_reset_leaves_queued_jobs_alone(hcn, fake_sbatch, capsys):
    run(hcn, "next", "--no-squeue")
    run(hcn, "submit", "--no-squeue")                         # R1 and P1 submitted: jobs 1001 and 1002
    capsys.readouterr()
    run(hcn, "reset", "R1", "--no-squeue")                    # squeue not asked: the state is unknown
    out = capsys.readouterr().out
    assert "squeue could not be asked" in out and "Nothing to reset" in out
    squeue = os.path.join(fake_sbatch, "squeue")              # squeue: both jobs are still queued
    open(squeue, "w").write("#!/bin/sh\necho 1001\necho 1002\n")
    os.chmod(squeue, 0o755)
    run(hcn, "reset", "R1", "--force")
    out = capsys.readouterr().out
    assert "R1 Step1 is queued in Slurm (job 1001)" in out and "Nothing to reset" in out    # never, even with --force
    d = Project.load(hcn).species_job_dir("R1")
    open(os.path.join(d, "output.out"), "w").write("ORCA finished by error termination in SCF\n")
    assert job_state(d) == "failed"                           # an error in the output, but Slurm still has the job
    run(hcn, "reset", "--force")
    out = capsys.readouterr().out
    assert "Nothing to reset" in out and os.path.isfile(os.path.join(d, "output.out"))
    run(hcn, "reset", "R1", "--force")
    assert "R1 Step1 is running in Slurm (job 1001)" in capsys.readouterr().out
    os.remove(os.path.join(d, "output.out"))
    run(hcn, "reset", "R1", "--force", "--no-squeue")         # unknown state + --force: reset
    assert "Reset 1 job(s)." in capsys.readouterr().out
    p = Project.load(hcn)
    assert job_state(p.species_job_dir("R1")) == "ready" and job_state(p.species_job_dir("P1")) == "submitted"
