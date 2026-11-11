"""chain.enabled: each finished job runs 'mechflow chain', which writes and submits what can start now."""
import glob
import os
import time

import fake_orca
import pytest
from mechflow.cli import CHAIN_LOCK, main
from mechflow.orca import job_state
from mechflow.project import Project
from mechflow.workflow import step_progress

SETUP = ["module load Python/3.11.6-foss-2023a", "source ~/mechflow-venv/bin/activate"]


def run(root, *args):
    return main(["-p", root] + list(args))


def enable_chain(root):
    with open(os.path.join(root, "project.yaml"), "a") as fh:
        fh.write("\nchain:\n  enabled: true\n  setup: [%s]\n" % ", ".join('"%s"' % s for s in SETUP))


def test_job_scripts_chain_only_when_asked(hcn):
    run(hcn, "next", "--no-squeue")
    plain = open(os.path.join(hcn, "calcs/Species/R1/Step1_Geo_Opt/submit.sl")).read()
    assert "mechflow chain" not in plain and plain.rstrip().endswith("${orca_exe} orca.inp > output.out")
    enable_chain(hcn)
    run(hcn, "next", "--refresh", "--no-squeue")
    sl = open(os.path.join(hcn, "calcs/Species/R1/Step1_Geo_Opt/submit.sl")).read()
    assert "${orca_exe} orca.inp > output.out || orca_status=$?" in sl
    assert 'grep -q "ORCA TERMINATED NORMALLY" output.out' in sl
    for line in SETUP + ["mechflow chain", "exit ${orca_status}"]:
        assert line in sl, line
    assert sl.index("module load ORCA") < sl.index(SETUP[0]) < sl.index("mechflow chain")


def test_chain_submits_what_can_start(hcn, fake_sbatch, capsys):
    enable_chain(hcn)
    run(hcn, "next", "--no-squeue")
    run(hcn, "submit", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)                        # Step 1 of R1 and P1 finished
    capsys.readouterr()
    run(hcn, "chain", "--no-squeue")                         # what the last Step 1 job runs
    out = capsys.readouterr().out
    assert "SCAN H2-N1" in out and "Submitted 1 job(s)." in out
    scan = os.path.join(hcn, "calcs/Steps/S1/Step2Ba_SCAN")
    assert job_state(scan) == "submitted"
    assert not os.path.exists(os.path.join(hcn, CHAIN_LOCK))
    run(hcn, "chain", "--no-squeue")                         # a second finished job: nothing new, nothing twice
    assert "No prepared jobs" in capsys.readouterr().out
    for _ in range(4):                                       # SCAN -> OptTS -> IRC -> Fwd/Bwd, each by chain
        fake_orca.finish_ready_jobs(hcn)
        run(hcn, "chain", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    p = Project.load(hcn)
    assert step_progress(p, p.steps["S1"])[0] == "DONE"


def test_chain_lock(hcn, fake_sbatch, capsys):
    run(hcn, "next", "--no-squeue")
    lock = os.path.join(hcn, CHAIN_LOCK)
    os.mkdir(lock)                                           # another chain is running
    with pytest.raises(SystemExit):
        run(hcn, "chain", "--no-squeue", "--wait", "0.2")
    assert "another 'mechflow chain'" in capsys.readouterr().err
    old = time.time() - 3600                                 # left over from a killed job: taken over
    os.utime(lock, (old, old))
    run(hcn, "chain", "--no-squeue", "--wait", "0.2")
    assert "Submitted 2 job(s)." in capsys.readouterr().out and not os.path.exists(lock)
