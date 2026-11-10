"""The whole procedure on the HCN example, with fabricated ORCA outputs (tests/fake_orca.py)."""
import glob
import os

import fake_orca
from mechflow.cli import main
from mechflow.orca import job_state
from mechflow.project import Project
from mechflow.workflow import step_progress


def run(root, *args):
    return main(["-p", root] + list(args))


def stage(root):
    p = Project.load(root)
    return step_progress(p, p.steps["S1"])


def test_full_procedure(hcn, capsys):
    run(hcn, "next", "--no-squeue")
    assert os.path.isfile(os.path.join(hcn, "calcs/Species/R1/Step1_Geo_Opt/orca.inp"))
    assert stage(hcn)[0] == "waiting"
    run(hcn, "submit", "--no-squeue")
    sh = open(os.path.join(hcn, "submit_ready.sh")).read()
    assert sh.count("sbatch --parsable submit.sl") == 2
    assert fake_orca.finish_ready_jobs(hcn) == 2
    assert stage(hcn)[2] == "prepare_scan"
    for expected in ("OptTS", "IRC", "Fwd/Bwd", "DONE"):
        run(hcn, "next", "--no-squeue")
        fake_orca.finish_ready_jobs(hcn)
        if expected != "DONE":
            assert stage(hcn)[0] == expected
    st, detail, action = stage(hcn)
    assert st == "DONE", detail
    assert "connects R1 (Backwards) and P1 (Forwards)" in detail
    scan = open(os.path.join(hcn, "calcs/Steps/S1/Step2Ba_SCAN/orca.inp")).read()
    assert "SCAN B 2 1 = 2.211, 0.991, 62 END" in scan
    run(hcn, "collect", "--no-squeue")
    out = capsys.readouterr().out
    assert "| M1 | P1 | TS1 |" in out
    for f in ("species.csv", "steps.csv", "profiles.csv", "mechanisms.csv", "summary.md"):
        assert os.path.isfile(os.path.join(hcn, "results", f))
    run(hcn, "map", "--no-squeue", "--status", "--energies", "--format", "dot")
    assert "TS1" in open(os.path.join(hcn, "maps", "network.dot")).read()
    assert "flowchart LR" in open(os.path.join(hcn, "maps", "network.md")).read()


def test_scan_restart_and_image_choice(hcn, capsys):
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn, die_after=12)          # the SCAN hits the time limit
    st, detail, _ = stage(hcn)
    assert st == "SCAN" and "restart-scan" in detail
    run(hcn, "restart-scan", "S1")
    scan2 = os.path.join(hcn, "calcs/Steps/S1/Step2Ba_SCAN/SCAN_2")
    assert job_state(scan2) == "ready"
    assert "SCAN B 2 1 = 1.991, 0.991, 51 END" in open(os.path.join(scan2, "orca.inp")).read()
    fake_orca.finish_ready_jobs(hcn)
    st, detail, action = stage(hcn)
    assert action == "prepare_optts", detail
    run(hcn, "scan-profile", "S1")
    assert "<== TS guess" in capsys.readouterr().out
    run(hcn, "next", "S1", "--no-squeue")
    first = glob.glob(os.path.join(hcn, "calcs/Steps/S1/Step3_Opt_TS/orca.*.xyz"))
    assert len(first) == 1
    run(hcn, "next", "S1", "--image", "40", "--force", "--no-squeue")
    assert os.path.isdir(os.path.join(hcn, "calcs/Steps/S1/Step3_Opt_TS_old1"))
    assert "image 40" in open(os.path.join(hcn, "calcs/Steps/S1/Step3_Opt_TS/TS_guess_from_SCAN.txt")).read()


def test_refresh_and_changed_species(hcn, capsys):
    run(hcn, "next", "--no-squeue")
    inp = os.path.join(hcn, "calcs/Species/R1/Step1_Geo_Opt/orca.inp")
    yaml_path = os.path.join(hcn, "project.yaml")
    text = open(yaml_path).read().replace('method: "!B3LYP DEF2-SVP D3BJ"', 'method: "!PBE0 DEF2-SVP D3BJ"')
    open(yaml_path, "w").write(text)
    run(hcn, "next", "--refresh", "--no-squeue")
    assert open(inp).read().startswith("!PBE0 DEF2-SVP D3BJ")
    # a changed xyz file is picked up while the job has not been submitted
    xyz = os.path.join(hcn, "species/R1.xyz")
    text = open(xyz).read().replace("-1.04980000", "-1.05000000")
    open(xyz, "w").write(text)
    run(hcn, "next", "--no-squeue")
    assert "updated" in capsys.readouterr().out


def test_accept_species_with_changed_bonds(hcn, capsys):
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    # pretend the R1 optimisation moved the H onto N
    opt = os.path.join(hcn, "calcs/Species/R1/Step1_Geo_Opt/orca.xyz")
    open(opt, "w").write("3\nE\nC 0 0 0\nN 1.169 0 0\nH 2.163 0.02 0\n")
    st, detail, _ = stage(hcn)
    assert st == "waiting" and "bonds changed" in detail
    run(hcn, "accept", "R1")
    assert stage(hcn)[2] == "prepare_scan"


def test_init_writes_a_valid_project(tmp_path, capsys):
    root = tmp_path / "new"
    os.makedirs(root / "species")
    src = os.path.join(os.path.dirname(os.path.dirname(__file__)), "examples", "HCN_isomerisation", "species")
    for f in ("R1.xyz", "P1.xyz"):
        (root / "species" / f).write_text(open(os.path.join(src, f)).read())
    main(["-p", str(root), "init"])
    p = Project.load(str(root))
    assert not p.errors and list(p.species) == ["R1", "P1"]
    assert p.settings["orca"]["solvent_block"].strip().endswith("%CPCM EPSILON 4.81 REFRAC 1.4459 END")


def test_single_points(hcn, capsys):
    """single_point.enabled: 'next' writes a single point (SP) for each optimised species and TS; 'collect' then
    uses E(SP) + [G - E(el)] of the DFT job, or the DFT G with --level dft."""
    yaml_path = os.path.join(hcn, "project.yaml")
    with open(yaml_path, "a") as fh:
        fh.write("\nsingle_point:\n  enabled: true\n"
                 "stage_overrides:\n  sp: {maxcore: 4000, mem_per_cpu: 5000MB}\n")
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    run(hcn, "next", "--no-squeue")                 # Step 1 done: the species' single points and the SCAN
    sp_r1 = os.path.join(hcn, "calcs/Species/R1/SP_Single_Point")
    inp = open(os.path.join(sp_r1, "orca.inp")).read()
    assert inp.startswith("!DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C TightPNO\n!RIJCOSX def2/J TightSCF defgrid2\n")
    dft = open(os.path.join(hcn, "calcs/Species/R1/Step1_Geo_Opt/orca.inp")).read()
    assert "%maxcore 4000" in inp and "OPT" not in inp
    assert [ln for ln in inp.splitlines() if "CPCM" in ln or "SMD" in ln] == \
        [ln for ln in dft.splitlines() if "CPCM" in ln or "SMD" in ln]          # the same solvent as the DFT jobs
    assert "--mem-per-cpu=5000MB" in open(os.path.join(sp_r1, "submit.sl")).read()
    assert not os.path.exists(os.path.join(hcn, "calcs/Steps/S1/SP_Single_Point"))      # no TS yet
    for _ in range(4):
        fake_orca.finish_ready_jobs(hcn)
        run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    assert job_state(os.path.join(hcn, "calcs/Steps/S1/SP_Single_Point")) == "finished"
    capsys.readouterr()
    run(hcn, "status", "--no-squeue")
    assert "SINGLE POINTS (DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C TightPNO): finished 3" in capsys.readouterr().out
    import csv
    from mechflow.energies import Energies
    p = Project.load(hcn)
    sp, dft = Energies(p), Energies(p, level="dft")
    assert sp.level == "sp" and dft.level == "dft"
    for sid in ("R1", "P1"):           # E(SP) + (G - E(el)) = (G - 0.1 + shift) + 0.1
        assert abs(sp.g_species(sid) - dft.g_species(sid) - fake_orca.SP_SHIFT) < 1e-7
    assert abs(sp.steps["S1"]["gibbs"] - dft.steps["S1"]["gibbs"] - fake_orca.SP_SHIFT) < 1e-7
    run(hcn, "collect", "--no-squeue")
    assert "Energies: E(DLPNO-CCSD(T)" in capsys.readouterr().out
    rows = list(csv.DictReader(open(os.path.join(hcn, "results", "species.csv"))))
    assert rows[0]["SP"] == "finished" and rows[0]["E_SP_Eh"] and rows[0]["G_used_Eh"]
    run(hcn, "collect", "--no-squeue", "--level", "dft")
    assert "Energies: DFT free energies" in capsys.readouterr().out


def test_single_points_off_by_default(hcn):
    run(hcn, "next", "--no-squeue")
    fake_orca.finish_ready_jobs(hcn)
    run(hcn, "next", "--no-squeue")
    assert not os.path.exists(os.path.join(hcn, "calcs/Species/R1/SP_Single_Point"))
    from mechflow.energies import Energies
    assert Energies(Project.load(hcn)).level == "dft"
