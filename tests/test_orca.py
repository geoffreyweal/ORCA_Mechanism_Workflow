"""The inputs must have the layout of the ORCA Mechanism Procedure templates."""
import os

from mechflow import orca
from mechflow.defaults import default_settings, deep_merge

from conftest import write

HEAD = """!B3LYP DEF2-TZVP D3BJ
{kw}
%SCF
    MaxIter 2000       # Here setting MaxIter to a very high number. Intended for systems that require sometimes 1000 iterations before converging (very rare).
    DIISMaxEq 5        # Default value is 5. A value of 15-40 necessary for difficult systems.
    directresetfreq 15 # Default value is 15. A value of 1 (very expensive) is sometimes required. A value between 1 and 15 may be more cost-effective.
END
# Solvent: chloroform (CHCl3), dielectric constant 4.81 and refractive index 1.4459 (20 C).
# Alternatives: "!CPCM(Chloroform)" on the method line, or SMD:  %CPCM SMD true SMDSOLVENT "CHLOROFORM" END
%CPCM EPSILON 4.81 REFRAC 1.4459 END
%PAL NPROCS 32 END
%maxcore 2000
"""


def test_opt_input():
    text = orca.render_input(default_settings(), "opt", "reactant.xyz", 1, 1)
    assert text == HEAD.format(kw="!OPT FREQ TightOPT TightSCF defgrid2") + "* xyzfile 1 1 reactant.xyz\n\n\n"


def test_scan_input():
    text = orca.render_input(default_settings(), "scan", "opt_product.xyz", 1, 1, scan=(11, 16, 3.245, 0.745, 126))
    assert text == HEAD.format(kw="!OPT NormalOPT TightSCF defgrid2 # Try TightOPT if you have convergence problems.") + (
        "%GEOM\n    MaxIter 200   # optimisation cycles allowed per SCAN point\n"
        "    SCAN B 11 16 = 3.245, 0.745, 126 END\nEND\n* xyzfile 1 1 opt_product.xyz\n\n\n")


def test_optts_and_irc_inputs():
    s = default_settings()
    assert orca.render_input(s, "optts", "orca.073.xyz", 1, 1) == HEAD.format(
        kw="!OptTS NumFreq TightOPT TightSCF defgrid2") + (
        "%GEOM\n    Calc_Hess true # Calculate Hessian in the beginning\n"
        "    NumHess true # Request numerical Hessian (analytical not available)\n"
        "    Recalc_Hess 1  # Recalculate the Hessian for every step\nEND\n* xyzfile 1 1 orca.073.xyz\n\n\n")
    irc = orca.render_input(s, "irc", "TS.xyz", 1, 1)
    assert irc.startswith(HEAD.format(kw="!IRC TightSCF defgrid2") + "%IRC\n    MaxIter 2000\n    InitHess calc_numfreq\n")
    assert irc.endswith("    PrintLevel 1\nEND\n* xyzfile 1 1 TS.xyz\n\n\n")


def test_settings_overrides_and_stage_overrides():
    s = deep_merge(default_settings(), {"resources": {"nprocs": 8}, "orca": {"solvent_block": ""},
                                        "stage_overrides": {"optts": {"time": "7-00:00", "nprocs": 16}}})
    assert "%PAL NPROCS 8 END" in orca.render_input(s, "opt", "a.xyz", 0, 1)
    assert "CPCM" not in orca.render_input(s, "opt", "a.xyz", 0, 1)
    assert "%PAL NPROCS 16 END" in orca.render_input(s, "optts", "a.xyz", 0, 1)
    sl = orca.render_submit(s, "optts", "job")
    assert "#SBATCH --time=7-00:00     # Walltime" in sl and "#SBATCH --ntasks=16" in sl
    assert "#SBATCH --time=3-00:00" in orca.render_submit(s, "opt", "job")


def test_submit_script_matches_template():
    sl = orca.render_submit(default_settings(), "opt", "A3_Step1_Reactants")
    assert sl == """#!/bin/bash -e
#SBATCH --job-name=A3_Step1_Reactants
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=2000MB
#SBATCH --partition=genoa,milan
#SBATCH --time=3-00:00     # Walltime
#SBATCH --output=slurm-%j.out      # %x and %j are replaced by job name and ID
#SBATCH --error=slurm-%j.err
#SBATCH --nodes=1 # OpenMPI can have problems with ORCA over multiple nodes sometimes, depending on your system.

# Load ORCA
module load ORCA/6.1.0-f.0-OpenMPI-4.1.5

# ORCA under MPI requires that it be called via its full absolute path
orca_exe=$(which orca)

# Don't use "srun" as ORCA does that itself when launching its MPI process.
${orca_exe} orca.inp > output.out
"""


def test_parse_output_and_job_state(tmp_path):
    from fake_orca import out_text
    d = str(tmp_path / "job")
    orca.write_job(d, default_settings(), "optts", "x.xyz", 0, 1, "job")
    assert orca.job_state(d) == "ready"
    write(os.path.join(d, ".submitted"), "12345\n")
    assert orca.job_state(d) == "submitted"
    write(os.path.join(d, "output.out"), out_text(-100.5, imag=-412.3, term=False))
    assert orca.job_state(d) == "running"
    write(os.path.join(d, "slurm-12345.err"), "slurmstepd: error: *** JOB 12345 CANCELLED DUE TO TIME LIMIT ***")
    assert orca.job_state(d) == "failed"
    write(os.path.join(d, "output.out"), out_text(-100.5, imag=-412.3))
    assert orca.job_state(d) == "finished"
    r = orca.parse_output(os.path.join(d, "output.out"))
    assert r["opt_converged"] and r["terminated_normally"]
    assert r["imag"] == [-412.3] and abs(r["gibbs"] + 100.5) < 1e-9


def test_rewrite_job_uses_current_settings(tmp_path):
    d = str(tmp_path / "job")
    orca.write_job(d, default_settings(), "opt", "x.xyz", 0, 1, "job")
    s = deep_merge(default_settings(), {"orca": {"method": "!PBE0 DEF2-SVP D3BJ"}})
    assert orca.rewrite_job(d, s)
    assert open(os.path.join(d, "orca.inp")).read().startswith("!PBE0 DEF2-SVP D3BJ\n")


GENOA = {"slurm": {"partition": "genoa", "module_purge": True, "print_binding": True,
                   "extra_sbatch": ["#SBATCH --cores-per-socket=21", "#SBATCH --threads-per-core=1",
                                    "#SBATCH --distribution=*:block:*", "#SBATCH --mem-bind=local"],
                   "extra_lines": ["export OMP_PROC_BIND=true", "export OMP_PLACES=cores"]}}


def test_submit_script_numa_settings():
    """module_purge and print_binding (off by default), with NeSI's NUMA settings for Genoa nodes."""
    plain = orca.render_submit(default_settings(), "opt", "x")
    assert "module purge" not in plain and "print_binding_stats" not in plain
    sl = orca.render_submit(deep_merge(default_settings(), GENOA), "opt", "x")
    lines = sl.splitlines()
    assert "#SBATCH --partition=genoa" in lines and "#SBATCH --distribution=*:block:*" in lines
    k = lines.index("module purge")
    first = [ln for ln in lines if ln.strip() and not ln.startswith("#")][0]
    assert first == "module purge"                           # the first command after the #SBATCH lines
    assert lines[k + 1:k + 4] == ["", "# Load ORCA", "module load ORCA/6.1.0-f.0-OpenMPI-4.1.5"]
    assert lines[k + 4:k + 6] == ["export OMP_PROC_BIND=true", "export OMP_PLACES=cores"]
    srun = [ln for ln in lines if ln.startswith("srun ")]
    assert srun == [orca.PRINT_BINDING]
    assert lines.index(srun[0]) < lines.index("orca_exe=$(which orca)")
    sbatch = [k for k, ln in enumerate(lines) if ln.startswith("#SBATCH")]
    assert max(sbatch) < lines.index("module purge")       # Slurm reads no #SBATCH line after the first command


def test_print_binding_runs_in_bash(tmp_path):
    """The binding line is valid bash, and prints what it should (srun, lscpu and taskset faked)."""
    import shutil
    import subprocess
    import pytest
    if not shutil.which("bash"):
        pytest.skip("no bash")
    fake = {"srun": '#!/bin/sh\nwhile [ "$1" != "bash" ]; do shift; done\nSLURM_PROCID=0 exec "$@"\n',
            "lscpu": "#!/bin/sh\necho 'NUMA node(s):        8'\n",
            "taskset": "#!/bin/sh\necho \"pid $3's current affinity list: 0-20\"\n",
            "nproc": "#!/bin/sh\necho 21\n"}
    for name, text in fake.items():
        write(str(tmp_path / name), text)
        os.chmod(str(tmp_path / name), 0o755)
    env = dict(os.environ, PATH="%s%s%s" % (tmp_path, os.pathsep, os.environ.get("PATH", "")))
    res = subprocess.run(["bash", "-e", "-c", orca.PRINT_BINDING], capture_output=True, text=True, env=env)
    assert res.returncode == 0, res.stderr
    assert "Task #0 is running on node" in res.stdout
    assert "NUMA node(s):        8" in res.stdout and "Task #0 has 21 CPUs, their core IDs are 0-20" in res.stdout
    env["PATH"] = str(tmp_path / "nothing")                  # no srun at all: the job goes on
    res = subprocess.run(["/bin/bash", "-e", "-c", orca.PRINT_BINDING + "; echo still-running"],
                         capture_output=True, text=True, env=env)
    assert "could not be printed" in res.stdout and "still-running" in res.stdout
