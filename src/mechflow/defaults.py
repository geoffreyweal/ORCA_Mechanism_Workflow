"""Default settings: the ORCA Mechanism Procedure templates, with chloroform as the solvent.

SETTINGS_YAML is written into every new project.yaml (``mechflow init``) and is also
the source of the defaults, so a project only has to list what it changes.
"""
import copy

import yaml

SETTINGS_YAML = '''\
# -------------------------------------------------------------------------------------
# ORCA settings. Every input mechflow writes is built from these lines, which are
# taken from the ORCA Mechanism Procedure templates.
# Decide on them BEFORE the first job runs: every calculation in a project must use
# the same method, basis set and solvent, or the energies cannot be compared.
# Inputs that already exist are never rewritten (see 'mechflow next --refresh').
# -------------------------------------------------------------------------------------
orca:
  method: "!B3LYP DEF2-TZVP D3BJ"
  keywords:
    opt: "!OPT FREQ TightOPT TightSCF defgrid2"           # Step 1 and Step 4B
    scan: "!OPT NormalOPT TightSCF defgrid2 # Try TightOPT if you have convergence problems."   # Step 2Ba
    optts: "!OptTS NumFreq TightOPT TightSCF defgrid2"     # Step 3
    irc: "!IRC TightSCF defgrid2"                           # Step 4A
    neb: "!NEB-CI TightSCF defgrid2"                        # Step 2Bb (with neb.enabled)
    nebend: "!OPT TightOPT TightSCF defgrid2"               # Step 2Bb: an NEB end of more than one molecule
  scf_block: |
    %SCF
        MaxIter 2000       # Here setting MaxIter to a very high number. Intended for systems that require sometimes 1000 iterations before converging (very rare).
        DIISMaxEq 5        # Default value is 5. A value of 15-40 necessary for difficult systems.
        directresetfreq 15 # Default value is 15. A value of 1 (very expensive) is sometimes required. A value between 1 and 15 may be more cost-effective.
    END
  solvent_block: |
    # Solvent: chloroform (CHCl3), dielectric constant 4.81 and refractive index 1.4459 (20 C).
    # Alternatives: "!CPCM(Chloroform)" on the method line, or SMD:  %CPCM SMD true SMDSOLVENT "CHLOROFORM" END
    %CPCM EPSILON 4.81 REFRAC 1.4459 END
  blocks:                 # extra input blocks for one kind of job: opt, scan, optts, irc, neb or nebend
    optts: |
      %GEOM
          Calc_Hess true # Calculate Hessian in the beginning
          NumHess true # Request numerical Hessian (analytical not available)
          Recalc_Hess 1  # Recalculate the Hessian for every step
      END
    irc: |
      %IRC
          MaxIter 2000
          InitHess calc_numfreq
          Direction both
          TolMAXG 2.e-3  # Max gradient (a.u.). Currently set to values equivalent to LooseOpt.
          TolRMSG 5.e-4  # RMS gradient (a.u.). Currently set to values equivalent to LooseOpt.
          PrintLevel 1
      END
  scan_maxiter: 200       # adds "MaxIter 200" to the SCAN %GEOM block (the procedure's fix for SCAN points
                          # that hit the 50-cycle limit); set to null to leave it out

# -------------------------------------------------------------------------------------
# Single points (SP, optional): a higher-level single-point energy on every optimised structure
# (each species of one molecule from Step 1, each transition state from Step 3), e.g. DLPNO-CCSD(T)
# on the DFT geometries. 'mechflow collect' then uses
#     G = E(single point) + [G - E(el)] of the DFT frequency calculation
# (use 'mechflow collect --level dft' for the DFT free energies alone).
# The single points use the scf_block and solvent_block above (the same solvent).
# -------------------------------------------------------------------------------------
single_point:
  enabled: false          # true: the project has single points; 'mechflow sp' runs them when you say so
  auto: false             # true: 'mechflow next' writes them by itself as soon as each structure is ready (and
                          # 'mechflow submit' submits them); false: only 'mechflow sp' writes and submits them
  method: "!DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C TightPNO"
  keywords: "!RIJCOSX def2/J TightSCF defgrid2"
  block: ""               # extra input for the single points, e.g. a %MDCI block
  complexes: false        # true: also species of more than one molecule (complexes, e.g. SCAN starts). They need
                          # none: a state's energy is the sum of its molecules' (a complex that is itself a state of
                          # a mechanism gets one anyway)

# -------------------------------------------------------------------------------------
# What to calculate. skip_steps: steps that 'mechflow next', 'submit', 'sp' and 'reset' (and so the chain)
# leave alone unless you name them, e.g. steps that make or use an intermediate that turned out not to
# exist. 'mechflow status' lists them as skipped; their species are still optimised.
# -------------------------------------------------------------------------------------
calculate:
  skip_steps: []

# -------------------------------------------------------------------------------------
# Imaginary frequencies where there should be none (a species or a Step 4B end: any; a TS: more than
# one). 'mechflow reset --imaginary' restarts such an optimisation from its last geometry pushed along
# the mode(s), read from orca.hess. With auto: true, 'mechflow next' (and so the chain) does that by
# itself, at most max_pushes times per job; after that the job waits for you.
# -------------------------------------------------------------------------------------
imaginary:
  auto: false
  max_pushes: 2
  amplitude: 0.3          # Angstrom: how far the atom that moves most in the mode(s) is pushed

# -------------------------------------------------------------------------------------
# Complexes: species of more than one molecule that are only the start of a SCAN (e.g. a substrate
# hydrogen-bonded to an acid). With optimise: false they get no Step 1 of their own: each SCAN starts
# from the complex built from its molecules' optimised structures (each molecule is a species of its
# own, its atoms listed in the same order), every molecule placed where the complex's xyz file has it.
# The first SCAN point relaxes the complex anyway, and no state's energy needs it (a state is the sum
# of its molecules). A complex that a mechanism has as a state is always optimised.
# -------------------------------------------------------------------------------------
complexes:
  optimise: true          # false: build the complexes from their optimised molecules (saves their Step 1)

# -------------------------------------------------------------------------------------
# Step 2 by NEB-CI instead of the SCAN (the procedure's Step 2Bb). Each step needs its two ends as xyz
# files at 'ends' ({step}: the step's ID; {side}: left or right, the two sides of its reaction), with
# the same atoms in the same order: the optimised structures (Step 1) of the molecules of that side,
# placed where they react (made e.g. with a script, and looked at first). An end of one molecule is
# used as it is; an end of more than one molecule is optimised first (Step2Bb_Ends; if its bonds
# change there, the step waits for you). The NEB-CI then runs between the two ends (Step2Bb_NEB) and
# its climbing image is the TS guess for the OptTS (Step 3).
# -------------------------------------------------------------------------------------
neb:
  enabled: false
  ends: "neb/{step}_{side}.xyz"
  images: 8               # NImages: the images between the two ends (ORCA's default is 8)
  block: |                # more lines for the %NEB block
    Interpolation IDPP

# -------------------------------------------------------------------------------------
# A chemist in the loop (optional). With enabled: true, Step 2 of a step (its SCAN, or its NEB with
# the end optimisations) waits until you approve where it starts from, and the OptTS waits until you
# approve the TS guess (the SCAN maximum or the NEB's climbing image): 'mechflow approve <steps>
# --start' or '--guess' (a file in the step's folder), or list the steps here.
# -------------------------------------------------------------------------------------
approval:
  enabled: false
  start: []
  guess: []

resources:
  nprocs: 32              # cores per job (%PAL NPROCS and --ntasks)
  maxcore: 2000           # MB per core (%maxcore)

slurm:
  orca_module: ORCA/6.1.0-f.0-OpenMPI-4.1.5
  partition: genoa,milan
  time: "3-00:00"
  mem_per_cpu: 2000MB
  module_purge: false     # true: 'module purge' just before the ORCA module is loaded (every job starts clean)
  print_binding: false    # true: before ORCA starts, print the cores each task was given and the node's NUMA
                          # domains (into slurm-<id>.out), to check that a job was not spread over NUMA domains
  extra_sbatch: []        # extra "#SBATCH ..." lines, e.g. ["#SBATCH --account=nesi12345"]
  extra_lines: []         # shell lines run after the module is loaded, before ORCA starts,
                          # e.g. ["export OMP_PROC_BIND=true", "export OMP_PLACES=cores"]

# Different resources for one kind of job (opt = Step 1 and Step 4B, scan, optts, irc,
# sp = the single points), e.g.
#   stage_overrides:
#     optts: {time: "5-00:00"}
#     sp: {maxcore: 4000, mem_per_cpu: 5000MB}     # DLPNO-CCSD(T) needs more memory per core
stage_overrides: {}

# -------------------------------------------------------------------------------------
# Chained jobs (optional). When a job finishes, its Slurm script runs 'mechflow chain' in the
# project: that writes and submits every stage that can start now ('mechflow next' + 'mechflow
# submit'), so Step 1 -> SCAN -> OptTS -> IRC -> Step 4B follow one another by themselves (the
# single points too with single_point.auto; otherwise 'mechflow sp' runs them). A result that
# needs you (CHECK, no SCAN maximum, a failed job) stops the chain there: deal with it, then run
# 'mechflow next' and 'mechflow submit' (or 'mechflow chain').
# -------------------------------------------------------------------------------------
chain:
  enabled: false
  setup: []               # shell lines that make the 'mechflow' command work inside a job, e.g.
                          #   ["module load Python/3.11.6-foss-2023a", "source ~/mechflow-venv/bin/activate"]

scan:
  step: 0.02              # distance between SCAN points in Angstrom (as in the procedure)
  stretch_to: 3.2         # default final distance when a bond is stretched (no 'to' given)
'''

PROJECT_DEFAULTS = {
    "name": None,              # used in Slurm job names (default: the project folder name)
    "title": "",
    "atom_numbering": 0,       # 0: atoms counted from 0 as in ORCA inputs; 1: counted from 1 (Avogadro, GaussView)
    "bond_tolerance": 1.25,    # bonded if distance <= tolerance x sum of covalent radii
}

ENERGY_DEFAULTS = {
    "reference": None,         # e.g. "R1 + R2 + 2 C1"; default: every R species once
    "balance_with": [],        # species added to balance each state, e.g. [P2, C1, C2]
    "temperature": 298.15,     # K, used for product ratios
    "units": "kJ/mol",         # or kcal/mol
    "standard_state": "1 atm", # ORCA's values; "1 M" adds RT ln(24.5) (7.9 kJ/mol at 298 K) per molecule
    "proton_transfers": "auto",  # barrier-free proton transfers in the network search: auto, path or steps
}

#: energies.proton_transfers: which proton transfers the network search may make without a barrier
PROTON_TRANSFER_MODES = ("auto", "path", "steps")

JOB_KINDS = ("opt", "scan", "optts", "irc", "sp", "neb", "nebend")
DFT_KINDS = ("opt", "scan", "optts", "irc", "neb", "nebend")   # the kinds that use orca.method and orca.keywords


def default_settings():
    d = yaml.safe_load(SETTINGS_YAML)
    d["project"] = copy.deepcopy(PROJECT_DEFAULTS)
    d["energies"] = copy.deepcopy(ENERGY_DEFAULTS)
    return d


def deep_merge(base, override):
    """Mappings are merged key by key; anything else in `override` replaces the value in `base`."""
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out
