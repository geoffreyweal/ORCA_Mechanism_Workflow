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
  blocks:                 # extra input blocks for one kind of job: opt, scan, optts or irc
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
# (each species from Step 1, each transition state from Step 3), e.g. DLPNO-CCSD(T)
# on the DFT geometries. 'mechflow collect' then uses
#     G = E(single point) + [G - E(el)] of the DFT frequency calculation
# (use 'mechflow collect --level dft' for the DFT free energies alone).
# The single points use the scf_block and solvent_block above (the same solvent).
# -------------------------------------------------------------------------------------
single_point:
  enabled: false          # true: 'mechflow next' also writes the Step 5 single points
  method: "!DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C TightPNO"
  keywords: "!RIJCOSX def2/J TightSCF defgrid2"
  block: ""               # extra input for the single points, e.g. a %MDCI block

resources:
  nprocs: 32              # cores per job (%PAL NPROCS and --ntasks)
  maxcore: 2000           # MB per core (%maxcore)

slurm:
  orca_module: ORCA/6.1.0-f.0-OpenMPI-4.1.5
  partition: genoa,milan
  time: "3-00:00"
  mem_per_cpu: 2000MB
  extra_sbatch: []        # extra "#SBATCH ..." lines, e.g. ["#SBATCH --account=nesi12345"]
  extra_lines: []         # shell lines run before ORCA starts, e.g. ["module purge"]

# Different resources for one kind of job (opt = Step 1 and Step 4B, scan, optts, irc,
# sp = the single points), e.g.
#   stage_overrides:
#     optts: {time: "5-00:00"}
#     sp: {maxcore: 4000, mem_per_cpu: 5000MB}     # DLPNO-CCSD(T) needs more memory per core
stage_overrides: {}

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

JOB_KINDS = ("opt", "scan", "optts", "irc", "sp")
DFT_KINDS = ("opt", "scan", "optts", "irc")       # the kinds that use orca.method and orca.keywords


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
