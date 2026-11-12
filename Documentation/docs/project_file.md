# project.yaml reference

A project is a folder with a `project.yaml` file and a `species/` folder. mechflow finds `project.yaml` in
the current folder or the nearest folder above it, so commands work from anywhere inside the project.
`-p FOLDER` points at another project.

```
my_reaction/
    project.yaml
    species/           your xyz files
    calcs/             written by mechflow: one folder per species and per step
    results/           written by 'mechflow collect'
    maps/              written by 'mechflow map'
    overview.md        written by 'mechflow validate'
```

Every setting has a default (the ORCA Mechanism Procedure templates, with chloroform as the solvent), so
a project only needs to list what it changes. `mechflow init` writes all of them out with comments.

## project

```yaml
project:
  name: epoxide_MeOH       # used in Slurm job names (default: the folder name)
  title: Propylene oxide + methanol, acid catalysed   # shown in overview.md, results and maps
  atom_numbering: 0        # 0: atoms counted from 0 as in ORCA inputs; 1: from 1 (Avogadro, GaussView)
  bond_tolerance: 1.25     # two atoms are bonded if closer than 1.25 x the sum of their covalent radii
```

## orca

```yaml
orca:
  method: "!B3LYP DEF2-TZVP D3BJ"
  keywords:
    opt:   "!OPT FREQ TightOPT TightSCF defgrid2"          # Step 1 and Step 4B
    scan:  "!OPT NormalOPT TightSCF defgrid2 # Try TightOPT if you have convergence problems."
    optts: "!OptTS NumFreq TightOPT TightSCF defgrid2"     # Step 3
    irc:   "!IRC TightSCF defgrid2"                         # Step 4A
  scf_block: |
    %SCF
        MaxIter 2000
        DIISMaxEq 5
        directresetfreq 15
    END
  solvent_block: |
    # Solvent: chloroform (CHCl3), dielectric constant 4.81 and refractive index 1.4459 (20 C).
    %CPCM EPSILON 4.81 REFRAC 1.4459 END
  blocks:                  # extra blocks for one kind of job: opt, scan, optts or irc
    optts: |
      %GEOM
          Calc_Hess true
          NumHess true
          Recalc_Hess 1
      END
    irc: |
      %IRC
          MaxIter 2000
          InitHess calc_numfreq
          Direction both
          TolMAXG 2.e-3
          TolRMSG 5.e-4
          PrintLevel 1
      END
  scan_maxiter: 200        # "MaxIter 200" in the SCAN %GEOM block; null to leave it out
```

Every input is: the `method` line, the keyword line for that kind of job, `scf_block`, `solvent_block`,
`%PAL NPROCS`, `%maxcore`, the SCAN block (SCANs only), the `blocks` entry for that kind of job, and
`* xyzfile charge mult file`. That is the layout of the procedure's templates.

!!! warning "Use the same settings for every calculation"
    Energies can only be compared if every job used the same method, basis set and solvent. Decide on
    these before the first job runs. Inputs already written are never changed automatically; use
    `mechflow next --refresh` to rewrite prepared jobs that have not been submitted yet.

!!! note "Comments in ORCA inputs"
    ORCA reads `#` as the start of a comment and a second `#` on the same line as its end, so keep to one
    `#` per line. `mechflow validate` warns about lines with two.

Other solvents: change the dielectric constant and refractive index, or use ORCA's names, e.g.
`solvent_block: "%CPCM SMD true SMDSOLVENT \"WATER\" END"`, or `CPCM(Water)` on the method line with
`solvent_block: ""`.

## single_point

An optional single point (SP) at a higher level on every optimised species and transition state, used by
`collect` as G = E(SP) + [G - E(el)] of the DFT frequency job (see [Energies](energies.md#higher-level-single-points-optional)).

```yaml
single_point:
  enabled: true
  method: "!DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C TightPNO"
  keywords: "!RIJCOSX def2/J TightSCF defgrid2"
  block: ""            # extra input, e.g. a %MDCI block

stage_overrides:
  sp: {maxcore: 4000, mem_per_cpu: 5000MB}     # DLPNO-CCSD(T) needs more memory per core
```

The single points use the same `scf_block` and `solvent_block` as the DFT jobs, so the solvent model is the same.
Their resources are those of the job kind `sp` (`stage_overrides.sp`). Off by default.

## resources, slurm and stage_overrides

```yaml
resources:
  nprocs: 32               # %PAL NPROCS and #SBATCH --ntasks
  maxcore: 2000            # %maxcore (MB per core)

slurm:
  orca_module: ORCA/6.1.0-f.0-OpenMPI-4.1.5
  partition: genoa,milan
  time: "3-00:00"
  mem_per_cpu: 2000MB
  module_purge: true                     # 'module purge' just before the ORCA module is loaded
  print_binding: false                   # print each task's cores and the node's NUMA domains before ORCA starts
  extra_sbatch: ["#SBATCH --account=nesi12345"]
  extra_lines: ["export OMP_PROC_BIND=true"]   # shell lines run after the module is loaded, before ORCA

stage_overrides:                         # different settings for one kind of job
  optts: {time: "7-00:00", nprocs: 64}
  opt:   {time: "1-00:00"}
  sp:    {maxcore: 4000, mem_per_cpu: 5000MB}   # the single points (single_point)
```

Individual species (in `species:`) and steps (`resources:` in a step) can override these too.

## scan

```yaml
scan:
  step: 0.02               # distance between SCAN points, Angstrom
  stretch_to: 3.2          # final distance when a step gives no 'to' and the bond is stretched
```

## species, steps, mechanisms

See [Species](species.md), [Steps](steps.md) and [Mechanisms](mechanisms.md).

```yaml
species:
  R1: hydrogen cyanide
  P1: {name: hydrogen isocyanide, charge: 0, mult: 1}

steps:
  S1:
    name: 1,2-hydrogen shift
    reaction: R1 -> P1
    scan: {from: R1, atoms: [2, 1], to: 1.00}

mechanisms:
  M1:
    name: 1,2-hydrogen shift
    path: [R1, S1, P1]
    product: P1              # optional: product label for comparisons (default: last species of the path)
    # reference: 2 R1        # optional: the mechanism's own reference state (default: energies.reference)
```

A mechanism may also have a `family` (see [Mechanisms](mechanisms.md)). `reference` is for a route that needs
more molecules than `energies.reference`, such as a second molecule of a reactant; see
[A mechanism's own reference](mechanisms.md#a-mechanisms-own-reference).

## energies

```yaml
energies:
  reference: R1 + R2 + C1          # G = 0 (default: every R species once)
  balance_with: [C1, C2]           # species that may be added to balance a state
  temperature: 298.15              # K, for the product ratios
  units: kJ/mol                    # or kcal/mol
  standard_state: 1 atm            # ORCA's values; "1 M" adds RT ln(24.5) = 7.9 kJ/mol per molecule
  proton_transfers: auto           # proton transfers without a barrier: auto, path or steps
```

See [Energies and selectivity](energies.md), and [Proton transfers without a barrier](energies.md#proton-transfers-without-a-barrier)
for `proton_transfers`.
