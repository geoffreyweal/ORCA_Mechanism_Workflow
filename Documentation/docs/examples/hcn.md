# Example: HCN → HNC

Folder: [`examples/HCN_isomerisation`](https://github.com/geoffreyweal/ORCA_Mechanism_Workflow/tree/main/examples/HCN_isomerisation)

A 1,2-hydrogen shift turns hydrogen cyanide (R1) into hydrogen isocyanide (P1) through one transition state
(TS1). The example uses B3LYP-D3BJ/def2-SVP in the gas phase on 4 cores, so every job runs in seconds to
minutes, even on a laptop (`mechflow submit --local`).

The [Quick start](../quick_start.md) goes through it command by command.

Things to notice in `project.yaml`:

- Only settings that differ from the defaults are given (`orca.method`, `orca.solvent_block`, `resources`,
  `slurm.time`).
- The step SCANs the H···N distance (atoms 2 and 1) from R1, where they are not bonded, so the SCAN brings
  them together, down to 1.00 Å.
- R1's starting structure is slightly bent. A perfectly linear HCN gives the SCAN no direction in which to
  bend.

![HCN map](../figures/HCN.svg)
