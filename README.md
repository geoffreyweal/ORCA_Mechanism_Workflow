# mechflow: the ORCA Mechanism Procedure for whole reaction networks

**mechflow** runs the [ORCA Mechanism Procedure](https://geoffreyweal.github.io/ORCA_Mechanism_Procedure/)
(SCAN → OptTS → IRC → validation → energies) for every step of a reaction network, starting from a folder
of xyz files.

**Documentation: <https://geoffreyweal.github.io/ORCA_Mechanism_Workflow/>**

You provide:

- `species/`: one xyz file per species. The file name is the ID: **R1, R2 …** reactants, **I1 …**
  intermediates, **P1 …** products, **C1 …** catalysts. Charge and multiplicity go on the comment line
  (`charge=1 mult=1`).
- `project.yaml`: the steps (**S1, S2 …**, each with transition state **TS1, TS2 …**), which two atoms to
  SCAN for each step, the mechanisms to compare, and the ORCA and Slurm settings (default: the procedure's
  templates, solvent chloroform).

mechflow writes every ORCA input and Slurm script, runs the procedure's checks on each output, and writes
the input of the next stage. It also checks that every Step 4B end point is the species you said the step
connects. Finally it collects energy profiles and product ratios, and draws network maps.

```bash
pip install git+https://github.com/geoffreyweal/ORCA_Mechanism_Workflow.git

cd my_reaction          # contains species/*.xyz
mechflow init           # writes project.yaml; add your steps and mechanisms
mechflow validate       # checks everything, writes overview.md
mechflow next           # writes every input that can be written now
mechflow submit         # submits every prepared job (--dry-run: only list them)
mechflow reset          # failed jobs: set them up to run again (old files kept in failed_attempts/)
# ... repeat  next / submit  until 'mechflow status' says DONE ...
mechflow collect        # energies, profiles, product ratios -> results/
mechflow map --status   # network map -> maps/
```

## Examples

| folder | what |
|---|---|
| [`examples/HCN_isomerisation`](examples/HCN_isomerisation) | HCN → HNC: one step, runs on a laptop in minutes |

## Requirements

Python ≥ 3.9 and PyYAML. ORCA and Slurm are needed to run the calculations. Graphviz is optional, for map
images.

## Licence

GNU General Public License v3.0 or later. See [LICENSE](LICENSE).
