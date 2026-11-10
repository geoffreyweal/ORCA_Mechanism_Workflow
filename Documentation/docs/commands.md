# Commands

Run every command inside a project folder (or any folder below it), or give `-p FOLDER`.
`mechflow COMMAND --help` lists the options of each command.

| command | what it does |
|---|---|
| `mechflow init` | writes `project.yaml` for the xyz files in `species/`, with every setting and commented examples |
| `mechflow validate [-v]` | checks the project and writes `overview.md` |
| `mechflow next` | writes every input that can be written now |
| `mechflow next S5 --image 57 [--force]` | OptTS of step S5 from SCAN point 57 instead of the automatic choice |
| `mechflow next --refresh` | first rewrites prepared, unsubmitted inputs from the current settings |
| `mechflow submit` | writes `submit_ready.sh`, which submits every prepared job and records its job id |
| `mechflow submit --run` | calls `sbatch` straight away |
| `mechflow submit --local` | writes `run_ready.sh` to run the jobs here without Slurm |
| `mechflow submit --limit 50` | at most 50 jobs |
| `mechflow status [--brief]` | the state of every species and step; `--brief` hides finished ones |
| `mechflow check FOLDER` | the procedure's checks on any ORCA job folder |
| `mechflow scan-profile S5` | the SCAN energy profile of a step, with the TS guess marked |
| `mechflow restart-scan S5` | restarts a SCAN that stopped |
| `mechflow accept ID ...` | accepts a species' Step 1 structure, or a step's Step 4B result, that was flagged CHECK |
| `mechflow collect` | Step 5: energies, profiles, product ratios and the best routes through the whole network → `results/` |
| `mechflow map` | network maps → `maps/` |
| `mechflow atoms I4 [i j]` | the numbered atoms of a species, or the distance between two of them |

## Choosing what to work on

`status`, `next`, `submit`, `collect` and `map` take filters:

| filter | selects |
|---|---|
| `--mech M1A,M6B` | these mechanisms (and the species and steps they need) |
| `--mech M7` | the family M7 (M7A and M7B) |
| `--product P1a` | the mechanisms that give P1a |
| `--step S2,S4` | these steps; `S2` also means `S2a` and `S2b` |
| `--step S1-S10` | a range of steps |
| `--species I4a,I4b` | these species (their Step 1 jobs only) |

For example, to run only the steps that decide a selectivity first:

```bash
mechflow next   --step S2,S4,S6,S7
mechflow submit --step S2,S4,S6,S7
```

`--no-squeue` stops `status`, `next` and `submit` from asking Slurm which jobs are queued. Without that
information, a job that disappeared from the queue without an output cannot be detected.

## collect options

| option | meaning |
|---|---|
| `--level dft/sp` | the energies to use: `dft` (G of the DFT frequency jobs) or `sp` (E of the single points + [G - E(el)] of the DFT jobs; the default when `single_point.enabled`); see [Energies](energies.md#higher-level-single-points-optional) |
| `--proton-transfers auto/path/steps` | which proton transfers the network search may make without a barrier, for this run only (default: `energies.proton_transfers`; see [Energies](energies.md#proton-transfers-without-a-barrier)) |

## map options

| option | meaning |
|---|---|
| `--status` | colour by progress: green done, orange CHECK, blue running, grey not run yet, red needs attention |
| `--energies` | add Gibbs energies where they are known |
| `--by-product` | one map per product |
| `--colour-mechanisms` | colour steps used by only one mechanism |
| `--no-names` | IDs only |
| `--direction TB` | top to bottom instead of left to right |
| `--format svg,pdf,png` | Graphviz formats to render (needs Graphviz) |
| `--name`, `--out` | file name and folder |
