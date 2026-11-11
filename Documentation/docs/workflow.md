# The workflow stages

`mechflow next` looks at every species and step, checks the stage that finished last, and writes the
input for the next stage wherever it can. It never overwrites an input that is already there. Running it
again is always safe.

## What `next` checks and writes

| finished stage | checks (from the procedure) | writes |
|---|---|---|
| (nothing yet) | | `calcs/Species/<ID>/Step1_Geo_Opt`: Step 1 optimisation of every species |
| Step 1 | ORCA terminated normally; optimisation converged; no imaginary frequency; bonds as in your xyz file | `calcs/Steps/<ID>/Step2Ba_SCAN` for every step that starts from this species |
| Step 2Ba SCAN | the energy has a maximum inside the SCAN (all SCAN parts combined) | `Step3_Opt_TS` from the highest maximum (`TS_guess_from_SCAN.txt` says which point) |
| Step 3 OptTS | converged; exactly one imaginary frequency | `Step4A_IRC` from the optimised TS |
| Step 4A IRC | both IRC end points written | `Step4B_Validate/Forwards` and `Backwards` |
| Step 4B | both converged with no imaginary frequency; the scanned pair is bonded at one end only; the two ends are the two sides of the step's reaction | the step is DONE (or CHECK) |

The last check compares connectivity. mechflow finds the separate molecules in each end-point structure
and compares their bonds with the species on each side of the reaction. Atom order does not matter;
stereochemistry and conformation are not compared.

## Folder layout

```
calcs/
    Species/R1/Step1_Geo_Opt/         orca.inp  submit.sl  start_structure.xyz  (output.out, orca.xyz ...)
    Steps/S1/
        step_info.txt                 what the step is and how it is scanned
        Step2Ba_SCAN/                 orca.inp  submit.sl  start_opt.xyz  (orca.001.xyz ...)
        Step3_Opt_TS/                 orca.inp  submit.sl  orca.037.xyz  TS_guess_from_SCAN.txt
        Step4A_IRC/                   orca.inp  submit.sl  TS.xyz  (orca_IRC_F.xyz, orca_IRC_B.xyz)
        Step4B_Validate/Forwards/     orca.inp  submit.sl  orca_IRC_F.xyz
        Step4B_Validate/Backwards/    orca.inp  submit.sl  orca_IRC_B.xyz
        SP_Single_Point/              orca.inp  submit.sl  opt_structure.xyz   (single point of the TS)
    Species/R1/SP_Single_Point/       orca.inp  submit.sl  opt_structure.xyz   (single point of the species)
```

The `SP_Single_Point` folders appear only with `single_point.enabled: true`: `next` writes the single point of
a species as soon as its Step 1 is usable, and that of a transition state as soon as its OptTS has exactly one
imaginary frequency (it does not wait for the IRC). `status` counts them under SINGLE POINTS.

Every job folder also has a small `.mechflow.json` (what the job is) and, once submitted, `.submitted`
(the Slurm job id).

## Job states

`mechflow status` reports each job as:

| state | meaning |
|---|---|
| absent | no input yet |
| ready | input written, not submitted |
| submitted | submitted, no output yet |
| running | output being written |
| finished | `ORCA TERMINATED NORMALLY` |
| failed | ORCA error, or Slurm reports a time limit, out of memory, or cancellation, or the job left the queue without finishing |

## When something needs you

These are your decisions, as in the procedure, and mechflow points them out:

| status message | what to do |
|---|---|
| `no energy maximum` | The SCAN is barrierless, or did not go far enough. Look at it with `mechflow scan-profile S5`. Then change `scan.to` or the atoms and run `mechflow reset S5 --stage SCAN --redo-done --from-scratch`, then `next` |
| `the SCAN stopped early` | `mechflow reset S5`: the procedure's restart (see below) |
| the wrong SCAN point chosen | `mechflow next S5 --image 57 --force`. The old OptTS folder is kept as `Step3_Opt_TS_old1` |
| `2 imaginary frequencies` | Look at the modes; restart OptTS from a better point (`--image`) or from the displaced structure |
| `CHECK` after Step 4B | The TS does not connect the two sides of the reaction you wrote. Look at both ends (`viewORCA opt`). If the result is right, run `mechflow accept S5`; if not, the step needs another SCAN coordinate |
| species `bonds changed during the optimisation` | The species is not stable as drawn. Fix the xyz file, or run `mechflow accept I2a` |

## Restarting a SCAN

When a SCAN runs out of time, `mechflow reset S5` (or `mechflow restart-scan S5`) follows the procedure's
restart instructions:

1. The SCAN files move into `Step2Ba_SCAN/SCAN_1`.
2. `SCAN_2` starts from the last converged SCAN point (`last_converged_SCAN_image_restart.xyz`).
3. It keeps the same spacing and final distance.

All `SCAN_n` parts are read together when looking for the maximum, so nothing else changes.

## Changing settings later

- Inputs already written are not changed when you edit `project.yaml`.
- `mechflow next --refresh` rewrites the inputs of jobs that have been prepared but not yet submitted.
- If you edit a species' xyz file before its Step 1 job is submitted, `next` picks up the change. After
  that, `mechflow reset <ID> --redo-done --from-scratch --downstream` starts the species again.
- To redo a stage, `mechflow reset` it (see [when jobs fail](cluster.md#when-jobs-fail)) and run
  `mechflow next`.
