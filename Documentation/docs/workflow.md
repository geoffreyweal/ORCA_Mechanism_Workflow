# The workflow stages

`mechflow next` looks at every species and step, checks the stage that finished last, and writes the
input for the next stage wherever it can. It never overwrites an input that is already there. Running it
again is always safe.

## What `next` checks and writes

| finished stage | checks (from the procedure) | writes |
|---|---|---|
| (nothing yet) | | `calcs/Species/<ID>/Step1_Geo_Opt`: Step 1 optimisation of every species (with `complexes.optimise: false`, not of the [complexes](project_file.md#complexes) built from their molecules) |
| Step 1 | ORCA terminated normally; optimisation converged; no imaginary frequency; bonds as in your xyz file | `calcs/Steps/<ID>/Step2Ba_SCAN` for every step that starts from this species (or from a complex built from it and the other molecules, once they are all optimised) |
| Step 2Ba SCAN | the energy has a maximum inside the SCAN (all SCAN parts combined) | `Step3_Opt_TS` from the highest maximum (`TS_guess_from_SCAN.txt` says which point) |
| Step 3 OptTS | converged; exactly one imaginary frequency | `Step4A_IRC` from the optimised TS |
| Step 4A IRC | both IRC end points written | `Step4B_Validate/Forwards` and `Backwards` (with `complexes.optimise: false`, only for an end of one molecule: an end of several molecules is checked as it is) |
| Step 4B | both converged with no imaginary frequency; the scanned pair is bonded at one end only; the two ends are the two sides of the step's reaction | the step is DONE (or CHECK) |

With [`neb.enabled`](project_file.md#neb), Step 2 is an NEB-CI (the procedure's Step 2Bb) instead of the SCAN:

| finished stage | checks | writes |
|---|---|---|
| Step 1 | as above, for every molecule on both sides of the step; the step's two end files (`neb.ends`) have the same atoms in the same order, and each is the molecules of its side | `Step2Bb_Ends/Left` and `Right` for an end of more than one molecule (an end of one molecule is used as it is) |
| Step 2Bb ends | converged; bonds as in the end file (if a proton went back, or a molecule added again, the step stops: CHECK) | `Step2Bb_NEB`: the NEB-CI between the two ends |
| Step 2Bb NEB | ORCA's `THE NEB OPTIMIZATION HAS CONVERGED` and a climbing image (`orca_NEB-CI_converged.xyz`) | `Step3_Opt_TS` from the climbing image (`TS_guess_from_NEB.txt`), then as above |

With [`approval.enabled`](project_file.md#approval) (a chemist in the loop), two of these wait for you. Step 2 (the
SCAN, or the NEB with its end optimisations) waits until you approve where it starts from, and the OptTS waits until
you approve the TS guess (the SCAN maximum, or the NEB's climbing image). `mechflow status` shows those steps as
`approve`, with what to look at. After `mechflow approve S5 --start` (or `--guess`), `next` and the chain go on as
usual. Approving a step's NEB start holds only while its end files stay the same.

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
        Step2Bb_Ends/Left/, Right/    orca.inp  submit.sl  left_end.xyz     (neb.enabled: an end of several molecules)
        Step2Bb_NEB/                  orca.inp  submit.sl  left_end.xyz  right_end.xyz  (orca_MEP_trj.xyz,
                                      orca_NEB-CI_converged.xyz, orca.interp.final ...)
        approved_start.txt            (approval.enabled) written by 'mechflow approve S1 --start'
        approved_guess.txt            (approval.enabled) written by 'mechflow approve S1 --guess'
        Step3_Opt_TS/                 orca.inp  submit.sl  orca.037.xyz  TS_guess_from_SCAN.txt
        Step4A_IRC/                   orca.inp  submit.sl  TS.xyz  (orca_IRC_F.xyz, orca_IRC_B.xyz)
        Step4B_Validate/Forwards/     orca.inp  submit.sl  orca_IRC_F.xyz
        Step4B_Validate/Backwards/    orca.inp  submit.sl  orca_IRC_B.xyz
        SP_Single_Point/              orca.inp  submit.sl  opt_structure.xyz   (single point of the TS)
    Species/R1/SP_Single_Point/       orca.inp  submit.sl  opt_structure.xyz   (single point of the species)
```

The `SP_Single_Point` folders appear only with `single_point.enabled: true`, when `mechflow sp` writes them (a
species once its Step 1 is usable, a transition state once its step is DONE), or, with `single_point.auto: true`,
when `next` does (a TS as soon as its OptTS has exactly one imaginary frequency). `status` counts them under SINGLE
POINTS.

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
| species or Step 4B: `imaginary frequency` | A saddle point, not a minimum (often a methyl or phenyl group stuck eclipsed). `mechflow reset --imaginary` restarts every such optimisation from its last geometry pushed along the imaginary mode(s) (0.3 Å for the atom that moves most; `--amplitude` changes it), then `mechflow submit`. With [`imaginary.auto: true`](project_file.md#imaginary), `mechflow next` does this by itself (twice at most) |
| OptTS: `2 imaginary frequencies` | Look at the modes. `mechflow reset S5 --imaginary` restarts the OptTS pushed along the extra mode(s) (every one but the largest, the reaction's); or restart it from a better SCAN point (`--image`) |
| `CHECK` after Step 4B | The TS does not connect the two sides of the reaction you wrote. Look at both ends (`viewORCA opt`). If the result is right, run `mechflow accept S5`; if not, the step needs another SCAN coordinate |
| `approve` (with `approval.enabled`) | Look at the step: where Step 2 starts (the SCAN start, or the two NEB ends and the path between them), or the TS guess (the SCAN profile, or the NEB path: `Step2Bb_NEB/orca_MEP_trj.xyz` and `orca.interp.final`). Then `mechflow approve S5 --start` or `--guess`, and `mechflow next` |
| NEB: `the right end changed during its optimisation` | The end of several molecules is not stable as placed (often a proton went back, or a leaving group added again). Look at `Step2Bb_Ends/Right/orca.xyz`. If it is right, `mechflow accept S5`; if not, place the molecules differently in the end file and `mechflow reset S5 --stage Ends --redo-done --from-scratch` |
| NEB: `did not converge` | `mechflow reset S5 --stage NEB` runs it on from its last path (`orca_MEP.allxyz`, as `Restart_ALLXYZFile`) |
| NEB: the climbing image is not the TS you want | `mechflow next S5 --image 12 --force` starts the OptTS from image 12 of the NEB path instead |
| species `bonds changed during the optimisation` | The species is not stable as drawn. Often a proton was hydrogen-bonded to a better base in the start structure and moved over: turn it away from that atom (or draw another conformer) in the xyz file and run `mechflow reset I2a --redo-done --from-scratch --downstream`. If the new structure is what you want, run `mechflow accept I2a` |

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
