# Running on a Slurm cluster

mechflow is made for clusters that use Slurm (such as Mahuika). It does not submit anything by itself:
`mechflow submit` writes a script that you run.

## The loop

```bash
cd ~/projects/my_reaction           # the folder with project.yaml
mechflow next                       # write every input that can be written now
mechflow submit                     # write submit_ready.sh
bash submit_ready.sh                # submit; each job's id goes into its folder (.submitted)
# ... later ...
mechflow status --brief             # what is running, what failed, what is ready
mechflow next                       # next stages for whatever finished
mechflow submit && bash submit_ready.sh
```

`submit_ready.sh` holds one line per job:

```bash
( cd "calcs/Steps/S5a/Step3_Opt_TS" && jid=$(sbatch --parsable submit.sl) && echo "$jid" > .submitted && echo "submitted $jid calcs/Steps/S5a/Step3_Opt_TS" )
```

`mechflow submit --run` does the same without the script. `--limit 50` submits at most 50 jobs, if your
cluster limits the number of queued jobs.

## How mechflow knows a job's state

- `.submitted` (written when the job is submitted) and `slurm-<id>.out` give the job id.
- `squeue -u $USER` (asked once per command) says whether that job is still queued or running.
- `output.out` says whether ORCA finished (`ORCA TERMINATED NORMALLY`) or stopped with an error.
- `slurm-<id>.err` and `.out` show time limits, memory problems and cancellations.

A job that left the queue without a finished output is reported as failed.

## Resources

The defaults (32 cores, 2000 MB per core, 3 days, `genoa,milan`) are the procedure's templates. Change
them for the whole project in `resources:` and `slurm:`, or for one kind of job:

```yaml
stage_overrides:
  optts: {time: "7-00:00"}        # OptTS with NumHess and Recalc_Hess 1 is the most expensive stage
  opt:   {time: "1-00:00"}
slurm:
  extra_sbatch: ["#SBATCH --account=nesi12345"]
```

Small species can use fewer cores (`P2: {name: water, nprocs: 4}` in `species:`).

!!! tip "OptTS cost"
    `Recalc_Hess 1` recalculates the numerical Hessian at every OptTS cycle. That is the most expensive
    part of the procedure. `Recalc_Hess 5` in `orca.blocks.optts` is a common compromise.

## Time limits and restarts

- **SCAN:** `mechflow restart-scan S5` (see [the workflow stages](workflow.md#restarting-a-scan)).
- **Other stages:** rename or delete the stage folder and run `mechflow next`. It writes the stage again
  from the previous stage's result. For a longer time limit, set `stage_overrides` first.

## Working on a laptop and the cluster

mechflow runs anywhere with Python. A common pattern:

1. Write `project.yaml` and run `mechflow validate` on your laptop.
2. Copy the project folder to the cluster and run the `next` / `submit` loop there.
3. Copy `results/` back, or run `mechflow map` on your laptop if Graphviz is not on the cluster.

`mechflow submit --local` writes `run_ready.sh` instead, which runs the jobs one after the other on the
current computer. This is useful for the small [HCN example](examples/hcn.md).
