# Running on a Slurm cluster

mechflow is made for clusters that use Slurm (such as Mahuika). It does not submit anything by itself:
`mechflow submit` submits the jobs that `mechflow next` has prepared.

## The loop

```bash
cd ~/projects/my_reaction           # the folder with project.yaml
mechflow next                       # write every input that can be written now
mechflow submit                     # submit them; each job's id goes into its folder (.submitted)
# ... later ...
mechflow status --brief             # what is running, what failed, what is ready
mechflow next                       # next stages for whatever finished
mechflow submit
```

`mechflow submit` runs `sbatch submit.sl` in each prepared job's folder and writes the job id into its
`.submitted` file. That file is made before `sbatch` runs (and removed again if `sbatch` refuses the job), so
a job is never submitted twice, even by two `mechflow submit` commands at once. A job that `sbatch` refuses
is listed with `sbatch`'s message and stays ready to submit. `mechflow submit --dry-run` lists the jobs
without submitting them, and `--limit 50` submits at most 50 jobs, if your cluster limits the number of
queued jobs.

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

## Jobs that start the next ones (chain)

With

```yaml
chain:
  enabled: true
  setup: ["module load Python/3.11.6-foss-2023a", "source ~/mechflow-venv/bin/activate"]
```

every job script ends, once ORCA has finished normally, by running `mechflow chain` in the project. That does what
you would do by hand: `mechflow next` works out which stages can start now (a SCAN once its start structure is
optimised, the OptTS from the SCAN maximum, the IRC once the TS has one imaginary frequency, Step 4B, the single
points) and writes them, and `mechflow submit` submits them. So after the first `mechflow next` and
`mechflow submit`, Step 1 → SCAN → OptTS → IRC → Step 4B follow one another by themselves.

- `setup` holds the shell lines that make the `mechflow` command work inside a job (it runs on a compute node).
- Only one `mechflow chain` works on a project at a time (a lock folder, `.mechflow_chain.lock`), so two jobs that
  finish together cannot write or submit the same stage twice.
- A result that needs you stops the chain at that point: a structure flagged CHECK, a SCAN with no maximum, a TS
  with the wrong number of imaginary frequencies, or a job that failed. `mechflow status` shows these; deal with
  them (`mechflow accept`, `mechflow reset`, ...), then run `mechflow next` and `mechflow submit` (or
  `mechflow chain`) once, and the chain carries on.
- What each chain did is at the end of the job's `slurm-<id>.out`.
- The cluster must let jobs submit jobs (`sbatch` on a compute node), as most do.
- Jobs written before `chain` was switched on do not chain: `mechflow next --refresh` rewrites the ones that have
  not been submitted.

## When jobs fail

`mechflow status` counts the failed jobs and says why each one stopped, when it can tell (time limit, out of
memory, cancelled, an ORCA error). `mechflow reset` sets them up to run again, and `mechflow submit` submits
them:

```bash
mechflow reset --dry-run      # what it would do
mechflow reset                # every failed job
mechflow submit
```

Nothing is deleted: the files of each attempt move into `failed_attempts/<date>/` inside the job's folder.

| stage | the next attempt |
|---|---|
| Step 1, OptTS, Forwards, Backwards | restarts from the last geometry of the attempt (`last_geometry_restart.xyz`) |
| SCAN | a new part from the last converged point, as in the procedure (see [restarting a SCAN](workflow.md#restarting-a-scan)) |
| IRC, single points | runs again from the same structure |

The inputs are written again from the current `project.yaml`, so a longer time limit or more memory set
there (for example in `stage_overrides`) applies to the next attempt. Jobs that are queued or running in
Slurm are never touched. Other options:

| option | what it does |
|---|---|
| `mechflow reset I5a S7` | only these species and steps (S7 includes S7a and S7b) |
| `--stage OptTS` | only these stages: Step1, SP, SCAN, OptTS, IRC, Fwd, Bwd |
| `--from-scratch` | start the stage again from the beginning: `mechflow next` writes it again from the previous stage's result |
| `--redo-done` | also redo jobs that finished (only for jobs you name) |
| `--downstream` | also reset, from scratch, the later stages that used their results |
| `--force` | also reset submitted jobs whose state is unknown because `squeue` could not be asked |

For example, `mechflow reset S5a --stage OptTS --redo-done --from-scratch --downstream` throws away the TS
of S5a and everything after it; `mechflow next` then writes the OptTS again from the SCAN maximum (or
`mechflow next S5a --image 57` from another SCAN point).

## Working on a laptop and the cluster

mechflow runs anywhere with Python. A common pattern:

1. Write `project.yaml` and run `mechflow validate` on your laptop.
2. Copy the project folder to the cluster and run the `next` / `submit` loop there.
3. Copy `results/` back, or run `mechflow map` on your laptop if Graphviz is not on the cluster.

`mechflow submit --local` writes `run_ready.sh` instead, which runs the jobs one after the other on the
current computer. This is useful for the small [HCN example](examples/hcn.md).
