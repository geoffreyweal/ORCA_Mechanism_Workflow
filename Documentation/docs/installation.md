# Installation

mechflow needs **Python 3.9 or newer** and **PyYAML**. No chemistry toolkit is needed: everything is
done from xyz files.

To run the calculations you also need [ORCA](https://www.faccts.de/orca/) (version 6 is assumed by the
default settings) and, normally, a cluster with the **Slurm** scheduler.

## Install with pip

```bash
pip install git+https://github.com/geoffreyweal/ORCA_Mechanism_Workflow.git
mechflow --version
```

Or from a copy of the repository, so that you can change the code:

```bash
git clone https://github.com/geoffreyweal/ORCA_Mechanism_Workflow.git
cd ORCA_Mechanism_Workflow
pip install -e .
```

## On an HPC cluster (e.g. Mahuika)

Install mechflow into a virtual environment once:

```bash
module load Python                       # any Python 3.9 or newer
python3 -m venv ~/venvs/mechflow
source ~/venvs/mechflow/bin/activate
pip install git+https://github.com/geoffreyweal/ORCA_Mechanism_Workflow.git
```

Then, in every new session:

```bash
module load Python
source ~/venvs/mechflow/bin/activate
```

mechflow only writes files, reads files and (if you ask it to) calls `sbatch` and `squeue`. It runs on the
login node in seconds, even for projects with hundreds of jobs.

## Graphviz (optional, for map images)

`mechflow map` always writes the map as a Graphviz `.dot` file and as Mermaid. To turn the `.dot` file
into SVG, PDF or PNG it needs Graphviz's `dot` program:

| where | how |
|---|---|
| macOS | `brew install graphviz` |
| Linux | `sudo apt install graphviz` (or `conda install graphviz`) |
| a cluster | `module load Graphviz` if your cluster has it |

Without Graphviz you can still view the Mermaid version (`maps/network.md`) on GitHub or in VS Code. You
can also copy `project.yaml` and `species/` to a computer that has Graphviz.

## For development

```bash
pip install -e ".[test,docs]"            # docs: Zensical, which needs Python 3.11 or newer
pytest                                   # the tests
cd Documentation && zensical serve       # this website on http://localhost:8000
```
