# Development

## Repository layout

```
src/mechflow/
    cli.py          the commands
    project.py      project.yaml and species/: species, steps, mechanisms, balancing, filters
    workflow.py     the procedure's stages: checks, and the input for the next stage
    orca.py         ORCA inputs, Slurm scripts, reading outputs, job states
    chem.py         xyz files, formulas, bonds, fragments, connectivity fingerprints
    energies.py     Step 5: energies, profiles, product ratios
    network.py      maps (Graphviz and Mermaid)
    report.py       overview.md
    defaults.py     default settings (the procedure's templates, chloroform)
tests/              pytest; fake_orca.py fabricates ORCA outputs, so the tests need no ORCA
examples/           complete example projects
Documentation/      this website (Zensical): zensical.toml and the pages in docs/
```

mechflow depends only on PyYAML. It uses no chemistry toolkit, so it installs anywhere Python does.

## Tests

```bash
pip install -e ".[test]"
pytest
```

The tests take the HCN example through every stage with fabricated outputs (including a SCAN restart and a
choice of SCAN image). They also check that the inputs match the procedure's templates, test the balancing on
a small acid-catalysed reaction, and check that every project in `examples/` loads without problems.

`tests/fake_orca.py` can also be used by hand, to try a project without ORCA:

```bash
cp -r my_reaction /tmp/try && cd /tmp/try
mechflow next
python -c "import sys; sys.path.insert(0, 'PATH/TO/REPO/tests'); import fake_orca; fake_orca.finish_ready_jobs('.')"
# repeat next / fake ...
```

## This website

The website is built with [Zensical](https://zensical.org) from `Documentation/zensical.toml` and the pages
in `Documentation/docs/`:

```bash
pip install -e ".[docs]"         # Zensical needs Python 3.11 or newer
cd Documentation
zensical serve                   # http://localhost:8000, updates as you edit
zensical build --strict          # writes site/; stops on broken links
```

Pushing to `main` builds the site and publishes it to GitHub Pages
(`.github/workflows/build_documentation.yml`); the tests run on every push (`.github/workflows/tests.yml`).
GitHub Pages must publish from GitHub Actions (*Settings → Pages → Build and deployment → Source: GitHub
Actions*).
