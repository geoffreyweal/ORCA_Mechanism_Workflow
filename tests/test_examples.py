"""Every project in examples/ loads without problems and follows the naming scheme."""
import os

import pytest

from mechflow import network
from mechflow.project import Project

EXAMPLES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")
FOLDERS = sorted(d for d in os.listdir(EXAMPLES) if os.path.isfile(os.path.join(EXAMPLES, d, "project.yaml")))


def test_examples_are_found():
    assert "HCN_isomerisation" in FOLDERS


@pytest.mark.parametrize("folder", FOLDERS)
def test_example_has_no_problems(folder):
    p = Project.load(os.path.join(EXAMPLES, folder))
    assert not p.errors, p.errors[:3]
    assert not p.warnings, p.warnings[:3]
    for sid in p.species:
        assert sid[0] in "RICP"                   # reactants, intermediates, catalysts, products
    for sid, st in p.steps.items():
        assert st.ts == "T" + sid                 # step Sn goes through TSn
    g = network.build(p, list(p.mechanisms.values()))
    assert "__start__" in g.nodes
    assert network.to_dot(p, g, folder).startswith("digraph")
