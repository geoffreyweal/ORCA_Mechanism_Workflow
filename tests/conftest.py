import os
import shutil
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
EXAMPLES = os.path.join(REPO, "examples")
sys.path.insert(0, HERE)          # for fake_orca


@pytest.fixture
def hcn(tmp_path):
    """A fresh copy of the HCN example project."""
    dst = tmp_path / "HCN"
    shutil.copytree(os.path.join(EXAMPLES, "HCN_isomerisation"), dst)
    return str(dst)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
