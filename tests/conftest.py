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


FAKE_SBATCH = """#!/bin/sh
# A fake sbatch for the tests: prints a new job id as 'sbatch --parsable' does (or fails, if asked to).
if [ -n "${FAKE_SBATCH_FAIL:-}" ]; then
    echo "sbatch: error: Batch job submission failed: Invalid account or account/partition combination specified" >&2
    exit 1
fi
count="${FAKE_SBATCH_DIR}/count"
n=$(( $(cat "$count" 2>/dev/null || echo 1000) + 1 ))
echo "$n" > "$count"
echo "$n;cluster"
"""


@pytest.fixture
def fake_sbatch(tmp_path, monkeypatch):
    """Put a fake sbatch first on PATH (job ids 1001, 1002 ...)."""
    d = tmp_path / "fake_bin"
    d.mkdir()
    exe = d / "sbatch"
    exe.write_text(FAKE_SBATCH)
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", "%s%s%s" % (d, os.pathsep, os.environ.get("PATH", "")))
    monkeypatch.setenv("FAKE_SBATCH_DIR", str(d))
    monkeypatch.delenv("FAKE_SBATCH_FAIL", raising=False)
    return str(d)
