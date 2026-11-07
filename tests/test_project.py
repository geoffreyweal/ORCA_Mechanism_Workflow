import os
from collections import Counter

import pytest

from mechflow.project import Project, ProjectError, match_any, parse_expression, split_reaction

from conftest import write

# A small acid-catalysed hydration: CH2O + H2O -> CH2(OH)2, catalysed by HF (C1) / F- (C2).
XYZ = {
    "R1": ("CH2O formaldehyde charge=0 mult=1",
           [("C", 0, 0, 0), ("O", 1.21, 0, 0), ("H", -0.55, 0.94, 0), ("H", -0.55, -0.94, 0)]),
    "R2": ("water charge=0 mult=1", [("O", 0, 0, 0), ("H", 0.96, 0, 0), ("H", -0.24, 0.93, 0)]),
    "C1": ("HF charge=0 mult=1", [("F", 0, 0, 0), ("H", 0.92, 0, 0)]),
    "C2": ("fluoride charge=-1 mult=1", [("F", 0, 0, 0)]),
    "I1": ("protonated formaldehyde charge=1 mult=1",
           [("C", 0, 0, 0), ("O", 1.25, 0, 0), ("H", -0.55, 0.94, 0), ("H", -0.55, -0.94, 0), ("H", 1.80, 0.80, 0)]),
    "I2": ("protonated methanediol charge=1 mult=1",
           [("C", 0, 0, 0), ("O", 1.40, 0, 0), ("O", -0.70, 1.30, 0), ("H", -0.40, -0.60, 0.89),
            ("H", -0.40, -0.60, -0.89), ("H", 1.75, 0.90, 0), ("H", -0.30, 2.15, 0.30), ("H", -1.65, 1.30, -0.20)]),
    "P1": ("methanediol charge=0 mult=1",
           [("C", 0, 0, 0), ("O", 1.40, 0, 0), ("O", -0.70, 1.21, 0), ("H", -0.40, -0.60, 0.89),
            ("H", -0.40, -0.60, -0.89), ("H", 1.75, 0.90, 0), ("H", -0.30, 2.05, 0.30)]),
}

PROJECT = """
project: {name: hydration}
species:
  R1: formaldehyde
  R2: water
  C1: HF
  C2: fluoride
  I1: protonated formaldehyde
  I2: protonated methanediol
  P1: methanediol
steps:
  S1:
    name: water adds to the carbon
    reaction: I1 + R2 -> I2
    scan: {from: I2, atoms: [0, 2], to: 3.2}
mechanisms:
  M1:
    path: [R1 + R2, I1, S1, I2, P1]
energies:
  reference: R1 + R2 + C1
  balance_with: [C1, C2]
"""


def make(tmp_path, text=PROJECT, skip=()):
    root = str(tmp_path / "proj")
    for sid, (comment, atoms) in XYZ.items():
        if sid in skip:
            continue
        lines = ["%d" % len(atoms), comment] + ["%s %.4f %.4f %.4f" % a for a in atoms]
        write(os.path.join(root, "species", sid + ".xyz"), "\n".join(lines) + "\n")
    write(os.path.join(root, "project.yaml"), text)
    return root


def test_expressions():
    assert parse_expression("R1 + R2 + 2 C1") == Counter({"R1": 1, "R2": 1, "C1": 2})
    assert parse_expression("2C1") == Counter({"C1": 2})
    left, right = split_reaction("I2 + R2 -> I3a")
    assert left == Counter({"I2": 1, "R2": 1}) and right == Counter({"I3a": 1})
    with pytest.raises(ProjectError):
        split_reaction("I1 I2")


def test_match_any():
    assert match_any("S5a", "S5") and match_any("S5", "S5") and not match_any("S50", "S5")
    assert match_any("S7b", "S1-S9") and not match_any("S10a", "S1-S9")
    assert match_any("M10A", "M10") and not match_any("M10A", "M1")
    assert match_any("I4b", "I2,I4b")


def test_hcn_example_loads():
    p = Project.load(os.path.join(os.path.dirname(os.path.dirname(__file__)), "examples", "HCN_isomerisation"))
    assert not p.errors and not p.warnings
    st = p.steps["S1"]
    assert st.ts == "TS1" and st.mode == "contract" and st.atoms == (2, 1)
    assert abs(st.d0 - 2.2106) < 1e-3
    assert [s.kind for s in p.mechanisms["M1"].states] == ["species", "ts", "species"]


def test_balancing_and_states(tmp_path):
    p = Project.load(make(tmp_path))
    assert not p.errors, p.errors
    st = p.steps["S1"]
    assert st.mode == "stretch" and st.start == "I2" and st.charge == 1
    states = p.mechanisms["M1"].states
    assert [(s.kind, s.label) for s in states] == [("species", "R1 + R2"), ("species", "I1"), ("ts", "TS1"),
                                                   ("species", "I2"), ("species", "P1")]
    assert states[0].full == Counter({"R1": 1, "R2": 1, "C1": 1})      # the reference, labelled as written
    assert states[1].added == Counter({"R2": 1, "C2": 1})        # R1 + C1 -> I1 + C2, water not used yet
    assert states[2].added == Counter({"C2": 1})                  # the TS has the atoms of I2
    assert states[4].added == Counter({"C1": 1})
    assert p.mechanisms["M1"].product == "P1"
    assert p.balance(Counter({"P1": 1})) == Counter({"C1": 1})
    with pytest.raises(ProjectError):
        p.balance(Counter({"I1": 2}))                             # more atoms than the reference


def test_step_direction_from_the_next_state(tmp_path):
    # P1 + I2 holds both sides of S2 (P1 + C1 -> I2 + C2). The state written next, 2 P1, shows that S2 runs
    # backwards here (I2 gives its proton to F-); forwards would protonate the second P1 instead.
    text = PROJECT.replace("  P1: methanediol", "  P1: methanediol\n  I3: methanediol HF complex").replace(
        "mechanisms:",
        "  S2:\n    name: HF protonates methanediol\n    reaction: P1 + C1 -> I2 + C2\n"
        "    scan: {from: I3, atoms: [1, 8], to: 1.0}\nmechanisms:").replace(
        "energies:",
        "  M2:\n    reference: 2 R1 + 2 R2 + 2 C1\n    path: [2 R1 + 2 R2, P1 + I2, S2, 2 P1]\nenergies:")
    hb = XYZ["P1"][1] + [("F", 3.787, -1.193, 0), ("H", 2.965, -0.782, 0)]     # O-H...F 1.75 A from O1

    def load(where, text):
        root = make(where, text)
        lines = ["%d" % len(hb), "methanediol HF complex charge=0 mult=1"] + ["%s %.4f %.4f %.4f" % a for a in hb]
        write(os.path.join(root, "species", "I3.xyz"), "\n".join(lines) + "\n")
        return Project.load(root)

    p = load(tmp_path, text)
    assert not p.errors, p.errors
    states = p.mechanisms["M2"].states
    assert [s.kind for s in states] == ["species", "species", "ts", "species"]
    assert states[1].species == Counter({"P1": 1, "I2": 1}) and states[2].step == "S2"
    assert states[3].species == Counter({"P1": 2})
    text = text.replace("path: [2 R1 + 2 R2, P1 + I2, S2, 2 P1]", "path: [2 R1 + 2 R2, P1 + I2, S2, 2 I2]")
    states = load(tmp_path / "fwd", text).mechanisms["M2"].states
    assert [s.kind for s in states] == ["species", "species", "ts", "species"]
    assert states[3].species == Counter({"I2": 2})


def test_errors_are_reported(tmp_path):
    bad = PROJECT.replace("reaction: I1 + R2 -> I2", "reaction: I1 -> I2")
    p = Project.load(make(tmp_path, bad))
    assert any("two sides" in e for e in p.errors)
    bad = PROJECT.replace("to: 3.2", "to: 1.0")
    p = Project.load(make(tmp_path / "b", bad))
    assert any("not longer" in e for e in p.errors)
    p = Project.load(make(tmp_path / "c", skip=("P1",)))
    assert any("P1" in e and "xyz" in e for e in p.errors)


def test_proton_transfers_setting(tmp_path):
    add = "balance_with: [C1, C2]\n  proton_transfers: %s"
    p = Project.load(make(tmp_path, PROJECT.replace("balance_with: [C1, C2]", add % "steps")))
    assert not p.errors and p.proton_transfers == "steps"
    p = Project.load(make(tmp_path / "b", PROJECT.replace("balance_with: [C1, C2]", add % "never")))
    assert any("energies.proton_transfers: 'never' is not one of auto, path, steps" in e for e in p.errors)
    assert Project.load(make(tmp_path / "c")).proton_transfers == "auto"


def test_charge_conflict_and_multiplicity(tmp_path):
    p = Project.load(make(tmp_path, PROJECT.replace("C2: fluoride", "C2: {name: fluoride, charge: 0}")))
    assert any("charge is 0 in project.yaml but -1" in e for e in p.errors)


def test_select(tmp_path):
    p = Project.load(make(tmp_path))
    species, steps, mechs = p.select(mech="M1")
    assert [s.id for s in steps] == ["S1"] and {s.id for s in species} >= {"I1", "I2", "R2", "C1", "C2", "P1"}
    species, steps, mechs = p.select(species="I1")
    assert [s.id for s in species] == ["I1"] and steps == []
