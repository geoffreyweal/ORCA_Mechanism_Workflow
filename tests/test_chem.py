import random

import pytest

from mechflow import chem

HCN = chem.Molecule(["C", "N", "H"], [(0, 0, 0), (1.153, 0, 0), (-1.066, 0, 0)])
HNC = chem.Molecule(["C", "N", "H"], [(0, 0, 0), (1.169, 0, 0), (2.163, 0, 0)])
WATER = chem.Molecule(["O", "H", "H"], [(0, 0, 0), (0.96, 0, 0), (-0.24, 0.93, 0)])


def test_formula_and_string():
    assert chem.formula_string(chem.formula_of(["C", "H", "H", "O"])) == "CH2O"
    assert chem.formula_string({"O": 1, "H": 2}) == "H2O"                 # no carbon: alphabetical
    assert chem.formula_string({"C": 12, "H": 10, "O": 4, "P": 1}, -1) == "C12H10O4P-"
    assert chem.formula_string({"C": 1}, 3) == "C(+3)"


def test_symbols_are_normalised():
    assert chem.normalise_symbol("CL") == "Cl"
    assert chem.normalise_symbol("c1") == "C"
    with pytest.raises(chem.XYZError):
        chem.normalise_symbol("Xx")


def test_comment_charge_mult():
    assert chem.comment_charge_mult("I2a oxonium ion  charge=1 mult=1 | key atoms") == (1, 1)
    assert chem.comment_charge_mult("charge = -1") == (-1, None)
    assert chem.comment_charge_mult("nothing here") == (None, None)


def test_multiplicity_check():
    assert chem.multiplicity_problem(["O", "H", "H"], 0, 1) is None
    assert "check the charge" in chem.multiplicity_problem(["O", "H", "H"], 1, 1)     # 9 electrons, singlet
    assert chem.multiplicity_problem(["O", "H", "H"], 1, 2) is None


def test_bonds_and_fragments():
    assert chem.bonds(HCN) == {(0, 1), (0, 2)}
    assert chem.bonds(HNC) == {(0, 1), (1, 2)}
    two = chem.Molecule(WATER.symbols + HCN.symbols, WATER.coords + [(x + 6, y, z) for x, y, z in HCN.coords])
    frags = chem.fragments(len(two), chem.bonds(two))
    assert sorted(len(f) for f in frags) == [3, 3]
    assert sorted(f for _, f in chem.fragment_fingerprints(two)) == ["CHN", "H2O"]


def test_fingerprint_ignores_atom_order_but_sees_connectivity():
    order = [2, 0, 1]
    shuffled = chem.Molecule([HCN.symbols[k] for k in order], [HCN.coords[k] for k in order])
    assert chem.molecule_fingerprint(shuffled) == chem.molecule_fingerprint(HCN)
    assert chem.molecule_fingerprint(HCN) != chem.molecule_fingerprint(HNC)


def test_fingerprint_bigger_molecule_any_order():
    # ethanol vs dimethyl ether: same formula, different connectivity
    etoh = chem.Molecule(["C", "C", "O", "H", "H", "H", "H", "H", "H"],
                         [(0, 0, 0), (1.52, 0, 0), (2.0, 1.35, 0), (-0.36, 1.03, 0), (-0.36, -0.51, 0.89),
                          (-0.36, -0.51, -0.89), (1.88, -0.51, 0.89), (1.88, -0.51, -0.89), (2.96, 1.35, 0)])
    dme = chem.Molecule(["C", "O", "C", "H", "H", "H", "H", "H", "H"],
                        [(0, 0, 0), (1.41, 0, 0), (1.88, 1.33, 0), (-0.36, 1.03, 0), (-0.36, -0.51, 0.89),
                         (-0.36, -0.51, -0.89), (2.97, 1.33, 0), (1.52, 1.84, 0.89), (1.52, 1.84, -0.89)])
    assert chem.formula_of(etoh.symbols) == chem.formula_of(dme.symbols)
    assert chem.molecule_fingerprint(etoh) != chem.molecule_fingerprint(dme)
    rnd = random.Random(3)
    for _ in range(5):
        order = list(range(len(etoh)))
        rnd.shuffle(order)
        perm = chem.Molecule([etoh.symbols[k] for k in order], [etoh.coords[k] for k in order])
        assert chem.molecule_fingerprint(perm) == chem.molecule_fingerprint(etoh)


def test_xyz_round_trip(tmp_path):
    p = tmp_path / "w.xyz"
    chem.write_xyz(str(p), WATER, "water charge=0 mult=1")
    m = chem.read_xyz(str(p))
    assert m.symbols == WATER.symbols and m.comment == "water charge=0 mult=1"
    assert chem.read_xyz_frames(str(p))[0].symbols == WATER.symbols


def test_bad_xyz(tmp_path):
    p = tmp_path / "bad.xyz"
    p.write_text("3\ncomment\nO 0 0 0\nH 1 0\n")
    with pytest.raises(chem.XYZError):
        chem.read_xyz(str(p))
