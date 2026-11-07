"""Best routes through the whole network (mechflow.routes), with made-up energies.

The project is a small network of C3H6O isomers (and a dimer), with HF / F- as the acid and its base:

    M1:  R1 -S1-> I1 = I2 -S2-> I3                    (I2 is I1 protonated on O; product P1 given by hand)
    M2:  R1 -S4-> I4 -S5-> I1 -S6-> I5 -S7-> P1       (S5 is written I1 -> I4, so M2 uses it backwards)
    M3:  2 R1 -S1-> I1 -S4-> I1 + I4 -S10-> P3        (its own reference: two R1; I1 is carried through S4)
    S8:  P2 -> I5                                     (in no mechanism)

I3 is P1 protonated on O, so I3 + F- and P1 + HF are joined by an automatic proton transfer.
"""
import csv
import os
from collections import Counter

from mechflow import energies, routes
from mechflow.energies import R_KJ, Energies
from mechflow.project import Project
from mechflow.routes import Network, proton_fingerprint

from conftest import write


def chain(comment, heavy, hs):
    """A chain of heavy atoms 1.45 A apart along x, with hs[i] H atoms 1.0 A from atom i."""
    atoms = [(el, 1.45 * i, 0.0, 0.0) for i, el in enumerate(heavy)]
    for i, n in enumerate(hs):
        atoms += [("H", 1.45 * i + dx, dy, dz) for dx, dy, dz in ((0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))[:n]]
    return comment, atoms


# H atoms on each atom of the chain C-C-C-O: the C-H counts tell the isomers apart
XYZ = {
    "R1": chain("charge=0", "CCCO", [3, 2, 1, 0]),
    "C1": chain("HF charge=0", "F", [1]),
    "C2": chain("fluoride charge=-1", "F", [0]),
    "I1": chain("charge=0", "CCCO", [3, 1, 1, 1]),
    "I2": chain("I1 protonated on O charge=1", "CCCO", [3, 1, 1, 2]),
    "I3": chain("P1 protonated on O charge=1", "CCCO", [2, 2, 2, 1]),
    "I4": chain("charge=0", "CCCO", [3, 2, 0, 1]),
    "I5": chain("charge=0", "CCCO", [2, 1, 2, 1]),
    "P1": chain("charge=0", "CCCO", [2, 2, 2, 0]),
    "P2": chain("charge=0", "CCCO", [1, 2, 2, 1]),
    "P3": chain("dimer charge=0", "CCCOCCCO", [3, 2, 1, 0, 3, 2, 1, 0]),
}

PROJECT = """
project: {name: network}
species: {R1: reactant, C1: HF, C2: fluoride, I1: enol, I2: protonated I1, I3: protonated P1, I4: isomer,
          I5: isomer, P1: first product, P2: second product, P3: dimer}
steps:
  S1: {reaction: R1 -> I1, scan: {from: R1, atoms: [0, 1], to: 3.0}}
  S2: {reaction: I2 -> I3, scan: {from: I2, atoms: [0, 1], to: 3.0}}
  S4: {reaction: R1 -> I4, scan: {from: R1, atoms: [0, 1], to: 3.0}}
  S5: {reaction: I1 -> I4, scan: {from: I1, atoms: [0, 1], to: 3.0}}
  S6: {reaction: I1 -> I5, scan: {from: I1, atoms: [0, 1], to: 3.0}}
  S7: {reaction: I5 -> P1, scan: {from: I5, atoms: [0, 1], to: 3.0}}
  S8: {reaction: P2 -> I5, scan: {from: P2, atoms: [0, 1], to: 3.0}}
  S10: {reaction: I1 + I4 -> P3, scan: {from: P3, atoms: [3, 4], to: 3.0}}
mechanisms:
  M1: {path: [R1, S1, I1, I2, S2, I3], product: P1}
  M2: {path: [R1, S4, I4, S5, I1, S6, I5, S7, P1]}
  M3: {reference: 2 R1 + 2 C1, path: [2 R1, S1, I1, S4, I1 + I4, S10, P3]}
energies:
  reference: R1 + C1
  balance_with: [C1, C2]
"""

# Gibbs energies in kJ/mol (R1, HF and F- at 0, so a state's G is the sum over its core species)
G = {"R1": 0.0, "C1": 0.0, "C2": 0.0, "I1": 10.0, "I2": 15.0, "I3": -5.0, "I4": 5.0, "I5": 0.0,
     "P1": -20.0, "P2": -30.0, "P3": -40.0}
TS = {"S1": 100.0, "S2": 40.0, "S4": 50.0, "S5": 45.0, "S6": 120.0, "S7": 30.0, "S8": 60.0, "S10": 70.0}


def load(tmp_path):
    root = str(tmp_path / "proj")
    for sid, (comment, atoms) in XYZ.items():
        lines = ["%d" % len(atoms), comment] + ["%s %.4f %.4f %.4f" % a for a in atoms]
        write(os.path.join(root, "species", sid + ".xyz"), "\n".join(lines) + "\n")
    write(os.path.join(root, "project.yaml"), PROJECT)
    p = Project.load(root)
    assert not p.errors and not p.warnings, p.errors + p.warnings
    return p


def fake_energies(p, species=None, ts=None):
    """An Energies object with the energies above (None = not computed) instead of ORCA outputs."""
    gs, gts = dict(G, **(species or {})), dict(TS, **(ts or {}))
    en = Energies.__new__(Energies)
    en.project, en.units, en.factor, en.T = p, "kJ/mol", 1.0, 298.15     # factor 1: G is given in kJ/mol
    en.RT, en.standard_state, en.per_molecule = R_KJ * en.T, "1 atm", 0.0
    en.species = {s: {"state": "finished", "ok": gs.get(s) is not None, "gibbs": gs.get(s), "e_el": None,
                      "note": ""} for s in p.species}
    en.steps = {s: {"stage": "DONE", "detail": "", "gibbs": gts.get(s), "imag": -300.0, "validated": True}
                for s in p.steps}
    return en


def by_label(net):
    return {n.label: n for n in net.nodes.values()}


def test_proton_fingerprints(tmp_path):
    p = load(tmp_path)
    fp = {s: proton_fingerprint(sp.mol) for s, sp in p.species.items()}
    assert fp["C1"] == fp["C2"]                               # HA and A-
    assert fp["I1"] == fp["I2"] and fp["P1"] == fp["I3"]      # the same molecule with one more H on O
    assert len({fp[s] for s in ("R1", "I1", "I4", "I5", "P1", "P2")}) == 6     # C-H counts differ


def test_states_and_automatic_proton_transfer(tmp_path):
    p = load(tmp_path)
    net = Network(p, fake_energies(p))
    nodes = by_label(net)
    assert set(nodes) == {"R1 + C1", "I1", "I2", "I3", "I4", "I5", "P1", "P2", "I1 + I4", "P3"}
    assert nodes["R1 + C1"].energy == 0.0 and nodes["I3"].full == Counter({"I3": 1, "C2": 1})
    pts = {frozenset((net.nodes[e.a].label, net.nodes[e.b].label)): e for e in net.edges if e.kind == "pt"}
    assert set(pts) == {frozenset(("I1", "I2")), frozenset(("I3", "P1"))}
    assert pts[frozenset(("I3", "P1"))].auto                  # I3 + F- and P1 + HF: found automatically
    assert not pts[frozenset(("I1", "I2"))].auto              # written in M1
    assert list(net.products) == ["P1", "P2", "P3"]


def test_best_route_combines_steps_of_two_mechanisms(tmp_path):
    p = load(tmp_path)
    res = Network(p, fake_energies(p)).routes_to("P1")
    r = res.best
    # M1 alone climbs to TS1 (100) and M2 alone to TS6 (120); M2's start and M1's end only reach TS4 (50)
    assert r.tokens == ["R1 + C1", "S4", "I4", "S5", "I1", "I2", "S2", "I3", "P1"]
    assert r.highest == ("TS4", 50.0) and r.eff == ("TS4", 50.0, "R1 + C1")
    assert r.mechanisms == "M2 (S4, S5) + M1 (S2)"
    assert r.auto_pt == ["I3 <=> P1"]
    assert [(a.highest, a.mechanisms) for a in res.alternatives] == [
        (("TS1", 100.0), "M1 (S1, S2)"), (("TS6", 120.0), "M2 (S4, S5, S6, S7)")]
    assert [x.highest[0] for x in res.routes] == ["TS4", "TS1", "TS6"]


def test_route_uses_steps_backwards(tmp_path):
    p = load(tmp_path)
    r = Network(p, fake_energies(p)).routes_to("P2").best
    # P2 is reached from P1 through S7 (M2 goes I5 -> P1) and S8 (written P2 -> I5)
    assert r.tokens == ["R1 + C1", "S4", "I4", "S5", "I1", "I2", "S2", "I3", "P1", "S7", "I5", "S8", "P2"]
    assert r.highest == ("TS8", 60.0)
    assert r.backwards == ["S5", "S7", "S8"]
    assert r.mechanisms == "M2 (S4, S5) + M1 (S2) + M2 (S7) + no mechanism (S8)"


def test_lowest_highest_point_wins(tmp_path):
    p = load(tmp_path)

    def best(**ts):
        return Network(p, fake_energies(p, ts=ts)).routes_to("P1").best
    # one step (S1) against two (S4 at 50, then S5 at 45): the lower highest point wins, however many steps
    assert best(S1=55.0).tokens[:3] == ["R1 + C1", "S4", "I4"]
    assert best(S1=48.0).tokens[:3] == ["R1 + C1", "S1", "I1"]
    # the same highest point (50): the next-highest point decides (TS2 at 40 after S1, TS5 at 45 after S4)
    r = best(S1=50.0)
    assert r.tokens[:3] == ["R1 + C1", "S1", "I1"] and r.highest == ("TS1", 50.0)


def test_unknown_energies_block_routes(tmp_path):
    p = load(tmp_path)
    net = Network(p, fake_energies(p, ts={"S4": None}))           # TS4 not converged
    assert net.routes_to("P1").best.highest == ("TS1", 100.0)
    assert net.routes_to("P3").best is None                       # I1 + I4 is only made through S4
    net = Network(p, fake_energies(p, species={"I5": None}))      # I5 not computed
    assert net.routes_to("P2").best is None and net.routes_to("P1").best is not None


def test_own_reference_and_spectator(tmp_path):
    p = load(tmp_path)
    ts4 = [s for s in p.mechanisms["M3"].states if s.step == "S4"][0]
    assert ts4.species == Counter({"R1": 1, "I1": 1}) and ts4.added == Counter({"C1": 2})   # I1 rides along
    net = Network(p, fake_energies(p))
    nodes = by_label(net)
    assert nodes["I1 + I4"].energy == 15.0 and nodes["P3"].energy == -40.0     # relative to 2 R1 + 2 C1
    assert [e.energy for e in net.edges if e.label == "TS4 + I1"] == [60.0]    # G(TS4) + G(I1)
    r = net.routes_to("P3").best
    assert r.tokens == ["R1 + C1", "S4", "I4", "S5", "I1", "S4", "I1 + I4", "S10", "P3"]
    assert r.highest == ("TS10", 70.0) and r.reference == "2 R1 + 2 C1"
    assert r.mechanisms == "M2 (S4, S5) + M3 (S4, S10)"


def test_collect_writes_network_routes(tmp_path, monkeypatch):
    p = load(tmp_path)
    monkeypatch.setattr(energies, "Energies", lambda project, queue=None: fake_energies(project))
    text = energies.collect(p)
    assert "## Mechanisms" in text and "## Products" in text
    part = text.split(routes.HEADING)[1]
    assert "| P1 | TS4 | 50.0 | 50.0 (TS4 from R1 + C1) | M2 (S4, S5) + M1 (S2) |" in part
    assert ("- P1: `[R1 + C1, S4, I4, S5, I1, I2, S2, I3, P1]` (S5 used backwards; automatic proton transfer "
            "I3 <=> P1)") in part
    assert "- P1 : P2 cannot be estimated" in part                  # P2 is made from P1
    assert "- P1 : P3 = 1.8e+05 : 1  (TS2 vs TS10: difference 30.0 kJ/mol)" in part
    assert "| P1 | TS1 | 100.0 | +50.0 | M1 (S1, S2) |" in part
    with open(os.path.join(p.root, "results", "network_routes.csv")) as fh:
        rows = list(csv.DictReader(fh))
    assert [(r["product"], r["rank"], r["bottleneck_TS"]) for r in rows if r["product"] == "P1"] == [
        ("P1", "1", "TS4"), ("P1", "2", "TS1"), ("P1", "3", "TS6")]
    assert rows[0]["route"] == "R1 + C1, S4, I4, S5, I1, I2, S2, I3, P1"


def test_collect_with_nothing_computed(tmp_path, monkeypatch):
    p = load(tmp_path)
    nothing = dict(species={s: None for s in G}, ts={s: None for s in TS})
    monkeypatch.setattr(energies, "Energies", lambda project, queue=None: fake_energies(project, **nothing))
    text = energies.collect(p)
    assert text.split(routes.HEADING)[1].strip().splitlines()[0].startswith("Nothing to report yet")
    with open(os.path.join(p.root, "results", "network_routes.csv")) as fh:
        assert len(list(csv.reader(fh))) == 1                       # the header only
