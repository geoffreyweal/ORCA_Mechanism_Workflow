"""Placing the molecules of a complex built from their optimised structures (complexes.optimise: false).

Each molecule keeps its optimised shape exactly (a rigid body). It is first fitted onto the place the complex's xyz
file has it (workflow.assemble); then every molecule but the largest is moved and turned, as a whole, until

* the contacts drawn in the complex's xyz file are as long as drawn: its hydrogen bonds between molecules (a polar
  H ... N/O/S up to HBOND_MAX), and any pair drawn closer than the minimum below (e.g. H(A) ... C of a C-protonation);
* no other atoms of different molecules are closer than CLASH_SCALE x the sum of their van der Waals radii (a polar
  H and an N/O/S: HBOND_MIN), nor close enough for mechflow to see a bond between them.

Without this second part, a molecule whose optimised shape differs from the one drawn in the complex (a phenyl ring
turned another way, a ring flipped, an OH pointing elsewhere) can be put on top of the other molecule. Nothing is
optimised: there is no energy but this simple penalty, and the first SCAN point relaxes the whole complex anyway.
"""
import math
import random
import zlib

from . import chem

VDW_RADIUS = {"H": 1.20, "He": 1.40, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "Si": 2.10, "P": 1.80,
              "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98}       # Bondi (B: Mantina); others 1.80
CLASH_SCALE = 0.75     # atoms of different molecules kept at least this x (sum of van der Waals radii) apart
HBOND_MIN = 1.70       # ... a polar H and an N, O or S at least this far (Angstrom)
HBOND_MAX = 2.60       # a polar H ... N/O/S drawn up to this long is a contact kept at its drawn length
POLAR = ("N", "O", "S")
SKIN = 2.5             # Angstrom: pairs further apart than their limit + SKIN are left out until the next list
TRIES = 16             # further starts (the first placement turned about its contacts) while the penalty stays high
GOOD = 0.05            # penalty below which a placement is taken without trying further starts
TOO_CLOSE = 0.8        # a built complex with atoms of different molecules closer than this x their minimum (apart
                       # from the drawn contacts) is not used: no SCAN starts from it


def _sub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _turn(v, axis, angle):
    """v turned by angle about the unit vector axis (Rodrigues)."""
    c, s = math.cos(angle), math.sin(angle)
    kv = _cross(axis, v)
    kd = axis[0] * v[0] + axis[1] * v[1] + axis[2] * v[2]
    return [v[i] * c + kv[i] * s + axis[i] * kd * (1 - c) for i in range(3)]


def polar_hydrogens(symbols, bonds):
    """The H atoms bonded to an N, O or S."""
    out = set()
    for i, j in bonds:
        for h, x in ((i, j), (j, i)):
            if symbols[h] == "H" and symbols[x] in POLAR:
                out.add(h)
    return out


def minimum(symbols, polar, i, j, tolerance=chem.BOND_TOLERANCE):
    """How close atoms i and j of different molecules may come."""
    if (i in polar and symbols[j] in POLAR) or (j in polar and symbols[i] in POLAR):
        d = HBOND_MIN
    else:
        d = CLASH_SCALE * (VDW_RADIUS.get(symbols[i], 1.8) + VDW_RADIUS.get(symbols[j], 1.8))
    return max(d, 1.1 * chem.bond_cutoff(symbols[i], symbols[j], tolerance))


def contacts(symbols, drawn, part_of, polar, tolerance=chem.BOND_TOLERANCE):
    """[(i, j, drawn length)]: the hydrogen bonds between molecules drawn in the complex, and the pairs drawn closer
    than their minimum; if there is none, the closest pair of atoms of different molecules."""
    n = len(symbols)
    out, closest = [], None
    for i in range(n):
        for j in range(i + 1, n):
            if part_of[i] == part_of[j]:
                continue
            d = math.dist(drawn[i], drawn[j])
            if closest is None or d < closest[2]:
                closest = (i, j, d)
            hb = (i in polar and symbols[j] in POLAR) or (j in polar and symbols[i] in POLAR)
            if (hb and d <= HBOND_MAX) or d < minimum(symbols, polar, i, j, tolerance):
                out.append((i, j, d))
    return out or ([closest] if closest else [])


class _Penalty:
    """Contacts held at their drawn lengths, other pairs kept beyond their minimum (three times as stiff: a contact
    stretched a little is better than atoms too close); with a neighbour list."""

    def __init__(self, symbols, polar, part_of, contact, tolerance, k=10.0, k_clash=30.0):
        self.k, self.k_clash = k, k_clash
        self.contact = contact
        self.part_of = part_of
        keep = {(min(i, j), max(i, j)) for i, j, _ in contact}
        n = len(symbols)
        self.pairs = [(i, j, minimum(symbols, polar, i, j, tolerance)) for i in range(n) for j in range(i + 1, n)
                      if part_of[i] != part_of[j] and (i, j) not in keep]
        self.near = None
        self.listed_at = None

    def neighbours(self, coords):
        self.near = [(i, j, r) for i, j, r in self.pairs if math.dist(coords[i], coords[j]) < r + SKIN]
        self.listed_at = [list(c) for c in coords]

    def __call__(self, coords):
        if self.near is None or max(math.dist(a, b) for a, b in zip(coords, self.listed_at)) > SKIN / 2:
            self.neighbours(coords)
        k, kc, e = self.k, self.k_clash, 0.0
        g = [[0.0, 0.0, 0.0] for _ in coords]
        for i, j, d0 in self.contact:
            v = _sub(coords[i], coords[j])
            d = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) or 1e-9
            e += k * (d - d0) ** 2
            f = 2 * k * (d - d0) / d
            for x in range(3):
                g[i][x] += f * v[x]
                g[j][x] -= f * v[x]
        for i, j, r in self.near:
            v = _sub(coords[i], coords[j])
            d2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
            if d2 >= r * r:
                continue
            d = math.sqrt(d2) or 1e-9
            e += kc * (r - d) ** 2
            f = -2 * kc * (r - d) / d
            for x in range(3):
                g[i][x] += f * v[x]
                g[j][x] -= f * v[x]
        return e, g


def _relax(coords, bodies, moving, penalty, steps=2000):
    """Move and turn each moving body as a whole down the penalty (steepest descent, step sizes adapted)."""
    coords = [list(c) for c in coords]
    e, g = penalty(coords)
    move, turn = 0.02, 0.01                   # Angstrom, radian
    for _ in range(steps):
        if e < 1e-6:
            break
        trial = [list(c) for c in coords]
        for b in moving:
            idx = bodies[b]
            centre = [sum(coords[a][x] for a in idx) / len(idx) for x in range(3)]
            force = [-sum(g[a][x] for a in idx) for x in range(3)]
            torque = [0.0, 0.0, 0.0]
            for a in idx:
                t = _cross(_sub(coords[a], centre), [-g[a][0], -g[a][1], -g[a][2]])
                torque = [torque[x] + t[x] for x in range(3)]
            fn = math.sqrt(sum(f * f for f in force))
            tn = math.sqrt(sum(t * t for t in torque))
            shift = [move * f / fn * min(fn, 1.0) for f in force] if fn > 1e-12 else [0.0, 0.0, 0.0]
            axis = [t / tn for t in torque] if tn > 1e-12 else [1.0, 0.0, 0.0]
            angle = turn * min(tn, 1.0)
            for a in idx:
                r = _turn(_sub(coords[a], centre), axis, angle)
                trial[a] = [centre[x] + r[x] + shift[x] for x in range(3)]
        e2, g2 = penalty(trial)
        if e2 < e:
            coords, e, g = trial, e2, g2
            move, turn = min(move * 1.2, 0.3), min(turn * 1.2, 0.15)
        else:
            move, turn = move * 0.5, turn * 0.5
            if move < 1e-6 and turn < 1e-6:
                break
    return coords, e


def place(symbols, drawn, start, bodies, bonds, tolerance=chem.BOND_TOLERANCE, seed="complex"):
    """Move the molecules of a complex (bodies: lists of atom numbers; start: their first placement) as rigid bodies so
    that the drawn contacts are as drawn and nothing else is too close. The largest molecule stays where it is.
    Deterministic for a given seed (a text, e.g. the complex's ID).

    Returns (coords, contacts as [(i, j, drawn length)], penalty left)."""
    part_of = {a: b for b, idx in enumerate(bodies) for a in idx}
    polar = polar_hydrogens(symbols, bonds)
    contact = contacts(symbols, drawn, part_of, polar, tolerance)
    if len(bodies) < 2:
        return [list(c) for c in start], contact, 0.0
    anchor = max(range(len(bodies)), key=lambda b: (len(bodies[b]), -b))
    moving = [b for b in range(len(bodies)) if b != anchor]
    penalty = _Penalty(symbols, polar, part_of, contact, tolerance)
    best = _relax(start, bodies, moving, penalty)
    rng = random.Random(zlib.crc32(str(seed).encode()))
    for _ in range(TRIES):
        if best[1] < GOOD:
            break
        trial = [list(c) for c in start]
        for b in moving:
            idx = bodies[b]
            mine = [i if part_of[i] == b else j for i, j, _ in contact if b in (part_of[i], part_of[j])] or idx
            pivot = [sum(start[a][x] for a in mine) / len(mine) for x in range(3)]
            axis = [rng.gauss(0.0, 1.0) for _ in range(3)]
            n = math.sqrt(sum(a * a for a in axis)) or 1.0
            axis = [a / n for a in axis]
            angle = rng.uniform(-math.pi / 2, math.pi / 2)
            for a in idx:
                r = _turn(_sub(start[a], pivot), axis, angle)
                trial[a] = [pivot[x] + r[x] for x in range(3)]
        res = _relax(trial, bodies, moving, penalty)
        if res[1] < best[1]:
            best = res
    return best[0], contact, best[1]


def too_close(symbols, coords, part_of, polar, contact, tolerance=chem.BOND_TOLERANCE):
    """(i, j, distance, minimum) of the pair of atoms of different molecules closest to each other for their minimum,
    if it is closer than TOO_CLOSE x its minimum (the drawn contacts aside); else None."""
    keep = {(min(i, j), max(i, j)) for i, j, _ in contact}
    worst = None
    n = len(symbols)
    for i in range(n):
        for j in range(i + 1, n):
            if part_of[i] == part_of[j] or (i, j) in keep:
                continue
            r = minimum(symbols, polar, i, j, tolerance)
            d = math.dist(coords[i], coords[j])
            if d < TOO_CLOSE * r and (worst is None or d / r < worst[2] / worst[3]):
                worst = (i, j, d, r)
    return worst
