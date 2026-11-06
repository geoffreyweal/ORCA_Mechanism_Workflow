"""Elements, xyz files, formulas and connectivity.

Everything here works from plain xyz files: no SMILES and no chemistry toolkit.
Bonds are assigned from distances (covalent radii x a tolerance), which is how
the program checks SCAN start structures and the end points of an IRC.
"""
import hashlib
import math
import re
from collections import Counter

SYMBOLS = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
]
ATOMIC_NUMBER = {s: z for z, s in enumerate(SYMBOLS) if z}

# Covalent radii in Angstrom (Cordero et al., Dalton Trans. 2008, 2832), elements 1-86.
_RADII = [
    0.00, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
    2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16,
    2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40,
    2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87,
    1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50,
]
COVALENT_RADIUS = {s: _RADII[z] for s, z in ATOMIC_NUMBER.items()}

#: two atoms are bonded when their distance is at most BOND_TOLERANCE x (sum of covalent radii)
BOND_TOLERANCE = 1.25


class XYZError(ValueError):
    pass


class Molecule:
    """Symbols, coordinates (Angstrom) and the comment line of one xyz frame."""

    def __init__(self, symbols, coords, comment=""):
        self.symbols = list(symbols)
        self.coords = [tuple(c) for c in coords]
        self.comment = comment

    def __len__(self):
        return len(self.symbols)

    @property
    def formula(self):
        return formula_of(self.symbols)


def normalise_symbol(s):
    s = re.sub(r"[^A-Za-z]", "", s)          # "C1" -> "C", "Cl2" -> "Cl"
    if not s:
        raise XYZError("empty element symbol")
    s = s[0].upper() + s[1:].lower()
    if s not in ATOMIC_NUMBER:
        raise XYZError("unknown element '%s'" % s)
    return s


def _parse_frame(lines, k, path):
    try:
        n = int(lines[k].split()[0])
    except (ValueError, IndexError):
        raise XYZError("%s, line %d: expected the number of atoms" % (path, k + 1))
    if len(lines) < k + 2 + n:
        raise XYZError("%s: expected %d atoms but the file ends early" % (path, n))
    comment = lines[k + 1]
    sym, xyz = [], []
    for m, ln in enumerate(lines[k + 2:k + 2 + n]):
        p = ln.split()
        if len(p) < 4:
            raise XYZError("%s, line %d: expected 'element x y z'" % (path, k + 3 + m))
        try:
            sym.append(normalise_symbol(p[0]))
            xyz.append((float(p[1]), float(p[2]), float(p[3])))
        except ValueError as exc:
            raise XYZError("%s, line %d: %s" % (path, k + 3 + m, exc))
    return Molecule(sym, xyz, comment), k + 2 + n


def read_xyz(path):
    """First frame of an xyz file."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    while lines and not lines[0].strip():
        lines.pop(0)
    if not lines:
        raise XYZError("%s is empty" % path)
    return _parse_frame(lines, 0, path)[0]


def read_xyz_frames(path):
    """All frames of a (multi-frame) xyz file, e.g. an ORCA trajectory."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    frames, k = [], 0
    while k < len(lines):
        if not lines[k].strip():
            k += 1
            continue
        mol, k = _parse_frame(lines, k, path)
        frames.append(mol)
    return frames


def write_xyz(path, mol, comment=None):
    with open(path, "w") as fh:
        fh.write("%d\n%s\n" % (len(mol), mol.comment if comment is None else comment))
        for s, (x, y, z) in zip(mol.symbols, mol.coords):
            fh.write("%-2s %14.8f %14.8f %14.8f\n" % (s, x, y, z))


_CHARGE_RE = re.compile(r"\bcharge\s*[=:]\s*([+-]?\d+)", re.I)
_MULT_RE = re.compile(r"\b(?:mult|multiplicity|spin)\s*[=:]\s*(\d+)", re.I)


def comment_charge_mult(comment):
    """charge / multiplicity written on an xyz comment line as 'charge=1 mult=1' (either may be absent)."""
    c = _CHARGE_RE.search(comment or "")
    m = _MULT_RE.search(comment or "")
    return (int(c.group(1)) if c else None, int(m.group(1)) if m else None)


def energy_from_comment(comment):
    """Energy written by ORCA on the comment line of SCAN images ('... E -1234.567890')."""
    m = re.search(r"\bE\s+(-?\d+\.\d+)", comment or "")
    return float(m.group(1)) if m else None


# ----------------------------------------------------------------------------- formulas

def formula_of(symbols):
    return Counter(symbols)


def add_formulas(*formulas):
    tot = Counter()
    for f in formulas:
        tot.update(f)
    return tot


def formula_string(formula, charge=0):
    """Hill order: C, H, then alphabetical (alphabetical throughout when there is no C)."""
    f = {k: v for k, v in formula.items() if v}
    keys = (["C"] + (["H"] if "H" in f else []) + sorted(k for k in f if k not in ("C", "H"))) if "C" in f else sorted(f)
    s = "".join("%s%s" % (k, f[k] if f[k] != 1 else "") for k in keys)
    if charge:
        s += ("+" if charge > 0 else "-") * abs(charge) if abs(charge) <= 2 else "(%+d)" % charge
    return s or "(nothing)"


def electron_count(symbols, charge):
    return sum(ATOMIC_NUMBER[s] for s in symbols) - charge


def multiplicity_problem(symbols, charge, mult):
    """Text if the charge / multiplicity combination is impossible, else None."""
    if mult < 1:
        return "multiplicity must be 1 or more"
    ne = electron_count(symbols, charge)
    if ne < 0:
        return "charge %+d leaves a negative number of electrons" % charge
    if (ne % 2) == (mult % 2):
        return ("%d electrons cannot have multiplicity %d (an %s electron count needs an %s multiplicity): "
                "check the charge" % (ne, mult, "even" if ne % 2 == 0 else "odd", "odd" if ne % 2 == 0 else "even"))
    if mult - 1 > ne:
        return "multiplicity %d needs at least %d unpaired electrons" % (mult, mult - 1)
    return None


# ----------------------------------------------------------------------------- geometry / bonds

def distance(coords, i, j):
    return math.dist(coords[i], coords[j])


def bond_cutoff(si, sj, tolerance=BOND_TOLERANCE):
    return tolerance * (COVALENT_RADIUS.get(si, 1.5) + COVALENT_RADIUS.get(sj, 1.5))


def is_bonded(mol, i, j, tolerance=BOND_TOLERANCE):
    return distance(mol.coords, i, j) <= bond_cutoff(mol.symbols[i], mol.symbols[j], tolerance)


def bonds(mol, tolerance=BOND_TOLERANCE):
    """Set of bonded atom pairs (i < j)."""
    out = set()
    n = len(mol)
    sym, xyz = mol.symbols, mol.coords
    maxr = max((COVALENT_RADIUS.get(s, 1.5) for s in sym), default=0.0)
    for i in range(n):
        xi, yi, zi = xyz[i]
        ri = COVALENT_RADIUS.get(sym[i], 1.5)
        for j in range(i + 1, n):
            dx = abs(xyz[j][0] - xi)
            if dx > tolerance * (ri + maxr):
                continue
            d = math.sqrt(dx * dx + (xyz[j][1] - yi) ** 2 + (xyz[j][2] - zi) ** 2)
            if d <= tolerance * (ri + COVALENT_RADIUS.get(sym[j], 1.5)):
                out.add((i, j))
    return out


def adjacency(n, bondset):
    adj = [[] for _ in range(n)]
    for i, j in bondset:
        adj[i].append(j)
        adj[j].append(i)
    return adj


def fragments(n, bondset):
    """Connected components (lists of atom indices), largest first."""
    adj = adjacency(n, bondset)
    seen, comps = [False] * n, []
    for s in range(n):
        if seen[s]:
            continue
        stack, comp = [s], []
        seen[s] = True
        while stack:
            a = stack.pop()
            comp.append(a)
            for b in adj[a]:
                if not seen[b]:
                    seen[b] = True
                    stack.append(b)
        comps.append(sorted(comp))
    comps.sort(key=lambda c: (-len(c), c[0]))
    return comps


def _h(text):
    return hashlib.sha1(text.encode()).hexdigest()[:16]


def graph_hash(symbols, adj, atoms):
    """Connectivity fingerprint of the atoms `atoms` (Weisfeiler-Lehman colour refinement).

    Identical for the same molecule whatever the atom order; different for different
    connectivity (stereochemistry and conformation are not seen)."""
    atoms = list(atoms)
    inside = set(atoms)
    labels = {a: symbols[a] for a in atoms}
    nclass = len(set(labels.values()))
    for _ in range(len(atoms) + 1):
        new = {a: _h(labels[a] + "(" + ",".join(sorted(labels[b] for b in adj[a] if b in inside)) + ")") for a in atoms}
        nnew = len(set(new.values()))
        labels = new
        if nnew == nclass:
            break
        nclass = nnew
    return _h(";".join(sorted(labels.values())))


def fragment_fingerprints(mol, tolerance=BOND_TOLERANCE):
    """Sorted list of (fingerprint, formula string) for every separate molecule in `mol`."""
    bs = bonds(mol, tolerance)
    adj = adjacency(len(mol), bs)
    out = []
    for comp in fragments(len(mol), bs):
        out.append((graph_hash(mol.symbols, adj, comp), formula_string(formula_of(mol.symbols[a] for a in comp))))
    return sorted(out)


def molecule_fingerprint(mol, tolerance=BOND_TOLERANCE):
    """One fingerprint for the whole structure (all fragments)."""
    return _h("+".join(h for h, _ in fragment_fingerprints(mol, tolerance)))
