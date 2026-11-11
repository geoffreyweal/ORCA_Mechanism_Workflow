"""A mechflow project: project.yaml plus a folder of species xyz files.

    my_reaction/
        project.yaml          settings, species names, steps, mechanisms, energies
        species/R1.xyz ...    one xyz file per species (the file name is the species ID)
        calcs/                written by mechflow (one folder per species and per step)
        results/  maps/       written by 'mechflow collect' and 'mechflow map'

Naming convention (recommended, used for map styles): R1, R2 ... reactants, I1 ... intermediates,
P1 ... products, C1 ... catalysts and other species used during the reaction, S1 ... steps
(step Sn goes through transition state TSn). Letters after the number (I2a / I2b) are allowed.
"""
import os
import re
from collections import Counter, OrderedDict

import yaml

from . import chem
from .defaults import DFT_KINDS, JOB_KINDS, PROTON_TRANSFER_MODES, deep_merge, default_settings

PROJECT_FILE = "project.yaml"
ID_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.\-]*$")
ROLE_OF_PREFIX = OrderedDict([("R", "reactant"), ("I", "intermediate"), ("P", "product"), ("C", "catalyst")])
RESOURCE_KEYS = ("nprocs", "maxcore", "time", "mem_per_cpu", "partition", "extra_sbatch", "extra_lines", "orca_module")


class ProjectError(Exception):
    pass


class BalanceError(ProjectError):
    pass


def find_project(start="."):
    """The folder holding project.yaml: `start` itself or the nearest parent."""
    d = os.path.abspath(start)
    if os.path.isfile(d):
        d = os.path.dirname(d)
    while True:
        if os.path.isfile(os.path.join(d, PROJECT_FILE)):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            return None
        d = parent


def natural_key(text):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", str(text))]


def species_sort_key(sid):
    order = list(ROLE_OF_PREFIX)
    k = order.index(sid[0]) if sid[:1] in order else len(order)
    return (k, natural_key(sid))


_TERM = re.compile(r"^\s*(?:(\d+)\s*\*?\s*)?([^\s+*]+)\s*$")


def parse_expression(text):
    """'R1 + R2 + 2 C1' -> Counter({'R1': 1, 'R2': 1, 'C1': 2})"""
    c = Counter()
    for part in str(text).split("+"):
        m = _TERM.match(part)
        if not m:
            raise ProjectError("cannot read '%s' (write species as 'I3 + R2' or '2 C1')" % text)
        c[m.group(2)] += int(m.group(1) or 1)
    return c


def expression(counter, order=None):
    items = sorted(counter.items(), key=lambda kv: species_sort_key(kv[0])) if order is None else \
        [(k, counter[k]) for k in order if counter.get(k)]
    return " + ".join(("%d %s" % (n, s)) if n > 1 else s for s, n in items if n > 0) or "(nothing)"


def split_reaction(text):
    parts = re.split(r"\s*(?:<=>|->|=>|→)\s*", str(text))
    if len(parts) != 2:
        raise ProjectError("reaction '%s' must look like 'I3 + R2 -> I4'" % text)
    return parse_expression(parts[0]), parse_expression(parts[1])


# ============================================================================== objects

class Species:
    def __init__(self, sid, path):
        self.id = sid
        self.path = path
        self.name = ""
        self.role = None
        self.mol = None
        self.charge = 0
        self.mult = 1
        self.overrides = {}
        self.used_by = []

    @property
    def formula(self):
        return chem.formula_of(self.mol.symbols) if self.mol else Counter()

    @property
    def formula_str(self):
        return chem.formula_string(self.formula, self.charge)

    @property
    def natoms(self):
        return len(self.mol) if self.mol else 0

    def __repr__(self):
        return "Species(%s)" % self.id


class Step:
    def __init__(self, sid):
        self.id = sid
        self.ts = "TS" + sid[1:] if sid.startswith("S") else "TS-" + sid
        self.name = ""
        self.left = Counter()
        self.right = Counter()
        self.start = None             # species whose optimised structure starts the SCAN
        self.atoms = None             # 0-based (ORCA) atom numbers of the scanned pair
        self.to = None                # final distance (Angstrom)
        self.scan_step = None
        self.mode = None              # 'stretch' (pair bonded at the start) or 'contract'
        self.d0 = None                # distance in the species xyz file
        self.charge = 0
        self.mult = 1
        self.overrides = {}
        self.used_by = []

    def reaction_text(self):
        return "%s -> %s" % (expression(self.left), expression(self.right))

    def __repr__(self):
        return "Step(%s)" % self.id


class State:
    """One point of an energy profile: species (and/or a TS) plus the species added to balance it."""

    def __init__(self, kind, species, step=None, label=None, source="path"):
        self.kind = kind                      # 'species' or 'ts'
        self.species = Counter(species)       # for a TS: its step's start species (same atoms as the TS)
        self.step = step
        self.label = label
        self.source = source                  # reference | path | step | product
        self.added = Counter()

    @property
    def full(self):
        return self.species + self.added

    def same_as(self, other):
        return self.kind == other.kind == "species" and self.full == other.full

    def __repr__(self):
        return "State(%s, %s)" % (self.kind, self.label)


class Mechanism:
    def __init__(self, mid):
        self.id = mid
        self.name = ""
        self.family = None
        self.tokens = []
        self.states = []
        self.product = None
        self.reference = None         # own reference state (e.g. two molecules of a reactant); None = project's
        self.notes = []


# ============================================================================== project

class Project:
    def __init__(self, root):
        self.root = os.path.abspath(root)
        self.errors = []
        self.warnings = []
        self.notes = []
        self.raw = {}
        self.settings = default_settings()
        self.species = OrderedDict()
        self.steps = OrderedDict()
        self.mechanisms = OrderedDict()
        self.reference = Counter()
        self.balance_with = []
        self.proton_transfers = "auto"
        self._balance_cache = {}

    # ------------------------------------------------------------------ paths
    @property
    def name(self):
        return self.settings["project"].get("name") or os.path.basename(self.root)

    @property
    def index_base(self):
        return int(self.settings["project"].get("atom_numbering", 0))

    @property
    def tolerance(self):
        return float(self.settings["project"].get("bond_tolerance", chem.BOND_TOLERANCE))

    def path(self, *parts):
        return os.path.join(self.root, *parts)

    def species_job_dir(self, sid):
        return self.path("calcs", "Species", sid, "Step1_Geo_Opt")

    def step_dir(self, sid):
        return self.path("calcs", "Steps", sid)

    def atom_label(self, mol, i):
        """'C12' style label in the project's atom numbering."""
        return "%s%d" % (mol.symbols[i], i + self.index_base)

    def role(self, sid):
        sp = self.species.get(sid)
        if sp is not None and sp.role:
            return sp.role
        return ROLE_OF_PREFIX.get(sid[:1], "species")

    # ------------------------------------------------------------------ loading
    @classmethod
    def load(cls, root):
        p = cls(root)
        p._load()
        return p

    def _load(self):
        fn = self.path(PROJECT_FILE)
        try:
            with open(fn) as fh:
                raw = yaml.safe_load(fh) or {}
        except OSError:
            raise ProjectError("no %s in %s" % (PROJECT_FILE, self.root))
        except yaml.YAMLError as exc:
            raise ProjectError("%s is not valid YAML:\n%s" % (fn, exc))
        if not isinstance(raw, dict):
            raise ProjectError("%s must be a mapping of sections (project, orca, species, steps, ...)" % fn)
        self.raw = raw
        known = {"project", "orca", "resources", "slurm", "stage_overrides", "scan", "species", "steps",
                 "mechanisms", "energies", "single_point", "chain"}
        for k in raw:
            if k not in known:
                self.warnings.append("project.yaml: unknown section '%s' (ignored)" % k)
        settings_part = {k: v for k, v in raw.items() if k not in ("species", "steps", "mechanisms")}
        self.settings = deep_merge(default_settings(), settings_part)
        self._check_settings()
        self._load_species(raw.get("species") or {})
        self._load_steps(raw.get("steps") or {})
        self._load_energies()
        self._load_mechanisms(raw.get("mechanisms") or {})

    def _check_settings(self):
        o = self.settings["orca"]
        for kind in DFT_KINDS:
            if not (o.get("keywords") or {}).get(kind):
                self.errors.append("orca.keywords.%s is missing" % kind)
        sp = self.settings.get("single_point") or {}
        if sp.get("enabled") and not str(sp.get("method") or "").strip():
            self.errors.append("single_point.method is missing (single_point.enabled is true)")
        texts = [o.get("method", ""), o.get("scf_block", ""), o.get("solvent_block", ""),
                 sp.get("method", ""), sp.get("keywords", ""), sp.get("block", "")]
        texts += list((o.get("keywords") or {}).values()) + list((o.get("blocks") or {}).values())
        for t in texts:
            for ln in str(t or "").splitlines():
                if ln.count("#") >= 2:
                    self.warnings.append("ORCA ends a comment at a second '#', so this line is read partly as input: %r"
                                         % ln.strip())
        ch = self.settings.get("chain") or {}
        if not isinstance(ch.get("setup") or [], list):
            self.errors.append("chain.setup must be a list of shell lines")
        if self.index_base not in (0, 1):
            self.errors.append("project.atom_numbering must be 0 (as in ORCA) or 1")
        for kind in (self.settings.get("stage_overrides") or {}):
            if kind not in JOB_KINDS:
                self.errors.append("stage_overrides: '%s' is not one of %s" % (kind, ", ".join(JOB_KINDS)))

    def _load_species(self, entries):
        sdir = self.path("species")
        files = {}
        if os.path.isdir(sdir):
            for fn in os.listdir(sdir):
                if fn.lower().endswith(".xyz") and not fn.startswith("."):
                    files[fn[:-4]] = os.path.join(sdir, fn)
        else:
            self.errors.append("no species/ folder: put one xyz file per species in %s" % sdir)
        if not isinstance(entries, dict):
            self.errors.append("'species' must map species IDs to names or settings")
            entries = {}
        ids = list(entries) + sorted((f for f in files if f not in entries), key=species_sort_key)
        for sid in ids:
            sid = str(sid)
            e = entries.get(sid)
            if e is None:
                e = {}
                if entries:
                    self.notes.append("species/%s.xyz is not listed in project.yaml (included anyway)" % sid)
            elif isinstance(e, str):
                e = {"name": e}
            elif not isinstance(e, dict):
                self.errors.append("species %s: expected a name or a mapping" % sid)
                continue
            if not ID_RE.match(sid):
                self.errors.append("species ID '%s': use letters, digits, '_', '-' or '.', starting with a letter "
                                   "(no '+' or spaces)" % sid)
                continue
            path = os.path.join(self.root, e["file"]) if e.get("file") else files.get(sid)
            sp = Species(sid, path)
            sp.name = str(e.get("name", "") or "")
            sp.role = e.get("role")
            sp.overrides = {k: v for k, v in e.items() if k in RESOURCE_KEYS}
            if not path or not os.path.isfile(path):
                self.errors.append("species %s: no xyz file (expected species/%s.xyz)" % (sid, sid))
                self.species[sid] = sp
                continue
            try:
                sp.mol = chem.read_xyz(path)
            except (chem.XYZError, OSError) as exc:
                self.errors.append("species %s: %s" % (sid, exc))
                self.species[sid] = sp
                continue
            cq, cm = chem.comment_charge_mult(sp.mol.comment)
            for what, yval, cval, default in (("charge", e.get("charge"), cq, 0), ("mult", e.get("mult"), cm, 1)):
                if yval is not None and cval is not None and int(yval) != cval:
                    self.errors.append("species %s: %s is %s in project.yaml but %s in the xyz comment line"
                                       % (sid, what, yval, cval))
                val = int(yval) if yval is not None else (cval if cval is not None else default)
                setattr(sp, what, val)
            prob = chem.multiplicity_problem(sp.mol.symbols, sp.charge, sp.mult)
            if prob:
                self.errors.append("species %s (%s): %s" % (sid, sp.formula_str, prob))
            self.species[sid] = sp
        # the same molecule twice?
        seen = {}
        for sp in self.species.values():
            if sp.mol is None:
                continue
            fp = (chem.molecule_fingerprint(sp.mol, self.tolerance), sp.charge)
            if fp in seen:
                self.warnings.append("species %s and %s have the same atoms and bonds (conformers or stereoisomers? "
                                     "otherwise one of them is not needed)" % (seen[fp], sp.id))
            else:
                seen[fp] = sp.id

    def side_formula(self, counter):
        f, q = Counter(), 0
        for sid, n in counter.items():
            sp = self.species.get(sid)
            if sp is None or sp.mol is None:
                raise ProjectError("unknown species '%s'" % sid)
            for el, c in sp.formula.items():
                f[el] += c * n
            q += sp.charge * n
        return f, q

    def _load_steps(self, entries):
        if not isinstance(entries, dict):
            self.errors.append("'steps' must map step IDs (S1, S2 ...) to their settings")
            return
        for sid, e in entries.items():
            sid = str(sid)
            st = Step(sid)
            self.steps[sid] = st
            if not ID_RE.match(sid):
                self.errors.append("step ID '%s' is not allowed (letters, digits, '_', '-', '.')" % sid)
                continue
            if not isinstance(e, dict):
                self.errors.append("step %s: expected 'name', 'reaction' and 'scan'" % sid)
                continue
            st.name = str(e.get("name", "") or "")
            if e.get("ts"):
                st.ts = str(e["ts"])
            st.overrides = {k: v for k, v in (e.get("resources") or {}).items()}
            try:
                st.left, st.right = split_reaction(e.get("reaction", ""))
            except ProjectError as exc:
                self.errors.append("step %s: %s" % (sid, exc))
                continue
            missing = [s for s in list(st.left) + list(st.right) if s not in self.species]
            if missing:
                self.errors.append("step %s: unknown species %s" % (sid, ", ".join(missing)))
                continue
            try:
                fl, ql = self.side_formula(st.left)
                fr, qr = self.side_formula(st.right)
            except ProjectError as exc:
                self.errors.append("step %s: %s" % (sid, exc))
                continue
            if fl != fr or ql != qr:
                self.errors.append("step %s: the two sides of '%s' differ: %s vs %s" % (
                    sid, st.reaction_text(), chem.formula_string(fl, ql), chem.formula_string(fr, qr)))
                continue
            scan = e.get("scan") or {}
            if not isinstance(scan, dict):
                self.errors.append("step %s: scan must look like {from: I4, atoms: [3, 9], to: 3.2}" % sid)
                continue
            start = scan.get("from")
            if start is None:
                singles = [c for c in (st.right, st.left) if sum(c.values()) == 1]
                start = list(singles[-1])[0] if singles else None
                if start is None:
                    self.errors.append("step %s: give scan.from (a species holding every atom of the step)" % sid)
                    continue
            start = str(start)
            if start not in self.species or self.species[start].mol is None:
                self.errors.append("step %s: scan.from '%s' is not a species" % (sid, start))
                continue
            sp = self.species[start]
            if (sp.formula, sp.charge) != (fl, ql):
                self.errors.append("step %s: scan.from %s (%s) does not hold every atom of the step (%s); "
                                   "use a species with all the atoms, e.g. the bonded adduct" % (
                                       sid, start, sp.formula_str, chem.formula_string(fl, ql)))
                continue
            st.start, st.charge, st.mult = start, sp.charge, sp.mult
            atoms = scan.get("atoms")
            if not (isinstance(atoms, (list, tuple)) and len(atoms) == 2 and all(isinstance(a, int) for a in atoms)):
                self.errors.append("step %s: scan.atoms must be two atom numbers, e.g. [3, 9]" % sid)
                continue
            i, j = (a - self.index_base for a in atoms)
            if not (0 <= i < sp.natoms and 0 <= j < sp.natoms) or i == j:
                self.errors.append("step %s: scan.atoms %s are not two different atoms of %s (%d atoms, numbered from %d)"
                                   % (sid, list(atoms), start, sp.natoms, self.index_base))
                continue
            st.atoms = (i, j)
            st.d0 = chem.distance(sp.mol.coords, i, j)
            st.mode = "stretch" if chem.is_bonded(sp.mol, i, j, self.tolerance) else "contract"
            to = scan.get("to")
            if to is None:
                if st.mode == "stretch":
                    to = float(self.settings["scan"].get("stretch_to", 3.2))
                else:
                    si, sj = sp.mol.symbols[i], sp.mol.symbols[j]
                    to = round(1.05 * (chem.COVALENT_RADIUS[si] + chem.COVALENT_RADIUS[sj]), 2)
                self.notes.append("step %s: no scan.to given, using %.2f A" % (sid, to))
            st.to = float(to)
            st.scan_step = float(scan.get("step") or self.settings["scan"]["step"])
            if st.mode == "stretch" and st.to < st.d0 + 0.3:
                self.errors.append("step %s: atoms %s are bonded in %s (%.2f A), so the SCAN stretches them, "
                                   "but scan.to = %.2f A is not longer" % (sid, list(atoms), start, st.d0, st.to))
            if st.mode == "contract" and st.to > st.d0 - 0.3:
                self.errors.append("step %s: atoms %s are not bonded in %s (%.2f A), so the SCAN brings them "
                                   "together, but scan.to = %.2f A is not shorter" % (sid, list(atoms), start, st.d0, st.to))

    def _load_energies(self):
        e = self.settings["energies"]
        try:
            if e.get("reference"):
                self.reference = parse_expression(e["reference"])
            else:
                self.reference = Counter({s: 1 for s in self.species if self.role(s) == "reactant"})
        except ProjectError as exc:
            self.errors.append("energies.reference: %s" % exc)
        bad = [s for s in self.reference if s not in self.species]
        if bad:
            self.errors.append("energies.reference: unknown species %s" % ", ".join(bad))
            self.reference = Counter()
        bw = e.get("balance_with") or []
        if isinstance(bw, str):
            bw = [x.strip() for x in bw.replace("+", ",").split(",") if x.strip()]
        self.balance_with = [str(s) for s in bw]
        bad = [s for s in self.balance_with if s not in self.species]
        if bad:
            self.errors.append("energies.balance_with: unknown species %s" % ", ".join(bad))
            self.balance_with = [s for s in self.balance_with if s in self.species]
        pt = str(e.get("proton_transfers") or "auto").strip().lower()
        if pt not in PROTON_TRANSFER_MODES:
            self.errors.append("energies.proton_transfers: '%s' is not one of %s" % (
                e.get("proton_transfers"), ", ".join(PROTON_TRANSFER_MODES)))
            pt = "auto"
        self.proton_transfers = pt

    # ------------------------------------------------------------------ balancing
    @property
    def auto_species(self):
        """Species that are added to states automatically (unreacted reference species and balancers)."""
        out = self.auto_for(self.reference)
        for m in self.mechanisms.values():
            for s in (m.reference or ()):
                if s not in out:
                    out.append(s)
        return out

    def auto_for(self, reference):
        """Species added automatically when balancing against `reference`."""
        out = list(reference)
        for s in self.balance_with:
            if s not in out:
                out.append(s)
        return out

    def balance(self, counter, reference=None):
        """Species to add to `counter` so that it has exactly the atoms and charge of the reference.

        Uses the reference species (still unreacted) and energies.balance_with. `reference` is the
        project's reference unless a mechanism has its own. Returns a Counter."""
        ref = self.reference if reference is None else reference
        key = (tuple(sorted(counter.items())), tuple(sorted(ref.items())))
        if key in self._balance_cache:
            res = self._balance_cache[key]
            if isinstance(res, Exception):
                raise res
            return Counter(res)
        try:
            res = self._solve_balance(counter, ref)
        except BalanceError as exc:
            self._balance_cache[key] = exc
            raise
        self._balance_cache[key] = res
        return Counter(res)

    def _solve_balance(self, counter, reference=None):
        ref = self.reference if reference is None else reference
        fref, qref = self.side_formula(ref)
        f, q = self.side_formula(counter)
        diff = {el: fref.get(el, 0) - f.get(el, 0) for el in set(fref) | set(f)}
        neg = {el: -v for el, v in diff.items() if v < 0}
        if neg:
            raise BalanceError("%s has more atoms than the reference %s (extra %s)" % (
                expression(counter), expression(ref), chem.formula_string(neg)))
        dq = qref - q
        cands = [s for s in self.auto_for(ref) if self.species[s].mol is not None]
        forms = [(s, self.species[s].formula, self.species[s].charge) for s in cands]
        sols = []

        def dfs(k, rem, rq, chosen):
            if len(sols) > 20:
                return
            if all(v == 0 for v in rem.values()):
                if rq == 0:
                    sols.append(Counter({s: n for s, n in chosen if n}))
                return
            if k == len(forms):
                return
            s, fs, qs = forms[k]
            if not fs:
                return dfs(k + 1, rem, rq, chosen)
            nmax = min(rem.get(el, 0) // c for el, c in fs.items())
            for n in range(nmax, -1, -1):
                rem2 = dict(rem)
                for el, c in fs.items():
                    rem2[el] -= n * c
                dfs(k + 1, rem2, rq - n * qs, chosen + [(s, n)])

        dfs(0, diff, dq, [])
        if not sols:
            raise BalanceError("%s cannot be balanced against the reference %s with %s (missing %s, charge %+d)" % (
                expression(counter), expression(ref), ", ".join(cands) or "nothing",
                chem.formula_string({k: v for k, v in diff.items() if v}), dq))
        best = min(sum(s.values()) for s in sols)
        bests = [s for s in sols if sum(s.values()) == best]
        uniq = []
        for s in bests:
            if s not in uniq:
                uniq.append(s)
        if len(uniq) > 1:
            raise BalanceError("%s can be balanced in more than one way (%s): write the state out in full" % (
                expression(counter), " or ".join(expression(u) for u in uniq)))
        return uniq[0]

    # ------------------------------------------------------------------ mechanisms
    def _load_mechanisms(self, entries):
        if not isinstance(entries, dict):
            self.errors.append("'mechanisms' must map mechanism IDs to a name and a path")
            return
        for mid, e in entries.items():
            mid = str(mid)
            m = Mechanism(mid)
            self.mechanisms[mid] = m
            if isinstance(e, list):
                e = {"path": e}
            if not isinstance(e, dict) or not isinstance(e.get("path"), list):
                self.errors.append("mechanism %s: needs a path, e.g. path: [R1 + R2, S1, I1, S2, P1]" % mid)
                continue
            m.name = str(e.get("name", "") or "")
            m.family = str(e["family"]) if e.get("family") else (
                mid[:-1] if re.match(r"^.*\d[A-Za-z]$", mid) else mid)
            m.tokens = [str(t) for t in e["path"]]
            if e.get("reference"):
                try:
                    ref = parse_expression(e["reference"])
                except ProjectError as exc:
                    self.errors.append("mechanism %s: reference: %s" % (mid, exc))
                    continue
                bad = [s for s in ref if s not in self.species]
                if bad:
                    self.errors.append("mechanism %s: reference: unknown species %s" % (mid, ", ".join(bad)))
                    continue
                m.reference = ref
            try:
                m.states = self.expand(m.tokens, mid, reference=m.reference)
            except ProjectError as exc:
                self.errors.append("mechanism %s: %s" % (mid, exc))
                m.states = []
                continue
            prod = e.get("product")
            if prod:
                m.product = str(prod)
            else:
                last = m.states[-1]
                main = [s for s in last.species if s not in self.auto_species] or list(last.species)
                m.product = expression(Counter({s: last.species[s] for s in main}))
            for s in m.states:
                if s.kind == "ts":
                    st = self.steps[s.step]
                    if mid not in st.used_by:
                        st.used_by.append(mid)
                for sid in s.full:
                    if sid in self.species and mid not in self.species[sid].used_by:
                        self.species[sid].used_by.append(mid)
            for st in self.steps.values():
                for sid in list(st.left) + list(st.right) + ([st.start] if st.start else []):
                    if sid in self.species and st.id not in self.species[sid].used_by:
                        self.species[sid].used_by.append(st.id)

    def expand(self, tokens, mid="", reference=None):
        """Turn a path (species and step IDs) into the full list of balanced states.

        A step is used forwards or backwards, whichever side matches the state before it. When both
        sides match (a spectator is the same species as one side of the step), the state written
        after the step decides. Species of the state before the step that the step does not use (e.g.
        a second molecule that has already reacted) are carried through the step unchanged (spectators)."""
        ref = self.reference if reference is None else reference
        if not ref:
            raise ProjectError("no reference state (set energies.reference)")
        auto = set(self.auto_for(ref))

        def core(c):
            return Counter({s: n for s, n in c.items() if s not in auto and n > 0})

        def contains(big, small):
            return all(big[s] >= n for s, n in small.items())

        def written_after(k):
            """Core of the species state written right after token k (None if a step or nothing follows)."""
            if k + 1 >= len(tokens) or tokens[k + 1] in self.steps:
                return None
            try:
                return core(parse_expression(tokens[k + 1]))
            except ProjectError:
                return None

        raw = [State("species", ref, label=expression(ref), source="reference")]
        for k, tok in enumerate(tokens):
            if tok in self.steps:
                st = self.steps[tok]
                if not st.start:
                    raise ProjectError("step %s has problems (see above)" % tok)
                prev = raw[-1]
                fwd = [st.left, st.right]
                spect = Counter()
                if prev.kind == "species":
                    pc, lc, rc = core(prev.species), core(st.left), core(st.right)
                    lin, rin = contains(pc, lc), contains(pc, rc)
                    nxt = written_after(k) if (lin and rin) else None
                    if lin and not rin:
                        spect = pc - lc
                    elif rin and not lin:
                        fwd, spect = [st.right, st.left], pc - rc
                    elif nxt is not None and ((pc - lc) + rc == nxt) != ((pc - rc) + lc == nxt):
                        if (pc - rc) + lc == nxt:                 # only the backward step gives the next state
                            fwd = [st.right, st.left]
                        spect = pc - core(fwd[0])
                    else:
                        try:
                            pf = prev.species + self.balance(prev.species, ref)
                            lf = st.left + self.balance(st.left, ref)
                            rf = st.right + self.balance(st.right, ref)
                            if pf == rf and pf != lf:
                                fwd = [st.right, st.left]
                        except BalanceError:
                            pass
                        if lin and rin:
                            spect = pc - core(fwd[0])
                a, b = fwd[0] + spect, fwd[1] + spect
                raw.append(State("species", a, label=expression(a), source="step"))
                raw.append(State("ts", Counter({st.start: 1}) + spect, step=tok, label=st.ts))
                raw.append(State("species", b, label=expression(b), source="step"))
            else:
                try:
                    c = parse_expression(tok)
                except ProjectError as exc:
                    raise ProjectError("path item '%s': %s" % (tok, exc))
                unknown = [s for s in c if s not in self.species]
                if unknown:
                    hint = " (not a species or step ID)" if len(c) == 1 else ""
                    raise ProjectError("path item '%s': unknown %s%s" % (tok, ", ".join(unknown), hint))
                raw.append(State("species", c, label=tok, source="path"))
        for s in raw:
            try:
                s.added = self.balance(s.species, ref)
            except BalanceError as exc:
                raise ProjectError("state '%s': %s" % (s.label, exc))
        # merge neighbours that are the same state (a path species written next to the step side it equals)
        out = []
        for s in raw:
            if out and s.same_as(out[-1]):
                if s.source == "path" and out[-1].source == "step":
                    out[-1] = s                      # keep the state as the user wrote it
                elif s.source == "path" and out[-1].source == "reference":
                    out[-1].label = s.label          # the starting state, labelled as written
                continue
            out.append(s)
        return out

    def state_main(self, state):
        """The species of a state that a map should show (without unreacted reactants and balancers)."""
        if state.kind == "ts":
            return Counter()
        main = Counter({s: n for s, n in state.species.items() if s not in self.auto_species})
        return main if main else Counter(state.species)

    # ------------------------------------------------------------------ selections
    def select(self, mech=None, steps=None, species=None, product=None):
        """Species and steps needed by the chosen mechanisms / steps / species (None = everything)."""
        mechs = list(self.mechanisms.values())
        if product:
            want = set(x.strip() for x in product.split(","))
            mechs = [m for m in mechs if m.product in want]
        if mech:
            mechs = [m for m in mechs if match_any(m.id, mech) or (m.family and match_any(m.family, mech, exact=True))]
        chosen_steps = list(self.steps.values())
        chosen_species = list(self.species.values())
        if mech or product:
            ids = set()
            sps = set()
            for m in mechs:
                for s in m.states:
                    if s.kind == "ts":
                        ids.add(s.step)
                    sps.update(s.full)
            chosen_steps = [st for st in chosen_steps if st.id in ids]
            for st in chosen_steps:
                sps.update(list(st.left) + list(st.right) + [st.start])
            chosen_species = [sp for sp in chosen_species if sp.id in sps]
        if steps:
            chosen_steps = [st for st in chosen_steps if match_any(st.id, steps)]
            sps = set()
            for st in chosen_steps:
                sps.update(list(st.left) + list(st.right) + [st.start])
            chosen_species = [sp for sp in chosen_species if sp.id in sps]
        if species:
            chosen_species = [sp for sp in chosen_species if match_any(sp.id, species)]
            if not steps and not mech and not product:
                chosen_steps = []                    # --species alone: only the Step 1 jobs of those species
        return chosen_species, chosen_steps, mechs

    def resolve_step(self, ref):
        ref = str(ref).rstrip("/")
        if ref in self.steps:
            return self.steps[ref]
        base = os.path.basename(os.path.abspath(ref))
        if base in self.steps:
            return self.steps[base]
        parts = os.path.abspath(ref).split(os.sep)
        if "Steps" in parts:
            k = len(parts) - 1 - parts[::-1].index("Steps")
            if k + 1 < len(parts) and parts[k + 1] in self.steps:
                return self.steps[parts[k + 1]]
        raise ProjectError("unknown step '%s'" % ref)


def match_any(ident, patterns, exact=False):
    """'S5' matches S5, S5a, S5b (not S50); 'S1-S9' is a range; patterns are comma separated."""
    for pat in str(patterns).split(","):
        pat = pat.strip()
        if not pat:
            continue
        m = re.match(r"^([A-Za-z]+)(\d+)\s*-\s*(?:[A-Za-z]+)?(\d+)$", pat)
        if m:
            mm = re.match(r"^([A-Za-z]+)(\d+)([^\d].*)?$", ident)
            if mm and mm.group(1) == m.group(1) and int(m.group(2)) <= int(mm.group(2)) <= int(m.group(3)):
                return True
            continue
        if ident == pat:
            return True
        if not exact and ident.startswith(pat) and not ident[len(pat)].isdigit():
            return True
    return False
