"""Best routes to each product through the whole reaction network (part of 'mechflow collect').

Every mechanism and every step of a project is put into one network:

    states       balanced states, named by their core species: the species that are not added
                 automatically (unreacted reference species and energies.balance_with), e.g. I5a or
                 I7a + I9a. The reference state has no core species.
    connections  - two neighbouring states of a mechanism path: through the TS between them, or a proton
                   transfer (no barrier) if there is none;
                 - every step, either way round, through its TS;
                 - an automatic proton transfer between two states that hold the same molecules apart from
                   where H atoms sit on N, O, S, P, Se or a halogen (see proton_fingerprint).

energies.proton_transfers (or 'mechflow collect --proton-transfers') says which barrier-free proton transfers
are used: 'auto' both kinds (the default), 'path' only those written in mechanism paths, 'steps' none, so
that every proton transfer has to go through a step and its TS (when the proton transfers are steps of
their own, e.g. X + HA -> XH+ + A-).

For each product, a minimax ("bottleneck") Dijkstra search finds the route from the reference whose highest
point (state or TS) is lowest. Routes with the same highest point are ranked by their next-highest point,
and so on. States and TSs without an energy cannot be passed.
"""
import csv
import heapq
import itertools
import math
import os
from collections import Counter, OrderedDict

from . import chem
from .defaults import PROTON_TRANSFER_MODES
from .project import ProjectError, State, expression, species_sort_key

#: H atoms on these elements are mobile protons; H atoms on any other element (C-H) are counted
MOBILE_H = ("N", "O", "S", "P", "Se", "F", "Cl", "Br", "I")
N_ROUTES = 10            # routes per product in network_routes.csv
N_ALTERNATIVES = 5       # alternative bottlenecks per product in summary.md
SAME = 0.05              # two energies (in the project's units) closer than this are the same
HEADING = "## Best routes through the whole network"


def proton_fingerprint(mol, tolerance=chem.BOND_TOLERANCE):
    """Connectivity fingerprint of a molecule that does not see where its mobile protons sit.

    A Weisfeiler-Lehman hash (chem.graph_hash) of the heavy atoms and the bonds between them; H atoms
    are left out. Each heavy atom is labelled by its element and, unless it is in MOBILE_H, by its number
    of H atoms. HA and A-, or an alcohol and its protonated form, therefore get the same fingerprint, a
    ketone and its enol (one C-H fewer) do not. Stereoisomers and conformers get the same fingerprint."""
    adj = chem.adjacency(len(mol), chem.bonds(mol, tolerance))
    heavy = [i for i, s in enumerate(mol.symbols) if s != "H"]
    if not heavy:
        return "H%d" % len(mol)
    labels = list(mol.symbols)
    for i in heavy:
        if labels[i] not in MOBILE_H:
            labels[i] += "H%d" % sum(1 for j in adj[i] if mol.symbols[j] == "H")
    return chem.graph_hash(labels, adj, heavy)


def _contains(big, small):
    return all(big.get(s, 0) >= n for s, n in small.items())


def _merge(key, points):
    """The energies met along a route, highest first. Comparing two of these tuples ranks routes by their
    highest point, then by their next-highest point, and so on; of two otherwise equal routes the one with
    fewer points comes first. Adding the same point to two routes keeps their order, so Dijkstra works."""
    return tuple(sorted(key + tuple(points), reverse=True))


def _back(pred, k):
    """Nodes and edges from node k back to the source of a search: [k, edge, node, ..., source]."""
    seq = [k]
    while pred[k] is not None:
        prev, edge = pred[k]
        seq += [edge, prev]
        k = prev
    return seq


def _fmt(v, f="%.1f"):
    return "" if v is None else f % v


def _count(n, word, many=None):
    return "%d %s" % (n, word if n == 1 else (many or word + "s"))


class Node:
    """A balanced state of the network, named by its core species."""

    def __init__(self, key, core, label):
        self.key = key
        self.core = core
        self.label = label
        self.variants = []         # (full state, reference, G, source) for every place the state was found
        self.full = None           # the balanced state whose energy and molecules are used
        self.ref = None            # ... and the reference it is balanced against
        self.energy = None
        self.cls = None            # proton-transfer class: fingerprints of every molecule in self.full
        self.mechs = []

    def __repr__(self):
        return "Node(%s)" % self.label


class Edge:
    """A connection: a step through its TS (kind 'ts'), or a proton transfer without a barrier ('pt')."""

    def __init__(self, index, a, b, step=None, label="PT", energy=None, left=None, auto=False):
        self.index = index
        self.a, self.b = a, b
        self.kind = "ts" if step else "pt"
        self.step = step
        self.label = label         # 'TS5a', or 'TS5a + I7a' when a spectator is carried through the step
        self.energy = energy
        self.left = left           # the node on the step's left-hand side (to tell the direction of use)
        self.auto = auto           # a proton transfer that no mechanism writes
        self.mechs = []            # mechanisms whose path contains this connection

    def __repr__(self):
        return "Edge(%s)" % (self.step or "PT")


class Route:
    """A route through the network: seq = [node key, edge, node key, edge, ..., node key]."""

    def __init__(self, net, seq):
        self.seq = seq
        nodes = net.nodes
        self.nodes = seq[::2]
        self.tokens, self.points, self.backwards, self.auto_pt, ts = [], [], [], [], []
        for i, x in enumerate(seq):
            if isinstance(x, Edge):
                if x.kind == "ts":
                    ts.append(x)
                    self.tokens.append(x.step)
                    self.points.append((x.label, x.energy, "ts"))
                    if x.left is not None and seq[i - 1] != x.left and x.step not in self.backwards:
                        self.backwards.append(x.step)
                elif x.auto:
                    self.auto_pt.append("%s <=> %s" % (nodes[seq[i - 1]].label, nodes[seq[i + 1]].label))
            else:
                self.tokens.append(nodes[x].label)
                self.points.append((nodes[x].label, nodes[x].energy, "state"))
        self.key = tuple(sorted((v for _, v, _ in self.points), reverse=True))
        # highest point, and the effective barrier as in energies.mechanism_profile
        self.highest = self.eff = lowest = None
        for label, v, kind in self.points:
            if self.highest is None or v > self.highest[1]:
                self.highest = (label, v)
            if kind == "ts" and lowest is not None and (self.eff is None or v - lowest[1] > self.eff[1]):
                self.eff = (label, v - lowest[1], lowest[0])
            if lowest is None or v < lowest[1]:
                lowest = (label, v)
        self.segments = net.segments(ts)
        refs = [nodes[k].ref for k in self.nodes if nodes[k].ref != net.project.reference]
        self.reference = expression(refs[0]) if refs else None

    @property
    def simple(self):
        return len(set(self.nodes)) == len(self.nodes)

    @property
    def path(self):
        return "[%s]" % ", ".join(self.tokens)

    @property
    def mechanisms(self):
        """'M2 (S4, S5) + M1 (S2)': the mechanisms whose steps the route uses, in order."""
        return " + ".join("%s (%s)" % (m or "no mechanism", ", ".join(steps)) for m, steps in self.segments)

    @property
    def notes(self):
        out = []
        if self.reference:
            out.append("needs reference %s" % self.reference)
        if self.backwards:
            out.append("%s used backwards" % ", ".join(self.backwards))
        if self.auto_pt:
            out.append("automatic proton transfer %s" % ", ".join(self.auto_pt))
        return "; ".join(out)


class ProductRoutes:
    """The routes found to one product."""

    def __init__(self, label, nodes):
        self.label = label
        self.nodes = nodes         # keys of the nodes that count as this product (e.g. P1a and 2 P1a)
        self.routes = []           # different routes, best first (routes[0] is the best route)
        self.alternatives = []     # the best route limited by each other highest point, lowest first

    @property
    def best(self):
        return self.routes[0] if self.routes else None


class Network:
    """Every mechanism and every step of a project as one graph of balanced states (see the module doc).

    `en` is an energies.Energies object (or anything with its g_state, g_counter and rel methods).
    `proton_transfers` ('auto', 'path' or 'steps') replaces the project's energies.proton_transfers."""

    def __init__(self, project, en, proton_transfers=None):
        self.project = p = project
        self.en = en
        self.mode = proton_transfers or getattr(p, "proton_transfers", "auto")
        if self.mode not in PROTON_TRANSFER_MODES:
            raise ProjectError("proton_transfers must be one of %s, not '%s'" % (
                ", ".join(PROTON_TRANSFER_MODES), self.mode))
        self.auto = set(p.auto_species)
        self.nodes = OrderedDict()         # core key -> Node
        self.edges = []
        self.adj = {}                      # core key -> [(edge, other core key)]
        self.notes = []
        self.unused_links = []             # (mechanism, state, state) joined in a path without a step ('steps')
        self._edge_of = {}
        self._fp = {}
        self.g_ref = en.g_counter(p.reference) if p.reference else None
        self.start = ()
        self._node(p.reference, p.reference, self._rel(self.g_ref, self.g_ref), "reference")
        self._add_mechanisms()
        self._add_steps()
        self._choose_states()
        if self.mode == "auto":
            self._add_proton_transfers()
        if self.unused_links:
            m, a, b = self.unused_links[0]
            self._note("%s in mechanism paths without a step %s not used, because proton transfers only go through "
                       "steps (proton_transfers: steps), e.g. %s: %s -> %s" % (
                           _count(len(self.unused_links), "link"), "is" if len(self.unused_links) == 1 else "are",
                           m, a, b))
        self.products = self._find_products()

    # ------------------------------------------------------------------ building
    def _rel(self, g, g0):
        v = self.en.rel(g, g0)
        return None if v is None else round(v, 6)

    def _note(self, text):
        if text not in self.notes:
            self.notes.append(text)

    def _core(self, counter):
        return Counter({s: n for s, n in counter.items() if s not in self.auto and n > 0})

    def _node(self, full, ref, energy, source):
        core = self._core(full)
        key = tuple(sorted(core.items()))
        node = self.nodes.get(key)
        if node is None:
            label = expression(core) if core else expression(self.project.reference)
            node = self.nodes[key] = Node(key, core, label)
            self.adj[key] = []
        node.variants.append((Counter(full), Counter(ref), energy, source))
        if source in self.project.mechanisms and source not in node.mechs:
            node.mechs.append(source)
        return node

    def _edge(self, x, y, step=None, energy=None, label="PT", mech=None, left=None, auto=False):
        if x.key == y.key:
            return None
        k = (frozenset((x.key, y.key)), step)
        e = self._edge_of.get(k)
        if e is None:
            e = Edge(len(self.edges), x.key, y.key, step, label, energy, left, auto)
            self.edges.append(e)
            self._edge_of[k] = e
            self.adj[x.key].append((e, y.key))
            self.adj[y.key].append((e, x.key))
        else:
            if e.left is None:
                e.left = left
            if e.energy is None:
                e.energy = energy
            elif energy is not None and abs(energy - e.energy) > SAME:
                self._note("%s between %s and %s: G differs between the places it appears (%.1f and %.1f); "
                           "%.1f is used" % (e.label, self.nodes[e.a].label, self.nodes[e.b].label,
                                             e.energy, energy, e.energy))
            e.auto = e.auto and auto
        if mech and mech not in e.mechs:
            e.mechs.append(mech)
        return e

    def _left(self, st, x, y):
        """The one of nodes x and y that holds the step's left-hand side (None if that cannot be told)."""
        lc, rc = self._core(st.left), self._core(st.right)
        xl, yr = _contains(x.core, lc), _contains(y.core, rc)
        yl, xr = _contains(y.core, lc), _contains(x.core, rc)
        if xl and yr and not (yl and xr):
            return x.key
        if yl and xr and not (xl and yr):
            return y.key
        return None

    def _add_mechanisms(self):
        """Neighbouring states of every mechanism path, with energies relative to the mechanism's reference."""
        p, en = self.project, self.en
        for m in p.mechanisms.values():
            if not m.states:
                continue
            if m.reference is not None and not _contains(m.reference, p.reference):
                self._note("mechanism %s is left out: its reference %s does not contain the project's "
                           "reference %s" % (m.id, expression(m.reference), expression(p.reference)))
                continue
            ref = p.reference if m.reference is None else m.reference
            g0 = en.g_state(m.states[0])
            prev, ts = None, None
            for s in m.states:
                v = self._rel(en.g_state(s), g0)
                if s.kind == "ts":
                    ts = (s, v)
                    continue
                node = self._node(s.full, ref, v, m.id)
                if prev is not None and ts is None and prev.key != node.key:
                    if self.mode == "steps":
                        self.unused_links.append((m.id, prev.label, node.label))
                    else:
                        self._edge(prev, node, mech=m.id)
                elif prev is not None:
                    st = p.steps[ts[0].step]
                    spect = self._core(ts[0].species - Counter({st.start: 1}))
                    label = st.ts + (" + " + expression(spect) if spect else "")
                    self._edge(prev, node, st.id, ts[1], label, m.id, self._left(st, prev, node))
                prev, ts = node, None

    def _add_steps(self):
        """Every step between its two sides, balanced against the project's reference."""
        p, en = self.project, self.en
        for st in p.steps.values():
            if not (st.start and st.left and st.right):
                continue
            sides = []
            for side in (st.left, st.right):
                try:
                    full = side + p.balance(side)
                except ProjectError:
                    sides.append(None)
                    continue
                sides.append(self._node(full, p.reference, self._rel(en.g_counter(full), self.g_ref),
                                        "step " + st.id))
            if None in sides:
                continue
            ts = State("ts", Counter({st.start: 1}), step=st.id, label=st.ts)
            try:
                ts.added = p.balance(ts.species)
            except ProjectError:
                continue
            self._edge(sides[0], sides[1], st.id, self._rel(en.g_state(ts), self.g_ref), st.ts, left=sides[0].key)

    def _choose_states(self):
        """Pick the balanced state of each node: preferably one balanced against the project's reference."""
        p, en = self.project, self.en
        for node in self.nodes.values():
            canon = next((v for v in node.variants if v[1] == p.reference), None)
            if canon is None:
                try:
                    full = node.core + p.balance(node.core)
                except ProjectError:
                    canon = node.variants[0]
                else:
                    canon = (full, Counter(p.reference), self._rel(en.g_counter(full), self.g_ref), "balanced")
                    node.variants.append(canon)
            node.full, node.ref = canon[0], canon[1]
            known = [v[2] for v in node.variants if v[2] is not None]
            node.energy = canon[2] if canon[2] is not None else (known[0] if known else None)
            if known and max(known) - min(known) > SAME:
                self._note("%s: G differs between the places it appears (%.1f to %.1f); %.1f is used"
                           % (node.label, min(known), max(known), node.energy))

    def fingerprint(self, sid):
        if sid not in self._fp:
            sp = self.project.species.get(sid)
            self._fp[sid] = None if sp is None or sp.mol is None else proton_fingerprint(sp.mol, self.project.tolerance)
        return self._fp[sid]

    def _add_proton_transfers(self):
        """Join states whose molecules have the same proton fingerprints (balancing species included)."""
        classes = OrderedDict()
        for node in self.nodes.values():
            fps = [self.fingerprint(s) for s, n in node.full.items() for _ in range(n)]
            if not fps or None in fps:
                continue
            node.cls = tuple(sorted(fps))
            classes.setdefault(node.cls, []).append(node)
        for members in classes.values():
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    self._edge(members[i], members[j], auto=True)

    def _find_products(self):
        """{product label: [node keys]}: states whose core species are all products, grouped by those species."""
        p = self.project
        groups = OrderedDict()
        for node in self.nodes.values():
            if node.core and all(p.role(s) == "product" for s in node.core):
                groups.setdefault(frozenset(node.core), []).append(node.key)
        out = OrderedDict()
        for g in sorted(groups, key=lambda g: sorted(species_sort_key(s) for s in g)):
            out[" + ".join(sorted(g, key=species_sort_key))] = groups[g]
        return out

    def segments(self, ts):
        """The TS edges of a route split into runs taken from one mechanism: [(mechanism or None, [steps])]."""
        out, i = [], 0
        while i < len(ts):
            best, n = None, 1
            for m in ts[i].mechs:
                k = 1
                while i + k < len(ts) and m in ts[i + k].mechs:
                    k += 1
                if best is None or k > n:
                    best, n = m, k
            if best is None:
                while i + n < len(ts) and not ts[i + n].mechs:
                    n += 1
            out.append((best, [e.step for e in ts[i:i + n]]))
            i += n
        return out

    # ------------------------------------------------------------------ searching
    def search(self, sources, stop=()):
        """Minimax Dijkstra from `sources` (node keys). Nodes in `stop` are reached but not left.

        Returns ({node: key}, {node: (previous node, edge) or None}); a key is the tuple of energies met on
        the best route to the node, highest first (see _merge)."""
        best, pred, done = {}, {}, set()
        heap, tie = [], itertools.count()
        for k in sources:
            e = self.nodes[k].energy
            if e is not None and k not in best:
                best[k], pred[k] = (e,), None
                heapq.heappush(heap, ((e,), next(tie), k))
        while heap:
            key, _, u = heapq.heappop(heap)
            if u in done:
                continue
            done.add(u)
            if u in stop:
                continue
            for edge, v in self.adj[u]:
                ev = self.nodes[v].energy
                if v in done or ev is None:
                    continue
                if edge.kind == "ts":
                    if edge.energy is None:
                        continue
                    nk = _merge(key, (edge.energy, ev))
                else:
                    nk = _merge(key, (ev,))
                if v not in best or nk < best[v]:
                    best[v], pred[v] = nk, (u, edge)
                    heapq.heappush(heap, (nk, next(tie), v))
        return best, pred

    def routes_to(self, label, n_routes=N_ROUTES, n_alternatives=N_ALTERNATIVES):
        """The best route to a product, the next-best different routes, and the alternative bottlenecks.

        For every TS edge u-v, the best route through it is the best route from the start to u, the
        edge, and the best route from v to the product (a search from the product). Routes that would
        pass a state twice are left out."""
        group = self.products[label]
        res = ProductRoutes(label, group)
        stop = set(group)
        to, pred_to = self.search([self.start], stop)
        reached = [k for k in group if k in to]
        if not reached:
            return res
        best = Route(self, _back(pred_to, min(reached, key=lambda k: to[k]))[::-1])
        res.routes.append(best)
        fro, pred_fro = self.search(group, {self.start})
        cands = []
        for e in self.edges:
            if e.kind != "ts" or e.energy is None:
                continue
            for u, v in ((e.a, e.b), (e.b, e.a)):
                if u in stop or u not in to or v not in fro:
                    continue
                cands.append((_merge(to[u], (e.energy,) + fro[v]), e.index, u, v, e))
        cands.sort(key=lambda c: c[:2])
        seen, tops = {tuple(best.tokens)}, {best.highest[0]}
        for _, _, u, v, e in cands:
            if len(res.routes) >= n_routes and len(res.alternatives) >= n_alternatives:
                break
            r = Route(self, _back(pred_to, u)[::-1] + [e] + _back(pred_fro, v))
            if not r.simple or tuple(r.tokens) in seen:
                continue
            seen.add(tuple(r.tokens))
            if len(res.routes) < n_routes:
                res.routes.append(r)
            if r.highest[0] not in tops and len(res.alternatives) < n_alternatives:
                tops.add(r.highest[0])
                res.alternatives.append(r)
        return res

    def all_routes(self):
        return [self.routes_to(label) for label in self.products]


def deciding(ra, rb):
    """Highest TS of each of two routes after the point where they part, as energies.deciding does for
    mechanisms. Returns ((label, G) or None, (label, G) or None, index of the first different item)."""
    k = 0
    while k < min(len(ra.seq), len(rb.seq)) and ra.seq[k] == rb.seq[k]:
        k += 1

    def top(r):
        best = None
        for x in r.seq[k:]:
            if isinstance(x, Edge) and x.kind == "ts" and (best is None or x.energy > best[1]):
                best = (x.label, x.energy)
        return best
    return top(ra), top(rb), k


def _ratio(en, ga, gb):
    return math.exp(max(min((gb - ga) / en.RT, 700.0), -700.0))


def section(project, en, outdir, proton_transfers=None):
    """Search the whole network, write network_routes.csv to `outdir` and return the summary.md lines."""
    u = en.units
    net = Network(project, en, proton_transfers)
    results = net.all_routes()
    with open(os.path.join(outdir, "network_routes.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["product", "rank", "bottleneck_TS", "G_highest_point_" + u, "effective_barrier_" + u,
                    "rate_limiting_TS", "from_state", "mechanisms", "notes", "route"])
        for res in results:
            for rank, r in enumerate(res.routes, 1):
                eff = r.eff or ("", None, "")
                w.writerow([res.label, rank, r.highest[0], _fmt(r.highest[1]), _fmt(eff[1]), eff[0], eff[2],
                            r.mechanisms, r.notes, ", ".join(r.tokens)])
    L = ["", HEADING, ""]
    found = [res for res in results if res.best]
    if not results:
        return L + ["No product states (species with role product) in the network yet."]
    if not found:
        return L + ["Nothing to report yet: no route from the reference to a product has an energy for every "
                    "state and TS on it."]
    n_ts = sum(1 for e in net.edges if e.kind == "ts")
    n_pt = sum(1 for e in net.edges if e.kind == "pt")
    n_auto = sum(1 for e in net.edges if e.auto)
    size = "%s and %s through a TS" % (_count(len(net.nodes), "state"), _count(n_ts, "connection"))
    if n_pt:
        size = "%s, %s through a TS and %s without a barrier (%d added automatically, between states with the " \
               "same molecules apart from where H atoms sit on N, O, S, P or a halogen)" % (
                   _count(len(net.nodes), "state"), _count(n_ts, "connection"), _count(n_pt, "proton transfer"),
                   n_auto)
    mode = {"auto": "",
            "path": " Proton transfers without a TS are used only where a mechanism path has one (proton_transfers: "
                    "path).",
            "steps": " Proton transfers only go through steps and their TSs (proton_transfers: steps)."}[net.mode]
    L += ["Every mechanism and every step joined into one network: %s. Steps can be used either way round.%s "
          "The best route to a product is the one whose highest point (state or TS) is lowest. G in %s relative "
          "to %s. Up to %d different routes per product are in network_routes.csv." % (
              size, mode, u, expression(project.reference), N_ROUTES), "",
          "| product | highest point | G(highest point) | effective barrier | steps from |", "|---|---|---|---|---|"]
    for res in results:
        r = res.best
        if r is None:
            L.append("| %s | - | no route with known energies yet | - | - |" % res.label)
            continue
        eff = "%s (%s from %s)" % (_fmt(r.eff[1]), r.eff[0], r.eff[2]) if r.eff else "-"
        L.append("| %s | %s | %s | %s | %s |" % (res.label, r.highest[0], _fmt(r.highest[1]), eff, r.mechanisms or "-"))
    L += ["", "Best routes, written as mechanism paths:", ""]
    for res in found:
        L.append("- %s: `%s`%s" % (res.label, res.best.path, (" (%s)" % res.best.notes) if res.best.notes else ""))
    if len(found) >= 2:
        L += ["", "Product ratios (Curtin-Hammett, %.2f K) from the highest TS of each best route after the point "
              "where the two routes part:" % en.T, ""]
        for i in range(len(found)):
            for j in range(i + 1, len(found)):
                a, b = found[i], found[j]
                ta, tb, _ = deciding(a.best, b.best)
                if ta is None or tb is None:
                    L.append("- %s : %s cannot be estimated (no transition state after the branch point)"
                             % (a.label, b.label))
                    continue
                L.append("- %s : %s = %.3g : 1  (%s vs %s: difference %s %s)" % (
                    a.label, b.label, _ratio(en, ta[1], tb[1]), ta[0], tb[0], _fmt(tb[1] - ta[1]), u))
    alts = [(res, r) for res in found for r in res.alternatives]
    if alts:
        L += ["", "Alternative bottlenecks: the best routes to each product whose highest point is another TS or "
              "state (up to %d per product, lowest first):" % N_ALTERNATIVES, "",
              "| product | highest point | G(highest point) | above the best route | steps from |",
              "|---|---|---|---|---|"]
        for res, r in alts:
            L.append("| %s | %s | %s | %+.1f | %s |" % (res.label, r.highest[0], _fmt(r.highest[1]),
                                                       r.highest[1] - res.best.highest[1], r.mechanisms or "-"))
    if net.notes:
        L += ["", "Notes on the network:", ""] + ["- " + n for n in net.notes[:10]]
        if len(net.notes) > 10:
            L.append("- ... and %d more" % (len(net.notes) - 10))
    return L
