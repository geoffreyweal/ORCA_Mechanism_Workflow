"""Network maps of the mechanisms: Graphviz (dot / svg / pdf / png) and Mermaid.

Nodes are the species a mechanism goes through (unreacted reactants and balancing species
are left out); a solid arrow is a step (its label is the step ID, e.g. S5a); a dashed line is a
proton-transfer equilibrium (two species written next to each other in a path).
"""
import html
import os
import shutil
import subprocess
from collections import OrderedDict

from .project import expression

STATUS_COLOURS = {          # fill (nodes) / line (edges)
    "done": "#2e7d32", "check": "#ef6c00", "active": "#1565c0", "failed": "#c62828",
    "ready": "#6d6d6d", "waiting": "#9e9e9e", "absent": "#bdbdbd",
}
NODE_FILL = {"done": "#c8e6c9", "check": "#ffe0b2", "active": "#bbdefb", "failed": "#ffcdd2",
             "ready": "#eeeeee", "waiting": "#ffffff", "absent": "#ffffff"}
MECH_COLOURS = ["#1b5e20", "#0d47a1", "#b71c1c", "#e65100", "#4a148c", "#006064", "#827717", "#3e2723",
                "#880e4f", "#263238", "#33691e", "#01579b"]


class Graph:
    def __init__(self):
        self.nodes = OrderedDict()     # key -> {"label", "species", "kind", "mechs", "state"}
        self.edges = OrderedDict()     # (a, b, step) -> {"mechs", "step", "note", "state"}


def _direction_notes(project, st, prev_state, shown_before, shown_after):
    """Reactants / balancing species a step uses (+X) or releases (-X) that the two nodes do not show."""
    forward = True
    try:
        lf = st.left + project.balance(st.left)
        rf = st.right + project.balance(st.right)
        if prev_state is not None and prev_state.full == rf and prev_state.full != lf:
            forward = False
    except Exception:
        pass
    a, b = (st.left, st.right) if forward else (st.right, st.left)
    auto = set(project.auto_species)
    used = [s for s in a if s in auto and s not in shown_before]
    made = [s for s in b if s in auto and s not in shown_after]
    return " ".join(["+" + s for s in used] + ["-" + s for s in made])


def build(project, mechs):
    g = Graph()
    for m in mechs:
        prev_key, prev_state, ts = None, None, None
        for s in m.states:
            if s.kind == "ts":
                ts = s
                continue
            if s.source == "reference":
                key = "__start__"
                label, species = s.label, list(s.species)
                kind = "start"
            else:
                main = project.state_main(s)
                key = expression(main)
                label, species = key, list(main)
                kind = "product" if species and all(project.role(x) == "product" for x in species) else "intermediate"
            node = g.nodes.setdefault(key, {"label": label, "species": species, "kind": kind, "mechs": [], "state": s})
            if m.id not in node["mechs"]:
                node["mechs"].append(m.id)
            if prev_key is not None and (prev_key != key or ts is not None):
                step = ts.step if ts is not None else None
                note = _direction_notes(project, project.steps[step], prev_state, g.nodes[prev_key]["species"],
                                        species) if step else ""
                e = g.edges.setdefault((prev_key, key, step), {"mechs": [], "step": step, "note": note,
                                                              "state": ts})
                if m.id not in e["mechs"]:
                    e["mechs"].append(m.id)
            prev_key, prev_state, ts = key, s, None
    return g


# ----------------------------------------------------------------------------- status / energies

def species_status(result):
    if result["ok"]:
        return "done"
    st = result["state"]
    if st == "finished" or st == "failed":
        return "failed"
    if st in ("submitted", "running"):
        return "active"
    return {"ready": "ready", "absent": "absent"}.get(st, "waiting")


def step_status(stage, detail, action):
    if stage == "DONE":
        return "done"
    if stage == "CHECK":
        return "check"
    if stage == "waiting":
        return "waiting"
    d = detail.lower()
    if any(w in d for w in ("fail", "no energy maximum", "imaginary frequencies", "did not converge", "stopped")):
        return "failed"
    if action or d == "ready":
        return "ready"
    return "active"


def _node_status(node, sp_status):
    states = [sp_status.get(s, "absent") for s in node["species"]]
    for k in ("failed", "check", "active", "ready", "waiting", "absent"):
        if k in states:
            return k
    return "done" if states else "absent"


# ----------------------------------------------------------------------------- writers

def _short(text, n=34):
    text = str(text)
    return text if len(text) <= n else text[:n - 1] + "…"


def to_dot(project, g, title="", status=None, energies=None, colour_mechs=False, names=True, direction="LR"):
    """status = (species_status dict, step_status dict); energies = Energies object (optional)."""
    sp_status, st_status = status if status else ({}, {})
    mech_list = sorted(set(m for e in g.edges.values() for m in e["mechs"]))
    mcol = {m: MECH_COLOURS[k % len(MECH_COLOURS)] for k, m in enumerate(mech_list)}
    legend = ("solid arrow: step Sn (transition state TSn) &#160; dashed: proton-transfer equilibrium &#160; "
              "+X / -X: reactant or balancing species used / released")
    if status:
        legend += ("<br/>colours: green done, orange check, blue running or submitted, grey not run yet, "
                   "red needs attention")
    if energies is not None:
        legend += "<br/>numbers: G in %s relative to the start" % energies.units
    out = ["digraph mechanisms {",
           '  graph [rankdir=%s, fontname="Helvetica", fontsize=16, nodesep=0.18, ranksep=0.45, labelloc=t, '
           'label=<<b>%s</b><br/><font point-size="10">%s</font><br/> >];' % (direction, html.escape(title), legend),
           '  node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#424242", fontname="Helvetica", '
           'fontsize=11, margin="0.08,0.04"];',
           '  edge [fontname="Helvetica", fontsize=9, color="#424242", arrowsize=0.7];']
    ids = {k: "n%d" % i for i, k in enumerate(g.nodes)}
    g0 = None
    if energies is not None and "__start__" in g.nodes:
        g0 = energies.g_state(g.nodes["__start__"]["state"])
    for k, n in g.nodes.items():
        lines = ["<b>%s</b>" % html.escape(n["label"])]
        if names and len(n["species"]) == 1 and n["kind"] != "start":
            nm = project.species[n["species"][0]].name
            if nm:
                lines.append('<font point-size="8">%s</font>' % html.escape(_short(nm)))
        if energies is not None:
            v = energies.rel(energies.g_state(n["state"]), g0)
            if v is not None:
                lines.append('<font point-size="9" color="#1a237e">%.1f</font>' % v)
        a = OrderedDict([("label", "<%s>" % "<br/>".join(lines))])
        if n["kind"] == "start":
            a.update(style='"rounded,filled,bold"', fillcolor='"#fff9c4"', penwidth="2")
        elif n["kind"] == "product":
            a.update(peripheries="2", fillcolor='"#e1f5fe"')
        if status and n["kind"] != "start":
            s = _node_status(n, sp_status)
            a.update(fillcolor='"%s"' % NODE_FILL[s], color='"%s"' % STATUS_COLOURS[s])
        attrs = ["%s=%s" % kv for kv in a.items()]
        out.append('  %s [%s];' % (ids[k], ", ".join(attrs)))
    for (a, b, step), e in g.edges.items():
        attrs = []
        if step:
            st = project.steps[step]
            lab = step if not e["note"] else "%s %s" % (step, e["note"])
            if energies is not None:
                v = energies.rel(energies.g_state(e["state"]), g0)
                if v is not None:
                    lab += "\\n%s %.1f" % (st.ts, v)
            attrs.append('label="%s"' % lab.replace('"', "'"))
            attrs.append('tooltip="%s: %s"' % (step, (st.name or st.reaction_text()).replace('"', "'")))
            if status:
                s = st_status.get(step, "absent")
                attrs.append('color="%s", fontcolor="%s", penwidth=%s' % (STATUS_COLOURS[s], STATUS_COLOURS[s],
                                                                          "2.2" if s in ("done", "check") else "1.3"))
            elif colour_mechs and len(e["mechs"]) == 1:
                attrs.append('color="%s", fontcolor="%s"' % (mcol[e["mechs"][0]], mcol[e["mechs"][0]]))
        else:
            attrs.append('style=dashed, dir=both, arrowsize=0.5, color="#757575", tooltip="proton transfer (equilibrium)"')
        out.append("  %s -> %s [%s];" % (ids[a], ids[b], ", ".join(attrs)))
    out.append("}")
    return "\n".join(out) + "\n"


def _mm(text):
    return str(text).replace('"', "#quot;").replace("<", "#lt;").replace(">", "#gt;")


def to_mermaid(project, g, status=None, names=True):
    sp_status, st_status = status if status else ({}, {})
    ids = {k: "n%d" % i for i, k in enumerate(g.nodes)}
    out = ["flowchart LR"]
    for k, n in g.nodes.items():
        label = _mm(n["label"])
        if names and len(n["species"]) == 1 and n["kind"] != "start":
            nm = project.species[n["species"][0]].name
            if nm:
                label += "<br/><small>%s</small>" % _mm(_short(nm))
        shape = {"start": '(["%s"])', "product": '[["%s"]]'}.get(n["kind"], '["%s"]')
        out.append("  %s%s" % (ids[k], shape % label))
    links = []
    for (a, b, step), e in g.edges.items():
        if step:
            lab = step if not e["note"] else "%s %s" % (step, e["note"])
            out.append('  %s -->|"%s"| %s' % (ids[a], _mm(lab), ids[b]))
        else:
            out.append("  %s -.- %s" % (ids[a], ids[b]))
        links.append(step)
    out.append("  classDef start fill:#fff9c4,stroke:#424242,stroke-width:2px")
    out.append("  classDef product fill:#e1f5fe,stroke:#424242,stroke-width:2px")
    starts = [ids[k] for k, n in g.nodes.items() if n["kind"] == "start"]
    prods = [ids[k] for k, n in g.nodes.items() if n["kind"] == "product"]
    if starts:
        out.append("  class %s start" % ",".join(starts))
    if prods:
        out.append("  class %s product" % ",".join(prods))
    if status:
        for s in ("done", "check", "active", "failed", "ready"):
            out.append("  classDef %s fill:%s,stroke:%s" % (s, NODE_FILL[s], STATUS_COLOURS[s]))
        for k, n in g.nodes.items():
            if n["kind"] == "intermediate":
                s = _node_status(n, sp_status)
                if s in ("done", "check", "active", "failed", "ready"):
                    out.append("  class %s %s" % (ids[k], s))
        for idx, step in enumerate(links):
            if step:
                s = st_status.get(step, "absent")
                out.append("  linkStyle %d stroke:%s,stroke-width:%s" % (idx, STATUS_COLOURS[s],
                                                                         "3px" if s in ("done", "check") else "1.5px"))
    return "\n".join(out) + "\n"


def render(dotfile, formats):
    """Render a .dot file with Graphviz. Returns (written files, problem or None)."""
    exe = shutil.which("dot")
    if not exe:
        return [], "Graphviz 'dot' not found: only the .dot and Mermaid files were written"
    done = []
    for fmt in formats:
        target = os.path.splitext(dotfile)[0] + "." + fmt
        res = subprocess.run([exe, "-T" + fmt, dotfile, "-o", target], capture_output=True, text=True)
        if res.returncode != 0:
            return done, "dot failed for %s: %s" % (fmt, res.stderr.strip()[:300])
        done.append(target)
    return done, None
