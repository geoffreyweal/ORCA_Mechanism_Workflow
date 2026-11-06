"""overview.md: every species, step and mechanism of a project in one readable page."""
from . import chem, network
from .project import expression


def plural(n, word, many=None):
    return "%d %s" % (n, word if n == 1 else (many or word + "s"))


def _md(text):
    return str(text).replace("|", "\\|")


def path_text(project, mech):
    """'R1 + R2 -> I1 =[S1]=> I2 <=> I3 ...' using the species that matter in each state."""
    parts, prev_species = [], False
    for s in mech.states:
        if s.kind == "ts":
            parts.append(" —%s→ " % project.steps[s.step].id)
            prev_species = False
            continue
        lab = s.label if s.source == "reference" else expression(project.state_main(s))
        if prev_species:
            parts.append(" ⇌ ")
        parts.append("**%s**" % lab if project.role(lab) == "product" else lab)
        prev_species = True
    return "".join(parts)


def overview(project, mermaid=True):
    p = project
    o = p.settings["orca"]
    solvent = [ln for ln in str(o.get("solvent_block") or "").splitlines() if ln.strip() and not ln.strip().startswith("#")]
    L = ["# %s" % (p.settings["project"].get("title") or p.name), "",
         "Project `%s`: %s, %s (one TS search each), %s." % (
             p.name, plural(len(p.species), "species", "species"), plural(len(p.steps), "step"),
             plural(len(p.mechanisms), "mechanism")), "",
         "- Method: `%s`" % o["method"].strip(),
         "- Solvent: `%s`" % (solvent[0].strip() if solvent else "none"),
         "- Atom numbers below count from %d%s." % (p.index_base, " (as in ORCA inputs)" if p.index_base == 0 else ""),
         "- Names: R reactants, I intermediates, P products, C catalysts / species used during the reaction, "
         "S steps; step Sn goes through transition state TSn.",
         "- Energy reference (G = 0): %s; states are balanced with %s." % (
             expression(p.reference), ", ".join(p.auto_species) or "nothing"), ""]
    L += ["## Species", "", "| ID | role | name | formula | charge | mult | atoms | used by |",
          "|---|---|---|---|---|---|---|---|"]
    for sid, sp in p.species.items():
        mechs = [u for u in sp.used_by if u in p.mechanisms]
        steps = [u for u in sp.used_by if u in p.steps]
        used = ", ".join(steps) if len(mechs) == 0 else "%s%s" % (
            plural(len(mechs), "mechanism"), (", " + ", ".join(steps)) if steps else "")
        L.append("| %s | %s | %s | %s | %s | %d | %d | %s |" % (
            sid, p.role(sid), _md(sp.name), chem.formula_string(sp.formula),
            "%+d" % sp.charge if sp.charge else "0", sp.mult, sp.natoms,
            _md(used) or "-"))
    L += ["", "## Steps", "",
          "SCAN: the distance between two atoms of the `from` structure is scanned from its value there "
          "(measured in the species xyz file below; re-measured after Step 1) to `to`.", "",
          "| step | TS | what happens | reaction | SCAN | used by |", "|---|---|---|---|---|---|"]
    for sid, st in p.steps.items():
        if st.atoms:
            mol = p.species[st.start].mol
            scan = "from %s: %s–%s %.2f → %.2f Å (%s)" % (
                st.start, p.atom_label(mol, st.atoms[0]), p.atom_label(mol, st.atoms[1]), st.d0, st.to,
                "stretch" if st.mode == "stretch" else "together")
        else:
            scan = "(not set)"
        L.append("| %s | %s | %s | %s | %s | %s |" % (sid, st.ts, _md(st.name), _md(st.reaction_text()), scan,
                                                     ", ".join(st.used_by) or "-"))
    L += ["", "## Mechanisms", "",
          "—Sn→ is a step (transition state TSn); ⇌ is a proton-transfer equilibrium. "
          "Unreacted reactants and balancing species are not shown.", ""]
    for mid, m in p.mechanisms.items():
        L.append("- **%s** %s (product **%s**)  " % (mid, _md(m.name), m.product))
        L.append("  %s" % path_text(p, m))
    if mermaid and p.mechanisms:
        g = network.build(p, list(p.mechanisms.values()))
        L += ["", "## Network map", "", "```mermaid", network.to_mermaid(p, g).rstrip(), "```"]
    return "\n".join(L) + "\n"
