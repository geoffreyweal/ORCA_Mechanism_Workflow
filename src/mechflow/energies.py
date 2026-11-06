"""Step 5: Gibbs energies of species and transition states, energy profiles and product ratios."""
import csv
import math
import os

from . import chem
from .orca import parse_output
from .project import expression
from .workflow import HARTREE_KCAL, HARTREE_KJ, OPTTS, species_result, step_progress

R_KJ = 8.314462618e-3        # kJ / (mol K)


class Energies:
    """Gibbs free energies (Eh) of every species and transition state that has finished."""

    def __init__(self, project, queue=None):
        self.project = project
        self.units = project.settings["energies"].get("units", "kJ/mol")
        self.factor = HARTREE_KCAL if "kcal" in self.units.lower() else HARTREE_KJ
        self.T = float(project.settings["energies"].get("temperature", 298.15))
        self.RT = R_KJ * self.T * (self.factor / HARTREE_KJ)
        # ORCA's G is for 1 atm. In solution a 1 mol/L standard state adds RT ln(RT/p) per molecule
        # (+7.9 kJ/mol at 298 K), which matters when two states have different numbers of molecules.
        ss = str(project.settings["energies"].get("standard_state", "1 atm")).replace(" ", "").lower()
        self.standard_state = "1 M" if ss in ("1m", "1mol/l", "1mol/dm3") else "1 atm"
        self.per_molecule = (R_KJ * self.T * math.log(0.0820573661 * self.T) / HARTREE_KJ
                             if self.standard_state == "1 M" else 0.0)
        self.species = {sid: species_result(project, sp, queue) for sid, sp in project.species.items()}
        self.steps = {}
        for sid, st in project.steps.items():
            stage, detail, _ = step_progress(project, st, queue)
            r3 = parse_output(os.path.join(project.step_dir(sid), OPTTS, "output.out"))
            g = r3["gibbs"] if (r3["opt_converged"] and len(r3["imag"]) == 1) else None
            self.steps[sid] = {"stage": stage, "detail": detail, "gibbs": g,
                               "imag": r3["imag"][0] if len(r3["imag"]) == 1 else None,
                               "validated": stage == "DONE"}

    def g_species(self, sid):
        r = self.species.get(sid)
        return (r["gibbs"] + self.per_molecule) if r and r["ok"] and r["gibbs"] is not None else None

    def g_counter(self, counter):
        tot = 0.0
        for sid, n in counter.items():
            g = self.g_species(sid)
            if g is None:
                return None
            tot += n * g
        return tot

    def g_state(self, state):
        if state.kind == "ts":
            g = self.steps[state.step]["gibbs"]
            extra = self.g_counter(state.added)
            return None if g is None or extra is None else g + self.per_molecule + extra
        return self.g_counter(state.full)

    def rel(self, a, b):
        return None if a is None or b is None else (a - b) * self.factor


def _fmt(v, f="%.1f"):
    return "" if v is None else f % v


def mechanism_profile(en, mech):
    """[(state, G_rel)], highest TS, effective barrier, complete, validated"""
    rows = []
    g0 = en.g_state(mech.states[0]) if mech.states else None
    for s in mech.states:
        rows.append((s, en.rel(en.g_state(s), g0)))
    complete = bool(rows) and all(v is not None for _, v in rows)
    validated = all(en.steps[s.step]["validated"] for s, _ in rows if s.kind == "ts")
    highest, eff = None, None
    if complete:
        lowest = None
        for s, v in rows:
            if s.kind == "ts":
                if highest is None or v > highest[1]:
                    highest = (s.label, v)
                if lowest is not None and (eff is None or v - lowest[1] > eff[1]):
                    eff = (s.label, v - lowest[1], lowest[0])
            if lowest is None or v < lowest[1]:
                lowest = (s.label, v)
    return rows, highest, eff, complete, validated


def _same_state(a, b):
    return a.kind == b.kind and a.step == b.step and a.full == b.full


def deciding(rows_a, rows_b):
    """Highest TS of each route AFTER the point where the two routes part (shared steps cannot decide
    between them). Returns ((label, G) or None, (label, G) or None, index of the first different state)."""
    k = 0
    while k < min(len(rows_a), len(rows_b)) and _same_state(rows_a[k][0], rows_b[k][0]):
        k += 1

    def top(rows):
        best = None
        for s, v in rows[k:]:
            if s.kind == "ts" and v is not None and (best is None or v > best[1]):
                best = (s.label, v)
        return best
    return top(rows_a), top(rows_b), k


def collect(project, mechs=None, outdir=None, queue=None):
    """Write results/*.csv and results/summary.md; return the summary text."""
    p = project
    en = Energies(p, queue)
    outdir = outdir or p.path("results")
    os.makedirs(outdir, exist_ok=True)
    u = en.units
    mechs = list(p.mechanisms.values()) if mechs is None else mechs

    with open(os.path.join(outdir, "species.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "role", "name", "formula", "charge", "mult", "Step1", "usable", "G_Eh", "E_el_Eh", "note"])
        for sid, sp in p.species.items():
            r = en.species[sid]
            w.writerow([sid, p.role(sid), sp.name, chem.formula_string(sp.formula), sp.charge, sp.mult, r["state"],
                        r["ok"], _fmt(r["gibbs"], "%.8f"), _fmt(r["e_el"], "%.8f"), r["note"]])

    with open(os.path.join(outdir, "steps.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "TS", "name", "reaction", "stage", "imag_cm-1", "G_TS_Eh",
                    "dG_act_forward_" + u, "dG_act_reverse_" + u, "dG_reaction_" + u, "used_by", "detail"])
        for sid, st in p.steps.items():
            e = en.steps[sid]
            gl, gr = en.g_counter(st.left), en.g_counter(st.right)
            gts = None if e["gibbs"] is None else e["gibbs"] + en.per_molecule
            w.writerow([sid, st.ts, st.name, st.reaction_text(), e["stage"], _fmt(e["imag"]), _fmt(e["gibbs"], "%.8f"),
                        _fmt(en.rel(gts, gl)), _fmt(en.rel(gts, gr)), _fmt(en.rel(gr, gl)),
                        " ".join(st.used_by), e["detail"]])

    results, profile_rows = {}, {}
    with open(os.path.join(outdir, "profiles.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mechanism", "product", "index", "kind", "state", "added_to_balance", "G_rel_" + u])
        for m in mechs:
            rows, highest, eff, complete, validated = mechanism_profile(en, m)
            results[m.id] = (m, highest, eff, complete, validated)
            profile_rows[m.id] = rows
            for k, (s, v) in enumerate(rows):
                w.writerow([m.id, m.product, k, s.kind, s.label, expression(s.added) if s.added else "", _fmt(v)])

    with open(os.path.join(outdir, "mechanisms.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mechanism", "family", "product", "highest_TS", "G_highest_TS_" + u, "rate_limiting_TS",
                    "effective_barrier_" + u, "from_state", "complete", "all_steps_validated"])
        for mid, (m, hi, eff, comp, val) in results.items():
            w.writerow([mid, m.family, m.product, hi[0] if hi else "", _fmt(hi[1] if hi else None),
                        eff[0] if eff else "", _fmt(eff[1] if eff else None), eff[2] if eff else "", comp, val])

    # best route to each product (Curtin-Hammett: the lowest highest-TS, all relative to the same reference)
    best = {}
    for mid, (m, hi, eff, comp, val) in results.items():
        if comp and hi and (m.product not in best or hi[1] < best[m.product][1]):
            best[m.product] = (mid, hi[1], hi[0], val)
    prods = sorted(set(m.product for m, *_ in results.values()))

    def ratio(ga, gb):
        x = max(min((gb - ga) / en.RT, 700.0), -700.0)
        return math.exp(x)

    L = ["# Results: %s" % (p.settings["project"].get("title") or p.name), "",
         "Gibbs free energies in %s relative to %s (G = 0), at %.2f K, standard state %s." % (
             u, expression(p.reference), en.T, en.standard_state),
         "Each state is balanced with %s where needed." % (", ".join(p.auto_species) or "nothing"), ""]
    L += ["## Mechanisms", "",
          "| mechanism | product | highest TS | G(highest TS) | rate-limiting TS | effective barrier | all steps validated |",
          "|---|---|---|---|---|---|---|"]
    for mid, (m, hi, eff, comp, val) in results.items():
        L.append("| %s | %s | %s | %s | %s | %s | %s |" % (
            mid, m.product, hi[0] if hi else "-", _fmt(hi[1]) if hi else ("incomplete" if not comp else "-"),
            eff[0] if eff else "-", ("%s (from %s)" % (_fmt(eff[1]), eff[2])) if eff else "-", "yes" if val else "no"))
    L += ["", "## Products", "", "Best (lowest highest-TS) complete mechanism to each product:", "",
          "| product | mechanism | highest TS | G(highest TS) |", "|---|---|---|---|"]
    for pr in prods:
        if pr in best:
            mid, g, ts, val = best[pr]
            L.append("| %s | %s | %s | %s%s |" % (pr, mid, ts, _fmt(g), "" if val else " (not all steps validated)"))
        else:
            L.append("| %s | - | - | no complete mechanism yet |" % pr)
    ready = [pr for pr in prods if pr in best]
    if len(ready) >= 2:
        L += ["", "Product ratios (Curtin-Hammett, %.2f K) from the highest TS of each best mechanism after the "
              "point where the two mechanisms part (steps they share cannot decide between them):" % en.T, ""]
        for a_ in range(len(ready)):
            for b_ in range(a_ + 1, len(ready)):
                a, b = ready[a_], ready[b_]
                ta, tb, _ = deciding(profile_rows[best[a][0]], profile_rows[best[b][0]])
                if ta is None or tb is None:
                    L.append("- %s : %s cannot be estimated (no transition state after the branch point)" % (a, b))
                    continue
                L.append("- %s : %s = %.3g : 1  (%s %s vs %s %s: difference %s %s)" % (
                    a, b, ratio(ta[1], tb[1]), best[a][0], ta[0], best[b][0], tb[0], _fmt(tb[1] - ta[1]), u))
    fams = {}
    for mid, (m, hi, eff, comp, val) in results.items():
        fams.setdefault(m.family, []).append((mid, m.product, comp))
    multi = {f: v for f, v in fams.items() if len(set(x[1] for x in v)) > 1}
    with open(os.path.join(outdir, "families.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["family", "mechanism_a", "product_a", "deciding_TS_a", "G_deciding_TS_a", "mechanism_b",
                    "product_b", "deciding_TS_b", "G_deciding_TS_b", "ddG_b_minus_a_" + u, "ratio_a_to_b"])
        if multi:
            L += ["", "## Within each mechanism family", "",
                  "Deciding TS = the highest TS after the two mechanisms part.", "",
                  "| family | mechanisms | deciding TS | G(deciding TS) | difference | ratio |",
                  "|---|---|---|---|---|---|"]
        for f, v in multi.items():
            v = sorted(v)
            for a_ in range(len(v)):
                for b_ in range(a_ + 1, len(v)):
                    (ma, pa, ca), (mb, pb, cb) = v[a_], v[b_]
                    ta, tb, _ = deciding(profile_rows[ma], profile_rows[mb])
                    ga = ta[1] if (ta and ca) else None
                    gb = tb[1] if (tb and cb) else None
                    dd = None if ga is None or gb is None else gb - ga
                    r = None if dd is None else ratio(ga, gb)
                    w.writerow([f, ma, pa, ta[0] if ta else "", _fmt(ga), mb, pb, tb[0] if tb else "", _fmt(gb),
                                _fmt(dd), _fmt(r, "%.3g")])
                    L.append("| %s | %s (%s) vs %s (%s) | %s / %s | %s / %s | %s | %s |" % (
                        f, ma, pa, mb, pb, ta[0] if ta else "-", tb[0] if tb else "-", _fmt(ga) or "-",
                        _fmt(gb) or "-", _fmt(dd) or "-",
                        ("%s : %s = %s : 1" % (pa, pb, _fmt(r, "%.3g"))) if r is not None else "-"))
    L += ["", "Files: species.csv, steps.csv, profiles.csv, mechanisms.csv, families.csv (this folder).", ""]
    text = "\n".join(L)
    with open(os.path.join(outdir, "summary.md"), "w") as fh:
        fh.write(text)
    return text
