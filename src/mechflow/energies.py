"""Step 5: Gibbs energies of species and transition states, energy profiles and product ratios."""
import csv
import math
import os
import traceback
from collections import Counter

from . import chem, routes
from .orca import parse_output
from .project import expression
from .workflow import (HARTREE_KCAL, HARTREE_KJ, OPTTS, sp_enabled, sp_energy, species_result, species_sp_dir,
                       step_progress, step_sp_dir)

LEVELS = ("dft", "sp")

R_KJ = 8.314462618e-3        # kJ / (mol K)
KB_H = 2.083661912e10        # Boltzmann constant / Planck constant, 1 / (K s)


def eyring(en, barrier):
    """Eyring rate constant k = (kB T / h) exp(-dG / RT), transmission coefficient 1, for a barrier in the
    project's units. In s-1 for a step that starts from one molecule; for two molecules it is second order
    (M-1 s-1 with standard_state: 1 M)."""
    if barrier is None:
        return None
    return KB_H * en.T * math.exp(max(min(-barrier / en.RT, 700.0), -700.0))


class Energies:
    """Gibbs free energies (Eh) of every species and transition state that has finished.

    level 'dft': G of the DFT frequency calculation. level 'sp' (the default when single_point.enabled):
    G = E(single point, SP) + [G - E(el)] of the DFT frequency calculation, e.g. DLPNO-CCSD(T)//DFT; a species
    or TS whose single point has not finished has no energy yet."""

    level = "dft"

    def __init__(self, project, queue=None, level=None):
        self.project = project
        self.level = level or ("sp" if sp_enabled(project) else "dft")
        if self.level not in LEVELS:
            raise ValueError("level must be one of %s" % ", ".join(LEVELS))
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
        for sid, r in self.species.items():
            r["sp_state"], r["e_sp"] = sp_energy(species_sp_dir(project, sid), queue)
            r["g_used"] = self._combine(r["gibbs"] if r["ok"] else None, r["g_minus_eel"], r["e_sp"])
        self.steps = {}
        for sid, st in project.steps.items():
            stage, detail, _ = step_progress(project, st, queue)
            r3 = parse_output(os.path.join(project.step_dir(sid), OPTTS, "output.out"))
            ts_ok = r3["opt_converged"] and len(r3["imag"]) == 1
            sp_state, e_sp = sp_energy(step_sp_dir(project, sid), queue)
            self.steps[sid] = {"stage": stage, "detail": detail, "gibbs_dft": r3["gibbs"] if ts_ok else None,
                               "gibbs": self._combine(r3["gibbs"] if ts_ok else None, r3["g_minus_eel"], e_sp),
                               "e_sp": e_sp, "sp_state": sp_state,
                               "imag": r3["imag"][0] if len(r3["imag"]) == 1 else None,
                               "validated": stage == "DONE"}

    def _combine(self, g_dft, g_minus_eel, e_sp):
        """the free energy used at this level (Eh), or None"""
        if g_dft is None:
            return None
        if self.level == "dft":
            return g_dft
        if e_sp is None or g_minus_eel is None:
            return None
        return e_sp + g_minus_eel

    @property
    def level_text(self):
        if self.level == "dft":
            return "DFT free energies (%s)" % self.project.settings["orca"]["method"].lstrip("!").strip()
        sp = self.project.settings.get("single_point") or {}
        return ("E(%s) + [G - E(el)] of %s" % (str(sp.get("method", "")).lstrip("!").strip(),
                                                self.project.settings["orca"]["method"].lstrip("!").strip()))

    def g_species(self, sid):
        r = self.species.get(sid)
        g = r.get("g_used", r.get("gibbs")) if r else None
        return (g + self.per_molecule) if r and r["ok"] and g is not None else None

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
            # the TS replaces its step's start species; other species in the state are spectators
            spect = state.species - Counter({self.project.steps[state.step].start: 1})
            extra = self.g_counter(state.added + spect)
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


def collect(project, mechs=None, outdir=None, queue=None, proton_transfers=None, level=None):
    """Write results/*.csv and results/summary.md; return the summary text.

    `proton_transfers` ('auto', 'path' or 'steps') replaces energies.proton_transfers for the network search.
    `level` ('dft' or 'sp') picks the energies (default: 'sp' when single_point.enabled, else 'dft')."""
    p = project
    en = Energies(p, queue) if level is None else Energies(p, queue, level)
    outdir = outdir or p.path("results")
    os.makedirs(outdir, exist_ok=True)
    u = en.units
    mechs = list(p.mechanisms.values()) if mechs is None else mechs

    with open(os.path.join(outdir, "species.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "role", "name", "formula", "charge", "mult", "Step1", "usable", "G_Eh", "E_el_Eh",
                    "G_minus_E_el_Eh", "SP", "E_SP_Eh", "G_used_Eh", "note"])
        for sid, sp in p.species.items():
            r = en.species[sid]
            w.writerow([sid, p.role(sid), sp.name, chem.formula_string(sp.formula), sp.charge, sp.mult, r["state"],
                        r["ok"], _fmt(r["gibbs"], "%.8f"), _fmt(r["e_el"], "%.8f"), _fmt(r.get("g_minus_eel"), "%.8f"),
                        r.get("sp_state", ""), _fmt(r.get("e_sp"), "%.8f"), _fmt(r.get("g_used", r["gibbs"]), "%.8f"),
                        r["note"]])

    with open(os.path.join(outdir, "steps.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "TS", "name", "reaction", "stage", "imag_cm-1", "G_TS_Eh", "E_SP_TS_Eh",
                    "dG_act_forward_" + u, "dG_act_reverse_" + u, "dG_reaction_" + u, "k_forward", "k_reverse",
                    "scan_start", "dG_act_from_start_" + u, "k_from_start_s-1", "used_by", "detail"])
        for sid, st in p.steps.items():
            e = en.steps[sid]
            gl, gr = en.g_counter(st.left), en.g_counter(st.right)
            gts = None if e["gibbs"] is None else e["gibbs"] + en.per_molecule
            fwd, rev = en.rel(gts, gl), en.rel(gts, gr)
            from_start = en.rel(gts, en.g_species(st.start)) if st.start else None
            w.writerow([sid, st.ts, st.name, st.reaction_text(), e["stage"], _fmt(e["imag"]), _fmt(e["gibbs"], "%.8f"),
                        _fmt(e.get("e_sp"), "%.8f"),
                        _fmt(fwd), _fmt(rev), _fmt(en.rel(gr, gl)), _fmt(eyring(en, fwd), "%.3g"),
                        _fmt(eyring(en, rev), "%.3g"), st.start or "", _fmt(from_start),
                        _fmt(eyring(en, from_start), "%.3g"), " ".join(st.used_by), e["detail"]])

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
    prods = sorted(set(m.product for m, *_ in results.values() if m.product))

    def ratio(ga, gb):
        x = max(min((gb - ga) / en.RT, 700.0), -700.0)
        return math.exp(x)

    L = ["# Results: %s" % (p.settings["project"].get("title") or p.name), "",
         "Gibbs free energies in %s relative to %s (G = 0), at %.2f K, standard state %s." % (
             u, expression(p.reference), en.T, en.standard_state),
         "Each state is balanced with %s where needed." % (", ".join(p.auto_species) or "nothing"),
         "Energies: %s." % en.level_text, ""]
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
    try:
        L += routes.section(p, en, outdir, proton_transfers)
    except Exception as exc:          # a problem in the network search must not stop summary.md being written
        traceback.print_exc()
        L += ["", routes.HEADING, "", "The search stopped with an error: %s" % exc]
    L += ["", "Files: species.csv, steps.csv, profiles.csv, mechanisms.csv, families.csv, network_routes.csv "
          "(this folder).", ""]
    text = "\n".join(L)
    with open(os.path.join(outdir, "summary.md"), "w") as fh:
        fh.write(text)
    return text
