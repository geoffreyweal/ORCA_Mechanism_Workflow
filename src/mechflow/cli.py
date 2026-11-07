"""mechflow command line.

    mechflow init                  write project.yaml for the xyz files in species/
    mechflow validate              check the project, write overview.md
    mechflow next                  write every input that can be written now (Step 1 ... Step 4B)
    mechflow submit                write submit_ready.sh for the prepared jobs (--run: sbatch now)
    mechflow status                where every species and step is up to
    mechflow collect               Step 5: energies, profiles, product ratios -> results/
    mechflow map                   network maps (Graphviz / Mermaid) -> maps/

Run it inside a project folder (or give -p FOLDER). See https://geoffreyweal.github.io/ORCA_Mechanism_Workflow/
"""
import argparse
import os
import subprocess
import sys
from collections import Counter

from . import __version__, chem, network
from .defaults import SETTINGS_YAML
from .energies import Energies, collect
from .orca import SUBMITTED, Queue, job_state, parse_output, slurm_problem
from .project import PROJECT_FILE, Project, ProjectError, expression, find_project, species_sort_key
from .report import overview, plural
from .workflow import (ACCEPTED, VALIDATE, job_folders, prepare, refresh, restart_scan, retire, scan_maximum,
                       scan_points, scan_profile_text, species_result, species_setup, step_progress)


def _die(msg, code=1):
    print("ERROR: " + msg, file=sys.stderr)
    sys.exit(code)


def _root(args):
    if args.project:
        if not os.path.isfile(os.path.join(args.project, PROJECT_FILE)):
            _die("no %s in %s" % (PROJECT_FILE, args.project))
        return os.path.abspath(args.project)
    root = find_project(os.getcwd())
    if not root:
        _die("no %s here or in any parent folder (run 'mechflow init' in your project folder, or use -p)" % PROJECT_FILE)
    return root


def _load(args, strict=True):
    try:
        p = Project.load(_root(args))
    except ProjectError as exc:
        _die(str(exc))
    if p.errors:
        if strict:
            for e in p.errors:
                print("ERROR: " + e, file=sys.stderr)
            _die("%d problem(s) in the project: fix them first (mechflow validate lists them)" % len(p.errors))
        print("WARNING: the project has %d problem(s); run 'mechflow validate'\n" % len(p.errors), file=sys.stderr)
    return p


def _select(p, args):
    return p.select(mech=getattr(args, "mech", None), steps=getattr(args, "step", None),
                    species=getattr(args, "species", None), product=getattr(args, "product", None))


def _rel(p, path):
    return os.path.relpath(path, p.root)


# ============================================================================== init / validate

def cmd_init(args):
    root = os.path.abspath(args.project or os.getcwd())
    target = os.path.join(root, PROJECT_FILE)
    if os.path.exists(target) and not args.force:
        _die("%s already exists (use --force to overwrite it)" % target)
    sdir = os.path.join(root, "species")
    os.makedirs(sdir, exist_ok=True)
    files = sorted((f[:-4] for f in os.listdir(sdir) if f.lower().endswith(".xyz")), key=species_sort_key)
    L = ["# mechflow project file (https://geoffreyweal.github.io/ORCA_Mechanism_Workflow/)", "#",
         "# Species are the xyz files in species/: the file name is the species ID. Recommended names:",
         "#   R1, R2 ... reactants   I1, I2 ... intermediates   P1, P2 ... products",
         "#   C1, C2 ... catalysts and other species used during the reaction",
         "#   S1, S2 ... steps (step Sn goes through transition state TSn)",
         "# Charge and multiplicity: on the xyz comment line as 'charge=1 mult=1', or here (default 0 and 1).", "",
         "project:",
         "  name: %s              # used in Slurm job names" % os.path.basename(root).replace(" ", "_"),
         '  title: ""',
         "  atom_numbering: 0      # 0: atoms counted from 0, as in ORCA inputs; 1: from 1 (Avogadro, GaussView)",
         "", SETTINGS_YAML.rstrip("\n"), "",
         "# -------------------------------------------------------------------------------------",
         "# Species (names are only labels)",
         "# -------------------------------------------------------------------------------------",
         "species:" + ("" if files else " {}")]
    bad = []
    for sid in files:
        try:
            mol = chem.read_xyz(os.path.join(sdir, sid + ".xyz"))
            q, m = chem.comment_charge_mult(mol.comment)
            info = "%s, charge %s, mult %s" % (chem.formula_string(mol.formula), q if q is not None else "0 (default)",
                                               m if m is not None else "1 (default)")
        except (chem.XYZError, OSError) as exc:
            info = "cannot read: %s" % exc
            bad.append(sid)
        L.append('  %s: ""%s# %s' % (sid, " " * max(1, 10 - len(sid)), info))
    L += ["", "# -------------------------------------------------------------------------------------",
          "# Steps: one TS search each (SCAN -> OptTS -> IRC -> Forwards/Backwards).",
          "#   reaction: the species on each side",
          "#   scan:     from  = a species holding every atom of the step (its optimised structure starts the SCAN)",
          "#             atoms = the two atoms whose distance is scanned (numbered as in atom_numbering)",
          "#             to    = final distance in Angstrom",
          "# Example:",
          "#   S1:",
          "#     name: N adds to the carbonyl carbon",
          "#     reaction: R1 + R2 -> I1",
          "#     scan: {from: I1, atoms: [0, 7], to: 3.2}",
          "# -------------------------------------------------------------------------------------",
          "steps: {}", "",
          "# -------------------------------------------------------------------------------------",
          "# Mechanisms: species and steps from the reactants to a product. Example:",
          "#   M1:",
          "#     name: addition then elimination",
          "#     path: [R1 + R2, S1, I1, I2, S2, P1]",
          "# -------------------------------------------------------------------------------------",
          "mechanisms: {}", "",
          "energies:",
          "  reference: %s      # G = 0 for every energy profile" % (
              " + ".join(f for f in files if f.startswith("R")) or "R1 + R2"),
          "  balance_with: []    # species added to balance states, e.g. [P2, C1, C2] (water, acid, base)",
          "  temperature: 298.15",
          "  units: kJ/mol", ""]
    with open(target, "w") as fh:
        fh.write("\n".join(L))
    print("Wrote %s with %d species%s." % (target, len(files), "" if files else " (put your xyz files in species/)"))
    if bad:
        print("Could not read: %s" % ", ".join(bad))
    print("Next: add the steps and mechanisms, then run 'mechflow validate'.")


def cmd_validate(args):
    p = _load(args, strict=False)
    for e in p.errors:
        print("ERROR:   " + e)
    for w in p.warnings:
        print("WARNING: " + w)
    if args.verbose:
        for n in p.notes:
            print("note:    " + n)
    unused_sp = [s for s, sp in p.species.items() if not sp.used_by]
    unused_st = [s for s, st in p.steps.items() if not st.used_by]
    if unused_sp:
        print("note:    species not used by any step or mechanism: %s" % ", ".join(unused_sp))
    if unused_st:
        print("note:    steps not used by any mechanism: %s" % ", ".join(unused_st))
    text = overview(p)
    with open(p.path("overview.md"), "w") as fh:
        fh.write(text)
    n_states = sum(len(m.states) for m in p.mechanisms.values())
    own = [m.id for m in p.mechanisms.values() if m.reference]
    print("\n%s, %s, %s (%s, balanced against %s%s)." % (
        plural(len(p.species), "species", "species"), plural(len(p.steps), "step"),
        plural(len(p.mechanisms), "mechanism"), plural(n_states, "state"), expression(p.reference),
        "; %s use their own reference" % ", ".join(own) if own else ""))
    print("Wrote overview.md")
    if p.errors:
        print("\n%d problem(s) to fix." % len(p.errors))
        sys.exit(1)
    print("No problems found." if not p.warnings else "No errors (%d warning(s))." % len(p.warnings))


# ============================================================================== progress

def cmd_status(args):
    p = _load(args, strict=False)
    species, steps, _ = _select(p, args)
    q = Queue(enabled=not args.no_squeue)
    counts = Counter()
    if not args.steps_only and species:
        print("SPECIES (Step 1 optimisation)")
        for sp in species:
            r = species_result(p, sp, q)
            label = "ok" if r["ok"] else ("CHECK" if r["state"] == "finished" else r["state"])
            counts["species " + label] += 1
            if args.brief and (r["ok"] or label == "absent"):
                continue
            print("  %-7s %-10s %s" % (sp.id, label, r["note"] or sp.name))
    if steps:
        print("\nSTEPS (SCAN -> OptTS -> IRC -> Forwards/Backwards)")
        for st in steps:
            stage, detail, action = step_progress(p, st, q)
            counts["steps " + stage] += 1
            if args.brief and stage in ("DONE", "waiting"):
                continue
            print("  %-7s %-8s %s%s" % (st.id, stage, detail, "   -> mechflow next" if action else ""))
    sp_c = ", ".join("%s %d" % (k.split(" ", 1)[1], v) for k, v in sorted(counts.items()) if k.startswith("species"))
    st_c = ", ".join("%s %d" % (k.split(" ", 1)[1], v) for k, v in sorted(counts.items()) if k.startswith("steps"))
    print("\nSpecies: %s\nSteps:   %s" % (sp_c or "-", st_c or "-"))
    if counts.get("species absent"):
        print("(%d species have no Step 1 input yet: run 'mechflow next')" % counts["species absent"])


def cmd_next(args):
    p = _load(args)
    q = Queue(enabled=not args.no_squeue)
    if args.target:
        try:
            st = p.resolve_step(args.target)
        except ProjectError as exc:
            _die(str(exc))
        if args.image is not None:
            t3 = os.path.join(p.step_dir(st.id), "Step3_Opt_TS")
            if job_state(t3, q) != "absent":
                if not args.force:
                    _die("%s already has a Step3_Opt_TS; add --force to replace it (the old one is kept as "
                         "Step3_Opt_TS_old1, ...)" % st.id)
                print("moved the old OptTS to %s" % _rel(p, retire(t3)))
            try:
                print("%s: %s" % (st.id, prepare(p, st, "prepare_optts", image=args.image)))
            except ValueError as exc:
                _die(str(exc))
            return
        stage, detail, action = step_progress(p, st, q)
        if not action:
            print("%s: %s - %s (nothing to prepare)" % (st.id, stage, detail))
            return
        print("%s: %s" % (st.id, prepare(p, st, action)))
        return
    species, steps, _ = _select(p, args)
    if args.refresh:
        n = refresh(p, job_folders(p, species, steps), q)
        print("rewrote %d prepared job(s) from the current settings" % n)
    n_sp = n_st = 0
    for sp in species:
        msg, warn = species_setup(p, sp, q)
        if msg:
            print("%-7s %s%s" % (sp.id, "WARNING: " if warn else "", msg))
            n_sp += 0 if warn else 1
    for st in steps:
        stage, detail, action = step_progress(p, st, q)
        if action:
            try:
                print("%-7s %s" % (st.id, prepare(p, st, action)))
                n_st += 1
            except (ValueError, OSError) as exc:
                print("%-7s could not prepare %s: %s" % (st.id, stage, exc))
    print("\nPrepared %d species and %d step stage(s). Next: mechflow submit" % (n_sp, n_st)
          if (n_sp or n_st) else "\nNothing new to prepare (see 'mechflow status').")


def cmd_submit(args):
    p = _load(args)
    species, steps, _ = _select(p, args)
    q = Queue(enabled=not args.no_squeue)
    ready = [f for f in job_folders(p, species, steps) if job_state(f, q) == "ready"]
    if args.limit:
        ready = ready[:args.limit]
    if not ready:
        print("No prepared jobs waiting to be submitted (run 'mechflow next' first?).")
        return
    if args.run:
        for f in ready:
            res = subprocess.run(["sbatch", "--parsable", "submit.sl"], cwd=f, capture_output=True, text=True)
            if res.returncode == 0:
                jid = res.stdout.strip().split(";")[0]
                with open(os.path.join(f, SUBMITTED), "w") as fh:
                    fh.write(jid + "\n")
                print("submitted %-10s %s" % (jid, _rel(p, f)))
            else:
                print("FAILED to submit %s: %s" % (_rel(p, f), res.stderr.strip()))
        return
    if args.local:
        path = p.path("run_ready.sh")
        with open(path, "w") as fh:
            fh.write("#!/bin/bash\n# %d prepared job(s), run one after the other on this computer (no Slurm).\n"
                     "# Written by 'mechflow submit --local'. Run it from %s\nset -u\norca_exe=$(which orca)\n"
                     % (len(ready), p.root))
            for f in ready:
                rel = _rel(p, f)
                fh.write('( cd "%s" && echo local > %s && "$orca_exe" orca.inp > output.out ; echo "finished %s" )\n'
                         % (rel, SUBMITTED, rel))
        os.chmod(path, 0o755)
        print("Wrote %s with %d job(s). Run it with:  bash run_ready.sh" % (_rel(p, path), len(ready)))
        return
    path = p.path("submit_ready.sh")
    with open(path, "w") as fh:
        fh.write("#!/bin/bash\n# %d prepared job(s). Written by 'mechflow submit'. Run it from %s\nset -u\n"
                 % (len(ready), p.root))
        for f in ready:
            rel = _rel(p, f)
            fh.write('( cd "%s" && jid=$(sbatch --parsable submit.sl) && echo "$jid" > %s && echo "submitted $jid %s" )\n'
                     % (rel, SUBMITTED, rel))
    os.chmod(path, 0o755)
    print("Wrote %s with %d job(s). Submit them with:  bash submit_ready.sh" % (_rel(p, path), len(ready)))


def cmd_check(args):
    d = os.path.abspath(args.folder)
    print("job state: %s" % job_state(d, Queue(enabled=not args.no_squeue)))
    out = os.path.join(d, "output.out")
    if not os.path.isfile(out):
        return
    r = parse_output(out)
    print("terminated normally:    %s" % r["terminated_normally"])
    print("optimisation converged: %s" % r["opt_converged"])
    if r["irc_converged"]:
        print("IRC directions converged: %d of 2" % r["irc_converged"])
    if r["freqs"]:
        print("imaginary frequencies:  %d %s" % (len(r["imag"]), ", ".join("%.2f" % f for f in r["imag"])))
    if r["gibbs"] is not None:
        print("Final Gibbs free energy:   %.8f Eh" % r["gibbs"])
    if r["e_el"] is not None:
        print("FINAL SINGLE POINT ENERGY: %.8f Eh" % r["e_el"])
    prob = slurm_problem(d)
    if prob:
        print("Slurm reports: %s" % prob)


def cmd_scan_profile(args):
    p = _load(args, strict=False)
    try:
        st = p.resolve_step(args.step_ref)
    except ProjectError as exc:
        _die(str(exc))
    pts = scan_points(p.step_dir(st.id), *st.atoms)
    mol = p.species[st.start].mol
    print("%s (%s): SCAN of %s-%s, %s" % (st.id, st.ts, p.atom_label(mol, st.atoms[0]), p.atom_label(mol, st.atoms[1]),
                                        st.name))
    print(scan_profile_text(pts, scan_maximum(pts)))


def cmd_restart_scan(args):
    p = _load(args)
    try:
        st = p.resolve_step(args.step_ref)
        print("%s: %s" % (st.id, restart_scan(p, st)))
    except (ProjectError, ValueError) as exc:
        _die(str(exc))


def cmd_accept(args):
    p = _load(args, strict=False)
    for ident in args.ids:
        if ident in p.species:
            d = p.species_job_dir(ident)
            what = "the Step 1 structure of %s (bond changes accepted)" % ident
        elif ident in p.steps:
            d = os.path.join(p.step_dir(ident), VALIDATE)
            what = "the Step 4B end points of %s" % ident
        else:
            _die("'%s' is not a species or step ID" % ident)
        if not os.path.isdir(d):
            _die("%s has nothing to accept yet (%s does not exist)" % (ident, _rel(p, d)))
        with open(os.path.join(d, ACCEPTED), "w") as fh:
            fh.write("accepted with 'mechflow accept'\n")
        print("accepted %s" % what)


def cmd_collect(args):
    p = _load(args, strict=False)
    _, _, mechs = _select(p, args)
    text = collect(p, mechs, queue=Queue(enabled=not args.no_squeue))
    print(text)
    print("Written to %s" % _rel(p, p.path("results")))


def cmd_map(args):
    p = _load(args, strict=False)
    _, _, mechs = _select(p, args)
    if not mechs:
        _die("no mechanisms selected")
    q = Queue(enabled=not args.no_squeue)
    status = None
    if args.status:
        sp_status = {sid: network.species_status(species_result(p, sp, q)) for sid, sp in p.species.items()}
        st_status = {sid: network.step_status(*step_progress(p, st, q)) for sid, st in p.steps.items()}
        status = (sp_status, st_status)
    energies = Energies(p, q) if args.energies else None
    outdir = p.path(args.out)
    os.makedirs(outdir, exist_ok=True)
    groups = [(args.name, mechs, "")]
    if args.by_product:
        prods = sorted(set(m.product for m in mechs))
        groups = [("%s_%s" % (args.name, pr.replace(" ", "").replace("+", "_")), [m for m in mechs if m.product == pr],
                   " giving %s" % pr) for pr in prods]
    formats = [f.strip() for f in args.format.split(",") if f.strip()]
    for name, ms, extra in groups:
        g = network.build(p, ms)
        title = "%s: %s%s" % (p.settings["project"].get("title") or p.name, plural(len(ms), "mechanism"), extra)
        dot = os.path.join(outdir, name + ".dot")
        with open(dot, "w") as fh:
            fh.write(network.to_dot(p, g, title, status, energies, args.colour_mechanisms, not args.no_names,
                                    args.direction))
        mmd = network.to_mermaid(p, g, status, not args.no_names)
        with open(os.path.join(outdir, name + ".mmd"), "w") as fh:
            fh.write(mmd)
        with open(os.path.join(outdir, name + ".md"), "w") as fh:
            fh.write("# %s\n\n```mermaid\n%s```\n" % (title, mmd))
        done, problem = network.render(dot, [f for f in formats if f not in ("dot", "mermaid", "md")])
        print("%s: %d nodes, %d connections -> %s" % (name, len(g.nodes), len(g.edges), ", ".join(
            _rel(p, x) for x in [dot, os.path.join(outdir, name + ".md")] + done)))
        if problem:
            print("  " + problem)


def cmd_atoms(args):
    p = _load(args, strict=False)
    if args.species not in p.species or p.species[args.species].mol is None:
        _die("unknown species '%s'" % args.species)
    sp = p.species[args.species]
    mol = sp.mol
    if args.optimised:
        r = species_result(p, sp)
        if not r["xyz"]:
            _die("%s has no optimised structure yet" % sp.id)
        mol = chem.read_xyz(r["xyz"])
    b = chem.bonds(mol, p.tolerance)
    adj = chem.adjacency(len(mol), b)
    if args.pair:
        i, j = (a - p.index_base for a in args.pair)
        if not (0 <= i < len(mol) and 0 <= j < len(mol)):
            _die("atoms are numbered %d to %d" % (p.index_base, len(mol) - 1 + p.index_base))
        d = chem.distance(mol.coords, i, j)
        bonded = chem.is_bonded(mol, i, j, p.tolerance)
        print("%s-%s in %s: %.3f A, %s" % (p.atom_label(mol, i), p.atom_label(mol, j), sp.id, d,
                                          "bonded" if bonded else "not bonded"))
        if bonded:
            print("A SCAN would stretch this bond, e.g.  scan: {from: %s, atoms: [%d, %d], to: %.1f}" % (
                sp.id, args.pair[0], args.pair[1], max(float(p.settings["scan"].get("stretch_to", 3.2)), d + 1.2)))
        else:
            to = 1.05 * (chem.COVALENT_RADIUS[mol.symbols[i]] + chem.COVALENT_RADIUS[mol.symbols[j]])
            print("A SCAN would bring them together, e.g.  scan: {from: %s, atoms: [%d, %d], to: %.2f}" % (
                sp.id, args.pair[0], args.pair[1], to))
        return
    print("%s %s  %s  charge %d  mult %d  (atoms numbered from %d)" % (sp.id, sp.name, sp.formula_str, sp.charge,
                                                                     sp.mult, p.index_base))
    for k in range(len(mol)):
        nb = ", ".join("%s (%.2f)" % (p.atom_label(mol, x), chem.distance(mol.coords, k, x)) for x in sorted(adj[k]))
        print("  %4d  %-2s  bonded to: %s" % (k + p.index_base, mol.symbols[k], nb or "-"))


# ============================================================================== main

def build_parser():
    ap = argparse.ArgumentParser(prog="mechflow", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version="mechflow " + __version__)
    ap.add_argument("-p", "--project", help="project folder (default: the folder holding project.yaml, here or above)")
    sub = ap.add_subparsers(dest="cmd", metavar="command")

    def filters(q):
        q.add_argument("--mech", help="mechanisms, e.g. M1A,M6B or a family M7 (comma separated)")
        q.add_argument("--product", help="mechanisms giving this product, e.g. P1a")
        q.add_argument("--step", help="steps, e.g. S5,S7a or a range S1-S10 (S5 includes S5a and S5b)")
        q.add_argument("--species", help="species, e.g. I4a,I4b")

    def squeue_flag(q):
        q.add_argument("--no-squeue", action="store_true", help="do not ask Slurm which jobs are in the queue")

    q = sub.add_parser("init", help="write project.yaml for the xyz files in species/")
    q.add_argument("--force", action="store_true", help="overwrite an existing project.yaml")
    q.set_defaults(func=cmd_init)

    q = sub.add_parser("validate", help="check the project and write overview.md")
    q.add_argument("-v", "--verbose", action="store_true", help="also show notes")
    q.set_defaults(func=cmd_validate)

    q = sub.add_parser("status", help="where every species and step is up to")
    filters(q)
    squeue_flag(q)
    q.add_argument("--brief", action="store_true", help="hide finished items and species not started")
    q.add_argument("--steps-only", action="store_true", help="only show the steps")
    q.set_defaults(func=cmd_status)

    q = sub.add_parser("next", help="write every input that can be written now")
    filters(q)
    squeue_flag(q)
    q.add_argument("target", nargs="?", help="only this step (ID or folder)")
    q.add_argument("--image", type=int, help="with a step: use this SCAN image as the TS guess")
    q.add_argument("--force", action="store_true", help="with --image: replace an existing Step3_Opt_TS")
    q.add_argument("--refresh", action="store_true",
                   help="first rewrite prepared (not yet submitted) inputs from the current settings")
    q.set_defaults(func=cmd_next)

    q = sub.add_parser("submit", help="submit prepared jobs (writes submit_ready.sh unless --run)")
    filters(q)
    squeue_flag(q)
    q.add_argument("--run", action="store_true", help="call sbatch now")
    q.add_argument("--local", action="store_true", help="write run_ready.sh to run the jobs here without Slurm")
    q.add_argument("--limit", type=int, help="at most this many jobs")
    q.set_defaults(func=cmd_submit)

    q = sub.add_parser("check", help="the procedure's checks on one ORCA job folder")
    q.add_argument("folder")
    squeue_flag(q)
    q.set_defaults(func=cmd_check)

    q = sub.add_parser("scan-profile", help="energy profile of a step's SCAN")
    q.add_argument("step_ref")
    q.set_defaults(func=cmd_scan_profile)

    q = sub.add_parser("restart-scan", help="restart a SCAN from its last converged point")
    q.add_argument("step_ref")
    q.set_defaults(func=cmd_restart_scan)

    q = sub.add_parser("accept", help="accept a Step 1 structure or Step 4B result that was flagged CHECK")
    q.add_argument("ids", nargs="+", help="species or step IDs")
    q.set_defaults(func=cmd_accept)

    q = sub.add_parser("collect", help="Step 5: energies, profiles and product ratios -> results/")
    filters(q)
    squeue_flag(q)
    q.set_defaults(func=cmd_collect)

    q = sub.add_parser("map", help="network map of the mechanisms -> maps/")
    filters(q)
    squeue_flag(q)
    q.add_argument("--status", action="store_true", help="colour by calculation progress")
    q.add_argument("--energies", action="store_true", help="add Gibbs energies where they are known")
    q.add_argument("--by-product", action="store_true", help="one map per product")
    q.add_argument("--colour-mechanisms", action="store_true", help="colour steps used by a single mechanism")
    q.add_argument("--no-names", action="store_true", help="IDs only (no species names)")
    q.add_argument("--format", default="svg,pdf", help="Graphviz formats to render (default svg,pdf)")
    q.add_argument("--direction", default="LR", choices=["LR", "TB"], help="left to right (default) or top to bottom")
    q.add_argument("--name", default="network", help="file name (default network)")
    q.add_argument("--out", default="maps", help="folder (default maps)")
    q.set_defaults(func=cmd_map)

    q = sub.add_parser("atoms", help="numbered atoms of a species, or the distance between two of them")
    q.add_argument("species")
    q.add_argument("pair", nargs="*", type=int, help="two atom numbers")
    q.add_argument("--optimised", action="store_true", help="use the optimised Step 1 structure")
    q.set_defaults(func=cmd_atoms)
    return ap


def main(argv=None):
    ap = build_parser()
    args = ap.parse_args(argv)
    if not getattr(args, "func", None):
        ap.print_help()
        return 0
    if getattr(args, "pair", None) and len(args.pair) != 2:
        _die("give two atom numbers (or none)")
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
