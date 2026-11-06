# Steps and SCANs

A **step** is one elementary reaction with one transition state. mechflow takes every step through the
procedure: SCAN → OptTS → IRC → Forwards/Backwards. Step **Sn** goes through transition state **TSn**.

```yaml
steps:
  S1a:
    name: methanol opens the protonated epoxide at the CHMe carbon
    reaction: I1 + R2 -> I2a
    scan: {from: I2a, atoms: [0, 3], to: 3.20}
```

| key | meaning |
|---|---|
| `name` | what happens (free text) |
| `reaction` | the species on each side, joined by `+`, e.g. `I1 + R2 -> I2a`. Both sides must have the same atoms and charge |
| `scan.from` | the species whose **optimised** structure starts the SCAN. It must contain every atom of the step |
| `scan.atoms` | the two atoms whose distance is scanned, numbered as in `atom_numbering` |
| `scan.to` | the distance (Å) the SCAN goes to |
| `scan.step` | distance between SCAN points (default `scan.step` in the settings, 0.02 Å) |
| `ts` | the transition state's name, if not TS + the step number (S1a → TS1a) |
| `resources` | resources for this step's jobs, e.g. `{time: "5-00:00", nprocs: 64}` |

## Which structure to start from

The SCAN needs one structure that holds **all the atoms of the step**. That is why `scan.from` is usually
the side of the reaction that is a single molecule:

- **Two molecules join** (`I1 + R2 -> I2a`): start from the adduct I2a and *stretch* the new bond
  (3.2 Å by default). This scans the reaction backwards, which needs no extra start structure and avoids
  placing two molecules by hand (the procedure's Step 2A).
- **A molecule falls apart** (`I3 -> I4 + P2`): start from I3 and stretch the bond that breaks.
- **A shift or a ring closure** (one molecule on each side): start from either side. If the two atoms are
  not bonded yet, the SCAN brings them together to about a bond length.

mechflow looks at the start structure to decide which: if the two atoms are bonded there, the SCAN
stretches the bond (and `to` must be longer than the bond); otherwise it brings them together (and `to`
must be shorter).

## Choosing the two atoms

Pick the bond that forms or breaks in the step. When several bonds change (a concerted step), pick the one
that changes most. `mechflow atoms` helps:

```console
$ mechflow atoms I2a 0 3
C0-O3 in I2a: 1.427 A, bonded
A SCAN would stretch this bond, e.g.  scan: {from: I2a, atoms: [0, 3], to: 3.2}
```

## What the SCAN input looks like

mechflow measures the distance in the optimised start structure and keeps the procedure's 0.02 Å
spacing, finishing at the point closest to `to`. For S1a (1.427 Å in the xyz file; the optimised value will
differ):

```
%GEOM
    MaxIter 200   # optimisation cycles allowed per SCAN point
    SCAN B 0 3 = 1.427, 3.207, 90 END
END
```

`MaxIter 200` is the procedure's fix for SCAN points that hit ORCA's 50-cycle limit. Set
`orca.scan_maxiter: null` to leave it out.

## Proton transfers

Proton transfers between a species and the acid or base are usually fast equilibria. Treat them as
connections between species in a [mechanism](mechanisms.md) (two species written next to each other)
rather than as steps. Only bond-making and bond-breaking steps at heavy atoms need a TS search. If a
proton transfer matters (an intramolecular shift, say), make it a step like any other.

## Isomeric series of steps

When the same step happens on two isomeric branches, give the two versions the same number with a letter,
e.g. `S1a` and `S1b`. Filters understand this: `--step S1` selects both, `--step S1-S9` a range.
