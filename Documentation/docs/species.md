# Species: the xyz files

A project starts from a folder called `species/` with **one xyz file per species**. Every reactant,
intermediate, product and catalyst that appears in a step or a mechanism needs one.

```
my_reaction/
    species/
        R1.xyz   R2.xyz        reactants
        C1.xyz   C2.xyz        catalyst (e.g. an acid HA and its base A-)
        I1.xyz   I2.xyz ...    intermediates
        P1.xyz   P2.xyz        products (and by-products such as water)
```

## Names

The file name (without `.xyz`) is the species **ID**, used everywhere else in the project.

| prefix | role | examples |
|---|---|---|
| R | reactant | `R1`, `R2` |
| I | intermediate | `I1`, `I2a`, `I2b` |
| P | product (including by-products) | `P1`, `P1a`, `P2` |
| C | catalyst, or anything used during the reaction and given back | `C1` (HA), `C2` (A⁻) |

The prefix only decides how a species is drawn in the [network maps](maps.md): R is the starting point,
and P nodes get a double border. IDs may use letters, digits, `_`, `-` and `.`, and must start with a letter.
They cannot contain `+` or spaces, because `+` joins species in reactions.

!!! tip "Isomeric series"
    When a reaction can give two isomers, a letter after the number keeps related species together: P1a /
    P1b for the two products, I2a / I2b for the same intermediate on the two branches, and S1a / S1b for the
    same step.

## Charge and multiplicity

Write them on the comment line (line 2) of the xyz file:

```
17
I2a  oxonium ion, attack at CHMe  charge=1 mult=1  |  key atoms: 0 C(H)Me, 1 CH2, 2 O(ring), 3 O(MeOH)
C      0.12345678   ...
```

`charge=` and `mult=` can be anywhere on the line; the rest is free text. If they are missing, the charge
is 0 and the multiplicity 1. They can also be given in `project.yaml` (see below). If both places give
a value, the two must agree.

mechflow checks every species: a charge that leaves an even number of electrons with an even multiplicity
(or odd with odd) is reported as an error. This catches most wrong charges.

## Listing the species in project.yaml

The `species:` section gives each species a name (just a label) and, if needed, other settings. The pages
under *Setting up a reaction* use one small reaction as their example: the acid-catalysed opening of
propylene oxide by methanol. Methanol can attack either carbon of the protonated epoxide, so there are two
products, P1a and P1b.

```yaml
species:
  R1: propylene oxide
  R2: {name: methanol, nprocs: 8}          # fewer cores for a small molecule
  C1: triflic acid (HA)
  C2: {name: triflate (A-), charge: -1}
  I1: protonated epoxide
  I2a: oxonium ion, attack at CHMe
  I2b: oxonium ion, attack at CH2
  P1a: 2-methoxypropan-1-ol
  P1b: 1-methoxypropan-2-ol
```

| key | meaning |
|---|---|
| `name` | a label shown in `overview.md`, `status` and the maps |
| `charge`, `mult` | instead of the xyz comment line |
| `file` | a different xyz file (path relative to the project folder) |
| `role` | `reactant`, `intermediate`, `product` or `catalyst`, if the ID prefix does not say it |
| `nprocs`, `maxcore`, `time`, `mem_per_cpu`, `partition` | resources for this species' Step 1 job only |

A species whose xyz file is in `species/` but which is not listed is still included, with the defaults.
`mechflow init` writes the `species:` section for you from the files in `species/`.

## Atom numbers

Steps refer to atoms by number. By default atoms are counted **from 0, as in ORCA input files**
(`atom_numbering: 0`). If you read numbers off a program that counts from 1 (Avogadro, GaussView,
Chemcraft), set `atom_numbering: 1` in the `project:` section and mechflow converts them for ORCA.

Different species do not need the same atom order: each step's SCAN atoms refer to the one species the
SCAN starts from. Putting the atoms that react first in every file makes the steps easier to read. In the
epoxide example, atom 0 is always the CHMe carbon, atom 1 the CH2 carbon, atom 2 the ring oxygen and
atom 3 the oxygen of methanol.

`mechflow atoms` prints the numbered atoms of a species and what each is bonded to:

```console
$ mechflow atoms R1
R1 hydrogen cyanide  CHN  charge 0  mult 1  (atoms numbered from 0)
     0  C   bonded to: N1 (1.15), H2 (1.07)
     1  N   bonded to: C0 (1.15)
     2  H   bonded to: C0 (1.07)
```

## Good starting structures

- Any reasonable 3D structure will do, for example from Avogadro, a force field, or CREST. Step 1 optimises it.
- **Conformers matter** for transition states. For the steps that decide a selectivity, it is worth using
  the lowest conformer from a conformer search.
- After Step 1, mechflow compares the bonds of the optimised structure with your xyz file. If a bond broke
  or formed (a proton jumped, a ring opened), the species is flagged CHECK. Look at it, then either fix the
  xyz file or run `mechflow accept <ID>` to keep the result.
- If two species have exactly the same atoms and bonds, `mechflow validate` warns you. They are either
  conformers or stereoisomers, or one of them is not needed.
