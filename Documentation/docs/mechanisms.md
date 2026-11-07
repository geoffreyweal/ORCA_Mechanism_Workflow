# Mechanisms

A **mechanism** is a path from the reactants to a product, written as the species and steps it goes
through. In the acid-catalysed opening of propylene oxide (R1) by methanol (R2), methanol can attack either
carbon of the protonated epoxide (I1), which gives two mechanisms:

```yaml
mechanisms:
  M1A:
    name: attack at the more substituted carbon
    path: [R1 + R2, I1, S1a, I2a, P1a]
  M1B:
    name: attack at the less substituted carbon
    path: [R1 + R2, I1, S1b, I2b, P1b]
```

| key | meaning |
|---|---|
| `path` | species (or `+`-joined groups of species) and step IDs, in order |
| `name` | free text |
| `product` | the product label used for comparisons (default: the last species of the path) |
| `family` | mechanisms of the same family are compared with each other (default: the ID without a final letter, so M1A and M1B are family M1) |
| `reference` | the mechanism's own reference state, for a route that needs more molecules than `energies.reference` (see [below](#a-mechanisms-own-reference)) |

## Reading a path

- A **step** (`S1a`) means: the species on one side of the step, its transition state TS1a, then the species
  on the other side. You do not have to write the step's own species around it, but you may.
- **Two species next to each other** (`R1 + R2, I1` or `I2a, P1a`) are joined by an equilibrium, usually a
  proton transfer. There is no transition state between them.
- mechflow works out which way round a step is used: if the path arrives at the step's right-hand side, the
  step is used backwards.

The path of M1A therefore reads:

R1 + R2 ⇌ I1 —S1a→ I2a ⇌ P1a

The epoxide is protonated, methanol opens it at the CHMe carbon (step S1a, through TS1a), and the oxonium ion
I2a gives its proton back to the base.

## States and balancing

For energies, every point of a path (a **state**) must contain the same atoms as the starting point, the
**reference**. You only write the species that matter; mechflow adds the rest:

```yaml
energies:
  reference: R1 + R2 + C1          # G = 0: both reactants and one acid molecule
  balance_with: [C1, C2]           # HA and A- may be added to any state
```

A number in front of a species counts it more than once, e.g. `R1 + R2 + 2 C1`.

To each state mechflow adds:

1. reactants of the reference that have not reacted yet (R2 in the states before methanol adds), and
2. the `balance_with` species (here the acid HA and its base A⁻; water, for example, in a condensation), as
   separate molecules.

It solves for whole numbers of these, using as few molecules as possible. In the example, the state I1
(the protonated epoxide) becomes `I1 + R2 + C2`: the acid C1 has given its proton to the epoxide, so HA
became A⁻, and methanol R2 has not reacted yet. TS1a and I2a become `TS1a + C2` and `I2a + C2`, and the
product state is `P1a + C1`, with the proton back on the acid.

`mechflow validate` reports any state that cannot be balanced. A state that needs atoms the reference does
not have, or that could be balanced in two ways, must be written out in full.

## Spectators

When a state holds two molecules and the next step uses only one of them, the other molecule is carried
through the step unchanged. It is a **spectator**. In a path such as

```yaml
path: [R1 + R2, ..., I7a + I9a, S5a, I7a + I11a, ...]     # step S5a is I9a -> I11a
```

I7a is a spectator of S5a. The three states are `I7a + I9a`, `TS5a + I7a` and `I7a + I11a`, and the energy
of the TS state includes G(I7a). You do not have to write anything for this: mechflow compares the step with
the state before it and carries the rest along.

## A mechanism's own reference

Every state of a mechanism must contain the atoms of the reference. A route that needs more molecules than
`energies.reference` has, for example a second molecule of a reactant, needs its own reference. In the
epoxide example, the product P1a can open a second protonated epoxide and give a dimer ether (P2):

```yaml
mechanisms:
  M2:
    name: the product opens a second protonated epoxide
    reference: 2 R1 + R2 + C1          # two epoxides, one methanol, one acid
    path: [R1 + R2, I1, S1a, I2a, P1a + I1, S3, I3, P2]
```

Until the second epoxide reacts, the balancing adds it as an unreacted R1: the state I1 becomes
`I1 + R1 + R2 + C2`. I2a then gives its proton to the second epoxide (`P1a + I1`), and the OH group of P1a
opens it in step S3 (`P1a + I1 -> I3`).

The extra R1 is in every state of M2 and in its reference, so it cancels: the energies of M2 are on the same
scale as those of the other mechanisms. For the [search through the whole
network](energies.md#best-routes-through-the-whole-network), a mechanism's own reference must contain the
project's reference (here `R1 + R2 + C1` and one more R1); a mechanism whose reference does not is left out of
that search.

## Products and families

`mechflow collect` compares mechanisms in two ways:

- **By product.** For each product it finds the mechanism with the lowest highest TS, then estimates
  the product ratio from the best mechanisms (see [Energies](energies.md)).
- **Within a family.** Mechanisms of one family that end in different products, such as M1A (→ P1a) and
  M1B (→ P1b), are compared directly.

It also joins every mechanism and every step into one network and finds the best route to each product
through it, combining steps of different mechanisms (see [Energies](energies.md#best-routes-through-the-whole-network)).

## Checks

`mechflow validate` checks every mechanism:

- every species and step exists,
- every state balances against the reference,
- the steps' two sides have the same atoms and charge.

It then writes `overview.md`, with every path written out as above.
