# Energies and selectivity (Step 5)

`mechflow collect` reads the Gibbs free energies (`Final Gibbs free energy` in each `output.out`) and
writes the results to `results/`:

| file | contents |
|---|---|
| `species.csv` | every species: Step 1 state, G, electronic energy, notes |
| `steps.csv` | every step: stage, imaginary frequency, G of the TS, forward and reverse barriers, reaction energy, rate constants, and the barrier from the SCAN start (see [Rate constants](#rate-constants)) |
| `profiles.csv` | every state of every mechanism, with the species added to balance it, and G relative to the reference |
| `mechanisms.csv` | for each mechanism: highest TS, rate-limiting TS, effective barrier, whether every step is validated |
| `families.csv` | mechanisms of one family that give different products, compared |
| `network_routes.csv` | the best routes to each product through the whole network, up to 10 per product (see [below](#best-routes-through-the-whole-network)) |
| `summary.md` | the tables above, readable (also printed on screen) |

## Where the numbers come from

- **Species:** G from Step 1, only if it converged with no imaginary frequency.
- **Transition states:** G from Step 3 (OptTS), only if it converged with exactly one imaginary frequency.
  A mechanism counts as validated only when all of its steps are DONE (Step 4B passed).
- **States:** the sum over the species in the state, including the species added to balance it (see
  [Mechanisms](mechanisms.md#states-and-balancing)). A TS state is the TS plus its balancing species.

Every state of every mechanism has the same atoms, so all energies are relative to the same zero: the
reference.

## Barriers

For each mechanism, mechflow reports:

- **highest TS:** the highest transition state along the path, relative to the reference;
- **effective barrier:** the largest climb from a state to a later transition state, i.e. the TS minus the
  lowest state before it. This is the barrier that limits the rate when an intermediate lies below the
  reference;
- **rate-limiting TS:** the TS of the effective barrier.

## Rate constants

`steps.csv` gives the barriers of every step, and their rate constants from the Eyring equation
(transmission coefficient 1):

\[
k = \frac{k_\mathrm{B}T}{h}\exp\left(-\frac{\Delta G^{\ddagger}}{RT}\right)
\]

| columns of `steps.csv` | meaning |
|---|---|
| `dG_act_forward_kJ/mol`, `k_forward` | from the left side of the step (as written in `reaction:`) to its TS |
| `dG_act_reverse_kJ/mol`, `k_reverse` | from the right side to the TS: the step used backwards |
| `scan_start`, `dG_act_from_start_kJ/mol`, `k_from_start_s-1` | from the SCAN start species (`scan.from`) to the TS |

- k is in s⁻¹ when the side a step starts from is one molecule. When that side is two molecules, k is a
  second-order rate constant: in M⁻¹ s⁻¹ with `standard_state: 1 M` (see [Standard state](#standard-state)),
  which is the right choice for reactions in solution.
- The barrier from the SCAN start matters for a step scanned from a complex of two molecules, for example a
  proton transfer `X + HA -> XH+ + A-` scanned from the hydrogen-bonded complex X···H–A. It is the barrier
  inside the complex (first order, s⁻¹), without the cost of bringing the two molecules together, which is in
  the forward barrier.
- At 298 K, k<sub>B</sub>T/h = 6.2 × 10¹² s⁻¹: a barrier of 50 kJ/mol gives about 10⁴ s⁻¹ and one of
  100 kJ/mol about 10⁻⁵ s⁻¹ (a half-life of about 10 hours).
- The barriers are in the units of `energies.units`.

## Product ratios

If the steps that lead to different products are slower than the equilibria before them (Curtin–Hammett
conditions), the product ratio follows from the difference in the free energies of the transition states
that decide between them:

\[
\frac{[\mathrm{P_a}]}{[\mathrm{P_b}]} = \exp\left(-\frac{G^{\ddagger}_\mathrm{a} - G^{\ddagger}_\mathrm{b}}{RT}\right)
\]

The deciding transition state of each mechanism is its highest TS **after the point where the two
mechanisms part**. Steps they share cannot decide between them. If two mechanisms share their first step,
and that step has the highest barrier of all, the ratio comes from the highest TSs after that step.

`collect` gives two kinds of ratio:

- **Between products:** for each product, the mechanism with the lowest highest TS, compared with the best
  mechanism to each other product.
- **Within each family:** e.g. M1A (→ P1a) against M1B (→ P1b). This shows what each mechanism alone would
  predict.

## Best routes through the whole network

The mechanisms in `project.yaml` are the routes you thought of. `collect` also joins every mechanism and
every step into one network and searches it for the most likely route to each product. A route found this
way can combine steps from different mechanisms, use any step in either direction, and use steps that are in
no mechanism at all.

### What is in the network

**States.** Every state of every mechanism, and both sides of every step. A state is named by its *core*
species: the species that are not added automatically (unreacted reference species and the `balance_with`
species). In the epoxide example of [Mechanisms](mechanisms.md), `I1` stands for the balanced state
`I1 + R2 + C2`, and `I1 + P1a` for a state with two molecules. The reference has no core species. A state
that appears in several places (in two mechanisms, or as the side of a step) is one state.

**Connections.**

1. Two neighbouring states of a mechanism path: through the TS between them, or by a proton transfer with
   no barrier when there is no TS (unless `proton_transfers: steps`, see below).
2. Every step, between its two sides, either way round, through its TS.
3. Automatic proton transfers (only with `proton_transfers: auto`, the default; see below).

All energies are relative to `energies.reference`. States of a mechanism with its own `reference:` (see
[Mechanisms](mechanisms.md#a-mechanisms-own-reference)) are measured from that reference instead. Its extra
molecules are in the state and in the reference, so they cancel, and every state has the same zero.

### Proton transfers without a barrier

`energies.proton_transfers` says which proton transfers the search may make without a barrier:

| value | proton transfers without a TS |
|---|---|
| `auto` (default) | between two neighbouring states of a mechanism path with no step between them, and the automatic ones (next section) |
| `path` | only between two neighbouring states of a mechanism path with no step between them |
| `steps` | none: every proton transfer has to be a step, with its own TS |

Use `steps` when you write every proton transfer as a step, e.g. `X + C1 -> XH+ + C2` with the acid (C1) and
its conjugate base (C2) as species of the step, scanned from the hydrogen-bonded complex. Otherwise the
connections without a barrier would let a route skip those steps and their barriers. With `steps`, two
neighbouring states of a mechanism path that have no step between them are not joined, and `summary.md` says
where that happens.

`mechflow collect --proton-transfers auto` (or `path`, `steps`) replaces the setting for one run, for example to
see how much the barriers of the proton transfers change the best routes.

### Automatic proton transfers

With `proton_transfers: auto`, two states are joined by a proton transfer when they hold the same molecules
apart from where H atoms sit on N, O, S, P, Se or a halogen. For example, the oxonium ion with A⁻ (`I2a + C2`) and the product with HA
(`P1a + C1`) are joined, whether or not a mechanism writes them next to each other.

mechflow compares the molecules by their connectivity, read from the xyz files in `species/`:

- the bonds between heavy atoms count; H atoms are left out;
- each carbon (and any other atom that is not N, O, S, P, Se or a halogen) also counts its H atoms;
- every molecule of the balanced state counts, including the balancing species. HA and A⁻ look the same,
  so a proton may move between a species and the acid or base.

Limits:

- A change in the number of H atoms on a carbon is not joined automatically: keto–enol tautomers and
  hydride shifts need a step, or two species written next to each other in a mechanism path.
- Stereoisomers and conformers have the same connectivity, so they look the same. Two diastereomers (a cis
  and a trans diol, say) written as separate species are joined as if a proton transfer turned one into the
  other.
- The connectivity is that of your xyz files (bonds from distances, as everywhere in mechflow). A bond that
  is too long or too short in the file gives the wrong connectivity.

`network_routes.csv` and `summary.md` list every automatic proton transfer a route uses, so you can check them.

### The best route

The best route to a product is the one whose **highest point** (state or TS) is lowest. Every state and TS
is relative to the same reference, so the highest point is the top of the climb from the reactants, and the
route with the lowest top is the easiest way to the product. If two routes share their highest point, the
next-highest point decides, and so on. When a route passes an intermediate far below the reference, its
effective barrier (also given) can matter more than its highest point.

- A state or TS without an energy (not computed yet, or not converged) cannot be passed. The best routes can
  change as more calculations finish.
- A product is any state whose core species are all products (role `product`). `P1a` and `2 P1a` count as
  the same product.
- A route may pass through another product, if one product turns into the other.

For each product, `summary.md` gives the highest point, the effective barrier along the route (as for
mechanisms: a TS minus the lowest state before it), the mechanisms whose steps the route uses, and the route
written as a mechanism path, with notes on steps used backwards and automatic proton transfers:

```
- P1a: `[R1 + R2 + C1, I1, S1a, I2a, P1a]`
- P2: `[R1 + R2 + C1, I1, S1a, I2a, I1 + P1a, S3, I3, P2]` (needs reference 2 R1 + R2 + C1)
```

Copy a route into `project.yaml` as a new mechanism to get its full energy profile and map. If the route says
`needs reference ...`, give the mechanism that `reference:`.

**Ratios** between products come from their best routes, in the same way as for mechanisms: the highest TS
of each route after the point where the two routes part.

**Alternatives.** For every TS, mechflow also finds the best route to each product that goes through it
(the best route from the reference to one side of the TS, then the TS, then the best route from the other
side to the product). Routes that would pass the same state twice are left out. From these:

- `network_routes.csv` lists up to 10 different routes per product, best first;
- `summary.md` lists up to five **alternative bottlenecks** per product: the best routes whose highest point
  is another TS (or state), and how far above the best route each one is.

| column of `network_routes.csv` | meaning |
|---|---|
| `product`, `rank` | the product; rank 1 is the best route |
| `bottleneck_TS` | the highest point of the route: a TS, or a state if it is higher than every TS on the route |
| `G_highest_point_kJ/mol` | G of the highest point, relative to the reference (in the units of `energies.units`) |
| `effective_barrier_kJ/mol`, `rate_limiting_TS`, `from_state` | the largest climb from a state to a later TS |
| `mechanisms` | the mechanisms the steps come from, in order, e.g. `M2 (S4, S5) + M1 (S2)` |
| `notes` | steps used backwards, automatic proton transfers, and the reference the route needs |
| `route` | the route as a mechanism path |

The search always uses the whole project, also when `collect` is given filters such as `--mech`. A step with
a second molecule carried through it (a [spectator](mechanisms.md#spectators)) is used only where a mechanism
path writes it.

## Standard state

ORCA's Gibbs energies are for an ideal gas at 1 atm. In solution, a 1 mol/L standard state adds
RT ln(RT/p) = 7.9 kJ/mol (298 K) to each separate molecule. This changes energy differences between states
that have different numbers of separate molecules. Set

```yaml
energies:
  standard_state: 1 M
```

to include it. Comparisons between mechanisms at the same point (e.g. two isomeric transition states)
are not affected.

## Things to keep in mind

- Charged species are treated as separate ions, with implicit solvent. Absolute energies of states with
  separated ions are therefore rough. Differences between isomeric states (A versus B) are much more
  reliable, because the same errors appear on both sides.
- Unless they are steps of their own (with `proton_transfers: steps`), proton transfers are treated as
  equilibria without barriers.
- The energies are only as good as the conformers. For the transition states that decide a selectivity,
  check the conformers.
