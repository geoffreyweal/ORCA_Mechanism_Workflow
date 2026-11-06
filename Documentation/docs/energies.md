# Energies and selectivity (Step 5)

`mechflow collect` reads the Gibbs free energies (`Final Gibbs free energy` in each `output.out`) and
writes the results to `results/`:

| file | contents |
|---|---|
| `species.csv` | every species: Step 1 state, G, electronic energy, notes |
| `steps.csv` | every step: stage, imaginary frequency, G of the TS, forward and reverse barriers, reaction energy |
| `profiles.csv` | every state of every mechanism, with the species added to balance it, and G relative to the reference |
| `mechanisms.csv` | for each mechanism: highest TS, rate-limiting TS, effective barrier, whether every step is validated |
| `families.csv` | mechanisms of one family that give different products, compared |
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
- Proton transfers are treated as equilibria without barriers.
- The energies are only as good as the conformers. For the transition states that decide a selectivity,
  check the conformers.
