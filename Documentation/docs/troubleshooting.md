# Troubleshooting

## The project

**`mechflow validate` says a state cannot be balanced.**
The species written in that state plus the reference's unreacted species and the `balance_with` species
cannot make the reference. Usually a species is missing from the path (e.g. a second molecule that has
already reacted), or `balance_with` lacks a species (water, the acid or the base). The message gives the
missing atoms and charge.

**"… electrons cannot have multiplicity …: check the charge".**
The charge on the xyz comment line (or in `project.yaml`) is probably wrong, e.g. a protonated species
written with `charge=0`.

**"species X and Y have the same atoms and bonds".**
Two files describe the same molecule: conformers, stereoisomers (mechflow does not compare
stereochemistry), or a duplicate.

**"scan.from … does not hold every atom of the step".**
A SCAN starts from one structure, so `from` must contain all the atoms of the step. Use the side of the
reaction that is one molecule (for an addition, the adduct), or make a species of the two molecules
together.

## SCAN

**No energy maximum.**
`mechflow scan-profile S5` shows the profile. Either the step is barrierless in this model (common for
ions combining), or the SCAN did not go far enough (increase `scan.to`), or another distance describes
the reaction better. After changing the step, `mechflow reset S5 --stage SCAN --redo-done --from-scratch`
and `mechflow next` start the SCAN again.

**The SCAN stopped before the end.**
`mechflow reset S5` (the procedure's restart, from the last converged point). If it stopped after passing a
maximum, `next` already offers OptTS from the points it has.

**The automatic TS guess is not the one I want.**
`mechflow next S5 --image 57 --force` (the image numbers are in `scan-profile`).

## OptTS, IRC and Step 4B

**More than one imaginary frequency.**
Start OptTS again from a different SCAN point (`--image`), or from the TS displaced along the second
imaginary mode.

**CHECK: "the Forwards end is neither … nor …".**
The IRC led somewhere else. Look at the end structures (`viewORCA opt` in `Step4B_Validate/Forwards`).
Often the TS is real but belongs to a different step (a proton moved as well). Add that step or species,
or SCAN another distance. If the result is what you wanted after all, run `mechflow accept S5`.

**CHECK: "… is bonded at BOTH ends".**
The TS is not for the bond you scanned. The scanned bond is the same at both ends, so the IRC connects
two other structures.

## Settings

**I changed project.yaml but the inputs did not change.**
Inputs that already exist are never rewritten by `next`. Use `mechflow next --refresh` for prepared jobs that
have not been submitted. A job reset with `mechflow reset` is written again from the current settings (and
remember that all calculations must use the same method, basis set and solvent).

**A comment in my settings breaks the ORCA input.**
ORCA ends a comment at a second `#` on the same line. Keep to one `#` per line; `mechflow validate`
warns about this.
