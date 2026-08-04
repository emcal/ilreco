# Generating shower profile tables

How to produce the `.dat` [shower profile](./profiles) for a calorimeter
the shipped tables do not describe. For what the table *is* and how the
reconstruction uses it, see the algorithm page.

## When you need your own

The table's offsets are in cell units, so it encodes the ratio of cell size
to Molière radius, the material, and the longitudinal averaging of the
setup it was made for. Regenerate when any of these changes:

- different **cell size / pitch** relative to the Molière radius (the
  shipped `prof_pwo.dat` is for 2.05 cm PbWO₄ cells, `prof_lg.dat` for
  lead glass) — even the same material with a different crystal width
  needs a new table;
- different **material** (X₀, R_M);
- significantly different **incidence geometry** (large angles, very
  different target distance, thick upstream material);
- a response chain that reshapes the *measured* lateral profile — the
  table must describe reconstructed cell energies, not bare deposits.

## File format

Text, one line per node:

```
i  j  amean  ad2c
```

- `i, j` = 100 × |dx|, 100 × |dy| — offsets in units of 0.01 cell, from 0
  to 500 (up to 5.00 cells); only the triangle `j <= i` is stored, the
  symmetric part is mirrored on load. 125 751 lines total.
- `amean` = mean energy fraction in the cell at that offset.
- `ad2c` = variance of that fraction.

## Recipe

Monte-Carlo (data-driven works the same way if you can select clean single
showers with a known impact point):

1. Simulate single electrons/photons over your energy range, impact points
   uniform across one cell, normal incidence (or your real geometry);
   apply your full response chain to get per-cell *measured* energies.
2. For every event and every cell within 5 cells of the true impact,
   compute the offset `(|dx|, |dy|)` in cell units and the cell's fraction
   of the event's total measured energy — **including zero-signal cells**
   (omitting them biases the tails).
3. Accumulate mean and variance of the fraction in 0.01-cell offset bins;
   symmetrize, smooth thin bins.
4. Write the `j <= i` triangle in the format above.

## Validation

This recipe was validated end-to-end: a table generated from an
independent Geant4 simulation reproduces the shipped `prof_pwo.dat` within
a few percent in the core, gives identical energy resolution, and — because
it is self-consistent with that simulation's response — even improves
position resolution by 15–20 % on the simulated data.

Check a fresh table before trusting it: single showers should give χ²/ndof
of order one, the reconstructed energy scale should be flat versus impact
position within a cell, and the position S-curve should be small and
symmetric.
