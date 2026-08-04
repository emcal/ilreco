# Island clustering and peak finding

The first two reconstruction stages turn a list of fired cells into shower
candidates. They are purely combinatorial — no profile tables involved —
and are identical under both reconstruction metrics.

## Islands

Cells that touch — sharing an edge **or a corner** (8-connectivity) — are
grouped into one *island* (raw cluster). An electromagnetic shower in a PbWO₄-type
calorimeter typically fires a 3×3 to 5×5 block of cells; two showers
landing within a few cells of each other merge into a single island.

There is no energy threshold at this stage: every stored cell
participates. Islands are found by ordering cells by their packed address
and merging address-adjacent runs, which makes the stage allocation-free
and fast (it was designed for online use).

## Peaks

Within an island, a cell is a **peak** if none of its 8 neighbors (within
the island) has a higher energy, and its energy exceeds the seed
threshold. Each peak is one shower candidate.

The seed threshold is energy-dependent: starting from the configured value
(default 10 MeV), it is scaled up with the island's total energy as
`seed × 7·ln(1 + E_island)` (the factor is 20 for lead glass). The scaling
suppresses fake peaks from shower fluctuations, which grow with shower
energy: a 100 MeV fluctuation two cells from the core of a 10 GeV shower
would otherwise become a shower candidate.

What happens next depends on the peak count:

- **one peak** — the island goes directly to the
  [single-shower fit](./chi2-method), and optionally to the
  [second-order separation](./two-gamma) ("is this really two merged
  showers?");
- **several peaks** — the island's cells are first divided between the
  peaks by iterative profile-based sharing
  ([overlapping showers](./two-gamma)), and each share is then fitted as
  a single shower.

## Fiducial classification

Each candidate is labeled by where its seed cell sits: interior (`type 0`),
at the border of a beam hole or dead region (`type 1`), or on the outer
boundary ring (`type 2`). With a cell-existence mask configured, the labels
derive from the actual geometry; candidates produced by a two-shower split
carry `+10`. The labels never change energies or positions — they let the
analysis decide which regions to trust.

## Practical limits

Fake peaks from fluctuations are the dominant artifact of this stage: with
the default thresholds, a few-GeV shower has a ~20 % chance that a
secondary fluctuation peak survives and produces a low-energy satellite
cluster. This rate is a property of the peak finder (identical in both
metrics); raise the seed threshold to trade satellite suppression against
efficiency for genuine soft particles.
