# The log-likelihood method

The historic position fit and cluster separation minimize a profile χ²: each
cell's measured energy fraction is compared to the *mean* fraction a shower
at the hypothesized position would deposit, weighted by a parametrized
variance. The log-likelihood (LL) method replaces this metric with the full
tabulated probability distribution of the fraction — instead of asking "how
far is this cell from the mean profile", it asks "how *probable* is this
cell's content for a shower at this position".

```python
import ilreco
calo = ilreco.Calorimeter(
    20, 20, profile="pwo",
    method="log_likelihood",
    ll_profiles="llprof_asset.bin",   # per-energy tables
    split=True)                       # optional two-shower separation
```

The cluster **energy** is unchanged — it stays the containment-corrected
moments estimate of the χ² path (by design; refining the energy by
likelihood maximization is a possible later extension).

## The tabulated distributions

For a particle hitting the calorimeter at a known point, the quantity
tabulated is the fraction `f = e_vis / E_incident` a cell collects, as a
function of the cell's offset from the impact point. Offsets are measured
in 0.01-cell steps; cells are square, so distributions are symmetric under
sign flips and axis swaps, and each offset is stored once as
`(offset_major, offset_minor)` with `minor <= major <= 500` — 125 751 nodes.
Cells below the readout threshold record `f = 0`; that probability is
tracked separately.

![node distributions](./ll-node-distributions.png)

*Measured `f` distributions at six offset nodes, 1–20 GeV. Near the shower
center the distribution is a narrow peak; away from it, a rapidly falling
near-exponential from `f = 0`. Both the shapes and the zero-probability `p0`
depend strongly on energy — hence per-energy tables, with the table nearest
in log-energy chosen per cluster.*

## What the library evaluates

Per node and energy the asset stores the negative log-likelihood
`nll(f) = -ln(p(f) / pmax)` on a shared log-`f` bin grid, where `p` is the
probability density of `f` and `pmax` its maximum — so a cell sitting at
the most probable fraction for the hypothesis costs zero, and improbable
cells cost more. A separate `nll(0)` covers cells that recorded nothing:
a zero cell is a *censored* observation ("somewhere below the readout
threshold"), and its likelihood is the measured zero probability spread
over that interval.

No fitting happens at run time — the tables are read as they are. How they
are produced from simulation is
[generating log-likelihood tables](./generating-ll-tables).

## Position fit

The likelihood of a position hypothesis is the sum of `nll` over the
cluster's cells (`f = e_cell / E_cluster`) and over the zero-signal
neighbor cells (`nll(0)`), with every offset rounded to the nearest table
node. The fit is a coarse-to-fine descent on the 0.01-cell node grid
(steps of 5, 2, 1 nodes); the final position is the likelihood-weighted
average of the best four nodes of the final 3×3 neighborhood. The
`chi2` field of the output carries `nll / ndof` for this method.

## Separating merged showers

With `split=True`, a single-peak cluster is tested against a two-shower
hypothesis (both methods use the same machinery, differing only in the
metric; the full comparison is on [overlapping showers](./two-gamma)):

1. first and second moments of the cluster give the principal axis;
2. scanning the energy asymmetry `a = (E1-E2)/E` seeds the pair — the
   second moment fixes the separation via `lam = E1 E2 / E^2 * d^2`;
3. coordinate descent refines the two positions;
4. the split is accepted if the per-dof metric improves by more than
   `split_threshold` **and** the pair passes the invariant-mass
   physicality cut (`zcal`-dependent, as in the historic code).

For the χ² metric each cell is compared to the sum of both predicted
showers with a combined variance. For the LL metric each cell is scored
once, against the shower predicting the larger share, after subtracting
the other shower's predicted contribution.

## Measured performance

20×20 PbWO₄ (2.05 cm cells), photons with the full response chain,
compared against the [χ² fit](./chi2-method) on identical events:

![single-photon comparison](./ll-single-comparison.png)

- **Single-photon positions**: 9–13 % smaller σ68 at 1–10 GeV, and the
  systematic mid-cell pull is up to 3× smaller. The largest gain is for
  photons near cell edges at low energy; the one exception is 0.5 GeV
  photons at the exact cell center (slightly worse).
- **Single-photon energies**: identical by construction (shared path).
- **Merged clusters**: at the same wrong-split rate the LL decision
  recognizes ~2× more genuinely merged two-photon clusters, and divides
  the cluster energy between the two photons markedly better — see
  [overlapping showers](./two-gamma) for the full study.
- **Speed**: without the split stage the LL fit is *faster* than the χ²
  fit (node-grid lookups vs repeated interpolation); with the split on it
  is currently slower because it tests every single-peak cluster.

## Known limitations (v1)

- Table nodes at offsets ≡ 50 (mod 100) are neighbor interpolations
  (crack-impact artifact of the ideal training geometry).
- Overlapping showers are scored with the simpler per-peak assignment
  above; a joint two-shower likelihood is a recorded future improvement.
- Multi-peak clusters: the per-peak SHARED cell energies produced by the
  iterative profile sharing are deterministic profile fractions, not
  single-shower samples — the tables do not describe them, and a measured
  2x position degradation resulted. The LL method therefore falls back to
  the chi2-metric position fit for shared clusters (single-peak clusters,
  the vast majority, use the full LL fit). A proper joint likelihood for
  overlaps would supersede this.
- The energy is not likelihood-refined (kept from the moments estimate).
- `split_threshold` defaults to 0.5 per dof pending tuning on simulated
  two-photon samples.
- Table selection is nearest-in-log-energy; interpolating `nll` between
  the two bracketing energy tables is a recorded improvement for the
  tuning phase.
- A table with missing nodes (incomplete production) must not be shipped:
  empty nodes evaluate as "no information" (nll = 0), which makes hole
  patterns attractive to the position fit. Leave such an energy out of the
  asset — clusters fall back to the nearest complete table in log-energy —
  until the production is complete.
