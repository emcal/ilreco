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

## Edge cases and failure behavior

What the method does when inputs are not the happy path. Everything here
is deliberate; "refused" means an error return, never a crash or exit.

**Loading and configuration**

- A missing, unreadable, truncated, or foreign table file is **refused at
  load** with a message (`ilreco_config_load_ll_profiles` returns −1, the
  Python constructor raises `ValueError`). A corrupt header, an
  unsupported version, and energy blocks out of ascending order are all
  detected. On failure the configuration keeps whatever tables it had
  before — a failed reload cannot leave it half-updated.
- Selecting `method="log_likelihood"` **without loaded tables is refused**
  — you cannot configure a calorimeter that would silently have nothing to
  evaluate. (Defensively, if the reconstruction ever ran in LL mode with
  no tables, it falls through to the χ² path rather than crash.)
- The χ² profile is still **required**: the energy path, the containment
  correction, and the sharing of multi-peak clusters always use it.

**Cluster energy vs the tabulated energies**

- The table **nearest in log-energy** to the cluster's reconstructed
  energy is used; there is no interpolation and no extrapolation. A
  cluster far outside the tabulated range is scored against the nearest
  table and the resolution advantage degrades gracefully toward (or
  below) the χ² fit — measured example: 20 GeV photons scored against a
  10 GeV table resolve to 1.4 mm versus 0.64 mm for χ². Ship tables
  covering your energy range.
- The selection uses the *reconstructed* energy, so a badly mismeasured
  cluster may be scored against a neighboring energy's table — a
  second-order effect, since neighboring tables differ smoothly.
- A **non-positive cluster energy** (pathological input) skips the LL fit
  entirely; the cluster keeps its moments position, as in the χ² path.

**Geometry and range**

- Cells (and silent neighbors) **beyond the tabulated 5-cell offset**
  contribute nothing — "no information", the same for every hypothesis.
- The position fit is **confined to ±1 cell around the moments seed**. A
  genuine single-shower correction is always sub-cell; the confinement
  also removes a pathological attractor — far from the cluster all cells
  would fall outside the tabulated range and a runaway hypothesis would
  look artificially "likely".
- Fractions outside the binned range are clamped to the first/last bin;
  cells with **zero recorded energy** are scored with the dedicated
  censored-zero value `nll(0)` (see above), including cells silent
  because of the readout threshold.
- The **cell-existence mask** works as in the χ² path: cells in a beam
  hole or outside the calorimeter are neither expected nor punished, and
  clusters at edges simply have fewer informative neighbors.

**Table content**

- A table with **missing (empty) nodes must never be shipped**: an empty
  node is indistinguishable from "beyond range" at run time and attracts
  the fit. Leave incomplete energies out of the asset — see
  [generating likelihood tables](./generating-ll-tables) for this and the
  other preparation-time pitfalls (crack nodes, zero-probability checks).

**Interpreting the output**

- The `chi2` field of a cluster carries `nll/ndof` under this method — a
  different scale from χ²/ndof. Any analysis cut on fit quality must be
  retuned when switching methods, and clusters reconstructed by the two
  methods must not be mixed in one quality cut.
- Multi-peak clusters report the χ²-metric fit (see limitations below),
  so within one event the `chi2` field can legitimately mix the two
  scales between single-peak and multi-peak clusters.

## Known limitations (v1) and planned improvements

- **Multi-peak clusters keep the χ² fit.** The per-peak *shared* cell
  energies produced by the iterative profile sharing are deterministic
  profile fractions, not single-shower samples — the tables do not
  describe them, and using them anyway degraded split-pair positions 2×.
  Single-peak clusters (the vast majority) use the full LL fit.
- **Overlapping showers are scored per peak** (partner-subtracted), not
  with a joint two-shower likelihood; the joint fit is the planned cure
  for both the shared-cluster fallback above and the residual split-pair
  position deficit.
- **The energy is not likelihood-refined** — it stays the moments +
  containment estimate (deliberate for now; refinement is a possible
  extension).
- **No interpolation between energy tables** (nearest-in-log-energy only);
  interpolating `nll` between the bracketing tables is a recorded
  improvement for the tuning phase.
- **Table nodes at offsets ≡ 50 (mod 100)** are neighbor interpolations —
  a training-geometry artifact, see
  [generating likelihood tables](./generating-ll-tables).
- **`split_threshold` default (0.5/dof) is provisional** pending tuning;
  choose it from your wrong-split budget
  ([overlapping showers](./two-gamma)).
