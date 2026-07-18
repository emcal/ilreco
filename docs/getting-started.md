# Getting started

## Quick start (C)

```c
#include <ilreco.h>

char err[256];
/* once per process: geometry + shower-profile tables, immutable afterwards */
ilreco_config *cfg = ilreco_config_create(/*n_cols=*/30, /*n_rows=*/30,
                                          "prof_pwo.dat", err, sizeof err);
if (!cfg) { fprintf(stderr, "%s\n", err); return 1; }

/* once per thread: all working memory, one allocation */
ilreco_workspace *ws = ilreco_workspace_create(cfg);

/* per event: zero allocations */
ilreco_hit hits[] = { {12, 14, 1.85}, {13, 14, 0.42}, {12, 15, 0.31} };
ilreco_cluster out[16];
int n = ilreco_reconstruct(ws, hits, 3, out, 16);
for (int i = 0; i < n && i < 16; ++i)
    printf("E=%.3f GeV at (%.2f, %.2f) cells, chi2=%.2f\n",
           out[i].e, out[i].x, out[i].y, out[i].chi2);

ilreco_workspace_destroy(ws);
ilreco_config_destroy(cfg);
```

Cells are 0-based `(col, row)`, energies in GeV, positions in cell units
(`x == col` means the center of that column). Clusters come back
energy-descending; the return value is the total found and `-1` flags
invalid input.

## Tutorial: an 8×8 calorimeter with a central 2×2 beam hole

Step-by-step example. The detector: an 8×8 crystal array with the central
2×2 block removed for the beam (the classic layout). The event: three
photons, landing on a hole-border cell, a normal interior cell, and an
outer-edge cell — every cell class appears, each with its own cluster in
the output.

![tutorial event](./tutorial-8x8.png)

### Step 1 — create the config (once per process)

```c
char err[256];
ilreco_config *cfg = ilreco_config_create(/*n_cols=*/8, /*n_rows=*/8,
                                          "prof_pwo.dat", err, sizeof err);
if (!cfg) { fprintf(stderr, "%s\n", err); return 1; }
```

The config holds the geometry and the shower-profile tables; after this
call it is immutable and shared by every thread.

### Step 2 — which cells are normal, edge, hole border, hole

You do **not** register any of this explicitly:

- **Normal (interior) cells** — the 20 white cells: a cluster seeded there
  is labeled `type = 0`.
- **Outer edge** — the boundary ring (col 0/7, row 0/7): no phantom
  neighbors are assumed beyond it; a cluster *seeded* there gets `type = 2`.
- **Hole border** — the 12 cells surrounding the central 2×2: seeds there
  get **`type = 1`**, automatically. The built-in classification encodes
  this layout (a central 2×2 hole). For a hole of any *other* shape or
  position, use the [cell-existence mask](./python-api#arbitrary-shapes)
  instead — labels then derive from the mask.
- **The hole is four cells you never feed.** (3,3), (4,3), (3,4), (4,4)
  exist in the geometry but no hit ever arrives for them, so clustering
  treats them as silent (zero-signal) neighbors. Energy deposited there is
  physically lost; the profile-based containment correction compensates on
  average, the same mechanism as at the outer edge. (With a cell mask the
  treatment is stricter: the cells are known absent and the correction is
  exact-by-construction rather than measured-zero.)

### Step 3 — create a workspace (once per thread)

```c
ilreco_workspace *ws = ilreco_workspace_create(cfg);
```

One allocation, sized for the 8×8, reused for every event this thread
processes. The workspace remembers which config it belongs to, so from
here on every call takes only `ws`. One workspace per thread, never shared
(see the [threading model](./threading)).

### Step 4 — put the event into the input array

Each fired cell becomes one `{col, row, energy}` entry — 0-based indices,
GeV, **no entries for the hole**. The array goes in **raw DAQ order**: you
do not sort it and you do not mark seeds — the library orders hits
internally and finds peak (seed) cells itself:

```c
ilreco_hit hits[15] = {
    {1, 1, 0.120}, {5, 4, 0.950}, {0, 6, 0.520},
    {2, 2, 0.090}, {6, 4, 0.140}, {2, 1, 0.780},
    {0, 5, 0.070}, {5, 5, 0.130}, {3, 1, 0.100},
    {1, 6, 0.080}, {5, 3, 0.110}, {2, 0, 0.070},
    {6, 5, 0.040}, {1, 7, 0.030}, {6, 3, 0.035},
};   /* hole cells (3,3) (4,3) (3,4) (4,4): no entries */

ilreco_cluster out[8];
int n = ilreco_reconstruct(ws, hits, 15, out, 8);
```

### Step 5 — interpret the output

Running the above gives `n = 3`; clusters arrive energy-descending, and
this event yields one of each type:

| | `e` [GeV] | `x`, `y` [cells] | `chi2` | `size` | `type` | interpretation |
|---|---|---|---|---|---|---|
| `out[0]` (A) | **1.4239** | 5.278, 4.047 | 1.38 | 6 | **1** | seeded on the **hole border**; E above the 1.405 GeV input sum (containment correction includes what leaked into the silent hole) |
| `out[1]` (B) | **1.1789** | 1.949, 1.051 | 1.79 | 5 | **0** | **normal** interior shower |
| `out[2]` (C) | **0.7450** | 0.325, 5.815 | 0.79 | 4 | **2** | seeded on the **outer edge**; strongest correction (0.745 vs Σ 0.700) — part of the shower left the detector |

Convert to millimeters for an array centered on the origin with pitch `p`:
`x_mm = (x - (n_cols-1)/2.0) * p` → photon A:
`(5.278 - 3.5) * 20.9 = +37.2 mm`.

### Step 6 — multiple clusters and merged showers

The return value `n` is the number of clusters found (it can exceed
`max_out`; `out[]` receives at most `max_out`, highest energies first).
Always loop:

```c
for (int i = 0; i < n && i < 8; ++i) { /* ... out[i] ... */ }
```

What `n` means physically — three regimes, on the same 8×8:

1. **Well-separated showers** (≥1 empty cell between them): distinct
   islands → one clean cluster each, nothing shared.
2. **Touching showers with distinct maxima** (seeds two cells apart, 1.0
   and 0.9 GeV with a 0.3 GeV valley): ONE island, but peak finding sees
   two local maxima and profile-fits both: `n = 2`,
   `E = 1.453 / 1.303` — and **both report `size = 7`**: every cell of the
   island belongs to both clusters, its energy *shared* between them
   according to the shower profile. Cluster sizes do not add up to the hit
   count in this regime.
3. **Genuinely merged** (maxima on adjacent cells): one maximum survives →
   `n = 1`, `E = 2.402` (the pair summed), and **`chi2 = 10.5`** against
   ~0.5–2 for a clean single shower. A one-cell separation is below the
   method's resolving power; cut on `chi2` to flag such candidates.
