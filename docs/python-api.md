# Python API

```bash
pip install ilreco
```

Table in, table out:

```python
import ilreco

calo = ilreco.Calorimeter(20, 20, profile="pwo")
clusters = calo.reconstruct([[12, 14, 1.85], [13, 14, 0.42], [12, 15, 0.31]])
# array([(0, 2.633, 12.361, 14.322, 2.07, 3, 0)],
#       dtype=[('event','<i8'),('e','<f8'),('x','<f8'),('y','<f8'),
#              ('chi2','<f8'),('size','<i4'),('type','<i4')])
```

## Input tables

`reconstruct` accepts anything numpy can view as a 2-D table — a numpy
array, a pandas DataFrame, a list of rows:

- **3 columns** → `col, row, e` — a single event;
- **4 columns** → `event_id, col, row, e` — many events. Rows of one event
  must be **contiguous** (grouped, as data out of any store or DataFrame
  groupby already is); an interleaved table raises a `ValueError`.

Cells are 0-based, energies in GeV. A DataFrame with recognizable column
names (`event`/`event_id`, `col`/`ix`, `row`/`iy`, `e`/`energy`) is mapped
by name, in any column order:

```python
import pandas as pd
frame = pd.read_parquet("hits.parquet")        # event, ix, iy, e columns
clusters = calo.reconstruct(frame[["event", "ix", "iy", "e"]])
```

## Output

Always a structured numpy array with one row per cluster,
energy-descending within each event:

| field | dtype | meaning |
|---|---|---|
| `event` | int64 | the input `event_id` (0 for a 3-column input) |
| `e` | float64 | energy [GeV], containment-corrected |
| `x`, `y` | float64 | 0-based cell units; `x == col` = center of that column |
| `chi2` | float64 | profile-fit χ²/ndof |
| `size` | int32 | hits in the cluster |
| `type` | int32 | 0 interior, 1 hole border / rim, 2 outer boundary |

Convert positions to millimeters for an array centered on the origin:

```python
x_mm = (clusters["x"] - 0.5 * (calo.n_cols - 1)) * pitch_mm
```

`pd.DataFrame(clusters)` turns the result into a DataFrame directly.

## Geometry: rectangles, holes, arbitrary shapes {#arbitrary-shapes}

```python
ilreco.Calorimeter(20, 20)                 # full rectangle, no hole
ilreco.Calorimeter(10, 10, hole_size=2)    # centered 2x2 beam hole
ilreco.Calorimeter(mask=exists_2d)         # any shape; grid = mask.shape
```

- `hole_size=` removes a centered square block via the cell mask: energy
  leaking into the hole is corrected for, border clusters get `type = 1`,
  and a hit arriving on a hole cell raises `ValueError` (it is a caller
  bug, not an event). The hole must be centrable — a 2×2 hole on a 5×5
  grid raises and asks for an explicit mask.
- `mask=` takes a 2-D boolean/integer array, nonzero = the cell exists —
  circular calorimeters, asymmetric beam holes (EIC B0), dead regions.

Tunables mirror the C API: `Calorimeter(..., seed_threshold_gev=0.02,
zcal_cm=732.0)`.

## Threads and `n_jobs`

```python
clusters = calo.reconstruct(table)              # n_jobs="auto" (default)
clusters = calo.reconstruct(table, n_jobs=1)    # force serial
clusters = calo.reconstruct(table, n_jobs=8)    # force 8 threads
```

- The C core releases the GIL and takes a workspace from an internal pool,
  so concurrent calls — yours or the internal workers — are safe on any
  Python, including free-threaded builds.
- `"auto"` stays serial for small jobs (a handful of events is not worth a
  thread pool) and otherwise splits the events at event boundaries over a
  `concurrent.futures.ThreadPoolExecutor`.
- Results are **identical for every `n_jobs` value** (asserted bitwise in
  the test suite): chunking only splits the work, and the output order is
  restored at concatenation.

## Validation

The binding is gated by the same frozen golden tables as the C library:
every event of every golden set must reproduce its recorded response
(`python/tests/`).
