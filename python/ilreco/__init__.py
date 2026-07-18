"""ilreco — island-clustering reconstruction for square-cell calorimeters.

Table in, table out:

    import ilreco
    calo = ilreco.Calorimeter(20, 20, profile="pwo")
    clusters = calo.reconstruct(hits)

``hits`` is anything numpy can view as a table (numpy array, pandas
DataFrame, list of rows):

* 3 columns -> ``col, row, e``  — a single event;
* 4 columns -> ``event_id, col, row, e`` — many events; rows of one event
  must be contiguous (grouped, as data from any store or DataFrame is).

The result is a numpy structured array with fields
``event, e, x, y, chi2, size, type`` (one row per cluster, energy-descending
within each event; a single-event input comes back with ``event = 0``).
Positions are in 0-based cell units — ``x == col`` means the center of that
column. See the documentation for conversion to millimeters.
"""
from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from ilreco._core import Calorimeter as _CoreCalorimeter

__all__ = ["Calorimeter", "CLUSTER_DTYPE"]
__version__ = "1.0.0"

CLUSTER_DTYPE = np.dtype([
    ("event", np.int64),
    ("e", np.float64),      # GeV, containment-corrected
    ("x", np.float64),      # 0-based cell units
    ("y", np.float64),
    ("chi2", np.float64),   # profile-fit chi2 / ndof
    ("size", np.int32),     # hits in the cluster
    ("type", np.int32),     # 0 interior, 1 hole border/rim, 2 outer boundary
])

# Column names accepted from DataFrames, in the meaning they map to.
_EVENT_NAMES = ("event", "event_id", "evt")
_COL_NAMES = ("col", "cols", "ix", "column")
_ROW_NAMES = ("row", "rows", "iy")
_ENERGY_NAMES = ("e", "energy", "e_gev", "edep")

# Serial reconstruction cost scale used by n_jobs="auto" (order of magnitude
# of the measured single-shower benchmark; only the threshold logic uses it).
_APPROX_SECONDS_PER_EVENT = 40e-6
_AUTO_PARALLEL_THRESHOLD_SECONDS = 20e-3
_MIN_EVENTS_PER_WORKER = 200


def _profile_path(profile: str | os.PathLike) -> str:
    """Resolve 'pwo'/'lg' to the shipped tables, or pass a file path through."""
    shortcuts = {"pwo": "prof_pwo.dat", "lg": "prof_lg.dat"}
    if isinstance(profile, str) and profile in shortcuts:
        packaged = Path(__file__).parent / "data" / shortcuts[profile]
        if packaged.exists():
            return str(packaged)
        repo_data = Path(__file__).parents[2] / "data" / shortcuts[profile]
        if repo_data.exists():
            return str(repo_data)
        raise FileNotFoundError(
            f"shipped profile '{profile}' not found (looked in {packaged} "
            f"and {repo_data})")
    return os.fspath(profile)


def _central_hole_mask(n_cols: int, n_rows: int, hole_size: int) -> np.ndarray:
    """All-cells-exist mask with a centered hole_size x hole_size hole.

    A centered square hole only exists when grid and hole parities match
    (even hole on even grid, odd hole on odd grid) in BOTH dimensions;
    anything else is an off-center hole and needs an explicit mask.
    """
    if hole_size < 1:
        raise ValueError("hole_size must be >= 1")
    if (n_cols - hole_size) % 2 != 0 or (n_rows - hole_size) % 2 != 0:
        raise ValueError(
            f"a {hole_size}x{hole_size} hole cannot be centered on a "
            f"{n_cols}x{n_rows} grid (parity mismatch) — pass an explicit "
            f"mask= instead")
    if hole_size >= min(n_cols, n_rows):
        raise ValueError("hole_size leaves no cells")
    mask = np.ones((n_rows, n_cols), dtype=np.uint8)
    col0 = (n_cols - hole_size) // 2
    row0 = (n_rows - hole_size) // 2
    mask[row0:row0 + hole_size, col0:col0 + hole_size] = 0
    return mask


def _as_table(hits) -> np.ndarray:
    """Coerce the input to a 2-D float64 table, using DataFrame column names
    when they are recognizable."""
    if hasattr(hits, "columns") and hasattr(hits, "to_numpy"):   # DataFrame
        names = [str(c).lower() for c in hits.columns]
        def pick(candidates):
            for candidate in candidates:
                if candidate in names:
                    return names.index(candidate)
            return None
        event_i = pick(_EVENT_NAMES)
        col_i = pick(_COL_NAMES)
        row_i = pick(_ROW_NAMES)
        energy_i = pick(_ENERGY_NAMES)
        if col_i is not None and row_i is not None and energy_i is not None:
            order = ([event_i] if event_i is not None else []) + \
                    [col_i, row_i, energy_i]
            hits = hits.iloc[:, order]
        hits = hits.to_numpy()
    table = np.asarray(hits, dtype=np.float64)
    if table.ndim != 2 or table.shape[1] not in (3, 4):
        raise ValueError(
            "hits must be a table with 3 columns (col, row, e) or 4 columns "
            f"(event_id, col, row, e); got shape {table.shape}")
    return table


def _event_starts(event_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Offsets of contiguous event blocks + the id of each block.

    Rows of one event must be contiguous; a repeated id after other events
    means the table is not grouped and is rejected.
    """
    change = np.flatnonzero(event_ids[1:] != event_ids[:-1]) + 1
    starts = np.concatenate(([0], change, [len(event_ids)]))
    block_ids = event_ids[starts[:-1]]
    if len(np.unique(block_ids)) != len(block_ids):
        raise ValueError(
            "rows of one event must be contiguous — group the table by "
            "event_id first")
    return starts.astype(np.int64), block_ids


class Calorimeter:
    """A calorimeter ready to reconstruct events.

    Geometry is given one of three ways:

    * ``Calorimeter(n_cols, n_rows, ...)`` — full rectangle, no hole;
    * ``Calorimeter(n_cols, n_rows, hole_size=2, ...)`` — centered square
      beam hole; hole cells are removed via the cell mask, so energy leaking
      into the hole is corrected for and clusters at its border get
      ``type = 1``. The hole must be centrable (parity match), otherwise a
      ValueError asks for an explicit mask;
    * ``Calorimeter(mask=mask_2d, ...)`` — arbitrary shape: 2-D boolean/int
      array, nonzero = cell exists; grid size is the mask shape.

    ``profile`` is ``"pwo"`` (2.05 cm PbWO4, default), ``"lg"`` (lead
    glass), or a path to your own table (see the profile documentation).
    """

    def __init__(self, n_cols: int | None = None, n_rows: int | None = None,
                 *, mask=None, hole_size: int | None = None,
                 profile: str | os.PathLike = "pwo",
                 seed_threshold_gev: float | None = None,
                 zcal_cm: float | None = None):
        if mask is not None:
            mask = np.asarray(mask)
            if mask.ndim != 2:
                raise ValueError("mask must be a 2-D array (rows x cols)")
            if hole_size is not None:
                raise ValueError("give either mask= or hole_size=, not both")
            if n_cols is not None or n_rows is not None:
                if (n_rows, n_cols) != mask.shape:
                    raise ValueError(
                        f"mask shape {mask.shape} contradicts n_cols/n_rows "
                        f"({n_cols}, {n_rows}); grid size is taken from the "
                        "mask — drop the explicit sizes")
            n_rows, n_cols = mask.shape
        elif n_cols is None or n_rows is None:
            raise ValueError("give n_cols and n_rows, or a mask=")

        self._core = _CoreCalorimeter(int(n_cols), int(n_rows),
                                      _profile_path(profile))

        if hole_size is not None:
            mask = _central_hole_mask(int(n_cols), int(n_rows), int(hole_size))
        if mask is not None:
            self._core.set_cell_mask(
                np.ascontiguousarray(mask != 0, dtype=np.uint8))
        else:
            # full rectangle: there is no hole, so the built-in central-2x2
            # labeling would mislabel interior cells
            self._core.set_hole_classification(0)
        if seed_threshold_gev is not None:
            self._core.set_seed_threshold(float(seed_threshold_gev))
        if zcal_cm is not None:
            self._core.set_zcal(float(zcal_cm))
        self.n_cols = int(n_cols)
        self.n_rows = int(n_rows)

    def reconstruct(self, hits, *, n_jobs="auto") -> np.ndarray:
        """Reconstruct one event (3-column table) or many (4-column table).

        ``n_jobs``: any falsy value or 1 = serial; an integer = that many
        threads; ``"auto"`` = serial for small jobs, otherwise splits the
        events over threads (the GIL is released inside each chunk, so the
        threads genuinely run in parallel).
        """
        table = _as_table(hits)
        if table.shape[1] == 3:
            event_ids = np.zeros(len(table), dtype=np.int64)
            cols, rows, energies = (table[:, 0], table[:, 1], table[:, 2])
        else:
            event_ids = table[:, 0].astype(np.int64)
            cols, rows, energies = (table[:, 1], table[:, 2], table[:, 3])

        cols = cols.astype(np.int32)
        rows = rows.astype(np.int32)
        energies = np.ascontiguousarray(energies)
        if len(table) == 0:
            return np.empty(0, dtype=CLUSTER_DTYPE)

        starts, block_ids = _event_starts(event_ids)
        n_events = len(block_ids)
        n_workers = self._resolve_n_jobs(n_jobs, n_events)

        if n_workers <= 1:
            chunk_results = [self._run_chunk(cols, rows, energies, starts, 0)]
        else:
            # split at event boundaries into one chunk per worker
            boundaries = np.linspace(0, n_events, n_workers + 1, dtype=np.int64)
            jobs = []
            for w in range(n_workers):
                first_event, last_event = boundaries[w], boundaries[w + 1]
                if first_event == last_event:
                    continue
                lo = starts[first_event]
                hi = starts[last_event]
                chunk_starts = starts[first_event:last_event + 1] - lo
                jobs.append((cols[lo:hi], rows[lo:hi], energies[lo:hi],
                             chunk_starts, int(first_event)))
            with ThreadPoolExecutor(max_workers=len(jobs)) as pool:
                chunk_results = list(pool.map(
                    lambda job: self._run_chunk(*job), jobs))

        n_clusters = sum(len(r[0]) for r in chunk_results)
        out = np.empty(n_clusters, dtype=CLUSTER_DTYPE)
        cursor = 0
        for local_event, e, x, y, chi2, size, cluster_type in chunk_results:
            block = out[cursor:cursor + len(e)]
            block["event"] = block_ids[local_event]
            block["e"] = e
            block["x"] = x
            block["y"] = y
            block["chi2"] = chi2
            block["size"] = size
            block["type"] = cluster_type
            cursor += len(e)
        return out

    def _run_chunk(self, cols, rows, energies, chunk_starts, first_event):
        """One serial C call over a contiguous slice of events."""
        (local_event, e, x, y, chi2, size, cluster_type) = \
            self._core.reconstruct_batch(cols, rows, energies, chunk_starts)
        return (local_event + first_event, e, x, y, chi2, size, cluster_type)

    def _resolve_n_jobs(self, n_jobs, n_events: int) -> int:
        if not n_jobs or n_jobs == 1:
            return 1
        if n_jobs == "auto":
            estimated_seconds = n_events * _APPROX_SECONDS_PER_EVENT * \
                max(1.0, (self.n_cols * self.n_rows) / 1000.0)
            if estimated_seconds < _AUTO_PARALLEL_THRESHOLD_SECONDS:
                return 1
            return max(1, min(os.cpu_count() or 1,
                              n_events // _MIN_EVENTS_PER_WORKER))
        return max(1, int(n_jobs))
