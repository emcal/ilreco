"""Shared fixtures for the ilreco Python binding tests.

The binding is validated against the SAME golden tables as the C library
(tests/data at the repository root): every (input_*.csv, expected_*.csv)
pair is replayed through the Python layer. The CSV files keep the historical
1-based columns/rows/positions; conversion to the 0-based API happens here,
mirroring tests/common/csv_io.h.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).parents[2]
GOLDEN_DIR = REPO_ROOT / "tests" / "data"
PROFILE_PWO = REPO_ROOT / "data" / "prof_pwo.dat"

GOLDEN_GRID = 34   # the golden fixtures were recorded on a 34x34 grid


def golden_set_names() -> list[str]:
    return sorted(p.name[len("input_"):-len(".csv")]
                  for p in GOLDEN_DIR.glob("input_*.csv"))


def read_input_table(name: str) -> np.ndarray:
    """Golden input as a 4-column (event, col, row, e) 0-based table."""
    rows = []
    with open(GOLDEN_DIR / f"input_{name}.csv") as file:
        for record in csv.DictReader(file):
            rows.append([int(record["event"]),
                         int(record["col"]) - 1,     # file is 1-based
                         int(record["row"]) - 1,
                         float(record["e"])])
    return np.array(rows, dtype=np.float64)


def read_expected(name: str) -> dict[int, list[dict]]:
    """Golden clusters per event, positions converted to 0-based."""
    events: dict[int, list[dict]] = {}
    with open(GOLDEN_DIR / f"expected_{name}.csv") as file:
        for record in csv.DictReader(file):
            event = int(record["event"])
            events.setdefault(event, [])
            if int(record["icl"]) >= 0:
                events[event].append({
                    "e": float(record["e"]),
                    "x": float(record["x"]) - 1.0,   # exact, see csv_io.h
                    "y": float(record["y"]) - 1.0,
                    "chi2": float(record["chi2"]),
                    "size": int(record["size"]),
                })
    return events


@pytest.fixture(scope="session")
def golden_calorimeter():
    import ilreco
    return ilreco.Calorimeter(GOLDEN_GRID, GOLDEN_GRID, profile=PROFILE_PWO)
