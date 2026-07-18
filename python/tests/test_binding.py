"""Tests of the ilreco Python binding.

The core physics gate is the golden replay: the binding must reproduce the
same frozen event responses that gate the C library. The rest covers the
binding's own contracts — table dispatch, grouping requirement, hole/mask
constructor logic, threading determinism, and error reporting.
"""
from __future__ import annotations

import numpy as np
import pytest

import ilreco
from conftest import (GOLDEN_GRID, PROFILE_PWO, golden_set_names,
                      read_expected, read_input_table)


@pytest.mark.parametrize("set_name", golden_set_names())
def test_golden_tables_replay(golden_calorimeter, set_name):
    table = read_input_table(set_name)
    expected = read_expected(set_name)
    clusters = golden_calorimeter.reconstruct(table, n_jobs=1)

    for event, expected_clusters in expected.items():
        found = clusters[clusters["event"] == event]
        assert len(found) == len(expected_clusters), \
            f"set {set_name} event {event}: cluster count"
        for got, want in zip(found, expected_clusters):
            assert got["e"] == pytest.approx(want["e"], rel=1e-6)
            assert got["x"] == pytest.approx(want["x"], abs=1e-5)
            assert got["y"] == pytest.approx(want["y"], abs=1e-5)
            assert got["size"] == want["size"]


def test_n_jobs_gives_identical_results(golden_calorimeter):
    table = np.concatenate([read_input_table(name)
                            for name in golden_set_names()[:3]])
    # re-key events so ids stay unique across the concatenated sets
    boundaries = np.flatnonzero(np.diff(table[:, 0]) != 0) + 1
    new_ids = np.zeros(len(table))
    new_ids[boundaries] = 1
    table[:, 0] = np.cumsum(new_ids)

    serial = golden_calorimeter.reconstruct(table, n_jobs=1)
    threaded = golden_calorimeter.reconstruct(table, n_jobs=8)
    automatic = golden_calorimeter.reconstruct(table, n_jobs="auto")
    assert np.array_equal(serial, threaded)
    assert np.array_equal(serial, automatic)


def test_single_event_three_columns():
    calo = ilreco.Calorimeter(20, 20, profile=PROFILE_PWO)
    clusters = calo.reconstruct([[12, 14, 1.85], [13, 14, 0.42], [12, 15, 0.31]])
    assert len(clusters) == 1
    assert clusters["event"][0] == 0
    assert round(clusters["x"][0]) == 12
    assert round(clusters["y"][0]) == 14
    assert clusters["e"][0] >= 1.85


def test_empty_input():
    calo = ilreco.Calorimeter(10, 10, profile=PROFILE_PWO)
    clusters = calo.reconstruct(np.empty((0, 3)))
    assert len(clusters) == 0
    assert clusters.dtype == ilreco.CLUSTER_DTYPE


def test_hole_size_labels_and_rejection():
    calo = ilreco.Calorimeter(8, 8, hole_size=2, profile=PROFILE_PWO)
    # shower on the hole border: labeled type 1, energy corrected above raw sum
    shower = [[5, 4, 0.95], [6, 4, 0.14], [5, 5, 0.13], [5, 3, 0.11]]
    clusters = calo.reconstruct(shower)
    assert clusters["type"][0] == 1
    assert clusters["e"][0] > sum(hit[2] for hit in shower)
    # a hit on a hole cell is a caller bug, reported loudly
    with pytest.raises(ValueError, match="masked-out cell"):
        calo.reconstruct([[3, 3, 1.0]])


def test_hole_size_parity_validation():
    with pytest.raises(ValueError, match="explicit mask"):
        ilreco.Calorimeter(5, 5, hole_size=2, profile=PROFILE_PWO)
    with pytest.raises(ValueError, match="explicit mask"):
        ilreco.Calorimeter(4, 4, hole_size=1, profile=PROFILE_PWO)
    # matching parity works
    ilreco.Calorimeter(5, 5, hole_size=1, profile=PROFILE_PWO)
    ilreco.Calorimeter(8, 8, hole_size=2, profile=PROFILE_PWO)


def test_arbitrary_mask_types():
    # circular 9x9 with an off-center hole, as in the C mask tests
    mask = np.zeros((9, 9), dtype=bool)
    for row in range(9):
        for col in range(9):
            if (col - 4) ** 2 + (row - 4) ** 2 <= 18:
                mask[row, col] = True
    mask[4, 5] = False   # beam hole
    mask[5, 5] = False
    calo = ilreco.Calorimeter(mask=mask, profile=PROFILE_PWO)
    assert (calo.n_cols, calo.n_rows) == (9, 9)
    assert calo.reconstruct([[4, 4, 1.0]])["type"][0] == 1   # next to hole
    assert calo.reconstruct([[4, 0, 1.0]])["type"][0] == 2   # bbox ring
    assert calo.reconstruct([[2, 4, 1.0]])["type"][0] == 0   # interior


def test_mask_and_sizes_must_agree():
    mask = np.ones((6, 7), dtype=np.uint8)
    calo = ilreco.Calorimeter(mask=mask, profile=PROFILE_PWO)
    assert (calo.n_cols, calo.n_rows) == (7, 6)
    with pytest.raises(ValueError, match="grid size is taken from the mask"):
        ilreco.Calorimeter(5, 5, mask=mask, profile=PROFILE_PWO)
    with pytest.raises(ValueError, match="not both"):
        ilreco.Calorimeter(8, 8, mask=np.ones((8, 8)), hole_size=2,
                           profile=PROFILE_PWO)


def test_ungrouped_events_rejected():
    calo = ilreco.Calorimeter(10, 10, profile=PROFILE_PWO)
    table = [[0, 4, 4, 1.0], [1, 6, 6, 1.0], [0, 5, 4, 0.2]]
    with pytest.raises(ValueError, match="contiguous"):
        calo.reconstruct(table)


def test_bad_shapes_rejected():
    calo = ilreco.Calorimeter(10, 10, profile=PROFILE_PWO)
    with pytest.raises(ValueError, match="3 columns"):
        calo.reconstruct(np.ones((4, 5)))
    with pytest.raises(ValueError, match="outside the 10x10 grid"):
        calo.reconstruct([[10, 0, 1.0]])


def test_dataframe_with_named_columns():
    pandas = pytest.importorskip("pandas")
    calo = ilreco.Calorimeter(10, 10, profile=PROFILE_PWO)
    frame = pandas.DataFrame({
        "e": [1.0, 0.3],
        "row": [4, 4],
        "event_id": [7, 7],
        "col": [4, 5],
    })   # deliberately scrambled column order: names must win
    clusters = calo.reconstruct(frame)
    assert clusters["event"][0] == 7
    assert 4.0 < clusters["x"][0] < 4.5


def test_output_positions_are_zero_based():
    calo = ilreco.Calorimeter(10, 10, profile=PROFILE_PWO)
    clusters = calo.reconstruct([[0, 0, 1.0]])
    assert abs(clusters["x"][0]) < 0.5
    assert abs(clusters["y"][0]) < 0.5
