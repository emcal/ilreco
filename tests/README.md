# ilreco test suite

Safety net built BEFORE any modernization of the library (0-based indexing, dynamic
allocation are planned): the goal is that any future code change either reproduces
these outputs or consciously re-baselines them. Library sources are untouched by
this suite; only CMake files were edited (version bumps + new targets).

Build & run (Catch2 v3 fetched automatically):

```bash
cmake -B build && cmake --build build -j
cd build && ctest --output-on-failure
```

## Categories

| dir | ctest name | what it protects |
|---|---|---|
| `unit/` | `unit` | API contracts: addressing round-trips (incl. corners), cluster counting, order-invariance, energy scaling/envelope, finiteness. Tests contracts, not internals — they should survive refactoring. Tags: `[contract]`, `[edge]`, `[documented]` (= pins observed behavior whose *intent* is unclear; check KNOWN_ISSUES.md before "fixing"). |
| `smoke/` | `smoke` | thousands of generated events (single/pairs/noise/whole-grid/state-leak check); asserts only crash-freedom and coarse sanity. First net against memory corruption. |
| `data/` | `data` | golden event-by-event regression tables (input_*.csv + expected_*.csv, 660 events, 364 KB) frozen from Geant4-simulated detector data: central showers at 1/5/20 GeV, TRUE grid corners/edges (leakage response), two-particle overlaps, pi- topologies, synthetic exotics. A mismatch = the physics output changed. |
| `benchmark/` | `benchmark` | `ilreco_bench [--reps N]` prints JSON (us/event by category); `baseline.json` is the reference. Compare before/after any change. |
| `tools/` | — | `ilreco_gen_golden input.csv expected.csv` re-freezes tables. Only run when a physics change is intended, and commit tables with that change. |

## Fixture provenance

Inputs were produced from the standalone Geant4 PbWO4 simulation (halld-matched
response chain; currently github-bound as a separate emcal repo, today in the
cal-fpga workspace — `scripts/750_make_ilreco_fixtures.py`), translated onto the
stock 34x34 grid; edge/corner sets deliberately shift showers so part of the energy
falls off-grid. Golden outputs = the library's own responses at freeze time
(tolerances 1e-6 rel on E, 1e-5 abs on position absorb FP noise only).

## Simulation-reconstruction (physics) testing

End-to-end resolution validation (simulate -> reconstruct -> fit sigma_E/E and
position curves -> compare against frozen coefficients + plots) lives with the
simulation for now: cal-fpga `scripts/700_simreco_regression.py`. Plan of record:
the Geant4 setup becomes its own repo under github.com/emcal and the sim-reco
regression a third repo wiring both together.

## CI

`.github/workflows/ci.yml`: build + ctest (blocking), AddressSanitizer+UBSan test
run and cppcheck/clang-tidy static analysis (non-blocking, informational — see
KNOWN_ISSUES.md for triaged findings), benchmark artifact upload.

## What is deliberately NOT unit-tested

Internal stages (peak finding, chi2 iteration, profile interpolation) are not
tested in isolation: they are not reachable through a stable public seam without
modifying sources, and pinning them would weld tests to the implementation.
They are covered end-to-end by the golden tables. When the modernization refactor
introduces proper seams, unit tests can move inward.
