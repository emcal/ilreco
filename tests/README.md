# ilreco test suite

The safety net around the library: any code change must either reproduce these
outputs exactly or consciously re-baseline them (and commit the new tables
together with the change that caused them).

Build & run (Catch2 v3 fetched automatically):

```bash
cmake -B build && cmake --build build -j
cd build && ctest --output-on-failure
```

## Infrastructure map — where everything comes from

| provides | what |
|---|---|
| `tests/CMakeLists.txt` | the four ctest targets, and two compile definitions visible in every test: `ILRECO_PROF_PWO` = absolute path of the shipped shower-profile file `data/prof_pwo.dat`; `ILRECO_TEST_DATA` = absolute path of the golden fixture directory `tests/data` |
| `common/ilreco_test_utils.h` | namespace `ilt::` — everything the tests share: `GRID_COLS`/`GRID_ROWS` (the 34×34 default test grid), `TestContext` (RAII config + workspace), `default_context()`, `run_event()` (reconstruct one event on the shared context), `finite()`, `Rng` (deterministic LCG), `synth_shower()` (generated shower events) |
| `common/csv_io.h` | golden-table CSV reading/writing; converts the frozen 1-based file format to the library's 0-based convention at the file boundary — the only place that conversion exists |
| Catch2 v3 | fetched by CMake `FetchContent` at configure time; test binaries link `Catch2::Catch2WithMain` (it supplies `main`) |

House rule: a test file includes only `ilreco.h` (via the helpers), Catch2
headers, and the two `common/` headers. Anything shared lives in `common/`;
nothing is defined in one test file and used from another.

## Categories

| dir | ctest name | what it protects |
|---|---|---|
| `unit/` | `unit` | API contracts: addressing round-trips (incl. corners), cluster counting, order-invariance, energy scaling/envelope, finiteness, error handling, cell-existence mask, multithreading (bitwise vs serial). Tests contracts, not internals — they survive refactoring. Tags: `[contract]`, `[edge]`, `[context]`, `[mask]`, `[threads]`, `[documented]` (= pins observed behavior whose *intent* is unclear; check KNOWN_ISSUES.md before "fixing"). |
| `smoke/` | `smoke` | thousands of generated events (single/pairs/noise/whole-grid/state-leak check); asserts only crash-freedom and coarse sanity. First net against memory corruption. |
| `data/` | `data` | golden event-by-event regression tables (`input_*.csv` + `expected_*.csv`, 660 events, 364 KB) frozen from Geant4-simulated detector data: central showers at 1/5/20 GeV, TRUE grid corners/edges (leakage response), two-particle overlaps, pi- topologies, synthetic exotics. The files keep their historical 1-based format (see `common/csv_io.h`). A mismatch = the physics output changed. |
| `benchmark/` | `benchmark` | `ilreco_bench [--reps N]` prints JSON (µs/event by category); `baseline.json` is the reference. Compare before/after any change. |
| `tools/` | — | `ilreco_gen_golden input.csv expected.csv` re-freezes tables. Only run when a physics change is intended, and commit tables with that change. |

## Fixture provenance

Inputs were produced from the standalone Geant4 PbWO4 simulation (halld-matched
response chain; currently github-bound as a separate emcal repo, today in the
cal-fpga workspace — `scripts/750_make_ilreco_fixtures.py`), translated onto the
34x34 test grid; edge/corner sets deliberately shift showers so part of the energy
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
tested in isolation: they are not reachable through a stable public seam, and
pinning them would weld tests to the implementation. They are covered
end-to-end by the golden tables.
