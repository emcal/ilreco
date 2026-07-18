# Development

## Repository layout

```
src/ilreco/     the library: ilreco.c, ilreco.h, ilreco.hpp, CMake export
data/           shipped shower-profile tables (prof_pwo.dat, prof_lg.dat)
tests/          C/C++ test suite (Catch2) + golden tables + benchmark
python/         Python binding: C extension, wrapper package, pytest suite
docs/           this documentation (VitePress)
```

## The safety net

Any code change must either reproduce the recorded outputs exactly or
consciously re-baseline them, committing the new tables together with the
change that caused them. The layers:

| suite | what it protects |
|---|---|
| `unit` | API contracts: addressing, cluster counting, order-invariance, energy envelope, error handling, cell mask, C++ wrapper, multithreading (bitwise vs serial) |
| `smoke` | thousands of generated events; crash-freedom and coarse sanity — the first net against memory corruption |
| `data` | golden event-by-event tables (660 events frozen from Geant4-simulated detector data); a mismatch means the physics output changed |
| `benchmark` | µs/event by category vs a committed `baseline.json` |
| python | the same golden tables replayed through `pip`-installed ilreco, plus binding contracts and `n_jobs` bitwise determinism |

```bash
cmake -B build && cmake --build build -j
cd build && ctest --output-on-failure     # unit / smoke / data / benchmark

pip install -e . && pytest python/tests   # the binding gate
```

`tests/README.md` carries the infrastructure map (where helpers and
compile definitions come from); `tests/KNOWN_ISSUES.md` records every
audited quirk and its status — read it before "fixing" a `[documented]`
test.

## API rules

- The C API in `ilreco.h` is the ground truth; `ilreco.hpp` and the Python
  layer add lifetimes and ergonomics, never semantics.
- Interface and stored integers are exact-width (`int32_t`).
- Every public function documents every parameter.
- The algorithm section of `ilreco.c` is intentionally
  expression-for-expression identical to the original PrimEx code:
  evaluation order and integer quantization are part of the contract,
  pinned bit-exact by the golden tables. Do not simplify arithmetic there
  without re-baselining.
- Every build of `ilreco.c` compiles with `-ffp-contract=off`: fused
  multiply-add contraction (the arm64 default) changes last-ulp results
  and can flip the ordering of near-degenerate clusters, breaking the
  cross-platform bit-identity that the golden tables enforce. If you
  embed the source directly into your own build, keep that flag.

## CI

GitHub Actions:

- **ci** — build + ctest on Linux and macOS (blocking), AddressSanitizer +
  UBSan run and cppcheck/clang-tidy (informational), benchmark artifact;
- **python** — build the wheel and run the pytest gate on Linux and macOS;
- **docs** — build this site and deploy it to GitHub Pages.
