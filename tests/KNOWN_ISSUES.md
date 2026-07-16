# Known issues / findings (test-suite audit, 2026-07-16)

Found while building the test suite. Per the phase-4 ground rule, NONE of these are
fixed yet — the library sources are frozen until the safety net (this test suite) is
in place and agreed baselines exist. Each item has a documenting test or note.

1. **`main()` lives inside the library source** (`ilreco.c:89`). Works only because
   the library is built SHARED; static linking or LTO would clash with any host
   program. Downstream builds (cal-fpga) compile with `-Dmain=<renamed>`.
2. **`elfract` units are inconsistent across branches**: two-gamma splits store
   GeV-weighted fractions (`id[j]*e1/(e1+e2)`, ilreco.c:405), the single-gamma branch
   stores `nint(1e4*id[j])` fixed-point counts (ilreco.c:412), and the multi-peak path
   mixes both (ilreco.c:665-672). A consumer cannot interpret the field uniformly.
3. **`CLUSTER_MIN_ENERGY` (0.1 GeV) is declared in ilreco.h and never used.**
   No such cut exists in the code — comments/docs implying it are misleading.
4. **`MIN_COUNTER_ENERGY` (1 MeV) is only used by the built-in test generator**
   (`read_event`), not enforced by the library. Callers must apply their own hit
   threshold. (Documented in tests/unit/test_edge_inputs.cpp.)
5. **Undocumented cluster-seed threshold**: `process_cluster` requires a seed cell
   above `minpk = 0.01 GeV`, scaled by `7*log(1+E_cluster)` for clusters with >= 3
   hits (ilreco.c:299,313). Isolated deposits below ~10 MeV vanish silently. Real
   contract, nowhere documented before. (Frozen in tests/data/input_synthetic.csv
   and the seed-threshold unit test.)
6. **`read_profile_data` calls `exit(1)`** on a corrupt profile file (ilreco.c:291) —
   a library should report failure, not kill the host process.
7. **Duplicate cell addresses are undefined behavior by design** (upstream must
   merge). The library survives them (documented test) but the output is unspecified.
8. **Global mutable state / not thread-safe**: profile tables and work buffers are
   file-scope statics; one instance per process, single-threaded only.
9. **Geometry hardcoded**: `_NCOL_/_NROW_ = 34`, `_OFFSET_ = 100`, 1-based indexing,
   and `_ZCAL_ = 732` (target distance, used in the two-gamma invariant-mass
   separation cut) are compile-time constants. This is the agreed refactor target
   (0-based indexing + one-time dynamic allocation from a config) — NOT NOW.
10. Positive finding: the full suite (unit/smoke/data/benchmark) runs CLEAN under
    AddressSanitizer + UndefinedBehaviorSanitizer — no leaks (all-static design),
    no UB triggered on any tested input, including the whole-grid and >_MAXLEN_
    stress cases.

## Post-refactor notes (2026-07-16, context-API modernization)

Fixed by the refactor: #1 (main() now behind ILRECO_BUILD_TEST_MAIN), #6 (context
API returns errors; legacy shim still exits by design), #8 (context API is
thread-safe; legacy shim remains single-threaded), #9 (geometry is runtime config
in the context API; legacy macros remain as defaults), and the latent fill_zero_
hits overflow (workspace sizes iazero at 8*madr+2). Unchanged by design: #2, #3,
#4, #5 (now tunable via ilreco_config_set_seed_threshold), #7, plus the dead
tgamma_cluster stub and the chisq2t_hyc typo (preserved, still dead code).

New notes:
- ~5% throughput cost vs the macro-configured original (40.8 vs 38.9 us/event
  single-shower): packed-address divisions by a runtime `offset` cannot be
  strength-reduced by the compiler. Future optimization: carry (col,row)
  unpacked through the pipeline (planned for the 0-based internal refactor).
- The legacy API's hidden context is created on first use and intentionally
  never freed (process-lifetime); LeakSanitizer treats it as still-reachable.
