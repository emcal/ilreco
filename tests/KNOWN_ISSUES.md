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

## Re-scan of the rewritten code (2026-07-17)

Fresh review of post-refactor ilreco.c (manual line-by-line; CI carries
cppcheck/clang-tidy; ASan/UBSan clean; all 4 ctest suites green). New findings,
none affecting physics output — documented, not fixed:

11. **`ilreco_reconstruct` does not validate `out`/`max_out`**: `out == NULL`
    with `max_out > 0` crashes; `max_out < 0` is accepted and returns nadcgam
    with nothing written (caller may read uninitialized memory). Hits pointer
    and geometry ARE validated.
12. **Inconsistent NULL handling in setters**: `ilreco_config_set_cell_mask`
    checks `cfg` for NULL and returns -1; `set_zcal` / `set_seed_threshold` /
    `set_hole_classification` dereference without checking.
13. **The library prints to stdout on internal conditions** (preserved from
    the original): "maximum number of clusters reached" (cluster_search also
    silently truncates to mcl-1), "WRN: lost maximum", "WARNING NEGATIVE CORR",
    "case 0 ch". A drop-in library should report via counters/return codes.
14. ~~`order[_MADCGAM_]` fixed-size local in `ilreco_reconstruct`~~ FIXED
    2026-07-17: moved into the workspace arena, sized by cfg->madcgam (golden
    tables unchanged). The `ipnpk/epk/xpk/ypk/fw[_MPK_]` locals stay on the
    stack deliberately: _MPK_ is an algorithm constant (max peaks per island,
    bounded by the npk break at _MPK_-2, verified), not geometry.
15. **Exported non-static helper symbols pollute the namespace**: `nint`,
    `ZBQLINI`, `ZBQLU01`, `ZBQLUAB`, `read_event`, `dump_clusters` — link-time
    collision risk for a drop-in library (should be static or prefixed).
16. **Packed-address `int` overflow for absurd grids**: `(col+1)*offset` can
    exceed INT_MAX around ~46000x46000 cells — theoretical, no guard.
17. Cosmetic, preserved: `chisq1_cluser` (typo in name), "does not exists" in
    the legacy error message; `read_profile_data` reports an unreadable file
    twice (access() check + ensure_legacy).

## Legacy removal (2026-07-17)

The legacy API (read_profile_data / cluster_search / process_cluster), the
built-in test event generator (read_event, ZBQL* random generator,
dump_clusters, the opt-in main()), the src/test/ilreco_dump tool, and the
legacy compile-time geometry macros were REMOVED — the context API is the
only interface. adcgam_t is now internal to ilreco.c. Interface and stored
integer types are exact-width (int32_t). Validation: full ctest green
(golden tables replay through the 0-based driver), ASan/UBSan clean, and the
cal-fpga snapshot protocol reports the rebuilt chain byte-identical to the
pre-removal snapshot (33 reco files, 14 number sets).

Status effect on the list above: #1 (main in library), #3 (CLUSTER_MIN_ENERGY)
and #4 (MIN_COUNTER_ENERGY) are gone with the removed code; #15 is resolved
(ZBQL*/read_event/dump_clusters removed, nint static); #12/#13 from the
re-scan and #2 (elfract units, now internal-only) remain open. The default
peak budget (madr0) is now 100*n_cols*n_rows — far above any physical
occupancy, equivalent on all validated data; the guard binds only for
pathological >= _MPK_-2-peak islands.

## Stage-3 note (cell-existence mask, 2026-07-16)

Without a mask, cells beyond a physical rim or inside a beam hole (any cell
that exists in the grid but never fires) are treated as measured-zero
neighbors: they enter fill_zero_hits and the chi2 as real measurements, so
clusters hugging a hole/rim are slightly under-corrected. This behavior is
kept bit-identical (golden tables) as the default. With
`ilreco_config_set_cell_mask` the missing cells are excluded from the
zero-neighbor treatment and the chi2, `type` labels derive from the mask, and
hits on nonexistent cells are rejected — see unit/test_cell_mask.cpp.
