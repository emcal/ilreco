// Shared helpers for the ilreco test suite — the ONLY place test
// infrastructure is defined (see tests/README.md, "Infrastructure map").
//
// Everything here speaks the library's own conventions (see ilreco.h):
// 0-based (col, row) cell indices, energies in GeV, cluster positions in
// 0-based cell units. Most tests run on one shared GRID_COLS x GRID_ROWS
// context — the golden fixtures in tests/data were recorded on that grid.
//
// ILRECO_PROF_PWO (used below) is a compile definition set in
// tests/CMakeLists.txt: the absolute path of the shipped shower-profile
// file data/prof_pwo.dat. Its sibling ILRECO_TEST_DATA (used by the golden
// tests) is the absolute path of the fixture directory tests/data.
#ifndef ILRECO_TEST_UTILS_H
#define ILRECO_TEST_UTILS_H

#include <ilreco.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace ilt {

constexpr int GRID_COLS = 34;
constexpr int GRID_ROWS = 34;

// One immutable config + one workspace. The default context below is shared
// by the single-threaded Catch2 tests; the multithreading test creates its
// own workspaces (one per thread), as real callers must.
struct TestContext {
    ilreco_config* config = nullptr;
    ilreco_workspace* workspace = nullptr;

    TestContext(int n_cols, int n_rows) {
        char error[256] = {0};
        config = ilreco_config_create(n_cols, n_rows, ILRECO_PROF_PWO,
                                      error, sizeof error);
        if (!config) {
            std::fprintf(stderr, "ilreco_config_create failed: %s\n", error);
            std::abort();
        }
        workspace = ilreco_workspace_create(config);
        if (!workspace) {
            std::fprintf(stderr, "ilreco_workspace_create failed\n");
            std::abort();
        }
    }
    ~TestContext() {
        ilreco_workspace_destroy(workspace);
        ilreco_config_destroy(config);
    }
    TestContext(const TestContext&) = delete;
    TestContext& operator=(const TestContext&) = delete;
};

inline TestContext& default_context() {
    static TestContext context(GRID_COLS, GRID_ROWS);
    return context;
}

// Reconstruct one event on the default context. Clusters come back exactly
// as the library returns them: energy-descending. Inputs the library rejects
// (out-of-grid cells etc.) throw — tests that PROBE rejection call
// ilreco_reconstruct directly instead.
inline std::vector<ilreco_cluster> run_event(const std::vector<ilreco_hit>& hits) {
    constexpr int MAX_CLUSTERS = 256;
    TestContext& context = default_context();
    std::vector<ilreco_cluster> clusters(MAX_CLUSTERS);
    const int n_found = ilreco_reconstruct(context.workspace,
                                           hits.data(), (int)hits.size(),
                                           clusters.data(), MAX_CLUSTERS);
    if (n_found < 0)
        throw std::runtime_error("ilreco_reconstruct rejected the test event");
    clusters.resize(std::min(n_found, MAX_CLUSTERS));
    return clusters;
}

inline bool finite(const ilreco_cluster& cluster) {
    return std::isfinite(cluster.e) && std::isfinite(cluster.x) &&
           std::isfinite(cluster.y);
}

// A deterministic little LCG so generated inputs are reproducible
// without <random>.
struct Rng {
    unsigned long long state;
    explicit Rng(unsigned long long seed) : state(seed) {}
    double uniform() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return double((state >> 11) & 0x1FFFFFFFFFFFFFULL) /
               double(0x20000000000000ULL);
    }
    int uniform_int(int lo, int hi) {  // inclusive
        return lo + int(uniform() * (hi - lo + 1));
    }
};

// Synthetic shower: energy-weighted blob around (center_col, center_row)
// with a rough exponential transverse profile — good enough for smoke and
// benchmark inputs (real-data fixtures live in tests/data/).
inline std::vector<ilreco_hit> synth_shower(double energy_gev, int center_col,
                                            int center_row, Rng& rng) {
    std::vector<ilreco_hit> hits;
    for (int dcol = -2; dcol <= 2; ++dcol) {
        for (int drow = -2; drow <= 2; ++drow) {
            const int col = center_col + dcol;
            const int row = center_row + drow;
            if (col < 0 || col >= GRID_COLS || row < 0 || row >= GRID_ROWS)
                continue;
            const double weight =
                std::exp(-1.6 * std::hypot(double(dcol), double(drow)));
            const double energy = energy_gev * weight * (0.8 + 0.4 * rng.uniform());
            if (energy > 0.005) hits.push_back({col, row, energy});
        }
    }
    return hits;
}

}  // namespace ilt

#endif
