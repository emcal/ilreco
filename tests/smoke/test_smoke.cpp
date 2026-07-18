// Smoke sweeps: large volumes of generated inputs, asserting only the coarse
// health contract — no crash, finite outputs, positions on the grid. This is
// the first net that catches memory corruption after a refactor.
//
// Shared helpers (ilt::run_event, TestContext, synth_shower, Rng) come from
// common/ilreco_test_utils.h; infrastructure map: tests/README.md.

#include <catch2/catch_test_macros.hpp>

#include <vector>

#include "ilreco_test_utils.h"

using ilt::run_event;
using ilt::GRID_COLS;
using ilt::GRID_ROWS;

namespace {

void check_sane(const std::vector<ilreco_cluster>& clusters) {
    for (const auto& cluster : clusters) {
        REQUIRE(ilt::finite(cluster));
        REQUIRE(cluster.e > 0.0);
        REQUIRE(cluster.x > -1.0);
        REQUIRE(cluster.x < GRID_COLS);
        REQUIRE(cluster.y > -1.0);
        REQUIRE(cluster.y < GRID_ROWS);
        REQUIRE(cluster.size >= 1);
    }
}

}  // namespace

TEST_CASE("1000 random single showers", "[smoke]") {
    ilt::Rng rng(1);
    for (int i = 0; i < 1000; ++i) {
        const auto hits = ilt::synth_shower(0.05 + 20.0 * rng.uniform(),
                                            rng.uniform_int(0, GRID_COLS - 1),
                                            rng.uniform_int(0, GRID_ROWS - 1), rng);
        if (hits.empty()) continue;
        check_sane(run_event(hits));
    }
}

TEST_CASE("500 random overlapping shower pairs", "[smoke]") {
    ilt::Rng rng(2);
    for (int i = 0; i < 500; ++i) {
        auto hits = ilt::synth_shower(1.0 + 10.0 * rng.uniform(),
                                      rng.uniform_int(1, GRID_COLS - 2),
                                      rng.uniform_int(1, GRID_ROWS - 2), rng);
        const auto second = ilt::synth_shower(1.0 + 10.0 * rng.uniform(),
                                              rng.uniform_int(1, GRID_COLS - 2),
                                              rng.uniform_int(1, GRID_ROWS - 2), rng);
        hits.insert(hits.end(), second.begin(), second.end());
        check_sane(run_event(hits));
    }
}

TEST_CASE("random sparse noise events", "[smoke]") {
    ilt::Rng rng(3);
    for (int i = 0; i < 500; ++i) {
        std::vector<ilreco_hit> hits;
        const int n_hits = rng.uniform_int(1, 40);
        for (int k = 0; k < n_hits; ++k)
            hits.push_back({rng.uniform_int(0, GRID_COLS - 1),
                            rng.uniform_int(0, GRID_ROWS - 1),
                            0.005 + 0.5 * rng.uniform()});
        check_sane(run_event(hits));
    }
}

TEST_CASE("whole grid lit at once", "[smoke]") {
    std::vector<ilreco_hit> hits;
    for (int col = 0; col < GRID_COLS; ++col)
        for (int row = 0; row < GRID_ROWS; ++row)
            hits.push_back({col, row, 0.01 + 0.001 * ((col * 7 + row * 13) % 17)});
    check_sane(run_event(hits));
}

TEST_CASE("workspace reuse leaves no state behind between events", "[smoke]") {
    // The same workspace is reused for every event, so its buffers still hold
    // whatever the previous event wrote there. Reconstruction must depend only
    // on the current input: reconstruct a small two-hit event, then an event
    // that lights every cell of the grid (overwriting every buffer with
    // non-trivial data), then the same two-hit event again. If any leftover
    // buffer content influenced the result, the two small-event results would
    // differ; they must be bitwise identical.
    const std::vector<ilreco_hit> small_event = {{9, 9, 1.0}, {10, 9, 0.4}};
    const auto before = run_event(small_event);

    std::vector<ilreco_hit> grid_filling_event;
    for (int col = 0; col < GRID_COLS; ++col)
        for (int row = 0; row < GRID_ROWS; ++row)
            grid_filling_event.push_back({col, row, 0.02});
    run_event(grid_filling_event);

    const auto after = run_event(small_event);
    REQUIRE(before.size() == after.size());
    for (size_t i = 0; i < before.size(); ++i) {
        REQUIRE(before[i].e == after[i].e);
        REQUIRE(before[i].x == after[i].x);
        REQUIRE(before[i].y == after[i].y);
    }
}
