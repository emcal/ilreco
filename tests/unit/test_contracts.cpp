// Contract tests for the ilreco public API. These pin down what callers can
// rely on — NOT the internal algorithm. Cells are 0-based (col, row), cluster
// positions in 0-based cell units, exactly as in ilreco.h.
//
// Shared helpers (ilt::run_event, TestContext, synth_shower, Rng) come from
// common/ilreco_test_utils.h; infrastructure map: tests/README.md.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>

#include "ilreco_test_utils.h"

using ilt::run_event;
using ilt::GRID_COLS;
using ilt::GRID_ROWS;

TEST_CASE("empty event yields zero clusters", "[unit][contract]") {
    REQUIRE(run_event({}).empty());
}

TEST_CASE("single hit reconstructs at its own cell", "[unit][contract]") {
    // corners, edge midpoints and center of the grid
    const int last_col = GRID_COLS - 1, last_row = GRID_ROWS - 1;
    const int mid_col = GRID_COLS / 2, mid_row = GRID_ROWS / 2;
    const int probe_cells[][2] = {
        {0, 0},        {0, last_row},       {last_col, 0}, {last_col, last_row},
        {mid_col, 0},  {0, mid_row},        {last_col, mid_row},
        {mid_col, last_row},                {mid_col, mid_row}};
    for (const auto& cell : probe_cells) {
        const auto clusters = run_event({{cell[0], cell[1], 1.0}});
        INFO("probe cell (" << cell[0] << "," << cell[1] << ")");
        REQUIRE(clusters.size() == 1);
        CHECK(std::lround(clusters[0].x) == cell[0]);
        CHECK(std::lround(clusters[0].y) == cell[1]);
        // profile containment correction is clamped to [0.8, 1.0]:
        CHECK(clusters[0].e >= 1.0 - 1e-9);
        CHECK(clusters[0].e <= 1.0 / 0.8 + 1e-9);
    }
}

TEST_CASE("two well-separated showers give two clusters at the right places",
          "[unit][contract]") {
    const std::vector<ilreco_hit> hits = {
        {5, 5, 2.0},   {4, 5, 0.3},    {6, 5, 0.3},    {5, 4, 0.3}, {5, 6, 0.3},
        {20, 20, 1.0}, {19, 20, 0.15}, {21, 20, 0.15}};
    const auto clusters = run_event(hits);
    REQUIRE(clusters.size() == 2);
    CHECK(std::lround(clusters[0].x) == 5);
    CHECK(std::lround(clusters[0].y) == 5);
    CHECK(std::lround(clusters[1].x) == 20);
    CHECK(std::lround(clusters[1].y) == 20);
    CHECK(clusters[0].e > clusters[1].e);   // energy-descending output
}

TEST_CASE("asymmetric adjacent pair: one cluster between cells, nearer the hot one",
          "[unit][contract]") {
    const auto clusters = run_event({{10, 10, 0.7}, {11, 10, 0.3}});
    REQUIRE(clusters.size() == 1);
    CHECK(clusters[0].x > 10.0);
    CHECK(clusters[0].x < 10.5);
    CHECK(std::lround(clusters[0].y) == 10);
}

TEST_CASE("input hit order does not change the result", "[unit][contract]") {
    std::vector<ilreco_hit> hits = {{7, 7, 1.5},  {8, 7, 0.4},   {7, 8, 0.3},
                                    {6, 7, 0.2},  {7, 6, 0.1},   {25, 12, 0.8},
                                    {26, 12, 0.2}};
    const auto original = run_event(hits);
    std::reverse(hits.begin(), hits.end());
    const auto reversed = run_event(hits);
    std::swap(hits[0], hits[3]);
    const auto shuffled = run_event(hits);
    REQUIRE(original.size() == reversed.size());
    REQUIRE(original.size() == shuffled.size());
    for (size_t i = 0; i < original.size(); ++i) {
        CHECK_THAT(original[i].e, Catch::Matchers::WithinRel(reversed[i].e, 1e-12));
        CHECK_THAT(original[i].x, Catch::Matchers::WithinAbs(reversed[i].x, 1e-12));
        CHECK_THAT(original[i].y, Catch::Matchers::WithinAbs(shuffled[i].y, 1e-12));
        CHECK_THAT(original[i].e, Catch::Matchers::WithinRel(shuffled[i].e, 1e-12));
    }
}

TEST_CASE("energy scaling: doubling all inputs doubles cluster energies",
          "[unit][contract]") {
    std::vector<ilreco_hit> hits = {{15, 15, 1.0},  {16, 15, 0.25}, {15, 16, 0.2},
                                    {14, 15, 0.15}, {15, 14, 0.1}};
    const auto nominal = run_event(hits);
    for (auto& hit : hits) hit.e *= 2.0;
    const auto doubled = run_event(hits);
    REQUIRE(nominal.size() == doubled.size());
    for (size_t i = 0; i < nominal.size(); ++i) {
        CHECK_THAT(doubled[i].e / nominal[i].e, Catch::Matchers::WithinAbs(2.0, 0.02));
        CHECK_THAT(doubled[i].x, Catch::Matchers::WithinAbs(nominal[i].x, 0.05));
        CHECK_THAT(doubled[i].y, Catch::Matchers::WithinAbs(nominal[i].y, 0.05));
    }
}

TEST_CASE("cluster energy stays within the containment-correction envelope",
          "[unit][contract]") {
    ilt::Rng rng(20260716);
    for (int trial = 0; trial < 50; ++trial) {
        // keep the shower fully inside the grid (synth_shower spans +-2 cells)
        const int center_col = rng.uniform_int(2, GRID_COLS - 3);
        const int center_row = rng.uniform_int(2, GRID_ROWS - 3);
        const auto hits = ilt::synth_shower(1.0 + 9.0 * rng.uniform(),
                                            center_col, center_row, rng);
        double raw_sum = 0;
        for (const auto& hit : hits) raw_sum += hit.e;
        double reco_sum = 0;
        for (const auto& cluster : run_event(hits)) reco_sum += cluster.e;
        INFO("shower at (" << center_col << "," << center_row << "), raw " << raw_sum);
        CHECK(reco_sum >= raw_sum - 1e-9);         // correction only adds
        CHECK(reco_sum <= raw_sum / 0.8 + 1e-9);   // clamp floor 0.8
    }
}

TEST_CASE("outputs are finite and inside the grid for physical inputs",
          "[unit][contract]") {
    ilt::Rng rng(424242);
    for (int trial = 0; trial < 100; ++trial) {
        const auto hits = ilt::synth_shower(0.1 + 19.9 * rng.uniform(),
                                            rng.uniform_int(0, GRID_COLS - 1),
                                            rng.uniform_int(0, GRID_ROWS - 1), rng);
        if (hits.empty()) continue;
        for (const auto& cluster : run_event(hits)) {
            CHECK(ilt::finite(cluster));
            CHECK(cluster.x > -1.0);
            CHECK(cluster.x < GRID_COLS);
            CHECK(cluster.y > -1.0);
            CHECK(cluster.y < GRID_ROWS);
            CHECK(cluster.e > 0.0);
        }
    }
}

TEST_CASE("independent contexts from the same profile give identical results",
          "[unit][contract]") {
    const std::vector<ilreco_hit> hits = {{12, 12, 1.0}, {13, 12, 0.4}};
    const auto from_shared = run_event(hits);

    ilt::TestContext fresh(GRID_COLS, GRID_ROWS);
    ilreco_cluster clusters[8];
    const int n_found = ilreco_reconstruct(fresh.workspace,
                                           hits.data(), (int)hits.size(),
                                           clusters, 8);
    REQUIRE(n_found == (int)from_shared.size());
    CHECK(clusters[0].e == from_shared[0].e);   // bitwise: same config, same math
    CHECK(clusters[0].x == from_shared[0].x);
    CHECK(clusters[0].y == from_shared[0].y);
}
