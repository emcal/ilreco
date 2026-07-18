// Edge-input behavior. Some of these document CURRENT behavior rather than a
// designed contract (marked [documented]); if a future change alters them,
// that may be intentional — check KNOWN_ISSUES.md before "fixing" the test.
//
// Shared helpers (ilt::run_event, TestContext, synth_shower, Rng) come from
// common/ilreco_test_utils.h; infrastructure map: tests/README.md.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "ilreco_test_utils.h"

using ilt::run_event;

TEST_CASE("corner shower: cluster reconstructed on the grid, energy corrected up",
          "[unit][edge]") {
    // 2x2 blob in the (0,0) corner — half the shower is off the calorimeter
    const auto clusters = run_event({{0, 0, 1.0}, {1, 0, 0.3}, {0, 1, 0.3},
                                     {1, 1, 0.1}});
    REQUIRE(clusters.size() == 1);
    CHECK(clusters[0].x >= 0.0);
    CHECK(clusters[0].y >= 0.0);
    CHECK(clusters[0].x < 1.0);
    CHECK(clusters[0].y < 1.0);
    // leakage correction should push e above the raw sum, bounded by the clamp
    CHECK(clusters[0].e >= 1.7 - 1e-9);
    CHECK(clusters[0].e <= 1.7 / 0.8 + 1e-9);
}

TEST_CASE("full edge row shower keeps position on the boundary", "[unit][edge]") {
    const auto clusters = run_event({{16, 0, 1.0}, {15, 0, 0.4}, {17, 0, 0.4},
                                     {16, 1, 0.25}});
    REQUIRE(clusters.size() == 1);
    CHECK(std::lround(clusters[0].x) == 16);
    CHECK(clusters[0].y >= 0.0);
    CHECK(clusters[0].y < 0.6);
}

TEST_CASE("cluster seed threshold: isolated deposits below ~10 MeV are dropped",
          "[unit][edge][documented]") {
    // DISCOVERED (test-suite finding, see KNOWN_ISSUES): a cluster is only
    // reconstructed if a seed cell exceeds min_seed = 0.01 GeV (scaled up
    // ~7*log(1+E_cluster) for >=3-hit clusters). Isolated hits below that
    // vanish entirely; the same hits DO get counted as members of a cluster
    // with a valid seed. Tunable via ilreco_config_set_seed_threshold().
    CHECK(run_event({{10, 10, 0.0005}}).empty());
    CHECK(run_event({{10, 10, 0.005}}).empty());
    const auto seeded = run_event({{10, 10, 0.011}});
    REQUIRE(seeded.size() == 1);
    CHECK(std::lround(seeded[0].x) == 10);
    // sub-threshold neighbor still contributes to a seeded cluster
    const auto with_neighbor = run_event({{10, 10, 1.0}, {11, 10, 0.0005}});
    REQUIRE(with_neighbor.size() == 1);
    CHECK(with_neighbor[0].size == 2);
}

TEST_CASE("zero-energy hit does not crash", "[unit][edge][documented]") {
    const auto clusters = run_event({{10, 10, 1.0}, {11, 10, 0.0}});
    REQUIRE(!clusters.empty());
    for (const auto& cluster : clusters) CHECK(ilt::finite(cluster));
}

TEST_CASE("duplicate cell address does not crash", "[unit][edge][documented]") {
    // Feeding the same cell twice is undefined by design (upstream should
    // merge); this documents that the library at least survives it.
    const auto clusters = run_event({{10, 10, 0.6}, {10, 10, 0.6}, {11, 10, 0.2}});
    REQUIRE(!clusters.empty());
    for (const auto& cluster : clusters) CHECK(ilt::finite(cluster));
}

TEST_CASE("two maxima in one island separate into two clusters", "[unit][edge]") {
    // peaks at columns 10 and 13 connected by a low valley
    const auto clusters = run_event({{10, 10, 1.0}, {11, 10, 0.25},
                                     {12, 10, 0.25}, {13, 10, 1.0}});
    REQUIRE(clusters.size() >= 2);
    CHECK(std::lround(clusters[0].x + clusters[1].x) == 23);  // 10 + 13
}

TEST_CASE("cluster larger than the per-object hit-list capacity does not crash",
          "[unit][edge][documented]") {
    // one connected blob of 12x12 = 144 cells, more than the 100-entry
    // per-object hit list the library keeps internally; the cluster itself
    // must still come out sane
    std::vector<ilreco_hit> hits;
    for (int col = 5; col < 17; ++col)
        for (int row = 5; row < 17; ++row)
            hits.push_back({col, row, col == 10 && row == 10 ? 2.0 : 0.05});
    const auto clusters = run_event(hits);
    REQUIRE(!clusters.empty());
    for (const auto& cluster : clusters) CHECK(ilt::finite(cluster));
}
