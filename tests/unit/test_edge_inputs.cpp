// Edge-input behavior. Some of these document CURRENT behavior rather than a
// designed contract (marked [documented]); if a future change alters them,
// that may be intentional — check KNOWN_ISSUES.md before "fixing" the test.

#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "ilreco_test_utils.h"

using ilt::Hit;
using ilt::run_event;

namespace {
struct Init {
    Init() { ilt::init_profile(); }
} init_once;
}  // namespace

TEST_CASE("corner shower: cluster reconstructed on the grid, energy corrected up",
          "[unit][edge]") {
    // 2x2 blob in the (1,1) corner — half the shower is off the calorimeter
    const auto out = run_event({{1, 1, 1.0}, {2, 1, 0.3}, {1, 2, 0.3}, {2, 2, 0.1}});
    REQUIRE(out.size() == 1);
    CHECK(out[0].x >= 1.0);
    CHECK(out[0].y >= 1.0);
    CHECK(out[0].x < 2.0);
    CHECK(out[0].y < 2.0);
    // leakage correction should push e above the raw sum, bounded by the clamp
    CHECK(out[0].e >= 1.7 - 1e-9);
    CHECK(out[0].e <= 1.7 / 0.8 + 1e-9);
}

TEST_CASE("full edge row shower keeps position on the boundary", "[unit][edge]") {
    const auto out = run_event({{17, 1, 1.0}, {16, 1, 0.4}, {18, 1, 0.4},
                                {17, 2, 0.25}});
    REQUIRE(out.size() == 1);
    CHECK(std::lround(out[0].x) == 17);
    CHECK(out[0].y >= 1.0);
    CHECK(out[0].y < 1.6);
}

TEST_CASE("cluster seed threshold: isolated deposits below ~10 MeV are dropped",
          "[unit][edge][documented]") {
    // DISCOVERED (test-suite finding, see KNOWN_ISSUES): process_cluster requires a
    // seed cell above minpk = 0.01 GeV (scaled up ~7*log(1+E_cluster) for >=3-hit
    // clusters). Isolated hits below that vanish entirely; the same hits DO get
    // counted as members of a cluster with a valid seed. MIN_COUNTER_ENERGY (1 MeV)
    // in the header is unrelated and unused by the library.
    CHECK(run_event({{10, 10, 0.0005}}).empty());
    CHECK(run_event({{10, 10, 0.005}}).empty());
    const auto seeded = run_event({{10, 10, 0.011}});
    REQUIRE(seeded.size() == 1);
    CHECK(std::lround(seeded[0].x) == 10);
    // sub-threshold neighbor still contributes to a seeded cluster
    const auto pair = run_event({{10, 10, 1.0}, {11, 10, 0.0005}});
    REQUIRE(pair.size() == 1);
    CHECK(pair[0].size == 2);
}

TEST_CASE("zero-energy hit does not crash", "[unit][edge][documented]") {
    const auto out = run_event({{10, 10, 1.0}, {11, 10, 0.0}});
    REQUIRE(!out.empty());
    for (const auto& c : out) CHECK(ilt::finite(c));
}

TEST_CASE("duplicate cell address does not crash", "[unit][edge][documented]") {
    // Feeding the same cell twice is undefined by design (upstream should merge);
    // this documents that the library at least survives it. See KNOWN_ISSUES #5.
    const auto out = run_event({{10, 10, 0.6}, {10, 10, 0.6}, {11, 10, 0.2}});
    REQUIRE(!out.empty());
    for (const auto& c : out) CHECK(ilt::finite(c));
}

TEST_CASE("two maxima in one island separate into two clusters", "[unit][edge]") {
    // peaks at columns 10 and 13 connected by a low valley
    const auto out = run_event({{10, 10, 1.0}, {11, 10, 0.25}, {12, 10, 0.25},
                                {13, 10, 1.0}});
    REQUIRE(out.size() >= 2);
    CHECK(std::lround(out[0].x + out[1].x) == 23);  // 10 + 13
}

TEST_CASE("cluster larger than _MAXLEN_ hits does not crash",
          "[unit][edge][documented]") {
    // one connected blob of 12x12 = 144 cells > _MAXLEN_ = 100 (element/elfract
    // arrays cap at _MAXLEN_; the cluster itself must still come out sane)
    std::vector<Hit> ev;
    for (int c = 5; c < 17; ++c)
        for (int r = 5; r < 17; ++r)
            ev.push_back({c, r, c == 10 && r == 10 ? 2.0 : 0.05});
    const auto out = run_event(ev);
    REQUIRE(!out.empty());
    for (const auto& c : out) CHECK(ilt::finite(c));
}
