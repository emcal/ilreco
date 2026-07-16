// Contract tests for the ilreco public API. These pin down what callers can
// rely on — NOT the internal algorithm. They should survive refactoring
// (0-based indexing / dynamic allocation are planned) as long as the physics
// contract holds; only the helper in tests/common encodes the calling
// convention and can be updated alongside such changes.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

#include "ilreco_test_utils.h"

using ilt::Hit;
using ilt::run_event;

namespace {
struct Init {
    Init() { ilt::init_profile(); }
} init_once;
}  // namespace

TEST_CASE("empty event yields zero clusters", "[unit][contract]") {
    REQUIRE(run_event({}).empty());
}

TEST_CASE("single hit reconstructs at its own cell", "[unit][contract]") {
    // corners, edge midpoints and center of the compiled-in grid
    const int C = _NCOL_, R = _NROW_;
    const int probes[][2] = {{1, 1},     {1, R},      {C, 1},    {C, R},
                             {C / 2, 1}, {1, R / 2},  {C, R / 2}, {C / 2, R},
                             {C / 2, R / 2}};
    for (const auto& p : probes) {
        const auto out = run_event({{p[0], p[1], 1.0}});
        INFO("probe cell (" << p[0] << "," << p[1] << ")");
        REQUIRE(out.size() == 1);
        CHECK(std::lround(out[0].x) == p[0]);
        CHECK(std::lround(out[0].y) == p[1]);
        // profile containment correction is clamped to [0.8, 1.0]:
        CHECK(out[0].e >= 1.0 - 1e-9);
        CHECK(out[0].e <= 1.0 / 0.8 + 1e-9);
    }
}

TEST_CASE("two well-separated showers give two clusters at the right places",
          "[unit][contract]") {
    std::vector<Hit> ev = {{5, 5, 2.0}, {4, 5, 0.3}, {6, 5, 0.3}, {5, 4, 0.3},
                           {5, 6, 0.3},
                           {20, 20, 1.0}, {19, 20, 0.15}, {21, 20, 0.15}};
    const auto out = run_event(ev);
    REQUIRE(out.size() == 2);
    CHECK(std::lround(out[0].x) == 5);
    CHECK(std::lround(out[0].y) == 5);
    CHECK(std::lround(out[1].x) == 20);
    CHECK(std::lround(out[1].y) == 20);
    CHECK(out[0].e > out[1].e);
}

TEST_CASE("asymmetric adjacent pair: one cluster between cells, nearer the hot one",
          "[unit][contract]") {
    const auto out = run_event({{10, 10, 0.7}, {11, 10, 0.3}});
    REQUIRE(out.size() == 1);
    CHECK(out[0].x > 10.0);
    CHECK(out[0].x < 10.5);
    CHECK(std::lround(out[0].y) == 10);
}

TEST_CASE("input hit order does not change the result", "[unit][contract]") {
    std::vector<Hit> ev = {{7, 7, 1.5}, {8, 7, 0.4}, {7, 8, 0.3}, {6, 7, 0.2},
                           {7, 6, 0.1}, {25, 12, 0.8}, {26, 12, 0.2}};
    const auto a = run_event(ev);
    std::reverse(ev.begin(), ev.end());
    const auto b = run_event(ev);
    std::swap(ev[0], ev[3]);
    const auto c = run_event(ev);
    REQUIRE(a.size() == b.size());
    REQUIRE(a.size() == c.size());
    for (size_t i = 0; i < a.size(); ++i) {
        CHECK_THAT(a[i].e, Catch::Matchers::WithinRel(b[i].e, 1e-12));
        CHECK_THAT(a[i].x, Catch::Matchers::WithinAbs(b[i].x, 1e-12));
        CHECK_THAT(a[i].y, Catch::Matchers::WithinAbs(c[i].y, 1e-12));
        CHECK_THAT(a[i].e, Catch::Matchers::WithinRel(c[i].e, 1e-12));
    }
}

TEST_CASE("energy scaling: doubling all inputs doubles cluster energies",
          "[unit][contract]") {
    std::vector<Hit> ev = {{15, 15, 1.0}, {16, 15, 0.25}, {15, 16, 0.2},
                           {14, 15, 0.15}, {15, 14, 0.1}};
    const auto a = run_event(ev);
    for (auto& h : ev) h.e *= 2.0;
    const auto b = run_event(ev);
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        CHECK_THAT(b[i].e / a[i].e, Catch::Matchers::WithinAbs(2.0, 0.02));
        CHECK_THAT(b[i].x, Catch::Matchers::WithinAbs(a[i].x, 0.05));
        CHECK_THAT(b[i].y, Catch::Matchers::WithinAbs(a[i].y, 0.05));
    }
}

TEST_CASE("cluster energy stays within the containment-correction envelope",
          "[unit][contract]") {
    ilt::Rng rng(20260716);
    for (int trial = 0; trial < 50; ++trial) {
        const int c = rng.uniform_int(3, _NCOL_ - 2);
        const int r = rng.uniform_int(3, _NROW_ - 2);
        const auto ev = ilt::synth_shower(1.0 + 9.0 * rng.uniform(), c, r, rng);
        double raw = 0;
        for (const auto& h : ev) raw += h.e;
        double reco = 0;
        for (const auto& cl : run_event(ev)) reco += cl.e;
        INFO("shower at (" << c << "," << r << "), raw " << raw);
        CHECK(reco >= raw - 1e-9);              // correction only adds
        CHECK(reco <= raw / 0.8 + 1e-9);        // clamp floor 0.8
    }
}

TEST_CASE("outputs are finite and inside the grid for physical inputs",
          "[unit][contract]") {
    ilt::Rng rng(424242);
    for (int trial = 0; trial < 100; ++trial) {
        const auto ev = ilt::synth_shower(0.1 + 19.9 * rng.uniform(),
                                          rng.uniform_int(1, _NCOL_),
                                          rng.uniform_int(1, _NROW_), rng);
        if (ev.empty()) continue;
        for (const auto& cl : run_event(ev)) {
            CHECK(ilt::finite(cl));
            CHECK(cl.x > 0.0);
            CHECK(cl.x < _NCOL_ + 1.0);
            CHECK(cl.y > 0.0);
            CHECK(cl.y < _NROW_ + 1.0);
            CHECK(cl.e > 0.0);
        }
    }
}

TEST_CASE("profile re-initialization does not change results",
          "[unit][contract]") {
    const std::vector<Hit> ev = {{12, 12, 1.0}, {13, 12, 0.4}};
    const auto a = run_event(ev);
    read_profile_data(ILRECO_PROF_PWO);
    const auto b = run_event(ev);
    REQUIRE(a.size() == b.size());
    CHECK_THAT(a[0].e, Catch::Matchers::WithinRel(b[0].e, 1e-12));
    CHECK_THAT(a[0].x, Catch::Matchers::WithinAbs(b[0].x, 1e-12));
}
