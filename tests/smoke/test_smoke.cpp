// Smoke sweeps: large volumes of generated inputs, asserting only the coarse
// health contract — no crash, finite outputs, positions on the grid. This is
// the first net that catches memory corruption after a refactor.

#include <catch2/catch_test_macros.hpp>

#include "ilreco_test_utils.h"

using ilt::Hit;
using ilt::run_event;

namespace {
struct Init {
    Init() { ilt::init_profile(); }
} init_once;

void check_sane(const std::vector<ilt::Cluster>& out) {
    for (const auto& c : out) {
        REQUIRE(ilt::finite(c));
        REQUIRE(c.e > 0.0);
        REQUIRE(c.x > 0.0);
        REQUIRE(c.x < _NCOL_ + 1.0);
        REQUIRE(c.y > 0.0);
        REQUIRE(c.y < _NROW_ + 1.0);
        REQUIRE(c.size >= 1);
    }
}
}  // namespace

TEST_CASE("1000 random single showers", "[smoke]") {
    ilt::Rng rng(1);
    for (int i = 0; i < 1000; ++i) {
        const auto ev = ilt::synth_shower(0.05 + 20.0 * rng.uniform(),
                                          rng.uniform_int(1, _NCOL_),
                                          rng.uniform_int(1, _NROW_), rng);
        if (ev.empty()) continue;
        check_sane(run_event(ev));
    }
}

TEST_CASE("500 random overlapping shower pairs", "[smoke]") {
    ilt::Rng rng(2);
    for (int i = 0; i < 500; ++i) {
        auto ev = ilt::synth_shower(1.0 + 10.0 * rng.uniform(),
                                    rng.uniform_int(2, _NCOL_ - 1),
                                    rng.uniform_int(2, _NROW_ - 1), rng);
        const auto b = ilt::synth_shower(1.0 + 10.0 * rng.uniform(),
                                         rng.uniform_int(2, _NCOL_ - 1),
                                         rng.uniform_int(2, _NROW_ - 1), rng);
        ev.insert(ev.end(), b.begin(), b.end());
        check_sane(run_event(ev));
    }
}

TEST_CASE("random sparse noise events", "[smoke]") {
    ilt::Rng rng(3);
    for (int i = 0; i < 500; ++i) {
        std::vector<Hit> ev;
        const int n = rng.uniform_int(1, 40);
        for (int k = 0; k < n; ++k)
            ev.push_back({rng.uniform_int(1, _NCOL_), rng.uniform_int(1, _NROW_),
                          0.005 + 0.5 * rng.uniform()});
        check_sane(run_event(ev));
    }
}

TEST_CASE("whole grid lit at once", "[smoke]") {
    std::vector<Hit> ev;
    for (int c = 1; c <= _NCOL_; ++c)
        for (int r = 1; r <= _NROW_; ++r)
            ev.push_back({c, r, 0.01 + 0.001 * ((c * 7 + r * 13) % 17)});
    check_sane(run_event(ev));
}

TEST_CASE("back-to-back events do not leak state", "[smoke]") {
    // identical event before and after a monster event must give identical output
    const std::vector<Hit> probe = {{9, 9, 1.0}, {10, 9, 0.4}};
    const auto before = run_event(probe);
    std::vector<Hit> monster;
    for (int c = 1; c <= _NCOL_; ++c)
        for (int r = 1; r <= _NROW_; ++r) monster.push_back({c, r, 0.02});
    run_event(monster);
    const auto after = run_event(probe);
    REQUIRE(before.size() == after.size());
    for (size_t i = 0; i < before.size(); ++i) {
        REQUIRE(before[i].e == after[i].e);
        REQUIRE(before[i].x == after[i].x);
        REQUIRE(before[i].y == after[i].y);
    }
}
