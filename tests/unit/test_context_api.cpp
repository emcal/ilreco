// Context-API tests: 0-based conventions, old-vs-new equivalence (bitwise),
// error handling, and true multithreading (N threads, shared const config,
// one workspace each — results must equal the serial reference exactly).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstring>
#include <thread>
#include <vector>

#include "ilreco_test_utils.h"

namespace {

struct Ctx {
    ilreco_config *cfg = nullptr;
    ilreco_workspace *ws = nullptr;
    Ctx(int nc, int nr, const char *profile = ILRECO_PROF_PWO) {
        char err[256] = {0};
        cfg = ilreco_config_create(nc, nr, profile, err, sizeof err);
        INFO(err);
        REQUIRE(cfg != nullptr);
        ws = ilreco_workspace_create(cfg);
        REQUIRE(ws != nullptr);
    }
    ~Ctx() {
        ilreco_workspace_destroy(ws);
        ilreco_config_destroy(cfg);
    }
};

std::vector<ilreco_cluster> run_new(const Ctx &c, const std::vector<ilreco_hit> &hits) {
    std::vector<ilreco_cluster> out(64);
    const int n = ilreco_reconstruct(c.cfg, c.ws, hits.data(), (int)hits.size(),
                                     out.data(), (int)out.size());
    REQUIRE(n >= 0);
    out.resize(std::min<int>(n, 64));
    return out;
}

struct Init {
    Init() { ilt::init_profile(); }   // legacy context for the equivalence runs
} init_once;

}  // namespace

TEST_CASE("context API: 0-based single-hit round trip on an arbitrary grid",
          "[unit][context]") {
    Ctx c(17, 23);   // deliberately non-square, non-default
    for (auto [col, row] : {std::pair{0, 0}, {16, 22}, {8, 11}, {0, 22}, {16, 0}}) {
        const auto out = run_new(c, {{col, row, 1.0}});
        INFO("cell (" << col << "," << row << ")");
        REQUIRE(out.size() == 1);
        CHECK(std::lround(out[0].x) == col);
        CHECK(std::lround(out[0].y) == row);
    }
}

TEST_CASE("context API rejects invalid input", "[unit][context]") {
    Ctx c(10, 10);
    ilreco_cluster out[4];
    const ilreco_hit bad1{10, 0, 1.0};   // col out of range
    const ilreco_hit bad2{0, -1, 1.0};
    CHECK(ilreco_reconstruct(c.cfg, c.ws, &bad1, 1, out, 4) == -1);
    CHECK(ilreco_reconstruct(c.cfg, c.ws, &bad2, 1, out, 4) == -1);
    CHECK(ilreco_reconstruct(c.cfg, c.ws, nullptr, 0, out, 4) == 0);
}

TEST_CASE("config creation errors are reported, not fatal", "[unit][context]") {
    char err[256] = {0};
    CHECK(ilreco_config_create(10, 10, "/nonexistent/profile.dat",
                               err, sizeof err) == nullptr);
    CHECK(std::strlen(err) > 0);
    CHECK(ilreco_config_create(0, 10, ILRECO_PROF_PWO, err, sizeof err) == nullptr);
}

TEST_CASE("old and new API agree bitwise on the legacy grid", "[unit][context]") {
    Ctx c(_NCOL_, _NROW_);
    ilt::Rng rng(505);
    for (int trial = 0; trial < 200; ++trial) {
        auto ev = ilt::synth_shower(0.5 + 15.0 * rng.uniform(),
                                    rng.uniform_int(1, _NCOL_),
                                    rng.uniform_int(1, _NROW_), rng);
        if (trial % 3 == 0) {   // add a second shower for multi-cluster coverage
            const auto b = ilt::synth_shower(0.5 + 10.0 * rng.uniform(),
                                             rng.uniform_int(1, _NCOL_),
                                             rng.uniform_int(1, _NROW_), rng);
            ev.insert(ev.end(), b.begin(), b.end());
        }
        const auto old_out = ilt::run_event(ev);   // legacy API, energy-sorted
        std::vector<ilreco_hit> hits;
        for (const auto &h : ev) hits.push_back({h.col - 1, h.row - 1, h.e});
        const auto new_out = run_new(c, hits);
        INFO("trial " << trial);
        REQUIRE(old_out.size() == new_out.size());
        for (size_t i = 0; i < old_out.size(); ++i) {
            CHECK(new_out[i].e == old_out[i].e);              // bitwise
            CHECK(new_out[i].x == old_out[i].x - 1.0);
            CHECK(new_out[i].y == old_out[i].y - 1.0);
            CHECK(new_out[i].size == old_out[i].size);
        }
    }
}

TEST_CASE("multithreading: shared config, workspace per thread, exact results",
          "[unit][context][threads]") {
    Ctx serial(_NCOL_, _NROW_);
    constexpr int N_THREADS = 4, N_EVENTS = 200;

    // pre-generate events + serial reference
    std::vector<std::vector<ilreco_hit>> events(N_EVENTS);
    std::vector<std::vector<ilreco_cluster>> ref(N_EVENTS);
    ilt::Rng rng(909);
    for (int i = 0; i < N_EVENTS; ++i) {
        auto ev = ilt::synth_shower(0.5 + 18.0 * rng.uniform(),
                                    rng.uniform_int(1, _NCOL_),
                                    rng.uniform_int(1, _NROW_), rng);
        for (const auto &h : ev) events[i].push_back({h.col - 1, h.row - 1, h.e});
        ref[i] = run_new(serial, events[i]);
    }

    std::vector<std::vector<std::vector<ilreco_cluster>>> results(N_THREADS);
    std::vector<std::thread> threads;
    for (int t = 0; t < N_THREADS; ++t) {
        threads.emplace_back([&, t]() {
            ilreco_workspace *ws = ilreco_workspace_create(serial.cfg);
            results[t].resize(N_EVENTS);
            for (int i = 0; i < N_EVENTS; ++i) {
                ilreco_cluster out[64];
                const int n = ilreco_reconstruct(serial.cfg, ws,
                                                 events[i].data(),
                                                 (int)events[i].size(), out, 64);
                results[t][i].assign(out, out + std::min(n, 64));
            }
            ilreco_workspace_destroy(ws);
        });
    }
    for (auto &th : threads) th.join();

    for (int t = 0; t < N_THREADS; ++t)
        for (int i = 0; i < N_EVENTS; ++i) {
            REQUIRE(results[t][i].size() == ref[i].size());
            for (size_t k = 0; k < ref[i].size(); ++k) {
                REQUIRE(results[t][i][k].e == ref[i][k].e);   // bitwise
                REQUIRE(results[t][i][k].x == ref[i][k].x);
                REQUIRE(results[t][i][k].y == ref[i][k].y);
            }
        }
}

TEST_CASE("context API handles a 100x100 grid (dynamic sizing)", "[unit][context]") {
    Ctx c(100, 100);
    const auto out = run_new(c, {{99, 99, 1.0}, {98, 99, 0.3}, {99, 98, 0.3},
                                 {2, 2, 2.0}, {3, 2, 0.5}});
    REQUIRE(out.size() == 2);
    CHECK(std::lround(out[0].x) == 2);
    CHECK(std::lround(out[1].x) == 99);
}
