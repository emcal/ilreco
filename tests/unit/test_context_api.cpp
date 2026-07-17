// Context-API tests: 0-based conventions, error handling, dynamic grid
// sizing, and true multithreading (N threads, one shared const config, one
// workspace per thread — results must equal the serial reference exactly).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstring>
#include <thread>
#include <vector>

#include "ilreco_test_utils.h"

namespace {

std::vector<ilreco_cluster> run_on(const ilt::TestContext& context,
                                   const std::vector<ilreco_hit>& hits) {
    std::vector<ilreco_cluster> clusters(64);
    const int n_found = ilreco_reconstruct(context.config, context.workspace,
                                           hits.data(), (int)hits.size(),
                                           clusters.data(), (int)clusters.size());
    REQUIRE(n_found >= 0);
    clusters.resize(std::min<int>(n_found, 64));
    return clusters;
}

}  // namespace

TEST_CASE("context API: 0-based single-hit round trip on an arbitrary grid",
          "[unit][context]") {
    ilt::TestContext context(17, 23);   // deliberately non-square, non-default
    for (auto [col, row] : {std::pair{0, 0}, {16, 22}, {8, 11}, {0, 22}, {16, 0}}) {
        const auto clusters = run_on(context, {{col, row, 1.0}});
        INFO("cell (" << col << "," << row << ")");
        REQUIRE(clusters.size() == 1);
        CHECK(std::lround(clusters[0].x) == col);
        CHECK(std::lround(clusters[0].y) == row);
    }
}

TEST_CASE("context API rejects invalid input", "[unit][context]") {
    ilt::TestContext context(10, 10);
    ilreco_cluster clusters[4];
    const ilreco_hit col_out_of_range{10, 0, 1.0};
    const ilreco_hit row_negative{0, -1, 1.0};
    CHECK(ilreco_reconstruct(context.config, context.workspace,
                             &col_out_of_range, 1, clusters, 4) == -1);
    CHECK(ilreco_reconstruct(context.config, context.workspace,
                             &row_negative, 1, clusters, 4) == -1);
    CHECK(ilreco_reconstruct(context.config, context.workspace,
                             nullptr, 0, clusters, 4) == 0);
}

TEST_CASE("config creation errors are reported, not fatal", "[unit][context]") {
    char error[256] = {0};
    CHECK(ilreco_config_create(10, 10, "/nonexistent/profile.dat",
                               error, sizeof error) == nullptr);
    CHECK(std::strlen(error) > 0);
    CHECK(ilreco_config_create(0, 10, ILRECO_PROF_PWO, error, sizeof error) == nullptr);
    CHECK(ilreco_config_create(10, 10, nullptr, error, sizeof error) == nullptr);
}

TEST_CASE("multithreading: shared config, workspace per thread, exact results",
          "[unit][context][threads]") {
    ilt::TestContext serial(ilt::GRID_COLS, ilt::GRID_ROWS);
    constexpr int N_THREADS = 4, N_EVENTS = 200;

    // pre-generate events + serial reference
    std::vector<std::vector<ilreco_hit>> events(N_EVENTS);
    std::vector<std::vector<ilreco_cluster>> reference(N_EVENTS);
    ilt::Rng rng(909);
    for (int i = 0; i < N_EVENTS; ++i) {
        events[i] = ilt::synth_shower(0.5 + 18.0 * rng.uniform(),
                                      rng.uniform_int(0, ilt::GRID_COLS - 1),
                                      rng.uniform_int(0, ilt::GRID_ROWS - 1), rng);
        reference[i] = run_on(serial, events[i]);
    }

    std::vector<std::vector<std::vector<ilreco_cluster>>> per_thread(N_THREADS);
    std::vector<std::thread> threads;
    for (int t = 0; t < N_THREADS; ++t) {
        threads.emplace_back([&, t]() {
            ilreco_workspace* workspace = ilreco_workspace_create(serial.config);
            per_thread[t].resize(N_EVENTS);
            for (int i = 0; i < N_EVENTS; ++i) {
                ilreco_cluster clusters[64];
                const int n_found =
                    ilreco_reconstruct(serial.config, workspace, events[i].data(),
                                       (int)events[i].size(), clusters, 64);
                per_thread[t][i].assign(clusters, clusters + std::min(n_found, 64));
            }
            ilreco_workspace_destroy(workspace);
        });
    }
    for (auto& thread : threads) thread.join();

    for (int t = 0; t < N_THREADS; ++t)
        for (int i = 0; i < N_EVENTS; ++i) {
            REQUIRE(per_thread[t][i].size() == reference[i].size());
            for (size_t k = 0; k < reference[i].size(); ++k) {
                REQUIRE(per_thread[t][i][k].e == reference[i][k].e);   // bitwise
                REQUIRE(per_thread[t][i][k].x == reference[i][k].x);
                REQUIRE(per_thread[t][i][k].y == reference[i][k].y);
            }
        }
}

TEST_CASE("context API handles a 100x100 grid (dynamic sizing)", "[unit][context]") {
    ilt::TestContext context(100, 100);
    const auto clusters = run_on(context, {{99, 99, 1.0}, {98, 99, 0.3},
                                           {99, 98, 0.3}, {2, 2, 2.0},
                                           {3, 2, 0.5}});
    REQUIRE(clusters.size() == 2);
    CHECK(std::lround(clusters[0].x) == 2);    // energy-descending: 2.5 GeV first
    CHECK(std::lround(clusters[1].x) == 99);
}
