// The C++ RAII wrapper (ilreco.hpp) must behave exactly like the C API it
// wraps: same results, exceptions instead of error codes, move-only
// ownership that never double-frees.
//
// Shared helpers (ilt::run_event, TestContext, synth_shower, Rng) come from
// common/ilreco_test_utils.h; infrastructure map: tests/README.md.

#include <catch2/catch_test_macros.hpp>

#include <ilreco.hpp>

#include <utility>
#include <vector>

#include "ilreco_test_utils.h"

TEST_CASE("C++ wrapper reproduces the C API results", "[unit][cpp]") {
    ilreco::Config config(ilt::GRID_COLS, ilt::GRID_ROWS, ILRECO_PROF_PWO);
    ilreco::Workspace workspace(config);

    const std::vector<ilreco_hit> hits = {{12, 14, 1.85}, {13, 14, 0.42},
                                          {12, 15, 0.31}};
    const auto from_wrapper = workspace.reconstruct(hits);
    const auto from_c_api = ilt::run_event(hits);

    REQUIRE(from_wrapper.size() == from_c_api.size());
    for (size_t i = 0; i < from_wrapper.size(); ++i) {
        CHECK(from_wrapper[i].e == from_c_api[i].e);      // bitwise
        CHECK(from_wrapper[i].x == from_c_api[i].x);
        CHECK(from_wrapper[i].y == from_c_api[i].y);
    }
}

TEST_CASE("C++ wrapper reports errors as exceptions", "[unit][cpp]") {
    CHECK_THROWS_AS(ilreco::Config(10, 10, "/nonexistent/profile.dat"),
                    std::runtime_error);

    ilreco::Config config(10, 10, ILRECO_PROF_PWO);
    ilreco::Workspace workspace(config);
    CHECK_THROWS_AS(workspace.reconstruct({{10, 0, 1.0}}),
                    std::invalid_argument);
}

TEST_CASE("C++ wrapper: cell mask and move semantics", "[unit][cpp]") {
    ilreco::Config config(9, 9, ILRECO_PROF_PWO);
    std::vector<unsigned char> all_cells_exist(9 * 9, 1);
    config.set_cell_mask(all_cells_exist);

    ilreco::Workspace workspace(config);
    // moving transfers ownership; the moved-to object keeps working
    ilreco::Workspace moved_workspace(std::move(workspace));
    const auto clusters = moved_workspace.reconstruct({{4, 4, 1.0}});
    REQUIRE(clusters.size() == 1);
    CHECK(clusters[0].type == 0);   // interior, per mask-derived labels
}
