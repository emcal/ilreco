// Cell-existence mask: arbitrary shapes (circular calorimeter, asymmetric
// hole — the EIC-B0-like case). Verifies input validation, mask-derived type
// labels, the physics effect (honest containment correction at real
// boundaries), and equivalence when the mask marks every cell as existing.
//
// Shared helpers (ilt::run_event, TestContext, synth_shower, Rng) come from
// common/ilreco_test_utils.h; infrastructure map: tests/README.md.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "ilreco_test_utils.h"

namespace {

constexpr int GRID_SIZE = 9;

// Circular-ish 9x9 detector with an off-center 1x2 beam hole at
// (col 5, row 4) and (col 5, row 5).
std::vector<unsigned char> circular_detector_mask() {
    std::vector<unsigned char> mask(GRID_SIZE * GRID_SIZE, 0);
    for (int row = 0; row < GRID_SIZE; ++row)
        for (int col = 0; col < GRID_SIZE; ++col) {
            const int dist2 = (col - 4) * (col - 4) + (row - 4) * (row - 4);
            if (dist2 <= 18) mask[row * GRID_SIZE + col] = 1;   // inside the circle
        }
    mask[4 * GRID_SIZE + 5] = 0;   // beam hole, cell (col 5, row 4)
    mask[5 * GRID_SIZE + 5] = 0;   // beam hole, cell (col 5, row 5)
    return mask;
}

// A 9x9 context, optionally with a cell mask applied at creation.
struct MaskedContext {
    ilreco_config* config = nullptr;
    ilreco_workspace* workspace = nullptr;

    explicit MaskedContext(const unsigned char* mask) {
        char error[256] = {0};
        config = ilreco_config_create(GRID_SIZE, GRID_SIZE, ILRECO_PROF_PWO,
                                      error, sizeof error);
        REQUIRE(config != nullptr);
        if (mask) REQUIRE(ilreco_config_set_cell_mask(config, mask) == 0);
        workspace = ilreco_workspace_create(config);
        REQUIRE(workspace != nullptr);
    }
    ~MaskedContext() {
        ilreco_workspace_destroy(workspace);
        ilreco_config_destroy(config);
    }

    std::vector<ilreco_cluster> run(const std::vector<ilreco_hit>& hits) {
        std::vector<ilreco_cluster> clusters(16);
        const int n_found = ilreco_reconstruct(workspace, hits.data(),
                                               (int)hits.size(), clusters.data(), 16);
        REQUIRE(n_found >= 0);
        clusters.resize(std::min(n_found, 16));
        return clusters;
    }
};

}  // namespace

TEST_CASE("mask: hits on missing cells are rejected", "[unit][mask]") {
    const auto mask = circular_detector_mask();
    MaskedContext detector(mask.data());
    ilreco_cluster clusters[4];
    const ilreco_hit on_hole_cell{5, 4, 1.0};
    const ilreco_hit outside_circle{0, 0, 1.0};
    CHECK(ilreco_reconstruct(detector.workspace,
                             &on_hole_cell, 1, clusters, 4) == -1);
    CHECK(ilreco_reconstruct(detector.workspace,
                             &outside_circle, 1, clusters, 4) == -1);
}

TEST_CASE("mask: type labels derive from the mask", "[unit][mask]") {
    const auto mask = circular_detector_mask();
    MaskedContext detector(mask.data());

    // seed left of the hole -> missing in-grid neighbor -> type 1
    auto clusters = detector.run({{4, 4, 1.0}});
    REQUIRE(clusters.size() == 1);
    CHECK(clusters[0].type == 1);

    // seed on the circle rim (missing in-grid neighbors) -> type 1
    clusters = detector.run({{1, 2, 1.0}});
    REQUIRE(clusters.size() == 1);
    CHECK(clusters[0].type == 1);

    // seed on the bounding-box ring with all in-grid neighbors present -> type 2
    clusters = detector.run({{4, 0, 1.0}});
    REQUIRE(clusters.size() == 1);
    CHECK(clusters[0].type == 2);

    // fully surrounded interior seed -> type 0
    clusters = detector.run({{2, 4, 1.0}});
    REQUIRE(clusters.size() == 1);
    CHECK(clusters[0].type == 0);
}

TEST_CASE("mask: energy next to a hole is corrected up (honest containment)",
          "[unit][mask]") {
    const auto mask = circular_detector_mask();
    MaskedContext masked_detector(mask.data());
    MaskedContext plain_detector(nullptr);

    // shower hugging the hole from the left
    const std::vector<ilreco_hit> shower = {{4, 4, 1.000}, {3, 4, 0.180},
                                            {4, 5, 0.150}, {4, 3, 0.140},
                                            {3, 3, 0.050}, {3, 5, 0.045}};
    const auto masked = masked_detector.run(shower);
    const auto plain = plain_detector.run(shower);
    REQUIRE(masked.size() == 1);
    REQUIRE(plain.size() == 1);
    // The masked config knows cells (5,4)/(5,5) are ABSENT, not measured-zero:
    // its containment correction must be at least as large, and the energy
    // must stay inside the correction envelope.
    CHECK(masked[0].e >= plain[0].e);
    const double raw_sum = 1.565;
    CHECK(masked[0].e >= raw_sum - 1e-9);
    CHECK(masked[0].e <= raw_sum / 0.8 + 1e-9);
}

TEST_CASE("mask: all-cells-exist mask matches the unmasked physics", "[unit][mask]") {
    const std::vector<unsigned char> all_exist(GRID_SIZE * GRID_SIZE, 1);
    MaskedContext masked_detector(all_exist.data());
    MaskedContext plain_detector(nullptr);

    ilt::Rng rng(303);
    for (int trial = 0; trial < 50; ++trial) {
        const auto shower = ilt::synth_shower(0.5 + 10.0 * rng.uniform(),
                                              rng.uniform_int(0, GRID_SIZE - 1),
                                              rng.uniform_int(0, GRID_SIZE - 1), rng);
        std::vector<ilreco_hit> hits;   // synth spans +-2 cells; clip to the 9x9
        for (const auto& hit : shower)
            if (hit.col < GRID_SIZE && hit.row < GRID_SIZE) hits.push_back(hit);
        if (hits.empty()) continue;

        const auto masked = masked_detector.run(hits);
        const auto plain = plain_detector.run(hits);
        REQUIRE(masked.size() == plain.size());
        for (size_t i = 0; i < masked.size(); ++i) {
            // identical energies/positions (the zero-neighbor treatment is
            // unchanged by an all-cells-exist mask); type may differ by design
            // (mask-derived labels vs the built-in beam-hole pattern)
            REQUIRE(masked[i].e == plain[i].e);
            REQUIRE(masked[i].x == plain[i].x);
            REQUIRE(masked[i].y == plain[i].y);
        }
    }
}
