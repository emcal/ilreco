// Cell-existence mask: arbitrary shapes (circular calorimeter, asymmetric
// hole — the EIC-B0-like case). Verifies input validation, mask-derived type
// labels, the physics effect (honest containment correction at real
// boundaries), and equivalence when the mask is all-ones.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "ilreco_test_utils.h"

namespace {

constexpr int N = 9;

// circular-ish 9x9 detector with an off-center 1x2 hole at (5,4),(5,5)
std::vector<unsigned char> b0_mask() {
    std::vector<unsigned char> m(N * N, 0);
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            if ((c - 4) * (c - 4) + (r - 4) * (r - 4) <= 18) m[r * N + c] = 1;
    m[4 * N + 5] = 0;   // hole (5,4)
    m[5 * N + 5] = 0;   // hole (5,5)
    return m;
}

struct Ctx {
    ilreco_config *cfg = nullptr;
    ilreco_workspace *ws = nullptr;
    explicit Ctx(const unsigned char *mask) {
        char err[256] = {0};
        cfg = ilreco_config_create(N, N, ILRECO_PROF_PWO, err, sizeof err);
        REQUIRE(cfg != nullptr);
        if (mask) REQUIRE(ilreco_config_set_cell_mask(cfg, mask) == 0);
        ws = ilreco_workspace_create(cfg);
        REQUIRE(ws != nullptr);
    }
    ~Ctx() {
        ilreco_workspace_destroy(ws);
        ilreco_config_destroy(cfg);
    }
    std::vector<ilreco_cluster> run(const std::vector<ilreco_hit> &hits) {
        std::vector<ilreco_cluster> out(16);
        const int n = ilreco_reconstruct(cfg, ws, hits.data(), (int)hits.size(),
                                         out.data(), 16);
        REQUIRE(n >= 0);
        out.resize(std::min(n, 16));
        return out;
    }
};

}  // namespace

TEST_CASE("mask: hits on missing cells are rejected", "[unit][mask]") {
    const auto m = b0_mask();
    Ctx c(m.data());
    ilreco_cluster out[4];
    const ilreco_hit on_hole{5, 4, 1.0};
    const ilreco_hit off_circle{0, 0, 1.0};
    CHECK(ilreco_reconstruct(c.cfg, c.ws, &on_hole, 1, out, 4) == -1);
    CHECK(ilreco_reconstruct(c.cfg, c.ws, &off_circle, 1, out, 4) == -1);
}

TEST_CASE("mask: type labels derive from the mask", "[unit][mask]") {
    const auto m = b0_mask();
    Ctx c(m.data());
    // seed left of the hole -> missing in-grid neighbor -> type 1
    auto out = c.run({{4, 4, 1.0}});
    REQUIRE(out.size() == 1);
    CHECK(out[0].type == 1);
    // seed on the circle rim (missing in-grid neighbors) -> type 1
    out = c.run({{1, 2, 1.0}});
    REQUIRE(out.size() == 1);
    CHECK(out[0].type == 1);
    // seed on the bounding-box ring with all in-grid neighbors present -> type 2
    out = c.run({{4, 0, 1.0}});
    REQUIRE(out.size() == 1);
    CHECK(out[0].type == 2);
    // fully surrounded interior seed -> type 0
    out = c.run({{2, 4, 1.0}});
    REQUIRE(out.size() == 1);
    CHECK(out[0].type == 0);
}

TEST_CASE("mask: energy next to a hole is corrected up (honest containment)",
          "[unit][mask]") {
    const auto m = b0_mask();
    Ctx masked(m.data());
    Ctx plain(nullptr);
    // shower hugging the hole from the left
    const std::vector<ilreco_hit> ev = {{4, 4, 1.000}, {3, 4, 0.180},
                                        {4, 5, 0.150}, {4, 3, 0.140},
                                        {3, 3, 0.050}, {3, 5, 0.045}};
    const auto em = masked.run(ev);
    const auto ep = plain.run(ev);
    REQUIRE(em.size() == 1);
    REQUIRE(ep.size() == 1);
    // masked config knows (5,4)/(5,5) are ABSENT (not measured-zero): the
    // containment correction must be at least as large, and the energies must
    // stay inside the correction envelope
    CHECK(em[0].e >= ep[0].e);
    const double raw = 1.565;
    CHECK(em[0].e >= raw - 1e-9);
    CHECK(em[0].e <= raw / 0.8 + 1e-9);
}

TEST_CASE("mask: all-ones mask matches the unmasked physics", "[unit][mask]") {
    std::vector<unsigned char> ones(N * N, 1);
    Ctx masked(ones.data());
    Ctx plain(nullptr);
    ilt::Rng rng(303);
    for (int t = 0; t < 50; ++t) {
        std::vector<ilreco_hit> ev;
        const auto sh = ilt::synth_shower(0.5 + 10.0 * rng.uniform(),
                                          rng.uniform_int(1, N),
                                          rng.uniform_int(1, N), rng);
        for (const auto &h : sh)
            if (h.col <= N && h.row <= N) ev.push_back({h.col - 1, h.row - 1, h.e});
        if (ev.empty()) continue;
        const auto a = masked.run(ev);
        const auto b = plain.run(ev);
        REQUIRE(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) {
            // identical energies/positions (fill_zero unchanged by all-ones mask);
            // type may differ by design (mask labels vs HyCal pattern)
            REQUIRE(a[i].e == b[i].e);
            REQUIRE(a[i].x == b[i].x);
            REQUIRE(a[i].y == b[i].y);
        }
    }
}
