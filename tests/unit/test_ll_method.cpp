// Log-likelihood method contract tests: ILLLPROF asset loading (validation,
// error reporting), method selection rules, and reconstruction behavior with
// synthetic LL tables (energy path shared with chi2, uninformative tables
// keep the moments position, method switch requires loaded tables).
//
// The physics performance of the method against real per-energy tables is
// validated in the Python comparison study, not here.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

#include "ilreco_test_utils.h"

namespace {

constexpr int kNprof = 500;
constexpr int kEntries = (kNprof + 1) * (kNprof + 2) / 2;

// Minimal valid ILLLPROF asset: one energy, n_bins bins, every node has
// nll(f) = bin_nll[k] and nll(0) = nll_zero.
std::string write_synthetic_asset(const std::vector<double>& bin_nll,
                                  float nll_zero, double energy = 5.0) {
    const std::string path =
        (std::filesystem::temp_directory_path() / "ilreco_ll_test.bin").string();
    FILE* fh = std::fopen(path.c_str(), "wb");
    REQUIRE(fh != nullptr);
    const uint32_t head[5] = {1u, 1u, (uint32_t)kNprof,
                              (uint32_t)bin_nll.size(), (uint32_t)kEntries};
    const double nll_scale = 32.0 / 65535.0;
    const double grid[3] = {std::log(1e-6), 0.023, nll_scale};
    std::fwrite("ILLLPROF", 1, 8, fh);
    std::fwrite(head, sizeof head[0], 5, fh);
    std::fwrite(grid, sizeof grid[0], 3, fh);
    std::fwrite(&energy, sizeof energy, 1, fh);
    std::vector<float> zeros(kEntries, nll_zero);
    std::fwrite(zeros.data(), sizeof(float), zeros.size(), fh);
    std::vector<uint16_t> row(bin_nll.size());
    for (size_t k = 0; k < bin_nll.size(); ++k)
        row[k] = (uint16_t)std::lround(bin_nll[k] / nll_scale);
    for (int i = 0; i < kEntries; ++i)
        std::fwrite(row.data(), sizeof(uint16_t), row.size(), fh);
    std::fclose(fh);
    return path;
}

}  // namespace

TEST_CASE("LL asset loader rejects garbage and truncation", "[unit][ll]") {
    ilt::TestContext context(10, 10);
    char error[256] = {0};

    SECTION("missing file") {
        CHECK(ilreco_config_load_ll_profiles(context.config,
                                             "/nonexistent/llprof.bin",
                                             error, sizeof error) == -1);
        CHECK(std::strlen(error) > 0);
    }
    SECTION("wrong magic") {
        const std::string path =
            (std::filesystem::temp_directory_path() / "ilreco_ll_bad.bin").string();
        FILE* fh = std::fopen(path.c_str(), "wb");
        std::fwrite("NOTAPROF11111111111111111111", 1, 28, fh);
        std::fclose(fh);
        CHECK(ilreco_config_load_ll_profiles(context.config, path.c_str(),
                                             error, sizeof error) == -1);
        CHECK(std::string(error).find("ILLLPROF") != std::string::npos);
    }
    SECTION("truncated energy block") {
        const std::string good = write_synthetic_asset({1.0, 2.0, 3.0}, 4.0f);
        std::filesystem::resize_file(good,
                                     std::filesystem::file_size(good) / 2);
        CHECK(ilreco_config_load_ll_profiles(context.config, good.c_str(),
                                             error, sizeof error) == -1);
        CHECK(std::string(error).find("truncated") != std::string::npos);
    }
}

TEST_CASE("method selection requires loaded LL tables", "[unit][ll]") {
    ilt::TestContext context(10, 10);
    CHECK(ilreco_config_set_method(context.config,
                                   ILRECO_METHOD_LOG_LIKELIHOOD) == -1);
    CHECK(ilreco_config_set_method(context.config, 42) == -1);
    CHECK(ilreco_config_set_method(context.config,
                                   ILRECO_METHOD_CHI_SQUARE) == 0);

    char error[256] = {0};
    const std::string asset = write_synthetic_asset({0.5, 1.0, 2.0}, 3.0f);
    REQUIRE(ilreco_config_load_ll_profiles(context.config, asset.c_str(),
                                           error, sizeof error) == 0);
    CHECK(ilreco_config_set_method(context.config,
                                   ILRECO_METHOD_LOG_LIKELIHOOD) == 0);
}

TEST_CASE("LL reconstruction shares the energy path with chi2 and keeps the "
          "moments position under uninformative tables", "[unit][ll]") {
    // flat nll => every position is equally likely => the fit must not move
    const std::vector<ilreco_hit> shower = {
        {7, 7, 0.10}, {8, 7, 0.55}, {9, 7, 0.12},
        {7, 8, 0.08}, {8, 8, 3.60}, {9, 8, 0.30},
        {7, 9, 0.02}, {8, 9, 0.25}, {9, 9, 0.05},
    };
    ilreco_cluster chi2_out[4];
    ilreco_cluster ll_out[4];

    ilt::TestContext chi2_context(20, 20);
    REQUIRE(ilreco_reconstruct(chi2_context.workspace, shower.data(),
                               (int)shower.size(), chi2_out, 4) == 1);

    ilt::TestContext ll_context(20, 20);
    char error[256] = {0};
    const std::string asset =
        write_synthetic_asset(std::vector<double>(64, 1.0), 1.0f);
    REQUIRE(ilreco_config_load_ll_profiles(ll_context.config, asset.c_str(),
                                           error, sizeof error) == 0);
    REQUIRE(ilreco_config_set_method(ll_context.config,
                                     ILRECO_METHOD_LOG_LIKELIHOOD) == 0);
    REQUIRE(ilreco_reconstruct(ll_context.workspace, shower.data(),
                               (int)shower.size(), ll_out, 4) == 1);

    // energy: identical moments + containment path
    CHECK_THAT(ll_out[0].e,
               Catch::Matchers::WithinRel(chi2_out[0].e, 1e-12));
    // flat likelihood: position stays within a node step of the moments seed;
    // the chi2 fit moves away from it, so only coarse agreement is expected
    CHECK_THAT(ll_out[0].x, Catch::Matchers::WithinAbs(chi2_out[0].x, 0.35));
    CHECK_THAT(ll_out[0].y, Catch::Matchers::WithinAbs(chi2_out[0].y, 0.35));
    CHECK(ll_out[0].chi2 >= 0.0);
}
