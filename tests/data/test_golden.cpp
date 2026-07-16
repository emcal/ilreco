// Data-logic regression: every (input_*.csv, expected_*.csv) pair under
// tests/data is replayed through the library and compared event by event.
// The tables were produced by tools/gen_golden.cpp from simulated detector
// data (see tests/README.md for provenance); tolerances absorb FP noise
// across compilers, nothing more. A mismatch here means the physics output
// changed — that is the point of this suite.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <filesystem>
#include <string>
#include <vector>

#include "csv_io.h"

namespace fs = std::filesystem;

namespace {
struct Init {
    Init() { ilt::init_profile(); }
} init_once;

std::vector<std::string> golden_sets() {
    std::vector<std::string> names;
    for (const auto& p : fs::directory_iterator(ILRECO_TEST_DATA)) {
        const auto fn = p.path().filename().string();
        if (fn.rfind("input_", 0) == 0 && p.path().extension() == ".csv")
            names.push_back(fn.substr(6, fn.size() - 10));
    }
    std::sort(names.begin(), names.end());
    return names;
}
}  // namespace

TEST_CASE("golden tables replay identically", "[data]") {
    const auto sets = golden_sets();
    if (sets.empty()) {
        WARN("no golden tables found under " ILRECO_TEST_DATA
             " — run scripts in cal-fpga (750) + ilreco_gen_golden to create them");
        return;
    }
    for (const auto& name : sets) {
        DYNAMIC_SECTION("set " << name) {
            const auto inputs =
                ilt::read_input_csv(std::string(ILRECO_TEST_DATA) + "/input_" + name + ".csv");
            const auto golden =
                ilt::read_expected_csv(std::string(ILRECO_TEST_DATA) + "/expected_" + name + ".csv");
            REQUIRE(!inputs.empty());
            REQUIRE(inputs.size() == golden.size());
            for (const auto& [evt, hits] : inputs) {
                const auto out = ilt::run_event(hits);
                const auto& exp = golden.at(evt);
                INFO("set " << name << " event " << evt);
                REQUIRE((int)out.size() == exp.ncl);
                for (size_t i = 0; i < out.size(); ++i) {
                    CHECK_THAT(out[i].e,
                               Catch::Matchers::WithinRel(exp.clusters[i].e, 1e-6));
                    CHECK_THAT(out[i].x,
                               Catch::Matchers::WithinAbs(exp.clusters[i].x, 1e-5));
                    CHECK_THAT(out[i].y,
                               Catch::Matchers::WithinAbs(exp.clusters[i].y, 1e-5));
                    CHECK(out[i].size == exp.clusters[i].size);
                }
            }
        }
    }
}
