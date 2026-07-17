// Data-logic regression: every (input_*.csv, expected_*.csv) pair under
// tests/data is replayed through the library and compared event by event.
// The tables were produced by tools/gen_golden.cpp from simulated detector
// data (see tests/README.md for provenance); tolerances absorb FP noise
// across compilers, nothing more. A mismatch here means the physics output
// changed — that is the point of this suite.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "csv_io.h"

namespace fs = std::filesystem;

namespace {

std::vector<std::string> golden_set_names() {
    std::vector<std::string> names;
    for (const auto& entry : fs::directory_iterator(ILRECO_TEST_DATA)) {
        const auto filename = entry.path().filename().string();
        if (filename.rfind("input_", 0) == 0 && entry.path().extension() == ".csv")
            names.push_back(filename.substr(6, filename.size() - 10));
    }
    std::sort(names.begin(), names.end());
    return names;
}

}  // namespace

TEST_CASE("golden tables replay identically", "[data]") {
    const auto set_names = golden_set_names();
    if (set_names.empty()) {
        WARN("no golden tables found under " ILRECO_TEST_DATA
             " — run scripts in cal-fpga (750) + ilreco_gen_golden to create them");
        return;
    }
    for (const auto& name : set_names) {
        DYNAMIC_SECTION("set " << name) {
            const auto inputs = ilt::read_input_csv(
                std::string(ILRECO_TEST_DATA) + "/input_" + name + ".csv");
            const auto golden = ilt::read_expected_csv(
                std::string(ILRECO_TEST_DATA) + "/expected_" + name + ".csv");
            REQUIRE(!inputs.empty());
            REQUIRE(inputs.size() == golden.size());
            for (const auto& [event, hits] : inputs) {
                const auto clusters = ilt::run_event(hits);
                const auto& expected = golden.at(event);
                INFO("set " << name << " event " << event);
                REQUIRE((int)clusters.size() == expected.n_clusters);
                for (size_t i = 0; i < clusters.size(); ++i) {
                    CHECK_THAT(clusters[i].e,
                               Catch::Matchers::WithinRel(expected.clusters[i].e, 1e-6));
                    CHECK_THAT(clusters[i].x,
                               Catch::Matchers::WithinAbs(expected.clusters[i].x, 1e-5));
                    CHECK_THAT(clusters[i].y,
                               Catch::Matchers::WithinAbs(expected.clusters[i].y, 1e-5));
                    CHECK(clusters[i].size == expected.clusters[i].size);
                }
            }
        }
    }
}
