// Speed benchmark: representative event categories timed through the full
// reconstruction chain. Prints a JSON object so CI can diff runs; compare
// against tests/benchmark/baseline.json (regenerate with --reps 2000 on a
// quiet machine when intentionally re-baselining).
//
//   ilreco_bench [--reps N]

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "ilreco_test_utils.h"

using Clock = std::chrono::steady_clock;
using ilt::GRID_COLS;
using ilt::GRID_ROWS;

namespace {

double time_us_per_event(const std::vector<std::vector<ilreco_hit>>& events, int reps) {
    for (const auto& event : events) ilt::run_event(event);   // warm-up
    const auto start = Clock::now();
    for (int rep = 0; rep < reps; ++rep)
        for (const auto& event : events) ilt::run_event(event);
    const double total_us =
        std::chrono::duration<double, std::micro>(Clock::now() - start).count();
    return total_us / (double(reps) * events.size());
}

std::vector<std::vector<ilreco_hit>> make_category(const char* kind, int n_events,
                                                   ilt::Rng& rng) {
    std::vector<std::vector<ilreco_hit>> events;
    for (int i = 0; i < n_events; ++i) {
        if (std::strcmp(kind, "single") == 0) {
            events.push_back(ilt::synth_shower(1.0 + 19.0 * rng.uniform(),
                                               rng.uniform_int(2, GRID_COLS - 3),
                                               rng.uniform_int(2, GRID_ROWS - 3), rng));
        } else if (std::strcmp(kind, "pair") == 0) {
            auto event = ilt::synth_shower(1.0 + 9.0 * rng.uniform(),
                                           rng.uniform_int(2, GRID_COLS - 3),
                                           rng.uniform_int(2, GRID_ROWS - 3), rng);
            const auto second = ilt::synth_shower(1.0 + 9.0 * rng.uniform(),
                                                  rng.uniform_int(2, GRID_COLS - 3),
                                                  rng.uniform_int(2, GRID_ROWS - 3), rng);
            event.insert(event.end(), second.begin(), second.end());
            events.push_back(event);
        } else {  // dense: 30% occupancy noise + 3 showers
            std::vector<ilreco_hit> event;
            for (int col = 0; col < GRID_COLS; ++col)
                for (int row = 0; row < GRID_ROWS; ++row)
                    if (rng.uniform() < 0.3)
                        event.push_back({col, row, 0.005 + 0.1 * rng.uniform()});
            for (int s = 0; s < 3; ++s) {
                const auto shower = ilt::synth_shower(5.0,
                                                      rng.uniform_int(2, GRID_COLS - 3),
                                                      rng.uniform_int(2, GRID_ROWS - 3),
                                                      rng);
                event.insert(event.end(), shower.begin(), shower.end());
            }
            events.push_back(event);
        }
    }
    return events;
}

}  // namespace

int main(int argc, char** argv) {
    int reps = 500;
    for (int i = 1; i + 1 < argc; i += 2)
        if (std::strcmp(argv[i], "--reps") == 0) reps = std::atoi(argv[i + 1]);

    ilt::Rng rng(777);

    std::printf("{\n  \"reps\": %d,\n  \"grid\": \"%dx%d\",\n", reps, GRID_COLS,
                GRID_ROWS);
    const char* kinds[] = {"single", "pair", "dense"};
    const int counts[] = {100, 100, 20};
    for (int k = 0; k < 3; ++k) {
        const auto events = make_category(kinds[k], counts[k], rng);
        double total_hits = 0;
        for (const auto& event : events) total_hits += event.size();
        const double us = time_us_per_event(events, reps);
        std::printf("  \"%s\": {\"us_per_event\": %.3f, \"events_per_sec\": %.0f, "
                    "\"mean_hits\": %.1f}%s\n",
                    kinds[k], us, 1e6 / us, total_hits / events.size(),
                    k < 2 ? "," : "");
    }
    std::printf("}\n");
    return 0;
}
