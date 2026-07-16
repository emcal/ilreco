// Speed benchmark: representative event categories timed through the full
// cluster_search + process_cluster chain. Prints a JSON object so CI can diff
// runs; compare against tests/benchmark/baseline.json (regenerate with
// --reps 2000 on a quiet machine when intentionally re-baselining).
//
//   ilreco_bench [--reps N]

#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "csv_io.h"

using Clock = std::chrono::steady_clock;

namespace {

double time_us_per_event(const std::vector<std::vector<ilt::Hit>>& events, int reps) {
    // warm-up
    for (const auto& ev : events) ilt::run_event(ev);
    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r)
        for (const auto& ev : events) ilt::run_event(ev);
    const double us =
        std::chrono::duration<double, std::micro>(Clock::now() - t0).count();
    return us / (double(reps) * events.size());
}

std::vector<std::vector<ilt::Hit>> make_category(const char* kind, int n, ilt::Rng& rng) {
    std::vector<std::vector<ilt::Hit>> events;
    for (int i = 0; i < n; ++i) {
        if (std::strcmp(kind, "single") == 0) {
            events.push_back(ilt::synth_shower(1.0 + 19.0 * rng.uniform(),
                                               rng.uniform_int(3, _NCOL_ - 2),
                                               rng.uniform_int(3, _NROW_ - 2), rng));
        } else if (std::strcmp(kind, "pair") == 0) {
            auto ev = ilt::synth_shower(1.0 + 9.0 * rng.uniform(),
                                        rng.uniform_int(3, _NCOL_ - 2),
                                        rng.uniform_int(3, _NROW_ - 2), rng);
            auto b = ilt::synth_shower(1.0 + 9.0 * rng.uniform(),
                                       rng.uniform_int(3, _NCOL_ - 2),
                                       rng.uniform_int(3, _NROW_ - 2), rng);
            ev.insert(ev.end(), b.begin(), b.end());
            events.push_back(ev);
        } else {  // dense: 30% occupancy noise + 3 showers
            std::vector<ilt::Hit> ev;
            for (int c = 1; c <= _NCOL_; ++c)
                for (int r = 1; r <= _NROW_; ++r)
                    if (rng.uniform() < 0.3) ev.push_back({c, r, 0.005 + 0.1 * rng.uniform()});
            for (int s = 0; s < 3; ++s) {
                auto sh = ilt::synth_shower(5.0, rng.uniform_int(3, _NCOL_ - 2),
                                            rng.uniform_int(3, _NROW_ - 2), rng);
                ev.insert(ev.end(), sh.begin(), sh.end());
            }
            events.push_back(ev);
        }
    }
    return events;
}

}  // namespace

int main(int argc, char** argv) {
    int reps = 500;
    for (int i = 1; i + 1 < argc; i += 2)
        if (std::strcmp(argv[i], "--reps") == 0) reps = std::atoi(argv[i + 1]);

    ilt::init_profile();
    ilt::Rng rng(777);

    std::printf("{\n  \"reps\": %d,\n  \"grid\": \"%dx%d\",\n", reps, _NCOL_, _NROW_);
    const char* kinds[] = {"single", "pair", "dense"};
    const int counts[] = {100, 100, 20};
    for (int k = 0; k < 3; ++k) {
        const auto events = make_category(kinds[k], counts[k], rng);
        double nhits = 0;
        for (const auto& ev : events) nhits += ev.size();
        const double us = time_us_per_event(events, reps);
        std::printf("  \"%s\": {\"us_per_event\": %.3f, \"events_per_sec\": %.0f, "
                    "\"mean_hits\": %.1f}%s\n",
                    kinds[k], us, 1e6 / us, nhits / events.size(), k < 2 ? "," : "");
    }
    std::printf("}\n");
    return 0;
}
