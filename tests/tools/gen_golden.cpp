// Golden-table generator: runs the CURRENT library over an input CSV and
// freezes the responses as the expected CSV. Run this ONLY when intentionally
// re-baselining (e.g. after an agreed physics-affecting change); commit the
// regenerated tables together with the change that caused them.
//
//   ilreco_gen_golden input_central.csv expected_central.csv
//
// Golden CSV I/O comes from common/csv_io.h (which converts the 1-based file
// format to the 0-based API); infrastructure map: tests/README.md.

#include <cstdio>
#include <map>
#include <vector>

#include "csv_io.h"

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: ilreco_gen_golden <input.csv> <expected.csv>\n");
        return 2;
    }
    const auto events = ilt::read_input_csv(argv[1]);
    if (events.empty()) {
        std::fprintf(stderr, "no events read from %s\n", argv[1]);
        return 3;
    }
    std::map<long, std::vector<ilreco_cluster>> results;
    for (const auto& [event, hits] : events) results[event] = ilt::run_event(hits);
    ilt::write_expected_csv(argv[2], results);
    std::fprintf(stderr, "%zu events -> %s\n", events.size(), argv[2]);
    return 0;
}
