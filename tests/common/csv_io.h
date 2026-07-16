// Minimal CSV I/O for the golden-table format shared by tests/data and
// tests/tools. Formats:
//   input:    event,col,row,e            (1-based col/row, e in GeV)
//   expected: event,icl,ncl,e,x,y,chi2,size   (icl = 0.. by descending e)
#ifndef ILRECO_TEST_CSV_IO_H
#define ILRECO_TEST_CSV_IO_H

#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "ilreco_test_utils.h"

namespace ilt {

inline std::map<long, std::vector<Hit>> read_input_csv(const std::string& path) {
    std::map<long, std::vector<Hit>> events;
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) return events;
    char line[256];
    if (!std::fgets(line, sizeof line, f)) { std::fclose(f); return events; }  // header
    long evt; int col, row; double e;
    while (std::fgets(line, sizeof line, f))
        if (std::sscanf(line, "%ld,%d,%d,%lf", &evt, &col, &row, &e) == 4)
            events[evt].push_back({col, row, e});
    std::fclose(f);
    return events;
}

struct GoldenRow {
    int ncl;
    std::vector<Cluster> clusters;
};

inline std::map<long, GoldenRow> read_expected_csv(const std::string& path) {
    std::map<long, GoldenRow> events;
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) return events;
    char line[256];
    if (!std::fgets(line, sizeof line, f)) { std::fclose(f); return events; }
    long evt; int icl, ncl, size; double e, x, y, chi2;
    while (std::fgets(line, sizeof line, f)) {
        if (std::sscanf(line, "%ld,%d,%d,%lf,%lf,%lf,%lf,%d",
                        &evt, &icl, &ncl, &e, &x, &y, &chi2, &size) == 8) {
            events[evt].ncl = ncl;
            if (icl >= 0) events[evt].clusters.push_back({e, x, y, chi2, size});
        }
    }
    std::fclose(f);
    return events;
}

inline void write_expected_csv(const std::string& path,
                               const std::map<long, std::vector<Cluster>>& results) {
    FILE* f = std::fopen(path.c_str(), "w");
    std::fprintf(f, "event,icl,ncl,e,x,y,chi2,size\n");
    for (const auto& [evt, cls] : results) {
        if (cls.empty()) {
            std::fprintf(f, "%ld,-1,0,0,0,0,0,0\n", evt);
            continue;
        }
        for (size_t i = 0; i < cls.size(); ++i)
            std::fprintf(f, "%ld,%zu,%zu,%.9g,%.9g,%.9g,%.9g,%d\n", evt, i,
                         cls.size(), cls[i].e, cls[i].x, cls[i].y, cls[i].chi2,
                         cls[i].size);
    }
    std::fclose(f);
}

}  // namespace ilt

#endif
