// CSV I/O for the golden-table format shared by tests/data and tests/tools.
//
// File formats (FROZEN — the fixtures predate the 0-based API and keep the
// historical 1-based columns/rows/positions; conversion to the library's
// 0-based convention happens here, at the file boundary, and nowhere else):
//   input:    event,col,row,e                 (col/row 1-based, e in GeV)
//   expected: event,icl,ncl,e,x,y,chi2,size   (x/y 1-based cell units,
//                                              icl = 0.. by descending e)
#ifndef ILRECO_TEST_CSV_IO_H
#define ILRECO_TEST_CSV_IO_H

#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "ilreco_test_utils.h"

namespace ilt {

inline std::map<long, std::vector<ilreco_hit>> read_input_csv(const std::string& path) {
    std::map<long, std::vector<ilreco_hit>> events;
    FILE* file = std::fopen(path.c_str(), "r");
    if (!file) return events;
    char line[256];
    if (!std::fgets(line, sizeof line, file)) { std::fclose(file); return events; }  // header
    long event;
    int col, row;
    double energy;
    while (std::fgets(line, sizeof line, file))
        if (std::sscanf(line, "%ld,%d,%d,%lf", &event, &col, &row, &energy) == 4)
            events[event].push_back({col - 1, row - 1, energy});  // file is 1-based
    std::fclose(file);
    return events;
}

struct GoldenEvent {
    int n_clusters;
    std::vector<ilreco_cluster> clusters;   // positions already 0-based
};

inline std::map<long, GoldenEvent> read_expected_csv(const std::string& path) {
    std::map<long, GoldenEvent> events;
    FILE* file = std::fopen(path.c_str(), "r");
    if (!file) return events;
    char line[256];
    if (!std::fgets(line, sizeof line, file)) { std::fclose(file); return events; }
    long event;
    int cluster_index, n_clusters, size;
    double e, x, y, chi2;
    while (std::fgets(line, sizeof line, file)) {
        if (std::sscanf(line, "%ld,%d,%d,%lf,%lf,%lf,%lf,%d", &event,
                        &cluster_index, &n_clusters, &e, &x, &y, &chi2, &size) == 8) {
            events[event].n_clusters = n_clusters;
            if (cluster_index >= 0)   // file positions are 1-based; -1.0 is exact
                events[event].clusters.push_back({e, x - 1.0, y - 1.0, chi2, size, 0});
        }
    }
    std::fclose(file);
    return events;
}

inline void write_expected_csv(const std::string& path,
                               const std::map<long, std::vector<ilreco_cluster>>& results) {
    FILE* file = std::fopen(path.c_str(), "w");
    std::fprintf(file, "event,icl,ncl,e,x,y,chi2,size\n");
    for (const auto& [event, clusters] : results) {
        if (clusters.empty()) {
            std::fprintf(file, "%ld,-1,0,0,0,0,0,0\n", event);
            continue;
        }
        for (size_t i = 0; i < clusters.size(); ++i)   // +1.0 back to the file convention
            std::fprintf(file, "%ld,%zu,%zu,%.9g,%.9g,%.9g,%.9g,%d\n", event, i,
                         clusters.size(), clusters[i].e, clusters[i].x + 1.0,
                         clusters[i].y + 1.0, clusters[i].chi2, clusters[i].size);
    }
    std::fclose(file);
}

}  // namespace ilt

#endif
