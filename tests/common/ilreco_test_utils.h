// Shared helpers for the ilreco test suite. This is the ONLY place the tests
// touch ilreco's raw calling convention, so the contract is written down once:
//
//   * arrays are used from index 1 (ia[1..nw], id[1..nw]; element 0 unused)
//   * address = column*_OFFSET_ + row, column/row 1-based, 1..(_NCOL_/_NROW_)
//   * energies in GeV
//   * cluster positions x[1]/y[1] are in column/row units (x == column means
//     the center of that column)
//   * read_profile_data(prof_pwo.dat) must be called once before processing
#ifndef ILRECO_TEST_UTILS_H
#define ILRECO_TEST_UTILS_H

#include <ilreco.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace ilt {

struct Hit {
    int col;     // 1-based
    int row;     // 1-based
    double e;    // GeV
};

struct Cluster {
    double e;
    double x;    // column units
    double y;    // row units
    double chi2;
    int size;
};

inline void init_profile() {
    static bool done = false;
    if (!done) {
        read_profile_data(ILRECO_PROF_PWO);
        done = true;
    }
}

inline int address(int col, int row) { return col * _OFFSET_ + row; }

// Run one event through cluster_search + process_cluster, exactly like the
// library's own driver loop. Returns clusters sorted by energy, descending.
inline std::vector<Cluster> run_event(const std::vector<Hit>& hits) {
    static int ia[_MADR_];
    static double id[_MADR_];
    static int lencl[_MCL_];
    static adcgam_t adcgam[_MADCGAM_];

    const int nw = static_cast<int>(hits.size());
    for (int i = 0; i < nw; ++i) {
        ia[i + 1] = address(hits[i].col, hits[i].row);
        id[i + 1] = hits[i].e;
    }

    int nadcgam = 0;
    const int ncl = cluster_search(nw, ia, id, lencl);
    for (int icl = 1, ipncl = 1; icl <= ncl && nadcgam < _MADCGAM_ - 2;
         ipncl += lencl[icl++])
        process_cluster(lencl[icl], &ia[ipncl - 1], &id[ipncl - 1], &nadcgam, adcgam);

    std::vector<Cluster> out;
    for (int g = 1; g <= nadcgam; ++g)
        out.push_back({adcgam[g].e, adcgam[g].x[1], adcgam[g].y[1],
                       adcgam[g].chi2, adcgam[g].size});
    std::sort(out.begin(), out.end(),
              [](const Cluster& a, const Cluster& b) { return a.e > b.e; });
    return out;
}

inline bool finite(const Cluster& c) {
    return std::isfinite(c.e) && std::isfinite(c.x) && std::isfinite(c.y);
}

// A deterministic little LCG so smoke inputs are reproducible without <random>.
struct Rng {
    unsigned long long s;
    explicit Rng(unsigned long long seed) : s(seed) {}
    double uniform() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return double((s >> 11) & 0x1FFFFFFFFFFFFFULL) / double(0x20000000000000ULL);
    }
    int uniform_int(int lo, int hi) {  // inclusive
        return lo + int(uniform() * (hi - lo + 1));
    }
};

// Synthetic shower: energy-weighted blob around (c, r) with rough transverse
// profile — good enough for smoke/benchmark inputs (real-data fixtures live in
// tests/data/).
inline std::vector<Hit> synth_shower(double e_gev, int c, int r, Rng& rng) {
    std::vector<Hit> hits;
    for (int dc = -2; dc <= 2; ++dc) {
        for (int dr = -2; dr <= 2; ++dr) {
            const int cc = c + dc, rr = r + dr;
            if (cc < 1 || cc > _NCOL_ || rr < 1 || rr > _NROW_) continue;
            const double w = std::exp(-1.6 * std::hypot(double(dc), double(dr)));
            const double e = e_gev * w * (0.8 + 0.4 * rng.uniform());
            if (e > 0.005) hits.push_back({cc, rr, e});
        }
    }
    return hits;
}

}  // namespace ilt

#endif
