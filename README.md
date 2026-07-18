# ilreco


[![sanitizers + static-analysis + tests](https://github.com/emcal/ilreco/actions/workflows/ci.yml/badge.svg)](https://github.com/emcal/ilreco/actions/workflows/ci.yml)
[![docs](https://github.com/emcal/ilreco/actions/workflows/docs.yml/badge.svg)](https://github.com/emcal/ilreco/actions/workflows/docs.yml)


Island-clustering reconstruction for square-cell calorimeters (e.g. PbWO4 / lead
glass): connected-region search, peak finding, measured-shower-profile energy
sharing between overlapping showers, χ²-refined positions and containment-
corrected energies. Authored by Ilya Larin.

Drop in, plain C11, zero dependencies, thread-safe, allocation-free
per event. With C++ RAII wrapper and python bindings for your convenience! 


![ilreco overview](docs/ilreco-overview.svg)

![energy resolution benchmark](docs/energy-resolution-benchmark.png)

*Physics benchmark: energy resolution obtained with ilreco on a Geant4-simulated
10x10 PbWO4 array (2.055x2.055x20 cm crystals, measurement-matched response
model), electrons at 1-20 GeV — reproducing the HallD measured PbWO4 model
(dashed). Produced by the simulation-reconstruction regression that gates every
change to this library.*


[== DOCUMENTATION ==](https://emcal.github.io/ilreco/) 

Install (pip / plain C / CMake), tutorial, threading and memory model (for C/C++), Python API, shower profiles...



## Example: an 8×8 calorimeter with a central 2×2 beam hole

To give you a glimplse. The detector: an 8×8 crystal array with the central
2×2 block removed for the beam (the classic layout). The event: three photons,
landing on a hole-border cell, a normal interior cell, and an outer-edge
cell — every cell class appears, each with its own cluster in the output.

![tutorial event](docs/tutorial-8x8.png)

## Code example 

For C++, mt, and other examples read 
[== DOCUMENTATION ==](https://emcal.github.io/ilreco/)

```c
#include <ilreco.h>

char err[256];
/* once per process: geometry + shower-profile tables, immutable afterwards */
ilreco_config *cfg = ilreco_config_create(/*n_cols=*/30, /*n_rows=*/30,
                                          "data/prof_pwo.dat", err, sizeof err);
if (!cfg) { fprintf(stderr, "%s\n", err); return 1; }

/* once per thread: all working memory, one allocation */
ilreco_workspace *ws = ilreco_workspace_create(cfg);

/* per event: zero allocations */
ilreco_hit hits[] = { {12, 14, 1.85}, {13, 14, 0.42}, {12, 15, 0.31} };  /* 0-based col,row; GeV */
ilreco_cluster out[16];
int n = ilreco_reconstruct(ws, hits, 3, out, 16);
for (int i = 0; i < n && i < 16; ++i)
    printf("E=%.3f GeV at (%.2f, %.2f) cells, chi2=%.2f\n",
           out[i].e, out[i].x, out[i].y, out[i].chi2);

ilreco_workspace_destroy(ws);
ilreco_config_destroy(cfg);
```

Or from Python (`pip install ilreco`):

```python
import ilreco
calo = ilreco.Calorimeter(30, 30, profile="pwo")
clusters = calo.reconstruct(hits_table)   # numpy table in, numpy table out
```

Build:

```bash
cmake -B build && cmake --build build -j     # library + test suite
cd build && ctest                            # unit / smoke / golden-data / benchmark
```


[== DOCUMENTATION ==](https://emcal.github.io/ilreco/)
