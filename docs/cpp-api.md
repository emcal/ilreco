# C++ API

`ilreco.hpp` is a header-only RAII wrapper over the C API: same physics,
same structs (`ilreco_hit`, `ilreco_cluster`), with lifetimes managed by
objects and errors reported as exceptions.

```cpp
#include <ilreco.hpp>

#include <cstdio>
#include <vector>

int main() {
    ilreco::Config config(20, 20, "prof_pwo.dat");   // throws on failure
    ilreco::Workspace workspace(config);             // once per thread

    const std::vector<ilreco_cluster> clusters =
        workspace.reconstruct({{12, 14, 1.85}, {13, 14, 0.42}, {12, 15, 0.31}});

    for (const auto& cluster : clusters)
        std::printf("E=%.3f GeV at (%.2f, %.2f), chi2=%.2f, type=%d\n",
                    cluster.e, cluster.x, cluster.y, cluster.chi2, cluster.type);
}
```

Both classes are move-only; destruction order is enforced by scope
(workspaces die before the config they came from — see the
[memory model](./memory-model)).

## Configuration

All tuning happens on `Config`, before the first `Workspace` is created:

```cpp
ilreco::Config config(15, 15, "prof_pwo.dat");

// arbitrary shapes: nonzero = cell exists, row-major n_rows*n_cols
std::vector<unsigned char> mask = build_my_detector_shape();
config.set_cell_mask(mask);

config.set_seed_threshold(0.02);   // GeV; default 0.01
config.set_zcal(732.0);            // cm; two-gamma separation cut
config.set_hole_classification(false);   // built-in central-2x2 labels
```

## Errors

| condition | behavior |
|---|---|
| profile file missing/corrupt, bad geometry | `Config` constructor throws `std::runtime_error` |
| allocation failure | `Workspace` constructor throws `std::runtime_error` |
| hit outside the grid or on a masked-out cell | `reconstruct` throws `std::invalid_argument` |

## The C API underneath

The wrapper adds nothing semantically — everything it does maps 1:1 onto
six C functions, documented in `ilreco.h`:

```c
ilreco_config *ilreco_config_create(int32_t n_cols, int32_t n_rows,
                                    const char *profile_path,
                                    char *errbuf, size_t errbuf_len);
void ilreco_config_destroy(ilreco_config *cfg);
int32_t ilreco_config_set_cell_mask(ilreco_config *cfg, const unsigned char *mask);
ilreco_workspace *ilreco_workspace_create(const ilreco_config *cfg);
void ilreco_workspace_destroy(ilreco_workspace *ws);
int32_t ilreco_reconstruct(ilreco_workspace *ws,
                           const ilreco_hit *hits, int32_t n_hits,
                           ilreco_cluster *out, int32_t max_out);
```

Use the C API directly when you need C linkage or want zero abstraction;
use the wrapper everywhere else.
