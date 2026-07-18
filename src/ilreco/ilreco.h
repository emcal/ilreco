#ifndef CALORIMETRY_STUDIES_ILRECO_H
#define CALORIMETRY_STUDIES_ILRECO_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * ilreco — island-clustering reconstruction for square-cell calorimeters.
 *
 * Usage model:
 *   - geometry and tunables are runtime configuration
 *   - all memory is allocated ONCE (config: shared, immutable after create;
 *     workspace: one per thread); event processing allocates nothing
 *   - thread-safe by construction: any number of threads may share one
 *     const config, each with its own workspace
 *   - 0-based cell indexing, errors returned (never exit())
 * See README.md for a step-by-step tutorial.
 * ==========================================================================*/

typedef struct ilreco_config    ilreco_config;    /* opaque, immutable, shared */
typedef struct ilreco_workspace ilreco_workspace; /* opaque, one per thread    */

typedef struct {
    int32_t col;   /* 0-based column, 0 .. n_cols-1 */
    int32_t row;   /* 0-based row,    0 .. n_rows-1 */
    double  e;     /* GeV */
} ilreco_hit;

typedef struct {
    double  e;        /* GeV (profile containment correction applied)         */
    double  x;        /* 0-based cell units: x == col means center of col     */
    double  y;        /* 0-based cell units: y == row means center of row     */
    double  chi2;     /* profile-fit chi2 / ndof                              */
    int32_t size;     /* number of hits in the cluster                        */
    int32_t type;     /* fiducial classification: 0 interior, 1 hole border,
                         2 outer boundary; +10 if a split gamma pair          */
} ilreco_cluster;

/* Create the shared configuration: geometry + shower-profile tables + tunables.
 * Allocates everything the config will ever need; immutable afterwards.
 *   n_cols, n_rows   grid size in cells
 *   profile_path     shower-profile table file (see README, "The shower
 *                    profile (weights) file")
 *   errbuf, errbuf_len   optional: receives the failure message
 * Returns NULL on failure with a message in errbuf (if given). */
ilreco_config *ilreco_config_create(int32_t n_cols, int32_t n_rows,
                                    const char *profile_path,
                                    char *errbuf, size_t errbuf_len);
void ilreco_config_destroy(ilreco_config *cfg);

/* Optional tuning — call before the config is shared across threads. */

/* Target-to-calorimeter distance [cm]; enters the two-gamma invariant-mass
 * cut that decides whether a chi2-improving split is physical (default 732,
 * PrimEx). */
void ilreco_config_set_zcal(ilreco_config *cfg, double zcal);

/* Cluster-seed threshold [GeV] (default 0.01): a cluster is reconstructed
 * only if one cell exceeds it (scaled up with cluster energy; see README
 * pitfall 1). */
void ilreco_config_set_seed_threshold(ilreco_config *cfg, double min_seed_gev);

/* Enable/disable the built-in central-2x2-beam-hole labeling (type = 1;
 * default enabled). Labels only — energies and positions are never
 * affected. Ignored when a cell mask is set. */
void ilreco_config_set_hole_classification(ilreco_config *cfg, int32_t enabled);

/* Cell-existence mask for non-rectangular shapes (circular calorimeters,
 * asymmetric beam holes, dead regions). mask[row*n_cols + col] != 0 means the
 * cell physically exists; the array is copied. Effects:
 *   - missing cells are excluded from the zero-signal neighbor treatment, so
 *     containment corrections and chi2 near real boundaries become honest
 *     (energy lost into a hole is corrected UP, not treated as measured-zero);
 *   - cluster `type` is derived from the mask: 1 = seed has a missing in-grid
 *     neighbor (hole border / rim), 2 = seed on the bounding-box ring, 0 =
 *     fully surrounded. The built-in central-hole pattern is ignored when a
 *     mask is set;
 *   - hits on masked-out cells are rejected by ilreco_reconstruct (-1).
 * Returns 0 on success, -1 on invalid input. */
int32_t ilreco_config_set_cell_mask(ilreco_config *cfg, const unsigned char *mask);

/* Per-thread workspace: ONE allocation (internal arena), reused for every
 * event. The workspace keeps a reference to cfg — after this call the config
 * is only needed to spawn further workspaces, and it must outlive every
 * workspace created from it. Never share a workspace between threads. */
ilreco_workspace *ilreco_workspace_create(const ilreco_config *cfg);
void ilreco_workspace_destroy(ilreco_workspace *ws);

/* Reconstruct one event on the workspace's calorimeter. Zero allocations.
 *   hits, n_hits   fired cells, any order, one entry per cell
 *   out, max_out   receives up to max_out clusters, energy-descending
 * Returns the number of clusters found (may exceed max_out), or -1 on
 * invalid input (hit outside the grid or on a masked-out cell, negative
 * n_hits, ...). */
int32_t ilreco_reconstruct(ilreco_workspace *ws,
                           const ilreco_hit *hits, int32_t n_hits,
                           ilreco_cluster *out, int32_t max_out);

#ifdef __cplusplus
}
#endif

#endif /* CALORIMETRY_STUDIES_ILRECO_H */
