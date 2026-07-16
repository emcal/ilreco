#ifndef CALORIMETRY_STUDIES_ILRECO_H
#define CALORIMETRY_STUDIES_ILRECO_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * ilreco — island clustering for square-cell calorimeters
 *
 * Two APIs live here:
 *
 *   1. The CONTEXT API (ilreco_config / ilreco_workspace) — use this.
 *      - geometry and tunables are runtime configuration
 *      - all memory is allocated ONCE (config: shared, immutable after create;
 *        workspace: one per thread); event processing allocates nothing
 *      - thread-safe by construction: any number of threads may share one
 *        const config, each with its own workspace
 *      - 0-based cell indexing, errors returned (never exit())
 *
 *   2. The LEGACY API (cluster_search/process_cluster/read_profile_data) —
 *      the original 1-based, macro-configured, single-threaded interface,
 *      kept source- and behavior-compatible (including its quirks; see
 *      tests/KNOWN_ISSUES.md). It runs on an internal default context sized
 *      by the compile-time macros below. New code should not use it.
 * ==========================================================================*/

/* ---- compile-time defaults (legacy API geometry; also new-API defaults) --- */
#define _MADR0_   10000       /* number of calorimeter modules                */
#define _MADR_OFFSET_ 800     /* number of calorimeter modules safety offset  */
#define _MADR_    (_MADR0_+_MADR_OFFSET_)
#define _MCL_       200       /* max number+1 of raw clusters                 */
#define _MADCGAM_   200       /* max number+1 of reconstructed "particles"    */
#define _NCOL_       34       /* legacy grid columns                          */
#define _NROW_       34       /* legacy grid rows                             */
#define _OFFSET_    100       /* legacy column stride in packed addresses     */
#define _MAXLEN_    100       /* max hits stored per reconstructed object     */

#define _MPK_        12       /* max number of peaks in raw cluster           */

#define _N_REC_METHODS_ 5     /* number of coord. reconstruction methods      */
#define MIN_COUNTER_ENERGY  1.e-3   /* NOTE: used only by the test generator  */
#define CLUSTER_MIN_ENERGY  0.1     /* NOTE: historical, not used anywhere    */
#define _ZCAL_      732.      /* target-to-calorimeter distance (two-gamma
                                 separation mass cut); config field in the
                                 context API                                  */

#define PWO_CALOR 1
#define LG_CALOR  2
#ifndef CALOR_MATERIAL
#define CALOR_MATERIAL PWO_CALOR
#endif

#define _N_PROFILE_POINTS_  500 /* number of 2d profile nods                  */

/* ---- reconstructed object, shared by both APIs (layout unchanged) -------- */
typedef struct {
    double e;
    double x[_N_REC_METHODS_];
    double y[_N_REC_METHODS_];
    double z[_N_REC_METHODS_];
    double chi2;
    int    size;              /* number of hits                               */
    int    type;              /* fiducial-region classification (+10 if split)*/
    int      id;              /* reconstruction-path tag                      */
    int    stat;              /* status of cluster                            */
    int    element[_MAXLEN_]; /* hit addresses (legacy packed; slots 1..)     */
    double elfract[_MAXLEN_]; /* hit weights — units differ by path, see
                                 tests/KNOWN_ISSUES.md                        */
} adcgam_t;

/* ==========================================================================
 * CONTEXT API
 * ==========================================================================*/

typedef struct ilreco_config    ilreco_config;    /* opaque, immutable, shared */
typedef struct ilreco_workspace ilreco_workspace; /* opaque, one per thread    */

typedef struct {
    int    col;   /* 0-based column, 0 .. n_cols-1 */
    int    row;   /* 0-based row,    0 .. n_rows-1 */
    double e;     /* GeV */
} ilreco_hit;

typedef struct {
    double e;         /* GeV (profile containment correction applied)         */
    double x, y;      /* 0-based cell units: x == col means center of col     */
    double chi2;
    int    size;      /* number of hits in the cluster                        */
    int    type;      /* fiducial classification (0 inside, 2 boundary, ...)  */
} ilreco_cluster;

/* Create the shared configuration: geometry + shower-profile tables + tunables.
 * Allocates everything the config will ever need; immutable afterwards.
 * Returns NULL on failure with a message in errbuf (if given). */
ilreco_config *ilreco_config_create(int n_cols, int n_rows,
                                    const char *profile_path,
                                    char *errbuf, size_t errbuf_len);
void ilreco_config_destroy(ilreco_config *cfg);

/* Optional tuning — call before the config is shared across threads. */
void ilreco_config_set_zcal(ilreco_config *cfg, double zcal);
void ilreco_config_set_seed_threshold(ilreco_config *cfg, double min_seed_gev);
void ilreco_config_set_hole_classification(ilreco_config *cfg, int enabled);

/* Cell-existence mask for non-rectangular shapes (circular calorimeters,
 * asymmetric beam holes, dead regions). mask[row*n_cols + col] != 0 means the
 * cell physically exists; the array is copied. Effects:
 *   - missing cells are excluded from the zero-signal neighbor treatment, so
 *     containment corrections and chi2 near real boundaries become honest
 *     (energy lost into a hole is corrected UP, not treated as measured-zero);
 *   - cluster `type` is derived from the mask: 1 = seed has a missing in-grid
 *     neighbor (hole border / rim), 2 = seed on the bounding-box ring, 0 =
 *     fully surrounded. The built-in HyCal hole pattern is ignored when a
 *     mask is set;
 *   - hits on masked-out cells are rejected by ilreco_reconstruct (-1).
 * Returns 0 on success, -1 on invalid input. */
int ilreco_config_set_cell_mask(ilreco_config *cfg, const unsigned char *mask);

/* Per-thread workspace: ONE allocation (internal arena), reused for every
 * event. Never share a workspace between threads. */
ilreco_workspace *ilreco_workspace_create(const ilreco_config *cfg);
void ilreco_workspace_destroy(ilreco_workspace *ws);

/* Reconstruct one event. Zero allocations. Returns the number of clusters
 * found (may exceed max_out; out receives at most max_out, energy-descending),
 * or -1 on invalid input (hit outside the grid, negative n_hits, ...). */
int ilreco_reconstruct(const ilreco_config *cfg, ilreco_workspace *ws,
                       const ilreco_hit *hits, int n_hits,
                       ilreco_cluster *out, int max_out);

/* ==========================================================================
 * LEGACY API — original interface, unchanged semantics (single-threaded;
 * grid fixed to _NCOL_ x _NROW_; addresses = col*_OFFSET_+row, 1-based;
 * arrays used from index 1).
 * ==========================================================================*/
void  read_profile_data(const char* prof_file_name);
int  cluster_search(int nw, int *ia, double *id, int *lencl);
void process_cluster(int nadc, int *ia, double *id, int *nadcgam, adcgam_t *adcgam);
int  read_event(int *nw, int *ia, double *id);
void dump_clusters(int nw, int *ia, double *id, int ncl, int *lencl, int nadcgam, adcgam_t *adcgam);

#ifdef __cplusplus
}
#endif

#endif /* CALORIMETRY_STUDIES_ILRECO_H */
