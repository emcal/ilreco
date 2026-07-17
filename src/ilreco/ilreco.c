//
//  ilreco — island-clustering reconstruction for square-cell calorimeters.
//  Public interface and its documentation: ilreco.h.
//
//  File layout:
//    1. contexts: ilreco_config (shared, immutable) / ilreco_workspace
//       (per-thread buffers) and their create/destroy/reconstruct entry points
//    2. the reconstruction algorithm
//
//  The algorithm section is intentionally identical, expression for
//  expression, to the original PrimEx/HyCal code: evaluation order and the
//  integer quantization steps are part of the contract, and the results are
//  pinned bit-exact by the golden tables in tests/data. Do not "clean up"
//  arithmetic there without re-baselining those tables.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ilreco.h"

/* ===================== internal constants ================================= */

#define _MCL_       200       /* cluster-length array capacity                */
#define _MADCGAM_   200       /* output-object array capacity                 */
#define _MAXLEN_    100       /* max hits stored per reconstructed object     */
#define _MPK_        12       /* max peaks per raw cluster                    */
#define _N_REC_METHODS_ 5     /* coordinate-reconstruction method slots       */
#define _ZCAL_      732.      /* default target-to-calorimeter distance [cm]  */
#define NPROF       500       /* profile table nodes per axis (0..NPROF)      */

#define PWO_CALOR 1
#define LG_CALOR  2
#ifndef CALOR_MATERIAL
#define CALOR_MATERIAL PWO_CALOR  /* selects the seed-threshold energy scaling */
#endif

/* internal representation of one reconstructed object */
typedef struct {
    double  e;                 /* energy [GeV], containment-corrected          */
    double  x[_N_REC_METHODS_];/* position per method; x[1] = standard, in
                                  1-based column units                         */
    double  y[_N_REC_METHODS_];/* same, row units                              */
    double  z[_N_REC_METHODS_];
    double  chi2;              /* profile-fit chi2 / ndof                      */
    int32_t size;              /* number of hits                               */
    int32_t type;              /* fiducial-region classification (+10 if split)*/
    int32_t id;                /* reconstruction-path tag                      */
    int32_t stat;              /* status of cluster                            */
    int32_t element[_MAXLEN_]; /* hit addresses (packed; slots 1..)            */
    double  elfract[_MAXLEN_]; /* hit weights — units differ by path, see
                                  tests/KNOWN_ISSUES.md                        */
} adcgam_t;

/* ============================ contexts ==================================== */

struct ilreco_config {
    int32_t ncol;          /* grid columns                                    */
    int32_t nrow;          /* grid rows                                       */
    int32_t offset;        /* packed 1-based cell address = col*offset + row;
                              offset > nrow+1 so addresses never collide      */
    int32_t madr;          /* capacity of the per-event hit/address buffers   */
    int32_t madr0;         /* peak-count budget: a cluster of nadc hits may
                              have at most madr0/nadc - 3 peaks, bounding the
                              peaks x cluster-length work per cluster         */
    int32_t mcl;           /* cluster-length array capacity: at most mcl-1
                              raw clusters per event                          */
    int32_t madcgam;       /* output-object array capacity: at most madcgam-2
                              reconstructed objects per event                 */
    double zcal;           /* target-to-calorimeter distance [cm] (two-gamma
                              invariant-mass separation cut)                  */
    double min_seed_gev;   /* cluster seed threshold [GeV] (historical 0.01)  */
    double seed_scale;     /* seed-threshold energy scaling factor            */
    int32_t hole_enabled;  /* label seeds around a central 2x2 beam hole
                              as type 1 (classic layout)                       */
    unsigned char *cellmask; /* NULL = all cells exist; else [nrow*ncol]
                                row-major, 0 = cell physically missing        */
    double *amean;         /* (NPROF+1)^2 shower-profile mean energy fraction */
    double *ad2c;          /* (NPROF+1)^2 shower-profile fraction variance    */
};

struct ilreco_workspace {
    int32_t stride;        /* row length of iwrk/idp/fwrk; == cfg->madr       */
    int32_t mcl;           /* copies of the config capacities, for checking   */
    int32_t madcgam;       /*   that workspace and config match               */
    /* per-event buffers (filled by ilreco_reconstruct) */
    int32_t *ia;           /* [madr]  packed cell addresses, used from 1      */
    double  *id;           /* [madr]  cell energies [GeV], parallel to ia     */
    int32_t *lencl;        /* [mcl]   raw-cluster lengths                     */
    adcgam_t *adcgam;      /* [madcgam] reconstructed objects                 */
    int32_t *order;        /* [madcgam] energy-descending sort scratch        */
    /* cluster_search scratch */
    int32_t *iwork;        /* [madr]  subcluster-glue address shuffle         */
    double  *dwork;        /* [madr]  subcluster-glue energy shuffle          */
    /* process_cluster / gamma_cluster scratch */
    int32_t *iazero;       /* [8*madr+2] zero-signal neighbor addresses       */
    int32_t *iwrk;         /* [(_MPK_+1) x stride] per-peak integer weights   */
    int32_t *idp;          /* [(_MPK_+1) x stride] per-peak deposit shares    */
    double  *fwrk;         /* [(_MPK_+1) x stride] per-peak float weights     */
    void *arena;           /* single allocation backing all of the above     */
};

#define AMEAN(cf, i, j) ((cf)->amean[(size_t)(i) * (NPROF + 1) + (j)])
#define AD2C(cf, i, j)  ((cf)->ad2c [(size_t)(i) * (NPROF + 1) + (j)])
#define IWRK(wk, p, i)  ((wk)->iwrk[(size_t)(p) * (wk)->stride + (i)])
#define IDP(wk, p, i)   ((wk)->idp [(size_t)(p) * (wk)->stride + (i)])
#define FWRK(wk, p, i)  ((wk)->fwrk[(size_t)(p) * (wk)->stride + (i)])

/* ===================== internal declarations ==============================
 *
 * Common parameters of the algorithm functions (original conventions, arrays
 * used from index 1):
 *   ia[]      packed 1-based cell addresses: ia[i] = col*cf->offset + row
 *   id[]      cell energies [GeV], parallel to ia[]
 *   nw/nadc   number of valid entries in ia[]/id[] (indices 1..nw)
 *   lencl[]   hit count of each raw cluster (lencl[1..ncl])
 *   iazero[]  addresses of existing zero-signal cells neighboring the
 *             cluster, nzero entries — they enter the containment correction
 *             and the chi2 as real zero measurements
 *   ipnpk[]   ia-indices of the peak (local-maximum) cells, 1..npk
 *   adcgam[]/nadcgam  output objects and their running count
 */

/* stage 1: group hits into connected clusters; reorders ia/id in place
 * (address-ascending, clusters contiguous); returns the cluster count */
static int32_t cluster_search_core(const ilreco_config *cf, ilreco_workspace *wk,
                                   int32_t nw, int32_t *ia, double *id,
                                   int32_t *lencl);

/* stage 2: split one raw cluster into reconstructed objects (peak finding,
 * then iterative profile-based energy sharing between the peaks) */
static void process_cluster_core(const ilreco_config *cf, ilreco_workspace *wk,
                                 int32_t nadc, int32_t *ia, double *id,
                                 int32_t *nadcgam, adcgam_t *adcgam);

/* insertion sort of ia (ascending) with id following */
static void order_hits(int32_t nw, int32_t *ia, double *id);

/* fiducial classification of a seed cell (1-based col ix, row iy):
 * 0 = interior, 1 = beam-hole border, 2 = outer boundary */
static int32_t peak_type(const ilreco_config *cf, int32_t ix, int32_t iy);

/* round to nearest integer, halves away from zero */
static int32_t nint(double x);

/* fit one cluster: energy/position via profile moments + chi2 descent.
 * In:  itype from peak_type; nadc/ia/id = the cluster's cells.
 * Out: e1/x1/y1 = the gamma (positions in 1-based cell units);
 *      chisq in = reference chi2 to beat, out = final chi2/ndof;
 *      e2/x2/y2 = second gamma if a physical two-gamma split was found,
 *      e2 = 0 otherwise. */
static void gamma_cluster(const ilreco_config *cf, ilreco_workspace *wk,
                          int32_t itype, int32_t nadc, int32_t *ia, double *id,
                          double *chisq, double *e1, double *x1, double *y1,
                          double *e2, double *x2, double *y2);

/* collect existing zero-signal neighbors of the cluster into iazero[] */
static void fill_zero_hits(const ilreco_config *cf, int32_t nadc, int32_t *ia,
                           int32_t *nzero, int32_t *iazero);

/* first moments: energy-weighted position + containment-corrected energy */
static void mom1_cluster(const ilreco_config *cf, int32_t nadc, int32_t *ia,
                         double *id, int32_t nzero, int32_t *iazero,
                         double *e, double *x, double *y);

/* first + second moments (only used by the dead two-gamma stub below) */
static void mom2_cluster(const ilreco_config *cf, int32_t nadc, int32_t *ia,
                         double *id, int32_t nzero, int32_t *iazero,
                         double *e, double *x, double *y,
                         double *xx, double *yy, double *xy);

/* chi2 of the (e, x, y) single-shower hypothesis against the profile */
static double chisq1_cluser(const ilreco_config *cf, int32_t nadc, int32_t *ia,
                            double *id, int32_t nzero, int32_t *iazero,
                            double e, double x, double y);

/* second-order two-gamma separation — intentionally an empty stub, see its
 * definition */
static void tgamma_cluster(const ilreco_config *cf, int32_t nadc, int32_t *ia,
                           double *id, int32_t nzero, int32_t *iazero,
                           double *chisq, double *ee, double *xx, double *yy,
                           double *e2, double *x2, double *y2);

/* bilinear interpolation of the profile tables at cell-unit offset (x, y) */
static double profile_mean(const ilreco_config *cf, double x, double y);
static double d2c(const ilreco_config *cf, double x, double y);

/* expected variance of the measured energy fraction at offset (dx, dy) */
static double sigma2(const ilreco_config *cf, double dx, double dy, double e);

/* emit the reconstructed object(s) for a single-peak / multi-peak cluster */
static void add_single_adcgam(const ilreco_config *cf, ilreco_workspace *wk,
                              int32_t nadc, int32_t *ia, double *id,
                              int32_t *ipnpk, int32_t *nadcgam, adcgam_t *adcgam);
static void add_many_adcgam(const ilreco_config *cf, ilreco_workspace *wk,
                            int32_t nadc, int32_t *ia, double *id, int32_t npk,
                            int32_t *ipnpk, int32_t *nadcgam, adcgam_t *adcgam);

/* ====================== config / workspace =============================== */

static int32_t load_profile(ilreco_config *cf, const char *path,
                            char *errbuf, size_t errlen) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        if (errbuf) snprintf(errbuf, errlen,
                             "profile file '%s' does not exist or access is denied", path);
        return -1;
    }
    for (int i = 0; i <= NPROF; ++i)
        for (int j = 0; j <= i; ++j) {
            int i1, j1; double f1, f2;
            if (fscanf(fp, "%d %d %lf %lf", &i1, &j1, &f1, &f2) != 4 ||
                i1 != i || j1 != j) {
                if (errbuf) snprintf(errbuf, errlen,
                                     "profile data file corruption in '%s' at (%d,%d)",
                                     path, i, j);
                fclose(fp);
                return -1;
            }
            AMEAN(cf, i, j) = f1;
            AMEAN(cf, j, i) = f1;
            AD2C(cf, i, j)  = f2;
            AD2C(cf, j, i)  = f2;
        }
    fclose(fp);
    return 0;
}

ilreco_config *ilreco_config_create(int32_t n_cols, int32_t n_rows,
                                    const char *profile_path,
                                    char *errbuf, size_t errbuf_len) {
    if (n_cols < 1 || n_rows < 1) {
        if (errbuf) snprintf(errbuf, errbuf_len, "invalid geometry %dx%d",
                             n_cols, n_rows);
        return NULL;
    }
    ilreco_config *cf = calloc(1, sizeof *cf);
    if (!cf) return NULL;
    cf->ncol = n_cols;
    cf->nrow = n_rows;
    cf->offset = n_rows + 2;
    cf->madr = (n_cols + 1) * cf->offset + n_rows + 2 + 64;
    /* peak budget scaled far above physical occupancies: the guard only ever
     * binds for pathological events with >= _MPK_-2 peaks in a huge island */
    const int64_t budget = 100LL * n_cols * n_rows;
    cf->madr0 = budget < INT32_MAX ? (int32_t)budget : INT32_MAX;
    cf->mcl = _MCL_;
    cf->madcgam = _MADCGAM_;
    cf->zcal = _ZCAL_;
    cf->min_seed_gev = 0.01;
#if CALOR_MATERIAL == LG_CALOR
    cf->seed_scale = 20.;
#else
    cf->seed_scale = 7.;
#endif
    cf->hole_enabled = 1;
    cf->amean = calloc((size_t)(NPROF + 1) * (NPROF + 1), sizeof(double));
    cf->ad2c  = calloc((size_t)(NPROF + 1) * (NPROF + 1), sizeof(double));
    if (!cf->amean || !cf->ad2c) {
        ilreco_config_destroy(cf);
        return NULL;
    }
    if (!profile_path) {
        if (errbuf) snprintf(errbuf, errbuf_len, "no shower-profile file given");
        ilreco_config_destroy(cf);
        return NULL;
    }
    if (load_profile(cf, profile_path, errbuf, errbuf_len) != 0) {
        ilreco_config_destroy(cf);
        return NULL;
    }
    return cf;
}

void ilreco_config_destroy(ilreco_config *cfg) {
    if (!cfg) return;
    free(cfg->cellmask);
    free(cfg->amean);
    free(cfg->ad2c);
    free(cfg);
}

/* 1-based column/row -> does the cell exist? (mask absent => yes) */
static int32_t cell_exists(const ilreco_config *cf, int32_t ix, int32_t iy) {
    if (!cf->cellmask) return 1;
    return cf->cellmask[(size_t)(iy - 1) * cf->ncol + (ix - 1)] != 0;
}

int32_t ilreco_config_set_cell_mask(ilreco_config *cfg, const unsigned char *mask) {
    if (!cfg || !mask) return -1;
    if (!cfg->cellmask) {
        cfg->cellmask = malloc((size_t)cfg->ncol * cfg->nrow);
        if (!cfg->cellmask) return -1;
    }
    memcpy(cfg->cellmask, mask, (size_t)cfg->ncol * cfg->nrow);
    return 0;
}

void ilreco_config_set_zcal(ilreco_config *cfg, double zcal) { cfg->zcal = zcal; }
void ilreco_config_set_seed_threshold(ilreco_config *cfg, double min_seed_gev) {
    cfg->min_seed_gev = min_seed_gev;
}
void ilreco_config_set_hole_classification(ilreco_config *cfg, int32_t enabled) {
    cfg->hole_enabled = enabled;
}

ilreco_workspace *ilreco_workspace_create(const ilreco_config *cfg) {
    if (!cfg) return NULL;
    ilreco_workspace *wk = calloc(1, sizeof *wk);
    if (!wk) return NULL;
    const size_t m = (size_t)cfg->madr;
    wk->stride = cfg->madr;
    wk->mcl = cfg->mcl;
    wk->madcgam = cfg->madcgam;
    const size_t n_int = m /*ia*/ + m /*iwork*/ + (size_t)cfg->mcl /*lencl*/
                       + (8 * m + 2) /*iazero*/
                       + (size_t)(_MPK_ + 1) * m /*iwrk*/
                       + (size_t)(_MPK_ + 1) * m /*idp*/
                       + (size_t)cfg->madcgam /*order*/;
    const size_t n_dbl = m /*id*/ + m /*dwork*/ + (size_t)(_MPK_ + 1) * m /*fwrk*/;
    const size_t bytes = n_dbl * sizeof(double) + n_int * sizeof(int32_t)
                       + (size_t)cfg->madcgam * sizeof(adcgam_t);
    char *p = malloc(bytes);
    if (!p) { free(wk); return NULL; }
    wk->arena = p;
    /* doubles first (strictest alignment), then adcgam_t, then ints */
    wk->id = (double *)p;                 p += m * sizeof(double);
    wk->dwork = (double *)p;              p += m * sizeof(double);
    wk->fwrk = (double *)p;               p += (size_t)(_MPK_ + 1) * m * sizeof(double);
    wk->adcgam = (adcgam_t *)p;           p += (size_t)cfg->madcgam * sizeof(adcgam_t);
    wk->ia = (int32_t *)p;                    p += m * sizeof(int32_t);
    wk->iwork = (int32_t *)p;                 p += m * sizeof(int32_t);
    wk->lencl = (int32_t *)p;                 p += (size_t)cfg->mcl * sizeof(int32_t);
    wk->iazero = (int32_t *)p;                p += (8 * m + 2) * sizeof(int32_t);
    wk->iwrk = (int32_t *)p;                  p += (size_t)(_MPK_ + 1) * m * sizeof(int32_t);
    wk->idp = (int32_t *)p;                   p += (size_t)(_MPK_ + 1) * m * sizeof(int32_t);
    wk->order = (int32_t *)p;
    return wk;
}

void ilreco_workspace_destroy(ilreco_workspace *ws) {
    if (!ws) return;
    free(ws->arena);
    free(ws);
}

int32_t ilreco_reconstruct(const ilreco_config *cfg, ilreco_workspace *ws,
                           const ilreco_hit *hits, int32_t n_hits,
                           ilreco_cluster *out, int32_t max_out) {
    if (!cfg || !ws || n_hits < 0 || (n_hits > 0 && !hits) ||
        ws->stride < cfg->madr || n_hits > cfg->madr - 2)
        return -1;
    for (int32_t i = 0; i < n_hits; ++i) {
        if (hits[i].col < 0 || hits[i].col >= cfg->ncol ||
            hits[i].row < 0 || hits[i].row >= cfg->nrow)
            return -1;
        if (!cell_exists(cfg, hits[i].col + 1, hits[i].row + 1))
            return -1;
        ws->ia[i + 1] = (hits[i].col + 1) * cfg->offset + (hits[i].row + 1);
        ws->id[i + 1] = hits[i].e;
    }

    int32_t nadcgam = 0;
    const int32_t ncl = cluster_search_core(cfg, ws, n_hits, ws->ia, ws->id,
                                            ws->lencl);
    for (int32_t icl = 1, ipncl = 1; icl <= ncl && nadcgam < cfg->madcgam - 2;
         ipncl += ws->lencl[icl++])
        process_cluster_core(cfg, ws, ws->lencl[icl], &ws->ia[ipncl - 1],
                             &ws->id[ipncl - 1], &nadcgam, ws->adcgam);

    /* order the output by energy, descending (stable insertion sort) */
    int32_t *order = ws->order;
    for (int32_t g = 1; g <= nadcgam; ++g) order[g - 1] = g;
    for (int32_t sorted_end = 1; sorted_end < nadcgam; ++sorted_end) {
        const int32_t key = order[sorted_end];
        int32_t slot = sorted_end - 1;
        while (slot >= 0 && ws->adcgam[order[slot]].e < ws->adcgam[key].e) {
            order[slot + 1] = order[slot];
            --slot;
        }
        order[slot + 1] = key;
    }
    const int32_t n_out = nadcgam < max_out ? nadcgam : max_out;
    for (int32_t k = 0; k < n_out; ++k) {
        const adcgam_t *g = &ws->adcgam[order[k]];
        out[k].e = g->e;
        out[k].x = g->x[1] - 1.0;   /* 1-based column units -> 0-based */
        out[k].y = g->y[1] - 1.0;
        out[k].chi2 = g->chi2;
        out[k].size = g->size;
        out[k].type = g->type;
    }
    return nadcgam;
}

/* ========================== algorithm ===================================== */

static int32_t cluster_search_core(const ilreco_config *cf, ilreco_workspace *wk,
                                   int32_t nw, int32_t *ia, double *id,
                                   int32_t *lencl) {
  if(nw<2) {
    lencl[1] = 1;
    return nw;
  }
  order_hits(nw,ia,id);       // addresses must be in increasing order

  int32_t *iwork = wk->iwork;  double *dwork = wk->dwork;
  int32_t ncl = 0, next  = 1, iak = 0;
  for(int32_t k = 2; k <= nw+1; ++k) {
    if(k<=nw) iak = ia[k];
    if(iak-ia[k-1]<=1 && k<=nw) continue;

    int32_t ib  = next; //  first word of the (sub)cluster
    int32_t ie  = k-1;  //  last  word of the (sub)cluster
    next    = k;    //  first word of the next (sub)cluster
    if(ncl>cf->mcl-2) {
      printf("maximum number of clusters reached\n");
      return cf->mcl-1;
    }
    lencl[++ncl]= next-ib;    // length of the (sub)cluster
    if(ncl==1)  continue;
//
//  glue subclusters
//
    int32_t ias     = ia[ib];
    int32_t iaf     = ia[ie];
    int32_t last    = ib-1;
    int32_t lastcl  = ncl-1;

    for(int32_t icl = lastcl; icl>0; --icl) {
      int32_t leng  = lencl[icl];
      if(ias-ia[last] > cf->offset) break;  //  no subclusters to glue
      for(int32_t i = last; i > last-leng; --i) {
        if(ias-ia[i]  >  cf->offset) break;
        if(iaf-ia[i]  >= cf->offset) {      //  subclusters to glue
          if(icl<ncl-1 && leng<cf->madr) {

            memmove(&iwork[1],&ia[last+1-leng],leng*sizeof(int32_t));
            memmove(&ia[last+1-leng],&ia[last+1],(ib-1-last)*sizeof(int32_t));
            memmove(&ia[ib-leng],&iwork[1],leng*sizeof(int32_t));
            memcpy(&dwork[1],&id[last+1-leng],leng*sizeof(double));
            memmove(&id[last+1-leng],&id[last+1],(ib-1-last)*sizeof(double));
            memcpy(&id[ib-leng],&dwork[1],leng*sizeof(double));

            for(int32_t j = icl; j < ncl-1; ++j) lencl[j] = lencl[j+1];
          }
          ib  -= leng;
          lencl[ncl-1]  = lencl[ncl] + leng;
          --ncl;
          break;
        }
      }
      last  -= leng;                    //  last word of tested subcluster
    }
  }
  return ncl;
}

static void order_hits(int32_t nw, int32_t *ia, double *id) {

  for(int32_t k = 2; k <= nw; ++k) {
    if(ia[k] > ia[k-1])  continue;
    int32_t     iat = ia[k];
    double  idt = id[k];
    for(int32_t i = k-1; i>=0; --i)
      if(i) {
        if(iat<ia[i]) {
          ia[i+1] = ia[i];
          id[i+1] = id[i];
        } else {
          ia[i+1] = iat;
          id[i+1] = idt;
          break;
        }
      } else {
        ia[1] = iat;
        id[1] = idt;
      }
  }

  return;
}

static void process_cluster_core(const ilreco_config *cf, ilreco_workspace *wk,
                                 int32_t nadc, int32_t *ia, double *id,
                                 int32_t *nadcgam, adcgam_t *adcgam) {
  order_hits(nadc,ia,id);     // addresses must be in increasing order

  double minpk = cf->min_seed_gev;
  int32_t ipnpk[_MPK_];

  if(nadc>=3) {
    double idsum = 0.;
    for(int32_t i = 1; i <= nadc; ++i) idsum += id[i];
    double ib = cf->seed_scale*log(1.+idsum);
    if(ib>1) minpk *= ib;
    minpk = nint(1.e2*minpk);
    minpk *= 1.e2;
    minpk = 1.e-4*minpk;
  }

  int32_t npk = 0;
  for(int32_t ic = 1; ic <= nadc; ++ic) {
    double iac = id[ic];
    if(iac<minpk) continue;
    int32_t ixy     = ia[ic];
    int32_t ixymax  = ixy + cf->offset + 1;
    int32_t ixymin  = ixy - cf->offset - 1;
    int32_t iyc     = ixy - (ixy/cf->offset)*cf->offset;

    int32_t in;
    in  = ic  + 1;
    int32_t skipflag = 0;
    while(in<=nadc && ia[in] <= ixymax) {
      int32_t iy  = ia[in] - (ia[in]/cf->offset)*cf->offset;
      if(abs(iy-iyc)<=1 && id[in]>=iac) {skipflag = 1; break;}
      ++in;
    }
    if(skipflag)  continue;
    in = ic - 1;
    while (in>=1 && ia[in] >= ixymin) {
      int32_t iy  = ia[in] - (ia[in]/cf->offset)*cf->offset;
      if(abs(iy-iyc)<=1 && id[in]>iac)  {skipflag = 1; break;}
      --in;
    }
    if(skipflag) continue;

    ++npk;          // peak found
    ipnpk[npk] = ic;
    if(npk >= _MPK_-2 || npk >= cf->madr0/nadc - 3) break;
  }
  if(!npk || *nadcgam >= cf->madcgam-3)  return;

  if(npk==1) add_single_adcgam(cf,wk,nadc,ia,id,ipnpk,nadcgam,adcgam);
  else     add_many_adcgam(cf,wk,nadc,ia,id,npk,ipnpk,nadcgam,adcgam);

  return;
}

static void add_single_adcgam(const ilreco_config *cf, ilreco_workspace *wk,
                              int32_t nadc, int32_t *ia, double *id, int32_t *ipnpk,
                              int32_t *nadcgam, adcgam_t *adcgam) {

  const double  chisq1  = 90.; // 3.0 value of chi2 for preliminary seperation

  int32_t n1  = ++(*nadcgam);
  int32_t ic  = ipnpk[1];
  int32_t ix  = ia[ic]/cf->offset;
  int32_t iy  = ia[ic]-ix*cf->offset;
  int32_t itype = peak_type(cf,ix,iy);

  double  e1 = 0., x1 = 0., y1 = 0., e2 = 0., x2 = 0., y2 = 0., chisq = chisq1;
  gamma_cluster(cf,wk,itype,nadc,ia,id,&chisq,&e1,&x1,&y1,&e2,&x2,&y2);

  adcgam[n1].e    = e1;
  adcgam[n1].x[1] = x1;
  adcgam[n1].y[1] = y1;
  adcgam[n1].chi2 = chisq;
  adcgam[n1].size = nadc;
  adcgam[n1].type = itype;
  adcgam[n1].id   = 0;
  adcgam[n1].stat = itype;

  if(e2>0. && *nadcgam < cf->madcgam-3) {
    int32_t n2        = ++(*nadcgam);
    adcgam[n2].e    = e2;
    adcgam[n2].x[1] = x2;
    adcgam[n2].y[1] = y2;
    adcgam[n2].chi2 = chisq;
    adcgam[n2].size = nadc;
    adcgam[n1].type = itype+10;
    adcgam[n2].type = itype+10;
    adcgam[n1].id   = 1;
    adcgam[n2].id   = 2;
    adcgam[n2].stat = itype;

    adcgam[n1].x[0] = 0.5*(x1-x2);
    adcgam[n1].y[0] = 0.5*(y1-y2);
    adcgam[n2].x[0] = 0.5*(x2-x1);
    adcgam[n2].y[0] = 0.5*(y2-y1);

    for(int32_t j = 1; j <= nadc; ++j) {
      if(j>=_MAXLEN_) break;
      adcgam[n1].element[j] = ia[j];
      adcgam[n2].element[j] = ia[j];
      adcgam[n1].elfract[j] = id[j]*e1/(e1+e2);
      adcgam[n2].elfract[j] = id[j]*e2/(e1+e2);
    }
  } else {
    for(int32_t j = 1; j <= nadc; ++j) {
      if(j>=_MAXLEN_) break;
      adcgam[n1].element[j] = ia[j];
      adcgam[n1].elfract[j] = nint(1.e4*id[j]);
    }
  }
  return;
}

static void add_many_adcgam(const ilreco_config *cf, ilreco_workspace *wk,
                            int32_t nadc, int32_t *ia, double *id, int32_t npk, int32_t *ipnpk,
                            int32_t *nadcgam, adcgam_t *adcgam) {

  const double  idelta  =  0.; //  min cell energy part to be a member of separated cluster
  const double  chisq1  = 90.; // 3.0 value of chi2 for preliminary seperation
  const double  chisq2  = 50.; // 0.8 value of chi2 for final seperation
  const int32_t     niter   =   6;
  int32_t ngam0 = *nadcgam;
  int32_t igmpk[_MPK_][3];

  double ratio = 1.;
  for(int32_t iter = 1; iter <= niter; ++iter) {
    for(int32_t i = 1; i <= nadc; ++i) {IWRK(wk,0,i) = 0; FWRK(wk,0,i) = 0.;}
    double epk[_MPK_], xpk[_MPK_], ypk[_MPK_];
    for(int32_t ipk = 1; ipk <= npk; ++ipk) {
      int32_t ic = ipnpk[ipk];
      if(iter!=1) ratio = FWRK(wk,ipk,ic) / FWRK(wk,npk+1,ic);
      double eg = id[ic] * ratio;
      int32_t ixypk = ia[ic];
      int32_t ixpk  = ixypk / cf->offset;
      int32_t iypk  = ixypk - ixpk * cf->offset;
      epk[ipk]  = eg;
      xpk[ipk]  = eg*ixpk;
      ypk[ipk]  = eg*iypk;
      for(int32_t in = ic + 1; in <= nadc; ++in) {
        int32_t ixy = ia[in];
        int32_t ix  = ixy / cf->offset;
        int32_t iy  = ixy - ix * cf->offset;
        if(ixy-ixypk>cf->offset+1) break;
        if(abs(iy-iypk)>1) continue;
        if(iter!=1) ratio = FWRK(wk,ipk,in) / FWRK(wk,npk+1,in);
        eg = id[in] * ratio;
        epk[ipk] += eg;
        xpk[ipk] += eg*ix;
        ypk[ipk] += eg*iy;
      }

      for(int32_t in = ic - 1; in>0; --in) {
        int32_t ixy = ia[in];
        int32_t ix  = ixy / cf->offset;
        int32_t iy  = ixy - ix * cf->offset;
        if(ixypk-ixy>cf->offset+1) break;
        if(abs(iy-iypk)>1) continue;
        if(iter!=1) ratio = FWRK(wk,ipk,in) / FWRK(wk,npk+1,in);
        eg = id[in] * ratio;
        epk[ipk] += eg;
        xpk[ipk] += eg*ix;
        ypk[ipk] += eg*iy;
      }

      if(epk[ipk]>0.) {
        xpk[ipk] /= epk[ipk];
        ypk[ipk] /= epk[ipk];
      } else {
        printf("WRN: lost maximum with peak energy %lf at ITER = %d\n", id[ic], iter);
      }

      for(int32_t i = 1; i <= nadc; ++i) {
        int32_t ixy = ia[i];
        int32_t ix  = ixy / cf->offset;
        int32_t iy  = ixy - ix * cf->offset;
        double dx = fabs(ix-xpk[ipk]);
        double dy = fabs(iy-ypk[ipk]);
        double a  = epk[ipk]*profile_mean(cf,dx,dy);
        IWRK(wk,ipk,i) = nint(a*1.e4);
        IWRK(wk,0,i)  += IWRK(wk,ipk,i);
        FWRK(wk,ipk,i) = a;
        FWRK(wk,0,i)  += FWRK(wk,ipk,i);
      }
    }

    for(int32_t i = 1; i <= nadc; ++i) {
      int32_t iwk = IWRK(wk,0,i);
      if(iwk<1) iwk = 1;
      IWRK(wk,npk+1,i) = iwk;
      if(FWRK(wk,0,i)>1.e-6)
        FWRK(wk,npk+1,i) = FWRK(wk,0,i);
      else
        FWRK(wk,npk+1,i) = 1.e-6;
    }
  }

  for(int32_t ipk = 1; ipk <= npk; ++ipk) {
    int32_t leng = 0;
    for(int32_t i = 1; i <= nadc; ++i) {
      if(FWRK(wk,0,i) <= 1.e-6) continue;
      int32_t ixy   = ia[i];
      double fe = 1.e4*id[i]*FWRK(wk,ipk,i)/FWRK(wk,0,i);
      if(fe<=idelta)  continue;
      ++leng;
      IWRK(wk,npk+1,leng) = ixy;
      IWRK(wk,npk+2,leng) = nint(fe);
      FWRK(wk,npk+1,leng) = ixy;
      FWRK(wk,npk+2,leng) = 1.e-4*fe;
    }

    if(*nadcgam >= cf->madcgam-2) return;
    igmpk[ipk][2] = 0;
    if(!leng) continue;
    int32_t n1  = ++(*nadcgam);
    int32_t ic  = ipnpk[ipk];
    int32_t ix  = ia[ic] / cf->offset;
    int32_t iy  = ia[ic] - ix * cf->offset;
    int32_t itype = peak_type(cf,ix,iy);

    double  e1 = 0., x1 = 0., y1 = 0., e2 = 0., x2 = 0., y2 = 0., chisq = chisq1;
    gamma_cluster(cf,wk,itype,leng,&IWRK(wk,npk+1,0),&FWRK(wk,npk+2,0),&chisq,&e1,&x1,&y1,&e2,&x2,&y2);

    adcgam[n1].e    = e1;
    adcgam[n1].x[1] = x1;
    adcgam[n1].y[1] = y1;
    adcgam[n1].chi2 = chisq;
    adcgam[n1].size = leng;
    adcgam[n1].type = itype;
    adcgam[n1].id   = 90;
    adcgam[n1].stat = itype;

    igmpk[ipk][1]   = n1;
    igmpk[ipk][2]   = n1;

    if(e2>0. && n1 < cf->madcgam-3) {
      int32_t n2        = ++(*nadcgam);
      adcgam[n2].e    = e2;
      adcgam[n2].x[1] = x2;
      adcgam[n2].y[1] = y2;
      adcgam[n2].chi2 = chisq;
      adcgam[n2].size = leng;
      adcgam[n2].type = itype;
      adcgam[n1].id   = 91;
      adcgam[n2].id   = 92;
      adcgam[n2].stat = itype;
      igmpk[ipk][2]   = n2;
      adcgam[n1].x[0] = 0.5*(x1-x2);
      adcgam[n1].y[0] = 0.5*(y1-y2);
      adcgam[n2].x[0] = 0.5*(x2-x1);
      adcgam[n2].y[0] = 0.5*(y2-y1);
    }
  }

  for(int32_t i = 1; i <= nadc; ++i) {IWRK(wk,0,i) = 0; IDP(wk,0,i) = 0; FWRK(wk,0,i) = 0.;}

  for(int32_t ipk = 1; ipk <= npk; ++ipk)
    for(int32_t i = 1; i <= nadc; ++i) {
      IWRK(wk,ipk,i) = 0; IDP(wk,ipk,i) = 0; FWRK(wk,ipk,i) = 0.;
      if(igmpk[ipk][2])
        for(int32_t ig = igmpk[ipk][1]; ig <= igmpk[ipk][2]; ++ig) {

          int32_t ixy = ia[i];
          int32_t ix  = ixy / cf->offset;
          int32_t iy  = ixy - ix * cf->offset;
          double dx = ix-adcgam[ig].x[1];
          double dy = iy-adcgam[ig].y[1];
          double a  = adcgam[ig].e*profile_mean(cf,dx,dy);
          int32_t   iia = nint(a*1.e4);
          IWRK(wk,ipk,i) += iia;
          IWRK(wk,0,i)   += iia;
          FWRK(wk,ipk,i) += a;
          FWRK(wk,0,i)   += a;
          IDP(wk,ipk,i)  += iia;
        }
    }

  for(int32_t i = 1; i <= nadc; ++i) {
    IDP(wk,0,i) = 0;
    for(int32_t ipk = 1; ipk <= npk; ++ipk)
      IDP(wk,0,i) += IDP(wk,ipk,i);
    int32_t ide = nint(id[i]*1.e4) - IDP(wk,0,i);
    if(!ide || FWRK(wk,0,i) == 0.) continue;

    double fw[_MPK_];
    for(int32_t ipk = 1; ipk <= npk; ++ipk)
      fw[ipk] = FWRK(wk,ipk,i)/FWRK(wk,0,i);

    int32_t idecorr = 0;
    for(int32_t ipk = 1; ipk <= npk; ++ipk) {
      double fia = ide*fw[ipk];
      if(FWRK(wk,ipk,i)+fia>0.) {
        FWRK(wk,ipk,i) += fia;
        FWRK(wk,0,i)   += fia;
      }
      int32_t iia = nint(fia);
      if(IWRK(wk,ipk,i)+iia>0)  {
        IWRK(wk,ipk,i) += iia;
        IWRK(wk,0,i)   += iia;
        idecorr        += iia;
      } else {
        if(IWRK(wk,ipk,i)+iia<0)
          printf("WARNING NEGATIVE CORR = %i %lf\n", ia[i], id[i]);
      }
    }
  }

  *nadcgam  = ngam0;          //  reanalize last (two) gamma(s)
  for(int32_t ipk = 1; ipk <= npk; ++ipk) {
    int32_t leng = 0;
    for(int32_t i = 1; i <= nadc; ++i) {
      if(IWRK(wk,0,i) <= 0) continue;
      double fe = 1.e4*id[i]*FWRK(wk,ipk,i)/FWRK(wk,0,i);
      if(fe<=idelta)  continue;
      ++leng;
      IWRK(wk,npk+1,leng) = ia[i];
      IWRK(wk,npk+2,leng) = nint(fe);
      FWRK(wk,npk+1,leng) = ia[i];
      FWRK(wk,npk+2,leng) = 1.e-4*fe;
    }
    if(*nadcgam >= cf->madcgam-2) return;
    if(!leng) continue;
    int32_t n1  = ++(*nadcgam);
    int32_t ic  = ipnpk[ipk];
    int32_t ix  = ia[ic] / cf->offset;
    int32_t iy  = ia[ic] - ix * cf->offset;
    int32_t itype = peak_type(cf,ix,iy);

    double  e1 = 0., x1 = 0., y1 = 0., e2 = 0., x2 = 0., y2 = 0., chisq = chisq2;
    gamma_cluster(cf,wk,itype,leng,&IWRK(wk,npk+1,0),&FWRK(wk,npk+2,0),&chisq,&e1,&x1,&y1,&e2,&x2,&y2);

    adcgam[n1].e    = e1;
    adcgam[n1].x[1] = x1;
    adcgam[n1].y[1] = y1;
    adcgam[n1].chi2 = chisq;
    adcgam[n1].size = leng;
    adcgam[n1].type = itype;
    adcgam[n1].id   = 10;
    adcgam[n1].stat = itype;

    if(e2>0. && n1 < cf->madcgam-3) {
      int32_t n2        = ++(*nadcgam);
      adcgam[n2].e    = e2;
      adcgam[n2].x[1] = x2;
      adcgam[n2].y[1] = y2;
      adcgam[n2].chi2 = chisq;
      adcgam[n2].size = leng;
      adcgam[n1].type = itype+10;
      adcgam[n2].type = itype+10;
      adcgam[n1].id   = 11;
      adcgam[n2].id   = 12;
      adcgam[n2].stat = itype;
      adcgam[n1].x[0] = 0.5*(x1-x2);
      adcgam[n1].y[0] = 0.5*(y1-y2);
      adcgam[n2].x[0] = 0.5*(x2-x1);
      adcgam[n2].y[0] = 0.5*(y2-y1);

      for(int32_t j = 1; j <= leng; ++j) {
        if(j>=_MAXLEN_) break;
        adcgam[n1].element[j] = IWRK(wk,npk+1,j);
        adcgam[n2].element[j] = IWRK(wk,npk+1,j);
        adcgam[n1].elfract[j] = nint(IWRK(wk,npk+2,j)*e1/(e1+e2));
        adcgam[n2].elfract[j] = nint(IWRK(wk,npk+2,j)*e2/(e1+e2));
      }
    } else {
      for(int32_t j = 1; j <= leng; ++j) {
        if(j>=_MAXLEN_) break;
        adcgam[n1].element[j] = IWRK(wk,npk+1,j);
        adcgam[n1].elfract[j] = IWRK(wk,npk+2,j);
      }
    }
  }

  return;
}

static int32_t peak_type(const ilreco_config *cf, int32_t ix, int32_t iy) {

  if(cf->cellmask) {
    /* mask authoritative: 1 = missing in-grid neighbor (hole border / rim),
       2 = bounding-box ring, 0 = fully surrounded */
    for(int32_t dx = -1; dx <= 1; ++dx)
      for(int32_t dy = -1; dy <= 1; ++dy) {
        if(!dx && !dy) continue;
        int32_t nx = ix + dx, ny = iy + dy;
        if(nx >= 1 && nx <= cf->ncol && ny >= 1 && ny <= cf->nrow &&
           !cell_exists(cf, nx, ny)) return 1;
      }
    if(ix == 1 || ix == cf->ncol || iy == 1 || iy == cf->nrow) return 2;
    return 0;
  }

  if(cf->hole_enabled) {
    if( (ix == cf->ncol/2-1 || ix == cf->ncol/2+2) &&
         iy >= cf->nrow/2-1 && iy <= cf->nrow/2+2) return 1;      //  hole

    if( (iy == cf->nrow/2-1 || iy == cf->nrow/2+2) &&
         ix >= cf->ncol/2-1 && ix <= cf->ncol/2+2) return 1;
  }

  if( ix == 1 || ix == cf->ncol || iy == 1 || iy == cf->nrow) return 2;   //  transition (outer boundary)

  return 0;         // inside matrix
}

static void gamma_cluster(const ilreco_config *cf, ilreco_workspace *wk,
                          int32_t itype, int32_t nadc, int32_t *ia, double *id,
                          double *chisq, double *e1, double *x1, double *y1,
                          double *e2, double *x2, double *y2)  {

  const double  xm2cut = 1.7e-3*cf->zcal;

  *e2 = 0.; *x2 = 0.; *y2 = 0.;
  int32_t nzero, *iazero = wk->iazero;
  fill_zero_hits(cf,nadc,ia,&nzero,iazero);          // make use of good but zero signal cells around clusters
  mom1_cluster(cf,nadc,ia,id,nzero,iazero,e1,x1,y1); // get initial e,x,y
  if(nadc <= 0) return;

  double chimem   = *chisq;   // memorize reference chi2 value
  double  chi0    = chisq1_cluser(cf,nadc,ia,id,nzero,iazero,*e1,*x1,*y1);  //  get initial chi2
  double  chisq0  = chi0;
  int32_t        ndof = nzero + nadc - 2;   //  number of degrees of freedom
  if(ndof<1) ndof = 1;

  *chisq = chi0/(double)ndof;
  double x0 = *x1, y0 = *y1;
  while(1) {        // iteration loop

    double const dxy = 0.05, stepmin = .002;
    double stepx, stepy;
    double chiright = chisq1_cluser(cf,nadc,ia,id,nzero,iazero,*e1,x0+dxy,y0);
    double chileft  = chisq1_cluser(cf,nadc,ia,id,nzero,iazero,*e1,x0-dxy,y0);
    double chiup    = chisq1_cluser(cf,nadc,ia,id,nzero,iazero,*e1,x0,y0+dxy);
    double chidown  = chisq1_cluser(cf,nadc,ia,id,nzero,iazero,*e1,x0,y0-dxy);

    if(chiright <= chi0 || chileft <= chi0) {
      stepx = dxy;
      if(chileft <= chiright) stepx *= -1.;
    } else {
      stepx = 0.;
      double  parx = 0.5*(chiright+chileft) - chi0;
      if(parx>0.) stepx = 0.25*dxy*(chileft-chiright)/parx;
    }
    if(chiup <= chi0 || chidown <= chi0) {
      stepy = dxy;
      if(chidown<=chiup) stepy *= -1.;
    } else {
      stepy = 0.;
      double pary = 0.5*(chiup+chidown) - chi0;
      if(pary>0.) stepy = 0.25*dxy*(chidown-chiup)/pary;
    }
    if(fabs(stepx)<stepmin && fabs(stepy)<stepmin) break;   // terminate: steps are too small already

    double chi00 = chisq1_cluser(cf,nadc,ia,id,nzero,iazero,*e1,x0+stepx,y0+stepy);
    if(chi0<chi00)  break;    // at min, terminate
    chi0 = chi00;
    x0  += stepx;
    y0  += stepy;
  }

  if(chi0<chisq0) { // if chi2 improved: store new values
    *x1 = x0;
    *y1 = y0;
    *chisq = chi0 / ndof;
  }

  if(*chisq <= chimem)  return;         //  no 2nr order separation of chi2 less than maximum single peak value

//  otherwise do 2nr order separation:

  double  chiold = *chisq, ee = 0., xx = 0., yy = 0.;
  tgamma_cluster(cf,nadc,ia,id,nzero,iazero,chisq,&ee,&xx,&yy,e2,x2,y2);

  if(*e2>0.)  {     //  if chi2 is improved, decide if the separation has physical meaning by inv. mass
    double  xm2 = sqrt(ee*(*e2))*hypot((xx-(*x2)),(yy-(*y2)));
    if(xm2>xm2cut)  {         //  if separation have physical meaning fix the parameters of first gamma
      *e1  = ee;
      *x1  = xx;
      *y1  = yy;
    } else {
      *e2  = 0.;              //  no physical meaning e2=0 -> second  gamma reset
      *chisq = chiold;
    }
  }

  return;
}

static void fill_zero_hits(const ilreco_config *cf, int32_t nadc, int32_t *ia,
                           int32_t *nzero, int32_t *iazero) {
  *nzero = 0;

  for(int32_t i = 1; i <= nadc; ++i) {
    int32_t ix = ia[i]/cf->offset;
    int32_t iy = ia[i] - ix * cf->offset;

    if(ix>1)  {
      if(cell_exists(cf,ix-1,iy))
        iazero[++(*nzero)] = iy + (ix-1) * cf->offset;                    //  left neib
      if(iy>1 && cell_exists(cf,ix-1,iy-1))
        iazero[++(*nzero)] = (iy-1) + (ix-1) * cf->offset;  //  bottom left neib
      if(iy<cf->nrow && cell_exists(cf,ix-1,iy+1))
        iazero[++(*nzero)] = (iy+1) + (ix-1) * cf->offset;  //  top left neib
    }
    if(ix<cf->ncol)  {
      if(cell_exists(cf,ix+1,iy))
        iazero[++(*nzero)] = iy + (ix+1) * cf->offset;                    //  right neib
      if(iy>1 && cell_exists(cf,ix+1,iy-1))
        iazero[++(*nzero)] = (iy-1) + (ix+1) * cf->offset;  //  bottom right neib
      if(iy<cf->nrow && cell_exists(cf,ix+1,iy+1))
        iazero[++(*nzero)] = (iy+1) + (ix+1) * cf->offset;  //  top right neib
    }
    if(iy>1 && cell_exists(cf,ix,iy-1))
      iazero[++(*nzero)] = iy-1 + ix * cf->offset;        //  bottom neib
    if(iy<cf->nrow && cell_exists(cf,ix,iy+1))
      iazero[++(*nzero)] = iy+1 + ix * cf->offset;        //  top neib
  }

  for(int32_t i = 1; i <= (*nzero); ++i)
    for(int32_t j = 1; j <= nadc; ++j)
      if(ia[j] == iazero[i]) iazero[i] = -1;

  for(int32_t i = 1; i <= (*nzero); ++i)
    if(iazero[i] != -1)
      for(int32_t j = i+1; j <= (*nzero); ++j)
        if(iazero[j]==iazero[i]) iazero[j] = -1;

  int32_t newzero = 0;
  for(int32_t i = 1; i <= (*nzero); ++i)
    if(iazero[i] != -1 /* && status_channel[?]==0 */)
      iazero[++newzero] = iazero[i];

  *nzero = newzero;
  return;
}

static void mom1_cluster(const ilreco_config *cf, int32_t nadc, int32_t *ia, double *id,
                         int32_t nzero, int32_t *iazero, double *e, double *x, double *y) {

  *e = *x = *y = 0.;
  for(int32_t i = 1; i <= nadc; ++i) {
    double a  = id[i];
    int32_t ix    = ia[i]/cf->offset;
    int32_t iy    = ia[i] - ix*cf->offset;
    *e  += a;
    *x  += a*ix;
    *y  += a*iy;
  }

  if(*e <= 0.) return;
  *x  /= *e;
  *y  /= *e;

  double corr = 0.;
  for(int32_t i = 1; i <= nadc; ++i) {
    double dx = (ia[i]/cf->offset)-(*x);
    double dy = (ia[i] - (ia[i]/cf->offset)*cf->offset)-(*y);
    corr  += profile_mean(cf,dx,dy);
  }
  for(int32_t i = 1; i <= nzero; ++i) {
    double dx = (iazero[i]/cf->offset)-(*x);
    double dy = (iazero[i] - (iazero[i]/cf->offset)*cf->offset)-(*y);
    corr  += profile_mean(cf,dx,dy);
  }

  corr /= 1.006;    // to be eliminated
  if(corr<0.8)  corr = 0.8;
  if(corr>1.0)  corr = 1.0;
  *e   /= corr;
  return;
}

static void mom2_cluster(const ilreco_config *cf, int32_t nadc, int32_t *ia, double *id,
                         int32_t nzero, int32_t *iazero, double *e, double *x, double *y,
                         double *xx, double *yy, double *xy) {
  *e = *x = *y = *xx = *yy = *xy = 0.;
  for(int32_t i = 1; i <= nadc; ++i) {
    double a  = id[i];
    int32_t ix    = ia[i]/cf->offset;
    int32_t iy    = ia[i] - ix*cf->offset;
    *e  += a;
    *x  += a*ix;
    *y  += a*iy;
  }

  if(*e <= 0.) return;
  *x  /= *e;
  *y  /= *e;

  for(int32_t i = 1; i <= nadc; ++i) {
    double a  = id[i];
    int32_t ix    = ia[i]/cf->offset;
    int32_t iy    = ia[i] - ix*cf->offset;
    *xx += a/(*e)*(ix-(*x))*(ix-(*x));
    *yy += a/(*e)*(iy-(*y))*(iy-(*y));
    *xy += a/(*e)*(ix-(*x))*(iy-(*y));
  }

  double corr = 0.;
  for(int32_t i = 1; i <= nadc; ++i) {
    double dx = (ia[i]/cf->offset)-(*x);
    double dy = (ia[i] - (ia[i]/cf->offset)*cf->offset)-(*y);
    corr  += profile_mean(cf,dx,dy);
  }
  for(int32_t i = 1; i <= nzero; ++i) {
    double dx = (iazero[i]/cf->offset)-(*x);
    double dy = (iazero[i] - (iazero[i]/cf->offset)*cf->offset)-(*y);
    corr  += profile_mean(cf,dx,dy);
  }

  corr /= 1.006;    // to be eliminated
  if(corr<0.8)  corr = 0.8;
  if(corr>1.0)  corr = 1.0;
  *e   /= corr;
  return;
}

static double chisq1_cluser(const ilreco_config *cf, int32_t nadc, int32_t *ia, double *id,
                            int32_t nzero, int32_t *iazero, double e, double x, double y)  {

  double chi2 = 0.;

  for(int32_t i = 1; i <= nadc; ++i) {
    int32_t ix    = ia[i]/cf->offset;
    int32_t iy    = ia[i] - ix*cf->offset;
    if(e) {
      if(fabs(x-ix)<=6.0 && fabs(y-iy)<=6.0)  {
        double f = (profile_mean(cf,x-ix,y-iy) - id[i]/e);
        chi2 += 1.e4*e*f*f/sigma2(cf,x-ix,y-iy,e);
      }
    } else {
      chi2 += id[i]*id[i] / 9.;
      printf(" case 0 ch\n"); //  should'nt come here
    }
  }

  for(int32_t i = 1; i <= nzero; ++i) {
    int32_t ix    = iazero[i]/cf->offset;
    int32_t iy    = iazero[i] - ix*cf->offset;
    if(e) {
      double f = profile_mean(cf,x-ix,y-iy);
      chi2 += 1.e4*e*f*f/sigma2(cf,x-ix,y-iy,e);
    } else {
      printf(" case 0 ch\n"); //  should'nt come here
    }
  }

  return  chi2;
}

static double sigma2(const ilreco_config *cf, double dx, double dy, double e) {

  if(dx*dx+dy*dy>25.) return 1.e2;

  double const alp = 0.816, bet1 = 32.1, bet2 = 17.2;

  double retval = 1.e2*(alp*profile_mean(cf,dx,dy) + (bet1+bet2*sqrt(e))*d2c(cf,dx,dy) + 2.e-3/e);

  return retval*pow(e,-0.166);
}

static double d2c(const ilreco_config *cf, double x, double y) {
  double ax = fabs(x*1.e2), ay = fabs(y*1.e2);
  int32_t i = (int32_t)ax, j = (int32_t)ay;
  double wx = ax-i, wy = ay-j;
  if(i<NPROF && j<NPROF)
    return AD2C(cf,i,j)     * (1.-wx) * (1.-wy) +
           AD2C(cf,i+1,j)   *     wx  * (1.-wy) +
           AD2C(cf,i,  j+1) * (1.-wx) *     wy  +
           AD2C(cf,i+1,j+1) *     wx  *     wy;
  return 1.;
}

/* NOTE: dead code path — tgamma_cluster below is an empty stub (as in the
 * original), so the second-order separation never runs; this helper also
 * carries the original's suspected typo (e1*f2 instead of e1*f1). Preserved
 * untouched; see tests/KNOWN_ISSUES.md. */
static double chisq2t_hyc(const ilreco_config *cf,
                          double ecell, double e1, double dx1, double dy1,
                          double e2, double dx2, double dy2, double f1, double f2) {
  double s;
  int32_t ic = !(!(e1))*10+!(!(e2));
  switch(ic) {
    case 11:
      s = e1*sigma2(cf,dx1,dy1,e1) + e2*sigma2(cf,dx2,dy2,e2);
      break;
    case 10:
      s = e1*sigma2(cf,dx1,dy1,e1);
      break;
    case  1:
      s = e2*sigma2(cf,dx2,dy2,e2);
      break;
    default:
      s = 9e-4;
  }
    s *= 1.e8;
    double d = e1*f2+e2*f2-ecell;
    return d*d/s;
}

static void tgamma_cluster(const ilreco_config *cf, int32_t nadc, int32_t *ia, double *id,
                           int32_t nzero, int32_t *iazero, double *chisq,
                           double *ee, double *xx, double *yy,
                           double *e2, double *x2, double *y2) {
  (void)cf; (void)nadc; (void)ia; (void)id; (void)nzero; (void)iazero;
  (void)chisq; (void)ee; (void)xx; (void)yy; (void)e2; (void)x2; (void)y2;
  (void)chisq2t_hyc; (void)mom2_cluster;
}

static double profile_mean(const ilreco_config *cf, double x, double y) {
  double ax = fabs(x*1.e2), ay = fabs(y*1.e2);
  int32_t i = (int32_t)ax, j = (int32_t)ay;
  double wx = ax-i, wy = ay-j;

  if(i<NPROF && j<NPROF)
    return AMEAN(cf,i,j)     * (1.-wx) * (1.-wy) +
           AMEAN(cf,i+1,j)   *     wx  * (1.-wy) +
           AMEAN(cf,i,j+1)   * (1.-wx) *     wy  +
           AMEAN(cf,i+1,j+1) *     wx  *     wy;
  return 0.;
}

/* round to nearest integer, halves away from zero */
static int32_t nint(double x) {
  int32_t n = (x<0.) ? -(int32_t)(0.5-x) : (int32_t)(0.5+x);
  return n;
}
