//
//  ilreco — island clustering for square-cell calorimeters.
//
//  2026 modernization: all mutable state moved into explicit contexts
//  (ilreco_config: immutable shared geometry/profile/tunables, allocated once;
//  ilreco_workspace: per-thread buffers, allocated once). The reconstruction
//  algorithm itself is UNCHANGED — every expression and its evaluation order
//  is preserved from the original (validated bit-for-bit by the golden test
//  tables in tests/data). The legacy 1-based macro-configured API is kept as
//  a thin shim over an internal default context.
//
//  gcc ilreco.c -O3 -Wall -std=c11 -lm
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "ilreco.h"

/* ============================ contexts ==================================== */

#define NPROF _N_PROFILE_POINTS_

struct ilreco_config {
    int ncol, nrow;        /* grid size (columns x rows)                      */
    int offset;            /* column stride in packed 1-based addresses       */
    int madr;              /* hit/address buffer capacity                     */
    int madr0;             /* nominal module count (peak-count guard)         */
    int mcl;               /* max raw clusters + 1                            */
    int madcgam;           /* max reconstructed objects + 1                   */
    double zcal;           /* target distance (two-gamma mass separation cut) */
    double min_seed_gev;   /* cluster seed threshold (historical 0.01)        */
    double seed_scale;     /* seed-threshold energy scaling (_IBVAL_)         */
    int hole_enabled;      /* classify central 2x2 beam-hole cells (HyCal)    */
    unsigned char *cellmask; /* NULL = all cells exist; else [nrow*ncol], 0=missing */
    double *amean;         /* (NPROF+1)^2 shower profile mean                 */
    double *ad2c;          /* (NPROF+1)^2 shower profile second moment        */
};

struct ilreco_workspace {
    int stride;            /* == cfg->madr at creation                        */
    int mcl, madcgam;
    /* event buffers (new API entry) */
    int *ia; double *id; int *lencl; adcgam_t *adcgam;
    /* cluster_search scratch */
    int *iwork; double *dwork;
    /* process_cluster / gamma_cluster scratch */
    int *iazero;           /* 8*stride + 2                                    */
    int *iwrk, *idp;       /* (_MPK_+1) * stride each                         */
    double *fwrk;          /* (_MPK_+1) * stride                              */
    void *arena;           /* single allocation backing all of the above      */
};

#define AMEAN(cf, i, j) ((cf)->amean[(size_t)(i) * (NPROF + 1) + (j)])
#define AD2C(cf, i, j)  ((cf)->ad2c [(size_t)(i) * (NPROF + 1) + (j)])
#define IWRK(wk, p, i)  ((wk)->iwrk[(size_t)(p) * (wk)->stride + (i)])
#define IDP(wk, p, i)   ((wk)->idp [(size_t)(p) * (wk)->stride + (i)])
#define FWRK(wk, p, i)  ((wk)->fwrk[(size_t)(p) * (wk)->stride + (i)])

/* ===================== internal declarations ============================== */

static int  cluster_search_core(const ilreco_config *cf, ilreco_workspace *wk,
                                int nw, int *ia, double *id, int *lencl);
static void process_cluster_core(const ilreco_config *cf, ilreco_workspace *wk,
                                 int nadc, int *ia, double *id,
                                 int *nadcgam, adcgam_t *adcgam);
static void order_hits(int nw, int *ia, double *id);
static int  peak_type(const ilreco_config *cf, int ix, int iy);
int   nint(double x);
static void gamma_cluster(const ilreco_config *cf, ilreco_workspace *wk,
                          int itype, int nadc, int *ia, double *id,
                          double *chisq, double *e1, double *x1, double *y1,
                          double *e2, double *x2, double *y2);
static void fill_zero_hits(const ilreco_config *cf, int nadc, int *ia,
                           int *nzero, int *iazero);
static void mom1_cluster(const ilreco_config *cf, int nadc, int *ia, double *id,
                         int nzero, int *iazero, double *e, double *x, double *y);
static void mom2_cluster(const ilreco_config *cf, int nadc, int *ia, double *id,
                         int nzero, int *iazero, double *e, double *x, double *y,
                         double *xx, double *yy, double *xy);
static double chisq1_cluser(const ilreco_config *cf, int nadc, int *ia, double *id,
                            int nzero, int *iazero, double e, double x, double y);
static void tgamma_cluster(const ilreco_config *cf, int nadc, int *ia, double *id,
                           int nzero, int *iazero, double *chisq,
                           double *ee, double *xx, double *yy,
                           double *e2, double *x2, double *y2);
static double profile_mean(const ilreco_config *cf, double x, double y);
static double sigma2(const ilreco_config *cf, double dx, double dy, double e);
static double d2c(const ilreco_config *cf, double x, double y);
static void add_single_adcgam(const ilreco_config *cf, ilreco_workspace *wk,
                              int nadc, int *ia, double *id, int *ipnpk,
                              int *nadcgam, adcgam_t *adcgam);
static void add_many_adcgam(const ilreco_config *cf, ilreco_workspace *wk,
                            int nadc, int *ia, double *id, int npk, int *ipnpk,
                            int *nadcgam, adcgam_t *adcgam);

/* random generator functions for the built-in test event generator */
double  ZBQLINI(int seed, double ZBQLIX[43+1]);
double  ZBQLUAB(double x1, double x2);
double  ZBQLU01();

static int debug_ncalls = 0;

/* ====================== config / workspace =============================== */

static int load_profile(ilreco_config *cf, const char *path,
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

static ilreco_config *config_create_internal(int n_cols, int n_rows, int offset,
                                             int madr, int madr0, const char *profile_path,
                                             char *errbuf, size_t errlen) {
    if (n_cols < 1 || n_rows < 1 || offset <= n_rows + 1) {
        if (errbuf) snprintf(errbuf, errlen, "invalid geometry %dx%d (offset %d)",
                             n_cols, n_rows, offset);
        return NULL;
    }
    ilreco_config *cf = calloc(1, sizeof *cf);
    if (!cf) return NULL;
    cf->ncol = n_cols;
    cf->nrow = n_rows;
    cf->offset = offset;
    cf->madr = madr;
    cf->madr0 = madr0;
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
    if (profile_path && load_profile(cf, profile_path, errbuf, errlen) != 0) {
        ilreco_config_destroy(cf);
        return NULL;
    }
    return cf;
}

ilreco_config *ilreco_config_create(int n_cols, int n_rows,
                                    const char *profile_path,
                                    char *errbuf, size_t errbuf_len) {
    const int offset = n_rows + 2;
    const int madr = (n_cols + 1) * offset + n_rows + 2 + 64;
    return config_create_internal(n_cols, n_rows, offset, madr, madr - 64,
                                  profile_path, errbuf, errbuf_len);
}

void ilreco_config_destroy(ilreco_config *cfg) {
    if (!cfg) return;
    free(cfg->cellmask);
    free(cfg->amean);
    free(cfg->ad2c);
    free(cfg);
}

/* 1-based column/row -> does the cell exist? (mask absent => yes) */
static int cell_exists(const ilreco_config *cf, int ix, int iy) {
    if (!cf->cellmask) return 1;
    return cf->cellmask[(size_t)(iy - 1) * cf->ncol + (ix - 1)] != 0;
}

int ilreco_config_set_cell_mask(ilreco_config *cfg, const unsigned char *mask) {
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
void ilreco_config_set_hole_classification(ilreco_config *cfg, int enabled) {
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
                       + (size_t)(_MPK_ + 1) * m /*idp*/;
    const size_t n_dbl = m /*id*/ + m /*dwork*/ + (size_t)(_MPK_ + 1) * m /*fwrk*/;
    const size_t bytes = n_dbl * sizeof(double) + n_int * sizeof(int)
                       + (size_t)cfg->madcgam * sizeof(adcgam_t);
    char *p = malloc(bytes);
    if (!p) { free(wk); return NULL; }
    wk->arena = p;
    /* doubles first (strictest alignment), then adcgam_t, then ints */
    wk->id = (double *)p;                 p += m * sizeof(double);
    wk->dwork = (double *)p;              p += m * sizeof(double);
    wk->fwrk = (double *)p;               p += (size_t)(_MPK_ + 1) * m * sizeof(double);
    wk->adcgam = (adcgam_t *)p;           p += (size_t)cfg->madcgam * sizeof(adcgam_t);
    wk->ia = (int *)p;                    p += m * sizeof(int);
    wk->iwork = (int *)p;                 p += m * sizeof(int);
    wk->lencl = (int *)p;                 p += (size_t)cfg->mcl * sizeof(int);
    wk->iazero = (int *)p;                p += (8 * m + 2) * sizeof(int);
    wk->iwrk = (int *)p;                  p += (size_t)(_MPK_ + 1) * m * sizeof(int);
    wk->idp = (int *)p;
    return wk;
}

void ilreco_workspace_destroy(ilreco_workspace *ws) {
    if (!ws) return;
    free(ws->arena);
    free(ws);
}

int ilreco_reconstruct(const ilreco_config *cfg, ilreco_workspace *ws,
                       const ilreco_hit *hits, int n_hits,
                       ilreco_cluster *out, int max_out) {
    if (!cfg || !ws || n_hits < 0 || (n_hits > 0 && !hits) ||
        ws->stride < cfg->madr || n_hits > cfg->madr - 2)
        return -1;
    for (int i = 0; i < n_hits; ++i) {
        if (hits[i].col < 0 || hits[i].col >= cfg->ncol ||
            hits[i].row < 0 || hits[i].row >= cfg->nrow)
            return -1;
        if (!cell_exists(cfg, hits[i].col + 1, hits[i].row + 1))
            return -1;
        ws->ia[i + 1] = (hits[i].col + 1) * cfg->offset + (hits[i].row + 1);
        ws->id[i + 1] = hits[i].e;
    }

    int nadcgam = 0;
    const int ncl = cluster_search_core(cfg, ws, n_hits, ws->ia, ws->id, ws->lencl);
    for (int icl = 1, ipncl = 1; icl <= ncl && nadcgam < cfg->madcgam - 2;
         ipncl += ws->lencl[icl++])
        process_cluster_core(cfg, ws, ws->lencl[icl], &ws->ia[ipncl - 1],
                             &ws->id[ipncl - 1], &nadcgam, ws->adcgam);

    /* energy-descending, stable */
    int order[_MADCGAM_];
    for (int g = 1; g <= nadcgam; ++g) order[g - 1] = g;
    for (int a = 1; a < nadcgam; ++a) {
        const int key = order[a];
        int b = a - 1;
        while (b >= 0 && ws->adcgam[order[b]].e < ws->adcgam[key].e) {
            order[b + 1] = order[b];
            --b;
        }
        order[b + 1] = key;
    }
    const int n_out = nadcgam < max_out ? nadcgam : max_out;
    for (int k = 0; k < n_out; ++k) {
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

/* ================== legacy API (default context shims) =================== */

static ilreco_config *g_legacy_cfg = NULL;
static ilreco_workspace *g_legacy_wk = NULL;

static int ensure_legacy(const char *profile_path_or_null) {
    if (g_legacy_cfg && profile_path_or_null) {           /* reload */
        ilreco_workspace_destroy(g_legacy_wk);
        ilreco_config_destroy(g_legacy_cfg);
        g_legacy_cfg = NULL;
        g_legacy_wk = NULL;
    }
    if (!g_legacy_cfg) {
        char err[256];
        g_legacy_cfg = config_create_internal(_NCOL_, _NROW_, _OFFSET_, _MADR_, _MADR0_,
                                              profile_path_or_null, err, sizeof err);
        if (!g_legacy_cfg) {
            /* legacy behavior: report and terminate on unreadable profile */
            fprintf(stderr, "(!)ERROR(!) %s\n", err);
            exit(1);
        }
    }
    if (!g_legacy_wk) g_legacy_wk = ilreco_workspace_create(g_legacy_cfg);
    return g_legacy_cfg && g_legacy_wk;
}

void read_profile_data(const char *prof_file_name) {
    if (access(prof_file_name, R_OK) == -1)
        fprintf(stderr, "(!)ERROR(!) profile file '%s' does not exists or access is denied\n",
                prof_file_name);
    ensure_legacy(prof_file_name);
}

int cluster_search(int nw, int *ia, double *id, int *lencl) {
    ensure_legacy(NULL);
    return cluster_search_core(g_legacy_cfg, g_legacy_wk, nw, ia, id, lencl);
}

void process_cluster(int nadc, int *ia, double *id, int *nadcgam, adcgam_t *adcgam) {
    ensure_legacy(NULL);
    process_cluster_core(g_legacy_cfg, g_legacy_wk, nadc, ia, id, nadcgam, adcgam);
}

/* ======================= algorithm (unchanged) ============================ */

static int cluster_search_core(const ilreco_config *cf, ilreco_workspace *wk,
                               int nw, int *ia, double *id, int *lencl) {
  if(nw<2) {
    lencl[1] = 1;
    return nw;
  }
  order_hits(nw,ia,id);       // addresses must be in increasing order

  int *iwork = wk->iwork;  double *dwork = wk->dwork;
  int ncl = 0, next  = 1, iak = 0;
  for(int k = 2; k <= nw+1; ++k) {
    if(k<=nw) iak = ia[k];
    if(iak-ia[k-1]<=1 && k<=nw) continue;

    int ib  = next; //  first word of the (sub)cluster
    int ie  = k-1;  //  last  word of the (sub)cluster
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
    int ias     = ia[ib];
    int iaf     = ia[ie];
    int last    = ib-1;
    int lastcl  = ncl-1;

    for(int icl = lastcl; icl>0; --icl) {
      int leng  = lencl[icl];
      if(ias-ia[last] > cf->offset) break;  //  no subclusters to glue
      for(int i = last; i > last-leng; --i) {
        if(ias-ia[i]  >  cf->offset) break;
        if(iaf-ia[i]  >= cf->offset) {      //  subclusters to glue
          if(icl<ncl-1 && leng<cf->madr) {

            memmove(&iwork[1],&ia[last+1-leng],leng*sizeof(int));
            memmove(&ia[last+1-leng],&ia[last+1],(ib-1-last)*sizeof(int));
            memmove(&ia[ib-leng],&iwork[1],leng*sizeof(int));
            memcpy(&dwork[1],&id[last+1-leng],leng*sizeof(double));
            memmove(&id[last+1-leng],&id[last+1],(ib-1-last)*sizeof(double));
            memcpy(&id[ib-leng],&dwork[1],leng*sizeof(double));

            for(int j = icl; j < ncl-1; ++j) lencl[j] = lencl[j+1];
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

static void order_hits(int nw, int *ia, double *id) {

  for(int k = 2; k <= nw; ++k) {
    if(ia[k] > ia[k-1])  continue;
    int     iat = ia[k];
    double  idt = id[k];
    for(int i = k-1; i>=0; --i)
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
                                 int nadc, int *ia, double *id,
                                 int *nadcgam, adcgam_t *adcgam) {
  order_hits(nadc,ia,id);     // addresses must be in increasing order

  double minpk = cf->min_seed_gev;
  int ipnpk[_MPK_];

  if(nadc>=3) {
    double idsum = 0.;
    for(int i = 1; i <= nadc; ++i) idsum += id[i];
    double ib = cf->seed_scale*log(1.+idsum);
    if(ib>1) minpk *= ib;
    minpk = nint(1.e2*minpk);
    minpk *= 1.e2;
    minpk = 1.e-4*minpk;
  }

  int npk = 0;
  for(int ic = 1; ic <= nadc; ++ic) {
    double iac = id[ic];
    if(iac<minpk) continue;
    int ixy     = ia[ic];
    int ixymax  = ixy + cf->offset + 1;
    int ixymin  = ixy - cf->offset - 1;
    int iyc     = ixy - (ixy/cf->offset)*cf->offset;

    int in;
    in  = ic  + 1;
    int skipflag = 0;
    while(in<=nadc && ia[in] <= ixymax) {
      int iy  = ia[in] - (ia[in]/cf->offset)*cf->offset;
      if(abs(iy-iyc)<=1 && id[in]>=iac) {skipflag = 1; break;}
      ++in;
    }
    if(skipflag)  continue;
    in = ic - 1;
    while (in>=1 && ia[in] >= ixymin) {
      int iy  = ia[in] - (ia[in]/cf->offset)*cf->offset;
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
                              int nadc, int *ia, double *id, int *ipnpk,
                              int *nadcgam, adcgam_t *adcgam) {

  const double  chisq1  = 90.; // 3.0 value of chi2 for preliminary seperation

  int n1  = ++(*nadcgam);
  int ic  = ipnpk[1];
  int ix  = ia[ic]/cf->offset;
  int iy  = ia[ic]-ix*cf->offset;
  int itype = peak_type(cf,ix,iy);

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
    int n2        = ++(*nadcgam);
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

    for(int j = 1; j <= nadc; ++j) {
      if(j>=_MAXLEN_) break;
      adcgam[n1].element[j] = ia[j];
      adcgam[n2].element[j] = ia[j];
      adcgam[n1].elfract[j] = id[j]*e1/(e1+e2);
      adcgam[n2].elfract[j] = id[j]*e2/(e1+e2);
    }
  } else {
    for(int j = 1; j <= nadc; ++j) {
      if(j>=_MAXLEN_) break;
      adcgam[n1].element[j] = ia[j];
      adcgam[n1].elfract[j] = nint(1.e4*id[j]);
    }
  }
  return;
}

static void add_many_adcgam(const ilreco_config *cf, ilreco_workspace *wk,
                            int nadc, int *ia, double *id, int npk, int *ipnpk,
                            int *nadcgam, adcgam_t *adcgam) {

  const double  idelta  =  0.; //  min cell energy part to be a member of separated cluster
  const double  chisq1  = 90.; // 3.0 value of chi2 for preliminary seperation
  const double  chisq2  = 50.; // 0.8 value of chi2 for final seperation
  const int     niter   =   6;
  int ngam0 = *nadcgam;
  int igmpk[_MPK_][3];

  double ratio = 1.;
  for(int iter = 1; iter <= niter; ++iter) {
    for(int i = 1; i <= nadc; ++i) {IWRK(wk,0,i) = 0; FWRK(wk,0,i) = 0.;}
    double epk[_MPK_], xpk[_MPK_], ypk[_MPK_];
    for(int ipk = 1; ipk <= npk; ++ipk) {
      int ic = ipnpk[ipk];
      if(iter!=1) ratio = FWRK(wk,ipk,ic) / FWRK(wk,npk+1,ic);
      double eg = id[ic] * ratio;
      int ixypk = ia[ic];
      int ixpk  = ixypk / cf->offset;
      int iypk  = ixypk - ixpk * cf->offset;
      epk[ipk]  = eg;
      xpk[ipk]  = eg*ixpk;
      ypk[ipk]  = eg*iypk;
      for(int in = ic + 1; in <= nadc; ++in) {
        int ixy = ia[in];
        int ix  = ixy / cf->offset;
        int iy  = ixy - ix * cf->offset;
        if(ixy-ixypk>cf->offset+1) break;
        if(abs(iy-iypk)>1) continue;
        if(iter!=1) ratio = FWRK(wk,ipk,in) / FWRK(wk,npk+1,in);
        eg = id[in] * ratio;
        epk[ipk] += eg;
        xpk[ipk] += eg*ix;
        ypk[ipk] += eg*iy;
      }

      for(int in = ic - 1; in>0; --in) {
        int ixy = ia[in];
        int ix  = ixy / cf->offset;
        int iy  = ixy - ix * cf->offset;
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

      for(int i = 1; i <= nadc; ++i) {
        int ixy = ia[i];
        int ix  = ixy / cf->offset;
        int iy  = ixy - ix * cf->offset;
        double dx = fabs(ix-xpk[ipk]);
        double dy = fabs(iy-ypk[ipk]);
        double a  = epk[ipk]*profile_mean(cf,dx,dy);
        IWRK(wk,ipk,i) = nint(a*1.e4);
        IWRK(wk,0,i)  += IWRK(wk,ipk,i);
        FWRK(wk,ipk,i) = a;
        FWRK(wk,0,i)  += FWRK(wk,ipk,i);
      }
    }

    for(int i = 1; i <= nadc; ++i) {
      int iwk = IWRK(wk,0,i);
      if(iwk<1) iwk = 1;
      IWRK(wk,npk+1,i) = iwk;
      if(FWRK(wk,0,i)>1.e-6)
        FWRK(wk,npk+1,i) = FWRK(wk,0,i);
      else
        FWRK(wk,npk+1,i) = 1.e-6;
    }
  }

  for(int ipk = 1; ipk <= npk; ++ipk) {
    int leng = 0;
    for(int i = 1; i <= nadc; ++i) {
      if(FWRK(wk,0,i) <= 1.e-6) continue;
      int ixy   = ia[i];
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
    int n1  = ++(*nadcgam);
    int ic  = ipnpk[ipk];
    int ix  = ia[ic] / cf->offset;
    int iy  = ia[ic] - ix * cf->offset;
    int itype = peak_type(cf,ix,iy);

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
      int n2        = ++(*nadcgam);
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

  for(int i = 1; i <= nadc; ++i) {IWRK(wk,0,i) = 0; IDP(wk,0,i) = 0; FWRK(wk,0,i) = 0.;}

  for(int ipk = 1; ipk <= npk; ++ipk)
    for(int i = 1; i <= nadc; ++i) {
      IWRK(wk,ipk,i) = 0; IDP(wk,ipk,i) = 0; FWRK(wk,ipk,i) = 0.;
      if(igmpk[ipk][2])
        for(int ig = igmpk[ipk][1]; ig <= igmpk[ipk][2]; ++ig) {

          int ixy = ia[i];
          int ix  = ixy / cf->offset;
          int iy  = ixy - ix * cf->offset;
          double dx = ix-adcgam[ig].x[1];
          double dy = iy-adcgam[ig].y[1];
          double a  = adcgam[ig].e*profile_mean(cf,dx,dy);
          int   iia = nint(a*1.e4);
          IWRK(wk,ipk,i) += iia;
          IWRK(wk,0,i)   += iia;
          FWRK(wk,ipk,i) += a;
          FWRK(wk,0,i)   += a;
          IDP(wk,ipk,i)  += iia;
        }
    }

  for(int i = 1; i <= nadc; ++i) {
    IDP(wk,0,i) = 0;
    for(int ipk = 1; ipk <= npk; ++ipk)
      IDP(wk,0,i) += IDP(wk,ipk,i);
    int ide = nint(id[i]*1.e4) - IDP(wk,0,i);
    if(!ide || FWRK(wk,0,i) == 0.) continue;

    double fw[_MPK_];
    for(int ipk = 1; ipk <= npk; ++ipk)
      fw[ipk] = FWRK(wk,ipk,i)/FWRK(wk,0,i);

    int idecorr = 0;
    for(int ipk = 1; ipk <= npk; ++ipk) {
      double fia = ide*fw[ipk];
      if(FWRK(wk,ipk,i)+fia>0.) {
        FWRK(wk,ipk,i) += fia;
        FWRK(wk,0,i)   += fia;
      }
      int iia = nint(fia);
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
  for(int ipk = 1; ipk <= npk; ++ipk) {
    int leng = 0;
    for(int i = 1; i <= nadc; ++i) {
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
    int n1  = ++(*nadcgam);
    int ic  = ipnpk[ipk];
    int ix  = ia[ic] / cf->offset;
    int iy  = ia[ic] - ix * cf->offset;
    int itype = peak_type(cf,ix,iy);

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
      int n2        = ++(*nadcgam);
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

      for(int j = 1; j <= leng; ++j) {
        if(j>=_MAXLEN_) break;
        adcgam[n1].element[j] = IWRK(wk,npk+1,j);
        adcgam[n2].element[j] = IWRK(wk,npk+1,j);
        adcgam[n1].elfract[j] = nint(IWRK(wk,npk+2,j)*e1/(e1+e2));
        adcgam[n2].elfract[j] = nint(IWRK(wk,npk+2,j)*e2/(e1+e2));
      }
    } else {
      for(int j = 1; j <= leng; ++j) {
        if(j>=_MAXLEN_) break;
        adcgam[n1].element[j] = IWRK(wk,npk+1,j);
        adcgam[n1].elfract[j] = IWRK(wk,npk+2,j);
      }
    }
  }

  return;
}

static int peak_type(const ilreco_config *cf, int ix, int iy) {

  if(cf->cellmask) {
    /* mask authoritative: 1 = missing in-grid neighbor (hole border / rim),
       2 = bounding-box ring, 0 = fully surrounded */
    for(int dx = -1; dx <= 1; ++dx)
      for(int dy = -1; dy <= 1; ++dy) {
        if(!dx && !dy) continue;
        int nx = ix + dx, ny = iy + dy;
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
                          int itype, int nadc, int *ia, double *id,
                          double *chisq, double *e1, double *x1, double *y1,
                          double *e2, double *x2, double *y2)  {

  const double  xm2cut = 1.7e-3*cf->zcal;

  *e2 = 0.; *x2 = 0.; *y2 = 0.;
  int nzero, *iazero = wk->iazero;
  fill_zero_hits(cf,nadc,ia,&nzero,iazero);          // make use of good but zero signal cells around clusters
  mom1_cluster(cf,nadc,ia,id,nzero,iazero,e1,x1,y1); // get initial e,x,y
  if(nadc <= 0) return;

  double chimem   = *chisq;   // memorize reference chi2 value
  double  chi0    = chisq1_cluser(cf,nadc,ia,id,nzero,iazero,*e1,*x1,*y1);  //  get initial chi2
  double  chisq0  = chi0;
  int        ndof = nzero + nadc - 2;   //  number of degrees of freedom
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

static void fill_zero_hits(const ilreco_config *cf, int nadc, int *ia,
                           int *nzero, int *iazero) {
  *nzero = 0;

  for(int i = 1; i <= nadc; ++i) {
    int ix = ia[i]/cf->offset;
    int iy = ia[i] - ix * cf->offset;

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

  for(int i = 1; i <= (*nzero); ++i)
    for(int j = 1; j <= nadc; ++j)
      if(ia[j] == iazero[i]) iazero[i] = -1;

  for(int i = 1; i <= (*nzero); ++i)
    if(iazero[i] != -1)
      for(int j = i+1; j <= (*nzero); ++j)
        if(iazero[j]==iazero[i]) iazero[j] = -1;

  int newzero = 0;
  for(int i = 1; i <= (*nzero); ++i)
    if(iazero[i] != -1 /* && status_channel[?]==0 */)
      iazero[++newzero] = iazero[i];

  *nzero = newzero;
  return;
}

static void mom1_cluster(const ilreco_config *cf, int nadc, int *ia, double *id,
                         int nzero, int *iazero, double *e, double *x, double *y) {

  *e = *x = *y = 0.;
  for(int i = 1; i <= nadc; ++i) {
    double a  = id[i];
    int ix    = ia[i]/cf->offset;
    int iy    = ia[i] - ix*cf->offset;
    *e  += a;
    *x  += a*ix;
    *y  += a*iy;
  }

  if(*e <= 0.) return;
  *x  /= *e;
  *y  /= *e;

  double corr = 0.;
  for(int i = 1; i <= nadc; ++i) {
    double dx = (ia[i]/cf->offset)-(*x);
    double dy = (ia[i] - (ia[i]/cf->offset)*cf->offset)-(*y);
    corr  += profile_mean(cf,dx,dy);
  }
  for(int i = 1; i <= nzero; ++i) {
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

static void mom2_cluster(const ilreco_config *cf, int nadc, int *ia, double *id,
                         int nzero, int *iazero, double *e, double *x, double *y,
                         double *xx, double *yy, double *xy) {
  *e = *x = *y = *xx = *yy = *xy = 0.;
  for(int i = 1; i <= nadc; ++i) {
    double a  = id[i];
    int ix    = ia[i]/cf->offset;
    int iy    = ia[i] - ix*cf->offset;
    *e  += a;
    *x  += a*ix;
    *y  += a*iy;
  }

  if(*e <= 0.) return;
  *x  /= *e;
  *y  /= *e;

  for(int i = 1; i <= nadc; ++i) {
    double a  = id[i];
    int ix    = ia[i]/cf->offset;
    int iy    = ia[i] - ix*cf->offset;
    *xx += a/(*e)*(ix-(*x))*(ix-(*x));
    *yy += a/(*e)*(iy-(*y))*(iy-(*y));
    *xy += a/(*e)*(ix-(*x))*(iy-(*y));
  }

  double corr = 0.;
  for(int i = 1; i <= nadc; ++i) {
    double dx = (ia[i]/cf->offset)-(*x);
    double dy = (ia[i] - (ia[i]/cf->offset)*cf->offset)-(*y);
    corr  += profile_mean(cf,dx,dy);
  }
  for(int i = 1; i <= nzero; ++i) {
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

static double chisq1_cluser(const ilreco_config *cf, int nadc, int *ia, double *id,
                            int nzero, int *iazero, double e, double x, double y)  {

  double chi2 = 0.;

  for(int i = 1; i <= nadc; ++i) {
    int ix    = ia[i]/cf->offset;
    int iy    = ia[i] - ix*cf->offset;
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

  for(int i = 1; i <= nzero; ++i) {
    int ix    = iazero[i]/cf->offset;
    int iy    = iazero[i] - ix*cf->offset;
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
  int i = (int)ax, j = (int)ay;
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
  int ic = !(!(e1))*10+!(!(e2));
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

static void tgamma_cluster(const ilreco_config *cf, int nadc, int *ia, double *id,
                           int nzero, int *iazero, double *chisq,
                           double *ee, double *xx, double *yy,
                           double *e2, double *x2, double *y2) {
  (void)cf; (void)nadc; (void)ia; (void)id; (void)nzero; (void)iazero;
  (void)chisq; (void)ee; (void)xx; (void)yy; (void)e2; (void)x2; (void)y2;
  (void)chisq2t_hyc; (void)mom2_cluster;
}

static double profile_mean(const ilreco_config *cf, double x, double y) {
  double ax = fabs(x*1.e2), ay = fabs(y*1.e2);
  int i = (int)ax, j = (int)ay;
  double wx = ax-i, wy = ay-j;

  if(i<NPROF && j<NPROF)
    return AMEAN(cf,i,j)     * (1.-wx) * (1.-wy) +
           AMEAN(cf,i+1,j)   *     wx  * (1.-wy) +
           AMEAN(cf,i,j+1)   * (1.-wx) *     wy  +
           AMEAN(cf,i+1,j+1) *     wx  *     wy;
  return 0.;
}

/* ================= legacy test scaffolding (unchanged) ==================== */

int read_event(int *nw, int *ia, double *id) {

  double  ech[_NCOL_+1][_NROW_+1];
  static int ievent = 0;
  ++ievent;

//  set event pattern for testing, units are [GeV]:

  memset(ech,0,sizeof(ech));

  const int col_shift[26] = {0, 0, 1,  1,  0, -1, -1, -1, 0, 1, 2,  2,  2,  1,  0, -1, -2, -2, -2, -2, -2, -1, 0, 1, 2, 2};
  const int row_shift[26] = {0, 0, 0, -1, -1, -1,  0,  1, 1, 1, 0, -1, -2, -2, -2, -2, -2, -1,  0,  1,  2,  2, 2, 2, 2, 1};

  int ncl = ZBQLUAB(1.01,10.99);
  for(int icl = 1; icl <= ncl; ++icl) {
    int dim   = ZBQLUAB(1.01,12.99);
    int ccol  = ZBQLUAB(1.01,(double)_NCOL_);
    int crow  = ZBQLUAB(1.01,(double)_NROW_);
    int ec    = ZBQLUAB(0.1,5.)*1.e4;
    ec /= 100;
    ec *= 100;
    ech[crow][ccol] = ec*1.e-4;
    for(int i = 2; i <= dim; ++i) {
      int col = ccol + col_shift[i];
      int row = crow + row_shift[i];
      if(col<1 || row<1 || col>=_NCOL_ || row>=_NROW_) continue;
      int e   = ec * exp(-ZBQLUAB(1.6,2.4)*hypot(col-ccol,row-crow));
      e /= 100;
      e *= 100;
      ech[row][col] = e*1.e-4;
    }
  }

  int n = 0;
  for(int i = 1; i <= _NCOL_; ++i) {
    for(int j = 1; j <= _NROW_; ++j) {
      if(ech[i][j]>MIN_COUNTER_ENERGY && n < _MADR_) {
        ++n;
        ia[n] = _OFFSET_*i+j;
        id[n] = ech[i][j];
      }
    }
  }
  *nw = n;

  return ievent;
}

void  dump_clusters(int nw, int *ia, double *id, int ncl, int *lencl, int nadcgam, adcgam_t *adcgam) {

  printf(" --- dump counters --- %11i\n",debug_ncalls);
  for(int i = 1; i <= nw; ++i) {
    printf("%3i %4i %6i\n", i, ia[i], (int)(1.e4*id[i]+0.5));
  }
  if(ncl) {
    printf(" --- dump clusters --- %11i\n",debug_ncalls);
    for(int i = 1; i <= ncl; ++i) {
      printf("%3i %4i\n", i, lencl[i]);
    }
  }
  if(nadcgam) {
    printf(" --- dump adcgam --- %11i\n",debug_ncalls);
    for(int i = 1; i <= nadcgam; ++i) {
      printf("%3i %9.6lf %9.6lf %9.6lf %17.6lf %9.6lf %9.6lf %3i %3i %3i %3i\n", i,
              adcgam[i].e, adcgam[i].x[1], adcgam[i].y[1], adcgam[i].chi2,
              adcgam[i].x[0], adcgam[i].y[0], adcgam[i].size, adcgam[i].type, adcgam[i].stat, adcgam[i].id);
    }
  }
  printf(" --- --- --- --- %11i\n",debug_ncalls);
  return;
}

#ifdef ILRECO_BUILD_TEST_MAIN
int main() {

#if   CALOR_MATERIAL==LG_CALOR
    char* profile_fname="prof_lg.dat";
#elif CALOR_MATERIAL==PWO_CALOR
    char* profile_fname="prof_pwo.dat";
#else
#error  Unknown calorimeter material
#endif

  read_profile_data(profile_fname);

  int nadcgam, ncl, nw;               // number of particles, clusters and hits
  int ia[_MADR_]; double id[_MADR_];  // arrays of adresses and energies
  int lencl[_MCL_];                   // array of clusters lengths
  adcgam_t adcgam[_MADCGAM_];         // final reconstructions storage

  while((debug_ncalls = read_event(&nw,ia,id))<=1000000) {      // event loop
    nadcgam = 0;                      // reset number of rec-d particles

    ncl = cluster_search(nw,ia,id,lencl);         //  1st stage cluster processing
    for(int icl = 1, ipncl = 1; icl <= ncl && nadcgam < _MADCGAM_-2; ipncl +=  lencl[icl++])
      process_cluster(lencl[icl],&ia[ipncl-1],&id[ipncl-1],&nadcgam,adcgam);    //  2nd stage cluster processing

    dump_clusters(nw,ia,id,ncl,lencl,nadcgam,adcgam);
  }

  return 0;
}
#endif /* ILRECO_BUILD_TEST_MAIN */

//
// initialize the random number generator
//

double ZBQLINI(int seed, double ZBQLIX[43+1]) {

  static int init = 0;
  if(init) {
    if(init==1) printf("***WARNING**** You have called routine ZBQLINI more than once. Ignoring any subsequent calls.\n");
    ++init;
    exit(1);
  } else  {init = 1;}

  double B = 4.294967291e9;
  if(!seed) seed = 1;

  ZBQLIX[1] = fmod((double)seed,B);

  for(int i = 1; i < 43; ++i) ZBQLIX[i+1] = fmod(ZBQLIX[i]*30269.,B);

  return B;
}


//
//  Returns a uniform random number between 0 & 1, using
//  a Marsaglia-Zaman type subtract-with-borrow generator
//

double  ZBQLU01() {

  static int init = 0;
  static double  ZBQLIX[43+1], B;

  if(!init) {
    double  zz[43+1];
    int iseed = 1;
    B = ZBQLINI(iseed,zz);
    memcpy(ZBQLIX,zz,sizeof(zz));
    init = 1;
  }

  static int curpos = 1, id22 = 22, id43 = 43;
  static double C = 0.;
  double x, B2 = B, BINV = 1./B;

  while(1) {
    x = ZBQLIX[id22] - ZBQLIX[id43] - C;
    if(x<0.) {x += B; C = 1.;} else {C = 0.;}
    ZBQLIX[id43] = x;
    --curpos; --id22; --id43;
    if(!curpos) {
      curpos = 43;
    } else {
      if(!id22) {id22 = 43;} else {if(!id43) id43 = 43;}
    }
    if(x<BINV) {B2 *= B;} else {break;}
  }

   return x/B2;
}

//
//  Returns a random number uniformly distributed on (x1,x2)
//  Even if x1 > x2, this will work as x2-x1 will then be -ve
//
double ZBQLUAB(double x1, double x2) {return x1+(x2-x1)*ZBQLU01();}

int nint(double x) {
  int n = (x<0.) ? -(int)(0.5-x) : (int)(0.5+x);
  return n;
}
