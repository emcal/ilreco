//
//  gcc il.c -O3 -Wall -std=c11 -o ilc -lm
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "ilreco.h"


void  read_profile_data();
int   read_event(int *nw, int *ia, double *id);
int   cluster_search(int nw, int *ia, double *id, int *lencl);
void  process_cluster(int len, int *ia, double *id, int *nadcgam, adcgam_t *adcgam);
void  dump_clusters(int nw, int *ia, double *id, int ncl, int *lencl, int nadcgam, adcgam_t *adcgam);
void  order_hits(int nw, int *ia, double *id);
int   peak_type(int ix, int iy);
int   nint(double x);
void  gamma_cluster(int itype, int nadc, int *ia, double *id,
                    double *chisq, double *e1, double *x1, double *y1, double *e2, double *x2, double *y2);
void  fill_zero_hits(int nadc, int *ia, int *nzero, int *iazero);
void  mom1_cluster(int nadc, int *ia, double *id, int nzero, int *iazero, double *e, double *x, double *y);
void  mom2_cluster(int nadc, int *ia, double *id, int nzero, int *iazero, double *e, double *x, double *y,
                                                                          double *xx, double *yy, double *xy);

double chisq1_cluser(int nadc, int *ia, double *id, int nzero, int *iazero, double e, double x, double y);
void  tgamma_cluster(int nadc, int *ia, double *id, int nzero, int *iazero, double *chisq,
                    double *ee, double *xx, double *yy, double *e2, double *x2, double *y2);

double  profile_mean(double x, double y);
double sigma2(double dx, double dy, double e);
double  d2c(double x, double y);

void add_single_adcgam(int nadc, int *ia, double *id, int *ipnpk, int *nadcgam, adcgam_t *adcgam);
void   add_many_adcgam(int nadc, int *ia, double *id, int npk, int *ipnpk, int *nadcgam, adcgam_t *adcgam);

static double amean[_N_PROFILE_POINTS_+1][_N_PROFILE_POINTS_+1], ad2c[_N_PROFILE_POINTS_+1][_N_PROFILE_POINTS_+1];


//
// random generator functions for debugging:
//
double  ZBQLINI(int seed, double ZBQLIX[43+1]);
double  ZBQLUAB(double x1, double x2);
double  ZBQLU01();

#define _IFDB_ if(0)
//  #define _IFDB_ if(debug_ncalls == 132)
static int debug_ncalls = 0;

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

//  dump_clusters(nw,ia,id,0,lencl,0,adcgam);

    ncl = cluster_search(nw,ia,id,lencl);         //  1st stage cluster processing
    for(int icl = 1, ipncl = 1; icl <= ncl && nadcgam < _MADCGAM_-2; ipncl +=  lencl[icl++])
      process_cluster(lencl[icl],&ia[ipncl-1],&id[ipncl-1],&nadcgam,adcgam);    //  2nd stage cluster processing

//  fill_clusters();

    dump_clusters(nw,ia,id,ncl,lencl,nadcgam,adcgam);
  }

  return 0;
}


int read_event(int *nw, int *ia, double *id) {

  double  ech[_NCOL_+1][_NROW_+1];
  static int ievent = 0;
  ++ievent;

//  set event pattern for testing, units are [GeV]:

  memset(ech,0,sizeof(ech));

/*
  ech [19][19] = 1.0e-2;
  ech [19][20] = 1.2e-1;
  ech [19][21] = 2.0e-2;

  ech [20][19] = 3.0e-1;
  ech [20][20] = 2.1e+0;
  ech [20][21] = 4.3e-1;

  ech [21][20] = 5.1e-2;
  ech [10][15] = 0.3e-2;
  ech [20][25] = 3.3e-1;

  ech [21][25] = 3.5e-1;
  ech [22][25] = 4.0e-1;
  ech [23][25] = 3.1e-1;
*/

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

int cluster_search(int nw, int *ia, double *id, int *lencl) {

  if(nw<2) {
    lencl[1] = 1;
    return nw;
  }
  order_hits(nw,ia,id);       // addresses must be in increasing order

  int iwork[_MADR_];  double dwork[_MADR_];
  int ncl = 0, next  = 1, iak = 0;
  for(int k = 2; k <= nw+1; ++k) {
    if(k<=nw) iak = ia[k];
    if(iak-ia[k-1]<=1 && k<=nw) continue;

    int ib  = next; //  first word of the (sub)cluster
    int ie  = k-1;  //  last  word of the (sub)cluster
    next    = k;    //  first word of the next (sub)cluster
    if(ncl>_MCL_-2) {
      printf("maximum number of clusters reached\n");
      return _MCL_-1;
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
      if(ias-ia[last] > _OFFSET_) break;  //  no subclusters to glue
      for(int i = last; i > last-leng; --i) {
        if(ias-ia[i]  >  _OFFSET_) break;
        if(iaf-ia[i]  >= _OFFSET_) {      //  subclusters to glue
          if(icl<ncl-1 && leng<_MADR_) {

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

void  order_hits(int nw, int *ia, double *id) {

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

void  read_profile_data(const char* prof_file_name) {
  FILE *fp;
    // R_OK - OK for reading
    if( access( prof_file_name, R_OK ) == -1 ) {
        fprintf(stderr, "(!)ERROR(!) profile file '%s' does not exists or access is denied");
    }

    fp = fopen(prof_file_name,"r");
  for(int i = 0; i<=_N_PROFILE_POINTS_; ++i)
    for(int j = 0; j <= i; ++j) {
      int i1, j1; double f1, f2;
      fscanf(fp,"%d %d %lf %lf",&i1,&j1,&f1,&f2);
      amean[i][j] = f1;
      amean[j][i] = f1;
       ad2c[i][j] = f2;
       ad2c[j][i] = f2;
       if(i1!=i || j1!=j) {printf("profile data file corruption\n"); exit(1);}
    }
  fclose(fp);
}

void  process_cluster(int nadc, int *ia, double *id, int *nadcgam, adcgam_t *adcgam) {
  order_hits(nadc,ia,id);     // addresses must be in increasing order

  double minpk = 0.01;
  int ipnpk[_MPK_];

  if(nadc>=3) {
    double idsum = 0.;
    for(int i = 1; i <= nadc; ++i) idsum += id[i];
#if   CALOR_MATERIAL==LG_CALOR
#define _IBVAL_ 20.
#elif CALOR_MATERIAL==PWO_CALOR
#define _IBVAL_  7.
#else
#error  Unknown calorimeter material
#endif
    double ib = _IBVAL_*log(1.+idsum);
    if(ib>1) minpk *= ib;
    minpk = nint(1.e2*minpk);
    minpk *= 1.e2;
    minpk = 1.e-4*minpk;
#ifndef PRELIM_VERSION
#endif
  }

  int npk = 0;
  for(int ic = 1; ic <= nadc; ++ic) {
    double iac = id[ic];
    if(iac<minpk) continue;
    int ixy     = ia[ic];
    int ixymax  = ixy + _OFFSET_ + 1;
    int ixymin  = ixy - _OFFSET_ - 1;
    int iyc     = ixy - (ixy/_OFFSET_)*_OFFSET_;

    int in;
    in  = ic  + 1;
    int skipflag = 0;
    while(in<=nadc && ia[in] <= ixymax) {
      int iy  = ia[in] - (ia[in]/_OFFSET_)*_OFFSET_;
      if(abs(iy-iyc)<=1 && id[in]>=iac) {skipflag = 1; break;}
      ++in;
    }
    if(skipflag) continue;

    in = ic - 1;
    while (in>=1 && ia[in] >= ixymin) {
      int iy  = ia[in] - (ia[in]/_OFFSET_)*_OFFSET_;
      if(abs(iy-iyc)<=1 && id[in]>iac)  {skipflag = 1; break;}
      --in;
    }
    if(skipflag) continue;

    ++npk;          // peak found
    ipnpk[npk] = ic;
    if(npk >= _MPK_-2 || npk >= _MADR0_/nadc - 3) break;
  }
  if(!npk || *nadcgam >= _MADCGAM_-3)  return;

  if(npk==1) add_single_adcgam(nadc,ia,id,ipnpk,nadcgam,adcgam);
  else     add_many_adcgam(nadc,ia,id,npk,ipnpk,nadcgam,adcgam);

  return;
}


void add_single_adcgam(int nadc, int *ia, double *id, int *ipnpk, int *nadcgam, adcgam_t *adcgam) {

  const double  chisq1  = 90.; // 3.0 value of chi2 for preliminary seperation

  int n1  = ++(*nadcgam);
  int ic  = ipnpk[1];
  int ix  = ia[ic]/_OFFSET_;
  int iy  = ia[ic]-ix*_OFFSET_;
  int itype = peak_type(ix,iy);

  double  e1 = 0., x1 = 0., y1 = 0., e2 = 0., x2 = 0., y2 = 0., chisq = chisq1;
  gamma_cluster(itype,nadc,ia,id,&chisq,&e1,&x1,&y1,&e2,&x2,&y2);

  adcgam[n1].e    = e1;
  adcgam[n1].x[1] = x1;
  adcgam[n1].y[1] = y1;
  adcgam[n1].chi2 = chisq;
  adcgam[n1].size = nadc;
  adcgam[n1].type = itype;
  adcgam[n1].id   = 0;
  adcgam[n1].stat = itype;

  if(e2>0. && *nadcgam < _MADCGAM_-3) {
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

void   add_many_adcgam(int nadc, int *ia, double *id, int npk, int *ipnpk, int *nadcgam, adcgam_t *adcgam) {

  int iwrk[_MPK_+1][_MADR_], idp[_MPK_+1][_MADR_]; double fwrk[_MPK_+1][_MADR_];

  const double  idelta  =  0.; //  min cell energy part to be a member of separated cluster
  const double  chisq1  = 90.; // 3.0 value of chi2 for preliminary seperation
  const double  chisq2  = 50.; // 0.8 value of chi2 for final seperation
  const int     niter   =   6;
  int ngam0 = *nadcgam;
  int igmpk[_MPK_][3];

  double ratio = 1.;
  for(int iter = 1; iter <= niter; ++iter) {
    for(int i = 1; i <= nadc; ++i) {iwrk[0][i] = 0; fwrk[0][i] = 0.;}
    double epk[_MPK_], xpk[_MPK_], ypk[_MPK_];
    for(int ipk = 1; ipk <= npk; ++ipk) {
      int ic = ipnpk[ipk];
      if(iter!=1) ratio = fwrk[ipk][ic] / fwrk[npk+1][ic];
      double eg = id[ic] * ratio;
      int ixypk = ia[ic];
      int ixpk  = ixypk / _OFFSET_;
      int iypk  = ixypk - ixpk * _OFFSET_;
      epk[ipk]  = eg;
      xpk[ipk]  = eg*ixpk;
      ypk[ipk]  = eg*iypk;
      for(int in = ic + 1; in <= nadc; ++in) {
        int ixy = ia[in];
        int ix  = ixy / _OFFSET_;
        int iy  = ixy - ix * _OFFSET_;
        if(ixy-ixypk>_OFFSET_+1) break;
        if(abs(iy-iypk)>1) continue;
        if(iter!=1) ratio = fwrk[ipk][in] / fwrk[npk+1][in];
        eg = id[in] * ratio;
        epk[ipk] += eg;
        xpk[ipk] += eg*ix;
        ypk[ipk] += eg*iy;
      }

      for(int in = ic - 1; in>0; --in) {
        int ixy = ia[in];
        int ix  = ixy / _OFFSET_;
        int iy  = ixy - ix * _OFFSET_;
        if(ixypk-ixy>_OFFSET_+1) break;
        if(abs(iy-iypk)>1) continue;
        if(iter!=1) ratio = fwrk[ipk][in] / fwrk[npk+1][in];
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
        int ix  = ixy / _OFFSET_;
        int iy  = ixy - ix * _OFFSET_;
        double dx = fabs(ix-xpk[ipk]);
        double dy = fabs(iy-ypk[ipk]);
        double a  = epk[ipk]*profile_mean(dx,dy);
        iwrk[ipk][i] = nint(a*1.e4);
        iwrk[0][i]  += iwrk[ipk][i];
        fwrk[ipk][i] = a;
        fwrk[0][i]  += fwrk[ipk][i];
      }
    }

    for(int i = 1; i <= nadc; ++i) {
      int iwk = iwrk[0][i];
      if(iwk<1) iwk = 1;
      iwrk[npk+1][i] = iwk;
      if(fwrk[0][i]>1.e-6)
        fwrk[npk+1][i] = fwrk[0][i];
      else
        fwrk[npk+1][i] = 1.e-6;
    }
  }

  for(int ipk = 1; ipk <= npk; ++ipk) {
    int leng = 0;
    for(int i = 1; i <= nadc; ++i) {
      if(fwrk[0][i] <= 1.e-6) continue;
      int ixy   = ia[i];
      double fe = 1.e4*id[i]*fwrk[ipk][i]/fwrk[0][i];
      if(fe<=idelta)  continue;
      ++leng;
      iwrk[npk+1][leng] = ixy;
      iwrk[npk+2][leng] = nint(fe);
      fwrk[npk+1][leng] = ixy;
      fwrk[npk+2][leng] = 1.e-4*fe;
    }

    if(*nadcgam >= _MADCGAM_-2) return;
    igmpk[ipk][2] = 0;
    if(!leng) continue;
    int n1  = ++(*nadcgam);
    int ic  = ipnpk[ipk];
    int ix  = ia[ic] / _OFFSET_;
    int iy  = ia[ic] - ix * _OFFSET_;
    int itype = peak_type(ix,iy);

    double  e1 = 0., x1 = 0., y1 = 0., e2 = 0., x2 = 0., y2 = 0., chisq = chisq1;
    gamma_cluster(itype,leng,&iwrk[npk+1][0],&fwrk[npk+2][0],&chisq,&e1,&x1,&y1,&e2,&x2,&y2);

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

    if(e2>0. && n1 < _MADCGAM_-3) {
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

  for(int i = 1; i <= nadc; ++i) {iwrk[0][i] = 0; idp[0][i] = 0; fwrk[0][i] = 0.;}

  for(int ipk = 1; ipk <= npk; ++ipk)
    for(int i = 1; i <= nadc; ++i) {
      iwrk[ipk][i] = 0; idp[ipk][i] = 0; fwrk[ipk][i] = 0.;
      if(igmpk[ipk][2])
        for(int ig = igmpk[ipk][1]; ig <= igmpk[ipk][2]; ++ig) {

          int ixy = ia[i];
          int ix  = ixy / _OFFSET_;
          int iy  = ixy - ix * _OFFSET_;
          double dx = ix-adcgam[ig].x[1];
          double dy = iy-adcgam[ig].y[1];
          double a  = adcgam[ig].e*profile_mean(dx,dy);
          int   iia = nint(a*1.e4);
          iwrk[ipk][i] += iia;
          iwrk[0][i]   += iia;
          fwrk[ipk][i] += a;
          fwrk[0][i]   += a;
          idp[ipk][i]  += iia;
        }
    }

  for(int i = 1; i <= nadc; ++i) {
    idp[0][i] = 0;
    for(int ipk = 1; ipk <= npk; ++ipk)
      idp[0][i] += idp[ipk][i];
    int ide = nint(id[i]*1.e4) - idp[0][i];
    if(!ide || fwrk[0][i] == 0.) continue;

    double fw[_MPK_];
    for(int ipk = 1; ipk <= npk; ++ipk)
      fw[ipk] = fwrk[ipk][i]/fwrk[0][i];

    int idecorr = 0;
    for(int ipk = 1; ipk <= npk; ++ipk) {
      double fia = ide*fw[ipk];
      if(fwrk[ipk][i]+fia>0.) {
        fwrk[ipk][i] += fia;
        fwrk[0][i]   += fia;
      }
      int iia = nint(fia);
      if(iwrk[ipk][i]+iia>0)  {
        iwrk[ipk][i] += iia;
        iwrk[0][i]   += iia;
        idecorr      += iia;
      } else {
        if(iwrk[ipk][i]+iia<0)
          printf("WARNING NEGATIVE CORR = %i %lf\n", ia[i], id[i]);
      }
    }
  }

  *nadcgam  = ngam0;          //  reanalize last (two) gamma(s)
  for(int ipk = 1; ipk <= npk; ++ipk) {
    int leng = 0;
    for(int i = 1; i <= nadc; ++i) {
      if(iwrk[0][i] <= 0) continue;
      double fe = 1.e4*id[i]*fwrk[ipk][i]/fwrk[0][i];
      if(fe<=idelta)  continue;
      ++leng;
      iwrk[npk+1][leng] = ia[i];
      iwrk[npk+2][leng] = nint(fe);
      fwrk[npk+1][leng] = ia[i];
      fwrk[npk+2][leng] = 1.e-4*fe;
    }
    if(*nadcgam >= _MADCGAM_-2) return;
    if(!leng) continue;
    int n1  = ++(*nadcgam);
    int ic  = ipnpk[ipk];
    int ix  = ia[ic] / _OFFSET_;
    int iy  = ia[ic] - ix * _OFFSET_;
    int itype = peak_type(ix,iy);

    double  e1 = 0., x1 = 0., y1 = 0., e2 = 0., x2 = 0., y2 = 0., chisq = chisq2;
    gamma_cluster(itype,leng,&iwrk[npk+1][0],&fwrk[npk+2][0],&chisq,&e1,&x1,&y1,&e2,&x2,&y2);

    adcgam[n1].e    = e1;
    adcgam[n1].x[1] = x1;
    adcgam[n1].y[1] = y1;
    adcgam[n1].chi2 = chisq;
    adcgam[n1].size = leng;
    adcgam[n1].type = itype;
    adcgam[n1].id   = 10;
    adcgam[n1].stat = itype;

    if(e2>0. && n1 < _MADCGAM_-3) {
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
        adcgam[n1].element[j] = iwrk[npk+1][j];
        adcgam[n2].element[j] = iwrk[npk+1][j];
        adcgam[n1].elfract[j] = nint(iwrk[npk+2][j]*e1/(e1+e2));
        adcgam[n2].elfract[j] = nint(iwrk[npk+2][j]*e2/(e1+e2));
      }
    } else {
      for(int j = 1; j <= leng; ++j) {
        if(j>=_MAXLEN_) break;
        adcgam[n1].element[j] = iwrk[npk+1][j];
        adcgam[n1].elfract[j] = iwrk[npk+2][j];
      }
    }
  }

  return;
}

int peak_type(int ix, int iy) {

  if( (ix == _NCOL_/2-1 || ix == _NCOL_/2+2) &&
       iy >= _NROW_/2-1 && iy <= _NROW_/2+2) return 1;      //  hole

  if( (iy == _NROW_/2-1 || iy == _NROW_/2+2) &&
       ix >= _NCOL_/2-1 && ix <= _NCOL_/2+2) return 1;

  if( ix == 1 || ix == _NCOL_ || iy == 1 || iy == _NROW_) return 2;   //  transition (outer boundary)

  return 0;         // inside matrix
}

void  gamma_cluster(int itype, int nadc, int *ia, double *id,
                    double *chisq, double *e1, double *x1, double *y1, double *e2, double *x2, double *y2)  {

  const double  xm2cut = 1.7e-3*_ZCAL_;

  *e2 = 0.; *x2 = 0.; *y2 = 0.;
  int nzero, iazero[_MADR_];
  fill_zero_hits(nadc,ia,&nzero,iazero);          // make use of good but zero signal cells around clusters
  mom1_cluster(nadc,ia,id,nzero,iazero,e1,x1,y1); // get initial e,x,y
  if(nadc <= 0) return;

  double chimem   = *chisq;   // memorize reference chi2 value
  double  chi0    = chisq1_cluser(nadc,ia,id,nzero,iazero,*e1,*x1,*y1);  //  get initial chi2
  double  chisq0  = chi0;
  int        ndof = nzero + nadc - 2;   //  number of degrees of freedom
  if(ndof<1) ndof = 1;

  *chisq = chi0/(double)ndof;
  double x0 = *x1, y0 = *y1;
  while(1) {        // iteration loop

    double const dxy = 0.05, stepmin = .002;
    double stepx, stepy;
    double chiright = chisq1_cluser(nadc,ia,id,nzero,iazero,*e1,x0+dxy,y0);
    double chileft  = chisq1_cluser(nadc,ia,id,nzero,iazero,*e1,x0-dxy,y0);
    double chiup    = chisq1_cluser(nadc,ia,id,nzero,iazero,*e1,x0,y0+dxy);
    double chidown  = chisq1_cluser(nadc,ia,id,nzero,iazero,*e1,x0,y0-dxy);

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

    double chi00 = chisq1_cluser(nadc,ia,id,nzero,iazero,*e1,x0+stepx,y0+stepy);
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
  tgamma_cluster(nadc,ia,id,nzero,iazero,chisq,&ee,&xx,&yy,e2,x2,y2);

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

void  fill_zero_hits(int nadc, int *ia, int *nzero, int *iazero) {
  *nzero = 0;

  for(int i = 1; i <= nadc; ++i) {
    int ix = ia[i]/_OFFSET_;
    int iy = ia[i] - ix * _OFFSET_;

    if(ix>1)  {
      iazero[++(*nzero)] = iy + (ix-1) * _OFFSET_;                    //  left neib
      if(iy>1)      iazero[++(*nzero)] = (iy-1) + (ix-1) * _OFFSET_;  //  bottom left neib
      if(iy<_NROW_) iazero[++(*nzero)] = (iy+1) + (ix-1) * _OFFSET_;  //  top left neib
    }
    if(ix<_NCOL_)  {
      iazero[++(*nzero)] = iy + (ix+1) * _OFFSET_;                    //  right neib
      if(iy>1)      iazero[++(*nzero)] = (iy-1) + (ix+1) * _OFFSET_;  //  bottom right neib
      if(iy<_NROW_) iazero[++(*nzero)] = (iy+1) + (ix+1) * _OFFSET_;  //  top right neib
    }
    if(iy>1)        iazero[++(*nzero)] = iy-1 + ix * _OFFSET_;        //  bottom neib
    if(iy<_NROW_)   iazero[++(*nzero)] = iy+1 + ix * _OFFSET_;        //  top neib
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

void  mom1_cluster(int nadc, int *ia, double *id, int nzero, int *iazero, double *e, double *x, double *y) {

  *e = *x = *y = 0.;
  for(int i = 1; i <= nadc; ++i) {
    double a  = id[i];
    int ix    = ia[i]/_OFFSET_;
    int iy    = ia[i] - ix*_OFFSET_;
    *e  += a;
    *x  += a*ix;
    *y  += a*iy;
  }

  if(*e <= 0.) return;
  *x  /= *e;
  *y  /= *e;

  double corr = 0.;
  for(int i = 1; i <= nadc; ++i) {
    double dx = (ia[i]/_OFFSET_)-(*x);
    double dy = (ia[i] - (ia[i]/_OFFSET_)*_OFFSET_)-(*y);
    corr  += profile_mean(dx,dy);
  }
  for(int i = 1; i <= nzero; ++i) {
    double dx = (iazero[i]/_OFFSET_)-(*x);
    double dy = (iazero[i] - (iazero[i]/_OFFSET_)*_OFFSET_)-(*y);
    corr  += profile_mean(dx,dy);
  }

  corr /= 1.006;    // to be eliminated
  if(corr<0.8)  corr = 0.8;
  if(corr>1.0)  corr = 1.0;
  *e   /= corr;
  return;
}

void  mom2_cluster(int nadc, int *ia, double *id, int nzero, int *iazero, double *e, double *x, double *y,
                                                                          double *xx, double *yy, double *xy) {
  *e = *x = *y = *xx = *yy = *xy = 0.;
  for(int i = 1; i <= nadc; ++i) {
    double a  = id[i];
    int ix    = ia[i]/_OFFSET_;
    int iy    = ia[i] - ix*_OFFSET_;
    *e  += a;
    *x  += a*ix;
    *y  += a*iy;
  }

  if(*e <= 0.) return;
  *x  /= *e;
  *y  /= *e;

  for(int i = 1; i <= nadc; ++i) {
    double a  = id[i];
    int ix    = ia[i]/_OFFSET_;
    int iy    = ia[i] - ix*_OFFSET_;
    *xx += a/(*e)*(ix-(*x))*(ix-(*x));
    *yy += a/(*e)*(iy-(*y))*(iy-(*y));
    *xy += a/(*e)*(ix-(*x))*(iy-(*y));
  }

  double corr = 0.;
  for(int i = 1; i <= nadc; ++i) {
    double dx = (ia[i]/_OFFSET_)-(*x);
    double dy = (ia[i] - (ia[i]/_OFFSET_)*_OFFSET_)-(*y);
    corr  += profile_mean(dx,dy);
  }
  for(int i = 1; i <= nzero; ++i) {
    double dx = (iazero[i]/_OFFSET_)-(*x);
    double dy = (iazero[i] - (iazero[i]/_OFFSET_)*_OFFSET_)-(*y);
    corr  += profile_mean(dx,dy);
  }

  corr /= 1.006;    // to be eliminated
  if(corr<0.8)  corr = 0.8;
  if(corr>1.0)  corr = 1.0;
  *e   /= corr;
  return;
}

double chisq1_cluser(int nadc, int *ia, double *id, int nzero, int *iazero, double e, double x, double y)  {

  double chi2 = 0.;

  for(int i = 1; i <= nadc; ++i) {
    int ix    = ia[i]/_OFFSET_;
    int iy    = ia[i] - ix*_OFFSET_;
    if(e) {
      if(fabs(x-ix)<=6.0 && fabs(y-iy)<=6.0)  {
        double f = (profile_mean(x-ix,y-iy) - id[i]/e);
        chi2 += 1.e4*e*f*f/sigma2(x-ix,y-iy,e);
      }
    } else {
      chi2 += id[i]*id[i] / 9.;
      printf(" case 0 ch\n"); //  should'nt come here
    }
  }

  for(int i = 1; i <= nzero; ++i) {
    int ix    = iazero[i]/_OFFSET_;
    int iy    = iazero[i] - ix*_OFFSET_;
    if(e) {
      double f = profile_mean(x-ix,y-iy);
      chi2 += 1.e4*e*f*f/sigma2(x-ix,y-iy,e);
    } else {
      printf(" case 0 ch\n"); //  should'nt come here
    }
  }

  return  chi2;
}

double sigma2(double dx, double dy, double e) {

  if(dx*dx+dy*dy>25.) return 1.e2;

  double const alp = 0.816, bet1 = 32.1, bet2 = 17.2;

  double retval = 1.e2*(alp*profile_mean(dx,dy) + (bet1+bet2*sqrt(e))*d2c(dx,dy) + 2.e-3/e);

  return retval*pow(e,-0.166);
}

double  d2c(double x, double y) {
  double ax = fabs(x*1.e2), ay = fabs(y*1.e2);
  int i = (int)ax, j = (int)ay;
  double wx = ax-i, wy = ay-j;
  if(i<_N_PROFILE_POINTS_ && j<_N_PROFILE_POINTS_)
    return ad2c[i][j]     * (1.-wx) * (1.-wy) +
           ad2c[i+1][j]   *     wx  * (1.-wy) +
           ad2c[i][  j+1] * (1.-wx) *     wy  +
           ad2c[i+1][j+1] *     wx  *     wy;
  return 1.;
}

double chisq2t_hyc( double ecell, double e1, double dx1, double dy1,
                                  double e2, double dx2, double dy2, double f1, double f2) {
  double s;
  int ic = !(!(e1))*10+!(!(e2));
  switch(ic) {
    case 11:
      s = e1*sigma2(dx1,dy1,e1) + e2*sigma2(dx2,dy2,e2);
      break;
    case 10:
      s = e1*sigma2(dx1,dy1,e1);
      break;
    case  1:
      s = e2*sigma2(dx2,dy2,e2);
      break;
    default:
      s = 9e-4;
  }
    s *= 1.e8;
    double d = e1*f2+e2*f2-ecell;
    return d*d/s;
}

void  tgamma_cluster(int nadc, int *ia, double *id, int nzero, int *iazero, double *chisq,
                    double *ee, double *xx, double *yy, double *e2, double *x2, double *y2) {}

double  profile_mean(double x, double y) {
  double ax = fabs(x*1.e2), ay = fabs(y*1.e2);
  int i = (int)ax, j = (int)ay;
  double wx = ax-i, wy = ay-j;

  if(i<_N_PROFILE_POINTS_ && j<_N_PROFILE_POINTS_)
    return amean[i][j]     * (1.-wx) * (1.-wy) +
           amean[i+1][j]   *     wx  * (1.-wy) +
           amean[i][j+1]   * (1.-wx) *     wy  +
           amean[i+1][j+1] *     wx  *     wy;
  return 0.;
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

//  ++debug_ncalls;

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
