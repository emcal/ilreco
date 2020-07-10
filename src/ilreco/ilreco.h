#ifndef CALORIMETRY_STUDIES_ILRECO_H
#define CALORIMETRY_STUDIES_ILRECO_H

#ifdef __cplusplus
extern "C" {
#endif

#define _MADR0_   10000       // number of calorimeter modules
#define _MADR_OFFSET_ 800     // number of calorimeter modules safety offset
#define _MADR_    (_MADR0_+_MADR_OFFSET_)         // number of calorimeter modules with some safety addition
#define _MCL_       200       // max number+1 of raw clusters
#define _MADCGAM_   200       // max number+1 of reconstructed "particles"
#define _NCOL_       34       // number of columns in calorimeter structure
#define _NROW_       34       // number of   rows  in calorimeter structure
#define _OFFSET_    100       // column offset in array addressing scheme
#define _MAXLEN_    100       // max length (number of hits) in final reconstruction object ("particle")

#define _MPK_        12       // max number of peaks in raw cluster

#define _N_REC_METHODS_ 5     // number of coord. recosntruction methods reserved
#define MIN_COUNTER_ENERGY  1.e-3
#define CLUSTER_MIN_ENERGY  0.1
#define _ZCAL_      732.      // target distance to calorimeter to define 2nd step separation cut

#define PWO_CALOR 1
#define LG_CALOR  2
#define CALOR_MATERIAL PWO_CALOR

#define _N_PROFILE_POINTS_  500 // number of 2d profile nods

typedef struct {
    double e;
    double x[_N_REC_METHODS_];
    double y[_N_REC_METHODS_];
    double z[_N_REC_METHODS_];
    double chi2;
    int    size;              // number of hits
    int    type;              // type of cluster (in most cases fiducial region of the calorimeter and if was merged from different subparts)
    int      id;              // guessed particle id
    int    stat;              //  status of cluster
    int    element[_MAXLEN_]; // link to hit array
    double elfract[_MAXLEN_]; // fraction of hit belonging to this cluster (important for overlapped clusters)
} adcgam_t;



void  read_profile_data(const char* prof_file_name);
int  cluster_search(int nw, int *ia, double *id, int *lencl);
void process_cluster(int nadc, int *ia, double *id, int *nadcgam, adcgam_t *adcgam);
int  read_event(int *nw, int *ia, double *id);
void dump_clusters(int nw, int *ia, double *id, int ncl, int *lencl, int nadcgam, adcgam_t *adcgam);

#ifdef __cplusplus
}
#endif

#endif //CALORIMETRY_STUDIES_ILRECO_H
