#include <iostream>

#include <ilreco.h>



int main() {
    static int debug_ncalls = 0;

    read_profile_data("prof_pwo.dat");

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


