# ilreco

Calorimeter reconstruction library

## Get & build

```bash
git clone https://github.com/emcal/ilreco.git
mkdir ilreco/build &&cd ilreco/build
cmake ../src
make          # build. One can just call `make install` 
make install  # (!) tries to install in system /bin /lib etc.
```

To install everything in the build directory instead of system lib:

```bash
# (kind of more development mode)
cmake -DCMAKE_INSTALL_PREFIX=`pwd` -DCMAKE_BUILD_TYPE=Debug -DFLAT_INSTALL=ON ../src
make
```

Special cmake flags:

```bash
-DFLAT_INSTALL=OFF   # ON  - all files in the same directory. 
                     # OFF - put binaries and libraries to .../bin .../lib etc.
```

## Possibilities

...

## Isle clustering algorithm

Algorithm first finds islands as connected areas in separate HyCal2
sectors.  Each island is subject to search for maxima.  
In case of many maxima are found each of them will be associated with the separate hit.
Each single hit is also subject to test if it could be split into two close hits (second step of separation). 
Clusters found in all HyCal sectors are subject to merge in case if they pass the close-enough test. 
The second step of separation was suppressed by high cut values because it can split single clusters 
with probability of few percent and real hits will be close enough to be subject of this step only 
in fraction of percent of all events.Thus the algorithm mostly seprates hits which produce different maxima.
The minimum distance between them to be resolved 


## Support

Author: Ilya Larin

