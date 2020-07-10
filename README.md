# ilreco

Calorimeter reconstruction library

## Get & build

```bash
git clone https://github.com/emcal/ilreco.git
mkdir ilreco/build
cd ilreco/build
cmake ../src
make          # build. One can just call `make install` 
make install  # (!) tries to install in system /bin /lib etc.


# To install everything in the build directory instead of system lib
# (kind of more development mode)
cmake -L -DCMAKE_INSTALL_PREFIX=`pwd` -DCMAKE_BUILD_TYPE=Debug -DFLAT_INSTALL=ON ../src
make

```

Author: Ilya Larin

