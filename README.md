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


## Support

Author: Ilya Larin

