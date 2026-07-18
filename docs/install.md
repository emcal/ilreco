# Install

Linux and macOS are supported. Pick the route that matches how you work.

## Python

```bash
pip install ilreco
```

Ships the compiled library, the Python API, and the PbWO4 / lead-glass
shower-profile tables. numpy is the only dependency.

```python
import ilreco
calo = ilreco.Calorimeter(20, 20, profile="pwo")
```

To build from a source checkout instead (requires CMake ≥ 3.18 and a C
compiler):

```bash
pip install ./ilreco
```

## Plain C — grab the files

The library is two files with no dependencies. Copy them into your project:

```
src/ilreco/ilreco.h     the public API
src/ilreco/ilreco.c     the implementation
src/ilreco/ilreco.hpp   optional C++ RAII wrapper
data/prof_pwo.dat       shower profile for 2.05 cm PbWO4 cells
data/prof_lg.dat        shower profile for lead glass
```

```bash
cc -O3 -std=c11 -c ilreco.c && cc your_program.c ilreco.o -lm
```

## CMake

Build and install:

```bash
git clone https://github.com/emcal/ilreco
cmake -B build ilreco
cmake --build build -j
cmake --install build            # add --prefix <dir> for a local install
```

This installs the library, both headers, the profile tables
(`<prefix>/share/ilreco/`) and CMake package files. Consume it with
`find_package`:

```cmake
find_package(ilreco REQUIRED)
target_link_libraries(your_target PRIVATE ilreco::ilreco)
```

Or embed it directly with FetchContent — the test suite is skipped
automatically when ilreco is not the top-level project:

```cmake
include(FetchContent)
FetchContent_Declare(ilreco
    GIT_REPOSITORY https://github.com/emcal/ilreco
    GIT_TAG        v1.0.0)
FetchContent_MakeAvailable(ilreco)
target_link_libraries(your_target PRIVATE ilreco::ilreco)
```
