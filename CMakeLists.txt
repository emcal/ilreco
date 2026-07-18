# ilreco — island-clustering reconstruction for square-cell calorimeters.
#
# Consumable three ways:
#   installed package:  find_package(ilreco)  ->  target ilreco::ilreco
#   embedded:           add_subdirectory / FetchContent -> same target
#   python wheel:       pip install .  (scikit-build-core drives this file
#                       with ILRECO_BUILD_PYTHON=ON)
cmake_minimum_required(VERSION 3.18)
project(ilreco VERSION 1.0.0 LANGUAGES C CXX)

option(ILRECO_BUILD_TESTS "Build the ilreco test suite" ${PROJECT_IS_TOP_LEVEL})
option(ILRECO_BUILD_PYTHON "Build the python extension module" OFF)

add_subdirectory(src/ilreco)

if(ILRECO_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(ILRECO_BUILD_PYTHON)
    find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
    # the extension compiles the library source directly: the wheel stays
    # self-contained (no shared-library rpath juggling)
    python_add_library(_core MODULE python/src/_core.c src/ilreco/ilreco.c
                       WITH_SOABI)
    target_include_directories(_core PRIVATE src/ilreco)
    # same cross-platform bit-identity requirement as the library target
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(_core PRIVATE -ffp-contract=off)
    endif()
    target_link_libraries(_core PRIVATE Python::NumPy)
    if(UNIX)
        target_link_libraries(_core PRIVATE m)
    endif()
    install(TARGETS _core DESTINATION ilreco)
    install(FILES data/prof_pwo.dat data/prof_lg.dat DESTINATION ilreco/data)
endif()
