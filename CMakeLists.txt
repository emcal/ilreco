# Top-level build: the ilreco library (src/) plus the test suite (tests/).
# The library sources are NOT modified by the test effort; see tests/README.md.
cmake_minimum_required(VERSION 3.14)
project(ilreco_dev LANGUAGES C CXX)

add_subdirectory(src)

option(ILRECO_BUILD_TESTS "Build the ilreco test suite" ON)
if(ILRECO_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
