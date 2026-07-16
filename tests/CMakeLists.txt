# ilreco test suite. Categories (see tests/README.md):
#   unit       — contract tests of the public API (Catch2)
#   smoke      — no-crash / sanity sweeps over generated inputs (Catch2)
#   data       — golden event-by-event regression tables (Catch2 driver + CSVs)
#   benchmark  — speed benchmark with JSON output (standalone binary)
#   tools      — golden-table generator
include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.7.1
)
FetchContent_MakeAvailable(Catch2)

set(ILRECO_PROF_PWO "${CMAKE_CURRENT_SOURCE_DIR}/../data/prof_pwo.dat")
set(ILRECO_TEST_DATA "${CMAKE_CURRENT_SOURCE_DIR}/data")

add_library(ilreco_test_common INTERFACE)
target_include_directories(ilreco_test_common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/common)
target_compile_definitions(ilreco_test_common INTERFACE
    ILRECO_PROF_PWO="${ILRECO_PROF_PWO}"
    ILRECO_TEST_DATA="${ILRECO_TEST_DATA}")
target_link_libraries(ilreco_test_common INTERFACE ilreco)

add_executable(ilreco_unit_tests
    unit/test_contracts.cpp
    unit/test_edge_inputs.cpp
    unit/test_context_api.cpp
    unit/test_cell_mask.cpp)
target_link_libraries(ilreco_unit_tests PRIVATE ilreco_test_common Catch2::Catch2WithMain)
find_package(Threads REQUIRED)
target_link_libraries(ilreco_unit_tests PRIVATE Threads::Threads)
add_test(NAME unit COMMAND ilreco_unit_tests)

add_executable(ilreco_smoke_tests smoke/test_smoke.cpp)
target_link_libraries(ilreco_smoke_tests PRIVATE ilreco_test_common Catch2::Catch2WithMain)
add_test(NAME smoke COMMAND ilreco_smoke_tests)

add_executable(ilreco_data_tests data/test_golden.cpp)
target_link_libraries(ilreco_data_tests PRIVATE ilreco_test_common Catch2::Catch2WithMain)
add_test(NAME data COMMAND ilreco_data_tests)

add_executable(ilreco_bench benchmark/bench_main.cpp)
target_link_libraries(ilreco_bench PRIVATE ilreco_test_common)
# CI-friendly quick pass (full runs: ./ilreco_bench --reps 2000)
add_test(NAME benchmark COMMAND ilreco_bench --reps 50)

add_executable(ilreco_gen_golden tools/gen_golden.cpp)
target_link_libraries(ilreco_gen_golden PRIVATE ilreco_test_common)
