cmake_minimum_required(VERSION 3.0)
project('il_emcal_reconstruction')

# Using flat installation (all files in the same directory)
option(FLAT_INSTALL "Using flat installation (all files in the same directory)" OFF)
if(${FLAT_INSTALL})
    message(STATUS "Using flat installation (all files in the same directory)")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)
else()
    message(STATUS "Using regular installation (.../bin, .../lib and others )")
endif()

add_subdirectory(ilreco)
add_subdirectory(test)


