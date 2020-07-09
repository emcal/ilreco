# This cmake file here is only to work with both plugins simultaneously within IDEs such as CLion or VSCode

cmake_minimum_required(VERSION 3.9)
project('il_reconstruction')

add_subdirectory(ilreco)
add_subdirectory(test)
