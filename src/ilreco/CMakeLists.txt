cmake_minimum_required(VERSION 3.9)

set(ILRECO_LIBNAME ilreco)

# => ADD Sources here:
set(ILRECO_SOURCES ilreco.c)

set(ILRECO_TARGET_NAME ${ILRECO_LIBNAME})

#   Most of the time you are not interested in lines below
# ==========================================================

project(${ILRECO_LIBNAME})

add_library(${ILRECO_TARGET_NAME} SHARED ${ILRECO_SOURCES})

target_include_directories(ilreco PUBLIC ${CMAKE_CURRENT_LIST_DIR})
# target_link_libraries(${PLUGIN_TARGET_NAME}  dl pthread)
#set_target_properties(${ILRECO_TARGET_NAME} PROPERTIES PREFIX "" OUTPUT_NAME ${PLUGIN_NAME} SUFFIX ".so")

install(TARGETS ${ILRECO_TARGET_NAME})