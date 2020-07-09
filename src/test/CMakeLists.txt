cmake_minimum_required(VERSION 3.9)

project(ilreco_dump)

# => ADD Sources here:
set(ILRECO_DUMP_SOURCES ilreco_dump.cpp)

#   Most of the time you are not interested in lines below
# ==========================================================



add_executable(ilreco_dump ${ILRECO_DUMP_SOURCES})

#target_include_directories(${ILRECO_TARGET_NAME} INTERNAL ${ilreco_INCLUDE_DIRS})
target_link_libraries(ilreco_dump ilreco)
#set_target_properties(${ILRECO_TARGET_NAME} PROPERTIES PREFIX "" OUTPUT_NAME ${PLUGIN_NAME} SUFFIX ".so")

install(TARGETS ${ILRECO_TARGET_NAME})