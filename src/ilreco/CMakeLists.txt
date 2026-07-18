# The ilreco library target + install/export rules (find_package support).
add_library(ilreco ilreco.c)
add_library(ilreco::ilreco ALIAS ilreco)

target_include_directories(ilreco PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<INSTALL_INTERFACE:include>)
if(UNIX)
    target_link_libraries(ilreco PUBLIC m)
endif()
set_target_properties(ilreco PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER "ilreco.h;ilreco.hpp")

include(GNUInstallDirs)
install(TARGETS ilreco EXPORT ilrecoTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${PROJECT_SOURCE_DIR}/data/prof_pwo.dat
              ${PROJECT_SOURCE_DIR}/data/prof_lg.dat
        DESTINATION ${CMAKE_INSTALL_DATADIR}/ilreco)

install(EXPORT ilrecoTargets
    NAMESPACE ilreco::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ilreco)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/ilrecoConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/ilrecoConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ilreco)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/ilrecoConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/ilrecoConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/ilrecoConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ilreco)
