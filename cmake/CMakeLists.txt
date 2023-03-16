cmake_minimum_required(VERSION 3.19)

project(ufcx VERSION 0.7.0 DESCRIPTION "UFCx interface header for finite element kernels"
  LANGUAGES C
  HOMEPAGE_URL https://github.com/fenics/ffcx)
include(GNUInstallDirs)

file(SHA1 ${PROJECT_SOURCE_DIR}/../ffcx/codegeneration/ufcx.h UFCX_HASH)
message("Test hash: ${UFCX_HASH}")

add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Prepare and install CMake target/config files
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY AnyNewerVersion)
configure_package_config_file("${PROJECT_NAME}Config.cmake.in" "${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)
install(EXPORT ${PROJECT_NAME}_Targets FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)
install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

# Install header file
install(FILES ${PROJECT_SOURCE_DIR}/../ffcx/codegeneration/ufcx.h TYPE INCLUDE)

# Configure and install pkgconfig file
configure_file(ufcx.pc.in ufcx.pc @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/ufcx.pc DESTINATION ${CMAKE_INSTALL_DATADIR}/pkgconfig)

