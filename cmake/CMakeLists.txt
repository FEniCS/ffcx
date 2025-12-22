cmake_minimum_required(VERSION 3.23)

# Development version: x.x.x.0
# Release version: x.x.x
project(ufcx VERSION 0.11.0.0 DESCRIPTION "UFCx interface header for finite element kernels"
  LANGUAGES C
  HOMEPAGE_URL https://github.com/fenics/ffcx)
include(GNUInstallDirs)

file(SHA1 ${PROJECT_SOURCE_DIR}/../ffcx/codegeneration/ufcx.h UFCX_HASH)
message("Test hash: ${UFCX_HASH}")

add_library(${PROJECT_NAME} INTERFACE)
target_compile_features(${PROJECT_NAME} INTERFACE c_std_17)
target_sources(
  ${PROJECT_NAME}
  INTERFACE FILE_SET
        public_headers
        TYPE
        HEADERS
        BASE_DIRS ../ffcx/codegeneration
        FILES ../ffcx/codegeneration/ufcx.h)

# Install header file
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}Targets
        FILE_SET public_headers)

# Install CMake Target files
install(EXPORT ${PROJECT_NAME}Targets
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

# Prepare and install CMake config files
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY AnyNewerVersion)
configure_package_config_file("${PROJECT_NAME}Config.cmake.in" "${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)
install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

# Configure and install pkgconfig file
configure_file(ufcx.pc.in ufcx.pc @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/ufcx.pc DESTINATION ${CMAKE_INSTALL_DATADIR}/pkgconfig)
