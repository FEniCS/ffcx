# Top level CMakeLists.txt file for UFC

# Require CMake 2.8.9 (CMake bug http://public.kitware.com/Bug/view.php?id=13216)
cmake_minimum_required(VERSION 2.8.9)

#------------------------------------------------------------------------------
# Set project name and version number

project(UFC)
set(UFC_VERSION_MAJOR "2")
set(UFC_VERSION_MINOR "2")
set(UFC_VERSION_MICRO "0+")
set(UFC_VERSION_STRING "${UFC_VERSION_MAJOR}.${UFC_VERSION_MINOR}.${UFC_VERSION_MICRO}")

# Set special link option, see `cmake --help-policy CMP0003`
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif()

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set location of our FindFoo.cmake modules
set(UFC_CMAKE_DIR "${UFC_SOURCE_DIR}/cmake" CACHE INTERNAL "")
set(CMAKE_MODULE_PATH "${UFC_CMAKE_DIR}/modules")

#------------------------------------------------------------------------------
# Options

option(UFC_ENABLE_PYTHON "Enable Python extensions." ON)

#------------------------------------------------------------------------------
# Run tests to find required packages

set(PythonInterp_FIND_VERSION 2)
find_package(PythonInterp QUIET)

# Find Python library corresponding to Python interpreter
if (UFC_ENABLE_PYTHON)
  find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT QUIET REQUIRED)
else()
  find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT QUIET)
endif()

# Find SWIG
find_package(SWIG QUIET)
include(UseSWIG)

# Look for Boost (shared_ptr support for Python wrappers)
set(Boost_ADDITIONAL_VERSIONS 1.43 1.43.0)
set(BOOST_ROOT $ENV{BOOST_DIR})
find_package(Boost 1.36 QUIET)

#------------------------------------------------------------------------------
# Get installation path for Python modules

# Get Python module path from distutils
if (PYTHONINTERP_FOUND)

  if (NOT DEFINED UFC_INSTALL_PYTHON_EXT_DIR)
    # Get path for platform-dependent Python modules (since we install a binary libary)
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
      OUTPUT_VARIABLE UFC_INSTALL_PYTHON_EXT_DIR
      )
    # Strip off CMAKE_INSTALL_PREFIX (is added later by CMake)
    string(REGEX REPLACE "${CMAKE_INSTALL_PREFIX}(/|\\\\)([^ ]*)" "\\2"
      UFC_INSTALL_PYTHON_EXT_DIR "${UFC_INSTALL_PYTHON_EXT_DIR}")
    set(UFC_INSTALL_PYTHON_EXT_DIR ${UFC_INSTALL_PYTHON_EXT_DIR}
      CACHE PATH "Python extension module installation directory.")
  endif()

  if (NOT DEFINED UFC_INSTALL_PYTHON_MODULE_DIR)
    # Get path for pure Python modules
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=False, prefix='${CMAKE_INSTALL_PREFIX}'))"
      OUTPUT_VARIABLE UFC_INSTALL_PYTHON_MODULE_DIR
      )
    # Strip off CMAKE_INSTALL_PREFIX (is added later by CMake)
    string(REGEX REPLACE "${CMAKE_INSTALL_PREFIX}(/|\\\\)([^ ]*)" "\\2"
      UFC_INSTALL_PYTHON_MODULE_DIR "${UFC_INSTALL_PYTHON_MODULE_DIR}")
    set(UFC_INSTALL_PYTHON_MODULE_DIR ${UFC_INSTALL_PYTHON_MODULE_DIR}
      CACHE PATH "Python module installation directory.")
  endif()
endif (PYTHONINTERP_FOUND)

#------------------------------------------------------------------------------
# Target names and installation directories

# Set UFC install sub-directories
if (LIB_INSTALL_DIR)
  set(UFC_LIB_DIR "${LIB_INSTALL_DIR}" CACHE PATH "Library installation directory.")
else()
  set(UFC_LIB_DIR "lib" CACHE PATH "Library installation directory.")
endif()
set(UFC_INCLUDE_DIR "include" CACHE PATH "C/C++ header installation directory.")
set(UFC_PKGCONFIG_DIR "${UFC_LIB_DIR}/pkgconfig" CACHE PATH "pkg-config file installation directory.")
set(UFC_CMAKE_CONFIG_DIR "share/ufc" CACHE PATH "CMake configuration file directory.")

#------------------------------------------------------------------------------
# Install ufc.h

set(UFC_H src/ufc/ufc.h src/ufc/ufc_geometry.h)
install(FILES ${UFC_H} DESTINATION ${UFC_INCLUDE_DIR} COMPONENT Development)

#------------------------------------------------------------------------------
# Build SWIG extension and install

if (UFC_ENABLE_PYTHON AND SWIG_FOUND AND Boost_FOUND)

  # Check SWIG version
  if (${SWIG_VERSION} LESS 2.0)
      message(ERROR " UFC requires SWIG version 2.0 or greater. You have version ${SWIG_VERSION}. Set UFC_ENABLE_PYTHON to False or install correct SWIG version.")
    endif()

  # Default to release build (can be overridden by user)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
  endif()

  # Make build directory for SWIG-generated C++ file
  FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/src/ufc")

  # In CMake 2.6 PYTHON_INCLUDE_DIRS was named PYTHON_INCLUDE_PATH
  if (NOT DEFINED PYTHON_INCLUDE_DIRS)
    set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_PATH})
  endif()

  # Include directories
  include_directories(${UFC_SOURCE_DIR}/src/ufc ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})

  # Set module name
  set(SWIG_MODULE_NAME ufc)

  # SWIG flags
  set(CMAKE_SWIG_FLAGS
    -module ${SWIG_MODULE_NAME}
    -shadow
    -modern
    -modernargs
    -fastdispatch
    -fvirtual
    -nosafecstrings
    -noproxydel
    -fastproxy
    -fastinit
    -fastunpack
    -fastquery
    -nobuildnone
    )

  # SWIG sources
  set(SWIG_SOURCES src/ufc/ufc.i)
  set(SWIG_MODULE_${SWIG_MODULE_NAME}_EXTRA_DEPS ${UFC_H})
  set_source_files_properties(${SWIG_SOURCES} PROPERTIES CPLUSPLUS ON)
  swig_add_module(${SWIG_MODULE_NAME} python ${SWIG_SOURCES})

  # Is this required?
  swig_link_libraries(ufc ${PYTHON_LIBRARIES})
  get_target_property(SWIG_MODULE_LOCATION ${SWIG_MODULE_ufc_REAL_NAME} LOCATION)
  message("SWIG_MODULE_LOCATION ${SWIG_MODULE_LOCATION}")

  # Install the swig file
  install(FILES src/ufc/ufc.i
    DESTINATION ${UFC_INCLUDE_DIR}/swig
    COMPONENT Development
    )

  # Install _ufc.so and ufc.py
  install(FILES
    ${SWIG_MODULE_LOCATION}
    ${CMAKE_CURRENT_BINARY_DIR}/ufc.py src/ufc/__init__.py
    DESTINATION ${UFC_INSTALL_PYTHON_EXT_DIR}/ufc
    COMPONENT Development
    )

endif()

#------------------------------------------------------------------------------
# Install Python utils (ufc_utils)

install(DIRECTORY
  src/utils/python/ufc_utils
  DESTINATION ${UFC_INSTALL_PYTHON_MODULE_DIR}
  USE_SOURCE_PERMISSIONS
  COMPONENT Runtime
  )

#------------------------------------------------------------------------------
# Generate UFCConfig.cmake file

configure_file(${UFC_CMAKE_DIR}/UFCConfig.cmake.in UFCConfig.cmake @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/UFCConfig.cmake
  DESTINATION ${UFC_CMAKE_CONFIG_DIR}
  COMPONENT Development
  )

#------------------------------------------------------------------------------
# Generate UFCConfigVersion.cmake file

configure_file(${UFC_CMAKE_DIR}/UFCConfigVersion.cmake.in
  UFCConfigVersion.cmake @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/UFCConfigVersion.cmake
  DESTINATION ${UFC_CMAKE_CONFIG_DIR}
  COMPONENT Development
  )

#------------------------------------------------------------------------------
# Generate UseUFC.cmake file

configure_file(${UFC_CMAKE_DIR}/UseUFC.cmake.in UseUFC.cmake @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/UseUFC.cmake
  DESTINATION ${UFC_CMAKE_CONFIG_DIR}
  COMPONENT Development
  )

#------------------------------------------------------------------------------
# Generate pkg-config config file (ufc-1.pc)

configure_file(${UFC_CMAKE_DIR}/ufc-1.pc.in ufc-1.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ufc-1.pc
  DESTINATION ${UFC_PKGCONFIG_DIR}
  COMPONENT Development
  )

#------------------------------------------------------------------------------
