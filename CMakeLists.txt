# Top level CMakeLists.txt file for UFC

#------------------------------------------------------------------------------

# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

# Set verbose output while testing CMake
set(CMAKE_VERBOSE_MAKEFILE 1)

#------------------------------------------------------------------------------
project(UFC)

# Set special link option, see `cmake --help-policy CMP0003`
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif()

#------------------------------------------------------------------------------
# Set UFC version number

set(UFC_VERSION_MAJOR "0")
set(UFC_VERSION_MINOR "9")
set(UFC_VERSION_MICRO "8")
set(UFC_VERSION "${UFC_VERSION_MAJOR}.${UFC_VERSION_MINOR}.${UFC_VERSION_MICRO}")
#add_definitions(-DPACKAGE_VERSION="${UFC_VERSION}")

# Compiler flags
# FIXME: Default to optimisation?

#------------------------------------------------------------------------------
# Options

option(UFC_ENABLE_PYTHON "Enable Python extensions." ON)

#------------------------------------------------------------------------------
# Run tests to find required packages

find_package(PythonInterp QUIET)
find_package(PythonLibs QUIET)
find_package(SWIG QUIET)
include(UseSWIG)

# FIXME: Look for Boost (shared_ptr support)

#------------------------------------------------------------------------------
# Get installation path for Python modules

# Get Python module path from distutils
if (PYTHONINTERP_FOUND)
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
    OUTPUT_VARIABLE PYTHON_INSTALL_PATH
    )

  # Set Python installation path
  set(UFC_INSTALL_PYTHON_MODULE_DIR ${PYTHON_INSTALL_PATH})
  set(UFC_INSTALL_PYTHON_EXT_DIR ${PYTHON_INSTALL_PATH})
endif (PYTHONINTERP_FOUND)

#------------------------------------------------------------------------------
# Target names and installation directories

# Set UFC install sub-directories
set(UFC_LIB_DIR "lib")
set(UFC_INCLUDE_DIR "include")
set(UFC_PKGCONFIG_DIR "lib/pkgconfig")
set(UFC_SHARE_DIR "share")
set(UFC_MAN_DIR "share/man")
#set(UFC_CMAKE_CONFIG_DIR "share/ufc")

#------------------------------------------------------------------------------
# Install ufc.h

set(UFC_H src/ufc/ufc.h)
install(FILES ${UFC_H} DESTINATION ${UFC_INCLUDE_DIR} COMPONENT Development)

#------------------------------------------------------------------------------
# Build SWIG extension and install

if (UFC_ENABLE_PYTHON AND SWIG_FOUND)

  set(SWIG_MODULE_NAME ufc)
  set(CMAKE_SWIG_FLAGS -module ${SWIG_MODULE_NAME})
  set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR})
  set(SWIG_SOURCES src/ufc/ufc.i)
  set_source_files_properties(${SWIG_SOURCES} PROPERTIES CPLUSPLUS ON)

  swig_add_module(${SWIG_MODULE_NAME} python ${SWIG_SOURCES})
  include_directories(${PYTHON_INCLUDE_DIRS})

  # Is this required?
  swig_link_libraries(ufc ${PYTHON_LIBRARIES})


  get_target_property(SWIG_MODULE_LOCATION ${SWIG_MODULE_ufc_REAL_NAME} LOCATION)

  # Install the swig file
  install(FILES src/ufc/ufc.i
    DESTINATION ${UFC_INCLUDE_DIR}/swig
    COMPONENT Development
    )

  # Install _ufc.so and ufc.py
  install(FILES ${SWIG_MODULE_LOCATION} ${CMAKE_CURRENT_BINARY_DIR}/ufc.py
    DESTINATION ${UFC_INSTALL_PYTHON_MODULE_DIR}
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
# Generate ufc-config.cmake file

# FIXME:

#------------------------------------------------------------------------------
