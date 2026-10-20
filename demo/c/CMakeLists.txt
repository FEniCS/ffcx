cmake_minimum_required(VERSION 3.21...4.4)

project(ffcx-c-demo LANGUAGES C)

find_package(ufcx REQUIRED)
find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Compile the UFL file with FFCx. This writes poisson.c and poisson.h
# into the build directory.
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/poisson.c ${CMAKE_CURRENT_BINARY_DIR}/poisson.h
  COMMAND Python3::Interpreter -m ffcx --dir ${CMAKE_CURRENT_BINARY_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/poisson.py
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/poisson.py
  COMMENT "Compiling poisson.py with FFCx"
  VERBATIM)

add_executable(poisson main.c ${CMAKE_CURRENT_BINARY_DIR}/poisson.c)
set_target_properties(poisson PROPERTIES
  C_STANDARD 17
  C_STANDARD_REQUIRED ON
  C_EXTENSIONS OFF)
target_include_directories(poisson PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(poisson PRIVATE ufcx::ufcx)

enable_testing()
add_test(NAME poisson COMMAND poisson)

# The generated code calls into the C maths library. It is part of the
# C runtime on some platforms, so link it only where it is a separate
# library.
find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
  target_link_libraries(poisson PRIVATE ${MATH_LIBRARY})
endif()
