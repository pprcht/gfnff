# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
#
# gfnff is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gfnff is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

cmake_minimum_required(VERSION 3.17)

project(
  "gfnff"
  LANGUAGES "Fortran" "C" "CXX"
  VERSION "0.0.1"
  DESCRIPTION "Standalone implementation of general the `Geometries, Frequencies, Non-covalent interactions Force-Field` (GFN-FF)"
)


# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# General configuration information
set(libs)
add_subdirectory("config")

###############################################################################
####################### SUBPROJECTS & DEPENDENCIES ############################
###############################################################################

# OpenMP parallelization
if(NOT TARGET "OpenMP::OpenMP_Fortran")
  if (WITH_OpenMP)
    find_package("OpenMP" REQUIRED)
    add_compile_definitions(WITH_OpenMP) 
  else()
    # Create dummy library
    add_library("OpenMP::OpenMP_Fortran" INTERFACE IMPORTED)
  endif()
endif()

# Linear algebra (contains BLAS)
if(NOT TARGET "LAPACK::LAPACK")
  find_package("LAPACK" REQUIRED)
endif()

# Fortran unit test interface (also used by other subprojects)
if(NOT TARGET "test-drive::test-drive" AND WITH_TESTS)
   find_package("test-drive" REQUIRED)
endif()

set(
  lib-deps
  "OpenMP::OpenMP_Fortran"
  "LAPACK::LAPACK"
  "$<$<VERSION_LESS:${CMAKE_VERSION},3.20>:BLAS::BLAS>"
)

# Custom compile arguments
IF(WITH_GBSA)
add_definitions(-DWITH_GBSA)
ENDIF()

# gfnff-lib sources
set(srcs)
add_subdirectory("src")


###############################################################################
############################### Static Library ################################
###############################################################################
add_library(
  "${PROJECT_NAME}-lib"
  "${srcs}"
)
set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  WINDOWS_EXPORT_ALL_SYMBOLS TRUE
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "${lib-deps}"
)
target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
  file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
endif()

###############################################################################
##################### Install and Export Instructions #########################
###############################################################################
# Export targets for other projects
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")
install(
  TARGETS
  "${PROJECT_NAME}"
  "${PROJECT_NAME}-lib"
  EXPORT
  "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)
# Package license files
install(
  FILES
  "COPYING"
  "COPYING.LESSER"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)


# build a minimal app that can read xyz files and perform singlepoints
if (build_exe)
add_subdirectory("app")
endif()

###############################################################################
################################ python bindings ##############################
###############################################################################
# if (PYTHON_BINDINGS)
# find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
# find_package(pybind11 REQUIRED)
# 
# # Add the pybind11 module
# #pybind11_add_module(_gfnff src/pygfnff.cpp)
# python_add_library(_gfnff src/pygfnff.cpp WITH_SOABI)
# 
# # Link against any existing libraries if necessary
# target_link_libraries(_gfnff 
#                       PRIVATE 
#                       "${PROJECT_NAME}-lib"
#                       pybind11::headers)
# 
# # This is passing in the version as a define
# target_compile_definitions(_gfnff PRIVATE VERSION_INFO=${PROJECT_VERSION})
# 
# # The install directory is the output (wheel) directory
# install(TARGETS _gfnff DESTINATION pygfnff)
# endif()

###############################################################################
################################## Testing ####################################
###############################################################################
# # add the testsuite
# if (WITH_TESTS)
#   enable_testing()
#   add_subdirectory("test")
# endif()

