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

cmake_minimum_required(VERSION 3.21)

project(
  "gfnff"
  LANGUAGES "Fortran" "C" "CXX"
  VERSION "0.2.0"
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
    # Do not use REQUIRED: on macOS the C/CXX compiler (Apple Clang) has no
    # native OpenMP support, so find_package(OpenMP REQUIRED) would fail for
    # those languages even though gfortran supports it fine.  We only need
    # Fortran-level OpenMP; check OpenMP_Fortran_FOUND after the search.
    find_package("OpenMP")
    if(OpenMP_Fortran_FOUND)
      add_compile_definitions(WITH_OpenMP)
    else()
      message(STATUS "OpenMP not found for Fortran; building without OpenMP parallelization")
      add_library("OpenMP::OpenMP_Fortran" INTERFACE IMPORTED)
    endif()
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

# TOML parameter file support via toml-f.
#
# Optional, like the meson build: without it the library still builds and the
# reader/writer report that support was not compiled in.  It is ON by default
# here because the shared library built from this file is what the Python
# bindings load, and a Python caller passing parametrisation="...toml" to a
# library built without toml-f gets a runtime failure with no way to fix it.
# An already-present toml-f (e.g. from a parent project) is used as-is;
# otherwise it is fetched at the same revision the meson wrap pins.
option(WITH_TOMLF "TOML parameter file support via toml-f" ON)
set(GFNFF_TOMLF_REV "28f4601fe0992ed17a6b97e230a61289f56bcd8e"
    CACHE STRING "toml-f revision to fetch; keep in step with subprojects/toml-f.wrap")
set(tomlf-dep)
if(WITH_TOMLF)
  if(NOT TARGET "toml-f::toml-f" AND NOT TARGET "toml-f")
    find_package("toml-f" QUIET)
    if(NOT toml-f_FOUND)
      include(FetchContent)
      set(WITH_TESTS_SAVED "${WITH_TESTS}")
      set(WITH_TESTS OFF CACHE BOOL "" FORCE)  # do not build toml-f's own tests
      FetchContent_Declare(
        "toml-f"
        GIT_REPOSITORY "https://github.com/toml-f/toml-f"
        GIT_TAG "${GFNFF_TOMLF_REV}"
        GIT_SHALLOW FALSE  # a pinned SHA cannot be fetched shallowly
      )
      FetchContent_MakeAvailable("toml-f")
      set(WITH_TESTS "${WITH_TESTS_SAVED}" CACHE BOOL "" FORCE)
    endif()
  endif()
  if(TARGET "toml-f::toml-f")
    set(tomlf-dep "toml-f::toml-f")
  elseif(TARGET "toml-f-lib")
    set(tomlf-dep "toml-f-lib")
  elseif(TARGET "toml-f")
    set(tomlf-dep "toml-f")
  else()
    message(FATAL_ERROR "WITH_TOMLF=ON but no toml-f target was produced")
  endif()
  add_compile_definitions(WITH_TOMLF)
endif()

set(
  lib-deps
  "OpenMP::OpenMP_Fortran"
  "LAPACK::LAPACK"
  "${tomlf-dep}"
  "$<$<VERSION_LESS:${CMAKE_VERSION},3.20>:BLAS::BLAS>"
)

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
  "LICENSE"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)


# build a minimal app that can read xyz files and perform singlepoints
if (build_exe)
add_subdirectory("app")
endif()

###############################################################################
################################ python bindings ##############################
###############################################################################
# Shared library for Python ctypes bindings.
# Enable with -DPYTHON_BINDINGS=ON (scikit-build-core passes this automatically).
# The static gfnff-lib target is kept for downstream Fortran/CMake consumers.
if(PYTHON_BINDINGS)
  add_library(
    "${PROJECT_NAME}-shared" SHARED
    "${srcs}"
  )
  set_target_properties(
    "${PROJECT_NAME}-shared"
    PROPERTIES
    OUTPUT_NAME "${PROJECT_NAME}"
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION_MAJOR}"
    WINDOWS_EXPORT_ALL_SYMBOLS TRUE
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include-shared"
  )
  target_link_libraries(
    "${PROJECT_NAME}-shared"
    PUBLIC
    "${lib-deps}"
  )
  target_include_directories(
    "${PROJECT_NAME}-shared"
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include-shared>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
  # Install the shared library so scikit-build-core bundles it into the wheel.
  install(
    TARGETS "${PROJECT_NAME}-shared"
    LIBRARY DESTINATION "gfnff"
    RUNTIME DESTINATION "gfnff"
  )
  # After every build, also copy the shared library into the Python source tree
  # so that editable installs work without LD_LIBRARY_PATH.
  add_custom_command(
    TARGET "${PROJECT_NAME}-shared" POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      $<TARGET_FILE:${PROJECT_NAME}-shared>
      "${CMAKE_CURRENT_SOURCE_DIR}/python/gfnff/$<TARGET_FILE_NAME:${PROJECT_NAME}-shared>"
    COMMENT "Copying libgfnff to python/gfnff/ for development use"
  )
endif()

###############################################################################
################################## Testing ####################################
###############################################################################
# add the testsuite
if (WITH_TESTS)
  enable_testing()
  add_subdirectory("test")
endif()

