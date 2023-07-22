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

cmake_minimum_required(VERSION 3.14)

project(
  "gfnff"
  LANGUAGES "Fortran" "C"
  VERSION "0.0.1"
)


# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# General configuration information
set(libs)
add_subdirectory("config")


#
# Libraries
#

# OpenMP dependency
if(NOT TARGET OpenMP::OpenMP_Fortran AND WITH_OpenMP)
  find_package(OpenMP REQUIRED)
  list(
    APPEND libs
    OpenMP::OpenMP_Fortran
  )
endif()

# BLAS and LAPACK
if(NOT TARGET BLAS::BLAS)
  find_package(BLAS REQUIRED)
  if(NOT TARGET BLAS::BLAS AND BLAS_FOUND)
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    target_link_libraries(BLAS::BLAS INTERFACE "${BLAS_LIBRARIES}")
    target_link_options(BLAS::BLAS INTERFACE "${BLAS_LINKER_FLAGS}")
  endif()
endif()
if(NOT TARGET LAPACK::LAPACK)
  find_package(LAPACK REQUIRED)
  if(NOT TARGET LAPACK::LAPACK AND LAPACK_FOUND)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    target_link_libraries(LAPACK::LAPACK INTERFACE "${LAPACK_LIBRARIES}")
    target_link_options(LAPACK::LAPACK INTERFACE "${LAPACK_LINKER_FLAGS}")
  endif()
endif()
list(
  APPEND libs
  LAPACK::LAPACK
  BLAS::BLAS
)


#
# gfnff sources
#
set(srcs)
add_subdirectory("src")


#
# Build the lib
#
add_library(
  "lib${PROJECT_NAME}.a"
  "${srcs}"
)
set_target_properties(
  "lib${PROJECT_NAME}.a"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  WINDOWS_EXPORT_ALL_SYMBOLS TRUE
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_include_directories(
  "lib${PROJECT_NAME}.a"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
  "lib${PROJECT_NAME}.a"
  PUBLIC
  "${libs}"
)

IF(WITH_GBSA)
add_definitions(-DWITH_GBSA)
ENDIF()

#
# Executables
#
IF(build_exe)
set(srcs2)
add_subdirectory("testprog")
add_executable(
  ${PROJECT_NAME}-tester
  "${srcs2}" 
)
set_target_properties(
  ${PROJECT_NAME}-tester
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  OUTPUT_NAME "${PROJECT_NAME}-tester"
)
target_link_libraries(
  ${PROJECT_NAME}-tester
  PRIVATE
  "lib${PROJECT_NAME}.a"
)
ENDIF()

#
# Export targets for other projects
#
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE   "lib${PROJECT_NAME}.a")
install(
  TARGETS
  "${PROJECT_NAME}"
  "lib${PROJECT_NAME}.a"
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


