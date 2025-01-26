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
#
# You should have received a copy of the GNU Lesser General Public License
# along with gfnff.  If not, see <https://www.gnu.org/licenses/>.
# Set the module path for CMake includes

#########################################################################################
#########################################################################################

# Add modules to the CMake build 
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)


# specify module installation directory
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/modules/"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
set(
  module-dir
  "${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}"
)
set(module-dir "${module-dir}" PARENT_SCOPE)


# Custom options
option(WITH_OpenMP "Enable OpenMP support"    TRUE)
option(WITH_GBSA   "Build with implicit solvation support" TRUE)
option(WITH_TESTS  "Enable test binary build" TRUE)
option(build_exe   "Build an executable"      FALSE)
# For python bindings install scikit-build and pybind11
option(PYTHON_BINDINGS  "Create Python bindings"  FALSE)

#########################################################################################
#########################################################################################

# Set build type as CMake does not provide defaults
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(
    CMAKE_BUILD_TYPE "RelWithDebInfo"
    CACHE STRING "Build type to be used."
    FORCE
  )
  message(
    STATUS
    "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified."
  )
endif()


### Options and defaults
include("${CMAKE_CURRENT_SOURCE_DIR}/modules/gfnff-utils.cmake")
set(fortran_minimal_versions "GNU;7.5" "Intel;19.0")
check_minimal_compiler_version("Fortran" "${fortran_minimal_versions}")

if(NOT DEFINED "${PROJECT_NAME}-dependency-method")
  set(
    "${PROJECT_NAME}-dependency-method"
    "subproject" "cmake" "pkgconf" "fetch"
  )
endif()

#########################################################################################
#########################################################################################

# Compiler settings for GNU and Intel Fortran compilers
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-g -O0 -fbacktrace -ffree-line-length-none -fbacktrace")
  set(bounds "-fbounds-check -ffpe-trap=invalid,zero,overflow")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-g -O2 -r8 -align array64byte -traceback")
  set(bounds "-check all -fpe0")
else()
  message(FATAL_ERROR "Please use an Intel or GNU compiler!")
endif()

# Apply the compiler flags
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}" PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}" PARENT_SCOPE)

#########################################################################################
#########################################################################################

