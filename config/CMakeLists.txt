
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)
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

option(WITH_OpenMP "Enable support for shared memory parallelisation with OpenMP" TRUE)

option(WITH_GBSA "Enable compilation with GBSA/ALPB implicit solvation routines?" TRUE)

option(build_exe "Build a test binary?" FALSE)


if(NOT DEFINED "${PROJECT_NAME}-dependency-method")
  set(
    "${PROJECT_NAME}-dependency-method"
    "subproject" "cmake" "pkgconf" "fetch"
  )
endif()


#
# To generate gfnff-config.cmake upon a cmake --install
#
include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMinorVersion
)
install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

#
# Compiler settings
#
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#  set(dialect "-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fbacktrace")
  set(dialect "-g -O0 -fbacktrace -ffree-line-length-none -fbacktrace")
  set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#  set(dialect "-axAVX2 -r8 -traceback")
  set(dialect "-g -O2 -r8 -align array64byte -traceback") 
  set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set(dialect "-Mbackslash -Mallocatable=03 -r8 -traceback")
endif()
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}" PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}" PARENT_SCOPE)
