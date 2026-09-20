# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/dftd3param.f90"
  "${dir}/dftd4param.f90"
  "${dir}/lebedev.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
