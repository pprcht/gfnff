# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/math.f90"
  "${dir}/types.f90"
  "${dir}/neighbor.f90"
  "${dir}/param_tables.f90"
  "${dir}/bond_potential.f90"
  "${dir}/param.f90"
  "${dir}/param_io.F90"
  "${dir}/geometry.f90"
  "${dir}/fragments.f90"
  "${dir}/strings.f90"
  "${dir}/search.f90"
  "${dir}/wsc.f90"
  "${dir}/latticepoint.f90"
  "${dir}/timer.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
