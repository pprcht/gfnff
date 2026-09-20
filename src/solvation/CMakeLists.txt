# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/type.f90"
  "${dir}/state.f90"
  "${dir}/input.F90"
  "${dir}/util.F90"
  "${dir}/vdwrad.f90"
  "${dir}/kernel.f90"
  "${dir}/born.f90"
  "${dir}/sasa.f90"
  "${dir}/gbsa.f90"
  "${dir}/model.f90"
  "${dir}/alpb.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
