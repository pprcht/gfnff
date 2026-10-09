# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/rab.f90"
  "${dir}/cn.f90"
  "${dir}/qm.f90"
  "${dir}/neighborlist.f90"
  "${dir}/rings.f90"
  "${dir}/predicates.f90"
  "${dir}/bondmat.f90"
  "${dir}/eeq.f90"
  "${dir}/hbset.f90"
  "${dir}/ini_phases.f90"
  "${dir}/ini.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
