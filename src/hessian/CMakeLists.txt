# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/hess_pair.f90"
  "${dir}/hess_analysis.f90"
  "${dir}/hess_rep.f90"
  "${dir}/hess_prim.f90"
  "${dir}/hess_cn.f90"
  "${dir}/hess_bonded.f90"
  "${dir}/hess_solv.f90"
  "${dir}/hess_es.f90"
  "${dir}/hess_disp.f90"
  "${dir}/hess_hb.f90"
  "${dir}/driver.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
