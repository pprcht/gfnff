# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/gdisp0.f90"
  "${dir}/eg_terms.f90"
  "${dir}/eg_es.f90"
  "${dir}/eg_hb.f90"
  "${dir}/eg_bonded.f90"
  "${dir}/eg_rep.f90"
  "${dir}/driver.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
