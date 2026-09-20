# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/restart.f90"
  "${dir}/setup.f90"
  "${dir}/interface.f90"
  "${dir}/c_interface.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
