# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
#
# See meson.build in this directory for what each subdirectory holds.
# test/tools/checkbuild.py asserts that this list and the meson one agree
# with each other and with the tree.

add_subdirectory("data")
add_subdirectory("core")
add_subdirectory("solvation")
add_subdirectory("topo")
add_subdirectory("engrad")
add_subdirectory("hessian")
add_subdirectory("api")

set(srcs ${srcs} PARENT_SCOPE)
