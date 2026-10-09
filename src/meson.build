# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
#
# Source layout:
#   data/       vendored constant tables (D3/D4 references, Lebedev grids)
#   core/       derived types, parameters, shared math and utilities
#   solvation/  ALPB/GBSA implicit solvation
#   topo/       topology construction: distances, CN, rings, lists, setup
#   engrad/     energy and gradient, one module per term
#   hessian/    second derivatives, one module per term
#   api/        public Fortran and C interfaces
#
# Both this file and CMakeLists.txt list every source by hand;
# test/tools/checkbuild.py asserts that the two agree with each other and
# with the tree, so the lists cannot drift apart again.

subdir('data')
subdir('core')
subdir('solvation')
subdir('topo')
subdir('engrad')
subdir('hessian')
subdir('api')
