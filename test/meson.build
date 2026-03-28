# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later

# ── Fortran test (test-drive) ──────────────────────────────────────────────────
gfnff_tester = executable(
  'gfnff-tester',
  sources : files(
    'coffeine.f90',
    'supermol.f90',
    'sio2.f90',
    'timer.f90',
    'main.f90',
    'test_gfnff.F90',
    'test_pbc_ini.f90',
    'test_pbc_sp.f90',
  ),
  dependencies : [gfnff_dep, testdrive_dep],
)

test('gfnff/gfnff',   gfnff_tester, args : ['gfnff'],   timeout : 120)
test('gfnff/pbc-ini', gfnff_tester, args : ['pbc-ini'], timeout : 120)
test('gfnff/pbc-sp',  gfnff_tester, args : ['pbc-sp'],  timeout : 120)

# ── C API test ────────────────────────────────────────────────────────────────
gfnff_c_tester = executable(
  'gfnff-c-tester',
  sources : files('main.c'),
  dependencies : gfnff_dep,
)

test('gfnff/gfnff-C-API', gfnff_c_tester, timeout : 120)

# ── C++ API test ──────────────────────────────────────────────────────────────
gfnff_cpp_tester = executable(
  'gfnff-cpp-tester',
  sources      : files('main.cpp'),
  dependencies : gfnff_dep,
  # No link_language override: meson sees .cpp sources and picks g++ automatically,
  # which pulls in libstdc++. libgfortran comes in transitively via gfnff_dep.
)

test('gfnff/gfnff-C++-API', gfnff_cpp_tester, timeout : 120)
