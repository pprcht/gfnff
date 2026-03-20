# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
#
# Included from the root meson.build via subdir('config').
# Sets: fc, cc, fc_id, cc_id, static_build

# ═══════════════════════════════════════════════════════════════════════════════
#  Compiler identification
# ═══════════════════════════════════════════════════════════════════════════════
fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')
fc_id = fc.get_id()
cc_id = cc.get_id()

if fc_id != cc_id
  warning('Fortran and C compilers are from different vendors (@0@ vs @1@). '
    .format(fc_id, cc_id)
    + 'See config/intel-fortran-gnu-c.ini for guidance on OpenMP runtime ownership.')
endif

if fc_id not in ['gcc', 'intel', 'intel-llvm']
  warning('Unsupported Fortran compiler "@0@". Proceeding, but results may be wrong.'.format(fc_id))
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Compiler flags
# ═══════════════════════════════════════════════════════════════════════════════
_is_debug = get_option('buildtype') in ['debug', 'debugoptimized']

# --- GNU (gfortran / gcc) -----------------------------------------------------
if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language : 'fortran',
  )
  if _is_debug
    add_project_arguments(
      '-fcheck=all',
      '-ffpe-trap=invalid,zero,overflow',
      '-fbounds-check',
      '-finit-real=snan',
      '-finit-integer=-999',
      '-Wall',
      '-Wextra',
      language : 'fortran',
    )
  endif
endif

# --- Intel classic (ifort / icc) -----------------------------------------------
if fc_id == 'intel'
  add_project_arguments(
    '-r8',
    '-align', 'array64byte',
    '-traceback',
    '-Wl,--allow-multiple-definition',
    language : 'fortran',
  )
  add_project_arguments('-DLINUX', language : 'c')
  if _is_debug
    add_project_arguments('-check', 'all', '-fpe0', language : 'fortran')
  endif
endif

# --- Intel LLVM (ifx / icx) ----------------------------------------------------
if fc_id == 'intel-llvm'
  add_project_arguments(
    '-r8',
    '-align', 'array64byte',
    '-traceback',
    language : 'fortran',
  )
  if _is_debug
    add_project_arguments('-check', 'all', '-fpe0', language : 'fortran')
  endif
endif

# Needed for some C files on Linux (float128 compatibility)
add_project_arguments('-D_Float128=__float128', language : 'c')

# ═══════════════════════════════════════════════════════════════════════════════
#  Static-binary link arguments
# ═══════════════════════════════════════════════════════════════════════════════
static_build = get_option('static')

if static_build
  if fc_id == 'gcc'
    add_project_link_arguments('-static', language : 'fortran')
    add_project_link_arguments('-static', language : 'c')
  elif fc_id == 'intel'
    add_project_link_arguments('-static-intel', '-qopenmp-link=static', language : 'fortran')
  elif fc_id == 'intel-llvm'
    add_project_link_arguments('-static-intel', '-qopenmp-link=static', language : 'fortran')
  endif
endif
