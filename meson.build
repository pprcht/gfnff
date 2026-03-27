# This file is part of gfnff.
# SPDX-Identifier: LGPL-3.0-or-later
#
# gfnff is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gfnff is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with gfnff.  If not, see <https://www.gnu.org/licenses/>.

project(
  'gfnff',
  ['fortran', 'c', 'cpp'],
  version         : '0.0.1',
  license         : 'LGPL-3.0-or-later',
  meson_version   : '>=0.59.0',
  default_options : [
    'buildtype=release',
    'default_library=static',  # gfnff is a dependency library; always build static
    'c_std=gnu11',
    'cpp_std=c++11',
  ],
)

# ── Compiler flags and static link args ───────────────────────────────────────
# Sets: fc, cc, fc_id, cc_id, static_build
subdir('config')

# inc_dirs: 'include/' source dir + build dir (for generated .mod files)
inc_dirs = include_directories('include', '.')

# ═══════════════════════════════════════════════════════════════════════════════
#  OpenMP
# ═══════════════════════════════════════════════════════════════════════════════
# ifx (intel-llvm) triggers a Meson bug in dependency('openmp') because its
# preprocessor output lacks the delimiter strings Meson searches for.
# Work around it by constructing the dependency manually for that compiler.
if fc_id == 'intel-llvm'
  if get_option('openmp')
    omp_dep = declare_dependency(
      compile_args : ['-qopenmp'],
      link_args    : ['-qopenmp'],
    )
  else
    omp_dep = declare_dependency()
  endif
else
  omp_dep = dependency('openmp',
    required : get_option('openmp'),
    language : 'fortran',
  )
endif

_omp_link_dep = (omp_dep.found() and not static_build) ? [omp_dep] : []

# ═══════════════════════════════════════════════════════════════════════════════
#  LAPACK / BLAS
# ═══════════════════════════════════════════════════════════════════════════════
# Same provider logic as crest — fc_id and static_build set by config/meson.build.
# yield:true on the lapack/openmp options means when built as a crest subproject
# these are automatically inherited from the parent, so no duplication occurs.

lapack_opt    = get_option('lapack')
with_mkl      = false
with_openblas = false
lapack_dep    = []
blas_dep      = []

_prefer_mkl = (fc_id in ['intel', 'intel-llvm']) and (lapack_opt == 'auto')

# ── MKL ───────────────────────────────────────────────────────────────────────
if lapack_opt == 'mkl' or _prefer_mkl
  _mkl = dependency('mkl-sdl', required : false, static : static_build)
  if not _mkl.found()
    _thread_lib = (fc_id == 'gcc') ? 'mkl_gnu_thread' : 'mkl_intel_thread'
    _omp_lib    = (fc_id == 'gcc') ? 'gomp'           : 'iomp5'
    _mkl_parts  = []
    foreach lib : ['mkl_intel_lp64', _thread_lib, 'mkl_core', _omp_lib, 'm', 'dl']
      _l = fc.find_library(lib, required : false, static : static_build)
      if _l.found()
        _mkl_parts += [_l]
      endif
    endforeach
    if _mkl_parts.length() > 0
      _mkl = declare_dependency(dependencies : _mkl_parts)
    endif
  endif
  if _mkl.found()
    lapack_dep = _mkl
    blas_dep   = _mkl
    with_mkl   = true
    message('LAPACK/BLAS provider: Intel MKL')
  elif lapack_opt == 'mkl'
    error('MKL requested but not found.')
  else
    lapack_opt = 'auto'
  endif
endif

# ── OpenBLAS ──────────────────────────────────────────────────────────────────
if lapack_opt in ['auto', 'openblas'] and not with_mkl
  _opa = dependency('openblas', required : false, static : static_build)
  if not _opa.found()
    _opa = dependency('blas', required : false, static : static_build)
  endif
  if _opa.found()
    # Check whether this OpenBLAS also bundles LAPACK (dsygvd present)
    _opa_has_lapack = fc.links(
      'external dsygvd; call dsygvd(); end',
      dependencies : _opa,
    )
    lapack_dep    = _opa_has_lapack ? _opa : dependency('lapack',
                      required : false, static : static_build)
    blas_dep      = _opa
    with_openblas = true
    message('LAPACK/BLAS provider: OpenBLAS')
  elif lapack_opt == 'openblas'
    error('OpenBLAS requested but not found.')
  endif
endif

# ── Netlib reference ──────────────────────────────────────────────────────────
if lapack_opt in ['auto', 'netlib'] and not with_mkl and not with_openblas
  _blas   = dependency('blas',   required : false, static : static_build)
  _lapack = dependency('lapack', required : false, static : static_build)
  if _blas.found() and _lapack.found()
    blas_dep   = _blas
    lapack_dep = _lapack
    message('LAPACK/BLAS provider: netlib reference')
  else
    _blas   = fc.find_library('blas',   required : false, static : static_build)
    _lapack = fc.find_library('lapack', required : false, static : static_build)
    if _blas.found() and _lapack.found()
      blas_dep   = declare_dependency(dependencies : [_blas])
      lapack_dep = declare_dependency(dependencies : [_lapack])
    elif lapack_opt == 'netlib'
      error('Netlib BLAS/LAPACK not found.')
    else
      error('No LAPACK/BLAS found. Install libopenblas-dev or set -Dlapack=custom.')
    endif
  endif
endif

# ── Custom ────────────────────────────────────────────────────────────────────
if lapack_opt == 'custom'
  _ll = get_option('lapack_libs')
  _bl = get_option('blas_libs')
  if _ll.length() == 0 and _bl.length() == 0
    error('lapack=custom requires lapack_libs and/or blas_libs.')
  endif
  _deps = []
  foreach lib : _ll + _bl
    _deps += [fc.find_library(lib, required : true, static : static_build)]
  endforeach
  lapack_dep = declare_dependency(dependencies : _deps)
  blas_dep   = lapack_dep
endif

if with_mkl
  add_project_arguments('-DWITH_MKL', language : ['c', 'fortran'])
endif
if with_openblas
  add_project_arguments('-DWITH_OPENBLAS', language : ['c', 'fortran'])
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Source collection
# ═══════════════════════════════════════════════════════════════════════════════
srcs = []
prog = []
subdir('src')

# ═══════════════════════════════════════════════════════════════════════════════
#  Library target
# ═══════════════════════════════════════════════════════════════════════════════
_lib_deps = [lapack_dep, blas_dep] + _omp_link_dep

gfnff_lib = library(
  meson.project_name(),
  sources              : srcs,
  include_directories  : inc_dirs,
  dependencies         : _lib_deps,
  install              : not (meson.is_subproject() and get_option('default_library') == 'static'),
)

# ── Exported dependency (what crest links against) ────────────────────────────
gfnff_dep = declare_dependency(
  link_with           : gfnff_lib,
  include_directories : [inc_dirs, gfnff_lib.private_dir_include()],
  dependencies        : _lib_deps,
)

# ═══════════════════════════════════════════════════════════════════════════════
#  Optional standalone app
# ═══════════════════════════════════════════════════════════════════════════════
if get_option('build_exe')
  subdir('app')
  executable(
    'gfnff',
    sources      : prog,
    dependencies : gfnff_dep,
    install      : true,
    link_language: 'fortran',
  )
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Unit tests
# ═══════════════════════════════════════════════════════════════════════════════
if get_option('tests')
  testdrive_dep = dependency('test-drive',
    version  : '>=0.4.0',
    fallback : ['test-drive', 'testdrive_dep'],
    required : true,
    static   : static_build,
  )
  subdir('test')
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Install
# ═══════════════════════════════════════════════════════════════════════════════
install_data(
  'LICENSE',
  install_dir : get_option('datadir') / 'licenses' / meson.project_name(),
)

# ═══════════════════════════════════════════════════════════════════════════════
#  Summary
# ═══════════════════════════════════════════════════════════════════════════════
summary({
  'Fortran compiler' : '@0@ @1@'.format(fc_id, fc.version()),
  'C compiler'       : '@0@ @1@'.format(cc_id, cc.version()),
  'Build type'       : get_option('buildtype'),
  'Static binary'    : static_build,
  'OpenMP'           : omp_dep.found(),
  'MKL'              : with_mkl,
  'OpenBLAS'         : with_openblas,
  'Standalone app'   : get_option('build_exe'),
  'Unit tests'       : get_option('tests'),
}, section : 'gfnff build configuration')
