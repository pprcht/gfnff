
fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')
fc_id = fc.get_id()
cc_id = cc.get_id()

if fc.get_id() != cc.get_id()
  warning('FC and CC are not from the same vendor')
endif

## ========================================= ##
## Compiler specific default arguments
## ========================================= ##  
if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    '-Wno-maybe-uninitialized',
    '-Wno-uninitialized',
    '-Wno-unused-variable',
    '-Wno-unused-dummy-argument',
    language: 'fortran',
  )
elif fc_id == 'intel'
  add_project_link_arguments(
    '-Wl,--allow-multiple-definition',
    language: 'fortran',
  )
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
  add_project_arguments(
    '-DLINUX',
    language: 'c',
  )
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
elif fc_id == 'flang'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    language: 'fortran',
  )
endif

add_project_arguments('-D_Float128=__float128', language: 'c')

## ========================================= ##
## OpenMP
## ========================================= ##
if get_option('openmp')
  omp_dep = dependency('openmp', required: fc.get_id() != 'intel')
  if not omp_dep.found()
    omp_dep = declare_dependency(
      compile_args: '-qopenmp',
      link_args: '-fopenmp',
    )
  endif
  exe_deps += omp_dep
endif

## ======================================== ##
## Linear Algebra Libraries
## ======================================== ##
la_backend = get_option('la_backend')
if la_backend == 'mkl' or la_backend == 'mkl-static'
  add_project_arguments('-DWITH_MKL', language: 'fortran')
  if la_backend == 'mkl-static'
    add_project_link_arguments('-static', language: 'fortran')
    add_project_link_arguments('-static', language: 'c') # icc will do linking
  endif

  if get_option('default_library') == 'shared'
    mkl_rt_dep = cc.find_library('mkl_rt', required: true)
    exe_deps += mkl_rt_dep
  else
    if fc.get_id() == 'gcc'
      libmkl_exe = [cc.find_library('mkl_gf_lp64')]
      if get_option('openmp')
        libmkl_exe += cc.find_library('mkl_gnu_thread')
      endif
    elif fc.get_id() == 'intel' or fc.get_id() == 'intel-cl'
      libmkl_exe = [cc.find_library('mkl_intel_lp64')]
      if get_option('openmp')
        libmkl_exe += cc.find_library('mkl_intel_thread')
      endif
    elif fc.get_id() == 'pgi' or fc.get_id() == 'nvidia_hpc'
      libmkl_exe = [cc.find_library('mkl_intel_lp64')]
      if get_option('openmp')
        libmkl_exe += cc.find_library('mkl_pgi_thread')
      endif
    endif
    if not get_option('openmp')
      libmkl_exe += cc.find_library('mkl_sequential')
    endif
    libmkl_exe += cc.find_library('mkl_core')
    exe_deps += libmkl_exe
  endif

elif la_backend == 'mkl-rt'
  add_project_arguments('-DWITH_MKL', language: 'fortran')

  mkl_rt_dep = fc.find_library('mkl_rt', required: true)
  if fc.get_id() == 'intel'
    exe_deps += fc.find_library('ifcore')
  endif

  exe_deps += mkl_rt_dep

elif la_backend == 'openblas'
  # search for OpenBLAS
  blas_dep = dependency('openblas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('openblas', required: true)
  endif
  exe_deps += blas_dep
  # some OpenBLAS versions can provide lapack, check if we can find dsygvd
  openblas_provides_lapack = fc.links(
    'external dsygvd; call dsygvd(); end',
    dependencies: blas_dep,
  )
  # otherwise we fall back to LAPACK
  if not openblas_provides_lapack
    lapack_dep = dependency('lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack', required: true)
    endif
    exe_deps += lapack_dep
  endif

elif la_backend == 'custom'
  foreach lib: get_option('custom_libraries')
    exe_deps += fc.find_library(lib)
  endforeach

else
  # Find BLAS (usually netlib, but in conda also OpenBLAS/MKL)
  blas_dep = dependency('blas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('blas', required: true)
  endif
  exe_deps += blas_dep
  # Find LAPACK (usually netlib, but in conda also MKL)
  lapack_dep = dependency('lapack', required: false)
  if not lapack_dep.found()
    lapack_dep = fc.find_library('lapack', required: true)
  endif
  exe_deps += lapack_dep
endif


## ========================================= ##
## External subprojects
## ========================================= ##  

if get_option('with_gbsa')
add_project_arguments('-DWITH_GBSA', language: 'fortran')
#subdir('gbsa')
endif

