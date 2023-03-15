set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/born.f90"
  "${dir}/bisectsearch.f90"
  "${dir}/cm5.f90"
  "${dir}/gbsa.f90"
  "${dir}/input.F90"
  "${dir}/kernel.f90"
  "${dir}/lebedev.f90"
  "${dir}/model.f90"
  "${dir}/sasa.f90"
  "${dir}/solvation_math_wrapper.f90"
  "${dir}/state.f90"
  "${dir}/type.f90"
  "${dir}/util.F90"
  "${dir}/vdwrad.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
