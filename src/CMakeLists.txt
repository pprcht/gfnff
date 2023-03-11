
IF(WITH_GBSA)
add_subdirectory("gbsa")
ENDIF()

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
   "${dir}/gfnff_math_wrapper.f90"
   "${dir}/gfnff_data_types.f90"
   "${dir}/gfnff_param.f90"
   "${dir}/dftd3param.f90"
   "${dir}/dftd4param.f90"
   "${dir}/gfnff_helpers.f90"
   "${dir}/gfnff_fraghess.f90"
   "${dir}/gfnff_rab.f90"
   "${dir}/gfnff_cn.f90"
   "${dir}/gfnff_gdisp0.f90"
   "${dir}/gfnff_qm.f90"
   "${dir}/gfnff_ini2.f90"
   "${dir}/gfnff_ini.f90"
   "${dir}/gfnff_alpb.F90"
   "${dir}/gfnff_engrad.f90"
   "${dir}/gfnff_restart.f90"
   "${dir}/gfnff_setup.f90"
   "${dir}/gfnff_interface.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
