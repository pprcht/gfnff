
set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs2
  "${dir}/testmol.f90"
  "${dir}/main.f90"
)

set(srcs2 ${srcs2} PARENT_SCOPE)
