set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND ffapp
  "${dir}/testmol.f90"
  "${dir}/main.f90"
)

add_executable(
  "${PROJECT_NAME}-exe"
  "${ffapp}"
)
set_target_properties(
  "${PROJECT_NAME}-exe"
  PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}"
)
target_compile_options(
  "${PROJECT_NAME}-exe"
  PRIVATE
  "$<$<Fortran_COMPILER_ID:GNU>:-fall-intrinsics>"
)
target_link_libraries(
  "${PROJECT_NAME}-exe"
  PRIVATE
  "${PROJECT_NAME}-lib"
)

install(
  TARGETS
  "${PROJECT_NAME}-exe"
  DESTINATION
  "${CMAKE_INSTALL_BINDIR}"
)

add_custom_command(TARGET ${PROJECT_NAME}-exe
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                   $<TARGET_FILE:${PROJECT_NAME}-exe>
                   ${CMAKE_BINARY_DIR}/$<TARGET_FILE_NAME:${PROJECT_NAME}-exe>)

