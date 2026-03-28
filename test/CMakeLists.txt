###########################################################################
# Fortran
###########################################################################

set(
  test-srcs
  "coffeine.f90"
  "supermol.f90"
  "sio2.f90"
  "timer.f90"
  "main.f90"
  "test_gfnff.F90"
  "test_pbc_ini.f90"
  "test_pbc_sp.f90"
)

add_executable(
  "${PROJECT_NAME}-tester"
  "${test-srcs}"
)
target_include_directories(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "${CMAKE_BINARY_DIR}/include"  # Path where the .mod files are located
)
target_link_libraries(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
  "test-drive::test-drive"
)

# Ensure tester is built after the library
add_dependencies("${PROJECT_NAME}-tester" "${PROJECT_NAME}-lib")

add_test("${PROJECT_NAME}/gfnff"   "${PROJECT_NAME}-tester" "gfnff")
add_test("${PROJECT_NAME}/pbc-ini" "${PROJECT_NAME}-tester" "pbc-ini")
add_test("${PROJECT_NAME}/pbc-sp"  "${PROJECT_NAME}-tester" "pbc-sp")


###########################################################################
# C
###########################################################################

set(
  tests
  "gfnff-C-API"
)
set(
  test-srcs-c
  "main.c"
)

add_executable(
  "${PROJECT_NAME}-c-tester"
  "${test-srcs-c}"
)
target_include_directories(
  "${PROJECT_NAME}-c-tester"
  PRIVATE
  "${CMAKE_BINARY_DIR}/include"  # Path where the .mod files are located
)
target_link_libraries(
  "${PROJECT_NAME}-c-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
)


foreach(t IN LISTS tests)
  add_test("${PROJECT_NAME}/${t}" "${PROJECT_NAME}-c-tester" "${t}")
endforeach()


###########################################################################
# C++
###########################################################################

# Specify that C++ is used and set the C++ standard
set(CMAKE_CXX_STANDARD 11)  # You can adjust this to 14, 17, or later as needed
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(
  cpptest
  "gfnff-C++-API"
)
# Add the C++ tester
set(test-srcs-cpp "main.cpp")

add_executable(
  "${PROJECT_NAME}-cpp-tester"
  "${test-srcs-cpp}"
)
set_target_properties("${PROJECT_NAME}-cpp-tester" PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(
  "${PROJECT_NAME}-cpp-tester"
  PRIVATE
  "${CMAKE_BINARY_DIR}/include"  # Include path for headers
)

target_link_libraries(
  "${PROJECT_NAME}-cpp-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
)

add_test("${PROJECT_NAME}/${cpptest}" "${PROJECT_NAME}-cpp-tester")

