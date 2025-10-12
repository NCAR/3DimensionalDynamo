cmake_minimum_required(VERSION 3.5)

project(
 "3D-dynamo"
 LANGUAGES "Fortran" "C"
)

find_package(MPI REQUIRED Fortran)
find_package(superlu REQUIRED C)

set(src)

add_subdirectory(driver)
add_subdirectory(apex)
add_subdirectory(source)

# Create target
add_executable(
  "${PROJECT_NAME}-exe"
  "${src}"
)

target_link_libraries("${PROJECT_NAME}-exe" PRIVATE superlu::superlu)

# Using set_target_properties for each target is usually preferred
set_target_properties(
  "${PROJECT_NAME}-exe"
  PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
