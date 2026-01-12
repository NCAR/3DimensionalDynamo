cmake_minimum_required(VERSION 3.5)

project(
 "dynamo_driver"
 LANGUAGES "Fortran" "C"
)

set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(MPI REQUIRED Fortran)
find_package(superlu REQUIRED C)

set(src)

add_subdirectory(driver)
add_subdirectory(apex)
add_subdirectory(source)

# Create target
add_executable(
  "${PROJECT_NAME}"
  "${src}"
)

target_include_directories("${PROJECT_NAME}" PRIVATE $ENV{NETCDF}/include)
target_link_directories("${PROJECT_NAME}" PRIVATE $ENV{NETCDF}/lib)
target_link_libraries("${PROJECT_NAME}" PRIVATE superlu::superlu netcdff)

# Using set_target_properties for each target is usually preferred
set_target_properties(
  "${PROJECT_NAME}"
  PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
