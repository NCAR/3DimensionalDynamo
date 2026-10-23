cmake_minimum_required(VERSION 3.18)

project(
 "dynamo_driver"
 LANGUAGES "Fortran" "C"
)

set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(superlu REQUIRED)

# netCDF-Fortran: try pkg-config first, then fall back to the NETCDF env var
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(NETCDFF QUIET IMPORTED_TARGET netcdf-fortran)
endif()
if(TARGET PkgConfig::NETCDFF)
  set(netcdf_lib PkgConfig::NETCDFF)
elseif(DEFINED ENV{NETCDF})
  add_library(netcdff_ext INTERFACE)
  target_include_directories(netcdff_ext INTERFACE $ENV{NETCDF}/include)
  target_link_directories(netcdff_ext INTERFACE $ENV{NETCDF}/lib)
  target_link_libraries(netcdff_ext INTERFACE netcdff netcdf)
  set(netcdf_lib netcdff_ext)
else()
  message(FATAL_ERROR "netCDF-Fortran not found: set PKG_CONFIG_PATH or NETCDF")
endif()

set(src)

add_subdirectory(driver)
add_subdirectory(apex)
add_subdirectory(source)

# Create target
add_executable(
  "${PROJECT_NAME}"
  "${src}"
)

target_link_libraries("${PROJECT_NAME}" PRIVATE MPI::MPI_Fortran superlu::superlu ${netcdf_lib})

# Using set_target_properties for each target is usually preferred
set_target_properties(
  "${PROJECT_NAME}"
  PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
