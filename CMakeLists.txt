cmake_minimum_required(VERSION 3.14...3.31)

project(dynamo_driver LANGUAGES Fortran C)

set(CMAKE_VERBOSE_MAKEFILE ON)

# Default to Release if the user did not specify a build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# Make cmake/ visible to find_package
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ── MPI ───────────────────────────────────────────────────────────────────────
find_package(MPI REQUIRED COMPONENTS Fortran)

# ── NetCDF ────────────────────────────────────────────────────────────────────
# On NCAR/Derecho the $NETCDF env var points to the installation prefix;
# on other systems set NETCDF_ROOT on the cmake command line:
#   cmake -DNETCDF_ROOT=/path/to/netcdf ..
if(DEFINED ENV{NETCDF} AND NOT NETCDF_ROOT)
  set(NETCDF_ROOT "$ENV{NETCDF}" CACHE PATH "NetCDF installation root")
endif()

# Try the CMake config package first (present in modern/Spack builds).
# Fall back to find_library + find_path for older or system installs.
find_package(netCDF CONFIG QUIET
  HINTS
    "${NETCDF_ROOT}/lib64/cmake"
    "${NETCDF_ROOT}/lib/cmake"
    "${NETCDF_ROOT}")

if(netCDF_FOUND)
  message(STATUS "Found netCDF via config package: ${netCDF_INSTALL_PREFIX}")
else()
  message(STATUS "netCDF config package not found -- falling back to find_library")

  find_library(NETCDF_C_LIBRARY
    NAMES netcdf
    HINTS "${NETCDF_ROOT}/lib" "${NETCDF_ROOT}/lib64"
    REQUIRED)
  find_path(NETCDF_C_INCLUDE_DIR
    NAMES netcdf.h
    HINTS "${NETCDF_ROOT}/include"
    REQUIRED)

  add_library(netCDF::netcdf UNKNOWN IMPORTED)
  set_target_properties(netCDF::netcdf PROPERTIES
    IMPORTED_LOCATION             "${NETCDF_C_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_C_INCLUDE_DIR}")

  message(STATUS "Found netCDF-C library: ${NETCDF_C_LIBRARY}")
endif()

# NetCDF-Fortran: always found separately (libnetcdff) since netCDFConfig.cmake
# covers only the C library.
find_library(NETCDF_Fortran_LIBRARY
  NAMES netcdff
  HINTS "${NETCDF_ROOT}/lib" "${NETCDF_ROOT}/lib64"
  REQUIRED)

find_path(NETCDF_Fortran_INCLUDE_DIR
  NAMES netcdf.mod
  HINTS "${NETCDF_ROOT}/include"
  REQUIRED)

add_library(netCDF::netcdff UNKNOWN IMPORTED)
set_target_properties(netCDF::netcdff PROPERTIES
  IMPORTED_LOCATION             "${NETCDF_Fortran_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_Fortran_INCLUDE_DIR}")

# ── SuperLU-DIST (no upstream CMake config; uses cmake/FindSuperLUDIST.cmake) ─
# On NCAR/Derecho the module sets PKG_CONFIG_PATH so pkg-config finds it
# automatically.  On other systems set SUPERLU_DIST_ROOT:
#   cmake -DSUPERLU_DIST_ROOT=/path/to/superlu_dist ..
find_package(SuperLUDIST REQUIRED)

# ── Sources ───────────────────────────────────────────────────────────────────
set(src)
add_subdirectory(driver)
add_subdirectory(apex)
add_subdirectory(source)

# ── Executable ────────────────────────────────────────────────────────────────
add_executable(${PROJECT_NAME} ${src})

target_link_libraries(${PROJECT_NAME} PRIVATE
  MPI::MPI_Fortran
  SuperLUDIST::superlu_dist_fortran
  SuperLUDIST::superlu_dist
  netCDF::netcdff
  netCDF::netcdf)

set_target_properties(${PROJECT_NAME} PROPERTIES
  OUTPUT_NAME            ${PROJECT_NAME}
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
