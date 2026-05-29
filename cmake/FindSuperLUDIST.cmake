# FindSuperLUDIST.cmake
#
# Finds SuperLU-DIST and creates imported targets:
#   SuperLUDIST::superlu_dist
#   SuperLUDIST::superlu_dist_fortran  (if present)
#
# Hints (set via cmake -D or environment before configuring):
#   SUPERLU_DIST_ROOT          — installation prefix
#   ENV{SUPERLU_DIST_ROOT}
#   ENV{NCAR_ROOT_SUPERLU_DIST}  — set by NCAR module system

include(FindPackageHandleStandardArgs)

# 1. Try pkg-config (works on NCAR after "module load superlu-dist",
#    and on any system where PKG_CONFIG_PATH is set)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_SUPERLU_DIST QUIET superlu_dist)
endif()

# 2. Resolve a root hint from environment if not already set as a cache var
if(NOT SUPERLU_DIST_ROOT)
  if(DEFINED ENV{SUPERLU_DIST_ROOT})
    set(SUPERLU_DIST_ROOT "$ENV{SUPERLU_DIST_ROOT}" CACHE PATH "SuperLU-DIST installation root")
  elseif(DEFINED ENV{NCAR_ROOT_SUPERLU_DIST})
    set(SUPERLU_DIST_ROOT "$ENV{NCAR_ROOT_SUPERLU_DIST}" CACHE PATH "SuperLU-DIST installation root")
  endif()
endif()

# 3. Find header
find_path(SUPERLU_DIST_INCLUDE_DIR
  NAMES superlu_defs.h
  HINTS
    ${PC_SUPERLU_DIST_INCLUDE_DIRS}
    "${SUPERLU_DIST_ROOT}/include")

# 4. Find libraries
find_library(SUPERLU_DIST_LIBRARY
  NAMES superlu_dist
  HINTS
    ${PC_SUPERLU_DIST_LIBRARY_DIRS}
    "${SUPERLU_DIST_ROOT}/lib"
    "${SUPERLU_DIST_ROOT}/lib64")

find_library(SUPERLU_DIST_FORTRAN_LIBRARY
  NAMES superlu_dist_fortran
  HINTS
    ${PC_SUPERLU_DIST_LIBRARY_DIRS}
    "${SUPERLU_DIST_ROOT}/lib"
    "${SUPERLU_DIST_ROOT}/lib64")

# 5. Standard result + version reporting
find_package_handle_standard_args(SuperLUDIST
  REQUIRED_VARS SUPERLU_DIST_LIBRARY SUPERLU_DIST_INCLUDE_DIR
  VERSION_VAR   PC_SUPERLU_DIST_VERSION)

# 6. Create imported targets
if(SuperLUDIST_FOUND AND NOT TARGET SuperLUDIST::superlu_dist)
  add_library(SuperLUDIST::superlu_dist UNKNOWN IMPORTED)
  set_target_properties(SuperLUDIST::superlu_dist PROPERTIES
    IMPORTED_LOCATION             "${SUPERLU_DIST_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${SUPERLU_DIST_INCLUDE_DIR}")

  if(SUPERLU_DIST_FORTRAN_LIBRARY)
    add_library(SuperLUDIST::superlu_dist_fortran UNKNOWN IMPORTED)
    set_target_properties(SuperLUDIST::superlu_dist_fortran PROPERTIES
      IMPORTED_LOCATION             "${SUPERLU_DIST_FORTRAN_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${SUPERLU_DIST_INCLUDE_DIR}")
  endif()
endif()

mark_as_advanced(SUPERLU_DIST_INCLUDE_DIR SUPERLU_DIST_LIBRARY SUPERLU_DIST_FORTRAN_LIBRARY)
