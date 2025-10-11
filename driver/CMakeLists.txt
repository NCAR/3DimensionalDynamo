set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND src
  "${dir}/main.F90"
)

set(src "${src}" PARENT_SCOPE)
