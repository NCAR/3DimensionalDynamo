set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND src
  "${dir}/dynamo_interface_mod.F90"
  "${dir}/alloc_mod.F90"
  "${dir}/calculate_terms_mod.F90"
  "${dir}/cons_mod.F90"
  "${dir}/fieldline_mod.F90"
  "${dir}/grid_mod.F90"
  "${dir}/init_mod.F90"
  "${dir}/mpi_mod.F90"
  "${dir}/params_mod.F90"
  "${dir}/prec.F90"
  "${dir}/solver_mod.F90"
  "${dir}/stencil_mod.F90"
  "${dir}/util_mod.F90"
  "${dir}/c_fortran_dgssv.c"
)

set(src "${src}" PARENT_SCOPE)
