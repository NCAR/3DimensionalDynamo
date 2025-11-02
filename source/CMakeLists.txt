set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND src
  "${dir}/dynamo_interface_mod.F90"
  "${dir}/alloc_module.F90"
  "${dir}/calculate_terms_module.F90"
  "${dir}/cons_module.F90"
  "${dir}/fieldline_module.F90"
  "${dir}/grid_module.F90"
  "${dir}/init_module.F90"
  "${dir}/mpi_module.F90"
  "${dir}/params_module.F90"
  "${dir}/prec.F90"
  "${dir}/solver_module.F90"
  "${dir}/stencil_module.F90"
  "${dir}/util_module.F90"
  "${dir}/c_fortran_dgssv.c"
)

set(src "${src}" PARENT_SCOPE)
