program main
  use iso_fortran_env, only: r8=>real64 ! double precision
  use dynamo_interface_mod
  use apex
  use mpi_module, only: mpi_rank, mpi_size, lat_size, lon_size
  use params_module,only: nmlat_h, nmlon, nhgt_fix, hgt_fix_r

  implicit none

  include 'mpif.h'

  integer, parameter :: edyn3d_nmlon = 12
  integer, parameter :: edyn3d_nhgt = 45
  integer, parameter :: edyn3d_nmlat_h = 90

  integer, parameter :: nglats = 90
  integer, parameter :: nglons = 180

  real(r8), parameter :: geomag_year = 2020.5

  character(len=*), parameter :: prefix = '3D-Dynamamo test: '

  integer :: ierr, job_size, my_rank
  print*,'BEGIN TEST...'

  call mpi_init(ierr)
  print*,'... mpi_init ierr: ',ierr, MPI_ERROR,' ERROR? ',ierr==MPI_ERROR

  call mpi_comm_size(MPI_COMM_WORLD, job_size, ierr)

  call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)

  print*,' my rank : ',my_rank,' of ',job_size

  write(*,*) prefix,'call dynamo_init1 ...'
  ! before apex init
  call dynamo_init1( &
       set_hilat_pot_in=.true., &
       set_hilat_fac_in=.false., &
       mpicom_atm=MPI_COMM_WORLD, npes_edyn3D=job_size, &
       edyn3d_nmlat_h=edyn3d_nmlat_h, edyn3d_nmlon=edyn3d_nmlon, edyn3d_nhgt=edyn3d_nhgt, &
       real_kind=r8 )

  write(*,*) prefix,'3D Edyn grid params nmlat_h,nmlon,nhgt_fix: ',nmlat_h,nmlon,nhgt_fix
  write(*,*) prefix,'3D Edyn mpi_rank, mpi_size: ',mpi_rank,mpi_size
  write(*,*) prefix,'3D Edyn lon_size, lat_size: ',lon_size,lat_size


  ! initialize APEX -- after dynamo_init1 and before dynamo_init2
  ! call mo_apex_init1( alts_in=hgt_fix_r*1.e-3_r8 ) ! m --> km

  write(*,*) prefix,'call APEX setup ...'
  call apex_setup(date=2010.2_r8 ,altmax=maxval(hgt_fix_r*1.e-3_r8))

  write(*,*) prefix,'call dynamo_init2 ...'
  ! after apex init
  call dynamo_init2()

  write(*,*) prefix,'DONE...'

  call mpi_finalize(ierr)
  print*,'... mpi_finalize ierr: ',ierr

  print*,'END TEST.'

end program main
