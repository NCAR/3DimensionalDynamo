program main
  use iso_fortran_env, only: r8=>real64 ! double precision
  use dynamo_interface_mod, only: dynamo_init1, dynamo_init2, dynamo_calc
  use apex, only: apex_setup

  use inputdata_mod, only: inputdata_init, inputdata_close, npflpts1, npflpts2, ntimes
  use inputdata_mod, only: inputdata_read_s1_fld, inputdata_read_s2_fld
  use inputdata_mod, only: un_s1_vid, un_s2_vid, vn_s1_vid, vn_s2_vid
  use inputdata_mod, only: sigma_hal_s1_vid, sigma_hal_s2_vid, sigma_ped_s1_vid, sigma_ped_s2_vid
  use mpi_module, only: mpi_rank, mpi_size, lat_size, lon_size
  use params_module,only: nmlat_h, nmlon, nhgt_fix, hgt_fix_r
  use indices_mod, only : indices_init
  use mpi_module, only: mlat0, mlat1, mlon0, mlon1


  implicit none

  include 'mpif.h'

  integer, parameter :: edyn3d_nmlon = 180
  integer, parameter :: edyn3d_nhgt = 54
  integer, parameter :: edyn3d_nmlat_h = 91

  integer, parameter :: nglats = 90
  integer, parameter :: nglons = 180

  real(r8), parameter :: geomag_year = 2020.5

  character(len=*), parameter :: prefix = '3D-Dynamamo test: '

  integer :: ierr, job_size, my_rank, itime

  real(r8), allocatable :: un_s1(:,:,:,:)
  real(r8), allocatable :: un_s2(:,:,:,:)
  real(r8), allocatable :: vn_s1(:,:,:,:)
  real(r8), allocatable :: vn_s2(:,:,:,:)

  real(r8), allocatable :: sigma_hal_s1(:,:,:,:)
  real(r8), allocatable :: sigma_hal_s2(:,:,:,:)
  real(r8), allocatable :: sigma_ped_s1(:,:,:,:)
  real(r8), allocatable :: sigma_ped_s2(:,:,:,:)

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

  write(*,*) prefix,'call inputdata_read ...'
  call inputdata_init('../data/FX2000_ne16pg3_wcmx_3Dedyn_test02.cam.h5i.0001-01-01-00000.nc')

  call indices_init(npflpts1,npflpts2)

  allocate( un_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  un_s1 = 0._r8
  allocate( vn_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  vn_s1 = 0._r8
  allocate( un_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  un_s2 = 0._r8
  allocate( vn_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  vn_s2 = 0._r8

  allocate( sigma_hal_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  sigma_hal_s1 = 0._r8
  allocate( sigma_ped_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  sigma_ped_s1 = 0._r8

  allocate( sigma_hal_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  sigma_hal_s2 = 0._r8
  allocate( sigma_ped_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  sigma_ped_s2 = 0._r8

  do itime = 1,ntimes
     call inputdata_read_s1_fld( varid=un_s1_vid, time_ndx=itime, fld=un_s1 )
     call inputdata_read_s1_fld( varid=vn_s1_vid, time_ndx=itime, fld=vn_s1 )
     call inputdata_read_s1_fld( varid=sigma_hal_s1_vid, time_ndx=itime, fld=sigma_hal_s1 )
     call inputdata_read_s1_fld( varid=sigma_ped_s1_vid, time_ndx=itime, fld=sigma_ped_s1 )
     call inputdata_read_s2_fld( varid=un_s2_vid, time_ndx=itime, fld=un_s2 )
     call inputdata_read_s2_fld( varid=vn_s2_vid, time_ndx=itime, fld=vn_s2 )
     call inputdata_read_s2_fld( varid=sigma_hal_s2_vid, time_ndx=itime, fld=sigma_hal_s2 )
     call inputdata_read_s2_fld( varid=sigma_ped_s2_vid, time_ndx=itime, fld=sigma_ped_s2 )

     write(*,*) prefix,' ... itime:',itime,' maxval( sigma_ped_s2 ): ',maxval( sigma_ped_s2 )
  end do


  call inputdata_close()
  write(*,*) prefix,'DONE...'

  call mpi_finalize(ierr)
  print*,'... mpi_finalize ierr: ',ierr

  print*,'END TEST.'

end program main
