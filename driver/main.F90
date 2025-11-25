program main
  use iso_fortran_env, only: r8=>real64 ! double precision
  use dynamo_interface_mod, only: dynamo_init1, dynamo_init2, dynamo_calc, dynamo_final
  use apex, only: apex_setup

  use inputdata_mod, only: inputdata_init, inputdata_close, npflpts1, npflpts2, ntimes
  use inputdata_mod, only: inputdata_read_s1_fld, inputdata_read_s2_fld, inputdata_read_2dp_fld
  use inputdata_mod, only: un_s1_vid, un_s2_vid, vn_s1_vid, vn_s2_vid
  use inputdata_mod, only: sigma_hal_s1_vid, sigma_hal_s2_vid, sigma_ped_s1_vid, sigma_ped_s2_vid
  use inputdata_mod, only: hilat_pot_vid
  use mpi_module, only: mpi_rank, mpi_size, lat_size, lon_size
  use params_module,only: nmlat_h, nmlon, nhgt_fix, hgt_fix_r
  use indices_mod, only : indices_init
  use mpi_module, only: mlat0, mlat1, mlon0, mlon1
  use mpi_module, only: mlond0, mlond1, mlatd0, mlatd1

  use outputdata_mod, only: outputdata_init, outputdata_close
  use outputdata_mod, only: outputdata_write_s1_fld
  use outputdata_mod, only: outputdata_write_s2_fld
  use outputdata_mod, only: outputdata_write_2dp_fld

  implicit none

  include 'mpif.h'

  !integer, parameter :: edyn3d_nmlon = 180
  !integer, parameter :: edyn3d_nhgt = 54
  !integer, parameter :: edyn3d_nmlat_h = 91

  integer, parameter :: edyn3d_nmlat_h = 45
  integer, parameter :: edyn3d_nmlon = 60
  integer, parameter :: edyn3d_nhgt = 26

! ionos_edyn3d_nmlat_h = 45
! ionos_edyn3d_nmlon = 60
! ionos_edyn3d_nhgt = 26

!  integer, parameter :: nglats = 90
!  integer, parameter :: nglons = 180

  real(r8), parameter :: geomag_year = 2000.50

  character(len=*), parameter :: prefix = '3D-Dynamamo test: '

  integer :: ierr, job_size, my_rank, itime

  real(r8), allocatable :: un_s1(:,:,:,:)
  real(r8), allocatable :: un_s2(:,:,:,:)
  real(r8), allocatable :: vn_s1(:,:,:,:)
  real(r8), allocatable :: vn_s2(:,:,:,:)

  real(r8), allocatable :: sigped_p(:,:,:,:)
  real(r8), allocatable :: sigma_hal_s1(:,:,:,:)
  real(r8), allocatable :: sigma_hal_s2(:,:,:,:)
  real(r8), allocatable :: sigma_ped_s1(:,:,:,:)
  real(r8), allocatable :: sigma_ped_s2(:,:,:,:)

  real(r8), allocatable :: efld1_s1(:,:,:)
  real(r8), allocatable :: efld2_s1(:,:,:)
  real(r8), allocatable :: efld1_s2(:,:,:)
  real(r8), allocatable :: efld2_s2(:,:,:)
  real(r8), allocatable :: ionvel1_s1(:,:,:)
  real(r8), allocatable :: ionvel2_s1(:,:,:)
  real(r8), allocatable :: ionvel1_s2(:,:,:)
  real(r8), allocatable :: ionvel2_s2(:,:,:)

  real(r8), allocatable :: ui_s1(:,:,:,:)
  real(r8), allocatable :: vi_s1(:,:,:,:)
  real(r8), allocatable :: wi_s1(:,:,:,:)

  real(r8), allocatable :: ui_s2(:,:,:,:)
  real(r8), allocatable :: vi_s2(:,:,:,:)
  real(r8), allocatable :: wi_s2(:,:,:,:)

  real(r8), allocatable :: pot_hl_p(:,:,:)
  real(r8), allocatable :: fac_hl_p(:,:,:)
  real(r8), allocatable :: elec_pot_p(:,:,:)
  real(r8), allocatable :: ped_cond_p(:,:,:)

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
  !call inputdata_init('../data/FX2000_ne16pg3_wcmx_3Dedyn_test02.cam.h5i.0001-01-01-00000.nc')
  !call outputdata_init('../data/FX2000_ne16pg3_wcmx_3Dedyn_test02.cam.h5i.0001-01-01-00000.nc','testout.nc')

  call inputdata_init('../data/FX2000_ne16pg3_wcmx_3Dedyn_test04.cam.h2i.0001-01-01-00000.nc')
  call outputdata_init('../data/FX2000_ne16pg3_wcmx_3Dedyn_test04.cam.h2i.0001-01-01-00000.nc','testout.nc')

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

  allocate( sigped_p(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  sigped_p = 0._r8

  allocate( ui_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  allocate( vi_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  allocate( wi_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  allocate( ui_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  allocate( vi_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )
  allocate( wi_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1) )

  allocate( efld1_s1(2,mlat0:mlat1,mlon0:mlon1) )
  efld1_s1 = 0._r8
  allocate( efld2_s1(2,mlat0:mlat1,mlon0:mlon1) )
  efld2_s1 = 0._r8
  allocate( efld1_s2(2,mlat0:mlat1,mlon0:mlon1) )
  allocate( efld2_s2(2,mlat0:mlat1,mlon0:mlon1) )
  allocate( ionvel1_s1(2,mlat0:mlat1,mlon0:mlon1) )
  allocate( ionvel2_s1(2,mlat0:mlat1,mlon0:mlon1) )
  allocate( ionvel1_s2(2,mlat0:mlat1,mlon0:mlon1) )
  allocate( ionvel2_s2(2,mlat0:mlat1,mlon0:mlon1) )

  allocate( pot_hl_p(2,mlatd0:mlatd1,mlond0:mlond1) )
  pot_hl_p = 0._r8
  allocate( fac_hl_p(2,mlatd0:mlatd1,mlond0:mlond1) )
  fac_hl_p = 0._r8
  allocate( elec_pot_p(2,mlat0:mlat1,mlon0:mlon1) )
  elec_pot_p = 0._r8
  allocate( ped_cond_p(2,mlat0:mlat1,mlon0:mlon1) )
  ped_cond_p = 0._r8

  do itime = 1,ntimes
     call inputdata_read_s1_fld( varid=un_s1_vid, time_ndx=itime, fld=un_s1 )
     call inputdata_read_s1_fld( varid=vn_s1_vid, time_ndx=itime, fld=vn_s1 )
     call inputdata_read_s1_fld( varid=sigma_hal_s1_vid, time_ndx=itime, fld=sigma_hal_s1 )
     call inputdata_read_s1_fld( varid=sigma_ped_s1_vid, time_ndx=itime, fld=sigma_ped_s1 )
     call inputdata_read_s2_fld( varid=un_s2_vid, time_ndx=itime, fld=un_s2 )
     call inputdata_read_s2_fld( varid=vn_s2_vid, time_ndx=itime, fld=vn_s2 )
     call inputdata_read_s2_fld( varid=sigma_hal_s2_vid, time_ndx=itime, fld=sigma_hal_s2 )
     call inputdata_read_s2_fld( varid=sigma_ped_s2_vid, time_ndx=itime, fld=sigma_ped_s2 )
     call inputdata_read_2dp_fld(hilat_pot_vid, itime, pot_hl_p)

     write(*,*) prefix,' ... itime:',itime,' maxval( sigma_ped_s2 ): ',maxval( sigma_ped_s2 )

     call outputdata_write_s1_fld( varname='un_s1', time_ndx=itime, fld=un_s1 )
     call outputdata_write_s2_fld( varname='un_s2', time_ndx=itime, fld=un_s2 )
     call outputdata_write_s1_fld( varname='vn_s1', time_ndx=itime, fld=vn_s1 )
     call outputdata_write_s2_fld( varname='vn_s2', time_ndx=itime, fld=vn_s2 )
     call outputdata_write_s1_fld( varname='sigma_hal_s1', time_ndx=itime, fld=sigma_hal_s1 )
     call outputdata_write_s1_fld( varname='sigma_ped_s1', time_ndx=itime, fld=sigma_ped_s1 )
     call outputdata_write_s2_fld( varname='sigma_hal_s2', time_ndx=itime, fld=sigma_hal_s2 )
     call outputdata_write_s2_fld( varname='sigma_ped_s2', time_ndx=itime, fld=sigma_ped_s2 )

     call dynamo_calc( &
          sigma_ped_s1, sigma_hal_s1, un_s1, vn_s1, &
          sigma_ped_s2, sigma_hal_s2, un_s2, vn_s2, &
          sigped_p, pot_hl_p, fac_hl_p, &
          ui_s1, vi_s1, wi_s1, ui_s2, vi_s2, wi_s2, &
          elec_pot_p, ped_cond_p, &
          efld1_s1, efld2_s1, efld1_s2, efld2_s2, &
          ionvel1_s1, ionvel2_s1, ionvel1_s2, ionvel2_s2)

     print*,'*** max min elec_pot_p: ',minval(elec_pot_p), maxval(elec_pot_p)

     call outputdata_write_2dp_fld('HILAT_POT',itime, pot_hl_p(:,mlat0:mlat1,mlon0:mlon1))
     call outputdata_write_2dp_fld('HILAT_FAC',itime, fac_hl_p(:,mlat0:mlat1,mlon0:mlon1))

     call outputdata_write_2dp_fld('ELECPOTEN',itime, elec_pot_p)

     call outputdata_write_s1_fld('IonU_s1', itime, ui_s1)
     call outputdata_write_s2_fld('IonU_s2', itime, ui_s2)
     call outputdata_write_s1_fld('IonV_s1', itime, vi_s1)
     call outputdata_write_s2_fld('IonV_s2', itime, vi_s2)
     call outputdata_write_s1_fld('IonW_s1', itime, wi_s1)
     call outputdata_write_s2_fld('IonW_s2', itime, wi_s2)

   !  call outputdata_write_2dp_fld('ED1s1',itime, efld1_s1)

     !efld1_s1 = 1._r8
    ! call outputdata_write_s1_fld('ED1s1', itime, efld1_s1)
    ! call outputdata_write_s1_fld('ED2s1', itime, efld2_s1)
    ! call outputdata_write_s2_fld('ED1s2', itime, efld1_s2)
    ! call outputdata_write_s2_fld('ED2s2', itime, efld2_s2)

    ! call outputdata_write_s1_fld('Ve1s1', itime, ionvel1_s1)
    ! call outputdata_write_s1_fld('Ve2s1', itime, ionvel2_s1)
    ! call outputdata_write_s2_fld('Ve1s2', itime, ionvel1_s2)
    ! call outputdata_write_s2_fld('Ve2s2', itime, ionvel2_s2)

  end do

  call dynamo_final()

  call inputdata_close()

  call outputdata_close()
  write(*,*) prefix,'DONE...'

  call mpi_finalize(ierr)
  print*,'... mpi_finalize ierr: ',ierr

  print*,'END TEST.'

end program main
