module inputdata_mod
  use iso_fortran_env, only: r8=>real64 ! double precision

  use mpi_module, only: mpi_rank

  use netcdf, only: NF90_NOWRITE, NF90_NOERR, nf90_open, nf90_inquire
  use netcdf, only: nf90_inquire_dimension, nf90_inq_dimid
  use netcdf, only: nf90_close, nf90_inq_varid, nf90_get_var
  use mpi_module, only: mlat0, mlat1, mlon0, mlon1
  use params_module, only: nhgt_fix, nmlatS2_h

  implicit none

  integer, protected :: ncid = -1
  integer, protected :: ntimes = -1
  integer, protected :: npflpts1 = -1
  integer, protected :: npflpts2 = -1
  integer, protected :: nlons1 = -1
  integer, protected :: nlons2 = -1

  integer, protected :: un_s1_vid = -1
  integer, protected :: vn_s1_vid = -1
  integer, protected :: un_s2_vid = -1
  integer, protected :: vn_s2_vid = -1
  integer, protected :: sigma_hal_s1_vid = -1
  integer, protected :: sigma_hal_s2_vid = -1
  integer, protected :: sigma_ped_s1_vid = -1
  integer, protected :: sigma_ped_s2_vid = -1

contains

  subroutine inputdata_close()
    call handle_error( nf90_close(ncid), 'ERROR: nf90_close' )
  end subroutine inputdata_close

  subroutine inputdata_init(filepath)

    character(len=*), intent(in) :: filepath

    integer :: status, ndims, did

    character(len=*), parameter :: prefix = 'inputdata_read: '
    character(len=32) :: dim_name

    write(*,*) prefix//' open: ',trim(filepath)

    call handle_error( nf90_open(path=filepath, mode=NF90_NOWRITE, ncid=ncid), prefix//'nf90_open ERROR: '//filepath )

    call handle_error(nf90_inquire(ncid, ndimensions=ndims), prefix//'nf90_inquire ndims ERROR')

    npflpts1 = get_dimsize( 'pflpts1' )
    npflpts2 = get_dimsize( 'pflpts2' )
    nlons1 = get_dimsize( 'lon_s1' )
    nlons2 = get_dimsize( 'lon_s2' )
    ntimes = get_dimsize( 'time' )

    call handle_error( nf90_inq_varid(ncid, 'un_s1', un_s1_vid), 'ERROR: nf90_inq_varid un_s1' )
    call handle_error( nf90_inq_varid(ncid, 'un_s2', un_s2_vid), 'ERROR: nf90_inq_varid un_s2' )
    call handle_error( nf90_inq_varid(ncid, 'vn_s1', vn_s1_vid), 'ERROR: nf90_inq_varid vn_s1' )
    call handle_error( nf90_inq_varid(ncid, 'vn_s2', vn_s2_vid), 'ERROR: nf90_inq_varid vn_s2' )

    call handle_error( nf90_inq_varid(ncid, 'sigma_hal_s1', sigma_hal_s1_vid), 'ERROR: nf90_inq_varid sigma_hal_s1' )
    call handle_error( nf90_inq_varid(ncid, 'sigma_hal_s2', sigma_hal_s2_vid), 'ERROR: nf90_inq_varid sigma_hal_s2' )
    call handle_error( nf90_inq_varid(ncid, 'sigma_ped_s1', sigma_ped_s1_vid), 'ERROR: nf90_inq_varid sigma_ped_s1' )
    call handle_error( nf90_inq_varid(ncid, 'sigma_ped_s2', sigma_ped_s2_vid), 'ERROR: nf90_inq_varid sigma_ped_s2' )

  contains

    function get_dimsize( dimname ) result(size)
      character(len=*), intent(in) :: dimname
      integer :: size
      integer :: did

      call handle_error( nf90_inq_dimid(ncid, dimname, did), 'ERROR: nf90_inq_dimid '//dimname)
      call handle_error( nf90_inquire_dimension(ncid, did, len=size), 'ERROR: nf90_inquire_dimension '//dimname)
    end function get_dimsize

  end subroutine inputdata_init


  subroutine handle_error( ierror, msg )
    integer, intent(in) :: ierror
    character(len=*), intent(in) :: msg

    if (ierror/=NF90_NOERR) then
       write(*,*) msg
       write(*,*) 'ERROR -- nf90 status number: ',ierror
       stop 1
    end if
  end subroutine handle_error

  subroutine inputdata_read_s1_fld( varid, time_ndx, fld )
    use indices_mod, only: flpts1ndx
    use dynamo_interface_mod, only: npts_s1

    integer, intent(in) :: varid, time_ndx

    real(r8), intent(out) :: fld(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: flddata(nlons1,npflpts1)
    integer :: cnt(3) ! nlons1,npflpts1,ntimes)
    integer :: strt(3)
    integer :: h,i,j,k

    fld = -huge(1._r8)

    cnt(1) = nlons1
    cnt(2) = npflpts1
    cnt(3) = 1

    strt(1) = 1
    strt(2) = 1
    strt(3) = time_ndx

    call handle_error( nf90_get_var(ncid, varid, flddata, start=strt, count=cnt), ' nf90_get_var error ')

    do h = 1,2
       do i = mlon0,mlon1
          do j = mlat0,mlat1
             do k = 1,npts_s1(j)
                fld(k,h,j,i) = flddata(i,flpts1ndx(k,h,j))
             end do
          end do
       end do
    end do
  end subroutine inputdata_read_s1_fld


  subroutine inputdata_read_s2_fld( varid, time_ndx, fld )
    use indices_mod, only: flpts2ndx
    use dynamo_interface_mod, only: npts_s2

    integer, intent(in) :: varid, time_ndx

    real(r8), intent(out) :: fld(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: flddata(nlons2,npflpts2)
    integer :: cnt(3) ! nlons2,npflpts2,ntimes)
    integer :: strt(3)
    integer :: h,i,j,k

    fld = -huge(1._r8)

    cnt(1) = nlons2
    cnt(2) = npflpts2
    cnt(3) = 1

    strt(1) = 1
    strt(2) = 1
    strt(3) = time_ndx

    call handle_error( nf90_get_var(ncid, varid, flddata, start=strt, count=cnt), ' nf90_get_var error ')

    do h = 1,2
       do i = mlon0,mlon1
          do j = mlat0,min(mlat1,nmlatS2_h)
             do k = 1,npts_s2(j)
                fld(k,h,j,i) = flddata(i,flpts2ndx(k,h,j))
             end do
          end do
       end do
    end do
  end subroutine inputdata_read_s2_fld


end module inputdata_mod
