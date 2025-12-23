module outputdata_mod
  use iso_fortran_env, only: r8=>real64 ! double precision
  use mpi_module, only: mpi_rank
  use mpi_module, only: mlat0, mlat1, mlon0, mlon1
  use params_module, only: nhgt_fix, nmlatS2_h, nmlat_T1, nmlat_T2
  use inputdata_mod, only: nlons1, npflpts1, nlons2, npflpts2
  use inputdata_mod, only: nmaglat, nmaglon
  use inputdata_mod, only: nmaglat_s, nmaglon_s

  use netcdf

  implicit none

  integer :: outfid = -huge(1)

contains

!~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
  subroutine outputdata_close()
    if (mpi_rank==0) then
       call handle_error( nf90_close(outfid), 'ERROR: outputdata nf90_close' )
    end if
  end subroutine outputdata_close

!~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
  subroutine outputdata_init(infile, outfile)

    character(len=*), intent(in) :: infile
    character(len=*), intent(in) :: outfile

    integer :: infid, xtype, vid_in, ival
    integer :: dimids_in(1), len
    integer, allocatable :: xint(:)
    real(r8), allocatable :: xdbl(:)
    integer :: ndims, nvars, vid, did,  dimsize, dimid,  varsize, varid, n, natts
    integer, allocatable :: dimids(:)

    character(len=32) :: dimname
    character(len=32) :: varname
    character(len=32) :: attname

    character(len=*), parameter :: prefix = 'outputdata_init: '

    if (mpi_rank==0) then

       call handle_error(nf90_open(path=infile, mode=NF90_NOWRITE, ncid=infid), prefix//'nf90_open ERROR: '//infile)

       call handle_error(nf90_inquire(infid, ndimensions=ndims, nvariables=nvars), prefix//'nf90_inquire ndims ERROR')

       call handle_error(nf90_create(outfile, NF90_CLOBBER, outfid),  prefix//'nf90_create file: '//outfile)

       do did = 1,ndims
          call handle_error( nf90_inquire_dimension(infid, did, name=dimname, len=dimsize), 'ERROR: nf90_inquire_dimension ')
          if (dimname=='time') dimsize = NF90_UNLIMITED
          call handle_error(nf90_def_dim(outfid, dimname, dimsize, dimid), prefix//'nf90_def_dim '//dimname)
       end do

       do vid = 1,nvars
          call handle_error(nf90_inquire_variable(infid, vid, name=varname, ndims=ndims, natts=natts, xtype=xtype), &
               'ERROR: nf90_inquire_variable name, ndims, natts')

          allocate(dimids(ndims))
          call handle_error(nf90_inquire_variable(infid, vid, dimids=dimids), 'ERROR: nf90_inquire_variable dimids ')

          call handle_error(nf90_def_var(outfid, varname, xtype, dimids, varid), prefix//'nf90_def_dim '//dimname)

          do n = 1, natts
             call handle_error(nf90_inq_attname( infid, vid, n, attname ), prefix//'nf90_inq_attname' )
             call handle_error(nf90_copy_att( infid, vid, trim(attname), outfid, varid), prefix//'nf90_copy_att '//trim(attname) )
          enddo
          deallocate(dimids)
       end do

       call handle_error(nf90_enddef(outfid), prefix//' ERROR: nf90_enddef' )

       ! copy coords, etc...
       do vid = 1,nvars
          call handle_error(nf90_inquire_variable(outfid, vid, name=varname, xtype=xtype, ndims=ndims ), &
               prefix//' ERROR: nf90_inquire_variable')

          if (ndims==0) then
             call handle_error( nf90_inq_varid(infid, varname, vid_in), prefix//'nf90_inq_varid' )
             if (xtype==nf90_int) then
                call handle_error( nf90_get_var( infid, vid_in, ival ), prefix//'nf90_get_var' )
                call handle_error( nf90_put_var( outfid, vid, ival ), prefix//'nf90_put_var' )
             end if
          endif


          if (ndims==1) then
             call handle_error( nf90_inq_varid(infid, varname, vid_in), prefix//'nf90_inq_varid' )
             call handle_error( nf90_inquire_variable(infid, vid_in, dimids=dimids_in), prefix//' nf90_inquire_variable' )
             call handle_error( nf90_inquire_dimension( infid, dimids_in(1), len=len ), prefix//' nf90_inquire_dimension' )
             if (xtype==nf90_int) then
                allocate(xint(len))
                call handle_error( nf90_get_var( infid, vid_in, xint ), prefix//'nf90_get_var' )
                call handle_error( nf90_put_var( outfid, vid, xint ), prefix//'nf90_put_var' )
                deallocate(xint)
             else if (xtype==nf90_double .or. xtype==nf90_real) then
                allocate(xdbl(len))
                call handle_error( nf90_get_var( infid, vid_in, xdbl ), prefix//'nf90_get_var' )
                call handle_error( nf90_put_var( outfid, vid, xdbl ), prefix//'nf90_put_var' )
                deallocate(xdbl)
             end if
          end if
       end do

       call handle_error(nf90_close(infid), 'ERROR: nf90_close' )

    endif

  end subroutine outputdata_init

!~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
  subroutine outputdata_write_s1_fld( varname, time_ndx, fld )
    use indices_mod, only: flpts1ndx
    use dynamo_interface_mod, only: npts_s1

    use mpi, only: MPI_INTEGER, MPI_REAL8, MPI_SUCCESS, MPI_SUM, MPI_COMM_WORLD

    character(len=*), intent(in) :: varname
    integer, intent(in) :: time_ndx

    real(r8), intent(in) :: fld(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: sndbuf(nlons1,npflpts1)
    real(r8) :: rcvbuf(nlons1,npflpts1)

    integer :: cnt(3) ! nlons1,npflpts1,ntimes)
    integer :: strt(3)
    integer :: h,i,j,k, len, ier, varid

    len = nlons1*npflpts1

    sndbuf = 0._r8
    rcvbuf = 0._r8

    cnt(1) = nlons1
    cnt(2) = npflpts1
    cnt(3) = 1

    strt(1) = 1
    strt(2) = 1
    strt(3) = time_ndx

    do h = 1,2
       do i = mlon0,mlon1
          do j = mlat0,mlat1
             do k = 1,npts_s1(j)
                sndbuf(i,flpts1ndx(k,h,j)) = fld(k,h,j,i)
             end do
          end do
       end do
    end do

    ! Gather to root by using scalable reduce method:
    call mpi_reduce(sndbuf, rcvbuf, len, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ier )

    if (mpi_rank==0) then
       call handle_error( nf90_inq_varid(outfid, varname, varid), 'ERROR: nf90_inq_varid '//varname )
       call handle_error( nf90_put_var(outfid, varid, rcvbuf, start=strt, count=cnt), ' nf90_put_var error ')
    end if

  end subroutine outputdata_write_s1_fld


!~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
  subroutine outputdata_write_s2_fld( varname, time_ndx, fld )
    use indices_mod, only: flpts2ndx
    use dynamo_interface_mod, only: npts_s2

    use mpi, only: MPI_INTEGER, MPI_REAL8, MPI_SUCCESS, MPI_SUM, MPI_COMM_WORLD

    character(len=*), intent(in) :: varname
    integer, intent(in) :: time_ndx

    real(r8), intent(in) :: fld(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: sndbuf(nlons2,npflpts2)
    real(r8) :: rcvbuf(nlons2,npflpts2)

    integer :: cnt(3) ! nlons2,npflpts2,ntimes)
    integer :: strt(3)
    integer :: h,i,j,k, len, ier, varid

    len = nlons2*npflpts2

    sndbuf = 0._r8
    rcvbuf = 0._r8

    cnt(1) = nlons2
    cnt(2) = npflpts2
    cnt(3) = 1

    strt(1) = 1
    strt(2) = 1
    strt(3) = time_ndx

    do h = 1,2
       do i = mlon0,mlon1
          do j = mlat0,min(mlat1,nmlatS2_h)
             do k = 1,npts_s2(j)
                sndbuf(i,flpts2ndx(k,h,j)) = fld(k,h,j,i)
             end do
          end do
       end do
    end do

    ! Gather to root by using scalable reduce method:
    call mpi_reduce(sndbuf, rcvbuf, len, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ier )

    if (mpi_rank==0) then
       call handle_error( nf90_inq_varid(outfid, varname, varid), 'ERROR: nf90_inq_varid '//varname )
       call handle_error( nf90_put_var(outfid, varid, rcvbuf, start=strt, count=cnt), ' nf90_put_var error ')
    end if

  end subroutine outputdata_write_s2_fld

!~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
  subroutine outputdata_write_2dp_fld( varname, time_ndx, fld )
    use mpi, only: MPI_INTEGER, MPI_REAL8, MPI_SUCCESS, MPI_SUM, MPI_COMM_WORLD

    character(len=*), intent(in) :: varname
    integer, intent(in) :: time_ndx
    real(r8), intent(in) :: fld(2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: sndbuf(nmaglon,nmaglat)
    real(r8) :: rcvbuf(nmaglon,nmaglat)

    integer :: cnt(3)
    integer :: strt(3)
    integer :: h,i,j, jj
    integer :: len, ier, varid

    len = nmaglon*nmaglat

    sndbuf = 0._r8
    rcvbuf = 0._r8

    cnt(1) = nmaglon
    cnt(2) = nmaglat
    cnt(3) = 1

    strt(1) = 1
    strt(2) = 1
    strt(3) = time_ndx


    do h = 1,2
       do i = mlon0,mlon1
          do j = mlat0,mlat1
             if (h==2) then
                jj = nmlat_T1 - j + 1
             else
                jj = j
             end if
             sndbuf(i,jj) = fld(h,j,i)
          end do
       end do
    end do

    ! Gather to root by using scalable reduce method:
    call mpi_reduce(sndbuf, rcvbuf, len, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ier )

    if (mpi_rank==0) then
       call handle_error( nf90_inq_varid(outfid, varname, varid), 'ERROR: nf90_inq_varid '//varname )
       call handle_error( nf90_put_var(outfid, varid, rcvbuf, start=strt, count=cnt), ' nf90_put_var error ')
    end if

  end subroutine outputdata_write_2dp_fld

!~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
  subroutine outputdata_write_2ds1_fld( varname, time_ndx, fld )
    use mpi, only: MPI_INTEGER, MPI_REAL8, MPI_SUCCESS, MPI_SUM, MPI_COMM_WORLD

    character(len=*), intent(in) :: varname
    integer, intent(in) :: time_ndx
    real(r8), intent(in) :: fld(2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: sndbuf(nmaglon_s,nmaglat)
    real(r8) :: rcvbuf(nmaglon_s,nmaglat)

    integer :: cnt(3)
    integer :: strt(3)
    integer :: h,i,j, jj
    integer :: len, ier, varid

    len = nmaglon_s*nmaglat

    sndbuf = 0._r8
    rcvbuf = 0._r8

    cnt(1) = nmaglon_s
    cnt(2) = nmaglat
    cnt(3) = 1

    strt(1) = 1
    strt(2) = 1
    strt(3) = time_ndx


    do h = 1,2
       do i = mlon0,mlon1
          do j = mlat0,mlat1
             if (h==2) then
                jj = nmlat_T1 - j + 1
             else
                jj = j
             end if
             sndbuf(i,jj) = fld(h,j,i)
          end do
       end do
    end do

    ! Gather to root by using scalable reduce method:
    call mpi_reduce(sndbuf, rcvbuf, len, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ier )

    if (mpi_rank==0) then
       call handle_error( nf90_inq_varid(outfid, varname, varid), 'ERROR: nf90_inq_varid '//varname )
       call handle_error( nf90_put_var(outfid, varid, rcvbuf, start=strt, count=cnt), ' nf90_put_var error ')
    end if

  end subroutine outputdata_write_2ds1_fld


!~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
  subroutine outputdata_write_2ds2_fld( varname, time_ndx, fld )
    use mpi, only: MPI_INTEGER, MPI_REAL8, MPI_SUCCESS, MPI_SUM, MPI_COMM_WORLD

    character(len=*), intent(in) :: varname
    integer, intent(in) :: time_ndx
    real(r8), intent(in) :: fld(2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: sndbuf(nmaglon,nmaglat_s)
    real(r8) :: rcvbuf(nmaglon,nmaglat_s)

    integer :: cnt(3)
    integer :: strt(3)
    integer :: h,i,j, jj
    integer :: len, ier, varid

    len = nmaglon*nmaglat_s

    sndbuf = 0._r8
    rcvbuf = 0._r8

    cnt(1) = nmaglon
    cnt(2) = nmaglat_s
    cnt(3) = 1

    strt(1) = 1
    strt(2) = 1
    strt(3) = time_ndx

    do h = 1,2
       do i = mlon0,mlon1
          do j = mlat0,min(mlat1,nmlatS2_h)
             if (h==2) then
                jj = nmlat_T2 - j + 1
             else
                jj = j
             end if
             sndbuf(i,jj) = fld(h,j,i)
          end do
       end do
    end do

    ! Gather to root by using scalable reduce method:
    call mpi_reduce(sndbuf, rcvbuf, len, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ier )

    if (mpi_rank==0) then
       call handle_error( nf90_inq_varid(outfid, varname, varid), 'ERROR: nf90_inq_varid '//varname )
       call handle_error( nf90_put_var(outfid, varid, rcvbuf, start=strt, count=cnt), ' nf90_put_var error ')
    end if

  end subroutine outputdata_write_2ds2_fld

!~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
  subroutine handle_error( ierror, msg )
    integer, intent(in) :: ierror
    character(len=*), intent(in) :: msg

    if (ierror/=NF90_NOERR) then
       write(*,*) msg
       write(*,*) 'ERROR -- nf90 status number: ',ierror
       stop 1
    end if
  end subroutine handle_error

end module outputdata_mod
