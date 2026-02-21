module mpi_mod

  use iso_fortran_env, only: real32,real64
  use prec, only: rp
  use mpi

  implicit none

  integer :: mpi_rp = -huge(1)

! magnetic decomposition
  integer :: mag_world = -huge(1)
  integer :: mag_size=0, mag_rank=-1, &
    mag_lat_size=0, mag_lat_rank=-huge(1), mag_nlat=0, mag_maxlat=0, &
    mag_lat0=1, mag_lat1=0, mag_latd0=1, mag_latd1=0, &
    mag_lon_size=0, mag_lon_rank=-huge(1), mag_nlon=0, mag_maxlon=0, &
    mag_lon0=1, mag_lon1=0, mag_lond0=1, mag_lond1=0

  integer, dimension(:), allocatable :: &
    mag_nlat_task, mag_lat0_task, mag_lat1_task, &
    mag_nlon_task, mag_lon0_task, mag_lon1_task

  interface gather_mag_lon ! gather magnetic fields in longitudes
    module procedure gather_mag_lon_0d, gather_mag_lon_1d, gather_mag_lon_3d
  endinterface gather_mag_lon

  interface gather_mag ! gather magnetic fields
    module procedure gather_mag_2d, gather_mag_3d, gather_mag_4d, gather_mag_5d
  endinterface gather_mag

  contains

!-----------------------------------------------------------------------
  subroutine setup_comm(mpi_comm_host, mag_size_in)

    integer, intent(in) :: mpi_comm_host
    integer, intent(in) :: mag_size_in

    integer :: color, ierror, npes_host

    if (rp == real32) then
      mpi_rp = MPI_REAL4
    elseif (rp == real64) then
      mpi_rp = MPI_REAL8
    else
      stop 'unknown real precision'
    endif

! Get the size and rank in the new communicator
    call mpi_comm_size(mpi_comm_host, npes_host, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_size', ierror)

    mag_size = min(npes_host, mag_size_in)

    call MPI_Comm_rank(mpi_comm_host, mag_rank, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_rank', ierror)

    color = mag_rank/mag_size
    call mpi_comm_split(mpi_comm_host, color, mag_rank, mag_world, ierror)

! factorize group size to the nearest two numbers
    do mag_lat_size = int(sqrt(real(mag_size, kind=rp))), 1, -1
      mag_lon_size = mag_size / mag_lat_size
      if (mag_lon_size*mag_lat_size == mag_size) exit
    enddo

! 2D index increases faster in latitudes
! (stack along latitudes first then longitudes)
! (lat_size=3)
!  8  9 10 11
!  4  5  6  7
!  0  1  2  3 (lon_size=4)
    mag_lat_rank = mag_rank / mag_lon_size
    mag_lon_rank = modulo(mag_rank, mag_lon_size)

  endsubroutine setup_comm

!-----------------------------------------------------------------------
  subroutine setup_mag_topology(nlat, nlon)
! setup magnetic decomposition and the connectivity matrix
! similar to geographic decomposition

    integer, intent(in) :: nlat, nlon

    integer :: i, j, rnk, rnki, rnkj

    mag_nlat = nlat
    mag_nlon = nlon

    allocate(mag_nlat_task(0:mag_lat_size-1))
    allocate(mag_nlon_task(0:mag_lon_size-1))
    allocate(mag_lat0_task(0:mag_size-1))
    allocate(mag_lat1_task(0:mag_size-1))
    allocate(mag_lon0_task(0:mag_size-1))
    allocate(mag_lon1_task(0:mag_size-1))

    mag_nlat_task = generate_minvar_list(mag_nlat, mag_lat_size)
    mag_maxlat = maxval(mag_nlat_task)
    mag_lat0 = 1
    do j = 0, mag_lat_rank-1
      mag_lat0 = mag_lat0 + mag_nlat_task(j)
    enddo
    mag_lat1 = mag_lat0 + mag_nlat_task(mag_lat_rank) - 1

    mag_nlon_task = generate_minvar_list(mag_nlon, mag_lon_size)
    mag_maxlon = maxval(mag_nlon_task)
    mag_lon0 = 1
    do i = 0, mag_lon_rank-1
      mag_lon0 = mag_lon0 + mag_nlon_task(i)
    enddo
    mag_lon1 = mag_lon0 + mag_nlon_task(mag_lon_rank) - 1

    mag_latd0 = mag_lat0 - 1
    mag_latd1 = mag_lat1 + 1
    mag_lond0 = mag_lon0 - 1
    mag_lond1 = mag_lon1 + 1

    do rnk = 0, mag_size-1
      rnkj = rnk / mag_lon_size
      rnki = modulo(rnk, mag_lon_size)

      mag_lat0_task(rnk) = 1
      do j = 0, rnkj-1
        mag_lat0_task(rnk) = mag_lat0_task(rnk) + mag_nlat_task(j)
      enddo
      mag_lat1_task(rnk) = mag_lat0_task(rnk) + mag_nlat_task(rnkj) - 1

      mag_lon0_task(rnk) = 1
      do i = 0, rnki-1
        mag_lon0_task(rnk) = mag_lon0_task(rnk) + mag_nlon_task(i)
      enddo
      mag_lon1_task(rnk) = mag_lon0_task(rnk) + mag_nlon_task(rnki) - 1
    enddo

  endsubroutine setup_mag_topology

!-----------------------------------------------------------------------
  subroutine sync_mag_lat_5d(var, l, m, n)
! longitude halo points are not included

    integer, intent(in) :: l, m, n
    real(kind=rp), dimension(l, m, n, mag_latd0:mag_latd1, mag_lon0:mag_lon1), intent(inout) :: var

    integer :: below, above, cnt, i, lc, mc, nc, ierror
    real(kind=rp), dimension(l, m, n, mag_maxlon) :: &
      send_to_below, send_to_above, recv_from_below, recv_from_above
    integer, dimension(4) :: request

! find the rank of adjacent processes
    if (mag_lat_rank == 0) then
      below = MPI_PROC_NULL
    else
      below = mag_rank - mag_lon_size
    endif

    if (mag_lat_rank == mag_lat_size-1) then
      above = MPI_PROC_NULL
    else
      above = mag_rank + mag_lon_size
    endif

    cnt = l * m * n * mag_maxlon

! load to work array
    do concurrent (i = 1:mag_lon1-mag_lon0+1, nc = 1:n, mc = 1:m, lc = 1:l)
      send_to_below(lc, mc, nc, i) = var(lc, mc, nc, mag_lat0, i+mag_lon0-1)
      send_to_above(lc, mc, nc, i) = var(lc, mc, nc, mag_lat1, i+mag_lon0-1)
   enddo

! sync in latitude
    call MPI_Isend(send_to_below, cnt, mpi_rp, &
      below, 0, mag_world, request(1), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

    call MPI_Isend(send_to_above, cnt, mpi_rp, &
      above, 1, mag_world, request(2), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

    call MPI_Irecv(recv_from_above, cnt, mpi_rp, &
      above, 0, mag_world, request(3), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)

    call MPI_Irecv(recv_from_below, cnt, mpi_rp, &
      below, 1, mag_world, request(4), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)

! wait for sync to complete
    call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

! unpack to model fields
    if (mag_lat_rank /= 0) then
      do concurrent (i = mag_lon0:mag_lon1, nc = 1:n, mc = 1:m, lc = 1:l)
        var(lc, mc, nc, mag_latd0, i) = recv_from_below(lc, mc, nc, i-mag_lon0+1)
      enddo
    endif

    if (mag_lat_rank /= mag_lat_size-1) then
      do concurrent (i = mag_lon0:mag_lon1, nc = 1:n, mc = 1:m, lc = 1:l)
        var(lc, mc, nc, mag_latd1, i) = recv_from_above(lc, mc, nc, i-mag_lon0+1)
      enddo
    endif

  endsubroutine sync_mag_lat_5d
!-----------------------------------------------------------------------
  subroutine sync_mag_lon_5d(var, l, m, n)

    integer, intent(in) :: l, m, n
    real(kind=rp), dimension(l, m, n, mag_latd0:mag_latd1, mag_lond0:mag_lond1), intent(inout) :: var

    integer :: j, lc, mc, nc, left, right, cnt, ierror
    real(kind=rp), dimension(l, m, n, mag_maxlat+4) :: &
      send_to_left, send_to_right, recv_from_left, recv_from_right
    integer, dimension(4) :: request

! find the rank of adjacent processes
    if (mag_lon_rank == 0) then
      left = mag_rank - 1 + mag_lon_size
    else
      left = mag_rank - 1
    endif

    if (mag_lon_rank == mag_lon_size-1) then
      right = mag_rank + 1 - mag_lon_size
    else
      right = mag_rank + 1
    endif

    cnt = l * m * n * (mag_maxlat + 4)

! load to work array
    do concurrent (j = 1:mag_latd1-mag_latd0+1, nc = 1:n, mc = 1:m, lc = 1:l)
      send_to_left(lc, mc, nc, j) = var(lc, mc, nc, j+mag_latd0-1, mag_lon0)
      send_to_right(lc, mc, nc, j) = var(lc, mc, nc, j+mag_latd0-1, mag_lon1)
    enddo

! sync in longitude
    call MPI_Isend(send_to_left, cnt, mpi_rp, &
      left, 0, mag_world, request(1), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

    call MPI_Isend(send_to_right, cnt, mpi_rp, &
      right, 1, mag_world, request(2), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

    call MPI_Irecv(recv_from_right, cnt, mpi_rp, &
      right, 0, mag_world, request(3), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)

    call MPI_Irecv(recv_from_left, cnt, mpi_rp, &
      left, 1, mag_world, request(4), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)

! wait for sync to complete
    call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

! unpack to model fields
    do concurrent (j = mag_latd0:mag_latd1, nc = 1:n, mc = 1:m, lc = 1:l)
      var(lc, mc, nc, j, mag_lond0) = recv_from_left(lc, mc, nc, j-mag_latd0+1)
      var(lc, mc, nc, j, mag_lond1) = recv_from_right(lc, mc, nc, j-mag_latd0+1)
    enddo

  endsubroutine sync_mag_lon_5d
!-----------------------------------------------------------------------
  function gather_mag_lon_0d(varin, lat_rank, root) result(varout)

    real(kind=rp), intent(in) :: varin
    integer, intent(in), optional :: lat_rank, root
    real(kind=rp), dimension(0:mag_lon_size-1) :: varout

    integer :: root_lat_rank, rnk, rnki, ierror
    integer, dimension(0:mag_lon_size*2-1) :: request

! every process in the longitude ring sends, root process receives
    if (present(root)) then
      root_lat_rank = root / mag_lon_size
      if (present(lat_rank)) then
        if (lat_rank /= root_lat_rank) stop 'incompatible lat_rank and root'
      endif

      if (mag_lat_rank == root_lat_rank) then
        call MPI_Isend(varin, 1, mpi_rp, root, 0, &
          mag_world, request(mag_lon_size), ierror)
        if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

        if (mag_rank == root) then
          do rnki = 0, mag_lon_size-1
            rnk = root_lat_rank*mag_lon_size + rnki
            call MPI_Irecv(varout(rnki), 1, mpi_rp, &
              rnk, 0, mag_world, request(rnki), ierror)
            if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
          enddo

          call MPI_Waitall(mag_lon_size+1, request(0:mag_lon_size), MPI_STATUSES_IGNORE, ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)
        else
          call MPI_Wait(request(mag_lon_size), MPI_STATUS_IGNORE, ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Wait', ierror)
        endif
      endif
    else

! every process in the longitude ring sends and receives
      if (present(lat_rank)) then
        if (mag_lat_rank == lat_rank) then
          do rnki = 0, mag_lon_size-1
            rnk = lat_rank*mag_lon_size + rnki

            call MPI_Isend(varin, 1, mpi_rp, rnk, 0, &
              mag_world, request(mag_lon_size+rnki), ierror)
            if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

            call MPI_Irecv(varout(rnki), 1, mpi_rp, &
              rnk, 0, mag_world, request(rnki), ierror)
            if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
          enddo

          call MPI_Waitall(mag_lon_size*2, request, MPI_STATUSES_IGNORE, ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)
        endif

! every process sends and receives
      else
        root_lat_rank = mag_rank / mag_lon_size

        do rnki = 0, mag_lon_size-1
          rnk = root_lat_rank*mag_lon_size + rnki

          call MPI_Isend(varin, 1, mpi_rp, rnk, 0, &
            mag_world, request(mag_lon_size+rnki), ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

          call MPI_Irecv(varout(rnki), 1, mpi_rp, &
            rnk, 0, mag_world, request(rnki), ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
        enddo

        call MPI_Waitall(mag_lon_size*2, request, MPI_STATUSES_IGNORE, ierror)
        if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)
      endif
    endif

  endfunction gather_mag_lon_0d
!-----------------------------------------------------------------------
  function gather_mag_lon_1d(varin, lat_rank, root) result(varout)

    real(kind=rp), dimension(mag_lon0:mag_lon1), intent(in) :: varin
    integer, intent(in), optional :: lat_rank, root
    real(kind=rp), dimension(mag_nlon) :: varout

    integer :: cnt, i, root_lat_rank, rnk, rnki, i0, i1, ierror
    real(kind=rp), dimension(mag_maxlon) :: sendbuf
    real(kind=rp), dimension(mag_maxlon, 0:mag_lon_size-1) :: recvbuf
    integer, dimension(0:mag_lon_size*2-1) :: request

! load to work array
    cnt = mag_maxlon
    do concurrent (i = 1:mag_lon1-mag_lon0+1)
      sendbuf(i) = varin(i+mag_lon0-1)
    enddo

! every process in the longitude ring sends, root process receives
    if (present(root)) then
      root_lat_rank = root / mag_lon_size
      if (present(lat_rank)) then
        if (lat_rank /= root_lat_rank) stop 'incompatible lat_rank and root'
      endif

      if (mag_lat_rank == root_lat_rank) then
        call MPI_Isend(sendbuf, cnt, mpi_rp, root, 0, &
          mag_world, request(mag_lon_size), ierror)
        if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

        if (mag_rank == root) then
          do rnki = 0, mag_lon_size-1
            rnk = root_lat_rank*mag_lon_size + rnki
            call MPI_Irecv(recvbuf(:, rnki), cnt, mpi_rp, &
              rnk, 0, mag_world, request(rnki), ierror)
            if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
          enddo

          call MPI_Waitall(mag_lon_size+1, request(0:mag_lon_size), MPI_STATUSES_IGNORE, ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

! reconstruct longitude rings
          do rnki = 0, mag_lon_size-1
            rnk = root_lat_rank*mag_lon_size + rnki
            i0 = mag_lon0_task(rnk)
            i1 = mag_lon1_task(rnk)
            do concurrent (i = i0:i1)
              varout(i) = recvbuf(i-i0+1, rnki)
            enddo
          enddo
        else
          call MPI_Wait(request(mag_lon_size), MPI_STATUS_IGNORE, ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Wait', ierror)
        endif
      endif
    else

! every process in the longitude ring sends and receives
      if (present(lat_rank)) then
        if (mag_lat_rank == lat_rank) then
          do rnki = 0, mag_lon_size-1
            rnk = lat_rank*mag_lon_size + rnki

            call MPI_Isend(sendbuf, cnt, mpi_rp, rnk, 0, &
              mag_world, request(mag_lon_size+rnki), ierror)
            if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

            call MPI_Irecv(recvbuf(:, rnki), cnt, mpi_rp, &
              rnk, 0, mag_world, request(rnki), ierror)
            if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
          enddo

          call MPI_Waitall(mag_lon_size*2, request, MPI_STATUSES_IGNORE, ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

! reconstruct longitude rings
          do rnki = 0, mag_lon_size-1
            rnk = lat_rank*mag_lon_size + rnki
            i0 = mag_lon0_task(rnk)
            i1 = mag_lon1_task(rnk)
            do concurrent (i = i0:i1)
              varout(i) = recvbuf(i-i0+1, rnki)
            enddo
          enddo
        endif

! every process sends and receives
      else
        root_lat_rank = mag_rank / mag_lon_size

        do rnki = 0, mag_lon_size-1
          rnk = root_lat_rank*mag_lon_size + rnki

          call MPI_Isend(sendbuf, 1, mpi_rp, rnk, 0, &
            mag_world, request(mag_lon_size+rnki), ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

          call MPI_Irecv(recvbuf(:, rnki), 1, mpi_rp, &
            rnk, 0, mag_world, request(rnki), ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
        enddo

        call MPI_Waitall(mag_lon_size*2, request, MPI_STATUSES_IGNORE, ierror)
        if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

! reconstruct longitude rings
        do rnki = 0, mag_lon_size-1
          rnk = root_lat_rank*mag_lon_size + rnki
          i0 = mag_lon0_task(rnk)
          i1 = mag_lon1_task(rnk)
          do concurrent (i = i0:i1)
            varout(i) = recvbuf(i-i0+1, rnki)
          enddo
        enddo
      endif
    endif

  endfunction gather_mag_lon_1d
!-----------------------------------------------------------------------
  function gather_mag_lon_3d(varin, m, n, lat_rank, root) result(varout)

    integer, intent(in) :: m, n
    real(kind=rp), dimension(m, n, mag_lon0:mag_lon1), intent(in) :: varin
    integer, intent(in), optional :: lat_rank, root
    real(kind=rp), dimension(m, n, mag_nlon) :: varout

    integer :: cnt, i, mc, nc, root_lat_rank, rnk, rnki, i0, i1, ierror
    real(kind=rp), dimension(m, n, mag_maxlon) :: sendbuf
    real(kind=rp), dimension(m, n, mag_maxlon, 0:mag_lon_size-1) :: recvbuf
    integer, dimension(0:mag_lon_size*2-1) :: request

! load to work array
    cnt = m * n * mag_maxlon
    do concurrent (i = 1:mag_lon1-mag_lon0+1, nc = 1:n, mc = 1:m)
      sendbuf(mc, nc, i) = varin(mc, nc, i+mag_lon0-1)
    enddo

! every process in the longitude ring sends, root process receives
    if (present(root)) then
      root_lat_rank = root / mag_lon_size
      if (present(lat_rank)) then
        if (lat_rank /= root_lat_rank) stop 'incompatible lat_rank and root'
      endif

      if (mag_lat_rank == root_lat_rank) then
        call MPI_Isend(sendbuf, cnt, mpi_rp, root, 0, &
          mag_world, request(mag_lon_size), ierror)
        if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

        if (mag_rank == root) then
          do rnki = 0, mag_lon_size-1
            rnk = root_lat_rank*mag_lon_size + rnki
            call MPI_Irecv(recvbuf(:, :, :, rnki), cnt, mpi_rp, &
              rnk, 0, mag_world, request(rnki), ierror)
            if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
          enddo

          call MPI_Waitall(mag_lon_size+1, request(0:mag_lon_size), MPI_STATUSES_IGNORE, ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

! reconstruct longitude rings
          do rnki = 0, mag_lon_size-1
            rnk = root_lat_rank*mag_lon_size + rnki
            i0 = mag_lon0_task(rnk)
            i1 = mag_lon1_task(rnk)
            do concurrent (i = i0:i1, nc = 1:n, mc = 1:m)
              varout(mc, nc, i) = recvbuf(mc, nc, i-i0+1, rnki)
            enddo
          enddo
        else
          call MPI_Wait(request(mag_lon_size), MPI_STATUS_IGNORE, ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Wait', ierror)
        endif
      endif
    else

! every process in the longitude ring sends and receives
      if (present(lat_rank)) then
        if (mag_lat_rank == lat_rank) then
          do rnki = 0, mag_lon_size-1
            rnk = lat_rank*mag_lon_size + rnki

            call MPI_Isend(sendbuf, cnt, mpi_rp, rnk, 0, &
              mag_world, request(mag_lon_size+rnki), ierror)
            if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

            call MPI_Irecv(recvbuf(:, :, :, rnki), cnt, mpi_rp, &
              rnk, 0, mag_world, request(rnki), ierror)
            if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
          enddo

          call MPI_Waitall(mag_lon_size*2, request, MPI_STATUSES_IGNORE, ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

! reconstruct longitude rings
          do rnki = 0, mag_lon_size-1
            rnk = lat_rank*mag_lon_size + rnki
            i0 = mag_lon0_task(rnk)
            i1 = mag_lon1_task(rnk)
            do concurrent (i = i0:i1, nc = 1:n, mc = 1:m)
              varout(mc, nc, i) = recvbuf(mc, nc, i-i0+1, rnki)
            enddo
          enddo
        endif

! every process sends and receives
      else
        root_lat_rank = mag_rank / mag_lon_size

        do rnki = 0, mag_lon_size-1
          rnk = root_lat_rank*mag_lon_size + rnki

          call MPI_Isend(sendbuf, 1, mpi_rp, rnk, 0, &
            mag_world, request(mag_lon_size+rnki), ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

          call MPI_Irecv(recvbuf(:, :, :, rnki), 1, mpi_rp, &
            rnk, 0, mag_world, request(rnki), ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
        enddo

        call MPI_Waitall(mag_lon_size*2, request, MPI_STATUSES_IGNORE, ierror)
        if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

! reconstruct longitude rings
        do rnki = 0, mag_lon_size-1
          rnk = root_lat_rank*mag_lon_size + rnki
          i0 = mag_lon0_task(rnk)
          i1 = mag_lon1_task(rnk)
          do concurrent (i = i0:i1, nc = 1:n, mc = 1:m)
            varout(mc, nc, i) = recvbuf(mc, nc, i-i0+1, rnki)
          enddo
        enddo
      endif
    endif

  endfunction gather_mag_lon_3d
!-----------------------------------------------------------------------
  function gather_mag_2d(varin, root) result(varout)

    integer, intent(in) :: root
    real(kind=rp), dimension(mag_lat0:mag_lat1, mag_lon0:mag_lon1), intent(in) :: varin
    real(kind=rp), dimension(mag_nlat, mag_nlon) :: varout

    integer :: i, j, cnt, rnk, i0, i1, j0, j1, ierror
    real(kind=rp), dimension(mag_maxlat, mag_maxlon) :: sendbuf
    real(kind=rp), dimension(mag_maxlat, mag_maxlon, 0:mag_size-1) :: recvbuf

    cnt = mag_maxlat * mag_maxlon

! load to work array
    do concurrent (i = 1:mag_lon1-mag_lon0+1, j = 1:mag_lat1-mag_lat0+1)
      sendbuf(j, i) = varin(j+mag_lat0-1, i+mag_lon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, mag_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierror)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, mag_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)
    endif

! reconstruct based on lat-lon decomposition
    if (root<0 .or. mag_rank==root) then
      do rnk = 0, mag_size-1
        i0 = mag_lon0_task(rnk)
        i1 = mag_lon1_task(rnk)
        j0 = mag_lat0_task(rnk)
        j1 = mag_lat1_task(rnk)
        do concurrent (i = i0:i1, j = j0:j1)
          varout(j, i) = recvbuf(j-j0+1, i-i0+1, rnk)
        enddo
      enddo
    endif

  endfunction gather_mag_2d
!-----------------------------------------------------------------------
  function gather_mag_3d(varin, n, root) result(varout)

    integer, intent(in) :: n, root
    real(kind=rp), dimension(n, mag_lat0:mag_lat1, mag_lon0:mag_lon1), intent(in) :: varin
    real(kind=rp), dimension(n, mag_nlat, mag_nlon) :: varout

    integer :: i, j, nc, cnt, rnk, i0, i1, j0, j1, ierror
    real(kind=rp), dimension(n, mag_maxlat, mag_maxlon) :: sendbuf
    real(kind=rp), dimension(n, mag_maxlat, mag_maxlon, 0:mag_size-1) :: recvbuf

    cnt = n * mag_maxlat * mag_maxlon

! load to work array
    do concurrent (i = 1:mag_lon1-mag_lon0+1, j = 1:mag_lat1-mag_lat0+1, nc = 1:n)
      sendbuf(nc, j, i) = varin(nc, j+mag_lat0-1, i+mag_lon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, mag_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierror)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, mag_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)
    endif

! reconstruct based on lat-lon decomposition
    if (root<0 .or. mag_rank==root) then
      do rnk = 0, mag_size-1
        i0 = mag_lon0_task(rnk)
        i1 = mag_lon1_task(rnk)
        j0 = mag_lat0_task(rnk)
        j1 = mag_lat1_task(rnk)
        do concurrent (i = i0:i1, j = j0:j1, nc = 1:n)
          varout(nc, j, i) = recvbuf(nc, j-j0+1, i-i0+1, rnk)
        enddo
      enddo
    endif

  endfunction gather_mag_3d
!-----------------------------------------------------------------------
  function gather_mag_4d(varin, m, n, root) result(varout)

    integer, intent(in) :: m, n, root
    real(kind=rp), dimension(m, n, mag_lat0:mag_lat1, mag_lon0:mag_lon1), intent(in) :: varin
    real(kind=rp), dimension(m, n, mag_nlat, mag_nlon) :: varout

    integer :: i, j, mc, nc, cnt, rnk, i0, i1, j0, j1, ierror
    real(kind=rp), dimension(m, n, mag_maxlat, mag_maxlon) :: sendbuf
    real(kind=rp), dimension(m, n, mag_maxlat, mag_maxlon, 0:mag_size-1) :: recvbuf

    cnt = m * n * mag_maxlat * mag_maxlon

! load to work array
    do concurrent (i = 1:mag_lon1-mag_lon0+1, j = 1:mag_lat1-mag_lat0+1, nc = 1:n, mc = 1:m)
      sendbuf(mc, nc, j, i) = varin(mc, nc, j+mag_lat0-1, i+mag_lon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, mag_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierror)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, mag_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)
    endif

! reconstruct based on lat-lon decomposition
    if (root<0 .or. mag_rank==root) then
      do rnk = 0, mag_size-1
        i0 = mag_lon0_task(rnk)
        i1 = mag_lon1_task(rnk)
        j0 = mag_lat0_task(rnk)
        j1 = mag_lat1_task(rnk)
        do concurrent (i = i0:i1, j = j0:j1, nc = 1:n, mc = 1:m)
          varout(mc, nc, j, i) = recvbuf(mc, nc, j-j0+1, i-i0+1, rnk)
        enddo
      enddo
    endif

  endfunction gather_mag_4d
!-----------------------------------------------------------------------
  function gather_mag_5d(varin, l, m, n, root) result(varout)

    integer, intent(in) :: l, m, n, root
    real(kind=rp), dimension(l, m, n, mag_lat0:mag_lat1, mag_lon0:mag_lon1), intent(in) :: varin
    real(kind=rp), dimension(l, m, n, mag_nlat, mag_nlon) :: varout

    integer :: i, j, lc, mc, nc, cnt, rnk, i0, i1, j0, j1, ierror
    real(kind=rp), dimension(l, m, n, mag_maxlat, mag_maxlon) :: sendbuf
    real(kind=rp), dimension(l, m, n, mag_maxlat, mag_maxlon, 0:mag_size-1) :: recvbuf

    cnt = l * m * n * mag_maxlat * mag_maxlon

! load to work array
    do concurrent (i = 1:mag_lon1-mag_lon0+1, j = 1:mag_lat1-mag_lat0+1, nc = 1:n, mc = 1:m, lc = 1:l)
      sendbuf(lc, mc, nc, j, i) = varin(lc, mc, nc, j+mag_lat0-1, i+mag_lon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, mag_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierror)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, mag_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)
    endif

! reconstruct based on lat-lon decomposition
    if (root<0 .or. mag_rank==root) then
      do rnk = 0, mag_size-1
        i0 = mag_lon0_task(rnk)
        i1 = mag_lon1_task(rnk)
        j0 = mag_lat0_task(rnk)
        j1 = mag_lat1_task(rnk)
        do concurrent (i = i0:i1, j = j0:j1, nc = 1:n, mc = 1:m, lc = 1:l)
          varout(lc, mc, nc, j, i) = recvbuf(lc, mc, nc, j-j0+1, i-i0+1, rnk)
        enddo
      enddo
    endif

  endfunction gather_mag_5d
!-----------------------------------------------------------------------
  subroutine bcast_3d(var, l, m, n, root)

    integer, intent(in) :: l, m, n, root
    real(kind=rp), dimension(l, m, n), intent(inout) :: var

    integer :: cnt, ierror
    real(kind=rp), dimension(l, m, n) :: buffer

    cnt = l * m * n

! load to work array
    buffer = var

    call MPI_Bcast(buffer, cnt, mpi_rp, root, mag_world, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Bcast', ierror)

! unpack to model fields
    var = buffer

  endsubroutine bcast_3d
!-----------------------------------------------------------------------
  function reduce_sum_1d(varin, n, root) result(varout)

    integer, intent(in) :: n, root
    real(kind=rp), dimension(n), intent(in) :: varin
    real(kind=rp), dimension(n) :: varout

    integer :: ierror

    if (root < 0) then
      call MPI_Allreduce(varin, varout, n, mpi_rp, &
        MPI_SUM, mag_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allreduce', ierror)
    else
      call MPI_Reduce(varin, varout, n, mpi_rp, &
        MPI_SUM, root, mag_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Reduce', ierror)
    endif

  endfunction reduce_sum_1d
!-----------------------------------------------------------------------
  subroutine finalize

    integer :: ierror

    call MPI_Comm_free(mag_world, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_free', ierror)

  endsubroutine finalize
!-----------------------------------------------------------------------
  subroutine handle_error(funcname, errorcode)

    character(len=*), intent(in) :: funcname
    integer, intent(in) :: errorcode

    character(len=MPI_MAX_ERROR_STRING) :: string
    integer :: resultlen, ierror

    call MPI_Error_string(errorcode, string, resultlen, ierror)
    write(6, "('MPI error encountered: ', a, ', when calling ', a, '. Finalizing...')") &
      trim(string), trim(funcname)
    call MPI_Finalize(ierror)

  endsubroutine handle_error
!-----------------------------------------------------------------------
  pure function generate_minvar_list(summation, nfactor) result(factor_list)
! generate a list of the given length summing up to the given number
! which has the minimum variances (all factors differ at most by 1)
! e.g., 5 out of 26 would be 5,5,5,5,6

    integer, intent(in) :: summation, nfactor
    integer, dimension(nfactor) :: factor_list

    integer :: minfactor, maxfactor, mincnt
    real(kind=rp) :: factor

    factor = real(summation, kind=rp) / nfactor
    minfactor = floor(factor)
    maxfactor = ceiling(factor)

    mincnt = maxfactor*nfactor - summation

    factor_list(1:mincnt) = minfactor
    factor_list(mincnt+1:nfactor) = maxfactor

  endfunction generate_minvar_list
!-----------------------------------------------------------------------
endmodule mpi_mod
