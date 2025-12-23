module indices_mod
  use iso_fortran_env, only: r8=>real64

  implicit none

  type indx
     integer :: h = -huge(1)
     integer :: j = -huge(1)
     integer :: k = -huge(1)
  end type indx

  type(indx), protected, allocatable :: s1_ndx(:)
  type(indx), protected, allocatable :: s2_ndx(:)

  integer, protected, allocatable :: flpts1ndx(:,:,:)
  integer, protected, allocatable :: flpts2ndx(:,:,:)

contains

  subroutine indices_init( nflpts_s1, nflpts_s2 )

    use params_module, only: nmlat_h, nmlatS2_h, nhgt_fix
    use dynamo_interface_mod, only: npts_s1, npts_s2

    integer, intent(in) :: nflpts_s1, nflpts_s2

    integer :: j,k,isn
    integer :: ncnt, k0, k1, dk

    character(len=*), parameter :: prefix = 'indices_mod->init: '

    allocate( s1_ndx(nflpts_s1) )
    allocate( s2_ndx(nflpts_s2) )

    allocate(  flpts1ndx(nhgt_fix,2,nmlat_h) )
    allocate(  flpts2ndx(nhgt_fix,2,nmlatS2_h) )

    flpts1ndx = -huge(1)
    flpts2ndx = -huge(1)

    ncnt = 0
    do j = 1,nmlat_h
       do isn = 1,2

          if (isn==1) then
             k0 = 1
             k1 = npts_s1(j)
             dk = 1
          else
             k0 = npts_s1(j)
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             ncnt = ncnt + 1
             flpts1ndx(k,isn,j) = ncnt
             s1_ndx(ncnt)%j = j
             s1_ndx(ncnt)%k = k
             s1_ndx(ncnt)%h = isn
          end do

       end do
    end do

    if (ncnt/=nflpts_s1) then
       write(*,*) prefix//'ERROR: ncnt/=nflpts_s1',ncnt,nflpts_s1
       stop 1
    end if

    ncnt = 0
    do j = 1,nmlatS2_h
       do isn = 1,2

          if (isn==1) then
             k0 = 1
             k1 = npts_s2(j)
             dk = 1
          else
             k0 = npts_s2(j)
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             ncnt = ncnt + 1
             flpts2ndx(k,isn,j) = ncnt
             s2_ndx(ncnt)%j = j
             s2_ndx(ncnt)%k = k
             s2_ndx(ncnt)%h = isn
          end do

       end do

    end do

    if (ncnt/=nflpts_s2) then
       write(*,*) prefix//'ERROR: ncnt/=nflpts_s2',ncnt,nflpts_s2
       stop 1
    end if

  end subroutine indices_init


end module
