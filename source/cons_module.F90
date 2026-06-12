module cons_module

  use prec,only: rp, sp
  use params_module,only:nmlat_h,nmlon

  implicit none

  real(kind=rp),parameter :: &
    re = 6.37122e6_rp, &      ! earth radius (m)
    pi = 4._rp*atan(1._rp), &
    rtd = 180._rp/pi, &       ! radians to degrees
    dtr = pi/180._rp, &       ! degrees to radians
    h0 = 8.e4_rp, &           ! reference height (m) for dynamo calculations
    r0 = re+h0, &
    ylatm_JT = 45._rp*dtr, &  ! transition latitude where potential becomes symmetric/asymmetric
    phi_pol = 0._rp, &        ! north pole potential
    fill_value = huge(0.0_sp) ! filling value for uninitialized fields

  logical :: read_fac = .false. ! whether FAC is used at high latitude

  integer :: &
    jlatm_JT, & ! latitude index corresponding to the transition latitude
    ndays, &    ! number of days in this year
    ntime       ! number of time steps in this run
  integer,dimension(12) :: days_in_month ! number of days in each month

! lower boundary condition Je2LB is defined by lower atmosphere model
  logical,parameter :: use_lbJ = .false.
  real(kind=rp),dimension(:,:,:),allocatable :: J3LB ! = 0 ! R current [A/m2]

! b_mult is [|Phi|/Delta(Phi)]*(R/L)^2, where Phi is a characteristic potential value,
! Delta(Phi) is a characteristic allowed interhemispheric potential difference,
! R is Earth radius, and L is a characteristic N-S length scale for Phi.
! It is assumed that b_mult is similar for middle and auroral latitudes.
  real(kind=rp),parameter :: b_mult = 1.e3_rp

! pccolatrad is the polar cap colatitude in radians, which for now is fixed.
! But it can be made variable w.r.t. time and magnetic longitude in the future.
  real(kind=rp),parameter :: pccolatrad = 0.25_rp ! 14 degrees
  real(kind=rp),parameter :: rho_pc = sin(pccolatrad)

end module cons_module
