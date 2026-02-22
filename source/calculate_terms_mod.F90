module calculate_terms_mod

  use prec,only:rp

  implicit none

  contains
!-----------------------------------------------------------------------
  pure function calculate_conductance( &
    mlatd0,mlatd1,mlond0,mlond1,nmlat_h, &
    npts,vmp,bmag,sig) result(zig)
! calculate field-line integrated conductance

    use params_mod,only:nhgt_fix
    use cons_mod,only:fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1,nmlat_h
    integer,dimension(nmlat_h),intent(in) :: npts
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: vmp,bmag,sig
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1) :: zig

    integer :: i,j,isn,k
    real(kind=rp) :: sumcond

    zig = fill_value

    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
      sumcond = 0
      do k = 1,npts(j)-1
        sumcond = sumcond+ &
          (sig(k+1,isn,j,i)+sig(k,isn,j,i))* &
          abs(vmp(k+1,isn,j,i)-vmp(k,isn,j,i))/ &
          (bmag(k+1,isn,j,i)+bmag(k,isn,j,i))
      enddo
      zig(isn,j,i) = sumcond
    enddo

  endfunction calculate_conductance
!-----------------------------------------------------------------------
  subroutine calculate_ue( &
    mlatd0,mlatd1,mlond0,mlond1,nmlat_h, &
    npts,un,vn,wn,d1,d2,ue1,ue2)
! calculate winds in the perpendicular directions

    use params_mod,only:nhgt_fix
    use cons_mod,only:fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1,nmlat_h
    integer,dimension(nmlat_h),intent(in) :: npts
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: un,vn,wn
    real(kind=rp),dimension(3,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: d1,d2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: ue1,ue2

    integer :: i,j,isn,k

    ue1 = fill_value
    ue2 = fill_value

    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
       do k = 1,npts(j)
        !  print*,'FVDBG....u,v,w : ',un(k,isn,j,i), vn(k,isn,j,i), wn(k,isn,j,i)
        !  print*,'FVDBG..d1(1..3): ',d1(1,k,isn,j,i),d1(2,k,isn,j,i),d1(3,k,isn,j,i)
        ue1(k,isn,j,i) = un(k,isn,j,i)*d1(1,k,isn,j,i)+ &
                         vn(k,isn,j,i)*d1(2,k,isn,j,i)+ &
                         wn(k,isn,j,i)*d1(3,k,isn,j,i)
        ue2(k,isn,j,i) = un(k,isn,j,i)*d2(1,k,isn,j,i)+ &
                         vn(k,isn,j,i)*d2(2,k,isn,j,i)+ &
                         wn(k,isn,j,i)*d2(3,k,isn,j,i)
      enddo
    enddo

  endsubroutine calculate_ue
!-----------------------------------------------------------------------
  function balance_fac_hl(mlatd0,mlatd1,mlond0,mlond1, &
    M3_p,fac_hl_in_p) result(fac_hl_out_p)
! balance the input high latitude FAC to make sure it is zero
! when integrated over the globe

! this is called when direct FAC is read in (dynamo_fac==.true.)
! and needs to be corrected (direct FAC input may be unbalanced)
! M3_p is only the bottom level

    use params_mod,only:nmlat_h
    use mpi_mod,only:reduce_sum_1d

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      M3_p,fac_hl_in_p
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1) :: fac_hl_out_p

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn
    real(kind=rp) :: downfac
    real(kind=rp),dimension(2) :: subsum,fullsum

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

! calculate the sum of upward and downward FACs
    subsum = 0

! exclude halo points to avoid double counting
    do i = mlon0,mlon1
      do j = max(mlat0,2),min(mlat1,nmlat_h) ! no pole
        do isn = 1,2
          if (fac_hl_in_p(isn,j,i) > 0) &
            subsum(1) = subsum(1)+fac_hl_in_p(isn,j,i)*M3_p(isn,j,i)
          if (fac_hl_in_p(isn,j,i) < 0) &
            subsum(2) = subsum(2)-fac_hl_in_p(isn,j,i)*M3_p(isn,j,i)
        enddo
      enddo
    enddo

    fullsum = reduce_sum_1d(subsum,2,-1)

! initialize FAC to zero (not fill value)
    fac_hl_out_p = 0

! if the integrated upward current is stronger than the downward current
! scale the upward current down globally
    if (fullsum(1) > fullsum(2)) then
      downfac = fullsum(2)/fullsum(1)
      do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=2 .and. j<=nmlat_h)
        if (fac_hl_in_p(isn,j,i) > 0) fac_hl_out_p(isn,j,i) = fac_hl_in_p(isn,j,i)*downfac
        if (fac_hl_in_p(isn,j,i) < 0) fac_hl_out_p(isn,j,i) = fac_hl_in_p(isn,j,i)
      enddo
    endif

! if the integrated upward current is weaker than the downward current
! scale the downward current down globally
    if (fullsum(1) < fullsum(2)) then
      downfac = fullsum(1)/fullsum(2)
      do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=2 .and. j<=nmlat_h)
        if (fac_hl_in_p(isn,j,i) < 0) fac_hl_out_p(isn,j,i) = fac_hl_in_p(isn,j,i)*downfac
        if (fac_hl_in_p(isn,j,i) > 0) fac_hl_out_p(isn,j,i) = fac_hl_in_p(isn,j,i)
      enddo
    endif

  endfunction balance_fac_hl
!-----------------------------------------------------------------------
  function balance_fac_hl_2(mlatd0,mlatd1,mlond0,mlond1, &
    zigP_p,M3_p,fac_hl_in_p) result(fac_hl_out_p)
! balance the input high latitude FAC to make sure it is zero
! when integrated in each hemisphere

! this is called when direct FAC is read in (dynamo_fac==.true.)
! and needs to be corrected (direct FAC input may be unbalanced)
! M3_p is only the bottom level

    use params_mod,only:nmlat_h
    use mpi_mod,only:reduce_sum_1d

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      zigP_p,M3_p,fac_hl_in_p
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1) :: fac_hl_out_p

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn
    real(kind=rp) :: fac
    real(kind=rp),dimension(2) :: sumfac,sumzigP,corr
    real(kind=rp),dimension(4) :: tmp_sub,tmp_full

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

! corr = - zigP*abs(Jmr)/sinI * [sum_i^N Jmr*area] / [sum_i^N abs(Jmr)*zigP/sinI*area]
! sinI is ignored by assuming it is close to 1 at high latitudes
    sumfac = 0
    sumzigP = 0

! exclude halo points to avoid double counting
    do i = mlon0,mlon1
      do j = max(mlat0,2),min(mlat1,nmlat_h) ! no pole
        do isn = 1,2
          fac = fac_hl_in_p(isn,j,i)*M3_p(isn,j,i)
          sumfac(isn) = sumfac(isn)+fac
          sumzigP(isn) = sumzigP(isn)+zigP_p(isn,j,i)*abs(fac)
        enddo
      enddo
    enddo

    tmp_sub(1) = sumfac(1)
    tmp_sub(2) = sumfac(2)
    tmp_sub(3) = sumzigP(1)
    tmp_sub(4) = sumzigP(2)
    tmp_full = reduce_sum_1d(tmp_sub,4,-1)
    sumfac(1) = tmp_full(1)
    sumfac(2) = tmp_full(2)
    sumzigP(1) = tmp_full(3)
    sumzigP(2) = tmp_full(4)

    corr = sumfac/sumzigP

! initialize FAC to zero (not fill value)
    fac_hl_out_p = 0

! correct fac_hl
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=2 .and. j<=nmlat_h)
      fac_hl_out_p(isn,j,i) = fac_hl_in_p(isn,j,i)- &
        zigP_p(isn,j,i)*abs(fac_hl_in_p(isn,j,i))*corr(isn)
    enddo

  endfunction balance_fac_hl_2
!-----------------------------------------------------------------------
  pure subroutine calculate_n( &
    mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
    D_s1,M1_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1, &
    D_s2,M2_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2, &
    N1p_s1,N1h_s1,N2p_s2,N2h_s2)

    use params_mod,only:nhgt_fix,nmlat_h,nmlatS2_h,ylonm,rho,rho_s
    use cons_mod,only:r0,fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_s1
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      D_s1,M1_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1, &
      D_s2,M2_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: &
      N1p_s1,N1h_s1,N2p_s2,N2h_s2

    integer :: i,j,isn,k
    real(kind=rp) :: dlonm,sigC

    N1p_s1 = fill_value
    N1h_s1 = fill_value
    N2p_s2 = fill_value
    N2h_s2 = fill_value

! assume equidistant longitudinal grid points
    dlonm = ylonm(2)-ylonm(1)

! calculate N coefficients at (i+0.5,j,k) for S1 points (no pole or equator)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=2 .and. j<=nmlat_h-1)
      do k = 1,npts_s1(j)

! Equation (6.11)
! N1P(i+0.5) = M1(i+0.5)*[sigP*d1^2](i+0.5)/R/rho(j)/(phi(i+1)-phi(i))
        N1p_s1(k,isn,j,i) = M1_s1(k,isn,j,i)* &
          sigP_s1(k,isn,j,i)*d1d1_s1(k,isn,j,i)/r0/rho(j)/dlonm

! Equation (6.12)
! N1H(i+0.5) = M1(i+0.5)*[sigH*D-sigP*d1*d2](i+0.5)*sqrt(1-0.75*rho(j)^2)/2/R/(rho(j+1)-rho(j-1))
        N1h_s1(k,isn,j,i) = M1_s1(k,isn,j,i)* &
          (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)- &
          sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))* &
          sqrt(1-3*rho(j)**2/4)/2/r0/(rho(j+1)-rho(j-1))
      enddo
    enddo

! Equations (6.16) and (6.17) for the equatorial volume at (i+0.5,j=J,k=1)
! N1P <- M1(i+0.5)*sigC(i+0.5)/R/rho(j)/(phi(i+1)-phi(i))
! sigC = sigP*d1^2+((sigH*D)^2-(sigP*d1*d2)^2)/(sigP*d2^2)
! N1H <- 0
    j = nmlat_h
    k = 1
    if (j>=mlatd0 .and. j<=mlatd1) then
      do concurrent (i = mlond0:mlond1, isn = 1:2)
        sigC = sigP_s1(k,isn,j,i)*d1d1_s1(k,isn,j,i)+ &
          ((sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i))**2- &
          (sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))**2)/ &
          (sigP_s1(k,isn,j,i)*d2d2_s1(k,isn,j,i))
        N1p_s1(k,isn,j,i) = M1_s1(k,isn,j,i)*sigC/r0/rho(j)/dlonm
        N1h_s1(k,isn,j,i) = 0
      enddo
    endif

! calculate N coefficients at (i,j+0.5,k) for S2 points
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlatS2_h)
      do k = 1,npts_s2(j)

! Equation (6.20)
! N2P(j+0.5) = M2(j+0.5)*[sigP*d2^2](j+0.5)*sqrt(1-0.75*rho(j+0.5)^2)/R/(rho(j+1)-rho(j))
        N2p_s2(k,isn,j,i) = M2_s2(k,isn,j,i)* &
          sigP_s2(k,isn,j,i)*d2d2_s2(k,isn,j,i)* &
          sqrt(1-3*rho_s(j)**2/4)/r0/(rho(j+1)-rho(j))

! Equation (6.21)
! N2H(j+0.5) = M2(j+0.5)*[sigH*D+sigP*d1*d2](j+0.5)/2/R/rho(j+0.5)/(phi(i+1)-phi(i-1)))
        N2h_s2(k,isn,j,i) = M2_s2(k,isn,j,i)* &
          (sigH_s2(k,isn,j,i)*D_s2(k,isn,j,i)+ &
          sigP_s2(k,isn,j,i)*d1d2_s2(k,isn,j,i))/ &
          4/r0/rho_s(j)/dlonm
      enddo
    enddo

  endsubroutine calculate_n
!-----------------------------------------------------------------------
  pure subroutine calculate_je( &
    mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2,J3LB, &
    D_s1,be3_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1,ue1_s1,ue2_s1, &
    D_s2,be3_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2,ue1_s2,ue2_s2, &
    Je1D_s1,Je2D_s2)
! calculate wind driven currents

    use params_mod,only:nhgt_fix,nmlat_h,nmlatS2_h
    use cons_mod,only:fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_s1
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: J3LB
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      D_s1,be3_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1,ue1_s1,ue2_s1, &
      D_s2,be3_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2,ue1_s2,ue2_s2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: Je1D_s1,Je2D_s2

    integer :: i,j,isn,k
    real(kind=rp) :: je2d

    Je1D_s1 = fill_value
    Je2D_s2 = fill_value

! Equation (3.8)
! calculate Je1D = sigP*d1^2*ue2*Be3+(sigH*D-sigP*d1*d2)*ue1*Be3 at (i+0.5,j,k) for S1 points
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
      do k = 1,npts_s1(j)
        Je1D_s1(k,isn,j,i) = &
          sigP_s1(k,isn,j,i)*d1d1_s1(k,isn,j,i)* &
          ue2_s1(k,isn,j,i)*be3_s1(k,isn,j,i)+ &
          (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)- &
          sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))* &
          ue1_s1(k,isn,j,i)*be3_s1(k,isn,j,i)
      enddo
    enddo

! Equation (3.16) for the equatorial volume at (i+0.5,j=J,k=1)
! Je1D <- Je1D - (sigH*D-sigP*d1*d2)/(sigP*d2^2)*(Je2LB-Je2D)
    j = nmlat_h
    k = 1
    if (j>=mlatd0 .and. j<=mlatd1) then
      do i = mlond0,mlond1
        do isn = 1,2

! calculate Je2D = (sigH*D+sigP*d1*d2)*ue2*Be3-sigP*d2^2*ue1*Be3 at the equator
! there is no S2 point at the equator therefore it needs to be calculated
          je2d = (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)+ &
            sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))* &
            ue2_s1(k,isn,j,i)*be3_s1(k,isn,j,i)- &
            sigP_s1(k,isn,j,i)*d2d2_s1(k,isn,j,i)* &
            ue1_s1(k,isn,j,i)*be3_s1(k,isn,j,i)

! Je2LB is the current from lower atmosphere, we assume Je2LB = -J3LB
! not exactly since they are half height level apart but should be close
          Je1D_s1(k,isn,j,i) = Je1D_s1(k,isn,j,i)- &
            (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)- &
            sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))/ &
            (sigP_s1(k,isn,j,i)*d2d2_s1(k,isn,j,i))* &
            (-J3LB(isn,j,i)-je2d)
        enddo
      enddo
    endif

! Equation (3.9)
! calculate Je2D = (sigH*D+sigP*d1*d2)*ue2*Be3-sigP*d2^2*ue1*Be3 at (i,j+0.5,k) for S2 points
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlatS2_h)
      do k = 1,npts_s2(j)
        Je2D_s2(k,isn,j,i) = &
          (sigH_s2(k,isn,j,i)*D_s2(k,isn,j,i)+ &
          sigP_s2(k,isn,j,i)*d1d2_s2(k,isn,j,i))* &
          ue2_s2(k,isn,j,i)*be3_s2(k,isn,j,i)- &
          sigP_s2(k,isn,j,i)*d2d2_s2(k,isn,j,i)* &
          ue1_s2(k,isn,j,i)*be3_s2(k,isn,j,i)
      enddo
    enddo

  endsubroutine calculate_je
!-----------------------------------------------------------------------
  pure subroutine calculate_ed(mlatd0,mlatd1,mlond0,mlond1, &
    pot_p,ed1_s1,ed2_s1,ed1_s2,ed2_s2)
! calculate electric fields Ed1, Ed2 at S1 and S2 points

    use params_mod,only:nmlat_h,nmlatS2_h,ylonm,rho,rho_s
    use cons_mod,only:r0,fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: pot_p
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: &
      ed1_s1,ed2_s1,ed1_s2,ed2_s2

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn
    real(kind=rp) :: dlonm

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    ed1_s1 = fill_value
    ed2_s1 = fill_value
    ed1_s2 = fill_value
    ed2_s2 = fill_value

! assume equidistant longitudinal grid points
    dlonm = ylonm(2)-ylonm(1)

! Equations (6.4) and (6.5) (no pole or equator)
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j>=2 .and. j<=nmlat_h-1)
      ed1_s1(isn,j,i) = (pot_p(isn,j,i)-pot_p(isn,j,i+1))/r0/rho(j)/dlonm
      ed2_s1(isn,j,i) = &
        (pot_p(isn,j-1,i)+pot_p(isn,j-1,i+1)- &
         pot_p(isn,j+1,i)-pot_p(isn,j+1,i+1))* &
        sqrt(1-3*rho(j)**2/4)/r0/2/(rho(j+1)-rho(j-1))
    enddo

! Equations (6.8) and (6.9)
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j <= nmlatS2_h)
      ed1_s2(isn,j,i) = &
        (pot_p(isn,j,i-1)+pot_p(isn,j+1,i-1)- &
         pot_p(isn,j,i+1)-pot_p(isn,j+1,i+1))/ &
        4/r0/rho_s(j)/dlonm
      ed2_s2(isn,j,i) = (pot_p(isn,j,i)-pot_p(isn,j+1,i))* &
        sqrt(1-3*rho_s(j)**2/4)/r0/(rho(j+1)-rho(j))
    enddo

  endsubroutine calculate_ed
!-----------------------------------------------------------------------
  pure subroutine calculate_ve( &
    mlatd0,mlatd1,mlond0,mlond1,latbeg,latend, &
    ed1,ed2,be3,ve1,ve2)
! calculate drift velocities ve1, ve2 at S1 and S2 points

! be3 is only the bottom level

    use cons_mod,only:fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1,latbeg,latend
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: ed1,ed2,be3
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: ve1,ve2

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    ve1 = fill_value
    ve2 = fill_value

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j>=latbeg .and. j<=latend)
      ve1(isn,j,i) =  ed2(isn,j,i)/be3(isn,j,i)
      ve2(isn,j,i) = -ed1(isn,j,i)/be3(isn,j,i)
    enddo

  endsubroutine calculate_ve
!-----------------------------------------------------------------------
  pure subroutine calculate_exyz( &
    mlatd0,mlatd1,mlond0,mlond1,latbeg,latend, &
    npts,ed1,ed2,d1,d2,ex,ey,ez)
! calculate electric fields in geographic coordinates (Ed1,2 -> Ex,y,z)

    use params_mod,only:nhgt_fix
    use cons_mod,only:fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1,latbeg,latend
    integer,dimension(latbeg:latend),intent(in) :: npts
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: ed1,ed2
    real(kind=rp),dimension(3,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: d1,d2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: ex,ey,ez

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn,k

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    ex = fill_value
    ey = fill_value
    ez = fill_value

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j>=latbeg .and. j<=latend)
      do k = 1,npts(j)
        ex(k,isn,j,i) = &
          ed1(isn,j,i)*d1(1,k,isn,j,i)+ &
          ed2(isn,j,i)*d2(1,k,isn,j,i)
        ey(k,isn,j,i) = &
          ed1(isn,j,i)*d1(2,k,isn,j,i)+ &
          ed2(isn,j,i)*d2(2,k,isn,j,i)
        ez(k,isn,j,i) = &
          ed1(isn,j,i)*d1(3,k,isn,j,i)+ &
          ed2(isn,j,i)*d2(3,k,isn,j,i)
      enddo
    enddo

  endsubroutine calculate_exyz
!-----------------------------------------------------------------------
  pure subroutine calculate_vxyz( &
    mlatd0,mlatd1,mlond0,mlond1,latbeg,latend, &
    npts,ve1,ve2,e1,e2,vx,vy,vz)
! calculate drift velocities in geographic coordinates (Ve1,2 -> Vx,y,z)

    use params_mod,only:nhgt_fix
    use cons_mod,only:fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1,latbeg,latend
    integer,dimension(latbeg:latend),intent(in) :: npts
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: ve1,ve2
    real(kind=rp),dimension(3,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: e1,e2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: vx,vy,vz

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn,k

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    vx = fill_value
    vy = fill_value
    vz = fill_value

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j>=latbeg .and. j<=latend)
      do k = 1,npts(j)
        vx(k,isn,j,i) = &
          ve1(isn,j,i)*e1(1,k,isn,j,i)+ &
          ve2(isn,j,i)*e2(1,k,isn,j,i)
        vy(k,isn,j,i) = &
          ve1(isn,j,i)*e1(2,k,isn,j,i)+ &
          ve2(isn,j,i)*e2(2,k,isn,j,i)
        vz(k,isn,j,i) = &
          ve1(isn,j,i)*e1(3,k,isn,j,i)+ &
          ve2(isn,j,i)*e2(3,k,isn,j,i)
      enddo
    enddo

  endsubroutine calculate_vxyz
!-----------------------------------------------------------------------
  subroutine calculate_current( &
    mlatd0,mlatd1,mlond0,mlond1, &
    npts_s1,npts_s2,npts_r,J3LB,pot_p, &
    M1_s1,N1p_s1,N1h_s1,Je1D_s1, &
    M2_s2,N2p_s2,N2h_s2,Je2D_s2, &
    M3_r,I1_s1,I2_s2,I3_r)
! calculate area integrated currents

    use params_mod,only:nhgt_fix,nhgt_fix_r,nmlat_h,nmlatS2_h,nmlon
    use cons_mod,only:fill_value
    use mpi_mod,only:gather_mag_lon,sync_mag_lat_5d,sync_mag_lon_5d

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_s1,npts_r
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: J3LB,pot_p
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      M1_s1,N1p_s1,N1h_s1,Je1D_s1,M2_s2,N2p_s2,N2h_s2,Je2D_s2
    real(kind=rp),dimension(nhgt_fix_r,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: M3_r
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: I1_s1,I2_s2
    real(kind=rp),dimension(nhgt_fix_r,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: I3_r

    integer :: mlat0,mlat1,mlon0,mlon1,isn,i,j,k,iconj
    real(kind=rp),dimension(nhgt_fix,2,mlond0:mlond1) :: Je1_sub
    real(kind=rp),dimension(nhgt_fix,2,nmlon) :: Je1_full
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1) :: &
      I1_1,I1_2,I1_3,I2_1,I2_2,I2_3
    real(kind=rp),dimension(2,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1) :: tmpI

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    I1_s1 = fill_value
    I2_s2 = fill_value
    I3_r = fill_value

! Equation (6.10)
! Equations (6.16) to (6.18) for the equatorial volume at (i+0.5,j=J,k=1)
! no need to separate the case since they have already been considered in N1P and Je1D
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j >= 2) ! no pole
      do k = 1,npts_s1(j)
        I1_1(k,isn,j,i) = N1p_s1(k,isn,j,i)*(pot_p(isn,j,i)-pot_p(isn,j,i+1))
        if (j == nmlat_h) then
          I1_2(k,isn,j,i) = 0
        else
          I1_2(k,isn,j,i) = -N1h_s1(k,isn,j,i)* &
            (pot_p(isn,j-1,i)+pot_p(isn,j-1,i+1)- &
             pot_p(isn,j+1,i)-pot_p(isn,j+1,i+1))
        endif
        I1_3(k,isn,j,i) = M1_s1(k,isn,j,i)*Je1D_s1(k,isn,j,i)
        I1_s1(k,isn,j,i) = I1_1(k,isn,j,i)+I1_2(k,isn,j,i)+I1_3(k,isn,j,i)
      enddo
    enddo

! pole value see Equation (6.15)
! Je1(i+0.5,j=1,k) = 0.5*[Je1(i+0.5,j=2,k) - Je1(i'+0.5,j=2,k)]
! i' is the conjugate longitude
! I1(i+0.5,j=1,k) = Je1(i+0.5,j=1,k)*M1(i+0.5,j=1,k)
! we do not have i' in the current process, so first gather all longitudes
    j = mlat0
    if (j == 1) then
      do concurrent (i = mlon0:mlon1, isn = 1:2, k = 1:npts_s1(j))
        Je1_sub(k,isn,i) = I1_s1(k,isn,j+1,i)/M1_s1(k,isn,j+1,i) ! Je1(i+0.5,2,k)
      enddo
      Je1_full = gather_mag_lon(Je1_sub(:,:,mlon0:mlon1),nhgt_fix,2,lat_rank=0)
      do concurrent (i = mlon0:mlon1, isn = 1:2, k = 1:npts_s1(j))
        if (i > nmlon/2) then
          iconj = i-nmlon/2
        else
          iconj = i+nmlon/2
        endif
        I1_s1(k,isn,j,i) = (Je1_full(k,isn,i)-Je1_full(k,isn,iconj))/2*M1_s1(k,isn,j,i)
      enddo
    endif

! Equation (6.19)
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j <= nmlatS2_h)
      do k = 1,npts_s2(j)
        I2_1(k,isn,j,i) = N2p_s2(k,isn,j,i)*(pot_p(isn,j,i)-pot_p(isn,j+1,i))
        I2_2(k,isn,j,i) = N2h_s2(k,isn,j,i)* &
          (pot_p(isn,j,i-1)+pot_p(isn,j+1,i-1)- &
           pot_p(isn,j,i+1)-pot_p(isn,j+1,i+1))
        I2_3(k,isn,j,i) = M2_s2(k,isn,j,i)*Je2D_s2(k,isn,j,i)
        I2_s2(k,isn,j,i) = I2_1(k,isn,j,i)+I2_2(k,isn,j,i)+I2_3(k,isn,j,i)
      enddo
    enddo

! I1 at i-1 and I2 at j-1 are used in I3 calculation, so sync before used
    tmpI(1,:,:,:,:) = I1_s1
    tmpI(2,:,:,:,:) = I2_s2
    call sync_mag_lat_5d(tmpI(:,:,:,:,mlon0:mlon1),2,nhgt_fix,2)
    call sync_mag_lon_5d(tmpI,2,nhgt_fix,2)
    I1_s1 = tmpI(1,:,:,:,:)
    I2_s2 = tmpI(2,:,:,:,:)

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2)

! lower boundary is the current from lower atmosphere J3LB*M3
      k = 1
      I3_r(k,isn,j,i) = J3LB(isn,j,i)*M3_r(k,isn,j,i)

! Equation (6.27)
      do k = 2,npts_r(j)
        I3_r(k,isn,j,i) = I3_r(k-1,isn,j,i)+ &
          I1_s1(k-1,isn,j,i-1)-I1_s1(k-1,isn,j,i)- &
          I2_s2(k-1,isn,j,i)
        if (j /= 1) I3_r(k,isn,j,i) = I3_r(k,isn,j,i)+I2_s2(k-1,isn,j-1,i)
      enddo
    enddo

  endsubroutine calculate_current
!-----------------------------------------------------------------------
endmodule calculate_terms_mod
