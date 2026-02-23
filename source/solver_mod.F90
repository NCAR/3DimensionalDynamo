module solver_mod

  use prec,only:rp

  implicit none

  integer :: nlonlat ! total number of grids to be solved

! convention for grid numbers: longitude increases first then latitude
! so external loop is on latitudes, internal loop is on longitudes

! south pole:
! isn = 1, j = 1, i = 1; then gidx = 1

! southern latitudes (no pole):
! isn = 1, 2 <= j <= nmlat_h; then gidx = (j-2)*nmlon+i+1

! northern high latitudes (no pole):
! isn = 2, 2 <= j <= jlatm_JT-1; then gidx = (nmlat_h+jlatm_JT-2-j)*nmlon+i+1

  contains
!-----------------------------------------------------------------------
subroutine linear_system( mlatd0,mlatd1,mlond0,mlond1, dynamo_pot,dynamo_fac, &
                          bij,pot_hl,fac_hl,src,coef,pot)
! construct linear system based on src and coef and solve in pot

! if potential is read in, pot_hl is used, fac_hl is output
! if FAC is read in, pot_hl is not used, only fac_hl is used
! dynamo_pot and dynamo_fac can not both be true

    use params_mod,only:nmlat_h,nmlon
    use mpi_mod,only:mag_lon_size,mag_rank, &
      gather_mag,gather_mag_lon,bcast_3d

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    logical,intent(in) :: dynamo_pot,dynamo_fac
    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1),intent(in) :: bij
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: pot_hl,src
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(inout) :: fac_hl
    real(kind=rp),dimension(9,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: pot

! if two hemispheres are uncoupled at high latitudes, set beta to zero
    real(kind=rp),parameter :: beta = 0
    integer,parameter :: root = 0
    integer :: mlat0,mlat1,mlon0,mlon1,i,j,j0,j1
    real(kind=rp) :: phi_np ! north pole potential
    real(kind=rp),dimension(nlonlat) :: pot_hl_f
    real(kind=rp),dimension(2,nmlat_h,nmlon) :: pot_hl_full,fac_hl_full
    real(kind=rp),dimension(2,nmlat_h,0:nmlon+1) :: fac_hl_2,pot_2

! the first row is the south pole (only longitude i=1 is effective)
    real(kind=rp) :: src_sp,coef_sp_9
    integer,dimension(0:nmlon) :: jcol1,jcol1_nobij
    real(kind=rp),dimension(0:nmlon) :: nzval1,nzval1_nobij

! LHS matrix and RHS vector in the local subdomain
! reuse these variables after constructing the global matrix
    integer,dimension(mlond0:mlond1) :: jcol1_l
    real(kind=rp),dimension(mlond0:mlond1) :: nzval1_l
    integer,dimension(2,(mlatd1-mlatd0-1)*(mlond1-mlond0-1)) :: gidx_l,rowcnt_l
    integer,dimension(12,2,(mlatd1-mlatd0-1)*(mlond1-mlond0-1)) :: jcol_l
    real(kind=rp),dimension(12,2,(mlatd1-mlatd0-1)*(mlond1-mlond0-1)) :: nzval_l
    real(kind=rp),dimension(2,(mlatd1-mlatd0-1)*(mlond1-mlond0-1)) :: rhs_l

! global matrix after combining all subdomains
    integer,dimension(nlonlat) :: rowcnt,rowcnt_nobij
    integer,dimension(12,2:nlonlat) :: jcol,jcol_nobij
    real(kind=rp),dimension(12,2:nlonlat) :: nzval,nzval_nobij

! RHS vector (dense)
    real(kind=rp),dimension(nlonlat) :: ioncur,magcur,rhs,sol

! LHS matrix in CSR format
    integer :: nnz
    integer,dimension(nlonlat+1) :: rowptr,rowptr_nobij
    integer,dimension(12*nlonlat) :: colind,colind_nobij
    real(kind=rp),dimension(12*nlonlat) :: values_csr,values_nobij_csr

! LHS matrix in CSC format (used by SuperLU)
    integer,dimension(nlonlat+1) :: colptr
    integer,dimension(12*nlonlat) :: rowind
    real(kind=rp),dimension(12*nlonlat) :: values_csc

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    if (root/mag_lon_size /= 0) stop 'root task must include poles'

    if (dynamo_pot) then ! if high latitude potential is read in
! set the north pole potential based on the external potential
      pot_hl_full = gather_mag(pot_hl(:,mlat0:mlat1,mlon0:mlon1),2,root)
      phi_np = pot_hl_full(2,1,1)
    else ! set the north pole potential to zero
      phi_np = 0
      pot_hl_full = phi_np
    endif

    if (dynamo_fac) then
      fac_hl_full = gather_mag(fac_hl(:,mlat0:mlat1,mlon0:mlon1),2,root)
    else
      fac_hl_full = 0
    endif

    src_sp = sum(gather_mag_lon(sum(src(1,mlat0,mlon0:mlon1)),root=root))
    coef_sp_9 = sum(gather_mag_lon(sum(coef(9,1,mlat0,mlon0:mlon1)),root=root))

! source term is ionospheric currents
    call construct_rhs( &
      mlatd0,mlatd1,mlond0,mlond1,src,src_sp+beta*phi_np, &
      sum(coef(6:8,2,mlat0+1,:),dim=1)*phi_np,gidx_l,rhs_l)
    ioncur = gather_rhs(root,(mlat1-mlat0+1)*(mlon1-mlon0+1),gidx_l,rhs_l)

! construct LHS matrix in CSR format
    call construct_lhs(mlatd0,mlatd1,mlond0,mlond1, &
      bij,coef,jcol1_l,nzval1_l,gidx_l,rowcnt_l,jcol_l,nzval_l)
    call gather_lhs(root,mlond0,mlond1,(mlat1-mlat0+1)*(mlon1-mlon0+1), &
      coef_sp_9-beta,jcol1_l,nzval1_l,gidx_l,rowcnt_l,jcol_l,nzval_l, &
      rowcnt,jcol1,nzval1,jcol,nzval)

! A. Maute 2026/01/17: when high latitude potential is used to derive FAC,
! do not include bij in the stencil (these are actually on the right hand side)
    if (dynamo_fac) then
      call construct_lhs_nobij(mlatd0,mlatd1,mlond0,mlond1, &
        coef,jcol1_l,nzval1_l,gidx_l,rowcnt_l,jcol_l,nzval_l)
      call gather_lhs(root,mlond0,mlond1,(mlat1-mlat0+1)*(mlon1-mlon0+1), &
        coef_sp_9,jcol1_l,nzval1_l,gidx_l,rowcnt_l,jcol_l,nzval_l, &
        rowcnt_nobij,jcol1_nobij,nzval1_nobij,jcol_nobij,nzval_nobij)
    endif

    if (mag_rank == root) then

! determine magnetospheric currents (FAC)
! initialize to zero (overwritten later)
      magcur = 0

! input is pot_hl, fac_hl is to be calculated (output)
      if (dynamo_pot) then

! A. Maute 2023/11/21: put the high latitude potential in X
! and then use LHS to calculate FAC on the right hand side
        pot_hl_f = ravel(pot_hl_full)

! construct LHS matrix in CSR format
        call reconstruct_lhs(rowcnt_nobij, &
          jcol1_nobij,nzval1_nobij,jcol_nobij,nzval_nobij, &
          rowptr_nobij,colind_nobij,values_nobij_csr)

! magcur = matmul(lhs, pot_hl)
        magcur = 0
        do i = 1,nlonlat
          do j = rowptr_nobij(i),rowptr_nobij(i+1)-1
            magcur(i) = magcur(i)+values_nobij_csr(j)*pot_hl_f(colind_nobij(j))
          enddo
        enddo

! no need for correction since it is from the divergence of horizontal current

! reconstruct 2D distribution of FAC
        fac_hl_2(:,:,1:nmlon) = unravel(magcur)
        fac_hl_2(2,1,:) = 0 ! north pole is not set in unravel

! add periodic points
        fac_hl_2(:,:,0) = fac_hl_2(:,:,nmlon)
        fac_hl_2(:,:,nmlon+1) = fac_hl_2(:,:,1)
      endif

! input is corrected fac_hl, pot_hl is not used
      if (dynamo_fac) magcur = ravel(fac_hl_full)

! RHS = ionospheric currents + magnetospheric currents (FAC)
      rhs = ioncur+magcur

! construct LHS matrix in CSR format
      call reconstruct_lhs(rowcnt,jcol1,nzval1, &
        jcol,nzval,rowptr,colind,values_csr)
      nnz = rowptr(nlonlat+1)-1

#ifdef USE_MKL
      sol = solve_mkl(nlonlat,nnz,rowptr,colind(1:nnz),values_csr(1:nnz),rhs)
#else
      call csr_to_csc(nlonlat,nlonlat,nnz, &
        rowptr,colind(1:nnz),values_csr(1:nnz), &
        colptr,rowind(1:nnz),values_csc(1:nnz))
      sol = solve_superlu(nlonlat,nnz,colptr,rowind(1:nnz),values_csc(1:nnz),rhs)
#endif

! reconstruct 2D distribution of potential based on the solution
      pot_2(:,:,1:nmlon) = unravel(sol)
      pot_2(2,1,:) = phi_np ! north pole is not set in unravel

! periodic points
      pot_2(:,:,0) = pot_2(:,:,nmlon)
      pot_2(:,:,nmlon+1) = pot_2(:,:,1)
    endif

    j0 = max(mlatd0,1)
    j1 = min(mlatd1,nmlat_h)

    if (dynamo_pot) then
      call bcast_3d(fac_hl_2,2,nmlat_h,nmlon+2,root)

      fac_hl(:,j0:j1,:) = fac_hl_2(:,j0:j1,mlond0:mlond1)
    endif

    call bcast_3d(pot_2,2,nmlat_h,nmlon+2,root)
    pot(:,j0:j1,:) = pot_2(:,j0:j1,mlond0:mlond1)

  endsubroutine linear_system
!-----------------------------------------------------------------------
  pure subroutine construct_lhs(mlatd0,mlatd1,mlond0,mlond1, &
    bij,coef,jcol1,nzval1,gidx,rowcnt,jcol,nzval)
! construct LHS matrix (CSR format) in the local subdomain
! each grid in the lat-lon decomposition corresponds to a row in the matrix
! grids/rows in southern and northern hemispheres are saved separately

! this is not a 9-point stencil

    use params_mod,only:nmlat_h,nmlon
    use cons_mod,only:jlatm_JT

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1),intent(in) :: bij
    real(kind=rp),dimension(9,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef

! rowcnt,jcol,nzval use local indexing, 1:(mlat1-mlat0+1)*(mlon1-mlon0+1)
! gidx stores the corresponding global indices, 1:nlonlat
    integer,dimension(2,(mlatd1-mlatd0-1)*(mlond1-mlond0-1)),intent(out) :: &
      gidx,rowcnt ! number of non-zero elements in each row

! the first row is the south pole (only longitude i=1 is effective)
    integer,dimension(mlond0:mlond1),intent(out) :: jcol1
    real(kind=rp),dimension(mlond0:mlond1),intent(out) :: nzval1

! other rows have at most 12 elements
    integer,dimension(12,2,(mlatd1-mlatd0-1)*(mlond1-mlond0-1)),intent(out) :: jcol
    real(kind=rp),dimension(12,2,(mlatd1-mlatd0-1)*(mlond1-mlond0-1)),intent(out) :: nzval

    integer :: mlat0,mlat1,mlon0,mlon1,cnt,i,j,im,ip
    real(kind=rp) :: c1,c2,c3,c4,c5,c6,c7,c8,c9

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

    gidx = 0
    rowcnt = 0
    jcol1 = 0
    nzval1 = 0
    jcol = 0
    nzval = 0

! start with southern hemisphere

! number of grids/rows in the southern hemisphere
    cnt = 0

    do j = mlat0,mlat1

! south pole has only one point (longitudes are squeezed)
! the volume is a combined cell of all longitudes (a prism with nmlon edges)
      if (j == 1) then
        do i = mlon0,mlon1 ! Sum_i C3S(i) PhiS(2,i)
          jcol1(i) = index_3to1(1,j+1,i)
          nzval1(i) = coef(3,1,j,i)
        enddo
      endif

! j=2 in the southern hemisphere is different from others
! south pole (j=1) has only one point, so coefficients are summed up
      if (j == 2) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(1,cnt) = index_3to1(1,j,i)
          rowcnt(1,cnt) = 8

! C1S     PhiS(j  ,i+1)
! C2S     PhiS(j+1,i+1)
! C3S     PhiS(j+1,i  )
! C4S     PhiS(j+1,i-1)
! C5S     PhiS(j  ,i-1)
! (C9S-b) PhiS(j  ,i  )
! b       PhiN(j  ,i  )

! PhiS(j-1,i-1), PhiS(j-1,i), PhiS(j-1,i+1) are south pole
! so their coefficients C6S, C7S, C8S are summed up
! which becomes (C6S+C7S+C8S) PhiS(j-1,i)
          c1 = coef(1,1,j,i)
          c2 = coef(2,1,j,i)
          c3 = coef(3,1,j,i)
          c4 = coef(4,1,j,i)
          c5 = coef(5,1,j,i)
          c7 = coef(6,1,j,i)+coef(7,1,j,i)+coef(8,1,j,i)
          c9 = coef(9,1,j,i)-bij(j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:8,1,cnt) = &
              (/index_3to1(1,j-1,i), &
                index_3to1(1,j  ,i),index_3to1(1,j  ,ip),index_3to1(1,j  ,im), &
                index_3to1(1,j+1,i),index_3to1(1,j+1,ip),index_3to1(1,j+1,im), &
                index_3to1(2,j  ,i)/)
            nzval(1:8,1,cnt) = &
              (/c7, &
                c9,c1,c5, &
                c3,c2,c4, &
                bij(j,i)/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:8,1,cnt) = &
              (/                                          index_3to1(1,j-1,i), &
                index_3to1(1,j  ,ip),index_3to1(1,j  ,im),index_3to1(1,j  ,i), &
                index_3to1(1,j+1,ip),index_3to1(1,j+1,im),index_3to1(1,j+1,i), &
                                                          index_3to1(2,j  ,i)/)
            nzval(1:8,1,cnt) = &
              (/      c7, &
                c1,c5,c9, &
                c2,c4,c3, &
                bij(j,i)/)
          else
            im = i-1
            ip = i+1
            jcol(1:8,1,cnt) = &
              (/                     index_3to1(1,j-1,i), &
                index_3to1(1,j  ,im),index_3to1(1,j  ,i),index_3to1(1,j  ,ip), &
                index_3to1(1,j+1,im),index_3to1(1,j+1,i),index_3to1(1,j+1,ip), &
                                     index_3to1(2,j  ,i)/)
            nzval(1:8,1,cnt) = &
              (/   c7, &
                c5,c9,c1, &
                c4,c3,c2, &
                bij(j,i)/)
          endif
        enddo
      endif

! southern high latitudes, Equation (7.22)
      if (3<=j .and. j<=jlatm_JT-1) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(1,cnt) = index_3to1(1,j,i)
          rowcnt(1,cnt) = 10

! C1S     PhiS(j  ,i+1)
! C2S     PhiS(j+1,i+1)
! C3S     PhiS(j+1,i  )
! C4S     PhiS(j+1,i-1)
! C5S     PhiS(j  ,i-1)
! C6S     PhiS(j-1,i-1)
! C7S     PhiS(j-1,i  )
! C8S     PhiS(j-1,i+1)
! (C9S-b) PhiS(j  ,i  )
! b       PhiN(j  ,i  )
          c1 = coef(1,1,j,i)
          c2 = coef(2,1,j,i)
          c3 = coef(3,1,j,i)
          c4 = coef(4,1,j,i)
          c5 = coef(5,1,j,i)
          c6 = coef(6,1,j,i)
          c7 = coef(7,1,j,i)
          c8 = coef(8,1,j,i)
          c9 = coef(9,1,j,i)-bij(j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:10,1,cnt) = &
              (/index_3to1(1,j-1,i),index_3to1(1,j-1,ip),index_3to1(1,j-1,im), &
                index_3to1(1,j  ,i),index_3to1(1,j  ,ip),index_3to1(1,j  ,im), &
                index_3to1(1,j+1,i),index_3to1(1,j+1,ip),index_3to1(1,j+1,im), &
                index_3to1(2,j  ,i)/)
            nzval(1:10,1,cnt) = &
              (/c7,c8,c6, &
                c9,c1,c5, &
                c3,c2,c4, &
                bij(j,i)/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:10,1,cnt) = &
              (/index_3to1(1,j-1,ip),index_3to1(1,j-1,im),index_3to1(1,j-1,i), &
                index_3to1(1,j  ,ip),index_3to1(1,j  ,im),index_3to1(1,j  ,i), &
                index_3to1(1,j+1,ip),index_3to1(1,j+1,im),index_3to1(1,j+1,i), &
                                                          index_3to1(2,j  ,i)/)
            nzval(1:10,1,cnt) = &
              (/c8,c6,c7, &
                c1,c5,c9, &
                c2,c4,c3, &
                bij(j,i)/)
          else
            im = i-1
            ip = i+1
            jcol(1:10,1,cnt) = &
              (/index_3to1(1,j-1,im),index_3to1(1,j-1,i),index_3to1(1,j-1,ip), &
                index_3to1(1,j  ,im),index_3to1(1,j  ,i),index_3to1(1,j  ,ip), &
                index_3to1(1,j+1,im),index_3to1(1,j+1,i),index_3to1(1,j+1,ip), &
                                     index_3to1(2,j  ,i)/)
            nzval(1:10,1,cnt) = &
              (/c6,c7,c8, &
                c5,c9,c1, &
                c4,c3,c2, &
                bij(j,i)/)
          endif
        enddo
      endif

! southern transition latitude, Equation (7.13)
      if (j == jlatm_JT) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(1,cnt) = index_3to1(1,j,i)
          rowcnt(1,cnt) = 12

! (C1S+C1N) PhiS(j  ,i+1)
! (C2S+C2N) PhiS(j+1,i+1)
! (C3S+C3N) PhiS(j+1,i  )
! (C4S+C4N) PhiS(j+1,i-1)
! (C5S+C5N) PhiS(j  ,i-1)
! C6S       PhiS(j-1,i-1)
! C7S       PhiS(j-1,i  )
! C8S       PhiS(j-1,i+1)
! (C9S+C9N) PhiS(j  ,i  )
! C6N       PhiN(j-1,i-1)
! C7N       PhiN(j-1,i  )
! C8N       PhiN(j-1,i+1)
          c1 = coef(1,1,j,i)+coef(1,2,j,i)
          c2 = coef(2,1,j,i)+coef(2,2,j,i)
          c3 = coef(3,1,j,i)+coef(3,2,j,i)
          c4 = coef(4,1,j,i)+coef(4,2,j,i)
          c5 = coef(5,1,j,i)+coef(5,2,j,i)
          c6 = coef(6,1,j,i)
          c7 = coef(7,1,j,i)
          c8 = coef(8,1,j,i)
          c9 = coef(9,1,j,i)+coef(9,2,j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:12,1,cnt) = &
              (/index_3to1(1,j-1,i),index_3to1(1,j-1,ip),index_3to1(1,j-1,im), &
                index_3to1(1,j  ,i),index_3to1(1,j  ,ip),index_3to1(1,j  ,im), &
                index_3to1(1,j+1,i),index_3to1(1,j+1,ip),index_3to1(1,j+1,im), &
                index_3to1(2,j-1,i),index_3to1(2,j-1,ip),index_3to1(2,j-1,im)/)
            nzval(1:12,1,cnt) = &
              (/c7,c8,c6, &
                c9,c1,c5, &
                c3,c2,c4, &
                coef(7,2,j,i),coef(8,2,j,i),coef(6,2,j,i)/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:12,1,cnt) = &
              (/index_3to1(1,j-1,ip),index_3to1(1,j-1,im),index_3to1(1,j-1,i), &
                index_3to1(1,j  ,ip),index_3to1(1,j  ,im),index_3to1(1,j  ,i), &
                index_3to1(1,j+1,ip),index_3to1(1,j+1,im),index_3to1(1,j+1,i), &
                index_3to1(2,j-1,ip),index_3to1(2,j-1,im),index_3to1(2,j-1,i)/)
            nzval(1:12,1,cnt) = &
              (/c8,c6,c7, &
                c1,c5,c9, &
                c2,c4,c3, &
                coef(8,2,j,i),coef(6,2,j,i),coef(7,2,j,i)/)
          else
            im = i-1
            ip = i+1
            jcol(1:12,1,cnt) = &
              (/index_3to1(1,j-1,im),index_3to1(1,j-1,i),index_3to1(1,j-1,ip), &
                index_3to1(1,j  ,im),index_3to1(1,j  ,i),index_3to1(1,j  ,ip), &
                index_3to1(1,j+1,im),index_3to1(1,j+1,i),index_3to1(1,j+1,ip), &
                index_3to1(2,j-1,im),index_3to1(2,j-1,i),index_3to1(2,j-1,ip)/)
            nzval(1:12,1,cnt) = &
              (/c6,c7,c8, &
                c5,c9,c1, &
                c4,c3,c2, &
                coef(6,2,j,i),coef(7,2,j,i),coef(8,2,j,i)/)
          endif
        enddo
      endif

! southern low latitudes, Equation (7.14)
      if (jlatm_JT+1<=j .and. j<=nmlat_h-1) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(1,cnt) = index_3to1(1,j,i)
          rowcnt(1,cnt) = 9

! (C1S+C1N) PhiS(j  ,i+1)
! (C2S+C2N) PhiS(j+1,i+1)
! (C3S+C3N) PhiS(j+1,i  )
! (C4S+C4N) PhiS(j+1,i-1)
! (C5S+C5N) PhiS(j  ,i-1)
! (C6S+C6N) PhiS(j-1,i-1)
! (C7S+C7N) PhiS(j-1,i  )
! (C8S+C8N) PhiS(j-1,i+1)
! (C9S+C9N) PhiS(j  ,i  )
          c1 = coef(1,1,j,i)+coef(1,2,j,i)
          c2 = coef(2,1,j,i)+coef(2,2,j,i)
          c3 = coef(3,1,j,i)+coef(3,2,j,i)
          c4 = coef(4,1,j,i)+coef(4,2,j,i)
          c5 = coef(5,1,j,i)+coef(5,2,j,i)
          c6 = coef(6,1,j,i)+coef(6,2,j,i)
          c7 = coef(7,1,j,i)+coef(7,2,j,i)
          c8 = coef(8,1,j,i)+coef(8,2,j,i)
          c9 = coef(9,1,j,i)+coef(9,2,j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:9,1,cnt) = &
              (/index_3to1(1,j-1,i),index_3to1(1,j-1,ip),index_3to1(1,j-1,im), &
                index_3to1(1,j  ,i),index_3to1(1,j  ,ip),index_3to1(1,j  ,im), &
                index_3to1(1,j+1,i),index_3to1(1,j+1,ip),index_3to1(1,j+1,im)/)
            nzval(1:9,1,cnt) = &
              (/c7,c8,c6, &
                c9,c1,c5, &
                c3,c2,c4/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:9,1,cnt) = &
              (/index_3to1(1,j-1,ip),index_3to1(1,j-1,im),index_3to1(1,j-1,i), &
                index_3to1(1,j  ,ip),index_3to1(1,j  ,im),index_3to1(1,j  ,i), &
                index_3to1(1,j+1,ip),index_3to1(1,j+1,im),index_3to1(1,j+1,i)/)
            nzval(1:9,1,cnt) = &
              (/c8,c6,c7, &
                c1,c5,c9, &
                c2,c4,c3/)
          else
            im = i-1
            ip = i+1
            jcol(1:9,1,cnt) = &
              (/index_3to1(1,j-1,im),index_3to1(1,j-1,i),index_3to1(1,j-1,ip), &
                index_3to1(1,j  ,im),index_3to1(1,j  ,i),index_3to1(1,j  ,ip), &
                index_3to1(1,j+1,im),index_3to1(1,j+1,i),index_3to1(1,j+1,ip)/)
            nzval(1:9,1,cnt) = &
              (/c6,c7,c8, &
                c5,c9,c1, &
                c4,c3,c2/)
          endif
        enddo
      endif

! equator, Equation (7.14)
      if (j == nmlat_h) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(1,cnt) = index_3to1(1,j,i)
          rowcnt(1,cnt) = 6

! (C1S+C1N) PhiS(j  ,i+1)
! (C5S+C5N) PhiS(j  ,i-1)
! (C6S+C6N) PhiS(j-1,i-1)
! (C7S+C7N) PhiS(j-1,i  )
! (C8S+C8N) PhiS(j-1,i+1)
! (C9S+C9N) PhiS(j  ,i  )

! there are no PhiS(j+1,i-1), PhiS(j+1,i), PhiS(j+1,i+1)
          c1 = coef(1,1,j,i)+coef(1,2,j,i)
          c5 = coef(5,1,j,i)+coef(5,2,j,i)
          c6 = coef(6,1,j,i)+coef(6,2,j,i)
          c7 = coef(7,1,j,i)+coef(7,2,j,i)
          c8 = coef(8,1,j,i)+coef(8,2,j,i)
          c9 = coef(9,1,j,i)+coef(9,2,j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:6,1,cnt) = &
              (/index_3to1(1,j-1,i),index_3to1(1,j-1,ip),index_3to1(1,j-1,im), &
                index_3to1(1,j  ,i),index_3to1(1,j  ,ip),index_3to1(1,j  ,im)/)
            nzval(1:6,1,cnt) = &
              (/c7,c8,c6, &
                c9,c1,c5/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:6,1,cnt) = &
              (/index_3to1(1,j-1,ip),index_3to1(1,j-1,im),index_3to1(1,j-1,i), &
                index_3to1(1,j  ,ip),index_3to1(1,j  ,im),index_3to1(1,j  ,i)/)
            nzval(1:6,1,cnt) = &
              (/c8,c6,c7, &
                c1,c5,c9/)
          else
            im = i-1
            ip = i+1
            jcol(1:6,1,cnt) = &
              (/index_3to1(1,j-1,im),index_3to1(1,j-1,i),index_3to1(1,j-1,ip), &
                index_3to1(1,j  ,im),index_3to1(1,j  ,i),index_3to1(1,j  ,ip)/)
            nzval(1:6,1,cnt) = &
              (/c6,c7,c8, &
                c5,c9,c1/)
          endif
        enddo
      endif
    enddo

! continue with northern hemisphere

! number of grids/rows in the northern hemisphere
    cnt = 0

    do j = mlat1,mlat0,-1

! skip northern low latitudes (duplicated from southern low latitudes)

! northern high latitudes, Equation (7.22)
      if (3<=j .and. j<=jlatm_JT-1) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(2,cnt) = index_3to1(2,j,i)
          rowcnt(2,cnt) = 10

! C1N     PhiN(j  ,i+1)
! C2N     PhiN(j+1,i+1)
! C3N     PhiN(j+1,i  )
! C4N     PhiN(j+1,i-1)
! C5N     PhiN(j  ,i-1)
! C6N     PhiN(j-1,i-1)
! C7N     PhiN(j-1,i  )
! C8N     PhiN(j-1,i+1)
! (C9N-b) PhiN(j  ,i  )
! b       PhiS(j  ,i  )
          c1 = coef(1,2,j,i)
          c2 = coef(2,2,j,i)
          c3 = coef(3,2,j,i)
          c4 = coef(4,2,j,i)
          c5 = coef(5,2,j,i)
          c6 = coef(6,2,j,i)
          c7 = coef(7,2,j,i)
          c8 = coef(8,2,j,i)
          c9 = coef(9,2,j,i)-bij(j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:10,2,cnt) = &
              (/index_3to1(1,j  ,i), &
                index_3to1(2,j+1,i),index_3to1(2,j+1,ip),index_3to1(2,j+1,im), &
                index_3to1(2,j  ,i),index_3to1(2,j  ,ip),index_3to1(2,j  ,im), &
                index_3to1(2,j-1,i),index_3to1(2,j-1,ip),index_3to1(2,j-1,im)/)
            nzval(1:10,2,cnt) = &
              (/bij(j,i), &
                c3,c2,c4, &
                c9,c1,c5, &
                c7,c8,c6/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:10,2,cnt) = &
              (/                                          index_3to1(1,j  ,i), &
                index_3to1(2,j+1,ip),index_3to1(2,j+1,im),index_3to1(2,j+1,i), &
                index_3to1(2,j  ,ip),index_3to1(2,j  ,im),index_3to1(2,j  ,i), &
                index_3to1(2,j-1,ip),index_3to1(2,j-1,im),index_3to1(2,j-1,i)/)
            nzval(1:10,2,cnt) = &
              (/bij(j,i), &
                c2,c4,c3, &
                c1,c5,c9, &
                c8,c6,c7/)
          else
            im = i-1
            ip = i+1
            jcol(1:10,2,cnt) = &
              (/                     index_3to1(1,j  ,i), &
                index_3to1(2,j+1,im),index_3to1(2,j+1,i),index_3to1(2,j+1,ip), &
                index_3to1(2,j  ,im),index_3to1(2,j  ,i),index_3to1(2,j  ,ip), &
                index_3to1(2,j-1,im),index_3to1(2,j-1,i),index_3to1(2,j-1,ip)/)
            nzval(1:10,2,cnt) = &
              (/bij(j,i), &
                c4,c3,c2, &
                c5,c9,c1, &
                c6,c7,c8/)
          endif
        enddo
      endif

! j=2 in the northern hemisphere is different from others
! north pole (j=1) is set to phi_np, so it is on the right hand side
      if (j == 2) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(2,cnt) = index_3to1(2,j,i)
          rowcnt(2,cnt) = 7

! C1N     PhiN(j  ,i+1)
! C2N     PhiN(j+1,i+1)
! C3N     PhiN(j+1,i  )
! C4N     PhiN(j+1,i-1)
! C5N     PhiN(j  ,i-1)
! (C9N-b) PhiN(j  ,i  )
! b       PhiS(j  ,i  )

! PhiN(j-1,i-1), PhiN(j-1,i), PhiN(j-1,i+1) are north pole
! so they are moved to right hand side
          c1 = coef(1,2,j,i)
          c2 = coef(2,2,j,i)
          c3 = coef(3,2,j,i)
          c4 = coef(4,2,j,i)
          c5 = coef(5,2,j,i)
          c9 = coef(9,2,j,i)-bij(j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:7,2,cnt) = &
              (/index_3to1(1,j  ,i), &
                index_3to1(2,j+1,i),index_3to1(2,j+1,ip),index_3to1(2,j+1,im), &
                index_3to1(2,j  ,i),index_3to1(2,j  ,ip),index_3to1(2,j  ,im)/)
            nzval(1:7,2,cnt) = &
              (/bij(j,i), &
                c3,c2,c4, &
                c9,c1,c5/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:7,2,cnt) = &
              (/                                          index_3to1(1,j  ,i), &
                index_3to1(2,j+1,ip),index_3to1(2,j+1,im),index_3to1(2,j+1,i), &
                index_3to1(2,j  ,ip),index_3to1(2,j  ,im),index_3to1(2,j  ,i)/)
            nzval(1:7,2,cnt) = &
              (/bij(j,i), &
                c2,c4,c3, &
                c1,c5,c9/)
          else
            im = i-1
            ip = i+1
            jcol(1:7,2,cnt) = &
              (/                     index_3to1(1,j  ,i), &
                index_3to1(2,j+1,im),index_3to1(2,j+1,i),index_3to1(2,j+1,ip), &
                index_3to1(2,j  ,im),index_3to1(2,j  ,i),index_3to1(2,j  ,ip)/)
            nzval(1:7,2,cnt) = &
              (/bij(j,i), &
                c4,c3,c2, &
                c5,c9,c1/)
          endif
        enddo
      endif
    enddo

  endsubroutine construct_lhs
!-----------------------------------------------------------------------
  pure subroutine construct_lhs_nobij( &
    mlatd0,mlatd1,mlond0,mlond1, &
    coef,jcol1,nzval1,gidx,rowcnt,jcol,nzval)
! construct LHS matrix (CSR format) in the local subdomain
! each grid in the lat-lon decomposition corresponds to a row in the matrix
! grids/rows in southern and northern hemispheres are saved separately

! this is not a 9-point stencil

    use params_mod,only:nmlat_h,nmlon
    use cons_mod,only:jlatm_JT

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(9,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef

! rowcnt,jcol,nzval use local indexing, 1:(mlat1-mlat0+1)*(mlon1-mlon0+1)
! gidx stores the corresponding global indices, 1:nlonlat
    integer,dimension(2,(mlatd1-mlatd0-1)*(mlond1-mlond0-1)),intent(out) :: &
      gidx,rowcnt ! number of non-zero elements in each row

! the first row is the south pole (only longitude i=1 is effective)
    integer,dimension(mlond0:mlond1),intent(out) :: jcol1
    real(kind=rp),dimension(mlond0:mlond1),intent(out) :: nzval1

! other rows have at most 12 elements
    integer,dimension(12,2,(mlatd1-mlatd0-1)*(mlond1-mlond0-1)),intent(out) :: jcol
    real(kind=rp),dimension(12,2,(mlatd1-mlatd0-1)*(mlond1-mlond0-1)),intent(out) :: nzval

    integer :: mlat0,mlat1,mlon0,mlon1,cnt,i,j,im,ip
    real(kind=rp) :: c1,c2,c3,c4,c5,c6,c7,c8,c9

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

    gidx = 0
    rowcnt = 0
    jcol1 = 0
    nzval1 = 0
    jcol = 0
    nzval = 0

! start with southern hemisphere

! number of grids/rows in the southern hemisphere
    cnt = 0

    do j = mlat0,mlat1

! south pole has only one point (longitudes are squeezed)
! the volume is a combined cell of all longitudes (a prism with nmlon edges)
      if (j == 1) then
        do i = mlon0,mlon1 ! Sum_i C3S(i) PhiS(2,i)
          jcol1(i) = index_3to1(1,j+1,i)
          nzval1(i) = coef(3,1,j,i)
        enddo
      endif

! j=2 in the southern hemisphere is different from others
! south pole (j=1) has only one point, so coefficients are summed up
      if (j == 2) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(1,cnt) = index_3to1(1,j,i)
          rowcnt(1,cnt) = 7

! C1S PhiS(j  ,i+1)
! C2S PhiS(j+1,i+1)
! C3S PhiS(j+1,i  )
! C4S PhiS(j+1,i-1)
! C5S PhiS(j  ,i-1)
! C9S PhiS(j  ,i  )

! PhiS(j-1,i-1), PhiS(j-1,i), PhiS(j-1,i+1) are south pole
! so their coefficients C6S, C7S, C8S are summed up
! which becomes (C6S+C7S+C8S) PhiS(j-1,i)
          c1 = coef(1,1,j,i)
          c2 = coef(2,1,j,i)
          c3 = coef(3,1,j,i)
          c4 = coef(4,1,j,i)
          c5 = coef(5,1,j,i)
          c7 = coef(6,1,j,i)+coef(7,1,j,i)+coef(8,1,j,i)
          c9 = coef(9,1,j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:7,1,cnt) = &
              (/index_3to1(1,j-1,i), &
                index_3to1(1,j  ,i),index_3to1(1,j  ,ip),index_3to1(1,j  ,im), &
                index_3to1(1,j+1,i),index_3to1(1,j+1,ip),index_3to1(1,j+1,im)/)
            nzval(1:7,1,cnt) = &
              (/c7, &
                c9,c1,c5, &
                c3,c2,c4/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:7,1,cnt) = &
              (/                                          index_3to1(1,j-1,i), &
                index_3to1(1,j  ,ip),index_3to1(1,j  ,im),index_3to1(1,j  ,i), &
                index_3to1(1,j+1,ip),index_3to1(1,j+1,im),index_3to1(1,j+1,i)/)
            nzval(1:7,1,cnt) = &
              (/      c7, &
                c1,c5,c9, &
                c2,c4,c3/)
          else
            im = i-1
            ip = i+1
            jcol(1:7,1,cnt) = &
              (/                     index_3to1(1,j-1,i), &
                index_3to1(1,j  ,im),index_3to1(1,j  ,i),index_3to1(1,j  ,ip), &
                index_3to1(1,j+1,im),index_3to1(1,j+1,i),index_3to1(1,j+1,ip)/)
            nzval(1:7,1,cnt) = &
              (/   c7, &
                c5,c9,c1, &
                c4,c3,c2/)
          endif
        enddo
      endif

! southern high latitudes, Equation (7.22)
      if (3<=j .and. j<=jlatm_JT-1) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(1,cnt) = index_3to1(1,j,i)
          rowcnt(1,cnt) = 9

! C1S PhiS(j  ,i+1)
! C2S PhiS(j+1,i+1)
! C3S PhiS(j+1,i  )
! C4S PhiS(j+1,i-1)
! C5S PhiS(j  ,i-1)
! C6S PhiS(j-1,i-1)
! C7S PhiS(j-1,i  )
! C8S PhiS(j-1,i+1)
! C9S PhiS(j  ,i  )
          c1 = coef(1,1,j,i)
          c2 = coef(2,1,j,i)
          c3 = coef(3,1,j,i)
          c4 = coef(4,1,j,i)
          c5 = coef(5,1,j,i)
          c6 = coef(6,1,j,i)
          c7 = coef(7,1,j,i)
          c8 = coef(8,1,j,i)
          c9 = coef(9,1,j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:9,1,cnt) = &
              (/index_3to1(1,j-1,i),index_3to1(1,j-1,ip),index_3to1(1,j-1,im), &
                index_3to1(1,j  ,i),index_3to1(1,j  ,ip),index_3to1(1,j  ,im), &
                index_3to1(1,j+1,i),index_3to1(1,j+1,ip),index_3to1(1,j+1,im)/)
            nzval(1:9,1,cnt) = &
              (/c7,c8,c6, &
                c9,c1,c5, &
                c3,c2,c4/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:9,1,cnt) = &
              (/index_3to1(1,j-1,ip),index_3to1(1,j-1,im),index_3to1(1,j-1,i), &
                index_3to1(1,j  ,ip),index_3to1(1,j  ,im),index_3to1(1,j  ,i), &
                index_3to1(1,j+1,ip),index_3to1(1,j+1,im),index_3to1(1,j+1,i)/)
            nzval(1:9,1,cnt) = &
              (/c8,c6,c7, &
                c1,c5,c9, &
                c2,c4,c3/)
          else
            im = i-1
            ip = i+1
            jcol(1:9,1,cnt) = &
              (/index_3to1(1,j-1,im),index_3to1(1,j-1,i),index_3to1(1,j-1,ip), &
                index_3to1(1,j  ,im),index_3to1(1,j  ,i),index_3to1(1,j  ,ip), &
                index_3to1(1,j+1,im),index_3to1(1,j+1,i),index_3to1(1,j+1,ip)/)
            nzval(1:9,1,cnt) = &
              (/c6,c7,c8, &
                c5,c9,c1, &
                c4,c3,c2/)
          endif
        enddo
      endif

! southern transition latitude, Equation (7.13)
      if (j == jlatm_JT) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(1,cnt) = index_3to1(1,j,i)
          rowcnt(1,cnt) = 12

! (C1S+C1N) PhiS(j  ,i+1)
! (C2S+C2N) PhiS(j+1,i+1)
! (C3S+C3N) PhiS(j+1,i  )
! (C4S+C4N) PhiS(j+1,i-1)
! (C5S+C5N) PhiS(j  ,i-1)
! C6S       PhiS(j-1,i-1)
! C7S       PhiS(j-1,i  )
! C8S       PhiS(j-1,i+1)
! (C9S+C9N) PhiS(j  ,i  )
! C6N       PhiN(j-1,i-1)
! C7N       PhiN(j-1,i  )
! C8N       PhiN(j-1,i+1)
          c1 = coef(1,1,j,i)+coef(1,2,j,i)
          c2 = coef(2,1,j,i)+coef(2,2,j,i)
          c3 = coef(3,1,j,i)+coef(3,2,j,i)
          c4 = coef(4,1,j,i)+coef(4,2,j,i)
          c5 = coef(5,1,j,i)+coef(5,2,j,i)
          c6 = coef(6,1,j,i)
          c7 = coef(7,1,j,i)
          c8 = coef(8,1,j,i)
          c9 = coef(9,1,j,i)+coef(9,2,j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:12,1,cnt) = &
              (/index_3to1(1,j-1,i),index_3to1(1,j-1,ip),index_3to1(1,j-1,im), &
                index_3to1(1,j  ,i),index_3to1(1,j  ,ip),index_3to1(1,j  ,im), &
                index_3to1(1,j+1,i),index_3to1(1,j+1,ip),index_3to1(1,j+1,im), &
                index_3to1(2,j-1,i),index_3to1(2,j-1,ip),index_3to1(2,j-1,im)/)
            nzval(1:12,1,cnt) = &
              (/c7,c8,c6, &
                c9,c1,c5, &
                c3,c2,c4, &
                coef(7,2,j,i),coef(8,2,j,i),coef(6,2,j,i)/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:12,1,cnt) = &
              (/index_3to1(1,j-1,ip),index_3to1(1,j-1,im),index_3to1(1,j-1,i), &
                index_3to1(1,j  ,ip),index_3to1(1,j  ,im),index_3to1(1,j  ,i), &
                index_3to1(1,j+1,ip),index_3to1(1,j+1,im),index_3to1(1,j+1,i), &
                index_3to1(2,j-1,ip),index_3to1(2,j-1,im),index_3to1(2,j-1,i)/)
            nzval(1:12,1,cnt) = &
              (/c8,c6,c7, &
                c1,c5,c9, &
                c2,c4,c3, &
                coef(8,2,j,i),coef(6,2,j,i),coef(7,2,j,i)/)
          else
            im = i-1
            ip = i+1
            jcol(1:12,1,cnt) = &
              (/index_3to1(1,j-1,im),index_3to1(1,j-1,i),index_3to1(1,j-1,ip), &
                index_3to1(1,j  ,im),index_3to1(1,j  ,i),index_3to1(1,j  ,ip), &
                index_3to1(1,j+1,im),index_3to1(1,j+1,i),index_3to1(1,j+1,ip), &
                index_3to1(2,j-1,im),index_3to1(2,j-1,i),index_3to1(2,j-1,ip)/)
            nzval(1:12,1,cnt) = &
              (/c6,c7,c8, &
                c5,c9,c1, &
                c4,c3,c2, &
                coef(6,2,j,i),coef(7,2,j,i),coef(8,2,j,i)/)
          endif
        enddo
      endif

! southern low latitudes, Equation (7.14)
      if (jlatm_JT+1<=j .and. j<=nmlat_h-1) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(1,cnt) = index_3to1(1,j,i)
          rowcnt(1,cnt) = 9

! (C1S+C1N) PhiS(j  ,i+1)
! (C2S+C2N) PhiS(j+1,i+1)
! (C3S+C3N) PhiS(j+1,i  )
! (C4S+C4N) PhiS(j+1,i-1)
! (C5S+C5N) PhiS(j  ,i-1)
! (C6S+C6N) PhiS(j-1,i-1)
! (C7S+C7N) PhiS(j-1,i  )
! (C8S+C8N) PhiS(j-1,i+1)
! (C9S+C9N) PhiS(j  ,i  )
          c1 = coef(1,1,j,i)+coef(1,2,j,i)
          c2 = coef(2,1,j,i)+coef(2,2,j,i)
          c3 = coef(3,1,j,i)+coef(3,2,j,i)
          c4 = coef(4,1,j,i)+coef(4,2,j,i)
          c5 = coef(5,1,j,i)+coef(5,2,j,i)
          c6 = coef(6,1,j,i)+coef(6,2,j,i)
          c7 = coef(7,1,j,i)+coef(7,2,j,i)
          c8 = coef(8,1,j,i)+coef(8,2,j,i)
          c9 = coef(9,1,j,i)+coef(9,2,j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:9,1,cnt) = &
              (/index_3to1(1,j-1,i),index_3to1(1,j-1,ip),index_3to1(1,j-1,im), &
                index_3to1(1,j  ,i),index_3to1(1,j  ,ip),index_3to1(1,j  ,im), &
                index_3to1(1,j+1,i),index_3to1(1,j+1,ip),index_3to1(1,j+1,im)/)
            nzval(1:9,1,cnt) = &
              (/c7,c8,c6, &
                c9,c1,c5, &
                c3,c2,c4/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:9,1,cnt) = &
              (/index_3to1(1,j-1,ip),index_3to1(1,j-1,im),index_3to1(1,j-1,i), &
                index_3to1(1,j  ,ip),index_3to1(1,j  ,im),index_3to1(1,j  ,i), &
                index_3to1(1,j+1,ip),index_3to1(1,j+1,im),index_3to1(1,j+1,i)/)
            nzval(1:9,1,cnt) = &
              (/c8,c6,c7, &
                c1,c5,c9, &
                c2,c4,c3/)
          else
            im = i-1
            ip = i+1
            jcol(1:9,1,cnt) = &
              (/index_3to1(1,j-1,im),index_3to1(1,j-1,i),index_3to1(1,j-1,ip), &
                index_3to1(1,j  ,im),index_3to1(1,j  ,i),index_3to1(1,j  ,ip), &
                index_3to1(1,j+1,im),index_3to1(1,j+1,i),index_3to1(1,j+1,ip)/)
            nzval(1:9,1,cnt) = &
              (/c6,c7,c8, &
                c5,c9,c1, &
                c4,c3,c2/)
          endif
        enddo
      endif

! equator, Equation (7.14)
      if (j == nmlat_h) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(1,cnt) = index_3to1(1,j,i)
          rowcnt(1,cnt) = 6

! (C1S+C1N) PhiS(j  ,i+1)
! (C5S+C5N) PhiS(j  ,i-1)
! (C6S+C6N) PhiS(j-1,i-1)
! (C7S+C7N) PhiS(j-1,i  )
! (C8S+C8N) PhiS(j-1,i+1)
! (C9S+C9N) PhiS(j  ,i  )

! there are no PhiS(j+1,i-1), PhiS(j+1,i), PhiS(j+1,i+1)
          c1 = coef(1,1,j,i)+coef(1,2,j,i)
          c5 = coef(5,1,j,i)+coef(5,2,j,i)
          c6 = coef(6,1,j,i)+coef(6,2,j,i)
          c7 = coef(7,1,j,i)+coef(7,2,j,i)
          c8 = coef(8,1,j,i)+coef(8,2,j,i)
          c9 = coef(9,1,j,i)+coef(9,2,j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:6,1,cnt) = &
              (/index_3to1(1,j-1,i),index_3to1(1,j-1,ip),index_3to1(1,j-1,im), &
                index_3to1(1,j  ,i),index_3to1(1,j  ,ip),index_3to1(1,j  ,im)/)
            nzval(1:6,1,cnt) = &
              (/c7,c8,c6, &
                c9,c1,c5/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:6,1,cnt) = &
              (/index_3to1(1,j-1,ip),index_3to1(1,j-1,im),index_3to1(1,j-1,i), &
                index_3to1(1,j  ,ip),index_3to1(1,j  ,im),index_3to1(1,j  ,i)/)
            nzval(1:6,1,cnt) = &
              (/c8,c6,c7, &
                c1,c5,c9/)
          else
            im = i-1
            ip = i+1
            jcol(1:6,1,cnt) = &
              (/index_3to1(1,j-1,im),index_3to1(1,j-1,i),index_3to1(1,j-1,ip), &
                index_3to1(1,j  ,im),index_3to1(1,j  ,i),index_3to1(1,j  ,ip)/)
            nzval(1:6,1,cnt) = &
              (/c6,c7,c8, &
                c5,c9,c1/)
          endif
        enddo
      endif
    enddo

! continue with northern hemisphere

! number of grids/rows in the northern hemisphere
    cnt = 0

    do j = mlat1,mlat0,-1

! skip northern low latitudes (duplicated from southern low latitudes)

! northern high latitudes, Equation (7.22)
      if (3<=j .and. j<=jlatm_JT-1) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(2,cnt) = index_3to1(2,j,i)
          rowcnt(2,cnt) = 9

! C1N PhiN(j  ,i+1)
! C2N PhiN(j+1,i+1)
! C3N PhiN(j+1,i  )
! C4N PhiN(j+1,i-1)
! C5N PhiN(j  ,i-1)
! C6N PhiN(j-1,i-1)
! C7N PhiN(j-1,i  )
! C8N PhiN(j-1,i+1)
! C9N PhiN(j  ,i  )
          c1 = coef(1,2,j,i)
          c2 = coef(2,2,j,i)
          c3 = coef(3,2,j,i)
          c4 = coef(4,2,j,i)
          c5 = coef(5,2,j,i)
          c6 = coef(6,2,j,i)
          c7 = coef(7,2,j,i)
          c8 = coef(8,2,j,i)
          c9 = coef(9,2,j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:9,2,cnt) = &
              (/index_3to1(2,j+1,i),index_3to1(2,j+1,ip),index_3to1(2,j+1,im), &
                index_3to1(2,j  ,i),index_3to1(2,j  ,ip),index_3to1(2,j  ,im), &
                index_3to1(2,j-1,i),index_3to1(2,j-1,ip),index_3to1(2,j-1,im)/)
            nzval(1:9,2,cnt) = &
              (/c3,c2,c4, &
                c9,c1,c5, &
                c7,c8,c6/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:9,2,cnt) = &
              (/index_3to1(2,j+1,ip),index_3to1(2,j+1,im),index_3to1(2,j+1,i), &
                index_3to1(2,j  ,ip),index_3to1(2,j  ,im),index_3to1(2,j  ,i), &
                index_3to1(2,j-1,ip),index_3to1(2,j-1,im),index_3to1(2,j-1,i)/)
            nzval(1:9,2,cnt) = &
              (/c2,c4,c3, &
                c1,c5,c9, &
                c8,c6,c7/)
          else
            im = i-1
            ip = i+1
            jcol(1:9,2,cnt) = &
              (/index_3to1(2,j+1,im),index_3to1(2,j+1,i),index_3to1(2,j+1,ip), &
                index_3to1(2,j  ,im),index_3to1(2,j  ,i),index_3to1(2,j  ,ip), &
                index_3to1(2,j-1,im),index_3to1(2,j-1,i),index_3to1(2,j-1,ip)/)
            nzval(1:9,2,cnt) = &
              (/c4,c3,c2, &
                c5,c9,c1, &
                c6,c7,c8/)
          endif
        enddo
      endif

! j=2 in the northern hemisphere is different from others
! north pole (j=1) is set to phi_np, so it is on the right hand side
      if (j == 2) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(2,cnt) = index_3to1(2,j,i)
          rowcnt(2,cnt) = 6

! C1N PhiN(j  ,i+1)
! C2N PhiN(j+1,i+1)
! C3N PhiN(j+1,i  )
! C4N PhiN(j+1,i-1)
! C5N PhiN(j  ,i-1)
! C9N PhiN(j  ,i  )

! PhiN(j-1,i-1), PhiN(j-1,i), PhiN(j-1,i+1) are north pole
! so they are moved to right hand side
          c1 = coef(1,2,j,i)
          c2 = coef(2,2,j,i)
          c3 = coef(3,2,j,i)
          c4 = coef(4,2,j,i)
          c5 = coef(5,2,j,i)
          c9 = coef(9,2,j,i)

          if (i == 1) then
            im = nmlon
            ip = i+1
            jcol(1:6,2,cnt) = &
              (/index_3to1(2,j+1,i),index_3to1(2,j+1,ip),index_3to1(2,j+1,im), &
                index_3to1(2,j  ,i),index_3to1(2,j  ,ip),index_3to1(2,j  ,im)/)
            nzval(1:6,2,cnt) = &
              (/c3,c2,c4, &
                c9,c1,c5/)
          elseif (i == nmlon) then
            im = i-1
            ip = 1
            jcol(1:6,2,cnt) = &
              (/index_3to1(2,j+1,ip),index_3to1(2,j+1,im),index_3to1(2,j+1,i), &
                index_3to1(2,j  ,ip),index_3to1(2,j  ,im),index_3to1(2,j  ,i)/)
            nzval(1:6,2,cnt) = &
              (/c2,c4,c3, &
                c1,c5,c9/)
          else
            im = i-1
            ip = i+1
            jcol(1:6,2,cnt) = &
              (/index_3to1(2,j+1,im),index_3to1(2,j+1,i),index_3to1(2,j+1,ip), &
                index_3to1(2,j  ,im),index_3to1(2,j  ,i),index_3to1(2,j  ,ip)/)
            nzval(1:6,2,cnt) = &
              (/c4,c3,c2, &
                c5,c9,c1/)
          endif
        enddo
      endif
    enddo

  endsubroutine construct_lhs_nobij
!-----------------------------------------------------------------------
  subroutine gather_lhs(root,mlond0,mlond1,n, &
    coef_sp_9,jcol1_l,nzval1_l, &
    gidx_l,rowcnt_l,jcol_l,nzval_l, &
    rowcnt,jcol1,nzval1,jcol,nzval)

    use params_mod,only:nmlon
    use mpi_mod,only:mag_maxlat,mag_maxlon,mag_size,mag_lon_size, &
      mag_world,mag_rank,mag_lat_rank,mag_lon0_task,mag_lon1_task, &
      mpi_rp,mag_nlat_task,mag_nlon_task,gather_mag_lon,handle_error
    use mpi

! LHS matrix in the local subdomain
    integer,intent(in) :: root,mlond0,mlond1,n
    real(kind=rp),intent(in) :: coef_sp_9
    integer,dimension(mlond0:mlond1),intent(in) :: jcol1_l
    real(kind=rp),dimension(mlond0:mlond1),intent(in) :: nzval1_l
    integer,dimension(2,n),intent(in) :: gidx_l,rowcnt_l
    integer,dimension(12,2,n),intent(in) :: jcol_l
    real(kind=rp),dimension(12,2,n),intent(in) :: nzval_l

! reconstructed global LHS matrix
    integer,dimension(nlonlat),intent(out) :: rowcnt
    integer,dimension(0:nmlon),intent(out) :: jcol1
    real(kind=rp),dimension(0:nmlon),intent(out) :: nzval1
    integer,dimension(12,2:nlonlat),intent(out) :: jcol
    real(kind=rp),dimension(12,2:nlonlat),intent(out) :: nzval

    integer :: mlon0,mlon1,lat_rank,cnt,i,j,isn,rnk,rnki,rnkj,i0,i1,ielem,idx,ierror
    integer,dimension(mag_maxlon) :: sendbuf1_i
    integer,dimension(mag_maxlon,0:mag_lon_size-1) :: recvbuf1_i
    integer,dimension(14,2,mag_maxlat*mag_maxlon) :: sendbuf_i
    integer,dimension(14,2,mag_maxlat*mag_maxlon,0:mag_size-1) :: recvbuf_i
    real(kind=rp),dimension(12,2,mag_maxlat*mag_maxlon) :: sendbuf_r
    real(kind=rp),dimension(12,2,mag_maxlat*mag_maxlon,0:mag_size-1) :: recvbuf_r
    integer,dimension(0:mag_lon_size) :: request

    mlon0 = mlond0+1
    mlon1 = mlond1-1
    lat_rank = root/mag_lon_size

    if (mag_lat_rank == lat_rank) then
      cnt = mag_maxlon
      do i = 1, mlon1-mlon0+1
        sendbuf1_i(i) = jcol1_l(i+mlon0-1)
      enddo

      call MPI_Isend(sendbuf1_i, cnt, MPI_INTEGER, &
        root, 0, mag_world, request(mag_lon_size), ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

      if (mag_rank == root) then
        do rnki = 0, mag_lon_size-1
          rnk = lat_rank*mag_lon_size + rnki

          call MPI_Irecv(recvbuf1_i(:, rnki), cnt, MPI_INTEGER, &
            rnk, 0, mag_world, request(rnki), ierror)
          if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
        enddo

        call MPI_Waitall(mag_lon_size+1, request, MPI_STATUSES_IGNORE, ierror)
        if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

        do rnki = 0, mag_lon_size-1
          rnk = lat_rank*mag_lon_size + rnki
          i0 = mag_lon0_task(rnk)
          i1 = mag_lon1_task(rnk)
          do i = i0,i1
            jcol1(i) = recvbuf1_i(i-i0+1, rnki)
          enddo
        enddo
      else
        call MPI_Wait(request(mag_lon_size), MPI_STATUS_IGNORE, ierror)
        if (ierror /= MPI_SUCCESS) call handle_error('MPI_Wait', ierror)
      endif
    endif

    nzval1(1:nmlon) = gather_mag_lon(nzval1_l(mlon0:mlon1),root=root)

    sendbuf_i = 0
    sendbuf_r = 0
    do concurrent (i = 1:n, isn = 1:2)
      sendbuf_i(1,isn,i) = gidx_l(isn,i)
      sendbuf_i(2,isn,i) = rowcnt_l(isn,i)
      sendbuf_i(3:14,isn,i) = jcol_l(:,isn,i)
      sendbuf_r(:,isn,i) = nzval_l(:,isn,i)
    enddo

    recvbuf_i = 0
    recvbuf_r = 0

    cnt = 14*2*mag_maxlat*mag_maxlon
    call MPI_Gather(sendbuf_i, cnt, MPI_INTEGER, &
      recvbuf_i, cnt, MPI_INTEGER, root, mag_world, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)

    cnt = 12*2*mag_maxlat*mag_maxlon
    call MPI_Gather(sendbuf_r, cnt, mpi_rp, &
      recvbuf_r, cnt, mpi_rp, root, mag_world, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)

    if (mag_rank == root) then
      do rnk = 0,mag_size-1
        rnkj = rnk/mag_lon_size
        rnki = modulo(rnk,mag_lon_size)
        do ielem = 1,mag_nlat_task(rnkj)*mag_nlon_task(rnki)
          do isn = 1,2
            idx = recvbuf_i(1,isn,ielem,rnk)
            if (2<=idx .and. idx<=nlonlat) then
              rowcnt(idx) = recvbuf_i(2,isn,ielem,rnk)
              jcol(:,idx) = recvbuf_i(3:14,isn,ielem,rnk)
              nzval(:,idx) = recvbuf_r(:,isn,ielem,rnk)
            endif
          enddo
        enddo
      enddo

! complete south pole

! south pole has only one point (longitudes are squeezed)
! the volume is a combined cell of all longitudes (a prism with nmlon edges)
      j = 1

! Equation (7.23) is applied to i=1
      idx = index_3to1(1,j,1)
      rowcnt(idx) = nmlon+1

! A. Richmond 2023/06/20: add b effects
! (Sum_i C9S(i) - beta) PhiS(1,1)
! for potential input, no b effect: Sum_i C9S(i) PhiS(1,1)
      jcol1(0) = idx
      nzval1(0) = coef_sp_9
    endif

  endsubroutine gather_lhs
!-----------------------------------------------------------------------
  pure subroutine reconstruct_lhs(rowcnt,jcol1,nzval1,jcol,nzval,rowptr,colind,values)

    use params_mod,only:nmlon

    integer,dimension(nlonlat),intent(in) :: rowcnt
    integer,dimension(nmlon+1),intent(in) :: jcol1
    real(kind=rp),dimension(nmlon+1),intent(in) :: nzval1
    integer,dimension(12,2:nlonlat),intent(in) :: jcol
    real(kind=rp),dimension(12,2:nlonlat),intent(in) :: nzval
    integer,dimension(nlonlat+1),intent(out) :: rowptr
    integer,dimension(12*nlonlat),intent(out) :: colind
    real(kind=rp),dimension(12*nlonlat),intent(out) :: values

    integer :: i,j

    rowptr(1) = 1
    do i = 2,nlonlat+1
      rowptr(i) = rowptr(i-1)+rowcnt(i-1)
    enddo

    colind = 0
    values = 0

    do i = 1,rowcnt(1)
      colind(i) = jcol1(i)
      values(i) = nzval1(i)
    enddo

    do j = 2,nlonlat
      do i = 1,rowcnt(j)
        colind(rowptr(j)+i-1) = jcol(i,j)
        values(rowptr(j)+i-1) = nzval(i,j)
      enddo
    enddo

  endsubroutine reconstruct_lhs
!-----------------------------------------------------------------------
  pure subroutine construct_rhs( &
    mlatd0,mlatd1,mlond0,mlond1, &
    src,src_sp,src_np,gidx,rhs)
! construct RHS vector
! this is different from ravel

    use params_mod,only:nmlat_h
    use cons_mod,only:jlatm_JT

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: src
    real(kind=rp),intent(in) :: src_sp
    real(kind=rp),dimension(mlond0:mlond1),intent(in) :: src_np
    integer,dimension(2,(mlatd1-mlatd0-1)*(mlond1-mlond0-1)),intent(out) :: gidx
    real(kind=rp),dimension(2,(mlatd1-mlatd0-1)*(mlond1-mlond0-1)),intent(out) :: rhs

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,cnt

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

    gidx = 0
    rhs = 0

! start with southern hemisphere

! number of grids/rows in the southern hemisphere
    cnt = 0

    do j = mlat0,mlat1

! south pole has only one point (longitudes are squeezed)
! the volume is a combined cell of all longitudes (a prism with nmlon edges)
      if (j==1 .and. mlon0==1) then ! Equation (7.23) is applied to i=1
        cnt = cnt+1
        gidx(1,cnt) = index_3to1(1,j,1)

! Sum_i S(1,i) + beta PhiN(1,1)
! north pole is set to phi_np, so it is on the right hand side
        rhs(1,cnt) = -src_sp
      endif

! southern high latitudes, Equation (7.22)
      if (2<=j .and. j<=jlatm_JT-1) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(1,cnt) = index_3to1(1,j,i)
          rhs(1,cnt) = -src(1,j,i)
        enddo
      endif

! southern low latitudes, Equation (7.14)
      if (jlatm_JT<=j .and. j<=nmlat_h) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(1,cnt) = index_3to1(1,j,i)
          rhs(1,cnt) = -(src(1,j,i)+src(2,j,i))
        enddo
      endif
    enddo

! continue with northern hemisphere

! number of grids/rows in the northern hemisphere
    cnt = 0

    do j = mlat1,mlat0,-1

! skip northern low latitudes (duplicated from southern low latitudes)

! northern high latitudes, Equation (7.22)
      if (3<=j .and. j<=jlatm_JT-1) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(2,cnt) = index_3to1(2,j,i)
          rhs(2,cnt) = -src(2,j,i)
        enddo
      endif

! j=2 in the northern hemisphere is different from others
! north pole (j=1) is set to phi_np, so it is on the right hand side
! PhiN(j-1,i-1), PhiN(j-1,i), PhiN(j-1,i+1) are north pole
! so they are moved to right hand side
      if (j == 2) then
        do i = mlon0,mlon1
          cnt = cnt+1
          gidx(2,cnt) = index_3to1(2,j,i)
          rhs(2,cnt) = -src(2,j,i)-src_np(i)
        enddo
      endif
    enddo

  endsubroutine construct_rhs
!-----------------------------------------------------------------------
  function gather_rhs(root,n,gidx_l,rhs_l) result(rhs)

    use mpi_mod,only:mag_maxlat,mag_maxlon,mag_size,mag_world,mpi_rp, &
      mag_rank,mag_lon_size,mag_nlat_task,mag_nlon_task,handle_error
    use mpi

! RHS in the local subdomain
    integer,intent(in) :: root,n
    integer,dimension(2,n),intent(in) :: gidx_l
    real(kind=rp),dimension(2,n),intent(in) :: rhs_l

! reconstructed global RHS vector
    real(kind=rp),dimension(nlonlat) :: rhs

    integer :: i,isn,cnt,rnk,rnki,rnkj,ielem,idx,ierror
    integer,dimension(2,mag_maxlat*mag_maxlon) :: sendbuf_i
    integer,dimension(2,mag_maxlat*mag_maxlon,0:mag_size-1) :: recvbuf_i
    real(kind=rp),dimension(2,mag_maxlat*mag_maxlon) :: sendbuf_r
    real(kind=rp),dimension(2,mag_maxlat*mag_maxlon,0:mag_size-1) :: recvbuf_r

    sendbuf_i = 0
    sendbuf_r = 0
    do concurrent (i = 1:n, isn = 1:2)
      sendbuf_i(isn,i) = gidx_l(isn,i)
      sendbuf_r(isn,i) = rhs_l(isn,i)
    enddo

    recvbuf_i = 0
    recvbuf_r = 0

    cnt = 2*mag_maxlat*mag_maxlon

    call MPI_Gather(sendbuf_i, cnt, MPI_INTEGER, &
      recvbuf_i, cnt, MPI_INTEGER, root, mag_world, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)

    call MPI_Gather(sendbuf_r, cnt, mpi_rp, &
      recvbuf_r, cnt, mpi_rp, root, mag_world, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)

    if (mag_rank == root) then
      do rnk = 0,mag_size-1
        rnkj = rnk/mag_lon_size
        rnki = modulo(rnk,mag_lon_size)
        do ielem = 1,mag_nlat_task(rnkj)*mag_nlon_task(rnki)
          do isn = 1,2
            idx = recvbuf_i(isn,ielem,rnk)
            if (1<=idx .and. idx<=nlonlat) rhs(idx) = recvbuf_r(isn,ielem,rnk)
          enddo
        enddo
      enddo
    endif

  endfunction gather_rhs
!-----------------------------------------------------------------------
  pure function index_3to1(isn,j,i) result(gidx)
! convert lon-lat-hemi triplet (isn,j,i) to linear index gidx
! this is the inverse of index_1to3

! isn,j,i need to satisfy 1<=isn<=2, 1<=j<=nmlat_h, 1<=i<=nmlon
! but those conditions are not enforced

    use params_mod,only:nmlat_h,nmlon
    use cons_mod,only:jlatm_JT

    integer,intent(in) :: isn,j,i
    integer :: gidx

    if (isn == 1) then
      if (j == 1) then
        gidx = 1
      else
        gidx = (j-2)*nmlon+i+1
      endif
    else
      if (j == 1) then
        gidx = 0
      elseif (j <= jlatm_JT-1) then
        gidx = (nmlat_h+jlatm_JT-2-j)*nmlon+i+1
      else
        gidx = (j-2)*nmlon+i+1
      endif
    endif

  endfunction index_3to1
!-----------------------------------------------------------------------
  pure subroutine index_1to3(gidx,isn,j,i)
! convert linear index gidx to lon-lat-hemi triplet (isn,j,i)
! this is the inverse of index_3to1

! gidx needs to be within [1, nlonlat]
! but it is not enforced

    use params_mod,only:nmlat_h,nmlon
    use cons_mod,only:jlatm_JT

    integer,intent(in) :: gidx
    integer,intent(out) :: isn,j,i

    if (gidx == 1) then
      isn = 1
      j = 1
      i = 1
    elseif (gidx <= (nmlat_h-1)*nmlon+1) then
      isn = 1
      j = (gidx-1)/nmlon + 2
      i = modulo(gidx-1,nmlon)
      if (i == 0) then
        j = j-1
        i = nmlon
      endif
    else
      isn = 2
      j = nmlat_h+jlatm_JT-2 - (gidx-1)/nmlon
      i = modulo(gidx-1,nmlon)
      if (i == 0) then
        j = j+1
        i = nmlon
      endif
    endif

  endsubroutine index_1to3
!-----------------------------------------------------------------------
  pure function ravel(fin) result(fout)
! reorder 2D fields (lat-lon) into 1D vector (RHS)
! this is the inverse of unravel (except for the north pole)

    use params_mod,only:nmlat_h,nmlon
    use cons_mod,only:jlatm_JT

    real(kind=rp),dimension(2,nmlat_h,nmlon),intent(in) :: fin
    real(kind=rp),dimension(nlonlat) :: fout

    integer :: i,j,gidx

! start with southern hemisphere

! south pole has only one point (longitudes are squeezed)
    j = 1
    gidx = index_3to1(1,j,1)
    fout(gidx) = fin(1,j,1)

! southern latitudes (no pole)
    do concurrent (i = 1:nmlon, j = 2:nmlat_h)
      gidx = index_3to1(1,j,i)
      fout(gidx) = fin(1,j,i)
    enddo

! continue with northern hemisphere

! skip northern low latitudes (duplicated from southern low latitudes)

! northern high latitudes (no pole)
    do concurrent (i = 1:nmlon, j = 2:jlatm_JT-1)
      gidx = index_3to1(2,j,i)
      fout(gidx) = fin(2,j,i)
    enddo

  endfunction ravel
!-----------------------------------------------------------------------
  pure function unravel(fin) result(fout)
! reorder 1D vector (RHS) into 2D fields (lat-lon)
! this is the inverse of ravel (except for the north pole)

    use params_mod,only:nmlat_h,nmlon
    use cons_mod,only:jlatm_JT

    real(kind=rp),dimension(nlonlat),intent(in) :: fin
    real(kind=rp),dimension(2,nmlat_h,nmlon) :: fout

    integer :: i,j,gidx

! start with southern hemisphere

! south pole has only one point (longitudes are squeezed)
    j = 1
    gidx = index_3to1(1,j,1)
    do i = 1,nmlon
      fout(1,j,i) = fin(gidx)
    enddo

! southern latitudes (no pole)
    do concurrent (i = 1:nmlon, j = 2:nmlat_h)
      gidx = index_3to1(1,j,i)
      fout(1,j,i) = fin(gidx)
    enddo

! continue with northern hemisphere

! northern low latitudes (duplicated from southern low latitudes)
    do concurrent (i = 1:nmlon, j = jlatm_JT:nmlat_h)
      gidx = index_3to1(1,j,i)
      fout(2,j,i) = fin(gidx)
    enddo

! northern high latitudes (no pole)
    do concurrent (i = 1:nmlon, j = 2:jlatm_JT-1)
      gidx = index_3to1(2,j,i)
      fout(2,j,i) = fin(gidx)
    enddo

! north pole is left unset (not in fin)

  endfunction unravel
!-----------------------------------------------------------------------
  pure subroutine csr_to_csc(nrow,ncol,nnz, &
    rowptr,colind,values_csr,colptr,rowind,values_csc)

    integer,intent(in) :: nrow,ncol,nnz
    integer,dimension(nrow+1),intent(in) :: rowptr
    integer,dimension(nnz),intent(in) :: colind
    real(kind=rp),dimension(nnz),intent(in) :: values_csr
    integer,dimension(ncol+1),intent(out) :: colptr
    integer,dimension(nnz),intent(out) :: rowind
    real(kind=rp),dimension(nnz),intent(out) :: values_csc

    integer :: n,i,j,newidx
    integer,dimension(ncol) :: colcnt,cnt

    colcnt = 0
    do n = 1,nnz
      j = colind(n)
      colcnt(j) = colcnt(j)+1
    enddo

    colptr(1) = 1
    do j = 2,ncol+1
      colptr(j) = colptr(j-1)+colcnt(j-1)
    enddo

    cnt = 0
    do i = 1,nrow
      do n = rowptr(i),rowptr(i+1)-1
        j = colind(n)
        newidx = colptr(j)+cnt(j)
        rowind(newidx) = i
        values_csc(newidx) = values_csr(n)
        cnt(j) = cnt(j)+1
      enddo
    enddo

  endsubroutine csr_to_csc
!-----------------------------------------------------------------------
#ifdef USE_MKL
  function solve_mkl(n,nnz,rowptr,colind,values,rhs) result(sol)

    include 'mkl_pardiso.fi'

    integer,intent(in) :: n,nnz
    integer,dimension(n+1),intent(in) :: rowptr
    integer,dimension(nnz),intent(in) :: colind
    real(kind=rp),dimension(nnz),intent(in) :: values
    real(kind=rp),dimension(n),intent(in) :: rhs
    real(kind=rp),dimension(n) :: sol

! for PARDISO sparse matrix solver
    integer,parameter :: &
      maxfct = 1, & ! Maximal number of factors in memory
      mnum = 1, &   ! The number of matrix (from 1 to maxfct) to solve
      mtype = 11, & ! Matrix type, real and non-symmetric matrix
      nrhs = 1, &   ! Number of right-hand sides that need to be solved for
      msglvl = 1    ! Message level information, prints statistical information
    integer :: &
      phase, & ! Controls the execution of the solver
      error    ! Error indicator
    integer,dimension(64) :: iparm
    integer,dimension(n) :: perm ! permutation vector
    double precision,dimension(nnz) :: a
    double precision,dimension(n) :: b,x
    type(MKL_PARDISO_HANDLE),dimension(64) :: pt

! initialize PARDISO with default parameters in accordance with the matrix type
    call pardisoinit(pt,mtype,iparm)

! setup some nonzero iparm elements in debug mode
#ifdef DEBUG

! 0: iparm(2) - iparm(64) are filled with default values
    iparm(1) = 1

! Report the number of non-zero elements in the factors
    iparm(18) = -1

! Report number of floating point operations (in 10^6 floating point operations)
! that are necessary to factor the matrix A
    iparm(19) = -1

! Matrix checker, 1 checks integer arrays rowptr and colind
    iparm(27) = 1
#endif

! fill permutation vector with zero (iparm(5)==0, permutation is ignored)
    perm = 0

    a = values
    b = rhs

! Analysis, numerical factorization, solve
    phase = 13

! use 32-bit integer version
! if the number of non-zero elements is on the order of 500 million or more
! then use pardiso_64 (64-bit integer version)
    call pardiso(pt, maxfct, mnum, mtype, phase, n, &
      a, rowptr, colind, perm, nrhs, iparm, msglvl, b, x, error)

    write(6,"('phase ',i4,' error ',i4)") phase,error

    sol = x

! Release all internal memory for all matrices
    phase = -1

    call pardiso(pt, maxfct, mnum, mtype, phase, n, &
      a, rowptr, colind, perm, nrhs, iparm, msglvl, b, x, error)

    write(6,"('phase ',i4,' error ',i4)") phase,error

  endfunction solve_mkl
#endif
!-----------------------------------------------------------------------
  function solve_superlu(n,nnz,colptr,rowind,values,rhs) result(sol)

    use iso_c_binding,only:c_int,c_long_long,c_double

    integer,intent(in) :: n,nnz
    integer,dimension(n+1),intent(in) :: colptr
    integer,dimension(nnz),intent(in) :: rowind
    real(kind=rp),dimension(nnz),intent(in) :: values
    real(kind=rp),dimension(n),intent(in) :: rhs
    real(kind=rp),dimension(n) :: sol

! for SuperLU sparse matrix solver
    integer(kind=c_int),parameter :: nrhs = 1
    integer(kind=c_int) :: iopt,n_cp,nnz_cp,ldb,info
    real(kind=c_double),dimension(nnz) :: values_cp
    integer(kind=c_int),dimension(nnz) :: rowind_cp
    integer(kind=c_int),dimension(n+1) :: colptr_cp
    real(kind=c_double),dimension(n) :: b
    integer(kind=c_long_long) :: f_factors

    interface
      subroutine c_fortran_dgssv(iopt,n,nnz,nrhs, &
        values,rowind,colptr,b,ldb,f_factors,info) &
        bind(c,name='c_fortran_dgssv_')
        use iso_c_binding,only:c_int,c_long_long,c_double
        integer(kind=c_int) :: iopt,n,nnz,nrhs,ldb,info
        real(kind=c_double),dimension(nnz) :: values
        integer(kind=c_int),dimension(nnz) :: rowind
        integer(kind=c_int),dimension(n+1) :: colptr
        real(kind=c_double),dimension(ldb) :: b
        integer(kind=c_long_long) :: f_factors
      endsubroutine c_fortran_dgssv
    endinterface

    n_cp = n
    nnz_cp = nnz
    ldb = n
    values_cp = values
    rowind_cp = rowind
    colptr_cp = colptr
    b = rhs

! first, factorize the matrix, the factors are stored in *f_factors* handle
    iopt = 1
    call c_fortran_dgssv(iopt, n_cp, nnz_cp, nrhs, &
      values_cp, rowind_cp, colptr_cp, b, ldb, f_factors, info)
    write(6,"('INFO from LU decomposition = ',i4)") info

! second, solve the system using the existing factors
    iopt = 2
    call c_fortran_dgssv(iopt, n_cp, nnz_cp, nrhs, &
      values_cp, rowind_cp, colptr_cp, b, ldb, f_factors, info)
    write(6,"('INFO from triangular solve = ',i4)") info

! last, free the storage allocated inside SuperLU
    iopt = 3
    call c_fortran_dgssv(iopt, n_cp, nnz_cp, nrhs, &
      values_cp, rowind_cp, colptr_cp, b, ldb, f_factors, info)
    write(6,"('INFO from freeing storage = ',i4)") info

    sol = b

  endfunction solve_superlu
!-----------------------------------------------------------------------
endmodule solver_mod
