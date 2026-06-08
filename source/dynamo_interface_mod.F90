!%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>
!
! Interface module for the 3-dimensional dynamo
!
!%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>
module dynamo_interface_mod

  use prec, only: rp

  use mpi_module, only: mlat0, mlat1, mlon0, mlon1 ! local mag field dims
  use mpi_module, only: mlond0, mlond1, mlatd0, mlatd1 ! for ghost pnts
  use mpi_module, only: sync_mlat_5d, sync_mlon_5d

  use cons_module, only: read_fac
  use params_module, only: nmlat_h, nmlatS2_h, nmlon, nhgt_fix, nhgt_fix_r
  use params_module, only: reproducible

  use mpi_module, only: mpi_init => init, setup_topology
  use mpi_module, only: un_mpi_rank, mpi_rank

  use grid_module, only: generate_mag_grid

  use init_module, only: init_cons, init_fieldline
  use init_module, only: get_apex, calculate_m

  use alloc_module, only: alloc_fieldline, alloc_fieldline_lite, dealloc_fieldline

  use fieldline_module, only: F_p,F_s1,F_s2,F_r,M3_p,M1_s1,M2_s2,M3_r
  use fieldline_module, only: be3_s1,be3_s2,bmag_p,vmp_p,D1_s1,D1_s2,d1d1_s1,d1d2_s1,d1d2_s2
  use fieldline_module, only: d2_s1,d2_s2,d2d2_s1,d2d2_s2,d_s1,d_s2,e1_s1,e1_s2,e2_s1,e2_s2
  use fieldline_module, only: bmag_s1, bmag_s2, vmp_s1, vmp_s2
  use fieldline_module, only: &
    npts_p, npts_s1, npts_s2, npts_r, &
    jmax_p, jmax_s1, jmax_s2, jmax_r, &
    size_p, size_s1, size_s2, size_r
  use fieldline_module, only: &
       glat_p, glon_p, glat_s1, glon_s1, glat_s2, glon_s2, glat_r, glon_r
  use fieldline_module, only: &
       qdlat_p, qdlat_s1, qdlat_s2, qdlat_r

  use calculate_terms_module, only: calculate_n, calculate_je, calculate_s
  use calculate_terms_module, only: calculate_ed, calculate_ve, calculate_vxyz
  use calculate_terms_module, only: calculate_conductance !, balance_fac_hl

  use stencil_module, only: calculate_coef, calculate_coef_ns2, calculate_coef_ns
  use stencil_module, only: calculate_bij

  use dist_solver_module, only: dist_linear_system, dist_solver_init

  implicit none

  private

  public :: dynamo_init1
  public :: dynamo_init2
  public :: dynamo_calc
  public :: dynamo_final

  public :: &
       npts_p, npts_s1, npts_s2, npts_r, &
       jmax_p, jmax_s1, jmax_s2, jmax_r, &
       size_p, size_s1, size_s2, size_r

  public :: &
       glat_p, glon_p, glat_s1, glon_s1, glat_s2, glon_s2, glat_r, glon_r

  public :: &
       qdlat_p, qdlat_s1, qdlat_s2, qdlat_r

  logical, parameter :: setbij = .true.

contains

  !------------------------------------------------------------------------------
  ! phase 1 of initialization -- called before APEX is initialized
  !------------------------------------------------------------------------------
  subroutine dynamo_init1( set_hilat_pot_in, set_hilat_fac_in, &
       mpicom_atm, npes_edyn3D, edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt, &
       edyn3d_slu_refactor_int,  edyn3d_slu_refactor_berr, reprod_solution, real_kind )

    logical, intent(in) :: set_hilat_pot_in, set_hilat_fac_in
    integer, intent(in) :: mpicom_atm, npes_edyn3D
    integer, intent(in) :: edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt
    integer, intent(in) :: edyn3d_slu_refactor_int
    real(kind=rp), intent(in) :: edyn3d_slu_refactor_BERR
    logical, optional, intent(in) :: reprod_solution
    integer, optional, intent(in) :: real_kind

    integer :: ierror
    character(len=*), parameter :: prefix = 'dynamo_init1: '

    if (present(real_kind)) then
       if (real_kind /= rp) then
          write(*,*) prefix, 'ERROR real_kind /= rp ...'
          write(*,*) ' rp: ',rp,' real_kind: ',real_kind
          stop 1
       end if
    end if

    read_fac = set_hilat_fac_in

    ! init mpi for 3D edynamo
    call mpi_init( mpicom_atm, npes_edyn3D )

    !init solver items
    call dist_solver_init(edyn3d_slu_refactor_int, edyn3d_slu_refactor_berr)

    ! set up magnetic latitude and longitude grids
    call generate_mag_grid(edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt)

    ! set up constants
    call init_cons()

    call alloc_fieldline_lite(ierror)
    if (ierror/=0) then
       ! call endrun(prefix//'alloc_fieldline_lite failed')
       write(*,*) prefix, 'alloc_fieldline_lite failed'
       stop ierror
    end if

    ! set up field-line grids
    call init_fieldline(npts_p,npts_s1,npts_s2,npts_r, &
         jmax_p,jmax_s1,jmax_s2,jmax_r, &
         size_p,size_s1,size_s2,size_r, &
         qdlat_p,qdlat_s1,qdlat_s2,qdlat_r)

    ! set up MPI decomposition
    call setup_topology(nmlat_h,nmlon)

    ! reproducible solution
    if (present(reprod_solution)) then
       reproducible = reprod_solution
    end if

  end subroutine dynamo_init1

  !------------------------------------------------------------------------------
  ! phase 2 of initialization -- called after APEX is initialized
  !------------------------------------------------------------------------------
  subroutine dynamo_init2()

    integer :: ierror
    character(len=*), parameter :: prefix = 'dynamo_init2: '

    active_tasks: if (mpi_rank >= 0) then

       ! allocate memory for fieldline data
       call alloc_fieldline(ierror)
       if (ierror/=0) then
          write(*,*) prefix, 'alloc_fieldline failed'
          stop ierror
       end if

       glat_p  = -huge(1._rp)
       glon_p  = -huge(1._rp)
       glat_s1 = -huge(1._rp)
       glon_s1 = -huge(1._rp)
       glat_s2 = -huge(1._rp)
       glon_s2 = -huge(1._rp)

       ! get apex coordinates and unit vectors
       call get_apex()

       ! calculate M coefficients - P,S1,S2,R
       call calculate_m(npts_p,npts_s1,npts_s2,npts_r, &
            F_p,F_s1,F_s2,F_r,M3_p,M1_s1,M2_s2,M3_r)

    end if active_tasks

  end subroutine dynamo_init2

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine dynamo_calc( &
       sigped_s1, sighal_s1, un_s1, vn_s1, &
       sigped_s2, sighal_s2, un_s2, vn_s2, &
       sigped_p, pot_hl_p, fac_hl_p, &
       ui_s1, vi_s1, wi_s1, ui_s2, vi_s2, wi_s2, &
       elec_pot_p, ped_cond_p, &
       efld1_s1, efld2_s1, efld1_s2, efld2_s2, &
       ionvel1_s1, ionvel2_s1, ionvel1_s2, ionvel2_s2, &
       bij_out )

    ! args
    real(rp), intent(in) :: sigped_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(rp), intent(in) :: sighal_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(rp), intent(in) :: un_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(rp), intent(in) :: vn_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(rp), intent(in) :: sigped_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(rp), intent(in) :: sighal_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(rp), intent(in) :: un_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(rp), intent(in) :: vn_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(rp), intent(in) :: sigped_p(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(rp), intent(in) :: pot_hl_p(2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp), intent(inout) :: fac_hl_p(2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp), intent(out) :: ui_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(rp), intent(out) :: vi_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(rp), intent(out) :: wi_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(rp), intent(out) :: ui_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(rp), intent(out) :: vi_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(rp), intent(out) :: wi_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(rp), optional, intent(out) :: elec_pot_p(2,mlat0:mlat1,mlon0:mlon1)
    real(rp), optional, intent(out) :: ped_cond_p(2,mlat0:mlat1,mlon0:mlon1)

    real(rp), optional, intent(out) :: efld1_s1(2,mlat0:mlat1,mlon0:mlon1)
    real(rp), optional, intent(out) :: efld2_s1(2,mlat0:mlat1,mlon0:mlon1)

    real(rp), optional, intent(out) :: efld1_s2(2,mlat0:mlat1,mlon0:mlon1)
    real(rp), optional, intent(out) :: efld2_s2(2,mlat0:mlat1,mlon0:mlon1)

    real(rp), optional, intent(out) :: ionvel1_s1(2,mlat0:mlat1,mlon0:mlon1)
    real(rp), optional, intent(out) :: ionvel2_s1(2,mlat0:mlat1,mlon0:mlon1)

    real(rp), optional, intent(out) :: ionvel1_s2(2,mlat0:mlat1,mlon0:mlon1)
    real(rp), optional, intent(out) :: ionvel2_s2(2,mlat0:mlat1,mlon0:mlon1)

    real(rp), optional, intent(out) :: bij_out(2,mlat0:mlat1,mlon0:mlon1)

    ! local vars

    real(rp) :: tmp_ghost(4,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: sigP_p(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: sigP_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: sigH_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: sigP_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: sigH_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: ntlU_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: ntlV_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: ntlU_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: ntlV_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: coef3d(9,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: coef2d(9,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: src3d(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: src2d(2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: N1p_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: N1h_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: N2p_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: N2h_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: Je1D_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: Je2D_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: S_p(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: coef(10,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: coef_ns2(10,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: coef_ns(10,2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: bij(mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: pot_p(2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: zigP_p(2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: zigP_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: zigH_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: zigP_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: zigH_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: ed1_s1(2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: ed2_s1(2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: ve1_s1(2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: ve2_s1(2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: ed1_s2(2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: ed2_s2(2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: ve1_s2(2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: ve2_s2(2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: vx_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: vy_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: vz_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: vx_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: vy_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp) :: vz_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)

    real(rp) :: fac_hl_loc(2,mlatd0:mlatd1,mlond0:mlond1)
    real(rp), parameter :: thres = 0.5_rp
    integer :: i, isn, j

    character(len=*), parameter :: subname = 'dynamo_calc'

    !only active tasks + extras for solve
    if (mpi_rank >=0) then

       ! exchange S1 ghost points
       tmp_ghost(1,:,:,mlat0:mlat1,mlon0:mlon1) = sigped_s1(:,:,mlat0:mlat1,mlon0:mlon1)
       tmp_ghost(2,:,:,mlat0:mlat1,mlon0:mlon1) = sighal_s1(:,:,mlat0:mlat1,mlon0:mlon1)
       tmp_ghost(3,:,:,mlat0:mlat1,mlon0:mlon1) = un_s1(:,:,mlat0:mlat1,mlon0:mlon1)
       tmp_ghost(4,:,:,mlat0:mlat1,mlon0:mlon1) = vn_s1(:,:,mlat0:mlat1,mlon0:mlon1)
       call sync_mlat_5d(tmp_ghost(:,:,:,:,mlon0:mlon1), 4, nhgt_fix, 2)
       call sync_mlon_5d(tmp_ghost, 4, nhgt_fix, 2)

       sigP_s1(:,:,:,:) = tmp_ghost(1,:,:,:,:)
       sigH_s1(:,:,:,:) = tmp_ghost(2,:,:,:,:)
       ntlU_s1(:,:,:,:) = tmp_ghost(3,:,:,:,:)
       ntlV_s1(:,:,:,:) = tmp_ghost(4,:,:,:,:)

       ! exchange S2 ghost points
       tmp_ghost(1,:,:,mlat0:mlat1,mlon0:mlon1) = sigped_s2(:,:,mlat0:mlat1,mlon0:mlon1)
       tmp_ghost(2,:,:,mlat0:mlat1,mlon0:mlon1) = sighal_s2(:,:,mlat0:mlat1,mlon0:mlon1)
       tmp_ghost(3,:,:,mlat0:mlat1,mlon0:mlon1) = un_s2(:,:,mlat0:mlat1,mlon0:mlon1)
       tmp_ghost(4,:,:,mlat0:mlat1,mlon0:mlon1) = vn_s2(:,:,mlat0:mlat1,mlon0:mlon1)
       call sync_mlat_5d(tmp_ghost(:,:,:,:,mlon0:mlon1), 4, nhgt_fix, 2)
       call sync_mlon_5d(tmp_ghost, 4, nhgt_fix, 2)
       sigP_s2(:,:,:,:) = tmp_ghost(1,:,:,:,:)
       sigH_s2(:,:,:,:) = tmp_ghost(2,:,:,:,:)
       ntlU_s2(:,:,:,:) = tmp_ghost(3,:,:,:,:)
       ntlV_s2(:,:,:,:) = tmp_ghost(4,:,:,:,:)

       ! calculate field-line integrated conductance - S1,S2
       call calculate_conductance( &
            mlatd0,mlatd1,mlond0,mlond1, &
            npts_s1,npts_s2, &
            vmp_s1,bmag_s1,sigP_s1,sigH_s1, &
            vmp_s2,bmag_s2,sigP_s2,sigH_s2, &
            zigP_s1,zigH_s1,zigP_s2,zigH_s2 )

       ! calculate N coefficients - S1,S2
       call calculate_n( &
            mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
            D_s1,M1_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1, &
            D_s2,M2_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2, &
            N1p_s1,N1h_s1,N2p_s2,N2h_s2)

       ! calculate JeD-coefficients (right hand side) - S1,S2
       call calculate_je( &
            mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
            D_s1,be3_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1,ntlU_s1,ntlV_s1, &
            D_s2,be3_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2,ntlU_s2,ntlV_s2, &
            d1_s1,d2_s1,d1_s2,d2_s2,Je1D_s1,Je2D_s2)

       ! calculate S (right hand side)
       S_p = calculate_s(mlatd0,mlatd1,mlond0,mlond1, &
            npts_p,M1_s1,Je1D_s1,M2_s2,Je2D_s2,M3_r)

       ! calculate height-dependent matrix coefficients
       coef = calculate_coef(mlatd0,mlatd1,mlond0,mlond1, &
            npts_p,S_p,N1p_s1,N1h_s1,N2p_s2,N2h_s2)

       ! add the coefficients in height to get coefficients for each hemisphere
       coef_ns2 = calculate_coef_ns2(mlatd0,mlatd1,mlond0,mlond1,coef)

       ! set the coefficient matrix in both hemispheres
       coef_ns = calculate_coef_ns(mlatd0,mlatd1,mlond0,mlond1,coef_ns2)

       ! set field-aligned conductance (b) matrix
       if (setbij) then
          bij = calculate_bij(mlatd0,mlatd1,mlond0,mlond1,coef_ns2)
       else
          bij = 0._rp
       endif

       fac_hl_p = 0._rp

       pot_p = -huge(1._rp)

    endif !just dynamp procs

    !now dynamo + extra needed for solve
    ! construct linear system and solve

    !If in the future you want to force a refactor (the numerics of the matrix
    !change), then call the following:
    !call dist_solver_force_refactor()

    call dist_linear_system(mlatd0,mlatd1,mlond0,mlond1, bij,pot_hl_p,fac_hl_p,coef_ns,pot_p)

    !just dynamo procs again
    if (mpi_rank >= 0) then

       !call edyn3d_hist_mlonlat_out('HILAT_FAC',fac_hl_p(1:2,mlat0:mlat1,mlon0:mlon1))
       !call edyn3d_hist_mlonlat_out('ELECPOTEN', pot_p(1:2,mlat0:mlat1,mlon0:mlon1))

       ! calculate electric fields
       call calculate_ed( &
            mlatd0,mlatd1,mlond0,mlond1, &
            pot_p,ed1_s1,ed2_s1,ed1_s2,ed2_s2)

       ! calculate drift velocities
       call calculate_ve( &
            mlatd0,mlatd1,mlond0,mlond1, &
            ed1_s1,ed2_s1,be3_s1(1,:,:,:), &
            ed1_s2,ed2_s2,be3_s2(1,:,:,:), &
            ve1_s1,ve2_s1,ve1_s2,ve2_s2)

       ! calculate drift velocities in geographic coordinates
       call calculate_vxyz( &
            mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
            ve1_s1,ve2_s1,e1_s1,e2_s1, &
            ve1_s2,ve2_s2,e1_s2,e2_s2, &
            vx_s1,vy_s1,vz_s1,vx_s2,vy_s2,vz_s2)
    end if

    ui_s1(1:nhgt_fix,1:2,mlat0:mlat1,mlon0:mlon1) = vx_s1(1:nhgt_fix,1:2,mlat0:mlat1,mlon0:mlon1)
    vi_s1(1:nhgt_fix,1:2,mlat0:mlat1,mlon0:mlon1) = vy_s1(1:nhgt_fix,1:2,mlat0:mlat1,mlon0:mlon1)
    wi_s1(1:nhgt_fix,1:2,mlat0:mlat1,mlon0:mlon1) = vz_s1(1:nhgt_fix,1:2,mlat0:mlat1,mlon0:mlon1)

    ui_s2(1:nhgt_fix,1:2,mlat0:mlat1,mlon0:mlon1) = vx_s2(1:nhgt_fix,1:2,mlat0:mlat1,mlon0:mlon1)
    vi_s2(1:nhgt_fix,1:2,mlat0:mlat1,mlon0:mlon1) = vy_s2(1:nhgt_fix,1:2,mlat0:mlat1,mlon0:mlon1)
    wi_s2(1:nhgt_fix,1:2,mlat0:mlat1,mlon0:mlon1) = vz_s2(1:nhgt_fix,1:2,mlat0:mlat1,mlon0:mlon1)

    if (present(elec_pot_p)) then
       elec_pot_p(1:2,mlat0:mlat1,mlon0:mlon1) = pot_p(1:2,mlat0:mlat1,mlon0:mlon1)
    end if

    if (present(efld1_s1)) then
       efld1_s1(1:2,mlat0:mlat1,mlon0:mlon1) = ed1_s1(1:2,mlat0:mlat1,mlon0:mlon1)
    end if
    if (present(efld2_s1)) then
       efld2_s1(1:2,mlat0:mlat1,mlon0:mlon1) = ed2_s1(1:2,mlat0:mlat1,mlon0:mlon1)
    end if
    if (present(efld1_s2)) then
       efld1_s2(1:2,mlat0:mlat1,mlon0:mlon1) = ed1_s2(1:2,mlat0:mlat1,mlon0:mlon1)
    end if
    if (present(efld2_s2)) then
       efld2_s2(1:2,mlat0:mlat1,mlon0:mlon1) = ed2_s2(1:2,mlat0:mlat1,mlon0:mlon1)
    end if

    if (present(ionvel1_s1)) then
       ionvel1_s1(1:2,mlat0:mlat1,mlon0:mlon1) = ve1_s1(1:2,mlat0:mlat1,mlon0:mlon1)
    end if
    if (present(ionvel2_s1)) then
       ionvel2_s1(1:2,mlat0:mlat1,mlon0:mlon1) = ve2_s1(1:2,mlat0:mlat1,mlon0:mlon1)
    end if
    if (present(ionvel1_s2)) then
       ionvel1_s2(1:2,mlat0:mlat1,mlon0:mlon1) = ve1_s2(1:2,mlat0:mlat1,mlon0:mlon1)
    end if
    if (present(ionvel2_s2)) then
       ionvel2_s2(1:2,mlat0:mlat1,mlon0:mlon1) = ve2_s2(1:2,mlat0:mlat1,mlon0:mlon1)
    end if

    if (present(bij_out)) then
       bij_out(1,mlat0:mlat1,mlon0:mlon1) = bij(mlat0:mlat1,mlon0:mlon1)
       bij_out(2,mlat0:mlat1,mlon0:mlon1) = bij(mlat0:mlat1,mlon0:mlon1)
    end if

  end subroutine dynamo_calc

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine dynamo_final()

    call dealloc_fieldline()

  end subroutine dynamo_final

end module dynamo_interface_mod
