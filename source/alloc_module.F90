module alloc_module
  use prec, only: rp
  use fieldline_module, only: npts_p, npts_s2, npts_s1, npts_r, qdlat_p, qdlat_s1, qdlat_s2, qdlat_r
  use fieldline_module, only: jmax_p, jmax_s1, jmax_s2, jmax_r, size_p, size_s1, size_s2, size_r
  use params_module, only: nmlat_h, nmlatS2_h, nhgt_fix, nhgt_fix_r

  implicit none

contains
!-----------------------------------------------------------------------
  subroutine alloc_fieldline_lite(ierr)
    integer, intent(out) :: ierr

    ierr = 0

    allocate( npts_p(nmlat_h), npts_s1(nmlat_h), npts_r(nmlat_h), stat=ierr)
    if (ierr /= 0) return
    allocate(npts_s2(nmlatS2_h), stat=ierr)
    if (ierr /= 0) return
    allocate(qdlat_p(nhgt_fix,2,nmlat_h),qdlat_s1(nhgt_fix,2,nmlat_h), stat=ierr)
    if (ierr /= 0) return
    allocate(qdlat_s2(nhgt_fix,2,nmlatS2_h), stat=ierr)
    if (ierr /= 0) return
    allocate(qdlat_r(nhgt_fix_r,2,nmlat_h), stat=ierr)
    if (ierr /= 0) return

    qdlat_p = -huge(1._rp)
    qdlat_s1 = -huge(1._rp)
    qdlat_r = -huge(1._rp)
    qdlat_s2 = -huge(1._rp)

    allocate(jmax_p(nhgt_fix),jmax_s1(nhgt_fix),jmax_s2(nhgt_fix),size_p(nhgt_fix),size_s1(nhgt_fix),size_s2(nhgt_fix), stat=ierr)
    if (ierr /= 0) return
    allocate(jmax_r(nhgt_fix_r),size_r(nhgt_fix_r), stat=ierr)
    if (ierr /= 0) return


  end subroutine alloc_fieldline_lite

!-----------------------------------------------------------------------
  subroutine alloc_fieldline(ierr)

! All grids include one latitude halo point on each side,
! but only internal grids (1<=j<=nmlat_h) are defined,
! the halo points at j==0 and j==nmlat_h+1 are not defined.

! Although S2 grids have equal latitudes with P/S1/R grids,
! the grid at j==nmlat_h is not defined.

    use mpi_module,only:mlond0,mlond1,mlatd0,mlatd1
    use cons_module, only: J3LB
    use fieldline_module

    integer, intent(out) :: ierr

    ierr = 0

    allocate(J3LB(2,nmlat_h,mlond0:mlond1), stat=ierr)
    if (ierr /= 0) return

    J3LB = 0._rp

    allocate(f3d(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1,32), stat=ierr)
    if (ierr /= 0) return
    allocate(f3d_r(nhgt_fix_r,2,mlatd0:mlatd1,mlond0:mlond1,6), stat=ierr)
    if (ierr /= 0) return
    allocate(uvec(3,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1,12), stat=ierr)
    if (ierr /= 0) return

    glat_p(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,1)
    glon_p(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,2)
    sinI_p(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,3)
    D_p(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,4)
    F_p(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,5)
    vmp_p(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,6)
    bmag_p(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,7)
    M3_p(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,8)

    glat_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,9)
    glon_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,10)
    sinI_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,11)
    D_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,12)
    F_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,13)
    vmp_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,14)
    bmag_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,15)
    be3_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,16)
    M1_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,17)
    d1d1_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,18)
    d1d2_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,19)
    d2d2_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,20)

! S2 grid only goes to nmlatS2_h
    glat_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,21)
    glon_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,22)
    sinI_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,23)
    D_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,24)
    F_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,25)
    vmp_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,26)
    bmag_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,27)
    be3_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,28)
    M2_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,29)
    d1d1_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,30)
    d1d2_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,31)
    d2d2_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,32)

    glat_r(1:,1:,mlatd0:,mlond0:) => f3d_r(:,:,:,:,1)
    glon_r(1:,1:,mlatd0:,mlond0:) => f3d_r(:,:,:,:,2)
    sinI_r(1:,1:,mlatd0:,mlond0:) => f3d_r(:,:,:,:,3)
    D_r(1:,1:,mlatd0:,mlond0:) => f3d_r(:,:,:,:,4)
    F_r(1:,1:,mlatd0:,mlond0:) => f3d_r(:,:,:,:,5)
    M3_r(1:,1:,mlatd0:,mlond0:) => f3d_r(:,:,:,:,6)

    d1_s1(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,1)
    d2_s1(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,2)
    d3_s1(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,3)
    e1_s1(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,4)
    e2_s1(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,5)
    e3_s1(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,6)

! S2 grid only goes to nmlatS2_h
    d1_s2(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,7)
    d2_s2(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,8)
    d3_s2(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,9)
    e1_s2(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,10)
    e2_s2(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,11)
    e3_s2(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,12)

  endsubroutine alloc_fieldline
!-----------------------------------------------------------------------

  subroutine dealloc_fieldline
    use fieldline_module

    glat_p=>null()
    glon_p=>null()
    sinI_p=>null()
    D_p=>null()
    F_p=>null()
    vmp_p=>null()
    bmag_p=>null()
    M3_p=>null()
    
    glat_s1=>null()
    glon_s1=>null()
    sinI_s1=>null()
    D_s1=>null()
    F_s1=>null()
    vmp_s1=>null()
    bmag_s1=>null()
    be3_s1=>null()
    M1_s1=>null()
    d1d1_s1=>null()
    d1d2_s1=>null()
    d2d2_s1=>null()
    
    glat_s2=>null()
    glon_s2=>null()
    sinI_s2=>null()
    D_s2=>null()
    F_s2=>null()
    vmp_s2=>null()
    bmag_s2=>null()
    be3_s2=>null()
    M2_s2=>null()
    d1d1_s2=>null()
    d1d2_s2=>null()
    d2d2_s2=>null()
    
    glat_r=>null()
    glon_r=>null()
    sinI_r=>null()
    D_r=>null()
    F_r=>null()
    M3_r=>null()
    
    d1_s1=>null()
    d2_s1=>null()
    d3_s1=>null()
    e1_s1=>null()
    e2_s1=>null()
    e3_s1=>null()
    d1_s2=>null()
    d2_s2=>null()
    d3_s2=>null()
    e1_s2=>null()
    e2_s2=>null()
    e3_s2=>null()
    
    if (allocated(f3d)) deallocate(f3d)
    if (allocated(f3d_r)) deallocate(f3d_r)
    if (allocated(uvec)) deallocate(uvec)

    if (allocated(npts_p)) deallocate(npts_p)
    if (allocated(npts_s1)) deallocate(npts_s1)
    if (allocated(npts_r)) deallocate(npts_r)
    if (allocated(npts_s2)) deallocate(npts_s2)
    if (allocated(qdlat_p)) deallocate(qdlat_p)
    if (allocated(qdlat_s1)) deallocate(qdlat_s1)
    if (allocated(qdlat_s2)) deallocate(qdlat_s2)
    if (allocated(qdlat_r)) deallocate(qdlat_r)
    if (allocated(jmax_p)) deallocate(jmax_p)
    if (allocated(jmax_s1)) deallocate(jmax_s1)
    if (allocated(jmax_s2)) deallocate(jmax_s2)
    if (allocated(size_p)) deallocate(size_p)
    if (allocated(size_s1)) deallocate(size_s1)
    if (allocated(size_s2)) deallocate(size_s2)
    if (allocated(jmax_r)) deallocate(jmax_r)
    if (allocated(size_r)) deallocate(size_r)
    
  end subroutine dealloc_fieldline
  
endmodule alloc_module
