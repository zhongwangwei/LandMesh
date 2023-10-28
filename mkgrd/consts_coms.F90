!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives

!----------------------------------------------------------------------------
! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University;
! Colorado State University Research Foundation ; ATMET, LLC

! This software is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option)
! any later version.

! This software is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.

! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA
! (http://www.gnu.org/licenses/gpl.html)
!----------------------------------------------------------------------------

!===============================================================================
Module consts_coms

    implicit none


    ! Single (real*4) and double (real*8) precision real kinds:
    integer, parameter :: r4 = selected_real_kind(6, 37)
    integer, parameter :: r8 = selected_real_kind(13, 300)

    ! 1 byte, 2 byte, 4 byte, and 8 byte integer (and logical) types
    integer, parameter :: i1 = selected_int_kind(2)
    integer, parameter :: i2 = selected_int_kind(4)
    integer, parameter :: i4 = selected_int_kind(8)
    integer, parameter :: i8 = selected_int_kind(16)

    real(r4), parameter :: spval_r4 = -1.e36_r4

    real, parameter :: pio180 = 3.1415927 / 180.   ! radians per degree
    real, parameter :: piu180 = 180. / 3.1415927   ! degrees per radian
    real, parameter :: pio2 = 3.1415927 * .5
    real, parameter :: pi2 = 3.1415927 * 2.

    real :: dlat                          ! Earth 1-latitude-degree dist [m]
    real :: omega                         ! Earth rotational angular velocity
    real :: omega2                        ! omega * 2
    real :: xscale                        ! erad scaling factor (for NCAR DCMIP test cases)
    real :: erad                          ! Earth radius [m]
    real :: erad2
    real :: erad4
    real :: eradsq
    real :: erador5
    real :: eradi
    real :: erad2sq

    integer, parameter :: maxremote = 30   ! Max # of remote send/recv processes
    integer, parameter :: pathlen = 256  ! Max length of character strings for file paths

    character(pathlen) :: tmpdir = "/tmp"
    character(64) :: expnme
    character(pathlen) :: base_dir
    character(pathlen) :: source_dir
    character(16) :: lcs
    real :: rinit = 0.0
    real(r8) :: rinit8 = 0.0_r8

    integer :: io6
    integer :: nxp
    integer :: mode
    integer :: openmp
    integer :: ndm_domain
    logical :: refine
    logical :: no_caculate_fraction
    integer :: nlons_source
    integer :: nlats_source
    integer :: nlons_dest
    integer :: nlats_dest
    integer :: maxlc
    real(r8) :: edgee(10)
    real(r8) :: edgew(10)
    real(r8) :: edges(10)
    real(r8) :: edgen(10)

    Type oname_vars
        !!    RUNTYPE/NAME
        character(64) :: expnme = 'OLAM test'
        character(pathlen) :: base_dir = ' /tmp'
        character(pathlen) :: source_dir = ' /tmp'
        character(16) :: lcs = 'igbp'
        !!    GRID SPECIFICATIONS
        integer :: nxp = 0
        integer :: openmp = 16
        integer :: mode = 6
        integer :: maxlc = 17

        logical :: refine = .FALSE.
        logical :: no_caculate_fraction   =    .FALSE.
        integer :: nlons_source           =    43200
        integer :: nlats_source           =    21600
        integer :: nlons_dest             =    43200
        integer :: nlats_dest             =    21600
        integer :: ndm_domain             =    1
        real(r8) :: edgee(10)             =    0.
        real(r8) :: edgew(10)             =    0.
        real(r8) :: edges(10)             =    0.
        real(r8) :: edgen(10)             =    0.

    End Type oname_vars
    type (oname_vars) :: nl
End Module
!===============================================================================

Module mem_grid
    use consts_coms, only : r8
    implicit none

    private :: r8

    integer :: &  ! N values for full-domain reference on any process

            nza           &  ! Vertical number of all points
            , nsw_max       &  ! Max # vert atm levels in IW column with sfc flux
            , nve2_max      &  ! Max # underground v[xyz]e2 levels in IW column
            , nma, nua, nva, nwa  ! Horiz number of all M,U,V,W points in full domain

    integer :: &

            mza           &  ! Vertical number of all points
            , mma, mua, mva, mwa  ! Horiz number of all M,U,V,W points in sub-process

    integer, allocatable, dimension(:) :: &

            lpm, lpv, lpw   &  ! Lowest prognosed M,V,W/T
            , lsw           &  ! number of W/T levels in contact with surface
            , lve2             ! number of v[xyz]e2 levels in contact with surface

    real, allocatable, dimension(:) :: &

            zm, zt        &  ! Z coordinate at M,T point
            , dzm, dzt       &  ! Delta Z at M,T point
            , dzim, dzit      &  ! Delta Z inverse (1/dz) at M,T point

            , zfacm, zfact     &  ! expansion factor of delta_x with height at M,T point
            , zfacim, zfacit    &  ! inverse of zfacm, zfact
            , zfacm2, zfacim2   &  ! expansion factor of arw with height, and its inverse

            , xem, yem, zem    &  ! XE,YE,ZE coordinates of M point
            , xev, yev, zev    &  ! XE,YE,ZE coordinates of V point
            , xew, yew, zew    &  ! XE,YE,ZE coordinates of W point

            , unx, uny, unz    &  ! U face normal unit vector components
            , vnx, vny, vnz    &  ! V face normal unit vector components
            , wnx, wny, wnz    &  ! W face normal unit vector components

            , dnu              &  ! dxy across U face
            , dniu             &  ! (1/dxy) across U face for turb U XY flx

            , dnv              &  ! dxy across V face for PGF
            , dniv             &  ! (1/dxy) across V face for turb V XY flx

            , arm0             &  ! Area of IM triangle at earth surface
            , arw0             &  ! Area of IW polygon at earth surface

            , glatm, glonm     &  ! Latitude/longitude at M point
            , glatv, glonv     &  ! Latitude/longitude at V point
            , glatw, glonw     &  ! Latitude/longitude at W point

            , topm, topw          ! Topography height at M,W point

    real, allocatable :: gravm(:) ! gravity at M levels
    real, allocatable :: gravt(:) ! gravity at T levels

    real, allocatable, dimension(:, :) :: &

            arv, arw            ! Aperture area of V,W face

    real(r8), allocatable, dimension(:, :) :: &

            volt                ! Volume of T cell

    real, allocatable, dimension(:, :) :: &
            volti, &          ! 1 / (Volume of T cell)
            volwi, &          ! 1 / (Volumes of adjacent T cells)
            volvi             ! 1 / (Volumes of adjacent T cells)

    integer :: impent(12)  ! Scratch array for storing 12 pentagonal IM indices

    integer, parameter :: nrows = 5
    integer :: mrows

    ! "Derived" variables computed at beginning of integration rather than
    ! at the MAKEGRID stage

    real, allocatable, dimension(:) :: &

            wnxo2, wnyo2, wnzo2, & ! W-face unit normals divided by 2
            vnxo2, vnyo2, vnzo2, & ! V-face unit normals divided by 2
            unxo2, unyo2, unzo2, & ! U-face unit normals divided by 2

            dzt_top, & ! distance between ZM(k) and ZT(k)
            dzt_bot, & ! distance between ZT(k) and ZM(k-1)

            dzit_top, & ! distance between ZM(k) and ZT(k)
            dzit_bot, & ! distance between ZT(k) and ZM(k-1)

            zwgt_top, zwgt_bot, & ! weights for interpolating T levels to W
            dzto2, dzto4, & ! dzt(k)    / 2, dzt(k)    / 4
            dztsqo2, dztsqo4, & ! dzt(k)**2 / 2, dzt(k)**2 / 4
            dztsqo6, dztsqo12, & ! dzt(k)**2 / 6, dzt(k)**2 / 12
            dzimsq, & ! dzim(k)**2

            voa0, & ! ratio of cell volume to bottom area w/o terrain

            gdz_belo, gdz_abov, & ! weights for hydrostatic integration

            gdzim, & ! gravm / dzm

            arw0i, & ! 1 / arw0
            dnivo2                  ! 1/(2dxy) across V face

    ! double precision weights for interpolating T levels to W

    real(r8), allocatable :: zwgt_top8(:), zwgt_bot8(:)
    real(r8), allocatable :: gdz_belo8(:), gdz_abov8(:)
    real(r8), allocatable :: gdz_wgtm8(:), gdz_wgtp8(:)
    real, allocatable :: gdz_wgtm (:), gdz_wgtp (:)

    real, allocatable, dimension(:, :) :: &

            gxps_coef, gyps_coef    ! combined weights for grad_t2d

    real, allocatable, dimension(:) :: &

            vxn_ew, vyn_ew, vzn_ew, & ! unit normals of zonal (east-west) direction in earth cartesian coordinates
            vxn_ns, vyn_ns, vzn_ns, & ! unit normals of meridional (north-south) direction in earth cartesian coordinates
            vcn_ew, vcn_ns            ! components of zonal and merdional vectors in the direction of VC

Contains

    !===============================================================================



    !===============================================================================

    subroutine alloc_xyzem(lma)

        implicit none

        integer, intent(in) :: lma

        allocate (xem(lma));  xem(1:lma) = 0.
        allocate (yem(lma));  yem(1:lma) = 0.
        allocate (zem(lma));  zem(1:lma) = 0.

    end subroutine alloc_xyzem

    !===============================================================================

    subroutine alloc_xyzew(lwa)

        implicit none

        integer, intent(in) :: lwa

        allocate (xew(lwa));  xew(1:lwa) = 0.
        allocate (yew(lwa));  yew(1:lwa) = 0.
        allocate (zew(lwa));  zew(1:lwa) = 0.

    end subroutine alloc_xyzew

    !===============================================================================

    subroutine alloc_grid1(lma, lva, lwa)

        use consts_coms, only : io6
        implicit none

        integer, intent(in) :: lma, lva, lwa

        ! Allocate and initialize arrays (xem, yem, zem are already allocated)

        allocate (lsw (lwa));  lsw (1:lwa) = 0
        allocate (lve2(lwa));  lve2(1:lwa) = 0

        allocate (xev(lva));  xev(1:lva) = 0.
        allocate (yev(lva));  yev(1:lva) = 0.
        allocate (zev(lva));  zev(1:lva) = 0.

        allocate (glatv(lva));  glatv(1:lva) = 0.
        allocate (glonv(lva));  glonv(1:lva) = 0.

        allocate (unx(lva));  unx(1:lva) = 0.
        allocate (uny(lva));  uny(1:lva) = 0.
        allocate (unz(lva));  unz(1:lva) = 0.

        allocate (vnx(lva));  vnx(1:lva) = 0.
        allocate (vny(lva));  vny(1:lva) = 0.
        allocate (vnz(lva));  vnz(1:lva) = 0.

        allocate (wnx(lwa));  wnx(1:lwa) = 0.
        allocate (wny(lwa));  wny(1:lwa) = 0.
        allocate (wnz(lwa));  wnz(1:lwa) = 0.

        allocate (dnu  (lva));  dnu (1:lva) = 0.
        allocate (dniu (lva));  dniu(1:lva) = 0.

        allocate (dnv  (lva));  dnv (1:lva) = 0.
        allocate (dniv (lva));  dniv(1:lva) = 0.

        allocate  (arw0(lwa));   arw0(1:lwa) = 0.
        allocate  (topw(lwa));   topw(1:lwa) = 0.
        allocate (glatw(lwa));  glatw(1:lwa) = 0.
        allocate (glonw(lwa));  glonw(1:lwa) = 0.

        allocate  (arm0(lma));   arm0(1:lma) = 0.
        allocate  (topm(lma));   topm(1:lma) = 0.
        allocate (glatm(lma));  glatm(1:lma) = 0.
        allocate (glonm(lma));  glonm(1:lma) = 0.

        allocate (vxn_ew(lwa)) ; vxn_ew = 0.
        allocate (vyn_ew(lwa)) ; vyn_ew = 0.
        allocate (vzn_ew(lwa)) ; vzn_ew = 0.

        allocate (vxn_ns(lwa)) ; vxn_ns = 0.
        allocate (vyn_ns(lwa)) ; vyn_ns = 0.
        allocate (vzn_ns(lwa)) ; vzn_ns = 0.

        allocate (vcn_ew(lva)) ; vcn_ew = 0.
        allocate (vcn_ns(lva)) ; vcn_ns = 0.

        write(io6, *) 'finishing alloc_grid1'

    end subroutine alloc_grid1

    !===============================================================================

    subroutine alloc_grid2(lma, lva, lwa)

        use consts_coms, only : r8
        use consts_coms, only : io6

        implicit none

        integer, intent(in) :: lma, lva, lwa

        ! Allocate  and initialize arrays

        write(io6, *) 'alloc_grid2 ', lma, lva, lwa

        allocate (lpv(lva)); lpv(1:lva) = 0
        allocate (lpm(lma)); lpm(1:lma) = 0
        allocate (lpw(lwa)); lpw(1:lwa) = 0

        allocate (arv  (mza, lva));  arv  (1:mza, 1:lva) = 0.
        allocate (arw  (mza, lwa));  arw  (1:mza, 1:lwa) = 0.
        allocate (volt (mza, lwa));  volt (1:mza, 1:lwa) = 0._r8

        write(io6, *) 'finishing alloc_grid2'

    end subroutine alloc_grid2


End Module mem_grid

Module mem_ijtabs

    use consts_coms, only : maxremote

    implicit none

    private :: maxremote

    integer, parameter :: mloops = 7 ! max # non-para DO loops for M,V,W pts

    integer, parameter :: nloops_m = mloops + maxremote ! # M DO loops incl para
    integer, parameter :: nloops_v = mloops + maxremote ! # V DO loops incl para
    integer, parameter :: nloops_w = mloops + maxremote ! # W DO loops incl para

    ! M, U, V, and W loop indices. The last 4 letters of these indices have the
    ! following meanings:
    !
    ! *_grid is for setting up grid parameters in the ctrlvols subroutines
    ! *_init is for initialization of ATM fields (in ohhi, olhi, fldsisan,
    !        hurricane_init, omic_init)
    ! *_prog is for points where primary quantities such as VMC or THIL are
    !        prognosed
    ! *_wadj is for U, V, or W points that are adjacent to W points where primary
    !        quantities are prognosed
    ! *_wstn is for all U, V, or W points in the stencil of a W point where
    !        primary quantities are prognosed
    ! *_lbcp is for lateral boundary points whose values are copied from an
    !        interior prognostic point
    ! *_vadj is for M points that are adjacent to a prognostic V point, while
    !        jtw_vadj is the representation of jtm_vadj in the pre-hex-grid stage
    !        of model grid initialization)

    integer, parameter :: jtm_grid = 1, jtu_grid = 1, jtv_grid = 1, jtw_grid = 1
    integer, parameter :: jtm_init = 2, jtu_init = 2, jtv_init = 2, jtw_init = 2
    integer, parameter :: jtm_prog = 3, jtu_prog = 3, jtv_prog = 3, jtw_prog = 3
    integer, parameter :: jtm_wadj = 4, jtu_wadj = 4, jtv_wadj = 4, jtw_wadj = 4
    integer, parameter :: jtm_wstn = 5, jtu_wstn = 5, jtv_wstn = 5, jtw_wstn = 5
    integer, parameter :: jtm_lbcp = 6, jtu_lbcp = 6, jtv_lbcp = 6, jtw_lbcp = 6
    integer, parameter :: jtm_vadj = 7, jtu_wall = 7, jtv_wall = 7, jtw_vadj = 7

    integer :: nstp  ! # of finest grid acoustic timesteps in coarse grid dtlong
    integer :: istp  ! Current timestep counter from 1 to nstp
    integer :: mrls  ! Number of active mesh refinement levels (MRLs)

    integer, allocatable :: mrl_begl(:)  ! MRL at beginning of long timestep
    integer, allocatable :: mrl_begr(:)  ! MRL at beginning of RK step
    integer, allocatable :: mrl_begs(:)  ! MRL at beginning of short timestep
    integer, allocatable :: mrl_ends(:)  ! MRL at end of short timestep
    integer, allocatable :: mrl_endr(:)  ! MRL at end of RK step
    integer, allocatable :: mrl_endl(:)  ! MRL at end of long timestep
    real, allocatable :: dtrk    (:)  ! MRL RK timestep factor

    integer, allocatable :: leafstep(:)  ! flag to run leaf on any sub-timestep

    Type itab_m_vars             ! data structure for M pts (individual rank)
        logical, allocatable :: loop(:) ! flag to perform each DO loop at this M pt

        integer :: npoly = 0       ! number of V/W neighbors of this M pt
        integer :: imp = 1         ! M point from which to copy this M pt's values
        integer :: imglobe = 1     ! global index of this M pt (in parallel case)
        integer :: mrlm = 0        ! mesh refinement level of this M pt
        integer :: mrlm_orig = 0   ! original MRL of this M pt (hex only)
        integer :: mrow = 0        ! Full row number outside nest
        integer :: ngr = 0         ! Grid number
        integer :: iv(3) = 1       ! array of V neighbors of this M pt
        integer :: iw(3) = 1       ! array of W neighbors of this M pt
    End Type itab_m_vars

    Type itab_v_vars             ! data structure for V pts (individual rank)
        logical, allocatable :: loop(:) ! flag to perform each DO loop at this V pt

        integer :: ivp = 1       ! V pt from which to copy this V pt's values
        integer :: irank = -1    ! rank of parallel process at this V pt
        integer :: ivglobe = 1   ! global index of this V pt (in parallel case)
        integer :: mrlv = 0      ! mesh refinement level of this V pt
        integer :: im(6) = 1     ! neighbor M pts of this V pt
        integer :: iw(4) = 1     ! neighbor W pts of this V pt
        integer :: iv(4) = 1     ! neighbor V pts

        real :: farw(2) = 0.     ! Interp of ARW to V control volume [VADV + VDIFF]

        real :: cosv(2) = 0.     ! cosine of angle between V and zonal dir (Voronoi)
        real :: sinv(2) = 0.     ! sine of angle between V and zonal dir (Voronoi)

        real :: dxps(2) = 0.     ! xps (eastward) displacement from neighbor W pts
        real :: dyps(2) = 0.     ! yps (northward) displacement from neighbor W pts
    End Type itab_v_vars

    Type itab_w_vars             ! data structure for W pts (individual rank)
        logical, allocatable :: loop(:) ! flag to perform each DO loop at this W pt

        integer :: npoly = 0     ! number of M/V neighbors of this W pt
        integer :: iwp = 1       ! W pt from which to copy this W pt's values
        integer :: irank = -1    ! rank of parallel process at this W pt
        integer :: iwglobe = 1   ! global index of this W pt (in parallel run)
        integer :: mrlw = 0      ! mesh refinement level of this W pt
        integer :: mrlw_orig = 0 ! original MRL of this W pt
        integer :: ngr = 0       ! Grid number
        integer :: im(7) = 1     ! neighbor M pts
        integer :: iv(7) = 1     ! neighbor V pts
        integer :: iw(7) = 1     ! neighbor W pts

        real :: dirv(7) = 0.     ! pos direction of V neighbors

        real :: farm(7) = 0.     ! Fraction of arw0 in each M point sector
        real :: farv(7) = 0.     ! Fraction of arw0 in each V point sector

        real :: gxps1(7) = 0.    ! gradient weight xe component for point 1
        real :: gyps1(7) = 0.    ! gradient weight ye component for point 1

        real :: gxps2(7) = 0.    ! gradient weight xe component for point 2
        real :: gyps2(7) = 0.    ! gradient weight ye component for point 2

        real :: unx_w = 0.       ! xe component of eastward unit normal vector
        real :: uny_w = 0.       ! ye component of eastward unit normal vector

        real :: vnx_w = 0.       ! xe component of northward unit normal vector
        real :: vny_w = 0.       ! ye component of northward unit normal vector
        real :: vnz_w = 0.       ! ze component of northward unit normal vector

        real :: ecvec_vx(7) = 0. ! factors converting V to earth cart. velocity
        real :: ecvec_vy(7) = 0. ! factors converting V to earth cart. velocity
        real :: ecvec_vz(7) = 0. ! factors converting V to earth cart. velocity

        integer :: iwnud(3) = 1  ! local nudpoly pts
        real :: fnud(3) = 0.  ! local nudpoly coeffs

        !------------------------------------------------------------------------------

        integer :: jsfc2 = 0 ! number of surface cells attached to this W column
        integer :: jland1 = 0 ! beginning land cell counter (if any land cells present)
        integer :: jland2 = 0 ! ending    land cell counter (if any land cells present)
        integer :: jlake1 = 0 ! beginning lake cell counter (if any lake cells present)
        integer :: jlake2 = 0 ! ending    lake cell counter (if any lake cells present)
        integer :: jsea1 = 0 ! beginning sea  cell counter (if any sea  cells present)
        integer :: jsea2 = 0 ! ending    sea  cell counter (if any sea  cells present)

        integer, allocatable :: iwsfc (:) ! local-rank indices of attached surface cells
        integer, allocatable :: jasfc (:) ! atm j index of attached surface cells

    End Type itab_w_vars

    Type itabg_m_vars            ! data structure for M pts (global)
        integer :: im_myrank = -1 ! local (parallel subdomain) index of this M pt
        integer :: irank = -1     ! rank of parallel process at this M pt
        integer :: im_myrank_imp = -1 ! local M point that corresponds to global imp
    End Type itabg_m_vars

    Type itabg_v_vars            ! data structure for V pts (global)
        integer :: iv_myrank = -1 ! local (parallel subdomain) index of this V pt
        integer :: irank = -1     ! rank of parallel process at this V pt
        integer :: iv_myrank_ivp = -1 ! local V point that corresponds to global ivp
    End Type itabg_v_vars

    Type itabg_w_vars            ! data structure for W pts (global)
        integer :: iw_myrank = -1 ! local (parallel subdomain) index of this W pt
        integer :: irank = -1     ! rank of parallel process at this W pt
        integer :: iw_myrank_iwp = -1 ! local W point that corresponds to global iwp
    End Type itabg_w_vars

    Type jtab_m_vars
        integer, allocatable :: im(:)
        integer, allocatable :: jend(:)
    End Type jtab_m_vars

    Type jtab_v_vars
        integer, allocatable :: iv(:)
        integer, allocatable :: jend(:)
    End Type jtab_v_vars

    Type jtab_w_vars
        integer, allocatable :: iw(:)
        integer, allocatable :: jend(:)
    End Type jtab_w_vars

    type (itab_m_vars), allocatable, target :: itab_m(:)
    type (itab_v_vars), allocatable, target :: itab_v(:)
    type (itab_w_vars), allocatable, target :: itab_w(:)

    type (itabg_m_vars), allocatable, target :: itabg_m(:)
    type (itabg_v_vars), allocatable, target :: itabg_v(:)
    type (itabg_w_vars), allocatable, target :: itabg_w(:)

    type (jtab_m_vars) :: jtab_m(nloops_m)
    type (jtab_v_vars) :: jtab_v(nloops_v)
    type (jtab_w_vars) :: jtab_w(nloops_w)

Contains

    !===============================================================================

    subroutine alloc_itabs(mma, mva, mwa, input)

        implicit none

        integer, intent(in) :: mma, mva, mwa, input
        integer :: im, iv, iw

        allocate (itab_m(mma))
        allocate (itab_v(mva))
        allocate (itab_w(mwa))

        do im = 1, mma
            allocate(itab_m(im)%loop(mloops))
            itab_m(im)%loop(1:mloops) = .false.
        enddo

        do iv = 1, mva
            allocate(itab_v(iv)%loop(mloops))
            itab_v(iv)%loop(1:mloops) = .false.
        enddo

        do iw = 1, mwa
            allocate(itab_w(iw)%loop(mloops))
            itab_w(iw)%loop(1:mloops) = .false.
        enddo

    end subroutine alloc_itabs


    !===============================================================================

    subroutine fill_jtabs(mma, mva, mwa, input)

        implicit none

        integer, intent(in) :: mma, mva, mwa, input

        integer :: iw, iv, im
        integer :: iloop, jend
        integer :: nlm, nlv, nlw

        nlm = mloops
        nlv = mloops
        nlw = mloops
        ! Allocate and zero-fill jtab%jend()

        do iloop = 1, nlm
            allocate (jtab_m(iloop)%jend(mrls))
            jtab_m(iloop)%jend(1:mrls) = 0
        enddo

        if (allocated(itab_v)) then
            do iloop = 1, nlv
                allocate (jtab_v(iloop)%jend(mrls))
                jtab_v(iloop)%jend(1:mrls) = 0
            enddo
        endif

        do iloop = 1, nlw
            allocate (jtab_w(iloop)%jend(mrls))
            jtab_w(iloop)%jend(1:mrls) = 0
        enddo


        ! Compute and store jtab%jend(1)
        do iloop = 1, nlm
            jtab_m(iloop)%jend(1) = 0
            do im = 2, mma
                if (itab_m(im)%loop(iloop)) then
                    jtab_m(iloop)%jend(1) = jtab_m(iloop)%jend(1) + 1
                endif
            enddo
            jtab_m(iloop)%jend(1) = max(1, jtab_m(iloop)%jend(1))
        enddo

        do iloop = 1, nlv
            jtab_v(iloop)%jend(1) = 0
            do iv = 2, mva
                if (itab_v(iv)%loop(iloop)) then
                    jtab_v(iloop)%jend(1) = jtab_v(iloop)%jend(1) + 1
                endif
            enddo
            jtab_v(iloop)%jend(1) = max(1, jtab_v(iloop)%jend(1))
        enddo

        do iloop = 1, nlw
            jtab_w(iloop)%jend(1) = 0
            do iw = 2, mwa
                if (itab_w(iw)%loop(iloop)) then
                    jtab_w(iloop)%jend(1) = jtab_w(iloop)%jend(1) + 1
                endif
            enddo
            jtab_w(iloop)%jend(1) = max(1, jtab_w(iloop)%jend(1))
        enddo

        ! Allocate and zero-fill JTAB_M%IM, JTAB_V%IV, JTAB_W%IW
        do iloop = 1, nlm
            jend = jtab_m(iloop)%jend(1)
            allocate (jtab_m(iloop)%im(jend))
            jtab_m(iloop)%im(1:jend) = 0
        enddo

        do iloop = 1, nlv
            jend = jtab_v(iloop)%jend(1)
            allocate (jtab_v(iloop)%iv(jend))
            jtab_v(iloop)%iv(1:jend) = 0
        enddo

        do iloop = 1, nlw
            jend = jtab_w(iloop)%jend(1)
            allocate (jtab_w(iloop)%iw(jend))
            jtab_w(iloop)%iw(1:jend) = 0
        enddo
    end subroutine fill_jtabs

End Module mem_ijtabs

Module mem_delaunay

    use mem_ijtabs, only : mloops

    implicit none

    private :: mloops

    Type itab_md_vars             ! data structure for M pts (individual rank)
        ! on the Delaunay mesh
        logical :: loop(mloops) = .false. ! flag to perform each DO loop at this M pt
        integer :: npoly = 0       ! number of V/W neighbors of this M pt
        integer :: imp = 1         ! M point from which to copy this M pt's values
        integer :: mrlm = 0        ! mesh refinement level of this M pt
        integer :: mrlm_orig = 0   ! original MRL of this M pt (hex only)
        integer :: ngr = 0         ! Grid number
        integer :: im(7) = 1       ! array of M neighbors of this M pt
        integer :: iu(7) = 1       ! array of U neighbors of this M pt
        integer :: iw(7) = 1       ! array of W neighbors of this M pt
    End Type itab_md_vars

    Type itab_ud_vars             ! data structure for U pts (individual rank)
        ! on the Delaunay mesh
        logical :: loop(mloops) = .false. ! flag to perform each DO loop at this M pt
        integer :: iup = 1       ! U pt from which to copy this U pt's values
        integer :: mrlu = 0      ! mesh refinement level of this U pt
        integer :: im(2) = 1     ! neighbor M pts of this U pt
        integer :: iu(12) = 1    ! neighbor U pts
        integer :: iw(6) = 1     ! neighbor W pts
    End Type itab_ud_vars

    Type itab_wd_vars             ! data structure for W pts (individual rank)
        ! on the Delaunay mesh
        logical :: loop(mloops) = .false. ! flag to perform each DO loop at this W pt
        integer :: npoly = 0     ! number of M/V neighbors of this W pt
        integer :: iwp = 1       ! W pt from which to copy this W pt's values
        integer :: mrlw = 0      ! mesh refinement level of this W pt
        integer :: mrlw_orig = 0 ! original MRL of this W pt
        integer :: mrow = 0      ! Full row number outside nest
        integer :: ngr = 0       ! Grid number
        integer :: im(3) = 1     ! neighbor M pts
        integer :: iu(3) = 1     ! neighbor U pts
        integer :: iw(9) = 1     ! neighbor W pts
    End Type itab_wd_vars

    Type nest_ud_vars        ! temporary U-pt data structure for spawning nested grids
        integer :: im = 0, iu = 0 ! new M/U pts attached to this U pt
    End Type nest_ud_vars

    Type nest_wd_vars        ! temporary W-pt data structure for spawning nested grids
        integer :: iu(3) = 0  ! new U pts attached to this W pt
        integer :: iw(3) = 0  ! new W pts attached to this W pt
    End Type nest_wd_vars

    type (itab_md_vars), allocatable :: itab_md(:)
    type (itab_ud_vars), allocatable :: itab_ud(:)
    type (itab_wd_vars), allocatable :: itab_wd(:)

    type (itab_md_vars), allocatable :: itab_md_copy(:)
    type (itab_ud_vars), allocatable :: itab_ud_copy(:)
    type (itab_wd_vars), allocatable :: itab_wd_copy(:)

    real, allocatable :: xemd(:), yemd(:), zemd(:)

    real, allocatable :: xemd_copy(:), yemd_copy(:), zemd_copy(:)

    integer :: nmd, nud, nwd

    integer :: nmd_copy = 0
    integer :: nud_copy = 0
    integer :: nwd_copy = 0

    integer, allocatable :: iwdorig(:), iwdorig_temp(:)

Contains

    !===============================================================================

    subroutine alloc_itabsd(mma, mua, mwa)

        use consts_coms, only : rinit
        use mem_ijtabs, only : mloops

        implicit none

        integer, intent(in) :: mma, mua, mwa

        allocate (itab_md(mma))
        allocate (itab_ud(mua))
        allocate (itab_wd(mwa))

        allocate(xemd(mma)) ; xemd = rinit
        allocate(yemd(mma)) ; yemd = rinit
        allocate(zemd(mma)) ; zemd = rinit

        xemd(1) = 0.
        yemd(1) = 0.
        zemd(1) = 0.

    end subroutine alloc_itabsd

    !===============================================================================

    subroutine copy_tri_grid()

        implicit none

        ! Save a copy of triangle structure of ATM grid in its current state of
        ! construction for subsequent independent local refinement of SURFACE grid.

        nmd_copy = nmd
        nud_copy = nud
        nwd_copy = nwd

        allocate (xemd_copy(nmd))
        allocate (yemd_copy(nmd))
        allocate (zemd_copy(nmd))

        allocate (itab_md_copy(nmd))
        allocate (itab_ud_copy(nud))
        allocate (itab_wd_copy(nwd))

        xemd_copy = xemd
        yemd_copy = yemd
        zemd_copy = zemd

        itab_md_copy = itab_md
        itab_ud_copy = itab_ud
        itab_wd_copy = itab_wd

    end subroutine copy_tri_grid

    !===============================================================================

    subroutine copyback_tri_grid()

        implicit none

        integer :: iw

        ! Save a copy of triangle structure of ATM grid in its current state of
        ! construction for subsequent independent local refinement of SURFACE grid.

        if (nmd_copy == 0 .or. nud_copy == 0 .or. nwd_copy == 0) then
            write(*, *) " Error in copyback_tri_grid:"
            stop ' Delaunay mesh was not previously saved'
        endif

        nmd = nmd_copy
        nud = nud_copy
        nwd = nwd_copy

        call move_alloc(xemd_copy, xemd)
        call move_alloc(yemd_copy, yemd)
        call move_alloc(zemd_copy, zemd)

        call move_alloc(itab_md_copy, itab_md)
        call move_alloc(itab_ud_copy, itab_ud)
        call move_alloc(itab_wd_copy, itab_wd)

        allocate(iwdorig(nwd))
        do iw = 1, nwd
            iwdorig(iw) = iw
        enddo

    end subroutine copyback_tri_grid

    !===============================================================================

End Module mem_delaunay

Module refine_vars

    use consts_coms, only : r8

    implicit none

    integer :: ndm_refine              =  3

    real(r8) :: edgee_rf(10)           =  0.
    real(r8) :: edgew_rf(10)           =  0.
    real(r8) :: edges_rf(10)           =  0.
    real(r8) :: edgen_rf(10)           =  0.

    integer :: max_iter               =  4
    integer :: max_sa_iter            =  500

    integer :: th_num_landtypes       =  12
    real(r8) :: th_area_mainland      =  0.6
    real(r8) :: th_lai_m              =  52.0 !LAI*10
    real(r8) :: th_lai_s              =  16.0 !LAI*10
    real(r8) :: th_slope_m            =  7.0
    real(r8) :: th_slope_s            =  5.0
    real(r8) :: th_k_s_m              =  300.0
    real(r8) :: th_k_s_s              =  160.0
    real(r8) :: th_k_solids_m         =  1.8
    real(r8) :: th_k_solids_s         =  0.4
    real(r8) :: th_tkdry_m            =  0.3
    real(r8) :: th_tkdry_s            =  0.3
    real(r8) :: th_tksatf_m           =  3.3
    real(r8) :: th_tksatf_s           =  3.3
    real(r8) :: th_tksatu_m           =  1.8
    real(r8) :: th_tksatu_s           =  1.8

    logical  :: refine_num_landtypes   =  .TRUE.
    logical  :: refine_area_mainland   =  .TRUE.
    logical  :: refine_lai_m           =  .TRUE.
    logical  :: refine_lai_s           =  .TRUE.
    logical  :: refine_slope_m         =  .TRUE.
    logical  :: refine_slope_s         =  .TRUE.
    logical  :: refine_k_s_m           =  .TRUE.
    logical  :: refine_k_s_s           =  .TRUE.
    logical  :: refine_k_solids_m      =  .TRUE.
    logical  :: refine_k_solids_s      =  .TRUE.
    logical  :: refine_tkdry_m         =  .TRUE.
    logical  :: refine_tkdry_s         =  .TRUE.
    logical  :: refine_tksatf_m        =  .TRUE.
    logical  :: refine_tksatf_s        =  .TRUE.
    logical  :: refine_tksatu_m        =  .TRUE.
    logical  :: refine_tksatu_s        =  .TRUE.


    Type threshold_vars
        integer :: ndm_refine
        
        real(r8) :: edgee_rf(10)
        real(r8) :: edgew_rf(10)
        real(r8) :: edges_rf(10)
        real(r8) :: edgen_rf(10)

        integer :: max_iter
        integer :: max_sa_iter
        integer :: th_num_landtypes

        logical  :: refine_num_landtypes  =  .TRUE.
        logical  :: refine_area_mainland  =  .TRUE.
        logical  :: refine_lai_m          =  .TRUE.
        logical  :: refine_lai_s          =  .TRUE.
        logical  :: refine_slope_m        =  .TRUE.
        logical  :: refine_slope_s        =  .TRUE.
        logical  :: refine_k_s_m          =  .TRUE.
        logical  :: refine_k_s_s          =  .TRUE.
        logical  :: refine_k_solids_m     =  .TRUE.
        logical  :: refine_k_solids_s     =  .TRUE.
        logical  :: refine_tkdry_m        =  .TRUE.
        logical  :: refine_tkdry_s        =  .TRUE.
        logical  :: refine_tksatf_m       =  .TRUE.
        logical  :: refine_tksatf_s       =  .TRUE.
        logical  :: refine_tksatu_m       =  .TRUE.
        logical  :: refine_tksatu_s       =  .TRUE.

        real(r8) :: th_area_mainland
        real(r8) :: th_lai_m
        real(r8) :: th_lai_s
        real(r8) :: th_slope_m
        real(r8) :: th_slope_s
        real(r8) :: th_k_s_m
        real(r8) :: th_k_s_s
        real(r8) :: th_k_solids_m
        real(r8) :: th_k_solids_s
        real(r8) :: th_tkdry_m
        real(r8) :: th_tkdry_s
        real(r8) :: th_tksatf_m
        real(r8) :: th_tksatf_s
        real(r8) :: th_tksatu_m
        real(r8) :: th_tksatu_s

    End Type threshold_vars
    type (threshold_vars) :: rl


End Module refine_vars
