!DESCRIPTION
!===========
! Voronoi geometry module
! This module contains subroutines for Voronoi diagram computation
! and geometric calculations
!
!REVISION HISTORY
!----------------
! 2025.06.11  Zhongwang Wei @ SYSU (revised version)
! 2025.06.10  Rui Zhang @ SYSU (original version)

module MOD_voronoi_geometry
    use consts_coms, only : io6, r8
    implicit none
    
    public :: voronoi, pcvt, grid_geometry_hex

contains

!-- Computes the Voronoi diagram from an existing Delaunay triangulation.
SUBROUTINE voronoi()

    use mem_ijtabs, only : mloops, itab_m, itab_v, itab_w, alloc_itabs

    use mem_delaunay, only : itab_md, itab_ud, itab_wd, &
            xemd, yemd, zemd, nmd, nud, nwd

    use mem_grid, only : nma, nua, nva, nwa, mma, mua, mva, mwa, &
            xem, yem, zem, xew, yew, zew, &
            alloc_xyzem, alloc_xyzew

    use consts_coms, only : pi2, erad

    implicit none

    integer :: im1, im2
    integer :: iw1, iw2, iw3, im, iw
    integer :: iwd, iv, iud, iud1, iud2, imd, npoly, j, j1
    real :: expansion

    ! Interchange grid dimensions

    nma = nwd
    nua = nud
    nva = nud
    nwa = nmd

    mma = nma
    mua = nua
    mva = nva
    mwa = nwa

    ! Allocate Voronoi set of itabs

    call alloc_itabs(nma, nva, nwa, 0)

    ! Allocate XEW,YEW,ZEW arrays, and fill their values from XEMD,YEMD,ZEMD, which
    ! still have the OLD nmad dimension which is the NEW nwa dimension

    call move_alloc(xemd, xew)
    call move_alloc(yemd, yew)
    call move_alloc(zemd, zew)

    ! Allocate XEM,YEM,ZEM to NEW nma dimension

    call alloc_xyzem(nma)

    ! Since XEM,YEM,ZEM have just been re-allocated, initialize their values to be
    ! barycenters of Delaunay triangles whose vertices are at XEW,YEW,ZEW

    do iwd = 2, nwd
        im = iwd

        ! Indices of 3 M points surrounding WD point

        if (any(itab_wd(iwd)%im(1:3) < 2)) cycle

        iw1 = itab_wd(iwd)%im(1)
        iw2 = itab_wd(iwd)%im(2)
        iw3 = itab_wd(iwd)%im(3)

        xem(im) = (xew(iw1) + xew(iw2) + xew(iw3)) / 3.
        yem(im) = (yew(iw1) + yew(iw2) + yew(iw3)) / 3.
        zem(im) = (zew(iw1) + zew(iw2) + zew(iw3)) / 3.

        ! push M point coordinates out to earth radius
        expansion = erad / sqrt(xem(im) ** 2  &
                + yem(im) ** 2  &
                + zem(im) ** 2)

        xem(im) = xem(im) * expansion
        yem(im) = yem(im) * expansion
        zem(im) = zem(im) * expansion

    enddo

    ! Loop over V points

    do iv = 2, nva
        iud = iv

        itab_v(iv)%loop(1:mloops) = itab_ud(iud)%loop(1:mloops)

        itab_v(iv)%ivp = itab_ud(iud)%iup
        itab_v(iv)%ivglobe = iv
        itab_v(iv)%mrlv = itab_ud(iud)%mrlu

        itab_v(iv)%im(1:6) = itab_ud(iud)%iw(1:6)
        itab_v(iv)%iw(1:2) = itab_ud(iud)%im(1:2)

        itab_v(iv)%iv(1:4) = itab_ud(iud)%iu(1:4)
        ! itab_v(iv)%iv(1:12) = itab_ud(iud)%iu(1:12)

        ! For periodic Cartesian hex domain, compute coordinates for outer M points

        im1 = itab_v(iv)%im(1)
        im2 = itab_v(iv)%im(2)
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        if (itab_wd(im1)%npoly < 3) then ! itab_m(im1)%npoly not filled yet
            xem(im1) = xew(iw1) + xew(iw2) - xem(im2)
            yem(im1) = yew(iw1) + yew(iw2) - yem(im2)
            zem(im1) = 0.
        elseif (itab_wd(im2)%npoly < 3) then ! itab_m(im2)%npoly not filled yet
            xem(im2) = xew(iw1) + xew(iw2) - xem(im1)
            yem(im2) = yew(iw1) + yew(iw2) - yem(im1)
            zem(im2) = 0.
        endif

        ! Extract information from IMD1 neighbor

        imd = itab_ud(iud)%im(1)
        npoly = itab_md(imd)%npoly

        do j = 1, npoly
            j1 = j + 1
            if (j == npoly) j1 = 1

            iud1 = itab_md(imd)%iu(j)
            iud2 = itab_md(imd)%iu(j1)

            ! IW(3) and IW(4) neighbors of IV

            if (iud2 == iv) then
                iw1 = itab_ud(iud1)%im(1)
                iw2 = itab_ud(iud1)%im(2)

                if (iw1 == imd) then
                    itab_v(iv)%iw(3) = iw2
                else
                    itab_v(iv)%iw(3) = iw1
                endif
            endif

            if (iud1 == iv) then
                iw1 = itab_ud(iud2)%im(1)
                iw2 = itab_ud(iud2)%im(2)

                if (iw1 == imd) then
                    itab_v(iv)%iw(4) = iw2
                else
                    itab_v(iv)%iw(4) = iw1
                endif
            endif

        enddo

    enddo

    ! Loop over W points

    do iw = 2, nwa
        imd = iw

        itab_w(iw)%loop(1:mloops) = itab_md(imd)%loop(1:mloops)

        itab_w(iw)%iwp = iw

        itab_w(iw)%npoly = itab_md(imd)%npoly
        itab_w(iw)%iwglobe = iw

        itab_w(iw)%mrlw = itab_md(imd)%mrlm
        itab_w(iw)%mrlw_orig = itab_md(imd)%mrlm_orig
        itab_w(iw)%ngr = itab_md(imd)%ngr

        npoly = itab_w(iw)%npoly

        ! Loop over IM/IV neighbors of IW

        do j = 1, itab_w(iw)%npoly
            im = itab_md(imd)%iw(j)
            iwd = im
            iv = itab_md(imd)%iu(j)

            iw1 = itab_v(iv)%iw(1)
            iw2 = itab_v(iv)%iw(2)

            itab_w(iw)%im(j) = im
            itab_w(iw)%iv(j) = iv

            if (iw1 == iw) then
                itab_w(iw)%iw(j) = iw2
                itab_w(iw)%dirv(j) = -1.
            else
                itab_w(iw)%iw(j) = iw1
                itab_w(iw)%dirv(j) = 1.
            endif

        enddo

    enddo

    ! Loop over M points

    do im = 2, nma
        iwd = im

        itab_m(im)%loop(1:mloops) = itab_wd(iwd)%loop(1:mloops)

        itab_m(im)%imp = im

        itab_m(im)%npoly = itab_wd(iwd)%npoly
        itab_m(im)%imglobe = im

        itab_m(im)%mrlm = itab_wd(iwd)%mrlw
        itab_m(im)%ngr = itab_wd(iwd)%ngr

        itab_m(im)%mrlm_orig = itab_wd(iwd)%mrlw_orig
        itab_m(im)%mrow = itab_wd(iwd)%mrow

        itab_m(im)%iv(1:3) = itab_wd(iwd)%iu(1:3)
        itab_m(im)%iw(1:3) = itab_wd(iwd)%im(1:3)
    enddo

    deallocate(itab_md, itab_ud, itab_wd)

END SUBROUTINE voronoi

!-- Performs PCVT (Lloyd's algorithm) iterations to optimize Voronoi cell centroids.
!-- This subroutine iteratively adjusts the locations of Voronoi cell centers (M points)
!-- to be the circumcenters of their defining Delaunay triangle vertices (W points).
!-- This process helps make the Voronoi cells more regular and centroidal.
!-- For global domains, it involves transformations between spherical and polar stereographic (PS) plane coordinates.
!-- Loops over all M points (Voronoi cell centers).
!-- For each M point, identifies its 3 surrounding W points (vertices of the original Delaunay triangle
!-- whose circumcenter would ideally be the M point, or vertices of the Voronoi cell in the dual).
!-- Transforms these W points to a PS plane tangent at the current M point's barycenter (`xebc, yebc, zebc`).
!-- Calculates the circumcenter (`xcc, ycc`) of these 3 W points in the PS plane.
!-- Transforms this PS circumcenter back to Earth coordinates (`dxe, dye, dze`) relative to the barycenter.
!-- Updates the M point's location to this new circumcenter.
!-- Projects the updated M point coordinates back onto the Earth's surface (normalizes to `erad`).
!-- OpenMP directives are used for parallelization.
!===============================================================================

SUBROUTINE pcvt()

    ! Iterative procedure for defining centroidal voronoi cells

    use mem_ijtabs, only : itab_m             ! Neighbor table for M points (itab_m(im)%iw(1:3) are its surrounding W points)
    use mem_grid, only : nma, xem, yem, zem, xew, yew, zew ! Grid dimensions and coordinate arrays for M and W points
                                             ! nma: number of M points (Voronoi cells)
                                             ! xem, yem, zem: coordinates of M points
                                             ! xew, yew, zew: coordinates of W points (Voronoi vertices)
    use consts_coms, only : erad, eradi      ! Earth radius and its inverse

    implicit none

    integer :: im                           ! Loop counter for M points
    integer :: iw1, iw2, iw3                ! Indices of the 3 W points surrounding an M point

    real :: raxis, raxisi                   ! Distance from Earth's axis, and its inverse
    real :: expansion                       ! Expansion factor to project points onto Earth's surface
    real :: sinwlat, coswlat, sinwlon, coswlon! Sin/cos of latitude/longitude for PS transformations
    real :: dxe, dye, dze                   ! Earth-centered delta coordinates (transformed from PS)
    real :: xebc, yebc, zebc                ! Barycentric coordinates of the M point (initial estimate)
    real :: x1, x2, x3, y1, y2, y3          ! PS plane coordinates of the 3 W points
    real :: dx12, dx13, dx23                ! Differences in PS plane x-coordinates
    real :: s1, s2, s3                      ! Sum of squares of PS coordinates (x^2 + y^2) for each W point
    real :: xcc, ycc                        ! PS plane coordinates of the calculated circumcenter

    ! Compute XEM,YEM,ZEM location as circumcentric coordinates of 3 W points.
    ! This establishes W cell as voronoi. ( refer to the goal of making M points circumcenters??)

    ! Loop over all M points

    !$omp parallel
    !$omp do private(iw1,iw2,iw3,xebc,yebc,zebc,raxis,raxisi, &
    !$omp            sinwlat,coswlat,sinwlon,coswlon,dxe,dye,dze,x1,y1, &
    !$omp            x2,y2,x3,y3,dx12,dx13,dx23,s1,s2,s3,ycc,xcc)
    do im = 2, nma ! Loop starts from 2, assuming point 1 might be a pole or special point

        ! Indices of 3 W points surrounding M point (from Voronoi construction, these were Delaunay triangle vertices)
        if (any(itab_m(im)%iw(1:3) < 2)) cycle ! Skip if any W point index is invalid (e.g., boundary cases)

        iw1 = itab_m(im)%iw(1)
        iw2 = itab_m(im)%iw(2)
        iw3 = itab_m(im)%iw(3)

        ! These (xem(im), yem(im), zem(im)) were initialized to be the barycenter of each Delaunay triangle
        ! (which became the initial Voronoi cell center 'im').
        xebc = xem(im)
        yebc = yem(im)
        zebc = zem(im)


        ! For global domain, transform from sphere to PS (Polar Stereographic) plane tangent at (xebc, yebc, zebc)
        raxis = sqrt(xebc ** 2 + yebc ** 2) ! Distance from Earth's rotation axis
        raxisi = 1.0_r8 / raxis             ! Inverse of raxis

        sinwlat = zebc * eradi              ! Sin(latitude of the barycenter)
        coswlat = raxis * eradi             ! Cos(latitude of the barycenter)

        sinwlon = yebc * raxisi             ! Sin(longitude of the barycenter)
        coswlon = xebc * raxisi             ! Cos(longitude of the barycenter)

        ! Transform 3 W points to PS coordinates relative to the tangent point (xebc, yebc, zebc)
        dxe = xew(iw1) - xebc
        dye = yew(iw1) - yebc
        dze = zew(iw1) - zebc
        call de_ps(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, x1, y1) ! de_ps: Earth delta to PS

        dxe = xew(iw2) - xebc
        dye = yew(iw2) - yebc
        dze = zew(iw2) - zebc
        call de_ps(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, x2, y2)

        dxe = xew(iw3) - xebc
        dye = yew(iw3) - yebc
        dze = zew(iw3) - zebc
        call de_ps(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, x3, y3)

        ! Compute intermediate quantities for circumcenter calculation in PS plane
        dx12 = x2 - x1
        dx13 = x3 - x1
        dx23 = x3 - x2

        s1 = x1**2 + y1**2
        s2 = x2**2 + y2**2
        s3 = x3**2 + y3**2

        ! Algebraic solution for circumcenter Y coordinate (ycc) in PS plane
        ycc = 0.5_r8 * (dx13 * s2 - dx12 * s3 - dx23 * s1) &
                / (dx13 * y2 - dx12 * y3 - dx23 * y1)

        ! Algebraic solution for circumcenter X coordinate (xcc) in PS plane
        if (abs(dx12) > abs(dx13)) then
            xcc = (s2 - s1 - ycc * 2.0_r8 * (y2 - y1)) / (2.0_r8 * dx12)
        else
            xcc = (s3 - s1 - ycc * 2.0_r8 * (y3 - y1)) / (2.0_r8 * dx13)
        endif

        ! For global domain, transform circumcenter from PS plane back to Earth delta coordinates (dxe, dye, dze)
        call ps_de(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, xcc, ycc) ! ps_de: PS to Earth delta

        ! Update M point coordinates to the new circumcenter (relative to original barycenter, then absolute)
        xem(im) = dxe + xebc
        yem(im) = dye + yebc
        zem(im) = dze + zebc

    enddo
    !$omp end do nowait

    ! Adjust each M point to the Earth's radius for global domain (project onto sphere)
    !$omp do private(expansion)
    do im = 2, nma

        expansion = erad / sqrt(xem(im) ** 2 &
                + yem(im) ** 2 &
                + zem(im) ** 2)

        xem(im) = xem(im) * expansion
        yem(im) = yem(im) * expansion
        zem(im) = zem(im) * expansion

    enddo
    !$omp end do nowait


    !$omp end parallel

END SUBROUTINE pcvt

!-- Computes various geometric properties of the hexagonal (Voronoi) grid cells.
!-- This subroutine calculates latitudes, longitudes, distances, areas, and normal vectors
!-- for the M (cell center), V (edge midpoint), and W (vertex) points of the Voronoi grid.
!-- It handles transformations to polar stereographic planes for gradient calculations and
!-- computes coefficients for converting velocities between Earth-Cartesian and cell-relative components.
!-- Key calculations include:
!-- Lat/lon of M, V, W points.
!-- `dnu`, `dnv`: Normal distances across U (M-M) and V (W-W) faces/edges.
!-- `unx, vnx`: Unit normal vector components for U and V faces.
!-- `arm0`, `arw0`: Areas associated with M and W points (Voronoi cell area and related kite areas).
!-- `quarter_kite`: Area of quarter-sections of kites formed by M, V, and two W points, used for area distribution.
!-- Gradient coefficients (`gxps1`, `gyps1`, etc.) on a polar stereographic plane.
!-- Coefficients (`ecvec_vx`, etc.) for reconstructing Earth-Cartesian velocity components from cell-face normal components.
!-- OpenMP directives are used for parallelization of loops.
!===============================================================================

subroutine grid_geometry_hex()

    use mem_ijtabs, only : itab_m, itab_v, itab_w ! Neighbor tables for M, V, W points
    use mem_grid, only : nma, nva, nwa, xev, yev, zev, xem, yem, zem, & ! Grid dimensions and coordinates
            xew, yew, zew, unx, uny, unz, wnx, wny, wnz, &             ! Unit vectors (U-face normal, W-point normal)
            vnx, vny, vnz, glonw, glatw, dnu, dniu, dnv, dniv, arw0, & ! V-face normals, W-point lat/lon, distances, areas
            arm0, glonm, glatm, glatv, glonv                           ! M-point areas, M-point lat/lon, V-point lat/lon
    use consts_coms, only : erad, piu180, pio2 ! Physical constants (Earth radius, pi/180, pi/2)
    use consts_coms, only : r8                 ! Real kind for precision

    implicit none

    integer :: im, iv, iw                         !< Loop counters for M, V, W points
    integer :: im1, im2                          !< Indices of M points
    integer :: iw1, iw2                           !< Indices of W points
    integer :: j, npoly                          !< Loop counter, number of polygon edges/vertices
    real :: expansion                            !< Factor to project points onto Earth's surface
    real :: raxis                                !< Distance from Earth's rotation axis
    real :: dvm1, dvm2                           !< Distances from V point to M1 and M2 points
    integer :: j1, j2                            !< Loop counters for polygon edges
    integer :: npoly1, npoly2, np                !< Number of polygon edges for neighbor cells, loop counter
    real :: xv, yv, frac, alpha                 !< PS coordinates, fractional distance, angle
    real :: xw1, xw2, yw1, yw2                  !< PS coordinates of W points
    integer :: iwp, ivp, imp                     !< Global indices for W, V, M points (for periodic boundaries)
    logical :: dops                              !< Flag for deciding whether to do PS calculations
    real :: quarter_kite(2, nva)                 !< Area of a quarter of the kite formed by M-V-W-V geometry
    integer :: lwork                             !< Workspace size for LAPACK's dgels
    integer :: info                               !< Status info from LAPACK's dgels

    real(r8) :: b(7), fo(7)                      !< Arrays for dgels: RHS vector, initial solution guess
    real(r8) :: vnx_ps(7), vny_ps(7), vnz_ps(7)   !< V-normal components projected onto PS plane
    real(r8) :: vrot_x(7), vrot_y(7)             !< Rotated V-normal components in PS plane
    real(r8), allocatable :: work(:)              !< Workspace array for dgels
    real(r8), allocatable :: a(:, :)              !< Coefficient matrix for dgels
    real(r8) :: wsize(1)                         !< Optimal workspace size output from dgels
    real(r8) :: vdotw, vmag, fact               !< Dot product, magnitude, scaling factor

    ! Loop over all M points (Voronoi cell centers)
    !$omp parallel

    !$omp do private(raxis)
    do im = 2, nma
        ! Latitude and longitude at M points
        raxis = sqrt(xem(im) ** 2 + yem(im) ** 2)  ! Distance from earth axis
        glatm(im) = atan2(zem(im), raxis) * piu180 ! Latitude of M point
        glonm(im) = atan2(yem(im), xem(im)) * piu180 ! Longitude of M point

        ! Fill global index (replaced later if this run is parallel and domain is decomposed)
        itab_m(im)%imglobe = im
    enddo
    !$omp end do

    ! Loop over all V points (midpoints of Voronoi cell edges)
    !$omp do private(im1,im2,iw1,iw2,expansion,raxis,dvm1,dvm2,frac)
    do iv = 2, nva
        ! Fill global index
        itab_v(iv)%ivglobe = iv

        ! M-point indices of two end points of V segment (edge connecting two M points, V is its midpoint)
        im1 = itab_v(iv)%im(1) ! First M point connected by this V edge
        im2 = itab_v(iv)%im(2) ! Second M point connected by this V edge

        ! W-point indices on either side of V segment (these are Voronoi vertices, V is on edge between them)
        iw1 = itab_v(iv)%iw(1) ! First W point defining the V edge
        iw2 = itab_v(iv)%iw(2) ! Second W point defining the V edge

        ! V point is midway between W points of Voronoi cells
        xev(iv) = 0.5_r8 * (xew(iw1) + xew(iw2))
        yev(iv) = 0.5_r8 * (yew(iw1) + yew(iw2))
        zev(iv) = 0.5_r8 * (zew(iw1) + zew(iw2))

        !  push V point coordinates out to earth radius
        expansion = erad / sqrt(xev(iv) ** 2 &
                + yev(iv) ** 2 &
                + zev(iv) ** 2)

        xev(iv) = xev(iv) * expansion
        yev(iv) = yev(iv) * expansion
        zev(iv) = zev(iv) * expansion


        ! Latitude and longitude at V point
        raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! Distance from earth axis
        glatv(iv) = atan2(zev(iv), raxis) * piu180 ! Latitude of V point
        glonv(iv) = atan2(yev(iv), xev(iv)) * piu180 ! Longitude of V point

        ! Normal distance across U face (face connecting im1 and im2, with normal along V direction)
        ! Unit vector components of U face normal (points from im1 to im2)
        dnu(iv) = sqrt((xem(im1) - xem(im2))**2 &
                + (yem(im1) - yem(im2))**2 &
                + (zem(im1) - zem(im2))**2)  ! Length of the edge between M1 and M2
        unx(iv) = (xem(im2) - xem(im1)) / dnu(iv) ! U-normal x-component
        uny(iv) = (yem(im2) - yem(im1)) / dnu(iv) ! U-normal y-component
        unz(iv) = (zem(im2) - zem(im1)) / dnu(iv) ! U-normal z-component
        !x dnu(iv) = erad2 * asin(dnu(iv) / erad2) ! (alternative geodesic distance)
        dniu(iv) = 1.0_r8 / dnu(iv) ! Inverse of dnu

        ! Normal distance across V face (face connecting iw1 and iw2, with normal along U direction)
        ! Unit vector components of V face normal (points from iw1 to iw2)
        dnv(iv) = sqrt((xew(iw1) - xew(iw2))**2 &
                + (yew(iw1) - yew(iw2))**2 &
                + (zew(iw1) - zew(iw2))**2)  ! Length of the edge between W1 and W2 (this is the V-edge)
        vnx(iv) = (xew(iw2) - xew(iw1)) / dnv(iv) ! V-normal x-component
        vny(iv) = (yew(iw2) - yew(iw1)) / dnv(iv) ! V-normal y-component
        vnz(iv) = (zew(iw2) - zew(iw1)) / dnv(iv) ! V-normal z-component
        !x dnv(iv) = erad2 * asin(dnv(iv) / erad2) ! (Commented out: alternative geodesic distance)
        dniv(iv) = 1.0_r8 / dnv(iv) ! Inverse of dnv

        ! Skip this V point if iw1 < 2 or iw2 < 2 (boundary condition)
        if (iw1 < 2 .or. iw2 < 2) cycle

        ! refinement level at V point is max of levels at surrounding W points
        itab_v(iv)%mrlv = max(itab_w(iw1)%mrlw, itab_w(iw2)%mrlw)

        ! Compute IM1 and IM2 values of quarter kite area,
        ! and add to ARM0 (area associated with M point) and ARW0 (area associated with W point) arrays.
        ! A kite is formed by M1-V-M2-V'. V is the midpoint of W1-W2.
        ! The area is split into four "quarter-kites" by the U and V edges.
        dvm1 = sqrt((xev(iv) - xem(im1))**2 &
                + (yev(iv) - yem(im1))**2 &
                + (zev(iv) - zem(im1))**2) ! Distance from V to M1
        dvm2 = sqrt((xev(iv) - xem(im2))**2 &
                + (yev(iv) - yem(im2))**2 &
                + (zev(iv) - zem(im2))**2) ! Distance from V to M2

        ! Fractional distance along V edge (W1-W2) where intersection with U edge (M1-M2) is located.
        ! For orthogonal U and V edges in a centroidal Voronoi diagram, V should be the intersection.
        frac = dvm1 * dniu(iv) ! Should be close to 0.5 if V is midpoint of M1-M2 and U,V are orthogonal.
                               ! More accurately, this is projection of M1-V onto M1-M2, normalized by length of M1-M2.

        if (im1 > 1 .and. im2 > 1 .and. (frac < .0001 .or. frac > .9999)) then
            write(io6, *) 'Non-intersecting U-V edges detected in grid geometry'
            write(io6, *) 'FRAC  = ', frac
            write(io6, *) 'IW1 = ', iw1, ' IW2 = ', iw2
            write(io6, *) 'IV    = ', iv
            write(io6, *) 'GLATV = ', glatv(iv)
            write(io6, *) 'GLONV = ', glonv(iv)

            write(io6, *) 'dnu(iv),dniu(iv) ', dnu(iv), dniu(iv)

            stop 'STOP U-V edges'
        endif

        ! Area of the two quarter-kites associated with V point, M1, and (W1 or W2), and V, M2, and (W1 or W2)
        ! Area = 0.5 * base * height; here, 0.5 * (dvm1 or dvm2) * (0.5 * dnv(iv)).
        quarter_kite(1, iv) = 0.25_r8 * dvm1 * dnv(iv) ! Area related to M1
        quarter_kite(2, iv) = 0.25_r8 * dvm2 * dnv(iv) ! Area related to M2
    enddo
    !$omp end do
    !$omp end parallel

    ! Accumulate areas for M and W points from quarter kites
    !dir$ novector ! Hint to compiler not to vectorize this loop (if applicable)
    do iv = 2, nva
        im1 = itab_v(iv)%im(1)
        im2 = itab_v(iv)%im(2)

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        ! Each M point gets area from two quarter-kites per V-edge it's part of.
        arm0(im1) = arm0(im1) + 2.0_r8 * quarter_kite(1, iv)
        arm0(im2) = arm0(im2) + 2.0_r8 * quarter_kite(2, iv)

        ! Each W point's area is sum of all quarter-kites that meet at it.
        ! The V-edge (iw1-iw2) contributes (quarter_kite(1,iv) + quarter_kite(2,iv)) to arw0 for both iw1 and iw2.
        arw0(iw1) = arw0(iw1) + quarter_kite(1, iv) + quarter_kite(2, iv)
        arw0(iw2) = arw0(iw2) + quarter_kite(1, iv) + quarter_kite(2, iv)
    enddo

    ! Lateral boundary copy of arw0 (for periodic domains, copy values from master point)
    do iw = 2, nwa
        iwp = itab_w(iw)%iwp ! iwp is the primary (master) index for this W point
        if (iw /= iwp) arw0(iw) = arw0(iwp)
    enddo

    ! Lateral boundary copy of arm0
    do im = 2, nma
        imp = itab_m(im)%imp ! imp is the primary (master) index for this M point
        if (im /= imp) arm0(im) = arm0(imp)
    enddo

    !$omp parallel
    !$omp do private(raxis)
    do iw = 2, nwa ! Loop over all W points (Voronoi vertices)
        ! Fill global index
        itab_w(iw)%iwglobe = iw

        ! Fill outward unit vector components and latitude and longitude of W point
        raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2) ! Distance from Earth axis
        glatw(iw) = atan2(zew(iw), raxis) * piu180 ! Latitude of W point
        glonw(iw) = atan2(yew(iw), xew(iw)) * piu180 ! Longitude of W point

        ! Normal vector at W point (points radially outward from Earth center)
        wnx(iw) = xew(iw) / erad
        wny(iw) = yew(iw) / erad
        wnz(iw) = zew(iw) / erad
    enddo
    !$omp end do

    !$omp single
    ! Print min/max atmospheric grid spacing (derived from arw0 - area of Voronoi cell)
    iw = minloc(arw0(2:), 1) + 1 ! Find index of min arw0 (excluding point 1)
    write(*, *)
    write(*, '(A,f0.4,A)')       " Minimum atmos grid spacing is ", 0.001_r8 * sqrt(arw0(iw)), " km"
    write(*, '(A,I0,2(A,f0.4))') " at iw = ", iw, ", Lat = ", glatw(iw), ", Lon = ", glonw(iw)

    iw = maxloc(arw0(2:), 1) + 1 ! Find index of max arw0
    write(*, *)
    write(*, '(A,f0.4,A)')       " Maximum atmos grid spacing is ", 0.001_r8 * sqrt(arw0(iw)), " km"
    write(*, '(A,I0,2(A,f0.4))') " at iw = ", iw, ", Lat = ", glatw(iw), ", Lon = ", glonw(iw)
    write(*, *)
    !$omp end single

    !$omp do private(npoly,j2,j1,iv,iw1,iw2,dops,npoly1,npoly2,np,&
    !$omp            im1,xw1,xw2,yw1,yw2,xv,yv,alpha)
    do iw = 2, nwa ! Loop over W points again for gradient and other calculations
        ! Number of polygon edges/vertices for the Voronoi cell centered at W point 'iw'
        npoly = itab_w(iw)%npoly

        ! Loop over all polygon edges of Voronoi cell 'iw'
        do j2 = 1, npoly
            j1 = j2 - 1
            if (j2 == 1) j1 = npoly ! Previous vertex index for the current edge

            iv = itab_w(iw)%iv(j2)   ! V point on the current edge (j2) of polygon iw
            iw2 = itab_w(iw)%iw(j2)  ! Next W point (vertex of polygon iw) along the edge
            iw1 = itab_w(iw)%iw(j1)  ! Previous W point (vertex of polygon iw) defining the edge start

            ! Fractional area of arw0(iw) that is occupied by M and V sectors.
            ! farm: fraction of area of cell 'iw' associated with M-point sectors.
            ! farv: fraction of area of cell 'iw' associated with V-point sectors.
            if (itab_v(iv)%iw(1) == iw1) then ! Check orientation of V-edge relative to W1
                itab_w(iw)%farm(j1) = itab_w(iw)%farm(j1) + quarter_kite(1, iv) / arw0(iw)
                itab_w(iw)%farm(j2) = itab_w(iw)%farm(j2) + quarter_kite(2, iv) / arw0(iw)
            else
                itab_w(iw)%farm(j1) = itab_w(iw)%farm(j1) + quarter_kite(2, iv) / arw0(iw)
                itab_w(iw)%farm(j2) = itab_w(iw)%farm(j2) + quarter_kite(1, iv) / arw0(iw)
            endif
            itab_w(iw)%farv(j2) = (quarter_kite(1, iv) + quarter_kite(2, iv)) / arw0(iw)

            !----------------------------------------
            ! NEW SECTION JULY 2011 (Polar Stereographic calculations for gradients)
            !----------------------------------------

            ! Special - skip gradient calculation if we are at the periodic
            ! domain border and iw1 and iw2 do not share a common M vertex (cell center).
            if (iw == itab_w(iw)%iwp) then ! If iw is its own master point (not a periodic copy)
                dops = .true.
            else ! iw is a periodic copy
                dops = .false.
                npoly1 = itab_w(iw1)%npoly ! Neighbors of iw1
                npoly2 = itab_w(iw2)%npoly ! Neighbors of iw2
                ! Check if iw1 and iw2 share any common M point (cell center) neighbor
                do np = 1, npoly1
                    im1 = itab_w(iw1)%im(np)
                    if (im1 == 1) cycle ! Skip potential pole point
                    if (any(itab_w(iw2)%im(1:npoly2) == itab_w(iw1)%im(np))) then
                        dops = .true. ! Found a common M neighbor, so proceed with PS calculation
                        exit
                    endif
                enddo
            endif

            if (dops) then

                ! Evaluate x,y coordinates of IW1 and IW2 points on polar stereographic plane
                ! tangent at the current W point 'iw'.
                call e_ps(xew(iw1), yew(iw1), zew(iw1), glatw(iw), glonw(iw), xw1, yw1) ! Transform iw1 to PS
                call e_ps(xew(iw2), yew(iw2), zew(iw2), glatw(iw), glonw(iw), xw2, yw2) ! Transform iw2 to PS
                call e_ps(xev(iv), yev(iv), zev(iv), glatw(iw), glonw(iw), xv, yv)     ! Transform V point to PS

                ! Coefficients for eastward and northward components of gradient using values at iw1, iw2
                ! These are part of a finite difference/volume formulation for gradients on the PS plane.
                itab_w(iw)%gxps1(j1) = yw2 / (xw1 * yw2 - xw2 * yw1)  ! Coeff for value at iw1 for x-gradient
                itab_w(iw)%gxps2(j1) = -yw1 / (xw1 * yw2 - xw2 * yw1) ! Coeff for value at iw2 for x-gradient
                itab_w(iw)%gyps1(j1) = -xw2 / (xw1 * yw2 - xw2 * yw1) ! Coeff for value at iw1 for y-gradient
                itab_w(iw)%gyps2(j1) = xw1 / (xw1 * yw2 - xw2 * yw1)  ! Coeff for value at iw2 for y-gradient
                !----------------------------------------

                ! Store PS projection of V-edge normal and its angle (alpha)
                if (itab_w(iw)%dirv(j2) < 0.) then ! dirv indicates direction of V-edge relative to W-polygon
                    alpha = atan2(yw2, xw2)   ! Angle of vector from iw to iw2 in PS plane
                    itab_v(iv)%cosv(1) = cos(alpha) ! Store cos/sin for this V-edge direction
                    itab_v(iv)%sinv(1) = sin(alpha)
                    itab_v(iv)%dxps(1) = xv ! Store PS coordinates of V point relative to W point iw
                    itab_v(iv)%dyps(1) = yv
                else
                    alpha = atan2(-yw2, -xw2) ! Angle of vector from iw to iw1 in PS plane (opposite direction)
                    itab_v(iv)%cosv(2) = cos(alpha)
                    itab_v(iv)%sinv(2) = sin(alpha)

                    itab_v(iv)%dxps(2) = xv
                    itab_v(iv)%dyps(2) = yv
                endif

            endif

            ! Earth-grid components of rotated polar stereographic easterly unit vector at W point 'iw'
            ! (tangent vector pointing East)
            itab_w(iw)%unx_w = -sin(glonw(iw))
            itab_w(iw)%uny_w = cos(glonw(iw))
            ! itab_w(iw)%unz_w = 0 (implicitly, as it's horizontal)

            ! Earth-grid components of rotated polar stereographic northerly unit vector at W point 'iw'
            ! (tangent vector pointing North)
            itab_w(iw)%vnx_w = -sin(glatw(iw)) * cos(glonw(iw))
            itab_w(iw)%vny_w = -sin(glatw(iw)) * sin(glonw(iw))
            itab_w(iw)%vnz_w = cos(glatw(iw))

            !----------------------------------------
            ! END NEW SECTION JULY 2011
            !----------------------------------------

        enddo

    enddo
    !$omp end do

    ! Loop over all V points

    !$omp do private(ivp,iw1,iw2)
    do iv = 2, nva
        ! Let's not do this section on the boundary cells (ivp might be different for periodic boundaries)
        ivp = itab_v(iv)%ivp ! Primary index of V point
        iw1 = itab_v(iv)%iw(1) ! First W point of the V-edge
        iw2 = itab_v(iv)%iw(2) ! Second W point of the V-edge

        ! FARW(1) and FARW(2) interpolation coefficients for ARW (area of W-cell) and VOLT (volume, not used here)
        ! (taking V control volume to be full DNU(IV) * DNV(IV) rectangle - area of kite associated with V-edge)
        ! This seems to be an area weighting factor for interpolating values from W-cells to V-edge.
        itab_v(iv)%farw(1) = 2.0_r8 * (quarter_kite(1, iv) + quarter_kite(2, iv)) / arw0(iw1)
        itab_v(iv)%farw(2) = 2.0_r8 * (quarter_kite(1, iv) + quarter_kite(2, iv)) / arw0(iw2)

        ! Update mrlv (material/refinement level at V point) based on surrounding W points
        itab_v(iv)%mrlv = max(itab_w(iw1)%mrlw, itab_w(iw2)%mrlw)

    enddo  ! IV
    !$omp end do

    ! Scale eastward and northward gradient components by farm (fractional area)
    !$omp do private(j)
    do iw = 2, nwa
        ! The gradient components are not computed at the lateral boundaries (if iw is a copy)
        if (iw /= itab_w(iw)%iwp) cycle

        do j = 1, itab_w(iw)%npoly

            itab_w(iw)%gxps1(j) = itab_w(iw)%gxps1(j) * itab_w(iw)%farm(j)
            itab_w(iw)%gyps1(j) = itab_w(iw)%gyps1(j) * itab_w(iw)%farm(j)

            itab_w(iw)%gxps2(j) = itab_w(iw)%gxps2(j) * itab_w(iw)%farm(j)
            itab_w(iw)%gyps2(j) = itab_w(iw)%gyps2(j) * itab_w(iw)%farm(j)

        enddo
    enddo
    !$omp end do

    ! Coefficients for converting earth-cartesian velocity to V (edge-normal) and W (vertex-centered) components.
    ! This section uses LAPACK's dgels to solve a least-squares problem to find optimal coefficients.
    !$omp do private(npoly, fo, a, b, work, info, j, iv, vdotw, vmag, fact, &
    !$omp            vnx_ps, vny_ps, vnz_ps, vrot_x, vrot_y, wsize, lwork)
    do iw = 2, nwa
        npoly = itab_w(iw)%npoly ! Number of V-edges (and M-neighbors) for this W-cell

        ! Default coefficients from Perot (2000), related to fractional area of V-sectors
        fo(1:npoly) = 2.0_r8 * itab_w(iw)%farv(1:npoly)

        ! Allocate matrix 'a' for dgels if not allocated or size is wrong
        if (allocated(a)) then
            if (size(a, 2) /= npoly) deallocate(a)
        endif
        if (.not. allocated(a)) allocate(a(3, npoly))

        do j = 1, npoly ! Loop over V-edges of the W-cell
            iv = itab_w(iw)%iv(j) ! Current V-edge index

            ! Compute the components of the V unit normals perpendicular to W-normal (wnx, wny, wnz)
            ! This projects V-normals onto the tangent plane at W.
            vdotw = vnx(iv) * wnx(iw) + vny(iv) * wny(iw) + vnz(iv) * wnz(iw) ! Dot product of V-normal and W-normal
            vnx_ps(j) = vnx(iv) - vdotw * wnx(iw) ! Component of V-normal perpendicular to W-normal
            vny_ps(j) = vny(iv) - vdotw * wny(iw)
            vnz_ps(j) = vnz(iv) - vdotw * wnz(iw)

            ! Normalize these new vectors (vnx_ps, etc.) to unit length
            vmag = sqrt(vnx_ps(j)**2 + vny_ps(j)**2 + vnz_ps(j)**2)

            vnx_ps(j) = vnx_ps(j) / vmag
            vny_ps(j) = vny_ps(j) / vmag
            vnz_ps(j) = vnz_ps(j) / vmag

            ! Rotate these new unit normals (now in tangent plane) to a coordinate system
            ! where the local z-axis is aligned with W-normal (wnx,wny,wnz), and x,y are local tangent plane coords.
            if (wnz(iw) >= 0.0_r8) then ! Handles cases near poles for rotation
                fact = (wny(iw) * vnx_ps(j) - wnx(iw) * vny_ps(j)) / (1.0_r8 + wnz(iw))

                vrot_x(j) = vnx_ps(j) * wnz(iw) - vnz_ps(j) * wnx(iw) + wny(iw) * fact ! x-component in local tangent plane
                vrot_y(j) = vny_ps(j) * wnz(iw) - vnz_ps(j) * wny(iw) - wnx(iw) * fact ! y-component in local tangent plane

            else
                fact = (wny(iw) * vnx_ps(j) - wnx(iw) * vny_ps(j)) / (1.0_r8 - wnz(iw))
                vrot_x(j) = -vnx_ps(j) * wnz(iw) + vnz_ps(j) * wnx(iw) + wny(iw) * fact
                vrot_y(j) = -vny_ps(j) * wnz(iw) + vnz_ps(j) * wny(iw) - wnx(iw) * fact

            endif

        enddo

        ! Set up the matrix 'a' for the least-squares problem (dgels)
        a(1, 1:npoly) = vrot_x(1:npoly) * vrot_x(1:npoly)
        a(2, 1:npoly) = vrot_y(1:npoly) * vrot_y(1:npoly)
        a(3, 1:npoly) = vrot_x(1:npoly) * vrot_y(1:npoly)

        ! Set up the RHS vector 'b' for dgels
        b(1) = 1.0_r8 - sum(fo(1:npoly) * a(1, :))
        b(2) = 1.0_r8 - sum(fo(1:npoly) * a(2, :))
        b(3) = - sum(fo(1:npoly) * a(3, :))

        ! Query optimal workspace size for dgels
        call dgels('N', 3, npoly, 1, a, 3, b, 7, wsize, -1, info)
        lwork = nint(wsize(1)) + 1

        ! Allocate workspace array 'work'
        if (allocated(work)) then
            if (size(work) < lwork) deallocate(work)
        endif
        if (.not. allocated(work)) allocate(work(lwork))

        ! Solve the least-squares system A*x = b for x (correction to fo)
        call dgels('N', 3, npoly, 1, a, 3, b, 7, work, size(work), info)

        ! Vector b (now solution x) is the correction to the initial coefficients fo
        b(1:npoly) = b(1:npoly) + fo(1:npoly) ! These are the optimized coefficients

        ! If solution is valid and coefficients are within reasonable bounds, use them
        if (info == 0 .and. all(b(1:npoly) > 0.05_r8) .and. all(b(1:npoly) < 0.7_r8)) then
            ! Calculate Earth-Cartesian vector components for velocity reconstruction at W point
            ! These 'ecvec' components are used to get W-point velocity from V-edge normal velocities.
            fact = sum(b(1:npoly) * vnx_ps(1:npoly) * vnx(itab_w(iw)%iv(1:npoly)))
            fact = (1.0_r8 - wnx(iw)**2) / max(fact, 1.e-30_r8) ! Avoid division by zero
            itab_w(iw)%ecvec_vx(1:npoly) = b(1:npoly) * vnx_ps(1:npoly) * fact

            fact = sum(b(1:npoly) * vny_ps(1:npoly) * vny(itab_w(iw)%iv(1:npoly)))
            fact = (1.0_r8 - wny(iw)**2) / max(fact, 1.e-30_r8)
            itab_w(iw)%ecvec_vy(1:npoly) = b(1:npoly) * vny_ps(1:npoly) * fact

            fact = sum(b(1:npoly) * vnz_ps(1:npoly) * vnz(itab_w(iw)%iv(1:npoly)))
            fact = (1.0_r8 - wnz(iw)**2) / max(fact, 1.e-30_r8)
            itab_w(iw)%ecvec_vz(1:npoly) = b(1:npoly) * vnz_ps(1:npoly) * fact
        else ! If solution is problematic, use default (Perot) coefficients
            write(io6, *) "Problem optimizing vector coefficients for iw = ", iw
            write(io6, *) glatw(iw), glonw(iw)
            write(io6, *) info
            write(io6, *) real(b (1:npoly))
            write(io6, *) real(fo(1:npoly))
            write(io6, *) "Using default coefficients."
            itab_w(iw)%ecvec_vx(1:npoly) = fo(1:npoly) * vnx_ps(1:npoly)
            itab_w(iw)%ecvec_vy(1:npoly) = fo(1:npoly) * vny_ps(1:npoly)
            itab_w(iw)%ecvec_vz(1:npoly) = fo(1:npoly) * vnz_ps(1:npoly)

        endif

    enddo
    !omp end do

    ! Deallocate workspace arrays
    if (allocated(a))    deallocate(a)
    if (allocated(work)) deallocate(work)
    !$omp end parallel
end subroutine grid_geometry_hex


end module MOD_voronoi_geometry

