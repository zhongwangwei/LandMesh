!-- @brief Grid initialization module
!-- @details This module contains subroutines for initializing physical constants
!-- and basic grid structure generation
module MOD_grid_initialization
    implicit none
    
    public :: init_consts, gridinit, gridfile_write

contains

!-- @brief Initializes physical constants, primarily related to Earth's geometry.
!-- @details Sets values for Earth's radius and derived quantities.
subroutine init_consts()
    use consts_coms ! Module containing constants like erad, pio180
    implicit none

    ! Standard (Earth) values
    erad = 6371.22e3          ! Earth radius [m]
    write(io6, *) "erad = ", erad, "in the subroutine init_consts() in the mkgrd.F90"

    ! Secondary values (derived from erad)
    erad2 = erad * 2.0_r8       ! Twice Earth's radius
    erad4 = erad * 4.0_r8       ! Four times Earth's radius
    eradsq = erad * erad        ! Earth's radius squared
    erador5 = erad / sqrt(5.0_r8) ! Earth's radius divided by sqrt(5)
    eradi = 1.0_r8 / erad       ! Inverse of Earth's radius
    erad2sq = erad2 * erad2     ! Square of twice Earth's radius
    dlat = erad * pio180        ! Approximate distance [m] per degree of latitude (pio180 = pi/180)

end subroutine init_consts

!-- @brief Initializes the base grid structure, typically an icosahedral grid.
!-- @details This subroutine performs several steps to set up the initial global grid:
!-- 1. Calls `icosahedron` to generate a base Delaunay triangular grid on a sphere.
!-- 2. Copies this Delaunay grid.
!-- 3. Calls `voronoi` to compute the Voronoi dual of the Delaunay grid.
!-- 4. Calls `pcvt` (Lloyd's algorithm/iterations) to optimize the Voronoi cells towards Centroidal Voronoi Tessellation.
!-- 5. Allocates memory for grid arrays (`alloc_grid1`).
!-- 6. Fills neighborhood information tables (`fill_jtabs`).
!-- 7. Computes grid geometry details (`grid_geometry_hex`).
!-- 8. Copies back the original Delaunay grid information.
!-- 9. Allocates further memory for grid arrays (`alloc_grid2`).
!-- 10. Writes the initial grid to a file (`gridfile_write`).
!===============================================================================
subroutine gridinit()

    use consts_coms, only : io6, nxp             ! Constants (output unit, grid resolution parameter)
    use mem_ijtabs, only : fill_jtabs            ! Subroutine to fill neighbor tables
    use mem_delaunay, only : copy_tri_grid, copyback_tri_grid, nmd, nud, nwd ! Delaunay grid variables and routines
                                                                            ! nmd: number of Delaunay cells (triangles)
                                                                            ! nud: number of Delaunay unique edges
                                                                            ! nwd: number of Delaunay vertices
    use mem_grid, only : nma, nua, nva, nwa, alloc_grid1, alloc_grid2 ! Voronoi/Hex grid variables and allocation routines
                                                                     ! nma: number of Voronoi cells (polygons/hexagons)
                                                                     ! nua: number of Voronoi unique edges (same as nva)
                                                                     ! nva: number of Voronoi unique edges
                                                                     ! nwa: number of Voronoi vertices (same as nmd)
    use MOD_voronoi_geometry, only : voronoi, pcvt, grid_geometry_hex
    implicit none

    ! Horizontal grid setup

    ! Now generate global atmospheric grid
    write(io6, '(/,a)') 'gridinit calling icosahedron'
    call icosahedron(nxp)  ! Generate global spherical Delaunay triangular grid; calls 2 allocs within
    write(io6, '(/,a)') 'gridinit after icosahedron'
    write(io6, '(a,i0)')    ' nmd (Delaunay triangles) = ', nmd
    write(io6, '(a,i0)')    ' nud (Delaunay unique edges) = ', nud
    write(io6, '(a,i0)')    ' nwd (Delaunay vertices) = ', nwd

    ! Store a temporary copy of the full Delaunay mesh
    ! to be used later to construct the surface grid (if needed, or for reference)
    call copy_tri_grid()

    ! Compute Voronoi dual and optimize it
    call voronoi() ! Computes the Voronoi diagram from the Delaunay triangulation
                   ! After this, nma (Voronoi cells) = nwd (Delaunay vertices), etc.
    call pcvt()    ! Performs Lloyd's algorithm (PCVT iterations) to relax the Voronoi mesh
                   ! towards a Centroidal Voronoi Tessellation, making cells more regular.
    write(io6, '(/,a)') 'gridinit after voronoi and pcvt'
    write(io6, '(a,i8)')   ' nma (Voronoi cells) = ', nma
    write(io6, '(a,i8)')   ' nua (Voronoi unique edges A) = ', nua
    write(io6, '(a,i8)')   ' nwa (Voronoi vertices) = ', nwa

    ! Allocate remaining GRID FOOTPRINT arrays for full domain (Voronoi cells)
    write(io6, '(/,a)') 'gridinit calling alloc_grid1 for full domain'
    call alloc_grid1(nma, nva, nwa) ! nva is typically used for edge count in alloc_grid1

    ! Initialize dtlm, dtsm, ndtrat, and nacoust (likely related to time stepping or acoustic properties, not directly used in grid gen)
    ! and compute the timestep schedule for all grid operations.
    write(io6, '(/,a)') 'gridinit calling fill_jtabs'
    write(io6,*), "Calling fill_jtabs with nma, nva, nwa:", nma, nva, nwa
    call fill_jtabs(nma, nva, nwa, 0) ! Fill neighbor tables for Voronoi cells and their vertices/edges

    ! Fill remaining GRID FOOTPRINT geometry for full domain
    write(io6, '(/,a)') 'gridinit calling grid_geometry_hex'
    call grid_geometry_hex() ! Computes geometric properties of the hexagonal (Voronoi) cells

    ! Copy back the Delaunay triangulation data (if it was modified or needed for output)
    call copyback_tri_grid()

    ! Allocate remaining unstructured grid geometry arrays

    write(io6, '(/,a)') 'gridinit calling alloc_grid2'
    call alloc_grid2(nma, nva, nwa) ! Allocate further arrays based on Voronoi grid dimensions

    ! Write GRIDFILE (and potentially SFCGRILE - surface grid file, though not explicitly shown here)
    write(io6, '(/,a)') 'gridinit calling gridfile_write'
    call gridfile_write() ! Writes the generated grid to a NetCDF file

    write(io6, '(/,a)') 'gridinit completed'

end subroutine gridinit

!-- @brief Writes the generated grid data to a NetCDF file.
!-- @details This subroutine takes the primary grid information (cell centers, vertex coordinates,
!-- and connectivity tables for both Delaunay triangles/polygons and Voronoi cells/vertices)
!-- and saves it into a NetCDF file. The filename includes NXP, step, and mode_grid.
!-- Note: The comments mention "sjx_points" (triangles/polygons from Delaunay-like structure)
!-- and "lbx_points" (vertices of these polygons, or centers of Voronoi cells).
!-- It appears to save the primary cell structure (e.g., triangles if `mode_grid` is `tri`,
!-- or polygons if `mode_grid` is `hex` after Voronoi generation).
!===============================================================================

subroutine gridfile_write()
    ! Does not calculate dismm and disww because they are not used.
    use netcdf
    use consts_coms, only : r8, pathlen, io6, file_dir, EXPNME, NXP, mode_grid, refine, step ! Global constants and parameters
    use mem_ijtabs, only : mloops, itab_m, itab_w  ! Neighbor tables for M (cell centers) and W (vertices) points
                                                 ! itab_m(im)%iw(1:3) are vertices of triangle 'im'
                                                 ! itab_w(iw)%im(1:7) are cells surrounding vertex 'iw' (up to 7 for hex)
    use mem_grid, only : nma, nwa, glatw, glonw, glatm, glonm ! Grid dimensions and coordinate arrays
                                                              ! nma: number of main cells (e.g., triangles or hexagons)
                                                              ! nwa: number of vertices
                                                              ! glonm, glatm: longitude/latitude of cell centers (M points)
                                                              ! glonw, glatw: longitude/latitude of cell vertices (W points)
    use MOD_utilities, only : Unstructured_Mesh_Save ! Subroutine to handle NetCDF writing

    implicit none

    ! This routine writes the grid variables to the gridfile.
    integer :: im, iw                           ! Loop counters for cells (im) and vertices (iw)
    integer :: sjx_points, lbx_points           ! Number of "sjx_points" (polygons/triangles) and "lbx_points" (vertices defining them)
                                                 ! Typically, sjx_points = nma, lbx_points = nwa for this context.
    real(r8), allocatable :: mp(:,:),wp(:,:)    ! Arrays to hold cell center (mp) and vertex (wp) coordinates (lon, lat)
    integer, allocatable :: ngrmw(:,:),ngrwm(:,:)! Connectivity arrays:
                                                 ! ngrmw: for each cell (M), its vertices (W)
                                                 ! ngrwm: for each vertex (W), its surrounding cells (M)
    integer, allocatable :: n_ngrwm(:)          ! Number of neighboring cells for each vertex (not explicitly used in Unstructured_Mesh_Save call shown)
    character(pathlen) :: lndname               ! Output NetCDF filename
    character(5) :: nxpc,stepc                  ! Character representation of NXP and step for filename

    ! Populate ngrmw: for each cell (M point, index im), list its vertices (W points)
    ! For a triangular mesh (mode_grid='tri'), each cell im has 3 vertices.
    ! For a hexagonal mesh (mode_grid='hex'), itab_m would store neighbors differently if it were primary.
    ! Here, it seems to assume a triangular primary structure (Delaunay like) where itab_m lists 3 vertices.
    allocate (ngrmw(3, nma)) ; ngrmw = 0
    do im = 1, nma
       ngrmw(1:3, im) = itab_m(im)%iw(1:3) ! Get the 3 vertices for triangle 'im'
    enddo

    ! Populate ngrwm: for each vertex (W point, index iw), list its surrounding cells (M points)
    ! A vertex can be shared by up to 7 cells in a hexagonal grid (less for boundaries/pentagons).
    allocate (ngrwm(7, nwa)) ; ngrwm = 0
    !@RuiZhang: should 7 be replace by itab_w(iw)%ngr ?
    do iw = 1, nwa
       ngrwm(1:itab_w(iw)%ngr, iw) = itab_w(iw)%im(1:itab_w(iw)%ngr) ! Get surrounding cells for vertex 'iw'
                                                                   ! itab_w(iw)%ngr is the actual number of neighbors
    enddo

    ! Prepare coordinate arrays for saving
    allocate(mp(nma,2)); mp(:,1) = GLONM; mp(:,2) = GLATM ! Cell centers: (lon, lat)
    allocate(wp(nwa,2)); wp(:,1) = GLONW; wp(:,2) = GLATW ! Vertices: (lon, lat)

    ! Construct the output filename
    write(nxpc, '(I4.4)')NXP    ! Format NXP (e.g., 0030)
    write(stepc, '(I2.2)') step ! Format step (e.g., 01)
    lndname = trim(file_dir)// 'gridfile/gridfile_NXP' // trim(nxpc) // '_'//trim(stepc)// '_'// trim(mode_grid)// '.nc4'

    ! Print information and set point counts for saving
    write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(io6, *) 'grid_write: opening file:', lndname
    sjx_points = nma ! Number of polygons/triangles (M points)
    lbx_points = nwa ! Number of vertices (W points)
    write(io6, *) 'sjx_points (cells):', sjx_points
    write(io6, *) 'lbx_points (vertices):', lbx_points
    write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    ! Save the unstructured mesh data to NetCDF
    ! Initial file without any refinement (step usually 1 here)
    CALL Unstructured_Mesh_Save(lndname, sjx_points, lbx_points, mp, wp, ngrmw, ngrwm)

    ! Deallocate temporary arrays
    deallocate(ngrwm, ngrmw, mp, wp)

END SUBROUTINE gridfile_write

end module MOD_grid_initialization 