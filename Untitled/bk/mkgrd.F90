!DESCRIPTION
!===========
!===============================================================================
!-- @brief Main program for unstructured mesh generation.
!-- @details This program generates unstructured meshes for land surface models
!-- (e.g., CoLM). It handles various mesh types (land, ocean, earth),
!-- grid generation modes (hexagonal, triangular, lon-lat, lambert),
!-- and optional mesh refinement.
!===============================================================================


!REVISION HISTORY
!----------------
! 2025.06.11  Zhongwang Wei @ SYSU
! 2025.06.10  Rui Zhang @ SYSU
! 2023.02.21  Zhongwang Wei @ SYSU
! 2021.12.02  Zhongwang Wei @ SYSU 
! 2020.10.01  Zhongwang Wei @ SYSU

! The original code is from OLAM
!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

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
program main
    use netcdf
    use consts_coms                                                  ! Module containing physical and mathematical constants
    use refine_vars                                                  ! Module containing variables related to mesh refinement
    use MOD_data_preprocess , only : data_preprocess                 ! Module for preprocessing input data
    use MOD_grid_preprocess , only : grid_preprocess                 ! Module for preprocessing grid data
    use MOD_Area_judge      , only : Area_judge, Area_judge_refine   ! Module for judging areas for refinement/masking
    use MOD_GetContain      , only : Get_Contain                     ! Module to determine if points are contained in areas
    use MOD_GetRef          , only : GetRef                          ! Module to get refinement flags
    use MOD_Refine          , only : refine_loop                     ! Module for the main mesh refinement loop
    use MOD_mask_postproc       , only : mask_postproc               ! Module for post-processing masks
    use MOD_namelist        , only : read_nl                         ! Module for reading namelist files
    use MOD_grid_initialization, only : init_consts, gridinit        ! Module for grid initialization
    use MOD_utilities       , only : mode4mesh_make                  ! Module for utility functions

    implicit none
    character(pathlen) :: nlfile = 'mkgrd.mnl'                       ! Namelist file name, default 'mkgrd.mnl'
    character(pathlen) :: finfolist                                  ! Path to save the namelist file
    character(pathlen) :: lndname                                    ! Land grid file name
    character(5)       :: stepc                                      ! Character representation of the current step
    character(5)       :: nxpc                                       ! Character representation of NXP (grid resolution parameter)
    logical            :: exit_loop                                  ! Flag to exit the refinement loop
    logical            :: fexists                                    ! Flag to check if a file exists
    integer :: i, ncid, dimID_sjx, sjx_points, length                ! Index, ncid, dimID_sjx, sjx_points, length
    io6 = 6                                                          ! If run is sequential, default choice is to set io6 to standard output unit 6.
    CALL getarg(1, nlfile)                                           ! Get the namelist file name from the command line
    call read_nl(nlfile)                                             ! Read settings from the namelist file
    
    ! Check if the mesh type is valid
    if ((mesh_type /= 'landmesh')   .and. &
        (mesh_type /= 'oceanmesh')  .and. &
        (mesh_type /= 'earthmesh')) then
        write(io6, *) "ERROR! mesh_type = ", mesh_type
        STOP "ERROR! mesh_type mush be landmesh/oceanmesh/earthmesh"
    end if

    step = 1                                                         ! Initialize current refinement step
    num_vertex = 1                                                   ! Initialize number of vertices
    ! Handle mask restart scenario
    if (mask_restart) then
        call init_consts()                                           ! Initialize constants
        refine = .false.    ! Disable refinement for mask restart
        step = max_iter + 1 ! Set step beyond max_iter to skip refinement loop
        if ((mesh_type == 'oceanmesh') .and. (.not. mask_patch_on)) then
            ! This is for the case of only adjusting mask_sea_ratio, if you want to patch, you need to do it from the later process
            !please add description here @RuiZhang
            write(io6, *) "Remask_restart start"
            CALL mask_postproc(mesh_type) 
            write(io6, '(A)') "--------------------------------"
            write(io6, '(A)') ""
            write(io6, '(A)') "!! Successfully Make Grid End !!"
            write(io6, '(A)') ""
            write(io6, '(A)') "--------------------------------"
            stop "Remask_restart finish" ! Stop execution after remasking
        end if
    end if
   
    ! Save a copy of the namelist file
    finfolist = trim(file_dir)//'result/namelist.save' ! Define path for saving namelist
    CALL execute_command_line('cp '//trim(nlfile)//' '//trim(finfolist)) ! Copy namelist to result directory

    ! Grid generation logic based on mode_grid type
    inquire(file = mode_file, exist = fexists) ! Check if the mode_file (initial grid file) exists
    if ((mode_grid == 'hex') .or. &    ! Hexagonal grid
        (mode_grid == 'tri')) then     ! Triangular grid
        if (fexists) then ! If mode_file exists, use it as a base
            inquire(file = mode_file_description, exist = fexists) ! Check for a description file (currently leads to error if exists)
            if (.not. fexists) then
                write(io6, *) "mode_file is exist!"
                ! Check the consistency of mode_file
                CALL CHECK(NF90_OPEN(trim(mode_file), nf90_nowrite, ncid))
                CALL CHECK(NF90_INQ_DIMID(ncid, "sjx_points", dimID_sjx)) ! Get dimension ID for triangle/polygon points
                CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_sjx, len = sjx_points)) ! Get number of points
                CALL CHECK(NF90_CLOSE(ncid)) ! Close the NetCDF file
                if (int((sjx_points-1)/20) /= int(nxp*nxp)) then ! Validate NXP consistency
                    write(io6, *) "nxp read from namelist diff from nxp in the mode_file"
                    stop
                    ! write(io6, *) "nxp turn to ", nxp ," keep the same value as the mode_file"
                end if
            else
                stop "ERROR! can not use now" ! Error if mode_file_description exists
            end if
            ! Prepare output grid file name and copy the mode_file
            write(nxpc, '(I4.4)') NXP
            write(stepc, '(I2.2)') step
            lndname = trim(file_dir)// 'gridfile/gridfile_NXP' // trim(nxpc) // '_'// trim(stepc) //'_'//trim(mode_grid)// '.nc4'
            CALL execute_command_line('cp '//trim(mode_file)//' '//trim(lndname)) ! Copy mode_file to the output grid file
            write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(io6, *) 'grid_write: opening file:', lndname
            write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            call init_consts() ! Initialize constants
        else ! If mode_file does not exist, generate grid from scratch
            ! Initialize, execute, and end olam run (conceptual steps for grid generation)
            call init_consts() ! Initialize constants
            call gridinit()    ! Initialize the grid structure (e.g., icosahedron)
        end if
        write(io6, *) 'grid preporces start'
        CALL grid_preprocess() ! Preprocess the generated or loaded grid
        write(io6, *) 'grid preporces finish'

    else if ((mode_grid == 'lonlat' ) .or. &  ! Longitude-Latitude grid
             (mode_grid == 'lambert')) then ! Lambert conformal conic projection grid
        if (fexists) then ! mode_file must exist for these types
            write(io6, *) 'mode4mesh_make start'
            inquire(file = mode_file, exist = fexists) ! Double check existence (already done)
            if (.not. fexists) then
                write(*, *) "The input file " // trim(mode_file) // " is missing."
                stop "Stopping model run."
            endif
            CALL mode4mesh_make(mode_file, mode_grid) ! Create mesh from the mode_file
            write(io6, *) 'mode4mesh_make complete'
            write(io6, *) ""
        else
            write(io6, *) "ERROR! mode_file must fexists when mode_grid as ", mode_grid
            stop
        end if

        if (refine) then ! If refinement was intended, disable it for lonlat/lambert as it's handled differently or not supported
            write(io6, *) "turn refine from TRUE to FALSE"
            refine = .false.
        else
            write(io6, *) "refine is FALSE"
        end if

    else if ((mode_grid == 'dbx') .or. &     ! Direct read from database/file (future use)
             (mode_grid == 'cubical')) then  ! Cubical sphere grid (future use)
        ! only nc/nc4 format in current stage
        stop "ERROR! mode_grid == dbx/cubical can not use now!"
        length = len_trim(mode_file)
        if ('nml' == mode_file(length-2:length)) stop 'ERROR! can not use now in the dbx'
        if (refine) then
            write(io6, *) "turn refine from TRUE to FALSE when choose dbx"
            refine = .false.
         end if
    else
        stop 'ERROR mode_grid !!!' ! Invalid mode_grid type
    end if
    write(io6, *) 'mode_grid set as ', mode_grid    


    write(io6, *) 'data preporcess start'
    ! Preset some necessary data such as landtype 
    CALL data_preprocess()
    write(io6, *) 'data preporcess complete'
    write(io6, *) ""

    write(io6, *) 'area judge start'
    ! Multiple boundaries options available in the DmArea(RfArea)
    CALL Area_judge() ! Determine domain area (DmArea) or refinement area (RfArea), handles mask-patch modifications
    write(io6, *) 'area judge complete'
    write(io6, *) ""

    ! Main refinement loop
    if (refine) then
        write(io6, *) 'make grid with refine mesh'
        write(io6, *) ""
 
        max_iter = max(max_iter_cal, max_iter_spc) ! Determine maximum iterations from calculated and specified values
        write(io6, *) "max_iter_spc = ", max_iter_spc ! Max iterations for specified refinement (read from namelist)
        write(io6, *) "max_iter_cal = ", max_iter_cal ! Max iterations for calculated (threshold-based) refinement (read from namelist)
        write(io6, *) "max_iter = ", max_iter
        if (max_iter <= 0) stop 'Error! max_iter must more than zero'

        ! Validate halo vs. max_transition_row settings
        do i = 1, max_iter, 1
            if (halo(i) < max_transition_row(i)) then
                write(io6, *) 'i = ', i
                write(io6, *) 'halo(i) = ', halo(i)
                write(io6, *) 'max_transition_row(i) = ', max_transition_row(i)
                stop "ERROR! halo must larger than max_transition_row!"
            end if    
        end do

        write(io6, *) 'start do-while'
        exit_loop = .false. ! Initialize exit_loop flag
        do while(step <= max_iter) ! Loop through refinement steps
            write(io6, *) 'step = ',step, 'in the refine-circle'
            write(io6, *) 'Get ref_sjx start'
            ! Get refinement flags (ref_sjx) for triangles/polygons
            ! Only calculate for newly-generated triangles or polygons

            ! Threshold-based refinement
            if (refine_setting == 'calculate' .or. refine_setting == 'mixed') then
                if (step <= max_iter_cal) then
                    CALL Area_judge_refine(0)      ! Judge refinement area (0 indicates threshold-based)
                    CALL Get_Contain(0)            ! Determine containment for threshold-based
                    CALL GetRef(0, exit_loop)      ! Get refinement flags for threshold-based
                end if
            end if

            ! Specified refinement (e.g., based on predefined regions)
            if (refine_setting == 'specified' .or. refine_setting == 'mixed') then
                if (step <= max_iter_spc) then
                    ! Update map for specified refinement
                    CALL Area_judge_refine(step)   ! Judge refinement area for specified (step indicates current iteration)
                    CALL Get_Contain(step)         ! Determine containment for specified
                    CALL GetRef(step, exit_loop)   ! Get refinement flags for specified
                end if
            end if
            write(io6, *) 'Get ref_sjx complete'

            write(io6, *) 'refine_loop start'
            CALL refine_loop(exit_loop) ! Perform the actual mesh refinement based on ref_sjx
            write(io6, *) 'refine_loop complete'
            write(io6, *) ""
           
            if (exit_loop) then ! If GetRef indicated no more cells to refine
                print *, 'Exiting loop due to ref_sjx equal to zero! &
                          turn refine from to True to False !'
                exit ! Exit the DO WHILE loop
            end if
           
            step = step + 1 ! Increment refinement step
        end do
        write(io6, *) 'finish do-while'
    else
        write(io6, *) 'make grid with basic mesh' ! No refinement performed
    end if

    ! Final processing steps
    refine = .false. ! Ensure refine is false after the loop (or if it was never true)
                     ! This means no more cells need further calculation/refinement.
    CALL Get_Contain(0) ! Final containment check (purpose might need more context, 0 suggests a general pass)
    CALL mask_postproc(mesh_type) ! Post-process masks (e.g., apply land/sea mask)

    ! Success message
    write(io6, '(A)') "--------------------------------"
    write(io6, '(A)') ""
    write(io6, '(A)') "!! Successfully Make Grid End !!"
    write(io6, '(A)') ""
    write(io6, '(A)') "--------------------------------"

end program main



!-- @brief Reads namelist settings for grid generation and refinement.
!-- @details This subroutine reads two namelist groups: `mkgrd` for general grid
!-- generation parameters and `mkrefine` for refinement-specific parameters.
!-- It initializes various global variables based on these settings, creates
!-- necessary output directories, and calls `Mask_make` to process mask files
!-- for domain, refinement, and patches.
!-- @param nlfile The path to the main namelist file.
subroutine read_nl(nlfile)
    use consts_coms ! Module for constants and common variables
    use refine_vars ! Module for refinement related variables
    implicit none

    character(*), intent(in) :: nlfile     ! Input: Path to the namelist file
    integer :: i, pos, iostat             ! Loop counters, position, I/O status
    logical :: fexists                    ! Flag to check file existence
    character(pathlen) :: path, fprefix, filename, lndname! Path, file prefix, filename, land name (various path manipulations)

    namelist /mkgrd/ nl    ! Namelist group for general grid parameters
    namelist /mkrefine/ rl ! Namelist group for refinement parameters
    ! OPEN THE NAMELIST FILE
    inquire(file = nlfile, exist = fexists) ! Check if the namelist file exists
    write(io6, *) nlfile
    if (.not. fexists) then
        write(*, *) "The namelist file " // trim(nlfile) // " is missing."
        stop "Stopping model run."
    endif
    open(iunit, status = 'OLD', file = nlfile) ! Open the namelist file
    ! READ GRID POINT, MODEL OPTIONS, AND PLOTTING INFORMATION FROM THE NAMELIST
    REWIND(iunit)
    !write(io6, *)nl
    read(iunit, nml = mkgrd)   ! Read the 'mkgrd' namelist
    close(iunit)
    write(*, nml = mkgrd)      ! Print the read 'mkgrd' parameters for verification
    write(io6, *) ""

    !----------------------------------------------------------
    ! Assign values from 'mkgrd' namelist to global variables.
    ! Variables in the following section either must not be changed on a history
    ! restart or changing them would be irrelevant.  Thus, they are only copied
    ! from the namelist if a history file is not being read.
    expnme               = nl%expnme           ! Experiment name
    nxp                  = nl%nxp              ! Grid resolution parameter (e.g., for icosahedral grids)
    GXR                  = nl%GXR              ! Grid expansion ratio (placeholder, needs more context)
    base_dir             = nl%base_dir         ! Base directory for input/output
    source_dir           = nl%source_dir       ! Source data directory
    mesh_type            = nl%mesh_type        ! Type of mesh: 'landmesh', 'oceanmesh', 'earthmesh'
    mode_grid            = nl%mode_grid        ! Grid generation mode: 'hex', 'tri', 'lonlat', 'lambert', etc.
    mode_file            = nl%mode_file        ! Path to initial grid file (if any)
    mode_file_description= nl%mode_file_description! Description for the mode file (optional)
    refine               = nl%refine           ! Logical flag to enable/disable mesh refinement
    lcs                  = nl%lcs              ! Logical flag for local coordinate system (placeholder)
    openmp               = nl%openmp           ! Logical flag to enable OpenMP parallelization
    mask_sea_ratio       = nl%mask_sea_ratio   ! Threshold for sea ratio in mask generation
    mask_restart         = nl%mask_restart     ! Logical flag for mask restart mode
    mask_domain_type     = nl%mask_domain_type ! Type of domain mask ('bbox', 'lambert', 'circle', 'close')
    mask_domain_fprefix  = nl%mask_domain_fprefix! File prefix for domain mask files
    mask_patch_on        = nl%mask_patch_on    ! Logical flag to enable patch masks
    mask_patch_type      = nl%mask_patch_type  ! Type of patch mask
    mask_patch_fprefix   = nl%mask_patch_fprefix ! File prefix for patch mask files
    file_dir             = trim(base_dir) // trim(expnme) // '/'! Construct full experiment directory path

    if (GXR < 0) stop "ERROR! GXR must >= 0" ! Validate GXR

    ! If mask_restart is enabled, handle patch mask and return (skip full grid generation)
    if (mask_restart) then
        if (mask_patch_on) CALL Mask_make('mask_patch', mask_patch_type, mask_patch_fprefix) ! Create patch mask if enabled
        return
    end if

    ! Create output directories for the current experiment
    CALL execute_command_line('rm -rf '//trim(file_dir)) ! Remove old experiment directory if it exists
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"contain/")   ! Directory for containment check files
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"gridfile/")  ! Directory for generated grid files
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"patchtype/") ! Directory for patch type files
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"result/")    ! Directory for final mesh result
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"tmpfile/")   ! Directory for temporary files
    CALL execute_command_line('rm *_filelist.txt') ! Remove any existing file list temporary files
    ! Create domain mask
    CALL Mask_make('mask_domain', mask_domain_type, mask_domain_fprefix)
    ! Create patch mask if enabled
    if (mask_patch_on) CALL Mask_make('mask_patch', mask_patch_type, mask_patch_fprefix)

    ! If refinement is enabled, read 'mkrefine' namelist and set refinement parameters
    if (refine) then
        CALL execute_command_line('mkdir -p '//trim(file_dir)//"threshold/") ! Directory for threshold-based refinement files
        open(iunit, status = 'OLD', file = nlfile) ! Re-open namelist file
        REWIND(iunit)
        read(iunit, nml = mkrefine) ! Read the 'mkrefine' namelist
        close(iunit)
        write(*, nml = mkrefine)    ! Print read 'mkrefine' parameters for verification

        ! Assign values from 'mkrefine' namelist
        weak_concav_eliminate = rl%weak_concav_eliminate! Flag for weak concavity elimination
        Istransition          = rl%Istransition         ! Flag for transition zone handling
        max_sa_iter           = rl%max_sa_iter          ! Max iterations for simulated annealing (or similar optimization)
        halo                  = rl%halo                 ! Halo size for refinement transitions
        max_transition_row    = rl%max_transition_row   ! Max rows in transition zone
        if (Istransition == .false. .and. mode_grid /= 'tri') STOP "ERROR! not Istransition can only use in the tri" ! Validate Istransition setting

        refine_spc            = rl%refine_spc           ! Flag for specified refinement
        refine_cal            = rl%refine_cal           ! Flag for calculated (threshold-based) refinement
        if (refine_spc) max_iter_spc = rl%max_iter_spc  ! Max iterations for specified refinement (default 0, read if refine_spc is true)
        if (refine_cal) max_iter_cal = rl%max_iter_cal  ! Max iterations for calculated refinement (default 0, read if refine_cal is true)

        ! Determine overall refinement setting based on refine_spc and refine_cal
        if ((refine_spc .eqv. .TRUE.) .and. (refine_cal .eqv. .TRUE.)) then
            refine_setting = 'mixed'    ! Both specified and calculated refinement
        else if ((refine_spc .eqv. .TRUE.)  .and. (refine_cal .eqv. .FALSE.)) then
            refine_setting = 'specified'! Only specified refinement
        else if ((refine_spc .eqv. .FALSE.) .and. (refine_cal .eqv. .TRUE.)) then
            refine_setting = 'calculate'! Only calculated refinement
        else
            stop "ERROR! MUst one of TRUE in the refine_spc and refine_cal when refine is TRUE"
        end if
        write(io6, *) "refine_setting = ", refine_setting

        ! Specified refinement settings
        if (refine_setting == 'specified' .or. refine_setting == 'mixed') then
            mask_refine_spc_type       = RL%mask_refine_spc_type   ! Type of mask for specified refinement
            mask_refine_spc_fprefix    = RL%mask_refine_spc_fprefix! File prefix for specified refinement mask files
            CALL Mask_make('mask_refine', mask_refine_spc_type, mask_refine_spc_fprefix) ! Create specified refinement mask
            if (mask_refine_ndm(max_iter_spc) == 0) then ! Validate mask creation for specified refinement
                write(io6, *) "max_iter_spc = ", max_iter_spc
                write(io6, *) "mask_refine_ndm(max_iter_spc) = ", mask_refine_ndm(max_iter_spc)
                stop "ERROR! mask_refine_ndm(max_iter_spc) must larger then one, please modify max_iter_spc"
            end if
        end if

        ! Threshold-based (calculated) refinement settings
        if (refine_setting == 'calculate' .or. refine_setting == 'mixed') then
            ! Land surface related refinement criteria
            if ((mesh_type == 'landmesh') .or. (mesh_type == 'earthmesh')) then
                refine_num_landtypes      = rl%refine_num_landtypes   ! Refine based on number of land types
                refine_area_mainland      = rl%refine_area_mainland   ! Refine based on mainland area
                ! Single-layer land surface variables for refinement decision
                refine_onelayer_Lnd( 1)   = rl%refine_lai_m       ! Mean LAI
                refine_onelayer_Lnd( 2)   = rl%refine_lai_s       ! Stddev LAI
                refine_onelayer_Lnd( 3)   = rl%refine_slope_m     ! Mean slope
                refine_onelayer_Lnd( 4)   = rl%refine_slope_s     ! Stddev slope
                ! Two-layer (e.g., soil) land surface variables for refinement
                refine_twolayer_Lnd( 1)   = rl%refine_k_s_m       ! Mean soil hydraulic conductivity
                refine_twolayer_Lnd( 2)   = rl%refine_k_s_s       ! Stddev soil hydraulic conductivity
                refine_twolayer_Lnd( 3)   = rl%refine_k_solids_m  ! Mean soil solids thermal conductivity
                refine_twolayer_Lnd( 4)   = rl%refine_k_solids_s  ! Stddev soil solids thermal conductivity
                refine_twolayer_Lnd( 5)   = rl%refine_tkdry_m     ! Mean dry soil thermal conductivity
                refine_twolayer_Lnd( 6)   = rl%refine_tkdry_s     ! Stddev dry soil thermal conductivity
                refine_twolayer_Lnd( 7)   = rl%refine_tksatf_m    ! Mean frozen saturated soil thermal conductivity
                refine_twolayer_Lnd( 8)   = rl%refine_tksatf_s    ! Stddev frozen saturated soil thermal conductivity
                refine_twolayer_Lnd( 9)   = rl%refine_tksatu_m    ! Mean unfrozen saturated soil thermal conductivity
                refine_twolayer_Lnd(10)   = rl%refine_tksatu_s    ! Stddev unfrozen saturated soil thermal conductivity

                ! Corresponding thresholds for the above land variables
                th_num_landtypes          = rl%th_num_landtypes
                th_area_mainland          = rl%th_area_mainland
                th_onelayer_Lnd( 1)       = rl%th_lai_m
                th_onelayer_Lnd( 2)       = rl%th_lai_s
                th_onelayer_Lnd( 3)       = rl%th_slope_m
                th_onelayer_Lnd( 4)       = rl%th_slope_s
                th_twolayer_Lnd( 1, 1:2)  = rl%th_k_s_m        ! Thresholds (e.g., min/max)
                th_twolayer_Lnd( 2, 1:2)  = rl%th_k_s_s
                th_twolayer_Lnd( 3, 1:2)  = rl%th_k_solids_m
                th_twolayer_Lnd( 4, 1:2)  = rl%th_k_solids_s
                th_twolayer_Lnd( 5, 1:2)  = rl%th_tkdry_m
                th_twolayer_Lnd( 6, 1:2)  = rl%th_tkdry_s
                th_twolayer_Lnd( 7, 1:2)  = rl%th_tksatf_m
                th_twolayer_Lnd( 8, 1:2)  = rl%th_tksatf_s
                th_twolayer_Lnd( 9, 1:2)  = rl%th_tksatu_m
                th_twolayer_Lnd(10, 1:2)  = rl%th_tksatu_s
            end if

            ! Ocean related refinement criteria
            if ((mesh_type == 'oceanmesh') .or. (mesh_type == 'earthmesh')) then
                refine_sea_ratio          = rl%refine_sea_ratio       ! Refine based on sea ratio
                refine_Rossby_radius      = rl%refine_Rossby_radius   ! Refine based on Rossby radius of deformation
                ! Single-layer ocean variables
                refine_onelayer_Ocn(1)    = rl%refine_sst_m         ! Mean Sea Surface Temperature
                refine_onelayer_Ocn(2)    = rl%refine_sst_s         ! Stddev SST
                refine_onelayer_Ocn(3)    = rl%refine_ssh_m         ! Mean Sea Surface Height
                refine_onelayer_Ocn(4)    = rl%refine_ssh_s         ! Stddev SSH
                refine_onelayer_Ocn(5)    = rl%refine_eke_m         ! Mean Eddy Kinetic Energy
                refine_onelayer_Ocn(6)    = rl%refine_eke_s         ! Stddev EKE
                refine_onelayer_Ocn(7)    = rl%refine_sea_slope_m   ! Mean sea surface slope
                refine_onelayer_Ocn(8)    = rl%refine_sea_slope_s   ! Stddev sea surface slope

                ! Corresponding thresholds for ocean variables
                th_sea_ratio              = rl%th_sea_ratio
                th_Rossby_radius          = rl%th_Rossby_radius
                th_onelayer_Ocn(1)        = rl%th_sst_m
                th_onelayer_Ocn(2)        = rl%th_sst_s
                th_onelayer_Ocn(3)        = rl%th_ssh_m
                th_onelayer_Ocn(4)        = rl%th_ssh_s
                th_onelayer_Ocn(5)        = rl%th_eke_m
                th_onelayer_Ocn(6)        = rl%th_eke_s
                th_onelayer_Ocn(7)        = rl%th_sea_slope_m
                th_onelayer_Ocn(8)        = rl%th_sea_slope_s
            end if

            ! Earth system (coupled model) related refinement criteria
            if (mesh_type == 'earthmesh') then
                refine_onelayer_Earth( 1) = rl%refine_typhoon_m   ! Mean typhoon intensity/track related variable
                refine_onelayer_Earth( 2) = rl%refine_typhoon_s   ! Stddev typhoon variable
                th_onelayer_Earth( 1)     = rl%th_typhoon_m
                th_onelayer_Earth( 2)     = rl%th_typhoon_s
            end if

            ! Validation: If calculated refinement is on, at least one specific criterion must be enabled.
            if (refine_setting == 'calculate' .or. refine_setting == 'mixed') then
                if (mesh_type == 'landmesh') then
                    if ((refine_num_landtypes .eqv. .false.) .and. &
                        (refine_area_mainland .eqv. .false.) .and. &
                        (all(refine_onelayer_Lnd  .eqv. .false.)).and. &
                        (all(refine_twolayer_Lnd  .eqv. .false.))) then
                        stop "Error! MUst one of TRUE in the refine_num_landtypes or &
                                refine_area_mainland or refine_onelayer_Lnd or refine_twolayer_Lnd &
                                when refine is TRUE and meshtype = landmesh"
                    end if

                else if (mesh_type == 'oceanmesh') then
                    if ((refine_sea_ratio .eqv. .false.) .and. &
                        (refine_Rossby_radius .eqv. .false.) .and. &
                        (all(refine_onelayer_Ocn .eqv. .false.))) then
                        stop "ERROR! MUst one of TRUE in the refine_sea_ratio or refine_onelayer_Ocn when refine is TRUE and meshtype = oceanmesh"
                    end if

                else if (mesh_type == 'earthmesh') then
                    if ((refine_num_landtypes .eqv. .false.) .and. &
                        (refine_area_mainland .eqv. .false.) .and. &
                        (refine_sea_ratio     .eqv. .false.) .and. &
                        (refine_Rossby_radius .eqv. .false.) .and. &
                        (all(refine_onelayer_Lnd  .eqv. .false.)) .and. &
                        (all(refine_twolayer_Lnd  .eqv. .false.)) .and. &
                        (all(refine_onelayer_Ocn  .eqv. .false.)) .and. &
                        (all(refine_onelayer_Earth .eqv. .false.))) then
                        write(io6, *) "refine_num_landtypes = ", refine_num_landtypes
                        write(io6, *) "refine_area_mainland = ", refine_area_mainland
                        write(io6, *) "refine_sea_ratio = ", refine_sea_ratio
                        write(io6, *) "refine_Rossby_radius = ", refine_Rossby_radius
                        stop "Error! MUst one of TRUE in the refine_sea_ratio or refine_Rossby_radius or &
                                refine_num_landtypes or &
                                refine_area_mainland or refine_onelayer_Lnd or refine_twolayer_Lnd or &
                                refine_onelayer_Ocn or refine_onelayer_Earth &
                                when refine is TRUE and meshtype = earthmesh"
                    end if
                end if
            end if
            mask_refine_cal_type       = RL%mask_refine_cal_type   ! Type of mask for calculated refinement
            mask_refine_cal_fprefix    = RL%mask_refine_cal_fprefix! File prefix for calculated refinement mask files
            CALL Mask_make('mask_refine', mask_refine_cal_type, mask_refine_cal_fprefix) ! Create calculated refinement mask
        end if

        ! Validation checks for consistency between refinement flags and threshold values (999. indicates unset threshold)
        ! onelayer_Lnd
        if ((refine_onelayer_Lnd( 1) .eqv. .true.) .and. (th_onelayer_Lnd( 1) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Lnd( 1) and   th_onelayer_Lnd( 1) "
        if ((refine_onelayer_Lnd( 2) .eqv. .true.) .and. (th_onelayer_Lnd( 2) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Lnd( 2) and   th_onelayer_Lnd( 2) "
        if ((refine_onelayer_Lnd( 3) .eqv. .true.) .and. (th_onelayer_Lnd( 3) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Lnd( 3) and   th_onelayer_Lnd( 3) "
        if ((refine_onelayer_Lnd( 4) .eqv. .true.) .and. (th_onelayer_Lnd( 4) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Lnd( 4) and   th_onelayer_Lnd( 4) "
        ! twolayer_Lnd
        if ((refine_twolayer_Lnd( 1) .eqv. .true.) .and. any(th_twolayer_Lnd( 1, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 1)  and      th_twolayer_Lnd( 1, 1:2) "
        if ((refine_twolayer_Lnd( 2) .eqv. .true.) .and. any(th_twolayer_Lnd( 2, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 2)  and      th_twolayer_Lnd( 2, 1:2) "
        if ((refine_twolayer_Lnd( 3) .eqv. .true.) .and. any(th_twolayer_Lnd( 3, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 3)  and      th_twolayer_Lnd( 3, 1:2) "
        if ((refine_twolayer_Lnd( 4) .eqv. .true.) .and. any(th_twolayer_Lnd( 4, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 4)  and      th_twolayer_Lnd( 4, 1:2) "
        if ((refine_twolayer_Lnd( 5) .eqv. .true.) .and. any(th_twolayer_Lnd( 5, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 5)  and      th_twolayer_Lnd( 5, 1:2) "
        if ((refine_twolayer_Lnd( 6) .eqv. .true.) .and. any(th_twolayer_Lnd( 6, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 6)  and      th_twolayer_Lnd( 6, 1:2) "
        if ((refine_twolayer_Lnd( 7) .eqv. .true.) .and. any(th_twolayer_Lnd( 7, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 7)  and      th_twolayer_Lnd( 7, 1:2) "
        if ((refine_twolayer_Lnd( 8) .eqv. .true.) .and. any(th_twolayer_Lnd( 8, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 8)  and      th_twolayer_Lnd( 8, 1:2) "
        if ((refine_twolayer_Lnd( 9) .eqv. .true.) .and. any(th_twolayer_Lnd( 9, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 9)  and      th_twolayer_Lnd( 9, 1:2) "
        if ((refine_twolayer_Lnd(10) .eqv. .true.) .and. any(th_twolayer_Lnd(10, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd(10)  and      th_twolayer_Lnd(10, 1:2) "

        ! onelayer_Ocn
        if ((refine_onelayer_Ocn( 1) .eqv. .true.) .and. (th_onelayer_Ocn( 1) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 1) and   th_onelayer_Ocn( 1) "
        if ((refine_onelayer_Ocn( 2) .eqv. .true.) .and. (th_onelayer_Ocn( 2) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 2) and   th_onelayer_Ocn( 2) "
        if ((refine_onelayer_Ocn( 3) .eqv. .true.) .and. (th_onelayer_Ocn( 3) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 3) and   th_onelayer_Ocn( 3) "
        if ((refine_onelayer_Ocn( 4) .eqv. .true.) .and. (th_onelayer_Ocn( 4) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 4) and   th_onelayer_Ocn( 4) "
        if ((refine_onelayer_Ocn( 5) .eqv. .true.) .and. (th_onelayer_Ocn( 5) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 5) and   th_onelayer_Ocn( 5) "
        if ((refine_onelayer_Ocn( 6) .eqv. .true.) .and. (th_onelayer_Ocn( 6) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 6) and   th_onelayer_Ocn( 6) "
        if ((refine_onelayer_Ocn( 7) .eqv. .true.) .and. (th_onelayer_Ocn( 7) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 7) and   th_onelayer_Ocn( 7) "
        if ((refine_onelayer_Ocn( 8) .eqv. .true.) .and. (th_onelayer_Ocn( 8) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 8) and   th_onelayer_Ocn( 8) "

        ! onelayer_Earth
        if ((refine_onelayer_Earth( 1) .eqv. .true.) .and. (th_onelayer_Earth( 1) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Earth( 1) and   th_onelayer_Earth( 1) "
        if ((refine_onelayer_Earth( 2) .eqv. .true.) .and. (th_onelayer_Earth( 2) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Earth( 2) and   th_onelayer_Earth( 2) "

    end if

end subroutine read_nl

!-- @brief Manages the creation of different types of masks (domain, refine, patch).
!-- @details This subroutine identifies mask files based on a prefix and type,
!-- then calls specific subroutines (`bbox_mask_make`, `lamb_mask_make`,
!-- `circle_mask_make`, `close_mask_make`) to process and create them.
!-- It lists files matching `mask_fprefix*`, then iterates through them.
!-- @param mask_select Type of mask being made ('mask_domain', 'mask_refine', 'mask_patch').
!-- @param type_select Geometric type of the mask ('bbox', 'lambert', 'circle', 'close').
!-- @param mask_fprefix Full path and prefix for the input mask definition files.
SUBROUTINE Mask_make(mask_select, type_select, mask_fprefix)
    ! Common preprocessing for domain/refine/patch masks.
    use consts_coms     ! Module for constants and common variables
    use refine_vars     ! Module for refinement related variables
    use netcdf          ! For NetCDF operations (though not directly used here, used by called subroutines)
    implicit none
    character(*), intent(in) :: mask_select   ! Input: Specifies the purpose of the mask ('mask_domain', 'mask_refine', 'mask_patch')
    character(*), intent(in) :: type_select   ! Input: Specifies the geometric type of the mask ('bbox', 'lambert', 'circle', 'close')
    character(*), intent(in) :: mask_fprefix  ! Input: File prefix (including path) for mask definition files
    integer :: pos, i, iostat                 ! Position, loop counter, I/O status
    logical :: fexists                      ! Flag to indicate if matching files are found
    character(pathlen) :: path, fprefix, filename, lndname, command! Path, file prefix (parsed), current filename, list filename, system command

    pos = 0 ! Initialize position of the last '/'
    ! Find the last '/' in mask_fprefix to separate path and prefix
    do i = len_trim(mask_fprefix), 1, -1
        if (mask_fprefix(i:i) == '/') then
            pos = i
            exit
        end if
    end do

    if (pos > 0) then ! If a path separator is found
        path = mask_fprefix(1:pos)       ! Extract path part
        fprefix = mask_fprefix(pos+1:)   ! Extract file prefix part
        print *, 'path: ',    trim(adjustl(path))
        print *, 'fprefix: ', trim(adjustl(fprefix))
        fexists = .false. ! Initialize found flag

        ! Use system command to list all files matching the prefix and store in a temporary file
        lndname = trim(mask_select) // '_filelist.txt' ! Name of the temporary file list
        command = 'ls ' // trim(mask_fprefix) // '* > ' // trim(lndname) ! ls command
        call execute_command_line(command) ! Execute the command

        ! Open the temporary file to read the list of filenames
        open(unit=111, file=lndname, status='old')
        REWIND(111)
        ! Iterate through the file list
        do while (.true.)
            read(111,'(A)', iostat=iostat) filename ! Read a filename
            write(io6, *) "iostat = ", iostat, "mask_select = ", mask_select
            if (iostat < 0) exit ! End of file or error
            if (iostat > 0) then
                print *, "File read error in Mask_make in mkgrd.F90"
                stop
            end if
            ! If the filename contains the specific prefix, set found to .true. and process
            if (index(trim(adjustl(filename)), trim(adjustl(fprefix)) ) > 0) then
                write(io6, *) "filename : ", trim(adjustl(filename))
                fexists = .true.
                ! Call the appropriate mask making subroutine based on type_select
                if (type_select == 'bbox') then
                    CALL bbox_mask_make(filename, mask_select)
                else if (type_select == 'lambert') then
                    CALL lamb_mask_make(filename, mask_select)
                else if (type_select == 'circle') then
                    CALL circle_mask_make(filename, mask_select) 
                else if (type_select == 'close') then
                    CALL close_mask_make(filename, mask_select)
                else
                    write(io6, *) "ERROR! ", type_select, " must be bbox, lambert, circle, close"
                    stop
                end if
            end if
        end do

        ! Close the file list
        close(111)

        ! Check if any matching files were found
        if (.not. fexists) stop 'No matching files found. Incorrect path format in mask_refine_fprefix.'
    else
        write(io6, *) "Incorrect path format in mask_fprefix: ", trim(mask_fprefix)
        stop
    end if

END SUBROUTINE Mask_make

!-- @brief Creates a bounding box (bbox) mask.
!-- @details Reads bbox definitions (west, east, north, south coordinates and a refinement degree)
!-- from either a Namelist file (.nml) or a NetCDF file (.nc, .nc4).
!-- It then saves this information into a standardized NetCDF file in the `tmpfile/` directory.
!-- The refinement degree is checked against `max_iter_spc`.
!-- @param inputfile Path to the input file (Namelist or NetCDF) defining the bbox.
!-- @param mask_select Type of mask being created ('mask_domain', 'mask_refine', 'mask_patch'),
!-- used for naming the output file and incrementing counters.
subroutine bbox_mask_make(inputfile, mask_select)
    ! Converts file type to NetCDF for easier reading; if already NetCDF, assigns directly.
    use netcdf
    use consts_coms, only : r8, pathlen, file_dir, mask_domain_ndm, mask_patch_ndm
    use refine_vars, only : mask_refine_ndm, max_iter_spc
    use MOD_file_preprocess, only : bbox_Mesh_Save ! Subroutine to save bbox data to NetCDF

    implicit none
    character(pathlen), intent(in) :: inputfile   ! Input: Path to the bbox definition file (.nml or .nc/.nc4)
    character(*), intent(in)       :: mask_select ! Input: Type of mask ('mask_domain', 'mask_refine', 'mask_patch')
    integer :: ncid, varid                        ! NetCDF file ID, variable ID
    character(pathlen) :: lndname, line           ! Output NetCDF filename, line read from nml
    logical :: fexists                            ! File existence flag (not used in this scope)
    integer :: i, bbox_num, refine_degree, length, io_stat! Loop counter, number of bboxes, refinement degree, string length, I/O status
    real(r8), allocatable :: bbox_points(:,:)     ! Array to store bbox coordinates [west, east, north, south]
    character(5) :: numc, refinec                ! Character strings for numbering and refinement degree in filename

    length = len_trim(inputfile) ! Get length of input filename

    ! Process based on input file type
    if ('nml' == inputfile(length-2:length)) then ! If input is a Namelist file
        ! Namelist format:
        ! First line: number of points (bboxes)
        ! Second line: refinement degree
        ! Subsequent lines: west, east, north, south for each bbox
        open(unit=10, file=inputfile, status='old', action='read')
        REWIND(10)
        ! Read bbox_num (number of bounding boxes)
        read(10, '(A)', iostat=io_stat) line  ! Read the whole line
        if (io_stat /= 0) then
            print *, "Error: Failed to read bbox_num"
            stop
        end if
        read(line(index(line, '=')+1:), *) bbox_num  ! Parse bbox_num from the line "bbox_num = X"

        ! Read refine_degree
        read(10, '(A)', iostat=io_stat) line  ! Read the whole line
        if (io_stat /= 0) then
            print *, "Error: Failed to read refine_degree for bbox" ! Changed from circle_refine to generic
            stop
        end if
        read(line(index(line, '=')+1:), *) refine_degree  ! Parse refine_degree
        if (refine_degree > max_iter_spc) then ! Check if refine_degree exceeds max allowed for specified refinement
            write(io6, *) "refine_degree > max_iter_spc :", refine_degree, ">", max_iter_spc
            close(10)
            return ! Do not process this mask further
        end if

        write(io6, *) "bbox_points must be 1.West 2.East 3.North and 4.South"
        allocate(bbox_points(bbox_num, 4)) ! Allocate array for bbox coordinates
        do i = 1, bbox_num, 1
            ! Must ensure bbox_points are in the range -180 to 180 (lon) and -90 to 90 (lat)
            read(10, *) bbox_points(i, 1), bbox_points(i, 2), bbox_points(i, 3), bbox_points(i, 4) ! Read W, E, N, S
            if (bbox_points(i, 1) > bbox_points(i, 2)) stop "ERROR! bbox_points(i, 1) > bbox_points(i, 2) (West > East)"
            if (bbox_points(i, 3) < bbox_points(i, 4)) stop "ERROR! bbox_points(i, 3) < bbox_points(i, 4) (North < South)"
        end do
        close(10)

        ! Generate output filename and save to NetCDF
        write(refinec, '(I1)') refine_degree ! Format refinement degree for filename
        if (mask_select == 'mask_domain') then
            mask_domain_ndm = mask_domain_ndm + 1 ! Increment domain mask counter
            write(numc, '(I2.2)') mask_domain_ndm ! Format counter for filename
        else if (mask_select == 'mask_refine') then
            mask_refine_ndm(refine_degree) = mask_refine_ndm(refine_degree) + 1 ! Increment refinement mask counter for this degree
            write(numc, '(I2.2)') mask_refine_ndm(refine_degree)
        else if (mask_select == 'mask_patch') then
            mask_patch_ndm = mask_patch_ndm + 1 ! Increment patch mask counter
            write(numc, '(I2.2)') mask_patch_ndm
        end if

        lndname = trim(file_dir)// 'tmpfile/'//trim(mask_select)// '_bbox_'//trim(refinec)//'_'//trim(numc)//'.nc4'
        CALL bbox_Mesh_Save(lndname, bbox_num, bbox_points) ! Save bbox data
        deallocate(bbox_points)

    else if (('.nc' == inputfile(length-2:length)) .or. &  ! If input is a NetCDF file
             ('nc4' == inputfile(length-2:length))) then
        CALL CHECK(NF90_OPEN(trim(inputfile), nf90_nowrite, ncid)) ! Open the NetCDF file
        CALL CHECK(NF90_INQ_VARID(ncid, 'bbox_refine', varid))    ! Get variable ID for 'bbox_refine'
        CALL CHECK(NF90_GET_VAR(ncid, varid, refine_degree))      ! Read refinement degree
        CALL CHECK(NF90_CLOSE(ncid))                              ! Close NetCDF file
        if (refine_degree > max_iter_spc) then
            write(io6, *) "refine_degree > max_iter_spc :", refine_degree, ">", max_iter_spc
            return ! Do not process further
        end if

        ! Generate output filename and copy the input NetCDF file
        write(refinec, '(I1)') refine_degree
        if (mask_select == 'mask_domain') then
            mask_domain_ndm = mask_domain_ndm + 1
            write(numc, '(I2.2)') mask_domain_ndm
        else if (mask_select == 'mask_refine') then
            mask_refine_ndm(refine_degree) = mask_refine_ndm(refine_degree) + 1
            write(numc, '(I2.2)') mask_refine_ndm(refine_degree)
        else if (mask_select == 'mask_patch') then
            mask_patch_ndm = mask_patch_ndm + 1
            write(numc, '(I2.2)') mask_patch_ndm
        end if

        lndname = trim(file_dir)// 'tmpfile/'//trim(mask_select)// '_bbox_'//trim(refinec)//'_'//trim(numc)//'.nc4'
        CALL execute_command_line('cp '//trim(inputfile)//' '//trim(lndname)) ! Copy input file to standardized name in tmpfile dir
    else
        stop "ERROR! must choose 'ncfile' or 'nml' please check!" ! Invalid file type
    end if
    write(io6, *) lndname ! Print the name of the created/copied NetCDF mask file

end subroutine bbox_mask_make

!-- @brief Creates a Lambert projection mask.
!-- @details Reads Lambert projection grid vertex data (lon_vert, lat_vert) from a NetCDF input file.
!-- It converts these 2D vertex coordinates into a 1D list of unique vertices (`lonlat_bound`)
!-- and defines the cell connectivity (`ngr_bound`) for quadrilateral cells.
!-- This structured mask data is then saved to a new NetCDF file in the `tmpfile/` directory.
!-- Currently, this is intended for use in threshold-based calculations (`refine_degree` is set to 0).
!-- Namelist input for Lambert masks is not supported.
!-- @param inputfile Path to the input NetCDF file defining the Lambert grid.
!-- @param mask_select Type of mask ('mask_domain', 'mask_refine', 'mask_patch').
subroutine lamb_mask_make(inputfile, mask_select)
! Currently only usable for NetCDF files and only in threshold calculations.
    use netcdf
    use consts_coms, only : r8, pathlen, file_dir, mask_domain_ndm, mask_patch_ndm
    use refine_vars, only : mask_refine_ndm, max_iter_spc
    use MOD_file_preprocess, only : Mode4_Mesh_Save ! Subroutine to save mesh data

    implicit none
    character(pathlen), intent(in) :: inputfile   ! Input: Path to the Lambert grid definition file (.nc or .nc4)
    character(*), intent(in)       :: mask_select ! Input: Type of mask ('mask_domain', 'mask_refine', 'mask_patch')
    character(pathlen)             :: lndname       ! Output NetCDF filename
    logical                        :: fexists       ! File existence flag (not used in this scope)
    integer :: i, j, idx, length                  ! Loop counters, index, string length
    integer :: nlon_cal, nlat_cal                 ! Calculated number of lon/lat points (not used here)
    integer :: ncid, dimID_lon, dimID_lat         ! NetCDF IDs for file and dimensions
    integer :: lon_points, lat_points, bound_points, mode_points, refine_degree! Grid dimensions, number of boundary/mode points, refinement degree
    integer, dimension(10)         :: varid         ! Array for NetCDF variable IDs
    real(r8), dimension(:, :), allocatable :: lon_vert, lat_vert, lonlat_bound! 2D arrays for vertex lon/lat and 1D combined bounds
    integer,  dimension(:, :), allocatable :: ngr_bound     ! Connectivity array for cell bounds
    integer,  dimension(:),    allocatable :: n_ngr         ! Number of neighbors per cell
    character(5)                   :: nxpc, stepc, numc, refinec! Character strings for file naming

    length = len_trim(inputfile) ! Get length of input filename

    if ('nml' == inputfile(length-2:length)) then ! Namelist input not supported for Lambert masks
        stop 'ERROR! nml can not use in lambert now'
    else if (('.nc' == inputfile(length-2:length))  .or. &  ! If input is NetCDF
             ('nc4' == inputfile(length-2:length))) then
        CALL CHECK(NF90_OPEN(inputfile, nf90_nowrite, ncid)) ! Open NetCDF file
        ! The following lines for reading 'lamb_refine' are commented out, implying refine_degree is fixed or handled differently.
        ! CALL CHECK(NF90_INQ_VARID(ncid, 'lamb_refine', varid))
        ! CALL CHECK(NF90_GET_VAR(ncid, varid, refine_degree))
        ! if (refine_degree > max_iter_spc) then
        !     write(io6, *) "refine_degree > max_iter_spc :", refine_degree, ">", max_iter_spc
        !     CALL CHECK(NF90_CLOSE(ncid))
        !     return ! Exit if refinement degree is too high
        ! end if
        CALL CHECK(NF90_INQ_DIMID(ncid, "xi_vert",  dimID_lon))      ! Get x-dimension ID (e.g., "xi_vert")
        CALL CHECK(NF90_INQ_DIMID(ncid, "eta_vert", dimID_lat))      ! Get y-dimension ID (e.g., "eta_vert")
        CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lon, len = lon_points)) ! Get number of points in x-dim
        CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lat, len = lat_points)) ! Get number of points in y-dim
        allocate(lon_vert(lon_points, lat_points))
        allocate(lat_vert(lon_points, lat_points))
        CALL CHECK(NF90_INQ_VARID(ncid, "lon_vert", varid(1)))     ! Get var ID for longitude vertices
        CALL CHECK(NF90_INQ_VARID(ncid, "lat_vert", varid(2)))     ! Get var ID for latitude vertices
        CALL CHECK(NF90_GET_VAR(ncid, varid(1), lon_vert))         ! Read longitude vertex data
        CALL CHECK(NF90_GET_VAR(ncid, varid(2), lat_vert))         ! Read latitude vertex data
        CALL CHECK(NF90_CLOSE(ncid)) ! Close NetCDF file
        ! Because the input is boundary values, the number of grid cells is one less.
        lon_points = lon_points - 1 ! Adjust to number of cells in x-direction
        lat_points = lat_points - 1 ! Adjust to number of cells in y-direction
        where(lon_vert > 180.)  lon_vert = lon_vert - 360. ! Normalize longitudes
    else
        stop "ERROR! must choose 'ncfile' or 'nml' please check!" ! Invalid file type
    end if

    ! Convert 2D vertex data to 1D boundary list and connectivity (similar to mode4mesh_make)
    write(io6, *) "lonlat_bound and ngr_bound calculate start"
    write(io6, *) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    bound_points = (lon_points + 1) * (lat_points + 1) + 1
    mode_points = lon_points * lat_points + 1
    allocate(lonlat_bound(bound_points, 2)); lonlat_bound = -999.
    allocate(ngr_bound(4, mode_points)); ngr_bound = 1
    allocate(n_ngr(mode_points)); n_ngr = 4

    idx = 1
    do j = 1, lat_points + 1, 1
        do i = 1, lon_points + 1, 1
            idx = idx + 1
            lonlat_bound(idx, :) = [lon_vert(i, j), lat_vert(i, j)]
        end do
    end do
    
    ! ngr_bound calculate, never start from 1 !
    idx = 1
    do j = 1, lat_points, 1
        do i = 1, lon_points, 1
            idx = idx + 1
            ngr_bound(:, idx) = [i + (j - 1) * (lon_points + 1),     &
                                 i + (j - 1) * (lon_points + 1) + 1, &
                                 i +  j      * (lon_points + 1) + 1, &
                                 i +  j      * (lon_points + 1)]
        end do
    end do
    ngr_bound = ngr_bound + 1 ! Adjust indices

    write(io6, *) "lon_points : ", lon_points
    write(io6, *) "lat_points : ", lat_points
    write(io6, *) "bound_points : ", bound_points
    write(io6, *) "mode_points  : ", mode_points
    write(io6, *) "minval(lonlat_bound) : ", minval(lonlat_bound(2:bound_points, 1))
    write(io6, *) "maxval(lonlat_bound) : ", maxval(lonlat_bound(2:bound_points, 1))
    write(io6, *) "maxval(lonlat_bound) : ", maxval(lonlat_bound(2:bound_points, 2))
    write(io6, *) "minval(lonlat_bound) : ", minval(lonlat_bound(2:bound_points, 2))
    write(io6, *) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    write(io6, *) "lonlat_bound and ngr_bound calculate finish"
    write(io6, *) ""

    refine_degree = 0 ! For threshold calculations, refine_degree is typically 0
    write(refinec, '(I1)') refine_degree
    if (mask_select == 'mask_domain') then
        mask_domain_ndm = mask_domain_ndm + 1
        write(numc, '(I2.2)') mask_domain_ndm
    else if (mask_select == 'mask_refine') then
        mask_refine_ndm(refine_degree) = mask_refine_ndm(refine_degree) + 1
        write(numc, '(I2.2)') mask_refine_ndm(refine_degree)
    else if (mask_select == 'mask_patch') then
        mask_patch_ndm = mask_patch_ndm + 1
        write(numc, '(I2.2)') mask_patch_ndm
    end if
    lndname = trim(file_dir)// 'tmpfile/'//trim(mask_select)// '_lambert_'//trim(refinec)//'_'//trim(numc)//'.nc4'
    write(io6, *) lndname
    CALL Mode4_Mesh_Save(lndname, bound_points, mode_points, lonlat_bound, ngr_bound, n_ngr) ! Save the processed mask data
    write(io6, *) "lamb_mask_make finish"

end subroutine lamb_mask_make

!-- @brief Creates a circular mask.
!-- @details Reads circle definitions (center longitude, center latitude, radius in km, and
!-- refinement degree) from either a Namelist file (.nml) or a NetCDF file (.nc, .nc4).
!-- It then saves this information into a standardized NetCDF file in the `tmpfile/` directory.
!-- The refinement degree is checked against `max_iter_spc`.
!-- @param inputfile Path to the input file (Namelist or NetCDF) defining the circle(s).
!-- @param mask_select Type of mask being created ('mask_domain', 'mask_refine', 'mask_patch').
subroutine circle_mask_make(inputfile, mask_select)
    ! Converts file type to NetCDF for easier reading; if already NetCDF, assigns directly.
    use netcdf
    use consts_coms, only : r8, pathlen, file_dir, mask_domain_ndm, mask_patch_ndm
    use refine_vars, only : mask_refine_ndm, max_iter_spc
    use MOD_file_preprocess, only : circle_Mesh_Save ! Subroutine to save circle data to NetCDF

    implicit none
    character(pathlen), intent(in) :: inputfile   ! Input: Path to the circle definition file (.nml or .nc/.nc4)
    character(*), intent(in)       :: mask_select ! Input: Type of mask ('mask_domain', 'mask_refine', 'mask_patch')
    integer :: ncid, varid                        ! NetCDF file ID, variable ID
    character(pathlen) :: lndname, line           ! Output NetCDF filename, line read from nml
    logical :: fexists                            ! File existence flag (not used in this scope)
    integer :: i, circle_num, refine_degree, length, io_stat! Loop counter, number of circles, refinement degree, string length, I/O status
    real(r8), allocatable :: circle_points(:,:), circle_radius(:)! Arrays for circle center (lon, lat) and radius (km)
    character(5) :: numc, refinec                ! Character strings for numbering and refinement degree in filename

    length = len_trim(inputfile) ! Get length of input filename

    ! Process based on input file type
    if ('nml' == inputfile(length-2:length)) then ! If input is a Namelist file
        ! Namelist format:
        ! First line: number of circles
        ! Second line: refinement degree
        ! Subsequent lines: center_lon, center_lat, radius_km for each circle
        open(unit=10, file=inputfile, status='old', action='read')
        REWIND(10)
        ! Read circle_num (number of circles)
        read(10, '(A)', iostat=io_stat) line  ! Read the whole line
        if (io_stat /= 0) then
            print *, "Error: Failed to read circle_num"
            close(10)
            stop
        end if
        read(line(index(line, '=')+1:), *) circle_num  ! Parse circle_num

        ! Read refine_degree
        read(10, '(A)', iostat=io_stat) line  ! Read the whole line
        if (io_stat /= 0) then
            print *, "Error: Failed to read circle_refine" ! Name of parameter in nml
            close(10)
            stop
        end if
        read(line(index(line, '=')+1:), *) refine_degree  ! Parse refine_degree (referred to as close_refine in comment, but circle_refine in read)
        if (refine_degree > max_iter_spc) then
            write(io6, *) "refine_degree > max_iter_spc :", refine_degree, ">", max_iter_spc
            close(10)
            return ! Do not process further
        end if

        allocate(circle_points(circle_num, 2)) ! Allocate for center coordinates
        allocate(circle_radius(circle_num))   ! Allocate for radii
        do i = 1, circle_num, 1
            read(10, *) circle_points(i,1), circle_points(i,2), circle_radius(i) ! Read lon, lat, radius
        end do
        close(10)

        ! Save to NetCDF and deallocate
        write(refinec, '(I1)') refine_degree
        if (mask_select == 'mask_domain') then
            mask_domain_ndm = mask_domain_ndm + 1
            write(numc, '(I2.2)') mask_domain_ndm
        else if (mask_select == 'mask_refine') then
            mask_refine_ndm(refine_degree) = mask_refine_ndm(refine_degree) + 1
            write(numc, '(I2.2)') mask_refine_ndm(refine_degree)
        else if (mask_select == 'mask_patch') then
            mask_patch_ndm = mask_patch_ndm + 1
            write(numc, '(I2.2)') mask_patch_ndm
        end if
        lndname = trim(file_dir)// 'tmpfile/'//trim(mask_select)// '_circle_'//trim(refinec)//'_'//trim(numc)//'.nc4'
        CALL circle_Mesh_Save(lndname, circle_num, circle_points, circle_radius)
        deallocate(circle_points, circle_radius)

    else if (('.nc' == inputfile(length-2:length)) .or. &  ! If input is a NetCDF file
             ('nc4' == inputfile(length-2:length))) then
        CALL CHECK(NF90_OPEN(trim(inputfile), nf90_nowrite, ncid)) ! Open NetCDF
        CALL CHECK(NF90_INQ_VARID(ncid, 'circle_refine', varid))  ! Variable name for refinement degree
        CALL CHECK(NF90_GET_VAR(ncid, varid, refine_degree))    ! Read refinement degree
        CALL CHECK(NF90_CLOSE(ncid))                            ! Close NetCDF
        if (refine_degree > max_iter_spc) then
             write(io6, *) "refine_degree > max_iter_spc :", refine_degree, ">", max_iter_spc
            return ! Do not process further
        end if

        ! Generate output filename and copy input NetCDF
        write(refinec, '(I1)') refine_degree
        if (mask_select == 'mask_domain') then
            mask_domain_ndm = mask_domain_ndm + 1
            write(numc, '(I2.2)') mask_domain_ndm
        else if (mask_select == 'mask_refine') then
            mask_refine_ndm(refine_degree) = mask_refine_ndm(refine_degree) + 1
            write(numc, '(I2.2)') mask_refine_ndm(refine_degree)
        else if (mask_select == 'mask_patch') then
            mask_patch_ndm = mask_patch_ndm + 1
            write(numc, '(I2.2)') mask_patch_ndm
        end if
        lndname = trim(file_dir)// 'tmpfile/'//trim(mask_select)// '_circle_'//trim(refinec)//'_'//trim(numc)//'.nc4'
        CALL execute_command_line('cp '//trim(inputfile)//' '//trim(lndname)) ! Copy to standardized name
    else
        stop "ERROR! must choose 'ncfile' or 'nml' please check!" ! Invalid file type
    end if
   write(io6, *) lndname ! Print name of created/copied mask file

end subroutine circle_mask_make

!-- @brief Creates a closed polygon mask.
!-- @details Reads a series of vertex coordinates (longitude, latitude) defining one or more
!-- closed polygons, along with a refinement degree, from either a Namelist file (.nml)
!-- or a NetCDF file (.nc, .nc4). It then saves this information into a standardized
!-- NetCDF file in the `tmpfile/` directory.
!-- The refinement degree is checked against `max_iter_spc`.
!-- @param inputfile Path to the input file (Namelist or NetCDF) defining the polygon(s).
!-- @param mask_select Type of mask being created ('mask_domain', 'mask_refine', 'mask_patch').
subroutine close_mask_make(inputfile, mask_select)
    ! Converts file type to NetCDF for easier reading; if already NetCDF, assigns directly.
    use netcdf
    use consts_coms, only : r8, pathlen, file_dir, mask_domain_ndm, mask_patch_ndm
    use refine_vars, only : mask_refine_ndm, max_iter_spc
    use MOD_file_preprocess, only : close_Mesh_Save ! Subroutine to save closed polygon data to NetCDF

    implicit none
    character(pathlen), intent(in) :: inputfile   ! Input: Path to the closed polygon definition file (.nml or .nc/.nc4)
    character(*), intent(in)       :: mask_select ! Input: Type of mask ('mask_domain', 'mask_refine', 'mask_patch')
    integer :: ncid, varid                        ! NetCDF file ID, variable ID
    character(pathlen) :: lndname, line           ! Output NetCDF filename, line read from nml
    logical :: fexists                            ! File existence flag (not used in this scope)
    integer :: i, close_num, refine_degree, length, io_stat! Loop counter, number of points in polygon, refinement degree, string length, I/O status
    real(r8), allocatable :: close_points(:,:)     ! Array to store polygon vertex coordinates (lon, lat)
    character(5) :: numc, refinec                ! Character strings for numbering and refinement degree in filename

    ! Check for self-intersection of line segments, then number and sort segments (comment implies functionality not shown in code)
    length = len_trim(inputfile) ! Get length of input filename

    ! Process based on input file type
    if ('nml' == inputfile(length-2:length)) then ! If input is a Namelist file
        ! Namelist format:
        ! First line: number of points in the polygon
        ! Second line: refinement degree
        ! Subsequent lines: lon, lat for each point
        open(unit=10, file=inputfile, status='old', action='read')
        REWIND(10)
        ! Read close_num (number of points in the polygon)
        read(10, '(A)', iostat=io_stat) line  ! Read the whole line
        if (io_stat /= 0) then
            print *, "Error: Failed to read close_num" ! Changed from circle_num
            close(10)
            stop
        end if
        read(line(index(line, '=')+1:), *) close_num  ! Parse close_num

        ! Read refine_degree
        read(10, '(A)', iostat=io_stat) line  ! Read the whole line
        if (io_stat /= 0) then
            print *, "Error: Failed to read close_refine"
            close(10)
            stop
        end if
        read(line(index(line, '=')+1:), *) refine_degree  ! Parse refine_degree
        if (refine_degree > max_iter_spc) then
            write(io6, *) "refine_degree > max_iter_spc :", refine_degree, ">", max_iter_spc
            close(10)
            return ! Do not process further
        end if

        allocate(close_points(close_num, 2)) ! Allocate for polygon vertex coordinates
        do i = 1, close_num, 1
            read(10, *) close_points(i, 1), close_points(i, 2) ! Read lon, lat
        end do
        close(10)

        ! Save to NetCDF and deallocate
        write(refinec, '(I1)') refine_degree
        if (mask_select == 'mask_domain') then
            mask_domain_ndm = mask_domain_ndm + 1
            write(numc, '(I2.2)') mask_domain_ndm
        else if (mask_select == 'mask_refine') then
            mask_refine_ndm(refine_degree) = mask_refine_ndm(refine_degree) + 1
            write(numc, '(I2.2)') mask_refine_ndm(refine_degree)
        else if (mask_select == 'mask_patch') then
            mask_patch_ndm = mask_patch_ndm + 1
            write(numc, '(I2.2)') mask_patch_ndm
        end if
        lndname = trim(file_dir)// 'tmpfile/'//trim(mask_select)// '_close_'//trim(refinec)//'_'//trim(numc)//'.nc4'
        CALL close_Mesh_Save(lndname, close_num, close_points)
        deallocate(close_points)

    else if (('.nc' == inputfile(length-2:length)) .or. &  ! If input is a NetCDF file
             ('nc4' == inputfile(length-2:length))) then
        CALL CHECK(NF90_OPEN(trim(inputfile), nf90_nowrite, ncid)) ! Open NetCDF
        CALL CHECK(NF90_INQ_VARID(ncid, 'close_refine', varid))   ! Variable name for refinement degree
        CALL CHECK(NF90_GET_VAR(ncid, varid, refine_degree))     ! Read refinement degree
        CALL CHECK(NF90_CLOSE(ncid))                             ! Close NetCDF
        if (refine_degree > max_iter_spc) then
            write(io6, *) "refine_degree > max_iter_spc :", refine_degree, ">", max_iter_spc
            return ! Do not process further
        end if

        ! Generate output filename and copy input NetCDF
        write(refinec, '(I1)') refine_degree
        if (mask_select == 'mask_domain') then
            mask_domain_ndm = mask_domain_ndm + 1
            write(numc, '(I2.2)') mask_domain_ndm
        else if (mask_select == 'mask_refine') then
            mask_refine_ndm(refine_degree) = mask_refine_ndm(refine_degree) + 1
            write(numc, '(I2.2)') mask_refine_ndm(refine_degree)
        else if (mask_select == 'mask_patch') then
            mask_patch_ndm = mask_patch_ndm + 1
            write(numc, '(I2.2)') mask_patch_ndm
        end if
        lndname = trim(file_dir)// 'tmpfile/'//trim(mask_select)// '_close_'//trim(refinec)//'_'//trim(numc)//'.nc4'
        CALL execute_command_line('cp '//trim(inputfile)//' '//trim(lndname)) ! Copy to standardized name
    else
        stop "ERROR! must choose 'ncfile' or 'nml' please check!" ! Invalid file type
    end if
    write(io6, *) lndname ! Print name of created/copied mask file

end subroutine close_mask_make

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
    use MOD_file_preprocess, only : Unstructured_Mesh_Save ! Subroutine to handle NetCDF writing

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

!-- @brief Checks the status of a NetCDF operation and stops if an error occurred.
!-- @param STATUS The status code returned by a NetCDF library call.
SUBROUTINE CHECK(STATUS)
    use netcdf ! For NF90_NOERR and NF90_STRERROR
    INTEGER, intent (in) :: STATUS ! Input: Status code from a NetCDF operation
    if  (STATUS .NE. NF90_NOERR) then ! NF90_NOERR (usually 0) indicates no error
        print *, NF90_STRERROR(STATUS) ! Print the error message corresponding to the status code
        stop 'NetCDF operation failed. Stopping.' ! Stop execution
    endif
END SUBROUTINE CHECK

!-- @brief Computes the Voronoi diagram from an existing Delaunay triangulation.
!-- @details This subroutine transforms a Delaunay triangulation (defined by `mem_delaunay` module variables)
!-- into its dual Voronoi diagram. The key idea is that Voronoi cell centers (M points in the new grid)
!-- are the circumcenters of the Delaunay triangles, and Voronoi vertices (W points in the new grid)
!-- are the vertices of the Delaunay triangles.
!-- - It reassigns dimensions: `nma` (Voronoi cells) becomes `nwd` (Delaunay vertices), etc.
!-- - Allocates memory for Voronoi neighbor tables (`itab_v`, `itab_w`, `itab_m`).
!-- - Copies Delaunay vertex coordinates (`xemd, yemd, zemd`) to Voronoi vertex coordinates (`xew, yew, zew`).
!-- - Calculates Voronoi cell center coordinates (`xem, yem, zem`) as barycenters (centroids) of Delaunay triangles.
!-- - Populates the Voronoi neighbor tables based on the Delaunay connectivity.
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

!-- @brief Performs PCVT (Lloyd's algorithm) iterations to optimize Voronoi cell centroids.
!-- @details This subroutine iteratively adjusts the locations of Voronoi cell centers (M points)
!-- to be the circumcenters of their defining Delaunay triangle vertices (W points).
!-- This process helps make the Voronoi cells more regular and centroidal.
!-- For global domains, it involves transformations between spherical and polar stereographic (PS) plane coordinates.
!-- - Loops over all M points (Voronoi cell centers).
!-- - For each M point, identifies its 3 surrounding W points (vertices of the original Delaunay triangle
!--   whose circumcenter would ideally be the M point, or vertices of the Voronoi cell in the dual).
!-- - Transforms these W points to a PS plane tangent at the current M point's barycenter (`xebc, yebc, zebc`).
!-- - Calculates the circumcenter (`xcc, ycc`) of these 3 W points in the PS plane.
!-- - Transforms this PS circumcenter back to Earth coordinates (`dxe, dye, dze`) relative to the barycenter.
!-- - Updates the M point's location to this new circumcenter.
!-- - Projects the updated M point coordinates back onto the Earth's surface (normalizes to `erad`).
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

!-- @brief Computes various geometric properties of the hexagonal (Voronoi) grid cells.
!-- @details This subroutine calculates latitudes, longitudes, distances, areas, and normal vectors
!-- for the M (cell center), V (edge midpoint), and W (vertex) points of the Voronoi grid.
!-- It handles transformations to polar stereographic planes for gradient calculations and
!-- computes coefficients for converting velocities between Earth-Cartesian and cell-relative components.
!-- Key calculations include:
!-- - Lat/lon of M, V, W points.
!-- - `dnu`, `dnv`: Normal distances across U (M-M) and V (W-W) faces/edges.
!-- - `unx, vnx`: Unit normal vector components for U and V faces.
!-- - `arm0`, `arw0`: Areas associated with M and W points (Voronoi cell area and related kite areas).
!-- - `quarter_kite`: Area of quarter-sections of kites formed by M, V, and two W points, used for area distribution.
!-- - Gradient coefficients (`gxps1`, `gyps1`, etc.) on a polar stereographic plane.
!-- - Coefficients (`ecvec_vx`, etc.) for reconstructing Earth-Cartesian velocity components from cell-face normal components.
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
  
