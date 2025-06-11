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
    use consts_coms                                                         ! Module containing physical and mathematical constants
    use refine_vars                                                         ! Module containing variables related to mesh refinement
    use MOD_data_preprocessing     , only : data_preprocess                 ! Module for preprocessing input data
    use MOD_grid_preprocessing     , only : grid_preprocess                 ! Module for preprocessing grid data
    use MOD_Area_judge             , only : Area_judge, Area_judge_refine   ! Module for judging areas for refinement/masking
    use MOD_GetContain             , only : Get_Contain                     ! Module to determine if points are contained in areas
    use MOD_GetRef                 , only : GetRef                          ! Module to get refinement flags
    use MOD_Refine                 , only : refine_loop                     ! Module for the main mesh refinement loop
    use MOD_mask_postprocessing    , only : mask_postproc                   ! Module for post-processing masks
    use MOD_namelist               , only : read_nl                         ! Module for reading namelist files
    use MOD_grid_initialization    , only : init_consts, gridinit           ! Module for grid initialization
    use MOD_utilities              , only : mode4mesh_make, CHECK           ! Module for utility functions

    implicit none
    character(pathlen) :: nlfile = 'mkgrd.mnl'                              ! Namelist file name, default 'mkgrd.mnl'
    character(pathlen) :: finfolist                                         ! Path to save the namelist file
    character(pathlen) :: lndname                                           ! Land grid file name
    character(5)       :: stepc                                             ! Character representation of the current step
    character(5)       :: nxpc                                              ! Character representation of NXP (grid resolution parameter)
    logical            :: exit_loop                                         ! Flag to exit the refinement loop
    logical            :: fexists                                           ! Flag to check if a file exists
    integer :: i, ncid, dimID_sjx, sjx_points, length                       ! Index, ncid, dimID_sjx, sjx_points, length

    ! Read namelist file
    CALL getarg(1, nlfile)                                                  ! Get the namelist file name from the command line
    call read_nl(nlfile)                                                    ! Read settings from the namelist file
    
    mesh_type = trim(mesh_type)
    mode_grid = trim(mode_grid)
    ! Check if the mesh type is valid
    if ((mesh_type /= 'landmesh')   .and. &
        (mesh_type /= 'oceanmesh')  .and. &
        (mesh_type /= 'earthmesh')) then
        write(io6, *) "ERROR! mesh_type = ", mesh_type
        STOP "ERROR! mesh_type mush be landmesh/oceanmesh/earthmesh"
    end if

    step = 1                                                                 ! Initialize current refinement step
    num_vertex = 1                                                           ! Initialize number of vertices
    ! Handle mask restart scenario
    if (mask_restart) then
        call init_consts()                                                   ! Initialize constants
        refine = .false.                                                     ! Disable refinement for mask restart
        step = max_iter + 1                                                  ! Set step beyond max_iter to skip refinement loop
        if ((mesh_type == 'oceanmesh') .and. (.not. mask_patch_on)) then
            ! This is for the case of only adjusting mask_sea_ratio, if you want to patch, you need to do it from the later process
            ! TODO: need show the example of mask_patch_on
            !please add description here @RuiZhang
            write(io6, *) "Remask_restart start"
            CALL mask_postproc(mesh_type) 
            write(io6, '(A)') "--------------------------------"
            write(io6, '(A)') ""
            write(io6, '(A)') "!! Successfully Make Grid End !!"
            write(io6, '(A)') ""
            write(io6, '(A)') "--------------------------------"
            stop "Remask_restart finish"                                     ! Stop execution after remasking
        end if
    end if
   
    ! Save a copy of the namelist file
    finfolist = trim(file_dir)//'result/namelist.save'                       ! Define path for saving namelist
    CALL execute_command_line('cp '//trim(nlfile)//' '//trim(finfolist))     ! Copy namelist to result directory

    ! Grid generation logic based on mode_grid type
    inquire(file = mode_file, exist = fexists)                               ! Check if the mode_file (initial grid file) exists
    if ((mode_grid == 'hex') .or. &                                          ! Hexagonal grid
        (mode_grid == 'tri')) then                                           ! Triangular grid
        if (fexists) then                                                    ! If mode_file exists, use it as a base
            inquire(file = mode_file_description, exist = fexists)           ! Check for a description file (currently leads to error if exists)
            if (.not. fexists) then
                write(io6, *) "mode_file is exist!"
                ! Check the consistency of mode_file
                CALL CHECK(NF90_OPEN(trim(mode_file), nf90_nowrite, ncid))
                CALL CHECK(NF90_INQ_DIMID(ncid, "sjx_points", dimID_sjx))    ! Get dimension ID for triangle/polygon points
                CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_sjx, len = sjx_points)) ! Get number of points
                CALL CHECK(NF90_CLOSE(ncid))                                 ! Close the NetCDF file
                if (int((sjx_points-1)/20) /= int(nxp*nxp)) then             ! Validate NXP consistency
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
            call init_consts()                                              ! Initialize constants
        else                                                                ! If mode_file does not exist, generate grid from scratch
            ! Initialize, execute, and end olam run (conceptual steps for grid generation)
            call init_consts()                                              ! Initialize constants
            call gridinit()                                                 ! Initialize the grid structure (e.g., icosahedron)
        end if
        write(io6, *) 'grid preporces start'
        CALL grid_preprocess()                                              ! Preprocess the generated or loaded grid
        write(io6, *) 'grid preporces finish'

    else if ((mode_grid == 'lonlat' ) .or. &                                ! Longitude-Latitude grid
             (mode_grid == 'lambert')) then                                 ! Lambert conformal conic projection grid
        if (fexists) then                                                   ! mode_file must exist for these types
            write(io6, *) 'mode4mesh_make start'
            inquire(file = mode_file, exist = fexists)                      ! Double check existence (already done)
            if (.not. fexists) then
                write(io6, *) "The input file " // trim(mode_file) // " is missing."
                stop "Stopping model run."
            endif
            CALL mode4mesh_make(mode_file, mode_grid)                       ! Create mesh from the mode_file
            write(io6, *) 'mode4mesh_make complete'
            write(io6, *) ""
        else
            write(io6, *) "ERROR! mode_file must fexists when mode_grid as ", mode_grid
            stop
        end if

        if (refine) then                                                    ! If refinement was intended, disable it for lonlat/lambert as it's handled differently or not supported
            write(io6, *) "turn refine off"
            refine = .false.
        else
            write(io6, *) "refine is on"
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
