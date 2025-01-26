!DESCRIPTION
!===========
! This program is the unstructure mesh generation tool for land surface models (e.g.,CoLM).

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
!* :SUBROUTINE:"colm_CaMa_init" :  Initialization of the coupler


!REVISION HISTORY
!----------------
! 2023.02.21  Zhongwang Wei @ SYSU
! 2021.12.02  Zhongwang Wei @ SYSU 
! 2020.10.01  Zhongwang Wei @ SYSU


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
    use consts_coms
    use refine_vars
    use MOD_GetContain       , only : Get_Contain
    use MOD_Refine           , only : refine_loop
    use MOD_data_preprocess  , only : data_preprocess
    use MOD_Area_judge       , only : Area_judge
    use MOD_GetThreshold     , only : GetThreshold
    use MOD_Threshold_Read   , only : Threshold_Read

    implicit none
    character(pathlen) :: nlfile = 'mkgrd.mnl'
    character(LEN=256) :: finfolist, lndname
    character(5) :: modec, stepc, nxpc
    logical :: exit_loop, fexists
    io6 = 6! If run is sequential, default choice is to set io6 to standard output unit 6.
    ! Initialize HDF5 library

    !call h5open_f(hdferr)
    CALL getarg(1, nlfile)
    ! get nml name from command line. For example: if execute ./mkgrd.x ../mkgrd.nml ! Add comment by Rui Zhang
    ! nlfile = ../mkgrd.nml ! Add comment by Rui Zhang
    ! Read Fortran namelist
    call read_nl(nlfile)

    file_dir = trim(base_dir) // trim(expnme) // '/'
    CALL execute_command_line('rm -rf '//trim(file_dir)) ! rm old filedir
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"contain/") ! 为啥这里要有一个这个
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"gridfile/")
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"patchtype/")
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"result/")
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"threshold/")
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"tmpfile/")

    finfolist = trim(file_dir)//'result/namelist.save' !
    CALL execute_command_line('cp '//trim(nlfile)//' '//trim(finfolist)) ! cp ../mkgrd.nml finfolist

    if ((mode == 6) .or. (mode == 3)) then
        print*, 'mode set as ', mode
    else if (mode == 4) then
        print*, 'mode set as ', mode, "mode4_gridtype = ", mode4_gridtype
        if (refine) then
            print*, "turn refine from TRUE to FALSE"
            refine = .false.
        else
            print*, "refine is FALSE"
        end if
        if ((mode4_gridtype /= 'lonlat' ) .and. &
            (mode4_gridtype /= 'lambert') .and. &
            (mode4_gridtype /= 'cubical')) stop "ERROR! mode4_gridtype mismatch" 
    else
        stop 'mode can be only set as 3 or 4 or 6'
    end if

    if (mode == 4) then
        write(io6, *) 'mode4mesh_make start'
        ! if mode4_gridtype == 'lambert' need change the DmArea Range
        CALL mode4mesh_make()
        write(io6, *) 'mode4mesh_make complete'
        print*, ""
    else
        inquire(file = mode_filedir, exist = fexists)
        if (fexists) then
            ! check the file in the mode_filedir
            write(modec, '(I1.1)') mode
            write(nxpc, '(I4.4)') NXP
            write(stepc, '(I2.2)') step
            lndname = trim(file_dir)// 'gridfile/gridfile_NXP' // trim(nxpc) // '_'//trim(stepc)// '.nc4'
            CALL execute_command_line('cp '//trim(mode_filedir)//' '//trim(lndname))
            write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(io6, *) 'grid_write: opening file:', lndname
            write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

            if (refine == .false.) CALL execute_command_line('cp '//trim(lndname)//' &
                '//trim(file_dir)//'result/gridfile_NXP'//trim(nxpc)//'_mode'//trim(modec)//'.nc4')
        else
            ! Initialize, execute, and end olam run

            call init_consts()

            call gridinit()
        end if
    end if

    write(io6, *) 'data preporcess start' 
    CALL data_preprocess() ! Pre set some necessary data such as cellarea and landtype
    write(io6, *) 'data preporcess complete'
    print*, ""

    write(io6, *) 'area judge start'
    CALL Area_judge() ! Determination of DmArea(RfArea)
    write(io6, *) 'area judge complete'
    print*, ""

    if (refine) then
        step = 0 
        if (max_iter <= step) stop 'Error! max_iter must more than zero'
        write(io6, *) 'make grid with refine mesh'
        print*, ""
 
        if (th_file_read) then
            if (max_iter == 1) then
                print*, "read thresholdfile directly without threshold calculate"
            else
                stop "ERROR if th_file_read == .true. max_iter only set as 1"
            end if
        else
            if ((refine_num_landtypes .eqv. .false.) .and. &
                (refine_area_mainland .eqv. .false.) .and. &
                (all(refine_onelayer  .eqv. .false.)).and. &
                (all(refine_twolayer  .eqv. .false.))) then
                stop "Error! MUst one of TRUE in the refine_num_landtypes or &
                      refine_area_mainland or refine_onelayer or refine_twolayer &
                      when refine is TRUE"
            end if
            write(io6, *) 'Threshold_Read start'
            CALL Threshold_Read() ! var, var_m_s,.... 
            write(io6, *) 'Threshold_Read complete'
            print*, ""
        end if

        ! After merge refine_sjx and refine_lbx into refine_ustrgrid
        write(io6, *) 'step =', step
        write(io6, *) 'start do-while'
        write(io6, *) 'max_iter =', max_iter
        exit_loop = .false.
        do while(step < max_iter)
           
            if (th_file_read) then
                print*, "read thresholdfile directly without threshold calculate"
                print*, "jump out Get_Contain()" 
            else
                write(io6, *) 'Get_Contain start'
                ! only calculate for newly-generated tri or polygon
                CALL Get_Contain()
                write(io6, *) 'Get_Contain complete'
            end if
            
            write(io6, *) 'GetThreshold start'
            CALL GetThreshold(exit_loop)
            if (exit_loop) then
                refine = .false.
                print *, 'Exiting loop due to ref_sjx equal to zero! &
                          turn refine from to True to False !'
                exit  ! 退出外部的 DO WHILE 循环
            end if
            write(io6, *) 'GetThreshold complete'

            write(io6, *) 'refine_loop start'
            CALL refine_loop(exit_loop)
            if (exit_loop) then
                refine = .false.
                print *, 'Exiting loop due to ref_sjx equal to zero! &
                          turn refine from to True to False !'
                exit  ! 退出外部的 DO WHILE 循环
            end if
            write(io6, *) 'refine_loop complete'

            step = step + 1
            write(io6, *) 'step=',step
        end do
        write(io6, *) 'finish do-while'
    else
        write(io6, *) 'make grid with basic mesh'
    end if
    ! calculate for newly-generated tri or polygon
    ! calculate for tri or polygon in domain area rather than refine area
    call Get_Contain() ! 不管细化与否都是要获取patchID的
    write(io6, '(A)') "--------------------------------"
    write(io6, '(A)') ""
    write(io6, '(A)') "!! Successfully Make Grid End !!"
    write(io6, '(A)') ""
    write(io6, '(A)') "--------------------------------"

end program main


subroutine read_nl(file)
    use consts_coms
    use refine_vars
    implicit none

    character(*), intent(in) :: file

    logical :: fexists
    namelist /mkgrd/ nl
    namelist /mkrefine/ rl
    ! OPEN THE NAMELIST FILE
    inquire(file = file, exist = fexists)
    print*, file
    if (.not. fexists) then
        write(*, *) "The namelist file " // trim(file) // " is missing."
        stop "Stopping model run."
    endif
    open(10, status = 'OLD', file = file)
    ! READ GRID POINT, MODEL OPTIONS, AND PLOTTING INFORMATION FROM THE NAMELIST
    REWIND(10)
    !print*,nl
    read(10, nml = mkgrd)
    close(10)
    write(*, nml = mkgrd)

    !----------------------------------------------------------
    ! Variables in the following section either must not be changed on a history
    ! restart or changing them would be irrelevant.  Thus, they are only copied
    ! from the namelist if a history file is not being read.
    expnme               = nl%expnme
    nxp                  = nl%nxp
    openmp               = nl%openmp
    refine               = nl%refine
    mode                 = nl%mode
    base_dir             = nl%base_dir
    source_dir           = nl%source_dir
    lcs                  = nl%lcs
    mode_filedir         = nl%mode_filedir
    if (mode == 4) then
        ndm_domain       = 1
        mode4_gridtype   = nl%mode4_gridtype
        mode4_datatype   = nl%mode4_datatype
    else ! only use for mode = 3 or 6
        ndm_domain           = nl%ndm_domain
        edgee(1:ndm_domain)  = nl%edgee(1:ndm_domain)
        edgew(1:ndm_domain)  = nl%edgew(1:ndm_domain)
        edges(1:ndm_domain)  = nl%edges(1:ndm_domain)
        edgen(1:ndm_domain)  = nl%edgen(1:ndm_domain)
    end if
    if (refine) then
        open(10, status = 'OLD', file = file)
        ! READ GRID POINT, MODEL OPTIONS, AND PLOTTING INFORMATION FROM THE NAMELIST
        REWIND(10)
        !print*,nl
        read(10, nml = mkrefine)
        close(10)
        write(*, nml = mkrefine)

        ndm_refine            = RL%ndm_refine
        edgee_rf(1:ndm_refine)= RL%edgee_rf(1:ndm_refine)
        edgew_rf(1:ndm_refine)= RL%edgew_rf(1:ndm_refine)
        edges_rf(1:ndm_refine)= RL%edges_rf(1:ndm_refine)
        edgen_rf(1:ndm_refine)= RL%edgen_rf(1:ndm_refine)

        max_iter              = rl%max_iter
        max_sa_iter           = rl%max_sa_iter
        th_file_read          = rl%th_file_read
        if (th_file_read) then
            th_filedir        = rl%th_filedir
            return
        end if
        th_num_landtypes      = rl%th_num_landtypes

        refine_num_landtypes  =  rl%refine_num_landtypes
        refine_area_mainland  =  rl%refine_area_mainland
        refine_onelayer( 1)   =  rl%refine_lai_m
        refine_onelayer( 2)   =  rl%refine_lai_s
        refine_onelayer( 3)   =  rl%refine_slope_m
        refine_onelayer( 4)   =  rl%refine_slope_s
        refine_twolayer( 1)   =  rl%refine_k_s_m
        refine_twolayer( 2)   =  rl%refine_k_s_s
        refine_twolayer( 3)   =  rl%refine_k_solids_m
        refine_twolayer( 4)   =  rl%refine_k_solids_s
        refine_twolayer( 5)   =  rl%refine_tkdry_m
        refine_twolayer( 6)   =  rl%refine_tkdry_s
        refine_twolayer( 7)   =  rl%refine_tksatf_m
        refine_twolayer( 9)   =  rl%refine_tksatf_s
        refine_twolayer( 9)   =  rl%refine_tksatu_m
        refine_twolayer(10)   =  rl%refine_tksatu_s

        th_num_landtypes      = rl%th_num_landtypes
        th_area_mainland      = rl%th_area_mainland
        th_onelayer( 1)       = rl%th_lai_m
        th_onelayer( 2)       = rl%th_lai_s
        th_onelayer( 3)       = rl%th_slope_m
        th_onelayer( 4)       = rl%th_slope_s
        th_twolayer( 1, 1:2)  = rl%th_k_s_m
        th_twolayer( 2, 1:2)  = rl%th_k_s_s
        th_twolayer( 3, 1:2)  = rl%th_k_solids_m
        th_twolayer( 4, 1:2)  = rl%th_k_solids_s
        th_twolayer( 5, 1:2)  = rl%th_tkdry_m
        th_twolayer( 6, 1:2)  = rl%th_tkdry_s
        th_twolayer( 7, 1:2)  = rl%th_tksatf_m
        th_twolayer( 8, 1:2)  = rl%th_tksatf_s
        th_twolayer( 9, 1:2)  = rl%th_tksatu_m
        th_twolayer(10, 1:2)  = rl%th_tksatu_s

        ! onelayer
        if ((refine_onelayer( 1) .eqv. .true.) .and. (th_onelayer( 1) == 999.) ) then
            stop "stop for mismatch between refine_onelayer( 1) and   th_onelayer( 1) "
        end if
        if ((refine_onelayer( 2) .eqv. .true.) .and. (th_onelayer( 2) == 999.) ) stop "stop for &
            mismatch between refine_onelayer( 2) and   th_onelayer( 2) "
        if ((refine_onelayer( 3) .eqv. .true.) .and. (th_onelayer( 3) == 999.) ) stop "stop for &
            mismatch between refine_onelayer( 3) and   th_onelayer( 3) "
        if ((refine_onelayer( 4) .eqv. .true.) .and. (th_onelayer( 4) == 999.) ) stop "stop for &
            mismatch between refine_onelayer( 4) and   th_onelayer( 4) "

        ! twolayer
        if ((refine_twolayer( 1) .eqv. .true.) .and. any(th_twolayer( 1, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer( 1)  and      th_twolayer( 1, 1:2) "
        if ((refine_twolayer( 2) .eqv. .true.) .and. any(th_twolayer( 2, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer( 2)  and      th_twolayer( 2, 1:2) "
        if ((refine_twolayer( 3) .eqv. .true.) .and. any(th_twolayer( 3, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer( 3)  and      th_twolayer( 3, 1:2) "
        if ((refine_twolayer( 4) .eqv. .true.) .and. any(th_twolayer( 4, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer( 4)  and      th_twolayer( 4, 1:2) "
        if ((refine_twolayer( 5) .eqv. .true.) .and. any(th_twolayer( 5, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer( 5)  and      th_twolayer( 5, 1:2) "
        if ((refine_twolayer( 6) .eqv. .true.) .and. any(th_twolayer( 6, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer( 6)  and      th_twolayer( 6, 1:2) "
        if ((refine_twolayer( 7) .eqv. .true.) .and. any(th_twolayer( 7, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer( 7)  and      th_twolayer( 7, 1:2) "
        if ((refine_twolayer( 8) .eqv. .true.) .and. any(th_twolayer( 8, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer( 8)  and      th_twolayer( 8, 1:2) "
        if ((refine_twolayer( 9) .eqv. .true.) .and. any(th_twolayer( 9, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer( 9)  and      th_twolayer( 9, 1:2) "
        if ((refine_twolayer(10) .eqv. .true.) .and. any(th_twolayer(10, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer(10)  and      th_twolayer(10, 1:2) "

    end if

end subroutine read_nl

subroutine mode4mesh_make()

    use consts_coms, only : r8, pathlen, io6, file_dir, EXPNME, NXP, mode, refine, mode4_gridtype, mode4_datatype, mode_filedir, edgee, edgew, edges, edgen
    use netcdf
    use lonlatmesh_coms, only : mesh
    USE refine_vars, only : step
    use MOD_file_preprocess, only : Mode4_Mesh_Save ! Add by Rui Zhang

    implicit none
    logical :: fexists
    integer :: i, j, idx
    integer :: nlon_cal, nlat_cal
    integer :: ncid, dimID_lon, dimID_lat
    integer :: lon_points, lat_points, bound_points, mode4_points
    integer, dimension(10) :: varid
    real(r8) :: lon_start, lat_start, lon_end, lat_end, lon_grid_interval, lat_grid_interval
    real(r8), dimension(:),  allocatable :: lon_center, lat_center, lon_bound, lat_bound
    real(r8), dimension(:, :), allocatable :: lon_vert, lat_vert, lonlat_bound
    integer,  dimension(:, :), allocatable :: ngr_bound
    character(pathlen) :: lndname
    character(5) :: nxpc, stepc
    character(10) :: mesh_type
    
    lndname = trim(mode_filedir)
    print*, lndname
    if (mode4_gridtype == 'lonlat') then
        if (mode4_datatype == 'ncfile') then
            !  如果从外部读入信息，要求读入的是网格的中心点
            CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
            CALL CHECK(NF90_INQ_DIMID(ncid, "lon", dimID_lon))
            CALL CHECK(NF90_INQ_DIMID(ncid, "lat", dimID_lat))
            CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lon, len = lon_points))
            CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lat, len = lat_points))
            allocate(lon_center(lon_points))
            allocate(lat_center(lat_points))
            CALL CHECK(NF90_INQ_VARID(ncid, "long", varid(1)))
            CALL CHECK(NF90_INQ_VARID(ncid, "lati", varid(2)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), lon_center))
            CALL CHECK(NF90_GET_VAR(ncid, varid(2), lat_center))
            CALL CHECK(NF90_CLOSE(ncid))
            
            lon_grid_interval = lon_center(2) - lon_center(1)
            lat_grid_interval = lat_center(2) - lat_center(1)
            allocate(lon_bound(lon_points + 1))
            allocate(lat_bound(lat_points + 1))
            lon_bound(1 : lon_points) = lon_center             - lon_grid_interval / 2.0
            lon_bound(1 + lon_points) = lon_center(lon_points) + lon_grid_interval / 2.0
            where(lon_center > 180.) lon_center = lon_center - 360. ! between -180. and 180.
            where(lon_bound  > 180.) lon_bound  = lon_bound  - 360. ! between -180. and 180.
            lat_bound(1 : lat_points) = lat_center             - lat_grid_interval / 2.0
            lat_bound(1 + lat_points) = lat_center(lat_points) + lat_grid_interval / 2.0
        else if (mode4_datatype == 'namelist') then
            namelist /lonlatmesh/ mesh
            ! OPEN THE NAMELIST FILE
            inquire(file = lndname, exist = fexists)
            if (.not. fexists) then
                write(*, *) "The namelist file " // trim(lndname) // " is missing."
                stop "Stopping model run."
            endif
            open(10, status = 'OLD', file = lndname)
            ! READ GRID POINT, MODEL OPTIONS, AND PLOTTING INFORMATION FROM THE NAMELIST
            REWIND(10)
            !print*,nl
            read(10, nml = lonlatmesh)
            close(10)
            write(*, nml = lonlatmesh)

            mesh_type              = mesh%mesh_type
            lon_start              = mesh%lon_start
            lon_end                = mesh%lon_end
            lon_grid_interval      = mesh%lon_grid_interval
            lon_points             = mesh%lon_points
            lat_start              = mesh%lat_start
            lat_end                = mesh%lat_end
            lat_grid_interval      = mesh%lat_grid_interval
            lat_points             = mesh%lat_points
            ! make sure Data self consistency
            nlon_cal = int( (lon_end - lon_start) / lon_grid_interval) + 1
            nlat_cal = int( (lat_end - lat_start) / lat_grid_interval) + 1
            if (nlon_cal /= lon_points) then
                print*, "nlon_cal = ", nlon_cal
                print*, "lon_points = ", lon_points
                stop 'ERROR nlon_cal /= lon_points'
            else if (nlat_cal /= lat_points) then
                print*, "nlat_cal = ", nlat_cal
                print*, "lat_points = ", lat_points
                stop 'ERROR nlat_cal /= lat_points'
            end if
            
            if (mesh_type == 'center') then
                allocate(lon_center(lon_points))
                allocate(lat_center(lat_points))
                allocate(lon_bound(lon_points + 1))
                allocate(lat_bound(lat_points + 1))
                lon_center = lon_start + lon_grid_interval * ([1:lon_points] - 1)
                lat_center = lat_start + lat_grid_interval * ([1:lat_points] - 1)
                lon_bound(1 : lon_points) = lon_center - lon_grid_interval / 2.0
                lon_bound(1 + lon_points) = lon_end    + lon_grid_interval / 2.0
                lat_bound(1 : lat_points) = lat_center - lat_grid_interval / 2.0
                lat_bound(1 + lat_points) = lat_end    + lat_grid_interval / 2.0
            else if (mesh_type == 'bound') then
                allocate(lon_bound(lon_points))
                allocate(lat_bound(lat_points))
                ! 因为读入是边界值，网格数量比边界个数少一个
                lon_points = lon_points - 1
                lat_points = lat_points - 1
                allocate(lon_center(lon_points))
                allocate(lat_center(lat_points))
                lon_bound(1 : lon_points) = lon_start + lon_grid_interval * ([1:lon_points] - 1)
                lon_bound(1 + lon_points) = lon_end
                lat_bound(1 : lat_points) = lat_start + lat_grid_interval * ([1:lat_points] - 1)
                lat_bound(1 + lat_points) = lat_end
                lon_center = (lon_bound(1 : lon_points) + lon_bound(2 : 1 + lon_points)) / 2.0
                lat_center = (lat_bound(1 : lat_points) + lat_bound(2 : 1 + lat_points)) / 2.0
            else
                stop "error mesh_type mismatch"
            end if
        else
            stop "error mode4_datatype must choose 'ncfile' or 'namelist' please check!"
        end if

        write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        print*, "lonlat bound range"
        print*, "lon_bound(1) = ", lon_bound(1)
        print*, "lon_bound(1 + lon_points) = ", lon_bound(1 + lon_points)
        print*, "lat_bound(1) = ", lat_bound(1)
        print*, "lat_bound(1 + lat_points) = ", lat_bound(1 + lat_points)
        write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        print*, ""
        
        write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        print*, "lon_points = ", lon_points
        print*, "lat_points = ", lat_points
        print*, "lonlat start end interval"
        print*, "lon_center start from :", lon_center(1)
        print*, "lon_center start end  :", lon_center(lon_points)
        print*, "lon_grid_interval     :", lon_grid_interval
        print*, "lat_center start from :", lat_center(1)
        print*, "lat_center start end  :", lat_center(lat_points)
        print*, "lat_grid_interval     :", lat_grid_interval
        write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        print*, ""
        where(lon_center > 180.)  lon_center = lon_center - 360. ! between -180. and 180.
        where(lon_bound  > 180.)  lon_bound  = lon_bound  - 360. ! between -180. and 180.
        NXP = lon_points / 5 ! compare with NXP in triangle or polygon   

        ! turn lon_bound/lat_bound(1D) to lon_vert.lat_vert(2D)
        print*, "turn lon_bound/lat_bound(1D) to lon_vert/lat_vert(2D) start"
        allocate(lon_vert(lon_points + 1, lat_points + 1)); lon_vert = 0.
        allocate(lat_vert(lon_points + 1, lat_points + 1)); lat_vert = 0.
        do j = 1, lat_points + 1, 1
            lon_vert(:, j) = lon_bound
        end do 
        do i = 1, lon_points + 1, 1
            lat_vert(i, :) = lat_bound 
        end do
        print*, "turn lon_bound/lat_bound(1D) to lon_vert.lat_vert(2D) finish"
        print*, ""

    else if (mode4_gridtype == 'lambert') then
        if (mode4_datatype == 'ncfile') then
            CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
            CALL CHECK(NF90_INQ_DIMID(ncid, "xi_vert",  dimID_lon))
            CALL CHECK(NF90_INQ_DIMID(ncid, "eta_vert", dimID_lat))
            CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lon, len = lon_points))
            CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lat, len = lat_points))
            allocate(lon_vert(lon_points, lat_points))
            allocate(lat_vert(lon_points, lat_points))
            CALL CHECK(NF90_INQ_VARID(ncid, "lon_vert", varid(1)))
            CALL CHECK(NF90_INQ_VARID(ncid, "lat_vert", varid(2)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), lon_vert))
            CALL CHECK(NF90_GET_VAR(ncid, varid(2), lat_vert))
            CALL CHECK(NF90_CLOSE(ncid))
            ! 因为读入是边界值，网格数量比边界个数少一个
            lon_points = lon_points - 1
            lat_points = lat_points - 1
        end if
    else if (mode4_gridtype == 'cubical') then
        stop "error lambert or cubical can not use now!"
    else
        stop "error mode4_gridtype must choose 'lonlat' or 'lambert' or 'cubical' please check!"  
    end if

    ! turn 1D to 2D
    print*, "lonlat_bound and ngr_bound calculate start"
    bound_points = (lon_points + 1) * (lat_points + 1) + 1
    mode4_points = lon_points * lat_points + 1
    allocate(lonlat_bound(bound_points, 2)); lonlat_bound = 0.
    allocate(ngr_bound(4, mode4_points)); ngr_bound = 1

    idx = 1
    do j = 1, lat_points + 1, 1
        do i = 1, lon_points + 1, 1
            idx = idx + 1
            lonlat_bound(idx, :) = [lon_vert(i, j), lat_vert(i, j)]
        end do
    end do
    print*, "lonlat_bound calculate finish"
    print*, ""

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
    ngr_bound = ngr_bound + 1 ! never start from 1 !
    print*, "lonlat_bound and ngr_bound calculate finish"

    ! set values for DmArea Range
    print*, "set values for DmArea Range start "
    edgew(1) = max(minval(lonlat_bound(:, 1)), -180.) 
    edgee(1) = min(maxval(lonlat_bound(:, 1)),  180.)
    edgen(1) = min(maxval(lonlat_bound(:, 2)),   90.)
    edges(1) = max(minval(lonlat_bound(:, 2)),  -90.)
    if (lon_points * lon_grid_interval >= 360.) then
        edgew(1) = -180.
        edgee(1) =  180.
    end if
    print*, "edgew = ", edgew(1)
    print*, "edgee = ", edgee(1)
    print*, "edgen = ", edgen(1)
    print*, "edges = ", edges(1)
    print*, "set values for DmArea Range finish"
    print*, ""

    write(nxpc, '(I4.4)') NXP
    write(stepc, '(I2.2)') step
    lndname = trim(file_dir)// 'gridfile/gridfile_NXP' // trim(nxpc) // '_'//trim(stepc)// '.nc4'
    write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(io6, *) 'grid_write: opening file:', lndname
    write(io6, *) 'mode4_points : ', mode4_points
    write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    ! initial file without any refine
    CALL Mode4_Mesh_Save(lndname, bound_points, mode4_points, lonlat_bound, ngr_bound)

    if (refine == .false.) call execute_command_line('cp '//trim(lndname)//' &
        '//trim(file_dir)//'result/gridfile_NXP'//trim(nxpc)//'_mode4.nc4')

    print*, "mode4mesh_make finish"

end subroutine mode4mesh_make

subroutine init_consts()
    use consts_coms
    implicit none


    ! Standard (Earth) values

    erad = 6371.22e3          ! Earth radius [km]
    
    ! Secondary values
    erad2 = erad * 2.
    erad4 = erad * 4.
    eradsq = erad * erad
    erador5 = erad / sqrt(5.)
    eradi = 1.0 / erad
    erad2sq = erad2 * erad2
    dlat = erad * pio180

end subroutine init_consts


!===============================================================================
subroutine gridinit()

    use consts_coms, only : io6, nxp

    use mem_ijtabs, only : fill_jtabs

    use mem_delaunay, only : copy_tri_grid, copyback_tri_grid, nmd, nud, nwd

    use mem_grid, only : nma, nua, nva, nwa, alloc_grid1, alloc_grid2

    implicit none


    ! Horizontal grid setup


    ! Now generate global atmospheric grid

    write(io6, '(/,a)') 'gridinit calling icosahedron'
    call icosahedron(nxp)  ! global spherical domain; calls 2 allocs
    write(io6, '(/,a)') 'gridinit after icosahedron'
    write(io6, '(a,i0)')    ' nmd = ', nmd
    write(io6, '(a,i0)')    ' nud = ', nud
    write(io6, '(a,i0)')    ' nwd = ', nwd


    ! Store a temporaty copy of the full Delaunay mesh
    ! to be used later to construct the surface grid
    call copy_tri_grid()

    call voronoi()
    call pcvt()
    write(io6, '(/,a)') 'gridinit after voronoi'
    write(io6, '(a,i8)')   ' nma = ', nma
    write(io6, '(a,i8)')   ' nua = ', nua
    write(io6, '(a,i8)')   ' nwa = ', nwa

    ! Allocate remaining GRID FOOTPRINT arrays for full domain

    write(io6, '(/,a)') 'gridinit calling alloc_grid1 for full domain'

    call alloc_grid1(nma, nva, nwa)

    ! Initialize dtlm, dtsm, ndtrat, and nacoust,
    ! and compute the timestep schedule for all grid operations.

    write(io6, '(/,a)') 'gridinit calling fill_jtabs'
    print*, nma, nva, nwa

    call fill_jtabs(nma, nva, nwa, 0)

    ! Fill remaining GRID FOOTPRINT geometry for full domain

    write(io6, '(/,a)') 'gridinit calling grid_geometry'

    call grid_geometry_hex()

    call copyback_tri_grid()

    ! Allocate remaining unstructured grid geometry arrays

    write(io6, '(/,a)') 'gridinit calling alloc_grid2'

    call alloc_grid2(nma, nva, nwa)

    ! Write GRIDFILE and SFCGRILE

    write(io6, '(/,a)') 'gridinit calling gridfile_write'
    call gridfile_write()

    write(io6, '(/,a)') 'gridinit completed'

end subroutine gridinit


!===============================================================================

subroutine gridfile_write()
    ! do not calcullate dismm and disww because no use 
    use netcdf
    USE refine_vars, only: step
    use consts_coms, only : r8, pathlen, io6, file_dir, EXPNME, NXP, mode, refine
    use mem_ijtabs, only : mloops, itab_m, itab_w
    use mem_grid, only : nma, nwa, glatw, glonw, glatm, glonm
    use MOD_file_preprocess, only : Unstructured_Mesh_Save ! Add by Rui Zhang
    implicit none

    ! This routine writes the grid variables to the gridfile.
    integer :: im, iw, sjx_points, lbx_points
    real(r8), allocatable :: mp(:,:),wp(:,:)
    integer, allocatable :: ngrmw(:,:),ngrwm(:,:), n_ngrwm(:)
    character(pathlen) :: lndname
    character(5) :: nxpc,stepc

    allocate (ngrmw(3, nma)) ; ngrmw = 0
    allocate (ngrwm(7, nwa)) ; ngrwm = 0
    do im = 1, nma
       ngrmw(1:3, im) = itab_m(im)%iw(1:3)
    enddo
    do iw = 1, nwa
       ngrwm(1:7, iw) = itab_w(iw)%im(1:7)
    enddo

    allocate(mp(nma,2)); mp(:,1) = GLONM; mp(:,2) = GLATM
    allocate(wp(nwa,2)); wp(:,1) = GLONW; wp(:,2) = GLATW

    ! Execute the command
    write(nxpc, '(I4.4)')NXP
    write(stepc, '(I2.2)') step
    ! Open gridfile
    lndname = trim(file_dir)// 'gridfile/gridfile_NXP' // trim(nxpc) // '_'//trim(stepc)// '.nc4'

    write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(io6, *) 'grid_write: opening file:', lndname
    write(io6, *) 'sjx_points:', nma
    write(io6, *) 'lbx_points:', nwa
    write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    sjx_points = nma
    lbx_points = nwa
    ! initial file without any refine
    CALL Unstructured_Mesh_Save(lndname, sjx_points, lbx_points, mp, wp, ngrmw, ngrwm)
    
    deallocate(ngrwm, ngrmw, mp, wp)

    if (refine == .false.) then
       if (mode == 3) then
          call execute_command_line('cp '//trim(lndname)//' '//trim(file_dir)//'result/gridfile_NXP'//trim(nxpc)//'_sjx.nc4')
       else if(mode == 6)then
          call execute_command_line('cp '//trim(lndname)//' '//trim(file_dir)//'result/gridfile_NXP'//trim(nxpc)//'_lbx.nc4')
       end if
    end if

END SUBROUTINE gridfile_write

SUBROUTINE CHECK(STATUS)
    use netcdf
    INTEGER, intent (in) :: STATUS
    if  (STATUS .NE. NF90_NOERR) then ! nf_noerr=0 表示没有错误
        print *, NF90_STRERROR(STATUS)
        stop 'stopped'
    endif
END SUBROUTINE CHECK

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

!===============================================================================

SUBROUTINE pcvt()

    ! Iterative procedure for defining centroidal voronoi cells

    use mem_ijtabs, only : itab_m
    use mem_grid, only : nma, xem, yem, zem, xew, yew, zew
    use consts_coms, only : erad, eradi

    implicit none

    integer :: im, iw1, iw2, iw3

    real :: raxis, raxisi, expansion
    real :: sinwlat, coswlat, sinwlon, coswlon
    real :: dxe, dye, dze
    real :: xebc, yebc, zebc
    real :: x1, x2, x3, y1, y2, y3
    real :: dx12, dx13, dx23
    real :: s1, s2, s3
    real :: xcc, ycc

    ! Compute XEM,YEM,ZEM location as circumcentric coordinates of 3 W points.
    ! This establishes W cell as voronoi.

    ! Loop over all M points

    !$omp parallel
    !$omp do private(iw1,iw2,iw3,xebc,yebc,zebc,raxis,raxisi, &
    !$omp            sinwlat,coswlat,sinwlon,coswlon,dxe,dye,dze,x1,y1, &
    !$omp            x2,y2,x3,y3,dx12,dx13,dx23,s1,s2,s3,ycc,xcc)
    do im = 2, nma

        ! Indices of 3 W points surrounding M point

        if (any(itab_m(im)%iw(1:3) < 2)) cycle

        iw1 = itab_m(im)%iw(1)
        iw2 = itab_m(im)%iw(2)
        iw3 = itab_m(im)%iw(3)

        ! These were initialized to be the barycenter of each triangle

        xebc = xem(im)
        yebc = yem(im)
        zebc = zem(im)


        ! For global domain, transform from sphere to PS plane

        raxis = sqrt(xebc ** 2 + yebc ** 2)
        raxisi = 1.0 / raxis

        sinwlat = zebc * eradi
        coswlat = raxis * eradi

        sinwlon = yebc * raxisi
        coswlon = xebc * raxisi

        ! Transform 3 W points to PS coordinates

        dxe = xew(iw1) - xebc
        dye = yew(iw1) - yebc
        dze = zew(iw1) - zebc
        call de_ps(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, x1, y1)

        dxe = xew(iw2) - xebc
        dye = yew(iw2) - yebc
        dze = zew(iw2) - zebc
        call de_ps(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, x2, y2)

        dxe = xew(iw3) - xebc
        dye = yew(iw3) - yebc
        dze = zew(iw3) - zebc
        call de_ps(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, x3, y3)

        ! Compute intermediate quanties

        dx12 = x2 - x1
        dx13 = x3 - x1
        dx23 = x3 - x2

        s1 = x1**2 + y1**2
        s2 = x2**2 + y2**2
        s3 = x3**2 + y3**2

        ! Algebraic solution for circumcenter Y coordinate

        ycc = .5 * (dx13 * s2 - dx12 * s3 - dx23 * s1) &
                / (dx13 * y2 - dx12 * y3 - dx23 * y1)

        ! Algebraic solution for circumcenter X coordinate

        if (abs(dx12) > abs(dx13)) then
            xcc = (s2 - s1 - ycc * 2. * (y2 - y1)) / (2. * dx12)
        else
            xcc = (s3 - s1 - ycc * 2. * (y3 - y1)) / (2. * dx13)
        endif

        ! For global domain, transform circumcenter from PS to earth coordinates

        call ps_de(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, xcc, ycc)

        xem(im) = dxe + xebc
        yem(im) = dye + yebc
        zem(im) = dze + zebc

    enddo
    !$omp end do nowait

    ! Adjust each M point to the Earth's radius for global domain



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

!===============================================================================

subroutine grid_geometry_hex()

    use mem_ijtabs, only : itab_m, itab_v, itab_w
    use mem_grid, only : nma, nva, nwa, xev, yev, zev, xem, yem, zem, &
            xew, yew, zew, unx, uny, unz, wnx, wny, wnz, &
            vnx, vny, vnz, glonw, glatw, dnu, dniu, dnv, dniv, arw0, &
            arm0, glonm, glatm, glatv, glonv
    use consts_coms, only : erad, piu180, pio2
    use consts_coms, only : r8

    implicit none

    integer :: im, iv, iw
    integer :: im1, im2
    integer :: iw1, iw2
    integer :: j, npoly
    real :: expansion
    real :: raxis
    real :: dvm1, dvm2
    integer :: j1, j2
    integer :: npoly1, npoly2, np
    real :: xv, yv, frac, alpha
    real :: xw1, xw2, yw1, yw2
    integer :: iwp, ivp, imp
    logical :: dops
    real :: quarter_kite(2, nva)
    integer :: lwork
    integer :: info

    real(r8) :: b(7), fo(7), vnx_ps(7), vny_ps(7), vnz_ps(7), vrot_x(7), vrot_y(7)
    real(r8), allocatable :: work(:)
    real(r8), allocatable :: a(:, :)
    real(r8) :: wsize(1), vdotw, vmag, fact

    ! Loop over all M points

    !$omp parallel

    !$omp do private(raxis)
    do im = 2, nma

        ! Latitude and longitude at M points
        raxis = sqrt(xem(im) ** 2 + yem(im) ** 2)  ! dist from earth axis
        glatm(im) = atan2(zem(im), raxis) * piu180
        glonm(im) = atan2(yem(im), xem(im)) * piu180

        ! Fill global index (replaced later if this run is parallel)
        itab_m(im)%imglobe = im

    enddo
    !$omp end do

    ! Loop over all V points

    !$omp do private(im1,im2,iw1,iw2,expansion,raxis,dvm1,dvm2,frac)
    do iv = 2, nva

        ! Fill global index (replaced later if this run is parallel)

        itab_v(iv)%ivglobe = iv

        ! M-point indices of two end points of V segment

        im1 = itab_v(iv)%im(1)
        im2 = itab_v(iv)%im(2)

        ! W-point indices on either side of V segment

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        ! V point is midway between W points of Voronoi cells

        xev(iv) = .5 * (xew(iw1) + xew(iw2))
        yev(iv) = .5 * (yew(iw1) + yew(iw2))
        zev(iv) = .5 * (zew(iw1) + zew(iw2))

        !  push V point coordinates out to earth radius
        expansion = erad / sqrt(xev(iv) ** 2 &
                + yev(iv) ** 2 &
                + zev(iv) ** 2)

        xev(iv) = xev(iv) * expansion
        yev(iv) = yev(iv) * expansion
        zev(iv) = zev(iv) * expansion


        ! Latitude and longitude at V point

        raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis
        glatv(iv) = atan2(zev(iv), raxis) * piu180
        glonv(iv) = atan2(yev(iv), xev(iv)) * piu180


        ! Normal distance across U face
        ! Unit vector components of U face
        !x Convert normal distance across U face to geodesic arc length

        dnu(iv) = sqrt((xem(im1) - xem(im2))**2 &
                + (yem(im1) - yem(im2))**2 &
                + (zem(im1) - zem(im2))**2)

        unx(iv) = (xem(im2) - xem(im1)) / dnu(iv)
        uny(iv) = (yem(im2) - yem(im1)) / dnu(iv)
        unz(iv) = (zem(im2) - zem(im1)) / dnu(iv)

        !x   dnu(iv) = erad2 * asin(dnu(iv) / erad2)
        dniu(iv) = 1. / dnu(iv)

        ! Normal distance across V face
        ! Unit vector components of V face
        !x Convert normal distance across V face to geodesic arc length

        dnv(iv) = sqrt((xew(iw1) - xew(iw2))**2 &
                + (yew(iw1) - yew(iw2))**2 &
                + (zew(iw1) - zew(iw2))**2)

        vnx(iv) = (xew(iw2) - xew(iw1)) / dnv(iv)
        vny(iv) = (yew(iw2) - yew(iw1)) / dnv(iv)
        vnz(iv) = (zew(iw2) - zew(iw1)) / dnv(iv)

        !x   dnv(iv) = erad2 * asin(dnv(iv) / erad2)
        dniv(iv) = 1. / dnv(iv)

        ! Skip this U point if iw1 < 2 or iw2 < 2

        if (iw1 < 2 .or. iw2 < 2) cycle

        itab_v(iv)%mrlv = max(itab_w(iw1)%mrlw, itab_w(iw2)%mrlw)

        ! Compute IM1 and IM2 values of quarter kite area,
        ! and add to ARM0 and ARW0 arrays

        dvm1 = sqrt((xev(iv) - xem(im1))**2 &
                + (yev(iv) - yem(im1))**2 &
                + (zev(iv) - zem(im1))**2)

        dvm2 = sqrt((xev(iv) - xem(im2))**2 &
                + (yev(iv) - yem(im2))**2 &
                + (zev(iv) - zem(im2))**2)

        ! Fractional distance along V edge where intersection with U edge is located

        frac = dvm1 * dniu(iv)

        if (im1 > 1 .and. im2 > 1 .and. (frac < .0001 .or. frac > .9999)) then
            print*, 'Non-intersecting U-V edges detected in grid geometry'
            print*, 'FRAC  = ', frac
            print*, 'IW1 = ', iw1, ' IW2 = ', iw2
            print*, 'IV    = ', iv
            print*, 'GLATV = ', glatv(iv)
            print*, 'GLONV = ', glonv(iv)

            print*, 'dnu(iv),dniu(iv) ', dnu(iv), dniu(iv)

            stop 'STOP U-V edges'
        endif

        quarter_kite(1, iv) = .25 * dvm1 * dnv(iv)
        quarter_kite(2, iv) = .25 * dvm2 * dnv(iv)

    enddo
    !$omp end do
    !$omp end parallel

    !dir$ novector
    do iv = 2, nva
        im1 = itab_v(iv)%im(1)
        im2 = itab_v(iv)%im(2)

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        arm0(im1) = arm0(im1) + 2. * quarter_kite(1, iv)
        arm0(im2) = arm0(im2) + 2. * quarter_kite(2, iv)

        arw0(iw1) = arw0(iw1) + quarter_kite(1, iv) + quarter_kite(2, iv)
        arw0(iw2) = arw0(iw2) + quarter_kite(1, iv) + quarter_kite(2, iv)
    enddo

    ! Lateral boundary copy of arw0

    do iw = 2, nwa
        iwp = itab_w(iw)%iwp
        if (iw /= iwp) arw0(iw) = arw0(iwp)
    enddo

    ! Lateral boundary copy of arm0

    do im = 2, nma
        imp = itab_m(im)%imp
        if (im /= imp) arm0(im) = arm0(imp)
    enddo

    !$omp parallel
    !$omp do private(raxis)
    do iw = 2, nwa

        ! Fill global index (replaced later if this run is parallel)

        itab_w(iw)%iwglobe = iw

        ! Fill outward unit vector components and latitude and longitude of W point

        raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)

        glatw(iw) = atan2(zew(iw), raxis) * piu180
        glonw(iw) = atan2(yew(iw), xew(iw)) * piu180

        wnx(iw) = xew(iw) / erad
        wny(iw) = yew(iw) / erad
        wnz(iw) = zew(iw) / erad
    enddo
    !$omp end do

    !$omp single
    iw = minloc(arw0(2:), 1) + 1
    write(*, *)
    write(*, '(A,f0.4,A)')       " Minimum atmos grid spacing is ", .001 * sqrt(arw0(iw)), " km"
    write(*, '(A,I0,2(A,f0.4))') " at iw = ", iw, ", Lat = ", glatw(iw), ", Lon = ", glonw(iw)

    iw = maxloc(arw0(2:), 1) + 1
    write(*, *)
    write(*, '(A,f0.4,A)')       " Maximum atmos grid spacing is ", .001 * sqrt(arw0(iw)), " km"
    write(*, '(A,I0,2(A,f0.4))') " at iw = ", iw, ", Lat = ", glatw(iw), ", Lon = ", glonw(iw)
    write(*, *)
    !$omp end single

    !$omp do private(npoly,j2,j1,iv,iw1,iw2,dops,npoly1,npoly2,np,&
    !$omp            im1,xw1,xw2,yw1,yw2,xv,yv,alpha)
    do iw = 2, nwa

        ! Number of polygon edges/vertices

        npoly = itab_w(iw)%npoly

        ! Loop over all polygon edges

        do j2 = 1, npoly
            j1 = j2 - 1
            if (j2 == 1) j1 = npoly

            iv = itab_w(iw)%iv(j2)
            iw2 = itab_w(iw)%iw(j2)
            iw1 = itab_w(iw)%iw(j1)

            ! Fractional area of arw0(iw) that is occupied by M and V sectors

            if (itab_v(iv)%iw(1) == iw1) then
                itab_w(iw)%farm(j1) = itab_w(iw)%farm(j1) + quarter_kite(1, iv) / arw0(iw)
                itab_w(iw)%farm(j2) = itab_w(iw)%farm(j2) + quarter_kite(2, iv) / arw0(iw)
            else
                itab_w(iw)%farm(j1) = itab_w(iw)%farm(j1) + quarter_kite(2, iv) / arw0(iw)
                itab_w(iw)%farm(j2) = itab_w(iw)%farm(j2) + quarter_kite(1, iv) / arw0(iw)
            endif
            itab_w(iw)%farv(j2) = (quarter_kite(1, iv) + quarter_kite(2, iv)) / arw0(iw)

            !----------------------------------------
            ! NEW SECTION JULY 2011
            !----------------------------------------

            ! Special - skip gradient calculation if we are at the periodic
            ! domain border and iw1 and iw2 do not share a common vertex

            if (iw == itab_w(iw)%iwp) then
                dops = .true.
            else
                dops = .false.
                npoly1 = itab_w(iw1)%npoly
                npoly2 = itab_w(iw2)%npoly

                do np = 1, npoly1
                    im1 = itab_w(iw1)%im(np)
                    if (im1 == 1) cycle
                    if (any(itab_w(iw2)%im(1:npoly2) == itab_w(iw1)%im(np))) then
                        dops = .true.
                        exit
                    endif
                enddo
            endif

            if (dops) then

                ! Evaluate x,y coordinates of IW1 and IW2 points on polar stereographic plane
                ! tangent at IW

                call e_ps(xew(iw1), yew(iw1), zew(iw1), glatw(iw), glonw(iw), xw1, yw1)
                call e_ps(xew(iw2), yew(iw2), zew(iw2), glatw(iw), glonw(iw), xw2, yw2)
                call e_ps(xev(iv), yev(iv), zev(iv), glatw(iw), glonw(iw), xv, yv)

                ! Coefficients for eastward and northward components of gradient

                itab_w(iw)%gxps1(j1) = yw2 / (xw1 * yw2 - xw2 * yw1)
                itab_w(iw)%gxps2(j1) = -yw1 / (xw1 * yw2 - xw2 * yw1)

                itab_w(iw)%gyps1(j1) = -xw2 / (xw1 * yw2 - xw2 * yw1)
                itab_w(iw)%gyps2(j1) = xw1 / (xw1 * yw2 - xw2 * yw1)

                !----------------------------------------

                if (itab_w(iw)%dirv(j2) < 0.) then
                    alpha = atan2(yw2, xw2)   ! VC(iv) direction counterclockwise from east

                    itab_v(iv)%cosv(1) = cos(alpha)
                    itab_v(iv)%sinv(1) = sin(alpha)

                    itab_v(iv)%dxps(1) = xv
                    itab_v(iv)%dyps(1) = yv
                else
                    alpha = atan2(-yw2, -xw2) ! VC(iv) direction counterclockwise from east

                    itab_v(iv)%cosv(2) = cos(alpha)
                    itab_v(iv)%sinv(2) = sin(alpha)

                    itab_v(iv)%dxps(2) = xv
                    itab_v(iv)%dyps(2) = yv
                endif

            endif

            ! Earth-grid components of rotated polar stereographic easterly
            ! (or cartesian positive x) horizontal unit vector

            itab_w(iw)%unx_w = -sin(glonw(iw))
            itab_w(iw)%uny_w = cos(glonw(iw))



            ! Earth-grid components of rotated polar stereographic northerly
            ! (or cartesian positive y) horizontal unit vector

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

        ! Let's not do this section on the boundary cells

        ivp = itab_v(iv)%ivp
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        ! FARW(1) and FARW(2) interpolation coefficients for ARW and VOLT
        ! (taking V control volume to be full DNU(IV) * DNV(IV) rectangle)

        itab_v(iv)%farw(1) = 2. * (quarter_kite(1, iv) + quarter_kite(2, iv)) / arw0(iw1)
        itab_v(iv)%farw(2) = 2. * (quarter_kite(1, iv) + quarter_kite(2, iv)) / arw0(iw2)

        itab_v(iv)%mrlv = max(itab_w(iw1)%mrlw, itab_w(iw2)%mrlw)

    enddo  ! IV
    !$omp end do



    ! Scale eastward and northward gradient components by farm

    !$omp do private(j)
    do iw = 2, nwa

        ! The gradient components are not computed at the lateral boundaries
        if (iw /= itab_w(iw)%iwp) cycle

        do j = 1, itab_w(iw)%npoly

            itab_w(iw)%gxps1(j) = itab_w(iw)%gxps1(j) * itab_w(iw)%farm(j)
            itab_w(iw)%gyps1(j) = itab_w(iw)%gyps1(j) * itab_w(iw)%farm(j)

            itab_w(iw)%gxps2(j) = itab_w(iw)%gxps2(j) * itab_w(iw)%farm(j)
            itab_w(iw)%gyps2(j) = itab_w(iw)%gyps2(j) * itab_w(iw)%farm(j)

        enddo
    enddo
    !$omp end do

    ! Coefficients for converting earth-cartesian velocity to V and W


    !$omp do private(npoly, fo, a, b, work, info, j, iv, vdotw, vmag, fact, &
    !$omp            vnx_ps, vny_ps, vnz_ps, vrot_x, vrot_y, wsize, lwork)
    do iw = 2, nwa

        npoly = itab_w(iw)%npoly

        ! Default coefficients from Perot
        fo(1:npoly) = 2.0_r8 * itab_w(iw)%farv(1:npoly)

        if (allocated(a)) then
            if (size(a, 2) /= npoly) deallocate(a)
        endif

        if (.not. allocated(a)) allocate(a(3, npoly))

        do j = 1, npoly
            iv = itab_w(iw)%iv(j)

            ! Compute the components of the V unit normals perpendicular to W

            vdotw = vnx(iv) * wnx(iw) + vny(iv) * wny(iw) + vnz(iv) * wnz(iw)

            vnx_ps(j) = vnx(iv) - vdotw * wnx(iw)
            vny_ps(j) = vny(iv) - vdotw * wny(iw)
            vnz_ps(j) = vnz(iv) - vdotw * wnz(iw)

            ! Normalize these new vectors to unit length

            vmag = sqrt(vnx_ps(j)**2 + vny_ps(j)**2 + vnz_ps(j)**2)

            vnx_ps(j) = vnx_ps(j) / vmag
            vny_ps(j) = vny_ps(j) / vmag
            vnz_ps(j) = vnz_ps(j) / vmag

            ! Rotate these new unit normals to a coordinate system with Z aligned with W

            if (wnz(iw) >= 0.0) then

                fact = (wny(iw) * vnx_ps(j) - wnx(iw) * vny_ps(j)) / (1.0 + wnz(iw))

                vrot_x(j) = vnx_ps(j) * wnz(iw) - vnz_ps(j) * wnx(iw) + wny(iw) * fact
                vrot_y(j) = vny_ps(j) * wnz(iw) - vnz_ps(j) * wny(iw) - wnx(iw) * fact

            else

                fact = (wny(iw) * vnx_ps(j) - wnx(iw) * vny_ps(j)) / (1._r8 - wnz(iw))

                vrot_x(j) = -vnx_ps(j) * wnz(iw) + vnz_ps(j) * wnx(iw) + wny(iw) * fact
                vrot_y(j) = -vny_ps(j) * wnz(iw) + vnz_ps(j) * wny(iw) - wnx(iw) * fact

            endif

        enddo

        a(1, 1:npoly) = vrot_x(1:npoly) * vrot_x(1:npoly)
        a(2, 1:npoly) = vrot_y(1:npoly) * vrot_y(1:npoly)
        a(3, 1:npoly) = vrot_x(1:npoly) * vrot_y(1:npoly)

        b(1) = 1._r8 - sum(fo(1:npoly) * a(1, :))
        b(2) = 1._r8 - sum(fo(1:npoly) * a(2, :))
        b(3) = - sum(fo(1:npoly) * a(3, :))

        call dgels('N', 3, npoly, 1, a, 3, b, 7, wsize, -1, info)
        lwork = nint(wsize(1)) + 1

        if (allocated(work)) then
            if (size(work) < lwork) deallocate(work)
        endif
        if (.not. allocated(work)) allocate(work(lwork))

        call dgels('N', 3, npoly, 1, a, 3, b, 7, work, size(work), info)

        ! Vector b is now the correction to the coefficients fo
        b(1:npoly) = b(1:npoly) + fo(1:npoly)

        if (info == 0 .and. all(b(1:npoly) > 0.05_r8) .and. all(b(1:npoly) < 0.7_r8)) then

            fact = sum(b(1:npoly) * vnx_ps(1:npoly) * vnx(itab_w(iw)%iv(1:npoly)))
            fact = (1._r8 - wnx(iw)**2) / max(fact, 1.e-30_r8)

            itab_w(iw)%ecvec_vx(1:npoly) = b(1:npoly) * vnx_ps(1:npoly) * fact

            fact = sum(b(1:npoly) * vny_ps(1:npoly) * vny(itab_w(iw)%iv(1:npoly)))
            fact = (1._r8 - wny(iw)**2) / max(fact, 1.e-30_r8)

            itab_w(iw)%ecvec_vy(1:npoly) = b(1:npoly) * vny_ps(1:npoly) * fact

            fact = sum(b(1:npoly) * vnz_ps(1:npoly) * vnz(itab_w(iw)%iv(1:npoly)))
            fact = (1._r8 - wnz(iw)**2) / max(fact, 1.e-30_r8)

            itab_w(iw)%ecvec_vz(1:npoly) = b(1:npoly) * vnz_ps(1:npoly) * fact

        else

            write(*, *) "Problem optimizing vector coefficients for iw = ", iw
            write(*, *) glatw(iw), glonw(iw)
            write(*, *) info
            write(*, *) real(b (1:npoly))
            write(*, *) real(fo(1:npoly))
            write(*, *) "Using default coefficients."

            itab_w(iw)%ecvec_vx(1:npoly) = fo(1:npoly) * vnx_ps(1:npoly)
            itab_w(iw)%ecvec_vy(1:npoly) = fo(1:npoly) * vny_ps(1:npoly)
            itab_w(iw)%ecvec_vz(1:npoly) = fo(1:npoly) * vnz_ps(1:npoly)

        endif

    enddo
    !omp end do

    if (allocated(a))    deallocate(a)
    if (allocated(work)) deallocate(work)
    !$omp end parallel
end subroutine grid_geometry_hex

