
! This module determines which grid cells of a high-resolution source grid fall within specified domain or refinement areas.
! It defines areas using various geometric shapes (bounding box, Lambert projection outline, circle, closed polygon).
! It also handles land/sea classification of grid cells based on land type data and can modify these classifications using patch files.
! The results (masks indicating inclusion in domain/refinement areas, and land/sea status) are stored and can be saved to NetCDF files.
module MOD_Area_judge
    USE consts_coms
    USE MOD_utilities, only : CHECK
    USE refine_vars, only: refine_setting, mask_refine_cal_type, mask_refine_spc_type, mask_refine_ndm
    USE MOD_utilities, only: bbox_mesh_read,  mode4_mesh_Read, circle_mesh_read, close_mesh_read
    USE MOD_data_preprocessing, only: nlons_Rf_select, nlats_Rf_select, landtypes_global, landtypes, lon_i, lat_i, lon_vertex, lat_vertex, Threshold_Read_Lnd, Threshold_Read_Ocn, Threshold_Read_Earth
    implicit none

   ! Mask indicating land (1) or sea (0) for each cell in the source global grid.
    integer, allocatable, public :: seaorland(:, :)
   ! Mask indicating if a source grid cell is within the primary domain (1) or not (0).
    integer, allocatable, public :: IsInDmArea_grid(:, :)
   ! Mask indicating if a source grid cell is within a specified refinement area (1) or not (0). Used for user-defined refinement regions.
    integer, allocatable, public :: IsInRfArea_grid(:, :)
   ! Mask indicating if a source grid cell is within a calculated refinement area (1) or not (0). Used for threshold-based refinement.
    integer, allocatable, public :: IsInRfArea_cal_grid(:, :)
   ! Number of selected longitudes and latitudes for the primary domain.
    integer, public :: nlons_Dm_select, nlats_Dm_select
   ! Min/max longitude/latitude indices (1-based) defining the bounding box of the primary domain within the source global grid.
    integer, public :: minlon_DmArea, maxlon_DmArea, maxlat_DmArea, minlat_DmArea
   ! Min/max longitude/latitude indices defining the bounding box of the specified refinement area.
    integer, public :: minlon_RfArea, maxlon_RfArea, maxlat_RfArea, minlat_RfArea
   ! Min/max longitude/latitude indices defining the bounding box of the calculated (threshold-based) refinement area.
    integer, public :: minlon_RfArea_cal, maxlon_RfArea_cal, maxlat_RfArea_cal, minlat_RfArea_cal
   ! Min/max longitude/latitude indices defining the bounding box of a patch area.
    integer, public :: minlon_PaArea, maxlon_PaArea, maxlat_PaArea, minlat_PaArea
    contains

   ! Determines the primary domain area and processes land/sea masks and patches.
    !
    !1. If not in `mask_restart` mode:
    !   - Calculates `IsInDmArea_grid` based on `mask_domain_type` (bbox, lambert, circle, close).
    !   - Determines `seaorland` mask based on `landtypes_global` for cells within the domain.
    !   - Saves `IsInDmArea_grid` and related metadata.
    !2. If `mask_patch_on` is true:
    !   - If `mask_restart` is true, reads previously saved domain and sea/land masks.
    !   - Calls `mask_patch_modify` to adjust the `seaorland` mask using patch data.
    !3. If refinement is enabled and is not 'specified' only:
    !   - Calculates `IsInRfArea_cal_grid` for threshold-based refinement based on `mask_refine_cal_type`.
    !   - Validates that the refinement area is within the domain area.
    !   - Reads necessary threshold datasets for land, ocean, or earth system variables.
    !   - Saves `IsInRfArea_cal_grid` and its metadata.
    SUBROUTINE Area_judge()
        IMPLICIT NONE
        integer :: i, j, sum_land_grid, iter        ! Loop counters, sum of land grid cells, iteration variable (refinement degree).
        integer :: numpatch_Dm                      ! Number of patches for the domain mask (optional output from IsInArea_*_Calculation).
        integer, allocatable :: IsInDmArea_select(:, :), seaorland_select(:, :) ! Temporary arrays for reading mask data in restart mode.
        character(16) :: type_select                ! String to select mask type ('mask_domain', 'mask_refine', 'mask_patch').
        character(LEN = 256) :: inputfile, outputfile ! File names for I/O.

        if (.not. mask_restart) then
            write(io6, *) 'IsInArea_grid_Calculation start'
            allocate(IsInDmArea_grid(nlons_source, nlats_source)); IsInDmArea_grid = 0 ! Initialize domain mask to 0 (outside).
            ! Initialize bounding box indices for the domain.
            minlon_DmArea = nlons_source; maxlon_DmArea = 1
            maxlat_DmArea = nlats_source; minlat_DmArea = 1
            iter = 0 ! Iter might represent a refinement level or step; 0 for base domain.
            type_select = 'mask_domain'

            ! Calculate IsInDmArea_grid based on the specified geometric type.
            if (mask_domain_type == 'bbox') then
                CALL IsInArea_bbox_Calculation(type_select, iter, mask_domain_ndm, IsInDmArea_grid, numpatch_Dm) 
            else if (mask_domain_type == 'lambert') then
                CALL IsInArea_lambert_Calculation(type_select, iter, mask_domain_ndm, IsInDmArea_grid, numpatch_Dm)
            else if (mask_domain_type == 'circle') then
                CALL IsInArea_circle_Calculation(type_select, iter, mask_domain_ndm, IsInDmArea_grid, numpatch_Dm)
            else if (mask_domain_type == 'close') then
                CALL IsInArea_close_Calculation(type_select, iter, mask_domain_ndm, IsInDmArea_grid, numpatch_Dm)
            end if
            print*, "numpatch_Dm (pixels in domain mask) = ", numpatch_Dm
            print*, "minlon_DmArea = ", minlon_DmArea
            print*, "maxlon_DmArea = ", maxlon_DmArea
            print*, "maxlat_DmArea = ", maxlat_DmArea
            print*, "minlat_DmArea = ", minlat_DmArea
            write(io6, *) 'IsInArea_grid_Calculation complete'

            ! Determine land/sea classification for cells within the domain.
            write(io6, *) 'sea or land judge start'
            allocate(seaorland(nlons_source, nlats_source)); seaorland = 0 ! 0: sea/no data, 1: land.
            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
            !$OMP PRIVATE(i, j)
            do j = maxlat_DmArea, minlat_DmArea, 1 ! Loop order affects inner loop first (Original comment).
                do i = minlon_DmArea, maxlon_DmArea, 1
                    if (IsInDmArea_grid(i, j) == 1) then ! If cell is within the domain.
                        ! Land/sea judgment for cells within the calculation area.
                        ! Only process land type 0 (water) vs non-0 (land).
                        ! maxlc (max land class) for glaciers or water bodies could be handled differently. (Original comment)
                        if (landtypes_global(i, j) /= 0) seaorland(i, j) = 1 ! 1 indicates land grid cell.
                    end if
                end do
            end do
            !$OMP END PARALLEL DO
	    !TODO: better to use count(seaorland == 1) here?
            !sum_land_grid = sum(seaorland)
	    sum_land_grid = count(seaorland == 1) ! Count land grid cells.
            print*,"Number of land grid cells within the calculation domain: ",sum_land_grid
            nlons_Dm_select = maxlon_DmArea - minlon_DmArea + 1 ! Number of longitudes in selected domain.
            nlats_Dm_select = minlat_DmArea - maxlat_DmArea + 1 ! Number of latitudes in selected domain.
            print*, "IsInDmArea_grid save start"
            outputfile = trim(file_dir) // 'result/IsInDmArea_grid.nc4' ! Final output file.
            CALL IsInArea_grid_Save(outputfile, IsInDmArea_grid, minlon_DmArea, maxlon_DmArea, maxlat_DmArea, minlat_DmArea)
            print*, "IsInDmArea_grid save finish"
            print*, ""
 
        end if

        ! Patching essentially modifies the seaorland mask.
        if (mask_patch_on) then
            write(io6, *) 'make grid with patch mesh in the MOD_Area_judge.F90'
            if (mask_restart) then ! If restarting, read previously saved domain and sea/land data.
                allocate(seaorland(nlons_source, nlats_source)); seaorland = 0
                allocate(IsInDmArea_grid(nlons_source, nlats_source)); IsInDmArea_grid = 0
		!TODO: it is correct of inputfile, original code did not define it!
                inputfile = trim(file_dir) // 'result/IsInDmArea_grid.nc4' ! Path to previously saved data.
                CALL IsInArea_grid_Read(inputfile, IsInDmArea_select, seaorland_select)
                seaorland(minlon_DmArea:maxlon_DmArea, maxlat_DmArea:minlat_DmArea) = seaorland_select(:,:)
                IsInDmArea_grid(minlon_DmArea:maxlon_DmArea, maxlat_DmArea:minlat_DmArea) = IsInDmArea_select(:,:)
		!TODO: check deallocate of all variable
                deallocate(IsInDmArea_select, seaorland_select)
            end if
            type_select = 'mask_patch'
            CALL mask_patch_modify(type_select, iter) ! iter is 0 if not in restart refine context.
        end if
        
        ! If refinement is disabled or only 'specified' refinement is on, return.
        if ((.not. refine) .or. (refine_setting == 'specified')) return

        ! Calculate refinement area for threshold-based refinement ('calculate' or 'mixed').
        allocate(IsInRfArea_cal_grid(nlons_source, nlats_source)); IsInRfArea_cal_grid = 0
        minlon_RfArea = nlons_source; maxlon_RfArea = 1 ! Initialize bounding box for refinement area.
        maxlat_RfArea = nlats_source; minlat_RfArea = 1
        type_select = 'mask_refine'
        iter = 0 ! For calculated refinement mask, iter (degree) is 0.
        if (mask_refine_cal_type == 'bbox') then
            CALL IsInArea_bbox_Calculation(type_select, iter, mask_refine_ndm(iter), IsInRfArea_cal_grid)
        else if (mask_refine_cal_type == 'lambert') then
            CALL IsInArea_lambert_Calculation(type_select, iter, mask_refine_ndm(iter), IsInRfArea_cal_grid)
        else if (mask_refine_cal_type == 'circle') then
            CALL IsInArea_circle_Calculation(type_select, iter, mask_refine_ndm(iter), IsInRfArea_cal_grid)
        else if (mask_refine_cal_type == 'close') then
            CALL IsInArea_close_Calculation(type_select, iter, mask_refine_ndm(iter), IsInRfArea_cal_grid)
        end if
        ! Store these calculated refinement area bounds.
        minlon_RfArea_cal = minlon_RfArea
        maxlon_RfArea_cal = maxlon_RfArea
        maxlat_RfArea_cal = maxlat_RfArea
        minlat_RfArea_cal = minlat_RfArea
        print*, ""
	!TODO:check here, i think it should be  *_cal
        print*, "refine_degree (for calculated mask) = ", iter
        print*, "minlon_RfArea_cal = ", minlon_RfArea_cal
        print*, "maxlon_RfArea_cal = ", maxlon_RfArea_cal
        print*, "maxlat_RfArea_cal = ", maxlat_RfArea_cal
        print*, "minlat_RfArea_cal = ", minlat_RfArea_cal
        print*, ""
        if (minlon_RfArea_cal > maxlon_RfArea_cal) stop "ERROR! minlon_RfArea_cal > maxlon_RfArea_cal"
        if (maxlat_RfArea_cal > minlat_RfArea_cal) stop "ERROR! maxlat_RfArea_cal > minlat_RfArea_cal" 
        
        ! Efficiency of this check? (Original comment)
	!TODO: check here.
        ! Ensure refinement area is within the domain area.
        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
        !$OMP PRIVATE(i, j)
        do j = maxlat_RfArea_cal, minlat_RfArea_cal, 1
            do i = minlon_RfArea_cal, maxlon_RfArea_cal, 1
	    !I think IsInRfArea_cal_grid(i, j) == 1 is correct
                if (IsInRfArea_cal_grid(i, j)== 1) then
                    if (IsInDmArea_grid(i, j) == 0) then ! Ensure refined cells are within the domain.
                        print*,"ERROR!!! the refine area exceed the domain area!!!"
                        stop ! Error! No intersection between refinement area and calculation domain!! (Original comment)
                    end if
                end if
            end do
        end do
        !$OMP END PARALLEL DO
        print*,"the refine area Completely locate at the domain area!!!" ! Refinement area is completely within the calculation domain. (Original comment)
        ! Read threshold data for the calculated refinement area.
        nlons_Rf_select = maxlon_RfArea_cal - minlon_RfArea_cal + 1
        nlats_Rf_select = minlat_RfArea_cal - maxlat_RfArea_cal + 1
        allocate(landtypes(nlons_Rf_select, nlats_Rf_select))
        landtypes = landtypes_global(minlon_RfArea_cal:maxlon_RfArea_cal, maxlat_RfArea_cal:minlat_RfArea_cal) ! Subset landtypes.
        
        if (mesh_type == 'landmesh') then
            CALL Threshold_Read_Lnd(minlon_RfArea_cal, maxlon_RfArea_cal, maxlat_RfArea_cal, minlat_RfArea_cal) 
        else if (mesh_type == 'oceanmesh') then
            CALL Threshold_Read_Ocn(minlon_RfArea_cal, maxlon_RfArea_cal, maxlat_RfArea_cal, minlat_RfArea_cal)
        else if (mesh_type == 'earthmesh') then
            CALL Threshold_Read_Lnd(minlon_RfArea_cal, maxlon_RfArea_cal, maxlat_RfArea_cal, minlat_RfArea_cal)
            CALL Threshold_Read_Ocn(minlon_RfArea_cal, maxlon_RfArea_cal, maxlat_RfArea_cal, minlat_RfArea_cal)
            CALL Threshold_Read_Earth(minlon_RfArea_cal, maxlon_RfArea_cal, maxlat_RfArea_cal, minlat_RfArea_cal) 
        end if

        print*, "IsInRfArea_cal_grid save start"
        outputfile = trim(file_dir) // 'result/IsInRfArea_cal_grid.nc4'
	!!! here is also weird, I do believe it should be *_cal
	!TODO:CHECK HERE 
        CALL IsInArea_grid_Save(outputfile, IsInRfArea_cal_grid, minlon_RfArea_cal, maxlon_RfArea_cal, maxlat_RfArea_cal, minlat_RfArea_cal)
        print*, "IsInRfArea_cal_grid save finish"
        print*, ""

    END SUBROUTINE Area_judge

   ! Determines refinement area for 'specified' refinement at a given iteration/degree.
    ! If `iter` is 0, it uses the pre-calculated `IsInRfArea_cal_grid`.
    !Otherwise, it calculates `IsInRfArea_grid` based on `mask_refine_spc_type` for the
    !specified refinement degree `iter`. It ensures this area is within the primary domain.
    ! iter Input: Current refinement iteration/degree. For iter=0, uses calculated mask.
    SUBROUTINE Area_judge_refine(iter)
        ! For 'specified' refinement, the base map needs to be updated each time. (Original Chinese comment)
        IMPLICIT NONE
        integer, intent(in) :: iter             ! Current refinement degree/iteration.
        integer :: i, j                         ! Loop counters.
        character(16) :: type_select           ! Mask type selector ('mask_refine').
        character(16) :: iterc              ! Character representation of iter for filename.
        character(256) :: outputfile           ! Output filename.

        if (iter == 0) then ! Use pre-calculated threshold-based refinement area.
            IsInRfArea_grid = IsInRfArea_cal_grid
            minlon_RfArea   = minlon_RfArea_cal
            maxlon_RfArea   = maxlon_RfArea_cal
            maxlat_RfArea   = maxlat_RfArea_cal
            minlat_RfArea   = minlat_RfArea_cal
            return
        end if

        ! For specified refinement (iter > 0), calculate IsInRfArea_grid for this degree.
        if (allocated(IsInRfArea_grid)) deallocate(IsInRfArea_grid)
        allocate(IsInRfArea_grid(nlons_source, nlats_source))
        IsInRfArea_grid = 0
        minlon_RfArea = nlons_source; maxlon_RfArea = 1
        maxlat_RfArea = nlats_source; minlat_RfArea = 1

        type_select = 'mask_refine'
        if (mask_refine_spc_type == 'bbox') then
            CALL IsInArea_bbox_Calculation(type_select, iter, mask_refine_ndm(iter), IsInRfArea_grid) 
        else if (mask_refine_spc_type == 'lambert') then
            CALL IsInArea_lambert_Calculation(type_select, iter, mask_refine_ndm(iter), IsInRfArea_grid)
        else if (mask_refine_spc_type == 'circle') then
            CALL IsInArea_circle_Calculation(type_select, iter, mask_refine_ndm(iter), IsInRfArea_grid)
        else if (mask_refine_spc_type == 'close') then
            CALL IsInArea_close_Calculation(type_select, iter, mask_refine_ndm(iter), IsInRfArea_grid)
        end if
        print*, "refine_degree = ", iter
        print*, "minlon_RfArea (specified) = ", minlon_RfArea
        print*, "maxlon_RfArea (specified) = ", maxlon_RfArea
        print*, "maxlat_RfArea (specified) = ", maxlat_RfArea
        print*, "minlat_RfArea (specified) = ", minlat_RfArea
        if (minlon_RfArea > maxlon_RfArea .AND. maxlon_RfArea /=1) stop "ERROR! minlon_RfArea > maxlon_RfArea" ! maxlon_RfArea=1 if no area found
        if (maxlat_RfArea > minlat_RfArea .AND. minlat_RfArea /=1) stop "ERROR! maxlat_RfArea > minlat_RfArea" ! minlat_RfArea=1 if no area found
        
        ! Check if refinement area is within domain.
        !THIS IS also confusing
        if (maxlon_RfArea /=1 .AND. minlat_RfArea /=1) then ! Proceed only if a valid refinement area was found
            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
            !$OMP PRIVATE(i, j)
            do j = maxlat_RfArea, minlat_RfArea, 1
                do i = minlon_RfArea, maxlon_RfArea, 1
                    if (IsInRfArea_grid(i, j) == 1) then
                        if (IsInDmArea_grid(i, j) == 0) then
                            print*,"ERROR!!! the refine area exceed the domain area!!!"
                            stop ! Error! No intersection between refinement area and calculation domain!!
                        end if
                    end if
                end do
            end do
            !$OMP END PARALLEL DO
        else
            print*, "Warning: No specified refinement area found for iteration/degree: ", iter
        endif
        print*,"the refine area Completely locate at the domain area!!!"
        print*, "IsInRfArea_grid_spc save start"
        write(iterc, '(I2.2)') iter
        outputfile = trim(file_dir) // 'result/IsInRfArea_grid_spc_'//trim(iterc)//'.nc4'
        CALL IsInArea_grid_Save(outputfile, IsInRfArea_grid, minlon_RfArea, maxlon_RfArea, maxlat_RfArea, minlat_RfArea)
        print*, "IsInRfArea_grid_spc save finish"
        print*, ""
    END SUBROUTINE Area_judge_refine

   ! Modifies the `seaorland` mask based on patch data.
    ! Calculates `IsInPaArea_grid` based on `mask_patch_type` and then sets
    !corresponding cells in `seaorland` to 0 (sea/water), effectively applying patches
    !that convert land to water.
    ! type_select Input: Type of mask ('mask_patch').
    ! iter Input: Refinement degree/iteration (0 for base patches).
    !@note Only needs to modify `seaorland`, not `IsInDmArea_grid`.
    SUBROUTINE mask_patch_modify(type_select, iter)
        
        IMPLICIT NONE
        character(16), intent(in) :: type_select     ! Should be 'mask_patch'.
        integer, intent(in) :: iter                  ! Iteration/degree for patch definition.
        integer, allocatable :: IsInPaArea_grid(:, :) ! Temporary mask for the current patch area.
        integer :: i, j                             ! Loop counters.

        allocate(IsInPaArea_grid(nlons_source, nlats_source)); IsInPaArea_grid = 0
        ! Initialize patch area bounding box.
        minlon_PaArea = nlons_source; maxlon_PaArea = 1
        maxlat_PaArea = nlats_source; minlat_PaArea = 1

        if (mask_patch_type == 'bbox') then
            CALL IsInArea_bbox_Calculation(type_select, iter, mask_patch_ndm, IsInPaArea_grid) 
        else if (mask_patch_type == 'lambert') then
            CALL IsInArea_lambert_Calculation(type_select, iter, mask_patch_ndm, IsInPaArea_grid)
        else if (mask_patch_type == 'circle') then
            CALL IsInArea_circle_Calculation(type_select, iter, mask_patch_ndm, IsInPaArea_grid)
        else if (mask_patch_type == 'close') then
            CALL IsInArea_close_Calculation(type_select, iter, mask_patch_ndm, IsInPaArea_grid)
        end if

        ! Modify seaorland mask.
        print*, "Modifying seaorland mask with patch start!"
        !TODO: CHECK HERE
        if (maxlon_PaArea /=1 .AND. minlat_PaArea /=1) then ! Proceed only if a valid patch area was found
            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
            !$OMP PRIVATE(i, j)
            do j = maxlat_PaArea, minlat_PaArea, 1
                do i = minlon_PaArea, maxlon_PaArea, 1
                    if (IsInPaArea_grid(i, j) == 1) seaorland(i, j) = 0 ! Convert land pixels to sea/water within the patch.
                end do
            end do
            !$OMP END PARALLEL DO
        endif
        deallocate(IsInPaArea_grid)
        print*, "Modifying seaorland mask with patch finish!"

    END SUBROUTINE mask_patch_modify

   ! Determines if source grid cells are within a bounding box defined area.
    ! Reads bounding box definitions (west, east, north, south) from temporary files
    !(previously created by `mkgrd.F90:bbox_mask_make`). For each box, it updates `IsInArea_grid`
    !for all source grid cells falling within that box. It also updates the overall min/max
    !lat/lon indices covering all processed boxes for the given `type_select`.
    ! type_select Input: Type of mask being processed ('mask_domain', 'mask_refine', 'mask_patch').
    ! iter Input: Refinement degree/iteration, used for filename.
    ! ndm Input: Number of mask definition files for this type and iteration.
    ! IsInArea_grid Input/Output: The main mask grid being updated.
    ! numpatch Optional Output: Total number of source grid cells marked as inside by these boxes.
    SUBROUTINE IsInArea_bbox_Calculation(type_select, iter, ndm, IsInArea_grid, numpatch)
        ! Determines if lat/lon grid cells are within the refinement/domain area.
        implicit none
        character(16), intent(in) :: type_select           ! Type of mask ('mask_domain', 'mask_refine', 'mask_patch').
        integer, intent(in) :: iter, ndm                    ! Refinement degree for filename, number of mask files.
        character(pathlen) :: lndname                       ! Path to mask definition file.
        real(r8) :: edgee_temp, edgew_temp, edges_temp, edgen_temp ! East, West, South, North edges of a bbox.
        integer :: n, i, bbox_num                           ! Loop counters, number of bboxes in a file.
        integer :: temp1, temp2, temp3, temp4               ! Min/max lon/lat indices for a bbox.
        real(r8), allocatable :: bbox_points(:,:)           ! Coordinates of bounding boxes [num_boxes, 4 (W,E,N,S)].
        integer, intent(inout) :: IsInArea_grid(:, :)       ! Output: Grid indicating cells within the area (1=inside).
        integer, intent(out), optional :: numpatch          ! Optional Output: Total count of grid cells marked as inside.
        character(5) :: numc, refinec                       ! Character representation of n and iter for filename.

        if (present(numpatch)) numpatch = 0
        write(refinec, '(I1)') iter ! Format iter (refinement degree) for filename.
        do n = 1, ndm, 1 ! Loop through all mask files for this type/degree.
            write(numc, '(I2.2)') n ! Format file number.
            lndname = trim(file_dir)// 'tmpfile/'//trim(type_select)//'_bbox_'//trim(refinec)//'_'//trim(numc)//'.nc4'
            print*, trim(lndname), ': reading'
            CALL bbox_Mesh_Read(lndname, bbox_num, bbox_points) ! Read bbox definitions.
            do i = 1, bbox_num, 1 ! For each bounding box in the file.
                edgew_temp = bbox_points(i, 1)
                edgee_temp = bbox_points(i, 2)
                edgen_temp = bbox_points(i, 3)
                edges_temp = bbox_points(i, 4)
                ! Determine min/max source grid indices covering this bbox and update global min/max for this mask type.
                CALL minmax_range_make(type_select, edgew_temp, edgee_temp, edgen_temp, edges_temp, temp1, temp2, temp3, temp4)
                ! Mark cells within this specific bbox as inside.
                IsInArea_grid(temp1:temp2,temp3:temp4) = 1
                if (present(numpatch)) numpatch = numpatch + (temp2 - temp1 + 1) * (temp4 - temp3 + 1)
            end do
            deallocate(bbox_points)
        end do

    END SUBROUTINE IsInArea_bbox_Calculation

   ! Determines if source grid cells are within a Lambert projection defined area.
    ! Reads Lambert projection polygon definitions from temporary files. For each polygon,
    !it identifies the source grid cells whose centers fall within that polygon using the
    !`is_point_in_convex_polygon` function and updates `IsInArea_grid`.
    !Handles dateline crossing for polygons.
    SUBROUTINE IsInArea_lambert_Calculation(type_select, iter, ndm, IsInArea_grid, numpatch)
        ! Determines if lat/lon grid cells are within the refinement/domain area defined by Lambert polygons.
        implicit none
        character(16), intent(in) :: type_select       ! Type of mask ('mask_domain', 'mask_refine', 'mask_patch').
        integer, intent(in) :: iter, ndm                ! Refinement degree, number of mask files.
        character(pathlen) :: lndname                   ! Path to mask definition file.
        integer :: i, j, k, ii, num_edges, n            ! Loop counters.
        integer :: temp1, temp2, temp3, temp4           ! Temporary storage for bounding box indices.
        integer :: minlon_source, maxlon_source, maxlat_source, minlat_source ! Source grid indices covering a polygon.
        integer :: bound_points, mode_points            ! Number of boundary (vertex) and mode (cell) points from Lambert file.
        logical :: is_inside_temp_points                ! Flag: true if a point is inside the current polygon.
        real(r8):: temp_points(7, 2)                    ! Coordinates of polygon vertices (max 7 assumed for `is_point_in_convex_polygon`).
        real(r8):: point_i(2), point_c(2)               ! Coordinates of current source grid cell center, temp for dateline check.
        real(r8):: edgee_temp, edgew_temp, edges_temp, edgen_temp ! Polygon bounding box.
        real(r8), dimension(:, :), allocatable :: lonlat_bound   ! Vertex coordinates of the Lambert polygon.
        integer,  dimension(:, :), allocatable :: ngr_bound    ! Connectivity of the Lambert polygon (not directly used for point-in-poly).
        integer,  dimension(:),    allocatable :: n_ngr        ! Number of vertices per cell in Lambert definition.
        integer,  dimension(:),    allocatable :: icl_points   ! Flag per polygon cell if it crosses dateline.
        integer,  dimension(:),    allocatable :: num_points   ! Number of source grid cells inside each Lambert polygon cell.
        integer, dimension(:,:), intent(inout) :: IsInArea_grid ! Output: Grid indicating cells within the area.
        integer, intent(out), optional :: numpatch             ! Optional Output: Total count of grid cells marked inside.
        character(5) :: numc, refinec                          ! Character representation of n and iter.

        num_edges = 4 ! Default number of edges for checking (seems to assume quadrilateral cells from Lambert data).
                      ! This might be overridden by n_ngr(k) later if processing individual cells from mode4_mesh_Read.
        temp_points = 0.0_r8
        write(refinec, '(I1)') iter
        do n = 1, ndm, 1 ! Loop through mask files.
            write(numc, '(I2.2)') n
            lndname = trim(file_dir)// 'tmpfile/'//trim(type_select)//'_lambert_'//trim(refinec)//'_'//trim(numc)//'.nc4'
            CALL Mode4_Mesh_Read(lndname, bound_points, mode_points, lonlat_bound, ngr_bound, n_ngr) ! Read Lambert polygon data.
            
            ! Determine overall bounding box of the Lambert defined area.
            edgew_temp = minval(lonlat_bound(2:bound_points, 1)) ! Skip first point which might be placeholder.
            edgee_temp = maxval(lonlat_bound(2:bound_points, 1))
            edgen_temp = maxval(lonlat_bound(2:bound_points, 2))
            edges_temp = minval(lonlat_bound(2:bound_points, 2))
            print*, "n = ", n, " (Lambert mask file)"
            print*, "Raw BBox: W=", edgew_temp, " E=", edgee_temp, " N=", edgen_temp, " S=", edges_temp

            ! Check for dateline crossing of the entire Lambert shape.
            point_c = 0.0_r8
            do i = 2, bound_points - 1, 1 ! Check segments of the overall boundary.
                if (point_c(1) < abs(lonlat_bound(i+1, 1) - lonlat_bound(i, 1))) then
                    point_c(1) = abs(lonlat_bound(i+1, 1) - lonlat_bound(i, 1)) ! Max longitude diff between adjacent vertices.
                end if
            end do
            point_c(2) = abs(edgee_temp-edgew_temp) ! Overall longitude span.
            if (point_c(1) > point_c(2)) then ! If max segment span > overall span, likely dateline crossing.
                print*, "cross 180! need modified global BBox"
                edgew_temp = -180.0_r8
                edgee_temp =  180.0_r8
                print*, "Adjusted BBox: W=", edgew_temp, " E=", edgee_temp, " N=", edgen_temp, " S=", edges_temp
            else
                print*, "not cross 180! not need modified global BBox"
            end if
            CALL minmax_range_make(type_select, edgew_temp, edgee_temp, edgen_temp, edges_temp) ! Update global bounds for this mask type.

            allocate(icl_points(mode_points)); icl_points = 0 ! Dateline crossing flag for each polygon cell.
            allocate(num_points(mode_points)); num_points = 0 ! Pixels inside each polygon cell.

            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
            !$OMP PRIVATE(i, j, k, ii, temp_points, minlon_source, maxlon_source, maxlat_source, minlat_source, point_i, is_inside_temp_points)
            do k = 2, mode_points, 1 ! Loop through each cell defined in the Lambert file.
                num_edges = n_ngr(k) ! Actual number of vertices for this specific polygon cell.
                if (num_edges > size(temp_points,1)) then
                    print*, "Error: num_edges ", num_edges, " > allocated size ", size(temp_points,1)
                    CYCLE ! Or stop
                endif
                temp_points(1:num_edges, :) = lonlat_bound(ngr_bound(1:num_edges, k), :) ! Get vertices of this polygon cell.
                
                if (maxval(temp_points(1:num_edges, 1)) - minval(temp_points(1:num_edges, 1)) > 180.0_r8) then
                    icl_points(k) = 1 ! Mark as crossing dateline.
                    call CheckCrossing(num_edges, temp_points(1:num_edges, :)) ! Adjust coordinates for consistent point-in-polygon test.
                end if

                ! Determine source grid index range for this polygon cell.
                CALL Source_Find(minval(temp_points(1:num_edges, 1)), lon_vertex, 'lon', minlon_source)
                CALL Source_Find(maxval(temp_points(1:num_edges, 1)), lon_vertex, 'lon', maxlon_source)
                CALL Source_Find(maxval(temp_points(1:num_edges, 2)), lat_vertex, 'lat', maxlat_source)
                CALL Source_Find(minval(temp_points(1:num_edges, 2)), lat_vertex, 'lat', minlat_source)
                minlon_source = max(1, minlon_source - 1) ! Expand search box slightly.
                maxlat_source = max(1, maxlat_source - 1)

                ! Point in polygon test for source grid cells within this bounding box.
                if (icl_points(k) == 1) then ! If polygon crossed dateline.
                    do i = minlon_source, maxlon_source - 1, 1
                        do j = maxlat_source, minlat_source - 1, 1
                            point_i = [lon_i(i), lat_i(j)]
                            is_inside_temp_points = is_point_in_convex_polygon(temp_points(1:num_edges,:), point_i, num_edges)
                            ! Remap longitude index 'ii' if it was adjusted due to dateline crossing.
                            if (i < nlons_source/2 + 1) then 
                                ii = int( i + nlons_source/2)
                            else
                                ii = int( i - nlons_source/2)
                            end if
                            if (is_inside_temp_points) then 
                                IsInArea_grid(ii, j) = 1
                                num_points(k) = num_points(k) + 1
                            end if
                        end do
                    end do
                else ! If polygon does not cross dateline.
                    do i = minlon_source, maxlon_source - 1, 1
                        do j = maxlat_source, minlat_source - 1, 1
                            point_i = [lon_i(i), lat_i(j)]
                            is_inside_temp_points = is_point_in_convex_polygon(temp_points(1:num_edges,:), point_i, num_edges)
                            if (is_inside_temp_points) then
                                IsInArea_grid(i, j) = 1
                                num_points(k) = num_points(k) + 1
                            end if 
                        end do
                    end do
                end if
            end do
            !$OMP END PARALLEL DO
            if (present(numpatch)) numpatch = sum(num_points)
            deallocate(icl_points, num_points, lonlat_bound, ngr_bound, n_ngr)
        end do

    END SUBROUTINE IsInArea_lambert_Calculation

   ! Determines if source grid cells are within circular defined areas.
    ! Reads circle definitions (center, radius) from temporary files. For each circle,
    !it identifies source grid cells whose centers fall within that circle using the
    !`is_point_in_circle` function and updates `IsInArea_grid`.
    !Handles cases where circles might cross poles or the dateline by expanding the search box.
    SUBROUTINE IsInArea_circle_Calculation(type_select, iter, ndm, IsInArea_grid, numpatch)
        ! Requires knowing the coordinates of intersections of the circle with various latitudes. (Original comment - not fully implemented this way, uses point-in-circle.)
        implicit none
        character(16), intent(in) :: type_select       ! Type of mask.
        integer, intent(in) :: iter, ndm                ! Refinement degree, number of mask files.
        character(pathlen) :: lndname                   ! Path to mask definition file.
        integer :: i, j, k, ii, circle_num, n           ! Loop counters.
        integer :: minlon_source, maxlon_source, maxlat_source, minlat_source ! Source grid indices covering a circle's bounding box.
        real(r8) :: point_i(2), point_c(2), radius_c, temp ! Current grid point, circle center, circle radius, temp for radius conversion.
        real(r8) :: edgee_temp, edgew_temp, edges_temp, edgen_temp ! Circle's bounding box.
        real(r8), allocatable :: circle_points(:,:), circle_radius(:) ! Circle definitions [num_circles, 2 (lon,lat)], [num_circles (radius_km)].
        integer,  dimension(:), allocatable :: num_points       ! Number of source grid cells inside each circle.
        integer,  dimension(:, :), allocatable :: num_points_ij  ! Temp for parallel sum, (lon_idx, lat_idx) = 1 if inside.
        integer, dimension(:,:), intent(inout) :: IsInArea_grid   ! Output: Grid indicating cells within circles.
        integer, intent(out), optional :: numpatch               ! Optional Output: Total count of grid cells marked inside.
        character(5) :: numc, refinec                            ! Character representation of n and iter.
        !TODO: do not use exist as a variable
	logical :: is_inside_circle, exist_cross                ! Flag for point-in-circle, flag if bbox crosses pole/dateline.
        
        write(refinec, '(I1)') iter
        do n = 1, ndm, 1 ! Loop through mask files.
            write(numc, '(I2.2)') n
            lndname = trim(file_dir)// 'tmpfile/'//trim(type_select)//'_circle_'//trim(refinec)//'_'//trim(numc)//'.nc4'
            CALL circle_Mesh_Read(lndname, circle_num, circle_points, circle_radius) ! Read circle definitions.
            allocate(num_points(circle_num)); num_points = 0

            print*, "circle_num = ", circle_num
            do k = 1, circle_num, 1 ! For each circle.
                point_c  = circle_points(k, :)
                radius_c = circle_radius(k) ! Radius in km.
                ! Determine bounding box for the circle.
                temp = pio180 * erad / 1000.0_r8 ! Conversion factor: km per degree latitude.
                edgew_temp = point_c(1) - radius_c/(temp*cos(point_c(2)*pio180)) ! Approx West edge.
                edgee_temp = point_c(1) + radius_c/(temp*cos(point_c(2)*pio180)) ! Approx East edge.
                edgen_temp = point_c(2) + radius_c/temp                         ! Approx North edge.
                edges_temp = point_c(2) - radius_c/temp                         ! Approx South edge.
                print*, "Circle center = ", point_c, " radius_km = ", radius_c
                print*, "Initial BBox: W=", edgew_temp, " E=", edgee_temp, " N=", edgen_temp, " S=", edges_temp
                exist_cross = .false.
                ! If circle's bounding box crosses 180 meridian or poles, expand BBox to global/polar cap.
                if ((edgee_temp > 180.0_r8) .or. (edgew_temp < -180.0_r8)) then
                    edgew_temp = -180.0_r8; edgee_temp =  180.0_r8; exist_cross = .true.
                end if
                if (edgen_temp > 90.0_r8 ) then
                    edges_temp = min(edges_temp, edgen_temp) ! Ensure S edge is still valid if N edge went over pole.
                    edgen_temp = 90.0_r8; exist_cross = .true.
                else if (edges_temp < -90.0_r8 ) then
                    edgen_temp = max(edges_temp, edgen_temp)
                    edges_temp = -90.0_r8; exist_cross = .true.
                end if
                if (exist_cross) then
                    print*, "Adjusted BBox due to pole/dateline crossing:"
                    print*, "W=", edgew_temp, " E=", edgee_temp, " N=", edgen_temp, " S=", edges_temp
                end if

                ! Determine source grid index range for this circle's bounding box.
                CALL Source_Find(edgew_temp, lon_vertex, 'lon', minlon_source)
                CALL Source_Find(edgee_temp, lon_vertex, 'lon', maxlon_source)
                CALL Source_Find(edgen_temp, lat_vertex, 'lat', maxlat_source)
                CALL Source_Find(edges_temp, lat_vertex, 'lat', minlat_source)
                CALL minmax_range_make(type_select, edgew_temp, edgee_temp, edgen_temp, edges_temp) ! Update global bounds for this mask type.
                print*, "Source grid index range: minlon=", minlon_source, " maxlon=", maxlon_source, " maxlat=", maxlat_source, " minlat=", minlat_source
                if (minlon_source >= maxlon_source .OR. maxlat_source >= minlat_source) then
                    print*, "Warning: Invalid index range for circle ", k, " - skipping."
                    cycle
                endif
                allocate(num_points_ij(maxlon_source-minlon_source+1, minlat_source-maxlat_source+1)); num_points_ij = 0

                !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
                !$OMP PRIVATE(i, j, point_i, is_inside_circle)
                do i = minlon_source, maxlon_source - 1, 1 ! Loop through source grid cells in bounding box.
                    do j = maxlat_source, minlat_source - 1, 1
                        if (IsInArea_grid(i, j) == 1) cycle ! Skip if already marked by a previous circle/shape.
                        point_i = [lon_i(i), lat_i(j)]
                        is_inside_circle = is_point_in_circle(point_i, point_c, radius_c)
                        if (is_inside_circle) then 
                            IsInArea_grid(i, j) = 1
                            num_points_ij(i-minlon_source+1, j-maxlat_source+1) = 1
                        end if
                    end do
                end do
                !$OMP END PARALLEL DO
                num_points(k) = sum(num_points_ij)
                deallocate(num_points_ij)
            end do
            if (present(numpatch)) numpatch = sum(num_points)
            deallocate(num_points, circle_points, circle_radius)
        end do

    END SUBROUTINE IsInArea_circle_Calculation

   ! Determines if source grid cells are within a closed polygon defined area.
    ! Reads closed polygon definitions (list of vertices) from temporary files.
    !For each polygon:
    !1. Checks for self-intersections in the polygon.
    !2. Determines its bounding box and handles dateline crossing for the bounding box.
    !3. For each latitude row spanning the polygon:
    !   - Calculates intersections of a horizontal ray with all polygon segments.
    !   - Sorts intersection longitudes.
    !   - Marks source grid cells between odd-even pairs of intersections as inside.
    !   - Handles dateline crossing for filling.
    SUBROUTINE IsInArea_close_Calculation(type_select, iter, ndm, IsInArea_grid, numpatch)
        ! Calculates containment for closed polygons using the ray casting method.
        ! Accounts for crossing the 180th meridian.
        IMPLICIT NONE
        character(16), intent(in) :: type_select       ! Type of mask.
        integer, intent(in) :: iter, ndm                ! Refinement degree, number of mask files.
        character(pathlen) :: lndname                   ! Path to mask definition file.
        integer :: i, j, k, ii, close_num, n, icl_point, num_intersect ! Loop counters, num vertices in polygon, dateline flag, intersection count.
        integer :: minlon_source, maxlon_source, maxlat_source, minlat_source ! Source grid indices for polygon bounding box.
        real(r8) :: point_i(2), point_c(2)              ! Current grid point, temp for dateline check.
        real(r8) :: lon1, lat1, lon2, lat2, lon_intersect ! Segment endpoint coordinates, intersection longitude.
        real(r8) :: edgee_temp, edgew_temp, edges_temp, edgen_temp ! Polygon bounding box.
        real(r8), allocatable :: close_points(:,:)      ! Coordinates of polygon vertices.
        real(r8), allocatable :: ray_segment_intersect_lon(:,:) ! Stores longitudes of ray intersections for each latitude row.
        integer,  allocatable :: ray_segment_intersect_num(:)   ! Number of intersections for each latitude row.
        integer, dimension(:,:), intent(inout) :: IsInArea_grid   ! Output: Grid indicating cells within polygons.
        integer, intent(out), optional :: numpatch               ! Optional Output: Total count of grid cells marked inside.
        character(5) :: numc, refinec                            ! Character representation of n and iter.
        logical :: intersect                                     ! Flag for segment intersection.

        if (present(numpatch)) numpatch = 0
        write(refinec, '(I1)') iter
        do n = 1, ndm, 1 ! Loop through mask files.
            write(numc, '(I2.2)') n
            ! Filename construction seems to assume iter=0 for 'close' type, or refinec should be used.
            !lndname = trim(file_dir)// 'tmpfile/'//trim(type_select)//'_close_0_'//trim(numc)//'.nc4'
            lndname = trim(file_dir)// 'tmpfile/'//trim(type_select)//'_close_'//trim(refinec)//'_'//trim(numc)//'.nc4'
            print*, trim(lndname)
            CALL close_Mesh_Read(lndname, close_num, close_points) ! Read polygon vertices.
            
            ! Check for self-intersections in the polygon.
            CALL check_self_intersection(close_num, close_points)
            print*, "No line segment self intersection found!"

            ! Determine polygon bounding box.
            edgew_temp = minval(close_points(:, 1))
            edgee_temp = maxval(close_points(:, 1))
            edgen_temp = maxval(close_points(:, 2))
            edges_temp = minval(close_points(:, 2))
            print*, "Polygon ", n, " BBox: W=", edgew_temp, " E=", edgee_temp, " N=", edgen_temp, " S=", edges_temp
            
            ! Check for dateline crossing of the polygon's bounding box.
            point_c = 0.0_r8
            do i = 1, close_num - 1, 1
                if (point_c(1) < abs(close_points(i+1, 1) - close_points(i, 1))) then
                    point_c(1) = abs(close_points(i+1, 1) - close_points(i, 1))
                end if
            end do
            point_c(2) = abs(edgee_temp-edgew_temp)
            icl_point = 0 ! Initialize dateline crossing flag.
            if (point_c(1) > point_c(2)) then ! If max segment longitude span > overall BBox longitude span.
                print*, "Polygon crosses 180 meridian! Modifying coordinates for ray casting."
                icl_point = 1
                edgew_temp = -180.0_r8; edgee_temp =  180.0_r8 ! Expand BBox.
                CALL CheckCrossing(close_num, close_points) ! Adjust polygon coordinates.
                print*, "Adjusted BBox: W=", edgew_temp, " E=", edgee_temp, " N=", edgen_temp, " S=", edges_temp
            else
                print*, "Polygon does not cross 180 meridian."
            end if
            CALL minmax_range_make(type_select, edgew_temp, edgee_temp, edgen_temp, edges_temp) ! Update global bounds.

            ! Determine source grid index range for this polygon.
            CALL Source_Find(edgen_temp, lat_vertex, 'lat', maxlat_source)
            CALL Source_Find(edges_temp, lat_vertex, 'lat', minlat_source)
            print*, "Lat index range for ray casting: maxlat_idx=", maxlat_source, " minlat_idx=", minlat_source 
            
            ! Get intersections of rays with polygon segments for each latitude row.
            allocate(ray_segment_intersect_lon(minlat_source - maxlat_source + 1, close_num*2)) ! Max intersections = num_segments * 2 (worst case).
            ray_segment_intersect_lon = -9999.0_r8 ! Initialize with a very small number.
            allocate(ray_segment_intersect_num(minlat_source - maxlat_source + 1)); ray_segment_intersect_num = 0

            do j = maxlat_source, minlat_source - 1, 1 ! For each latitude row in the polygon's lat range.
                point_i = [-200.0_r8, lat_i(j)] ! Ray starting point (far west of current latitude).
                do i = 1, close_num, 1 ! For each segment of the polygon.
                    lon1 = close_points(i, 1)
                    lat1 = close_points(i, 2)
                    lon2 = close_points(mod(i, close_num) + 1, 1) ! Next vertex (cyclic).
                    lat2 = close_points(mod(i, close_num) + 1, 2)
                    CALL ray_segment_intersect(point_i, lat1, lon1, lat2, lon2, lon_intersect)
		    ! TODO: HERE is also need to check, I am not fully understand but I think it should be like this
		    ! /= is not good
                    if (lon_intersect > point_i(1) + 1.0e-6_r8) then ! If intersection found to the right of ray start.
                           !if (lon_intersect /= point_i(1)) then 
                        ray_segment_intersect_num(j - maxlat_source + 1) = ray_segment_intersect_num(j - maxlat_source + 1) + 1
                        if (ray_segment_intersect_num(j - maxlat_source + 1) > size(ray_segment_intersect_lon,2)) then
                            print*, "Error: Exceeded allocated size for ray_segment_intersect_lon"
                            stop
                        endif
                        ray_segment_intersect_lon(j - maxlat_source + 1,   ray_segment_intersect_num(j - maxlat_source + 1)) = lon_intersect
                    end if
                end do
            end do

            ! Sort intersection longitudes for each latitude row.
            CALL bubble_sort(ray_segment_intersect_num, ray_segment_intersect_lon)

            ! Fill IsInArea_grid based on sorted intersections (odd-even rule).
            if (icl_point == 1) then ! Polygon crossed dateline, handle remapping of indices.
                do j = maxlat_source, minlat_source - 1, 1 ! For each latitude row.
                    num_intersect = int(ray_segment_intersect_num(j - maxlat_source + 1)/2) ! Number of pairs of intersections.
                    do k = 1, num_intersect, 1
                        edgew_temp = ray_segment_intersect_lon(j - maxlat_source + 1, 2*k-1) ! Start of inside segment.
                        edgee_temp = ray_segment_intersect_lon(j - maxlat_source + 1, 2*k)   ! End of inside segment.
                        CALL Source_Find(edgew_temp, lon_vertex, 'lon', minlon_source)
                        CALL Source_Find(edgee_temp, lon_vertex, 'lon', maxlon_source)
                        if (present(numpatch)) numpatch = numpatch + (maxlon_source - minlon_source) !same result, but more logical
                        do i = minlon_source, maxlon_source - 1, 1 ! For longitude indices within this segment.
                            ! Remap longitude index 'ii' due to dateline crossing adjustment of polygon.
                            if (i < nlons_source/2 + 1) then 
                                ii = int( i + nlons_source/2)
                            else
                                ii = int( i - nlons_source/2)
                            end if
                            IsInArea_grid(ii, j) = 1
                        end do
                    end do
                end do
            else ! Polygon does not cross dateline.
                do j = maxlat_source, minlat_source - 1, 1
                    num_intersect = int(ray_segment_intersect_num(j - maxlat_source + 1)/2)
                    do k = 1, num_intersect, 1
                        edgew_temp = ray_segment_intersect_lon(j - maxlat_source + 1, 2*k-1)
                        edgee_temp = ray_segment_intersect_lon(j - maxlat_source + 1, 2*k)
                        CALL Source_Find(edgew_temp, lon_vertex, 'lon', minlon_source)
                        CALL Source_Find(edgee_temp, lon_vertex, 'lon', maxlon_source)
                        if (present(numpatch)) numpatch = numpatch + (maxlon_source - minlon_source)
                        IsInArea_grid(minlon_source:maxlon_source - 1, j) = 1
                    end do
                end do
            end if
            deallocate(ray_segment_intersect_lon, ray_segment_intersect_num, close_points)
        end do

    END SUBROUTINE IsInArea_close_Calculation

   ! Calculates the intersection point of a horizontal ray with a line segment.
    ! Given a ray starting at `point_i` and extending horizontally to the right (positive longitude),
    !and a line segment defined by (lon1, lat1) and (lon2, lat2), this finds the longitude
    !of intersection if the ray intersects the segment.
    ! point_i Input: Start point of the ray [lon_p, lat_p].
    ! lat1, lon1 Input: Coordinates of the first endpoint of the segment.
    ! lat2, lon2 Input: Coordinates of the second endpoint of the segment.
    ! lon_intersect Output: Longitude of intersection. If no intersection, it's set to `point_i(1)`.
    SUBROUTINE ray_segment_intersect(point_i, lat1, lon1, lat2, lon2, lon_intersect)
        implicit none
        real(r8), intent(in)  :: point_i(2)     ! Ray start point [lon_p, lat_p].
        real(r8), intent(in)  :: lat1, lon1, lat2, lon2 ! Segment endpoints.
        real(r8), intent(out) :: lon_intersect      ! Output: Longitude of intersection or ray start lon if no intersection.
        !local
	real(r8) :: lat_p, lon_p,  m              ! Ray coordinates, slope of the segment.

        lon_p = point_i(1); lat_p = point_i(2) ! Ray's starting longitude and latitude.

        ! If the segment is horizontal (lat1 == lat2), it cannot intersect a horizontal ray unless collinear.
        ! For simplicity in ray casting, horizontal segments are often ignored or handled by specific rules.
        ! Here, it's treated as no intersection for standard ray casting (returns ray's own longitude).
	!do not use ==, Allow for slight numerical error
        if (abs(lat1 - lat2) < 1e-9_r8) then ! Use tolerance for float comparison
            lon_intersect = lon_p
            return
        end if

        ! If the ray's latitude is outside the segment's latitude range, no intersection.
        if ((lat1 > lat_p .and. lat2 > lat_p) .or. (lat1 < lat_p .and. lat2 < lat_p)) then
            lon_intersect = lon_p
            return
        end if
        
	! TODO: CHECK here, potential bugs here in original code
        ! Calculate slope of the segment. If vertical, intersection lon is lon1.
        if (abs(lon2 - lon1) < 1e-9_r8) then
            lon_intersect = lon1
        else
            m = (lat2 - lat1) / (lon2 - lon1) ! Slope.
            if (abs(m) < 1e-9_r8) then ! Segment is horizontal (already checked, but as a safeguard for near-horizontal)
                lon_intersect = lon_p
                return
            endif
            lon_intersect = lon1 + (lat_p - lat1) / m ! X-coordinate of intersection.
        endif
        
        ! Check if intersection point lies within the segment's longitude bounds.
        ! This is crucial. If not, it's not a valid intersection with the segment.
        if (.not. ((lon_intersect >= min(lon1,lon2) - 1e-6_r8) .and. (lon_intersect <= max(lon1,lon2) + 1e-6_r8))) then
             lon_intersect = lon_p ! Not within segment's longitude span.
        endif
        ! Additional check for horizontal ray: intersection must be to the right of ray origin
        if (lon_intersect < lon_p - 1e-6_r8) then ! Allow for slight numerical error
            lon_intersect = lon_p
        endif

    END SUBROUTINE ray_segment_intersect

   ! Sorts intersection longitudes for each latitude row using bubble sort.
    ! ray_segment_intersect_num Input: Array containing the number of intersections for each latitude row.
    ! ray_segment_intersect_lon Input/Output: 2D array of intersection longitudes [lat_row_idx, intersection_idx]. Sorted in place.
    SUBROUTINE bubble_sort(ray_segment_intersect_num, ray_segment_intersect_lon)

        implicit none
        integer,  intent(in) :: ray_segment_intersect_num(:)    ! Number of intersections per latitude row.
        integer :: i, j, k, num                                 ! Loop counters.
        real(r8) :: temp                                        ! Temporary variable for swapping.
        real(r8), intent(inout) :: ray_segment_intersect_lon(:, :) ! Intersection longitudes (sorted in place).

        do k = 1, size(ray_segment_intersect_num), 1 ! For each latitude row.
            num = ray_segment_intersect_num(k) ! Number of intersections in this row.
            if (num <= 1) cycle ! No need to sort if 0 or 1 intersection.
            do i = 1, num - 1 ! Standard bubble sort.
                do j = i + 1, num
                    if (ray_segment_intersect_lon(k, i) > ray_segment_intersect_lon(k, j)) then
                        temp = ray_segment_intersect_lon(k, i)
                        ray_segment_intersect_lon(k, i) = ray_segment_intersect_lon(k, j)
                        ray_segment_intersect_lon(k, j) = temp
                    end if
                end do
            end do
        end do

    END SUBROUTINE bubble_sort

   ! Checks if a polygon (defined by `close_points`) has any self-intersections.
    ! It iterates through all pairs of non-adjacent segments of the polygon
    !and calls `segments_intersect` to check if they cross. Stops if an intersection is found.
    !@note The conversion to Cartesian coordinates commented out suggests this might be simplified
    !or that `segments_intersect` works with geographic coordinates if they are locally planar.
    SUBROUTINE check_self_intersection(close_num, close_points)
        ! Uses planar geometry line segment intersection algorithm by checking cross products.
        implicit none
        integer, intent(in) :: close_num                ! Number of vertices in the polygon.
        real(r8), allocatable, intent(in) :: close_points(:,:) ! Coordinates of polygon vertices.
        integer :: n, i, j                              ! Loop counters.
        real(r8), allocatable :: lon(:), lat(:)          ! Temporary arrays for vertex coordinates.
        logical :: intersect                             ! Flag, true if segments intersect.

        allocate(lon(close_num+1)); lon = 0.0_r8
        allocate(lat(close_num+1)); lat = 0.0_r8
        ! Read longitudes and latitudes. (Original comment)
        ! Convert to Cartesian coordinates (Original comment - this block is commented out in provided code).
        ! do i = 1, close_num, 1
        !     ! Convert to Cartesian.
        !     ! lon(i) = erad * cos(close_points(i, 2)*pio180) * cos(close_points(i, 1)*pio180) ! pio180 = pi/180.0
        !     ! lat(i) = erad * cos(close_points(i, 2)*pio180) * sin(close_points(i, 1)*pio180)
        ! end do
        lon(1:close_num) = close_points(1:close_num, 1)
        lat(1:close_num) = close_points(1:close_num, 2)
        lon(1+close_num) = lon(1)! Close the polygon for segment iteration.
        lat(1+close_num) = lat(1)

        ! Check if line segments self-intersect.
        do i = 1, close_num - 2, 1 ! For each segment (i, i+1).
            do j = i + 2, close_num, 1 ! For each non-adjacent segment (j, j+1), excluding the segment itself and its immediate neighbor.
                                       ! The last segment is (close_num, 1) which is (close_num, close_num+1) due to cyclic copy.
                !TODO: check here, if they share same vertex
		if (j == close_num .and. i == 1) cycle ! Avoid checking segment (1,2) against (num,1) if they share vertex 1.
                intersect = segments_intersect(lat(i), lon(i), lat(i+1), lon(i+1), &
                                        lat(j), lon(j), lat(j+1), lon(j+1))
                if (intersect) then
                    print*, "Segment i: (", lon(i), ",", lat(i), ") to (", lon(i+1), ",", lat(i+1), ")"
                    print*, "Segment j: (", lon(j), ",", lat(j), ") to (", lon(j+1), ",", lat(j+1), ")"
                    stop "ERROR! Segments i and Segments j intersect."
                end if
            end do
        end do
        deallocate(lon, lat)
    END SUBROUTINE check_self_intersection

   ! Checks if two line segments intersect using the cross-product method.
    ! Calculates four cross products to determine if the segments straddle each other.
    !Assumes segments are (lon1,lat1)-(lon2,lat2) and (lon3,lat3)-(lon4,lat4).
    !@return .TRUE. if segments intersect, .FALSE. otherwise.
    LOGICAL function segments_intersect(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4)
        implicit none
        real(r8), intent(in) :: lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4 ! Coordinates of segment endpoints.
        real(r8) :: cp1, cp2, cp3, cp4                                      ! Cross products.

        ! Calculate cross products.
        ! cp1: (P2-P1) x (P3-P1)
        cp1 = cross_product2(lon2 - lon1, lat2 - lat1, lon3 - lon1, lat3 - lat1)
        ! cp2: (P2-P1) x (P4-P1)
        cp2 = cross_product2(lon2 - lon1, lat2 - lat1, lon4 - lon1, lat4 - lat1)
        ! cp3: (P4-P3) x (P1-P3)
        cp3 = cross_product2(lon4 - lon3, lat4 - lat3, lon1 - lon3, lat1 - lat3)
        ! cp4: (P4-P3) x (P2-P3)
        cp4 = cross_product2(lon4 - lon3, lat4 - lat3, lon2 - lon3, lat2 - lat3)

        ! Check if segments intersect. Intersection occurs if orientations change (signs of cross products differ).
        if ((cp1 * cp2 < 0.0_r8) .and. (cp3 * cp4 < 0.0_r8)) then
            segments_intersect = .true.
        else
            segments_intersect = .false.
        end if
    end function segments_intersect

   ! Calculates the 2D cross product of two vectors (x1,y1) and (x2,y2).
    !@return The scalar value (z-component) of the 3D cross product (x1*y2 - x2*y1).
    Real(r8) function cross_product2(x1, y1, x2, y2)
        implicit none
        real(r8), intent(in)  :: x1, y1, x2, y2   ! Components of the two vectors.
        cross_product2 = x1 * y2 - x2 * y1
        
    end function cross_product2

   ! Updates the min/max lat/lon indices for a given mask type based on a new bounding box.
    ! This subroutine takes a bounding box (edgew_temp, etc.) and determines the
    !corresponding source grid indices (temp1-temp4). It then updates the global
    !min/max indices for the specified `type_select` (domain, refine, or patch)
    !if the new box extends the current overall bounds.
    ! type_select Input: Type of mask ('mask_domain', 'mask_refine', 'mask_patch').
    ! edgew_temp,edgee_temp,edgen_temp,edges_temp Input: West,East,North,South boundaries of the new box.
    ! temp1,temp2,temp3,temp4 Optional Output: Source grid indices [minlon,maxlon,maxlat,minlat] for the input box.
    SUBROUTINE minmax_range_make(type_select, edgew_temp, edgee_temp, edgen_temp, edges_temp, temp1, temp2, temp3, temp4)
        ! Determines the minimum rectangle that can contain all data.
        implicit none
        character(16), intent(in) :: type_select  ! Type of mask being processed.
        real(r8), intent(in) :: edgew_temp, edgee_temp, edgen_temp, edges_temp ! Input bounding box.
        integer, intent(out), optional :: temp1, temp2, temp3, temp4  ! Optional Output: Indices for the input box.
        integer :: temp1_local, temp2_local, temp3_local, temp4_local  ! Local variables for calculated indices.

        ! Call Source_Find to get source grid indices corresponding to the bounding box coordinates.
        CALL Source_Find(edgew_temp, lon_vertex, 'lon', temp1_local)  ! minlon_source index
        CALL Source_Find(edgee_temp, lon_vertex, 'lon', temp2_local)  ! maxlon_source index
        CALL Source_Find(edgen_temp, lat_vertex, 'lat', temp3_local)  ! maxlat_source index (smallest index for northernmost)
        CALL Source_Find(edges_temp, lat_vertex, 'lat', temp4_local)  ! minlat_source index (largest index for southernmost)

        ! Adjust indices (Source_Find might return index of vertex line *at or after* the coordinate).
        ! For maxlon_DmArea/maxlon_RfArea, we want the cell index *before* the eastern edge, so -1 or -2.
        ! For minlat_DmArea/minlat_RfArea, we want cell index *before* the southern edge.
        temp2_local = temp2_local - 2 ! Adjust max longitude index.
        temp4_local = temp4_local - 2 ! Adjust min latitude index (max index value).
        if (temp2_local == nlons_source-1) temp2_local = temp2_local + 1 ! Boundary adjustment.
        if (temp4_local == nlats_source-1) temp4_local = temp4_local + 1 ! Boundary adjustment.

        ! Update global min/max indices for the specified mask type.
        if (type_select == 'mask_domain') then
            if (temp1_local < minlon_DmArea) minlon_DmArea = temp1_local
            if (temp2_local > maxlon_DmArea) maxlon_DmArea = temp2_local
            if (temp3_local < maxlat_DmArea) maxlat_DmArea = temp3_local
            if (temp4_local > minlat_DmArea) minlat_DmArea = temp4_local
        else if (type_select == 'mask_refine') then
            if (temp1_local < minlon_RfArea) minlon_RfArea = temp1_local
            if (temp2_local > maxlon_RfArea) maxlon_RfArea = temp2_local
            if (temp3_local < maxlat_RfArea) maxlat_RfArea = temp3_local
            if (temp4_local > minlat_RfArea) minlat_RfArea = temp4_local
        else if (type_select == 'mask_patch') then
            if (temp1_local < minlon_PaArea) minlon_PaArea = temp1_local
            if (temp2_local > maxlon_PaArea) maxlon_PaArea = temp2_local
            if (temp3_local < maxlat_PaArea) maxlat_PaArea = temp3_local
            if (temp4_local > minlat_PaArea) minlat_PaArea = temp4_local
        end if

        ! If output arguments are present, assign the calculated local indices.
        if (present(temp1)) temp1 = temp1_local
        if (present(temp2)) temp2 = temp2_local
        if (present(temp3)) temp3 = temp3_local
        if (present(temp4)) temp4 = temp4_local
    END SUBROUTINE minmax_range_make

   ! Finds the index in a sorted coordinate array (`seq_lonlat`) that corresponds to a given coordinate `temp`.
    ! It performs a search within a slightly expanded range around an estimated position.
    !For longitude, it finds `i` such that `temp <= seq_lonlat(i)`.
    !For latitude, it finds `i` such that `temp >= seq_lonlat(i)` (assuming latitude array is N to S).
    ! temp Input: The coordinate value (lon or lat) to find.
    ! seq_lonlat Input: Sorted array of global grid vertex longitudes or latitudes.
    ! str1 Input: String "lon" or "lat" to indicate coordinate type.
    ! source Output: The found 1-based index in `seq_lonlat`.
    SUBROUTINE Source_Find(temp, seq_lonlat, str1, source)
        implicit none
        real(r8), intent(in) :: temp                 ! Coordinate value to find.
        real(r8), dimension(:), intent(in) :: seq_lonlat ! Sorted array of vertex coordinates.
        character(LEN = 3), intent(in) :: str1       ! "lon" or "lat".
        integer :: i, gridnum_perdegree, minsource, maxsource ! Loop counter, grid cells per degree, search range.
        integer, intent(out) :: source                ! Output: Found index.

        gridnum_perdegree = 120 ! Default grid cells per degree (e.g., for 0.25 deg resolution, 4 cells/deg, but this seems different).
                               ! This variable seems related to an expected resolution to define search window.
        if (lcs == 'igbp') gridnum_perdegree = 240 ! Higher for IGBP?
        !TODO: check here whether the revision is correct or not
        if (trim(str1) == 'lon') then ! For longitude sequence. Longitude from -180 to 180.
            ! Estimate search range.
            minsource = int( ( temp - (-180.0_r8) ) * gridnum_perdegree ) ! This formula is problematic.
                                                                      ! Should be based on array index and value relationship.
                                                                      ! Assuming seq_lonlat(1) = -180.
            maxsource = int( ( temp - (-180.0_r8) ) * gridnum_perdegree ) + 2*gridnum_perdegree ! Broaden range.
            minsource = max(1,              minsource-10) ! Ensure within array bounds, with some buffer.
            maxsource = min(1+nlons_source, maxsource+10)
            do i = minsource, maxsource, 1
                if (temp <= seq_lonlat(i)) then ! Find first vertex line at or to the east of temp.
                    source = i ! For minlon_source, this 'i' is the index of the eastern boundary of cell i-1.
                               ! So cell i-1 contains temp if temp > seq_lonlat(i-1).
                               ! If searching for cell index, usually source = i-1.
                               ! If searching for vertex index, source = i.
                    return
                end if
            end do
            ! Error if not found (e.g., temp is outside range of seq_lonlat after adjustments).
            print*, "temp <= seq_lonlat(i) is not exist!!!!!!!!!!!!!!!"
            print*, "temp = ", temp, "str1 = ", str1
            print*, "minsource = ", minsource, "maxsource = ", maxsource
            if (minsource <= size(seq_lonlat) .and. minsource > 0) print*, "seq_lonlat(minsource) = ", seq_lonlat(minsource)
            if (maxsource <= size(seq_lonlat) .and. maxsource > 0) print*, "seq_lonlat(maxsource) = ", seq_lonlat(maxsource)
            source = i ! Returns last index if not found.
        else ! For latitude sequence. Latitude from 90 to -90.
            !TODO:check here!
	    ! which is correct???!!
            minsource = ( 90 - ceiling(temp)) * gridnum_perdegree
            maxsource = ( 90 -   floor(temp)) * gridnum_perdegree
            !minsource = int( ( 90.0_r8 - temp ) * gridnum_perdegree )
            !maxsource = int( ( 90.0_r8 - temp ) * gridnum_perdegree ) + 2*gridnum_perdegree
            minsource = max(1,              minsource-10)
            maxsource = min(1+nlats_source, maxsource+10)
            do i = minsource, maxsource, 1
                if (temp >= seq_lonlat(i)) then ! Find first vertex line at or to the south of temp.
                                               ! seq_lonlat for latitude is likely 90, 90-dy, ... -90.
                                               ! So seq_lonlat(i) > seq_lonlat(i+1).
                                               ! If temp >= seq_lonlat(i), then temp is in cell i-1 (if index is N->S).
                    source = i
                    return
                end if
            end do
            print*, "temp >= seq_lonlat(i) is not exist!!!!!!!!!!!!!!!"
            print*, "temp = ", temp, "str1 = ", str1
            print*, "minsource = ", minsource, "maxsource = ", maxsource
            if (minsource <= size(seq_lonlat) .and. minsource > 0) print*, "seq_lonlat(minsource) = ", seq_lonlat(minsource)
            if (maxsource <= size(seq_lonlat) .and. maxsource > 0) print*, "seq_lonlat(maxsource) = ", seq_lonlat(maxsource)
            source = i
        end if
    END SUBROUTINE Source_Find

   ! Adjusts longitudes by +/- 180 degrees if they cross the dateline.
    ! This is a specific coordinate transformation used when a polygon/set of points
    !is known to cross the 180-degree meridian. It shifts points with negative longitudes
    !by +180 and positive by -180. This is NOT a general normalization to [-180, 180]
    !but a specific transformation likely for a particular geometric algorithm
    !that expects points to be on one side of a shifted meridian.
    SUBROUTINE CheckCrossing(num_edges, points)
        implicit none
        integer, intent(in) :: num_edges                      ! Number of points.
        real(r8), dimension(num_edges, 2), intent(inout) :: points ! Coordinates [(lon,lat)].
        integer :: j                                          ! Loop counter.
        do j = 1, num_edges, 1
            if(points(j, 1) < 0.0_r8)then
                points(j, 1) = points(j, 1) + 180.0_r8
            else
                points(j, 1) = points(j, 1) - 180.0_r8
            end if
        end do
    END SUBROUTINE CheckCrossing

   ! Checks if a point is inside a circle on a sphere.
    ! Calculates the great-circle distance between `point_i` and circle center `point_c`
    !using the Haversine formula.
    ! point_i Input: Coordinates [lon, lat] of the point to check.
    ! point_c Input: Coordinates [lon, lat] of the circle's center.
    ! center_radius Input: Radius of the circle in kilometers.
    !@return .TRUE. if `point_i` is inside or on the boundary of the circle, .FALSE. otherwise.
    LOGICAL FUNCTION is_point_in_circle(point_i, point_c, center_radius)
        ! Determines if a point is inside a circle based on distance.
        real(r8), intent(in) :: point_i(2), point_c(2)  ! Point to check, circle center [lon, lat] in degrees.
        real(r8), intent(in) :: center_radius         ! Circle radius in kilometers.
        real(r8) :: distance                          ! Calculated distance.

        is_point_in_circle = .true.
        distance = haversine(point_i, point_c) ! Calculate distance in km.
        if (distance > center_radius) is_point_in_circle = .false.

    END FUNCTION is_point_in_circle

   ! Calculates the great-circle distance between two points on a sphere using the Haversine formula.
    ! point_i Input: Coordinates [lon, lat] of the first point in degrees.
    ! point_c Input: Coordinates [lon, lat] of the second point in degrees.
    !@return Distance in kilometers.
    REAL(r8) Function haversine(point_i, point_c)

        implicit none
        ! Haversine formula, measures in radians!! Medium distance: suitable for hundreds to thousands of kilometers.
        ! Reference: https://zhuanlan.zhihu.com/p/658990378 (Chinese)
        real(r8), intent(in) :: point_i(2), point_c(2) ! Coordinates [lon,lat] in degrees.
        real(r8) :: px1, py1, px2, py2, v             ! Coordinates in radians, intermediate variable.
        px1 = point_i(1) * pio180 ! lon1 in radians.
        py1 = point_i(2) * pio180 ! lat1 in radians.
        px2 = point_c(1) * pio180 ! lon2 in radians.
        py2 = point_c(2) * pio180 ! lat2 in radians.
        
        ! v = sqrt( sin(py1/2-py2/2)**2 + cos(py2)*cos(py1)*sin(px1/2-px2/2)**2 )
        ! haversine = erad * 2 * asin(v)
        ! new
        v = sin((py1-py2)/2.0_r8)**2 + cos(py2)*cos(py1)*sin((px1-px2)/2.0_r8)**2
        haversine = erad / 1000.0_r8 * 2.0_r8 * atan2(sqrt(v), sqrt(1.0_r8-v)) ! Convert m to km. 2*atan2 is 2*asin(sqrt(v)) for small v.
        return

    END Function haversine

   ! Checks if a point is inside a convex polygon using the cross-product method.
    ! Iterates through the edges of the polygon. For each edge (p1, p2), it calculates
    !the cross product of vectors (p2-p1) and (ndm_point-p1). If the sign of the cross product
    !is consistent for all edges (all positive or all negative, depending on vertex order),
    !the point is inside. If signs change, the point is outside.
    ! polygon Input: Array of vertex coordinates [(lon,lat)] for the polygon. Assumed max 7 vertices.
    ! ndm_point Input: Coordinates [lon, lat] of the point to check.
    ! num_edges Input: Actual number of vertices/edges in the polygon.
    !@return .TRUE. if the point is inside or on the boundary, .FALSE. otherwise.
    LOGICAL FUNCTION is_point_in_convex_polygon(polygon, ndm_point, num_edges)
        ! Determines if a point is inside a convex polygon.
        real(r8), intent(in) :: polygon(7, 2)    ! Polygon vertices (max 7 allowed by fixed dimension).
        real(r8), intent(in) :: ndm_point(2)     ! Point to check [lon, lat].
        integer, intent(in) :: num_edges         ! Actual number of edges in the polygon.
        real(r8) :: p1(2), p2(2)                 ! Current edge endpoints.
        integer :: i, next_index                ! Loop counter, index of next vertex.
        real(r8) :: prev_cross, cross            ! Previous and current cross product values.
        real(r8), dimension(:,:), allocatable :: polygon_select ! Local copy of relevant polygon vertices.

        allocate(polygon_select(num_edges,2))
        polygon_select = polygon(1:num_edges, :) ! Use only the actual number of edges.
        prev_cross = 0.0_r8
        is_point_in_convex_polygon = .true. ! Assume inside until proven otherwise.

        do i = 1, num_edges
            next_index = mod(i, num_edges) + 1  ! Next vertex, handles wrap-around for the last edge.
            p1 = polygon_select(i, :)           ! Current vertex.
            p2 = polygon_select(next_index, :)  ! Next vertex.
            cross = cross_product(p1, p2, ndm_point) ! Cross product (P2-P1) x (Point-P1).
            if (abs(cross) > 1e-9_r8) then ! If not collinear (within tolerance).
                if (abs(prev_cross) < 1e-9_r8) then ! If this is the first non-zero cross product.
                    prev_cross = cross
                else if ((prev_cross > 0.0_r8) .and. (cross < 0.0_r8) .or. &
                         (prev_cross < 0.0_r8) .and. (cross > 0.0_r8)) then ! If sign changes.
                    is_point_in_convex_polygon = .false. ! Point is outside.
                    deallocate(polygon_select)
                    return
                end if
            end if
        end do
        deallocate(polygon_select)
    END FUNCTION is_point_in_convex_polygon

   ! Calculates the 2D cross product of two vectors P1P2 and P1P3.
    ! (P2 - P1) x (P3 - P1) = (p2x-p1x)*(p3y-p1y) - (p2y-p1y)*(p3x-p1x).
    !The sign indicates whether P3 is to the left or right of vector P1P2.
    ! p1 Input: Coordinates of point P1.
    ! p2 Input: Coordinates of point P2.
    ! p3 Input: Coordinates of point P3.
    !@return Scalar result of the 2D cross product.
    REAL FUNCTION cross_product(p1, p2, p3)
        ! Calculates cross product: (p2 - p1) × (p3 - p1).
        real(r8), intent(in) :: p1(2), p2(2), p3(2)
        cross_product = (p2(1) - p1(1)) * (p3(2) - p1(2)) - (p2(2) - p1(2)) * (p3(1) - p1(1))
    END FUNCTION cross_product
    
  ! Reads domain area and land/sea mask data from a NetCDF file.
   ! This is typically used in `mask_restart` mode to load previously processed masks.
   ! inputfile Input: Path to the NetCDF file.
   ! IsInDmArea_select Output: Mask read from file indicating cells in domain.
   ! seaorland_select Output: Land/sea mask read from file.
   SUBROUTINE IsInArea_grid_Read(inputfile, IsInDmArea_select, seaorland_select)
      USE NETCDF  
      IMPLICIT NONE
      character(len = 256), intent(in) :: inputfile     ! Input: File name.
      integer :: ncID, dimID_lon, dimID_lat, varid(12)   ! NetCDF IDs.
      integer :: nlons_select_read, nlats_select_read    ! Dimensions read from file.
      integer, allocatable, intent(out) :: IsInDmArea_select(:, :), seaorland_select(:, :) ! Output arrays.
      
      print*, trim(inputfile)
      CALL CHECK(NF90_OPEN(trim(inputfile), nf90_nowrite, ncid))
      CALL CHECK(NF90_INQ_DIMID(ncid, "nlons_select", dimID_lon))
      CALL CHECK(NF90_INQ_DIMID(ncid, "nlats_select", dimID_lat))
      CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, DimID_lon, len = nlons_select_read))
      CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, DimID_lat, len = nlats_select_read))

      allocate(IsInDmArea_select(nlons_select_read,  nlats_select_read))
      allocate(seaorland_select(nlons_select_read, nlats_select_read))
      !TODO:check here
      CALL CHECK(NF90_INQ_VARID(ncid, 'IsInArea_select', varid(1)))
      !CALL CHECK(NF90_INQ_VARID(ncid, 'IsInDmArea_select', varid(1)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'seaorland_select',  varid(2)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'minlon_DmArea', varid(3)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'maxlon_DmArea', varid(4)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'maxlat_DmArea', varid(5)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'minlat_DmArea', varid(6)))

      CALL CHECK(NF90_GET_VAR(ncID, varid(1), IsInDmArea_select))
      CALL CHECK(NF90_GET_VAR(ncID, varid(2), seaorland_select))
      CALL CHECK(NF90_GET_VAR(ncID, varid(3), minlon_DmArea))
      CALL CHECK(NF90_GET_VAR(ncID, varid(4), maxlon_DmArea))
      CALL CHECK(NF90_GET_VAR(ncID, varid(5), maxlat_DmArea))
      CALL CHECK(NF90_GET_VAR(ncID, varid(6), minlat_DmArea))
      CALL CHECK(NF90_CLOSE(ncID))
      ! Update global select dimensions based on what was read.
      nlons_Dm_select = nlons_select_read
      nlats_Dm_select = nlats_select_read
   END SUBROUTINE IsInArea_grid_Read

   ! Saves domain area and land/sea mask data to a NetCDF file.
   ! Saves a subset of `IsInArea_grid` and `seaorland` corresponding to the
   ! calculated bounding box (`minlon_Area`, etc.) of the active area.
   ! Also saves the longitude/latitude coordinates and bounding box indices.
   ! outputfile Input: Path for the output NetCDF file.
   ! IsInArea_grid Input: Full domain mask.
   ! minlon_Area,maxlon_Area,maxlat_Area,minlat_Area Input: Bounding box indices of the area to save.
   SUBROUTINE IsInArea_grid_Save(outputfile, IsInArea_grid, minlon_Area, maxlon_Area, maxlat_Area, minlat_Area)
      USE NETCDF 
      IMPLICIT NONE
      character(len = 256), intent(in) :: outputfile      ! Input: File name for saving.
      integer, intent(in) :: minlon_Area, maxlon_Area, maxlat_Area, minlat_Area ! Input: Bounding box of the area.
      integer, allocatable, intent(in) :: IsInArea_grid(:, :) ! Input: Mask grid to be saved (subset of).
      integer :: ncID, dimID_lon, dimID_lat, varid(8)        ! NetCDF IDs.
      integer :: nlons_select_save, nlats_select_save          ! Dimensions of the subset to save.
      integer, allocatable :: lons_select(:), lats_select(:)  ! Indices for subsetting.
      integer, allocatable :: IsInArea_select(:, :), seaorland_select(:, :) ! Subsetted arrays for saving.
      real(r8), allocatable :: longitude(:), latitude(:)    ! Coordinate values for the subset.
      integer :: i  ! Loop variable for array construction
      
      ! Define the subset of indices to save based on the active area's bounding box.
      if (maxlon_Area < minlon_Area .or. minlat_Area < maxlat_Area) then
          print*, "Warning: Invalid area specified for IsInArea_grid_Save. No data to save."
          print*, "minlon,maxlon,maxlat,minlat: ", minlon_Area, maxlon_Area, maxlat_Area, minlat_Area
          return
      endif
      allocate(lons_select(maxlon_Area - minlon_Area + 1))
      allocate(lats_select(minlat_Area - maxlat_Area + 1))
      do i = 1, maxlon_Area - minlon_Area + 1
          lons_select(i) = minlon_Area + i - 1
      end do
      do i = 1, minlat_Area - maxlat_Area + 1
          lats_select(i) = maxlat_Area + i - 1
      end do
      nlons_select_save = size(lons_select)
      nlats_select_save = size(lats_select)
      
      ! Allocate and populate arrays for the subset.
      allocate(longitude(nlons_select_save))
      allocate(latitude(nlats_select_save))
      allocate(IsInArea_select(nlons_select_save, nlats_select_save))
      allocate(seaorland_select(nlons_select_save, nlats_select_save))
      longitude = lon_i(lons_select)
      latitude  = lat_i(lats_select)  
      IsInArea_select  = IsInArea_grid(lons_select, lats_select) 
      seaorland_select = seaorland(lons_select, lats_select) 

      print*, trim(outputfile)
      CALL CHECK(NF90_CREATE(trim(outputfile), ior(nf90_clobber, nf90_netcdf4), ncID))
      CALL CHECK(NF90_DEF_DIM(ncID, "nlons_select", nlons_select_save, dimID_lon))
      CALL CHECK(NF90_DEF_DIM(ncID, "nlats_select", nlats_select_save, dimID_lat))
      CALL CHECK(NF90_DEF_VAR(ncID, "IsInArea_select",  NF90_INT,   (/ dimID_lon, dimID_lat /), varid(1)))
      CALL CHECK(NF90_DEF_VAR(ncID, "seaorland_select", NF90_INT,   (/ dimID_lon, dimID_lat /), varid(2)))
      CALL CHECK(NF90_DEF_VAR(ncID, "minlon_DmArea",    NF90_INT,   varid(3))) ! Save scalar metadata.
      CALL CHECK(NF90_DEF_VAR(ncID, "maxlon_DmArea",    NF90_INT,   varid(4)))
      CALL CHECK(NF90_DEF_VAR(ncID, "maxlat_DmArea",    NF90_INT,   varid(5)))
      CALL CHECK(NF90_DEF_VAR(ncID, "minlat_DmArea",    NF90_INT,   varid(6)))
      CALL CHECK(NF90_DEF_VAR(ncID, "longitude",        NF90_DOUBLE, (/ dimID_lon /), varid(7))) ! Changed to NF90_DOUBLE
      CALL CHECK(NF90_DEF_VAR(ncID, "latitude",         NF90_DOUBLE, (/ dimID_lat /), varid(8))) ! Changed to NF90_DOUBLE
      CALL CHECK(NF90_ENDDEF(ncID))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(1), IsInArea_select))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(2), seaorland_select))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(3), minlon_Area))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(4), maxlon_Area))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(5), maxlat_Area))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(6), minlat_Area))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(7), longitude))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(8), latitude))
      CALL CHECK(NF90_CLOSE(ncID))
      deallocate(lons_select, lats_select, IsInArea_select, seaorland_select, longitude, latitude)
   
   END SUBROUTINE IsInArea_grid_Save

END Module MOD_Area_judge
