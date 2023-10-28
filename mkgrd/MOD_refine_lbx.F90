! --------------------------------------
! Preliminary refinement of polygonal substructure mesh (single layer)
! *** Initial data
! ***_new Updated data
! ***_f Final data
! --------------------------------------

module MOD_refine_lbx
    USE consts_coms
    USE refine_vars
    USE NETCDF
    USE MOD_Get_Distance
    implicit none
Contains

    SUBROUTINE refine_lbx()

        implicit none
        
        integer :: i, j, k, L, m, n, x, y, z
        integer :: w1, w2, w3, w4, w5, w6, m1, m2, m3, m4      ! Serial number of triangles and polygons
        integer :: row, col, error, ustr_points, edges
        integer :: num_ref                                     ! Refine the number of triangles each time 
        integer :: num_all                                     ! The unstructured grid contains the total number of structured grids
        integer :: refed(1000)                                 ! Record the number of triangles that have been refined at each refinement
        integer :: iter                                        ! Number of mesh refinement
        integer :: sa_iter                                     ! The number of iterations of mesh quality adjustment
        integer :: num_sjx, num_dbx                            ! The number of triangles and polygons after refinement
        integer :: sjx_points, lbx_points                      ! Number of triangles and polygons (read from the original file)
        integer :: nmp(1000), nwp(1000)                        ! Record the number of m and w points after each refinement
        integer :: icl(3)                                      ! Record the number of m and w points after each refinement
        !integer :: maxlc                                       ! Land type maximum number, 17 or 24

        ! nc文件读写相关
        integer :: ncid, varid(10), ncvarid(13), upDimID, iunit
        integer :: spDImID, lpDimID, twDimID, thDimID, seDimID, sxDimID, dimID_sjx, dimID_lbx

        integer :: n_wbx, n_lbx, n_qbx                         ! Number of pentagons, hexagons and heptagons

        real(r8), allocatable :: wp(:, :), mp(:, :)            ! Initial data of center point of triangular and polygon mesh
        real(r8), allocatable :: mp_new(:, :), wp_new(:, :)    ! Update data at the center point of triangular and polygon mesh
        real(r8), allocatable :: mp_f(:, :), wp_f(:, :)        ! Triangle, polygon mesh center point final data 
        real(r8), allocatable :: mp_f_tmp(:, :), wp_f_tmp(:, :) 
        real(r8), allocatable :: ref_lbx(:, :)                 ! Refinement of triangles that form polygons
        real(r8), allocatable :: dismm(:,:),disww(:,:)

        real(r8), allocatable :: f_mainarea(:, :), slope(:, :), lai(:, :)
        real(r8), allocatable :: k_s(:, :, :), k_sl(:, :, :),tkdry(:, :, :)
        real(r8), allocatable :: tksatf(:, :, :), tksatu(:, :, :)

        ! Mesh quality adjustment related
	     ! Extr_*** refers to the minimum Angle and maximum Angle of the grid
	     ! Eavg_*** refers to the average of the minimum Angle and maximum Angle of the grid
	     ! Savg_*** refers to the standard deviation between the Angle of the mesh and the Angle of the regular polygon
	     ! less30 refers to the number of angles less than 30 degrees in a triangular grid        
	     real(r8), allocatable :: less30(:), length(:, :, :), angle(:, :, :)
        real(r8), allocatable :: Savg_sjx(:), Extr_sjx(:, :), Eavg_sjx(:, :)
        real(r8), allocatable :: angle_wbx(:, :, :), angle_lbx(:, :, :), angle_qbx(:, :, :)
        real(r8), allocatable :: Savg_wbx(:), Savg_lbx(:), Savg_qbx(:)
        real(r8), allocatable :: Eavg_wbx(:, :), Eavg_lbx(:, :), Eavg_qbx(:, :)
        real(r8), allocatable :: Extr_wbx(:, :), Extr_lbx(:, :), Extr_qbx(:, :)

        real(r8), allocatable :: MoveDis(:, :)                 ! Record the adjustment distance of point w in the x and y directions

        real(r8) :: sjx(3, 2), newsjx(3, 2), dbx(7, 2), pi
        real(r8) :: rx, ry, fra
        integer, allocatable :: n_landtypes(:) 

        integer, allocatable :: ref(:)          !  Initial refinement of triangular mesh
	     integer, allocatable :: ref_l(:, :)     !  Polygon mesh refinement
	     integer, allocatable :: ref_th(:, :)    !  Initial triangular mesh threshold refinement
	     integer, allocatable :: ref_tr(:, :)    !  The final triangular mesh threshold refinement
	     integer, allocatable :: ref_pl(:, :)    !  Record what thresholds the polygon is refined by

	     integer, allocatable :: ngrmw(:, :), ngrwm(:, :)          !  Index of w/m points adjacent to m/w points (before refinement)
	     integer, allocatable :: ngrmw_new(:, :), ngrwm_new(:, :)  !  Index of w/m points adjacent to m/w points (after refinement)
	     integer, allocatable :: ngrmw_f(:, :), ngrwm_f(:, :)      !  Index of w/m points adjacent to m/w points (final)
        integer, allocatable :: ngrmw_f_tmp(:, :), ngrwm_f_tmp(:, :)
        integer, allocatable :: ngrmm(:, :)        !  Index of m points adjacent to m points (before thinning)
        integer, allocatable :: ngrmm_new(:, :)    !  Index of m points adjacent to m points (after refinement)
        integer, allocatable :: mrl(:)             !  Triangular mesh refinement degree (before refinement)
        integer, allocatable :: mrl_new(:)         !  Triangular mesh refinement degree (after refinement)
        integer, allocatable :: mrl_f(:)           !  Triangular mesh refinement (final)
        integer, allocatable :: ngr_mrl(:, :)      !  mrl of adjacent triangular mesh points of triangular mesh (before thinning)
        integer, allocatable :: ngr_mrl_new(:, :)  !  mrl of adjacent triangular mesh points of triangular mesh (after refinement)
        integer, allocatable :: mp_ref(:)          !  Record the index of the triangle being refined
        character(LEN = 256) :: lndname, nxpc

	logical :: isexist                !  Determine whether there are repeated w points after refinement

	! To prevent the collision of the intersection belt, one is divided into four
	logical :: iterA                  !  When iteration B and iteration C pass at the same time, iteration A passes
	logical :: iterB                  !  Judging from the triangular grid
	logical :: iterC                  !  Judging from the polygonal mesh
	logical :: End_SpringAjustment    !  Determine whether the mesh quality adjustment is complete
        logical,allocatable :: IsInRfArea(:) ! Determine whether the triangular mesh is in the refinement area 
        real(r8) :: tmpa(2),tmpb(2),tmpc(2),tmp_angle(3),tmp_length(3),tmp_angle7(7)
        character(LEN = 20) :: p_name(6) = (/"GLONW", "GLATW", "GLONM", "GLATM", "itab_w%im", "itab_m%iw"/)

        pi = 3.1415926535

        iterA = .false.
        iterB = .false.
        iterC = .false.

        !-------------------------------------------
        ! read unstructure mesh
        !-------------------------------------------
        print*, "start to read unstructure mesh data"
        print*, ""

        write(nxpc, '(I3.3)') NXP
        lndname = trim(base_dir) // trim(EXPNME) // '/makegrid/gridfile/gridfile_NXP' // trim(nxpc) // '.nc4'
        print*,lndname
        CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid))
        CALL CHECK(NF90_INQ_VARID(ncid, p_name(1), varid(1)))
        CALL CHECK(NF90_INQ_VARID(ncid, p_name(2), varid(2)))
        CALL CHECK(NF90_INQ_VARID(ncid, p_name(3), varid(3)))
        CALL CHECK(NF90_INQ_VARID(ncid, p_name(4), varid(4)))
        CALL CHECK(NF90_INQ_VARID(ncid, p_name(5), varid(5)))
        CALL CHECK(NF90_INQ_VARID(ncid, p_name(6), varid(6)))
        CALL CHECK(NF90_INQ_DIMID(ncid, "sjx_points", dimID_sjx))
        CALL CHECK(NF90_INQ_DIMID(ncid, "lbx_points", dimID_lbx))
        CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lbx, len = lbx_points))
        CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_sjx, len = sjx_points))

        print*, "sjx_points = ", sjx_points
        print*, "lbx_points = ", lbx_points

        allocate(wp(lbx_points, 2))            !  Initial data for center point of polygon mesh Initial data for center point of Polygon mesh
	allocate(wp_new(lbx_points * 4, 2))    !  Update data at center point of polygon mesh Update data at center point of Polygon mesh
	allocate(mp(sjx_points, 2))            !  Initial data at the center point of the triangular grid Initial data at the center point of the triangular grid
	allocate(mp_new(sjx_points * 4, 2))    !  The center point of the triangular grid updates the data The center point of the triangular grid updates the data
	allocate(ngrwm(7, lbx_points))         !  wp initial index table of adjacent mp points wp Initial index table of adjacent mp points
	allocate(ngrwm_new(7, lbx_points * 4)) !  wp's adjacent mp points update the index table wp's adjacent mp points update the index table
	allocate(ngrmw(3, sjx_points))         !  The initial index table of adjacent wp points of mp The initial index table of adjacent wp points of mp
	allocate(ngrmw_new(3, sjx_points * 4)) !  The adjacent wp points of mp update the index table The adjacent wp points of mp update the index table
        wp = 0.
        mp = 0.
        wp_new = 0
        mp_new = 0
        ngrmw = 1
        ngrwm = 1
        ngrmw_new = 1
        ngrwm_new = 1

        CALL CHECK(NF90_GET_VAR(ncid, varid(1), wp(:, 1)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(2), wp(:, 2)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(3), mp(:, 1)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(4), mp(:, 2)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(5), ngrwm(1:7, :)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(6), ngrmw(1:3, :)))
        CALL CHECK(NF90_CLOSE(ncid))

        print*, "The unstructured grid data reading have done "
        print*, ""
        print*, "In total, triangular mesh number: ", sjx_points, "polygon mesh number: ", lbx_points
        print*, ""

        allocate(mrl(sjx_points))                  !  Triangle mesh refinement degree Triangle mesh refinement degree
	allocate(mrl_new(sjx_points * 4))          !  Updated mrl data Updated MRL data
	allocate(ngr_mrl(3, sjx_points))           !  mrl of adjacent triangular mesh points of a triangular mesh mrl of adjacent triangular mesh points of a triangular mesh
	allocate(ngr_mrl_new(3, sjx_points * 4))   !  Updated ngr_mrl data (Updated data)
	allocate(ngrmm(3, sjx_points))             !  mp adjacent mp initial index table mp adjacent mp Initial Index Table
	allocate(ngrmm_new(3, sjx_points * 4))     !  Updated ngrmm data Updated NGRMM Data
        mrl = 1
        mrl_new = 1
        ngr_mrl = 1
        ngr_mrl_new = 1
        ngrmm = 1
        ngrmm_new = 1

        do i = 1, sjx_points, 1
            if(mp(i, 1) == -180.)then
                mp(i, 1) = 180.
            end if
        end do

        do i = 1, lbx_points, 1
            if(wp(i, 1) == -180.)then
                wp(i, 1) = 180.
            end if
        end do

        !----------------------------------------------------------------
        ! ngrmm calculation
        !----------------------------------------------------------------
        do i = 1, sjx_points, 1
            if(ngrmw(3, i) == 1)then
                mrl(i) = 0
            end if
        end do

        do i = 1, sjx_points, 1
            do j = 1, sjx_points, 1
                if(i /= j)then
                    error = IsNgrmm(ngrmw(1:3, i), ngrmw(1:3, j))
                    if(error /= 0)then
                        ngrmm(error, i) = j
                    end if
                end if
            end do
        end do

        ! save ngrmm
        !lndname = '/tera01/fanhw21/olam/Ouhe/output/MAKEGRID/tool/ngrmm.nc4'
        !print*,lndname
        !CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber,nf90_netcdf4), ncID))
        !CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points"  , sjx_points , spDimID))
        !CALL CHECK(NF90_DEF_DIM(ncID, "dim"         , 3          , thDimID))
        !CALL CHECK(NF90_DEF_VAR(ncID, "ngrmm", NF90_INT, (/ thDimID, spDimID /), ncVarID(1)))
        !CALL CHECK(NF90_ENDDEF(ncID))
        !CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1), ngrmm))
        !CALL CHECK(NF90_CLOSE(ncID))

        do i = 1, sjx_points, 1
            do j = 1, 3, 1
                if(ngrmm(j, i) /= 1)then
                    if(mrl(ngrmm(j, i)) == 0)then
                        ngr_mrl(j, i) = 0
                    end if
                end if
            end do
        end do

        ! Allocate threshold dependent array memory (Allocate threshold dependent array memory)
        allocate(n_landtypes(sjx_points))
        allocate(f_mainarea(sjx_points, 2))
        allocate(slope(sjx_points, 2))
        allocate(lai(sjx_points, 2))
        allocate(k_s(sjx_points, 2, 2))
        allocate(k_sl(sjx_points, 2, 2))
        allocate(tkdry(sjx_points, 2, 2))
        allocate(tksatf(sjx_points, 2, 2))
        allocate(tksatu(sjx_points, 2, 2))
        n_landtypes = 0
        f_mainarea = 0.
        slope = 0.
        lai = 0.
        k_s = 0.
        k_sl = 0.
        tkdry = 0.
        tksatf = 0.
        tksatu = 0.

        !--------------------------------------------------
        ! Reading threshold data (latitude and longitude grid) (Reading threshold data (latitude and longitude grid))
        !--------------------------------------------------
        lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/threshold/threshold_" // trim(nxpc) // ".nc4"
        print*, lndname
        CALL CHECK(NF90_OPEN(lndname, ior(nf90_clobber, nf90_netcdf4), ncid))
        CALL CHECK(NF90_INQ_VARID(ncid, "num_landtypes", varid(1)))
        CALL CHECK(NF90_INQ_VARID(ncid, "fraction_mainarea", varid(2)))
        CALL CHECK(NF90_INQ_VARID(ncid, "p_slope", varid(3)))
        CALL CHECK(NF90_INQ_VARID(ncid, "p_lai", varid(4)))
        CALL CHECK(NF90_INQ_VARID(ncid, "p_k_s", varid(5)))
        CALL CHECK(NF90_INQ_VARID(ncid, "p_k_sl", varid(6)))
        CALL CHECK(NF90_INQ_VARID(ncid, "p_tkdry", varid(7)))
        CALL CHECK(NF90_INQ_VARID(ncid, "p_tksatf", varid(8)))
        CALL CHECK(NF90_INQ_VARID(ncid, "p_tksatu", varid(9)))

        CALL CHECK(NF90_GET_VAR(ncid, varid(1), n_landtypes))
        CALL CHECK(NF90_GET_VAR(ncid, varid(2), f_mainarea))
        CALL CHECK(NF90_GET_VAR(ncid, varid(3), slope))
        CALL CHECK(NF90_GET_VAR(ncid, varid(4), lai))
        CALL CHECK(NF90_GET_VAR(ncid, varid(5), k_s))
        CALL CHECK(NF90_GET_VAR(ncid, varid(6), K_sl))
        CALL CHECK(NF90_GET_VAR(ncid, varid(7), tkdry))
        CALL CHECK(NF90_GET_VAR(ncid, varid(8), tksatf))
        CALL CHECK(NF90_GET_VAR(ncid, varid(9), tksatu))
        CALL CHECK(NF90_CLOSE(ncid))

        allocate(ref(sjx_points))             !  Initial refinement of triangular mesh Initial refinement of triangular mesh
	allocate(ref_l(lbx_points, 7))        !  Polygon mesh refinement
	allocate(ref_lbx(lbx_points, 8))
	allocate(ref_th(sjx_points, 16))      !  Initial triangular mesh threshold refinement Initial triangular mesh threshold refinement
	allocate(ref_tr(sjx_points * 4, 16))  !  The final triangular mesh threshold refinement The final triangular mesh threshold refinement        
	ref = 0
        ref_l = 0
        ref_lbx = 0.
        ref_th = 0
        ref_tr = 0

        allocate(IsInRfArea(sjx_points))
        IsInRfArea = .false.
        CALL IsInRefineArea(IsInRfArea,sjx_points,lbx_points,ngrmw,wp)

        !if(lcs == "igbp")then
        !   maxlc = 17
        !else
        !   maxlc = 24
        !end if

        !--------------------------------------------------
        ! 1.1 Used to record triangles that require preliminary refinement (Records require preliminary refinement of triangles)
        !--------------------------------------------------
        do i = 1, sjx_points, 1
            !if((wp(ngrmw(i,1),2)>85.).or.(wp(ngrmw(i,2),2)>85.).or.(wp(ngrmw(i,3),2)>85.))then
            !   cycle
            !end if
            if(IsInRfArea(i) == .false.)then
               cycle
            end if
            
            if((f_mainarea(i, 1) == 0.).or.(f_mainarea(i, 1) == maxlc))then
                cycle
            end if

            if (refine_num_landtypes == .True. .and. n_landtypes(i)>th_num_landtypes) then
                ref(i) = 1
                ref_th(i, 1) = 1
            end if

            if (refine_area_mainland == .True. .and. f_mainarea(i, 2) < th_area_mainland) then
                ref(i) = 1
                ref_th(i, 2) = 1
            end if

            if (refine_lai_m == .True. .and. lai(i, 1) > th_lai_m) then
                ref(i) = 1
                ref_th(i, 3) = 1
            end if

            if (refine_lai_s == .True. .and. lai(i, 2) > th_lai_s) then
                ref(i) = 1
                ref_th(i, 4) = 1
            end if

            if (refine_slope_m == .True. .and. slope(i, 1) > th_slope_m) then
                ref(i) = 1
                ref_th(i, 5) = 1
            end if

            if (refine_slope_s == .True. .and. slope(i, 2) > th_slope_s) then
                ref(i) = 1
                ref_th(i, 6) = 1
            end if

            if (refine_k_s_m == .True. .and. ((k_s(i, 1, 1) > th_k_s_m).or.(k_s(i, 1, 2) > th_k_s_m))) then
                ref(i) = 1
                ref_th(i, 7) = 1
            end if

            if (refine_k_s_s == .True. .and. ((k_s(i, 2, 1) > th_k_s_s).or.(k_s(i, 2, 2) > th_k_s_s))) then
                ref(i) = 1
                ref_th(i, 8) = 1
            end if

            if (refine_k_solids_m == .True. .and. ((k_sl(i, 1, 1) > th_k_solids_m).or.(k_sl(i, 1, 2) > th_k_solids_m))) then
                ref(i) = 1
                ref_th(i, 9) = 1
            end if

            if (refine_k_solids_s == .True. .and. ((k_sl(i, 2, 1) > th_k_solids_s).or.(k_sl(i, 2, 2) > th_k_solids_s))) then
                ref(i) = 1
                ref_th(i, 10) = 1
            end if

            if (refine_tkdry_m == .True. .and. ((tkdry(i, 1, 1) > th_tkdry_m).or.(tkdry(i, 1, 2) > th_tkdry_m))) then
                ref(i) = 1
                ref_th(i, 11) = 1
            end if

            if (refine_tkdry_s == .True. .and. ((tkdry(i, 2, 1) > th_tkdry_s).or.(tkdry(i, 2, 2) > th_tkdry_s))) then
                ref(i) = 1
                ref_th(i, 12) = 1
            end if

            if (refine_tksatf_m == .True. .and. ((tksatf(i, 1, 1) > th_tksatf_m).or.(tksatf(i, 1, 2) > th_tksatf_m))) then
                ref(i) = 1
                ref_th(i, 13) = 1
            end if

            if (refine_tksatf_s == .True. .and. ((tksatf(i, 2, 1) > th_tksatf_s).or.(tksatf(i, 2, 2) > th_tksatf_s))) then
                ref(i) = 1
                ref_th(i, 14) = 1
            end if

            if (refine_tksatu_m == .True. .and. ((tksatu(i, 1, 1) > th_tksatu_m).or.(tksatu(i, 1, 2) > th_tksatu_m))) then
                ref(i) = 1
                ref_th(i, 15) = 1
            end if

            if (refine_tksatu_s == .True. .and. ((tksatu(i, 2, 1) > th_tksatu_s).or.(tksatu(i, 2, 2) > th_tksatu_s))) then
                ref(i) = 1
                ref_th(i, 16) = 1
            end if

        end do

        iter = 1                                 ! Number of Iterations (Number of iterations)
        num_ref = INT(sum(ref))                  ! Number of triangles per refinement (Refine the number of triangles each time)
        nmp(1) = sjx_points + 4 * num_ref        ! Record the number of triangles after each iteration (Record the number of triangles after each iteration)
        nwp(1) = lbx_points + 3 * num_ref        ! Record the number of polygons after each iteration (Record the number of polygons after each iteration)

        mrl_new(1:sjx_points) = mrl
        ngr_mrl_new(1:3, 1:sjx_points) = ngr_mrl
        mp_new(1:sjx_points, 1:2) = mp
        wp_new(1:lbx_points, 1:2) = wp
        ngrmw_new(1:3, 1:sjx_points) = ngrmw(1:3, 1:sjx_points)
        ngrmm_new(1:3, 1:sjx_points) = ngrmm(1:3, 1:sjx_points)
        ngrwm_new(1:7, 1:lbx_points) = ngrwm(1:7, 1:lbx_points)

        deallocate(n_landtypes)
        deallocate(f_mainarea)
        deallocate(slope)
        deallocate(lai)
        deallocate(k_s)
        deallocate(k_sl)
        deallocate(tkdry)
        deallocate(tksatf)
        deallocate(tksatu)

        !--------------------------------------------------
        ! 1.2 Preliminary refinement (one into four) 【 Preliminary refinement (divided into four) 】
        !--------------------------------------------------
        print*, "Start to refine (Start preliminary refinement)"
        print*, "iter =", iter, "num =", num_ref
        refed = 0
        do i = 1, sjx_points, 1
            if(ref(i) == 1)then
                icl = 0
                sjx = 0.
                newsjx = 0.
                sjx(1, 1:2) = wp(ngrmw(1, i), 1:2)
                sjx(2, 1:2) = wp(ngrmw(2, i), 1:2)
                sjx(3, 1:2) = wp(ngrmw(3, i), 1:2)

                icl(1) = IsCrossLine(sjx(2, 1), sjx(3, 1))
                icl(2) = IsCrossLine(sjx(1, 1), sjx(3, 1))
                icl(3) = IsCrossLine(sjx(1, 1), sjx(2, 1))

                if(INT(sum(icl)) > 0)then
                    do j = 1, 3, 1
                        if(sjx(j, 1) < 0.)then
                            sjx(j, 1) = sjx(j, 1) + 360.
                        end if
                    end do
                end if

                newsjx(1, 1:2) = (sjx(2, 1:2) + sjx(3, 1:2)) / 2.
                newsjx(2, 1:2) = (sjx(1, 1:2) + sjx(3, 1:2)) / 2.
                newsjx(3, 1:2) = (sjx(1, 1:2) + sjx(2, 1:2)) / 2.

                m1 = sjx_points + refed(1) * 4 + 1
                m2 = sjx_points + refed(1) * 4 + 2
                m3 = sjx_points + refed(1) * 4 + 3
                m4 = sjx_points + refed(1) * 4 + 4

                w1 = lbx_points + refed(1) * 3 + 1
                w2 = lbx_points + refed(1) * 3 + 2
                w3 = lbx_points + refed(1) * 3 + 3

                mp_new(m1, 1:2) = (sjx(1, :) + newsjx(2, :) + newsjx(3, :)) / 3.
                mp_new(m2, 1:2) = (sjx(2, :) + newsjx(1, :) + newsjx(3, :)) / 3.
                mp_new(m3, 1:2) = (sjx(3, :) + newsjx(1, :) + newsjx(2, :)) / 3.
                mp_new(m4, 1:2) = (newsjx(3, :) + newsjx(1, :) + newsjx(2, :)) / 3.

                wp_new(w1, 1:2) = newsjx(1, 1:2)
                wp_new(w2, 1:2) = newsjx(2, 1:2)
                wp_new(w3, 1:2) = newsjx(3, 1:2)

                do k = 1, 3, 1
                    ngrmm_new(k, i) = 1
                    ngrmw_new(k, i) = 1
                    ngrmw_new(1, sjx_points + refed(1) * 4 + k) = ngrmw(k, i)
                    ngrmm_new(1, sjx_points + refed(1) * 4 + k) = m4
                end do

                ngrmw_new(2, m1) = w3
                ngrmw_new(3, m1) = w2
                ngrmw_new(2, m2) = w1
                ngrmw_new(3, m2) = w3
                ngrmw_new(2, m3) = w2
                ngrmw_new(3, m3) = w1
                ngrmw_new(1, m4) = w1
                ngrmw_new(2, m4) = w2
                ngrmw_new(3, m4) = w3
                ngrmm_new(1, m4) = m1
                ngrmm_new(2, m4) = m2
                ngrmm_new(3, m4) = m3

                if(INT(sum(icl)) > 0)then
                    do k = 1, 4, 1
                        call CheckLon(mp_new(sjx_points + refed(1) * 4 + k, 1))
                        call CheckLon(wp_new(lbx_points + refed(1) * 3 + k, 1))
                    end do
                end if

                mrl_new(i) = 4 ! Indicates that the triangular grid is evenly divided into four parts (Indicates that the triangular grid is evenly divided into four parts)
                mrl_new(sjx_points + refed(iter) * 4 + 1:sjx_points + refed(1) * 4 + 4) = 4

                do j = 1, 3, 1
                    do k = 1, 3, 1
                        if(ngrmm(k, ngrmm(j, i)) == i)then
                            ngr_mrl_new(k, ngrmm(j, i)) = 4
                        end if
                    end do
                end do
                ngr_mrl_new(1:3, m4) = 4

                ref_lbx(8, ngrmw_new(1, i)) = 1
                ref_lbx(8, ngrmw_new(2, i)) = 1
                ref_lbx(8, ngrmw_new(3, i)) = 1

                do j = 1, 16, 1
                    if(ref_th(i, j) == 1)then
                        ref_tr(m1, j) = 1
                        ref_tr(m2, j) = 1
                        ref_tr(m3, j) = 1
                    end if
                end do

                refed(1) = refed(1) + 1
            end if
        end do
        print*, "itered_num =", refed(1)
        print*, "Refining step 1 is complete (Refining step 1 completed)"

        !--------------------------------------------------
        ! 1.3 Store the initial grid data and the new grid data after preliminary refinement
        !--------------------------------------------------
        lndname = trim(base_dir) // trim(EXPNME)  // "/makegrid/tmp/gridfile_NXP" // trim(nxpc) // "_01_ori.nc4"
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", sjx_points, spDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", lbx_points, lpDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 3, thDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_b", 7, seDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONM", NF90_FLOAT, (/ spDimID /), ncVarID(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATM", NF90_FLOAT, (/ spDimID /), ncVarID(2)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONW", NF90_FLOAT, (/ lpDimID /), ncVarID(3)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATW", NF90_FLOAT, (/ lpDimID /), ncVarID(4)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_m%iw", NF90_INT, (/ thDimID, spDimID /), ncVarID(5)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_w%im", NF90_INT, (/ seDimID, lpDimID /), ncVarID(6)))
        CALL CHECK(NF90_ENDDEF(ncID))

        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1), mp(:, 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(2), mp(:, 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(3), wp(:, 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(4), wp(:, 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5), ngrmw(:, :)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6), ngrwm(:, :)))
        CALL CHECK(NF90_CLOSE(ncID))

        lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/gridfile_NXP" // trim(nxpc) // "_01_1.nc4"
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", nmp(1), spDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", nwp(1), lpDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 3, thDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_b", 7, seDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONM", NF90_FLOAT, (/ spDimID /), ncVarID(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATM", NF90_FLOAT, (/ spDimID /), ncVarID(2)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONW", NF90_FLOAT, (/ lpDimID /), ncVarID(3)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATW", NF90_FLOAT, (/ lpDimID /), ncVarID(4)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_m%iw", NF90_INT, (/ thDimID, spDimID /), ncVarID(5)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_w%im", NF90_INT, (/ seDimID, lpDimID /), ncVarID(6)))
        CALL CHECK(NF90_ENDDEF(ncID))

        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1), mp_new(1:nmp(1), 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(2), mp_new(1:nmp(1), 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(3), wp_new(1:nwp(1), 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(4), wp_new(1:nwp(1), 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5), ngrmw_new(1:3, 1:nmp(1))))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6), ngrwm_new(1:7, 1:nwp(1))))
        CALL CHECK(NF90_CLOSE(ncID))

        !stop

        !--------------------------------------------------
        ! 2.1 Perform iterations (to prevent conflicts in the refinement intersection belt, one is divided into four)
        !--------------------------------------------------
        print*, "iterA start"
        do while(iterA == .false.)

            iterA = .true.    ! Determine whether iteration B and iteration C meet the conditions
            iterB = .false.   ! Judging from the triangular grid
            iterC = .false.   ! Judging from the polygonal mesh

            print*, "iterB start"
            do while(iterB == .false.)

                ref = 0

                ! When an unrefined triangle has two or more refined adjacent triangles, the triangle is refined
                do i = 1, sjx_points, 1
                    if(mrl_new(i) == 1)then
                        if((sum(ngr_mrl_new(1:3, i)) == 12.).or.(sum(ngr_mrl_new(1:3, i)) == 9.))then
                            ref(i) = 1
                        end if
                    end if
                end do

                num_ref = INT(sum(ref))

                if(num_ref == 0)then
                    iterB = .true.
                else
                    iterA = .false.
                    iter = iter + 1
                    nmp(iter) = nmp(iter - 1) + 4 * num_ref
                    nwp(iter) = nwp(iter - 1) + 3 * num_ref

                    print*, "iter =", iter, "num =", num_ref

                    do i = 1, sjx_points, 1
                        if((ref(i) == 1).and.(mrl_new(i) == 1))then
                            icl = 0
                            sjx(1, 1:2) = wp(ngrmw(1, i), 1:2)
                            sjx(2, 1:2) = wp(ngrmw(2, i), 1:2)
                            sjx(3, 1:2) = wp(ngrmw(3, i), 1:2)

                            icl(1) = IsCrossLine(sjx(2, 1), sjx(3, 1))
                            icl(2) = IsCrossLine(sjx(1, 1), sjx(3, 1))
                            icl(3) = IsCrossLine(sjx(1, 1), sjx(2, 1))

                            if(INT(sum(icl)) > 0.5)then
                                do j = 1, 3, 1
                                    if(sjx(j, 1) < 0)then
                                        sjx(j, 1) = sjx(j, 1) + 360.
                                    end if
                                end do
                            end if

                            newsjx(1, 1:2) = (sjx(2, 1:2) + sjx(3, 1:2)) / 2.
                            newsjx(2, 1:2) = (sjx(1, 1:2) + sjx(3, 1:2)) / 2.
                            newsjx(3, 1:2) = (sjx(1, 1:2) + sjx(2, 1:2)) / 2.

                            m1 = nmp(iter - 1) + refed(iter) * 4 + 1
                            m2 = nmp(iter - 1) + refed(iter) * 4 + 2
                            m3 = nmp(iter - 1) + refed(iter) * 4 + 3
                            m4 = nmp(iter - 1) + refed(iter) * 4 + 4

                            mp_new(m1, 1:2) = (sjx(1, :) + newsjx(2, :) + newsjx(3, :)) / 3.
                            mp_new(m2, 1:2) = (sjx(2, :) + newsjx(1, :) + newsjx(3, :)) / 3.
                            mp_new(m3, 1:2) = (sjx(3, :) + newsjx(1, :) + newsjx(2, :)) / 3.
                            mp_new(m4, 1:2) = (newsjx(3, :) + newsjx(1, :) + newsjx(2, :)) / 3.

                            w1 = nwp(iter - 1) + refed(iter) * 3 + 1
                            w2 = nwp(iter - 1) + refed(iter) * 3 + 2
                            w3 = nwp(iter - 1) + refed(iter) * 3 + 3

                            wp_new(w1, 1:2) = newsjx(1, 1:2)
                            wp_new(w2, 1:2) = newsjx(2, 1:2)
                            wp_new(w3, 1:2) = newsjx(3, 1:2)

                            do k = 1, 3, 1
                                ngrmm_new(k, i) = 1
                                ngrmw_new(k, i) = 1
                                ngrmw_new(1, nmp(iter - 1) + refed(iter) * 4 + k) = ngrmw(k, i)
                                ngrmm_new(1, nmp(iter - 1) + refed(iter) * 4 + k) = m4
                            end do

                            ngrmw_new(2, m1) = w3
                            ngrmw_new(3, m1) = w2
                            ngrmw_new(2, m2) = w1
                            ngrmw_new(3, m2) = w3
                            ngrmw_new(2, m3) = w2
                            ngrmw_new(3, m3) = w1
                            ngrmw_new(1, m4) = w1
                            ngrmw_new(2, m4) = w2
                            ngrmw_new(3, m4) = w3

                            if(INT(sum(icl)) > 0)then
                                do k = 1, 4, 1
                                    call CheckLon(mp_new(nmp(iter - 1) + refed(iter) * 4 + k, 1))
                                    call CheckLon(wp_new(nwp(iter - 1) + refed(iter) * 3 + k, 1))
                                end do
                            end if

                            mrl_new(i) = 4
                            mrl_new(nmp(iter - 1) + refed(iter) * 4 + 1:nmp(iter - 1) + refed(iter) * 4 + 4) = 4

                            do j = 1, 3, 1
                                do k = 1, 3, 1
                                    if(ngrmm(k, ngrmm(j, i)) == i)then
                                        ngr_mrl_new(k, ngrmm(j, i)) = 4
                                    end if
                                end do
                            end do
                            ngr_mrl_new(1:3, m4) = 4

                            !do j = 1,7,1
                            !   do k = 1,3,1
                            !      if(ngrwm_new(ngrmw(i,k),j) == i)then
                            !         ngrwm_new(ngrmw(i,k),j) = nmp(iter-1)+refed(iter)*3+k
                            !      end if
                            !   end do
                            !end do

                            ref_lbx(ngrmw(1, i), 8) = 1
                            ref_lbx(ngrmw(2, i), 8) = 1
                            ref_lbx(ngrmw(3, i), 8) = 1

                            refed(iter) = refed(iter) + 1

                        end if
                    end do
                end if
            end do ! iter B
            print*, "iterB end"

            print*, "iterC start"
            do while(iterC == .false.)

                ref = 0

                do i = 1, lbx_points, 1

                    edges = 0    

                    do j = 1, 7, 1   ! Calculate the number of sides of a polygon made of unrefined triangles
                        if(ngrwm(j, i) /= 1)then
                            m1 = ngrwm(j, i)
                            if(sum(ngr_mrl_new(:, m1)) == 6.)then
                                ref_lbx(i, j) = 1  ! Indicates that the triangles that make up part of the polygon have been refined
                            end if
                            edges = edges + 1
                        end if
                    end do

                    if(ref_lbx(i, 8) == 0)then     ! The triangles that make up the polygon have not been refined
                        if(edges == 5)then        
                            do j = 1, 5, 1
                                m1 = ngrwm(j, i)
                                m2 = ngrwm(j + 1, i)
                                if(j == 5)then
                                    m2 = ngrwm(1, i)
                                end if
                                if((sum(ngr_mrl_new(:, m1)) == 6).and.(sum(ngr_mrl_new(:, m2)) == 6.))then
                                    ref_lbx(i, j) = 0.5
                                    if(j < 5)then
                                        ref_lbx(i, j + 1) = 0.5
                                    else
                                        ref_lbx(i, 1) = 0.5
                                    end if
                                end if
                            end do
                        else if(edges == 6)then   
                            do j = 1, 6, 1
                                m1 = ngrwm(j, i)
                                m2 = ngrwm(j + 1, i)
                                if(j == 6)then
                                    m2 = ngrwm(1, i)
                                end if
                                if((sum(ngr_mrl_new(:, m1)) == 6).and.(sum(ngr_mrl_new(:, m2)) == 6.))then
                                    ref_lbx(i, j) = 0.5
                                    if(j < 6)then
                                        ref_lbx(i, j + 1) = 0.5
                                    else
                                        ref_lbx(i, 1) = 0.5
                                    end if
                                end if
                            end do
                        end if

                        if((edges == 5).and.(sum(ref_lbx(i, 1:7)) > 2.))then
                            do j = 1, 7, 1
                                if((ref_lbx(i, j) /= 0).and.(mrl_new(ngrwm(j, i)) == 1))then
                                    ref(ngrwm(j, i)) = 1
                                end if
                            end do
                        else if((edges == 6).and.(sum(ref_lbx(i, 1:7)) > 1.))then
                            do j = 1, 7, 1
                                if((ref_lbx(i, j) /= 0).and.(mrl_new(ngrwm(j, i)) == 1))then
                                    ref(ngrwm(j, i)) = 1
                                end if
                            end do
                        end if

                    else if(ref_lbx(i, 8) /= 0)then
                        if(edges == 5)then
                            if(sum(mrl_new(ngrwm(1:5, i))) > 10)then
                                do j = 1, 5, 1
                                    if(mrl_new(ngrwm(j, i)) == 1)then
                                        ref(ngrwm(j, i)) = 1
                                    end if
                                end do
                            end if
                        else if(edges == 6)then
                            if(sum(mrl_new(ngrwm(1:6, i))) == 12)then
                                do j = 1, 3, 1
                                    if((mrl_new(ngrwm(j, i)) == 4).and.(mrl_new(ngrwm(j + 3, i)) == 4))then
                                        if((mrl_new(ngrwm(j + 1, i)) == 1).and.(mrl_new(ngrwm(j + 1, i)) == 1))then
                                            ref(ngrwm(j + 1, i)) = 1
                                            ref(ngrwm(j + 2, i)) = 1
                                        end if
                                    end if
                                end do
                            else if(sum(mrl_new(ngrwm(1:6, i))) == 9)then
                                k = 0
                                do j = 1, 6, 1
                                    l = ngrwm(j, i)
                                    if(mrl_new(l) == 1)then
                                        if(sum(mrl_new(ngrmm(1:3, l))) == 6)then
                                            k = k + 1
                                        end if
                                    end if
                                end do
                                if(k > 3)then
                                    do j = 1, 6, 1
                                        l = ngrwm(j, i)
                                        if(mrl_new(l) == 1)then
                                            ref(l) = 1
                                        end if
                                    end do
                                end if
                            end if
                        end if
                    end if

                end do

                num_ref = INT(sum(ref))

                if(num_ref == 0)then
                    iterC = .true.
                else
                    iterA = .false.
                    iter = iter + 1
                    nmp(iter) = nmp(iter - 1) + 4 * num_ref
                    nwp(iter) = nwp(iter - 1) + 3 * num_ref

                    print*, "iter =", iter, "num =", num_ref

                    do i = 1, sjx_points, 1
                        if((ref(i) == 1).and.(mrl_new(i) == 1))then
                            icl = 0
                            sjx(1, 1:2) = wp(ngrmw(1, i), 1:2)
                            sjx(2, 1:2) = wp(ngrmw(2, i), 1:2)
                            sjx(3, 1:2) = wp(ngrmw(3, i), 1:2)

                            icl(1) = IsCrossLine(sjx(2, 1), sjx(3, 1))
                            icl(2) = IsCrossLine(sjx(1, 1), sjx(3, 1))
                            icl(3) = IsCrossLine(sjx(1, 1), sjx(2, 1))

                            if(INT(sum(icl)) > 0.5)then
                                do j = 1, 3, 1
                                    if(sjx(j, 1) < 0)then
                                        sjx(j, 1) = sjx(j, 1) + 360.
                                    end if
                                end do
                            end if

                            newsjx(1, 1:2) = (sjx(2, 1:2) + sjx(3, 1:2)) / 2.
                            newsjx(2, 1:2) = (sjx(1, 1:2) + sjx(3, 1:2)) / 2.
                            newsjx(3, 1:2) = (sjx(1, 1:2) + sjx(2, 1:2)) / 2.

                            m1 = nmp(iter - 1) + refed(iter) * 4 + 1
                            m2 = nmp(iter - 1) + refed(iter) * 4 + 2
                            m3 = nmp(iter - 1) + refed(iter) * 4 + 3
                            m4 = nmp(iter - 1) + refed(iter) * 4 + 4

                            mp_new(m1, 1:2) = (sjx(1, :) + newsjx(2, :) + newsjx(3, :)) / 3.
                            mp_new(m2, 1:2) = (sjx(2, :) + newsjx(1, :) + newsjx(3, :)) / 3.
                            mp_new(m3, 1:2) = (sjx(3, :) + newsjx(1, :) + newsjx(2, :)) / 3.
                            mp_new(m4, 1:2) = (newsjx(3, :) + newsjx(1, :) + newsjx(2, :)) / 3.

                            w1 = nwp(iter - 1) + refed(iter) * 3 + 1
                            w2 = nwp(iter - 1) + refed(iter) * 3 + 2
                            w3 = nwp(iter - 1) + refed(iter) * 3 + 3

                            wp_new(w1, 1:2) = newsjx(1, 1:2)
                            wp_new(w2, 1:2) = newsjx(2, 1:2)
                            wp_new(w3, 1:2) = newsjx(3, 1:2)

                            do k = 1, 3, 1
                                ngrmm_new(k, i) = 1
                                ngrmw_new(k, i) = 1
                                ngrmw_new(1, nmp(iter - 1) + refed(iter) * 4 + k) = ngrmw(k, i)
                                ngrmm_new(1, nmp(iter - 1) + refed(iter) * 4 + k) = m4
                            end do

                            ngrmw_new(2, m1) = w3
                            ngrmw_new(3, m1) = w2
                            ngrmw_new(2, m2) = w1
                            ngrmw_new(3, m2) = w3
                            ngrmw_new(2, m3) = w2
                            ngrmw_new(3, m3) = w1
                            ngrmw_new(1, m4) = w1
                            ngrmw_new(2, m4) = w2
                            ngrmw_new(3, m4) = w3

                            if(INT(sum(icl)) > 0)then
                                do k = 1, 4, 1
                                    call CheckLon(mp_new(nmp(iter - 1) + refed(iter) * 4 + k, 1))
                                    call CheckLon(wp_new(nwp(iter - 1) + refed(iter) * 3 + k, 1))
                                end do
                            end if

                            mrl_new(i) = 4
                            mrl_new(nmp(iter - 1) + refed(iter) * 4 + 1:nmp(iter - 1) + refed(iter) * 4 + 4) = 4

                            do j = 1, 3, 1
                                do k = 1, 3, 1
                                    if(ngrmm(k, ngrmm(j, i)) == i)then
                                        ngr_mrl_new(k, ngrmm(j, i)) = 4
                                    end if
                                end do
                            end do
                            ngr_mrl_new(1:3, m4) = 4

                            !do j = 1,7,1
                            !   do k = 1,3,1
                            !      if(ngrwm_new(ngrmw(i,k),j) == i)then
                            !         ngrwm_new(ngrmw(i,k),j) = nmp(iter-1)+refed(iter)*3+k
                            !      end if
                            !   end do
                            !end do

                            ref_lbx(ngrmw(1, i), 8) = 1
                            ref_lbx(ngrmw(2, i), 8) = 1
                            ref_lbx(ngrmw(3, i), 8) = 1

                            refed(iter) = refed(iter) + 1

                        end if
                    end do
                end if
            end do ! iterC
            print*, "iterC end"
        end do ! iterA
        print*, "iterA end"
        print*, "Refining step 2 is complete"

        !--------------------------------------------------
        ! 2.2 Store the new grid data refined in step 2
        !--------------------------------------------------
        lndname = trim(base_dir) // trim(EXPNME)  // "/makegrid/tmp/gridfile_NXP" // trim(nxpc) // "_01_2.nc4"
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", nmp(iter), spDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", nwp(iter), lpDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 3, thDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_b", 7, seDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONM", NF90_FLOAT, (/ spDimID /), ncVarID(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATM", NF90_FLOAT, (/ spDimID /), ncVarID(2)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONW", NF90_FLOAT, (/ lpDimID /), ncVarID(3)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATW", NF90_FLOAT, (/ lpDimID /), ncVarID(4)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_m%iw", NF90_INT, (/ thDimID, spDimID /), ncVarID(5)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_w%im", NF90_INT, (/ seDimID, lpDimID /), ncVarID(6)))
        CALL CHECK(NF90_ENDDEF(ncID))

        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1), mp_new(1:nmp(iter), 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(2), mp_new(1:nmp(iter), 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(3), wp_new(1:nwp(iter), 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(4), wp_new(1:nwp(iter), 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5), ngrmw_new(1:3, 1:nmp(iter))))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6), ngrwm_new(1:7, 1:nwp(iter))))
        CALL CHECK(NF90_CLOSE(ncID))

        !stop

        !--------------------------------------------------
        ! 3.1 Look for weak pits
        !--------------------------------------------------
        ref = 0
        iter = iter + 1
        do i = 1, sjx_points, 1
            if(mrl_new(i) == 1)then
                if(sum(ngr_mrl_new(1:3, i)) == 6.)then
                    do j = 1, 3, 1
                        if(mrl_new(ngrmm(j, i)) == 1)then
                            if(sum(ngr_mrl_new(1:3, ngrmm(j, i))) == 6.)then
                                do k = 1, 3, 1
                                    if(ngr_mrl_new(k, i) == 4)then
                                        m = k
                                    end if
                                    if(ngr_mrl_new(k, ngrmm_new(j, i)) == 4)then
                                        n = k
                                    end if
                                end do
                                do k = 1, 3, 1
                                    do l = 1, 3, 1
                                        if(ngrmw(k, ngrmm(m, i)) == ngrmw(l, ngrmm(n, ngrmm(j, i))))then
                                            ref(i) = 1
                                        end if
                                    end do
                                end do
                            end if
                        end if
                    end do
                end if
            end if
        end do

        num_ref = INT(sum(ref))
        print*, "Begin to refine the weak concave points"

        nmp(iter) = nmp(iter - 1) + 2 * num_ref
        nwp(iter) = nwp(iter - 1) + num_ref

        !--------------------------------------------------
        ! 3.2 Refine the weak concave points
        !--------------------------------------------------
        print*, "iter =", iter, "num =", num_ref

        do i = 1, sjx_points, 1
            if(ref(i) == 1)then
                j = 0
                do x = 1, 3, 1
                    if(ref(ngrmm(x, i)) == 1)then
                        j = ngrmm(x, i)
                    end if
                end do
                if(j == 0)then
                    cycle
                end if

                ref(i) = 0
                ref(j) = 0

                do x = 1, 3, 1
                    if(ngr_mrl_new(x, i) == 4)then
                        m = ngrmm(x, i)
                    end if
                    if(ngr_mrl_new(x, j) == 4)then
                        n = ngrmm(x, j)
                    end if
                end do

                do x = 1, 3, 1
                    do y = 1, 3, 1

                        if(ngrmw(x, i) == ngrmw(y, m))then
                            do z = 1, 3, 1
                                if(ngrmw(x, i) == ngrmw(z, j))then
                                    w2 = ngrmw(x, i)
                                end if
                            end do
                            if((ngrmw(x, i)/=ngrmw(1, j)).and.(ngrmw(x, i)/=ngrmw(2, j)).and.(ngrmw(x, i)/=ngrmw(3, j)))then
                                w1 = ngrmw(x, i)
                            end if
                        end if

                        if(ngrmw(x, j) == ngrmw(y, n))then
                            if((ngrmw(x, j)/=ngrmw(1, i)).and.(ngrmw(x, j)/=ngrmw(2, i)).and.(ngrmw(x, j)/=ngrmw(3, i)))then
                                w3 = ngrmw(x, j)
                            end if
                        end if

                        if(ngrmw(x, i) == ngrmw(y, j))then
                            if((ngrmw(x, i)/=ngrmw(1, m)).and.(ngrmw(x, i)/=ngrmw(2, m)).and.(ngrmw(x, i)/=ngrmw(3, m)))then
                                w4 = ngrmw(x, i)
                            end if
                        end if

                    end do
                end do

                mrl_new(i) = 0
                mrl_new(j) = 0

                w5 = nwp(iter - 1) + refed(iter) * 2 + 1
                w6 = nwp(iter - 1) + refed(iter) * 2 + 2

                m1 = nmp(iter - 1) + refed(iter) * 4 + 1
                m2 = nmp(iter - 1) + refed(iter) * 4 + 2
                m3 = nmp(iter - 1) + refed(iter) * 4 + 3
                m4 = nmp(iter - 1) + refed(iter) * 4 + 4

                icl = 0

                icl(1) = IsCrossLine(wp_new(w1, 1), wp_new(w2, 1))
                icl(2) = IsCrossLine(wp_new(w2, 1), wp_new(w3, 1))
                icl(3) = IsCrossLine(wp_new(w3, 1), wp_new(w4, 1))

                if(INT(sum(icl)) > 0.)then
                    call MoveLon(wp_new(w1, 1))
                    call MoveLon(wp_new(w2, 1))
                    call MoveLon(wp_new(w3, 1))
                    call MoveLon(wp_new(w4, 1))
                end if

                wp_new(w5, 1:2) = (wp_new(w1, 1:2) + wp_new(w2, 1:2)) / 2.
                wp_new(w6, 1:2) = (wp_new(w2, 1:2) + wp_new(w3, 1:2)) / 2.

                mp_new(m1, 1:2) = (wp_new(w1, 1:2) + wp_new(w4, 1:2) + wp_new(w5, 1:2)) / 3.
                mp_new(m2, 1:2) = (wp_new(w2, 1:2) + wp_new(w5, 1:2) + wp_new(w6, 1:2)) / 3.
                mp_new(m3, 1:2) = (wp_new(w3, 1:2) + wp_new(w4, 1:2) + wp_new(w6, 1:2)) / 3.
                mp_new(m4, 1:2) = (wp_new(w4, 1:2) + wp_new(w5, 1:2) + wp_new(w6, 1:2)) / 3.

                ngrmw_new(1, m1) = w4
                ngrmw_new(2, m1) = w1
                ngrmw_new(3, m1) = w5
                ngrmw_new(1, m2) = w2
                ngrmw_new(2, m2) = w5
                ngrmw_new(3, m2) = w6
                ngrmw_new(1, m3) = w4
                ngrmw_new(2, m3) = w6
                ngrmw_new(3, m3) = w3
                ngrmw_new(1, m4) = w4
                ngrmw_new(2, m4) = w5
                ngrmw_new(3, m4) = w6

                ngrmm_new(:, i) = 1
                ngrmm_new(:, j) = 1
                ngrmw_new(:, i) = 1
                ngrmw_new(:, j) = 1

                if(sum(icl) > 0)then
                    call CheckLon(wp_new(w1, 1))
                    call CheckLon(wp_new(w2, 1))
                    call CheckLon(wp_new(w3, 1))
                    call CheckLon(wp_new(w4, 1))
                    call CheckLon(wp_new(w5, 1))
                    call CheckLon(wp_new(w6, 1))
                    call CheckLon(mp_new(m1, 1))
                    call CheckLon(mp_new(m2, 1))
                    call CheckLon(mp_new(m3, 1))
                    call CheckLon(mp_new(m4, 1))
                end if

                refed(iter) = refed(iter) + 1
            end if
        end do

        print*, "The third step is to refine the total number of triangles", refed(iter) * 2
        print*, "Refining step 3 is complete"

        ref = 0
        iter = iter + 1

        !--------------------------------------------------
        ! 3.3 Store the new grid data refined in step 3
        !--------------------------------------------------
        lndname = trim(base_dir) // trim(EXPNME)  // "/makegrid/tmp/gridfile_NXP" // trim(nxpc) // "_01_3.nc4"
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", nmp(iter - 1), spDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", nwp(iter - 1), lpDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 3, thDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_b", 7, seDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONM", NF90_FLOAT, (/ spDimID /), ncVarID(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATM", NF90_FLOAT, (/ spDimID /), ncVarID(2)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONW", NF90_FLOAT, (/ lpDimID /), ncVarID(3)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATW", NF90_FLOAT, (/ lpDimID /), ncVarID(4)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_m%iw", NF90_INT, (/ thDimID, spDimID /), ncVarID(5)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_w%im", NF90_INT, (/ seDimID, lpDimID /), ncVarID(6)))
        CALL CHECK(NF90_ENDDEF(ncID))

        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1), mp_new(1:nmp(iter - 1), 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(2), mp_new(1:nmp(iter - 1), 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(3), wp_new(1:nwp(iter - 1), 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(4), wp_new(1:nwp(iter - 1), 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5), ngrmw_new(1:3, 1:nmp(iter - 1))))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6), ngrwm_new(1:7, 1:nwp(iter - 1))))
        CALL CHECK(NF90_CLOSE(ncID))

        !stop

        !--------------------------------------------------
        ! 4.1 Records are adjacent to only one refined triangle
        !--------------------------------------------------
        ref = 0
        do i = 1, sjx_points, 1
            if(mrl_new(i) == 1)then
                if(sum(ngr_mrl_new(:, i)) == 6)then
                    ref(i) = 1
                end if
            end if
        end do

        num_ref = INT(sum(ref))

        print*, "Begin to refine adjacent triangles Triangles with a preliminary refinement number of 1"
        print*, "iter =", iter, "num =", num_ref

        nmp(iter) = nmp(iter - 1) + 2 * num_ref
        nwp(iter) = nwp(iter - 1) + num_ref

        !--------------------------------------------------
        ! 4.2 Thinning Adjacent triangle with only one thinning triangle (split in two)
        !--------------------------------------------------
        do i = 1, sjx_points, 1
            if(ref(i) == 1)then

                do j = 1, 3, 1
                    if(mrl_new(ngrmm(j, i)) == 4)then
                        m2 = ngrmm(j, i)
                    end if
                end do

                do j = 1, 3, 1
                    if((ngrmw_new(j, i) /= ngrmw(1, m2)).and.(ngrmw_new(j, i) /= ngrmw(2, m2)).and.&
                            (ngrmw_new(j, i) /= ngrmw(3, m2)))then

                        if(j == 1)then
                            w1 = ngrmw_new(1, i)
                            w2 = ngrmw_new(2, i)
                            w3 = ngrmw_new(3, i)
                        else if(j == 2)then
                            w1 = ngrmw_new(2, i)
                            w2 = ngrmw_new(1, i)
                            w3 = ngrmw_new(3, i)
                        else
                            w1 = ngrmw_new(3, i)
                            w2 = ngrmw_new(1, i)
                            w3 = ngrmw_new(2, i)
                        end if

                        sjx(1, 1:2) = wp(w1, 1:2)
                        sjx(2, 1:2) = wp(w2, 1:2)
                        sjx(3, 1:2) = wp(w3, 1:2)

                    end if
                end do

                icl = 0

                icl(1) = IsCrossLine(sjx(2, 1), sjx(3, 1))
                icl(2) = IsCrossLine(sjx(1, 1), sjx(3, 1))
                icl(3) = IsCrossLine(sjx(1, 1), sjx(2, 1))

                if(INT(sum(icl)) > 0.)then
                    do j = 1, 3, 1
                        if(sjx(j, 1) < 0.)then
                            sjx(j, 1) = sjx(j, 1) + 360.
                        end if
                    end do
                end if

                newsjx(1, :) = (sjx(2, :) + sjx(3, :)) / 2.

                m1 = nmp(iter - 1) + refed(iter) * 2 + 1
                m2 = nmp(iter - 1) + refed(iter) * 2 + 2

                mp_new(m1, :) = (sjx(1, :) + newsjx(1, :) + sjx(2, :)) / 3.
                mp_new(m2, :) = (sjx(1, :) + newsjx(1, :) + sjx(3, :)) / 3.

                w4 = nwp(iter - 1) + refed(iter) + 1
                wp_new(w4, :) = newsjx(1, :)

                ngrmm_new(:, i) = 1
                ngrmw_new(:, i) = 1

                ngrmw_new(1, m1) = w1
                ngrmw_new(2, m1) = w2
                ngrmw_new(3, m1) = w4
                ngrmw_new(1, m2) = w1
                ngrmw_new(2, m2) = w3
                ngrmw_new(3, m2) = w4

                if(INT(sum(icl)) > 0.)then
                    do k = 1, 2, 1
                        call CheckLon(mp_new(nmp(iter - 1) + refed(iter) * 2 + k, 1))
                    end do
                    call CheckLon(wp_new(nwp(iter - 1) + refed(iter) + 1, 1))
                end if

                mrl_new(i) = 0
                ngr_mrl_new(:, i) = 0

                mrl_new(nmp(iter - 1) + refed(iter) * 2 + 1) = 2
                mrl_new(nmp(iter - 1) + refed(iter) * 2 + 2) = 2

                refed(iter) = refed(iter) + 1

            end if

        end do

        print*, "itered_num =", refed(iter)
        print*, "Refine step 4 complete"
        iter = iter + 1

        !--------------------------------------------------
        ! 4.3 Store the new grid data refined in step 4
        !--------------------------------------------------
        lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/gridfile_NXP" // trim(nxpc) // "_01_4.nc4"
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", nmp(iter - 1), spDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", nwp(iter - 1), lpDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 3, thDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_b", 7, seDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONM", NF90_FLOAT, (/ spDimID /), ncVarID(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATM", NF90_FLOAT, (/ spDimID /), ncVarID(2)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONW", NF90_FLOAT, (/ lpDimID /), ncVarID(3)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATW", NF90_FLOAT, (/ lpDimID /), ncVarID(4)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_m%iw", NF90_INT, (/ thDimID, spDimID /), ncVarID(5)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_w%im", NF90_INT, (/ seDimID, lpDimID /), ncVarID(6)))
        CALL CHECK(NF90_ENDDEF(ncID))

        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1), mp_new(1:nmp(iter - 1), 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(2), mp_new(1:nmp(iter - 1), 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(3), wp_new(1:nwp(iter - 1), 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(4), wp_new(1:nwp(iter - 1), 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5), ngrmw_new(1:3, 1:nmp(iter - 1))))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6), ngrwm_new(1:7, 1:nwp(iter - 1))))
        CALL CHECK(NF90_CLOSE(ncID))

        !stop
        !--------------------------------------------------
        ! 5.1 Calculate and store non-duplicate w points
        !--------------------------------------------------
        print*, "Start making polygonal mesh"
        num_dbx = 0
        allocate(wp_f(nwp(iter - 1), 5))
        wp_f(:, 1:2) = 9999.
        wp_f(:, 3:5) = 0.

        do i = 1, nwp(iter - 1), 1
            isexist = .false.
            do j = 1, nwp(iter - 1), 1
                if((wp_f(j, 1) == wp_new(i, 1)).and.(wp_f(j, 2) == wp_new(i, 2)))then
                    isexist = .true.
                end if
            end do
            if(isexist == .false.)then
                num_dbx = num_dbx + 1
                wp_f(num_dbx, 1) = wp_new(i, 1)
                wp_f(num_dbx, 2) = wp_new(i, 2)
            end if
        end do

        !--------------------------------------------------
        ! 5.2 Recalculate and store ngrmw
        !--------------------------------------------------
        print*, "Recalculate and store ngrmw"
        do i = 1, nmp(iter - 1), 1
            do j = 1, 3, 1
                if(ngrmw_new(j, i) /= 1)then
                    do k = 1, num_dbx, 1
                        if((wp_f(k, 1) == wp_new(ngrmw_new(j, i), 1)).and.(wp_f(k, 2) == wp_new(ngrmw_new(j, i), 2)))then
                            ngrmw_new(j, i) = k
                        end if
                    end do
                end if
            end do
        end do
        print*, "Calculation complete"

        !--------------------------------------------------
        ! 5.3 Recalculate and store ngrwm
        !--------------------------------------------------
        print*, "Recalculate and store ngrwm"
        allocate(ngrwm_f(8, nwp(iter - 1)))
        ngrwm_f(8, :) = 0
        ngrwm_f(1:7, :) = 1

        do i = 1, nmp(iter - 1), 1
            do j = 1, 3, 1
                k = ngrmw_new(j, i)
                if((k /= 0).and.(k /= 1))then

                    ngrwm_f(8, k) = ngrwm_f(8, k) + 1
                    if(ngrwm_f(8, k) > 7)then
                        !ngrwm_f(k,8) = 7
                        cycle
                    end if

                    ngrwm_f(ngrwm_f(8, k), k) = i

                end if
            end do
        end do
        print*, "The calculation is complete and the sorting begins"
        print*, minval(ngrwm_f(8, :)), maxval(ngrwm_f(8, :))

        !do i = 1,num_dbx,1
        !   if((ngrwm_f(i,9) > 8).or.(ngrwm_f(i,9) < 5))then
        !      print*,i,ngrwm_f(i,9)
        !   end if
        !end do
        !stop

        do i = 1, nmp(iter - 1), 1
            if(mp_new(i, 1) == -180.)then
                mp_new(i, 1) = 180.
            end if
        end do

        do i = 1, num_dbx, 1
            if(wp_f(i, 1) == -180.)then
                wp_f(i, 1) = 180.
            end if
        end do

        !do i = 1,num_dbx,1
        !   do j = 1,7,1
        !      if(ngrwm(j,i) == 1)then
        !         ngrwm(j,i) = 0
        !      end if
        !   end do
        !end do

        ! Sort the ngrwm (form polygons in order)
        CALL GetSort(ngrwm_f, nwp(iter - 1), mp_new, sjx_points, num_dbx, nmp(iter - 1))
        !CALL GetSort(ngrwm_f(:, 1:num_dbx), mp_new(1:nmp(iter - 1), 1:2), num_dbx, nmp(iter - 1))

        print*, "Sort completion"

        allocate(ref_pl(num_dbx, 16))
        ref_pl = 0

        do i = 1, num_dbx, 1
            do j = 1, ngrwm_f(8, i), 1
                do k = 1, 16, 1
                    if(ref_tr(ngrwm_f(j, i), k) == 1)then
                        ref_pl(i, k) = 1
                    end if
                end do
            end do
        end do

        !--------------------------------------------------
        ! 5.4 Remove refined triangles (m points)
        !--------------------------------------------------
        print*, "Start removing the refined triangle"
        num_ref = 0

        !do i = 1,sjx_points,1
        !   do j = 1,3,1
        !      if(ngrmw(j,i) == 1)then
        !         ngrmw(j,i) = 0
        !      end if
        !   end do
        !end do

        do i = 2, sjx_points, 1
            if(mrl_new(i) /= 1)then
                ngrmw_new(:, i) = 1
                num_ref = num_ref + 1
                cycle
            end if

            if((ngrmw_new(1, i) == 1).or.(ngrmw_new(2, i) == 1).or.&
                    (ngrmw_new(3, i) == 1))then
                ngrmw_new(:, i) = 1
                num_ref = num_ref + 1
            end if
        end do

        print*, "The number of triangles that have been refined is", num_ref
        allocate(mp_ref(num_ref))

        num_ref = 0

        do i = 2, sjx_points, 1
            if(mrl_new(i) /= 1)then
                num_ref = num_ref + 1
                mp_ref(num_ref) = i
                cycle
            end if

            if(ngrmw_new(1, i) == 1)then
                num_ref = num_ref + 1
                mp_ref(num_ref) = i
            end if
        end do

        num_sjx = nmp(iter - 1) - num_ref

        allocate(mp_f(num_sjx, 2))
        allocate(mrl_f(num_sjx))
        allocate(ngrmw_f(3, num_sjx))

        mp_f = 0.
        mrl_f = 0
        ngrmw_f = 1

        mrl_new(1) = 1

        i = 0

        do j = 1, nmp(iter - 1), 1
            if((j <= sjx_points).and.(mrl_new(j) /= 1))then
                cycle
            end if

            i = i + 1
            mp_f(i, 1:2) = mp_new(j, 1:2)
            mrl_f(i) = mrl_new(j)
            ngrmw_f(1:3, i) = ngrmw_new(1:3, j)
        end do

        do i = 1, num_dbx, 1
            do j = 1, 7, 1
                do k = 1, num_ref, 1
                    if(ngrwm_f(j, i) < mp_ref(k))then
                        ngrwm_f(j, i) = ngrwm_f(j, i) - (k - 1)
                        exit
                    end if
                    if(k == num_ref)then
                        ngrwm_f(j, i) = ngrwm_f(j, i) - k
                    end if
                end do
            end do
        end do

        !print*,mp_ref
        print*, minval(ngrwm_f(1:7, 1:num_dbx)), minval(ngrmw_f)
        print*, "The number of triangles before removing triangles that have been refined is", nmp(iter - 1)
        print*, "The number of triangles after removing the refined triangles is", num_sjx

        !--------------------------------------------------
        ! 5.5 Storage grid data
        !--------------------------------------------------
        lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/gridfile_NXP" // trim(nxpc) // "_01_5.nc4"
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", num_sjx, spDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", num_dbx, lpDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 3, thDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_b", 7, seDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_c", 16, sxDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONM", NF90_FLOAT, (/ spDimID /), ncVarID(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATM", NF90_FLOAT, (/ spDimID /), ncVarID(2)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONW", NF90_FLOAT, (/ lpDimID /), ncVarID(3)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATW", NF90_FLOAT, (/ lpDimID /), ncVarID(4)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_m%iw", NF90_INT, (/ thDimID, spDimID /), ncVarID(5)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_w%im", NF90_INT, (/ seDimID, lpDimID /), ncVarID(6)))
        CALL CHECK(NF90_DEF_VAR(ncID, "ref_pl", NF90_INT, (/ lpDimID, sxDimID /), ncVarID(7)))
        CALL CHECK(NF90_ENDDEF(ncID))

        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1), mp_f(1:num_sjx, 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(2), mp_f(1:num_sjx, 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(3), wp_f(1:num_dbx, 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(4), wp_f(1:num_dbx, 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5), ngrmw_f(1:3, 1:num_sjx)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6), ngrwm_f(1:7, 1:num_dbx)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(7), ref_pl))
        CALL CHECK(NF90_CLOSE(ncID))

        !stop
        !--------------------------------------------------
        ! 6.1 Calculate the side length and Angle of the triangular mesh
        !--------------------------------------------------
        sa_iter = 0
        allocate(length(0:max_sa_iter, num_sjx, 3))
        allocate(angle(0:max_sa_iter, num_sjx, 3))
        length(:, :, :) = 0.
        angle(:, :, :) = 0.

        do i = 2, num_sjx, 1
            tmpa = wp_f(ngrmw_f(1, i), 1:2)
            tmpb = wp_f(ngrmw_f(2, i), 1:2)
            tmpc = wp_f(ngrmw_f(3, i), 1:2)
            tmp_angle = angle(sa_iter, i, :)
            tmp_length = length(sa_iter, i, :)
            CALL GetTriangleLength(tmp_length, tmpa, tmpb, tmpc)
            length(sa_iter, i, :) = tmp_length
            !CALL GetTriangleLength(length(sa_iter, i, :), wp_f(ngrmw_f(1, i), 1:2), &
            !        wp_f(ngrmw_f(2, i), 1:2), wp_f(ngrmw_f(3, i), 1:2))
            CALL GetTriangleAngle(tmp_length, tmp_angle)
            angle(sa_iter, i, :) = tmp_angle
            !CALL GetTriangleAngle(length(sa_iter, i, :), angle(sa_iter, i, :))
        end do

        print*, "Triangle side length:", minval(length(sa_iter, 2:num_sjx, :)), maxval(length(sa_iter, 2:num_sjx, :))
        print*, "Triangle Angle:", minval(angle(sa_iter, 2:num_sjx, :)), maxval(angle(sa_iter, 2:num_sjx, :))
        !stop

	! --------------------------------------------------
	! 6.2 Calculate the quality of triangular mesh before elastic adjustment
	! --------------------------------------------------
	! Extr_*** refers to the minimum Angle and maximum Angle of the grid
	! Eavg_*** refers to the average of the minimum Angle and maximum Angle of the grid
	! Savg_*** refers to the standard deviation between the Angle of the mesh and the Angle of the regular polygon
	! less30 refers to the number of angles less than 30 degrees in a triangular grid        

	sa_iter = 0

        allocate(Eavg_sjx(0:max_sa_iter, 2))
        allocate(Extr_sjx(0:max_sa_iter, 2))
        allocate(Savg_sjx(0:max_sa_iter))
        allocate(less30(0:max_sa_iter))
        Extr_sjx = 0.
        Eavg_sjx = 0.
        Savg_sjx = 0.
        less30 = 0.
        ! real Gmin,Savg,Aavg,less30

        Extr_sjx(sa_iter, 1) = minval(angle(sa_iter, 2:num_sjx, :))
        Extr_sjx(sa_iter, 2) = maxval(angle(sa_iter, 2:num_sjx, :))

        do i = 2, num_sjx, 1

            Eavg_sjx(sa_iter, 1) = Eavg_sjx(sa_iter, 1) + minval(angle(sa_iter, i, :))
            Eavg_sjx(sa_iter, 2) = Eavg_sjx(sa_iter, 2) + maxval(angle(sa_iter, i, :))

            do j = 1, 3, 1
                Savg_sjx(sa_iter) = Savg_sjx(sa_iter) + (angle(sa_iter, i, j) - 60.)**2
            end do

            if(minval(angle(sa_iter, i, :)) < 30.)then
                less30(sa_iter) = less30(sa_iter) + 1.
            end if

        end do

        Eavg_sjx(sa_iter, :) = Eavg_sjx(sa_iter, :) / (num_sjx - 1)
        Savg_sjx(sa_iter) = sqrt(Savg_sjx(sa_iter) / ((num_sjx - 1) * 3))
        less30(sa_iter) = less30(sa_iter) / (num_sjx - 1)

        !--------------------------------------------------
        ! 6.3 Adjust all m points to the center of the triangle grid        
	!--------------------------------------------------
        do i = 2, num_sjx, 1
            w1 = ngrmw_f(1, i)
            w2 = ngrmw_f(2, i)
            w3 = ngrmw_f(3, i)
            icl(1) = 0
            if((w1/=1).and.(w2/=1).and.(w3/=1))then
                !      CALL MedianToCircum(wp_new(ngrmw_new(i,1),1:2),wp_new(ngrmw_new(i,2),1:2),&
                !            wp_new(ngrmw_new(i,3),1:2),mp_new(i,1:2))
                if((abs(wp_f(w1, 1) - wp_f(w2, 1)) > 180.).or.(abs(wp_f(w1, 1) - wp_f(w3, 1)) > 180.))then
                    icl(1) = 1
                    Call MoveLon(wp_f(w1, 1))
                    Call MoveLon(wp_f(w2, 1))
                    Call MoveLon(wp_f(w3, 1))
                end if

                mp_f(i, 1:2) = (wp_f(w1, 1:2) + wp_f(w2, 1:2) + wp_f(w3, 1:2)) / 3.

                if(icl(1) == 1)then
                    Call CheckLon(wp_f(w1, 1))
                    Call CheckLon(wp_f(w2, 1))
                    Call CheckLon(wp_f(w3, 1))
                    Call CheckLon(mp_f(i, 1))
                end if

            end if
        end do

        !--------------------------------------------------
        ! 6.4 Calculate the quality of polygon mesh before elastic adjustment
        !--------------------------------------------------
        n_wbx = 0
        n_lbx = 0
        n_qbx = 0

        do i = 1, num_dbx, 1
            edges = 7

            do j = 1, 7, 1
                if(ngrwm_f(j, i) == 1)then
                    edges = edges - 1
                end if
            end do

            if(edges == 5)then
                n_wbx = n_wbx + 1
            else if(edges == 6)then
                n_lbx = n_lbx + 1
            else if(edges == 7)then
                n_qbx = n_qbx + 1
            end if

            ngrwm_f(8, i) = edges

        end do

	! Extr_*** refers to the minimum Angle and maximum Angle of the grid
	! Eavg_*** refers to the average of the minimum Angle and maximum Angle of the grid
	! Savg_*** refers to the standard deviation between the Angle of the mesh and the Angle of the regular polygon        
	allocate(angle_wbx(0:max_sa_iter, n_wbx, 7))
        allocate(angle_lbx(0:max_sa_iter, n_lbx, 7))
        allocate(angle_qbx(0:max_sa_iter, n_qbx, 7))
        allocate(Eavg_wbx(0:max_sa_iter, 2))
        allocate(Eavg_lbx(0:max_sa_iter, 2))
        allocate(Eavg_qbx(0:max_sa_iter, 2))
        allocate(Savg_wbx(0:max_sa_iter))
        allocate(Savg_lbx(0:max_sa_iter))
        allocate(Savg_qbx(0:max_sa_iter))
        allocate(Extr_wbx(0:max_sa_iter, 2))
        allocate(Extr_lbx(0:max_sa_iter, 2))
        allocate(Extr_qbx(0:max_sa_iter, 2))
        angle_wbx = 0.
        angle_lbx = 0.
        angle_qbx = 0.
        Eavg_wbx = 0.
        Eavg_lbx = 0.
        Eavg_qbx = 0.
        Savg_wbx = 0.
        Savg_lbx = 0.
        Savg_qbx = 0.
        Extr_wbx = 0.
        Extr_lbx = 0.
        Extr_qbx = 0.

        w1 = 0
        w2 = 0
        w3 = 0

        do i = 2, num_dbx, 1
            if(ngrwm_f(8, i) < 5)then
                cycle
            end if

            do j = 1, ngrwm_f(8, i), 1
                dbx(j, :) = mp_f(ngrwm_f(j, i), :)
            end do

            if(ngrwm_f(8, i) == 5)then
                w1 = w1 + 1
                tmp_angle7 = angle_wbx(sa_iter, w1, :)
                CALL GetPolygonAngle(tmp_angle7, dbx, 5)
                angle_wbx(sa_iter, w1, :) = tmp_angle7
                !CALL GetPolygonAngle(angle_wbx(sa_iter, w1, :), dbx, 5)
            else if(ngrwm_f(8, i) == 6)then
                w2 = w2 + 1
                tmp_angle7 = angle_lbx(sa_iter, w2, :)
                CALL GetPolygonAngle(tmp_angle7, dbx, 6)
                angle_lbx(sa_iter, w2, :) = tmp_angle7
                !CALL GetPolygonAngle(angle_lbx(sa_iter, w2, :), dbx, 6)
            else if(ngrwm_f(8, i) == 7)then
                w3 = w3 + 1
                tmp_angle7 = angle_qbx(sa_iter, w3, :)
                CALL GetPolygonAngle(tmp_angle7, dbx, 7)
                angle_qbx(sa_iter, w3, :) = tmp_angle7
                !CALL GetPolygonAngle(angle_qbx(sa_iter, w3, :), dbx, 7)
            end if
        end do

        Extr_wbx(sa_iter, 1) = minval(angle_wbx(sa_iter, :, 1:5))
        Extr_wbx(sa_iter, 2) = maxval(angle_wbx(sa_iter, :, 1:5))
        Extr_lbx(sa_iter, 1) = minval(angle_lbx(sa_iter, :, 1:6))
        Extr_lbx(sa_iter, 2) = maxval(angle_lbx(sa_iter, :, 1:6))
        Extr_qbx(sa_iter, 1) = minval(angle_qbx(sa_iter, :, 1:7))
        Extr_qbx(sa_iter, 2) = maxval(angle_qbx(sa_iter, :, 1:7))

        do i = 1, n_wbx, 1
            Eavg_wbx(sa_iter, 1) = Eavg_wbx(sa_iter, 1) + minval(angle_wbx(sa_iter, i, 1:5))
            Eavg_wbx(sa_iter, 2) = Eavg_wbx(sa_iter, 2) + maxval(angle_wbx(sa_iter, i, 1:5))

            do j = 1, 5, 1
                Savg_wbx(sa_iter) = Savg_wbx(sa_iter) + (angle_wbx(sa_iter, i, j) - 108.)**2
            end do
        end do

        do i = 1, n_lbx, 1
            Eavg_lbx(sa_iter, 1) = Eavg_lbx(sa_iter, 1) + minval(angle_lbx(sa_iter, i, 1:6))
            Eavg_lbx(sa_iter, 2) = Eavg_lbx(sa_iter, 2) + maxval(angle_lbx(sa_iter, i, 1:6))

            do j = 1, 6, 1
                Savg_lbx(sa_iter) = Savg_lbx(sa_iter) + (angle_lbx(sa_iter, i, j) - 120.)**2
            end do
        end do

        do i = 1, n_qbx, 1
            Eavg_qbx(sa_iter, 1) = Eavg_qbx(sa_iter, 1) + minval(angle_qbx(sa_iter, i, 1:7))
            Eavg_qbx(sa_iter, 2) = Eavg_qbx(sa_iter, 2) + maxval(angle_qbx(sa_iter, i, 1:7))

            do j = 1, 7, 1
                Savg_qbx(sa_iter) = Savg_qbx(sa_iter) + (angle_qbx(sa_iter, i, j) - 129.)**2
            end do
        end do

        Eavg_wbx(sa_iter, :) = Eavg_wbx(sa_iter, :) / n_wbx
        Eavg_lbx(sa_iter, :) = Eavg_lbx(sa_iter, :) / n_lbx
        Eavg_qbx(sa_iter, :) = Eavg_qbx(sa_iter, :) / n_qbx
        Savg_wbx(sa_iter) = sqrt(Savg_wbx(sa_iter) / (n_wbx * 5.))
        Savg_lbx(sa_iter) = sqrt(Savg_lbx(sa_iter) / (n_lbx * 6.))
        Savg_qbx(sa_iter) = sqrt(Savg_qbx(sa_iter) / (n_qbx * 7.))

        print*, "The pentagon grid information is as follows:"
         print*, "Number of grids ", n_wbx
         print *, "minimum angle", minval (angle_wbx (sa_iter, :, 1)),"maximum angle ", maxval (angle_wbx (sa_iter, :, 1))
         print*, "Hexagonal grid information is as follows:"
         print*, "Number of grids ", n_lbx
         print *, " minimum angle", minval (angle_lbx (sa_iter, :, 1:6)),"maximum angle ", maxval (angle_lbx (sa_iter, :, 1:6))
         print*, "Heptagon grid information is as follows:"
         print*, "Number of grids ", n_qbx
         Print *, "minimum angle", minval (angle_qbx (sa_iter, :, 1:7)),"maximum angle ", maxval (angle_qbx (sa_iter, :, 1:7))
        !stop

        !--------------------------------------------------
        ! 6.3 Elastic Adjustment & Calculate the quality of the triangular mesh after each adjustment
        !--------------------------------------------------
        print*, "Initial elastic adjustment"
        allocate(MoveDis(num_dbx, 2))

        sa_iter = 1

        End_SpringAjustment = .false.

        do while(End_SpringAjustment == .false.)

            MoveDis = 0.
            End_SpringAjustment = .true.
            k = 0

            do i = sjx_points - num_ref + 1, num_sjx, 1

                do j = 1, 3, 1
                    if(angle(sa_iter - 1, i, j) > 90.)then
                        k = k + 1
                        fra = (1 - (72. / angle(sa_iter - 1, i, j))) / 100.
                        if(j == 1)then
                            w1 = ngrmw_f(2, i)
                            w2 = ngrmw_f(3, i)
                        else if(j == 2)then
                            w1 = ngrmw_f(1, i)
                            w2 = ngrmw_f(3, i)
                        else if(j == 3)then
                            w1 = ngrmw_f(1, i)
                            w2 = ngrmw_f(2, i)
                        end if
                        rx = wp_f(w2, 1) - wp_f(w1, 1)
                        ry = wp_f(w2, 2) - wp_f(w1, 2)

                        call CheckLon(rx)

                        MoveDis(w1, 1) = MoveDis(w1, 1) + rx * fra / 2.
                        MoveDis(w1, 2) = MoveDis(w1, 2) + ry * fra / 2.
                        MoveDis(w2, 1) = MoveDis(w2, 1) - rx * fra / 2.
                        MoveDis(w2, 2) = MoveDis(w2, 2) - ry * fra / 2.

                    end if
                end do
            end do

            do i = 1, num_dbx, 1
                wp_f(i, 1:2) = wp_f(i, 1:2) + MoveDis(i, 1:2)
                call CheckLon(wp_f(i, 1))
            end do

            do i = 2, num_sjx, 1
                tmpa = wp_f(ngrmw_f(1, i), 1:2)
                tmpb = wp_f(ngrmw_f(2, i), 1:2)
                tmpc = wp_f(ngrmw_f(3, i), 1:2)
                tmp_angle = angle(sa_iter, i, :)
                tmp_length = length(sa_iter, i, :)
                CALL GetTriangleLength(tmp_length, tmpa, tmpb, tmpc)      
                length(sa_iter, i, :) = tmp_length
                !CALL GetTriangleLength(length(sa_iter, i, :), wp_f(ngrmw_f(1, i), 1:2), &
                !        wp_f(ngrmw_f(2, i), 1:2), wp_f(ngrmw_f(3, i), 1:2))
                CALL GetTriangleAngle(tmp_length, tmp_angle)
                angle(sa_iter, i, :) = tmp_angle
                !CALL GetTriangleAngle(length(sa_iter, i, :), angle(sa_iter, i, :))
            end do

            Extr_sjx(sa_iter, 1) = minval(angle(sa_iter, 2:num_sjx, :))
            Extr_sjx(sa_iter, 2) = maxval(angle(sa_iter, 2:num_sjx, :))

            do i = 2, num_sjx, 1

                Eavg_sjx(sa_iter, 1) = Eavg_sjx(sa_iter, 1) + minval(angle(sa_iter, i, :))
                Eavg_sjx(sa_iter, 2) = Eavg_sjx(sa_iter, 2) + maxval(angle(sa_iter, i, :))

                do j = 1, 3, 1
                    Savg_sjx(sa_iter) = Savg_sjx(sa_iter) + (angle(sa_iter, i, j) - 60.)**2
                end do

                if(minval(angle(sa_iter, i, :)) < 30.)then
                    less30(sa_iter) = less30(sa_iter) + 1.
                end if

            end do

            Eavg_sjx(sa_iter, :) = Eavg_sjx(sa_iter, :) / (num_sjx - 1)
            Savg_sjx(sa_iter) = sqrt(Savg_sjx(sa_iter) / ((num_sjx - 1) * 3))
            less30(sa_iter) = less30(sa_iter) / (num_sjx - 1)

            if((k /= 0).and.(sa_iter < max_sa_iter))then
                End_SpringAjustment = .false.
            end if

            do i = 2, num_sjx, 1
                w1 = ngrmw_f(1, i)
                w2 = ngrmw_f(2, i)
                w3 = ngrmw_f(3, i)
                icl(1) = 0
                if((w1/=1).and.(w2/=1).and.(w3/=1))then
                    if((abs(wp_f(w1, 1) - wp_f(w2, 1)) > 180.).or.(abs(wp_f(w1, 1) - wp_f(w3, 1)) > 180.))then
                        icl(1) = 1
                        Call MoveLon(wp_f(w1, 1))
                        Call MoveLon(wp_f(w2, 1))
                        Call MoveLon(wp_f(w3, 1))
                    end if

                    mp_f(i, 1:2) = (wp_f(w1, 1:2) + wp_f(w2, 1:2) + wp_f(w3, 1:2)) / 3.

                    if(icl(1) == 1)then
                        Call CheckLon(wp_f(w1, 1))
                        Call CheckLon(wp_f(w2, 1))
                        Call CheckLon(wp_f(w3, 1))
                        Call CheckLon(mp_f(i, 1))
                    end if

                end if
            end do

            w1 = 0
            w2 = 0
            w3 = 0

            do i = 2, num_dbx, 1
                if(ngrwm_f(8, i) < 5)then
                    cycle
                end if

                do j = 1, ngrwm_f(8, i), 1
                    dbx(j, :) = mp_f(ngrwm_f(j, i), :)
                end do

                if(ngrwm_f(8, i) == 5)then
                    w1 = w1 + 1
                    tmp_angle7 = angle_wbx(sa_iter, w1, :)
                    CALL GetPolygonAngle(tmp_angle7, dbx, 5)
                    angle_wbx(sa_iter, w1, :) = tmp_angle7
                    !CALL GetPolygonAngle(angle_wbx(sa_iter, w1, :), dbx, 5)
                else if(ngrwm_f(8, i) == 6)then
                    w2 = w2 + 1
                    tmp_angle7 = angle_lbx(sa_iter, w2, :)
                    CALL GetPolygonAngle(tmp_angle7, dbx, 6)
                    angle_lbx(sa_iter, w2, :) = tmp_angle7
                    !CALL GetPolygonAngle(angle_lbx(sa_iter, w2, :), dbx, 6)
                else if(ngrwm_f(8, i) == 7)then
                    w3 = w3 + 1
                    tmp_angle7 = angle_qbx(sa_iter, w3, :)
                    CALL GetPolygonAngle(tmp_angle7, dbx, 7)
                    angle_qbx(sa_iter, w3, :) = tmp_angle7
                    !CALL GetPolygonAngle(angle_qbx(sa_iter, w3, :), dbx, 7)
                end if
            end do

            Extr_wbx(sa_iter, 1) = minval(angle_wbx(sa_iter, :, 1:5))
            Extr_wbx(sa_iter, 2) = maxval(angle_wbx(sa_iter, :, 1:5))
            Extr_lbx(sa_iter, 1) = minval(angle_lbx(sa_iter, :, 1:6))
            Extr_lbx(sa_iter, 2) = maxval(angle_lbx(sa_iter, :, 1:6))
            Extr_qbx(sa_iter, 1) = minval(angle_qbx(sa_iter, :, 1:7))
            Extr_qbx(sa_iter, 2) = maxval(angle_qbx(sa_iter, :, 1:7))

            do i = 1, n_wbx, 1
                Eavg_wbx(sa_iter, 1) = Eavg_wbx(sa_iter, 1) + minval(angle_wbx(sa_iter, i, 1:5))
                Eavg_wbx(sa_iter, 2) = Eavg_wbx(sa_iter, 2) + maxval(angle_wbx(sa_iter, i, 1:5))

                do j = 1, 5, 1
                    Savg_wbx(sa_iter) = Savg_wbx(sa_iter) + (angle_wbx(sa_iter, i, j) - 108.)**2
                end do
            end do

            do i = 1, n_lbx, 1
                Eavg_lbx(sa_iter, 1) = Eavg_lbx(sa_iter, 1) + minval(angle_lbx(sa_iter, i, 1:6))
                Eavg_lbx(sa_iter, 2) = Eavg_lbx(sa_iter, 2) + maxval(angle_lbx(sa_iter, i, 1:6))

                do j = 1, 6, 1
                    Savg_lbx(sa_iter) = Savg_lbx(sa_iter) + (angle_lbx(sa_iter, i, j) - 120.)**2
                end do
            end do

            do i = 1, n_qbx, 1
                Eavg_qbx(sa_iter, 1) = Eavg_qbx(sa_iter, 1) + minval(angle_qbx(sa_iter, i, 1:7))
                Eavg_qbx(sa_iter, 2) = Eavg_qbx(sa_iter, 2) + maxval(angle_qbx(sa_iter, i, 1:7))

                do j = 1, 7, 1
                    Savg_qbx(sa_iter) = Savg_qbx(sa_iter) + (angle_qbx(sa_iter, i, j) - 129.)**2
                end do
            end do

            Eavg_wbx(sa_iter, :) = Eavg_wbx(sa_iter, :) / n_wbx
            Eavg_lbx(sa_iter, :) = Eavg_lbx(sa_iter, :) / n_lbx
            Eavg_qbx(sa_iter, :) = Eavg_qbx(sa_iter, :) / n_qbx
            Savg_wbx(sa_iter) = sqrt(Savg_wbx(sa_iter) / (n_wbx * 5))
            Savg_lbx(sa_iter) = sqrt(Savg_lbx(sa_iter) / (n_lbx * 6))
            Savg_qbx(sa_iter) = sqrt(Savg_qbx(sa_iter) / (n_qbx * 7))

            sa_iter = sa_iter + 1

        end do

        !--------------------------------------------------
        ! 6.4 Store grid quality data
        !--------------------------------------------------
        if(max_iter == 1)then
            lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/result/quality_NXP" // trim(nxpc) // "_lbx.nc4"
        else
            lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/quality_NXP" // trim(nxpc) // "_1.nc4"
        end if

        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "num_iter", max_sa_iter + 1, lpDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 2, twDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "less30", NF90_FLOAT, (/ lpDimID /), ncVarID(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Eavg_sjx", NF90_FLOAT, (/ lpDimID, twDimID /), ncVarID(2)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Savg_sjx", NF90_FLOAT, (/ lpDimID /), ncVarID(3)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Extr_sjx", NF90_FLOAT, (/ lpDimID, twDimID /), ncVarID(4)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Eavg_wbx", NF90_FLOAT, (/ lpDimID, twDimID /), ncVarID(5)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Eavg_lbx", NF90_FLOAT, (/ lpDimID, twDimID /), ncVarID(6)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Eavg_qbx", NF90_FLOAT, (/ lpDimID, twDimID /), ncVarID(7)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Savg_wbx", NF90_FLOAT, (/ lpDimID /), ncVarID(8)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Savg_lbx", NF90_FLOAT, (/ lpDimID /), ncVarID(9)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Savg_qbx", NF90_FLOAT, (/ lpDimID /), ncVarID(10)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Extr_wbx", NF90_FLOAT, (/ lpDimID, twDimID /), ncVarID(11)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Extr_lbx", NF90_FLOAT, (/ lpDimID, twDimID /), ncVarID(12)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Extr_qbx", NF90_FLOAT, (/ lpDimID, twDimID /), ncVarID(13)))

        CALL CHECK(NF90_ENDDEF(ncID))

        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1), less30(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(2), Eavg_sjx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(3), Savg_sjx(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(4), Extr_sjx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5), Eavg_wbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6), Eavg_lbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(7), Eavg_qbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(8), Savg_wbx(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(9), Savg_lbx(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(10), Savg_qbx(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(11), Extr_wbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(12), Extr_lbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(13), Extr_qbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_CLOSE(ncID))
        !stop

        !--------------------------------------------------
        ! 6.5 Store the final grid data
        !--------------------------------------------------

        allocate(mp_f_tmp(1:num_sjx,1:2)); mp_f_tmp = mp_f(1:num_sjx,1:2)
        allocate(wp_f_tmp(1:num_dbx,1:2)); wp_f_tmp = wp_f(1:num_dbx,1:2)
        allocate(ngrmw_f_tmp(1:3, 1:num_sjx)); ngrmw_f_tmp = ngrmw_f(1:3, 1:num_sjx)
        allocate(ngrwm_f_tmp(1:7, 1:num_dbx)); ngrwm_f_tmp = ngrwm_f(1:7, 1:num_dbx)
        allocate(dismm(num_dbx,7)); dismm = 0.
        allocate(disww(num_sjx,3)); disww = 0.
        Call Get_Dis(mp_f_tmp,wp_f_tmp,ngrmw_f_tmp,ngrwm_f_tmp,dismm,disww,num_sjx,num_dbx)

        if(max_iter == 1)then
            lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/result/gridfile_NXP" // trim(nxpc) // "_lbx.nc4"
        else
            lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/gridfile_NXP" // trim(nxpc) // "_01_6.nc4"
        end if

        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", num_sjx, spDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", num_dbx, lpDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 3, thDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_b", 7, seDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_c", 16, sxDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONM", NF90_FLOAT, (/ spDimID /), ncVarID(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATM", NF90_FLOAT, (/ spDimID /), ncVarID(2)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLONW", NF90_FLOAT, (/ lpDimID /), ncVarID(3)))
        CALL CHECK(NF90_DEF_VAR(ncID, "GLATW", NF90_FLOAT, (/ lpDimID /), ncVarID(4)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_m%iw", NF90_INT, (/ thDimID, spDimID /), ncVarID(5)))
        CALL CHECK(NF90_DEF_VAR(ncID, "itab_w%im", NF90_INT, (/ seDimID, lpDimID /), ncVarID(6)))
        CALL CHECK(NF90_DEF_VAR(ncID, "ref_pl", NF90_INT, (/ lpDimID, sxDimID /), ncVarID(7)))
        CALL CHECK(NF90_DEF_VAR(ncID, "dis_w%iw", NF90_FLOAT, (/ spDimID, thDimID /), ncVarID(8)))
        CALL CHECK(NF90_DEF_VAR(ncID, "dis_m%im", NF90_FLOAT, (/ lpDimID, seDimID /), ncVarID(9)))
        CALL CHECK(NF90_ENDDEF(ncID))

        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1), mp_f(1:num_sjx, 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(2), mp_f(1:num_sjx, 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(3), wp_f(1:num_dbx, 1)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(4), wp_f(1:num_dbx, 2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5), ngrmw_f(1:3, 1:num_sjx)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6), ngrwm_f(1:7, 1:num_dbx)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(7), ref_pl))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(8), disww))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(9), dismm))
        CALL CHECK(NF90_CLOSE(ncID))

        deallocate(dismm)
        deallocate(disww)
        deallocate(mp_f_tmp)
        deallocate(wp_f_tmp)
        deallocate(ngrmw_f_tmp)
        deallocate(ngrwm_f_tmp)

        !stop
    END SUBROUTINE refine_lbx

    SUBROUTINE CHECK(STATUS)
        INTEGER, intent (in) :: STATUS
        if  (STATUS /= NF90_NOERR) then 
            print *, NF90_STRERROR(STATUS)
            stop 'stopped'
        endif
    END SUBROUTINE CHECK


    INTEGER FUNCTION IsCrossLine(x1, x2)

        implicit none

        real(r8), intent(in) :: x1, x2
        IsCrossLine = 0

        if(abs(x1 - x2) > 180.)then
            IsCrossLine = 1
        end if

    END FUNCTION IsCrossLine


    INTEGER FUNCTION IsNgrmm(a, b)

        IMPLICIT NONE

        integer, intent(in) :: a(3), b(3)

        IsNgrmm = 0

        if((a(3) == 1).or.(b(3) == 1))then
            return
        end if

        if((a(1)==b(1)).or.(a(1)==b(2)).or.(a(1)==b(3)))then
            if((a(2)==b(1)).or.(a(2)==b(2)).or.(a(2)==b(3)))then
                IsNgrmm = 3
                return
            else if((a(3)==b(1)).or.(a(3)==b(2)).or.(a(3)==b(3)))then
                IsNgrmm = 2
                return
            else
                return
            end if
        else if((a(2)==b(1)).or.(a(2)==b(2)).or.(a(2)==b(3)))then
            if((a(3)==b(1)).or.(a(3)==b(2)).or.(a(3)==b(3)))then
                IsNgrmm = 1
                return
            else
                return
            end if
        end if

    END FUNCTION IsNgrmm


    SUBROUTINE CheckLon(x)

        implicit none

        real(r8) :: x

        if(x > 180.)then
            x = x - 360.
        end if

        if(x == -180.)then
            x = 180.
        end if

        if(x < -180.)then
            x = x + 360.
        end if

    END SUBROUTINE CheckLon


    SUBROUTINE MoveLon(x)

        implicit none

        real(r8) :: x

        if(x < 0.)then
            x = x + 360.
        end if

    END SUBROUTINE MoveLon


    SUBROUTINE GetSort(ngr, num, mp, sjx_points, nwp, nmp)

        implicit none

        integer, intent(in) :: nwp, nmp, num, sjx_points
        integer :: i, j, k
        logical :: icl

        integer, dimension(8, num), intent(out) :: ngr
        real(r8), dimension(sjx_points*4, 2), intent(in) :: mp
        !real(r8),dimension(nwp,2),intent(in) :: wp
        real(r8) :: pi, temp, center(2)

        real(r8), allocatable :: angle(:), points(:, :)

        pi = 3.1415926535
        temp = 0.
        allocate(points(7, 2))
        allocate(angle(7))

        do i = 1, nwp, 1
            if(ngr(8, i) > 4)then

                icl = .false.
                points = 0.
                angle = 0.
                center = 0.

                do j = 1, ngr(8, i), 1
                    points(j, 1:2) = mp(ngr(j, i), 1:2)
                end do

                do j = 1, ngr(8, i) - 1, 1
                    if(abs(points(j, 1) - points(j + 1, 1)) > 180.)then
                        icl = .true.
                    end if
                end do

                if(icl == .true.)then
                    do j = 1, ngr(8, i), 1
                        if(points(j, 1) < 0.)then
                            points(j, 1) = points(j, 1) + 360.
                        end if
                    end do
                end if

                do j = 1, ngr(8, i), 1
                    center(1:2) = center(1:2) + points(j, 1:2)
                end do
                center(1:2) = center(1:2) / ngr(8, i)

                do j = 1, ngr(8, i), 1
                    points(j, 1:2) = points(j, 1:2) - center(1:2)
                    if(points(j, 2) > 0.)then
                        if(points(j, 1) == 0.)then
                            angle(j) = pi / 2.
                        else
                            angle(j) = atan(points(j, 2) / points(j, 1))
                            if(points(j, 1) < 0.)then
                                angle(j) = angle(j) + pi
                            end if
                        end if
                    else if(points(j, 2) < 0.)then
                        if(points(j, 1) == 0.)then
                            angle(j) = 1.5 * pi
                        else if(points(j, 1) < 0.)then
                            angle(j) = atan(points(j, 2) / points(j, 1)) + pi
                        else
                            angle(j) = atan(points(j, 2) / points(j, 1)) + 2. * pi
                        end if
                    else
                        if(points(j, 1) > 0.)then
                            angle(j) = 0.
                        else if(points(j, 1) < 0.)then
                            angle(j) = pi
                        end if
                    end if
                end do

                do j = 1, ngr(8, i) - 1, 1
                    do k = j + 1, ngr(8, i), 1
                        if(angle(j) > angle(k))then
                            temp = angle(j)
                            angle(j) = angle(k)
                            angle(k) = temp
                            temp = ngr(j, i)
                            ngr(j, i) = ngr(k, i)
                            ngr(k, i) = int(temp)
                        end if
                    end do
                end do

            end if
        end do

        deallocate(angle)
        deallocate(points)

    END SUBROUTINE GetSort


    ! Calculate the outer center (intersection of vertical bisectors)p of the triangle abc
    SUBROUTINE MedianToCircum(a, b, c, p)

        implicit none

        real(r8), intent(in) :: a(2), b(2), c(2)
        real(r8), intent(out) :: p(2)
        real(r8) :: k1, k2, b1, b2, points(3, 2)
        real(r8) :: m1, m2, n1, n2

        integer :: i
        logical :: iscross

        iscross = .false.

        points(1, :) = a(:)
        points(2, :) = b(:)
        points(3, :) = c(:)

        k1 = (points(2, 2) - points(1, 2)) / (points(2, 1) - points(1, 1))
        k2 = (points(2, 2) - points(3, 2)) / (points(2, 1) - points(3, 1))
        b1 = points(2, 2) - k1 * points(2, 1)
        b2 = points(2, 2) - k2 * points(2, 1)
        m1 = -1 / k1
        m2 = -1 / k2
        n1 = (points(2, 2) + points(1, 2)) / 2 - m1 * (points(2, 1) + points(1, 1)) / 2
        n2 = (points(2, 2) + points(3, 2)) / 2 - m2 * (points(2, 1) + points(3, 1)) / 2

        if(points(1, 2)==points(2, 2))then
            p(1) = (points(2, 1) + points(1, 1)) / 2
        else if(points(3, 2)==points(2, 2))then
            p(1) = (points(2, 1) + points(3, 1)) / 2
        else if((points(1, 1)/=points(2, 1)).and.(points(3, 1)/=points(2, 1)))then
            p(1) = (n2 - n1) / (m1 - m2)
            p(2) = (m1 * n2 - m2 * n1) / (m1 - m2)
        else if(points(1, 1)==points(2, 1))then
            p(2) = (points(2, 2) + points(1, 2)) / 2
            p(1) = (p(2) - n2) / m2
        else if(points(3, 1)==points(2, 1))then
            p(2) = (points(2, 2) + points(3, 2)) / 2
            p(1) = (p(2) - n1) / m1
        end if

    END SUBROUTINE MedianToCircum


    SUBROUTINE GetTriangleLength(length, a, b, c)

        implicit none

        integer :: i
        real(r8) :: R, pi, px(3), py(3), v(3), sjx(3, 2)
        real(r8), intent(in) :: a(2), b(2), c(2)
        real(r8), intent(out) :: length(3)

        R = 6371.
        pi = 3.1415926535

        sjx(1, :) = a(:)
        sjx(2, :) = b(:)
        sjx(3, :) = c(:)

        if((abs(a(1) - b(1))>180.).or.(abs(b(1) - c(1))>180.))then
            do i = 1, 3, 1
                if(sjx(i, 1)<0.)then
                    sjx(i, 1) = sjx(i, 1) + 360.
                end if
            end do
        end if

        px(1) = sjx(1, 1) * pi / 180.
        px(2) = sjx(2, 1) * pi / 180.
        px(3) = sjx(3, 1) * pi / 180.
        py(1) = sjx(1, 2) * pi / 180.
        py(2) = sjx(2, 2) * pi / 180.
        py(3) = sjx(3, 2) * pi / 180.

        v(1) = sin(py(2)) * sin(py(3)) + cos(py(2)) * cos(py(3)) * cos(px(2) - px(3))
        v(2) = sin(py(1)) * sin(py(3)) + cos(py(1)) * cos(py(3)) * cos(px(1) - px(3))
        v(3) = sin(py(1)) * sin(py(2)) + cos(py(1)) * cos(py(2)) * cos(px(1) - px(2))

        length(1) = R * acos(v(1)) * 1000
        length(2) = R * acos(v(2)) * 1000
        length(3) = R * acos(v(3)) * 1000

    END SUBROUTINE GetTriangleLength


    SUBROUTINE GetTriangleAngle(l, angle)

        implicit none

        real(r8) :: pi
        real(r8), intent(in) :: l(3)
        real(r8), intent(out) :: angle(3)

        pi = 3.1415926535

        angle(1) = acos((l(2) * l(2) + l(3) * l(3) - l(1) * l(1)) / (2 * l(2) * l(3))) * 180. / pi
        angle(2) = acos((l(1) * l(1) + l(3) * l(3) - l(2) * l(2)) / (2 * l(1) * l(3))) * 180. / pi
        angle(3) = acos((l(2) * l(2) + l(1) * l(1) - l(3) * l(3)) / (2 * l(2) * l(1))) * 180. / pi

    END SUBROUTINE GetTriangleAngle


    SUBROUTINE GetPolygonAngle(angle, dbx, edges)

        implicit none

        integer :: i
        integer, intent(in) :: edges
        real(r8), intent(in) :: dbx(7, 2)
        real(r8), intent(out) :: angle(7)
        real(r8) :: l(3), a(2), b(2), c(2), pi

        pi = 3.1415926535
        l = 0.

        do i = 1, edges, 1
            if(i <= (edges - 2))then
                a(:) = dbx(i, :)
                b(:) = dbx(i + 1, :)
                c(:) = dbx(i + 2, :)
            else if(i == (edges - 1))then
                a(:) = dbx(i, :)
                b(:) = dbx(i + 1, :)
                c(:) = dbx(1, :)
            else if(i == edges)then
                a(:) = dbx(i, :)
                b(:) = dbx(1, :)
                c(:) = dbx(2, :)
            end if

            CALL GetTriangleLength(l, a, b, c)

            angle(i) = acos((l(1) * l(1) + l(3) * l(3) - l(2) * l(2)) / (2. * l(1) * l(3))) * 180. / pi

        end do

    END SUBROUTINE GetPolygonAngle


    SUBROUTINE IsInRefineArea(IsInRfArea,sjx_points,lbx_points,ngrmw,wp)

      implicit none

      integer :: i,j,n
      integer,intent(in) :: sjx_points,lbx_points
      integer,dimension(3,sjx_points),intent(in) :: ngrmw
      real(r8),dimension(lbx_points,2),intent(in) :: wp
      real(r8) :: sjx(3,2)

      logical,dimension(sjx_points),intent(out) :: IsInRfArea

      sjx = 0.
      IsInRfArea = .false.

      do i = 1,sjx_points,1
         do j = 1,3,1
            sjx(j,1:2) = wp(ngrmw(j,i),1:2)
         end do

         do j = 1,3,1
            do n = 1,ndm_refine,1
               if((sjx(j,1) >= edgew_rf(n)).and.(sjx(j,1) <= edgee_rf(n)).and.&
                     (sjx(j,2) >= edges_rf(n)).and.(sjx(j,2) <= edgen_rf(n)))then
                  IsInRfArea(i) = .true.
                  cycle
               end if
            end do
         end do

      end do

     END SUBROUTINE IsInRefineArea

END module MOD_refine_lbx
