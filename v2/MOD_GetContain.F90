module MOD_GetContain
    USE consts_coms, only: r8, mode, nxp, refine, file_dir, edgee, edgew, edges, edgen, ndm_domain, openmp, nlons_source, nlats_source, erad, pio180, no_caculate_fraction
    USE refine_vars, only: step, max_iter, edgee_rf, edgew_rf, edges_rf, edgen_rf, ndm_refine
    USE netcdf
    USE MOD_file_preprocess, only: Unstructured_Mesh_Read, Contain_Save, Contain_Read
    USE MOD_data_preprocess, only: lon_i, lat_i, lon_vertex, lat_vertex
    USE MOD_Area_judge, only: seaorland, IsInRfArea_grid, IsInDmArea_grid
    implicit none
    logical, dimension(:), allocatable, public :: IsInRfArea_sjx
    integer,  dimension(:, :), allocatable, public :: patchtypes
    contains
    ! 包含关系的计算与patchID的获取
    SUBROUTINE  Get_Contain()

        integer :: ustr_points, sjx_points, lbx_points            ! 三角形网格数与多边形网格数
        integer :: numpatch_new
        integer :: varid, ncid, Dim_ustrID, Dim_aID
        integer,  dimension(:), allocatable :: n_ngrmw, n_ngrwm, ustr_n_ngr
        integer,  dimension(:, :), allocatable:: ngrmw, ngrwm, ustr_ngr ! in the subroutine as ngrwm (sjx--ngrmw lbx--ngrwm)
        real(r8), dimension(:, :), allocatable:: mp, wp, ustr_vertex ! in the subroutine as wp or mp
        integer,  dimension(:, :), allocatable :: ustr_id_new, ustr_ii_new
        character(LEN = 256) :: lndname, outputfile
        character(LEN = 5) :: nxpc, stepc, modec
        logical,  dimension(:), allocatable :: IsInDmArea_ustr
        
        write(nxpc, '(I3.3)') NXP ! 相当于是integer转字符串character
        write(modec, '(I1.1)') mode
        write(stepc, '(I2.2)') step 
        print*, "开始读取非结构网格数据 in the refine area MOD_GetContain.F90"
        if (refine == .true. .and. step < max_iter) then ! 
            print*, "refine on-going"
            lndname = trim(file_dir) // 'gridfile/gridfile_NXP' // trim(nxpc) // '_'//trim(stepc)//'.nc4'
            print*, lndname
            CALL Unstructured_Mesh_Read(lndname, sjx_points, lbx_points, mp, wp, ngrmw, ngrwm, n_ngrwm)
            print*, "mesh read : lndname, sjx_points, lbx_points, mp, wp, ngrmw, ngrwm, n_ngrwm"
            print*, "非结构网格数据读取完成 in the refine area MOD_GetContain.F90"
            print*, ""
            print*, "共有", sjx_points, "个三角形网格和", lbx_points, "个多边形网格"
            print*, ""

            allocate(n_ngrmw(sjx_points)); n_ngrmw = 3; n_ngrmw(1) = -1
            ! 获取当前迭代次数下，细化区域内需要计算的三角形/多边形
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            print*, "IsInArea_sjx_Calculation start"
            allocate(IsInRfArea_sjx(sjx_points)); IsInRfArea_sjx = .false.
            CALL IsInArea_ustr_Calculation(edgee_rf, edgew_rf, edges_rf, edgen_rf, ndm_refine, sjx_points, wp, ngrmw, n_ngrmw, IsInRfArea_sjx) 
            print*, "IsInArea_sjx_Calculation finish"

            print*, "Contain_Calculation start"
            ! 但是在阈值细化里面就是都是当作三角形来看，而且采用的是细化区域的信息 
            allocate(patchtypes(nlons_source, nlats_source)); patchtypes = 0
            CALL Contain_Calculation(sjx_points, wp, ngrmw, n_ngrmw, IsInRfArea_grid, IsInRfArea_sjx, ustr_id_new, ustr_ii_new, numpatch_new)
            print*, "Contain_Calculation finish"

            print*, "三角形网格包含数组计算完成，开始进行存储 in the GetContain.F90"
            lndname = trim(file_dir) // 'contain/contain_refine_NXP' //trim(nxpc)//'_'//trim(stepc)//'_mp.nc4'
            print*, lndname
            print*, ""
            CALL Contain_Save(lndname, sjx_points, numpatch_new, ustr_id_new, ustr_ii_new)
            print*, "非结构网格存储完成 in the GetContain.F90"
            print*, ""
            ! deallocate
            deallocate(n_ngrmw); deallocate(n_ngrwm);
            deallocate(ustr_id_new); deallocate(ustr_ii_new);
            deallocate(patchtypes)

        else
            print*,"refine finish and contain calculate in the domain region"
            print*, "开始读取非结构网格数据 in the GetContain.F90"
            print*, ""
            lndname = trim(file_dir) //  "gridfile/gridfile_NXP"//trim(nxpc)//"_"//trim(stepc)//".nc4"
            print*, lndname
            CALL Unstructured_Mesh_Read(lndname, sjx_points, lbx_points, mp, wp, ngrmw, ngrwm, n_ngrwm)
            outputfile = trim(file_dir)//"result/gridfile_NXP"//trim(nxpc)// "_"//trim(modec) //".nc4"
            CALL execute_command_line('cp '//trim(lndname)//' '//trim(outputfile))
            print*, "非结构网格数据读取完成 in the GetContain.F90"
            print*, ""
            print*, "共有", sjx_points, "个三角形网格和", lbx_points, "个多边形网格"
            print*, ""
            
            if(mode == 3)then
                print*, "开始计算三角形网格与经纬度网格包含关系数组大小......"
                ustr_points = sjx_points
                allocate(ustr_vertex(lbx_points, 2)); ustr_vertex = wp ! use wp rather than mp ???
                allocate(ustr_ngr(3, ustr_points));   ustr_ngr = ngrmw
                allocate(ustr_n_ngr(ustr_points)); ustr_n_ngr = 3; ustr_n_ngr(1) = -1
            else if(mode == 6)then
                print*, "开始计算六边形网格与经纬度网格包含关系数组大小......"
                ustr_points = lbx_points
                allocate(ustr_vertex(sjx_points, 2)); ustr_vertex = mp ! use mp rather than wp ??
                allocate(ustr_ngr(7, ustr_points));   ustr_ngr = ngrwm
                allocate(ustr_n_ngr(ustr_points));    ustr_n_ngr = n_ngrwm
                deallocate(n_ngrwm)
            end if
          
            print*, "IsInArea_ustr_Calculation start"
            ! 获取包含关系计算区域内需要计算的三角形/多边形
            allocate(IsInDmArea_ustr(ustr_points)); IsInDmArea_ustr = .false.
            call IsInArea_ustr_Calculation(edgee, edgew, edges, edgen, ndm_domain, ustr_points, ustr_vertex, ustr_ngr, ustr_n_ngr, IsInDmArea_ustr)
            print*, "IsInArea_ustr_Calculation finish"

            print*, "Contain_Calculation start"
            ! 这里应该是包含关系区域的才要计算
            allocate(patchtypes(nlons_source, nlats_source)); patchtypes = 0
            CALL Contain_Calculation(ustr_points, ustr_vertex, ustr_ngr, ustr_n_ngr, IsInDmArea_grid, IsInDmArea_ustr, ustr_id_new, ustr_ii_new, numpatch_new)
            print*, "包含关系数组大小计算完成......"
            print*, "Contain_Calculation finish"

            lndname = trim(file_dir) // 'contain/contain_domain_NXP'//trim(nxpc)//'_mode'//trim(modec)//'.nc4'
            print*, lndname
            print*, "非结构网格包含数组开始进行存储 in the GetContain.F90"
            CALL Contain_Save(lndname, ustr_points, numpatch_new, ustr_id_new, ustr_ii_new)
            print*, "非结构网格包含数组存储完成 in the GetContain.F90"
            print*, ""

            print*, "PatchID_save start"
            ! 开始计算mpi模式所需的patchID文件
            outputfile = trim(file_dir) // 'patchtype/patchtype_NXP'//trim(nxpc)//'_mode'//trim(modec)//'.nc4' ! 这是最终的>    输出结果
            CALL PatchID_Save(outputfile, patchtypes)
            print*, "PatchID_save finish"

            deallocate(ustr_vertex); deallocate(ustr_ngr); deallocate(ustr_n_ngr);
            deallocate(IsInDmArea_ustr)
            deallocate(ustr_id_new, ustr_ii_new)
            deallocate(patchtypes)
        end if

    END SUBROUTINE Get_Contain

    SUBROUTINE IsInArea_ustr_Calculation(edgee_temp, edgew_temp, edges_temp, edgen_temp, ndm, ustr_points, ustr_vertex, ustr_ngr, ustr_n_ngr, IsInArea_ustr)
        ! 判断非结构网格顶点是否位于细化区域/包含关系计算区域内容，只要有一个顶点在就可以了
        implicit none
     
        integer :: i, j, k, n
        real(r8) :: ndm_points(4, 2), ustr(7, 2), ndm_point(2)
        real(r8), dimension(:), intent(in) :: edgee_temp, edgew_temp, edges_temp, edgen_temp
        integer,  intent(in) :: ndm, ustr_points
        integer,  dimension(:,:),intent(in) :: ustr_ngr ! 三角形/多边形的顶点编号
        integer,  dimension(:),  intent(in) :: ustr_n_ngr 
        real(r8), dimension(:,:), allocatable, intent(in) :: ustr_vertex ! 三角形/多边形顶点的的经，纬度
        integer :: num_dbx, num
        logical,  dimension(:)  , allocatable, intent(inout) :: IsInArea_ustr
        logical :: is_inside_temp(4)
        ! how to skip to outer layer
        ustr = 0.
        num = 0
        do n = 1, ndm, 1  ! ndm_domain
            ndm_points(1, 1) = edgee_temp(n); ndm_points(1, 2) = edgen_temp(n)
            ndm_points(2, 1) = edgew_temp(n); ndm_points(2, 2) = edgen_temp(n)
            ndm_points(3, 1) = edgew_temp(n); ndm_points(3, 2) = edges_temp(n)
            ndm_points(4, 1) = edgee_temp(n); ndm_points(4, 2) = edges_temp(n)
            ! the first : determind the vertex of ustr in the domain or refine area
            do i = 2, ustr_points, 1 ! 第一个都是空三角形/多边形
                if (IsInArea_ustr(i) == .true.) cycle
                num_dbx = ustr_n_ngr(i)
                ustr(1:num_dbx, :) = ustr_vertex(ustr_ngr(1:num_dbx, i), :)! 直接矩阵赋值，不要循环
                do j = 1, num_dbx, 1
                    if( (ustr(j, 1) > edgew_temp(n)).and.(ustr(j, 1) < edgee_temp(n)) .and. &
                        (ustr(j, 2) > edges_temp(n)).and.(ustr(j, 2) < edgen_temp(n)) )then
                        IsInArea_ustr(i) = .true. ! 只有一个在就好了
                        num = num + 1
                        exit
                    end if
                end do
            end do

            ! the second : determind the vertex of domain/refine area in the ustr
            do i = 2, ustr_points, 1 ! 第一个都是空三角形/多边形
                if (IsInArea_ustr(i) == .true.) cycle
                num_dbx = ustr_n_ngr(i)
                ustr(1:num_dbx, :) = ustr_vertex(ustr_ngr(1:num_dbx, i), :)! 直接矩阵赋值，不要循环
                is_inside_temp = .false.
                do k = 1, 4, 1
                    ndm_point = ndm_points(k, :)
                    is_inside_temp(k) = is_point_in_convex_polygon(ustr, ndm_point, num_dbx)
                    if (is_inside_temp(k)) then
                        if (maxval(ustr(:, 1))-minval(ustr(:, 1)) > 180.) then
                            exit ! make sure ustr not cross 180 lontitude
                        end if
                        IsInArea_ustr(i) = .true.
                        num = num + 1
                        print*, "is_inside_temp is TRUE i = ",i, "k = ", k
                        print*, "ustr = ",ustr
                        exit
                    end if
                end do
            end do

        end do
        print*, "需要计算包含关系的三角形/多边形个数（含重复网格）为 = ", num
    END SUBROUTINE IsInArea_ustr_Calculation

    SUBROUTINE Contain_Calculation(ustr_points, ustr_vertex, ustr_ngr, ustr_n_ngr, IsInArea_grid, IsInArea_ustr, ustr_id_new, ustr_ii_new, numpatch_new)
        implicit none

        integer,  intent(in) :: ustr_points ! 三角形/多边形个数
        real(r8), dimension(:,:),  allocatable, intent(in) :: ustr_vertex
        integer,  dimension(:,:),  intent(in) :: ustr_ngr
        integer,  dimension(:),    intent(in) :: ustr_n_ngr !
        logical,  dimension(:,:),  intent(in) :: IsInArea_grid 
        integer,  dimension(:),    allocatable :: icl_points
        integer :: id, i, j, k, num_edges, numpatch, numpatch_select, ad1, ii
        integer :: maxlon_source, minlon_source, maxlat_source, minlat_source
        real(r8) :: ustr_move(7, 2), maxlat_m, minlat_m, point_i(2)
        character(len = 256) :: lndname
        character(LEN = 5) :: nxpc, stepc
        integer,  dimension(:),    allocatable :: ustr_id_reset
        integer,  dimension(:, :), allocatable :: Min_matrix_index
        integer,  dimension(:, :), allocatable :: ustr_id, ustr_ii
        logical,  dimension(:),    allocatable, intent(inout) :: IsInArea_ustr
        integer,  dimension(:, :), allocatable, intent(out) :: ustr_id_new, ustr_ii_new
        integer, intent(out) :: numpatch_new
        logical :: is_inside_ustr_move

        print*,"numpatch range calculate start"
        allocate(ustr_id(ustr_points, 2)); ustr_id = 0
        allocate(icl_points(ustr_points)); icl_points = 0
        allocate(Min_matrix_index(4, ustr_points)); Min_matrix_index = 0
        ! 这个部分赋值的情况下，没有被赋值的位置，数值又是多少呢？
        numpatch = 0
        maxlat_m = maxval(ustr_vertex(:, 2))
        minlat_m = minval(ustr_vertex(:, 2)) 
        print*, "minlat_m = ", minlat_m, "maxlat_m = ", maxlat_m 
        do k = 2, ustr_points, 1           ! 循环遍历三角形网格
            ! 在阈值计算的时候判断是否在细化区域内，在包含关系计算时候判断是否在包含区域内
            if (IsInArea_ustr(k) == .false.) cycle
            ustr_move = 0.
            num_edges = ustr_n_ngr(k) ! 主要是针对多边形，但是非极点三角形也是使用的
            ustr_move(1:num_edges, :) = ustr_vertex(ustr_ngr(1:num_edges, k), 1:2)! 直接矩阵赋值，不要循环
            if (num_edges == 3) then !!!!!!!! 这是针对mode = 3的情况
                if ((ustr_move(1, 2) == minlat_m) .or. (ustr_move(1, 2) ==  maxlat_m)) then
                    print*, "k = ", k, "north or south pole in the triangle jump out"
                    IsInArea_ustr(k) = .false. 
                    cycle
                end if
            else if (num_edges == 5) then !!!!!!!! 这是针对mode = 6的情况
                if ( (minval(ustr_move(:, 2)) == minlat_m) .or. (maxval(ustr_move(:, 2)) ==  maxlat_m) ) then
                    if (all(ustr_move(1:num_edges, 2) == minlat_m) .or. all(ustr_move(1:num_edges, 2) == maxlat_m)) then
                        print*, "k = ", k, "north or south pole in the dbx become polygon jump out"
                        IsInArea_ustr(k) = .false.
                        cycle
                    end if
                end if
            end if

            if (maxval(ustr_move(1:num_edges, 1)) - minval(ustr_move(1:num_edges, 1)) > 180.) then
                icl_points(k) = 1
                call CheckCrossing(num_edges, ustr_move)
            end if

            ! lon from -180 to 180   lat from 90 to -90
            CALL Source_Find(minval(ustr_move(1:num_edges, 1)), lon_vertex, 'lon', Min_matrix_index(1, k))! minlon_source
            CALL Source_Find(maxval(ustr_move(1:num_edges, 1)), lon_vertex, 'lon', Min_matrix_index(2, k))! maxlon_source
            CALL Source_Find(maxval(ustr_move(1:num_edges, 2)), lat_vertex, 'lat', Min_matrix_index(3, k))! maxlat_source
            CALL Source_Find(minval(ustr_move(1:num_edges, 2)), lat_vertex, 'lat', Min_matrix_index(4, k))! minlat_source 纬度大值反而是小的索引值
            Min_matrix_index(1, k) = max(1, Min_matrix_index(1, k) - 1)
            Min_matrix_index(3, k) = max(1, Min_matrix_index(3, k) - 1)
            ! calculate the size of ustrgrid need and adjust IsInArea_ustr
            numpatch_select = (Min_matrix_index(2, k) - Min_matrix_index(1, k)) * (Min_matrix_index(4, k) - Min_matrix_index(3, k))
            ad1 = 0
            if (icl_points(k)) then
                ! seaorland2 use for ustr_move which cross for 180 lontitude
                print*," cross 180 lontitude ustr_move"
                do i = Min_matrix_index(1, k), Min_matrix_index(2, k)-1, 1
                    if (i < nlons_source/2 + 1) then
                        ii = int(i + nlons_source/2)
                    else
                        ii = int(i - nlons_source/2)
                    end if
                    ad1 = ad1 + sum(seaorland(ii, Min_matrix_index(3, k):Min_matrix_index(4, k)-1))
                end do
            else
                ad1 = sum(seaorland(Min_matrix_index(1, k):Min_matrix_index(2, k)-1, Min_matrix_index(3, k):Min_matrix_index(4, k)-1))
            end if

            if (ad1 == 0) then
                ! print*, "k = ", k, "numpatch_select are all the sea"
                IsInArea_ustr(k) = .false.
                cycle
            end if
            ustr_id(k, 1) = min(ad1, numpatch_select)
            numpatch = numpatch + ustr_id(k, 1)
            
        end do

        ustr_id(1, 2) = 1
        do i = 2, ustr_points, 1
           ustr_id(i, 2) = ustr_id(i - 1, 2) + ustr_id(i - 1, 1)
        end do
        allocate(ustr_id_reset(ustr_points)); ustr_id_reset = ustr_id(:, 1)
        ustr_id(:, 1) = 0
        print*,"numpatch range calculate finish"

        if (mode == 3) then
            print*, "开始计算三角形网格与经纬度网格包含关系......"
        else if (mode == 6) then
            print*, "开始计算多边形网格与经纬度网格包含关系......"
        end if
        allocate(ustr_ii(numpatch, 2)); ustr_ii = 0

        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
        !$OMP PRIVATE(i, j, k, num_edges, ustr_move, minlon_source, maxlon_source, maxlat_source, minlat_source, point_i, is_inside_ustr_move)
        do k = 2, ustr_points, 1           ! 第一个三角形/多边形是不是都是空的？？？
            if (IsInArea_ustr(k) == .false.) cycle
            if (icl_points(k)) cycle
            ustr_move = 0.
            num_edges = ustr_n_ngr(k)
            ustr_move(1:num_edges, :) = ustr_vertex(ustr_ngr(1:num_edges, k), 1:2)! 直接矩阵赋值
            ! 获取三角形最小外接矩阵的经纬度索引值
            minlon_source = Min_matrix_index(1, k)
            maxlon_source = Min_matrix_index(2, k)
            maxlat_source = Min_matrix_index(3, k)
            minlat_source = Min_matrix_index(4, k)
            ! from 90 to -90 so range from maxlat_source to minlat_source
            do i = minlon_source, maxlon_source-1, 1
                do j = maxlat_source, minlat_source-1, 1 
                    if (IsInArea_grid(i, j) .eqv. .false.) cycle
                    ! 海洋网格跳过
                    if (seaorland(i, j) == 0) cycle

                    point_i = [lon_i(i), lat_i(j)]
                    is_inside_ustr_move = is_point_in_convex_polygon(ustr_move, point_i, num_edges)
                    if (is_inside_ustr_move) then
                        ustr_id(k, 1) = ustr_id(k, 1) + 1
                        id = ustr_id(k, 2) + ustr_id(k, 1) - 1
                        ustr_ii(id, 1) = i
                        ustr_ii(id, 2) = j
                        patchtypes(i, j) = k
                    end if
                end do     
            end do
            if (ustr_id(k, 1) == 0) then
                IsInArea_ustr(k) = 0
                ! print*, "k = ", k, "IsInArea_ustr(k) turn true to false"
            else
                ! print*, "k = ", k, "ustr_id finish" ! ustr_id(k, 1) = 0 is possible
            end if

        end do
        !$OMP END PARALLEL DO

        ! 用于计算上一步中出现跨域180经线的事情
        if (sum(icl_points) /= 0) then
            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
            !$OMP PRIVATE(i, j, k, num_edges, ustr_move, minlon_source, maxlon_source, maxlat_source, minlat_source, point_i, ii, is_inside_ustr_move)
            do k = 2, ustr_points, 1           ! 第一个三角形/多边形是不是都是空的？？？
                if (IsInArea_ustr(k) == .false.) cycle
                if (icl_points(k) == 0) cycle
                ustr_move = 0.
                num_edges = ustr_n_ngr(k)
                ustr_move(1:num_edges, :) = ustr_vertex(ustr_ngr(1:num_edges, k), 1:2)! 直接矩阵赋值
                call CheckCrossing(num_edges, ustr_move)
                ! 获取三角形最小外接矩阵的经纬度索引值
                minlon_source = Min_matrix_index(1, k)
                maxlon_source = Min_matrix_index(2, k)
                maxlat_source = Min_matrix_index(3, k)
                minlat_source = Min_matrix_index(4, k)
                ! force on cross 180 lontitude
                do i = minlon_source, maxlon_source-1, 1
                    do j = maxlat_source, minlat_source-1, 1
                        if (i < nlons_source/2 + 1) then
                            ii = int( i + nlons_source/2)
                        else
                            ii = int( i - nlons_source/2)
                        end if 
                        if (IsInArea_grid(ii, j) .eqv. .false.) cycle 
                        if (seaorland(ii, j) == 0) cycle
                        
                        point_i = [lon_i(i), lat_i(j)]
                        is_inside_ustr_move = is_point_in_convex_polygon(ustr_move, point_i, num_edges)
                        if (is_inside_ustr_move) then
                            ustr_id(k, 1) = ustr_id(k, 1) + 1
                            id = ustr_id(k, 2) + ustr_id(k, 1) - 1
                            ustr_ii(id, 1) = ii
                            ustr_ii(id, 2) = j
                            patchtypes(ii, j) = k
                        end if
                    end do
                end do
                if (ustr_id(k, 1) == 0) then
                    IsInArea_ustr(k) = 0
                    ! print*, "k = ", k, "IsInArea_ustr(k) turn true to false"
                else
                    ! print*, "k = ", k, "cross 180 lontitude ustr_id finish"
                end if
            end do
            !$OMP END PARALLEL DO
        end if
        deallocate(icl_points)
        ustr_move = 0.
        ! check for data range if sufficient
        do k = 2, ustr_points, 1           ! 循环遍历三角形网格
            if (ustr_id(k, 1) > ustr_id_reset(k) ) then
                print*, "k = ", k, "ustr_id(k, 1) = ", ustr_id(k, 1), "ustr_id_reset(k) = ",ustr_id_reset(k)
                num_edges = ustr_n_ngr(k)
                ustr_move(1:num_edges, :) = ustr_vertex(ustr_ngr(1:num_edges, k), 1:2)
                print*, "num_edges = ", num_edges
                print*, "ustr_move = ", ustr_move
                ! stop "stop for data range over reset in the MOD_Getcontain.F90 Line 499 please "
            end if
        end do

        if (mode == 3) then
            print*, "三角形网格包含数组计算完成......"
        else if (mode == 6) then
            print*, "多边形网格包含数组计算完成......"
        end if

         print*, "Contain_Renew start"
        ! 完成ustr_id ustr_ii数据更新! 将由于多预留空间导致的空值删去! 并重新计算序号
        allocate(ustr_id_new(ustr_points, 2))
        ustr_id_new = 0
        ustr_id_new(:, 1) = ustr_id(:, 1)
        ustr_id_new(1, 2) = 1
    
        ! 更新ustr_id_new(i, 2) 
        do i = 2, ustr_points, 1
            ustr_id_new(i, 2) = ustr_id_new(i - 1, 2) + ustr_id_new(i - 1, 1)
        end do

        numpatch_new = INT(sum(ustr_id_new(:, 1)))
        allocate(ustr_ii_new(numpatch_new, 2)); ustr_ii_new = 0
        print*,"非结构网格实际包含经纬度网格总数（不含重复网格）numpatch_new = ",numpatch_new
        do i = 2, ustr_points, 1
            if (ustr_id_new(i, 1) == 0) cycle
            ustr_ii_new(ustr_id_new(i, 2) : ustr_id_new(i, 2) + ustr_id_new(i, 1) - 1, :) &
            =   ustr_ii(    ustr_id(i, 2) : ustr_id(i, 2)     + ustr_id_new(i, 1) - 1, :)
            ! if (mp_same(i, 1) == 0) then
            !     ustr_ii_new(ustr_id_new(i, 2) : ustr_id_new(i, 2) + ustr_id_new(i, 1) - 1, :) &
            !     =   ustr_ii(    ustr_id(i, 2) : ustr_id(i, 2)     + ustr_id_new(i, 1) - 1, :)
            ! else ! 还要获取新的数据 ustr_id_LS 和 ustr_ii_LS
            !     ustr_ii_new(ustr_id_new(i, 2) :ustr_id_new(i, 2)  + ustr_id_new(i, 1) - 1, :) &
            !     =   mp_ii_LS(mp_id_LS(i, 2):mp_id_LS(i, 2)+ mp_id_LS(i, 1) - 1, :)
            ! end if
        end do
        
        ! deallocate(mp_same); 
        ! if (step >= 1 .and. step < max_iter .and. refine == .true.) then
        !     deallocate(mp_ii_LS); deallocate(mp_id_LS)
        ! end if
        deallocate(ustr_id); deallocate(ustr_ii);
        print*, "Contain_Renew finish"

    END SUBROUTINE Contain_Calculation

    SUBROUTINE CheckCrossing(num_edges, points)

        implicit none
        
        integer, intent(in) :: num_edges
        real(r8), dimension(num_edges, 2), intent(inout) :: points
        integer :: j
        do j = 1, num_edges, 1
            if(points(j, 1) < 0.)then
                points(j, 1) = points(j, 1) + 180.
            else
                points(j, 1) = points(j, 1) - 180.
            end if
        end do
    
    END SUBROUTINE CheckCrossing

    SUBROUTINE Source_Find(temp, seq_lonlat, str1, source)
        ! 找到第一个满足条件的索引值，要求进来的时候就不会出现跨越180的情况
        implicit none
    
        real(r8), intent(in) :: temp
        real(r8), dimension(:), intent(in) :: seq_lonlat
        character(LEN = 3), intent(in) :: str1
        integer :: i, gridnum_perdegree, minsource, maxsource
        integer, intent(out) :: source
    
        gridnum_perdegree = 120
        ! min(minsource) >= 1
        if( trim(str1) == 'lon')then ! 针对lon序列! 经度度是从-180到180 !注意出现跨越180经线的问题
            minsource = ( floor(temp)   - (-180)) * gridnum_perdegree 
            maxsource = ( ceiling(temp) - (-180)) * gridnum_perdegree
            minsource = max(1,              minsource-10)
            maxsource = min(1+nlons_source, maxsource+10)
            do i = minsource, maxsource, 1 ! because  
                if( temp <= seq_lonlat(i) )then ! 满足条件就退出，等号是针对极点出现的情况
                    source = i ! 找下/左边界,上/右边界都是一样的
                    return ! exit
                end if
            end do
            print*, "temp <= seq_lonlat(i) is not exist!!!!!!!!!!!!!!!"
            print*, "temp = ", temp
            print*, "minsource = ", minsource, "maxsource = ", maxsource
            print*, "seq_lonlat(minsource) = ", seq_lonlat(minsource)
            print*, "seq_lonlat(maxsource) = ", seq_lonlat(maxsource)
            source = i
        else
            minsource = ( 90 - ceiling(temp)) * gridnum_perdegree 
            maxsource = ( 90 -   floor(temp)) * gridnum_perdegree
            minsource = max(1,              minsource-10)
            maxsource = min(1+nlats_source, maxsource+10)
            do i = minsource, maxsource, 1 
                if( temp >= seq_lonlat(i) )then ! 满足条件就退出，等号是针对极点出现的情况
                    source = i ! 找下/左边界,上/右边界都是一样的
                    return ! exit
                end if
            end do
            print*, "temp >= seq_lonlat(i) is not exist!!!!!!!!!!!!!!!"
            print*, "temp = ", temp
            print*, "minsource = ", minsource, "maxsource = ", maxsource
            print*, "seq_lonlat(minsource) = ", seq_lonlat(minsource)
            print*, "seq_lonlat(maxsource) = ", seq_lonlat(maxsource)
            source = i
        end if
    
    END SUBROUTINE Source_Find


    LOGICAL FUNCTION is_point_in_convex_polygon(polygon, ndm_point, num_dbx)
        ! 判断点是否在凸多边形内
        real(r8), intent(in) :: polygon(7, 2), ndm_point(2)
        real(r8) :: p1(2), p2(2)
        integer, intent(in) :: num_dbx
        integer :: i, next_index
        real(r8) :: prev_cross, cross
        real(r8), dimension(:,:), allocatable :: polygon_select
        allocate(polygon_select(num_dbx,2))
        polygon_select = polygon(1:num_dbx, :)
        prev_cross = 0.
        is_point_in_convex_polygon = .true.

        do i = 1, num_dbx
            next_index = mod(i, num_dbx) + 1  ! 处理循环索引
            p1 = polygon_select(i, :)
            p2 = polygon_select(next_index, :)
            cross = cross_product(p1, p2, ndm_point)
            if (cross /= 0.) then
                if (prev_cross == 0.) then
                    prev_cross = cross
                else if ((prev_cross > 0.) /= (cross > 0.)) then
                    is_point_in_convex_polygon = .false.
                    deallocate(polygon_select)
                    return
                end if
            end if
        end do
        deallocate(polygon_select)

    END FUNCTION is_point_in_convex_polygon

    REAL FUNCTION cross_product(p1, p2, p3)
        ! 计算叉积: (p2 - p1) × (p3 - p1)
        real(8), intent(in) :: p1(2), p2(2), p3(2)
        cross_product = (p2(1) - p1(1)) * (p3(2) - p1(2)) - (p2(2) - p1(2)) * (p3(1) - p1(1))

    END FUNCTION cross_product

    SUBROUTINE PatchID_Save(outputfile, patchtypes)

        use NETCDF
    
        IMPLICIT NONE
    
        character(len = 256), intent(in) :: outputfile
        integer, dimension(:,:), intent(in) :: patchtypes !  记录各经纬度网格所在非结构网格序号
        integer :: ncID, dimID_lon, dimID_lat, varid(7)
        real(r8), allocatable :: lon_w(:), lon_e(:), lat_n(:), lat_s(:), longitude(:), latitude(:)

        ! 根据目标经纬度分辨率确定变量内存分配
        allocate(lon_e(nlons_source))
        allocate(lon_w(nlons_source))
        allocate(lat_n(nlats_source))
        allocate(lat_s(nlats_source))
        allocate(longitude(nlons_source))
        allocate(latitude(nlats_source))
        lon_w      = lon_vertex(1:nlons_source)
        lon_e      = lon_vertex(2:nlons_source+1)
        lat_n      = lat_vertex(1:nlats_source)
        lat_s      = lat_vertex(2:nlats_source+1)
        longitude  = (lon_w + lon_e) / 2.
        latitude   = (lat_n + lat_s) / 2.

        print*, trim(outputfile)
        CALL CHECK(NF90_CREATE(trim(outputfile), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "nlon", nlons_source, dimID_lon))
        CALL CHECK(NF90_DEF_DIM(ncID, "nlat", nlats_source, dimID_lat))
        CALL CHECK(NF90_DEF_VAR(ncID, "elmindex", NF90_INT, (/ dimID_lon, dimID_lat /), varid(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "lon_w", NF90_DOUBLE, (/ dimID_lon /), varid(2)))
        CALL CHECK(NF90_DEF_VAR(ncID, "lon_e", NF90_DOUBLE, (/ dimID_lon /), varid(3)))
        CALL CHECK(NF90_DEF_VAR(ncID, "lat_n", NF90_DOUBLE, (/ dimID_lat /), varid(4)))
        CALL CHECK(NF90_DEF_VAR(ncID, "lat_s", NF90_DOUBLE, (/ dimID_lat /), varid(5)))
        CALL CHECK(NF90_DEF_VAR(ncID, "longitude", NF90_DOUBLE, (/ dimID_lon /), varid(6)))
        CALL CHECK(NF90_DEF_VAR(ncID, "latitude", NF90_DOUBLE, (/ dimID_lat /), varid(7)))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, varid(1), patchtypes(:, :)))
        CALL CHECK(NF90_PUT_VAR(ncID, varid(2), lon_w))
        CALL CHECK(NF90_PUT_VAR(ncID, varid(3), lon_e))
        CALL CHECK(NF90_PUT_VAR(ncID, varid(4), lat_n))
        CALL CHECK(NF90_PUT_VAR(ncID, varid(5), lat_s))
        CALL CHECK(NF90_PUT_VAR(ncID, varid(6), longitude))
        CALL CHECK(NF90_PUT_VAR(ncID, varid(7), latitude))
        CALL CHECK(NF90_CLOSE(ncID))
        deallocate(lon_w, lon_e, lat_n, lat_s, longitude, latitude)
    
    END SUBROUTINE PatchID_Save

END module MOD_GetContain
