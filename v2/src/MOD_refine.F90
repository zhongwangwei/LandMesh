module MOD_refine
    USE consts_coms
    USE refine_vars
    USE NETCDF
    USE MOD_GetContain, only: CheckCrossing
    USE MOD_GetRef, only : ref_colnum, ref_sjx, ref_th
    use MOD_utilities, only : Unstructured_Mesh_Save, Unstructured_Mesh_Read ! Add by Rui Zhang
    use MOD_grid_preprocessing, only : SpringAjustment_refine, TriMeshQuality, PolyMeshQuality, GetTriangleLength, GetAngle 
    use MOD_utilities, only : CHECK
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! add by Rui Zhang !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! integer, allocatable, public :: ref_tr(:, :)    ! 最终三角形网格阈值细化情况
    ! integer, allocatable, public :: ref_pl(:, :)    ! 记录多边形是因何种阈值细化
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! add by Rui Zhang !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Contains

    SUBROUTINE refine_loop(exit_loop)

        implicit none
        logical, intent(inout) :: exit_loop ! use for exit refine
        integer :: set_dis, set_dis_in                    ! halo(step), max_transition_row(step) 
        integer :: TransitionRow_iter     
        integer :: spDimID, lpDimID, dimaID, dimbID
        integer :: varid(10), ncvarid(2), ncid, numDimID
        integer :: m, m0, m1, w, w0, w1
        integer :: m2, k1, k2, m11, m22, w11, w22, kk
        integer :: ii, i, jj, j, k, n, num_ngrmm, mm, l, num_edges
        integer :: mi, mk, ik, wj
        integer :: num_ref, num_halo                           ! 细化三角形数, halo三角形个数
        integer :: num_mp(800), num_wp(800)                    ! 记录每次细化后的m，w点数量
        integer :: iter                                        ! 网格细化次数
        integer :: num_sjx, num_dbx                            ! 细化后三角形、多边形数量
        integer :: sjx_points, lbx_points
        integer :: num_closed_curve ! 闭合边界总数
        integer :: num_bdy_refine_segment ! 分段总数
        integer :: num_ref_weak_concav, num_weak_concav_segment, num_weak_concav_pair
        integer :: tran_degree, num_end
        logical :: isexist
        integer,  allocatable :: mrl_renew(:)
        integer,  allocatable :: ref_select(:)
        integer,  allocatable :: mp_dis(:, :)                  ! 记录三角形网格m点与其他m点的距离
        real(r8), allocatable :: mp(:, :), wp(:, :)            ! 三角形、多边形网格中心点起始数据(上一步细化的结果)
        real(r8), allocatable :: mp_new(:, :), wp_new(:, :)    ! 三角形、多边形网格中心点更新数据
        real(r8), allocatable :: mp_f(:, :), wp_f(:, :)        ! 三角形、多边形网格中心点最终数据 
        real(r8), allocatable :: ref_lbx(:, :)                 ! 构成多边形的三角形细化情况，存放顶点信息
        integer,  allocatable :: n_ngrwm(:) ! 
        integer,  allocatable :: ngrmw(:, :), ngrwm(:, :) ! 用zero and one 表示顶点是否存在
        integer,  allocatable :: ngrmw_new(:, :), ngrwm_new(:, :)  ! m/w点相邻的w/m点索引(细化后)
        integer,  allocatable :: ngrmw_f(:, :), ngrwm_f(:, :)      ! m/w点相邻的w/m点索引(最终)
        integer,  allocatable :: n_ngrwm_f(:) ! n_ngrwm define in the MOD_utilities.F90
        integer,  allocatable :: ngrmm(:, :)        ! m点相邻的m点索引(细化前)
        integer,  allocatable :: mrl_new(:)         ! 三角形网格细化程度(细化后) 
        integer,  allocatable :: ngr_mrl_new(:, :)  ! 三角形网格的相邻三角形网格点的mrl(细化后)
        integer,  allocatable :: weak_concav_pair(:,:)
        integer,  allocatable :: weak_concav_segment(:,:), weak_concav_segment_old(:,:) ! 记录弱凹左右两侧所属的分段编号，以及进行LOP的tran数
        integer,  allocatable :: n_weak_concav_segment(:), n_weak_concav_segment_old(:)
        integer,  allocatable :: ref_sjx_segment(:) ! 用于记录需要细化/LOP变换的三角形
        integer,  allocatable :: ref_sjx_segment_temp(:, :), n_ref_sjx_segment_temp(:) ! 用于存储与LOP变换相关的三角形编号与个数
        integer,  allocatable :: close_curve(:,:), n_close_curve(:) ! 闭合曲线点位存储
        integer,  allocatable :: bdy_refine_segment(:,:), bdy_refine_segment_old(:,:) ! 存储分段中待细化三角形的编号
        integer,  allocatable :: n_bdy_refine_segment(:), n_bdy_refine_segment_old(:) ! 存储分段中待细化三角形个数
        integer,  allocatable :: sjx_child(:,:) ! 用于存储过渡细化中去除的父三角形与生成的子三角形的关系
        character(LEN = 256) :: lndname
        character(LEN = 5) :: nxpc, stepc, TransitionRow_iterc
        logical :: iterA                  ! 当迭代B与迭代C同时一次性通过，迭代A通过
        logical :: iterB                  ! 从三角形网格进行判断
        logical :: iterC                  ! 从多边形网格进行判断
        logical :: iterD                  ! 弱凹点判断
        logical :: isreverse

        ! read unstructure mesh
        write(io6, *)   "start to read unstructure mesh data in the Module MOD_refine in Line 55"
        ! 读取未细化初始网格数据
        write(nxpc, '(I4.4)') NXP
        write(stepc, '(I2.2)') step
        lndname = trim(file_dir) // 'gridfile/gridfile_NXP' // trim(nxpc) // '_'//trim(stepc)// '_' // trim(mode_grid) // '.nc4'
        write(io6, *)   lndname
        call Unstructured_Mesh_Read(lndname, sjx_points, lbx_points, mp, wp, ngrmw, ngrwm, n_ngrwm)
        write(io6, *)   "The unstructured grid data reading have done "
        write(io6, *)   ""
        write(io6, *)   "In total, triangular mesh number: ", sjx_points, "polygon mesh number: ", lbx_points
        write(io6, *)   ""

        iter = 1                                 ! 本次细化中的迭代次数
        num_mp(iter) = sjx_points ! 后面涉及sjx_points与num_mp，直接用num_mp(1)代替sjx_points
        num_wp(iter) = lbx_points ! 后面涉及lbx_points与num_wwp，直接用num_wp(1)代替lbx_points
        ! 存放原始数据，就是直接对读入数据cp就好了
        CALL execute_command_line('cp '//trim(lndname)//' '//trim(trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc) // "_"// trim(stepc) //"_ori.nc4"))

        write(io6, *)   "limit -180 >> 180 90 >> -90 start"
        do i = num_vertex + 1, sjx_points, 1 ! mp is the center of sjx
            if (mp(i, 1) == -180.) mp(i, 1) = 180.
        end do
        do i = num_center + 1, lbx_points, 1 ! wp is the center of lbx
            if (wp(i, 1) == -180.) wp(i, 1) = 180.
        end do
        write(io6, *)   "limit -180 >> 180 90 >> -90 finish"

        !----------------------------------------------------------------
        ! 完成ngrmm ,mrl_new, ngr_mrl_new 数据的初始化与更新
        ! 初始化，认为三角形都存在，赋值为1
        !----------------------------------------------------------------
        ! Triangle mesh refinement degree (三角形网格细化程度/方式) 1为不细化，2,4为细化, 0为三角形不存在
        allocate(mrl_new(sjx_points));                  mrl_new = 1     ! 反映三角形自身细化与否（0或者1），以及细化方法（2或者4）
        ! mp adjacent mp initial index table (m点相邻m点的初始索引表)
        allocate(ngrmm(3, sjx_points));             ngrmm = 1   ! 反映三角形的相邻三角形的邻域关系(0表示没有相邻三角形)
        ! 在FHW代码中，lbx的第一次细化的时候ngrmm初始为 1 但是在二次细化以及后面初始为0
        ! mrl of adjacent triangular mesh points of a triangular mesh (三角形网格的相邻三角形网格的mrl)
        allocate(ngr_mrl_new(3, sjx_points));           ngr_mrl_new = 1 ! 反映三角形的相邻三角形的细化程度（0,1,2,4）

        ! 更新ngrmm，获取三角形的邻域编号信息，邻域的指向关系实际上是相互的，可否用于减小计算？？
        write(io6, *)   "ngrmm renew start"
        do mi = num_vertex + 1, sjx_points, 1
            do j = 1, 3, 1
                wj = ngrmw(j, mi)
                n = n_ngrwm(wj)
                do k = 1, n, 1
                    mk = ngrwm(k, wj)
                    if (mi == mk) cycle ! 避免mi和mk是同一个三角形
                    ik = IsNgrmm(ngrmw(1:3, mi), ngrmw(1:3, mk))
                    ! i号m点的第error个w顶点的对边的另一个相邻m点索引为j
                    if (ik /= 0) ngrmm(ik, mi) = mk
                end do
            end do
        end do
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc)  // "_" // trim(stepc) // "_ngrmm.nc4"
        write(io6, *)   lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", sjx_points, spDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 3, DimaID))
        CALL CHECK(NF90_DEF_VAR(ncID, "ngrmm", NF90_INT, (/ DimaID, spDimID /), ncVarID(1)))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(1), ngrmm))
        CALL CHECK(NF90_CLOSE(ncID))
        write(io6, *)   "ngrmm renew finish"

        ! ref_lbx 前七列分别表示组成这个多边形的三角形是否被细化 分别是0, 0.5, 1
        ! ref_lbx 第八列表示该多边形是否存在被细化的情况! 最后一个位置存放细化情况，分为0, 1, 4
        allocate(ref_lbx(lbx_points, 8));      ref_lbx = 0. 
        !!!!!!!!!!!!!!!!!!!!!!!!!! add by Rui Zhang !!!!!!!!!!!!!!!!!!!!!!!!
        ! The final triangular mesh threshold refinement (最终三角形网格阈值细化情况)
        ! allocate(ref_tr(sjx_points * 4, ref_colnum)); ref_tr = 0 ! ref_tr第二列的大小应该与ref_th的大小一致
        !!!!!!!!!!!!!!!!!!!!!!!!!! add by Rui Zhang !!!!!!!!!!!!!!!!!!!!!!!!

        if (step > 1 .or. any(n_ngrwm(2:lbx_points)<5)) then
            write(io6, *)   "根据相邻距离筛选七边形附近点 或者 虚假多多边形"
            set_dis = halo(step)
            ref_sjx(1:num_vertex) = 0 ! 进一步保证N+1级细化出现在N级细化内
            write(io6, *)   "before sum(ref_sjx) = ", sum(ref_sjx)
            num_ref = 0
            allocate(ref_select(sum(ref_sjx)))
            do i = num_vertex + 1, sjx_points, 1
                if (ref_sjx(i) == 0) cycle
                num_ref = num_ref + 1
                ref_select(num_ref) = i
            end do
            
            write(io6, *)   "mp_dis calculate start"
            allocate(mp_dis(num_ref, sjx_points-num_vertex+1)); mp_dis = 9
            CALL GetTriangleDis(num_vertex, set_dis, num_ref, sjx_points, ngrmw, ngrwm, n_ngrwm, ref_select, mp_dis)
            write(io6, *)   "mp_dis calculate finish"

            do ii = 1, num_ref, 1
                i = ref_select(ii)
                do jj = 1, sjx_points-num_vertex+1 , 1
                    if (mp_dis(ii, jj) > set_dis) cycle ! 超过距离阈值跳过
                    j = jj + num_vertex - 1
                    do k = 1, 3, 1
                        ! 要求不可以遇到七边形和小于五边形（也就是加密边界）
                        if ((n_ngrwm(ngrmw(k, j)) /= 5) .and. &
                            (n_ngrwm(ngrmw(k, j)) /= 6)) then
                            ref_sjx(i) = 0
                            ! write(io6, *)   "i =", i
                        end if
                    end do
                end do
            end do
            num_ref = INT(sum(ref_sjx))                  ! 需要细化的三角形个数


            write(io6, *)   "去除孤立细化三角形前，需要细化的三角形:", num_ref
            do i = num_vertex + 1, sjx_points, 1
                if (ref_sjx(i) /= 1) cycle ! 跳过不需要细化的三角形
                if (sum(ref_sjx(ngrmm(:, i))) > 1) cycle
                ref_sjx(i) = 0
                num_ref = num_ref - 1
            end do

            write(io6, *)   "去除孤立细化三角形后，需要细化的三角形：", num_ref
            if (num_ref == 0) then
                exit_loop = .true.   
                return
            end if

            lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc)  // "_" // trim(stepc) // "_mp_dis.nc4"
            write(io6, *)   lndname
            CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
            CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", sjx_points-num_vertex+1, spDimID))
            CALL CHECK(NF90_DEF_DIM(ncID, "num_ref", num_ref, numDimID))
            CALL CHECK(NF90_DEF_VAR(ncID, "mp_dis", NF90_INT, (/ numDimID, spDimID /), ncVarID(1)))
            CALL CHECK(NF90_ENDDEF(ncID))
            CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(1), mp_dis))
            CALL CHECK(NF90_CLOSE(ncID))
            deallocate(ref_select, mp_dis) 
        end if

        !--------------------------------------------------
        ! 1.2 Preliminary refinement (one into four) 【初步细化（一分为四）】初步细化就是阈值细化
        !--------------------------------------------------
        
        allocate(ref_sjx_segment(sjx_points)); ref_sjx_segment = 0
        write(io6, *)   "Start to refine (开始初步细化)"
        write(io6, *)   "iter =", iter, "num =", num_ref ! iter：迭代次数 num_ref：需要细化的三角形个数
        iter = iter + 1 ! 相对于FHW的代码，这是加一之后的结果
        ! 先把预估的位置占出来
        num_mp(iter) = num_mp(iter - 1) + 4 * num_ref        ! 记录每次迭代后三角形数，每细化（一分为四）一个三角形，增加4个小三角形的中心点，
        num_wp(iter) = num_wp(iter - 1) + 3 * num_ref        ! 记录每次迭代后多边形数，每细化（一分为四）一个三角形，增加三个中点作为三角形的顶点
        ref_sjx_segment(num_vertex+1:sjx_points) = ref_sjx_segment(num_vertex+1:sjx_points) + ref_sjx(num_vertex+1:sjx_points)
        CALL OnedivideFour_connection(iter, sjx_points, ngrmw, ngrmm, ref_lbx, mrl_new, ngr_mrl_new)
        
        !--------------------------------------------------
        ! 2.1 进行迭代（防止细化交汇带出现冲突，一分为四）! 也是采用一分四算法 iterB 和 iterC 用三角形与多边形的角度去防止细化带的出现
        !--------------------------------------------------

        set_dis_in = max_transition_row(step) ! 过渡行行数设置
        write(io6, *)   "iterA start"
        iterA = .false.
        do while(iterA .eqv. .false.) ! 当iterA为true时，该步骤完成
            iterA = .true.    ! 判断迭代B和迭代C是否都已满足条件
            iterB = .false.   ! 从三角形网格进行判断
            iterC = .false.   ! 从多边形网格进行判断
            iterD = .false.   ! 从弱凹点处进行判断

            write(io6, *)   "iterB start" ! 从三角形网格进行判断
            do while (iterB .eqv. .false.)
                call iterB_judge(set_dis_in, sjx_points, mrl_new, ngrmm)
                num_ref = INT(sum(ref_sjx)) ! 获取需要细化的三角形个数
                if (num_ref == 0) then
                    iterB = .true.
                else
                    write(io6, *)   "iter =", iter, "num =", num_ref
                    iterA = .false.
                    iter = iter + 1
                    num_mp(iter) = num_mp(iter - 1) + 4 * num_ref
                    num_wp(iter) = num_wp(iter - 1) + 3 * num_ref
                    ref_sjx_segment(num_vertex+1:sjx_points) = ref_sjx_segment(num_vertex+1:sjx_points) + ref_sjx(num_vertex+1:sjx_points)
                    CALL OnedivideFour_connection(iter, sjx_points, ngrmw, ngrmm, ref_lbx, mrl_new, ngr_mrl_new)
                end if
            end do ! iterB
            write(io6, *)   "iterB end"

            write(io6, *)   "iterC start"! 从多边形形网格进行判断
            do while(iterC .eqv. .false.)

                call iterC_judge(1, sjx_points, lbx_points, ngrmm, ngrmw, ngrwm, n_ngrwm, mrl_new, ngr_mrl_new, ref_lbx)
                num_ref = INT(sum(ref_sjx)) ! 获取需要细化的三角形个数
                if (num_ref == 0) then
                    iterC = .true.
                else
                    write(io6, *)   "iter =", iter, "num =", num_ref
                    iterA = .false.
                    iter = iter + 1
                    num_mp(iter) = num_mp(iter - 1) + 4 * num_ref
                    num_wp(iter) = num_wp(iter - 1) + 3 * num_ref
                    ref_sjx_segment(num_vertex+1:sjx_points) = ref_sjx_segment(num_vertex+1:sjx_points) + ref_sjx(num_vertex+1:sjx_points)
                    CALL OnedivideFour_connection(iter, sjx_points, ngrmw, ngrmm, ref_lbx, mrl_new, ngr_mrl_new)
                end if

            end do ! iterC
            write(io6, *)   "iterC end"

            if (iterA .eqv. .false.) cycle
            write(io6, *)   "iterD start" ! 从三角形网格进行判断
            do while(iterD .eqv. .false.)
                ! 寻找弱凹点，有没有一种可能第一次细化就没有了弱凹点呢？
                ! New(五边形没有弱凹点，七边形也没有，只有六边形有)
                call iterD_judge(lbx_points, ngrwm, n_ngrwm, mrl_new)
                num_ref = INT(sum(ref_sjx)) ! 获取需要细化的三角形个数
                if (num_ref == 0) then
                    write(io6, *)  "no 弱凹点 in Line 277"
                    iterD = .true.
                else
                    if (.not. weak_concav_eliminate) exit
                    write(io6, *)   "iter =", iter, "num =", num_ref
                    iterA = .false.
                    iter = iter + 1
                    num_mp(iter) = num_mp(iter - 1) + 4 * num_ref
                    num_wp(iter) = num_wp(iter - 1) + 3 * num_ref
                    ref_sjx_segment(num_vertex+1:sjx_points) = ref_sjx_segment(num_vertex+1:sjx_points) + ref_sjx(num_vertex+1:sjx_points)
                    CALL OnedivideFour_connection(iter, sjx_points, ngrmw, ngrmm, ref_lbx, mrl_new, ngr_mrl_new)
                end if
            end do ! iterD
            write(io6, *)   "iterD end"
            if (.not. weak_concav_eliminate) exit

        end do ! iterA
        write(io6, *)   "iterA end"
        write(io6, *)   "细化第二步完成"
        write(io6, *)   ""

        !--------------------------------------------------
        ! 2.2 对一分四细化的网格处理并储存
        !--------------------------------------------------
        where (ref_sjx_segment(num_vertex+1:sjx_points) > 1) ref_sjx_segment(num_vertex+1:sjx_points) = 1
        num_halo = INT(sum(ref_sjx_segment)) ! 获取需要细化的三角形个数

        ! 确定halo所需要的网格个数，并确认最终*_new的数组长度
        CALL Array_length_calculation(set_dis_in, sjx_points, lbx_points, ngrwm, n_ngrwm, mrl_new, num_halo)
        write(io6, *)   "细化的三角形+外围halo三角形个数：", num_halo
        !----------------------------------------------------------------
        ! 完成mp_new, wp_new, ngrmw_new, ngrwm_new 数据的初始化与更新
        !----------------------------------------------------------------
        ! num_ref + num_halo 本质上是为三角形设计的，但是六边形增加的个数少于三角形增加个数，所以也可用于多边形
        allocate(mp_new(sjx_points + (num_halo) * 4, 2))    ; mp_new = 0.   ! The center point of the triangular grid updates the data (三角形网格中心点更新数据)
        allocate(wp_new(lbx_points + (num_halo) * 4, 2))    ; wp_new = 0.   ! Update data at center point of polygon mesh (多边形网格中心点更新数据)
        allocate(ngrmw_new(3, sjx_points + (num_halo) * 4)) ; ngrmw_new = 1 ! The adjacent wp points of mp update the index table (mp的相邻wp点更新索引表)，
        allocate(ngrwm_new(7, lbx_points + (num_halo) * 4)) ; ngrwm_new = 1 ! wp's adjacent mp points update the index table (wp的相邻mp点更新索引表)
        mp_new(1:sjx_points, :) = mp ! 三角形中心的的经，纬度
        wp_new(1:lbx_points, :) = wp ! 多边形中心的经，纬度
        ngrmw_new(:, 1:sjx_points) = ngrmw(1:3, 1:sjx_points) ! 三角形顶点的的经，纬度
        ngrwm_new(:, 1:lbx_points) = ngrwm(1:7, 1:lbx_points) ! 多边形顶点的的经，纬度

        CALL OnedivideFour_renew(iter, ngrmw, ref_sjx_segment, num_mp, num_wp, mp_new, wp_new, ngrmw_new)
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc) // "_"// trim(stepc) //"_2.nc4"
        write(io6, *)   lndname
        call Unstructured_Mesh_Save(lndname, num_mp(iter), num_wp(iter), mp_new, wp_new, ngrmw_new, ngrwm_new)
        deallocate(ref_sjx_segment)

        !if (step > 1) Istransition = .false.
        if (Istransition) then

            ! 创建ngrvv_refine数组，实现细化区域的外边界连接关系与分成若干闭合曲线
            CALL bdy_connection_refine(sjx_points, lbx_points, ngrmw, mrl_new, ngrwm, n_ngrwm, wp, num_closed_curve, close_curve, n_close_curve)
            write(io6, *)   "bdy_connection_refine finish"
            write(io6, *)   ""

            ! 创建bdy_refine_segment实现数据分段（这里并没有考虑弱凹的情况）
            CALL bdy_refine_segment_make(set_dis_in, num_closed_curve, close_curve, n_close_curve, ngrwm, n_ngrwm, mrl_new, bdy_refine_segment, n_bdy_refine_segment, num_bdy_refine_segment)
            write(io6, *)   "分段个数num_bdy_refine_segment为：", num_bdy_refine_segment
            write(io6, *)   "bdy_refine_segment_make finish"
            write(io6, *)   ""
            
            if (.not. weak_concav_eliminate) then ! 消除弱凹就肯定是拓展为六边形
                if (num_ref == 0) then
                    weak_concav_eliminate = .TRUE. ! 如果没有弱凹则直接视为弱凹拓展的情况
                else
                    if (set_dis_in == 1) stop "ERROR! set_dis_in must larger than one when weak concav isexist"
                    ! 在bdy_refine_segment分段中提起出弱凹分段
                    num_ref_weak_concav = num_ref
                    CALL weak_concav_segment_make(set_dis_in, num_bdy_refine_segment, num_ref_weak_concav, num_weak_concav_segment, num_weak_concav_pair, ngrmw, bdy_refine_segment, n_bdy_refine_segment, weak_concav_segment, n_weak_concav_segment, weak_concav_pair)
                    write(io6, *)   "num_ref_weak_concav = ", num_ref_weak_concav
                    write(io6, *)   "num_weak_concav_segment = ", num_weak_concav_segment
                    write(io6, *)   "num_weak_concav_pair = ", num_weak_concav_pair
                    write(io6, *)   "weak_concav_segment_make finish"
                    write(io6, *)   ""
                end if
            end if

            allocate(sjx_child(2, num_mp(1))); sjx_child = 0
            TransitionRow_iter = 1
            do while(TransitionRow_iter <= set_dis_in)
                write(io6, *)   "TransitionRow_iter = ", TransitionRow_iter, "in Line 364"
                write(TransitionRow_iterc, '(I1)') TransitionRow_iter

                !--------------------------------------------------
                ! 4.1 记录需要正向一分二的三角形编号与个数
                !-------------------------------------------------- 
                ref_sjx = 0
                ! 针对强凹
                if (num_bdy_refine_segment /= 0) then
                    allocate(bdy_refine_segment_old(set_dis_in, num_bdy_refine_segment)); bdy_refine_segment_old = bdy_refine_segment
                    allocate(n_bdy_refine_segment_old(num_bdy_refine_segment)); n_bdy_refine_segment_old = n_bdy_refine_segment
                    do i = 1, num_bdy_refine_segment, 1
                        if (n_bdy_refine_segment(i) == 0) cycle ! 跳过不符合的过渡等级 
                        do j = 1, set_dis_in, 1
                            if (bdy_refine_segment(j, i) == 1) exit ! 三角形不存在就跳过
                            ref_sjx(bdy_refine_segment(j, i)) = 1
                        end do
                        n_bdy_refine_segment(i) = n_bdy_refine_segment(i) - 1 ! 完成分段长度的缩减，这个放在这里好不好，或者放在其他地方呢？
                    end do
                end if

                ! 针对弱凹（融合非两端都是1和两端都是1两种情况）
                if (num_ref_weak_concav /= 0) then
                    allocate(weak_concav_segment_old(set_dis_in, num_ref_weak_concav)); weak_concav_segment_old = weak_concav_segment
                    allocate(n_weak_concav_segment_old(num_ref_weak_concav)); n_weak_concav_segment_old = n_weak_concav_segment
                    do i = 1, num_ref_weak_concav, 1
                        if (n_weak_concav_segment(i) == 0) cycle ! 跳过不符合的过渡等级 
                        do j = 1, set_dis_in, 1
                            if (weak_concav_segment(j, i) == 1) exit ! 三角形不存在就跳过
                            ref_sjx(weak_concav_segment(j, i)) = 1
                        end do
                        n_weak_concav_segment(i) = n_weak_concav_segment(i) - 1 ! 完成分段长度的缩减，这个放在这里好不好，或者放在其他地方呢？
                    end do
                end if

                num_ref = INT(sum(ref_sjx))
                if (num_ref == 0) then
                    stop "ERROR! impossible for NO 相邻三角形中只有一个三角形被细化"
                else
                    write(io6, *)   "开始细化相邻三角形经过初步细化数为1的三角形 in Line 402"
                    write(io6, *)   "iter =", iter, "num =", num_ref
                    iter = iter + 1
                    isreverse = .false. ! 正向一分为二
                    num_mp(iter) = num_mp(iter - 1) + 2 * num_ref
                    num_wp(iter) = num_wp(iter - 1) + num_ref
                    call OnedivideTwo(iter, isreverse, ngrmw, ngrmm, ngrwm, num_mp, num_wp, mp_new, wp_new, mrl_new, ngrmw_new, sjx_child)
                    do i = num_vertex + 1, num_mp(1), 1 ! 放在外面更新，要不然容易出错
                        if (ref_sjx(i) == 1) mrl_new(i) = 4
                    end do
                    lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc) // "_"// trim(stepc) //"_4_"//trim(TransitionRow_iterc) // ".nc4"
                    write(io6, *)   lndname
                    call Unstructured_Mesh_Save(lndname, num_mp(iter), num_wp(iter), mp_new, wp_new, ngrmw_new, ngrwm_new)
                end if
                write(io6, *)   "细化第四步完成"
                write(io6, *)   ""

                ! 这个需要小心，如果出现左右两侧分段长度的弱凹的话，这就会产生八边形，所有要合理的避免这种情况
                TransitionRow_iter = TransitionRow_iter + 1
                if (TransitionRow_iter > set_dis_in) cycle ! 跳过后面内容
                if (TransitionRow_iter == 3) then
                    num_weak_concav_pair = 0 
                    write(io6, *)  "num_weak_concav_pair turn to zero"
                end if

                !--------------------------------------------------
                ! 4.4 记录强凹与弱凹中需要反向一分二的三角形，并完成bdy_refine_segment的更新
                !--------------------------------------------------    
                ! 这里的TransitionRow_iter已经加过1了

                ref_sjx = 0
                ! 要求在强凹中，只对同一个分段的细化三角形确定需要反向一分二的三角形，借助ngrmm
                CALL ref_sjx_isreverse_judge(set_dis_in, num_bdy_refine_segment, ngrmm, mrl_new, bdy_refine_segment, n_bdy_refine_segment)
                CALL ref_sjx_isreverse_judge(set_dis_in, num_weak_concav_segment, ngrmm, mrl_new, weak_concav_segment, n_weak_concav_segment)
                ! 要求在弱凹中，这里只针对弱凹两端都是1的情况
                if (num_weak_concav_pair /= 0) then
                    write(io6, *)   "weak_concav_special start"
                    CALL weak_concav_special(num_weak_concav_pair, num_ref_weak_concav, ngrmm, ngrmw, mrl_new, weak_concav_pair, weak_concav_segment, n_weak_concav_segment)
                    write(io6, *)   "weak_concav_special finish"
                    write(io6, *)   ""
                end if

                num_ref = INT(sum(ref_sjx))
                if (num_ref == 0) then
                    write(io6, *)  "NO 相邻三角形之间的反向一分二细化"
                else
                    write(io6, *)   "开始相邻细化三角形之间的反向一分二细化"
                    write(io6, *)   "iter =", iter, "num =", num_ref
                    iter = iter + 1
                    isreverse = .true. ! 反向一分二
                    num_mp(iter) = num_mp(iter - 1) + 2 * num_ref
                    num_wp(iter) = num_wp(iter - 1) + num_ref
                    call OnedivideTwo(iter, isreverse, ngrmw, ngrmm, ngrwm, num_mp, num_wp, mp_new, wp_new, mrl_new, ngrmw_new, sjx_child)
                    do i = num_vertex + 1, num_mp(1), 1 ! 放在外面更新，要不然容易出错
                        if (ref_sjx(i) == 1) mrl_new(i) = 4
                    end do
                end if
                write(io6, *)   "相邻细化三角形之间的反向一分二细化完成 in Line 458"
                write(io6, *)   ""
                lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc) // "_"// trim(stepc) //"_41_"//trim(TransitionRow_iterc) // ".nc4"
                write(io6, *)   lndname
                call Unstructured_Mesh_Save(lndname, num_mp(iter), num_wp(iter), mp_new, wp_new, ngrmw_new, ngrwm_new)

                !!!!!!!!!!!!!!!!!!!!!!!!!
                ! write(io6, *)   "!!!!!!!!!!!!!!!!!! 临时修改！！！！！！！！！！！"
                ! TransitionRow_iter = TransitionRow_iter + 1
                ! if (TransitionRow_iter > set_dis_in) cycle ! 跳过后面内容
                !!!!!!!!!!!!!!!!!!!!!!!!!!!

                !--------------------------------------------------
                ! 4.5 记录需要对角变换的三角形并交换(针对强凹与弱凹有不同的处理方式)
                !-------------------------------------------------- 
                ! 需要确认需要对角变换的三角形的个数
                num_ref = 0 ! 再根据num_ref的大小确定ref_sjx_segment
                num_end = 4*(set_dis_in-1)
                write(io6, *)   "max value of num_end = ", num_end, "ref_sjx_segment_temp"
                write(io6, *)   ""
                allocate(ref_sjx_segment_temp(num_end, num_bdy_refine_segment+num_ref_weak_concav)); ref_sjx_segment_temp = 1 ! 三角形初始编号为1
                allocate(n_ref_sjx_segment_temp(num_bdy_refine_segment+num_ref_weak_concav)); n_ref_sjx_segment_temp = 0 ! 含有三角形的个数，初始为0
                n_ref_sjx_segment_temp(1:num_bdy_refine_segment) = n_bdy_refine_segment ! 已经在前面进行过减一的处理 ！！！！！！

                CALL sharp_concav_lop_judge(set_dis_in, num_ref, num_bdy_refine_segment, mrl_new, ngrmm, ngrmw_new, sjx_child, bdy_refine_segment, bdy_refine_segment_old, n_bdy_refine_segment, ref_sjx_segment_temp, n_ref_sjx_segment_temp)
                
                if (.not. weak_concav_eliminate) CALL weak_concav_lop_judge(set_dis_in, num_ref, num_bdy_refine_segment, num_ref_weak_concav, num_weak_concav_segment, num_weak_concav_pair, &
                                                mrl_new, ngrmm, ngrmw_new, sjx_child, weak_concav_segment, weak_concav_segment_old, n_weak_concav_segment, weak_concav_pair, ref_sjx_segment_temp, n_ref_sjx_segment_temp)

                if (num_ref == 0) then
                    write(io6, *)  "不需要对角交换"
                else
                    allocate(ref_sjx_segment(num_ref)); ref_sjx_segment = 1! 获取细化三角形的索引编号
                    m = 0 ! 用于推进
                    do i = 1, num_bdy_refine_segment + num_ref_weak_concav, 1
                        if (n_ref_sjx_segment_temp(i) == 0) cycle ! 跳过不存在的三角形
                        ref_sjx_segment(m+1:m+n_ref_sjx_segment_temp(i)) = ref_sjx_segment_temp(1:n_ref_sjx_segment_temp(i), i) ! 第一个是个数
                        m = m + n_ref_sjx_segment_temp(i)
                    end do

                    write(io6, *)   "开始对角变换"
                    write(io6, *)   "iter =", iter, "num =", num_ref
                    iter = iter + 1
                    ! 每次去掉两个三角形，而生成两个新三角形，不认为有新的多边形生成
                    num_mp(iter) = num_mp(iter - 1) + num_ref
                    num_wp(iter) = num_wp(iter - 1)
                    call Delaunay_Lop(iter, num_ref, num_mp, num_wp, mp_new, wp_new, ngrmw_new, ref_sjx_segment)  
                    deallocate(ref_sjx_segment)
                end if

                lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc) // "_"// trim(stepc) //"_LOP_"//trim(TransitionRow_iterc) // ".nc4"
                write(io6, *)   lndname
                call Unstructured_Mesh_Save(lndname, num_mp(iter), num_wp(iter), mp_new, wp_new, ngrmw_new, ngrwm_new)
                write(io6, *)   "对角变换完成 in Line 509"
                write(io6, *)   ""

                deallocate(bdy_refine_segment_old, n_bdy_refine_segment_old)
                if (allocated(weak_concav_segment_old)) deallocate(weak_concav_segment_old)
                if (allocated(n_weak_concav_segment_old)) deallocate(n_weak_concav_segment_old)
                deallocate(ref_sjx_segment_temp, n_ref_sjx_segment_temp)
                sjx_child = 0

                ! 将前面弱凹相关的需要一分二的三角形激活
                if (num_weak_concav_pair /= 0) then
                    do i = num_ref_weak_concav-num_weak_concav_pair+1, num_ref_weak_concav, 1 ! 有值但是个数为0，激活
                        if (n_weak_concav_segment(i) /= 0) cycle
                        if (weak_concav_segment(1, i) == 1) cycle 
                        n_weak_concav_segment(i) = 1
                    end do
                    write(io6, *)   "弱凹相关的需要一分二的三角形激活 完成"
                end if

                if (sum(n_weak_concav_segment) == 0) cycle
                if (sum(n_weak_concav_segment) == 0) then
                    write(io6, *)   "num_weak_concav_segment turn to zero"
                    num_weak_concav_segment = 0
                end if

                if (num_weak_concav_pair + num_weak_concav_segment == 0) cycle
                if (num_weak_concav_pair + num_weak_concav_segment == 0) then
                    write(io6, *)   " weak_concav_eliminate turn to TRUE"
                    weak_concav_eliminate = .TRUE.
                end if

            end do
            write(io6, *)   "过渡构建finish"
        end if
        write(io6, *)   "细化后共有", num_wp(iter), "个多边形网格"
        write(io6, *)   "细化后共有", num_mp(iter), "个三角形网格" 

        !--------------------------------------------------
        ! 5.5 存储网格数据
        !--------------------------------------------------
        CALL NGR_RENEW(iter, num_mp, num_wp, mp_new, wp_new, ngrmw_new, num_sjx, num_dbx, mp_f, wp_f, ngrmw_f, ngrwm_f, n_ngrwm_f)
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc) // "_"// trim(stepc) //"_5.nc4"
        write(io6, *)   lndname
        call Unstructured_Mesh_Save(lndname, num_sjx, num_dbx, mp_f, wp_f, ngrmw_f, ngrwm_f, n_ngrwm_f)
        
        !--------------------------------------------------
        ! 6.1 弹性调整
        !--------------------------------------------------
        write(io6, *)   "SpringAjustment_refine start"
        ! call SpringAjustment_refine(sjx_points, num_sjx, num_dbx, num_ref, ngrmw_f, ngrwm_f, n_ngrwm_f, mp_f, wp_f)
        write(io6, *)   "SpringAjustment_refine finish"

        !--------------------------------------------------
        ! 6.7 存储最终网格数据
        !--------------------------------------------------
        ! 规定180°经线圈上的经度为180° (限制经度范围)
        do i = num_vertex + 1, num_sjx, 1 ! mp is the center of sjx
            if (mp_f(i, 1) == -180.) mp_f(i, 1) = 180.
        end do
        do i = num_center + 1,  num_dbx, 1 ! wp is the center of lbx
            if (wp_f(i, 1) == -180.) wp_f(i, 1) = 180.
        end do
        write(stepc, '(I2.2)') step + 1
        CALL GetSort(n_ngrwm_f, num_dbx, mp_f, ngrwm_f)
        lndname = trim(file_dir) // 'gridfile/gridfile_NXP' // trim(nxpc) // '_'//trim(stepc)// '_' // trim(mode_grid) // '.nc4'
        write(io6, *)   lndname
        CALL Unstructured_Mesh_Save(lndname, num_sjx, num_dbx, mp_f, wp_f, ngrmw_f, ngrwm_f, n_ngrwm_f)
        num_mp_step(step) = sjx_points - (num_mp(iter) - num_sjx)
        deallocate(ref_sjx)

        ! deallocate
        deallocate(mp, wp, mp_new, wp_new, mp_f, wp_f)
        deallocate(ref_lbx)
        deallocate(n_ngrwm)
        deallocate(ngrmw, ngrwm, ngrmw_new, ngrwm_new, ngrmw_f, ngrwm_f, n_ngrwm_f)
        deallocate(ngrmm)
        deallocate(mrl_new, ngr_mrl_new)
        if (allocated(sjx_child)) deallocate(sjx_child)
        if (allocated(ref_th)) deallocate(ref_th)
        ! deallocate(ref_th, ref_tr, ref_pl)

    END SUBROUTINE refine_loop

    INTEGER FUNCTION IsNgrmm(a, b)
        ! 三角形中心点的相邻关系，0为不相邻，1，2，3为相邻的边不同
        IMPLICIT NONE
        integer, intent(in) :: a(3), b(3)! 三角形三个顶点的编号

        IsNgrmm = 0
        if ((a(1)==b(1)).or.(a(1)==b(2)).or.(a(1)==b(3))) then
            if ((a(2)==b(1)).or.(a(2)==b(2)).or.(a(2)==b(3))) then
                IsNgrmm = 3 ! 说明是顶点3的对边相邻
                return
            else if ((a(3)==b(1)).or.(a(3)==b(2)).or.(a(3)==b(3))) then
                IsNgrmm = 2 ! 说明是顶点2的对边相邻
                return
            end if
        else if ((a(2)==b(1)).or.(a(2)==b(2)).or.(a(2)==b(3))) then
            if ((a(3)==b(1)).or.(a(3)==b(2)).or.(a(3)==b(3))) then
                IsNgrmm = 1 ! 说明是顶点1的对边相邻
                return
            end if
        end if

    END FUNCTION IsNgrmm
    
    SUBROUTINE GetTriangleDis(num_vertex, set_dis, num_ref, sjx_points, ngrmw, ngrwm, n_ngrwm, ref_select, mp_dis)
        ! 这是利用三角形顶点进行的三角形距离计算，而且只计算新生成的三角形
        implicit none
        integer, intent(in) :: num_vertex, set_dis, num_ref, sjx_points
        integer, dimension(:,:), intent(in) :: ngrmw, ngrwm
        integer, dimension(:),   intent(in) :: n_ngrwm
        integer, dimension(:),   intent(in) :: ref_select
        integer :: ii, i, jj, j, k, w, l, m 
        integer :: ncid, spDimID, ncvarid, numDimID
        integer, dimension(:,:), allocatable, intent(inout) :: mp_dis
        write(io6, *)   "shape(mp_dis) = ", shape(mp_dis)
        write(io6, *)   "num_ref = ", num_ref
        write(io6, *)   "num_vertex = ", num_vertex
        write(io6, *)   "sjx_points = ", sjx_points
        !!$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !!$OMP PRIVATE(ii,i,jj,j,w,k)
        do ii = 1, num_ref, 1
            i = ref_select(ii)
            mp_dis(ii, i-num_vertex+1) = 0
            do j = 1, 3, 1
                w = ngrmw(j, i)
                do k = 1, n_ngrwm(w), 1
                    if (ngrwm(k, w) < num_vertex + 1) cycle 
                    if (mp_dis(ii, ngrwm(k, w)-num_vertex+1) == 9) mp_dis(ii, ngrwm(k, w)-num_vertex+1) = 1
                end do
            end do
        end do
        !!$OMP END PARALLEL DO
        if (set_dis == 1) return

        do m = 1, set_dis - 1, 1
            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
            !$OMP PRIVATE(ii,i,jj,j,k,w,l)
            do ii = 1, num_ref, 1
                i = ref_select(ii)
                do jj = 1, sjx_points-num_vertex+1 , 1
                    if (mp_dis(ii, jj) /= m) cycle
                    j = jj + num_vertex - 1
                    do k = 1, 3, 1
                        w = ngrmw(k, j)
                        do l = 1, n_ngrwm(w), 1
                            if (ngrwm(l, w) < num_vertex + 1) cycle 
                            if (mp_dis(ii, ngrwm(l, w)-num_vertex+1) == 9) mp_dis(ii, ngrwm(l, w)-num_vertex+1) = m + 1
                        end do
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
        end do

    END SUBROUTINE GetTriangleDis

    ! 主要是根据mrl_new进行判断
    SUBROUTINE iterB_judge(set_dis_in, sjx_points, mrl_new, ngrmm)
        ! 要继承之前已经算好的距离矩阵，减少计算量
        implicit none
        integer, intent(in)    :: set_dis_in, sjx_points
        integer, allocatable, intent(in) :: mrl_new(:), ngrmm(:,:)
        integer :: i, j, k, m1, m2, m3, hhh(5)
        integer, allocatable :: mrl_in(:), mrl_bk(:)
        hhh = [1,2,3,1,2] 
        allocate(mrl_in(sjx_points)); mrl_in = 0 ! 用于标记三角形是否需要进行细化
        
        do i = num_vertex + 1, sjx_points, 1
            if (mrl_new(i) /= 4) cycle
            do j = 1, 3, 1
                k = ngrmm(j, i)
                if (mrl_new(k) == 4) cycle
                mrl_in(ngrmm(j, i)) = mrl_in(ngrmm(j, i)) + 2
            end do
        end do

        k = 1
        do while(k < set_dis_in)
            allocate(mrl_bk(sjx_points)); mrl_bk = mrl_in ! 用于do-while循环迭代
            do i = num_vertex + 1, sjx_points, 1
                if (mrl_new(i) == 4) cycle
                if (mrl_in(i) /= 0) cycle
                if (sum(mrl_in(ngrmm(:, i))) /= 4) cycle
                do j = 1, 3, 1
                    m1 = ngrmm(hhh(j),   i)
                    m2 = ngrmm(hhh(j+1), i)
                    m3 = ngrmm(hhh(j+2), i)
                    if ((mrl_in(m1) == 2) .and. &
                        (mrl_in(m2) == 2)) then
                        mrl_bk(i) = mrl_bk(i) + 2
                        mrl_bk(m3)= mrl_bk(m3)+ 2
                        exit ! jump 
                    end if
                end do
            end do 
            k = k + 1
            mrl_in = mrl_bk
            deallocate(mrl_bk)
        end do 

        ref_sjx = 0
        do i = num_vertex + 1, sjx_points, 1
            if (mrl_new(i) == 4) cycle 
            if (mrl_in(i)  >= 4) ref_sjx(i) = 1
        end do
        deallocate(mrl_in)

    END SUBROUTINE iterB_judge
    
    SUBROUTINE iterC_judge(set_dis_in, sjx_points, lbx_points, ngrmm, ngrmw, ngrwm, n_ngrwm, mrl_new, ngr_mrl_new, ref_lbx)
        ! 分为两部分，一部分是与ref_lbx_in（也就是过渡行）无关的内容
        implicit none

        integer,  intent(in) :: set_dis_in, sjx_points, lbx_points
        integer,  allocatable, intent(in) :: ngrmm(:,:), ngrmw(:,:), ngrwm(:,:), n_ngrwm(:), mrl_new(:), ngr_mrl_new(:,:)
        real(r8), allocatable, intent(in) :: ref_lbx(:,:)
        integer :: num_edges, i, j, k, m1, m2, m3, num, hhh(5)
        integer,  allocatable :: mrl_in(:), mrl_bk(:), ngr_mrl_in(:,:), ref_lbx_sjx(:, :), ref_lbx_sjx_bk(:, :)
        real(r8), allocatable :: ref_lbx_in(:,:)
        real(r8) :: num_ref_lbx(8)
        logical :: fexist

        hhh = [1,2,3,1,2]
        allocate(mrl_in(sjx_points)); mrl_in = mrl_new
        allocate(ngr_mrl_in(3, sjx_points)); ngr_mrl_in = ngr_mrl_new
        allocate(ref_lbx_in(lbx_points, 8)); ref_lbx_in = ref_lbx
        allocate(ref_lbx_sjx(sjx_points, 3)); ref_lbx_sjx = 0
        ref_sjx = 0


        ! have no idea about ref_lbx_in(1:num_edges)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i = num_center + 1,  lbx_points, 1
            ! 尼玛就是这个遗漏了,导致一直出错！！
            num_edges = n_ngrwm(i) ! 用于统计多边形网格边数
            if (ref_lbx_in(i, 8) == 0) cycle ! 表示i号w点的相邻三角形均未被细化       
            
            ! 当多边形中心点相邻的三角形存在细化，即ref_lbx(i, 8) /= 0 可能存在弱凹点
            if (num_edges == 5) then ! 五边形没有弱凹点！！！
                ! 只可能是连续两个或者三个三角形被细化
                if (sum( mrl_new(ngrwm(1:num_edges, i)) ) > 10 ) then
                    do j = 1, num_edges, 1
                        if (mrl_new(ngrwm(j, i)) == 1) ref_sjx(ngrwm(j, i)) = 1
                    end do 
                end if

            else if (num_edges == 6) then 
                ! 可能1：是连续两个或者三个，四个三角形被细化，不需要处理
                ! 可能2：两个对角三角形被细化，中间都相隔两个没有被细化的三角形
                ! 可能3：只有一个三角形被细化
                if (sum(mrl_new(ngrwm(1:num_edges, i))) == 12) then! 两个三角形被细化 存在两种情况，相邻，隔两个（相对，同顶点）
                    do j = 1, 3, 1
                        if ((mrl_new(ngrwm(j, i)) == 4) .and. &
                            (mrl_new(ngrwm(j + 3, i)) == 4)) then! 两个被细化三角形是相对位置的
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RuiZhang !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            ! 我建议细化三角形存在对角情况的时候把剩下的三角形都细化了！！！！！！
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RuiZhang !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            ! new
                            ! do k = 1, num_edges, 1
                            !     if (mrl_new(ngrwm(j, i)) == 1) ref_sjx(ngrwm(j, i)) = 1
                            ! end do
                            ! old
                            if ((mrl_new(ngrwm(j + 1, i)) == 1) .and. &
                                (mrl_new(ngrwm(j + 2, i)) == 1)) then
                                ref_sjx(ngrwm(j+1:j+2, i)) = 1
                            end if
                        end if ! 细化一侧的三角形，将另一侧视作弱凹
                    end do
                end if

            end if ! num_edges == 6 或者 5 的问题

        end do ! i = num_center + 1, lbx_points, 1 循环
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! consider about ref_lbx_in(1:num_edges)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i = num_center + 1, lbx_points, 1
            num_edges = n_ngrwm(i) ! 用于统计多边形网格边数
            do j = 1, num_edges, 1
                ! 把ref_lbx_in改为只对射入的三角形有效果！
                m1 = ngrwm(j, i) ! 获取该相邻三角形的编号
                if (mrl_in(m1) == 4) cycle ! 如果这个三角形是细化三角形跳过
                if (sum(ngr_mrl_in(:, m1)) == 6) then ! 确保被细化的三角形在外部而不是多边形内部
                    fexist = .false. 
                    do k = 1, 3, 1
                        if (mrl_in(ngrmm(k, m1)) == 4) exit
                    end do
                    if (all(ngrwm(1:num_edges, i) /= ngrmm(k, m1))) fexist = .true.
                    if (fexist) then
                        ref_lbx_in(i, j) = 1
                        ref_lbx_sjx(m1, 1:3) = [1, i, j]
                    end if
                end if
            end do
        end do

        ! 进行过渡行的操作,重点在于修改ref_lbx_in
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        k = 1
        do while (k < set_dis_in)
            allocate(ref_lbx_sjx_bk(sjx_points, 3)); ref_lbx_sjx_bk = ref_lbx_sjx
            do i = num_vertex + 1, sjx_points, 1
                if (ref_lbx_sjx(i, 1) /= 0) cycle
                if (sum(ref_lbx_sjx(ngrmm(:, i), 1)) /= 2) cycle
                do j = 1, 3, 1
                    m1 = ngrmm(hhh(j),   i)
                    m2 = ngrmm(hhh(j+1), i)
                    m3 = ngrmm(hhh(j+2), i)
                    if ((ref_lbx_sjx(m1, 1) == 1) .and. &
                        (ref_lbx_sjx(m2, 2) == 1)) then
                        ref_lbx_in(ref_lbx_sjx_bk(m1, 2), ref_lbx_sjx_bk(m1, 3)) = 0
                        ref_lbx_in(ref_lbx_sjx_bk(m2, 2), ref_lbx_sjx_bk(m2, 3)) = 0
                        ref_lbx_sjx_bk(m1, :) = 0
                        ref_lbx_sjx_bk(m2, :) = 0
                        exit
                    end if
                end do

                ! ref_lbx_sjx_bk(m3, 1:3) 进行赋值
                do j = 1, 3, 1
                    k = ngrmw(j, m3) ! 获取三角形对应的六边形，编号为k
                    num_edges = n_ngrwm(k)
                    if (any(ngrwm(1:num_edges, k) == m3)) exit
                end do
                do j = 1, num_edges, 1
                    if (ngrwm(j, k) == m3) then
                        ref_lbx_in(k, j) = 1
                        ref_lbx_sjx_bk(m1, 1:3) = [1, k, j]
                    end if
                end do
            end do
            k = k + 1
            ref_lbx_sjx = ref_lbx_sjx_bk
            deallocate(ref_lbx_sjx)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        do i = num_center + 1, lbx_points, 1
            if (ref_lbx_in(i, 8) /= 0) cycle ! 表示i号w点的相邻三角形均未被细化
            num_edges = n_ngrwm(i) ! 用于统计多边形网格边数       
            if ((num_edges == 5) .or. (num_edges == 6)) then ! 对于edges 等于 5 或者 6 的情况 处理都是一样的   
                num_ref_lbx = ref_lbx_in(i, :)
                do j = 1, num_edges, 1
                    m1 = j
                    m2 = mod(j, num_edges) + 1
                    if ((ref_lbx_in(i, m1) == 1) .and. &
                        (ref_lbx_in(i, m2) == 1)) then ! 针对两个相邻的三角形，因此它们组成弱凹
                        num_ref_lbx(m1) = 0.5 
                        num_ref_lbx(m2) = 0.5
                    end if! 两个0.5合计1，表示要增加1条边
                end do
                if (sum(num_ref_lbx(1:num_edges)) + num_edges > 7.) then 
                    do j = 1, num_edges, 1 !New 应该不会循环到7
                        if ((ref_lbx_in(i, j) /= 0) .and. (mrl_new(ngrwm(j, i)) == 1)) then
                            ref_sjx(ngrwm(j, i)) = 1 ! 说明该三角形要细化（一分四那种）
                        end if
                    end do
                end if
            end if   ! if((num_edges == 5) .or. (num_edges == 6))then 
        end do


        do i = num_center + 1, lbx_points, 1
            if (ref_lbx_in(i, 8) == 0) cycle !    
            num_edges = n_ngrwm(i) ! 用于统计多边形网格边数
            if (num_edges == 6) then
                ! 这应该是为了避免七边形靠的太近
                ! 但是为什么只考虑六边形只有一个细化三角形的情况的呢。
                ! 如果含有两个细化三角形或者三个细化三角形呢？？
                if (sum(mrl_new(ngrwm(1:num_edges, i))) == 9) then ! 当i点相邻三角形中只有一个被细化
                    num = 2 + sum(ref_lbx_in(i, 1:num_edges))
                    ! 建议只要num > 2 即含有细化三角形的多边形不可以成为七边形
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RuiZhang !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! 细化的三角形不可以出现在七边形中
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RuiZhang !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! New
                    ! if (num >= 3) then! 因此6+2大于7时，细化所有相邻三角形，避免出现边数大于7的情况
                    ! Old
                    if (num > 3) then! 因此6+2大于7时，细化所有相邻三角形，避免出现边数大于7的情况
                        do j = 1, num_edges, 1
                            m1 = ngrwm(j, i)
                            if (mrl_new(m1) == 1) ref_sjx(m1) = 1
                        end do
                    end if
                end if
            end if
        end do ! i = num_center + 1, lbx_points, 1 循环
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RuiZhang !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 还需要考虑同一个三角形的三个顶点不可以都是七边形
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RuiZhang !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! do j = num_vertex + 1, sjx_points, 1
        !     num =  0
        !     do k = 1, 3, 1
        !         num_edges = n_ngrwm(ngrmw(k, j))
        !         if (INT(sum(ref_lbx_in(ngrmw(k, j), 1:num_edges)) + num_edges) == 7) num = num + 1
        !     end do
        !     ! if (num == 3) write(io6, *)   "同一个三角形的三个顶点都是七边形"
        !     ! if (num == 3) ref_sjx(j) = 1
        ! end do

        deallocate(mrl_in, ngr_mrl_in, ref_lbx_in)
        
    END SUBROUTINE iterC_judge

    SUBROUTINE iterD_judge(lbx_points, ngrwm, n_ngrwm, mrl_new)
        ! 需要弱凹点的算法，前提弱凹点只出现在六边形中
        implicit none

        integer, intent(in)    :: lbx_points
        integer, allocatable, intent(in) :: ngrwm(:,:), n_ngrwm(:), mrl_new(:)
        integer :: i, j, num_edges

        ref_sjx = 0
        do i = num_center + 1, lbx_points, 1
            num_edges = n_ngrwm(i)
            if (num_edges /= 6) cycle
            if (sum( mrl_new(ngrwm(1:num_edges, i)) ) /= 18) cycle 
            do j = 1, num_edges, 1
                if (mrl_new(ngrwm(j, i)) == 1) then
                    ref_sjx(ngrwm(j, i)) = 1
                    ! write(io6, *)   "i = ", ngrwm(j, i)
                    !!!!!!!!!!!!!!!!!!! add by RuiZhang !!!!!!!!!!!!!!!!!
                    ! ref_th(ngrwm(j, i), :) = -1 ! 表示因为过渡才产生的细化
                    !!!!!!!!!!!!!!!!!!! add by RuiZhang !!!!!!!!!!!!!!!!!
                end if
            end do 
        end do

    END SUBROUTINE iterD_judge

    SUBROUTINE OnedivideFour_connection(iter, sjx_points, ngrmw, ngrmm, ref_lbx, mrl_new, ngr_mrl_new)

        IMPLICIT NONE
        integer :: i, j, k
        integer,  intent(in) :: iter, sjx_points
        integer,  dimension(:, :), allocatable, intent(in) :: ngrmw, ngrmm
        real(r8), dimension(:, :), allocatable, intent(inout) :: ref_lbx
        integer,  dimension(:),    allocatable, intent(inout) :: mrl_new
        integer,  dimension(:, :), allocatable, intent(inout) :: ngr_mrl_new

        ! 需要建立refed_iter 与第几个三角形的映射关系
        do i = num_vertex + 1, sjx_points, 1
            if ((ref_sjx(i) == 0) .or. (mrl_new(i) /= 1)) cycle ! 若三角形需要细化而且还没别细化
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 连接关系调整 更新ref_lbx, mrl_new, ngr_mrl_new!!!!!!!!!!!!!!!!!!!!!!!!!
            ref_lbx(ngrmw(1:3, i), 8) = 1 ! 作用在iterC_judge

            ! 更新三角形的自身细化状态mrl_new(分为两个部分，一个是自身三角形，一个是细化生成的三角形)
            mrl_new(i) = 4 ! 原三角形网格被平均分为四份 ! 作用在iterB/C/D_judge

            ! 更新三角形的领域细化状态ngr_mrl_new ! 作用在iterC_judge
            do j = 1, 3, 1! 因为自身被细化了，所以需要说明自身的邻域三角形，指向自身的时候为4
                do k = 1, 3, 1
                    if (ngrmm(k, ngrmm(j, i)) == i) ngr_mrl_new(k, ngrmm(j, i)) = 4
                end do
            end do
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 连接关系调整 更新ref_lbx, mrl_new, ngr_mrl_new!!!!!!!!!!!!!!!!!!!!!!!!!
        end do

    END SUBROUTINE OnedivideFour_connection

    SUBROUTINE OnedivideFour_renew(iter, ngrmw, ref_sjx_segment, num_mp, num_wp, mp_new, wp_new, ngrmw_new)

        IMPLICIT NONE
        integer :: i, j, k, refed_iter
        integer :: icl, m0, w0      ! 三角形和多边形中心点序号起始索引
        real(r8) :: sjx(3, 2), newsjx(4, 2), newdbx(3, 2)
        integer,  intent(in) :: iter
        integer,  dimension(:, :), allocatable, intent(in) :: ngrmw
        integer,  dimension(:),    allocatable, intent(in) :: ref_sjx_segment
        integer,  dimension(:), intent(inout) :: num_mp, num_wp
        real(r8), dimension(:, :), allocatable, intent(inout) :: mp_new, wp_new
        integer,  dimension(:, :), allocatable, intent(inout) :: ngrmw_new

        sjx = 0.; newdbx = 0.; newsjx = 0.
        refed_iter = 0
        ! 需要建立refed_iter 与第几个三角形的映射关系
        do i = num_vertex + 1, num_mp(1), 1
            if (ref_sjx_segment(i) == 0) cycle ! 若三角形需要细化而且还没别细化
            ! write(io6, *)   "i = ", i, "OnedivideFour_renew"
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 新的三角形与多边形顶点坐标, ngrmw_new !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            icl = 0
            sjx = wp_new(ngrmw(:, i), :) 
            if (maxval(sjx(:, 1)) - minval(sjx(:, 1)) > 180.) then ! need to modify sjx first
                icl = 1
                call CheckCrossing(3, sjx)
            end if

            ! 生成新的多边形
            newdbx(1, 1:2) = (sjx(2, 1:2) + sjx(3, 1:2)) / 2. ! 顶点1对边中点为第一个三角形顶点
            newdbx(2, 1:2) = (sjx(1, 1:2) + sjx(3, 1:2)) / 2. !  顺序逆时针，
            newdbx(3, 1:2) = (sjx(1, 1:2) + sjx(2, 1:2)) / 2. ! 计算新生成的中间小三角形顶点（原三角形中点） 
            ! 生成新的三角形
            newsjx(1, 1:2) = (sjx(1, :) + newdbx(2, :) + newdbx(3, :)) / 3.
            newsjx(2, 1:2) = (sjx(2, :) + newdbx(1, :) + newdbx(3, :)) / 3.
            newsjx(3, 1:2) = (sjx(3, :) + newdbx(1, :) + newdbx(2, :)) / 3.
            newsjx(4, 1:2) = (newdbx(3, :) + newdbx(1, :) + newdbx(2, :)) / 3.

            if (icl) then ! 经度跨越修正! 将新生成m点、w点大于180°的经度减小360° 
                call CheckCrossing(4, newsjx)
                call CheckCrossing(3, newdbx)
            end if

            m0 = num_mp(1) + refed_iter * 4 ! 新三角形中心点编号基准
            w0 = num_wp(1) + refed_iter * 3 ! 新多边形中心点编号基准
            mp_new(m0 + 1 : m0 + 4, 1:2) = newsjx(:, 1:2)! 新三角形中心点经纬度(四个)
            wp_new(w0 + 1 : w0 + 3, 1:2) = newdbx(:, 1:2)! 新多边形中心点经纬度(三个)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 新的三角形与多边形顶点坐标 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! 增加第一，二，三个新三角形的剩下两个顶点的编号信息,ngrmw_new作用在ngr_renew中
            ngrmw_new(2, m0 + 1) = w0 + 3 !w3
            ngrmw_new(3, m0 + 1) = w0 + 2 !w2
            ngrmw_new(2, m0 + 2) = w0 + 1 !w1
            ngrmw_new(3, m0 + 2) = w0 + 3 !w3
            ngrmw_new(2, m0 + 3) = w0 + 2 !w2
            ngrmw_new(3, m0 + 3) = w0 + 1 !w1
            ! 增加第四个小三角形顶点编号信息
            ngrmw_new(1, m0 + 4) = w0 + 1
            ngrmw_new(2, m0 + 4) = w0 + 2
            ngrmw_new(3, m0 + 4) = w0 + 3
            ! 更新原三角形与新三角形的ngrmw_new
            do k = 1, 3, 1
                ngrmw_new(k, i) = 1 ! 1. 去掉原三角形的顶点编号信息
                ngrmw_new(1, m0 + k) = ngrmw(k, i)! 2. 增加新三角形的顶点编号信息
            end do
            refed_iter = refed_iter + 1
    

            !!!!!!!!!!!!!!!!!!!!! add by RuiZhang !!!!!!!!!!!!!!!!!!!!
            ! ref_tr(m0 + 1, :) = ref_th(i, 1:ref_colnum)
            ! ref_tr(m0 + 2, :) = ref_th(i, 1:ref_colnum)
            ! ref_tr(m0 + 3, :) = ref_th(i, 1:ref_colnum)
            !!!!!!!!!!!!!!!!!!!!! add by RuiZhang !!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !write(io6, *)   "sjx = ", sjx
            !write(io6, *)   "newsjx = ", newsjx
            !write(io6, *)   "newdbx = ", newdbx
            !write(io6, *)   "m0 = ", m0
            !write(io6, *)   "w0 = ", w0
            !write(io6, *)   "mp_new(m1:m4, 1:2) = ", mp_new(m0+1:m0+4, 1:2)
            !write(io6, *)   "wp_new(w1:w3, 1:2) = ", wp_new(w0+1:w0+3, 1:2)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end do

        ! 针对可能压着180经线的情况
        call crossline_check(iter, mp_new, wp_new, num_mp, num_wp)

    END SUBROUTINE OnedivideFour_renew

    SUBROUTINE Array_length_calculation(set_dis_in, sjx_points, lbx_points, ngrwm, n_ngrwm, mrl_new, num_halo)
        ! 确定halo所需要的网格个数，并确认最终*_new的数组长度
        implicit none
        integer, intent(in) :: set_dis_in, sjx_points, lbx_points
        integer,  allocatable, intent(in) :: ngrwm(:,:), n_ngrwm(:), mrl_new(:)
        integer, intent(inout) :: num_halo
        logical :: fexist
        integer :: i, j, k, m, num_edges, num
        integer,  allocatable  :: isbdy_refine(:), isbdy_array(:), mrl_in(:)
        integer :: ncid, spDimID, lpDimID, dimaID, dimbID, ncvarid(4)
        character(LEN = 256) :: lndname
        character(LEN = 5) :: nxpc, stepc
        write(nxpc, '(I4.4)') NXP
        write(stepc, '(I2.2)') step

        ! 如果该顶点（即该六边形是存在而且不完整，则认为是边界点，非边界点位必须满秩）
        m = 0
        allocate(mrl_in(sjx_points)); mrl_in = mrl_new ! 替代mrl_new便于后续在子程序内的赋值操作
        allocate(isbdy_array(lbx_points))
        isbdy_array = 0 ! 1表示是边界点位 ! 只适用于第一次循环中
        do i = num_center + 1, lbx_points, 1
            num_edges = n_ngrwm(i)
            num = sum(mrl_in(ngrwm(1:num_edges, i)))
            if (num == num_edges) cycle ! 多边形内没有细化三角形，跳过
            if (num == num_edges * 4) cycle ! 多边形内全是细化三角形，跳过
            isbdy_array(i) = 1
        end do
        allocate(isbdy_refine(lbx_points)); isbdy_refine = isbdy_array

        do while(m < set_dis_in)

            ! 这里应该只需要对mrl_in处理就好了
            do i = num_center + 1, lbx_points, 1
                if (isbdy_array(i) /= 1) cycle
                num_edges = n_ngrwm(i) ! 获取相连的三角形个数
                do j = 1, num_edges ,1
                    k = ngrwm(j, i) ! 获取对应的center网格编号
                    if (mrl_in(k) == 4) cycle ! 跳过细化三角形
                    mrl_in(k) = 4 ! 找到非细化网格标记为细化网格，便于下一次循环
                    num_halo = num_halo + 1 ! 统计num_halo个数
                end do
            end do

            isbdy_array = 0 ! 1表示是边界点位 ! 只适用于第一次循环中
            do i = num_center + 1, lbx_points, 1
                num_edges = n_ngrwm(i)
                num = sum(mrl_in(ngrwm(1:num_edges, i)))
                if (num == num_edges) cycle ! 多边形内没有细化三角形，跳过
                if (num == num_edges * 4) cycle ! 多边形内全是细化三角形，跳过
                isbdy_array(i) = 1
            end do
            m = m + 1
        end do

        ! 数据保存isbdy_refine, isbdy_array, mrl_new, mrl_in
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc)  // "_" // trim(stepc) // "_isbdy.nc4"
        write(io6, *)   lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", sjx_points, spDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", lbx_points, lpDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "isbdy_refine", NF90_INT, (/ lpDimID /), ncVarID(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "isbdy_array", NF90_INT, (/ lpDimID /), ncVarID(2)))
        CALL CHECK(NF90_DEF_VAR(ncID, "mrl_new", NF90_INT, (/ spDimID /), ncVarID(3)))
        CALL CHECK(NF90_DEF_VAR(ncID, "mrl_in", NF90_INT, (/ spDimID /), ncVarID(4)))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(1), isbdy_refine))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(2), isbdy_array))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(3), mrl_new))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(4), mrl_in))
        CALL CHECK(NF90_CLOSE(ncID))
        deallocate(isbdy_refine, isbdy_array, mrl_in)

    END SUBROUTINE Array_length_calculation

    SUBROUTINE bdy_connection_refine(sjx_points, lbx_points, ngrmw, mrl_new, ngrwm, n_ngrwm, wp, num_closed_curve, close_curve, n_close_curve)
        ! 获取细化细化边界的连接信息
        implicit none
        integer, intent(in) :: sjx_points, lbx_points
        integer, allocatable, intent(in) :: ngrmw(:,:), mrl_new(:)
        integer, allocatable, intent(in) :: ngrwm(:,:), n_ngrwm(:)
        real(r8), intent(in) :: wp(:,:)
        logical :: fexist
        integer :: i, j, k, m, n, im, ik
        integer :: ncid, spDimID, lpDimID, dimaID, dimbID, ncvarid(2)
        integer :: num_edges1, num1, num_edges2, num2
        integer :: bdy_num_in, num_bdy_long, bdy_num_in_save
        integer, allocatable :: ngrvv_refine(:, :), n_ngrvv_refine(:)
        integer, allocatable :: bdy_order(:)
        integer, allocatable :: bdy_ngr(:, :)
        real(r8),allocatable :: vertex_lonlat(:,:)
        integer, intent(out) :: num_closed_curve
        integer, allocatable, intent(out) :: close_curve(:,:), n_close_curve(:)

        ! 建立细化网格边界vertex-vertex的连接关系数组ngrvv_refine
        ! 通过理论分析，一个海洋网格边界理论上最多有四个走向（两真实，两虚假）
        allocate(ngrvv_refine(4, lbx_points));  ngrvv_refine = 1
        allocate(n_ngrvv_refine(lbx_points)); n_ngrvv_refine = 0
        bdy_num_in = 1 ! 空出第一个位置
        do i = num_vertex + 1, sjx_points, 1
            if (mrl_new(i) == 1) cycle ! 跳过非细化三角形
            do j = 1, 3, 1
                ! im, ik分别是两个相连的顶点编号，要求全是边界点位才进行下一步计算
                ! 如果这两个顶点对于的vertex形状不亏缺则说明三角形不在细化边界上
                im = ngrmw(j, i) ! 获取三角形顶点编号
                num_edges1 = n_ngrwm(im) ! 获取该顶点含有的三角形个数
                num1 = sum(mrl_new(ngrwm(1:num_edges1, im)))
                if (num1 == num_edges1 * 4) cycle ! 多边形内全是细化三角形，跳过

                ik = ngrmw(mod(j, 3) + 1, i) ! 获取相邻三角形顶点编号
                num_edges2 = n_ngrwm(ik) ! 获取该顶点含有的三角形个数
                num2 = sum(mrl_new(ngrwm(1:num_edges2, ik)))
                if (num2 == num_edges2 * 4) cycle ! 多边形内全是细化三角形，跳过

                bdy_num_in = bdy_num_in + 1
                n_ngrvv_refine(im) = n_ngrvv_refine(im) + 1
                n_ngrvv_refine(ik) = n_ngrvv_refine(ik) + 1
                ngrvv_refine(n_ngrvv_refine(im), im) = ik
                ngrvv_refine(n_ngrvv_refine(ik), ik) = im
                exit
            end do
        end do
        bdy_num_in_save = bdy_num_in
        write(io6, *)   "bdy_num_in_save = ", bdy_num_in_save

        ! check for n_ngrvv_refine and adjust ifneed
        allocate(vertex_lonlat(5, 2))
        do i = num_center + 1, lbx_points, 1
            ! 等于3的情况要成对出现，
            if (n_ngrvv_refine(i) == 0) cycle
            if (n_ngrvv_refine(i) == 2) cycle
            write(io6, *)   "i = ", i, "n_ngrvv_refine(i) = ", n_ngrvv_refine(i)
            write(io6, *)   "ngrvv_refine(1:4, i) = ", ngrvv_refine(1:4, i)
            write(io6, *)   ""
            ! write(io6, *)   "n_ngrvv_refine(i) must be 0 or 2! adjust scheme writing now!"
            ! stop "ERROR! n_ngrvv_refine(i) must be 0 or 2! adjust scheme writing now!"
            ! 将具有多个连接关系的ngrvv_refine数组进行调整，根据角度大小关系
            ! 而且还要小心角度相同的情况
            vertex_lonlat(1:n_ngrvv_refine(i), :) = wp(ngrvv_refine(1:n_ngrvv_refine(i), i), :)
            vertex_lonlat(5, :) = wp(i, :)
            CALL ngrvv_refine_adjust(ngrvv_refine(1:4, i), n_ngrvv_refine(i), vertex_lonlat)
            !stop "ERROR! n_ngrvv_refine(i) must be 0 or 2! adjust scheme writing now!"
        end do
        write(io6, *)   "ngrvv_refine_adjust finish"
        deallocate(vertex_lonlat)

        ! 对于每一个边界vertex，真正有效的连接关系与网格形状无关，应该有却只有两个连接的vertex
        allocate(bdy_order(bdy_num_in)); bdy_order = 1
        allocate(bdy_ngr(2, lbx_points))
        bdy_ngr = ngrvv_refine(1:2, :)
        bdy_num_in = 1
        do i = num_center + 1, lbx_points, 1
            if (n_ngrvv_refine(i) /= 2) cycle
            bdy_num_in = bdy_num_in + 1 !获取边界顶点信息
            bdy_order(bdy_num_in) = i
        end do
        write(io6, *)   "bdy_order"
        write(io6, *)   bdy_order
        if (bdy_num_in_save /= bdy_num_in) then
            write(io6, *)  "bdy_num_in_save = ", bdy_num_in_save
            write(io6, *)  "bdy_num_in = ", bdy_num_in
            stop "ERROR! bdy_num_in_save /= bdy_num_in"
        end if

        ! 获取num_closed_curve, num_bdy_long信息
        CALL bdy_connection_refine_closed_curve(bdy_num_in, bdy_order, bdy_ngr, num_closed_curve, num_bdy_long)
        
        ! 重新遍历，获取bdy_queue信息并保留！！！！
        write(io6, *)   "get close_curve and n_close_curve start"
        allocate(close_curve(num_bdy_long, num_closed_curve)); close_curve = 1
        allocate(n_close_curve(num_closed_curve)); n_close_curve = 0
        CALL bdy_connection_refine_closed_curve(bdy_num_in, bdy_order, bdy_ngr, num_closed_curve, num_bdy_long, close_curve, n_close_curve)
        write(io6, *)   "get close_curve and n_close_curve finish"

    END SUBROUTINE bdy_connection_refine

    SUBROUTINE ngrvv_refine_adjust(ngrvv_select, n_ngrvv_select, vertex_lonlat)

        IMPLICIT NONE
        integer, intent(inout) :: ngrvv_select(4), n_ngrvv_select
        real(r8), intent(in) :: vertex_lonlat(5, 2)
        integer :: n, i, j, k
        real(r8) :: sjx(3,2), angle_temp(3)
        real(r8), allocatable :: angle(:)
        integer,  allocatable :: angle_group(:,:)

        ! 形成射线计算角度
        n = INT(n_ngrvv_select*(n_ngrvv_select-1)/2) ! 
        allocate(angle(n))
        allocate(angle_group(n, 2)) !说明是那两个数
        k = 0
        do i = 1, n_ngrvv_select-1, 1
            do j = i+1, n_ngrvv_select, 1
                sjx(1, :) = vertex_lonlat(i, :)
                sjx(2, :) = vertex_lonlat(5, :)
                sjx(3, :) = vertex_lonlat(j, :)
                CALL GetAngle(3, sjx, angle_temp)! 计算多边形内角和
                k = k + 1
                angle_group(k, :) = [ngrvv_select(i), ngrvv_select(j)]
                angle(k) = angle_temp(2)
            end do
        end do

        ! 获取最大的角度对于的数据，修改结果并返回
        k = 1
        do i = 2, n, 1
            if (angle(k) > angle(1)) then
                k = i
                angle(1) =  angle(k)
            end if
        end do
        ngrvv_select(1:2) = angle_group(k, 2)
        ngrvv_select(3:4) = 0
        n_ngrvv_select = 2
        deallocate(angle, angle_group)
        
    END SUBROUTINE ngrvv_refine_adjust

    SUBROUTINE bdy_connection_refine_closed_curve(bdy_num_in, bdy_order, bdy_ngr, num_closed_curve, num_bdy_long, close_curve, n_close_curve)
        ! 进行vertex的前后连接, 形成闭合曲线
        integer, intent(in) :: bdy_num_in ! 细化边界点位总数（第一个位置是空出来的）
        integer, allocatable, intent(in) :: bdy_order(:), bdy_ngr(:, :)
        integer :: i, j
        integer :: num_points, bdy_end, ngr_select
        integer, allocatable :: bdy_queue(:), bdy_alternate(:) ! 
        integer, intent(out) :: num_closed_curve, num_bdy_long ! 闭合曲线个数，闭合曲线最大长度
        integer, allocatable, intent(inout), optional :: close_curve(:,:), n_close_curve(:) ! 闭合曲线点位与闭合曲线长度
        
        write(io6, *)   "bdy_num_in = ", bdy_num_in, "in Line 1288 SUBROUTINE bdy_connection_refine_closed_curve"
        allocate(bdy_queue(bdy_num_in))
        allocate(bdy_alternate(bdy_num_in)); bdy_alternate = 1 ! 1表示可以使用，0表示已经使用
        num_closed_curve = 0 ! 记录闭合曲线个数
        num_bdy_long = 0 ! 记录闭合曲线最大长度
        do while(sum(bdy_alternate) > 1)
            ! 寻找闭合曲线的起点
            num_points = 1 ! 数据初始化
            num_closed_curve = num_closed_curve + 1
            bdy_queue = 1 ! 数据初始化
            do j = 2, bdy_num_in, 1
                if (bdy_alternate(j) == 1) then ! the start of queue
                    bdy_queue(num_points) = bdy_order(j)
                    bdy_alternate(j) = 0
                    exit
                end if
            end do

            ! 开始进行vertex连接使其成为闭合曲线
            bdy_end    = bdy_ngr(2, bdy_order(j)) ! the end of queue
            ngr_select = bdy_ngr(1, bdy_order(j)) ! 获取编号, 还需要知道这个编号对应的顺序编号
            ! write(io6, *)   "start from : ", bdy_order(j)
            ! write(io6, *)   "end at : ", bdy_end
            ! write(io6, *)   "ngr_select = ", ngr_select

            do while(ngr_select /= bdy_end)
                num_points = num_points + 1
                bdy_queue(num_points)  = ngr_select
                do j = 2, bdy_num_in, 1 ! ngr_select 实际在bdy_order中的位置
                    if (bdy_order(j) == ngr_select) exit
                end do
                ! write(io6, *)   "j = ", j
                bdy_alternate(j) = 0
                do i = 1, 2, 1
                    if (bdy_ngr(i, ngr_select) == bdy_queue(num_points-1)) cycle
                    ngr_select = bdy_ngr(i, ngr_select)
                    exit ! aviod ngr_select change twice!
                end do
            end do
            write(io6, *)   ""

            num_points = num_points + 1
            bdy_queue(num_points)  = bdy_end
            do j = 2, bdy_num_in, 1
                if (bdy_order(j) == bdy_end) exit
            end do
            bdy_alternate(j) = 0
            if (num_points < 3) stop "ERROR! num_points < 3 !"

            ! bdy_queue是有顺序的，不可以随便处理，需要谨慎！
            ! 如果存在某一个参数，则执行这部分内容
            if (present(n_close_curve)) then
                n_close_curve(num_closed_curve) = num_points + 1
                close_curve(1:num_points, num_closed_curve) = bdy_queue(1:num_points)
                num_points = num_points + 1
                write(io6, *)   "num_closed_curve = ", num_closed_curve, "闭合曲线(尾部空1)点位个数为 ：", num_points
            end if

            ! num_bdy_long 更新
            if (num_points > num_bdy_long) num_bdy_long = num_points ! num_points变为最长
        end do ! do while(sum(bdy_alternate) > 1)

        if (.not. present(n_close_curve)) then
            num_bdy_long = num_bdy_long + 1
            write(io6, *)   "num_bdy_long = ", num_bdy_long, "start from two"
        end if

    END SUBROUTINE bdy_connection_refine_closed_curve
                
    SUBROUTINE bdy_refine_segment_make(set_dis_in, num_closed_curve, close_curve, n_close_curve, ngrwm, n_ngrwm, mrl_new, bdy_refine_segment, n_bdy_refine_segment, num_bdy_refine_segment)
        ! 根据mrl和ngrwm一起判断是弱凹（四个细化三角形）还是强凹（两个细化三角形），直线（三个细化三角形）和转折（一个细化三角形）
        IMPLICIT NONE
        integer, intent(in) :: set_dis_in, num_closed_curve
        integer, allocatable, intent(in) :: ngrwm(:,:), n_ngrwm(:), mrl_new(:)
        integer, allocatable, intent(inout) :: close_curve(:,:), n_close_curve(:)
        integer, allocatable, intent(out) :: bdy_refine_segment(:,:) ! 存储细化三角形分组情况
        integer, allocatable, intent(out) :: n_bdy_refine_segment(:) ! 存储每一个分段中三角形个数
        integer, intent(out) :: num_bdy_refine_segment
        integer :: i, j, k, w, num_edges, num, num_segement
        integer :: m, m1, m2, m3, num_edges1, num_edges2
        logical :: isexist
        integer, allocatable :: bdy_closed_curve(:) ! 临时存储闭合曲线，便于索引
        integer, allocatable :: segement_start_end(:, :) ! 临时存储闭合曲线分段信息起点与终点
        integer, allocatable :: bdy_refine_segment_temp(:, :), n_bdy_refine_segment_temp(:)
        ! 需要对边界相对顺序进行适当的调整，要求闭合曲线的起点左右两侧不在直线上

        allocate(bdy_refine_segment_temp(set_dis_in, sum(n_close_curve))); bdy_refine_segment_temp = 1 ! 初始三角形编号
        allocate(n_bdy_refine_segment_temp(sum(n_close_curve))); n_bdy_refine_segment_temp = 0 ! 标记每一分段含有的非细化三角形的个数
        num_bdy_refine_segment = 0 ! 标记分段个数
        do i = 1, num_closed_curve, 1
            if (set_dis_in == 1) then
                allocate(segement_start_end(n_close_curve(i), 2)); segement_start_end = 0 ! 初始为0
                segement_start_end(1:n_close_curve(i)-1, 1) = [1:n_close_curve(i)-1]
                segement_start_end(1:n_close_curve(i)-1, 2) = [2:n_close_curve(i)]
                close_curve(n_close_curve(i), i) = close_curve(1, i)
            else
                ! 需要对边界相对顺序进行适当的调整，要求闭合曲线的起点不是直线，闭合曲线
                do j = 1, n_close_curve(i)-1, 1 ! 因为最后一个数据还没有赋值
                    k = close_curve(j, i) ! 获取对应的三角形坐标编号
                    num_edges = n_ngrwm(k) ! 获取边界点位对应的三角形个数
                    num = sum(mrl_new(ngrwm(1:num_edges, k)))
                    if (INT((num - num_edges)/3) /= 3) exit ! 找到第一个非直线点位就好
                end do

                if (j /= 1) then
                    ! 所有要在移动前去掉首尾重复数据，再移动后再重新考虑首尾重复
                    ! 根据获取的j调整close_curve和bdy_closed_curve,数据范围是1到n_close_curve(i)-1
                    allocate(bdy_closed_curve(n_close_curve(i)))
                    bdy_closed_curve = close_curve(1:n_close_curve(i), i)
                    write(io6, *)   "j = ", j, "need to modify order of close_curve(1:n_close_curve(i), i)"
                    close_curve(n_close_curve(i)-j:n_close_curve(i)-1, i) = bdy_closed_curve(1:j)
                    close_curve(1:n_close_curve(i)-j-1, i) = bdy_closed_curve(j+1:n_close_curve(i)-1)
                    deallocate(bdy_closed_curve)
                end if
                close_curve(n_close_curve(i), i) = close_curve(1, i)

                allocate(segement_start_end(n_close_curve(i), 2)); segement_start_end = 0 ! 初始为0
                m = 1
                segement_start_end(m, 1) = 1
                do j = 2, n_close_curve(i)-1, 1 ! 因为第一个和最后一个数是一样的
                    k = close_curve(j, i) ! 获取对应的三角形坐标编号
                    num_edges = n_ngrwm(k) ! 获取边界点位对应的三角形个数
                    num = sum(mrl_new(ngrwm(1:num_edges, k)))
                    if (INT((num - num_edges)/3) == 3) cycle
                    segement_start_end(m, 2) = j
                    segement_start_end(j, 1) = j
                    m = j ! 标记分段起点
                end do
                segement_start_end(m, 2) = n_close_curve(i)

                ! 根据set_dis_in进一步划分
                do j = 1, n_close_curve(i)-1, 1
                    if (segement_start_end(j, 1) == 0) cycle ! 跳过
                    num = segement_start_end(j, 2) - segement_start_end(j, 1) ! 获取分段含有细化三角形的个数
                    if (num <= set_dis_in) cycle ! 如果比set_dis_in短或者一样长，跳过
                    num_segement = INT((num+1)/real(set_dis_in)) ! 确定分段个数, 至少有两个，不适合set_dis_in=1的情况
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    if (mod(num+1, set_dis_in) /= 0) num_segement = num_segement + 1
                    if (mod(num,   set_dis_in) == 0) num_segement = num_segement - 1 ! 不知道用num还是num+1
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    write(io6, *)   "num = ", num, "num_segement = ", num_segement

                    segement_start_end(j+num_segement-1, 2) = segement_start_end(j, 2)
                    do k = 1, num_segement - 1, 1
                        segement_start_end(j+k-1, 2) = segement_start_end(j+k-1, 1) + set_dis_in
                        segement_start_end(j+k,   1) = segement_start_end(j+k-1, 2)
                    end do
                    if (set_dis_in < 3) cycle
                    k = num_segement-1


                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 后面再考虑如何合理分配 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! 如果最后一个区间内的细化三角形个数小于最小数量的限制则重新分配
                    if ((segement_start_end(j+k, 2) - segement_start_end(j+k, 1)) < INT((set_dis_in+1)/2)) then ! 长度4，最大区间3，分为2+2
                        write(io6, *)   "refine sjx in the segement :", segement_start_end(j+k, 2) - segement_start_end(j+k, 1)
                        write(io6, *)   "j+k = ", j+k
                        write(io6, *)   "segement_start_end(j+k, 2) = ", segement_start_end(j+k, 2)
                        write(io6, *)   "segement_start_end(j+k, 1) = ", segement_start_end(j+k, 1)
                        write(io6, *)   "不满足最小区间要求",INT((set_dis_in+1)/2),"，数组修改"
                        write(io6, *)   ""
                        segement_start_end(j+k, 1) =  segement_start_end(j+k, 2) - INT((set_dis_in+1)/2) ! 修改分段起点
                        segement_start_end(j+k-1, 2) = segement_start_end(j+k, 1) ! 修改上一个分段的终点
                    end if
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 后面再考虑如何合理分配 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                end do
            end if
            ! write(io6, *)   "segement_start_end = "
            ! write(io6, *)   segement_start_end

            ! 根据segement_start_end找到对应的三角形编号并存入bdy_refine_segment_temp
            ! 这里num不应该出现负数，而且最后一个num为1才对，说明前面segement_start_end有误
            do j = 1, n_close_curve(i)-1, 1
                if (segement_start_end(j, 1) == 0) cycle ! 跳过
                num = segement_start_end(j, 2) - segement_start_end(j, 1) ! 获取分段含有细化三角形的个数
                num_bdy_refine_segment = num_bdy_refine_segment + 1 ! 向前推进
                n_bdy_refine_segment_temp(num_bdy_refine_segment) = num
                do k = 1, num, 1
                    m1 = close_curve(segement_start_end(j, 1) + k - 1, i) ! 获取边界点位编号1
                    m2 = close_curve(segement_start_end(j, 1) + k, i) ! 获取边界点位编号2
                    num_edges1 = n_ngrwm(m1) ! 获取边界点位1对应的三角形个数
                    isexist = .false.
                    do m = 1, num_edges1, 1
                        if (mrl_new(ngrwm(m, m1)) == 4) cycle ! 跳过已经细化的三角形
                        m3 = ngrwm(m, m1) ! 获取非细化的三角形编号
                        num_edges2 = n_ngrwm(m2) ! 获取边界点位1对应的三角形个数
                        do w = 1, num_edges2, 1
                            if (mrl_new(ngrwm(w, m2)) == 4) cycle ! 跳过已经细化的三角形
                            if (m3 == ngrwm(w, m2)) then
                                isexist = .true.
                                exit
                            end if
                        end do
                        if (isexist) exit
                    end do
                    bdy_refine_segment_temp(k, num_bdy_refine_segment) = m3
                end do
            end do
            deallocate(segement_start_end)
        end do

        allocate(bdy_refine_segment(set_dis_in, num_bdy_refine_segment))
        allocate(n_bdy_refine_segment(num_bdy_refine_segment))
        bdy_refine_segment = bdy_refine_segment_temp(:, 1:num_bdy_refine_segment)
        n_bdy_refine_segment = n_bdy_refine_segment_temp(1:num_bdy_refine_segment)
        deallocate(bdy_refine_segment_temp, n_bdy_refine_segment_temp)
        ! write(io6, *)   "bdy_refine_segment = "
        ! write(io6, *)   bdy_refine_segment

    END SUBROUTINE bdy_refine_segment_make

    SUBROUTINE weak_concav_segment_make(set_dis_in, num_bdy_refine_segment, num_ref_weak_concav, num_weak_concav_segment, num_weak_concav_pair, ngrmw, bdy_refine_segment, n_bdy_refine_segment, weak_concav_segment, n_weak_concav_segment, weak_concav_pair)
        ! 实现弱凹三角形的标记，弱凹三角形所在分段的标记，弱凹三角形分段中均只有一个三角形的标记与处理
        IMPLICIT NONE
        integer, intent(in) :: set_dis_in, num_ref_weak_concav, num_bdy_refine_segment
        integer, intent(out) :: num_weak_concav_segment, num_weak_concav_pair
        integer, allocatable, intent(in) :: ngrmw(:,:)
        integer, allocatable, intent(inout) :: bdy_refine_segment(:,:) ! 存储细化三角形分组情况
        integer, allocatable, intent(inout) :: n_bdy_refine_segment(:) ! 存储每一个分段中三角形个数
        integer, allocatable, intent(out) :: weak_concav_segment(:,:) ! 存储细化三角形分组情况
        integer, allocatable, intent(out) :: n_weak_concav_segment(:) ! 存储每一个分段中三角形个数
        integer, allocatable, intent(out) :: weak_concav_pair(:, :)
        integer :: i, j, m, w, ik
        integer :: num_bdy_refine_segment_temp
        integer, allocatable :: bdy_refine_segment_temp(:, :), n_bdy_refine_segment_temp(:)
        integer, allocatable :: weak_concav_segment_temp(:,:), n_weak_concav_segment_temp(:)
        integer, allocatable :: weak_concav_pair_temp(:)

        ! 先专门存储好，因为最好的数据范围可能会发生变化
        allocate(bdy_refine_segment_temp(set_dis_in, num_bdy_refine_segment)); bdy_refine_segment_temp = bdy_refine_segment
        allocate(n_bdy_refine_segment_temp(num_bdy_refine_segment)); n_bdy_refine_segment_temp = n_bdy_refine_segment
        num_bdy_refine_segment_temp = num_bdy_refine_segment

        allocate(weak_concav_segment_temp(set_dis_in, num_ref_weak_concav)); weak_concav_segment_temp = 1 ! 弱凹三角形所在分段的三角形编号，初始化为1 
        allocate(n_weak_concav_segment_temp(num_ref_weak_concav)); n_weak_concav_segment_temp = 0 ! 计算分段中三角形个数，初始化为0
        allocate(weak_concav_pair_temp(num_ref_weak_concav)); weak_concav_pair_temp = 1
        num_weak_concav_segment = 0 ! 用于记录弱凹左右两侧长度相同，而且不为1的情况
        num_weak_concav_pair = 0 ! 用于记录弱凹左右两侧长度为1的情况

        do i = 1, num_bdy_refine_segment, 1
            j = mod(i, num_bdy_refine_segment) + 1 ! 获取相邻下一分段信息
            m = bdy_refine_segment(n_bdy_refine_segment(i), i) ! 获取前一个分段最后一个三角形
            w = bdy_refine_segment(1, j) ! 获取后一个分段第一个三角形
            ik = IsNgrmm(ngrmw(1:3, m), ngrmw(1:3, w))
            if (ik == 0) cycle ! 说明三角形不是对偶弱凹

            if (n_bdy_refine_segment(i) == n_bdy_refine_segment(j)) then ! 分为长度为1和长度不为1两种情况
                if (n_bdy_refine_segment(i) == 1) then
                    weak_concav_pair_temp(num_weak_concav_pair + 1) = m ! 记录弱凹三角形编号   
                    weak_concav_pair_temp(num_weak_concav_pair + 2) = w ! 记录弱凹三角形编号  
                    num_weak_concav_pair = num_weak_concav_pair + 2 ! 针对长度为1的情况进行特殊处理
                else
                    weak_concav_segment_temp(:, num_weak_concav_segment+1) = bdy_refine_segment(:, i)
                    weak_concav_segment_temp(:, num_weak_concav_segment+2) = bdy_refine_segment(:, j)
                    n_weak_concav_segment_temp(num_weak_concav_segment+1:num_weak_concav_segment+2) = n_bdy_refine_segment(i)
                    num_weak_concav_segment = num_weak_concav_segment + 2 ! 针对两侧长度一致的处理
                end if
                bdy_refine_segment_temp(:, [i,j]) = 1
                n_bdy_refine_segment_temp([i,j]) = 0

            else ! 长度不一致，分为有一侧长度为1，或者长度都不为1两种情况，但是都会出现一行的弱凹三角形
                weak_concav_segment_temp(1, num_weak_concav_segment+1) = bdy_refine_segment(n_bdy_refine_segment(i), i)
                weak_concav_segment_temp(1, num_weak_concav_segment+2) = bdy_refine_segment(1, j)
                n_weak_concav_segment_temp(num_weak_concav_segment+1:num_weak_concav_segment+2) = 1
                num_weak_concav_segment = num_weak_concav_segment + 2 ! 针对两侧长度一致的处理

                if ((n_bdy_refine_segment(i) /= 1) .and. (n_bdy_refine_segment(j) /= 1)) then
                    if (n_bdy_refine_segment(i) > n_bdy_refine_segment(j)) then
                        bdy_refine_segment_temp(n_bdy_refine_segment_temp(i), i) = 1
                        n_bdy_refine_segment_temp(i) = n_bdy_refine_segment_temp(i) - 1
                    else
                        n_bdy_refine_segment_temp(j) = n_bdy_refine_segment_temp(j) - 1
                        bdy_refine_segment_temp(1:n_bdy_refine_segment_temp(j), i) = bdy_refine_segment(2:n_bdy_refine_segment(j), i)
                        bdy_refine_segment_temp(n_bdy_refine_segment(j), i) = 1
                    end if
                end if
            end if
        end do
        if (num_ref_weak_concav /= (num_weak_concav_segment + num_weak_concav_pair)) then
            stop "ERROR! num_ref_weak_concav /= (num_weak_concav_segment + num_weak_concav_pair) in SUBROUTINE weak_concav_segment_make"
        end if
        bdy_refine_segment = bdy_refine_segment_temp
        n_bdy_refine_segment = n_bdy_refine_segment_temp

        allocate(weak_concav_segment(set_dis_in, num_ref_weak_concav))
        allocate(n_weak_concav_segment(num_ref_weak_concav))
        weak_concav_segment = weak_concav_segment_temp
        n_weak_concav_segment = n_weak_concav_segment_temp

        if (num_weak_concav_pair /= 0) then
            allocate(weak_concav_pair(2, num_weak_concav_pair))
            weak_concav_pair(1, 1:num_weak_concav_pair) = weak_concav_pair_temp(1:num_weak_concav_pair)
            weak_concav_segment(1, num_weak_concav_segment+1:num_ref_weak_concav) = weak_concav_pair(1, 1:num_weak_concav_pair)
            n_weak_concav_segment(num_weak_concav_segment+1:num_ref_weak_concav) = 1
        end if
        deallocate(bdy_refine_segment_temp, n_bdy_refine_segment_temp, weak_concav_segment_temp, n_weak_concav_segment_temp, weak_concav_pair_temp)

    END SUBROUTINE weak_concav_segment_make

    SUBROUTINE OnedivideTwo(iter, isreverse, ngrmw, ngrmm, ngrwm, num_mp, num_wp, mp_new, wp_new, mrl_new, ngrmw_new, sjx_child)
        ! 一分二算法（针对细化向非细化的过渡）
        IMPLICIT NONE
        ! 内部自变量
        integer :: i, j, k, icl, num_ref, refed_iter 
        integer :: m1, m2, w1, w2, w3, w4     
        real(r8) :: sjx(3, 2), tempa(1,2),tempb(1,2), tempc(1,2), hhh(5)! 
        ! 外部读入变量
        integer,  intent(in) :: iter
        logical,  intent(in) :: isreverse
        integer,  dimension(:, :), allocatable, intent(in) :: ngrwm, ngrmm, ngrmw
        integer,  dimension(:), intent(inout) :: num_mp, num_wp
        real(r8), dimension(:, :), allocatable, intent(inout) :: mp_new, wp_new
        integer,  dimension(:),    allocatable, intent(inout) :: mrl_new
        integer,  dimension(:, :), allocatable, intent(inout) :: ngrmw_new
        integer,  dimension(:, :), allocatable, intent(inout) :: sjx_child

        ! 开始细化
        hhh = [1,2,3,1,2]
        refed_iter = 0
        do i = num_vertex + 1, num_mp(1), 1
            if (ref_sjx(i) == 0) cycle
            k = 0
            if (.not. isreverse) then ! 找到邻域中唯一一个被细化的三角形
                do j = 1, 3, 1
                    if (mrl_new(ngrmm(j, i)) == 4) k = ngrmm(j, i)
                end do
            else ! 找到邻域中唯一一个没有被细化的三角形(反向一分二)
                do j = 1, 3, 1
                    if (mrl_new(ngrmm(j, i)) == 1) k = ngrmm(j, i)
                end do
            end if
            if (k == 0) stop "k==0 in Line 1714"

            do j = 1, 3, 1
                if( (ngrmw_new(j, i) /= ngrmw(1, k)) .and. & 
                    (ngrmw_new(j, i) /= ngrmw(2, k)) .and. &
                    (ngrmw_new(j, i) /= ngrmw(3, k)) )then ! 找到三角形i与细化三角公共边不相交的顶点
                    ! 根据公共边的位置，获取顶点编号，不建议用w1到w3应该w1和m1都要新增加的含义，
                    ! 但是这里只有w4是新增加的
                    w1 = ngrmw_new(hhh(j), i)
                    w2 = ngrmw_new(hhh(j+1), i)
                    w3 = ngrmw_new(hhh(j+2), i)
                    sjx(1, 1:2) = wp_new(w1, 1:2)
                    sjx(2, 1:2) = wp_new(w2, 1:2)
                    sjx(3, 1:2) = wp_new(w3, 1:2)
                    !!!!!!!!!!!!!!!!!!!!!!!! 
                    ! write(io6, *)   "w1 = ", w1
                    ! write(io6, *)   "w2 = ", w2
                    ! write(io6, *)   "w3 = ", w3
                    !!!!!!!!!!!!!!!!!!!!!!!!
                end if       
            end do
            icl = 0
            if (maxval(sjx(:, 1)) - minval(sjx(:, 1)) > 180.) then
                icl = 1
                call CheckCrossing(3, sjx)
            end if
            tempc(1, :) = (sjx(2, :) + sjx(3, :)) / 2.! 获取公共边中点的经纬度

            m1 = num_mp(iter - 1) + refed_iter * 2 + 1
            m2 = num_mp(iter - 1) + refed_iter * 2 + 2
            tempa(1, :) = (sjx(1, :) + tempc(1, :) + sjx(2, :)) / 3.! 新增加两个三角形的中心点经纬度
            tempb(1, :) = (sjx(1, :) + tempc(1, :) + sjx(3, :)) / 3.

            w4 = num_wp(iter - 1) + refed_iter + 1
            
            ngrmw_new(1, m1) = w1! 增加第一个新三角形的顶点信息
            ngrmw_new(2, m1) = w2
            ngrmw_new(3, m1) = w4
            ngrmw_new(1, m2) = w1! 增加第二个新三角形的顶点信息
            ngrmw_new(2, m2) = w3
            ngrmw_new(3, m2) = w4

            if (icl) then! 经度复原
                Call CheckCrossing(1, tempa)
                Call CheckCrossing(1, tempb)
                Call CheckCrossing(1, tempc)
            end if
            mp_new(m1, 1:2) = tempa(1, :)
            mp_new(m2, 1:2) = tempb(1, :)
            wp_new(w4, 1:2) = tempc(1, :)
            ngrmw_new(:, i) = 1
            !!!!!!!!!!!!!!!!!!!!!!!!!!! add by Rui Zhang!!!!!!!!!!!!!!!!!!!
            ! ref_tr(m1 : m2, :) = -1 
            !!!!!!!!!!!!!!!!!!!!!!!!!!! add by Rui Zhang!!!!!!!!!!!!!!!!!!!
            refed_iter = refed_iter + 1
            sjx_child(:, i) = [m1, m2]
        end do
        call crossline_check(iter, mp_new, wp_new, num_mp, num_wp)

    END SUBROUTINE OnedivideTwo

    SUBROUTINE ref_sjx_isreverse_judge(set_dis_in, num_segment, ngrmm, mrl_new, segment, n_segment)
        ! 专门处理三角形所在分段中需要反向一分二而且确定下一轮需要正向一分二的三角形(强凹与弱凹都适应)
        IMPLICIT NONE
        integer :: set_dis_in, num_segment
        integer, allocatable, intent(in) :: ngrmm(:,:), mrl_new(:)
        integer, allocatable, intent(inout) :: segment(:,:), n_segment(:)
        integer :: i, j, m0, w0, m, w, m1, w1
        logical :: isexist
        integer, allocatable :: segment_select(:)

        allocate(segment_select(set_dis_in))
        do i = 1, num_segment, 1
            if (n_segment(i) == 0) cycle ! 跳过不符合的过渡等级 
            segment_select = segment(:, i) ! 临时存储
            segment(:, i) = 1 ! 重新初始化，存储下一轮需要细化的三角形
            do j = 1, set_dis_in-1, 1
                if (segment_select(j+1) == 1) exit ! 三角形不存在就直接结束，退出
                m0 = segment_select(j)
                w0 = segment_select(j+1)
                isexist = .false.
                do m = 1, 3, 1
                    m1 = ngrmm(m, m0) ! 获取对应的相邻三角形
                    do w = 1, 3, 1
                        w1 = ngrmm(w, w0) ! 获取对应的相邻三角形
                        if (m1 == w1) then
                            isexist = .true.
                            exit
                        end if
                    end do
                    if (isexist) exit
                end do
                ref_sjx(m1) = 1 ! 获取需要反向一分为二的三角形

                ! 获取下一轮需要细化的三角形，确保只有一个三角形被细化，所以一分二的细化标记放在外部进行
                do m = 1, 3, 1
                    if (mrl_new(ngrmm(m, m1)) == 4) cycle
                    segment(j, i) = ngrmm(m, m1)
                end do
            end do
        end do
        deallocate(segment_select)

    END SUBROUTINE ref_sjx_isreverse_judge

    SUBROUTINE weak_concav_special(num_weak_concav_pair, num_ref_weak_concav, ngrmm, ngrmw, mrl_new, weak_concav_pair, weak_concav_segment, n_weak_concav_segment)
        ! 专门处理弱凹三角形所在分段中三角形个数为1的情况，确定弱凹所对应的配对三角形并确定下一次迭代中需要正向一分二的三角形
        IMPLICIT NONE
        integer :: num_weak_concav_pair, num_ref_weak_concav
        integer, allocatable, intent(in) :: ngrmm(:,:), ngrmw(:,:)
        integer, allocatable, intent(in) :: n_weak_concav_segment(:)
        integer, allocatable, intent(inout) :: mrl_new(:)
        integer, allocatable, intent(inout) :: weak_concav_pair(:, :)
        integer, allocatable, intent(inout) :: weak_concav_segment(:, :)
        integer :: i, j, k, m1, m, mm, mk
        integer, allocatable :: mrl_renew(:)

        allocate(mrl_renew(num_weak_concav_pair)); mrl_renew = 1
        do k = 1, num_weak_concav_pair, 1
            i = weak_concav_pair(1, k) ! 获取弱凹三角形
            ! 获取对偶弱凹三角形
            if (mod(k, 2) == 0) then
                j = weak_concav_pair(1, k-1)
            else
                j = weak_concav_pair(1, k+1)
            end if

            do m = 1, 3, 1
                if (mrl_new(ngrmm(m, i)) == 4) cycle
                m1 = ngrmm(m, i) ! 获取弱凹中指向外侧的三角形
                exit 
            end do
            weak_concav_pair(2, k) = m1 ! 说明该弱凹指向外侧的三角形需要LOP变换
            ref_sjx(m1) = 1 ! 该三角形细化反向一分二 ! 这个三角形可能在强凹区已经确认要一分二细化

            ! 获取下一轮需要细化的三角形，有两个非细化的三角形，将其中一个标记为已经细化
            do m = 1, 3, 1
                if (mrl_new(ngrmm(m, m1)) == 4) cycle ! 说明指向了自身
                mk = ngrmm(m, m1) ! 指向非细化的三角形
                if( (ngrmw(1, mk) /= ngrmw(1, j)) .and. & 
                    (ngrmw(2, mk) /= ngrmw(2, j)) .and. &
                    (ngrmw(3, mk) /= ngrmw(3, j)) ) then ! 找到三角形i与细化三角公共边不相交的顶点
                    ! 将这个mk记录下来，最后统一更新
                    mrl_renew(k) = mk
                else
                    mm = num_ref_weak_concav - num_weak_concav_pair + k
                    weak_concav_segment(1, mm) = mk
                end if
            end do
        end do

        do k = 1, num_weak_concav_pair, 1
            mrl_new(mrl_renew(k)) = 4
        end do
        deallocate(mrl_renew)

    END SUBROUTINE weak_concav_special

    SUBROUTINE sharp_concav_lop_judge(set_dis_in, num_ref, num_bdy_refine_segment, mrl_new, ngrmm, ngrmw_new, sjx_child, bdy_refine_segment, bdy_refine_segment_old, n_bdy_refine_segment, ref_sjx_segment_temp, n_ref_sjx_segment_temp)
        ! 利用bdy_refine_segment_old 和 bdy_refine_segment,这个才是为强凹设计的
        IMPLICIT NONE
        integer, intent(in) :: set_dis_in, num_bdy_refine_segment
        integer, intent(inout) :: num_ref
        integer, allocatable, intent(in) :: mrl_new(:), ngrmm(:, :), ngrmw_new(:, :)
        integer, allocatable, intent(in) :: sjx_child(:, :)
        integer, allocatable, intent(in) :: bdy_refine_segment(:,:), bdy_refine_segment_old(:,:)
        integer, allocatable, intent(in) :: n_bdy_refine_segment(:)
        integer, allocatable, intent(inout) :: ref_sjx_segment_temp(:,:), n_ref_sjx_segment_temp(:)
        integer :: i, j, k, w0, w1, m, m1, m2, m11, w11, m22, w22, k2 
        integer :: num_end, tran_degree
        logical :: isexist

        do i = 1, num_bdy_refine_segment, 1
            tran_degree = n_ref_sjx_segment_temp(i) + 1 ! 获取当前过渡行等级，临时加一，便于后续操作
            if (tran_degree == 1) cycle ! 跳过不符合的过渡等级 
            do j = 1, tran_degree-1, 1 ! tran_degree == 1 说明是本轮原分段中还有两个三角形
                ! 开始按顺序存储两两配对的对边三角形
                m1 = bdy_refine_segment_old(j, i)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                w0 = bdy_refine_segment(j, i)
                ! w0是下一轮需要正向一分二的三角形，他的对偶三角形才是我们需要的
                do m = 1, 3, 1
                    if (mrl_new(ngrmm(m, w0)) == 1) cycle
                    w1 = ngrmm(m, w0) ! 这才是反向一分二的三角形
                    exit
                end do
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                m2 = bdy_refine_segment_old(j+1, i)

                CALL m1w1_to_m11w11(m1, w1, sjx_child, ngrmw_new, m11, w11) ! 获取m11, w11并赋值

                ! 获取w22, m22并赋值
                w22 = sjx_child(1, w1)
                if (w22 == w11) w22 = sjx_child(2, w1)
                do k2 = 1, 2, 1
                    m22 = sjx_child(k2, m2)
                    if (IsNgrmm(ngrmw_new(1:3, w22), ngrmw_new(1:3, m22)) /= 0) exit
                end do
                ref_sjx_segment_temp(4*j-3:4*j, i) = [m11, w11, w22, m22] ! 每次放四个三角形，两两配对
            end do
            num_end = 4*(tran_degree-1) ! 定制化处理
            n_ref_sjx_segment_temp(i) = INT(tran_degree/2) * 4 ! 获取num_ref的长度
            num_ref = num_ref + n_ref_sjx_segment_temp(i)
            if (tran_degree == 2) cycle
            ! write(io6, *)   "n_ref_sjx_segment_temp(i) = ", n_ref_sjx_segment_temp(i)
            do k = 1, n_ref_sjx_segment_temp(i), 4
                ! 在相邻位置获取另一端的数据
                ! write(io6, *)   "k = ", k
                ref_sjx_segment_temp(k+2:k+3, i) =  ref_sjx_segment_temp(num_end-k:num_end-k+1, i)
            end do
            ! write(io6, *)   "after ref_sjx_segment_temp(1:num_end, i) = ", ref_sjx_segment_temp(1:num_end, i)
            ! write(io6, *)  ""
        end do

    END SUBROUTINE sharp_concav_lop_judge

    SUBROUTINE m1w1_to_m11w11(m1, w1, sjx_child, ngrmw_new, m11, w11)

        IMPLICIT NONE
        integer, intent(in) :: m1, w1
        integer, allocatable, intent(in) :: sjx_child(:,:), ngrmw_new(:,:)
        integer, intent(out) :: m11, w11
        integer :: k1, k2
        logical :: isexist

        isexist = .false.
        do k1 = 1, 2, 1
            do k2 = 1, 2, 1
                m11 = sjx_child(k1, m1)
                w11 = sjx_child(k2, w1)
                if (IsNgrmm(ngrmw_new(1:3, w11), ngrmw_new(1:3, m11)) /= 0) then
                    isexist = .true.
                    exit
                end if
            end do
            if (isexist) exit
        end do
        if (isexist .eqv. .false.) stop "ERROR! isexist .eqv. .false. in SUBROUTINE m1w1_to_m11w11"

    END SUBROUTINE m1w1_to_m11w11

    SUBROUTINE weak_concav_lop_judge(set_dis_in, num_ref, num_bdy_refine_segment, num_ref_weak_concav, num_weak_concav_segment, num_weak_concav_pair, mrl_new, ngrmm, ngrmw_new, sjx_child, &
                                    weak_concav_segment, weak_concav_segment_old, n_weak_concav_segment, weak_concav_pair, ref_sjx_segment_temp, n_ref_sjx_segment_temp)
        ! 专门针对弱凹三角形的处理，需要分为两种情况去讨论，这部分代码还需要进一步修改
        IMPLICIT NONE
        integer, intent(in) :: set_dis_in, num_bdy_refine_segment
        integer, intent(in) :: num_ref_weak_concav, num_weak_concav_segment, num_weak_concav_pair
        integer, intent(inout) :: num_ref
        integer, allocatable, intent(in) :: mrl_new(:), ngrmm(:, :), ngrmw_new(:, :)
        integer, allocatable, intent(in) :: sjx_child(:, :)
        integer, allocatable, intent(inout) :: weak_concav_segment(:,:)
        integer, allocatable, intent(in) :: weak_concav_segment_old(:,:)
        integer, allocatable, intent(in) :: n_weak_concav_segment(:)
        integer, allocatable, intent(in) :: weak_concav_pair(:,:)
        integer, allocatable, intent(inout) :: ref_sjx_segment_temp(:,:), n_ref_sjx_segment_temp(:)
        integer :: i, j, k, w0, w1, m, m1, m11, w11, kk
        integer :: num_end

        ! 针对弱凹而且左右两侧长度均以1的情况
        if (num_weak_concav_pair /= 0) then
            do i = 1, num_weak_concav_pair, 1
                m1 = weak_concav_pair(1, i)
                w1 = weak_concav_pair(2, i)
                CALL m1w1_to_m11w11(m1, w1, sjx_child, ngrmw_new, m11, w11) ! 获取m11, w11并赋值
                m = num_bdy_refine_segment+num_weak_concav_segment+i
                n_ref_sjx_segment_temp(m) = 2 ! 获取num_ref的长度
                num_ref = num_ref + n_ref_sjx_segment_temp(m)
                ref_sjx_segment_temp(1:2, m) = [m11, w11]
            end do
            num_end = num_weak_concav_segment
        else
            num_end = num_ref_weak_concav
        end if

        ! 这个可以考虑修改为多对多的弱凹细化处理哈哈哈哈(存在两个方向的问题，一个是分段内部，一个是分段与分段之间)
        if (num_weak_concav_segment /= 0) then
            do i = 1, num_end, 1
                if (weak_concav_segment(1, i) == 1) cycle ! 跳过已经不存在的三角形
                ! write(io6, *)   "i = ", i, "in Line 1889"
                m = i + num_bdy_refine_segment
                kk = 0
                ! 分段之间
                if (mod(i, 2) /= 0) then ! 将去除八边形的LOP变换的三角形放在弱凹左侧
                    m1 = weak_concav_segment_old(n_weak_concav_segment(i)+1, i) ! 弱凹左侧
                    w1 = weak_concav_segment_old(1, i+1) ! 弱凹右侧
                    CALL m1w1_to_m11w11(m1, w1, sjx_child, ngrmw_new, m11, w11) ! 获取m11, w11并赋值
                    n_ref_sjx_segment_temp(m) = 2 ! 获取num_ref的长度
                    num_ref = num_ref + 2
                    ref_sjx_segment_temp(kk+1:kk+2, m) = [m11, w11]
                    kk = kk + 2
                    if (n_weak_concav_segment(i) == 0) then
                        weak_concav_segment(:, i:i+1) = 1
                        cycle
                    end if
                end if
                ! if (n_weak_concav_segment(i) == 0) cycle
                write(io6, *)   "n_weak_concav_segment(i) = ", n_weak_concav_segment(i), "说明存在两侧长度大于1的弱凹 in Line 1917"

                ! 分段内部
                do j = 1, n_weak_concav_segment(i), 1 ! 已经进行了减一操作
                    m1 = weak_concav_segment_old(j-mod(i, 2)+1, i) ! 理论上左右两侧都适用
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    w0 = weak_concav_segment(j, i)
                    ! w0是下一轮需要正向一分二的三角形，他的对偶三角形才是我们需要的
                    do m = 1, 3, 1
                        if (mrl_new(ngrmm(m, w0)) == 1) cycle
                        w1 = ngrmm(m, w0) ! 这才是反向一分二的三角形
                        exit
                    end do
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    CALL m1w1_to_m11w11(m1, w1, sjx_child, ngrmw_new, m11, w11) ! 获取m11, w11并赋值
                    n_ref_sjx_segment_temp(m) = n_ref_sjx_segment_temp(m) + 2 ! 获取num_ref的长度
                    num_ref = num_ref + 2
                    ref_sjx_segment_temp(kk+1:kk+2, m) = [m11, w11]
                    kk = kk + 2
                end do
            end do
        end if

    END SUBROUTINE weak_concav_lop_judge

    SUBROUTINE Delaunay_Lop(iter, num_ref, num_mp, num_wp, mp_new, wp_new, ngrmw_new, ref_sjx_segment)            
        ! 对角变换
        IMPLICIT NONE
        ! 内部自变量
        integer :: i, j, k, x, icl, refed_iter
        integer :: m, m1, m2
        integer :: w, w1, w2, w3, w4      ! 三角形和多边形中心点序号起始索引
        real(r8) :: newdbx(4, 2), newsjx(2, 2)  ! 数组大小是不一样的
        ! 外部读入变量
        integer,  intent(in) :: iter, num_ref
        integer,  dimension(:), intent(inout) :: num_mp, num_wp
        real(r8), dimension(:, :), allocatable, intent(inout) :: mp_new, wp_new
        integer,  dimension(:, :), allocatable, intent(inout) :: ngrmw_new
        integer,  dimension(:),    allocatable, intent(in) :: ref_sjx_segment

        ! 开始细化弱凹点 : two adjacent triangle in a polygon need to refine 
        refed_iter = 0
        do k = 1, num_ref/2, 1
            ! if (mod(k, 2) == 1) cycle ! 只部分执行看看效果
            i = ref_sjx_segment(2*k-1)
            j = ref_sjx_segment(2*k)
            if (i==0 .or. j==0) then
                write(io6, *)   "i = ", i, "j = ", j, "in Line 1971 SUBROUTINE Delaunay_Lop"
                cycle ! 不应该出现zero，暂时不知道为什么
            end if

            do x = 1, 3, 1
                ! 判断顶点w1到w3的位置 判断顶点是否在三角形j的顶点上
                ! 只有当所有这三个条件同时满足时，整个条件语句的结果才为 true
                if ((ngrmw_new(x, i) /= ngrmw_new(1, j)) .and. &
                    (ngrmw_new(x, i) /= ngrmw_new(2, j)) .and. & 
                    (ngrmw_new(x, i) /= ngrmw_new(3, j)) ) then 
                    w1 = ngrmw_new(x, i) ! 判断为真，说明x是三角形i is relative to triangle j 的游离的顶点
                end if
            end do

            do x = 1, 3, 1
                if (w1 /= ngrmw_new(x, i)) then 
                    w2 = ngrmw_new(x, i)
                    exit ! 找到一个就退出，
                end if
            end do

            do x = 1, 3, 1
                ! 只有当所有条件同时满足时，整个条件语句的结果才为 true
                if ((w1 /= ngrmw_new(x, i)) .and. (w2 /= ngrmw_new(x, i)) ) w4 = ngrmw_new(x, i)
            end do

            do x = 1, 3, 1
                ! 判断顶点w3的位置
                if ((ngrmw_new(x, j) /= ngrmw_new(1, i)) .and. &
                    (ngrmw_new(x, j) /= ngrmw_new(2, i)) .and. & 
                    (ngrmw_new(x, j) /= ngrmw_new(3, i)) ) then 
                    w3 = ngrmw_new(x, j) ! 判断为真，说明是三角形j is relative to triangle i的游离的顶点
                end if
            end do
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            m1 = num_mp(iter - 1) + refed_iter * 2 + 1
            m2 = num_mp(iter - 1) + refed_iter * 2 + 2

            ! m1
            ngrmw_new(1, m1) = w1
            ngrmw_new(2, m1) = w2
            ngrmw_new(3, m1) = w3
            ! m2
            ngrmw_new(1, m2) = w1
            ngrmw_new(2, m2) = w4
            ngrmw_new(3, m2) = w3

            ! 赋值，尽量不适用原来的数据
            newdbx(1, :) = wp_new(w1, :)
            newdbx(2, :) = wp_new(w2, :)
            newdbx(3, :) = wp_new(w3, :)
            newdbx(4, :) = wp_new(w4, :)
            icl = 0
            if (maxval(newdbx(1:4, 1)) - minval(newdbx(1:4, 1)) > 180.) then ! need to modify sjx first
                icl = 1! 判断是否存在跨越180经线，有则修正，并返回状态变量icl
                call CheckCrossing(4, newdbx)
            end if

            ! 两个相连三角形对边交换形成新的三角形
            newsjx(1, 1:2) = (newdbx(1, 1:2) + newdbx(2, 1:2) + newdbx(3, 1:2)) / 3.
            newsjx(2, 1:2) = (newdbx(1, 1:2) + newdbx(4, 1:2) + newdbx(3, 1:2)) / 3.
            ! 经度修正(对新生成的m或者w点的经度进行矫正)
            if (icl) call CheckCrossing(2, newsjx)
            mp_new(m1:m2,:) = newsjx
            ! 将旧三角形的信息去掉
            ngrmw_new(:, i) = 1
            ngrmw_new(:, j) = 1
            refed_iter = refed_iter + 1
        end do
        call crossline_check(iter, mp_new, wp_new, num_mp, num_wp)

    END SUBROUTINE Delaunay_Lop

    SUBROUTINE crossline_check(iter, mp_new, wp_new, num_mp, num_wp)

        IMPLICIT NONE
        integer :: i
        integer,  intent(in) :: iter
        real(r8), dimension(:, :), allocatable, intent(inout) :: mp_new, wp_new
        integer,  dimension(:), intent(in) :: num_mp, num_wp
        do i = num_mp(iter - 1) + 1, num_mp(iter), 1 ! mp_new is the center of sjx
            if (mp_new(i, 1) == -180.) mp_new(i, 1) = 180.
        end do
        do i = num_wp(iter - 1) + 1, num_wp(iter), 1 ! wp_new is the center of lbx
            if (wp_new(i, 1) == -180.) wp_new(i, 1) = 180.
        end do

    END SUBROUTINE crossline_check

    SUBROUTINE NGR_RENEW(iter, num_mp, num_wp, mp_new, wp_new, ngrmw_new, num_sjx, num_dbx, mp_f, wp_f, ngrmw_f, ngrwm_f, n_ngrwm_f)
        ! 更新 mp_f, wp_f, ngrmw_f, ngrwm_f, n_ngrwm_f
        implicit none
        integer, intent(in) :: iter
        integer,  dimension(:), intent(in) :: num_mp, num_wp
        real(r8), dimension(:, :), allocatable, intent(in) :: mp_new, wp_new
        integer,  dimension(:, :), allocatable, intent(in) :: ngrmw_new
        integer, intent(out) :: num_sjx, num_dbx
        real(r8), dimension(:, :), allocatable, intent(inout) :: mp_f, wp_f
        integer,  dimension(:, :), allocatable, intent(inout) :: ngrmw_f, ngrwm_f
        integer,  dimension(:),    allocatable, intent(inout) :: n_ngrwm_f
        integer :: i, j, k, num_ref
        integer :: ncid, spDimID, lpDimID, dimaID, dimbID, ncvarid(2) 
        logical :: isexist                ! 判断细化后是否存在重复w点
        integer,  allocatable :: vertex_mapping(:)  ! 新旧顶点之间的映射关系
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! add by Rui Zhang !!!!!!!!!!!!!!!!
        ! integer,  allocatable :: ref_tr_f(:, :)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! add by Rui Zhang !!!!!!!!!!!!!!!!
        character(LEN = 256) :: lndname
        character(LEN = 5) :: nxpc, stepc
        write(nxpc, '(I4.4)') NXP
        write(stepc, '(I2.2)') step
        ! m点 增加，不会重复；减少（只在原来的范围内）m点具有唯一性
        ! w点 增加，会重复，不同编号对应同一个点位；不会减少 w点没有唯一性


        write(io6, *)   "wp_f start"
        allocate(wp_f(num_wp(iter), 2)); wp_f(:, 1:2) = 9999. ! 初始化
        allocate(vertex_mapping(num_wp(iter))); vertex_mapping = 0
        num_dbx = num_wp(1)
        wp_f(1:num_wp(1), :) = wp_new(1:num_wp(1), :)
        vertex_mapping(1:num_wp(1)) = [1:num_wp(1)]
        do i = num_wp(1) + 1, num_wp(iter), 1
            isexist = .false.
            do j = num_wp(1) + 1, num_dbx + 1, 1 ! 新点是否与原始顶点有映射关系 
                if ((wp_f(j, 1) == wp_new(i, 1)) .and. &
                    (wp_f(j, 2) == wp_new(i, 2)) ) then
                    isexist = .true.
                    exit
                end if
            end do
            if (isexist .eqv. .false.) then
                num_dbx = num_dbx + 1
                wp_f(num_dbx, 1:2) = wp_new(i, 1:2)
                vertex_mapping(i) = num_dbx
            else
                vertex_mapping(i) = j
            end if
        end do
        write(io6, *)   "细化前共有", num_wp(1), "个多边形网格"
        write(io6, *)   "细化后共有", num_wp(iter), "个多边形网格"
        write(io6, *)   "去除重复点后，还剩", num_dbx, "个多边形网格"
        write(io6, *)   "wp_f finish"
        write(io6, *)   ""
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc)  // "_" // trim(stepc) // "_wp_f.nc4"
        write(io6, *)   lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
        CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", num_wp(iter), lpDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 2, DimaID))
        CALL CHECK(NF90_DEF_VAR(ncID, "vertex_mapping", NF90_INT, (/ lpDimID /), ncVarID(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "wp_f", NF90_FLOAT, (/ lpDimID, DimaID /), ncVarID(2)))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(1), vertex_mapping))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(2), wp_f))
        CALL CHECK(NF90_CLOSE(ncID))
        write(io6, *)   "max(vertex_mapping) = ", maxval(vertex_mapping)
        if (maxval(vertex_mapping) /= num_dbx) stop "maxval(vertex_mapping) /= num_dbx"

        ! 更新 mp_f ! 如果ngrmw_new不存在，则这个三角形不存在 
        write(io6, *)   "重新计算ngrmw_new and mp_f，并储存mp_f" 
        ! m点的特点是会增加，也是减少（只在原来的范围内），但是不会重复！！！！！！！！！！！
        ! 统计初始三角形中被细化的个数，此时的ngrmw_new 还没有进行重新编号
        write(io6, *)   "mp_f start"
        num_ref = 0
        ! 这里可能需要修改，因为对角变换的时候会删去新生成的三角形
        ! do i = num_vertex + 1, num_mp(1), 1
        do i = num_vertex + 1, num_mp(iter), 1
            if (ngrmw_new(1, i) == 1) num_ref = num_ref + 1 ! 当三角顶点不存在的时候，跳过
        end do

        num_sjx = num_mp(iter) - num_ref ! 获取三角形总数
        allocate(mp_f(num_sjx, 2)); mp_f = 0. ! 经纬度
        !!!!!!!!!!!!!!!!!!!!!!!! add by RuiZhang !!!!!!!!!!!!!!!!
        ! allocate(ref_tr_f(num_sjx, ref_colnum)); ref_tr_f = 0
        !!!!!!!!!!!!!!!!!!!!!!!! add by RuiZhang !!!!!!!!!!!!!!!!
        allocate(ngrmw_f(3, num_sjx)); ngrmw_f = 1 ! 顶点编号
        ngrmw_f(:, 1:num_vertex) = ngrmw_new(:, 1:num_vertex)
        mp_f(1:num_vertex, :) = mp_new(1:num_vertex, :)
        k = num_vertex
        do i = num_vertex + 1, num_mp(iter), 1 ! 因为细化的时候会把旧的三角形去除，所有还是从2还是比较稳妥
            if (ngrmw_new(1, i) == 1) cycle ! 跳过原来三角形中被细化的那部分三角形，
            k = k + 1 ! 累加，用于进位
            mp_f(k, :) = mp_new(i, :)
            !!!!!!!!!!!!!!!!!!!!!!!! add by RuiZhang !!!!!!!!!!!!!!!!
            ! ref_tr_f(k, :) = ref_tr(i, :)
            !!!!!!!!!!!!!!!!!!!!!!!! add by RuiZhang !!!!!!!!!!!!!!!!
            ngrmw_f(:, k) = ngrmw_new(:, i)
        end do
        write(io6, *)   "细化前共有", num_mp(1), "个三角形网格"
        write(io6, *)   "细化后共有", num_mp(iter), "个三角形网格"
        write(io6, *)   "去除重复点后，还剩", num_sjx, "个三角形网格"
        write(io6, *)   "mp_f finish"
        write(io6, *)   ""
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc)  // "_" // trim(stepc) // "_mp_f.nc4"
        write(io6, *)   lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", num_sjx, spDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 2, DimaID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_b", 3, DimbID))
        CALL CHECK(NF90_DEF_VAR(ncID, "mp_f", NF90_FLOAT, (/ spDimID, DimaID /), ncVarID(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "ngrmw_f_orial", NF90_INT, (/ DimbID, spDimID /), ncVarID(2)))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(1), mp_f))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(2), ngrmw_f))
        CALL CHECK(NF90_CLOSE(ncID))
        
        ! 更新ngrmw_f和ngrwm_f
        ! allocate(ngrwm_f(7, num_dbx)); ngrwm_f   = 1 ! 记录相邻三角形编号，初始化为1
        allocate(ngrwm_f(10, num_dbx)); ngrwm_f   = 1 ! 记录相邻三角形编号，初始化为1
        allocate(n_ngrwm_f(num_dbx));  n_ngrwm_f = 0 ! 记录相邻三角形
        do i = 2, num_sjx, 1 ! 三角形总数（含不存在的三角形），这个必须从2开始
            do j = 1, 3, 1
                ngrmw_f(j, i) = vertex_mapping(ngrmw_f(j, i))
                k = ngrmw_f(j, i) ! 获取多边形中心点（或三角形顶点）编号信息
                n_ngrwm_f(k) = n_ngrwm_f(k) + 1 ! 累加顶点个数
                ngrwm_f(n_ngrwm_f(k), k) = i
            end do
        end do
        CALL GetSort(n_ngrwm_f, num_dbx, mp_f, ngrwm_f)

    END SUBROUTINE NGR_RENEW

    SUBROUTINE GetSort(n_ngr, num_dbx, mp_f, ngr)

        implicit none
        integer, dimension(:), allocatable, intent(in) :: n_ngr
        integer, intent(in) :: num_dbx
        real(r8), dimension(:,:),allocatable, intent(in) :: mp_f
        integer, dimension(:,:), allocatable, intent(inout) :: ngr ! 顶点w的相邻点m的逆时针排序输出
        integer :: i, j, k, ref_temp, num_inter
        real(r8) :: angle_temp, center(2)
        real(r8), allocatable :: angle(:), points(:, :)
        ! 这里的points和angle就要根据边数确认了，因为我们是采用向量化运算的
        do i = num_center + 1, num_dbx, 1
            num_inter = n_ngr(i) ! 获取单个多边形的顶点总数 ! use for sort points 
            allocate(angle(num_inter)); angle = 0.
            allocate(points(num_inter, 2)); points = mp_f(ngr(1:num_inter, i), 1:2)
            center  = 0.
            ! 是否跨越180经线，存在跨越
            if (maxval(points(:,1))-minval(points(:,1))>180.) call CheckCrossing(num_inter, points)

            do j = 1, num_inter, 1
                center = center + points(j, :)
            end do
            center = center / num_inter

            do j = 1, num_inter, 1
                points(j, :) = points(j, :) - center
                angle(j) = atan2(points(j, 2),  points(j, 1))
            end do

            ! 针对南北极的顶点处的纬度坐标一致的特殊情况
            if (num_inter == 5) then 
                if (maxval(points(:, 2) - minval(points(:, 2))) < 0.00001) angle = points(:, 1)
            end if
            
            !  冒泡排序
            do j = 1, num_inter - 1, 1
                do k = j + 1, num_inter, 1
                    if (angle(j) > angle(k)) then
                        angle_temp = angle(j)
                        angle(j) = angle(k)
                        angle(k) = angle_temp! 角度交换
                        ref_temp = ngr(j, i)
                        ngr(j, i) = ngr(k, i)
                        ngr(k, i) = ref_temp! 编号交换
                    end if
                end do
            end do
            deallocate(angle, points)
        end do
    
    END SUBROUTINE GetSort

END module MOD_refine
