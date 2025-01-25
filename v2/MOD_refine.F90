module MOD_refine
    USE consts_coms
    USE refine_vars
    USE NETCDF
    USE MOD_GetContain, only: CheckCrossing
    USE MOD_GetThreshold, only : ref_sjx
    use MOD_file_preprocess, only : Unstructured_Mesh_Save, Unstructured_Mesh_Read ! Add by Rui Zhang
    implicit none
    Contains

    SUBROUTINE refine_loop(exit_loop)

        implicit none
        logical, intent(inout) :: exit_loop ! use for exit refine
        integer :: set_dis                   ! 细化过渡        
        integer :: spDimID, lpDimID, dimaID, dimbID, ncvarid(2) 
        integer :: i, j, k, m, n, num_ngrmm, mm, l
        integer :: mi, mk, ik, wj
        integer :: num_ref                                     ! 每次细化三角形数 
        integer :: num_mp(200), num_wp(200)                        ! 记录每次细化后的m，w点数量
        integer :: iter                                        ! 网格细化次数
        integer :: num_sjx, num_dbx                            ! 细化后三角形、多边形数量
        integer :: ncid, varid(10)
        integer :: sjx_points, lbx_points, num_inqbx
        integer,  allocatable :: sjx_inqbx(:)                  ! 在七边形中的三角形的编号，num_inqbx是个数
        integer,  allocatable :: mp_dis(:, :)                  ! 记录三角形网格m点与其他m点的距离
        real(r8), allocatable :: mp(:, :), wp(:, :)            ! 三角形、多边形网格中心点起始数据(上一步细化的结果)
        real(r8), allocatable :: mp_new(:, :), wp_new(:, :)    ! 三角形、多边形网格中心点更新数据
        real(r8), allocatable :: mp_f(:, :), wp_f(:, :)        ! 三角形、多边形网格中心点最终数据 
        real(r8), allocatable :: ref_lbx(:, :)                 ! 构成多边形的三角形细化情况，存放顶点信息
        integer,  allocatable :: n_ngrwm(:) ! 
        integer,  allocatable :: ngrmw(:, :), ngrwm(:, :) ! 用zero and one 表示顶点是否存在
        integer,  allocatable :: ngrmw_new(:, :), ngrwm_new(:, :)  ! m/w点相邻的w/m点索引(细化后)
        integer,  allocatable :: ngrmw_f(:, :), ngrwm_f(:, :)      ! m/w点相邻的w/m点索引(最终)
        integer,  allocatable :: n_ngrwm_f(:) ! n_ngrwm define in the MOD_file_preprocess.F90
        integer,  allocatable :: ngrmm(:, :)        ! m点相邻的m点索引(细化前)
        integer,  allocatable :: ngrmm_new(:, :)    ! m点相邻的m点索引(细化后)
        integer,  allocatable :: mrl_new(:)         ! 三角形网格细化程度(细化后) 
        integer,  allocatable :: ngr_mrl_new(:, :)  ! 三角形网格的相邻三角形网格点的mrl(细化后)
        character(LEN = 256) :: lndname
        character(LEN = 5) :: nxpc, stepc
        logical :: iterA                  ! 当迭代B与迭代C同时一次性通过，迭代A通过
        logical :: iterB                  ! 从三角形网格进行判断
        logical :: iterC                  ! 从多边形网格进行判断
        ! logical :: isexist                ! 判断细化后是否存在重复w点
        ! read unstructure mesh
        print*, "start to read unstructure mesh data"
        ! 读取未细化初始网格数据
        write(nxpc, '(I4.4)') NXP
        write(stepc, '(I2.2)') step
        lndname = trim(file_dir) // 'gridfile/gridfile_NXP' // trim(nxpc) // '_'//trim(stepc)//'.nc4'
        ! lndname = "/stu01/zhangr23/makegrid/cases/case3_landmeshv4_20250107_nxp144_global/tmpfile/gridfile_NXP144_00_ori.nc4"
        print*,lndname
        ! 只有在refine期间，可能有一些图形不存在！！！！！！！！
        call Unstructured_Mesh_Read(lndname, sjx_points, lbx_points, mp, wp, ngrmw, ngrwm, n_ngrwm)
        print*, "The unstructured grid data reading have done "
        print*, ""
        print*, "In total, triangular mesh number: ", sjx_points, "polygon mesh number: ", lbx_points
        print*, ""
        !----------------------------------------------------------------
        ! 完成mp, wp, ngrmw, ngrwm, mp_new, wp_new, ngrmw_new, ngrwm_new 数据的初始化与更新
        !----------------------------------------------------------------
        print*, "limit -180 >> 180 90 >> -90 start"
        do i = 2, sjx_points, 1 ! mp is the center of sjx
            if (mp(i, 1) == -180.) mp(i, 1) = 180.
        end do
        do i = 2, lbx_points, 1 ! wp is the center of lbx
            if (wp(i, 1) == -180.) wp(i, 1) = 180.
        end do
        print*, "limit -180 >> 180 90 >> -90 finish"

        ! 扩大四倍，是因为就算全部都要细化，也就是四倍 可以修改为根据实际大小确定，但是难度很大
        ! 对于连接关系ngr**（ngrmm,ngrmw,ngrwm）而言，为1就是无连接的含义
        allocate(mp_new(sjx_points * 4, 2))    ; mp_new = 0.   ! The center point of the triangular grid updates the data (三角形网格中心点更新数据)
        allocate(wp_new(lbx_points * 4, 2))    ; wp_new = 0.   ! Update data at center point of polygon mesh (多边形网格中心点更新数据)
        allocate(ngrmw_new(3, sjx_points * 4)) ; ngrmw_new = 1 ! The adjacent wp points of mp update the index table (mp的相邻wp点更新索引表)，
        allocate(ngrwm_new(7, lbx_points * 4)) ; ngrwm_new = 1 ! wp's adjacent mp points update the index table (wp的相邻mp点更新索引表)
        mp_new(1:sjx_points, :) = mp ! 三角形中心的的经，纬度
        wp_new(1:lbx_points, :) = wp ! 多边形中心的经，纬度
        ngrmw_new(:, 1:sjx_points) = ngrmw(1:3, 1:sjx_points) ! 三角形顶点的的经，纬度
        ngrwm_new(:, 1:lbx_points) = ngrwm(1:7, 1:lbx_points) ! 多边形顶点的的经，纬度

        !----------------------------------------------------------------
        ! 完成ngrmm ,mrl_new, ngrmm_new, ngr_mrl_new 数据的初始化与更新
        ! 初始化，认为三角形都存在，赋值为1
        !----------------------------------------------------------------
        ! Triangle mesh refinement degree (三角形网格细化程度/方式) 1为不细化，2,4为细化, 0为三角形不存在
        allocate(mrl_new(sjx_points));                  mrl_new = 1     ! 反映三角形自身细化与否（0或者1），以及细化方法（2或者4）
        ! mp adjacent mp initial index table (m点相邻m点的初始索引表)
        allocate(ngrmm(3, sjx_points));             ngrmm = 1   ! 反映三角形的相邻三角形的邻域关系(0表示没有相邻三角形)
        allocate(ngrmm_new(3, sjx_points)); ngrmm_new = ngrmm       ! Updated ngrmm data (更新后数据)
        ! 在FHW代码中，lbx的第一次细化的时候ngrmm初始为 1 但是在二次细化以及后面初始为0
        ! mrl of adjacent triangular mesh points of a triangular mesh (三角形网格的相邻三角形网格的mrl)
        allocate(ngr_mrl_new(3, sjx_points));           ngr_mrl_new = 1 ! 反映三角形的相邻三角形的细化程度（0,1,2,4）

        ! 更新ngrmm，获取三角形的邻域编号信息，邻域的指向关系实际上是相互的，可否用于减小计算？？
        print*, "ngrmm renew start"
        do mi = 2, sjx_points, 1
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
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", sjx_points, spDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 3, DimaID))
        CALL CHECK(NF90_DEF_VAR(ncID, "ngrmm", NF90_INT, (/ DimaID, spDimID /), ncVarID(1)))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(1), ngrmm))
        CALL CHECK(NF90_CLOSE(ncID))
        print*, "ngrmm renew finish"

        ! ref_lbx 前七列分别表示组成这个多边形的三角形是否被细化 分别是0, 0.5, 1
        ! ref_lbx 第八列表示该多边形是否存在被细化的情况! 最后一个位置存放细化情况，分为0, 1, 4
        allocate(ref_lbx(lbx_points, 8));      ref_lbx = 0. 
        if (step >= 1) then
            ! 只有经过一轮细化之后才有可能产生七边形,在七边形的set_dis之内的三角形都不可以细化，反过来思考就好了
            print*, "根据相邻距离筛选七边形附近点"
            set_dis = 5
            print*, "before sum(ref_sjx) = ", sum(ref_sjx)
            ! 可以弄一个GetTriangleDisv1和GetTriangleDisv2去对比效果，速度也看看
            ! Old
            allocate(mp_dis(sjx_points, sjx_points)); mp_dis = 0
            CALL GetTriangleDis(sjx_points, ngrmw, ngrwm, n_ngrwm, ref_sjx, mp_dis)
            do i = 2, sjx_points, 1
                if (ref_sjx(i) == 0) cycle
                do j = 2, sjx_points, 1
                    if (mp_dis(i, j) <= set_dis) then
                        do k = 1, 3, 1
                            if (n_ngrwm(ngrmw(k, j)) == 7) then
                                ref_sjx(i) = 0
                                ! print*, "i =", i
                            end if
                        end do
                    end if
                end do
            end do
            deallocate(mp_dis) 
            ! CALL GetTriangleDis2(set_dis, sjx_points, lbx_points, ngrmw, ngrwm, n_ngrwm, num_inqbx, sjx_inqbx, mp_dis)
            ! New
            ! do i = 1, num_inqbx, 1
            !     do j = 2, sjx_points, 1
            !         if (mp_dis(i, j) <= set_dis) then
            !             do k = 1, 3, 1
            !                 if (n_ngrwm(ngrmw(k, j)) == 7) then
            !                     ref_sjx(i) = 0
            !                     ! print*, "i =", i
            !                 end if
            !             end do
            !         end if
            !     end do
            ! end do

            print*, "after sum(ref_sjx) = ", sum(ref_sjx)
            if (INT(sum(ref_sjx)) == 0) then
                exit_loop = .true.   
                return
            end if

        end if

        !--------------------------------------------------
        ! 1.2 Preliminary refinement (one into four) 【初步细化（一分为四）】初步细化就是阈值细化
        !--------------------------------------------------
        iter = 1                                 ! 本次细化中的迭代次数
        num_mp(iter) = sjx_points ! 后面涉及sjx_points与num_mp，直接用num_mp(1)代替sjx_points
        num_wp(iter) = lbx_points ! 后面涉及lbx_points与num_wwp，直接用num_wp(1)代替lbx_points
        ! 存放原始数据，就是直接对读入数据cp就好了
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc) // "_"// trim(stepc) //"_ori.nc4"
        print*, lndname
        call Unstructured_Mesh_Save(lndname, num_mp(iter), num_wp(iter), mp_new, wp_new, ngrmw_new, ngrwm_new)
        
        
        num_ref = INT(sum(ref_sjx))                  ! 需要细化的三角形个数
        print*, "Start to refine (开始初步细化)"
        print*, "iter =", iter, "num =", num_ref ! iter：迭代次数 num_ref：需要细化的三角形个数
        iter = iter + 1 ! 相对于FHW的代码，这是加一之后的结果
        ! 先把预估的位置占出来
        num_mp(iter) = num_mp(iter - 1) + 4 * num_ref        ! 记录每次迭代后三角形数，每细化（一分为四）一个三角形，增加4个小三角形的中心点，
        num_wp(iter) = num_wp(iter - 1) + 3 * num_ref        ! 记录每次迭代后多边形数，每细化（一分为四）一个三角形，增加三个中点作为三角形的顶点
        call OnedivideFour(iter, ngrmw, ngrmm, ref_lbx, ngrmw_new, ngrmm_new, mp_new, wp_new, mrl_new, ngr_mrl_new, num_mp, num_wp)
        print*, "Refining step 1 is complete (细化第一步完成)"
        !-------------------------------------------------
        ! 1.3 储存初始网格数据和初步细化后新网格数据
        !--------------------------------------------------
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc) // "_"// trim(stepc) //"_1.nc4"
        print*, lndname
        call Unstructured_Mesh_Save(lndname, num_mp(iter), num_wp(iter), mp_new, wp_new, ngrmw_new, ngrwm_new)
        
        !--------------------------------------------------
        ! 2.1 进行迭代（防止细化交汇带出现冲突，一分为四）! 也是采用一分四算法 iterB 和 iterC 用三角形与多边形的角度去防止细化带的出现
        !--------------------------------------------------
        print*, "iterA start"
        iterA = .false.
        do while(iterA .eqv. .false.) ! 当iterA为true时，该步骤完成

            iterA = .true.    ! 判断迭代B和迭代C是否都已满足条件
            iterB = .false.   ! 从三角形网格进行判断
            iterC = .false.   ! 从多边形网格进行判断

            print*, "iterB start" ! 从三角形网格进行判断
            do while(iterB .eqv. .false.)

                call iterB_judge(sjx_points, mrl_new, ngr_mrl_new)
                num_ref = INT(sum(ref_sjx)) ! 获取需要细化的三角形个数
                if (num_ref == 0) then
                    iterB = .true.
                else
                    print*, "iter =", iter, "num =", num_ref
                    iterA = .false.
                    iter = iter + 1
                    num_mp(iter) = num_mp(iter - 1) + 4 * num_ref
                    num_wp(iter) = num_wp(iter - 1) + 3 * num_ref
                    CALL OnedivideFour(iter, ngrmw, ngrmm, ref_lbx, ngrmw_new, ngrmm_new, mp_new, wp_new, mrl_new, ngr_mrl_new, num_mp, num_wp)
                end if
            end do ! iterB
            print*, "iterB end"

            print*, "iterC start"! 从多边形形网格进行判断
            do while(iterC .eqv. .false.)

                call iterC_judge(lbx_points, ngrmm, ngrwm, n_ngrwm, mrl_new, ngr_mrl_new, ref_lbx)
                num_ref = INT(sum(ref_sjx)) ! 获取需要细化的三角形个数
                if (num_ref == 0) then
                    iterC = .true.
                else
                    print*, "iter =", iter, "num =", num_ref
                    iterA = .false.
                    iter = iter + 1
                    num_mp(iter) = num_mp(iter - 1) + 4 * num_ref
                    num_wp(iter) = num_wp(iter - 1) + 3 * num_ref
                    CALL OnedivideFour(iter, ngrmw, ngrmm, ref_lbx, ngrmw_new, ngrmm_new, mp_new, wp_new, mrl_new, ngr_mrl_new, num_mp, num_wp)
                end if
            end do ! iterC
            print*, "iterC end"

        end do ! iterA
        print*, "iterA end"
        print*, "细化第二步完成"
        ! stop "stop for test"
        !--------------------------------------------------
        ! 2.2 储存第二步细化后新网格数据
        !--------------------------------------------------
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc) // "_"// trim(stepc) //"_2.nc4"
        print*, lndname
        call Unstructured_Mesh_Save(lndname, num_mp(iter), num_wp(iter), mp_new, wp_new, ngrmw_new, ngrwm_new)
        
        !--------------------------------------------------
        ! 3.1 寻找并细化弱凹点
        !--------------------------------------------------
        ! 寻找弱凹点
        ref_sjx = 0
        do i = 2, sjx_points, 1
            if ((mrl_new(i) == 1) .and. (sum(ngr_mrl_new(1:3, i)) == 6)) then
                do j = 1, 3, 1
                    mm = ngrmm(j, i) 
                    if ((mrl_new(mm) == 1) .and. (sum(ngr_mrl_new(1:3, mm)) == 6)) then
                        do k = 1, 3, 1
                            if (ngr_mrl_new(k, i) == 4)  m = k
                            if (ngr_mrl_new(k, mm) == 4) n = k
                        end do
                        do k = 1, 3, 1
                            do l = 1, 3, 1
                                ! 这是共同顶点的判断,只能用ngrmw，因为ngrmw_new的信息已经被更新了
                                if (ngrmw(k, ngrmm(m, i)) == ngrmw(l, ngrmm(n, mm))) then
                                    ref_sjx(i) = 1
                                    ! print*, "i = ", i
                                end if
                            end do
                        end do
                    end if
                end do
            end if
        end do
        num_ref = INT(sum(ref_sjx))
        if (num_ref == 0) then
            print*,"no 弱凹点"
        else
            print*, "开始细化弱凹点"
            print*, "iter =", iter, "num =", num_ref
            iter = iter + 1
            num_mp(iter) = num_mp(iter - 1) + 2 * num_ref
            num_wp(iter) = num_wp(iter - 1) + num_ref         
            call TwodivideFour(iter, ngrmw, ngrmm, ngrwm, num_mp, num_wp, mp_new, wp_new, mrl_new, ngrmw_new, ngrmm_new, ngr_mrl_new)
        end if
        print*, "细化第三步完成"

        !--------------------------------------------------
        ! 3.3 储存第三步细化后新网格数据
        !--------------------------------------------------
        ! 这部分完全可以放在一个子例行程序里面去，然后放在最后面
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc) // "_"// trim(stepc) //"_3.nc4"
        print*, lndname
        call Unstructured_Mesh_Save(lndname, num_mp(iter), num_wp(iter), mp_new, wp_new, ngrmw_new, ngrwm_new)
        
        !--------------------------------------------------
        ! 4.1 记录相邻三角形中只有一个三角形但被细化的的情况
        !--------------------------------------------------        
        ref_sjx = 0
        do i = 2, sjx_points, 1
            if (mrl_new(i) == 1) then
                if (sum(ngr_mrl_new(:, i)) == 6) ref_sjx(i) = 1
            end if
        end do
        num_ref = INT(sum(ref_sjx))
        if (num_ref == 0) then
            print*,"NO 相邻三角形经过初步细化数为1的三角形"
        else
            print*, "开始细化相邻三角形经过初步细化数为1的三角形"
            print*, "iter =", iter, "num =", num_ref
            iter = iter + 1
            num_mp(iter) = num_mp(iter - 1) + 2 * num_ref
            num_wp(iter) = num_wp(iter - 1) + num_ref
            call OnedivideTwo(iter, ngrmw, ngrmm, ngrwm, num_mp, num_wp, mp_new, wp_new, mrl_new, ngrmw_new, ngrmm_new, ngr_mrl_new)
        end if
        print*, "细化第四步完成"
 
        !--------------------------------------------------
        ! 4.3 储存第四步细化后新网格数据
        !--------------------------------------------------
        ! 这部分完全可以放在一个子例行程序里面去，然后放在最后面，就是存入的路径可能有点不同而已
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc) // "_"// trim(stepc) //"_4.nc4"
        print*, lndname
        call Unstructured_Mesh_Save(lndname, num_mp(iter), num_wp(iter), mp_new, wp_new, ngrmw_new, ngrwm_new)
        
        print*, "细化后共有", num_wp(iter), "个多边形网格"
        print*, "细化后共有", num_mp(iter), "个三角形网格" 

        !--------------------------------------------------
        ! 5.5 存储网格数据
        !--------------------------------------------------
        CALL NGR_RENEW(iter, num_mp, num_wp, mp_new, wp_new, ngrmw_new, num_sjx, num_dbx, mp_f, wp_f, ngrmw_f, ngrwm_f, n_ngrwm_f)
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc) // "_"// trim(stepc) //"_5.nc4"
        print*, lndname
        call Unstructured_Mesh_Save(lndname, num_sjx, num_dbx, mp_f, wp_f, ngrmw_f, ngrwm_f)
        
        !--------------------------------------------------
        ! 6.1 弹性调整
        !--------------------------------------------------
        call SpringAjustment(sjx_points, num_sjx, num_dbx, num_ref, ngrmw_f, ngrwm_f, n_ngrwm_f, mp_f, wp_f)
        
        !--------------------------------------------------
        ! 6.7 存储最终网格数据
        !--------------------------------------------------
        ! 规定180°经线圈上的经度为180° (限制经度范围)
        do i = 2, num_sjx, 1 ! mp is the center of sjx
            if (mp_f(i, 1) == -180.) mp_f(i, 1) = 180.
        end do
        do i = 2, num_dbx, 1 ! wp is the center of lbx
            if (wp_f(i, 1) == -180.) wp_f(i, 1) = 180.
        end do
        write(stepc, '(I2.2)') step+1
        CALL GetSort(n_ngrwm_f, num_dbx, mp_f, ngrwm_f)
        lndname = trim(file_dir) // 'gridfile/gridfile_NXP' // trim(nxpc) // '_'//trim(stepc)//'.nc4'
        print*, lndname
        CALL Unstructured_Mesh_Save(lndname, num_sjx, num_dbx, mp_f, wp_f, ngrmw_f, ngrwm_f)
        deallocate(ref_sjx)

        ! deallocate
        deallocate(mp, wp, mp_new, wp_new, mp_f, wp_f)
        deallocate(ref_lbx)
        deallocate(n_ngrwm)
        deallocate(ngrmw, ngrwm, ngrmw_new, ngrwm_new, ngrmw_f, ngrwm_f, n_ngrwm_f)
        deallocate(ngrmm, ngrmm_new)
        deallocate(mrl_new)
        deallocate(ngr_mrl_new)

    END SUBROUTINE refine_loop

    SUBROUTINE iterB_judge(sjx_points, mrl_new, ngr_mrl_new)

        implicit none

        integer, intent(in)    :: sjx_points
        integer, allocatable, intent(in) :: mrl_new(:), ngr_mrl_new(:,:)
        integer :: i
        
        ref_sjx = 0 ! 记录三角形是否需要被细化，初始化
        ! 当一个未细化三角形有两个或两个以上细化相邻三角形时，细化该三角形
        do i = 2, sjx_points, 1
            if (mrl_new(i) == 1) then ! 自身未被细化
                if (sum(ngr_mrl_new(:, i)) > 7.) ref_sjx(i) = 1 ! 三个细化就是12 ，两个细化就是8
            end if
        end do

    END SUBROUTINE iterB_judge
    
    SUBROUTINE iterC_judge(lbx_points, ngrmm, ngrwm, n_ngrwm, mrl_new, ngr_mrl_new, ref_lbx)
        
        implicit none

        integer,  intent(in)    :: lbx_points
        integer,  allocatable, intent(in) :: ngrmm(:,:), ngrwm(:,:), n_ngrwm(:), mrl_new(:), ngr_mrl_new(:,:)
        real(r8), allocatable, intent(inout) :: ref_lbx(:,:)
        integer :: num_edges, i, j, m1, m2, num
 
        ref_sjx = 0

            do i = 2, lbx_points, 1
               
                num_edges = n_ngrwm(i) ! 用于统计多边形网格边数
                do j = 1, num_edges, 1
                    m1 = ngrwm(j, i) ! 获取该相邻三角形的编号
                    ! 1+1+4，即相邻三角形中一个细化两个未细化 标记这个三角形，是因为未来需要被细化
                    ! 但是没有要求自身是否已经被细化，所以这个判断很奇怪
                    if (sum(ngr_mrl_new(:, m1)) == 6) ref_lbx(i, j) = 1
                end do

                ! 只考虑五边形与六边形，因为有可能在第二/三次细化的时候起手就有七边形
                if (ref_lbx(i, 8) == 0) then ! 表示i号w点的相邻三角形均未被细化       
                    if ((num_edges == 5) .or. (num_edges == 6)) then ! 对于edges 等于 5 或者 6 的情况 处理都是一样的   
                        do j = 1, num_edges, 1
                            m1 = j
                            m2 = mod(j, num_edges) + 1
                            ! 可能存在来连续三个三角形都是有外接的细化三角形会如何？？？？？？？？？？
                            if ((sum(ngr_mrl_new(:, ngrwm(m1, i))) == 6) .and. &
                                (sum(ngr_mrl_new(:, ngrwm(m2, i))) == 6)) then ! 针对两个相邻的三角形，因此它们组成弱凹
                                ! 这是重新赋值的过程，因为如果出现相邻的情况可以减少一条边的增量
                                ref_lbx(i, m1) = 0.5 
                                ref_lbx(i, m2) = 0.5
                            end if! 两个0.5合计1，表示要增加1条边
                        end do

                        if (sum(ref_lbx(i, 1:num_edges)) + num_edges > 7.) then 
                            do j = 1, num_edges, 1 !New 应该不会循环到7
                                if ((ref_lbx(i, j) /= 0) .and. (mrl_new(ngrwm(j, i)) == 1)) then
                                    ref_sjx(ngrwm(j, i)) = 1 ! 说明该三角形要细化（一分四那种）
                                    ! print*, "ngrwm(j, i) = ", ngrwm(j, i), "ref_lbx==0"
                                end if
                            end do
                        end if
                    end if   ! if((num_edges == 5) .or. (num_edges == 6))then 

                else if (ref_lbx(i, 8) /= 0) then! 当多边形中心点相邻的三角形存在细化，即ref_lbx(i, 8) /= 0 可能存在弱凹点

                    if (num_edges == 5) then ! 五边形没有弱凹点！！！
                        ! 只可能是连续两个或者三个三角形被细化
                        if (sum( mrl_new(ngrwm(1:num_edges, i)) ) > 10 ) then
                            do j = 1, num_edges, 1
                                if (mrl_new(ngrwm(j, i)) == 1) ref_sjx(ngrwm(j, i)) = 1
                                ! if (mrl_new(ngrwm(j, i)) == 1) print*, "ngrwm(j, i) = ", ngrwm(j, i), "ref_lbx/=0"
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
                                    if ((mrl_new(ngrwm(j + 1, i)) == 1) .and. &
                                        (mrl_new(ngrwm(j + 2, i)) == 1)) then
                                        ref_sjx(ngrwm(j+1:j+2, i)) = 1
                                        ! print*, "ngrwm(j+1, i) = ", ngrwm(j+1, i), "ruo"
                                        ! print*, "ngrwm(j+2, i) = ", ngrwm(j+2, i), "ruo"
                                    end if
                                end if ! 细化一侧的三角形，将另一侧视作弱凹
                            end do

                        else if (sum(mrl_new(ngrwm(1:num_edges, i))) == 9) then ! 当i点相邻三角形中只有一个被细化
                            num = 0
                            do j = 1, num_edges, 1
                                m1 = ngrwm(j, i) ! 循环获取相邻三角形的编号
                                if (mrl_new(m1) == 1) then ! 如果该三角形本身还没有细化
                                    if (sum(mrl_new(ngrmm(1:3, m1))) == 6) then ! 如果该三角形相邻三角形只有一个被细化
                                        num = num + 1
                                    end if
                                end if
                            end do
                            if (num > 3) then! 因此6+2大于7时，细化所有相邻三角形，避免出现边数大于7的情况
                                do j = 1, num_edges, 1
                                    m1 = ngrwm(j, i)
                                    if (mrl_new(m1) == 1) ref_sjx(m1) = 1
                                    ! if (mrl_new(m1) == 1) print*, "m1 = ", m1
                                end do
                            end if
                        end if

                    end if ! num_edges == 6 或者 5 的问题

                end if ! ref_lbx(i, 8) 与 0 的关系

            end do ! i = 2, lbx_points, 1 循环

    END SUBROUTINE iterC_judge
 
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

    SUBROUTINE CheckLon(x)

        implicit none

        real(r8), intent(out) :: x

        if (x > 180.) then
            x = x - 360.
        else if (x < -180.) then
            x = x + 360.
        end if

    END SUBROUTINE CheckLon


    SUBROUTINE GetTriangleLength(sjx, length)

        implicit none
    
        real(r8), dimension(2), intent(in) :: sjx(3, 2)
        real(r8) :: a(2), b(2), c(2), v(3)
        real(r8), intent(inout) :: length(3)
        !  需要考虑跨越180经线的情况，而且是返回三个边长
        length = 0.
        a = sjx(1, :) * pio180 ! 顶点对应的边
        b = sjx(2, :) * pio180
        c = sjx(3, :) * pio180
        ! 半正矢公式, 以弧度制度量!!
        v(1) = sqrt( sin(c(2)/2-b(2)/2)**2 + cos(b(2))*cos(c(2))*sin(c(1)/2-b(1)/2)**2 )
        v(2) = sqrt( sin(c(2)/2-a(2)/2)**2 + cos(a(2))*cos(c(2))*sin(c(1)/2-a(1)/2)**2 )
        v(3) = sqrt( sin(a(2)/2-b(2)/2)**2 + cos(b(2))*cos(a(2))*sin(a(1)/2-b(1)/2)**2 )
        length(:) = erad * 2 * asin(v)

    END SUBROUTINE GetTriangleLength


    SUBROUTINE GetAngle(num_edges, dbx, angle_dbx) 
        ! 适用于三角形与五六七边形的球面内角计算形式  
        implicit none
    
        integer, intent(in) :: num_edges
        real(r8), intent(in) :: dbx(num_edges, 2) 
        integer :: i
        real(r8), allocatable :: dbx_cycle(:, :)
        real(r8) :: length(3), sjx(3, 2)
        real(r8), intent(out) :: angle_dbx(num_edges) ! 最后传入的边数一致
        ! 形成闭合, 设置为通用的闭合形式
        sjx = 0.
        allocate(dbx_cycle(num_edges + 2, 2))
        dbx_cycle(1, :)             = dbx(num_edges, :)
        dbx_cycle(2:num_edges+1, :) = dbx(1:num_edges, :)
        dbx_cycle(num_edges+2, :)   = dbx(1, :)
        do i = 1, num_edges, 1
            sjx = dbx_cycle(i:i + 2, :)
            ! 先计算长度，再计算角度
            CALL GetTriangleLength(sjx, length)
            ! 计算球面角，并实现弧度转角度，给出三点，但是角度到底是哪一个呢？？
            angle_dbx(i) = (length(1) * length(1) + length(3) * length(3) - length(2) * length(2)) &
                            / (2. * length(1) * length(3))
            angle_dbx(i) = acos(angle_dbx(i)) * piu180
        end do
        deallocate(dbx_cycle)
    
    END SUBROUTINE GetAngle

    SUBROUTINE SpringAjustment(sjx_points, num_sjx, num_dbx, num_ref, ngrmw_f, ngrwm_f, n_ngrwm_f, mp_f, wp_f)

        IMPLICIT NONE
        ! Extr_*** 指网格的最小角度与最大角度
        ! Eavg_*** 指网格的最小角度与最大角度平均值
        ! Savg_*** 指网格角度与正多边形角度的标准差
        ! less30   指三角形网格中小于30度角的数量占比
        integer,  intent(in) :: sjx_points, num_sjx, num_dbx, num_ref
        integer,  allocatable, intent(in) :: ngrmw_f(:, :), ngrwm_f(:, :), n_ngrwm_f(:)
        character(LEN = 256) :: lndname
        character(LEN = 5) :: stepc, nxpc
        integer :: w1, w2, w3, icl
        integer :: i, j, k, sa_iter, hhh(4) 
        integer :: lpDimID, twDimID, ncVarID(13), ncid
        integer :: num_wbx, num_lbx, num_qbx
        real(r8) :: rx, ry, fra, sjx(3, 2)
        real(r8) :: Extr_sjx_temp(2), Eavg_sjx_temp(2), Savg_sjx_temp, less30_temp
        real(r8) :: Extr_wbx_temp(2), Eavg_wbx_temp(2), Savg_wbx_temp
        real(r8) :: Extr_lbx_temp(2), Eavg_lbx_temp(2), Savg_lbx_temp
        real(r8) :: Extr_qbx_temp(2), Eavg_qbx_temp(2), Savg_qbx_temp
        logical :: End_SpringAjustment 
        real(r8), allocatable :: MoveDis(:, :)  ! 记录w点在x、y方向上的调整距离
        real(r8), allocatable :: length_sjx(:, :, :), length_sjx_temp(:, :)
        real(r8), allocatable :: angle_sjx(:, :, :), Extr_sjx(:, :), Eavg_sjx(:, :), Savg_sjx(:), less30(:)
        real(r8), allocatable :: angle_wbx(:, :, :), Extr_wbx(:, :), Eavg_wbx(:, :), Savg_wbx(:)
        real(r8), allocatable :: angle_lbx(:, :, :), Extr_lbx(:, :), Eavg_lbx(:, :), Savg_lbx(:)
        real(r8), allocatable :: angle_qbx(:, :, :), Extr_qbx(:, :), Eavg_qbx(:, :), Savg_qbx(:)
        real(r8), allocatable :: angle_sjx_temp(:, :), angle_wbx_temp(:, :), angle_lbx_temp(:, :), angle_qbx_temp(:, :) 
        real(r8), allocatable, intent(inout) :: wp_f(:, :)
        real(r8), allocatable, intent(inout) :: mp_f(:, :)

        write(stepc, '(I2.2)') step
        write(nxpc,  '(I4.4)') nxp

        ! obtain num_wbx, num_lbx, num_qbx
        num_wbx = 0; num_lbx = 0; num_qbx = 0; 
        do i = 2, num_dbx, 1
            if (n_ngrwm_f(i) == 5) then
                num_wbx = num_wbx + 1
            else if (n_ngrwm_f(i) == 6) then
                num_lbx = num_lbx + 1
            else if (n_ngrwm_f(i) == 7) then
                num_qbx = num_qbx + 1
            end if
        end do
        print*, "num_wbx = ", num_wbx
        print*, "num_lbx = ", num_lbx
        print*, "num_qbx = ", num_qbx 

        ! for sjx
        allocate(length_sjx(0:max_sa_iter, num_sjx, 3)); length_sjx(:, :, :) = 0.
        allocate(angle_sjx(0:max_sa_iter, num_sjx, 3));  angle_sjx(:, :, :)  = 0.
        allocate(Eavg_sjx(0:max_sa_iter, 2)); Extr_sjx = 0.
        allocate(Extr_sjx(0:max_sa_iter, 2)); Eavg_sjx = 0.
        allocate(Savg_sjx(0:max_sa_iter));    Savg_sjx = 0.
        allocate(less30(0:max_sa_iter)); less30 = 0.
        allocate(length_sjx_temp(num_sjx, 3));          length_sjx_temp(:, :) = 0.
        allocate(angle_sjx_temp(num_sjx, 3));           angle_sjx_temp(:, :)  = 0.
        ! for wbx
        allocate(angle_wbx(0:max_sa_iter, num_wbx, 5));  angle_wbx(:, :, :)  = 0.
        allocate(Eavg_wbx(0:max_sa_iter, 2)); Extr_wbx = 0.
        allocate(Extr_wbx(0:max_sa_iter, 2)); Eavg_wbx = 0.
        allocate(Savg_wbx(0:max_sa_iter));    Savg_wbx = 0.
        allocate(angle_wbx_temp(num_wbx, 5));           angle_wbx_temp(:, :)  = 0.
        ! for lbx
        allocate(angle_lbx(0:max_sa_iter, num_lbx, 6));  angle_lbx(:, :, :)  = 0.
        allocate(Eavg_lbx(0:max_sa_iter, 2)); Extr_lbx = 0.
        allocate(Extr_lbx(0:max_sa_iter, 2)); Eavg_lbx = 0.
        allocate(Savg_lbx(0:max_sa_iter));    Savg_lbx = 0.
        allocate(angle_lbx_temp(num_lbx, 6));           angle_lbx_temp(:, :)  = 0.
        ! for qbx
        allocate(angle_qbx(0:max_sa_iter, num_qbx, 7));  angle_qbx(:, :, :)  = 0.
        allocate(Eavg_qbx(0:max_sa_iter, 2)); Extr_qbx = 0.
        allocate(Extr_qbx(0:max_sa_iter, 2)); Eavg_qbx = 0.
        allocate(Savg_qbx(0:max_sa_iter));    Savg_qbx = 0.
        allocate(angle_qbx_temp(num_qbx, 7));           angle_qbx_temp(:, :)  = 0.
        Extr_sjx_temp = 0.; Eavg_sjx_temp = 0.; Savg_sjx_temp = 0.; less30_temp = 0.
        Extr_wbx_temp = 0.; Eavg_wbx_temp = 0.; Savg_wbx_temp = 0.;
        Extr_lbx_temp = 0.; Eavg_lbx_temp = 0.; Savg_lbx_temp = 0.;
        Extr_qbx_temp = 0.; Eavg_qbx_temp = 0.; Savg_qbx_temp = 0.;

        ! 弹性调整&计算每次调整后的三角形网格质量
        print*, "开始弹性调整"
        allocate(MoveDis(num_dbx, 2)) ! 记录w点在x、y方向上的调整距离
        sa_iter = 0
        k = -1 ! 初始化
        End_SpringAjustment = .false.
        hhh = [2,3,1,2] ! 用于简化计算
        
        do while(End_SpringAjustment .eqv. .false.)! 问题在于每次都要对全局进行调整，而且调整的方法是固定的
            MoveDis = 0. ! 调整距离初始化为0
            End_SpringAjustment = .true.
            ! 每次调整后，再计算三角形，五六七边形的网格质量
            Call TriMeshQuality(    num_sjx, wp_f, ngrmw_f, length_sjx_temp, angle_sjx_temp, Extr_sjx_temp, Eavg_sjx_temp, Savg_sjx_temp, less30_temp)! 三角形网格质量
            Call PolyMeshQuality(5, num_dbx, mp_f, ngrwm_f, n_ngrwm_f,       angle_wbx_temp, Extr_wbx_temp, Eavg_wbx_temp, Savg_wbx_temp)! 五边形网格质量
            Call PolyMeshQuality(6, num_dbx, mp_f, ngrwm_f, n_ngrwm_f,       angle_lbx_temp, Extr_lbx_temp, Eavg_lbx_temp, Savg_lbx_temp)! 六边形网格质量
            Call PolyMeshQuality(7, num_dbx, mp_f, ngrwm_f, n_ngrwm_f,       angle_qbx_temp, Extr_qbx_temp, Eavg_qbx_temp, Savg_qbx_temp)! 七边形网格质量
            
            ! tri
            length_sjx(sa_iter, :, :) = length_sjx_temp
            angle_sjx(sa_iter, :, :)  = angle_sjx_temp
            Extr_sjx(sa_iter, :)  = Extr_sjx_temp
            Eavg_sjx(sa_iter, :)  = Eavg_sjx_temp
            Savg_sjx(sa_iter)     = Savg_sjx_temp
            less30(sa_iter)       = less30_temp
            ! wbx
            angle_wbx(sa_iter, :, :)  = angle_wbx_temp
            Extr_wbx(sa_iter, :)  = Extr_wbx_temp
            Eavg_wbx(sa_iter, :)  = Eavg_wbx_temp
            Savg_wbx(sa_iter)     = Savg_wbx_temp
            ! lbx
            angle_lbx(sa_iter, :, :)  = angle_lbx_temp
            Extr_lbx(sa_iter, :)  = Extr_lbx_temp
            Eavg_lbx(sa_iter, :)  = Eavg_lbx_temp
            Savg_lbx(sa_iter)     = Savg_lbx_temp
            ! qbx
            angle_qbx(sa_iter, :, :)  = angle_qbx_temp
            Extr_qbx(sa_iter, :)  = Extr_qbx_temp
            Eavg_qbx(sa_iter, :)  = Eavg_qbx_temp
            Savg_qbx(sa_iter)     = Savg_qbx_temp

            print*, "第", sa_iter, "次弹性调整完成，调整角度个数为", k
            if (mod(sa_iter, 10) == 0) then
                print*, "三角形边长：", minval(length_sjx_temp(2:num_sjx, :)), maxval(length_sjx_temp(2:num_sjx, :))
                print*, "三角形角度：", minval(angle_sjx_temp(2:num_sjx, :)), maxval(angle_sjx_temp(2:num_sjx, :))
                print*, "五边形角度：", minval(angle_wbx_temp), maxval(angle_wbx_temp)
                print*, "六边形角度：", minval(angle_lbx_temp), maxval(angle_lbx_temp)
                print*, "七边形角度：", minval(angle_qbx_temp), maxval(angle_qbx_temp)
                print*, ""
            end if

            k = 0 ! 记录钝角三角形个数（需要调整的角度的个数）
            do i = sjx_points - num_ref + 1, num_sjx, 1! 为啥是从这个范围开始呀？？？
                do j = 1, 3, 1
                    if(angle_sjx(sa_iter, i, j) > 90.)then ! 出现钝角三角形
                        k = k + 1
                        fra = ( 1 - ( 72. / angle_sjx(sa_iter, i, j) ) ) / 100.
                        ! 获取三角形另外两个顶点的坐标，建议是逆时针获取
                        w1 = ngrmw_f(hhh(j),   i)
                        w2 = ngrmw_f(hhh(j+1), i)
                        rx = wp_f(w2, 1) - wp_f(w1, 1)
                        ry = wp_f(w2, 2) - wp_f(w1, 2)
                        
                        call CheckLon(rx)! 这是一种距离调整方式，但是只能记录上一步的调整结果
                        
                        MoveDis(w1, 1) = MoveDis(w1, 1) + rx * fra / 2.
                        MoveDis(w1, 2) = MoveDis(w1, 2) + ry * fra / 2.
                        MoveDis(w2, 1) = MoveDis(w2, 1) - rx * fra / 2.
                        MoveDis(w2, 2) = MoveDis(w2, 2) - ry * fra / 2.
                        exit ! 因为钝角只会出现一次
                    end if
                end do
            end do

            ! 全部弹性调整后，再调整三角形顶点的经纬度
            do i = 2, num_dbx, 1
                wp_f(i, 1:2) = wp_f(i, 1:2) + MoveDis(i, 1:2)
                call CheckLon(wp_f(i, 1))
            end do

            ! 调整所有m点至三角形网格重心     
            do i = 2, num_sjx, 1
                sjx = wp_f(ngrmw_f(:, i), :)
                if (maxval(sjx(:, 1)) - minval(sjx(:, 1)) > 180.) then
                    call CheckCrossing(3, sjx)
                    mp_f(i, 1:2) = (sjx(1, 1:2) + sjx(2, 1:2) + sjx(3, 1:2)) / 3.
                    if (mp_f(i, 1) < 0.) then
                        mp_f(i, 1) = mp_f(i, 1) + 180.
                    else
                        mp_f(i, 1) = mp_f(i, 1) - 180.
                    end if
                else
                    mp_f(i, 1:2) = (sjx(1, 1:2) + sjx(2, 1:2) + sjx(3, 1:2)) / 3.
                end if
            end do

            ! 不出现钝角或者达到最大调整次数的时候就退出while循环！
            if ((k /= 0) .and. (sa_iter < max_sa_iter)) End_SpringAjustment = .false.
            sa_iter = sa_iter + 1
        end do ! while(End_SpringAjustment .eqv. .false.)

        deallocate(MoveDis, length_sjx, length_sjx_temp)
        deallocate(angle_sjx, angle_wbx, angle_lbx, angle_qbx)
        deallocate(angle_sjx_temp, angle_wbx_temp, angle_lbx_temp, angle_qbx_temp)

        lndname = trim(file_dir) // "result/quality_NXP" // trim(nxpc) // '_' // trim(stepc) // "_lbx.nc4"
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "num_iter", max_sa_iter + 1, lpDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 2, twDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "less30",   NF90_FLOAT, (/ lpDimID /), ncVarID(1)))
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

        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1),  less30(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(2),  Eavg_sjx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(3),  Savg_sjx(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(4),  Extr_sjx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5),  Eavg_wbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6),  Eavg_lbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(7),  Eavg_qbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(8),  Savg_wbx(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(9),  Savg_lbx(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(10), Savg_qbx(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(11), Extr_wbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(12), Extr_lbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(13), Extr_qbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_CLOSE(ncID))
        deallocate(Extr_sjx, Eavg_sjx, Savg_sjx, less30)
        deallocate(Extr_wbx, Eavg_wbx, Savg_wbx)
        deallocate(Extr_lbx, Eavg_lbx, Savg_lbx)
        deallocate(Extr_qbx, Eavg_qbx, Savg_qbx)
        print*, "弹性调整结束"

    END SUBROUTINE SpringAjustment

    SUBROUTINE crossline_check(iter, mp_new, wp_new, num_mp, num_wp)

        IMPLICIT NONE
        integer :: i
        integer,  intent(in) :: iter
        real(r8), dimension(:, :), allocatable, intent(inout) :: mp_new, wp_new
        integer,  dimension(:), intent(inout) :: num_mp, num_wp
        do i = num_mp(iter - 1) + 1, num_mp(iter), 1 ! mp_new is the center of sjx
            if (mp_new(i, 1) == -180.) mp_new(i, 1) = 180.
        end do
        do i = num_wp(iter - 1) + 1, num_wp(iter), 1 ! wp_new is the center of lbx
            if (wp_new(i, 1) == -180.) wp_new(i, 1) = 180.
        end do

    END SUBROUTINE crossline_check
    
    SUBROUTINE OnedivideFour(iter, ngrmw, ngrmm, ref_lbx, ngrmw_new, ngrmm_new, mp_new, wp_new, mrl_new, ngr_mrl_new, num_mp, num_wp)

        IMPLICIT NONE
        integer :: i, j, k, refed_iter
        integer :: icl, m0, w0      ! 三角形和多边形中心点序号起始索引
        real(r8) :: sjx(3, 2), newsjx(4, 2), newdbx(3, 2)
        integer,  intent(in) :: iter
        integer,  dimension(:, :), allocatable, intent(in) :: ngrmw, ngrmm
        real(r8), dimension(:, :), allocatable, intent(inout) :: mp_new, wp_new, ref_lbx
        integer,  dimension(:), intent(inout) :: num_mp, num_wp
        integer,  dimension(:),    allocatable, intent(inout) :: mrl_new
        integer,  dimension(:, :), allocatable, intent(inout) :: ngrmw_new, ngrmm_new, ngr_mrl_new

        sjx = 0.; newdbx = 0.; newsjx = 0.
        refed_iter = 0
        ! 需要建立refed_iter 与第几个三角形的映射关系
        do i = 2, num_mp(1), 1
            if ((ref_sjx(i) == 0) .or. (mrl_new(i) /= 1)) cycle ! 若三角形需要细化而且还没别细化
            ! print*, "i = ", i, "OnedivideFour"
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

            m0 = num_mp(iter - 1) + refed_iter * 4 ! 新三角形中心点编号基准
            w0 = num_wp(iter - 1) + refed_iter * 3 ! 新多边形中心点编号基准
            mp_new(m0 + 1 : m0 + 4, 1:2) = newsjx(:, 1:2)! 新三角形中心点经纬度(四个)
            wp_new(w0 + 1 : w0 + 3, 1:2) = newdbx(:, 1:2)! 新多边形中心点经纬度(三个)

            ! 连接关系调整 更新ngrmm_new,ngrmw_new,mrl_new,ngr_mrl_new,ref_lbx
            ref_lbx(ngrmw(1:3, i), 8) = 1

            ! 增加第一，二，三个新三角形的剩下两个顶点的编号信息
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

            ! 更新三角形的自身细化状态mrl_new(分为两个部分，一个是自身三角形，一个是细化生成的三角形)
            mrl_new(i) = 4 ! 原三角形网格被平均分为四份
            ! mrl_new(m0 + 1 : m0 + 4) = 4! 新三角形是由一分四算法得到的

            ! 更新三角形的领域细化状态ngr_mrl_new
            ! 要求ngrmm_new还没有被更新！！！！！！！！！！！
            do j = 1, 3, 1! 因为自身被细化了，所以需要说明自身的邻域三角形，指向自身的时候为4
                do k = 1, 3, 1
                    if (ngrmm(k, ngrmm(j, i)) == i) ngr_mrl_new(k, ngrmm(j, i)) = 4
                end do
            end do

            ! 更新原三角形与新三角形的ngrmm_new,ngrmw_new
            ! 注意前后的逻辑联系，特别是要将三角形进行取消的操作，例如ngrmm_new(k, i) = 1
            do k = 1, 3, 1
                ngrmm_new(k, i) = 1 ! 1. 去掉原三角形的邻域编号信息
                ngrmw_new(k, i) = 1 ! 1. 去掉原三角形的顶点编号信息
                ngrmw_new(1, m0 + k) = ngrmw(k, i)! 2. 增加新三角形的顶点编号信息
            end do

            refed_iter = refed_iter + 1
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !print*, "sjx = ", sjx
            !print*, "newsjx = ", newsjx
            !print*, "newdbx = ", newdbx
            !print*, "m0 = ", m0
            !print*, "w0 = ", w0
            !print*, "mp_new(m1:m4, 1:2) = ", mp_new(m0+1:m0+4, 1:2)
            !print*, "wp_new(w1:w3, 1:2) = ", wp_new(w0+1:w0+3, 1:2)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end do
        call crossline_check(iter, mp_new, wp_new, num_mp, num_wp)
    END SUBROUTINE OnedivideFour
    
    SUBROUTINE TwodivideFour(iter, ngrmw, ngrmm, ngrwm, num_mp, num_wp, mp_new, wp_new, mrl_new, ngrmw_new, ngrmm_new, ngr_mrl_new)            
        ! 只能用ngrmw, ngrmm, ngrwm ，因为*_new的信息已经被更新了
        ! 二分四算法（针对弱凹点）
        ! 弱凹点满足：
        ! 1.在三角形自己与相邻三个三角形只有一个被细化的三角形
        ! 2.组成弱凹点的两个三角形以及全部相邻三角形有且只有两个被细化的三角形
        ! 组成弱凹点的两个三角形 : i and j
        ! 有且只有两个被细化的三角形 : m and n
        ! i adjacent with j and m
        ! j adjacent with i and n
        IMPLICIT NONE

        ! 内部自变量
        integer :: i, j, x, m, icl, refed_iter
        integer :: m1, m2, m3, m4, w, w1, w2, w3, w4, w5, w6      ! 三角形和多边形中心点序号起始索引
        real(r8) :: newdbx(6, 2), newsjx(4, 2)  ! 数组大小是不一样的
        ! 外部读入变量
        integer,  dimension(:, :), allocatable, intent(in) :: ngrmw, ngrmm, ngrwm
        integer,  intent(in) :: iter
        integer,  dimension(:), intent(inout) :: num_mp, num_wp
        real(r8), dimension(:, :), allocatable, intent(inout) :: mp_new, wp_new
        integer,  dimension(:),    allocatable, intent(inout) :: mrl_new
        integer,  dimension(:, :), allocatable, intent(inout) :: ngrmw_new, ngrmm_new, ngr_mrl_new 
        
        ! 开始细化弱凹点 : two adjacent triangle in a polygon need to refine 
        refed_iter = 0
        do i = 2, num_mp(1), 1
            if (ref_sjx(i) == 0) cycle ! 如果三角形i不需要细化
            j = 0 ! 获取编号
            do x = 1, 3, 1
                ! 寻找弱凹点里面配对的三角形
                if (ref_sjx(ngrmm(x, i)) == 1) j = ngrmm(x, i)
                if (ngr_mrl_new(x, i) == 4) m = ngrmm(x, i)
            end do

            if (j == 0) cycle ! 说明不是弱凹点，在这里不细化
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            ! 建议把w1和w2，以及w3和w4的顺序换一换，保证逆时针顺序 或者直接把原本w5，w6的位置预留出来也不是不可以
            do x = 1, 3, 1
                ! 判断顶点w1到w3的位置 判断顶点是否在三角形j的顶点上
                ! 只有当所有这三个条件同时满足时，整个条件语句的结果才为 true
                if ((ngrmw(x, i) /= ngrmw(1, j)) .and. &
                    (ngrmw(x, i) /= ngrmw(2, j)) .and. & 
                    (ngrmw(x, i) /= ngrmw(3, j)) ) then 
                    w1 = ngrmw(x, i) ! 判断为真，说明x是三角形i relative to triangle j 的游离的顶点
                else ! 判断为faluse，说明x是三角形j的公共顶点之一
                    w = ngrmw(x, i) ! 判断公共顶点w是否也是三角形m的顶点
                    if ((w /= ngrmw(1, m)) .and. &
                        (w /= ngrmw(2, m)) .and. & 
                        (w /= ngrmw(3, m)) ) then ! 判断为正真，说明
                        w4 = w ! 三角形i和j的公共顶点，但是不在三角形m上
                    else
                        w2 = w ! 三角形i，j和m的公共顶点
                    end if
                end if
                ! 判断顶点w4的位置
                if ((ngrmw(x, j) /= ngrmw(1, i)) .and. &
                    (ngrmw(x, j) /= ngrmw(2, i)) .and. & 
                    (ngrmw(x, j) /= ngrmw(3, i)) ) then 
                    w3 = ngrmw(x, j) ! 判断为真，说明是三角形j relative to triangle i的游离的顶点
                end if
            end do
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! w5和w6，m1到m4是新增加的信息
            w5 = num_wp(iter - 1) + refed_iter * 2 + 1
            w6 = num_wp(iter - 1) + refed_iter * 2 + 2

            m1 = num_mp(iter - 1) + refed_iter * 4 + 1
            m2 = num_mp(iter - 1) + refed_iter * 4 + 2
            m3 = num_mp(iter - 1) + refed_iter * 4 + 3
            m4 = num_mp(iter - 1) + refed_iter * 4 + 4
            !!!!!!!!!!!!!! 
            ! print*, "w1 = ", w1
            ! print*, "w2 = ", w2
            ! print*, "w3 = ", w3
            ! print*, "w4 = ", w4
            ! print*, "w5 = ", w5
            ! print*, "w6 = ", w6
            !!!!!!!!!!!!!!
            ! m1
            ngrmw_new(1, m1) = w4
            ngrmw_new(2, m1) = w1
            ngrmw_new(3, m1) = w5
            ! m2
            ngrmw_new(1, m2) = w2
            ngrmw_new(2, m2) = w5
            ngrmw_new(3, m2) = w6
            ! m3
            ngrmw_new(1, m3) = w4
            ngrmw_new(2, m3) = w6
            ngrmw_new(3, m3) = w3
            ! m4
            ngrmw_new(1, m4) = w4
            ngrmw_new(2, m4) = w5
            ngrmw_new(3, m4) = w6

            ! 赋值，尽量不适用原来的数据
            newdbx(1, :) = wp_new(w1, :) ! m1
            newdbx(2, :) = wp_new(w2, :) ! m2
            newdbx(3, :) = wp_new(w3, :) ! m3
            newdbx(4, :) = wp_new(w4, :) ! m4
            icl = 0
            if (maxval(newdbx(1:4, 1)) - minval(newdbx(1:4, 1)) > 180.) then ! need to modify sjx first
                icl = 1! 判断是否存在跨越180经线，有则修正，并返回状态变量icl
                call CheckCrossing(4, newdbx)
            end if

            newdbx(5, 1:2) = (newdbx(1, 1:2) + newdbx(2, 1:2)) / 2. ! w1和w2连线中点为w5
            newdbx(6, 1:2) = (newdbx(2, 1:2) + newdbx(3, 1:2)) / 2. ! w2和w3连线中点为w6
            ! 新三角形有无编号的顺序要求呢？
            ! 两个组成弱凹点的三角形变为四个三角形
            newsjx(1, 1:2) = (newdbx(1, 1:2) + newdbx(4, 1:2) + newdbx(5, 1:2)) / 3.
            newsjx(2, 1:2) = (newdbx(2, 1:2) + newdbx(5, 1:2) + newdbx(6, 1:2)) / 3.
            newsjx(3, 1:2) = (newdbx(3, 1:2) + newdbx(4, 1:2) + newdbx(6, 1:2)) / 3.
            newsjx(4, 1:2) = (newdbx(4, 1:2) + newdbx(5, 1:2) + newdbx(6, 1:2)) / 3.
            ! 经度修正(对新生成的m或者w点的经度进行矫正)
            if (icl) then
                call CheckCrossing(6, newdbx)
                call CheckCrossing(4, newsjx)
            end if
            wp_new(w5,:) = newdbx(5, :)
            wp_new(w6,:) = newdbx(6, :)
            mp_new(m1:m4,:) = newsjx

            ref_sjx(i) = 0 ! 赋值为0，避免重复细化
            ref_sjx(j) = 0 ! 赋值为0，避免重复细化
            mrl_new(i) = 0 ! 赋值为0，避免重复细化
            mrl_new(j) = 0 ! 赋值为0，避免重复细化

            ! 将全部组成的弱凹点的两个三角形的邻域与顶点信息去掉
            ngrmm_new(:, i) = 1
            ngrmm_new(:, j) = 1
            ngrmw_new(:, i) = 1
            ngrmw_new(:, j) = 1
            refed_iter = refed_iter + 1
        end do
        call crossline_check(iter, mp_new, wp_new, num_mp, num_wp)
    END SUBROUTINE TwodivideFour
    
    SUBROUTINE OnedivideTwo(iter, ngrmw, ngrmm, ngrwm, num_mp, num_wp, mp_new, wp_new, mrl_new, ngrmw_new, ngrmm_new, ngr_mrl_new)
        ! 一分二算法（针对细化向非细化的过渡）
        IMPLICIT NONE
        ! 内部自变量
        integer :: i, j, k, icl, num_ref, refed_iter 
        integer :: m1, m2, w1, w2, w3, w4     
        real(r8) :: sjx(3, 2), tempa(1,2),tempb(1,2), tempc(1,2), hhh(5)! 
        ! 外部读入变量
        integer,  intent(in) :: iter
        integer,  dimension(:, :), allocatable, intent(in) :: ngrwm, ngrmm, ngrmw
        integer,  dimension(:), intent(inout) :: num_mp, num_wp
        real(r8), dimension(:, :), allocatable, intent(inout) :: mp_new, wp_new
        integer,  dimension(:),    allocatable, intent(inout) :: mrl_new
        integer,  dimension(:, :), allocatable, intent(inout) :: ngrmw_new, ngrmm_new, ngr_mrl_new

        ! 开始细化
        hhh = [1,2,3,1,2]
        refed_iter = 0
        do i = 2, num_mp(1), 1
            if (ref_sjx(i) == 0) cycle
            do j = 1, 3, 1
                ! 找到邻域中唯一一个被细化的三角形
                if (mrl_new(ngrmm(j, i)) == 4) k = ngrmm(j, i)
            end do
            do j = 1, 3, 1
                if( (ngrmw_new(j, i) /= ngrmw(1, k)) .and. & 
                    (ngrmw_new(j, i) /= ngrmw(2, k)) .and. &
                    (ngrmw_new(j, i) /= ngrmw(3, k)) )then ! 找到三角形i与细化三角公共边不相交的顶点
                    ! 根据公共边的位置，获取顶点编号，不建议用w1到w3应该w1和m1都要新增加的含义，但是这里只有w4是新增加的
                    w1 = ngrmw_new(hhh(j), i)
                    w2 = ngrmw_new(hhh(j+1), i)
                    w3 = ngrmw_new(hhh(j+2), i)
                    sjx(1, 1:2) = wp_new(w1, 1:2)
                    sjx(2, 1:2) = wp_new(w2, 1:2)
                    sjx(3, 1:2) = wp_new(w3, 1:2)
                    !!!!!!!!!!!!!!!!!!!!!!!! 
                    ! print*, "w1 = ", w1
                    ! print*, "w2 = ", w2
                    ! print*, "w3 = ", w3
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
            mrl_new(i) = 0
            ngr_mrl_new(:, i) = 0
            ! mrl_new(m1:m2) = 2! 记录细化方式
            ngrmm_new(:, i) = 1! 去掉原来三角形的中心点与顶点信息
            ngrmw_new(:, i) = 1
            refed_iter = refed_iter + 1
        end do
        call crossline_check(iter, mp_new, wp_new, num_mp, num_wp)
    END SUBROUTINE OnedivideTwo

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
        do i = 2, num_dbx, 1
            num_inter = n_ngr(i) ! 获取单个多边形的顶点总数 ! use for sort points 
            allocate(angle(num_inter)); angle = 0.
            allocate(points(num_inter, 2)); points = mp_f(ngr(1:num_inter, i), 1:2)
            center  = 0.
            ! 是否跨越180经线，存在跨越
            if(maxval(points(:,1))-minval(points(:,1))>180.) call CheckCrossing(num_inter, points)

            do j = 1, num_inter, 1
                center = center + points(j, :)
            end do
            center = center / num_inter

            do j = 1, num_inter, 1
                points(j, :) = points(j, :) - center
                angle(j) = atan2(points(j, 2),  points(j, 1))
            end do
            if (all((angle(2:num_inter) - angle(1:num_inter-1)) >= 0.0)) then
                deallocate(angle, points)
                cycle
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


    SUBROUTINE TriMeshQuality(num_sjx, wp_f, ngrmw_f, length_sjx_temp, angle_sjx_temp, Extr_sjx_temp, Eavg_sjx_temp, Savg_sjx_temp, less30_temp)

        implicit none ! 定义的顺序最好是只读入的参数，内部参数，需要传出的参数
        integer, intent(in) :: num_sjx
        real(r8), dimension(:, :), intent(in) :: wp_f
        integer,  dimension(:, :), intent(in) :: ngrmw_f
        integer :: i
        real(r8) :: length_temp(3), angle_temp(3), sjx(3, 2) !angle_regular
        real(r8), dimension(:, :), intent(inout) :: length_sjx_temp, angle_sjx_temp
        real(r8), intent(inout) :: Extr_sjx_temp(2), Eavg_sjx_temp(2), Savg_sjx_temp, less30_temp
        ! 注意跨越180的情况
        do i = 2, num_sjx, 1 ! 从2开始，因为第一个不存在
            sjx = wp_f(ngrmw_f(:, i), 1:2)
            if (maxval(sjx(:, 1)) - minval(sjx(:, 1)) > 180.) call CheckCrossing(3, sjx)
            CALL GetTriangleLength(sjx, length_temp)
            CALL GetAngle(3, sjx, angle_temp)! 计算多边形内角和
            length_sjx_temp(i, :) = length_temp 
            angle_sjx_temp(i, :) = angle_temp
            Eavg_sjx_temp(1) = Eavg_sjx_temp(1) + minval(angle_temp)
            Eavg_sjx_temp(2) = Eavg_sjx_temp(2) + maxval(angle_temp)
            Savg_sjx_temp    = Savg_sjx_temp    + sum((angle_temp - 60.)**2)
            if (minval(angle_temp) < 30.) less30_temp = less30_temp + 1
        end do
        Extr_sjx_temp(1) = minval(angle_sjx_temp(2:num_sjx, :))
        Extr_sjx_temp(2) = maxval(angle_sjx_temp(2:num_sjx, :))
        Eavg_sjx_temp(:) = Eavg_sjx_temp(:) / (num_sjx - 1)
        Savg_sjx_temp    = sqrt(Savg_sjx_temp / (3 * num_sjx - 3))
        less30_temp      = less30_temp / (num_sjx - 1)

    END SUBROUTINE TriMeshQuality

    SUBROUTINE PolyMeshQuality(num_edges, num_dbx, mp_f, ngrwm_f, n_ngrwm_f, angle_dbx_temp, Extr_dbx_temp, Eavg_dbx_temp, Savg_dbx_temp)

        implicit none ! 定义的顺序最好是只读入的参数，内部参数，需要传出的参数
        integer, intent(in) :: num_edges, num_dbx
        real(r8), dimension(:, :), intent(in) :: mp_f
        integer,  dimension(:, :), intent(in) :: ngrwm_f
        integer,  dimension(:), allocatable, intent(in) :: n_ngrwm_f
        integer :: i, j
        real(r8) :: angle_regular
        real(r8), allocatable :: angle_temp(:), dbx(:,:)
        real(r8), dimension(:,:), intent(inout) :: angle_dbx_temp
        real(r8), intent(inout) :: Extr_dbx_temp(2), Eavg_dbx_temp(2), Savg_dbx_temp
        allocate(angle_temp(num_edges)); angle_temp = 0.
        allocate(dbx(num_edges, 2)); dbx = 0.
        j = 0
        angle_regular = (num_edges-2) * 180. / num_edges !正多边形内角度数公式
        do i = 2, num_dbx, 1 ! 为啥这里从2开始，而不是从1开始
            if (n_ngrwm_f(i) /= num_edges) cycle
            j = j + 1 ! 用于计数
            dbx = mp_f(ngrwm_f(1:num_edges, i), :)
            if (maxval(dbx(:, 1)) - minval(dbx(:, 1)) > 180.) call CheckCrossing(num_edges, dbx)
            CALL GetAngle(num_edges, dbx, angle_temp)! 计算多边形内角和
            angle_dbx_temp(j, :) = angle_temp
            Eavg_dbx_temp(1) = Eavg_dbx_temp(1) + minval(angle_temp)
            Eavg_dbx_temp(2) = Eavg_dbx_temp(2) + maxval(angle_temp)
            Savg_dbx_temp    = Savg_dbx_temp    + sum((angle_temp - angle_regular)**2)
        end do
        Extr_dbx_temp(1) = minval(angle_dbx_temp)
        Extr_dbx_temp(2) = maxval(angle_dbx_temp)
        Eavg_dbx_temp(:) = Eavg_dbx_temp(:) / j
        Savg_dbx_temp    = sqrt(Savg_dbx_temp / (num_edges * j))

    END SUBROUTINE PolyMeshQuality

    SUBROUTINE GetTriangleDis(sjx_points, ngrmw, ngrwm, n_ngrwm, ref_sjx, mp_dis)
        implicit none
        integer, intent(in) :: sjx_points
        integer, dimension(:,:), intent(in) :: ngrmw, ngrwm
        integer, dimension(:),   intent(in) :: n_ngrwm
        integer, dimension(:),   intent(in) :: ref_sjx
        integer :: i, j, k, w, num(5),l, m, ncid, spDimID, ncvarid
        integer, dimension(:,:), intent(out) :: mp_dis
        character(LEN = 256) :: lndname
        character(LEN = 5) :: nxpc, stepc
        write(nxpc, '(I4.4)') NXP
        write(stepc, '(I2.2)') step
        mp_dis = 9
        num = 0

        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i,j,k,w)
        do i = 2, sjx_points, 1
            if (ref_sjx(i) == 0) cycle
            mp_dis(i, i) = 0
            do j = 1, 3, 1
                w = ngrmw(j, i)
                do k = 1, n_ngrwm(w), 1
                    if (mp_dis(i, ngrwm(k, w)) == 9) then
                        mp_dis(i, ngrwm(k, w)) = 1
                        mp_dis(ngrwm(k, w), i) = 1
                    end if
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        do m = 1, 4, 1
            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
            !$OMP PRIVATE(i,j,k,l,w)
            do i = 2, sjx_points, 1
                if (ref_sjx(i) == 0) cycle
                do j = 2, sjx_points, 1
                    if (mp_dis(i, j) == m) then
                        do k = 1, 3, 1
                            w = ngrmw(k, j)
                            do l = 1, n_ngrwm(w), 1
                                if (mp_dis(i, ngrwm(l, w)) == 9) then
                                    mp_dis(i, ngrwm(l, w)) = m + 1
                                    mp_dis(ngrwm(l, w), i) = m + 1
                                end if
                            end do
                        end do
                    end if
                end do
            end do
            !$OMP END PARALLEL DO
        end do

        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc)  // "_" // trim(stepc) // "_mp_dis.nc4"
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", sjx_points, spDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "mp_dis", NF90_INT, (/ spDimID, spDimID /), ncVarID))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid, mp_dis))
        CALL CHECK(NF90_CLOSE(ncID))
        
    END SUBROUTINE GetTriangleDis

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
        integer,  allocatable :: hh(:)          ! 计数三角形个数
        integer,  allocatable :: mp_ref(:)          ! 记录被细化三角形的索引
        character(LEN = 256) :: lndname
        character(LEN = 5) :: nxpc, stepc
        write(nxpc, '(I4.4)') NXP
        write(stepc, '(I2.2)') step
        ! m点 增加，不会重复；减少（只在原来的范围内）m点具有唯一性
        ! w点 增加，会重复，不同编号对应同一个点位；不会减少 w点没有唯一性


        print*, "wp_f start"
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
        print*, "细化前共有", num_wp(1), "个多边形网格"
        print*, "细化后共有", num_wp(iter), "个多边形网格"
        print*, "去除重复点后，还剩", num_dbx, "个多边形网格"
        print*, "wp_f finfish"
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc)  // "_" // trim(stepc) // "_wp_f.nc4"
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
        CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", num_wp(iter), lpDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 2, DimaID))
        CALL CHECK(NF90_DEF_VAR(ncID, "vertex_mapping", NF90_INT, (/ lpDimID /), ncVarID(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "wp_f", NF90_FLOAT, (/ lpDimID, DimaID /), ncVarID(2)))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(1), vertex_mapping))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(2), wp_f))
        CALL CHECK(NF90_CLOSE(ncID))
        print*, "max(vertex_mapping) = ", maxval(vertex_mapping)
         if (maxval(vertex_mapping) /= num_dbx) stop "maxval(vertex_mapping) /= num_dbx"

        ! 更新 mp_f ! 如果ngrmw_new不存在，则这个三角形不存在 
        print*, "重新计算ngrmw_new and mp_f，并储存mp_f" 
        ! m点的特点是会增加，也是减少（只在原来的范围内），但是不会重复！！！！！！！！！！！
        ! 统计初始三角形中被细化的个数，此时的ngrmw_new 还没有进行重新编号
        print*, "mp_f start"
        num_ref = 0
        allocate(hh(num_mp(1))); hh = 1
        do i = 2, num_mp(1), 1
            if (ngrmw_new(1, i) == 1) then ! 当三角顶点不存在的时候，跳过
                num_ref = num_ref + 1
                hh(num_ref) = i
            end if
        end do
        allocate(mp_ref(num_ref)); mp_ref = hh(1:num_ref) ! 记录被细化三角形的索引
        deallocate(hh)

        num_sjx = num_mp(iter) - num_ref ! 获取三角形总数
        allocate(mp_f(num_sjx, 2)); mp_f = 0. ! 经纬度
        allocate(ngrmw_f(3, num_sjx)); ngrmw_f = 1 ! 顶点编号
        k = 1
        do i = 2, num_mp(iter), 1
            if (ngrmw_new(1, i) == 1) cycle ! 跳过原来三角形中被细化的那部分三角形，
            k = k + 1 ! 累加，用于进位
            mp_f(k, :) = mp_new(i, :)
            ngrmw_f(:, k) = ngrmw_new(:, i)
        end do
        print*, "细化前共有", num_mp(1), "个三角形网格"
        print*, "细化后共有", num_mp(iter), "个三角形网格"
        print*, "去除重复点后，还剩", num_sjx, "个三角形网格"
        print*, "mp_f finfish"
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc)  // "_" // trim(stepc) // "_mp_f.nc4"
        print*, lndname
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
        allocate(ngrwm_f(7, num_dbx)); ngrwm_f   = 1 ! 记录相邻三角形编号，初始化为1
        allocate(n_ngrwm_f(num_dbx));  n_ngrwm_f = 0 ! 记录相邻三角形
        do i = 2, num_sjx, 1 ! 三角形总数（含不存在的三角形）
            do j = 1, 3, 1
                ngrmw_f(j, i) = vertex_mapping(ngrmw_f(j, i))
                k = ngrmw_f(j, i) ! 获取多边形中心点（或三角形顶点）编号信息
                n_ngrwm_f(k) = n_ngrwm_f(k) + 1 ! 累加顶点个数
                ngrwm_f(n_ngrwm_f(k), k) = i
            end do
        end do
        CALL GetSort(n_ngrwm_f, num_dbx, mp_f, ngrwm_f)
        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc)  // "_" // trim(stepc) // "_ngrwm_f.nc4"
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
        CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", num_dbx, lpDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_b", 7, DimbID))
        CALL CHECK(NF90_DEF_VAR(ncID, "ngrwm_f", NF90_INT, (/ DimbID, lpDimID /), ncVarID(1)))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(1), ngrwm_f))
        CALL CHECK(NF90_CLOSE(ncID))

        lndname = trim(file_dir) // "tmpfile/gridfile_NXP" // trim(nxpc)  // "_" // trim(stepc) // "_ngrmw_f.nc4"
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", num_sjx, spDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dim_a", 3, DimaID))
        CALL CHECK(NF90_DEF_VAR(ncID, "ngrmw_f", NF90_INT, (/ DimbID, spDimID /), ncVarID(1)))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncvarid(1), ngrmw_f))
        CALL CHECK(NF90_CLOSE(ncID))

    END SUBROUTINE NGR_RENEW

END module MOD_refine
