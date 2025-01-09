!--------------------------------------
! 多边形分结构网格初步细化(单层)
! ***       初始数据
! ***_new   更新后数据
! ***_f     最终数据
!--------------------------------------

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
        integer :: w1, w2, w3, w4, w5, w6, m1, m2, m3, m4      ! 三角形和多边形中心点序号
        integer :: row, col, error,edges 
        integer :: ustr_points                                 ! 非结构网格数量
        integer :: num_ref                                     ! 每次细化三角形数 
        integer :: num_all                                     ! 非结构网格包含结构网格总数
        integer :: refed(1000)                                 ! 记录细化中每步（或每次迭代）细化的三角形数
        integer :: iter                                        ! 网格细化次数
        integer :: sa_iter                                     ! 网格质量调整迭代次数
        integer :: num_sjx, num_dbx                            ! 细化后三角形、多边形数量
        integer :: sjx_points, lbx_points                      ! 三角形与多边形数量(读取的原文件)
        integer :: nmp(1000), nwp(1000)                        ! 记录每次细化后的m，w点数量
        integer :: icl(3)                                      ! 判断三角形网格是否穿过180°经线
        !integer :: maxlc                                       ! 土地类型最大编号，17 or 24

        ! nc文件读写相关参数
        integer :: ncid, varid(10), ncvarid(13), upDimID, iunit
        integer :: spDImID, lpDimID, twDimID, thDimID, seDimID, sxDimID, dimID_sjx, dimID_lbx

        integer :: n_wbx, n_lbx, n_qbx                         ! 五边形、六边形、七边形数量

        real(r8), allocatable :: wp(:, :), mp(:, :)            ! 三角形、多边形网格中心点初始数据
        real(r8), allocatable :: mp_new(:, :), wp_new(:, :)    ! 三角形、多边形网格中心点更新数据
        real(r8), allocatable :: mp_f(:, :), wp_f(:, :)        ! 三角形、多边形网格中心点最终数据 
        real(r8), allocatable :: mp_f_tmp(:, :), wp_f_tmp(:, :)! mp_f/wp_f数组的缓存数据 
        real(r8), allocatable :: ref_lbx(:, :)                 ! 构成多边形的三角形细化情况
        real(r8), allocatable :: dismm(:,:),disww(:,:)         ! 相邻m/w点距离
        
        ! 阈值相关数组
        real(r8), allocatable :: f_mainarea(:, :), slope(:, :), lai(:, :)
        real(r8), allocatable :: k_s(:, :, :), k_sl(:, :, :),tkdry(:, :, :)
        real(r8), allocatable :: tksatf(:, :, :), tksatu(:, :, :)

        ! 网格质量调整相关
        ! Extr_*** 指网格的最小角度与最大角度
        ! Eavg_*** 指网格的最小角度与最大角度平均值
        ! Savg_*** 指网格角度与正多边形角度的标准差
        ! less30   指三角形网格中小于30度角的数量
        real(r8), allocatable :: less30(:), length(:, :, :), angle(:, :, :)
        real(r8), allocatable :: Savg_sjx(:), Extr_sjx(:, :), Eavg_sjx(:, :)
        real(r8), allocatable :: angle_wbx(:, :, :), angle_lbx(:, :, :), angle_qbx(:, :, :)
        real(r8), allocatable :: Savg_wbx(:), Savg_lbx(:), Savg_qbx(:)
        real(r8), allocatable :: Eavg_wbx(:, :), Eavg_lbx(:, :), Eavg_qbx(:, :)
        real(r8), allocatable :: Extr_wbx(:, :), Extr_lbx(:, :), Extr_qbx(:, :)

        real(r8), allocatable :: MoveDis(:, :)                 ! 记录w点在x、y方向上的调整距离

        real(r8) :: sjx(3, 2), newsjx(3, 2), dbx(7, 2)         ! 记录三角形/多边形顶点
        real(r8) :: pi
        real(r8) :: rx, ry, fra                 ! 网格质量调整相关
        integer, allocatable :: n_landtypes(:)  ! 土地利用类型

        integer, allocatable :: ref(:)          ! 三角形网格初始细化情况
        integer, allocatable :: ref_l(:, :)     ! 多边形网格细化情况
        integer, allocatable :: ref_th(:, :)    ! 初始三角形网格阈值细化情况
        integer, allocatable :: ref_tr(:, :)    ! 最终三角形网格阈值细化情况
        integer, allocatable :: ref_pl(:, :)    ! 记录多边形是因何种阈值细化

        integer, allocatable :: ngrmw(:, :), ngrwm(:, :)          ! m/w点相邻的w/m点索引(细化前)
        integer, allocatable :: ngrmw_new(:, :), ngrwm_new(:, :)  ! m/w点相邻的w/m点索引(细化后)
        integer, allocatable :: ngrmw_f(:, :), ngrwm_f(:, :)      ! m/w点相邻的w/m点索引(最终)
        integer, allocatable :: ngrmw_f_tmp(:, :), ngrwm_f_tmp(:, :)
        integer, allocatable :: ngrmm(:, :)        ! m点相邻的m点索引(细化前)
        integer, allocatable :: ngrmm_new(:, :)    ! m点相邻的m点索引(细化后)
        integer, allocatable :: mrl(:)             ! 三角形网格细化程度(细化前)
        integer, allocatable :: mrl_new(:)         ! 三角形网格细化程度(细化后) 
        integer, allocatable :: mrl_f(:)           ! 三角形网格细化程度(最终)
        integer, allocatable :: ngr_mrl(:, :)      ! 三角形网格的相邻三角形网格点的mrl(细化前)
        integer, allocatable :: ngr_mrl_new(:, :)  ! 三角形网格的相邻三角形网格点的mrl(细化后)
        integer, allocatable :: mp_ref(:)          ! 记录被细化三角形的索引

        character(LEN = 256) :: lndname, nxpc

        logical :: isexist                ! 判断细化后是否存在重复w点
        
        ! 防止细化交汇带出现冲突，一分为四
        logical :: iterA                  ! 当一轮迭代中迭代B与迭代C通过且无细化，迭代A通过
        logical :: iterB                  ! 从三角形网格进行判断
        logical :: iterC                  ! 从多边形网格进行判断
        logical :: End_SpringAjustment    ! 判断网格质量调整是否结束

        logical,allocatable :: IsInRfArea(:) ! 判断三角形网格是否位于细化区域内 
        
        ! 便于函数调用的缓存数组
        real(r8) :: tmpa(2),tmpb(2),tmpc(2),tmp_angle(3),tmp_length(3),tmp_angle7(7)
        
        ! nc文件读取参数名
        !character(LEN = 20) :: p_name(6) = (/"GLONW", "GLATW", "GLONM", "GLATM", "itab_w%im", "itab_m%iw"/)
        character(LEN = 20),dimension(6) :: p_name

        p_name = [character(len=20) :: "GLONW", "GLATW", "GLONM", "GLATM", "itab_w%im", "itab_m%iw"]
        pi = 3.1415926535

        iterA = .false.
        iterB = .false.
        iterC = .false.

        !-------------------------------------------
        ! read unstructure mesh
        !-------------------------------------------
        print*, "start to read unstructure mesh data"
        print*, ""
         
        ! 读取未细化初始网格数据
        write(nxpc, '(I4.4)') NXP
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

        allocate(wp(lbx_points, 2))            ! Initial data for center point of polygon mesh (多边形网格中心点初始数据)
        allocate(wp_new(lbx_points * 4, 2))    ! Update data at center point of polygon mesh (多边形网格中心点更新数据)
        allocate(mp(sjx_points, 2))            ! Initial data at the center point of the triangular grid (三角形网格中心点初始数据)
        allocate(mp_new(sjx_points * 4, 2))    ! The center point of the triangular grid updates the data (三角形网格中心点更新数据)
        allocate(ngrwm(7, lbx_points))         ! wp initial index table of adjacent mp points (wp的相邻mp点初始索引表)
        allocate(ngrwm_new(7, lbx_points * 4)) ! wp's adjacent mp points update the index table (wp的相邻mp点更新索引表)
        allocate(ngrmw(3, sjx_points))         ! The initial index table of adjacent wp points of mp (mp的相邻wp点初始索引表)
        allocate(ngrmw_new(3, sjx_points * 4)) ! The adjacent wp points of mp update the index table (mp的相邻wp点更新索引表)

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

        allocate(mrl(sjx_points))                  ! Triangle mesh refinement degree (三角形网格细化程度/方式)
        allocate(mrl_new(sjx_points * 4))          ! Updated mrl data (更新后数据)
        allocate(ngr_mrl(3, sjx_points))           ! mrl of adjacent triangular mesh points of a triangular mesh (三角形网格的相邻三角形网格的mrl)
        allocate(ngr_mrl_new(3, sjx_points * 4))   ! Updated ngr_mrl data (更新后数据)
        allocate(ngrmm(3, sjx_points))             ! mp adjacent mp initial index table (m点相邻m点的初始索引表)
        allocate(ngrmm_new(3, sjx_points * 4))     ! Updated ngrmm data (更新后数据)

        mrl = 1
        mrl_new = 1
        ngr_mrl = 1
        ngr_mrl_new = 1
        ngrmm = 1
        ngrmm_new = 1

        ! 规定180°经线圈上的经度为180°
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
            if(ngrmw(3, i) == 1)then ! 三角形不存在，mrl设置为0
                mrl(i) = 0
            end if
        end do

        do i = 1, sjx_points, 1
            do j = 1, sjx_points, 1
                if(i /= j)then   ! 若它们不是同一个m点
                    error = IsNgrmm(ngrmw(1:3, i), ngrmw(1:3, j))
                    if(error /= 0)then    ! error表示两个三角形公共边是i的第几个w顶点的对边
                        ngrmm(error, i) = j  ! i号m点的第error个w顶点的对边的另一个相邻m点索引为j
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

        ! 更新相邻三角形的mrl情况（暂时只有0或1）
        ! 0表示三角形不存在
        ! 1表示是未细化的初始三角形
        do i = 1, sjx_points, 1
            do j = 1, 3, 1
                if(ngrmm(j, i) /= 1)then
                    if(mrl(ngrmm(j, i)) == 0)then
                        ngr_mrl(j, i) = 0
                    end if
                end if
            end do
        end do

        ! Allocate threshold dependent array memory (分配阈值相关数组内存)
        ! 第二维记录均值与标准差
        ! 第三维记录土壤的第一层与第二层
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
        ! Reading threshold data (latitude and longitude grid) (读取初始网格的阈值数据（经纬度网格）)
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

        allocate(ref(sjx_points))             ! 记录三角形网格是否需要细化
        allocate(ref_l(lbx_points, 7))        ! Polygon mesh refinement (多边形网格细化情况)
        allocate(ref_lbx(lbx_points, 8))      ! 用于细化时记录信息
        allocate(ref_th(sjx_points, 16))      ! 记录三角形网格因何种阈值细化
        allocate(ref_tr(sjx_points * 4, 16))  ! The final triangular mesh threshold refinement (最终三角形网格阈值细化情况)
        ref = 0
        ref_l = 0
        ref_lbx = 0.
        ref_th = 0
        ref_tr = 0

        ! 判断哪些三角形网格位于细化区域内
        allocate(IsInRfArea(sjx_points))
        IsInRfArea = .false.
        CALL IsInRefineArea(IsInRfArea,sjx_points,lbx_points,ngrmw,wp)

        !if(lcs == "igbp")then
        !   maxlc = 17
        !else
        !   maxlc = 24
        !end if

        !--------------------------------------------------
        ! 1.1 Used to record triangles that  == ire preliminary refinement （记录需要初步细化的三角形）
        !--------------------------------------------------
        do i = 1, sjx_points, 1
            !if((wp(ngrmw(i,1),2)>85.).or.(wp(ngrmw(i,2),2)>85.).or.(wp(ngrmw(i,3),2)>85.))then
            !   cycle
            !end if
            if(IsInRfArea(i) .eqv. .false.)then ! 当三角形网格不位于细化区域内
               cycle
            end if
            
            ! 当主导土地类型为海洋
            if((f_mainarea(i, 1) == 0.).or.(f_mainarea(i, 1) == maxlc))then
                cycle
            end if

            if (refine_num_landtypes .eqv. .true. .and. n_landtypes(i)>th_num_landtypes) then
                ref(i) = 1
                ref_th(i, 1) = 1
            end if

            if (refine_area_mainland .eqv. .true. .and. f_mainarea(i, 2) < th_area_mainland) then
                ref(i) = 1
                ref_th(i, 2) = 1
            end if

            if (refine_lai_m .eqv. .true. .and. lai(i, 1) > th_lai_m) then
                ref(i) = 1
                ref_th(i, 3) = 1
            end if

            if (refine_lai_s .eqv. .true. .and. lai(i, 2) > th_lai_s) then
                ref(i) = 1
                ref_th(i, 4) = 1
            end if

            if (refine_slope_m .eqv. .true. .and. slope(i, 1) > th_slope_m) then
                ref(i) = 1
                ref_th(i, 5) = 1
            end if

            if (refine_slope_s .eqv. .true. .and. slope(i, 2) > th_slope_s) then
                ref(i) = 1
                ref_th(i, 6) = 1
            end if

            if (refine_k_s_m .eqv. .true. .and. ((k_s(i, 1, 1) > th_k_s_m).or.(k_s(i, 1, 2) > th_k_s_m))) then
                ref(i) = 1
                ref_th(i, 7) = 1
            end if

            if (refine_k_s_s .eqv. .true. .and. ((k_s(i, 2, 1) > th_k_s_s).or.(k_s(i, 2, 2) > th_k_s_s))) then
                ref(i) = 1
                ref_th(i, 8) = 1
            end if

            if (refine_k_solids_m .eqv. .true. .and. ((k_sl(i, 1, 1) > th_k_solids_m).or.(k_sl(i, 1, 2) > th_k_solids_m))) then
                ref(i) = 1
                ref_th(i, 9) = 1
            end if

            if (refine_k_solids_s .eqv. .true. .and. ((k_sl(i, 2, 1) > th_k_solids_s).or.(k_sl(i, 2, 2) > th_k_solids_s))) then
                ref(i) = 1
                ref_th(i, 10) = 1
            end if

            if (refine_tkdry_m .eqv. .true. .and. ((tkdry(i, 1, 1) > th_tkdry_m).or.(tkdry(i, 1, 2) > th_tkdry_m))) then
                ref(i) = 1
                ref_th(i, 11) = 1
            end if

            if (refine_tkdry_s .eqv. .true. .and. ((tkdry(i, 2, 1) > th_tkdry_s).or.(tkdry(i, 2, 2) > th_tkdry_s))) then
                ref(i) = 1
                ref_th(i, 12) = 1
            end if

            if (refine_tksatf_m .eqv. .true. .and. ((tksatf(i, 1, 1) > th_tksatf_m).or.(tksatf(i, 1, 2) > th_tksatf_m))) then
                ref(i) = 1
                ref_th(i, 13) = 1
            end if

            if (refine_tksatf_s .eqv. .true. .and. ((tksatf(i, 2, 1) > th_tksatf_s).or.(tksatf(i, 2, 2) > th_tksatf_s))) then
                ref(i) = 1
                ref_th(i, 14) = 1
            end if

            if (refine_tksatu_m .eqv. .true. .and. ((tksatu(i, 1, 1) > th_tksatu_m).or.(tksatu(i, 1, 2) > th_tksatu_m))) then
                ref(i) = 1
                ref_th(i, 15) = 1
            end if

            if (refine_tksatu_s .eqv. .true. .and. ((tksatu(i, 2, 1) > th_tksatu_s).or.(tksatu(i, 2, 2) > th_tksatu_s))) then
                ref(i) = 1
                ref_th(i, 16) = 1
            end if

        end do

        iter = 1                                 ! 本次细化中的迭代次数
        num_ref = INT(sum(ref))                  ! 需要细化的三角形数
        nmp(1) = sjx_points + 4 * num_ref        ! 记录每次迭代后三角形数
        nwp(1) = lbx_points + 3 * num_ref        ! 记录每次迭代后多边形数
        ! 每细化（一分为四）一个三角形，增加4个m点和3个w点

        ! 将初始数据导入内存更大的新数组中
        mrl_new(1:sjx_points) = mrl
        ngr_mrl_new(1:3, 1:sjx_points) = ngr_mrl
        mp_new(1:sjx_points, 1:2) = mp
        wp_new(1:lbx_points, 1:2) = wp
        ngrmw_new(1:3, 1:sjx_points) = ngrmw(1:3, 1:sjx_points)
        ngrmm_new(1:3, 1:sjx_points) = ngrmm(1:3, 1:sjx_points)
        ngrwm_new(1:7, 1:lbx_points) = ngrwm(1:7, 1:lbx_points)

        ! 清除阈值数组内存
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
        ! 1.2 Preliminary refinement (one into four) 【初步细化（一分为四）】
        !--------------------------------------------------
        print*, "Start to refine (开始初步细化)"
        print*, "iter =", iter, "num =", num_ref
        refed = 0
        do i = 1, sjx_points, 1
            if(ref(i) == 1)then ! 若三角形需要细化
                icl = 0
                sjx = 0.
                newsjx = 0.

                ! 读取三角形顶点经纬度
                sjx(1, 1:2) = wp(ngrmw(1, i), 1:2)
                sjx(2, 1:2) = wp(ngrmw(2, i), 1:2)
                sjx(3, 1:2) = wp(ngrmw(3, i), 1:2)

                ! 判断是否越过180°经线圈
                icl(1) = IsCrossLine(sjx(2, 1), sjx(3, 1))
                icl(2) = IsCrossLine(sjx(1, 1), sjx(3, 1))
                icl(3) = IsCrossLine(sjx(1, 1), sjx(2, 1))

                ! 若越线，则将小于0°经度加上360°
                if(INT(sum(icl)) > 0)then
                    do j = 1, 3, 1
                        if(sjx(j, 1) < 0.)then
                            sjx(j, 1) = sjx(j, 1) + 360.
                        end if
                    end do
                end if

                ! 计算新生成的中间小三角形顶点（原三角形中点）
                newsjx(1, 1:2) = (sjx(2, 1:2) + sjx(3, 1:2)) / 2.
                newsjx(2, 1:2) = (sjx(1, 1:2) + sjx(3, 1:2)) / 2.
                newsjx(3, 1:2) = (sjx(1, 1:2) + sjx(2, 1:2)) / 2.

                ! 为新生成的m点与w点编序号
                m1 = sjx_points + refed(1) * 4 + 1
                m2 = sjx_points + refed(1) * 4 + 2
                m3 = sjx_points + refed(1) * 4 + 3
                m4 = sjx_points + refed(1) * 4 + 4

                w1 = lbx_points + refed(1) * 3 + 1
                w2 = lbx_points + refed(1) * 3 + 2
                w3 = lbx_points + refed(1) * 3 + 3

                ! 计算新生成的m点和w点经纬度
                mp_new(m1, 1:2) = (sjx(1, :) + newsjx(2, :) + newsjx(3, :)) / 3.
                mp_new(m2, 1:2) = (sjx(2, :) + newsjx(1, :) + newsjx(3, :)) / 3.
                mp_new(m3, 1:2) = (sjx(3, :) + newsjx(1, :) + newsjx(2, :)) / 3.
                mp_new(m4, 1:2) = (newsjx(3, :) + newsjx(1, :) + newsjx(2, :)) / 3.

                wp_new(w1, 1:2) = newsjx(1, 1:2)
                wp_new(w2, 1:2) = newsjx(2, 1:2)
                wp_new(w3, 1:2) = newsjx(3, 1:2)

                ! 更新原三角形与新三角形的ngrmm,ngrmw索引
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

                ! 将新生成m点、w点大于180°的经度减小360°
                if(INT(sum(icl)) > 0)then
                    do k = 1, 4, 1
                        call CheckLon(mp_new(sjx_points + refed(1) * 4 + k, 1))
                        call CheckLon(wp_new(lbx_points + refed(1) * 3 + k, 1))
                    end do
                end if

                mrl_new(i) = 4 ! mrl==4 三角形网格被平均分为四份
                ! 将由被平均分为四份生成的小三角形mrl也记录为4
                mrl_new(sjx_points + refed(iter) * 4 + 1:sjx_points + refed(1) * 4 + 4) = 4

                ! 根据更新的相邻m点关系更新相邻mrl关系
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
                !         ngrwm_new(ngrmw(i,k),j) = sjx_points+refed(iter)*4+k
                !      end if
                !   end do
                !end do

                ! 记录w点相邻三角形有被细化情况发生
                ref_lbx(8, ngrmw_new(1, i)) = 1
                ref_lbx(8, ngrmw_new(2, i)) = 1
                ref_lbx(8, ngrmw_new(3, i)) = 1

                ! 记录小三角形是因何种阈值被细化产生的
                do j = 1, 16, 1
                    if(ref_th(i, j) == 1)then
                        ref_tr(m1, j) = 1
                        ref_tr(m2, j) = 1
                        ref_tr(m3, j) = 1
                    end if
                end do

                ! 第一次细化的操作数（或迭代次数）加一
                refed(1) = refed(1) + 1
            end if
        end do
        print*, "itered_num =", refed(1)
        print*, "Refining step 1 is complete (细化第一步完成)"

        !--------------------------------------------------
        ! 1.3 储存初始网格数据和初步细化后新网格数据
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
        ! 2.1 进行迭代（防止细化交汇带出现冲突，一分为四）
        !--------------------------------------------------
        print*, "iterA start"
        do while(iterA .eqv. .false.) ! 当iterA为true时，该步骤完成

            iterA = .true.    ! 判断迭代B和迭代C是否都已满足条件
            iterB = .false.   ! 从三角形网格进行判断
            iterC = .false.   ! 从多边形网格进行判断

            print*, "iterB start"
            do while(iterB .eqv. .false.)

                ref = 0 ! 记录三角形是否需要被细化，初始化

                ! 当一个未细化三角形的相邻三角形有两个或两个以上已经被细化时，细化该三角形
                do i = 1, sjx_points, 1
                    if(mrl_new(i) == 1)then  ! 用相邻三角形mrl的和来判断相邻三角形细化数量
                        if((sum(ngr_mrl_new(1:3, i)) == 12.).or.(sum(ngr_mrl_new(1:3, i)) == 9.))then
                            ref(i) = 1
                        end if
                    end if
                end do

                ! 需要细化的三角形总数
                num_ref = INT(sum(ref))

                ! 当产生新的细化三角形时，iterA不通过
                if(num_ref == 0)then
                    iterB = .true.
                else
                    iterA = .false.
                    iter = iter + 1

                    ! 将记录的这些三角形一分为四（类似以上操作）
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
            do while(iterC .eqv. .false.)

                ref = 0

                do i = 1, lbx_points, 1

                    edges = 0      ! 多边形网格边数

                    do j = 1, 7, 1   ! 计算未细化三角形构成的多边形的边数
                        if(ngrwm(j, i) /= 1)then
                            m1 = ngrwm(j, i)
                            if(sum(ngr_mrl_new(:, m1)) == 6.)then ! 1+1+4，即相邻一个细化两个未细化（图a）
                                ref_lbx(i, j) = 1  ! 表示i号w点的第j个相邻三角形将要被细化（一分为二）（图e）
                            end if
                            edges = edges + 1 ! 记录多边形边数
                        end if
                    end do

                    if(ref_lbx(i, 8) == 0)then     ! 表示i号w点的相邻三角形均未被细化
                        if(edges == 5)then         ! 五边形(图b)
                            do j = 1, 5, 1
                                m1 = ngrwm(j, i)
                                m2 = ngrwm(j + 1, i)
                                if(j == 5)then
                                    m2 = ngrwm(1, i)
                                end if ! m1，m2为多边形相邻顶点
                                ! m1,m2均为1+1+4结构，表示这两个三角形都有两个相邻三角形未被细化
                                ! 因此它们组成弱凹(图c)
                                if((sum(ngr_mrl_new(:, m1)) == 6).and.(sum(ngr_mrl_new(:, m2)) == 6.))then
                                    ref_lbx(i, j) = 0.5 ! 将弱凹三角形中心点记录
                                    ! 两个0.5合计1，表示要增加1条边（图d）
                                    if(j < 5)then
                                        ref_lbx(i, j + 1) = 0.5
                                    else
                                        ref_lbx(i, 1) = 0.5
                                    end if
                                end if
                            end do
                        else if(edges == 6)then    ! 六边形（和上述五边形情况类似）
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

                        ! sum(ref_lbx(i, 1:7))表示要增加的边数
                        ! 当增加后原边数+增加边数大于7，则细化整个多边形
                        if((edges == 5).and.(sum(ref_lbx(i, 1:7)) > 2.))then ! 5+2=7
                            do j = 1, 7, 1
                                if((ref_lbx(i, j) /= 0).and.(mrl_new(ngrwm(j, i)) == 1))then
                                    ref(ngrwm(j, i)) = 1
                                end if
                            end do
                        else if((edges == 6).and.(sum(ref_lbx(i, 1:7)) > 1.))then ! 6+1=7
                            do j = 1, 7, 1
                                if((ref_lbx(i, j) /= 0).and.(mrl_new(ngrwm(j, i)) == 1))then
                                    ref(ngrwm(j, i)) = 1
                                end if
                            end do
                        end if

                    else if(ref_lbx(i, 8) /= 0)then ! 当i点相邻三角形存在细化
                        if(edges == 5)then
                            ! 当i点相邻三角形至少有两个被细化 1+1+1+4+4
                            ! 则细化i点所有相邻三角形（这里有问题）
                            if(sum(mrl_new(ngrwm(1:5, i))) > 10)then 
                                do j = 1, 5, 1
                                    if(mrl_new(ngrwm(j, i)) == 1)then
                                        ref(ngrwm(j, i)) = 1
                                    end if
                                end do
                            end if
                        else if(edges == 6)then
                            ! 当i点相邻三角形有两个被细化 1+1+1+4+4+1
                            if(sum(mrl_new(ngrwm(1:6, i))) == 12)then
                                do j = 1, 3, 1
                                    ! 如果两个被细化三角形是相对的
                                    ! 细化一侧的三角形，将另一侧视作弱凹
                                    ! 弱凹处理将会减少1条边，若两边都视作弱凹
                                    ! 那么6-2=4，将会出现四边形（图f）
                                    if((mrl_new(ngrwm(j, i)) == 4).and.(mrl_new(ngrwm(j + 3, i)) == 4))then
                                        if((mrl_new(ngrwm(j + 1, i)) == 1).and.(mrl_new(ngrwm(j + 2, i)) == 1))then
                                            ref(ngrwm(j + 1, i)) = 1
                                            ref(ngrwm(j + 2, i)) = 1
                                        end if
                                    end if
                                end do
                            ! 当i点相邻三角形有一个被细化 1+1+1+4+1+1
                            else if(sum(mrl_new(ngrwm(1:6, i))) == 9)then
                                k = 0
                                do j = 1, 6, 1
                                    l = ngrwm(j, i)
                                    if(mrl_new(l) == 1)then
                                        ! 1+1+4=6，记录相邻三角形的相邻三角形有细化的情况
                                        if(sum(mrl_new(ngrmm(1:3, l))) == 6)then
                                            k = k + 1
                                        end if
                                    end if
                                end do
                                ! 由于本来存在一个已细化的相邻三角形
                                ! 2+2=4>3
                                ! 每个六边形外的相邻三角形要新增一条线
                                ! 因此6+2大于7时，细化所有相邻三角形（图g）
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
                else ! 存在三角形需要进行细化
                    iterA = .false. ! 需要循环进行循环B和循环C，下面过程与前面类似
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
        print*, "细化第二步完成"

        !--------------------------------------------------
        ! 2.2 储存第二步细化后新网格数据
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
        ! 3.1 寻找弱凹点
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
        print*, "开始细化弱凹点"

        nmp(iter) = nmp(iter - 1) + 2 * num_ref
        nwp(iter) = nwp(iter - 1) + num_ref

        !--------------------------------------------------
        ! 3.2 细化弱凹点
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

        print*, "第三步细化三角形总数为", refed(iter) * 2
        print*, "细化第三步完成"

        ref = 0
        iter = iter + 1

        !--------------------------------------------------
        ! 3.3 储存第三步细化后新网格数据
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
        ! 4.1 记录相邻只有一个细化的三角形
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

        print*, "开始细化相邻三角形经过初步细化数为1的三角形"
        print*, "iter =", iter, "num =", num_ref

        nmp(iter) = nmp(iter - 1) + 2 * num_ref
        nwp(iter) = nwp(iter - 1) + num_ref

        !--------------------------------------------------
        ! 4.2 细化相邻只有一个细化的三角形（一分为二）
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
        print*, "细化第四步完成"
        iter = iter + 1

        !--------------------------------------------------
        ! 4.3 储存第四步细化后新网格数据
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
        ! 5.1 计算并储存非重复w点
        !--------------------------------------------------
        print*, "开始制作多边形网格"
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
            if(isexist .eqv. .false.)then
                num_dbx = num_dbx + 1
                wp_f(num_dbx, 1) = wp_new(i, 1)
                wp_f(num_dbx, 2) = wp_new(i, 2)
            end if
        end do

        print*, "细化前共有", lbx_points, "个多边形网格"
        print*, "细化后共有", nwp(iter - 1), "个多边形网格"
        print*, "去除重复点后，还剩", num_dbx, "个多边形网格"

        !--------------------------------------------------
        ! 5.2 重新计算并储存ngrmw
        !--------------------------------------------------
        print*, "重新计算并储存ngrmw"
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
        print*, "计算完成"

        !--------------------------------------------------
        ! 5.3 重新计算并储存ngrwm
        !--------------------------------------------------
        print*, "重新计算并储存ngrwm"
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
        print*, "计算完成，开始排序"
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

        ! 对ngrwm进行排序（按顺序构成多边形）
        CALL GetSort(ngrwm_f, nwp(iter - 1), mp_new, sjx_points, num_dbx, nmp(iter - 1))
        !CALL GetSort(ngrwm_f(:, 1:num_dbx), mp_new(1:nmp(iter - 1), 1:2), num_dbx, nmp(iter - 1))

        print*, "排序完成"

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
        ! 5.4 去除已细化三角形（m点）
        !--------------------------------------------------
        print*, "开始去除已细化三角形"
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

        print*, "已被细化的三角形个数为", num_ref
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
        print*, "去除已被细化三角形前三角形个数为", nmp(iter - 1)
        print*, "去除已被细化三角形后三角形个数为", num_sjx

        !--------------------------------------------------
        ! 5.5 存储网格数据
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
        ! 6.1 计算三角形网格边长与角度
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

        print*, "三角形边长：", minval(length(sa_iter, 2:num_sjx, :)), maxval(length(sa_iter, 2:num_sjx, :))
        print*, "三角形角度：", minval(angle(sa_iter, 2:num_sjx, :)), maxval(angle(sa_iter, 2:num_sjx, :))
        !stop

        !--------------------------------------------------
        ! 6.2 计算弹性调整前三角形网格质量
        !--------------------------------------------------
        ! Extr_*** 指网格的最小角度与最大角度
        ! Eavg_*** 指网格的最小角度与最大角度平均值
        ! Savg_*** 指网格角度与正多边形角度的标准差
        ! less30   指三角形网格中小于30度角的数量
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
        ! 6.3 调整所有m点至三角形网格重心
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
        ! 6.4 计算弹性调整前多边形网格质量
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

        ! Extr_*** 指网格的最小角度与最大角度
        ! Eavg_*** 指网格的最小角度与最大角度平均值
        ! Savg_*** 指网格角度与正多边形角度的标准差
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

        print*, "五边形网格信息如下："
        print*, "网格数量", n_wbx
        print*, "最小角度", minval(angle_wbx(sa_iter, :, 1:5)), "最大角度", maxval(angle_wbx(sa_iter, :, 1:5))
        print*, "六边形网格信息如下："
        print*, "网格数量", n_lbx
        print*, "最小角度", minval(angle_lbx(sa_iter, :, 1:6)), "最大角度", maxval(angle_lbx(sa_iter, :, 1:6))
        print*, "七边形网格信息如下："
        print*, "网格数量", n_qbx
        print*, "最小角度", minval(angle_qbx(sa_iter, :, 1:7)), "最大角度", maxval(angle_qbx(sa_iter, :, 1:7))

        !stop

        !--------------------------------------------------
        ! 6.3 弹性调整&计算每次调整后的三角形网格质量
        !--------------------------------------------------
        print*, "开始弹性调整"
        allocate(MoveDis(num_dbx, 2))

        sa_iter = 1

        End_SpringAjustment = .false.

        do while(End_SpringAjustment .eqv. .false.)

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

            print*, "第", sa_iter, "次弹性调整完成，调整角度个数为", k

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
        ! 6.4 存储网格质量数据
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
        ! 6.5 存储最终网格数据
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
        if  (STATUS /= NF90_NOERR) then ! nf_noerr=0 表示没有错误
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

                if(icl .eqv. .true.)then
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


    ! 计算三角形abc的外心(垂直平分线交点)p
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
