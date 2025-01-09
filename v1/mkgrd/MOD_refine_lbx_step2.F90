! 多边形网格第2,3...次细化
! 1.计算上次最终输出网格文件包含关系
! 2.进行细化


module MOD_refine_lbx_step2
    USE consts_coms
    USE refine_vars
    USE NETCDF
    USE MOD_Get_Distance
    implicit none


contains

    subroutine refine_lbx_step2()

        implicit none

        integer :: i, j, k, l, m, n, x, y, z
        integer :: w1, w2, w3, w4, w5, w6, m1, m2, m3, m4      ! 三角形和多边形中心点序号
        integer :: row, col, edges, num_i, path
        integer :: num_ref                                     ! 每次细化三角形数
        integer :: num_all                                     ! 非结构网格包含结构网格总数
        integer :: num_null                                    ! 需要处理的空值数量
        integer :: refed(1000)                                 ! 记录每次细化时已经细化的三角形数
        integer :: iter                                        ! 网格细化次数
        integer :: sa_iter                                     ! 网格质量调整迭代次数
        integer :: num_sjx, num_dbx                            ! 细化后三角形、多边形数量
        integer :: nmp(1000), nwp(1000)                        ! 记录每次细化后的m，w点数量
        integer :: icl(3)                                      ! 判断三角形网格是否穿过180°经线
        !integer :: maxlc                                       ! 最大土地类型编号，17 or 24

        ! nc文件读写相关
        integer :: ncid, varid(10), ncvarid(13), upDimID, iunit
        integer :: spDImID, lpDimID, twDimID, thDimID, seDimID, sxDimID, dimID_sjx, dimID_lbx

        integer :: n_wbx, n_lbx, n_qbx                         ! 五边形、六边形、七边形数量
        integer :: sjx_points, lbx_points                      ! 三角形与多边形数量(上次细化后的文件)
        integer :: sjx_points_ori, lbx_points_ori              ! 三角形与多边形数量(未细化的文件)
        
        integer :: sum_sea,sum_land          ! 陆地和海洋网格数量
        integer :: maxid(1)
        integer :: numpatch                  ! 各非结构网格包含经纬度网格总数
        

        integer :: id(2, 3)                  ! 记录上次细化与未细化网格m点的相邻w点 
        integer :: set_dis                   ! 细化过渡行间隔
        integer :: step                      ! 表示第几次的细化迭代,从第二次开始

        real(r8) :: dx, dy, isinply, sjx(3, 2), newsjx(3, 2), dbx(7, 2)
        real(r8) :: fra                      ! 网格质量调整位移比例
        real(r8) :: rx, ry                   ! 网格质量调整位移分量

        real(r8), allocatable :: wp(:, :), mp(:, :)             ! 三角形、多边形网格中心点数据(上次细化)
        real(r8), allocatable :: wp_ori(:, :), mp_ori(:, :)     ! 三角形、多边形网格中心点数据(未细化)
        real(r8), allocatable :: wp_new(:, :), mp_new(:, :)     ! 三角形、多边形网格中心点数据(细化后)
        real(r8), allocatable :: wp_f(:, :), mp_f(:, :)         ! 三角形、多边形网格中心点数据(最终结果)
        real(r8), allocatable :: mp_f_tmp(:, :), wp_f_tmp(:, :)
        real(r8), allocatable :: mp_ii(:, :)                    ! 三角形网格与多边形网格包含的经纬度网格(上次细化)
        real(r8), allocatable :: mp_ii_ori(:, :)                ! 三角形网格与多边形网格包含的经纬度网格(未细化)
        real(r8), allocatable :: mp_ii_new(:, :)                ! 三角形网格与多边形网格包含的经纬度网格(细化后)

        integer, allocatable :: mp_id(:, :)                    ! 三角形网格包含经纬度网格数量与在mp_ii中的起始位置(上次细化)          
        integer, allocatable :: mp_id_ori(:, :)                ! 三角形网格包含经纬度网格数量与在mp_ii中的起始位置(未细化)   
        integer, allocatable :: mp_id_new(:, :)                ! 三角形网格包含经纬度网格数量与在mp_ii中的起始位置(细化后)   
        integer, allocatable :: mp_same(:, :)                  ! 记录未细化与细化后网格相同的m点 
        integer, allocatable :: mp_dis(:, :)                   ! 记录三角形网格m点与其他m点的距离

        integer, allocatable :: ngrmw(:, :), ngrwm(:, :)          ! m/w点相邻的w/m点索引(上次细化)
        integer, allocatable :: ngrmw_ori(:, :), ngrwm_ori(:, :)  ! m/w点相邻的w/m点索引(未细化)
        integer, allocatable :: ngrmw_new(:, :), ngrwm_new(:, :)  ! m/w点相邻的w/m点索引(细化后)
        integer, allocatable :: ngrwm_f(:, :), ngrmw_f(:, :)      ! m/w点相邻的w/m点索引(最终结果)
        integer, allocatable :: ngrmw_f_tmp(:, :), ngrwm_f_tmp(:, :)
        integer, allocatable :: ngrmm(:, :)                       ! m点相邻的m点索引(上次细化)
        integer, allocatable :: ngrmm_new(:, :)                   ! m点相邻的m点索引(细化后)

        integer, allocatable :: ref(:)          ! 三角形网格初始细化情况
        integer, allocatable :: ref_l(:, :)     ! 多边形网格细化情况
        integer, allocatable :: ref_th(:, :)    ! 初始三角形网格阈值细化情况
        integer, allocatable :: ref_tr(:, :)    ! 最终三角形网格阈值细化情况
        integer, allocatable :: ref_pl(:, :)    ! 记录多边形是因何种阈值细化
        integer, allocatable :: ref_jw(:, :)    ! 记录需要细化的经纬度网格
        real(r8), allocatable :: ref_lbx(:, :)
        real(r8), allocatable :: dismm(:,:),disww(:,:)

        integer, allocatable :: mrl(:)             ! 三角形网格细化程度(上次细化)
        integer, allocatable :: mrl_new(:)         ! 三角形网格细化程度(细化后)
        integer, allocatable :: mrl_f(:)           ! 三角形网格细化程度(最终)
        integer, allocatable :: ngr_mrl(:, :)      ! 三角形网格的相邻三角形网格点的mrl(上次细化)
        integer, allocatable :: ngr_mrl_new(:, :)  ! 三角形网格的相邻三角形网格点的mrl(细化后)
        integer, allocatable :: mp_ref(:)          ! 记录被细化三角形的索引


        integer, allocatable :: nl(:)              ! 非结构网格是否包含该种土地类型的标志
        real(r8), allocatable :: landtypes(:, :)    ! 经纬度网格土地类型
        integer, allocatable :: n_landtypes(:)     ! 三角形网格内包含的土地类型数量

        real(r8), allocatable :: lon_i(:), lat_i(:)               ! 经纬度网格中心点
        integer, allocatable :: seaorland(:, :)                   ! 判断经纬度网格为陆地或海洋

        ! 阈值数组
        real(r8),allocatable :: slope_max(:,:),slope_avg(:,:),p_slope(:,:)
        real(r8),allocatable :: area(:),area_fine_gridcell(:,:)
        real(r8),allocatable :: fraction_mainarea(:,:)
        real(r8),allocatable :: lai(:,:),p_lai(:,:)
        real(r8),allocatable :: k_s(:,:,:),p_k_s(:,:,:)
        real(r8),allocatable :: k_sl(:,:,:),p_k_sl(:,:,:)
        real(r8),allocatable :: tkdry(:,:,:),p_tkdry(:,:,:)
        real(r8),allocatable :: tksatf(:,:,:),p_tksatf(:,:,:)
        real(r8),allocatable :: tksatu(:,:,:),p_tksatu(:,:,:)


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

        character(LEN = 256) :: lndname, steps, nxpc
        logical,allocatable :: IsInRfArea(:)
        real(r8) :: tmpa(2),tmpb(2),tmpc(2),tmp_angle(3),tmp_length(3),tmp_angle7(7)
        !character(LEN = 20) :: p_name(6) = (/"GLONW", "GLATW", "GLONM", "GLATM", "itab_w%im", "itab_m%iw"/)
        character(LEN = 20),dimension(6) :: p_name

        logical :: ispart                ! 计算包含关系时，判断小网格是否被大网格完全包含
        logical :: isexist               ! 计算重复w点
        logical :: iter_is_end           ! 判断多重细化是否完成

        ! 防止细化交汇带出现冲突，一分为四
        logical :: iterA                  ! 当迭代B与迭代C同时一次性通过，迭代A通过
        logical :: iterB                  ! 从三角形网格进行判断
        logical :: iterC                  ! 从多边形网格进行判断
        logical :: End_SpringAjustment    ! 判断网格质量调整是否结束

        p_name = [character(len=20) :: "GLONW", "GLATW", "GLONM", "GLATM", "itab_w%im", "itab_m%iw"]
        set_dis = 5          ! 细化过渡行间隔，建议4-5

        step = 2             ! 表示第几次的细化迭代,默认从第二次开始

        iter_is_end = .false.

        write(nxpc, '(I4.4)') NXP

        do while(iter_is_end .eqv. .false.)

            if(step >= max_iter)then
                iter_is_end = .true.
            end if

            !-----------------------------------------------
            ! 1.读取原网格与细化网格文件
            !-----------------------------------------------
            write(steps, '(i2.2)')step - 1
            lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/gridfile_NXP" // trim(nxpc)  // "_" // trim(steps) // "_6.nc4"
            write(steps, '(i2.2)')step
            print*, lndname
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

            allocate(wp(lbx_points, 2))          ! 多边形网格中心点初始数据
            allocate(mp(sjx_points, 3))          ! 三角形网格中心点初始数据
            allocate(mp_id(sjx_points, 2))
            allocate(mp_id_new(sjx_points, 2))
            allocate(wp_new(lbx_points * 5, 2))
            allocate(mp_new(sjx_points * 5, 4))
            allocate(mp_same(sjx_points, 2))
            allocate(ngrwm(8, lbx_points))       ! wp的相邻mp点初始索引表
            allocate(ngrwm_new(7, lbx_points * 4)) ! wp的相邻mp点更新索引表
            allocate(ngrmw(4, sjx_points))       ! mp的相邻wp点初始索引表
            allocate(ngrmm(3, sjx_points))
            allocate(mrl(sjx_points))              ! 三角形网格细化程度
            allocate(mrl_new(sjx_points * 4))        ! 更新后数据
            allocate(ngr_mrl(3, sjx_points))        ! 三角形网格的相邻三角形网格点的mrl
            allocate(ngr_mrl_new(3, sjx_points * 4))  ! 更新后数据
            allocate(ngrmm_new(3, sjx_points * 4))    ! 更新后数据
            allocate(ngrmw_new(3, sjx_points * 5))

            wp = 0
            mp = 0
            mp_same = 0
            mp_id = 0
            mp_id_new = 0
            ngrmw = 0
            ngrwm = 0
            ngrmm = 0
            ngrwm(8, :) = 7
            ngrmw(4, :) = 3
            mrl = 1
            mrl_new = 1
            ngr_mrl = 1
            ngr_mrl_new = 1
            ngrmm_new = 1
            mp_new = 0.
            wp_new = 0.
            ngrwm_new = 1

            CALL CHECK(NF90_GET_VAR(ncid, varid(1), wp(:, 1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(2), wp(:, 2)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(3), mp(:, 1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(4), mp(:, 2)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(5), ngrwm(1:7, :)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(6), ngrmw(1:3, :)))
            CALL CHECK(NF90_CLOSE(ncid))

            lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/gridfile/gridfile_NXP" // trim(nxpc)  // ".nc4"
            print*, lndname
            CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid))
            CALL CHECK(NF90_INQ_VARID(ncid, p_name(1), varid(1)))
            CALL CHECK(NF90_INQ_VARID(ncid, p_name(2), varid(2)))
            CALL CHECK(NF90_INQ_VARID(ncid, p_name(3), varid(3)))
            CALL CHECK(NF90_INQ_VARID(ncid, p_name(4), varid(4)))
            CALL CHECK(NF90_INQ_VARID(ncid, p_name(5), varid(5)))
            CALL CHECK(NF90_INQ_VARID(ncid, p_name(6), varid(6)))
            CALL CHECK(NF90_INQ_DIMID(ncid, "sjx_points", dimID_sjx))
            CALL CHECK(NF90_INQ_DIMID(ncid, "lbx_points", dimID_lbx))
            CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lbx, len = lbx_points_ori))
            CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_sjx, len = sjx_points_ori))

            print*, "sjx_points_ori = ", sjx_points_ori
            print*, "lbx_points_ori = ", lbx_points_ori

            allocate(wp_ori(lbx_points_ori, 2))          ! 多边形网格中心点初始数据
            allocate(mp_ori(sjx_points_ori, 2))          ! 三角形网格中心点初始数据
            allocate(mp_id_ori(sjx_points_ori, 2))
            allocate(ngrwm_ori(8, lbx_points_ori))       ! wp的相邻mp点初始索引表
            allocate(ngrmw_ori(4, sjx_points_ori))       ! mp的相邻wp点初始索引表
            wp_ori = 0
            mp_ori = 0
            mp_id_ori = 0
            ngrmw_ori = 0
            ngrwm_ori = 0
            ngrmw_ori(4, :) = 3
            ngrwm_ori(8, :) = 7

            CALL CHECK(NF90_GET_VAR(ncid, varid(1), wp_ori(:, 1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(2), wp_ori(:, 2)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(3), mp_ori(:, 1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(4), mp_ori(:, 2)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(5), ngrwm_ori(1:7, :)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(6), ngrmw_ori(1:3, :)))
            CALL CHECK(NF90_CLOSE(ncid))

            ! 计算三角形网格与多边形网格
            do i = 1, sjx_points, 1
                do j = 1, 3, 1
                    if(ngrmw(j, i) == 1)then
                        ngrmw(4, i) = ngrmw(4, i) - 1
                    end if
                end do
            end do

            do i = 1, lbx_points, 1
                do j = 1, 7, 1
                    if(ngrwm(j, i) == 1)then
                        ngrwm(8, i) = ngrwm(8, i) - 1
                    end if
                end do
            end do

            do i = 1, sjx_points_ori, 1
                do j = 1, 3, 1
                    if(ngrmw_ori(j, i) == 1)then
                        ngrmw_ori(4, i) = ngrmw_ori(4, i) - 1
                    end if
                end do
            end do

            do i = 1, lbx_points_ori, 1
                do j = 1, 7, 1
                    if(ngrwm_ori(j, i) == 1)then
                        ngrwm_ori(8, i) = ngrwm_ori(8, i) - 1
                    end if
                end do
            end do

            do i = 1, sjx_points, 1
                do j = 1, sjx_points, 1
                    if(i /= j)then
                        k = IsNgrmm(ngrmw(1:3, i), ngrmw(1:3, j))
                        if(k /= 0)then
                            ngrmm(k, i) = j
                        end if
                    end if
                end do
            end do

            allocate(mp_dis(sjx_points, sjx_points))
            mp_dis = 0

            CALL GetTriangleDis(ngrwm, lbx_points, ngrmw, sjx_points, mp_dis)
            !stop
            !print*,maxval(ngrmw(4,:)),minval(ngrmw(4,:))
            !print*,maxval(ngrwm(8,:)),minval(ngrwm(8,:))
            !print*,maxval(ngrmw_ori(4,:)),minval(ngrmw_ori(4,:))
            !print*,maxval(ngrwm_ori(8,:)),minval(ngrwm_ori(8,:))
            !print*,minval(ngrmw(1:3,:)),maxval(ngrmw(1:3,:))

            lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/contain/initial/mp_id.bin"
            !lndname = "/stu01/fanhw21/olam/Ouhe/output/MKSRFDATA/grid/sjx/nxp"//trim(nxp)//"/mp_id.bin"
            print*, lndname
            OPEN(iunit, file = trim(lndname), form = 'unformatted', status = 'old')
            READ(iunit) mp_id_ori(:, 1:2)
            close(iunit)

            num_all = INT(sum(mp_id_ori(:, 1)))
            allocate(mp_ii_ori(num_all, 4))

            lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/contain/initial/mp_ii.bin"
            !lndname = "/stu01/fanhw21/olam/Ouhe/output/MKSRFDATA/grid/sjx/nxp"//trim(nxp)//"/mp_ii.bin"
            print*, lndname
            OPEN(iunit, file = trim(lndname), form = 'unformatted', status = 'old')
            READ(iunit) mp_ii_ori
            close(iunit)

            !-----------------------------------------------
            ! 2.计算网格包含关系
            !-----------------------------------------------
            print*, "开始定义初始数据网格......"
            print*, ""
            dx = 360. / nlons_source      ! = 1./120.
            dy = 180. / nlats_source      ! = 1./120.
            sum_sea = 0                      ! 海洋非结构网格数
            sum_land = 0                     ! 陆地经纬度网格数
            allocate(lon_i(nlons_source))
            allocate(lat_i(nlats_source))

            do i = 1, nlons_source, 1
                lon_i(i) = -180. + (2 * i - 1) * dx / 2.
            end do

            do i = 1, nlats_source, 1
                lat_i(i) = 90. - (2 * i - 1) * dy / 2.
            end do

            allocate(area_fine_gridcell(nlons_source, nlats_source))
            area_fine_gridcell(:, :) = 0.
            print*, "开始计算经纬度网格面积......"
            call cellarea(area_fine_gridcell)
            print*, "经纬度网格面积计算完成"

            allocate(landtypes(nlons_source, nlats_source))
            allocate(seaorland(nlons_source, nlats_source))
            landtypes = 0.
            seaorland = 0

            if(lcs == "igbp")then
               lndname = trim(source_dir) // '/landtype_igbp_update.nc'
               !maxlc = 17
            else
               lndname = trim(source_dir) // '/landtype_usgs_update.nc'
               !maxlc = 24
            end if
            print*,lndname
            CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
            CALL CHECK(NF90_INQ_VARID(ncid, "landtype", varid(1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), landtypes))
            CALL CHECK(NF90_CLOSE(ncid))
            print*,"landtypes",minval(landtypes),maxval(landtypes)

            !!$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1)&
            !!$OMP PRIVATE(i,j)
            do i = 1, nlons_source, 1
                do j = 1, nlats_source, 1
                    if(landtypes(i, j) /= 0.)then
                        seaorland(i, j) = 1
                        sum_land = sum_land + 1
                    else
                        sum_sea = sum_sea + 1
                    end if
                end do
            end do
            !!$OMP END PARALLEL DO

            print*, "海洋网格个数为", sum_sea, "，陆地网格个数为", sum_land
            print*, ""

            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1)&
            !$OMP PRIVATE(i,j,k,l,id)
            do i = 1, sjx_points, 1
                if(ngrmw(4, i) == 3)then
                    do k = 1, 3, 1
                        id(1, k) = ngrmw(k, i)
                    end do
                    do j = 1, sjx_points_ori, 1
                        if(ngrmw_ori(4, j) == 3)then
                            do l = 1, 3, 1
                                id(2, l) = ngrmw_ori(l, j)
                            end do

                            if((wp(id(1, 1), 1) == wp_ori(id(2, 1), 1)).and.(wp(id(1, 1), 2) == wp_ori(id(2, 1), 2)).and.&
                                    (wp(id(1, 2), 1) == wp_ori(id(2, 2), 1)).and.(wp(id(1, 2), 2) == wp_ori(id(2, 2), 2)).and.&
                                    (wp(id(1, 3), 1) == wp_ori(id(2, 3), 1)).and.(wp(id(1, 3), 2) == wp_ori(id(2, 3), 2)))then
                                mp_same(i, 1) = 1
                                mp_same(i, 2) = j

                            end if
                        end if
                    end do
                end if
            end do
            !$OMP END PARALLEL DO

            do i = 1, sjx_points, 1
                if(mp_same(i, 1) /= 0)then
                    if(mp_same(i, 2) == 0)then
                        print*, "error"
                        stop
                    end if
                end if
            end do

            allocate(ref_jw(nlons_source, nlats_source))
            ref_jw = 0

            do i = 1, sjx_points, 1
                if((mp_same(i, 1) == 0).and.(ngrmw(4, i) == 3))then
                    tmpa = wp(ngrmw(1, i), 1:2)
                    tmpb = wp(ngrmw(2, i), 1:2)
                    tmpc = wp(ngrmw(3, i), 1:2)
                    CALL GetRefineJW(tmpa, tmpb, tmpc, ref_jw)
                    !CALL GetRefineJW(wp(ngrmw(1, i), 1:2), wp(ngrmw(2, i), 1:2), wp(ngrmw(3, i), 1:2), ref_jw)
                end if
            end do

            print*, sum(ref_jw)
            !stop
            ! --------------------------------------------------------
            ! 获取非结构网格中经纬度网格的数量及所占面积比例
            ! 1.计算数组大小
            ! 2.分配内存
            ! 3.计算包含关系
            ! -------------------------------------------------------

            ! 首先计算数组大小
            print*, "开始计算三角形网格与经纬度网格包含关系数组大小......"
            print*, ""
            num_i = 0
            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
            !$OMP PRIVATE(i,j,k,l,sjx,isinply,ispart)
            do i = 1, nlons_source, 1
                num_i = num_i + 1
                !print*,num_i
                do j = 1, nlats_source, 1         ! 循环遍历初始网格单元
                    ispart = .false.
                    if((seaorland(i, j) == 0).or.(ref_jw(i, j) == 0))then
                        cycle
                    end if          ! 只遍历海洋网格
                    do k = 1, sjx_points, 1           ! 循环遍历三角形网格
                        if(ispart .eqv. .true.)then
                            exit
                        end if
                        if(mp_same(k, 1) == 1)then
                            cycle
                        end if
                        if(ngrmw(4, k) == 3)then          ! 判断m点邻域是否有三个被记录的w点
                            sjx = 0
                            do l = 1, 3, 1
                                sjx(l, 1:2) = wp(ngrmw(l, k), 1:2)
                            end do

                            isinply = IsInUstrGrid(sjx, lon_i(i), lat_i(j))

                            if(isinply == 1.)then
                                ispart = .true.
                            end if

                            if(isinply > 0.)then            ! 记录每个非结构网格包含经纬度网格数量
                                mp_id(k, 1) = mp_id(k, 1) + 1
                            end if
                        end if
                    end do
                end do
            end do
            !$OMP END PARALLEL DO

            print*, sum(mp_id(:, 1))

            ! 分配数组内存
            do i = 1, sjx_points, 1
                if(mp_same(i, 1) == 1)then
                    j = mp_same(i, 2)
                    mp_id(i, 1) = mp_id_ori(j, 1)
                end if
            end do

            numpatch = INT(sum(mp_id(:, 1))) * 2
            print*, "非结构网格包含经纬度网格总数为", numpatch / 2.
            allocate(mp_ii(numpatch, 4))
            mp_ii = 0.

            ! mp_id(:,1)记录mp或wp数组在mp_ii数组中的位置
            mp_id(1, 2) = 1

            do i = 2, sjx_points, 1
                mp_id(i, 2) = mp_id(i - 1, 2) + mp_id(i - 1, 1) * 2
            end do

            mp_id(:, 1) = 0
            !stop

            num_i = 0
            ! 计算包含关系
            print*, "开始计算三角形网格与经纬度网格包含关系......"
            print*, ""
            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
            !$OMP PRIVATE(i,j,k,l,sjx,isinply,ispart,path)
            do i = 1, nlons_source, 1
                num_i = num_i + 1
                !print*, num_i
                do j = 1, nlats_source, 1         ! 循环遍历初始网格单元
                    ispart = .false.
                    if((seaorland(i, j) == 0).or.(ref_jw(i, j) == 0))then
                        cycle
                    end if          ! 只遍历海洋网格
                    do k = 1, sjx_points, 1           ! 循环遍历三角形网格
                        if(ispart .eqv. .true.)then
                            exit
                        end if
                        if(mp_same(k, 1) == 1)then
                            cycle
                        end if
                        if(ngrmw(4, k) == 3)then          ! 判断m点邻域是否有三个被记录的w点
                            sjx = 0
                            do l = 1, 3, 1
                                sjx(l, 1:2) = wp(ngrmw(l, k), 1:2)
                            end do

                            isinply = IsInUstrGrid(sjx, lon_i(i), lat_i(j))
                            if(isinply == 1.)then            ! 完全包含
                                mp_id(k, 1) = mp_id(k, 1) + 1
                                mp(k, 3) = mp(k, 3) + area_fine_gridcell(i, j)

                                path = mp_id(k, 2) + mp_id(k, 1) - 1
                                mp_ii(path, 1) = i
                                mp_ii(path, 2) = j
                                mp_ii(path, 3) = isinply
                                mp_ii(path, 4) = area_fine_gridcell(i, j)
                                ispart = .true.               ! 若完全包含，则跳出循环
                            else if(isinply > 0.)then       ! 相交（部分包含）
                                mp_id(k, 1) = mp_id(k, 1) + 1
                                mp(k, 3) = mp(k, 3) + isinply * area_fine_gridcell(i, j)

                                path = mp_id(k, 2) + mp_id(k, 1) - 1
                                mp_ii(path, 1) = i
                                mp_ii(path, 2) = j
                                mp_ii(path, 3) = isinply
                                mp_ii(path, 4) = area_fine_gridcell(i, j) * isinply
                            end if
                        end if
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
            !stop
            print*, sum(mp_id(:, 1))

            do i = 1, sjx_points, 1
                num_null = 0
                if(mp_same(i, 1) == 0)then
                    do j = mp_id(i, 2), mp_id(i, 2) + mp_id(i, 1) - 1, 1
                        if(j > numpatch)then
                            print*, "j is bigger than numpatch!!!"
                            stop
                        else if(j < 1)then
                            print*, "j is smaller than 1!!!"
                            stop
                        end if
                        if(mp_ii(j, 3) == 0.)then
                            num_null = num_null + 1
                        end if
                    end do
                else
                    k = mp_same(i, 2)
                    if(k > sjx_points_ori)then
                        print*, "k is big than sjx_points_ori!!!"
                        stop
                    else if(k < 1)then
                        print*, "k is smaller than 1!!!"
                        stop
                    end if
                    do j = mp_id_ori(k, 2), mp_id_ori(k, 2) + mp_id_ori(k, 1) - 1, 1
                        if(j > num_all)then
                            print*, "j is big than num_all!!!"
                            stop
                        else if(j < 1)then
                            print*, "j is smaller than 1 in ori!!!"
                            stop
                        end if
                        if(mp_ii_ori(j, 3) == 0.)then
                            num_null = num_null + 1
                        end if
                    end do
                end if
                mp_id(i, 1) = mp_id(i, 1) - num_null
            end do

            mp_id_new = 0
            mp_id_new(:, 1) = mp_id(:, 1)
            mp_id_new(1, 2) = 1

            do i = 2, sjx_points, 1
                mp_id_new(i, 2) = mp_id_new(i - 1, 2) + mp_id_new(i - 1, 1)
            end do

            numpatch = INT(sum(mp_id_new(:, 1)))
            print*, numpatch
            allocate(mp_ii_new(numpatch, 4))
            mp_ii_new = 0

            do i = 1, sjx_points, 1
                do j = 1, 4, 1
                    if(mp_id_new(i, 1) /= 0)then
                        if(mp_same(i, 1) == 0)then
                            mp_ii_new(mp_id_new(i, 2):mp_id_new(i, 2) + mp_id_new(i, 1) - 1, j) = &
                                    mp_ii(mp_id(i, 2):mp_id(i, 2) + mp_id_new(i, 1) - 1, j)
                        else
                            mp_ii_new(mp_id_new(i, 2):mp_id_new(i, 2) + mp_id_new(i, 1) - 1, j) = &
                                    mp_ii_ori(mp_id_ori(i, 2):mp_id_ori(i, 2) + mp_id(i, 1) - 1, j)
                        end if
                    end if
                end do
            end do

            !stop

            !-----------------------------------------------
            ! 3.读取并计算阈值
            !-----------------------------------------------

            lndname = trim(source_dir) // 'slope_max.nc'
            print*, lndname
            allocate(slope_max(nlons_source, nlats_source))          ! 经纬度网格八邻域上坡度最大值
            slope_max = 0.
            CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
            CALL CHECK(NF90_INQ_VARID(ncid, "slope_max", varid(1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), slope_max))
            CALL CHECK(NF90_CLOSE(ncid))
            print*, "slope", minval(slope_max), maxval(slope_max)

            lndname = trim(source_dir) // 'slope_avg.nc'
            print*, lndname
            allocate(slope_avg(nlons_source, nlats_source))          ! 经纬度网格八邻域上坡度均值
            slope_avg = 0.
            CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
            CALL CHECK(NF90_INQ_VARID(ncid, "slope_avg", varid(1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), slope_avg))
            CALL CHECK(NF90_CLOSE(ncid))
            print*, "slope_avg", minval(slope_avg), maxval(slope_avg)

            lndname = trim(source_dir) // 'lai.nc'
            print*, lndname
            allocate(lai(nlons_source, nlats_source))
            lai = 0.
            CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
            CALL CHECK(NF90_INQ_VARID(ncid, "lai", varid(1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), lai))
            CALL CHECK(NF90_CLOSE(ncid))
            print*, "lai", minval(lai), maxval(lai)

            lndname = trim(source_dir) // 'k_s.nc'
            print*, lndname
            allocate(k_s(nlons_source, nlats_source, 2))
            k_s = 0.
            CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
            CALL CHECK(NF90_INQ_VARID(ncid, "k_s_l1", varid(1)))
            CALL CHECK(NF90_INQ_VARID(ncid, "k_s_l2", varid(2)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), k_s(:, :, 1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(2), k_s(:, :, 2)))
            CALL CHECK(NF90_CLOSE(ncid))
            print*, "k_s_l1", minval(k_s(:, :, 1)), maxval(k_s(:, :, 1))
            print*, "k_s_l2", minval(k_s(:, :, 2)), maxval(k_s(:, :, 2))

            lndname = trim(source_dir) // 'k_solids.nc'
            print*, lndname
            allocate(k_sl(nlons_source, nlats_source, 2))
            k_sl = 0.
            CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
            CALL CHECK(NF90_INQ_VARID(ncid, "k_solids_l1", varid(1)))
            CALL CHECK(NF90_INQ_VARID(ncid, "k_solids_l2", varid(2)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), k_sl(:, :, 1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(2), k_sl(:, :, 2)))
            CALL CHECK(NF90_CLOSE(ncid))
            print*, "k_solids_l1", minval(k_sl(:, :, 1)), maxval(k_sl(:, :, 1))
            print*, "k_solids_l2", minval(k_sl(:, :, 2)), maxval(k_sl(:, :, 2))

            lndname = trim(source_dir) // 'tkdry.nc'
            print*, lndname
            allocate(tkdry(nlons_source, nlats_source, 2))
            tkdry = 0.
            CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
            CALL CHECK(NF90_INQ_VARID(ncid, "tkdry_l1", varid(1)))
            CALL CHECK(NF90_INQ_VARID(ncid, "tkdry_l2", varid(2)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), tkdry(:, :, 1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(2), tkdry(:, :, 2)))
            CALL CHECK(NF90_CLOSE(ncid))
            print*, "tkdry_l1", minval(tkdry(:, :, 1)), maxval(tkdry(:, :, 1))
            print*, "tkdry_l2", minval(tkdry(:, :, 2)), maxval(tkdry(:, :, 2))

            lndname = trim(source_dir) // 'tksatf.nc'
            print*, lndname
            allocate(tksatf(nlons_source, nlats_source, 2))
            tksatf = 0.
            CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
            CALL CHECK(NF90_INQ_VARID(ncid, "tksatf_l1", varid(1)))
            CALL CHECK(NF90_INQ_VARID(ncid, "tksatf_l2", varid(2)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), tksatf(:, :, 1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(2), tksatf(:, :, 2)))
            CALL CHECK(NF90_CLOSE(ncid))
            print*, "tksatf_l1", minval(tksatf(:, :, 1)), maxval(tksatf(:, :, 1))
            print*, "tksatf_l2", minval(tksatf(:, :, 2)), maxval(tksatf(:, :, 2))

            lndname = trim(source_dir) // 'tksatu.nc'
            print*, lndname
            allocate(tksatu(nlons_source, nlats_source, 2))
            tksatu = 0.
            CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
            CALL CHECK(NF90_INQ_VARID(ncid, "tksatu_l1", varid(1)))
            CALL CHECK(NF90_INQ_VARID(ncid, "tksatu_l2", varid(2)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), tksatu(:, :, 1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(2), tksatu(:, :, 2)))
            CALL CHECK(NF90_CLOSE(ncid))
            print*, "tksatu_l1", minval(tksatu(:, :, 1)), maxval(tksatu(:, :, 1))
            print*, "tksatu_l2", minval(tksatu(:, :, 2)), maxval(tksatu(:, :, 2))

            !print*,"开始计算经纬度网格面积......"
            !call cellarea(area_fine_gridcell)
            !print*,"经纬度网格面积计算完成"

            !lndname = "/tera01/fanhw21/olam/Ouhe/output/MAKEGRID/tool/output/area.nc4"
            !print*,lndname
            !CALL CHECK(NF90_CREATE(trim(lndname), nf90_clobber, ncID))
            !CALL CHECK(NF90_DEF_DIM(ncID, "lon_points", nlons_source, spDimID))
            !CALL CHECK(NF90_DEF_DIM(ncID, "lat_points", nlats_source, lpDimID))
            !CALL CHECK(NF90_DEF_VAR(ncID, "area"    , NF90_float, (/ spDimID, lpDimID /), ncVarID(1)))
            !CALL CHECK(NF90_ENDDEF(ncID))
            !CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1), area_fine_gridcell))
            !CALL CHECK(NF90_CLOSE(ncID))

            !print*,area_fine_gridcell
            !stop

            allocate(nl(0:maxlc))                           ! 非结构网格是否包含该种土地类型的标志
            allocate(area(0:maxlc))                         ! 非结构网格内各土地类型面积
            allocate(fraction_mainarea(sjx_points, 2))    ! 非结构网格内主要土地类型编号与比例
            allocate(p_slope(sjx_points, 3))            ! 非结构网格内经纬度网格坡度最大值
            allocate(p_lai(sjx_points, 3))                ! 非结构网格内lai均值、最大值、网格数
            allocate(p_k_s(sjx_points, 3, 2))
            allocate(p_k_sl(sjx_points, 3, 2))
            allocate(p_tkdry(sjx_points, 3, 2))
            allocate(p_tksatf(sjx_points, 3, 2))
            allocate(p_tksatu(sjx_points, 3, 2))
            allocate(n_landtypes(sjx_points))
            allocate(ref_tr(sjx_points * 20, 16))
            fraction_mainarea = 0.
            p_slope = 0.
            n_landtypes = 0
            ref_jw = 0
            p_lai = 0.
            p_k_s = 0.
            p_k_sl = 0.
            p_tkdry = 0.
            p_tksatf = 0.
            p_tksatu = 0.
            ref_tr = 0

            !--------------------------------------------------------------------------
            ! 4.初步计算未迭代网格的阈值文件
            !--------------------------------------------------------------------------

            print*, "开始计算阈值文件"
            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC) &
            !$OMP PRIVATE(i,j,row,col,nl,L,area,maxid)
            do i = 1, sjx_points, 1
                nl = 0
                area = 0.
                maxid = 0
                if(mp_id_new(i, 1) == 0)then
                    cycle
                end if
                do j = 0, mp_id_new(i, 1) - 1, 1
                    col = mp_ii_new(mp_id_new(i, 2) + j, 1)
                    row = mp_ii_new(mp_id_new(i, 2) + j, 2)
                    if(col == 0. .or. row == 0.)then
                        cycle
                    end if
                    L = int(landtypes(col, row))

                    if((L /= 0).and.(L /= maxlc))then

                        nl(L) = 1
                        area(L) = area(L) + mp_ii_new(mp_id_new(i, 2) + j, 4)

                        p_lai(i, 1) = p_lai(i, 1) + lai(col, row)
                        p_lai(i, 3) = p_lai(i, 3) + 1

                        p_k_s(i, 1, :) = p_k_s(i, 1, :) + k_s(col, row, :)
                        p_k_s(i, 3, :) = p_k_s(i, 3, :) + 1

                        p_k_sl(i, 1, :) = p_k_sl(i, 1, :) + k_sl(col, row, :)
                        p_k_sl(i, 3, :) = p_k_sl(i, 3, :) + 1

                        p_tkdry(i, 1, :) = p_tkdry(i, 1, :) + tkdry(col, row, :)
                        p_tkdry(i, 3, :) = p_tkdry(i, 3, :) + 1

                        p_tksatf(i, 1, :) = p_tksatf(i, 1, :) + tksatf(col, row, :)
                        p_tksatf(i, 3, :) = p_tksatf(i, 3, :) + 1

                        p_tksatu(i, 1, :) = p_tksatu(i, 1, :) + tksatu(col, row, :)
                        p_tksatu(i, 3, :) = p_tksatu(i, 3, :) + 1

                        p_slope(i, 1) = p_slope(i, 1) + slope_avg(col, row)
                        p_slope(i, 3) = p_slope(i, 3) + 1

                    end if

                end do

                n_landtypes(i) = INT(sum(nl))
                maxid = maxloc(area) - 1
                fraction_mainarea(i, 1) = maxid(1)
                fraction_mainarea(i, 2) = area(maxid(1)) / mp(i, 3)

                if(fraction_mainarea(i, 2) > 1.)then
                    fraction_mainarea(i, 2) = 1.
                end if

                if(p_lai(i, 3) /= 0.)then
                    p_lai(i, 1) = p_lai(i, 1) / p_lai(i, 3)
                end if

                if(p_k_s(i, 3, 1) /= 0.)then
                    p_k_s(i, 1, 1) = p_k_s(i, 1, 1) / p_k_s(i, 3, 1)
                end if

                if(p_k_s(i, 3, 2) /= 0.)then
                    p_k_s(i, 1, 2) = p_k_s(i, 1, 2) / p_k_s(i, 3, 2)
                end if

                if(p_k_sl(i, 3, 1) /= 0.)then
                    p_k_sl(i, 1, 1) = p_k_sl(i, 1, 1) / p_k_sl(i, 3, 1)
                end if

                if(p_k_sl(i, 3, 2) /= 0.)then
                    p_k_sl(i, 1, 2) = p_k_sl(i, 1, 2) / p_k_sl(i, 3, 2)
                end if

                if(p_tkdry(i, 3, 1) /= 0.)then
                    p_tkdry(i, 1, 1) = p_tkdry(i, 1, 1) / p_tkdry(i, 3, 1)
                end if

                if(p_tkdry(i, 3, 2) /= 0.)then
                    p_tkdry(i, 1, 2) = p_tkdry(i, 1, 2) / p_tkdry(i, 3, 2)
                end if

                if(p_tksatf(i, 3, 1) /= 0.)then
                    p_tksatf(i, 1, 1) = p_tksatf(i, 1, 1) / p_tksatf(i, 3, 1)
                end if

                if(p_tksatf(i, 3, 2) /= 0.)then
                    p_tksatf(i, 1, 2) = p_tksatf(i, 1, 2) / p_tksatf(i, 3, 2)
                end if

                if(p_tksatu(i, 3, 1) /= 0.)then
                    p_tksatu(i, 1, 1) = p_tksatu(i, 1, 1) / p_tksatu(i, 3, 1)
                end if

                if(p_tksatu(i, 3, 2) /= 0.)then
                    p_tksatu(i, 1, 2) = p_tksatu(i, 1, 2) / p_tksatu(i, 3, 2)
                end if

                if(p_slope(i, 3) /= 0.)then
                    p_slope(i, 1) = p_slope(i, 1) / p_slope(i, 3)
                end if

            end do
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
            !$OMP PRIVATE(i,j,L,row,col)
            do i = 1, sjx_points, 1
                if(mp_id_new(i, 1) == 0)then
                    cycle
                end if
                do j = 0, mp_id_new(i, 1) - 1, 1
                    col = mp_ii_new(mp_id_new(i, 2) + j, 1)
                    row = mp_ii_new(mp_id_new(i, 2) + j, 2)
                    if(col == 0. .or. row == 0.)then
                        cycle
                    end if
                    L = int(landtypes(col, row))

                    if((L /= 0).and.(L /= maxlc))then

                        p_lai(i, 2) = p_lai(i, 2) + (lai(col, row) - p_lai(i, 1)) * (lai(col, row) - p_lai(i, 1))
                        p_k_s(i, 2, 1) = p_k_s(i, 2, 1) + (k_s(col, row, 1) - p_k_s(i, 1, 1)) * (k_s(col, row, 1) - p_k_s(i, 1, 1))
                        p_k_s(i, 2, 2) = p_k_s(i, 2, 2) + (k_s(col, row, 2) - p_k_s(i, 1, 2)) * (k_s(col, row, 2) - p_k_s(i, 1, 2))
                        p_k_sl(i, 2, 1) = p_k_sl(i, 2, 1) + (k_sl(col, row, 1) - p_k_sl(i, 1, 1)) * (k_sl(col, row, 1) - p_k_sl(i, 1, 1))
                        p_k_sl(i, 2, 2) = p_k_sl(i, 2, 2) + (k_sl(col, row, 2) - p_k_sl(i, 1, 2)) * (k_sl(col, row, 2) - p_k_sl(i, 1, 2))
                        p_tkdry(i, 2, 1) = p_tkdry(i, 2, 1) + (tkdry(col, row, 1) - p_tkdry(i, 1, 1)) * (tkdry(col, row, 1) - p_tkdry(i, 1, 1))
                        p_tkdry(i, 2, 2) = p_tkdry(i, 2, 2) + (tkdry(col, row, 2) - p_tkdry(i, 1, 2)) * (tkdry(col, row, 2) - p_tkdry(i, 1, 2))
                        p_tksatf(i, 2, 1) = p_tksatf(i, 2, 1) + (tksatf(col, row, 1) - p_tksatf(i, 1, 1)) * (tksatf(col, row, 1) - p_tksatf(i, 1, 1))
                        p_tksatf(i, 2, 2) = p_tksatf(i, 2, 2) + (tksatf(col, row, 2) - p_tksatf(i, 1, 2)) * (tksatf(col, row, 2) - p_tksatf(i, 1, 2))
                        p_tksatu(i, 2, 1) = p_tksatu(i, 2, 1) + (tksatu(col, row, 1) - p_tksatu(i, 1, 1)) * (tksatu(col, row, 1) - p_tksatu(i, 1, 1))
                        p_tksatu(i, 2, 2) = p_tksatu(i, 2, 2) + (tksatu(col, row, 2) - p_tksatu(i, 1, 2)) * (tksatu(col, row, 2) - p_tksatu(i, 1, 2))
                        p_slope(i, 2) = p_slope(i, 2) + (slope_avg(col, row) - p_slope(i, 1)) * (slope_avg(col, row) - p_slope(i, 1))

                    end if
                end do

                p_lai(i, 2) = sqrt(p_lai(i, 2) / p_lai(i, 3))
                p_k_s(i, 2, 1) = sqrt(p_k_s(i, 2, 1) / p_k_s(i, 3, 1))
                p_k_s(i, 2, 2) = sqrt(p_k_s(i, 2, 2) / p_k_s(i, 3, 2))
                p_k_sl(i, 2, 1) = sqrt(p_k_sl(i, 2, 1) / p_k_sl(i, 3, 1))
                p_k_sl(i, 2, 2) = sqrt(p_k_sl(i, 2, 2) / p_k_sl(i, 3, 2))
                p_tkdry(i, 2, 1) = sqrt(p_tkdry(i, 2, 1) / p_tkdry(i, 3, 1))
                p_tkdry(i, 2, 2) = sqrt(p_tkdry(i, 2, 2) / p_tkdry(i, 3, 2))
                p_tksatf(i, 2, 1) = sqrt(p_tksatf(i, 2, 1) / p_tksatf(i, 3, 1))
                p_tksatf(i, 2, 2) = sqrt(p_tksatf(i, 2, 2) / p_tksatf(i, 3, 2))
                p_tksatu(i, 2, 1) = sqrt(p_tksatu(i, 2, 1) / p_tksatu(i, 3, 1))
                p_tksatu(i, 2, 2) = sqrt(p_tksatu(i, 2, 2) / p_tksatu(i, 3, 2))
                p_slope(i, 2) = sqrt(p_slope(i, 2) / p_slope(i, 3))

            end do
            !$OMP END PARALLEL DO

            print*, "阈值文件计算完成"

            print*,"存储网格阈值文件"
   
            write(steps,'(i2.2)')step
            lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/threshold/threshold_NXP" // trim(nxpc) // "_" // trim(adjustl(steps)) // ".nc4"
            print*,lndname
            CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber,nf90_netcdf4), ncID))
            CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", sjx_points, spDimID))
            CALL CHECK(NF90_DEF_DIM(ncID, "dima"        , 2       , twDimID))
            CALL CHECK(NF90_DEF_VAR(ncID, "num_landtypes"    , NF90_INT  , (/ spDimID /)         , VarID(1)))
            CALL CHECK(NF90_DEF_VAR(ncID, "fraction_mainarea", NF90_float, (/ spDimID, twDimID /), VarID(2)))
            CALL CHECK(NF90_DEF_VAR(ncID, "p_slope"          , NF90_float, (/ spDimID, twDimID /), VarID(3)))
            CALL CHECK(NF90_DEF_VAR(ncID, "p_lai"            , NF90_float, (/ spDimID, twDimID /), VarID(4)))
            CALL CHECK(NF90_DEF_VAR(ncID, "p_k_s"            , NF90_float, (/ spDimID, twDimID, twDimID /), VarID(5)))
            CALL CHECK(NF90_DEF_VAR(ncID, "p_k_sl"           , NF90_float, (/ spDimID, twDimID, twDimID /), VarID(6)))
            CALL CHECK(NF90_DEF_VAR(ncID, "p_tkdry"          , NF90_float, (/ spDimID, twDimID, twDimID /), VarID(7)))
            CALL CHECK(NF90_DEF_VAR(ncID, "p_tksatf"         , NF90_float, (/ spDimID, twDimID, twDimID /), VarID(8)))
            CALL CHECK(NF90_DEF_VAR(ncID, "p_tksatu"         , NF90_float, (/ spDimID, twDimID, twDimID /), VarID(9)))
            CALL CHECK(NF90_ENDDEF(ncID))
            CALL CHECK(NF90_PUT_VAR(ncID, VarID(1), n_landtypes))
            CALL CHECK(NF90_PUT_VAR(ncID, VarID(2), fraction_mainarea))
            CALL CHECK(NF90_PUT_VAR(ncID, VarID(3), p_slope(:,1:2)))
            CALL CHECK(NF90_PUT_VAR(ncID, VarID(4), p_lai(:,1:2)))
            CALL CHECK(NF90_PUT_VAR(ncID, VarID(5), p_k_s(:,1:2,:)))
            CALL CHECK(NF90_PUT_VAR(ncID, VarID(6), p_k_sl(:,1:2,:)))
            CALL CHECK(NF90_PUT_VAR(ncID, VarID(7), p_tkdry(:,1:2,:)))
            CALL CHECK(NF90_PUT_VAR(ncID, VarID(8), p_tksatf(:,1:2,:)))
            CALL CHECK(NF90_PUT_VAR(ncID, VarID(9), p_tksatu(:,1:2,:)))
            CALL CHECK(NF90_CLOSE(ncID))


            !-----------------------------------------------
            ! 5.判断初步迭代网格
            !-----------------------------------------------
            print*, "判断初步迭代网格"

            allocate(ref(sjx_points))           ! 三角形网格初始细化情况
            allocate(ref_th(sjx_points, 16))     ! 初始三角形网格阈值细化情况
            allocate(ref_lbx(lbx_points, 8))
            allocate(ref_l(lbx_points, 7))       ! 多边形网格细化情况

            ref = 0
            ref_l = 0
            ref_th = 0
            ref_lbx = 0.

            allocate(IsInRfArea(sjx_points))
            IsInRfArea = .false.
            CALL IsInRefineArea(IsInRfArea,sjx_points,lbx_points,ngrmw,wp)

            do i = 1, sjx_points, 1
                !if((wp(ngrmw(i,1),2)>85.).or.(wp(ngrmw(i,2),2)>85.).or.(wp(ngrmw(i,3),2)>85.))then
                !   cycle
                !end if
                if(IsInRfArea(i) .eqv. .false.)then
                  cycle
                end if

                if((fraction_mainarea(i, 1) == 0.).or.(fraction_mainarea(i, 1) == maxlc))then
                    cycle
                end if

                if (refine_num_landtypes .eqv. .true. .and. (n_landtypes(i) > th_num_landtypes)) then
                    ref(i) = 1
                    ref_th(i, 1) = 1
                end if

                if (refine_area_mainland .eqv. .true. .and. (fraction_mainarea(i, 2) < th_area_mainland)) then
                    ref(i) = 1
                    ref_th(i, 2) = 1
                end if

                if (refine_lai_m .eqv. .true. .and. (p_lai(i, 1) > th_lai_m)) then
                    ref(i) = 1
                    ref_th(i, 3) = 1
                end if

                if (refine_lai_s .eqv. .true. .and. (p_lai(i, 2) > th_lai_s)) then
                    ref(i) = 1
                    ref_th(i, 4) = 1
                end if

                if (refine_slope_m .eqv. .true. .and. (p_slope(i, 1) > th_slope_m)) then
                    ref(i) = 1
                    ref_th(i, 5) = 1
                end if

                if (refine_slope_s .eqv. .true. .and. (p_slope(i, 2) > th_slope_s)) then
                    ref(i) = 1
                    ref_th(i, 6) = 1
                end if

                if (refine_k_s_m .eqv. .true. .and. ((p_k_s(i, 1, 1) > th_k_s_m).or.(p_k_s(i, 1, 2) > th_k_s_m))) then
                    ref(i) = 1
                    ref_th(i, 7) = 1
                end if

                if (refine_k_s_s .eqv. .true. .and. ((p_k_s(i, 2, 1) > th_k_s_s).or.(p_k_s(i, 2, 2) > th_k_s_s))) then
                    ref(i) = 1
                    ref_th(i, 8) = 1
                end if

                if (refine_k_solids_m .eqv. .true. .and. ((p_k_sl(i, 1, 1) > th_k_solids_m).or.(p_k_sl(i, 1, 2) > th_k_solids_m))) then
                    ref(i) = 1
                    ref_th(i, 9) = 1
                end if

                if (refine_k_solids_s .eqv. .true. .and. ((p_k_sl(i, 2, 1) > th_k_solids_s).or.(p_k_sl(i, 2, 2) > th_k_solids_s))) then
                    ref(i) = 1
                    ref_th(i, 10) = 1
                end if

                if (refine_tkdry_m .eqv. .true. .and. ((p_tkdry(i, 1, 1) > th_tkdry_m).or.(p_tkdry(i, 1, 2) > th_tkdry_m))) then
                    ref(i) = 1
                    ref_th(i, 11) = 1
                end if

                if (refine_tkdry_s .eqv. .true. .and. ((p_tkdry(i, 2, 1) > th_tkdry_s).or.(p_tkdry(i, 2, 2) > th_tkdry_s))) then
                    ref(i) = 1
                    ref_th(i, 12) = 1
                end if

                if (refine_tksatf_m .eqv. .true. .and. ((p_tksatf(i, 1, 1) > th_tksatf_m).or.(p_tksatf(i, 1, 2) > th_tksatf_m))) then
                    ref(i) = 1
                    ref_th(i, 13) = 1
                end if

                if (refine_tksatf_s .eqv. .true. .and. ((p_tksatf(i, 2, 1) > th_tksatf_s).or.(p_tksatf(i, 2, 2) > th_tksatf_s))) then
                    ref(i) = 1
                    ref_th(i, 14) = 1
                end if

                if (refine_tksatu_m .eqv. .true. .and. ((p_tksatu(i, 1, 1) > th_tksatu_m).or.(p_tksatu(i, 1, 2) > th_tksatu_m))) then
                    ref(i) = 1
                    ref_th(i, 15) = 1
                end if
                if (refine_tksatu_s .eqv. .true. .and. ((p_tksatu(i, 2, 1) > th_tksatu_s).or.(p_tksatu(i, 2, 2) > th_tksatu_s))) then
                    ref(i) = 1
                    ref_th(i, 16) = 1
                end if

            end do

            !----------------------------------
            ! 6.根据相邻距离筛选七边形附近点
            !----------------------------------
            print*, "根据相邻距离筛选七边形附近点"

            do i = 1, sjx_points, 1
                if(ref(i) == 1)then
                    do j = 1, sjx_points, 1
                        if(mp_dis(i, j) <= set_dis)then
                            do k = 1, 3, 1
                                if(ngrwm(8, ngrmw(k, j)) == 7)then
                                    ref(i) = 0
                                end if
                            end do
                        end if
                    end do
                end if
            end do

            !   do i = 1,sjx_points,1
            !      if(ref(i) == 0)then
            !         cycle
            !      end if
            !      do j = 1,3,1
            !         id(1,j) = ngrmm(j,i)
            !         do k = 1,3,1
            !            if(ngrwm(8,ngrmw(k,id(1,j))) == 7)then
            !               ref(i) = 0
            !               exit
            !            end if
            !         end do
            !
            !         if(ref(i) == 0)then
            !            exit
            !         end if
            !
            !         do k = 1,3,1
            !            id(2,k) = ngrmm(k,id(1,j))
            !            do l = 1,3,1
            !               if(ngrwm(8,ngrmw(k,id(2,k))) == 7)then
            !                  ref(i) = 0
            !                  exit
            !               end if
            !            end do
            !         end do
            !      end do
            !   end do

            iter = 1                         ! 迭代次数
            num_ref = INT(sum(ref))          ! 每次细化三角形数

            if(num_ref == 0)then
                iter_is_end = .true.
            end if

            nmp(1) = sjx_points + 4 * num_ref  ! 记录每次迭代后三角形数
            nwp(1) = lbx_points + 3 * num_ref  ! 记录每次迭代后多边形数

            mrl_new(1:sjx_points) = mrl
            ngr_mrl_new(1:3, 1:sjx_points) = ngr_mrl
            mp_new(1:sjx_points, 1:2) = mp(1:sjx_points, 1:2)
            wp_new(1:lbx_points, 1:2) = wp(1:lbx_points, 1:2)
            ngrmw_new(1:3, 1:sjx_points) = ngrmw(1:3, 1:sjx_points)
            ngrmm_new(1:3, 1:sjx_points) = ngrmm(1:3, 1:sjx_points)
            ngrwm_new(1:7, 1:lbx_points) = ngrwm(1:7, 1:lbx_points)

            deallocate(n_landtypes)
            deallocate(fraction_mainarea)
            deallocate(p_slope)
            deallocate(p_lai)
            deallocate(p_k_s)
            deallocate(p_k_sl)
            deallocate(p_tkdry)
            deallocate(p_tksatf)
            deallocate(p_tksatu)


            !-----------------------------------------------
            ! 7.初步细化（一分为四）
            !-----------------------------------------------
            print*, "开始初步细化"

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

                    mrl_new(i) = 4 ! 表示将三角形网格平均分为四份
                    mrl_new(sjx_points + refed(iter) * 4 + 1:sjx_points + refed(1) * 4 + 4) = 4

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
            print*, "细化第一步完成"

            !--------------------------------------------------
            ! 8.储存初始网格数据和初步细化后新网格数据
            !--------------------------------------------------
            lndname = trim(base_dir) // trim(EXPNME)  // "/makegrid/tmp/gridfile_NXP" // trim(nxpc)  // "_" // trim(steps) // "_ori.nc4"
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
            CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5), ngrmw(1:3, :)))
            CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6), ngrwm(1:7, :)))
            CALL CHECK(NF90_CLOSE(ncID))

            lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/gridfile_NXP" // trim(nxpc)  // "_" // trim(steps) // "_1.nc4"
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

            !-----------------------------------------------
            ! 9.进行迭代（防止细化交汇带出现冲突，一分为四）
            !-----------------------------------------------
            iterA = .false.
            iterB = .false.
            iterC = .false.

            print*, "iterA start"
            do while(iterA .eqv. .false.)

                iterA = .true.    ! 判断迭代B和迭代C是否都已满足条件
                iterB = .false.   ! 从三角形网格进行判断
                iterC = .false.   ! 从多边形网格进行判断

                print*, "iterB start"
                do while(iterB .eqv. .false.)

                    ref = 0

                    ! 当一个未细化三角形有两个或两个以上细化相邻三角形时，细化该三角形
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
                do while(iterC .eqv. .false.)

                    ref = 0

                    do i = 1, lbx_points, 1

                        edges = 0      ! 多边形网格边数

                        do j = 1, 7, 1   ! 计算未细化三角形构成的多边形的边数
                            if(ngrwm(j, i) /= 1)then
                                m1 = ngrwm(j, i)
                                if(sum(ngr_mrl_new(:, m1)) == 6.)then
                                    ref_lbx(i, j) = 1  ! 表示构成该多边形的部分三角形已被细化
                                end if
                                edges = edges + 1
                            end if
                        end do

                        if(ref_lbx(i, 8) == 0)then     ! 表示构成该多边形的部分三角形均未被细化
                            if(edges == 5)then         ! 五边形
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
                            else if(edges == 6)then    ! 六边形
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
            print*, "细化第二步完成"


            !-----------------------------------------------
            ! 10.储存第二步细化后新网格数据
            !-----------------------------------------------
            lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/gridfile_NXP" // trim(nxpc)  // "_" // trim(steps) // "_2.nc4"
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


            !-----------------------------------------------
            ! 11.寻找弱凹点
            !-----------------------------------------------
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


            !-----------------------------------------------
            ! 12.细化弱凹点
            !-----------------------------------------------
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

            !-----------------------------------------------
            ! 13.储存第三步细化后新网格数据
            !-----------------------------------------------
            lndname = trim(base_dir) // trim(EXPNME)  // "/makegrid/tmp/gridfile_NXP" // trim(nxpc)  // "_" // trim(steps) // "_3.nc4"
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


            !-----------------------------------------------
            ! 14.记录相邻只有一个细化的三角形
            !-----------------------------------------------
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


            !-----------------------------------------------
            ! 15.细化相邻只有一个细化的三角形（一分为二）
            !-----------------------------------------------
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

                    mp_new(m1, 1:2) = (sjx(1, :) + newsjx(1, :) + sjx(2, :)) / 3.
                    mp_new(m2, 1:2) = (sjx(1, :) + newsjx(1, :) + sjx(3, :)) / 3.

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
            ! 16.储存第四步细化后新网格数据
            !--------------------------------------------------
            lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/gridfile_NXP" // trim(nxpc)  // "_" // trim(steps) // "_4.nc4"
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


            !--------------------------------------------------
            ! 17.计算并储存非重复w点
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
            ! 18.重新计算并储存ngrmw
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
            ! 19.重新计算并储存ngrwm
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
            ! 20.去除已细化三角形（m点）
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
            ! 21.存储网格数据
            !--------------------------------------------------
            lndname = trim(base_dir) // trim(EXPNME)  // "/makegrid/tmp/gridfile_NXP" // trim(nxpc)  // "_" // trim(steps) // "_5.nc4"
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
            ! 22.计算三角形网格边长与角度
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
            ! 23.计算弹性调整前三角形网格质量
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
            ! 24.调整所有m点至三角形网格重心
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
            ! 25.计算弹性调整前多边形网格质量
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
            ! 26.弹性调整&计算每次调整后的三角形网格质量
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
            ! 27.存储网格质量数据
            !--------------------------------------------------
            if(iter_is_end .eqv. .true.)then
                !lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/result/quality.nc4"
                lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/result/quality_NXP" // trim(nxpc)  // "_lbx.nc4"
            else
                lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/quality_NXP" // trim(nxpc)  // "_" // trim(steps) // ".nc4"
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
            ! 28.存储最终网格数据
            !--------------------------------------------------
            
            allocate(mp_f_tmp(1:num_sjx,1:2)); mp_f_tmp = mp_f(1:num_sjx,1:2)
            allocate(wp_f_tmp(1:num_dbx,1:2)); wp_f_tmp = wp_f(1:num_dbx,1:2)
            allocate(ngrmw_f_tmp(1:3, 1:num_sjx)); ngrmw_f_tmp = ngrmw_f(1:3, 1:num_sjx)
            allocate(ngrwm_f_tmp(1:7, 1:num_dbx)); ngrwm_f_tmp = ngrwm_f(1:7, 1:num_dbx)
            allocate(dismm(num_dbx,7)); dismm = 0.
            allocate(disww(num_sjx,3)); disww = 0.
            Call Get_Dis(mp_f_tmp,wp_f_tmp,ngrmw_f_tmp,ngrwm_f_tmp,dismm,disww,num_sjx,num_dbx)

            if(iter_is_end .eqv. .true.)then
               !lndname = trim(base_dir) // trim(EXPNME) //"/makegrid/result/gridfile.nc4"
               lndname = trim(base_dir) // trim(EXPNME) //"/makegrid/result/gridfile_NXP" // trim(nxpc)  // "_lbx.nc4"
            else
               lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/gridfile_NXP" // trim(nxpc)  // "_" // trim(steps) // "_6.nc4"
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

            step = step + 1

            deallocate(wp)
            deallocate(mp)
            deallocate(mp_id)
            deallocate(mp_id_new)
            deallocate(wp_new)
            deallocate(mp_new)
            deallocate(mp_same)
            deallocate(ngrwm)
            deallocate(ngrwm_new)
            deallocate(ngrmw)
            deallocate(ngrmm)
            deallocate(mrl)
            deallocate(mrl_new)
            deallocate(ngr_mrl)
            deallocate(ngr_mrl_new)
            deallocate(ngrmm_new)
            deallocate(ngrmw_new)
            deallocate(wp_ori)
            deallocate(mp_ori)
            deallocate(mp_id_ori)
            deallocate(ngrwm_ori)
            deallocate(ngrmw_ori)
            deallocate(mp_dis)
            deallocate(mp_ii_ori)
            deallocate(lon_i)
            deallocate(lat_i)
            deallocate(area_fine_gridcell)
            deallocate(landtypes)
            deallocate(seaorland)
            deallocate(slope_max)
            deallocate(slope_avg)
            deallocate(ref_jw)
            deallocate(mp_ii)
            deallocate(mp_ii_new)
            deallocate(lai)
            deallocate(k_s)
            deallocate(k_sl)
            deallocate(tkdry)
            deallocate(tksatf)
            deallocate(tksatu)
            deallocate(nl)
            deallocate(area)
            deallocate(ref_tr)
            deallocate(ref)
            deallocate(ref_th)
            deallocate(ref_lbx)
            deallocate(ref_l)
            deallocate(wp_f)
            deallocate(ngrwm_f)
            deallocate(ref_pl)
            deallocate(mp_ref)
            deallocate(mp_f)
            deallocate(mrl_f)
            deallocate(ngrmw_f)
            deallocate(length)
            deallocate(angle)
            deallocate(Eavg_sjx)
            deallocate(Extr_sjx)
            deallocate(Savg_sjx)
            deallocate(less30)
            deallocate(angle_wbx)
            deallocate(angle_lbx)
            deallocate(angle_qbx)
            deallocate(Eavg_wbx)
            deallocate(Eavg_lbx)
            deallocate(Eavg_qbx)
            deallocate(Savg_wbx)
            deallocate(Savg_lbx)
            deallocate(Savg_qbx)
            deallocate(Extr_wbx)
            deallocate(Extr_lbx)
            deallocate(Extr_qbx)
            deallocate(MoveDis)
            deallocate(IsInRfArea)
            deallocate(dismm)
            deallocate(disww)
            deallocate(mp_f_tmp)
            deallocate(wp_f_tmp)
            deallocate(ngrmw_f_tmp)
            deallocate(ngrwm_f_tmp)

        end do

    end subroutine refine_lbx_step2

    subroutine CHECK(STATUS)
        INTEGER, intent (in) :: STATUS
        if  (STATUS /= NF90_NOERR) then ! nf_noerr=0 表示没有错误
            print *, NF90_STRERROR(STATUS)
            stop 'stopped'
        endif
    end subroutine CHECK


    INTEGER FUNCTION IsCrossLine(x1, x2)

        implicit none

        real(r8), intent(in) :: x1, x2
        IsCrossLine = 0

        if(abs(x1 - x2) > 180.)then
            IsCrossLine = 1
        end if

    END FUNCTION IsCrossLine


    SUBROUTINE CheckLon(x)

        implicit none

        real(r8), intent(out) :: x

        if(x > 180.)then
            x = x - 360.
        else if(x < -180.)then
            x = x + 360.
        end if

    END SUBROUTINE CheckLon

    SUBROUTINE GetRefineJW(a, b, c, ref_jw)

        implicit none

        real(r8), intent(out) :: a(2), b(2), c(2)
        real(r8) :: lon(nlons_source), lat(nlats_source), dx, dy
        real(r8) :: minlon, minlat, maxlon, maxlat
        integer, dimension(nlons_source, nlats_source), intent(out) :: ref_jw
        integer :: i, j, icl(3)

        dx = 360. / nlons_source
        dy = 180. / nlats_source

        do i = 1, nlons_source, 1
            lon(i) = -180. + (2 * i - 1) * dx / 2.
        end do

        do j = 1, nlats_source, 1
            lat(j) = 90. - (2 * j - 1) * dy / 2.
        end do

        icl = 0
        icl(1) = IsCrossLine(a(1), b(1))
        icl(2) = IsCrossLine(a(1), c(1))
        icl(3) = IsCrossLine(b(1), c(1))

        if(INT(sum(icl)) > 0)then
            if(a(1) < 0.)then
                a(1) = a(1) + 360.
            end if
            if(b(1) < 0.)then
                b(1) = b(1) + 360.
            end if
            if(c(1) < 0.)then
                c(1) = c(1) + 360.
            end if
            do i = 1, nlons_source, 1
                if(lon(i) < 0.)then
                    lon(i) = lon(i) + 360.
                end if
            end do
        end if

        minlon = min(a(1), b(1), c(1))
        minlat = min(a(2), b(2), c(2))
        maxlon = max(a(1), b(1), c(1))
        maxlat = max(a(2), b(2), c(2))

        do i = 1, nlons_source, 1
            if((lon(i)>minlon).and.(lon(i)<maxlon))then
                do j = 1, nlats_source, 1
                    if((lat(j)>minlat).and.(lat(j)<maxlat))then
                        ref_jw(i, j) = 1
                    end if
                end do
            end if
        end do

        if(INT(sum(icl)) > 0)then
            CALL CheckLon(a(1))
            CALL CheckLon(b(1))
            CALL CheckLon(c(1))
        end if

    END SUBROUTINE GetRefineJW


    SUBROUTINE CellArea(area)

        implicit none

        real(r8), intent(out) :: area(nlons_source, nlats_source)
        integer :: i, j
        real(r8) :: re, pi, deg2rad, global, dx, dy, error
        real(r8) :: lats(nlats_source), latn(nlats_source)
        real(r8) :: lonw(nlons_source), lone(nlons_source)

        re = 6.37122e6 * 0.001                    ! kilometer
        pi = 4. * atan(1.)
        deg2rad = pi / 180.
        global = 0.

        dx = 360. / nlons_source
        dy = 180. / nlats_source

        do i = 1, nlons_source, 1
            lone(i) = -180. + i * dx
            lonw(i) = -180. + (i - 1) * dx
        end do

        do i = 1, nlats_source, 1
            latn(i) = 90. - (i - 1) * dy
            lats(i) = 90. - i * dy
        end do

        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1)&
        !$OMP PRIVATE(i,j,dx,dy)
        do j = 1, nlats_source, 1
            do i = 1, nlons_source, 1
                if(lone(i)<lonw(i))then   ! west edge is more western than data line
                    ! 西部边缘处于日期线西方
                    dx = (lone(i) - lonw(i) + 360.0) * deg2rad
                else
                    dx = (lone(i) - lonw(i)) * deg2rad
                endif
                if(latn(j)>lats(j)) then          ! north to south grid
                    dy = sin(latn(j) * deg2rad) - sin(lats(j) * deg2rad)
                else                              ! south to north grid
                    dy = sin(lats(j) * deg2rad) - sin(latn(j) * deg2rad)
                end if
                area(i, j) = dx * dy * re * re
                ! 弧长公式解求面积
            end do
        end do
        !$OMP END PARALLEL DO

        global = sum(area(:, :))

        ! 确保网格单元的总面积与其边缘定义的网格面积相同
        dx = (180. - (-180.)) * deg2rad
        dy = sin(90. * deg2rad) - sin(-90. * deg2rad)
        error = dx * dy * re * re
        if(abs(global - error) / error > 1.0e-7) then
            print*, 'CELLAREA error: correct area is ', error, &
                    ' but summed area of grid cells is ', global
        end if

        return

    END SUBROUTINE CellArea


    REAL FUNCTION IsInUstrGrid(ustr, lon, lat)

        implicit none

        integer :: inc(4), i, j, iscross_l(2), stat, num_inter
        integer, allocatable :: iscross_g(:)                   ! 判断非结构网格线段是否穿过经纬度网格
        real(r8), intent(in) :: lon, lat
        real(r8) :: minlat, maxlat, dx, dy, minlon, maxlon
        real(r8), intent(in) :: ustr(3, 2)     ! 非结构网格单元顶点
        real(r8) :: ustr_move(3, 2)                ! 非结构网格经度移动后的点
        real(r8) :: point(4, 2)                      ! 经纬度网格顶点
        real(r8) :: center_point(2)                 ! 非结构网格单元中心点
        real(r8) :: interarea_points(20, 2)          ! 两网格相交区域顶点
        real(r8), allocatable :: inter_points(:, :, :)           ! 两种网格的交点
        real(r8), allocatable :: area(:)   ! 由非结构网格相邻两点与任一经纬度网格顶点组成的三角形面积
        real(r8), dimension(5) :: area_i   ! 由经纬度网格相邻两点与任一非结构网格顶点组成的三角形面积
        real(r8) :: tmpa(2),tmpb(2),tmpc(2),tmpd(3,2)

        !MoveLons              ! 移动越过180°经线的经度
        !IsCrossGrid           ! 判断非结构网格线段是否与经纬度网格相交
        !IsCrossLine2          ! 判断网格是否越过180°经线
        !SortPoints            ! 对多边形各个顶点进行排序
        !GetTriangleArea          ! 获取三角形网格面积
        !GetAreaPercent           ! 获取非结构网格占经纬度网格比例

        maxlat = 0.
        minlat = 0.
        maxlon = 0.
        minlon = 0.

        iscross_l = 0
        iscross_g = 0
        inc = 0                 ! 判断经纬度网格顶点在非结构网格中的数量
        IsInUstrGrid = 0        ! 判断经纬度网格与非结构网格位置关系
        center_point = 0        ! 经纬度网格中心点
        num_inter = 0           ! 两种网格重合面积多边形顶点数量
        interarea_points = 0.   ! 两种网格重叠多边形顶点
        area_i = 0.
        point = 0.

        dx = 360. / nlons_source
        dy = 180. / nlats_source

        ! 计算经纬度网格顶点坐标
        point(1, 1) = lon + dx / 2.
        point(1, 2) = lat + dy / 2.
        point(2, 1) = lon - dx / 2.
        point(2, 2) = lat + dy / 2.
        point(3, 1) = lon - dx / 2.
        point(3, 2) = lat - dy / 2.
        point(4, 1) = lon + dx / 2.
        point(4, 2) = lat - dy / 2.

        ! 确保经纬度网格纬度绝对值不超过90°
        do i = 1, 4, 1
            if(point(i, 2) > 90.)then
                point(i, 2) = 90.
            else if(point(i, 2) < -90.)then
                point(i, 2) = -90.
            end if
        end do

        !------------------------------------------------------------------
        ! 根据纬度初步筛选
        !------------------------------------------------------------------

        maxlat = maxval(ustr(1:3, 2))
        minlat = minval(ustr(1:3, 2))

        if((maxlat > 85.).or.(minlat < -85.))then
            IsInUstrGrid = -1
            return
        end if

        if((point(4, 2) > maxlat).or.(point(1, 2) < minlat))then
            IsInUstrGrid = -1
            return
        end if

        ustr_move = ustr

        allocate(inter_points(3, 3, 2))
        allocate(area(4))
        allocate(iscross_g(3))
        inter_points = 0
        iscross_g = 0
        area = 0

        !-------------------------------------------------------------------------------------
        ! 判断两网格是否越过±180°经线
        !-------------------------------------------------------------------------------------
        if((point(1, 1) > 180.).or.(point(2, 1) < -180.))then
            iscross_l(1) = 1
        end if

        iscross_l(2) = IsCrossLine2(ustr_move(:, 1), 3)

        if(iscross_l(2) == 1)then
            iscross_l(1) = 1
        end if

        !-------------------------------------------------------------------------------------
        ! 根据上述判断移动网格点经度
        !-------------------------------------------------------------------------------------
        if(iscross_l(1) == 1)then
            call MoveLons(point(:, 1), 4)
        end if
        if(iscross_l(2) == 1)then
            call MoveLons(ustr_move(:, 1), 3)
        end if

        !-------------------------------------------------------------------------------------
        ! 根据经度筛选网格
        !-------------------------------------------------------------------------------------
        minlon = minval(ustr_move(1:3, 1))
        maxlon = maxval(ustr_move(1:3, 1))

        if((point(2, 1) > maxlon).or.(point(1, 1) < minlon))then
            IsInUstrGrid = -1
            return
        end if

        !-------------------------------------------------------------------------------------
        ! 判断非结构网格是否被经纬度网格包含
        !-------------------------------------------------------------------------------------
        if((point(1, 1) > maxlon).and.(point(2, 1) < minlon).and.&
                (point(1, 2) > maxlat).and.(point(3, 2) < minlat))then
            IsInUstrGrid = 0.9
            return
        end if

        !-------------------------------------------------------------------------------------
        ! 开始判断两个网格间的位置关系，并记录关键点
        !-------------------------------------------------------------------------------------
        tmpa = ustr_move(1, 1:2)
        tmpb = ustr_move(2, 1:2)
        tmpd = inter_points(1, 1:3, 1:2)
        iscross_g(1) = IsCrossGrid(point, tmpa, tmpb, tmpd)
        !iscross_g(1) = IsCrossGrid(point, ustr_move(1, :), ustr_move(2, :), inter_points(1, :, :))
        tmpa = ustr_move(2, 1:2)
        tmpb = ustr_move(3, 1:2)
        tmpd = inter_points(2, 1:3, 1:2)
        iscross_g(2) = IsCrossGrid(point, tmpa, tmpb, tmpd)
        !iscross_g(2) = IsCrossGrid(point, ustr_move(2, :), ustr_move(3, :), inter_points(2, :, :))
        tmpa = ustr_move(1, 1:2)
        tmpb = ustr_move(3, 1:2)
        tmpd = inter_points(3, 1:3, 1:2)
        iscross_g(3) = IsCrossGrid(point, tmpa, tmpb, tmpd)
        !iscross_g(3) = IsCrossGrid(point, ustr_move(1, :), ustr_move(3, :), inter_points(3, :, :))

        ! 计算经纬度网格位于非结构网格中的顶点数目
        ! 若为4，则包含，否则相交
        do i = 1, 4, 1
            area = 0.
            tmpa = ustr_move(1, 1:2)
            tmpb = ustr_move(2, 1:2)
            tmpc = point(i, 1:2)
            area(1) = GetTriangleArea(tmpa, tmpb, tmpc)
            !area(1) = GetTriangleArea(ustr_move(1, :), ustr_move(2, :), point(i, :))
            tmpa = ustr_move(2, 1:2)
            tmpb = ustr_move(3, 1:2)
            tmpc = point(i, 1:2)
            area(2) = GetTriangleArea(tmpa, tmpb, tmpc)
            !area(2) = GetTriangleArea(ustr_move(2, :), ustr_move(3, :), point(i, :))
            tmpa = ustr_move(1, 1:2)
            tmpb = ustr_move(3, 1:2)
            tmpc = point(i, 1:2)
            area(3) = GetTriangleArea(tmpa, tmpb, tmpc)
            !area(3) = GetTriangleArea(ustr_move(1, :), ustr_move(3, :), point(i, :))
            tmpa = ustr_move(1, 1:2)
            tmpb = ustr_move(2, 1:2)
            tmpc = ustr_move(3, 1:2)
            area(4) = GetTriangleArea(tmpa, tmpb, tmpc)
            !area(4) = GetTriangleArea(ustr_move(1, :), ustr_move(2, :), ustr_move(3, :))

            if(abs(area(1) + area(2) + area(3) - area(4)) < 0.000005)then
                inc(i) = 1
            end if
        end do

        !print*,"inc",inc

        !----------------------------------------------------------------------------------------------------
        ! 计算非结构网格与经纬度网格重叠面积
        ! 1.寻找非结构网格位于经纬度网格中的点
        ! 2.寻找经纬度网格位于非结构网格中的点
        ! 3.寻找非结构网格与经纬度网格的交点
        ! 4.将它们进行排序，组合成凸多边形（或三角形）
        ! 5.计算该图形面积
        !---------------------------------------------------------------------------------------------------
        if (no_caculate_fraction) then
            if((sum(inc) > 0).or.(sum(iscross_g) > 0))then
                IsInUstrGrid = 1.
                return
            else
                IsInUstrGrid = 0
                return
            end if
        end if

        if(sum(inc) == 4)then   ! 若包含
            IsInUstrGrid = 1
            return
        else if((sum(inc) == 0).and.(sum(iscross_g) == 0))then  ! 若不包含也不相交
            IsInUstrGrid = -1
            return
        else
            do i = 1, 4, 1
                if(inc(i) == 1)then
                    num_inter = num_inter + 1
                    interarea_points(num_inter, :) = point(i, :)
                end if
            end do

            center_point = 0.
            do i = 1, 4, 1
                center_point = center_point + point(i, :)
            end do
            center_point = center_point / 4

            do i = 1, 3, 1
                area_i = 0.
                do j = 1, 3, 1
                    tmpa = ustr_move(i, 1:2)
                    tmpb = point(j, 1:2)
                    tmpc = point(j + 1, 1:2)
                    area_i(j) = GetTriangleArea(tmpa, tmpb, tmpc)
                    !area_i(j) = GetTriangleArea(ustr_move(i, :), point(j, :), point(j + 1, :))
                    tmpb = point(j, 1:2)
                    tmpc = point(j + 1, 1:2)
                    area_i(5) = area_i(5) + GetTriangleArea(center_point, tmpb, tmpc)
                    !area_i(5) = area_i(5) + GetTriangleArea(center_point, point(j, :), point(j + 1, :))
                end do

                tmpa = ustr_move(i, 1:2)
                tmpb = point(1, 1:2)
                tmpc = point(4, 1:2)
                area_i(4) = GetTriangleArea(tmpa, tmpb, tmpc)
                !area_i(4) = GetTriangleArea(ustr_move(i, :), point(1, :), point(4, :))
                tmpb = point(1, 1:2)
                tmpc = point(4, 1:2)
                area_i(5) = area_i(5) + GetTriangleArea(center_point, tmpb, tmpc)
                !area_i(5) = area_i(5) + GetTriangleArea(center_point, point(1, :), point(4, :))
                if(abs(sum(area_i(1:4)) - area_i(5)) < 0.00006)then
                    num_inter = num_inter + 1
                    interarea_points(num_inter, :) = ustr_move(i, :)
                end if
            end do

            do i = 1, 3, 1
                if(inter_points(i, 3, 1) /= 0)then
                    num_inter = num_inter + 1
                    interarea_points(num_inter:num_inter + inter_points(i, 3, 1) - 1, :) = &
                            inter_points(i, 1:inter_points(i, 3, 1), :)
                    num_inter = num_inter + inter_points(i, 3, 1) - 1
                end if
            end do

            if(num_inter == 0)then 
               IsInUstrGrid = 0
               return
            end if

            CALL SortPoints(interarea_points, num_inter)

            IsInUstrGrid = GetAreaPercent(interarea_points, num_inter, point)

        end if

    END FUNCTION IsInUstrGrid


    ! 判断网格是否越过189°与-180°经线
    INTEGER FUNCTION IsCrossLine2(lons, num)

        implicit none

        integer :: num, i, j
        real(r8), dimension(num) :: lons

        IsCrossLine2 = 0

        do i = 1, num - 1, 1
            do j = i + 1, num, 1
                if(abs(lons(j) - lons(i)) > 180.)then
                    IsCrossLine2 = 1
                    return
                end if
            end do
        end do

    END FUNCTION IsCrossLine2


    ! 调整网格点经度
    SUBROUTINE MoveLons(lons, num)         ! lor = left or right

        implicit none

        integer :: i
        integer, intent(in) :: num
        real(r8), dimension(num) :: lons

        do i = 1, num, 1
            if(lons(i) < 0.)then
                lons(i) = lons(i) + 360.
            end if
        end do

    END SUBROUTINE MoveLons


    ! 判断非结构网格线段是否穿过经纬度网格
    INTEGER FUNCTION IsCrossGrid(point, a, b, inter_point)

        implicit none

        real(r8), dimension(2) :: a, b      ! 直线端点
        real(r8), dimension(2) :: x, y
        real(r8), dimension(4, 2) :: point
        real(r8), dimension(3, 2) :: inter_point
        real(r8) :: x1, x2, y1, y2, m, n, num

        IsCrossGrid = 0
        inter_point = 0
        num = 0

        x(1) = max(a(1), b(1))
        x(2) = min(a(1), b(1))
        y(1) = max(a(2), b(2))
        y(2) = min(a(2), b(2))

        if(a(1) == b(1))then
            if(a(1)>point(2, 1).and.a(1)<point(1, 1))then
                if((y(1)>point(1, 2)).and.((y(2)>point(4, 2))).and.(y(2)<point(1, 2)))then
                    inter_point(1, 1) = a(1)
                    inter_point(1, 2) = point(1, 2)
                    inter_point(3, 1) = 1
                    IsCrossGrid = 1
                    return
                else if((y(1)>point(1, 2)).and.((y(2)<point(4, 2))))then
                    inter_point(1, 1) = a(1)
                    inter_point(1, 2) = point(1, 2)
                    inter_point(2, 1) = a(1)
                    inter_point(2, 2) = point(4, 2)
                    inter_point(3, 1) = 2
                    IsCrossGrid = 2
                    return
                else if((y(1)>point(4, 2)).and.(y(1)<point(1, 2)).and.(y(2)<point(4, 2)))then
                    inter_point(1, 1) = a(1)
                    inter_point(1, 2) = point(4, 2)
                    inter_point(3, 1) = 1
                    IsCrossGrid = 1
                    return
                end if
            end if
            IsCrossGrid = 0
            return
        else if(a(2) == b(2))then
            if(a(2)>point(4, 2).and.a(2)<point(1, 2))then
                if((x(1)>point(1, 1)).and.((x(2)>point(2, 1))).and.(x(2)<point(1, 1)))then
                    inter_point(1, 1) = point(1, 1)
                    inter_point(1, 2) = a(2)
                    inter_point(3, 1) = 1
                    IsCrossGrid = 1
                    return
                else if((x(1)>point(1, 1)).and.((x(2)<point(2, 1))))then
                    inter_point(1, 1) = point(1, 1)
                    inter_point(1, 2) = a(2)
                    inter_point(2, 1) = point(2, 1)
                    inter_point(2, 2) = a(2)
                    inter_point(3, 1) = 2
                    IsCrossGrid = 2
                    return
                else if((x(1)>point(2, 1)).and.(x(1)<point(1, 1)).and.(x(2)<point(2, 1)))then
                    inter_point(1, 1) = point(2, 1)
                    inter_point(1, 2) = a(2)
                    inter_point(3, 1) = 1
                    IsCrossGrid = 1
                    return
                end if
            end if
            IsCrossGrid = 0
            return
        else
            m = (a(2) - b(2)) / (a(1) - b(1))
            n = a(2) - m * a(1)
            y1 = m * point(1, 1) + n
            y2 = m * point(2, 1) + n
            x1 = (point(1, 2) - n) / m
            x2 = (point(4, 2) - n) / m
            if((y1>y(2)).and.(y1<y(1)).and.(y1>point(4, 2)).and.(y1<point(1, 2)))then
                IsCrossGrid = IsCrossGrid + 1
                inter_point(IsCrossGrid, 1) = point(1, 1)
                inter_point(IsCrossGrid, 2) = y1
            end if
            if((y2>y(2)).and.(y2<y(1)).and.(y2>point(4, 2)).and.(y2<point(1, 2)))then
                IsCrossGrid = IsCrossGrid + 1
                inter_point(IsCrossGrid, 1) = point(2, 1)
                inter_point(IsCrossGrid, 2) = y2
            end if
            if((x1>x(2)).and.(x1<x(1)).and.(x1>point(2, 1)).and.(x1<point(1, 1)))then
                IsCrossGrid = IsCrossGrid + 1
                inter_point(IsCrossGrid, 1) = x1
                inter_point(IsCrossGrid, 2) = point(1, 2)
            end if
            if((x2>x(2)).and.(x2<x(1)).and.(x2>point(2, 1)).and.(x2<point(1, 1)))then
                IsCrossGrid = IsCrossGrid + 1
                inter_point(IsCrossGrid, 1) = x2
                inter_point(IsCrossGrid, 2) = point(4, 2)
            end if

            if(IsCrossGrid > 2)then
                IsCrossGrid = 0
            end if

            inter_point(3, 1) = IsCrossGrid
        end if

    END FUNCTION IsCrossGrid


    ! 计算三角形面积（按经纬度，并非实际面积）
    REAL FUNCTION GetTriangleArea(a, b, c)

        implicit none

        real(r8), dimension(2), intent(in) :: a, b, c
        real(r8) :: aa, bb, cc, p

        GetTriangleArea = 0

        aa = sqrt((c(1) - b(1)) * (c(1) - b(1)) + (c(2) - b(2)) * (c(2) - b(2)))
        bb = sqrt((c(1) - a(1)) * (c(1) - a(1)) + (c(2) - a(2)) * (c(2) - a(2)))
        cc = sqrt((a(1) - b(1)) * (a(1) - b(1)) + (a(2) - b(2)) * (a(2) - b(2)))

        p = (aa + bb + cc) / 2

        GetTriangleArea = sqrt(p * (p - aa) * (p - bb) * (p - cc))

    END FUNCTION GetTriangleArea


    ! 获取非结构网格包含经纬度网格的比例
    REAL FUNCTION GetAreaPercent(inter_point, num, point)

        implicit none

        integer, intent(in) :: num
        real(r8), dimension(20, 2), intent(in) :: inter_point
        real(r8), dimension(4, 2), intent(in) :: point
        real(r8), dimension(2) :: center_point
        integer :: i
        real(r8) :: inter_area, tmpa(2), tmpb(2)

        GetAreaPercent = 0
        inter_area = 0
        center_point = 0

        do i = 1, num, 1
            center_point = center_point + inter_point(i, :)
        end do
        center_point = center_point / num

        do i = 1, num - 1, 1
            tmpa = inter_point(i, 1:2)
            tmpb = inter_point(i + 1, 1:2)
            inter_area = inter_area + GetTriangleArea(center_point, tmpa, tmpb)
            !inter_area = inter_area + GetTriangleArea(center_point, inter_point(i, :), inter_point(i + 1, :))
        end do

        tmpa = inter_point(1, 1:2)
        tmpb = inter_point(num, 1:2)
        inter_area = inter_area + GetTriangleArea(center_point, tmpa, tmpb)
        !inter_area = inter_area + GetTriangleArea(center_point, inter_point(1, :), inter_point(num, :))

        GetAreaPercent = inter_area / (abs((point(1, 1) - point(2, 1)) * (point(1, 2) - point(4, 2))))
        ! 面积占比为重合三角形除以经纬度网格面积

        !if(GetAreaPercent >= 1)then
        !        GetAreaPercent = 0
        !end if

    END FUNCTION GetAreaPercent


    ! 将点排序成多边形
    SUBROUTINE SortPoints(points, num)

        implicit none

        integer :: i, j, x
        real(r8) :: angle_x, pi
        integer, intent(in) :: num
        integer, dimension(num) :: sort_i
        real(r8), dimension(20, 2),intent(out) :: points
        real(r8), dimension(num, 2) :: points_i
        real(r8), dimension(2) :: center_point
        real(r8), dimension(num) :: angle

        center_point = 0.

        pi = 3.1415926535

        do i = 1, num, 1
            center_point = center_point + points(i, :)
            !print*,points(i,:)
            sort_i(i) = i
        end do
        center_point = center_point / num

        !print*,"center_point",center_point

        do i = 1, num, 1
            points_i(i, :) = points(i, :) - center_point
            if(points_i(i, 2) >= 0)then
                if(points_i(i, 1) == 0)then
                    angle(i) = pi / 2
                else
                    angle(i) = atan(points_i(i, 2) / points_i(i, 1))
                    if(points_i(i, 1) < 0)then
                        angle(i) = angle(i) + pi
                    end if
                end if
            else
                if(points_i(i, 1) == 0)then
                    angle(i) = 1.5 * pi
                else if(points_i(i, 1) < 0)then
                    angle(i) = atan(points_i(i, 2) / points_i(i, 1)) + pi
                else
                    angle(i) = atan(points_i(i, 2) / points_i(i, 1)) + 2 * pi
                end if
            end if
        end do
        do i = 1, num - 1, 1
            do j = i + 1, num, 1
                if(angle(j) < angle(i))then
                    angle_x = angle(j)
                    angle(j) = angle(i)
                    angle(i) = angle_x
                    x = sort_i(j)
                    sort_i(j) = sort_i(i)
                    sort_i(i) = x
                end if
            end do
        end do

        !print*,"angle2",angle

        do i = 1, num, 1
            points(i, :) = points_i(sort_i(i), :) + center_point
        end do

    END SUBROUTINE SortPoints


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
        real(r8), dimension(sjx_points*5, 4), intent(in) :: mp
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


    SUBROUTINE GetTriangleDis(ngrwm, lbx_points, ngrmw, sjx_points, mp_dis)

        implicit none

        integer :: i, j, k, w, num(5),l
        integer, intent(in) :: sjx_points, lbx_points
        integer, dimension(8, lbx_points), intent(in) :: ngrwm
        integer, dimension(4, sjx_points), intent(in) :: ngrmw
        integer, dimension(sjx_points, sjx_points), intent(out) :: mp_dis

        mp_dis = 9
        num = 0

        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i,j,k,w)
        do i = 1, sjx_points, 1
            mp_dis(i, i) = 0
            if(ngrmw(4, i) == 3)then
                do j = 1, 3, 1
                    w = ngrmw(j, i)
                    do k = 1, ngrwm(8, w), 1
                        if(mp_dis(i, ngrwm(k, w)) == 9)then
                            mp_dis(i, ngrwm(k, w)) = 1
                            mp_dis(ngrwm(k, w), i) = 1
                            num(1) = num(1) + 1
                        end if
                    end do
                end do
            end if
        end do
        !$OMP END PARALLEL DO
        !print*,num(1)
        !print*,sum(mp_dis)
        !num(1) = int((sjx_points*sjx_points*9-sjx_points*9-sum(mp_dis)) / 8)

        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i,j,k,w)
        do i = 1, sjx_points, 1
            do j = 1, sjx_points, 1
                if(mp_dis(i, j) == 1)then
                    do k = 1, 3, 1
                        w = ngrmw(k, j)
                        do l = 1, ngrwm(8, w), 1
                            if(mp_dis(i, ngrwm(l, w)) == 9)then
                                mp_dis(i, ngrwm(l, w)) = 2
                                mp_dis(ngrwm(l, w), i) = 2
                                num(2) = num(2) + 1
                            end if
                        end do
                    end do
                end if
            end do
        end do
        !$OMP END PARALLEL DO
        !print*,num(2)
        !print*,sum(mp_dis)

        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i,j,k,w)
        do i = 1, sjx_points, 1
            do j = 1, sjx_points, 1
                if(mp_dis(i, j) == 2)then
                    do k = 1, 3, 1
                        w = ngrmw(k, j)
                        do l = 1, ngrwm(8, w), 1
                            if(mp_dis(i, ngrwm(l, w)) == 9)then
                                mp_dis(i, ngrwm(l, w)) = 3
                                mp_dis(ngrwm(l, w), i) = 3
                                num(3) = num(3) + 1
                            end if
                        end do
                    end do
                end if
            end do
        end do
        !!$OMP END PARALLEL DO
        !print*,num(3)
        !print*,sum(mp_dis)

        !print*,"4"
        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i,j,k,w)
        do i = 1, sjx_points, 1
            do j = 1, sjx_points, 1
                if(mp_dis(i, j) == 3)then
                    do k = 1, 3, 1
                        w = ngrmw(k, j)
                        do l = 1, ngrwm(8, w), 1
                            if(mp_dis(i, ngrwm(l, w)) == 9)then
                                mp_dis(i, ngrwm(l, w)) = 4
                                mp_dis(ngrwm(l, w), i) = 4
                                num(4) = num(4) + 1
                            end if
                        end do
                    end do
                end if
            end do
        end do
        !$OMP END PARALLEL DO
        !print*,num(4)
        !print*,sum(mp_dis)

        !print*,"5"
        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i,j,k,w)
        do i = 1, sjx_points, 1
            do j = 1, sjx_points, 1
                if(mp_dis(i, j) == 4)then
                    do k = 1, 3, 1
                        w = ngrmw(k, j)
                        do l = 1, ngrwm(8, w), 1
                            if(mp_dis(i, ngrwm(l, w)) == 9)then
                                mp_dis(i, ngrwm(l, w)) = 5
                                mp_dis(ngrwm(l, w), i) = 5
                                num(5) = num(5) + 1
                            end if
                        end do
                    end do
                end if
            end do
        end do
        !$OMP END PARALLEL DO
        !print*,num(5)
        !print*,sum(mp_dis)

    END SUBROUTINE GetTriangleDis


    SUBROUTINE IsInRefineArea(IsInRfArea,sjx_points,lbx_points,ngrmw,wp)

      implicit none

      integer :: i,j,n
      integer,intent(in) :: sjx_points,lbx_points
      integer,dimension(4,sjx_points),intent(in) :: ngrmw
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


END module MOD_refine_lbx_step2
