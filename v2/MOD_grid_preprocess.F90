Module MOD_grid_preprocess
    USE NETCDF
    USE consts_coms
    USE refine_vars
    USE MOD_file_preprocess, only : Unstructured_Mesh_Read
    implicit none

    Contains

    SUBROUTINE grid_preprocess

        IMPLICIT NONE
        integer :: sjx_points, lbx_points, GXR_iter
        real(r8), allocatable :: mp(:, :), wp(:, :)
        integer,  allocatable :: ngrmw(:, :), ngrwm(:, :), n_ngrwm(:)
        character(LEN = 256) :: lndname, inputfile
        character(LEN = 5) :: nxpc, stepc, GXRC

        ! read unstructure mesh
        print*, "start to read unstructure mesh data in the Module MOD_grid_preprocess in Line 16"
        ! 读取未细化初始网格数据
        write(nxpc, '(I4.4)') NXP
        write(stepc, '(I2.2)') step
        lndname = trim(file_dir) // 'gridfile/gridfile_NXP' // trim(nxpc) // '_'//trim(stepc)// '_' // trim(mode_grid) // '.nc4'
        print*,lndname
        call Unstructured_Mesh_Read(lndname, sjx_points, lbx_points, mp, wp, ngrmw, ngrwm, n_ngrwm)
        print*, "The unstructured grid data reading have done "
        print*, ""
        print*, "In total, triangular mesh number: ", sjx_points, "polygon mesh number: ", lbx_points
        print*, ""

        GXR_iter = 0
        print*, "SpringAjustment_global start"
        call SpringAjustment_global(GXR_iter, sjx_points, lbx_points, ngrmw, ngrwm, n_ngrwm, mp, wp)
        print*, "SpringAjustment_global finish"

        if (GXR == 0) return
        ! 如果不为0，则进行循环，将网格分辨率加倍，目前是对全局，后面会改为只对局部网格生效
        write(GXRc, '(I1.1)') GXR_iter
        inputfile = trim(file_dir) // 'gridfile/gridfile_NXP' // trim(nxpc) // '_'//trim(stepc) // '_'//trim(GXRc)// '_' // trim(mode_grid) // '.nc4'
        CALL execute_command_line('cp '//trim(lndname)//' '//trim(inputfile))

        ! 这个应该是网格加倍之后采用的东西
        print*, "SpringAjustment_global start"
        call SpringAjustment_global(GXR_iter, sjx_points, lbx_points, ngrmw, ngrwm, n_ngrwm, mp, wp)
        print*, "SpringAjustment_global finish"

    END SUBROUTINE grid_preprocess

    ! 这里只根据三角形的角度判断是否要调整，是否需要根据多边形角度判断是否需要调整呢？
    SUBROUTINE SpringAjustment_global(GXR_iter, num_sjx, num_dbx, ngrmw_f, ngrwm_f, n_ngrwm_f, mp_f, wp_f)

        IMPLICIT NONE
        ! Extr_*** 指网格的最小角度与最大角度
        ! Eavg_*** 指网格的最小角度与最大角度平均值
        ! Savg_*** 指网格角度与正多边形角度的标准差
        ! less30   指三角形网格中小于30度角的数量占比
        integer,  intent(in) :: GXR_iter, num_sjx, num_dbx
        integer,  allocatable, intent(in) :: ngrmw_f(:, :), ngrwm_f(:, :), n_ngrwm_f(:)
        character(LEN = 256) :: lndname
        character(LEN = 5) :: stepc, nxpc
        integer :: w1, w2, w3, icl
        integer :: i, j, k, sa_iter, hhh(4) 
        integer :: lpDimID, twDimID, ncVarID(13), ncid
        integer :: num_wbx, num_lbx
        real(r8) :: rx, ry, fra, sjx(3, 2), mean_length_sjx
        real(r8) :: Extr_sjx_temp(2), Eavg_sjx_temp(2), Savg_sjx_temp, less30_temp
        real(r8) :: Extr_wbx_temp(2), Eavg_wbx_temp(2), Savg_wbx_temp
        real(r8) :: Extr_lbx_temp(2), Eavg_lbx_temp(2), Savg_lbx_temp
        logical :: End_SpringAjustment 
        integer,  allocatable :: adjust_sjx_flag(:), adjust_dbx_flag(:)
        real(r8), allocatable :: MoveDis(:, :)  ! 记录w点在x、y方向上的调整距离
        real(r8), allocatable :: fra_sjx(:, :)  ! 记录每次三角形顶点移动的幅度
        real(r8), allocatable :: length_sjx_temp(:, :) ! length_sjx(:, :, :)
        real(r8), allocatable :: angle_sjx(:, :, :), Extr_sjx(:, :), Eavg_sjx(:, :), Savg_sjx(:), less30(:)
        real(r8), allocatable :: angle_wbx(:, :, :), Extr_wbx(:, :), Eavg_wbx(:, :), Savg_wbx(:)
        real(r8), allocatable :: angle_lbx(:, :, :), Extr_lbx(:, :), Eavg_lbx(:, :), Savg_lbx(:)
        real(r8), allocatable :: angle_sjx_temp(:, :), angle_wbx_temp(:, :), angle_lbx_temp(:, :) 
        real(r8), allocatable, intent(inout) :: wp_f(:, :)
        real(r8), allocatable, intent(inout) :: mp_f(:, :)

        write(stepc, '(I2.2)') step
        write(nxpc,  '(I4.4)') nxp

        ! obtain num_wbx, num_lbx
        num_wbx = 0; num_lbx = 0;
        do i = 2, num_dbx, 1
            if (n_ngrwm_f(i) == 5) then
                num_wbx = num_wbx + 1
            else if (n_ngrwm_f(i) == 6) then
                num_lbx = num_lbx + 1
            else
                print*, "n_ngrwm_f(i) = ", n_ngrwm_f(i)
                stop "ERROR! n_ngrwm_f(i) must be 5 or 6"
            end if
        end do
        print*, "num_wbx = ", num_wbx
        print*, "num_lbx = ", num_lbx

        ! for sjx
        ! allocate(length_sjx(0:max_sa_iter, num_sjx, 3)); length_sjx(:, :, :) = 0.
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
        ! 这几个变量在这里赋值后，在while循环中会重复赋值，但是不会重新清零！
        Extr_sjx_temp = 0.; Eavg_sjx_temp = 0.; Savg_sjx_temp = 0.; less30_temp = 0.
        Extr_wbx_temp = 0.; Eavg_wbx_temp = 0.; Savg_wbx_temp = 0.;
        Extr_lbx_temp = 0.; Eavg_lbx_temp = 0.; Savg_lbx_temp = 0.;

        ! 弹性调整&计算每次调整后的三角形网格质量
        print*, "开始弹性调整"
        allocate(fra_sjx(num_sjx, 3))
        allocate(MoveDis(num_dbx, 2)) ! 记录w点在x、y方向上的调整距离
        allocate(adjust_sjx_flag(num_sjx));  adjust_sjx_flag = 1 ! 三角形初始化的时候都需要使用，不可跳过
        allocate(adjust_dbx_flag(num_dbx));  adjust_dbx_flag = 1 ! 多边形初始化的时候都需要使用，不可跳过
        sa_iter = 0
        k = -1 ! 初始化
        End_SpringAjustment = .false.
        hhh = [2,3,1,2] ! 用于简化计算
        
        do while(End_SpringAjustment .eqv. .false.)! 问题在于每次都要对全局进行调整，而且调整的方法是固定的
            End_SpringAjustment = .true.
            ! print*, "meshquality calculate start, sa_iter = ", sa_iter
            ! 每次调整后，再计算三角形，五六七边形的网格质量
            Call TriMeshQuality(    num_sjx, wp_f, ngrmw_f, adjust_sjx_flag, length_sjx_temp, angle_sjx_temp, Extr_sjx_temp, Eavg_sjx_temp, Savg_sjx_temp, less30_temp)! 三角形网格质量
            ! tri
            ! length_sjx(sa_iter, :, :) = length_sjx_temp
            angle_sjx(sa_iter, :, :)  = angle_sjx_temp
            ! print*, "meshquality calculate finish, sa_iter = ", sa_iter
            print*, "第", sa_iter, "次弹性调整完成，调整角度个数为", k
            if (mod(sa_iter, 10) == 0) then
                Call PolyMeshQuality(5, num_dbx, mp_f, ngrwm_f, n_ngrwm_f, adjust_dbx_flag, angle_wbx_temp, Extr_wbx_temp, Eavg_wbx_temp, Savg_wbx_temp)! 五边形网格质量
                Call PolyMeshQuality(6, num_dbx, mp_f, ngrwm_f, n_ngrwm_f, adjust_dbx_flag, angle_lbx_temp, Extr_lbx_temp, Eavg_lbx_temp, Savg_lbx_temp)! 六边形网格质量
    
                ! tri
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

                print*, "三角形边长：", minval(length_sjx_temp(2:num_sjx, :)), maxval(length_sjx_temp(2:num_sjx, :))
                print*, "三角形角度：", minval(angle_sjx_temp(2:num_sjx, :)), maxval(angle_sjx_temp(2:num_sjx, :))
                print*, "五边形角度：", minval(angle_wbx_temp), maxval(angle_wbx_temp)
                print*, "六边形角度：", minval(angle_lbx_temp), maxval(angle_lbx_temp)
                print*, ""
            end if
            if (GXR_iter == 0) cycle ! 直接跳出循环

            k = 0 ! 记录钝角三角形个数（需要调整的角度的个数）
            MoveDis = 0. ! 调整距离初始化为0
            adjust_sjx_flag = 0
            adjust_dbx_flag = 0
            fra_sjx = 0.
            ! New
            ! 根据三角形边长与平均边长的关系调整
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do i = 2, num_sjx, 1
                if (all(angle_sjx_temp(i, :) < Extr_sjx_GXR0(2))) cycle ! 跳过不存在钝角的三角形   
                if (all(angle_sjx_temp(i, :) > Extr_sjx_GXR0(1))) cycle ! 跳过不存在钝角的三角形 
                ! 找到最小角度小于40度，最大角度大于80度的三角形
                k = k + 1
                ! 获取三角形边长的平均长度
                mean_length_sjx = sum(length_sjx_temp(i, :))/3
                do j = 1, 3, 1
                    fra = (length_sjx_temp(i, j) - mean_length_sjx) / (length_sjx_temp(i, j) + mean_length_sjx) 
                    fra = fra * 0.1
                    fra_sjx(k, j) = fra
                    ! 获取三角形另外两个顶点的坐标，建议是逆时针获取
                    w1 = ngrmw_f(hhh(j),   i)
                    w2 = ngrmw_f(hhh(j+1), i)
                    rx = wp_f(w2, 1) - wp_f(w1, 1)
                    ry = wp_f(w2, 2) - wp_f(w1, 2)
                    call CheckLon(rx)! 这是一种距离调整方式，但是只能记录上一步的调整结果
                    adjust_dbx_flag(w1)= 1
                    adjust_dbx_flag(w2)= 1
                    MoveDis(w1, 1) = MoveDis(w1, 1) + rx * fra / 2.
                    MoveDis(w1, 2) = MoveDis(w1, 2) + ry * fra / 2.
                    MoveDis(w2, 1) = MoveDis(w2, 1) - rx * fra / 2.
                    MoveDis(w2, 2) = MoveDis(w2, 2) - ry * fra / 2.
                end do
            end do

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            print*, "new scheme, k = ", k
            if (k == 0) cycle
            print*, "max(fra_sjx) = ", maxval(fra_sjx(1:k, 1)), maxval(fra_sjx(1:k, 2)), maxval(fra_sjx(1:k, 3))
            print*, "min(fra_sjx) = ", minval(fra_sjx(1:k, 1)), minval(fra_sjx(1:k, 2)), minval(fra_sjx(1:k, 3))

            ! 全部弹性调整后，再调整三角形顶点的经纬度
            do i = 2, num_dbx, 1
                if (adjust_dbx_flag(i) == 0) cycle ! jump out if not change
                wp_f(i, 1:2) = wp_f(i, 1:2) + MoveDis(i, 1:2)
                call CheckLon(wp_f(i, 1))
            end do

            ! 调整所有m点至三角形网格重心，但是应该不是所有三角形都需要调整的
            ! do not use all range    
            do i = 2, num_sjx, 1
                ! 判断三角形的顶点是否涉及这些多边形
                if (sum(adjust_dbx_flag(ngrmw_f(:, i))) == 0) cycle
                adjust_sjx_flag(i) = 1 
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
            ! print*, "SpringAdjust finish sa_iter = ", sa_iter
            ! 不出现钝角或者达到最大调整次数的时候就退出while循环！
            if ((k /= 0) .and. (sa_iter < max_sa_iter)) End_SpringAjustment = .false.
            sa_iter = sa_iter + 1
        end do ! while(End_SpringAjustment .eqv. .false.)

        deallocate(MoveDis, length_sjx_temp)
        deallocate(angle_sjx, angle_wbx, angle_lbx)
        deallocate(angle_sjx_temp, angle_wbx_temp, angle_lbx_temp)
        
        lndname = trim(file_dir) // "result/quality_NXP" // trim(nxpc) // '_' // trim(stepc) // "_global.nc4"
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
        CALL CHECK(NF90_DEF_VAR(ncID, "Savg_wbx", NF90_FLOAT, (/ lpDimID /), ncVarID(7)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Savg_lbx", NF90_FLOAT, (/ lpDimID /), ncVarID(8)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Extr_wbx", NF90_FLOAT, (/ lpDimID, twDimID /), ncVarID(9)))
        CALL CHECK(NF90_DEF_VAR(ncID, "Extr_lbx", NF90_FLOAT, (/ lpDimID, twDimID /), ncVarID(10)))

        CALL CHECK(NF90_ENDDEF(ncID))

        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1),  less30(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(2),  Eavg_sjx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(3),  Savg_sjx(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(4),  Extr_sjx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5),  Eavg_wbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6),  Eavg_lbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(7),  Savg_wbx(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(8),  Savg_lbx(0:max_sa_iter)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(9), Extr_wbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(10), Extr_lbx(0:max_sa_iter, 1:2)))
        CALL CHECK(NF90_CLOSE(ncID))
        deallocate(Extr_sjx, Eavg_sjx, Savg_sjx, less30)
        deallocate(Extr_wbx, Eavg_wbx, Savg_wbx)
        deallocate(Extr_lbx, Eavg_lbx, Savg_lbx)
        print*, "弹性调整结束"
        if (GXR_iter == 0) Extr_sjx_GXR0 = Extr_sjx_temp

    END SUBROUTINE SpringAjustment_global

    SUBROUTINE SpringAjustment_refine(sjx_points, num_sjx, num_dbx, num_ref, ngrmw_f, ngrwm_f, n_ngrwm_f, mp_f, wp_f)

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
        real(r8) :: rx, ry, fra, sjx(3, 2), mean_length_sjx
        real(r8) :: Extr_sjx_temp(2), Eavg_sjx_temp(2), Savg_sjx_temp, less30_temp
        real(r8) :: Extr_wbx_temp(2), Eavg_wbx_temp(2), Savg_wbx_temp
        real(r8) :: Extr_lbx_temp(2), Eavg_lbx_temp(2), Savg_lbx_temp
        real(r8) :: Extr_qbx_temp(2), Eavg_qbx_temp(2), Savg_qbx_temp
        logical :: End_SpringAjustment 
        integer,  allocatable :: adjust_sjx_flag(:), adjust_dbx_flag(:)
        real(r8), allocatable :: MoveDis(:, :)  ! 记录w点在x、y方向上的调整距离
        real(r8), allocatable :: fra_sjx(:, :)  ! 记录每次三角形顶点移动的幅度
        real(r8), allocatable :: length_sjx_temp(:, :) ! length_sjx(:, :, :)
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
        if (num_dbx-1 /= num_wbx + num_lbx + num_qbx) print*, "存在小于五边形的多边形"

        ! for sjx
        ! allocate(length_sjx(0:max_sa_iter, num_sjx, 3)); length_sjx(:, :, :) = 0.
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
        if (num_qbx /= 0) then
            allocate(angle_qbx(0:max_sa_iter, num_qbx, 7));  angle_qbx(:, :, :)  = 0.
            allocate(Eavg_qbx(0:max_sa_iter, 2)); Extr_qbx = 0.
            allocate(Extr_qbx(0:max_sa_iter, 2)); Eavg_qbx = 0.
            allocate(Savg_qbx(0:max_sa_iter));    Savg_qbx = 0.
            allocate(angle_qbx_temp(num_qbx, 7));           angle_qbx_temp(:, :)  = 0.
        end if
        ! 这几个变量在这里赋值后，在while循环中会重复赋值，但是不会重新清零！
        Extr_sjx_temp = 0.; Eavg_sjx_temp = 0.; Savg_sjx_temp = 0.; less30_temp = 0.
        Extr_wbx_temp = 0.; Eavg_wbx_temp = 0.; Savg_wbx_temp = 0.;
        Extr_lbx_temp = 0.; Eavg_lbx_temp = 0.; Savg_lbx_temp = 0.;
        Extr_qbx_temp = 0.; Eavg_qbx_temp = 0.; Savg_qbx_temp = 0.;

        ! 弹性调整&计算每次调整后的三角形网格质量
        print*, "开始弹性调整"
        allocate(fra_sjx(num_sjx-(sjx_points - num_ref)+1, 3))
        allocate(MoveDis(num_dbx, 2)) ! 记录w点在x、y方向上的调整距离
        allocate(adjust_sjx_flag(num_sjx));  adjust_sjx_flag = 1 ! 三角形初始化的时候都需要使用，不可跳过
        allocate(adjust_dbx_flag(num_dbx));  adjust_dbx_flag = 1 ! 多边形初始化的时候都需要使用，不可跳过
        sa_iter = 0
        k = -1 ! 初始化
        End_SpringAjustment = .false.
        hhh = [2,3,1,2] ! 用于简化计算
        
        do while(End_SpringAjustment .eqv. .false.)! 问题在于每次都要对全局进行调整，而且调整的方法是固定的
            End_SpringAjustment = .true.
            ! print*, "meshquality calculate start, sa_iter = ", sa_iter
            ! 每次调整后，再计算三角形，五六七边形的网格质量
            adjust_sjx_flag = 1
            Call TriMeshQuality(    num_sjx, wp_f, ngrmw_f, adjust_sjx_flag, length_sjx_temp, angle_sjx_temp, Extr_sjx_temp, Eavg_sjx_temp, Savg_sjx_temp, less30_temp)! 三角形网格质量
            ! tri
            ! length_sjx(sa_iter, :, :) = length_sjx_temp
            angle_sjx(sa_iter, :, :)  = angle_sjx_temp
            ! print*, "meshquality calculate finish, sa_iter = ", sa_iter
            print*, "第", sa_iter, "次弹性调整完成，调整角度个数为", k
            if (mod(sa_iter, 10) == 0) then
                adjust_dbx_flag = 1
                Call PolyMeshQuality(5, num_dbx, mp_f, ngrwm_f, n_ngrwm_f, adjust_dbx_flag, angle_wbx_temp, Extr_wbx_temp, Eavg_wbx_temp, Savg_wbx_temp)! 五边形网格质量
                Call PolyMeshQuality(6, num_dbx, mp_f, ngrwm_f, n_ngrwm_f, adjust_dbx_flag, angle_lbx_temp, Extr_lbx_temp, Eavg_lbx_temp, Savg_lbx_temp)! 六边形网格质量
                CALL PolyMeshQuality(7, num_dbx, mp_f, ngrwm_f, n_ngrwm_f, adjust_dbx_flag, angle_qbx_temp, Extr_qbx_temp, Eavg_qbx_temp, Savg_qbx_temp)! 七边形网格质量
                ! tri
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

                print*, "三角形边长：", minval(length_sjx_temp(2:num_sjx, :)), maxval(length_sjx_temp(2:num_sjx, :))
                print*, "三角形角度：", minval(angle_sjx_temp(2:num_sjx, :)), maxval(angle_sjx_temp(2:num_sjx, :))
                print*, "五边形角度：", minval(angle_wbx_temp), maxval(angle_wbx_temp)
                print*, "六边形角度：", minval(angle_lbx_temp), maxval(angle_lbx_temp)
                print*, "七边形角度：", minval(angle_qbx_temp), maxval(angle_qbx_temp)
                print*, ""
            end if

            ! print*, "springAdjust start, sa_iter = ", sa_iter
            if (.true.) then
                k = 0 ! 记录钝角三角形个数（需要调整的角度的个数）
                MoveDis = 0. ! 调整距离初始化为0
                adjust_sjx_flag = 1
                adjust_dbx_flag = 1
                fra_sjx = 0.
                ! New
                ! 根据三角形边长与平均边长的关系调整
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                do i = sjx_points - num_ref + 1, num_sjx, 1
                    if (all(angle_sjx_temp(i, :) < 75.)) cycle ! 跳过不存在钝角的三角形   
                    if (all(angle_sjx_temp(i, :) > 45.)) cycle ! 跳过不存在钝角的三角形 
                    ! 找到最小角度小于45度，最大角度大于75度的三角形
                    k = k + 1
                    ! 获取三角形边长的平均长度
                    mean_length_sjx = sum(length_sjx_temp(i, :))/3
                    do j = 1, 3, 1
                        fra = (length_sjx_temp(i, j) - mean_length_sjx) / (length_sjx_temp(i, j) + mean_length_sjx) 
                        fra = fra * 0.1
                        fra_sjx(k, j) = fra
                        ! 获取三角形另外两个顶点的坐标，建议是逆时针获取
                        w1 = ngrmw_f(hhh(j),   i)
                        w2 = ngrmw_f(hhh(j+1), i)
                        rx = wp_f(w2, 1) - wp_f(w1, 1)
                        ry = wp_f(w2, 2) - wp_f(w1, 2)
                        call CheckLon(rx)! 这是一种距离调整方式，但是只能记录上一步的调整结果
                        adjust_dbx_flag(w1) = 1
                        adjust_dbx_flag(w2) = 1
                        MoveDis(w1, 1) = MoveDis(w1, 1) + rx * fra / 2.
                        MoveDis(w1, 2) = MoveDis(w1, 2) + ry * fra / 2.
                        MoveDis(w2, 1) = MoveDis(w2, 1) - rx * fra / 2.
                        MoveDis(w2, 2) = MoveDis(w2, 2) - ry * fra / 2.
                    end do
                end do

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                print*, "new scheme, k = ", k
                if (k == 0) cycle
                print*, "max(fra_sjx) = ", maxval(fra_sjx(1:k, 1)), maxval(fra_sjx(1:k, 2)), maxval(fra_sjx(1:k, 3))
                print*, "min(fra_sjx) = ", minval(fra_sjx(1:k, 1)), minval(fra_sjx(1:k, 2)), minval(fra_sjx(1:k, 3))
            end if

            print*, ""
            ! 全部弹性调整后，再调整三角形顶点的经纬度
            do i = 2, num_dbx, 1
                if (adjust_dbx_flag(i) == 0) cycle ! jump out if not change
                wp_f(i, 1:2) = wp_f(i, 1:2) + MoveDis(i, 1:2)
                call CheckLon(wp_f(i, 1))
            end do

            ! 调整所有m点至三角形网格重心，但是应该不是所有三角形都需要调整的
            ! do not use all range    
            do i = 2, num_sjx, 1
                ! 判断三角形的顶点是否涉及这些多边形
                if (sum(adjust_dbx_flag(ngrmw_f(:, i))) == 0) cycle
                adjust_sjx_flag(i) = 1 
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
            ! print*, "SpringAdjust finish sa_iter = ", sa_iter
            ! 不出现钝角或者达到最大调整次数的时候就退出while循环！
            if ((k /= 0) .and. (sa_iter < max_sa_iter)) End_SpringAjustment = .false.
            sa_iter = sa_iter + 1
        end do ! while(End_SpringAjustment .eqv. .false.)
        if ((k == 0) .and. (sa_iter < max_sa_iter))  print*, "没有需要调整的三角形，退出"
        if ((k /= 0) .and. (sa_iter >= max_sa_iter)) print*, "达到角度调整次数上限，退出"

        lndname = trim(file_dir) // "result/quality_NXP" // trim(nxpc) // '_' // trim(stepc) // "_refine.nc4"
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
        print*, "最后的网格质量如下："
        adjust_sjx_flag = 1
        adjust_dbx_flag = 1
        Call TriMeshQuality(    num_sjx, wp_f, ngrmw_f, adjust_sjx_flag, length_sjx_temp, angle_sjx_temp, Extr_sjx_temp, Eavg_sjx_temp, Savg_sjx_temp, less30_temp)! 三角形网格质量            
        Call PolyMeshQuality(5, num_dbx, mp_f, ngrwm_f, n_ngrwm_f, adjust_dbx_flag, angle_wbx_temp, Extr_wbx_temp, Eavg_wbx_temp, Savg_wbx_temp)! 五边形网格质量
        Call PolyMeshQuality(6, num_dbx, mp_f, ngrwm_f, n_ngrwm_f, adjust_dbx_flag, angle_lbx_temp, Extr_lbx_temp, Eavg_lbx_temp, Savg_lbx_temp)! 六边形网格质量
        Call PolyMeshQuality(7, num_dbx, mp_f, ngrwm_f, n_ngrwm_f, adjust_dbx_flag, angle_qbx_temp, Extr_qbx_temp, Eavg_qbx_temp, Savg_qbx_temp)! 七边形网格质量 
        print*, "三角形边长：", minval(length_sjx_temp(2:num_sjx, :)), maxval(length_sjx_temp(2:num_sjx, :))
        print*, "三角形角度：", minval(angle_sjx_temp(2:num_sjx, :)), maxval(angle_sjx_temp(2:num_sjx, :))
        print*, "五边形角度：", minval(angle_wbx_temp), maxval(angle_wbx_temp)
        print*, "六边形角度：", minval(angle_lbx_temp), maxval(angle_lbx_temp)
        print*, "七边形角度：", minval(angle_qbx_temp), maxval(angle_qbx_temp)
        print*, ""

        deallocate(MoveDis, length_sjx_temp)
        deallocate(angle_sjx, angle_wbx, angle_lbx)
        deallocate(angle_sjx_temp, angle_wbx_temp, angle_lbx_temp)
        deallocate(angle_qbx, angle_qbx_temp)

    END SUBROUTINE SpringAjustment_refine

    
    SUBROUTINE TriMeshQuality(num_sjx, wp_f, ngrmw_f, adjust_sjx_flag, length_sjx_temp, angle_sjx_temp, Extr_sjx_temp, Eavg_sjx_temp, Savg_sjx_temp, less30_temp)

        implicit none ! 定义的顺序最好是只读入的参数，内部参数，需要传出的参数
        integer, intent(in) :: num_sjx
        real(r8), dimension(:, :), intent(in) :: wp_f
        integer,  dimension(:, :), intent(in) :: ngrmw_f
        integer,  dimension(:), intent(in) :: adjust_sjx_flag
        integer :: i
        real(r8) :: length_temp(3), angle_temp(3), sjx(3, 2) !angle_regular
        real(r8), dimension(:, :), intent(inout) :: length_sjx_temp, angle_sjx_temp
        real(r8), intent(inout) :: Extr_sjx_temp(2), Eavg_sjx_temp(2), Savg_sjx_temp, less30_temp
        ! 注意跨越180的情况
        do i = 2, num_sjx, 1 ! 从2开始，因为第一个不存在
            if (adjust_sjx_flag(i) /= 0) then ! 赋值为1的时候才需要更新数据
                sjx = wp_f(ngrmw_f(:, i), 1:2)
                if (maxval(sjx(:, 1)) - minval(sjx(:, 1)) > 180.) call CheckCrossing(3, sjx)
                CALL GetTriangleLength(sjx, length_temp)
                CALL GetAngle(3, sjx, angle_temp)! 计算多边形内角和
                length_sjx_temp(i, :) = length_temp 
                angle_sjx_temp(i, :) = angle_temp
            else
                length_temp = length_sjx_temp(i, :)
                angle_temp  = angle_sjx_temp(i, :)
            end if
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

    SUBROUTINE PolyMeshQuality(num_edges, num_dbx, mp_f, ngrwm_f, n_ngrwm_f, adjust_dbx_flag, angle_dbx_temp, Extr_dbx_temp, Eavg_dbx_temp, Savg_dbx_temp)

        implicit none ! 定义的顺序最好是只读入的参数，内部参数，需要传出的参数
        integer, intent(in) :: num_edges, num_dbx
        real(r8), dimension(:, :), intent(in) :: mp_f
        integer,  dimension(:, :), intent(in) :: ngrwm_f
        integer,  dimension(:), allocatable, intent(in) :: n_ngrwm_f
        integer,  dimension(:), allocatable, intent(in) :: adjust_dbx_flag
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
            if (adjust_dbx_flag(i) /= 0) then ! 赋值为1的时候才需要更新数据
                dbx = mp_f(ngrwm_f(1:num_edges, i), :)
                if (maxval(dbx(:, 1)) - minval(dbx(:, 1)) > 180.) call CheckCrossing(num_edges, dbx)
                CALL GetAngle(num_edges, dbx, angle_temp)! 计算多边形内角和
                angle_dbx_temp(j, :) = angle_temp
            else ! jump out if not change
                angle_temp = angle_dbx_temp(j, :)
            end if
            Eavg_dbx_temp(1) = Eavg_dbx_temp(1) + minval(angle_temp)
            Eavg_dbx_temp(2) = Eavg_dbx_temp(2) + maxval(angle_temp)
            Savg_dbx_temp    = Savg_dbx_temp    + sum((angle_temp - angle_regular)**2)
        end do
        Extr_dbx_temp(1) = minval(angle_dbx_temp)
        Extr_dbx_temp(2) = maxval(angle_dbx_temp)
        Eavg_dbx_temp(:) = Eavg_dbx_temp(:) / j
        Savg_dbx_temp    = sqrt(Savg_dbx_temp / (num_edges * j))

    END SUBROUTINE PolyMeshQuality   

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
        ! 半正矢公式, 以弧度制度量!! 中距离：适用于几百到几千公里，
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

END Module MOD_grid_preprocess
