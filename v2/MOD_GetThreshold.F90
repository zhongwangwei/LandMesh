MODULE MOD_GetThreshold
  
    USE consts_coms, only: r8, mode, nxp, file_dir, openmp, maxlc, nlons_source, nlats_source
    USE refine_vars 
    USE netcdf
    USE MOD_data_preprocess, only: landtypes ! 土地覆盖类型编号0为海洋
    USE MOD_file_preprocess, only: Contain_Read, Unstructured_Threshold_Read
    USE MOD_GetContain,      only: IsInRfArea_sjx, patchtypes
    USE MOD_Threshold_Read,  only: RL_onelayer, RL_twolayer, input2d, input3d

    implicit none
    integer, dimension(:),    allocatable, public :: ref_sjx ! 需要细化的三角形标记 1和0
    integer, dimension(:, :), allocatable, public :: ref_th
    type :: var_data2d
        real(r8), allocatable :: varms(:)
    end type
    type(var_data2d), allocatable, public :: output2d(:)
    type :: var_data3d
        real(r8), allocatable :: varms(:, :)
    end type
    type(var_data3d), allocatable, public :: output3d(:)
    integer, public :: num_swithes
    
    !interface mean_std_cal
    !    module procedure mean_std_cal2d
    !    module procedure mean_std_cal3d
    !end interface mean_std_cal


    contains

    SUBROUTINE GetThreshold(exit_loop)
        logical, intent(inout) :: exit_loop ! use for exit refine
        integer :: sjx_points, numpatch, spDimID, twoDimID, DimID, ncid, varid(18), flag_temps
        integer,  dimension(:, :), allocatable :: mp_id, mp_ii
        integer,  dimension(:),    allocatable :: n_landtypes, nlaa, p_num
        real(r8), dimension(:),  allocatable :: f_mainarea
        real(r8), allocatable :: var2d_temp(:, :), var3d_temp(:, :, :)
        character(LEN = 5) :: nxpc, stepc
        character(LEN = 256) :: lndname
        integer :: i, j, row, col, L, num_onelayer, num_twolayer
        
        if (th_file_read) then
            lndname = th_filedir
            CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid))! 1. NF90_OPEN 打开文件
            CALL CHECK(NF90_INQ_DIMID(ncid, "sjx_points", spDimID))!
            CALL CHECK(NF90_INQ_DIMID(ncid, "num_swithes", DimID))!
            CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, spDimID, len = sjx_points))
            CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, DimID, len = num_swithes))
            allocate(ref_sjx(sjx_points)); ref_sjx = 0
            allocate(ref_th(sjx_points, num_swithes)); ref_th = 0
            CALL CHECK(NF90_INQ_VARID(ncid, 'ref_th', varid(1)))
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), ref_th))
            CALL CHECK(NF90_CLOSE(ncid))
            do i = 1, sjx_points, 1
                if (sum(ref_th(i, :)) > 0) ref_sjx(i) = 1 
            end do
            print*, "需要细化的三角形个数: ",INT(sum(ref_sjx))
            deallocate(ref_th)
            return
        end if

        num_swithes = 0
        write(nxpc, '(I4.4)') NXP
        write(stepc, '(I2.2)') step
        print*, "开始读取非结构网格数据 in the GetThreshold.F90"
        lndname = trim(file_dir) // 'contain/contain_refine_NXP' // trim(nxpc) //'_'// trim(stepc) //'_mp.nc4'
        print*, lndname
        call Contain_Read(lndname, sjx_points, numpatch, mp_id, mp_ii) ! 保存文件
        print*, "非结构网格数据读取完成 in the GetThreshold.F90"
        print*, ""

        allocate(ref_sjx(sjx_points)); ref_sjx = 0 ! will deallocate in the MOD_refine_lbx.F90
        allocate(ref_th(sjx_points,2+size(refine_onelayer)+size(refine_twolayer))); ref_th = 0

        ! 阈值文件的计算(可以保证进来GetThreshold就一定可以有阈值计算)
        if (refine_num_landtypes) then
            print*, "threshold_Calculation for refine_num_landtypes"
            allocate(n_landtypes(sjx_points)); n_landtypes = 0
            allocate(nlaa(0:maxlc)) ! zero is ocean
            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
            !$OMP PRIVATE(i, j, row, col, L, nlaa)
            do i = 2, sjx_points, 1
                ! 当三角形网格不位于细化区域内，or number_select is all in the sea跳过
                if (IsInRfArea_sjx(i) .eqv. .false.) cycle
                nlaa = 0  ! 单个非结构网格是否包含该种土地类型的标志（是一个数组），初值全部为0
                do j = 0, mp_id(i, 1)-1, 1
                    row = mp_ii(mp_id(i, 2)+j, 1) ! 获取非结构网格包含经纬度网格的经度序号
                    col = mp_ii(mp_id(i, 2)+j, 2) ! 获取非结构网格包含经纬度网格的纬度序号
                    L = landtypes(row, col) ! 获取这个经纬度网格的土地类型的编号 not ocean
                    if (L /= maxlc) nlaa(L) = 1
                end do
                n_landtypes(i) = sum(nlaa) ! 获取单个非结构网格中含有的土地类型数量总和
                if (n_landtypes(i) > th_num_landtypes) ref_sjx(i) = 1
            end do
            !$OMP END PARALLEL DO
            print*, "n_landtypes",minval(n_landtypes), maxval(n_landtypes)
            num_swithes  = num_swithes + 1
            ref_th(:, num_swithes) = ref_sjx
            deallocate(nlaa)
            print*, "threshold_Calculation for refine_num_landtypes : done"
        end if

        if (refine_area_mainland) then
            print*, "threshold_Calculation for refine_area_mainland"
            allocate(f_mainarea(sjx_points)); f_mainarea = 0.
            allocate(nlaa(0:maxlc)) ! zero is ocean
            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
            !$OMP PRIVATE(i, j, row, col, L, nlaa)
            do i = 2, sjx_points, 1
                if (IsInRfArea_sjx(i) .eqv. .false.) cycle
                nlaa = 0. ! 单个非结构网格内各土地类型面积，（是一个数组），初值为0
                do j = 0, mp_id(i, 1) - 1, 1
                    row = mp_ii(mp_id(i, 2) + j, 1) ! 获取非结构网格包含经纬度网格的经度序号
                    col = mp_ii(mp_id(i, 2) + j, 2) ! 获取非结构网格包含经纬度网格的纬度序号
                    L = landtypes(row, col) ! 获取这个经纬度网格的土地类型的编号 not ocean
                    if (L /= maxlc) nlaa(L) = nlaa(L) + 1 ! 该类型的土地面积累加
                end do
                f_mainarea(i) = min( maxval(nlaa)/mp_id(i, 1), 1) ! 非结构网格内主要土地类型占比
                if (f_mainarea(i) < th_area_mainland) ref_sjx(i) = 1
            end do
            !$OMP END PARALLEL DO
            print*, "f_mainarea",minval(f_mainarea), maxval(f_mainarea)
            num_swithes  = num_swithes + 1
            ref_th(:, num_swithes) = ref_sjx
            deallocate(nlaa)
            print*, "threshold_Calculation for refine_area_mainland : done"
        end if

        ! 改写为循环的形式， 而且先计算均值，再计算标准差, 单层与双层分开
        if (any(refine_onelayer .eqv. .true.) .or. &
            any(refine_twolayer .eqv. .true.)) then

            print*, "threshold_Calculation for refine_onelayer or refine_twolayer"
            allocate(p_num(sjx_points))

            ! refine_onelayer
            if (any(refine_onelayer .eqv. .true.)) then
                num_onelayer = 0 ! 单层数据计数
                allocate(output2d(size(refine_onelayer)/2))
                allocate(var2d_temp(nlons_source, nlats_source)); var2d_temp = 0.
                do i = 1, size(refine_onelayer)/2, 1 !
                    if ((refine_onelayer(2*i-1) .eqv. .false.) .and. (refine_onelayer(2*i) .eqv. .false.)) cycle ! 跳过
                    p_num = 0
                    flag_temps = 0
                    if (refine_onelayer(2*i)) flag_temps = 1
                    var2d_temp = input2d(i)%var2d
                    CALL mean_std_cal2d(i, sjx_points, flag_temps, mp_id, mp_ii, var2d_temp, num_onelayer, p_num)
                end do
                deallocate(var2d_temp)
            end if

            ! refine_twolayer
            if (any(refine_twolayer .eqv. .true.)) then
                num_twolayer = 0 ! 单层数据计数
                allocate(output3d(size(refine_twolayer)/2))
                allocate(var3d_temp(2, nlons_source, nlats_source)); var3d_temp = 0.
                do i = 1, size(refine_twolayer)/2, 1 !
                    if ((refine_twolayer(2*i-1) .eqv. .false.) .and. (refine_twolayer(2*i) .eqv. .false.)) cycle ! 跳过
                    p_num = 0
                    flag_temps = 0
                    if (refine_twolayer(2*i)) flag_temps = 1
                    var3d_temp = input3d(i)%var3d
                    CALL mean_std_cal3d(i, sjx_points, flag_temps, mp_id, mp_ii, var3d_temp, num_twolayer, p_num)
                end do
                deallocate(var3d_temp)
            end if

            print*, "threshold_Calculation for refine_onelayer or refine_twolayer : done"
        end if
        print*, "threshold_Calculation finish"

        ! 阈值计算结果的输出
        lndname = trim(file_dir) // "threshold/threshold_NXP" // trim(nxpc) //'_'// trim(stepc) // ".nc4"
        print*, lndname
        CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber,nf90_netcdf4), ncid))
        CALL CHECK(NF90_DEF_DIM(ncid, "sjx_points", sjx_points, spDimID))
        CALL CHECK(NF90_DEF_DIM(ncid, "dima", 2, twoDimID))
        CALL CHECK(NF90_DEF_DIM(ncid, "num_swithes", num_swithes, DimID))
        num_swithes = 0
        ! NF90_DEF_VAR
        if (refine_num_landtypes) then
            num_swithes = num_swithes + 1
            CALL CHECK(NF90_DEF_VAR(ncid, "n_landtypes", NF90_INT, (/ spDimID /), varid(num_swithes)))
        end if

        if (refine_area_mainland) then
            num_swithes = num_swithes + 1
            CALL CHECK(NF90_DEF_VAR(ncid, "f_mainarea", NF90_INT, (/ spDimID /), varid(num_swithes)))
        end if

        if (any(refine_onelayer .eqv. .true.) .or. &
            any(refine_twolayer .eqv. .true.)) then
            do i = 1, size(refine_onelayer)/2, 1
                if (refine_onelayer(2*i-1)) then
                    num_swithes = num_swithes + 1
                    CALL CHECK(NF90_DEF_VAR(ncid, trim(RL_onelayer(i))//"_m", NF90_float , (/ spDimID /), varid(num_swithes)))
                end if
                if (refine_onelayer(2*i)) then
                    num_swithes = num_swithes + 1
                    CALL CHECK(NF90_DEF_VAR(ncid, trim(RL_onelayer(i))//"_s", NF90_float , (/ spDimID /), varid(num_swithes)))
                end if
            end do
            do i = 1, size(refine_twolayer)/2, 1
                if (refine_twolayer(2*i-1)) then
                    num_swithes = num_swithes + 1
                    CALL CHECK(NF90_DEF_VAR(ncid, trim(RL_twolayer(i))//"_m", NF90_float , (/ spDimID, twoDimID /), varid(num_swithes)))
                end if
                if (refine_twolayer(2*i)) then
                    num_swithes = num_swithes + 1
                    CALL CHECK(NF90_DEF_VAR(ncid, trim(RL_twolayer(i))//"_s", NF90_float , (/ spDimID, twoDimID /), varid(num_swithes)))
                end if 
            end do
            num_swithes = num_swithes + 1
            CALL CHECK(NF90_DEF_VAR(ncid, "p_num", NF90_INT, (/ spDimID /), varid(num_swithes)))
        end if

        CALL CHECK(NF90_DEF_VAR(ncid, "ref_th", NF90_INT, (/ spDimID, DimID/), varid(num_swithes+1)))
        CALL CHECK(NF90_ENDDEF(ncid))

        ! NF90_PUT_VAR
        num_swithes = 0
        num_onelayer = 0 ! 单层数据计数
        num_twolayer = 0 ! 双层数据计数
        if (refine_num_landtypes) then
            num_swithes = num_swithes + 1
            CALL CHECK(NF90_PUT_VAR(ncid, varid(num_swithes), n_landtypes))
            deallocate(n_landtypes)
        end if

        if (refine_area_mainland) then
            num_swithes = num_swithes + 1
            CALL CHECK(NF90_PUT_VAR(ncid, varid(num_swithes), f_mainarea))
            deallocate(f_mainarea)
        end if

        if (any(refine_onelayer .eqv. .true.) .or. &
            any(refine_twolayer .eqv. .true.)) then
            do i = 1, size(refine_onelayer), 1
                if (refine_onelayer(i)) then
                    num_swithes = num_swithes + 1
                    num_onelayer = num_onelayer + 1
                    CALL CHECK(NF90_PUT_VAR(ncid, varid(num_swithes), output2d(num_onelayer)%varms))
                end if
            end do

            do i = 1, size(refine_twolayer), 1
                if (refine_twolayer(i)) then
                    num_swithes = num_swithes + 1
                    num_twolayer = num_twolayer + 1
                    CALL CHECK(NF90_PUT_VAR(ncid, varid(num_swithes), output3d(num_twolayer)%varms))
                end if
            end do
            num_swithes = num_swithes + 1
            CALL CHECK(NF90_PUT_VAR(ncid, varid(num_swithes), p_num))
            deallocate(p_num)
        end if

        CALL CHECK(NF90_PUT_VAR(ncid, varid(num_swithes+1), ref_th(:, 1:num_swithes-1)))
        CALL CHECK(NF90_CLOSE(ncid))

        if (allocated(output2d)) deallocate(output2d)
        if (allocated(output3d)) deallocate(output3d)
        deallocate(IsInRfArea_sjx) ! allocate in the MOD_GetContain.F90
        deallocate(ref_th)
        print*, "需要细化的三角形个数: ",INT(sum(ref_sjx))
        if (INT(sum(ref_sjx)) == 0) then
            exit_loop = .true.
            return
        end if

    END SUBROUTINE GetThreshold

    SUBROUTINE mean_std_cal2d(ii, sjx_points, flag_temps, mp_id, mp_ii, var2d, num_onelayer, p_num)
        integer, intent(in) :: ii, sjx_points, flag_temps
        integer,  dimension(:, :), intent(in) :: mp_id, mp_ii
        real(r8), dimension(:, :), intent(in) :: var2d
        
        real(r8), dimension(:),   allocatable :: tempm, temps! 均值计算/标准差计算的中间变量
        integer :: i, j, row, col, L
        integer,  intent(inout) :: num_onelayer
        integer,  dimension(:), intent(inout) :: p_num

        allocate(tempm(sjx_points)); tempm = 0.
        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i, j, row, col, L)
        do i = 2, sjx_points, 1
            if (IsInRfArea_sjx(i) .eqv. .false.) cycle
            do j = 0, mp_id(i, 1)-1, 1
                row = mp_ii(mp_id(i, 2)+j, 1) ! 获取非结构网格包含经纬度网格的经度序号
                col = mp_ii(mp_id(i, 2)+j, 2) ! 获取非结构网格包含经纬度网格的纬度序号
                if (L /= maxlc) then 
                    tempm(i) = tempm(i) + var2d(row, col)
                    p_num(i) = p_num(i) + 1
                end if
            end do
            tempm(i) = tempm(i) / p_num(i) !求均值
        end do
        !$OMP END PARALLEL DO

        if (refine_onelayer(2*ii-1)) then
            print*, trim(RL_onelayer(ii))//"_m", minval(tempm, mask=(tempm /= 0.)), maxval(tempm, mask=(tempm /= 0.)) 
            num_onelayer = num_onelayer + 1
            allocate(output2d(num_onelayer)%varms(sjx_points))
            output2d(num_onelayer)%varms = tempm
            do i = 2, sjx_points, 1
                if (tempm(i) > th_onelayer(2*ii-1)) ref_sjx(i) = 1
            end do
            num_swithes  = num_swithes + 1
            ref_th(:, num_swithes) = ref_sjx
        end if

        if (flag_temps == 0) then
            deallocate(tempm)
            return
        end if

        allocate(temps(sjx_points)); temps = 0.
        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i, j, row, col, L)
        do i = 2, sjx_points, 1
            if (IsInRfArea_sjx(i) .eqv. .false.) cycle
            do j = 0, mp_id(i, 1) - 1, 1
                row = mp_ii(mp_id(i, 2) + j, 1)
                col = mp_ii(mp_id(i, 2) + j, 2)
                L   = landtypes(row, col)
                if (L /= maxlc) temps(i) = temps(i) + (var2d(row, col)-tempm(i))**2
            end do
            ! 求标准差
            temps(i) = sqrt(temps(i) / p_num(i))
            if (temps(i) > th_onelayer(2*ii)) ref_sjx(i) = 1
        end do
        !$OMP END PARALLEL DO

        print*, trim(RL_onelayer(ii))//"_s", minval(temps, mask=(temps /= 0.)), maxval(temps, mask=(temps /= 0.))
        num_onelayer = num_onelayer + 1
        allocate(output2d(num_onelayer)%varms(sjx_points))
        output2d(num_onelayer)%varms = temps
        deallocate(temps)
        num_swithes  = num_swithes + 1
        ref_th(:, num_swithes) = ref_sjx

    END SUBROUTINE mean_std_cal2d

    SUBROUTINE mean_std_cal3d(ii, sjx_points, flag_temps, mp_id, mp_ii, var3d, num_twolayer, p_num)
        integer, intent(in) :: ii, sjx_points, flag_temps
        integer,  dimension(:, :), intent(in) :: mp_id, mp_ii
        real(r8), dimension(:, :, :), intent(in) :: var3d
        real(r8), dimension(:, :), allocatable :: tempm, temps! 均值计算/标准差计算的中间变量
        integer :: i, j, L, row, col
        integer,  intent(inout) :: num_twolayer
        integer,  dimension(:), intent(inout) :: p_num

        allocate(tempm(sjx_points, 2)); tempm = 0.
        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i, j, row, col, L)
        do i = 2, sjx_points, 1
            if (IsInRfArea_sjx(i) .eqv. .false.) cycle
            do j = 0, mp_id(i, 1)-1, 1
                row = mp_ii(mp_id(i, 2)+j, 1) ! 获取非结构网格包含经纬度网格的经度序号
                col = mp_ii(mp_id(i, 2)+j, 2) ! 获取非结构网格包含经纬度网格的纬度序号
                if (L /= maxlc) then
                    tempm(i, :) = tempm(i, :) + var3d(:, row, col)
                    p_num(i) = p_num(i) + 1
                end if
            end do
            tempm(i, :) = tempm(i, :) / p_num(i) !求均值
        end do
        !$OMP END PARALLEL DO

        if (refine_twolayer(2*ii-1)) then 
            num_twolayer = num_twolayer + 1
            allocate(output3d(num_twolayer)%varms(sjx_points, 2))
            output3d(num_twolayer)%varms = tempm
            print*, trim(RL_twolayer(ii))//"_m", minval(tempm, mask=(tempm /= 0.)), maxval(tempm, mask=(tempm /= 0.))
            do i = 2, sjx_points, 1
                if ((tempm(i, 1) > th_twolayer(2*ii-1, 1)) .or. (tempm(i, 2) > th_twolayer(2*ii, 2))) ref_sjx(i) = 1
            end do
            num_swithes = num_swithes + 1
            ref_th(:, num_swithes) = ref_sjx
        end if

        if (flag_temps == 0) then
            deallocate(tempm)
            return
        end if

        allocate(temps(sjx_points, 2)); temps = 0.
        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i, j, row, col, L)
        do i = 2, sjx_points, 1
            if (IsInRfArea_sjx(i) .eqv. .false.) cycle
            do j = 0, mp_id(i, 1) - 1, 1
                row = mp_ii(mp_id(i, 2) + j, 1)
                col = mp_ii(mp_id(i, 2) + j, 2)
                if (L /= maxlc) temps(i, :) = temps(i, :) + (var3d(:, row, col)-tempm(i, :))**2
            end do
            ! 求标准差
            temps(i, :) = sqrt(temps(i, :) / p_num(i))
            if ((temps(i, 1) > th_twolayer(2*ii-1, 1)) .or. (temps(i, 2) > th_twolayer(2*ii, 2))) ref_sjx(i) = 1
        end do
        !$OMP END PARALLEL DO

        num_twolayer = num_twolayer + 1
        allocate(output3d(num_twolayer)%varms(sjx_points, 2))
        output3d(num_twolayer)%varms = temps
        deallocate(temps)
        num_swithes = num_swithes + 1
        ref_th(:, num_swithes) = ref_sjx

    END SUBROUTINE mean_std_cal3d

END MODULE MOD_GetThreshold



