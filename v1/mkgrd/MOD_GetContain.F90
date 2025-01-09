!制作三角形网格包含数据

module MOD_GetContain

    use consts_coms
    use netcdf

    implicit none

contains
    subroutine  Get_Contain()

        real(r8), allocatable :: lon_i(:), lat_i(:)            ! 经纬度网格中心点
        real(r8), allocatable :: area_fine_gridcell(:, :)      ! 经纬度网格各单元面积
        real(r8), allocatable :: mp(:, :), wp(:, :)            ! 记录三角形网格与多边形网格中心点与面积
        real(r8), allocatable :: mp_ii(:, :)                   ! 记录三角形网格包含的经纬度网格(处理前)
        real(r8), allocatable :: mp_ii_new(:, :)               ! 记录三角形网格包含的经纬度网格(处理后)
        real(r8) :: maxlat, minlat                             ! 输入非结构网格的纬度极值       
        real(r8) :: isinply                                    ! 计算得出的包含比例
        real(r8) :: sjx(7, 2), dx, dy
        integer, allocatable :: mp_id(:, :)                    ! 记录三角形网格包含多边形网格数量与在mp_ii中的起始位置(处理前)
        integer, allocatable :: mp_id_new(:, :)                ! 记录三角形网格包含多边形网格数量与在mp_ii中的起始位置(处理后)
        integer, allocatable :: ngrmw(:, :), ngrwm(:, :)       ! 三角形网格与多边形网格中心点的邻域数组
        real(r8), allocatable :: landtypes(:, :)                ! 土地类型
        integer, allocatable :: seaorland(:, :)                ! 判断经纬度网格为陆地或海洋
        integer :: sjx_points, lbx_points                      ! 三角形网格数与多边形网格数
        integer :: ncid, varid(6), dimID_sjx, dimID_lbx        ! nc文件存储相关变量
        integer :: i, j, k, l, sum_land, sum_sea
        integer :: num_i, id, numpatch
        integer :: num_null                                    ! 需要处理的空值数量
        character(LEN = 256) :: lndname, nxpc,flnm
        character(LEN = 20),dimension(6) :: p_name
        logical :: ispart                                      ! 表示该经纬度网格是否被完全包含
        logical,allocatable :: IsInRfArea(:,:)                 ! 表示该经纬度网格是否在细化区域中
        
        isinply = 0.
        p_name = [character(len=20) :: "GLONW", "GLATW", "GLONM", "GLATM", "itab_w%im", "itab_m%iw"]

        print*, "开始读取非结构网格数据......"
        print*, ""
        write(nxpc, '(I4.4)')NXP
        flnm = trim(base_dir) // trim(EXPNME) // '/makegrid/gridfile/gridfile_NXP' // trim(nxpc) // '.nc4'
        print*, flnm
        CALL CHECK(NF90_OPEN(trim(flnm), nf90_nowrite, ncid))
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

        allocate(wp(lbx_points, 3))          ! 经度，纬度，面积
        allocate(mp(sjx_points, 3))
        allocate(mp_id(sjx_points, 2))    ! 非结构网格包含经纬度网格数量，在mp_ii位置
        allocate(ngrwm(8, lbx_points))       ! 邻域数组，多出来的一维表示相邻点数
        allocate(ngrmw(4, sjx_points))
        wp = 0.
        mp = 0.
        mp_id = 0
        ngrmw = 0
        ngrwm = 0

        CALL CHECK(NF90_GET_VAR(ncid, varid(1), wp(:, 1)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(2), wp(:, 2)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(3), mp(:, 1)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(4), mp(:, 2)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(5), ngrwm(1:7, :)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(6), ngrmw(1:3, :)))
        CALL CHECK(NF90_CLOSE(ncid))

        print*, "非结构网格数据读取完成"
        print*, ""
        print*, "共有", sjx_points, "个三角形网格和", lbx_points, "个多边形网格"
        print*, ""

        ! 获取每个m(w)点的相邻w(m)点数量
        call GetNgrNum(sjx_points, lbx_points, ngrmw, ngrwm)

        ! 地表覆盖类型数组
        allocate(landtypes(nlons_source, nlats_source))

        if(lcs == "igbp")then
           lndname = trim(source_dir) // 'landtype_igbp_update.nc'
        else
           lndname = trim(source_dir) // 'landtype_usgs_update.nc'
        end if

        print*, lndname
        CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
        CALL CHECK(NF90_INQ_VARID(ncid, "landtype", varid(1)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(1), landtypes))
        CALL CHECK(NF90_CLOSE(ncid))
        print*, "landtypes", minval(landtypes), maxval(landtypes)

        ! 地表网格分辨率
        dx = 360. / nlons_source
        dy = 180. / nlats_source

        allocate(lon_i(nlons_source))
        allocate(lat_i(nlats_source))

        do i = 1, nlons_source, 1
            lon_i(i) = -180. + (2 * i - 1) * dx / 2.
        end do

        do i = 1, nlats_source, 1
            lat_i(i) = 90. - (2 * i - 1) * dy / 2.
        end do

        allocate(seaorland(nlons_source, nlats_source))
        seaorland = 0
        sum_sea = 0
        sum_land = 0
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
        print*, "海洋网格个数为", sum_sea, "，陆地网格个数为", sum_land
        deallocate(landtypes)

        allocate(area_fine_gridcell(nlons_source, nlats_source))
        area_fine_gridcell(:, :) = 0.
        print*, "开始计算经纬度网格面积......"
        call cellarea(area_fine_gridcell, nlons_source, nlats_source)
        print*, "经纬度网格面积计算完成"

        maxlat = maxval(mp(:, 2))
        minlat = minval(mp(:, 2))

        num_i = 0

        ! 判断经纬度网格是否位于细化区域内
        allocate(IsInRfArea(nlons_source,nlats_source))
        IsInRfArea = .false.
        CALL IsInRefineArea(IsInRfArea,lon_i,lat_i,dx,dy)

        ! --------------------------------------------------------
        ! 获取三角形网格中经纬度网格的数量及所占面积比例
        ! 1.计算数组大小
        ! 2.分配内存
        ! 3.计算包含关系
        ! --------------------------------------------------------


        ! 首先计算数组大小
        print*, "开始计算三角形网格与经纬度网格包含关系数组大小......"
        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
        !$OMP PRIVATE(i,j,k,l,sjx,isinply,ispart)
        do i = 1, nlons_source, 1
            num_i = num_i + 1
            !print*, num_i
            do j = 1, nlats_source, 1                    ! 循环遍历初始网格单元
               if(IsInRfArea(i,j) .eqv. .false.)then
                  cycle
               end if
               ispart = .false.  ! 表示是否被完全包含
               if(seaorland(i, j) == 0)then
                  cycle
               end if                           ! 只遍历海洋网格
               do k = 1, sjx_points, 1           ! 循环遍历三角形网格
                  if(ispart .eqv. .true.)then
                     exit
                  end if
                  if(ngrmw(4, k) == 3)then       ! 判断m点邻域是否有三个被记录的w点
                     sjx = 0.
                     do l = 1, 3, 1
                        sjx(l, 1:2) = wp(ngrmw(l, k), 1:2)  ! 记录三角形网格顶点信息
                     end do

                     isinply = IsInUstrGrid(sjx, lon_i(i), lat_i(j), 3, maxlat, minlat)

                     if(isinply == 1.)then      ! 表示完全包含
                        ispart = .true.
                     end if

                     if(isinply > 0.)then       ! 记录每个非结构网格包含经纬度网格数量
                        mp_id(k, 1) = mp_id(k, 1) + 1
                     end if
                  end if
               end do
            end do
         end do
         !$OMP END PARALLEL DO

        print*, "包含关系数组大小计算完成......"

        ! 分配数组内存
        ! 为了防止在第二次重复计算时与第一次计算产生偏差
        ! 给每个非结构网格多预留一倍的空间
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
        num_i = 0

      print*, "开始计算三角形网格与经纬度网格包含关系......"
      !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
      !$OMP PRIVATE(i,j,k,l,ispart,sjx,isinply,id)
      do i = 1, nlons_source, 1
            num_i = num_i + 1
            !print*, num_i
            do j = 1, nlats_source, 1
               if(IsInRfArea(i,j) .eqv. .false.)then
                  cycle
               end if
               ispart = .false.
               if(seaorland(i, j) == 0)then
                  cycle
               end if
               do k = 1, sjx_points, 1
                  if(ispart .eqv. .true.)then
                     exit
                  end if
                  if(ngrmw(4, k) == 3)then
                     sjx = 0.
                     do l = 1, 3, 1
                        sjx(l, 1:2) = wp(ngrmw(l, k), 1:2)
                     end do

                     isinply = IsInUstrGrid(sjx, lon_i(i), lat_i(j), 3, maxlat, minlat)

                     if(isinply == 1.)then
                        mp_id(k, 1) = mp_id(k, 1) + 1
                        mp(k, 3) = mp(k, 3) + area_fine_gridcell(i, j)     ! 非结构网格面积

                        id = mp_id(k, 2) + mp_id(k, 1) - 1
                        mp_ii(id, 1) = i  ! 经度索引
                        mp_ii(id, 2) = j  ! 纬度索引
                        mp_ii(id, 3) = isinply  ! 包含比例
                        mp_ii(id, 4) = area_fine_gridcell(i, j)   ! 包含面积

                        ispart = .true.
                     else if(isinply > 0.)then
                        mp_id(k, 1) = mp_id(k, 1) + 1
                        mp(k, 3) = mp(k, 3) + area_fine_gridcell(i, j) * isinply

                        id = mp_id(k, 2) + mp_id(k, 1) - 1
                        mp_ii(id, 1) = i
                        mp_ii(id, 2) = j
                        mp_ii(id, 3) = isinply
                        mp_ii(id, 4) = area_fine_gridcell(i, j) * isinply

                     end if
                  end if
               end do
            end do
      end do
      !$OMP END PARALLEL DO

        ! 将由于多预留空间导致的空值删去
        ! 并重新计算序号
        do i = 1, sjx_points, 1
            num_null = 0
            do j = mp_id(i, 2), mp_id(i, 2) + mp_id(i, 1) - 1, 1
                if(mp_ii(j, 3) == 0.)then
                    num_null = num_null + 1
                end if
            end do
            mp_id(i, 1) = mp_id(i, 1) - num_null
        end do

        ! 分配储存消除空值后信息的数组内存
        allocate(mp_id_new(sjx_points, 2))
        mp_id_new = 0
        mp_id_new(:, 1) = mp_id(:, 1)
        mp_id_new(1, 2) = 1

        do i = 2, sjx_points, 1
            mp_id_new(i, 2) = mp_id_new(i - 1, 2) + mp_id_new(i - 1, 1)
        end do

        numpatch = INT(sum(mp_id_new(:, 1)))
        allocate(mp_ii_new(numpatch, 4))
        mp_ii_new = 0

        do i = 1, sjx_points, 1
            do j = 1, 4, 1
                if(mp_id_new(i, 1) /= 0)then
                    mp_ii_new(mp_id_new(i, 2):mp_id_new(i, 2) + mp_id_new(i, 1) - 1, j) = &
                            mp_ii(mp_id(i, 2):mp_id(i, 2) + mp_id_new(i, 1) - 1, j)
                end if
            end do
        end do

        print*, "三角形网格包含数组计算完成，开始进行存储......"
        print*, ""
        flnm = trim(base_dir) // trim(EXPNME) // '/makegrid/contain/initial' 
        call SaveFile(mp, mp_id_new, mp_ii_new, sjx_points, numpatch, flnm)
        print*, "非结构网格存储完成"
        print*, ""


    end SUBROUTINE Get_Contain

    SUBROUTINE CHECK(STATUS)
        INTEGER, intent (in) :: STATUS
        if  (STATUS .NE. NF90_NOERR) then ! nf_noerr=0 表示没有错误
            print *, NF90_STRERROR(STATUS)
            stop 'stopped'
        endif
    END SUBROUTINE CHECK


    ! 网格邻域计算处理，计算m(w)点相邻的w(m)点数量
    SUBROUTINE GetNgrNum(sjx_points, lbx_points, ngrmw, ngrwm)

        implicit None

        integer :: i, j, flag, sjx_points, lbx_points
        integer, dimension(4, sjx_points) :: ngrmw
        integer, dimension(8, lbx_points) :: ngrwm

        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i,j,flag)
        do i = 1, sjx_points, 1
            flag = 0
            do j = 1, 3, 1
                if(ngrmw(j, i) /= 1)then
                    flag = flag + 1
                end if
            end do  ! 1号m点和w点均为零值
            ngrmw(4, i) = flag
        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i,j,flag)
        do i = 1, lbx_points, 1
            flag = 0
            do j = 1, 7, 1
                if(ngrwm(j, i) /= 1)then
                    flag = flag + 1
                end if
            end do
            ngrwm(8, i) = flag
        end do
        !$OMP END PARALLEL DO

    END SUBROUTINE GetNgrNum


    ! 计算经纬度网格面积
    SUBROUTINE CellArea(area, nlons, nlats)

        implicit none

        integer :: i, j
        integer, intent(in) :: nlons, nlats
        real(r8), intent(out) :: area(nlons, nlats)
        real(r8) :: re, pi, deg2rad, global, dx, dy, error
        real(r8) :: lats(nlats), latn(nlats)
        real(r8) :: lonw(nlons), lone(nlons)

        re = 6.37122e6 * 0.001                    ! kilometer
        pi = 4. * atan(1.)
        deg2rad = pi / 180.
        global = 0.

        dx = 360. / nlons
        dy = 180. / nlats

        do i = 1, nlons, 1
            lone(i) = -180. + i * dx
            lonw(i) = -180. + (i - 1) * dx
        end do

        do i = 1, nlats, 1
            latn(i) = 90. - (i - 1) * dy
            lats(i) = 90. - i * dy
        end do

        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1)&
        !$OMP PRIVATE(i,j,dx,dy)
        do j = 1, nlats, 1
            do i = 1, nlons, 1
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


    SUBROUTINE Savefile(ustr, ustr_id, ustr_ii, num, num_ii, flnm)

        use NETCDF

        IMPLICIT NONE

        integer :: ncID, idDimID, pidDimID, infoDimID, ncVarID, iunit
        integer, intent(in) :: num, num_ii
        character(len = 256) :: outputfile(6), flnm
        real(r8), dimension(num, 3), intent(in) :: ustr
        real(r8), dimension(num_ii, 4), intent(in) :: ustr_ii
        integer, dimension(num, 2), intent(in) :: ustr_id

        outputfile(1) = trim(flnm)//"/mp.nc4"
        outputfile(2) = trim(flnm)//"/mp_ii.nc4"
        outputfile(3) = trim(flnm)//"/mp_id.nc4"
        outputfile(4) = trim(flnm)//"/mp.bin"
        outputfile(5) = trim(flnm)//"/mp_ii.bin"
        outputfile(6) = trim(flnm)//"/mp_id.bin"

        print*, outputfile(1)
        CALL CHECK(NF90_CREATE(trim(outputfile(1)), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", num, idDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dima", 3, infoDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "info", NF90_FLOAT, (/ idDimID, infoDimID /), ncVarID))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID, ustr))
        CALL CHECK(NF90_CLOSE(ncID))

        print*, outputfile(2)
        CALL CHECK(NF90_CREATE(trim(outputfile(2)), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", num_ii, idDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dima", 4, infoDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "info", NF90_FLOAT, (/ idDimID, infoDimID /), ncVarID))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID, ustr_ii))
        CALL CHECK(NF90_CLOSE(ncID))

        print*, outputfile(3)
        CALL CHECK(NF90_CREATE(trim(outputfile(3)), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", num, idDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dima", 2, infoDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "info", NF90_INT, (/ idDimID, infoDimID /), ncVarID))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID, ustr_id))
        CALL CHECK(NF90_CLOSE(ncID))

        iunit = 100
        print*, outputfile(4)
        open(iunit, file = trim(outputfile(4)), form = 'unformatted', status = 'unknown')
        write(iunit) ustr
        close(iunit)

        iunit = 101
        print*, outputfile(5)
        open(iunit, file = trim(outputfile(5)), form = 'unformatted', status = 'unknown')
        write(iunit) ustr_ii
        close(iunit)

        iunit = 102
        print*, outputfile(6)
        open(iunit, file = trim(outputfile(6)), form = 'unformatted', status = 'unknown')
        write(iunit) ustr_id
        close(iunit)

    END SUBROUTINE SaveFile


    ! 计算非结构网格与经纬度网格的包含关系
    REAL FUNCTION IsInUstrGrid(ustr, lon, lat, num_points, maxlat_m, minlat_m)

        implicit none

        integer :: inc(4), i, j, iscross_l(2), ispole
        integer, intent(in) :: num_points
        integer :: num_points_i, sjxorlbx_i, num_inter
        real(r8), intent(in) :: lon, lat, maxlat_m, minlat_m
        real(r8), dimension(7, 2), intent(in) :: ustr  ! 非结构网格单元顶点
        real(r8), allocatable :: ustr_move(:, :)             ! 非结构网格经度移动后的点
        real(r8), dimension(4, 2) :: point                   ! 经纬度网格顶点
        real(r8), dimension(2) :: center_point              ! 非结构网格单元中心点
        integer, allocatable :: iscross_g(:)                ! 判断非结构网格线段是否穿过经纬度网格
        real(r8), dimension(20, 2) :: interarea_points       ! 两网格相交区域顶点
        real(r8), allocatable :: inter_points(:, :, :)        ! 两种网格的交点
        real(r8), allocatable :: area(:)   ! 由非结构网格相邻两点与任一经纬度网格顶点组成的三角形面积
        real(r8), dimension(5) :: area_i   ! 由经纬度网格相邻两点与任一非结构网格顶点组成的三角形面积
        real(r8) :: dx, dy, maxlat, minlat, maxlon, minlon
        real(r8) :: tmpa(2),tmpb(2),tmpc(2),tmpd(3,2)

        inc = 0                 ! 判断经纬度网格顶点在非结构网格中的数量
        IsInUstrGrid = 0        ! 判断经纬度网格与非结构网格位置关系
        center_point = 0.       ! 经纬度网格中心点
        ispole = 0              ! 判断非结构网格是否包含极点
        interarea_points = 0.   ! 两种网格重叠多边形顶点
        iscross_l = 0
        point = 0.
        num_inter = 0
        area_i = 0.

        maxlat = 0.
        maxlon = 0.
        minlat = 0.
        minlon = 0.

        if(mode == 3)then
           sjxorlbx_i = 3
        else
           sjxorlbx_i = 7
        end if
      
        num_points_i = num_points

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

        !print*,point
        maxlat = maxval(ustr(1:num_points_i, 2))
        minlat = minval(ustr(1:num_points_i, 2))

        ! 确保经纬度网格纬度绝对值不超过90°
        do i = 1, 4, 1
            if(point(i, 2) > 90.)then
                point(i, 2) = 90.
            else if(point(i, 2) < -90.)then
                point(i, 2) = -90.
            end if
        end do

        ! 根据纬度筛选网格
        if((point(3, 2) > maxlat).or.(point(1, 2) < minlat))then
            IsInUstrGrid = -1
            return
        end if

        !------------------------------------------------------------------
        ! 处理三角形网格最顶部与最底部的网格（极点）
        !------------------------------------------------------------------
        if(sjxorlbx_i == 3)then
            if(ustr(1, 2) == 90.)then
                ispole = 1
                sjxorlbx_i = 7
                num_points_i = 4
                allocate(ustr_move(num_points_i, 2))
                ustr_move(3:4, :) = ustr(2:3, :)
                ustr_move(2, 1) = ustr(2, 1)
                ustr_move(1, 1) = ustr(3, 1)
                ustr_move(1:2, 2) = 90.
            else if(ustr(1, 2) == -90.)then
                ispole = -1
                sjxorlbx_i = 7
                num_points_i = 4
                allocate(ustr_move(num_points_i, 2))
                ustr_move(3:4, :) = ustr(2:3, :)
                ustr_move(2, 1) = ustr(2, 1)
                ustr_move(1, 1) = ustr(3, 1)
                ustr_move(1:2, 2) = -90.
            else
                allocate(ustr_move(num_points_i, 2))
                ustr_move = ustr
            end if
        else
            if(minval(ustr(1:num_points_i, 2)) == minlat_m)then
                ispole = 2
            else if(maxval(ustr(1:num_points_i, 2)) == maxlat_m)then
                ispole = -2
            else
                allocate(ustr_move(num_points_i, 2))
                ustr_move(:, :) = ustr(1:num_points_i, :)
            end if
        end if

        !------------------------------------------------------------------
        ! 处理多边形网格最顶部与最底部的网格（极点）
        !------------------------------------------------------------------
        if(ispole == 2)then
            if(point(1, 2) > maxlat_m)then
                if(point(3, 2) > maxlat_m)then
                    IsInUstrGrid = 1
                    return
                else
                    IsInUstrGrid = (point(1, 2) - maxlat_m) / dy
                    return
                end if
            end if
            return
        else if(ispole == -2)then
            if(point(3, 2) < minlat_m)then
                if(point(1, 2) < minlat_m)then
                    IsInUstrGrid = 1
                    return
                else
                    IsInUstrGrid = (minlat_m - point(3, 2)) / dy
                    return
                end if
            end if
            return
        end if

        allocate(inter_points(num_points_i, 3, 2))
        allocate(area(num_points_i + 1))
        allocate(iscross_g(num_points_i))
        inter_points = 0
        iscross_g = 0
        area = 0

        !-------------------------------------------------------------------------------------
        ! 判断两网格是否越过±180°经线
        !-------------------------------------------------------------------------------------
        if((point(1, 1) > 180.).or.(point(2, 1) < -180.))then
            iscross_l(1) = 1
        end if

        iscross_l(2) = IsCrossLine(ustr_move(:, 1), num_points_i)

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
            call MoveLons(ustr_move(:, 1), num_points_i)
        end if

        do i = 1, num_points_i, 1
            center_point = center_point + ustr_move(i, :)
        end do
        center_point = center_point / num_points_i

        !-------------------------------------------------------------------------------------
        ! 根据经度筛选网格
        !-------------------------------------------------------------------------------------
        minlon = minval(ustr_move(1:num_points_i, 1))
        maxlon = maxval(ustr_move(1:num_points_i, 1))

        if((point(2, 1) > maxlon).or.(point(1, 1) < minlon))then
            IsInUstrGrid = -1
            return
        end if

        !-------------------------------------------------------------------------------------
        ! 开始判断两个网格间的位置关系，并记录关键点
        !-------------------------------------------------------------------------------------
        ! 首先判断两个网格相交
        ! 判断非结构网格的边与经纬度网格相交关系
        do i = 1, num_points_i - 1, 1
            tmpa = ustr_move(i, 1:2)
            tmpb = ustr_move(i + 1, 1:2)
            tmpd = inter_points(i, 1:3, 1:2)

            iscross_g(i) = IsCrossGrid(point, tmpa, tmpb, tmpd)

            !iscross_g(i) = IsCrossGrid(point, ustr_move(i, :), ustr_move(i + 1, :), &
            !        inter_points(i, :, :))
        end do

        tmpa = ustr_move(1, 1:2)
        tmpb = ustr_move(num_points_i, 1:2)
        tmpd = inter_points(num_points_i, 1:3, 1:2)

        iscross_g(num_points_i) = IsCrossGrid(point, tmpa, tmpb, tmpd)

        !iscross_g(num_points_i) = IsCrossGrid(point, ustr_move(1, :), ustr_move(num_points_i, :), &
        !        inter_points(num_points_i, :, :))

        ! 计算经纬度网格位于非结构网格中的顶点数目
        ! 若为4，则包含，若为1~3，则相交
        do i = 1, 4, 1
            area = 0.
            if(sjxorlbx_i == 3)then
                tmpa = ustr_move(1, 1:2)
                tmpb = ustr_move(2, 1:2)
                tmpc = point(i, 1:2)
                area(1) = GetTriangleArea(tmpa, tmpb, tmpc)

                tmpa = ustr_move(2, 1:2)
                tmpb = ustr_move(3, 1:2)
                tmpc = point(i, 1:2)
                area(2) = GetTriangleArea(tmpa, tmpb, tmpc)

                tmpa = ustr_move(1, 1:2)
                tmpb = ustr_move(3, 1:2)
                tmpc = point(i, 1:2)
                area(3) = GetTriangleArea(tmpa, tmpb, tmpc)

                tmpa = ustr_move(1, 1:2)
                tmpb = ustr_move(2, 1:2)
                tmpc = ustr_move(3, 1:2)
                area(4) = GetTriangleArea(tmpa, tmpb, tmpc)

                !area(1) = GetTriangleArea(ustr_move(1, :), ustr_move(2, :), point(i, :))
                !area(2) = GetTriangleArea(ustr_move(2, :), ustr_move(3, :), point(i, :))
                !area(3) = GetTriangleArea(ustr_move(1, :), ustr_move(3, :), point(i, :))
                !area(4) = GetTriangleArea(ustr_move(1, :), ustr_move(2, :), ustr_move(3, :))

                ! 若该点与三角形网格三个顶点连线构成的三个三角形面积和与三角形网格面积相同
                ! 则该点位于三角形网格内
                if(abs(area(1) + area(2) + area(3) - area(4)) < 0.00005)then
                    inc(i) = 1
                end if
            else if(sjxorlbx_i == 7)then
                do j = 1, num_points_i - 1, 1
                    tmpa = ustr_move(j, 1:2)
                    tmpb = ustr_move(j + 1, 1:2)
                    tmpc = point(i, 1:2)
                    area(j) = GetTriangleArea(tmpa, tmpb, tmpc)

                    tmpa = ustr_move(j, 1:2)
                    tmpb = ustr_move(j + 1, 1:2)
                    area(num_points_i + 1) = area(num_points_i + 1) + GetTriangleArea(tmpa, tmpb, center_point)
                    !area(j) = GetTriangleArea(ustr_move(j, :), ustr_move(j + 1, :), point(i, :))
                    !area(num_points_i + 1) = area(num_points_i + 1) + GetTriangleArea(ustr_move(j, :), &
                    !        ustr_move(j + 1, :), center_point)
                end do

                tmpa = ustr_move(1, 1:2)
                tmpb = ustr_move(num_points_i, 1:2)
                tmpc = point(i, 1:2)
                area(num_points_i) = GetTriangleArea(tmpa, tmpb, tmpc)

                tmpa = ustr_move(1, 1:2)
                tmpb = ustr_move(num_points_i, 1:2)
                area(num_points_i + 1) = area(num_points_i + 1) + GetTriangleArea(tmpa, tmpb, center_point)

                !area(num_points_i) = GetTriangleArea(ustr_move(1, :), ustr_move(num_points_i, :), point(i, :))
                !area(num_points_i + 1) = area(num_points_i + 1) + GetTriangleArea(ustr_move(1, :), &
                !        ustr_move(num_points_i, :), center_point)

                if(abs(sum(area(1:num_points_i)) - area(num_points_i + 1)) < 0.00005)then
                    inc(i) = 1
                end if
            end if
        end do

        ! 当不计算比例时，只要相交就视为完全包含
        if (no_caculate_fraction) then
            if((sum(inc) > 0).or.(sum(iscross_g) > 0))then
               IsInUstrGrid = 1.
               return
            else
               IsInUstrGrid = 0
               return
            end if
        endif

        if(sum(inc) == 4)then   ! 若包含
            IsInUstrGrid = 1.
            return
        else if((sum(inc) == 0).and.(sum(iscross_g) == 0))then  ! 若不包含也不相交
            IsInUstrGrid = 0.
            return
        else
            ! 计算并记录非结构网格顶点是否位于经纬度网格内
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
            center_point = center_point / 4.

            do i = 1, num_points_i, 1
                area_i = 0.
                do j = 1, 3, 1
                    tmpa = ustr_move(i, 1:2)
                    tmpb = point(j, 1:2)
                    tmpc = point(j + 1, 1:2)
                    area_i(j) = GetTriangleArea(tmpa, tmpb, tmpc)

                    tmpb = point(j, 1:2)
                    tmpc = point(j + 1, 1:2)
                    area_i(5) = area_i(5) + GetTriangleArea(center_point, tmpb, tmpc)

                    !area_i(j) = GetTriangleArea(ustr_move(i, :), point(j, :), point(j + 1, :))
                    !area_i(5) = area_i(5) + GetTriangleArea(center_point, point(j, :), point(j + 1, :))
                end do

                tmpa = ustr_move(i, 1:2)
                tmpb = point(1, 1:2)
                tmpc = point(4, 1:2)
                area_i(4) = GetTriangleArea(tmpa, tmpb, tmpc)

                tmpb = point(1, 1:2)
                tmpc = point(4, 1:2)
                area_i(5) = area_i(5) + GetTriangleArea(center_point, tmpb, tmpc)

                !area_i(4) = GetTriangleArea(ustr_move(i, :), point(1, :), point(4, :))
                !area_i(5) = area_i(5) + GetTriangleArea(center_point, point(1, :), point(4, :))
                if(abs(sum(area_i(1:4)) - area_i(5)) < 0.0006)then
                    num_inter = num_inter + 1
                    interarea_points(num_inter, :) = ustr_move(i, :)
                end if
            end do

            do i = 1, num_points_i, 1
                if(inter_points(i, 3, 1) /= 0)then
                    num_inter = num_inter + 1
                    interarea_points(num_inter:num_inter + inter_points(i, 3, 1) - 1, :) = &
                            inter_points(i, 1:inter_points(i, 3, 1), :)
                    num_inter = num_inter + inter_points(i, 3, 1) - 1
                end if
            end do

            ! num_iter表示两种网格顶点位于对方内部与交点的总个数
            if(num_inter == 0)then
               IsInUstrGrid = 0.
               return
            end if

            ! 将两种网格顶点位于对方内部与交点进行排序
            call SortPoints(interarea_points, num_inter)

            ! 计算相交图形面积占经纬度网格面积比例
            IsInUstrGrid = GetAreaPercent(interarea_points, num_inter, point)

        end if

    END FUNCTION IsInUstrGrid


    ! 判断网格是否越过189°与-180°经线
    integer function IsCrossLine(lons, num)

        implicit none

        integer, intent(in) :: num
        integer :: i, j
        real(r8), dimension(num), intent(in) :: lons

        IsCrossLine = 0

        do i = 1, num - 1, 1
            do j = i + 1, num, 1
                if(abs(lons(j) - lons(i)) > 180.)then
                    IsCrossLine = 1
                    return
                end if
            end do
        end do

    end function IsCrossLine


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
    integer function IsCrossGrid(point, a, b, inter_point)

        implicit none

        real(r8), dimension(2), intent(in) :: a, b      ! 直线端点
        real(r8), dimension(2) :: x, y
        real(r8), dimension(4, 2), intent(in) :: point
        real(r8), dimension(3, 2), intent(out) :: inter_point
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

    end function IsCrossGrid


    ! 计算三角形面积（按经纬度，并非实际面积）
    real function GetTriangleArea(a, b, c)

        implicit none

        real(r8), dimension(2) :: a, b, c
        real(r8) :: aa, bb, cc, p

        GetTriangleArea = 0

        aa = sqrt((c(1) - b(1)) * (c(1) - b(1)) + (c(2) - b(2)) * (c(2) - b(2)))
        bb = sqrt((c(1) - a(1)) * (c(1) - a(1)) + (c(2) - a(2)) * (c(2) - a(2)))
        cc = sqrt((a(1) - b(1)) * (a(1) - b(1)) + (a(2) - b(2)) * (a(2) - b(2)))

        p = (aa + bb + cc) / 2

        GetTriangleArea = sqrt(p * (p - aa) * (p - bb) * (p - cc))

    end function GetTriangleArea


    ! 获取非结构网格包含经纬度网格的比例
    REAL FUNCTION GetAreaPercent(inter_point, num, point)


        implicit none

        integer, intent(in) :: num
        real(r8), dimension(20, 2), intent(in) :: inter_point
        real(r8), dimension(4, 2), intent(in) :: point

        real(r8), dimension(2) :: center_point
        integer :: i
        real(r8) :: inter_area,tmpa(2),tmpb(2)

        GetAreaPercent = 0.
        inter_area = 0.
        center_point = 0.

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

        if(GetAreaPercent >= 1)then
            GetAreaPercent = 1.
        end if

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
            sort_i(i) = i
        end do
        center_point = center_point / num

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

        do i = 1, num, 1
            points(i, :) = points_i(sort_i(i), :) + center_point
        end do

    END SUBROUTINE SortPoints


    ! 判断经纬度网格是否位于细化区域
    SUBROUTINE IsInRefineArea(IsInRfArea,lon_i,lat_i,dx,dy)

      use refine_vars

      implicit none

      integer :: i,j,n
      real(r8),intent(in) :: dx,dy
      real(r8),dimension(nlons_source),intent(in) :: lon_i
      real(r8),dimension(nlats_source),intent(in) :: lat_i
      real(r8),dimension(nlons_source) :: lone,lonw
      real(r8),dimension(nlats_source) :: latn,lats
      logical,dimension(nlons_source,nlats_source),intent(out) :: IsInRfArea


      do i = 1,nlons_source,1
         lone(i) = lon_i(i) + dx / 2.
         lonw(i) = lon_i(i) - dx / 2.
      end do

      do j = 1,nlats_source,1
         latn(j) = lat_i(j) + dy / 2.
         lats(j) = lat_i(j) - dy / 2.
      end do

      IsInRfArea = .false.

      do n = 1,ndm_refine,1
         do i = 1,nlons_source,1
            do j = 1,nlats_source,1
               if((lone(i) < edgee_rf(n)).and.(lonw(i) > edgew_rf(n)).and.&
                        (latn(j) < edgen_rf(n)).and.(lats(j) > edges_rf(n)))then
                  IsInRfArea(i,j) = .true.
               end if
            end do
         end do
      end do

   END SUBROUTINE IsInRefineArea

END module MOD_GetContain
