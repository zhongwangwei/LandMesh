module MOD_data_preprocess

   USE consts_coms, only : io6, r8, pio180, erad, source_dir, nlons_source, nlats_source, lcs
   USE netcdf
   implicit none
   integer,  allocatable, public :: landtypes(:, :)
   real(r8), allocatable, public :: lon_i(:), lat_i(:)
   real(r8), allocatable, public :: lon_vertex(:), lat_vertex(:)  ! 经纬度顶点信息
   real(r8), allocatable, public :: area_fine_gridcell(:,:)       ! 经纬度网格各单元面积 
   contains

   SUBROUTINE data_preprocess()

      IMPLICIT NONE
      real(r8) :: dx, dy
      integer :: ncid, varid
      character(LEN = 256) :: lndname
      ! 地表网格分辨率
      dx = 360. / nlons_source
      dy = 180. / nlats_source
      ! 经纬度网格中心点经纬度值
      allocate(lon_i(nlons_source)); allocate(lat_i(nlats_source))
      lon_i = -180. + (2 * [1:nlons_source] - 1) * dx / 2. ! [] mean array 
      lat_i =   90. - (2 * [1:nlats_source] - 1) * dy / 2.

      allocate(lon_vertex(1+nlons_source)); allocate(lat_vertex(1+nlats_source))
      ! lon_vertex combined lone and lonw from -180 to 180
      lon_vertex(2:nlons_source+1) = lon_i + dx / 2.
      lon_vertex(1) = -180. ! lon_vertex(nlons_source+1) 与 lon_vertex(1) 会不会冲突呢？要小心了
      lon_vertex(nlons_source+1) = 180.
      ! lat_vertex combined latn and lats from 90 to -90
      lat_vertex(2:nlats_source+1) = lat_i - dy / 2.
      lat_vertex(1) = 90. ! 也可以考虑去掉
      lat_vertex(nlats_source+1) = -90.
       
      ! 经纬度网格面积计算
      allocate(area_fine_gridcell(nlons_source, nlats_source)); area_fine_gridcell(:, :) = 0.
      write(io6, *) 'cellarea calculate start'
      call cellarea(area_fine_gridcell) ! 有很大的改进空间
      write(io6, *) 'cellarea calculate complete'

      ! landtypes数据读入
      allocate(landtypes(nlons_source, nlats_source))
      if(lcs == "igbp")then
         lndname = trim(source_dir) // 'landtype_igbp_update.nc'
      else
         lndname = trim(source_dir) // 'landtype_usgs_update.nc'
      end if
      print*, lndname
      CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
      CALL CHECK(NF90_INQ_VARID(ncid, "landtype", varid))
      CALL CHECK(NF90_GET_VAR(ncid, varid, landtypes))
      CALL CHECK(NF90_CLOSE(ncid))
      print*, "landtypes", minval(landtypes), maxval(landtypes)

   END SUBROUTINE data_preprocess

    ! 计算经纬度网格面积
    SUBROUTINE CellArea(area)

        implicit none

        integer :: i, j
        real(r8), intent(inout) :: area(:, :)
        real(r8) :: dx, dy, global_area, error_area
        
        dx = (lon_vertex(2) - lon_vertex(1)) * pio180
        !$OMP PARALLEL DO NUM_THREADS(96) SCHEDULE(DYNAMIC,1)&
        !$OMP PRIVATE(j, dy)
        do j = 1, nlats_source, 1
            dy = sin(lat_vertex(j) * pio180) - sin(lat_vertex(j+1) * pio180)
            area(1, j) = dx * dy * erad * erad / 1e6
        end do
        !$OMP END PARALLEL DO
        do i = 1, nlons_source, 1
            area(i, :) = area(1, :)
        end do
        global_area = sum(area(:, :))
        print*, "global area = ",global_area
        ! stop "stop for area calculate test"
        ! 确保网格单元的总面积与其边缘定义的网格面积相同
        dx = (180. - (-180.)) * pio180
        dy = sin(90. * pio180) - sin(-90. * pio180)
        ! error = dx * dy * re * re
        error_area = dx * dy * erad * erad / 1e6
        if(abs(global_area - error_area) / error_area > 1.0e-7) then
            print*, 'CELLAREA error: correct area is ', error_area, &
                    ' but summed area of grid cells is ', global_area
        end if

        return

    END SUBROUTINE CellArea

END Module MOD_data_preprocess
