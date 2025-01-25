module MOD_data_preprocess

   USE consts_coms, only : io6, r8, source_dir, nlons_source, nlats_source, lcs, maxlc
   USE netcdf
   implicit none
   integer,  allocatable, public :: landtypes(:, :)
   real(r8), allocatable, public :: lon_i(:), lat_i(:)
   real(r8), allocatable, public :: lon_vertex(:), lat_vertex(:)  ! 经纬度顶点信息
   contains

   SUBROUTINE data_preprocess()

      IMPLICIT NONE
      real(r8) :: dx, dy
      integer :: ncid, varid, dimID_lon, dimID_lat, ierr
      character(LEN = 10) :: lon_name, lat_name
      character(LEN = 256) :: lndname
      real(r8), allocatable :: landtypes_read(:, :)

      ! landtypes数据读入
      print*, "landtypes read start"
      source_dir = trim(source_dir)//trim(lcs) // '/'  
      if (lcs == "igbp") then 
         lndname = trim(source_dir) // 'landtype_igbp_update.nc'
         CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
         CALL CHECK(NF90_INQ_DIMID(ncid, "lon", dimID_lon))
         CALL CHECK(NF90_INQ_DIMID(ncid, "lat", dimID_lat))
      else
         lndname = trim(source_dir) // 'landtype_usgs_update.nc'
         CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid))
         CALL CHECK(NF90_INQ_DIMID(ncid, "longitude", dimID_lon))
         CALL CHECK(NF90_INQ_DIMID(ncid, "latitude", dimID_lat))
      end if 
      CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lon, len = nlons_source))
      CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lat, len = nlats_source))
      print*, lndname
      print*, "nlons_source = ", nlons_source
      print*, "nlats_source = ", nlats_source
      allocate(landtypes_read(nlons_source, nlats_source)); landtypes_read = 0.
      allocate(landtypes(nlons_source, nlats_source)); landtypes = 0
      CALL CHECK(NF90_INQ_VARID(ncid, "landtype", varid))
      CALL CHECK(NF90_GET_VAR(ncid, varid, landtypes_read))
      CALL CHECK(NF90_CLOSE(ncid))
      landtypes = int(landtypes_read)
      deallocate(landtypes_read)
      print*, "landtypes", minval(landtypes), maxval(landtypes)
      maxlc = maxval(landtypes)
      print*, "landtypes read finish"

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


   END SUBROUTINE data_preprocess

END Module MOD_data_preprocess
