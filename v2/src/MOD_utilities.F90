!DESCRIPTION
!===========
! This module contains utility subroutines like CHECK for NetCDF operations
! and mode4mesh_make for mesh creation from mode files
!
!REVISION HISTORY
!----------------
! 2025.06.11  Zhongwang Wei @ SYSU (revised version)
! 2025.06.10  Rui Zhang @ SYSU (original version)

module MOD_utilities
    use netcdf ! For NF90_NOERR and NF90_STRERROR
    use consts_coms
    implicit none
    public :: CHECK, mode4mesh_make,Mode4_Mesh_Save
contains

SUBROUTINE CHECK(STATUS)
    INTEGER, intent (in) :: STATUS ! Input: Status code from a NetCDF operation
    if  (STATUS .NE. NF90_NOERR) then ! NF90_NOERR (usually 0) indicates no error
        print *, NF90_STRERROR(STATUS) ! Print the error message corresponding to the status code
        stop 'NetCDF operation failed. Stopping.' ! Stop execution
    endif
END SUBROUTINE CHECK

! Creates a mesh from a mode file (lonlat or lambert).
subroutine mode4mesh_make(inputfile, grid_select)
    use lonlatmesh_coms, only : mesh          
    implicit none
    character(pathlen), intent(in) :: inputfile 
    character(*), intent(in)       :: grid_select   
    character(pathlen)             :: lndname       
    logical                        :: fexists       
    integer                        :: i, j, idx, length
    integer                        :: nlon_cal, nlat_cal
    integer                        :: ncid, dimID_lon, dimID_lat
    integer                        :: lon_points, lat_points, bound_points, mode_points
    integer, dimension(10)         :: varid         
    real(r8)                       :: lon_start, lat_start, lon_end, lat_end, lon_grid_interval, lat_grid_interval
    real(r8), dimension(:),  allocatable :: lon_center, lat_center, lon_bound, lat_bound
    real(r8), dimension(:, :), allocatable :: lon_vert, lat_vert, lonlat_bound
    integer,  dimension(:, :), allocatable :: ngr_bound     
    integer,  dimension(:),    allocatable :: n_ngr         
    character(5)                   :: nxpc, stepc, numc  
    character(16)                  :: definition    

    write(io6, *) "file_dir : ", file_dir  
    lndname = inputfile             
    length = len_trim(lndname)      
    
    if (trim(adjustl(grid_select)) == 'lonlat') then 
        if (('.nc' == lndname(length-2:length)) .or. &  
            ('nc4' == lndname(length-2:length))) then
            stop 'ERROR! can not use now in the lonlat' 

        else if ('nml' == lndname(length-2:length)) then 
            namelist /lonlatmesh/ mesh 
            
            open(10, status = 'OLD', file = lndname) 
            REWIND(10)
            read(10, nml = lonlatmesh) 
            close(10)
            write(*, nml = lonlatmesh) 

            definition             = mesh%definition
            lon_start              = mesh%lon_start
            lon_end                = mesh%lon_end
            lon_grid_interval      = mesh%lon_grid_interval
            lon_points             = mesh%lon_points
            lat_start              = mesh%lat_start
            lat_end                = mesh%lat_end
            lat_grid_interval      = mesh%lat_grid_interval
            lat_points             = mesh%lat_points
            
            nlon_cal = int( (lon_end - lon_start) / lon_grid_interval) + 1 
            nlat_cal = int( (lat_end - lat_start) / lat_grid_interval) + 1 
            if (nlon_cal /= lon_points) then
                write(io6, *) "nlon_cal = ", nlon_cal
                write(io6, *) "lon_points = ", lon_points
                stop 'ERROR nlon_cal /= lon_points'
            else if (nlat_cal /= lat_points) then
                write(io6, *) "nlat_cal = ", nlat_cal
                write(io6, *) "lat_points = ", lat_points
                stop 'ERROR nlat_cal /= lat_points'
            end if
            
            if (definition == 'center') then 
                allocate(lon_center(lon_points))
                allocate(lat_center(lat_points))
                allocate(lon_bound(lon_points + 1))
                allocate(lat_bound(lat_points + 1))
                lon_center = lon_start + lon_grid_interval * ([1:lon_points] - 1) 
                lat_center = lat_start + lat_grid_interval * ([1:lat_points] - 1) 
                lon_bound(1 : lon_points) = lon_center - lon_grid_interval / 2.0  
                lon_bound(1 + lon_points) = lon_end    + lon_grid_interval / 2.0  
                lat_bound(1 : lat_points) = lat_center - lat_grid_interval / 2.0  
                lat_bound(1 + lat_points) = lat_end    + lat_grid_interval / 2.0  
            else if (definition == 'bound') then 
                allocate(lon_bound(lon_points)) 
                allocate(lat_bound(lat_points)) 
                
                lon_points = lon_points - 1 
                lat_points = lat_points - 1 
                allocate(lon_center(lon_points))
                allocate(lat_center(lat_points))
                lon_bound(1 : lon_points) = lon_start + lon_grid_interval * ([1:lon_points] - 1) 
                lon_bound(1 + lon_points) = lon_end   
                lat_bound(1 : lat_points) = lat_start + lat_grid_interval * ([1:lat_points] - 1)  
                lat_bound(1 + lat_points) = lat_end    
                lon_center = (lon_bound(1 : lon_points) + lon_bound(2 : 1 + lon_points)) / 2.0 
                lat_center = (lat_bound(1 : lat_points) + lat_bound(2 : 1 + lat_points)) / 2.0 
            else
                stop "error definition mismatch" 
            end if
        else
            stop "error datatype must choose 'ncfile' or 'namelist' please check!" 
        end if
    
        write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(io6, *) "lonlat bound range"
        write(io6, *) "lon_bound(1) = ", lon_bound(1)
        write(io6, *) "lon_bound(1 + lon_points) = ", lon_bound(1 + lon_points)
        write(io6, *) "lat_bound(1) = ", lat_bound(1)
        write(io6, *) "lat_bound(1 + lat_points) = ", lat_bound(1 + lat_points)
        write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(io6, *) ""
        
        write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(io6, *) "lon_points = ", lon_points
        write(io6, *) "lat_points = ", lat_points
        write(io6, *) "lonlat start end interval"
        write(io6, *) "lon_center start from :", lon_center(1)
        write(io6, *) "lon_center start end  :", lon_center(lon_points)
        write(io6, *) "lon_grid_interval     :", lon_grid_interval
        write(io6, *) "lat_center start from :", lat_center(1)
        write(io6, *) "lat_center start end  :", lat_center(lat_points)
        write(io6, *) "lat_grid_interval     :", lat_grid_interval
        write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(io6, *) ""
        
        where(lon_center > 180.)  lon_center = lon_center - 360. 
        where(lon_bound  > 180.)  lon_bound  = lon_bound  - 360. 
        
        NXP = int( abs(360. / lon_grid_interval) / 5 )

        write(io6, *) "turn lon_bound/lat_bound(1D) to lon_vert/lat_vert(2D) start"
        allocate(lon_vert(lon_points + 1, lat_points + 1)); lon_vert = 0.
        allocate(lat_vert(lon_points + 1, lat_points + 1)); lat_vert = 0.
        do j = 1, lat_points + 1, 1
            lon_vert(:, j) = lon_bound 
        end do
        do i = 1, lon_points + 1, 1
            lat_vert(i, :) = lat_bound 
        end do
        write(io6, *) "turn lon_bound/lat_bound(1D) to lon_vert/lat_vert(2D) finish"
        write(io6, *) ""

    else if (trim(adjustl(grid_select)) == 'lambert') then 
        if (('.nc' == lndname(length-2:length))  .or. &  
            ('nc4' == lndname(length-2:length))) then
            CALL CHECK(NF90_OPEN(lndname, nf90_nowrite, ncid)) 
            CALL CHECK(NF90_INQ_DIMID(ncid, "xi_vert",  dimID_lon)) 
            CALL CHECK(NF90_INQ_DIMID(ncid, "eta_vert", dimID_lat)) 
            CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lon, len = lon_points)) 
            CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lat, len = lat_points)) 
            allocate(lon_vert(lon_points, lat_points)) 
            allocate(lat_vert(lon_points, lat_points))
            CALL CHECK(NF90_INQ_VARID(ncid, "lon_vert", varid(1))) 
            CALL CHECK(NF90_INQ_VARID(ncid, "lat_vert", varid(2))) 
            CALL CHECK(NF90_GET_VAR(ncid, varid(1), lon_vert))     
            CALL CHECK(NF90_GET_VAR(ncid, varid(2), lat_vert))     
            CALL CHECK(NF90_CLOSE(ncid)) 
            
            lon_points = lon_points - 1 
            lat_points = lat_points - 1 
            where(lon_vert > 180.)  lon_vert = lon_vert - 360. 
        else if ('nml' == lndname(length-2:length)) then 
            stop 'ERROR! nml can not use in lambert now' 
        end if

    else if (trim(adjustl(grid_select)) == 'cubical') then 
        stop "error lambert or cubical can not use now!" 
    else
        stop "error grid_select must choose 'lonlat' or 'lambert' or 'cubical' please check!" 
    end if

    write(io6, *) "lonlat_bound and ngr_bound calculate start"
    write(io6, *) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    
    bound_points = (lon_points + 1) * (lat_points + 1) + 1 
    mode_points  = lon_points * lat_points + 1             
    allocate(lonlat_bound(bound_points, 2)); lonlat_bound = -999. 
    allocate(ngr_bound(4, mode_points));     ngr_bound = 1        
    allocate(n_ngr(mode_points));            n_ngr = 4            

    idx = 1 
    do j = 1, lat_points + 1, 1 
        do i = 1, lon_points + 1, 1 
            idx = idx + 1
            lonlat_bound(idx, :) = [lon_vert(i, j), lat_vert(i, j)] 
        end do
    end do
    
    idx = 1 
    do j = 1, lat_points, 1 
        do i = 1, lon_points, 1 
            idx = idx + 1
            ngr_bound(:, idx) = [i + (j - 1) * (lon_points + 1),     & 
                                 i + (j - 1) * (lon_points + 1) + 1, & 
                                 i +  j      * (lon_points + 1) + 1, & 
                                 i +  j      * (lon_points + 1)]       
        end do
    end do
    ngr_bound = ngr_bound + 1 

    write(io6, *) "lon_points : ", lon_points
    write(io6, *) "lat_points : ", lat_points
    write(io6, *) "bound_points : ", bound_points
    write(io6, *) "mode_points  : ", mode_points
    write(io6, *) "minval(lonlat_bound) : ", minval(lonlat_bound(2:bound_points, 1))
    write(io6, *) "maxval(lonlat_bound) : ", maxval(lonlat_bound(2:bound_points, 1))
    write(io6, *) "maxval(lonlat_bound) : ", maxval(lonlat_bound(2:bound_points, 2))
    write(io6, *) "minval(lonlat_bound) : ", minval(lonlat_bound(2:bound_points, 2))
    write(io6, *) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    write(io6, *) "lonlat_bound and ngr_bound calculate finish"
    write(io6, *) ""

    write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(io6, *) 'grid_write: opening file:', lndname 
    write(io6, *) 'mode_points : ', mode_points
    write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    
    write(nxpc, '(I4.4)') NXP    
    write(stepc, '(I2.2)') step  
    lndname = trim(file_dir)// 'gridfile/gridfile_NXP' // trim(nxpc) // '_'// trim(stepc) //'_'//trim(grid_select)//'.nc4'
    write(io6, *) lndname
    CALL Mode4_Mesh_Save(lndname, bound_points, mode_points, lonlat_bound, ngr_bound, n_ngr)
    write(io6, *) "mode4mesh_make finish"

end subroutine mode4mesh_make

SUBROUTINE Mode4_Mesh_Save(lndname, bound_points, mode_points, lonlat_bound, ngr_bound, n_ngr)
    IMPLICIT NONE
    character(pathlen), intent(in) :: lndname
    integer :: ncid, dimID_bound, dimID_points, dimID_two, dimID_four, varid(3)
    integer, intent(in) :: bound_points, mode_points
    real(r8), dimension(:, :), allocatable, intent(in) :: lonlat_bound
    integer,  dimension(:, :), allocatable, intent(in) :: ngr_bound
    integer,  dimension(:),    allocatable, intent(in) :: n_ngr

    CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
    CALL CHECK(NF90_DEF_DIM(ncid, "bound_points", bound_points, DimID_bound))
    CALL CHECK(NF90_DEF_DIM(ncid, "mode_points", mode_points, DimID_points))
    CALL CHECK(NF90_DEF_DIM(ncid, "two", 2, DimID_two))
    CALL CHECK(NF90_DEF_DIM(ncid, "four", 4, DimID_four))
    CALL CHECK(NF90_DEF_VAR(ncid, "lonlat_bound", NF90_FLOAT, (/ DimID_bound, DimID_two /), varid(1)))
    CALL CHECK(NF90_DEF_VAR(ncid, "ngr_bound", NF90_INT, (/ DimID_four, DimID_points /), varid(2)))
    CALL CHECK(NF90_DEF_VAR(ncid, "n_ngr", NF90_INT, (/ DimID_points /), varid(3)))
    CALL CHECK(NF90_ENDDEF(ncid))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(1), lonlat_bound))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(2), ngr_bound))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(3), n_ngr))
    CALL CHECK(NF90_CLOSE(ncid))
    write(io6, *)   "mode4 mesh save finish"

END SUBROUTINE Mode4_Mesh_Save

SUBROUTINE Mode4_Mesh_Read(lndname, bound_points, mode_points, lonlat_bound, ngr_bound, n_ngr)
    IMPLICIT NONE
    integer :: i, ncid, dimID_bound, dimID_points
    integer, dimension(3) :: varid
    character(pathlen), intent(in) :: lndname
    integer,  intent(out) :: bound_points, mode_points
    real(r8), dimension(:, :), allocatable, intent(out) :: lonlat_bound
    integer,  dimension(:, :), allocatable, intent(out) :: ngr_bound
    integer,  dimension(:),    allocatable, intent(out) :: n_ngr
    CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid))! 1. NF90_OPEN 打开文件 
    CALL CHECK(NF90_INQ_DIMID(ncid, "bound_points", dimID_bound))
    CALL CHECK(NF90_INQ_DIMID(ncid, "mode_points", dimID_points))! 2. NF90_INQ_DIMID获取维度ID 
    CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_bound,  len = bound_points))
    CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_points, len = mode_points))! 3. NF90_INQUIRE_DIMENSION获取各维度的长度
    allocate(lonlat_bound(bound_points, 2))
    allocate(ngr_bound(4, mode_points))
    allocate(n_ngr(mode_points))
    CALL CHECK(NF90_INQ_VARID(ncid, 'lonlat_bound', varid(1)))
    CALL CHECK(NF90_INQ_VARID(ncid, 'ngr_bound', varid(2)))
    CALL CHECK(NF90_INQ_VARID(ncid, 'n_ngr', varid(3)))
    CALL CHECK(NF90_GET_VAR(ncid, varid(1), lonlat_bound))
    CALL CHECK(NF90_GET_VAR(ncid, varid(2), ngr_bound))
    CALL CHECK(NF90_GET_VAR(ncid, varid(3), n_ngr))
    CALL CHECK(NF90_CLOSE(ncid))! 7. NF90_CLOSE关闭文件

 END SUBROUTINE Mode4_Mesh_Read


 SUBROUTINE Unstructured_Mesh_Read(lndname, sjx_points, lbx_points, mp, wp, ngrmw, ngrwm, n_ngrwm)
    IMPLICIT NONE
    integer :: i, ncid, dimID_sjx, dimID_lbx
    integer, dimension(7) :: varid
    character(pathlen), intent(in) :: lndname
    integer,  intent(out) :: sjx_points, lbx_points
    integer,  dimension(:), allocatable, intent(out) :: n_ngrwm
    integer,  dimension(:, :), allocatable, intent(out) :: ngrmw, ngrwm
    real(r8), dimension(:, :), allocatable, intent(out) :: mp, wp

    CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid))! 1. NF90_OPEN 打开文件 
    CALL CHECK(NF90_INQ_DIMID(ncid, "sjx_points", dimID_sjx))! 2. NF90_INQ_DIMID获取维度ID 
    CALL CHECK(NF90_INQ_DIMID(ncid, "lbx_points", dimID_lbx))! 三个参数依次是：文件编号，NC文件中的维度名，维度ID
    CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_sjx, len = sjx_points))! 3. NF90_INQUIRE_DIMENSION获取各维度的长度
    CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_lbx, len = lbx_points))! 三个参数依次是：文件编号，维度ID，维度名
    allocate(mp(sjx_points, 2))
    allocate(wp(lbx_points, 2))            
    allocate(ngrmw(3, sjx_points))
    allocate(ngrwm(7, lbx_points))      
    allocate(n_ngrwm( lbx_points)) ! 4. ALLOCATE分配动态数组
    

    CALL CHECK(NF90_INQ_VARID(ncid, 'GLONM', varid(1)))
    CALL CHECK(NF90_INQ_VARID(ncid, 'GLATM', varid(2)))
    CALL CHECK(NF90_INQ_VARID(ncid, 'GLONW', varid(3)))
    CALL CHECK(NF90_INQ_VARID(ncid, 'GLATW', varid(4)))
    CALL CHECK(NF90_INQ_VARID(ncid, 'itab_m%iw', varid(5)))
    CALL CHECK(NF90_INQ_VARID(ncid, 'itab_w%im', varid(6)))
    CALL CHECK(NF90_INQ_VARID(ncid, 'n_ngrwm', varid(7)))

    CALL CHECK(NF90_GET_VAR(ncid, varid(1), mp(:, 1)))! 6. NF90_GET_VAR：文件编号，变量ID，程序中的变量名
    CALL CHECK(NF90_GET_VAR(ncid, varid(2), mp(:, 2)))
    CALL CHECK(NF90_GET_VAR(ncid, varid(3), wp(:, 1)))
    CALL CHECK(NF90_GET_VAR(ncid, varid(4), wp(:, 2)))
    CALL CHECK(NF90_GET_VAR(ncid, varid(5), ngrmw))
    CALL CHECK(NF90_GET_VAR(ncid, varid(6), ngrwm))
    CALL CHECK(NF90_GET_VAR(ncid, varid(7), n_ngrwm))
    CALL CHECK(NF90_CLOSE(ncid))! 7. NF90_CLOSE关闭文件

 END SUBROUTINE Unstructured_Mesh_Read

 SUBROUTINE Unstructured_Mesh_Save(lndname, sjx_points, lbx_points, mp, wp, ngrmw, ngrwm, n_ngrwm)

    IMPLICIT NONE
    integer :: i, j, ncid, spDimID, lpDimID, thDimID, seDimID
    integer, dimension(7) :: varid
    character(pathlen), intent(in) :: lndname
    integer,  intent(in) :: sjx_points, lbx_points
    real(r8), dimension(:, :), allocatable, intent(in) :: mp, wp
    integer,  dimension(:, :), allocatable, intent(in) :: ngrmw, ngrwm
    integer,  dimension(:),    intent(in), optional :: n_ngrwm
    
    CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
    CALL CHECK(NF90_DEF_DIM(ncid, "sjx_points", sjx_points, spDimID))
    CALL CHECK(NF90_DEF_DIM(ncid, "lbx_points", lbx_points, lpDimID))
    CALL CHECK(NF90_DEF_DIM(ncid, "dimb", 3, thDimID))
    if (present(n_ngrwm)) then
       CALL CHECK(NF90_DEF_DIM(ncid, "dimc", 10, seDimID))
       write(io6, *)   "maxval(n_ngrwm) = ", maxval(n_ngrwm)
    else
       CALL CHECK(NF90_DEF_DIM(ncid, "dimc", 7, seDimID))
    end if

    CALL CHECK(NF90_DEF_VAR(ncid, "GLONM", NF90_FLOAT, (/ spDimID /), varid(1)))
    CALL CHECK(NF90_DEF_VAR(ncid, "GLATM", NF90_FLOAT, (/ spDimID /), varid(2)))
    CALL CHECK(NF90_DEF_VAR(ncid, "GLONW", NF90_FLOAT, (/ lpDimID /), varid(3)))
    CALL CHECK(NF90_DEF_VAR(ncid, "GLATW", NF90_FLOAT, (/ lpDimID /), varid(4)))
    CALL CHECK(NF90_DEF_VAR(ncid, "itab_m%iw", NF90_INT, (/ thDimID, spDimID /), varid(5)))
    CALL CHECK(NF90_DEF_VAR(ncid, "itab_w%im", NF90_INT, (/ seDimID, lpDimID /), varid(6)))
    if (present(n_ngrwm)) CALL CHECK(NF90_DEF_VAR(ncid, "n_ngrwm", NF90_INT, (/ lpDimID /), varid(7)))
    CALL CHECK(NF90_ENDDEF(ncid))

    CALL CHECK(NF90_PUT_VAR(ncid, varid(1), mp(1:sjx_points, 1)))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(2), mp(1:sjx_points, 2)))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(3), wp(1:lbx_points, 1)))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(4), wp(1:lbx_points, 2)))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(5), ngrmw(:, 1:sjx_points)))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(6), ngrwm(:, 1:lbx_points)))
    if (present(n_ngrwm)) CALL CHECK(NF90_PUT_VAR(ncid, varid(7), n_ngrwm))
    CALL CHECK(NF90_CLOSE(ncid))

 END SUBROUTINE Unstructured_Mesh_Save


 SUBROUTINE bbox_Mesh_Read(lndname, bbox_num, bbox_points)
    ! 读取bboxmesh数据
    IMPLICIT NONE
    character(pathlen), intent(in) :: lndname
    integer,  intent(out) :: bbox_num
    real(r8), allocatable, intent(out) :: bbox_points(:,:)
    integer :: ncid, DimID_num, varid(2)

    CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid))! 1. NF90_OPEN 打开文件 
    CALL CHECK(NF90_INQ_DIMID(ncid, "bbox_num", DimID_num))! 2. NF90_INQ_DIMID获取维度ID 
    CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, DimID_num, len = bbox_num))! 3. NF90_INQUIRE_DIMENSION获取各维度的长度
    allocate(bbox_points(bbox_num, 4)); bbox_points = 0.
    CALL CHECK(NF90_INQ_VARID(ncid, 'bbox_points',   varid(1)))
    CALL CHECK(NF90_GET_VAR(ncid, varid(1), bbox_points))
    CALL CHECK(NF90_CLOSE(ncid))! 7. NF90_bbox关闭文件
    write(io6, *)   "bbox mesh read finish"

 END SUBROUTINE bbox_Mesh_Read

 SUBROUTINE bbox_Mesh_Save(lndname, bbox_num, bbox_points)
    ! 存储bboxmesh数据
    IMPLICIT NONE
    character(pathlen), intent(in) :: lndname
    integer,  intent(in) :: bbox_num
    real(r8), allocatable, intent(in) :: bbox_points(:,:)
    integer :: ncid, DimID_num, DimID_four, varid(2)

    CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
    CALL CHECK(NF90_DEF_DIM(ncid, "bbox_num", bbox_num, DimID_num))
    CALL CHECK(NF90_DEF_DIM(ncid, "four", 4, DimID_four))
    CALL CHECK(NF90_DEF_VAR(ncid, "bbox_points", NF90_FLOAT, (/ DimID_num, DimID_four /), varid(1)))
    CALL CHECK(NF90_ENDDEF(ncid))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(1), bbox_points))
    CALL CHECK(NF90_CLOSE(ncid))
    write(io6, *)   "bbox mesh save finish"

 END SUBROUTINE bbox_Mesh_Save


 SUBROUTINE circle_Mesh_Read(lndname, circle_num, circle_points, circle_radius)
    ! 读取circlemesh数据
    IMPLICIT NONE
    character(pathlen), intent(in) :: lndname
    integer,  intent(out) :: circle_num
    real(r8), allocatable, intent(out) :: circle_points(:,:), circle_radius(:)
    integer :: ncid, DimID_num, varid(3)

    CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid))! 1. NF90_OPEN 打开文件 
    CALL CHECK(NF90_INQ_DIMID(ncid, "circle_num", DimID_num))! 2. NF90_INQ_DIMID获取维度ID 
    CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, DimID_num, len = circle_num))! 3. NF90_INQUIRE_DIMENSION获取各维度的长度
    
    allocate(circle_points(circle_num, 2)); circle_points = 0.
    allocate(circle_radius(circle_num));    circle_radius = 0.

    CALL CHECK(NF90_INQ_VARID(ncid, 'circle_points',   varid(1)))
    CALL CHECK(NF90_INQ_VARID(ncid, 'circle_radius',   varid(2)))
    CALL CHECK(NF90_GET_VAR(ncid, varid(1), circle_points))
    CALL CHECK(NF90_GET_VAR(ncid, varid(2), circle_radius))
    CALL CHECK(NF90_CLOSE(ncid))! 7. NF90_CLOSE关闭文件
    write(io6, *)   "circle mesh read finish"

 END SUBROUTINE circle_Mesh_Read

 SUBROUTINE circle_Mesh_Save(lndname, circle_num, circle_points, circle_radius)
    ! 存储circlemesh数据
    IMPLICIT NONE
    character(pathlen), intent(in) :: lndname
    integer,  intent(in) :: circle_num
    real(r8), allocatable, intent(in) :: circle_points(:,:), circle_radius(:)
    integer :: ncid, DimID_num, DimID_two, varid(3)

    CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
    CALL CHECK(NF90_DEF_DIM(ncid, "circle_num", circle_num, DimID_num))
    CALL CHECK(NF90_DEF_DIM(ncid, "two", 2, DimID_two))
    CALL CHECK(NF90_DEF_VAR(ncid, "circle_points", NF90_FLOAT, (/ DimID_num, DimID_two /), varid(1)))
    CALL CHECK(NF90_DEF_VAR(ncid, "circle_radius", NF90_FLOAT, (/ DimID_num/), varid(2)))
    CALL CHECK(NF90_ENDDEF(ncid))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(1), circle_points))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(2), circle_radius))
    CALL CHECK(NF90_CLOSE(ncid))
    write(io6, *)   "circle mesh save finish"

 END SUBROUTINE circle_Mesh_Save

 SUBROUTINE close_Mesh_Read(lndname, close_num, close_points)
    ! 读取closemesh数据
    IMPLICIT NONE
    character(pathlen), intent(in) :: lndname
    integer,  intent(out) :: close_num
    real(r8), allocatable, intent(out) :: close_points(:,:)
    integer :: ncid, DimID_num, varid(2)

    CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid))! 1. NF90_OPEN 打开文件 
    CALL CHECK(NF90_INQ_DIMID(ncid, "close_num", DimID_num))! 2. NF90_INQ_DIMID获取维度ID 
    CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, DimID_num, len = close_num))! 3. NF90_INQUIRE_DIMENSION获取各维度的长度
    allocate(close_points(close_num, 2)); close_points = 0.

    CALL CHECK(NF90_INQ_VARID(ncid, 'close_points',   varid(1)))
    CALL CHECK(NF90_GET_VAR(ncid, varid(1), close_points))
    CALL CHECK(NF90_CLOSE(ncid))! 7. NF90_CLOSE关闭文件
    write(io6, *)   "close mesh read finish"

 END SUBROUTINE close_Mesh_Read

 SUBROUTINE close_Mesh_Save(lndname, close_num, close_points)
    ! 存储closemesh数据
    IMPLICIT NONE
    character(pathlen), intent(in) :: lndname
    integer,  intent(in) :: close_num
    real(r8), allocatable, intent(in) :: close_points(:,:)
    integer :: ncid, DimID_num, DimID_two, varid(2)

    CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
    CALL CHECK(NF90_DEF_DIM(ncid, "close_num", close_num, DimID_num))
    CALL CHECK(NF90_DEF_DIM(ncid, "two", 2, DimID_two))
    CALL CHECK(NF90_DEF_VAR(ncid, "close_points", NF90_FLOAT, (/ DimID_num, DimID_two /), varid(1)))
    CALL CHECK(NF90_ENDDEF(ncid))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(1), close_points))
    CALL CHECK(NF90_CLOSE(ncid))
    write(io6, *)   "close mesh save finish"

 END SUBROUTINE close_Mesh_Save

 ! 这个是应该只有三角形的情况，还需要检查一些
 SUBROUTINE Contain_Read(lndname, num_ustr, num_ii, ustr_id, ustr_ii, IsInArea_ustr)

    IMPLICIT NONE

    character(pathlen), intent(in) :: lndname
    integer,  intent(out) :: num_ustr, num_ii
    integer,  dimension(:,:), allocatable, intent(out) :: ustr_id, ustr_ii
    integer,  dimension(:),   allocatable, intent(out), optional :: IsInArea_ustr
    integer :: ncid, Dim_ustrID, Dim_iiID, Dim_aID, Dim_bID
    integer :: dim_a_len, dim_b_len
    integer :: varid(3)

    CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid))! 1. NF90_OPEN 打开文件 
    CALL CHECK(NF90_INQ_DIMID(ncid, "num_ustr", Dim_ustrID))! 2. NF90_INQ_DIMID获取维度ID 
    CALL CHECK(NF90_INQ_DIMID(ncid, "num_ii", Dim_iiID))! 三个参数依次是：文件编号，NC文件中的维度名，维度ID
    CALL CHECK(NF90_INQ_DIMID(ncid, "dim_a", Dim_aID))! 三个参数依次是：文件编号，NC文件中的维度名，维度ID
    CALL CHECK(NF90_INQ_DIMID(ncid, "dim_b", Dim_bID))! 三个参数依次是：文件编号，NC文件中的维度名，维度ID
    CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, Dim_ustrID, len = num_ustr))! 3. NF90_INQUIRE_DIMENSION获取各维度的长度
    CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, Dim_iiID, len = num_ii))! 三个参数依次是：文件编号，维度ID，维度名
    CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, Dim_aID, len = dim_a_len))! 三个参数依次是：文件编号，维度ID，维度名
    CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, Dim_bID, len = dim_b_len))! 三个参数依次是：文件编号，维度ID，维度名
    allocate(ustr_id(num_ustr, dim_a_len)) ! 4. ALLOCATE分配动态数组         
    allocate(ustr_ii(num_ii, dim_b_len))     
    CALL CHECK(NF90_INQ_VARID(ncid, 'ustr_id',   varid(1)))
    CALL CHECK(NF90_INQ_VARID(ncid, 'ustr_ii',   varid(2)))
    CALL CHECK(NF90_GET_VAR(ncid, varid(1), ustr_id))
    CALL CHECK(NF90_GET_VAR(ncid, varid(2), ustr_ii))
    if (present(IsInArea_ustr)) then
       allocate(IsInArea_ustr(num_ustr))
       CALL CHECK(NF90_INQ_VARID(ncid, 'IsInArea_ustr',   varid(3)))
       CALL CHECK(NF90_GET_VAR(ncid, varid(3), IsInArea_ustr))
    end if
    CALL CHECK(NF90_CLOSE(ncid))! 7. NF90_CLOSE关闭文件

 END SUBROUTINE Contain_Read

 SUBROUTINE Contain_Save(lndname, num_ustr, num_ii, ustr_id, ustr_ii, IsInArea_ustr)
    IMPLICIT NONE

    character(pathlen), intent(in) :: lndname
    integer, intent(in) :: num_ustr, num_ii
    integer,  dimension(:,:), allocatable, intent(in) :: ustr_id, ustr_ii
    integer,  dimension(:), allocatable, intent(in), optional :: IsInArea_ustr
    integer :: ncid, Dim_ustrID, Dim_iiID, Dim_aID, Dim_bID
    integer :: dim_a_len, dim_b_len
    integer :: varid(3)

    CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
    CALL CHECK(NF90_DEF_DIM(ncid, "num_ustr", num_ustr, Dim_ustrID))
    CALL CHECK(NF90_DEF_DIM(ncid, "num_ii", num_ii, Dim_iiID))
    dim_a_len = size(ustr_id, 2) ! 因为陆地网格有两列，海洋网格有三列, 海陆网格有四列
    dim_b_len = size(ustr_ii, 2) ! 因为陆地和海洋网格有两列，但是海陆网格有四列
    CALL CHECK(NF90_DEF_DIM(ncid, "dim_a", dim_a_len, Dim_aID)) ! *_id
    CALL CHECK(NF90_DEF_DIM(ncid, "dim_b", dim_b_len, Dim_bID))
    CALL CHECK(NF90_DEF_VAR(ncid, 'ustr_id', NF90_INT, (/ Dim_ustrID, Dim_aID /), varid(1)))
    CALL CHECK(NF90_DEF_VAR(ncid, 'ustr_ii', NF90_INT, (/ Dim_iiID, Dim_bID /), varid(2)))
    CALL CHECK(NF90_DEF_VAR(ncid, 'IsInArea_ustr', NF90_INT, (/ Dim_ustrID /), varid(3)))
    CALL CHECK(NF90_ENDDEF(ncid))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(1), ustr_id))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(2), ustr_ii))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(3), IsInArea_ustr))
    CALL CHECK(NF90_CLOSE(ncid))

 END SUBROUTINE Contain_Save

 SUBROUTINE earthmesh_info_save(lndname, num_step, num_ustr, num_step_f, refine_degree_f, seaorland_ustr_f)
    IMPLICIT NONE

    character(pathlen), intent(in) :: lndname
    integer, intent(in) :: num_step, num_ustr
    integer,  dimension(:), allocatable, intent(in) :: num_step_f, refine_degree_f, seaorland_ustr_f
    integer :: ncid, Dim_stepID, Dim_ustrID
    integer :: varid(3)

    CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
    CALL CHECK(NF90_DEF_DIM(ncid, "num_step", num_step, Dim_stepID))
    CALL CHECK(NF90_DEF_DIM(ncid, "num_ustr", num_ustr, Dim_ustrID))
    CALL CHECK(NF90_DEF_VAR(ncid, 'num_step_f', NF90_INT, (/ Dim_stepID /), varid(1)))
    CALL CHECK(NF90_DEF_VAR(ncid, 'refine_degree_f', NF90_INT, (/ Dim_ustrID /), varid(2)))
    CALL CHECK(NF90_DEF_VAR(ncid, 'seaorland_ustr_f', NF90_INT, (/ Dim_ustrID /), varid(3)))
    CALL CHECK(NF90_ENDDEF(ncid))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(1), num_step_f))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(2), refine_degree_f))
    CALL CHECK(NF90_PUT_VAR(ncid, varid(3), seaorland_ustr_f))
    CALL CHECK(NF90_CLOSE(ncid))

 END SUBROUTINE earthmesh_info_save

end module MOD_utilities 


