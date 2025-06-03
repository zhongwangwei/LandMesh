module MOD_file_preprocess
   ! 可以考虑把labdtype数据存在这里也是挺好的，还有一个问题
   USE consts_coms, only: r8, pathlen
   USE netcdf
   implicit none

   contains
   
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
      print*, "mode4 mesh save finish"

   END SUBROUTINE Mode4_Mesh_Save

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
         print*, "maxval(n_ngrwm) = ", maxval(n_ngrwm)
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
      print*, "bbox mesh read finish"

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
      print*, "bbox mesh save finish"

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
      print*, "circle mesh read finish"

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
      print*, "circle mesh save finish"

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
      print*, "close mesh read finish"

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
      print*, "close mesh save finish"

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

END Module MOD_file_preprocess
