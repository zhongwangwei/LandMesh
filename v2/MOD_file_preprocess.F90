module MOD_file_preprocess
   ! 可以考虑把labdtype数据存在这里也是挺好的，还有一个问题
   USE consts_coms, only: r8
   USE netcdf
   implicit none

   contains

   ! only use for threshold values for FHW NXP144 one step refine
   SUBROUTINE Unstructured_Threshold_Read(lndname, sjx_points, f_mainarea, num_landtypes, p_lai)
      IMPLICIT NONE
      character(len = 256), intent(in) :: lndname
      integer :: ncid, dimID_sjx
      integer, dimension(3) :: varid
      integer,  intent(out) :: sjx_points
      integer,  dimension(:), allocatable, intent(out) :: num_landtypes
      real(r8), dimension(:, :), allocatable, intent(out) :: f_mainarea, p_lai

      CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid))! 1. NF90_OPEN 打开文件 
      CALL CHECK(NF90_INQ_DIMID(ncid, "sjx_points", dimID_sjx))! 2. NF90_INQ_DIMID获取维度ID 
      CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_sjx, len = sjx_points))! 3. NF90_INQUIRE_DIMENSION获取各维度的长度
      allocate(f_mainarea(sjx_points, 2))
      allocate(num_landtypes(sjx_points))
      allocate(p_lai(sjx_points, 2))
      CALL CHECK(NF90_INQ_VARID(ncid, 'fraction_mainarea', varid(1)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'num_landtypes', varid(2)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'p_lai', varid(3)))
      CALL CHECK(NF90_GET_VAR(ncid, varid(1), f_mainarea))
      CALL CHECK(NF90_GET_VAR(ncid, varid(2), num_landtypes))! 6. NF90_GET_VAR：文件编号，变量ID，程序中的变量名
      CALL CHECK(NF90_GET_VAR(ncid, varid(3), p_lai))
      CALL CHECK(NF90_CLOSE(ncid))! 7. NF90_CLOSE关闭文件
  
   END SUBROUTINE Unstructured_Threshold_Read

   SUBROUTINE Mode4_Mesh_Read(lndname, bound_points, mode4_points, lonlat_bound, ngr_bound)!, lonlat_bound
      IMPLICIT NONE
      integer :: i, ncid, dimID_bound, dimID_points
      integer, dimension(4) :: varid
      character(len = 256), intent(in) :: lndname
      integer,  intent(out) :: bound_points, mode4_points
      real(r8), dimension(:, :), allocatable, intent(out) :: lonlat_bound
      integer,  dimension(:, :), allocatable, intent(out) :: ngr_bound
      CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid))! 1. NF90_OPEN 打开文件 
      CALL CHECK(NF90_INQ_DIMID(ncid, "bound_points", dimID_bound))
      CALL CHECK(NF90_INQ_DIMID(ncid, "mode4_points", dimID_points))! 2. NF90_INQ_DIMID获取维度ID 
      CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_bound,  len = bound_points))
      CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, dimID_points, len = mode4_points))! 3. NF90_INQUIRE_DIMENSION获取各维度的长度
      allocate(lonlat_bound(bound_points, 2))
      allocate(ngr_bound(4, mode4_points))
      CALL CHECK(NF90_INQ_VARID(ncid, 'lonlat_bound', varid(1)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'ngr_bound', varid(2)))
      CALL CHECK(NF90_GET_VAR(ncid, varid(1), lonlat_bound))
      CALL CHECK(NF90_GET_VAR(ncid, varid(2), ngr_bound))
      CALL CHECK(NF90_CLOSE(ncid))! 7. NF90_CLOSE关闭文件
  
   END SUBROUTINE Mode4_Mesh_Read

   SUBROUTINE Mode4_Mesh_Save(lndname, bound_points, mode4_points, lonlat_bound, ngr_bound)

      IMPLICIT NONE
      
      character(len = 256), intent(in) :: lndname
      integer :: ncid, dimID_bound, dimID_points, dimID_two, dimID_four, varid(2)
      integer, intent(in) :: bound_points, mode4_points
      real(r8), dimension(:, :), allocatable, intent(in) :: lonlat_bound
      integer,  dimension(:, :), allocatable, intent(in) :: ngr_bound
 
      CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
      CALL CHECK(NF90_DEF_DIM(ncid, "bound_points", bound_points, DimID_bound))
      CALL CHECK(NF90_DEF_DIM(ncid, "mode4_points", mode4_points, DimID_points))
      CALL CHECK(NF90_DEF_DIM(ncid, "two", 2, DimID_two))
      CALL CHECK(NF90_DEF_DIM(ncid, "four", 4, DimID_four))
      CALL CHECK(NF90_DEF_VAR(ncid, "lonlat_bound", NF90_FLOAT, (/ DimID_bound, DimID_two /), varid(1)))
      CALL CHECK(NF90_DEF_VAR(ncid, "ngr_bound", NF90_INT, (/ DimID_four, DimID_points /), varid(2)))
      CALL CHECK(NF90_ENDDEF(ncid))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(1), lonlat_bound))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(2), ngr_bound))
      CALL CHECK(NF90_CLOSE(ncid))
      print*, "mode4 mesh save finish"

   END SUBROUTINE Mode4_Mesh_Save

   SUBROUTINE Unstructured_Mesh_Read(lndname, sjx_points, lbx_points, mp, wp, ngrmw, ngrwm, n_ngrwm)
      !p_name(8) = (/"GLONM", "GLATM", "GLONW", "GLATW", "itab_w%iw", "itab_m%im", "n_ngrwm"/)
      IMPLICIT NONE
      integer :: i, ncid, dimID_sjx, dimID_lbx
      integer, dimension(7) :: varid
      character(len = 256), intent(in) :: lndname
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
      ! allocate(n_ngrmw( sjx_points))   
      allocate(n_ngrwm( lbx_points)) ! 4. ALLOCATE分配动态数组
      

      CALL CHECK(NF90_INQ_VARID(ncid, 'GLONM', varid(1)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'GLATM', varid(2)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'GLONW', varid(3)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'GLATW', varid(4)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'itab_m%iw', varid(5)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'itab_w%im', varid(6)))
      ! CALL CHECK(NF90_INQ_VARID(ncid, 'n_ngrmw', varid(7)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'n_ngrwm', varid(7)))

      CALL CHECK(NF90_GET_VAR(ncid, varid(1), mp(:, 1)))! 6. NF90_GET_VAR：文件编号，变量ID，程序中的变量名
      CALL CHECK(NF90_GET_VAR(ncid, varid(2), mp(:, 2)))
      CALL CHECK(NF90_GET_VAR(ncid, varid(3), wp(:, 1)))
      CALL CHECK(NF90_GET_VAR(ncid, varid(4), wp(:, 2)))
      CALL CHECK(NF90_GET_VAR(ncid, varid(5), ngrmw))
      CALL CHECK(NF90_GET_VAR(ncid, varid(6), ngrwm))
      ! CALL CHECK(NF90_GET_VAR(ncid, varid(7), n_ngrmw))
      CALL CHECK(NF90_GET_VAR(ncid, varid(7), n_ngrwm))
      CALL CHECK(NF90_CLOSE(ncid))! 7. NF90_CLOSE关闭文件
  
   END SUBROUTINE Unstructured_Mesh_Read

   SUBROUTINE Unstructured_Mesh_Save(lndname, sjx_points, lbx_points, mp, wp, ngrmw, ngrwm)

      IMPLICIT NONE
      integer :: i, j, ncid, spDimID, lpDimID, thDimID, seDimID
      integer, dimension(7) :: varid
      integer,  dimension(:), allocatable :: n_ngrwm
      character(len = 256), intent(in) :: lndname
      integer,  intent(in) :: sjx_points, lbx_points
      integer,  dimension(:, :), allocatable, intent(in) :: ngrmw, ngrwm
      real(r8), dimension(:, :), allocatable, intent(in) :: mp, wp

      ! allocate(n_ngrmw(sjx_points)); n_ngrmw = 3
      allocate(n_ngrwm(lbx_points)); n_ngrwm = 1
      !do i = 2, sjx_points, 1
      !    do j = 1, 3, 1
      !        if (ngrmw(j, i) == 1) n_ngrmw = n_ngrmw - 1
      !    end do
      !end do
      do i = 2, lbx_points, 1 ! 多边形
         if (ngrwm(6, i) == 1) then
             n_ngrwm(i) = 5
         else if (ngrwm(7, i) == 1) then
             n_ngrwm(i) = 6
         else
             n_ngrwm(i) = 7
         end if
      end do

      CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
      CALL CHECK(NF90_DEF_DIM(ncid, "sjx_points", sjx_points, spDimID))
      CALL CHECK(NF90_DEF_DIM(ncid, "lbx_points", lbx_points, lpDimID))
      CALL CHECK(NF90_DEF_DIM(ncid, "dimb", 3, thDimID))
      CALL CHECK(NF90_DEF_DIM(ncid, "dimc", 7, seDimID))

      CALL CHECK(NF90_DEF_VAR(ncid, "GLONM", NF90_FLOAT, (/ spDimID /), varid(1)))
      CALL CHECK(NF90_DEF_VAR(ncid, "GLATM", NF90_FLOAT, (/ spDimID /), varid(2)))
      CALL CHECK(NF90_DEF_VAR(ncid, "GLONW", NF90_FLOAT, (/ lpDimID /), varid(3)))
      CALL CHECK(NF90_DEF_VAR(ncid, "GLATW", NF90_FLOAT, (/ lpDimID /), varid(4)))
      CALL CHECK(NF90_DEF_VAR(ncid, "itab_m%iw", NF90_INT, (/ thDimID, spDimID /), varid(5)))
      CALL CHECK(NF90_DEF_VAR(ncid, "itab_w%im", NF90_INT, (/ seDimID, lpDimID /), varid(6)))
      ! CALL CHECK(NF90_DEF_VAR(ncid, "n_ngrmw", NF90_INT, (/ spDimID /), varid(7)))
      CALL CHECK(NF90_DEF_VAR(ncid, "n_ngrwm", NF90_INT, (/ lpDimID /), varid(7)))
      CALL CHECK(NF90_ENDDEF(ncid))

      CALL CHECK(NF90_PUT_VAR(ncid, varid(1), mp(1:sjx_points, 1)))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(2), mp(1:sjx_points, 2)))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(3), wp(1:lbx_points, 1)))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(4), wp(1:lbx_points, 2)))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(5), ngrmw(:, 1:sjx_points)))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(6), ngrwm(:, 1:lbx_points)))
      ! CALL CHECK(NF90_PUT_VAR(ncid, varid(7), n_ngrmw))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(7), n_ngrwm))
      CALL CHECK(NF90_CLOSE(ncid))
      deallocate(n_ngrwm) 
   END SUBROUTINE Unstructured_Mesh_Save

   ! 这个是应该只有三角形的情况，还需要检查一些
   SUBROUTINE Contain_Read(lndname, num_ustr, num_ii, ustr_id, ustr_ii)

      IMPLICIT NONE

      character(len = 256), intent(in) :: lndname
      integer,  intent(inout) :: num_ustr, num_ii
      integer,  dimension(:,:), allocatable, intent(out) :: ustr_id, ustr_ii
      integer :: ncid, Dim_ustrID, Dim_iiID, Dim_aID
      integer :: varid(3)

      CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid))! 1. NF90_OPEN 打开文件 
      CALL CHECK(NF90_INQ_DIMID(ncid, "num_ustr", Dim_ustrID))! 2. NF90_INQ_DIMID获取维度ID 
      CALL CHECK(NF90_INQ_DIMID(ncid, "num_ii", Dim_iiID))! 三个参数依次是：文件编号，NC文件中的维度名，维度ID
      CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, Dim_ustrID, len = num_ustr))! 3. NF90_INQUIRE_DIMENSION获取各维度的长度
      CALL CHECK(NF90_INQUIRE_DIMENSION(ncid, Dim_iiID, len = num_ii))! 三个参数依次是：文件编号，维度ID，维度名
      
      allocate(ustr_id(num_ustr, 2)) ! 4. ALLOCATE分配动态数组         
      allocate(ustr_ii(num_ii, 2))     

      CALL CHECK(NF90_INQ_VARID(ncid, 'ustr_id',   varid(1)))
      CALL CHECK(NF90_INQ_VARID(ncid, 'ustr_ii',   varid(2)))

      CALL CHECK(NF90_GET_VAR(ncid, varid(1), ustr_id))
      CALL CHECK(NF90_GET_VAR(ncid, varid(2), ustr_ii))
      CALL CHECK(NF90_CLOSE(ncid))! 7. NF90_CLOSE关闭文件

   END SUBROUTINE Contain_Read

   SUBROUTINE Contain_Save(lndname, num_ustr, num_ii, ustr_id, ustr_ii, Min_matrix_index)
      ! 建议*_id, *_area, *_ii都放在一个文件里面就好了
      IMPLICIT NONE

      character(len = 256), intent(in) :: lndname
      integer, intent(in) :: num_ustr, num_ii
      integer,  dimension(:,:), allocatable, intent(in) :: ustr_id, ustr_ii, Min_matrix_index
      integer :: ncid, Dim_ustrID, Dim_iiID, Dim_aID, Dim_bID
      integer :: varid(3)

      CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
      CALL CHECK(NF90_DEF_DIM(ncid, "num_ustr", num_ustr, Dim_ustrID))
      CALL CHECK(NF90_DEF_DIM(ncid, "num_ii", num_ii, Dim_iiID))
      CALL CHECK(NF90_DEF_DIM(ncid, "dim_a", 2, Dim_aID)) ! *_id
      CALL CHECK(NF90_DEF_DIM(ncid, "dim_b", 4, Dim_bID))
      CALL CHECK(NF90_DEF_VAR(ncid, 'ustr_id', NF90_INT, (/ Dim_ustrID, Dim_aID /), varid(1)))
      CALL CHECK(NF90_DEF_VAR(ncid, 'ustr_ii', NF90_INT, (/ Dim_iiID, Dim_aID /), varid(2)))
      CALL CHECK(NF90_DEF_VAR(ncid, 'Min_matrix_index', NF90_INT, (/ Dim_bID, Dim_ustrID /), varid(3)))
      CALL CHECK(NF90_ENDDEF(ncid))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(1), ustr_id))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(2), ustr_ii))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(3), Min_matrix_index(:, 1:num_ustr)))
      CALL CHECK(NF90_CLOSE(ncid))
      
   END SUBROUTINE Contain_Save

   SUBROUTINE Mode4_Contain_Save(lndname, num_ustr, ustr_id, Min_matrix_index)
      IMPLICIT NONE

      character(len = 256), intent(in) :: lndname
      integer, intent(in) :: num_ustr
      integer,  dimension(:,:), allocatable, intent(in) :: ustr_id, Min_matrix_index
      integer :: ncid, Dim_ustrID, Dim_iiID, Dim_aID, Dim_bID
      integer :: varid(2)

      CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber, nf90_netcdf4), ncid))
      CALL CHECK(NF90_DEF_DIM(ncid, "num_ustr", num_ustr, Dim_ustrID))
      CALL CHECK(NF90_DEF_DIM(ncid, "dim_a", 2, Dim_aID)) ! *_id
      CALL CHECK(NF90_DEF_DIM(ncid, "dim_b", 4, Dim_bID))
      CALL CHECK(NF90_DEF_VAR(ncid, 'ustr_id', NF90_INT, (/ Dim_ustrID, Dim_aID /), varid(1)))
      CALL CHECK(NF90_DEF_VAR(ncid, 'Min_matrix_index', NF90_INT, (/ Dim_bID, Dim_ustrID /), varid(2)))
      CALL CHECK(NF90_ENDDEF(ncid))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(1), ustr_id))
      CALL CHECK(NF90_PUT_VAR(ncid, varid(2), Min_matrix_index(:, 1:num_ustr)))
      CALL CHECK(NF90_CLOSE(ncid))

   END SUBROUTINE Mode4_Contain_Save


END Module MOD_file_preprocess
