! 1. Triangular mesh multiple refinement
! 2. Calculate the refined inclusion relationship
module MOD_refine_sjx

   use netcdf
   USE consts_coms
   USE refine_vars

   implicit none

   contains
   SUBROUTINE refine_sjx ()
   
   implicit none
   
   integer :: num_ref                     ! The number of triangular meshes to be refined
   integer :: num_all                     ! The unstructured grid contains the total number of structured grids
   integer :: num_refjw                   ! The number of latitude and longitude grids to be refined
   integer :: refed(1000)                 ! Record the number of triangles that have been refined at each refinement
   integer :: nmp(-1:100),nwp(-1:100)     ! Record the number of m and w points after each refinement 
   integer :: iter                        ! 迭代次数 
   integer :: icl(3)                      ! Record the number of m and w points after each refinement
   !integer :: maxlc                       ! Land type maximum number, 17 or 24

   integer :: i,j,k,L,m,n
   integer :: w1,w2,w3,m1,m2,m3,m4,row,col
   integer :: sum_sea,sum_land            ! Number of land and sea grids

   integer :: sjx_points,lbx_points
   integer :: maxid(1)
   
   integer :: ncid,varid(10)
   integer :: spDimID,iunit,lpDimID,twDimID,thDimID,nmDimID,fvDimID
   integer :: idDimID,infoDimID,seDimID,foDimID,sxDimID,dimID_sjx,dimID_lbx

   real(r8),allocatable :: lat_i(:),lon_i(:)          ! Center point of latitude and longitude grid
   real(r8),allocatable :: wp(:,:),mp(:,:)            ! Store longitude, latitude, area, threshold refinement (before refinement) 
   real(r8),allocatable :: mp_new(:,:),wp_new(:,:)    ! Store longitude, latitude, area, threshold refinement (after refinement) 
   real(r8),allocatable :: mp_i(:,:,:)                ! Record the number of unstructured grids including latitude and longitude, including proportion, including area (in the process of refinement)
   real(r8),allocatable :: mp_ii(:,:)                 ! Record the number of unstructured grids including latitude and longitude, including proportion, including area (in the process of refinement)
   real(r8),allocatable :: mp_ii_new(:,:)             ! Record the number of unstructured grids including latitude and longitude, including proportion, including area (after refinement)

   real(r8),allocatable :: slope_max(:,:),slope_avg(:,:),p_slope(:,:)
   real(r8),allocatable :: area(:),area_fine_gridcell(:,:)
   real(r8),allocatable :: fraction_mainarea(:,:)
   real(r8),allocatable :: lai(:,:),p_lai(:,:)
   real(r8),allocatable :: k_s(:,:,:),p_k_s(:,:,:)
   real(r8),allocatable :: k_sl(:,:,:),p_k_sl(:,:,:)
   real(r8),allocatable :: tkdry(:,:,:),p_tkdry(:,:,:)
   real(r8),allocatable :: tksatf(:,:,:),p_tksatf(:,:,:)
   real(r8),allocatable :: tksatu(:,:,:),p_tksatu(:,:,:)

   real(r8) :: sjx(3,2),newsjx(3,2)
   real(r8) :: isinply           ! Inclusion relationship between unstructured grid and latitude and longitude grid (proportion)
   real(r8) :: dx,dy
   real(r8) :: tmpa(2),tmpb(2),tmpc(2)
   
   integer,allocatable :: n_landtypes(:)        ! The number of land types contained within the triangular grid
   real(r8),allocatable :: landtypes(:,:)        ! Latitude and longitude grid land types
   integer,allocatable :: ref(:)                ! Record whether the triangular mesh needs to be refined
   integer,allocatable :: ref_jw(:,:)           ! Records whether the latitude and longitude grid needs to calculate inclusion relationships
   integer,allocatable :: ref_tr(:,:)           ! Record what thresholds the triangular mesh is refined with
   integer,allocatable :: ngrmw(:,:)            ! Index of w points adjacent to m points (before thinning)
   integer,allocatable :: ngrmw_new(:,:)        ! Index of w points adjacent to m points (after refinement)
   integer,allocatable :: mp_id(:,:)            ! Contains the number of blocks, starting id in mp_ii (before thinning)
   integer,allocatable :: mp_id_new(:,:)        ! Contains the number of blocks, starting id in mp_ii (after refinement)
   integer,allocatable :: nla(:)                ! Whether the unstructured grid contains an indication of the land type 
   integer,allocatable :: seaorland(:,:)        ! Determine whether the latitude and longitude grid is ocean or land
   integer,allocatable :: ngrwm_new(:,:)        ! Index of m points adjacent to w points (not calculated)

   character(LEN=256) :: lndname,it,dir_output,nxpc
   
   logical :: isover       ! Determine whether the refinement is complete
   logical :: ispart       ! Determine whether the latitude and longitude mesh is fully contained by the unstructured mesh
   logical,allocatable :: IsInRfArea(:)         ! Determine whether the triangular mesh is in the refinement area
   
   character(LEN=20):: p_name(6)=(/"GLONW","GLATW","GLONM","GLATM","itab_w%im","itab_m%iw"/)

   isover = .false.

   dir_output = trim(base_dir) // trim(EXPNME)

!--------------------------------------------------------------------------
! 1.Read the data and allocate the array memory
!--------------------------------------------------------------------------
   print*,"Start reading unstructured grid data......"
   print*,""


   write(nxpc, '(I3.3)') NXP
   lndname = trim(base_dir) // trim(EXPNME) // '/makegrid/gridfile/gridfile_NXP' // trim(nxpc) // '.nc4'
   print*,lndname

   CALL CHECK(NF90_OPEN(trim(lndname),nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,p_name(1),varid(1)))
   CALL CHECK(NF90_INQ_VARID(ncid,p_name(2),varid(2)))
   CALL CHECK(NF90_INQ_VARID(ncid,p_name(3),varid(3)))
   CALL CHECK(NF90_INQ_VARID(ncid,p_name(4),varid(4)))
   CALL CHECK(NF90_INQ_VARID(ncid,p_name(5),varid(5)))
   CALL CHECK(NF90_INQ_VARID(ncid,p_name(6),varid(6)))
   CALL CHECK(NF90_INQ_DIMID(ncid,"sjx_points",dimID_sjx))
   CALL CHECK(NF90_INQ_DIMID(ncid,"lbx_points",dimID_lbx))
   CALL CHECK(NF90_INQUIRE_DIMENSION(ncid,dimID_lbx,len=lbx_points))
   CALL CHECK(NF90_INQUIRE_DIMENSION(ncid,dimID_sjx,len=sjx_points))

   print*,"sjx_points = ",sjx_points
   print*,"lbx_points = ",lbx_points

   nmp(-1) = 0
   nwp(-1) = 0
   nwp(0) = lbx_points
   nmp(0) = sjx_points
   allocate(wp(lbx_points,2))    
   allocate(wp_new(lbx_points*5,2))
   allocate(mp(sjx_points,4))          ! Three dimensions respectively store longitude, latitude, area, threshold refinement
   allocate(mp_id(sjx_points,2))       ! Contains the number of blocks, starting id in mp_ii
   allocate(mp_new(sjx_points*5,4))
   allocate(mp_id_new(sjx_points*5,2))
   allocate(mp_i(sjx_points*5,25000,4))
   allocate(ngrmw(3,sjx_points))
   allocate(ngrmw_new(3,sjx_points*5))
   allocate(ngrwm_new(7,lbx_points*5))
   allocate(ref(sjx_points*5))
   wp = 0.
   mp = 0.
   mp_i = 0.
   mp_new = 0.
   wp_new = 0.
   ngrmw = 0
   ngrmw_new = 0
   ngrwm_new = 0
   mp_id_new = 0
   mp_id = 0
   ref = 0

   CALL CHECK(NF90_GET_VAR(ncid,varid(1),wp(:,1)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(2),wp(:,2)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(3),mp(:,1)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(4),mp(:,2)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(6),ngrmw(1:3,:)))
   CALL CHECK(NF90_CLOSE(ncid))

   wp_new(1:lbx_points,1:2) = wp

   do i = 1,sjx_points,1
      do j = 1,3,1
         ngrmw_new(j,i) = ngrmw(j,i)
      end do
   end do

   allocate(landtypes(nlons_source,nlats_source))
   allocate(slope_max(nlons_source,nlats_source))          ! The maximum slope value in the eight neighborhood of the latitude and longitude grid
   allocate(slope_avg(nlons_source,nlats_source))          ! The mean slope value in the eight neighborhood of the latitude and longitude grid
   landtypes = 0.
   slope_max = 0.
   slope_avg = 0.

   lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/contain/initial/mp.bin"
   print*,lndname
   OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
   READ(iunit) mp(:,1:3)
   close(iunit)
   mp_new(1:sjx_points,1:3) = mp(:,1:3)

   lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/contain/initial/mp_id.bin"
   print*,lndname
   OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
   READ(iunit) mp_id(:,1:2)
   close(iunit)
   mp_id_new(1:sjx_points,1:2) = mp_id(:,1:2)

   num_all = INT(sum(mp_id(:,1)))
   allocate(mp_ii(num_all,4))
   allocate(mp_ii_new(num_all+1:num_all*10,4))

   lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/contain/initial/mp_ii.bin"
   print*,lndname
   OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
   READ(iunit) mp_ii
   close(iunit)

   if(lcs == "igbp")then
      lndname = trim(source_dir) // 'landtype_igbp_update.nc'
      !maxlc = 17
   else 
      lndname = trim(source_dir) // 'landtype_usgs_update.nc'
      !maxlc = 24
   end if

   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"landtype",varid(1)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(1),landtypes))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"landtypes",minval(landtypes),maxval(landtypes)
   
   dy = 180. / nlats_source
   dx = 360. / nlons_source

   allocate(lon_i(nlons_source))
   allocate(lat_i(nlats_source))

   do i = 1,nlons_source,1
      lon_i(i) = -180. + (2 * i - 1) * dx / 2.
   end do

   do i = 1,nlats_source,1
      lat_i(i) = 90. - (2 * i - 1) * dy / 2.
   end do

   allocate(seaorland(nlons_source, nlats_source))
   sum_sea = 0
   sum_land = 0
   seaorland = 0

   do i = 1,nlons_source,1
      do j = 1,nlats_source,1
         if(landtypes(i,j) /= 0.)then
            seaorland(i,j) = 1
            sum_land = sum_land + 1
         else
            sum_sea = sum_sea + 1
         end if
      end do
   end do

   lndname = trim(source_dir) // 'slope_max.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"slope_max",varid(1)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(1),slope_max))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"slope",minval(slope_max),maxval(slope_max)

   lndname = trim(source_dir) // 'slope_avg.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"slope_avg",varid(1)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(1),slope_avg))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"slope_avg",minval(slope_avg),maxval(slope_avg)

   allocate(lai(nlons_source,nlats_source))
   lai = 0.
   lndname = trim(source_dir) // 'lai.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"lai",varid(1)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(1),lai))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"lai",minval(lai),maxval(lai)

   allocate(k_s(nlons_source,nlats_source,2))
   allocate(k_sl(nlons_source,nlats_source,2))
   allocate(tkdry(nlons_source,nlats_source,2))
   allocate(tksatf(nlons_source,nlats_source,2))
   allocate(tksatu(nlons_source,nlats_source,2))
   k_s = 0.
   k_sl = 0.
   tkdry = 0.
   tksatf = 0.
   tksatu = 0.

   lndname = trim(source_dir) //'k_s.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"k_s_l1",varid(1)))
   CALL CHECK(NF90_INQ_VARID(ncid,"k_s_l2",varid(2)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(1),k_s(:,:,1)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(2),k_s(:,:,2)))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"k_s_l1",minval(k_s(:,:,1)),maxval(k_s(:,:,1))
   print*,"k_s_l2",minval(k_s(:,:,2)),maxval(k_s(:,:,2))
   
   lndname = trim(source_dir) // 'k_solids.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"k_solids_l1",varid(1)))
   CALL CHECK(NF90_INQ_VARID(ncid,"k_solids_l2",varid(2)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(1),k_sl(:,:,1)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(2),k_sl(:,:,2)))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"k_solids_l1",minval(k_sl(:,:,1)),maxval(k_sl(:,:,1))
   print*,"k_solids_l2",minval(k_sl(:,:,2)),maxval(k_sl(:,:,2))

   lndname = trim(source_dir) // 'tkdry.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"tkdry_l1",varid(1)))
   CALL CHECK(NF90_INQ_VARID(ncid,"tkdry_l2",varid(2)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(1),tkdry(:,:,1)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(2),tkdry(:,:,2)))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"tkdry_l1",minval(tkdry(:,:,1)),maxval(tkdry(:,:,1))
   print*,"tkdry_l2",minval(tkdry(:,:,2)),maxval(tkdry(:,:,2))

   lndname = trim(source_dir) // 'tksatf.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"tksatf_l1",varid(1)))
   CALL CHECK(NF90_INQ_VARID(ncid,"tksatf_l2",varid(2)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(1),tksatf(:,:,1)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(2),tksatf(:,:,2)))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"tksatf_l1",minval(tksatf(:,:,1)),maxval(tksatf(:,:,1))
   print*,"tksatf_l2",minval(tksatf(:,:,2)),maxval(tksatf(:,:,2))

   lndname = trim(source_dir) // 'tksatu.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"tksatu_l1",varid(1)))
   CALL CHECK(NF90_INQ_VARID(ncid,"tksatu_l2",varid(2)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(1),tksatu(:,:,1)))
   CALL CHECK(NF90_GET_VAR(ncid,varid(2),tksatu(:,:,2)))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"tksatu_l1",minval(tksatu(:,:,1)),maxval(tksatu(:,:,1))
   print*,"tksatu_l2",minval(tksatu(:,:,2)),maxval(tksatu(:,:,2))

   print*,"Start calculating the area of the latitude and longitude grid......"
   allocate(area_fine_gridcell(nlons_source,nlats_source))
   call cellarea(area_fine_gridcell)
   print*,"The latitude and longitude grid area is calculated"

   allocate(nla(0:maxlc))                           ! Whether the unstructured grid contains an indication of the land type
   allocate(area(0:maxlc))                         ! Area of each land type in unstructured grid
   allocate(fraction_mainarea(sjx_points,2))    ! Number and proportion of main land types in unstructured grid
   allocate(p_slope(sjx_points,3))            ! Maximum slope of latitude and longitude grid in unstructured grid
   allocate(p_lai(sjx_points,3))                ! Mean value, maximum value and number of lai in unstructured grids
   allocate(p_k_s(sjx_points,3,2))
   allocate(p_k_sl(sjx_points,3,2))
   allocate(p_tkdry(sjx_points,3,2))
   allocate(p_tksatf(sjx_points,3,2))
   allocate(p_tksatu(sjx_points,3,2))
   allocate(n_landtypes(sjx_points))
   allocate(ref_jw(nlons_source,nlats_source))
   allocate(ref_tr(sjx_points*20,16))
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
! 2.The threshold file for the uniterated grid is preliminarily calculated
!--------------------------------------------------------------------------

   print*,"开始计算阈值文件"
!$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,row,col,nla,L,area,maxid)
   do i = 1,sjx_points,1
      nla = 0
      area = 0.
      maxid = 0
      if(mp_id_new(i,1) == 0)then
         cycle
      end if
      do j = 0,mp_id_new(i,1)-1,1
         col = mp_ii(mp_id_new(i,2)+j,1)
         row = mp_ii(mp_id_new(i,2)+j,2)
         if(col == 0. .or. row == 0.)then
            cycle
         end if
         L = int(landtypes(col,row))

         if((L /= 0).and.(L /= maxlc))then

            nla(L) = 1
            area(L) = area(L) + mp_ii(mp_id_new(i,2)+j,4)

            p_lai(i,1) = p_lai(i,1) + lai(col,row)
            p_lai(i,3) = p_lai(i,3) + 1

            p_k_s(i,1,:) = p_k_s(i,1,:) + k_s(col,row,:)
            p_k_s(i,3,:) = p_k_s(i,3,:) + 1

            p_k_sl(i,1,:) = p_k_sl(i,1,:) + k_sl(col,row,:)
            p_k_sl(i,3,:) = p_k_sl(i,3,:) + 1

            p_tkdry(i,1,:) = p_tkdry(i,1,:) + tkdry(col,row,:)
            p_tkdry(i,3,:) = p_tkdry(i,3,:) + 1

            p_tksatf(i,1,:) = p_tksatf(i,1,:) + tksatf(col,row,:)
            p_tksatf(i,3,:) = p_tksatf(i,3,:) + 1

            p_tksatu(i,1,:) = p_tksatu(i,1,:) + tksatu(col,row,:)
            p_tksatu(i,3,:) = p_tksatu(i,3,:) + 1

            p_slope(i,1) = p_slope(i,1) + slope_avg(col,row)
            p_slope(i,3) = p_slope(i,3) + 1

         end if

      end do

      n_landtypes(i) = INT(sum(nla))
      maxid = maxloc(area) - 1
      fraction_mainarea(i,1) = maxid(1)
      fraction_mainarea(i,2) = area(maxid(1)) / mp_new(i,3)

      if(fraction_mainarea(i,2) > 1.)then
         fraction_mainarea(i,2) = 1.
      end if

      if(p_lai(i,3) /= 0.)then
         p_lai(i,1) = p_lai(i,1) / p_lai(i,3)
      end if

      if(p_k_s(i,3,1) /= 0.)then
         p_k_s(i,1,1) = p_k_s(i,1,1) / p_k_s(i,3,1)
      end if

      if(p_k_s(i,3,2) /= 0.)then
         p_k_s(i,1,2) = p_k_s(i,1,2) / p_k_s(i,3,2)
      end if

      if(p_k_sl(i,3,1) /= 0.)then
         p_k_sl(i,1,1) = p_k_sl(i,1,1) / p_k_sl(i,3,1)
      end if

      if(p_k_sl(i,3,2) /= 0.)then
         p_k_sl(i,1,2) = p_k_sl(i,1,2) / p_k_sl(i,3,2)
      end if

      if(p_tkdry(i,3,1) /= 0.)then
         p_tkdry(i,1,1) = p_tkdry(i,1,1) / p_tkdry(i,3,1)
      end if

      if(p_tkdry(i,3,2) /= 0.)then
         p_tkdry(i,1,2) = p_tkdry(i,1,2) / p_tkdry(i,3,2)
      end if

      if(p_tksatf(i,3,1) /= 0.)then
         p_tksatf(i,1,1) = p_tksatf(i,1,1) / p_tksatf(i,3,1)
      end if

      if(p_tksatf(i,3,2) /= 0.)then
         p_tksatf(i,1,2) = p_tksatf(i,1,2) / p_tksatf(i,3,2)
      end if

      if(p_tksatu(i,3,1) /= 0.)then
         p_tksatu(i,1,1) = p_tksatu(i,1,1) / p_tksatu(i,3,1)
      end if

      if(p_tksatu(i,3,2) /= 0.)then
         p_tksatu(i,1,2) = p_tksatu(i,1,2) / p_tksatu(i,3,2)
      end if

      if(p_slope(i,3) /= 0.)then
         p_slope(i,1) = p_slope(i,1) / p_slope(i,3)
      end if
 
   end do
!$OMP END PARALLEL DO


!$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,L,row,col)
   
   do i = 1,sjx_points,1
      if(mp_id_new(i,1) == 0)then
         cycle
      end if
      do j = 0,mp_id_new(i,1)-1,1
         col = mp_ii(mp_id_new(i,2)+j,1)
         row = mp_ii(mp_id_new(i,2)+j,2)
         if(col == 0. .or. row == 0.)then
            cycle
         end if
         L = int(landtypes(col,row))

         if((L /= 0).and.(L /= maxlc))then

            p_lai(i,2)      = p_lai(i,2)      + (lai(col,row)-p_lai(i,1))*(lai(col,row)-p_lai(i,1))
            p_k_s(i,2,1)    = p_k_s(i,2,1)    + (k_s(col,row,1)-p_k_s(i,1,1))*(k_s(col,row,1)-p_k_s(i,1,1))
            p_k_s(i,2,2)    = p_k_s(i,2,2)    + (k_s(col,row,2)-p_k_s(i,1,2))*(k_s(col,row,2)-p_k_s(i,1,2))
            p_k_sl(i,2,1)   = p_k_sl(i,2,1)   + (k_sl(col,row,1)-p_k_sl(i,1,1))*(k_sl(col,row,1)-p_k_sl(i,1,1))
            p_k_sl(i,2,2)   = p_k_sl(i,2,2)   + (k_sl(col,row,2)-p_k_sl(i,1,2))*(k_sl(col,row,2)-p_k_sl(i,1,2))
            p_tkdry(i,2,1)  = p_tkdry(i,2,1)  + (tkdry(col,row,1)-p_tkdry(i,1,1))*(tkdry(col,row,1)-p_tkdry(i,1,1))
            p_tkdry(i,2,2)  = p_tkdry(i,2,2)  + (tkdry(col,row,2)-p_tkdry(i,1,2))*(tkdry(col,row,2)-p_tkdry(i,1,2))
            p_tksatf(i,2,1) = p_tksatf(i,2,1) + (tksatf(col,row,1)-p_tksatf(i,1,1))*(tksatf(col,row,1)-p_tksatf(i,1,1))
            p_tksatf(i,2,2) = p_tksatf(i,2,2) + (tksatf(col,row,2)-p_tksatf(i,1,2))*(tksatf(col,row,2)-p_tksatf(i,1,2))
            p_tksatu(i,2,1) = p_tksatu(i,2,1) + (tksatu(col,row,1)-p_tksatu(i,1,1))*(tksatu(col,row,1)-p_tksatu(i,1,1))
            p_tksatu(i,2,2) = p_tksatu(i,2,2) + (tksatu(col,row,2)-p_tksatu(i,1,2))*(tksatu(col,row,2)-p_tksatu(i,1,2))
            p_slope(i,2)    = p_slope(i,2)    + (slope_avg(col,row)-p_slope(i,1))*(slope_avg(col,row)-p_slope(i,1))

         end if
      end do

      p_lai(i,2)      = sqrt(p_lai(i,2)      / p_lai(i,3))
      p_k_s(i,2,1)    = sqrt(p_k_s(i,2,1)    / p_k_s(i,3,1))
      p_k_s(i,2,2)    = sqrt(p_k_s(i,2,2)    / p_k_s(i,3,2))
      p_k_sl(i,2,1)   = sqrt(p_k_sl(i,2,1)   / p_k_sl(i,3,1))
      p_k_sl(i,2,2)   = sqrt(p_k_sl(i,2,2)   / p_k_sl(i,3,2))
      p_tkdry(i,2,1)  = sqrt(p_tkdry(i,2,1)  / p_tkdry(i,3,1))
      p_tkdry(i,2,2)  = sqrt(p_tkdry(i,2,2)  / p_tkdry(i,3,2))
      p_tksatf(i,2,1) = sqrt(p_tksatf(i,2,1) / p_tksatf(i,3,1))
      p_tksatf(i,2,2) = sqrt(p_tksatf(i,2,2) / p_tksatf(i,3,2))
      p_tksatu(i,2,1) = sqrt(p_tksatu(i,2,1) / p_tksatu(i,3,1))
      p_tksatu(i,2,2) = sqrt(p_tksatu(i,2,2) / p_tksatu(i,3,2))
      p_slope(i,2)    = sqrt(p_slope(i,2)    / p_slope(i,3))

   end do
!$OMP END PARALLEL DO

   print*,"The threshold file is calculated"

!--------------------------------------------------------------------------
! 3.Stores threshold files and grid information for uniterated grids
!--------------------------------------------------------------------------

   print*,"Stores the unrefined grid threshold file"
   iter = 0
   write(it,'(i2.2)')iter
   lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/threshold/threshold_NXP" // trim(nxpc) // "_" // trim(adjustl(it)) // ".nc4"
   print*,lndname
   CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber,nf90_netcdf4), ncID))
   CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", sjx_points, spDimID))
   CALL CHECK(NF90_DEF_DIM(ncID, "dima"        , 2       , nmDimID))
   CALL CHECK(NF90_DEF_VAR(ncID, "num_landtypes"    , NF90_INT  , (/ spDimID /)         , VarID(1)))
   CALL CHECK(NF90_DEF_VAR(ncID, "fraction_mainarea", NF90_float, (/ spDimID, nmDimID /), VarID(2)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_slope"          , NF90_float, (/ spDimID, nmDimID /), VarID(3)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_lai"            , NF90_float, (/ spDimID, nmDimID /), VarID(4)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_k_s"            , NF90_float, (/ spDimID, nmDimID, nmDimID /), VarID(5)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_k_sl"           , NF90_float, (/ spDimID, nmDimID, nmDimID /), VarID(6)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_tkdry"          , NF90_float, (/ spDimID, nmDimID, nmDimID /), VarID(7)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_tksatf"         , NF90_float, (/ spDimID, nmDimID, nmDimID /), VarID(8)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_tksatu"         , NF90_float, (/ spDimID, nmDimID, nmDimID /), VarID(9)))
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

   print*,"n_landtypes",minval(n_landtypes),maxval(n_landtypes)
   print*,"mainarea_id",minval(fraction_mainarea(:,1)),maxval(fraction_mainarea(:,1))
   print*,"fraction_mainarea",minval(fraction_mainarea(:,2)),maxval(fraction_mainarea(:,2))
   print*,"Mean:"
   print*,"slope",minval(p_slope(:,1)),maxval(p_slope(:,1))
   print*,"lai",minval(p_lai(:,1)),maxval(p_lai(:,1))
   print*,"k_s_l1",minval(p_k_s(:,1,1)),maxval(p_k_s(:,1,1))
   print*,"k_s_l2",minval(p_k_s(:,1,2)),maxval(p_k_s(:,1,2))
   print*,"k_solids_l1",minval(p_k_sl(:,1,1)),maxval(p_k_sl(:,1,1))
   print*,"k_solids_l2",minval(p_k_sl(:,1,2)),maxval(p_k_sl(:,1,2))
   print*,"tkdry_l1",minval(p_tkdry(:,1,1)),maxval(p_tkdry(:,1,1))
   print*,"tkdry_l2",minval(p_tkdry(:,1,2)),maxval(p_tkdry(:,1,2))
   print*,"tksatf_l1",minval(p_tksatf(:,1,1)),maxval(p_tksatf(:,1,1))
   print*,"tksatf_l2",minval(p_tksatf(:,1,2)),maxval(p_tksatf(:,1,2))
   print*,"tksatu_l1",minval(p_tksatu(:,1,1)),maxval(p_tksatu(:,1,1))
   print*,"tksatu_l2",minval(p_tksatu(:,1,2)),maxval(p_tksatu(:,1,2))
   print*,"Standard Deviation:"
   print*,"slope",minval(p_slope(:,2)),maxval(p_slope(:,2))
   print*,"lai",minval(p_lai(:,2)),maxval(p_lai(:,2))
   print*,"k_s_l1",minval(p_k_s(:,2,1)),maxval(p_k_s(:,2,1))
   print*,"k_s_l2",minval(p_k_s(:,2,2)),maxval(p_k_s(:,2,2))
   print*,"k_solids_l1",minval(p_k_sl(:,2,1)),maxval(p_k_sl(:,2,1))
   print*,"k_solids_l2",minval(p_k_sl(:,2,2)),maxval(p_k_sl(:,2,2))
   print*,"tkdry_l1",minval(p_tkdry(:,2,1)),maxval(p_tkdry(:,2,1))
   print*,"tkdry_l2",minval(p_tkdry(:,2,2)),maxval(p_tkdry(:,2,2))
   print*,"tksatf_l1",minval(p_tksatf(:,2,1)),maxval(p_tksatf(:,2,1))
   print*,"tksatf_l2",minval(p_tksatf(:,2,2)),maxval(p_tksatf(:,2,2))
   print*,"tksatu_l1",minval(p_tksatu(:,2,1)),maxval(p_tksatu(:,2,1))
   print*,"tksatu_l2",minval(p_tksatu(:,2,2)),maxval(p_tksatu(:,2,2))

   lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/gridfile_NXP" // trim(nxpc) // "_" // trim(adjustl(it)) // ".nc4"
   print*,lndname
   CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber,nf90_netcdf4), ncID))
   CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", sjx_points , spDimID))
   CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", lbx_points , lpDimID))
   CALL CHECK(NF90_DEF_DIM(ncID, "dim_a"     , 3          , thDimID))
   CALL CHECK(NF90_DEF_DIM(ncID, "dim_c"     , 7          , seDimID))
   CALL CHECK(NF90_DEF_VAR(ncID, "GLONM"     , NF90_FLOAT, (/ spDimID /), varid(1)))
   CALL CHECK(NF90_DEF_VAR(ncID, "GLATM"     , NF90_FLOAT, (/ spDimID /), varid(2)))
   CALL CHECK(NF90_DEF_VAR(ncID, "GLONW"     , NF90_FLOAT, (/ lpDimID /), varid(3)))
   CALL CHECK(NF90_DEF_VAR(ncID, "GLATW"     , NF90_FLOAT, (/ lpDimID /), varid(4)))
   CALL CHECK(NF90_DEF_VAR(ncID, "itab_m%iw" , NF90_INT  , (/ thDimID, spDimID /), varid(5)))
   CALL CHECK(NF90_DEF_VAR(ncID, "itab_w%im" , NF90_INT  , (/ seDimID, lpDimID /), varid(6)))
   CALL CHECK(NF90_ENDDEF(ncID))

   CALL CHECK(NF90_PUT_VAR(ncID, varid(1), mp_new(1:sjx_points,1)))
   CALL CHECK(NF90_PUT_VAR(ncID, varid(2), mp_new(1:sjx_points,2)))
   CALL CHECK(NF90_PUT_VAR(ncID, varid(3), wp_new(1:lbx_points,1)))
   CALL CHECK(NF90_PUT_VAR(ncID, varid(4), wp_new(1:lbx_points,2)))
   CALL CHECK(NF90_PUT_VAR(ncID, varid(5), ngrmw_new(1:3,1:sjx_points)))
   CALL CHECK(NF90_PUT_VAR(ncID, varid(6), ngrwm_new(1:7,1:lbx_points)))
   CALL CHECK(NF90_CLOSE(ncID))

   ref = 0

!--------------------------------------------------------------------------
! 4.The initial iteration grid is selected according to the threshold options
!--------------------------------------------------------------------------
   allocate(IsInRfArea(sjx_points))
   IsInRfArea = .false.
   CALL IsInRefineArea(IsInRfArea,sjx_points,lbx_points,ngrmw_new,wp_new,sjx_points)

   do i = 1,sjx_points,1
      if(IsInRfArea(i) == .false.)then
         cycle
      end if

      if((fraction_mainarea(i,1) == 0.).or.(fraction_mainarea(i,1) == maxlc))then
         cycle
      end if

         if (refine_num_landtypes == .True. .and. n_landtypes(i)>th_num_landtypes) then

            ref(i) = 1
            ref_tr(i,1) = 1
         end if

         if (refine_area_mainland == .True. .and. fraction_mainarea(i, 2) < th_area_mainland) then
            ref(i) = 1
            ref_tr(i,2) = 1
         end if
         if (refine_lai_m == .True. .and. (p_lai(i, 1) > th_lai_m)) then
            ref(i) = 1
            ref_tr(i,3) = 1
         end if

         if (refine_lai_s == .True. .and. (p_lai(i, 2) > th_lai_s)) then
            ref(i) = 1
            ref_tr(i,4) = 1
         end if
         
         if (refine_slope_m == .True. .and. (p_slope(i, 1) > th_slope_m)) then

            ref(i) = 1
            ref_tr(i,5) = 1
         end if

         if (refine_slope_s == .True. .and. (p_slope(i, 2) > th_slope_s)) then

            ref(i) = 1
            ref_tr(i,6) = 1
         end if

         if (refine_k_s_m == .True. .and. ((p_k_s(i, 1, 1) > th_k_s_m).or.(p_k_s(i, 1, 2) > th_k_s_m))) then
            ref(i) = 1
            ref_tr(i,7) = 1
         end if

         if (refine_k_s_s == .True. .and. ((p_k_s(i, 2, 1) > th_k_s_s).or.(p_k_s(i, 2, 2) > th_k_s_s))) then
            ref(i) = 1
            ref_tr(i,8) = 1
         end if

         if (refine_k_solids_m == .True. .and. ((p_k_sl(i, 1, 1) > th_k_solids_m).or.(p_k_sl(i, 1, 2) > th_k_solids_m))) then
            ref(i) = 1
            ref_tr(i,9) = 1
         end if
         
         if (refine_k_solids_s == .True. .and. ((p_k_sl(i, 2, 1) > th_k_solids_s).or.(p_k_sl(i, 2, 2) > th_k_solids_s))) then
            ref(i) = 1
            ref_tr(i,10) = 1
         end if

         if (refine_tkdry_m == .True. .and. ((p_tkdry(i, 1, 1) > th_tkdry_m).or.(p_tkdry(i, 1, 2) > th_tkdry_m))) then
            ref(i) = 1
            ref_tr(i,11) = 1
         end if

         if (refine_tkdry_s == .True. .and. ((p_tkdry(i, 2, 1) > th_tkdry_s).or.(p_tkdry(i, 2, 2) > th_tkdry_s))) then
            ref(i) = 1
            ref_tr(i,12) = 1
         end if

         if (refine_tksatf_m == .True. .and. ((p_tksatf(i, 1, 1) > th_tksatf_m).or.(p_tksatf(i, 1, 2) > th_tksatf_m))) then
            ref(i) = 1
            ref_tr(i,13) = 1
         end if

         if (refine_tksatf_s == .True. .and. ((p_tksatf(i, 2, 1) > th_tksatf_s).or.(p_tksatf(i, 2, 2) > th_tksatf_s))) then
            ref(i) = 1
            ref_tr(i,14) = 1
         end if
         
         if (refine_tksatu_m == .True. .and. ((p_tksatu(i, 1, 1) > th_tksatu_m).or.(p_tksatu(i, 1, 2) > th_tksatu_m))) then
            ref(i) = 1
            ref_tr(i,15) = 1
         end if

         if (refine_tksatu_s == .True. .and. ((p_tksatu(i, 2, 1) > th_tksatu_s).or.(p_tksatu(i, 2, 2) > th_tksatu_s))) then

            ref(i) = 1
            ref_tr(i,16) = 1
         end if

   end do

   num_ref = INT(sum(ref))

   if(num_ref == 0)then
      isover = .true.
   end if

   refed = 0

!--------------------------------------------------------------------------
! 5.Start iteration, iteration process is:
!  5.1 Determine whether the iteration continues according to flag
!  5.2. Refinement
!  5.3. Record the refined grid file
!  5.4. Calculate the mesh inclusion relationship after refinement
!  5.5 Calculate and record the refined threshold file
!  5.6. Calculate the number of meshes to be refined according to the threshold
!  5.7. Set the flag according to the number of grids
!--------------------------------------------------------------------------

!------------------------------------
! 5.1 Determine whether the iteration continues according to the flag
!------------------------------------
   do while(isover == .false.)
      iter = iter + 1
      nmp(iter) = nmp(iter-1) + num_ref * 4
      nwp(iter) = nwp(iter-1) + num_ref * 3

!------------------------------------
! 5.2 refine
!------------------------------------
      do i = nmp(iter-2)+1,nmp(iter-1),1
         if(ref(i) == 1)then
            icl = 0
            sjx = 0.
            sjx(1,1:2) = wp_new(ngrmw_new(1,i),1:2)
            sjx(2,1:2) = wp_new(ngrmw_new(2,i),1:2)
            sjx(3,1:2) = wp_new(ngrmw_new(3,i),1:2)

            icl(1) = IsCrossLine(sjx(2,1),sjx(3,1))
            icl(2) = IsCrossLine(sjx(1,1),sjx(3,1))
            icl(3) = IsCrossLine(sjx(1,1),sjx(2,1))

            if(sum(icl) > 0.)then
               do j = 1,3,1
                  if(sjx(j,1) < 0.)then
                     sjx(j,1) = sjx(j,1) + 360.
                  end if
               end do
            end if

            newsjx(1,1:2) = (sjx(2,1:2) + sjx(3,1:2)) / 2.
            newsjx(2,1:2) = (sjx(1,1:2) + sjx(3,1:2)) / 2.
            newsjx(3,1:2) = (sjx(1,1:2) + sjx(2,1:2)) / 2.

            m1 = nmp(iter-1)+refed(iter)*4+1
            m2 = nmp(iter-1)+refed(iter)*4+2
            m3 = nmp(iter-1)+refed(iter)*4+3
            m4 = nmp(iter-1)+refed(iter)*4+4

            w1 = nwp(iter-1)+refed(iter)*3+1
            w2 = nwp(iter-1)+refed(iter)*3+2
            w3 = nwp(iter-1)+refed(iter)*3+3

            mp_new(m1,1:2) = (sjx(1,1:2) + newsjx(2,1:2) + newsjx(3,1:2)) / 3.
            mp_new(m2,1:2) = (sjx(2,1:2) + newsjx(1,1:2) + newsjx(3,1:2)) / 3.
            mp_new(m3,1:2) = (sjx(3,1:2) + newsjx(1,1:2) + newsjx(2,1:2)) / 3.
            mp_new(m4,1:2) = (sjx(1,1:2) + sjx(2,1:2) + sjx(3,1:2)) / 3.

            wp_new(w1,1:2) = newsjx(1,1:2)
            wp_new(w2,1:2) = newsjx(2,1:2)
            wp_new(w3,1:2) = newsjx(3,1:2)

            ngrmw_new(1,m1) = ngrmw_new(1,i) 
            ngrmw_new(2,m1) = w3
            ngrmw_new(3,m1) = w2
            ngrmw_new(1,m2) = ngrmw_new(2,i)
            ngrmw_new(2,m2) = w1
            ngrmw_new(3,m2) = w3
            ngrmw_new(1,m3) = ngrmw_new(3,i)
            ngrmw_new(2,m3) = w2
            ngrmw_new(3,m3) = w1
            ngrmw_new(1,m4) = w1
            ngrmw_new(2,m4) = w2
            ngrmw_new(3,m4) = w3

            mp_new(i,4) = -1
            mp_id_new(i,1) = 0
            mp_new(i,3) = 0.
            ngrmw_new(:,i) = 0

            CALL CheckLon(mp_new(m1,1))
            CALL CheckLon(mp_new(m2,1))
            CALL CheckLon(mp_new(m3,1))
            CALL CheckLon(mp_new(m4,1))
            CALL CheckLon(wp_new(w1,1))
            CALL CheckLon(wp_new(w2,1))
            CALL CheckLon(wp_new(w3,1))

            do j = 1,16,1
               if(ref_tr(i,j) == 1)then
                  ref_tr(m1,j) = 1
                  ref_tr(m2,j) = 1
                  ref_tr(m3,j) = 1
                  ref_tr(m4,j) = 1
               end if
            end do

            refed(iter) = refed(iter) + 1
         end if
      end do

      print*,"The number of refined triangles is",refed(iter)

!------------------------------------
! 5.3 Record refinement after grid file
!------------------------------------
      write(it,'(i2.2)')iter
      lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/tmp/gridfile_NXP" // trim(nxpc) // "_" // trim(adjustl(it)) // ".nc4"
      print*,lndname
      CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber,nf90_netcdf4), ncID))
      CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", nmp(iter)  , spDimID))
      CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", nwp(iter)  , lpDimID))
      CALL CHECK(NF90_DEF_DIM(ncID, "dim_a"     , 3          , thDimID))
      CALL CHECK(NF90_DEF_DIM(ncID, "dim_d"     , 16         , sxDimID))
      CALL CHECK(NF90_DEF_DIM(ncID, "dim_c"     , 7          , seDimID))
      CALL CHECK(NF90_DEF_VAR(ncID, "GLONM"     , NF90_FLOAT, (/ spDimID /), varid(1)))
      CALL CHECK(NF90_DEF_VAR(ncID, "GLATM"     , NF90_FLOAT, (/ spDimID /), varid(2)))
      CALL CHECK(NF90_DEF_VAR(ncID, "GLONW"     , NF90_FLOAT, (/ lpDimID /), varid(3)))
      CALL CHECK(NF90_DEF_VAR(ncID, "GLATW"     , NF90_FLOAT, (/ lpDimID /), varid(4)))
      CALL CHECK(NF90_DEF_VAR(ncID, "itab_m%iw" , NF90_INT  , (/ thDimID, spDimID /), varid(5)))
      CALL CHECK(NF90_DEF_VAR(ncID, "itab_w%im" , NF90_INT  , (/ seDimID, lpDimID /), varid(6)))
      CALL CHECK(NF90_DEF_VAR(ncID, "ref_tr"    , NF90_INT  , (/ spDimID, sxDimID /), varid(7)))
      CALL CHECK(NF90_ENDDEF(ncID))

      CALL CHECK(NF90_PUT_VAR(ncID, varid(1), mp_new(1:nmp(iter),1)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(2), mp_new(1:nmp(iter),2)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(3), wp_new(1:nwp(iter),1)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(4), wp_new(1:nwp(iter),2)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(5), ngrmw_new(1:3,1:nmp(iter))))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(6), ngrwm_new(1:7,1:nwp(iter))))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(7), ref_tr(1:nmp(iter),1:16)))
      CALL CHECK(NF90_CLOSE(ncID))

      !stop
      ref_jw = 0
      do j = nmp(iter-1)+1,nmp(iter),1
         tmpa(1:2) = wp_new(ngrmw_new(1,j),1:2)
         tmpb(1:2) = wp_new(ngrmw_new(2,j),1:2)
         tmpc(1:2) = wp_new(ngrmw_new(3,j),1:2)
         CALL GetRefineJW(tmpa, tmpb, tmpc, ref_jw)
         !CALL GetRefineJW(wp_new(ngrmw_new(1,j),1:2),wp_new(ngrmw_new(2,j),1:2),wp_new(ngrmw_new(3,j),1:2),ref_jw)
      end do

      num_refjw = INT(sum(ref_jw))
      print*,"The number of latitude and longitude grids to be adjusted is",num_refjw

!------------------------------------
! 5.4 The mesh inclusion relationship is calculated after refinement
!------------------------------------

      mp_i = 0.

!$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1)&
!$OMP PRIVATE(i,j,k,ispart,sjx,isinply)
      do i = 1,nlons_source,1
         do j = 1,nlats_source,1
            if((seaorland(i,j) == 0).or.(ref_jw(i,j) == 0))then
               cycle
            end if
            ispart = .false.
            do k = nmp(iter-1)+1,nmp(iter),1
               if(ispart == .true.)then
                  exit
               end if
               sjx = 0.
               sjx(1,1:2) = wp_new(ngrmw_new(1,k),1:2)
               sjx(2,1:2) = wp_new(ngrmw_new(2,k),1:2)
               sjx(3,1:2) = wp_new(ngrmw_new(3,k),1:2)
               isinply = IsInUstrGrid(sjx,lon_i(i),lat_i(j))

               if(isinply == 1)then
                  mp_new(k,3) = mp_new(k,3) + area_fine_gridcell(i,j)
                  mp_id_new(k,1) = mp_id_new(k,1) + 1
                  mp_i(k,mp_id_new(k,1),1) = i
                  mp_i(k,mp_id_new(k,1),2) = j
                  mp_i(k,mp_id_new(k,1),3) = isinply
                  mp_i(k,mp_id_new(k,1),4) = area_fine_gridcell(i,j)
                  ispart = .true.
                  !print*,i,j,k,isinply
               else if(isinply > 0)then
                  mp_new(k,3) = mp_new(k,3) + area_fine_gridcell(i,j)*isinply
                  mp_id_new(k,1) = mp_id_new(k,1) + 1
                  mp_i(k,mp_id_new(k,1),1) = i
                  mp_i(k,mp_id_new(k,1),2) = j
                  mp_i(k,mp_id_new(k,1),3) = isinply
                  mp_i(k,mp_id_new(k,1),4) = area_fine_gridcell(i,j)*isinply
                  !print*,i,j,k,isinply
               end if

            end do
         end do
      end do
!$OMP END PARALLEL DO

      print*,sum(mp_id_new(nmp(iter-1) + 1:nmp(iter),1))
      do i = nmp(iter-1) + 1,nmp(iter),1
         mp_id_new(i,2) = mp_id_new(i-1,2) + mp_id_new(i-1,1)
      end do

      do i = nmp(iter-1) + 1,nmp(iter),1
         if(mp_id_new(i,1) /= 0.)then
            do j = 1,4,1
               mp_ii_new(mp_id_new(i,2):(mp_id_new(i,2)+mp_id_new(i,1)-1),j) = mp_i(i,1:mp_id_new(i,1),j)
            end do
         end if
      end do

!------------------------------------
! 5.5 Calculate and record the refined threshold file
!------------------------------------

      deallocate(fraction_mainarea)
      deallocate(n_landtypes)
      deallocate(p_slope)
      deallocate(p_lai)
      deallocate(p_k_s)
      deallocate(p_k_sl)
      deallocate(p_tkdry)
      deallocate(p_tksatf)
      deallocate(p_tksatu)

      allocate(fraction_mainarea(nmp(iter),2))   ! Number and proportion of main land types in unstructured grid
      allocate(n_landtypes(nmp(iter)))
      allocate(p_slope(nmp(iter),3))
      allocate(p_lai(nmp(iter),3))
      allocate(p_k_s(nmp(iter),3,2))
      allocate(p_k_sl(nmp(iter),3,2))
      allocate(p_tkdry(nmp(iter),3,2))
      allocate(p_tksatf(nmp(iter),3,2))
      allocate(p_tksatu(nmp(iter),3,2))
      fraction_mainarea = 0.
      n_landtypes = 0
      p_slope = 0.
      p_lai = 0.
      p_k_s = 0.
      p_k_sl = 0.
      p_tkdry = 0.
      p_tksatf = 0.
      p_tksatu = 0.
      
!$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,row,col,nla,L,area,maxid)
      do i = 1,sjx_points,1
         nla = 0
         area = 0.
         maxid = 0
         if(mp_id_new(i,1) == 0)then
            cycle
         end if
         do j = 0,mp_id_new(i,1)-1,1
            col = mp_ii(mp_id_new(i,2)+j,1)
            row = mp_ii(mp_id_new(i,2)+j,2)
            if((col == 0.).or.(row == 0.))then
               cycle
            end if
            L = int(landtypes(col,row))

            if((L /= 0).and.(L /= maxlc))then

               nla(L) = 1
               area(L) = area(L) + mp_ii(mp_id_new(i,2)+j,4)

               p_lai(i,1) = p_lai(i,1) + lai(col,row)
               p_lai(i,3) = p_lai(i,3) + 1

               p_k_s(i,1,:) = p_k_s(i,1,:) + k_s(col,row,:)
               p_k_s(i,3,:) = p_k_s(i,3,:) + 1

               p_k_sl(i,1,:) = p_k_sl(i,1,:) + k_sl(col,row,:)
               p_k_sl(i,3,:) = p_k_sl(i,3,:) + 1

               p_tkdry(i,1,:) = p_tkdry(i,1,:) + tkdry(col,row,:)
               p_tkdry(i,3,:) = p_tkdry(i,3,:) + 1

               p_tksatf(i,1,:) = p_tksatf(i,1,:) + tksatf(col,row,:)
               p_tksatf(i,3,:) = p_tksatf(i,3,:) + 1

               p_tksatu(i,1,:) = p_tksatu(i,1,:) + tksatu(col,row,:)
               p_tksatu(i,3,:) = p_tksatu(i,3,:) + 1

               p_slope(i,1) = p_slope(i,1) + slope_avg(col,row)
               p_slope(i,3) = p_slope(i,3) + 1

            end if
         end do

         n_landtypes(i) = INT(sum(nla))
         maxid = maxloc(area) - 1
         fraction_mainarea(i,1) = maxid(1)
         fraction_mainarea(i,2) = area(maxid(1)) / mp_new(i,3)

         if(fraction_mainarea(i,2) > 1.)then
            fraction_mainarea(i,2) = 1.
         end if

         if(p_lai(i,3) /= 0.)then
            p_lai(i,1) = p_lai(i,1) / p_lai(i,3)
         end if

         if(p_k_s(i,3,1) /= 0.)then
            p_k_s(i,1,1) = p_k_s(i,1,1) / p_k_s(i,3,1)
         end if

         if(p_k_s(i,3,2) /= 0.)then
            p_k_s(i,1,2) = p_k_s(i,1,2) / p_k_s(i,3,2)
         end if

         if(p_k_sl(i,3,1) /= 0.)then
            p_k_sl(i,1,1) = p_k_sl(i,1,1) / p_k_sl(i,3,1)
         end if

         if(p_k_sl(i,3,2) /= 0.)then
            p_k_sl(i,1,2) = p_k_sl(i,1,2) / p_k_sl(i,3,2)
         end if

         if(p_tkdry(i,3,1) /= 0.)then
            p_tkdry(i,1,1) = p_tkdry(i,1,1) / p_tkdry(i,3,1)
         end if

         if(p_tkdry(i,3,2) /= 0.)then
            p_tkdry(i,1,2) = p_tkdry(i,1,2) / p_tkdry(i,3,2)
         end if

         if(p_tksatf(i,3,1) /= 0.)then
            p_tksatf(i,1,1) = p_tksatf(i,1,1) / p_tksatf(i,3,1)
         end if

         if(p_tksatf(i,3,2) /= 0.)then
            p_tksatf(i,1,2) = p_tksatf(i,1,2) / p_tksatf(i,3,2)
         end if

         if(p_tksatu(i,3,1) /= 0.)then
            p_tksatu(i,1,1) = p_tksatu(i,1,1) / p_tksatu(i,3,1)
         end if

         if(p_tksatu(i,3,2) /= 0.)then
            p_tksatu(i,1,2) = p_tksatu(i,1,2) / p_tksatu(i,3,2)
         end if

         if(p_slope(i,3) /= 0.)then
            p_slope(i,1) = p_slope(i,1) / p_slope(i,3)
         end if

      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,row,col,nla,L,area,maxid)
      do i = sjx_points+1,nmp(iter),1
         nla = 0
         area = 0.
         maxid = 0
         if(mp_id_new(i,1) == 0)then
            cycle
         end if
         do j = 0,mp_id_new(i,1)-1,1
            col = mp_ii_new(mp_id_new(i,2)+j,1)
            row = mp_ii_new(mp_id_new(i,2)+j,2)
            if((col == 0.).or.(row == 0.))then
               cycle
            end if
            L = int(landtypes(col,row))

            if((L /= 0).and.(L /= maxlc))then

               nla(L) = 1
               area(L) = area(L) + mp_ii_new(mp_id_new(i,2)+j,4)

               p_lai(i,1) = p_lai(i,1) + lai(col,row)
               p_lai(i,3) = p_lai(i,3) + 1

               p_k_s(i,1,:) = p_k_s(i,1,:) + k_s(col,row,:)
               p_k_s(i,3,:) = p_k_s(i,3,:) + 1

               p_k_sl(i,1,:) = p_k_sl(i,1,:) + k_sl(col,row,:)
               p_k_sl(i,3,:) = p_k_sl(i,3,:) + 1

               p_tkdry(i,1,:) = p_tkdry(i,1,:) + tkdry(col,row,:)
               p_tkdry(i,3,:) = p_tkdry(i,3,:) + 1

               p_tksatf(i,1,:) = p_tksatf(i,1,:) + tksatf(col,row,:)
               p_tksatf(i,3,:) = p_tksatf(i,3,:) + 1

               p_tksatu(i,1,:) = p_tksatu(i,1,:) + tksatu(col,row,:)
               p_tksatu(i,3,:) = p_tksatu(i,3,:) + 1

               p_slope(i,1) = p_slope(i,1) + slope_avg(col,row)
               p_slope(i,3) = p_slope(i,3) + 1

            end if
         end do

         n_landtypes(i) = INT(sum(nla))
         maxid = maxloc(area) - 1
         fraction_mainarea(i,1) = maxid(1)
         fraction_mainarea(i,2) = area(maxid(1)) / mp_new(i,3)

         if(fraction_mainarea(i,2) > 1.)then
            fraction_mainarea(i,2) = 1.
         end if

         if(p_lai(i,3) /= 0.)then
            p_lai(i,1) = p_lai(i,1) / p_lai(i,3)
         end if

         if(p_k_s(i,3,1) /= 0.)then
            p_k_s(i,1,1) = p_k_s(i,1,1) / p_k_s(i,3,1)
         end if

         if(p_k_s(i,3,2) /= 0.)then
            p_k_s(i,1,2) = p_k_s(i,1,2) / p_k_s(i,3,2)
         end if

         if(p_k_sl(i,3,1) /= 0.)then
            p_k_sl(i,1,1) = p_k_sl(i,1,1) / p_k_sl(i,3,1)
         end if

         if(p_k_sl(i,3,2) /= 0.)then
            p_k_sl(i,1,2) = p_k_sl(i,1,2) / p_k_sl(i,3,2)
         end if

         if(p_tkdry(i,3,1) /= 0.)then
            p_tkdry(i,1,1) = p_tkdry(i,1,1) / p_tkdry(i,3,1)
         end if

         if(p_tkdry(i,3,2) /= 0.)then
            p_tkdry(i,1,2) = p_tkdry(i,1,2) / p_tkdry(i,3,2)
         end if

         if(p_tksatf(i,3,1) /= 0.)then
            p_tksatf(i,1,1) = p_tksatf(i,1,1) / p_tksatf(i,3,1)
         end if

         if(p_tksatf(i,3,2) /= 0.)then
            p_tksatf(i,1,2) = p_tksatf(i,1,2) / p_tksatf(i,3,2)
         end if

         if(p_tksatu(i,3,1) /= 0.)then
            p_tksatu(i,1,1) = p_tksatu(i,1,1) / p_tksatu(i,3,1)
         end if

         if(p_tksatu(i,3,2) /= 0.)then
            p_tksatu(i,1,2) = p_tksatu(i,1,2) / p_tksatu(i,3,2)
         end if

         if(p_slope(i,3) /= 0.)then
            p_slope(i,1) = p_slope(i,1) / p_slope(i,3)
         end if

      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,L,row,col)
      do i = 1,sjx_points,1
         if(mp_id_new(i,1) == 0)then
            cycle
         end if
         do j = 0,mp_id_new(i,1)-1,1
            col = mp_ii(mp_id_new(i,2)+j,1)
            row = mp_ii(mp_id_new(i,2)+j,2)
            if(col == 0. .or. row == 0.)then
               cycle
            end if

            L = int(landtypes(col,row))

            if((L /= 0).and.(L /= maxlc))then

               p_lai(i,2)      = p_lai(i,2)      + (lai(col,row)-p_lai(i,1))*(lai(col,row)-p_lai(i,1))
               p_k_s(i,2,1)    = p_k_s(i,2,1)    + (k_s(col,row,1)-p_k_s(i,1,1))*(k_s(col,row,1)-p_k_s(i,1,1))
               p_k_s(i,2,2)    = p_k_s(i,2,2)    + (k_s(col,row,2)-p_k_s(i,1,2))*(k_s(col,row,2)-p_k_s(i,1,2))
               p_k_sl(i,2,1)   = p_k_sl(i,2,1)   + (k_sl(col,row,1)-p_k_sl(i,1,1))*(k_sl(col,row,1)-p_k_sl(i,1,1))
               p_k_sl(i,2,2)   = p_k_sl(i,2,2)   + (k_sl(col,row,2)-p_k_sl(i,1,2))*(k_sl(col,row,2)-p_k_sl(i,1,2))
               p_tkdry(i,2,1)  = p_tkdry(i,2,1)  + (tkdry(col,row,1)-p_tkdry(i,1,1))*(tkdry(col,row,1)-p_tkdry(i,1,1))
               p_tkdry(i,2,2)  = p_tkdry(i,2,2)  + (tkdry(col,row,2)-p_tkdry(i,1,2))*(tkdry(col,row,2)-p_tkdry(i,1,2))
               p_tksatf(i,2,1) = p_tksatf(i,2,1) + (tksatf(col,row,1)-p_tksatf(i,1,1))*(tksatf(col,row,1)-p_tksatf(i,1,1))
               p_tksatf(i,2,2) = p_tksatf(i,2,2) + (tksatf(col,row,2)-p_tksatf(i,1,2))*(tksatf(col,row,2)-p_tksatf(i,1,2))
               p_tksatu(i,2,1) = p_tksatu(i,2,1) + (tksatu(col,row,1)-p_tksatu(i,1,1))*(tksatu(col,row,1)-p_tksatu(i,1,1))
               p_tksatu(i,2,2) = p_tksatu(i,2,2) + (tksatu(col,row,2)-p_tksatu(i,1,2))*(tksatu(col,row,2)-p_tksatu(i,1,2))
               p_slope(i,2)    = p_slope(i,2)    + (slope_avg(col,row)-p_slope(i,1))*(slope_avg(col,row)-p_slope(i,1))

            end if
         end do

         p_lai(i,2)      = sqrt(p_lai(i,2)      / p_lai(i,3))
         p_k_s(i,2,1)    = sqrt(p_k_s(i,2,1)    / p_k_s(i,3,1))
         p_k_s(i,2,2)    = sqrt(p_k_s(i,2,2)    / p_k_s(i,3,2))
         p_k_sl(i,2,1)   = sqrt(p_k_sl(i,2,1)   / p_k_sl(i,3,1))
         p_k_sl(i,2,2)   = sqrt(p_k_sl(i,2,2)   / p_k_sl(i,3,2))
         p_tkdry(i,2,1)  = sqrt(p_tkdry(i,2,1)  / p_tkdry(i,3,1))
         p_tkdry(i,2,2)  = sqrt(p_tkdry(i,2,2)  / p_tkdry(i,3,2))
         p_tksatf(i,2,1) = sqrt(p_tksatf(i,2,1) / p_tksatf(i,3,1))
         p_tksatf(i,2,2) = sqrt(p_tksatf(i,2,2) / p_tksatf(i,3,2))
         p_tksatu(i,2,1) = sqrt(p_tksatu(i,2,1) / p_tksatu(i,3,1))
         p_tksatu(i,2,2) = sqrt(p_tksatu(i,2,2) / p_tksatu(i,3,2))
         p_slope(i,2)    = sqrt(p_slope(i,2)    / p_slope(i,3))

      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,L,row,col)
      do i = sjx_points + 1,nmp(iter),1
         if(mp_id_new(i,1) == 0)then
            cycle
         end if
         do j = 0,mp_id_new(i,1)-1,1
            col = mp_ii_new(mp_id_new(i,2)+j,1)
            row = mp_ii_new(mp_id_new(i,2)+j,2)
            if(col == 0. .or. row == 0.)then
               cycle
            end if

            L = int(landtypes(col,row))

            if((L /= 0).and.(L /= maxlc))then

               p_lai(i,2)      = p_lai(i,2)      + (lai(col,row)-p_lai(i,1))*(lai(col,row)-p_lai(i,1))
               p_k_s(i,2,1)    = p_k_s(i,2,1)    + (k_s(col,row,1)-p_k_s(i,1,1))*(k_s(col,row,1)-p_k_s(i,1,1))
               p_k_s(i,2,2)    = p_k_s(i,2,2)    + (k_s(col,row,2)-p_k_s(i,1,2))*(k_s(col,row,2)-p_k_s(i,1,2))
               p_k_sl(i,2,1)   = p_k_sl(i,2,1)   + (k_sl(col,row,1)-p_k_sl(i,1,1))*(k_sl(col,row,1)-p_k_sl(i,1,1))
               p_k_sl(i,2,2)   = p_k_sl(i,2,2)   + (k_sl(col,row,2)-p_k_sl(i,1,2))*(k_sl(col,row,2)-p_k_sl(i,1,2))
               p_tkdry(i,2,1)  = p_tkdry(i,2,1)  + (tkdry(col,row,1)-p_tkdry(i,1,1))*(tkdry(col,row,1)-p_tkdry(i,1,1))
               p_tkdry(i,2,2)  = p_tkdry(i,2,2)  + (tkdry(col,row,2)-p_tkdry(i,1,2))*(tkdry(col,row,2)-p_tkdry(i,1,2))
               p_tksatf(i,2,1) = p_tksatf(i,2,1) + (tksatf(col,row,1)-p_tksatf(i,1,1))*(tksatf(col,row,1)-p_tksatf(i,1,1))
               p_tksatf(i,2,2) = p_tksatf(i,2,2) + (tksatf(col,row,2)-p_tksatf(i,1,2))*(tksatf(col,row,2)-p_tksatf(i,1,2))
               p_tksatu(i,2,1) = p_tksatu(i,2,1) + (tksatu(col,row,1)-p_tksatu(i,1,1))*(tksatu(col,row,1)-p_tksatu(i,1,1))
               p_tksatu(i,2,2) = p_tksatu(i,2,2) + (tksatu(col,row,2)-p_tksatu(i,1,2))*(tksatu(col,row,2)-p_tksatu(i,1,2))
               p_slope(i,2)    = p_slope(i,2)    + (slope_avg(col,row)-p_slope(i,1))*(slope_avg(col,row)-p_slope(i,1))

            end if
         end do

         p_lai(i,2)      = sqrt(p_lai(i,2)      / p_lai(i,3))
         p_k_s(i,2,1)    = sqrt(p_k_s(i,2,1)    / p_k_s(i,3,1))
         p_k_s(i,2,2)    = sqrt(p_k_s(i,2,2)    / p_k_s(i,3,2))
         p_k_sl(i,2,1)   = sqrt(p_k_sl(i,2,1)   / p_k_sl(i,3,1))
         p_k_sl(i,2,2)   = sqrt(p_k_sl(i,2,2)   / p_k_sl(i,3,2))
         p_tkdry(i,2,1)  = sqrt(p_tkdry(i,2,1)  / p_tkdry(i,3,1))
         p_tkdry(i,2,2)  = sqrt(p_tkdry(i,2,2)  / p_tkdry(i,3,2))
         p_tksatf(i,2,1) = sqrt(p_tksatf(i,2,1) / p_tksatf(i,3,1))
         p_tksatf(i,2,2) = sqrt(p_tksatf(i,2,2) / p_tksatf(i,3,2))
         p_tksatu(i,2,1) = sqrt(p_tksatu(i,2,1) / p_tksatu(i,3,1))
         p_tksatu(i,2,2) = sqrt(p_tksatu(i,2,2) / p_tksatu(i,3,2))
         p_slope(i,2)    = sqrt(p_slope(i,2)    / p_slope(i,3))

      end do
!$OMP END PARALLEL DO

      write(it,'(i2.2)')iter
      lndname =  trim(base_dir) // trim(EXPNME) // "/makegrid/threshold/threshold_NXP" // trim(nxpc) // "_" // trim(adjustl(it)) // ".nc4"
      print*,lndname
      CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber,nf90_netcdf4), ncID))
      CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", nmp(iter), spDimID))
      CALL CHECK(NF90_DEF_DIM(ncID, "dima"      , 2        , nmDimID))
      CALL CHECK(NF90_DEF_VAR(ncID, "num_landtypes"    , NF90_INT  , (/ spDimID /)         , varid(1)))
      CALL CHECK(NF90_DEF_VAR(ncID, "fraction_mainarea", NF90_float, (/ spDimID, nmDimID /), varid(2)))
      CALL CHECK(NF90_DEF_VAR(ncID, "p_slope"          , NF90_float, (/ spDimID, nmDimID /), varid(3)))
      CALL CHECK(NF90_DEF_VAR(ncID, "p_lai"            , NF90_float, (/ spDimID, nmDimID /), varid(4)))
      CALL CHECK(NF90_DEF_VAR(ncID, "p_k_s"            , NF90_float, (/ spDimID, nmDimID, nmDimID /), varid(5)))
      CALL CHECK(NF90_DEF_VAR(ncID, "p_k_sl"           , NF90_float, (/ spDimID, nmDimID, nmDimID /), varid(6)))
      CALL CHECK(NF90_DEF_VAR(ncID, "p_tkdry"          , NF90_float, (/ spDimID, nmDimID, nmDimID /), varid(7)))
      CALL CHECK(NF90_DEF_VAR(ncID, "p_tksatf"         , NF90_float, (/ spDimID, nmDimID, nmDimID /), varid(8)))
      CALL CHECK(NF90_DEF_VAR(ncID, "p_tksatu"         , NF90_float, (/ spDimID, nmDimID, nmDimID /), varid(9)))
      CALL CHECK(NF90_ENDDEF(ncID))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(1), n_landtypes))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(2), fraction_mainarea))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(3), p_slope(:,1:2)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(4), p_lai(:,1:2)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(5), p_k_s(:,1:2,:)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(6), p_k_sl(:,1:2,:)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(7), p_tkdry(:,1:2,:)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(8), p_tksatf(:,1:2,:)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(9), p_tksatu(:,1:2,:)))
      CALL CHECK(NF90_CLOSE(ncID))

      print*,"n_landtypes",minval(n_landtypes),maxval(n_landtypes)
      print*,"mainarea_id",minval(fraction_mainarea(:,1)),maxval(fraction_mainarea(:,1))
      print*,"fraction_mainarea",minval(fraction_mainarea(:,2)),maxval(fraction_mainarea(:,2))
      print*,"Mean:"
      print*,"slope",minval(p_slope(:,1)),maxval(p_slope(:,1))
      print*,"lai",minval(p_lai(:,1)),maxval(p_lai(:,1))
      print*,"k_s_l1",minval(p_k_s(:,1,1)),maxval(p_k_s(:,1,1))
      print*,"k_s_l2",minval(p_k_s(:,1,2)),maxval(p_k_s(:,1,2))
      print*,"k_solids_l1",minval(p_k_sl(:,1,1)),maxval(p_k_sl(:,1,1))
      print*,"k_solids_l2",minval(p_k_sl(:,1,2)),maxval(p_k_sl(:,1,2))
      print*,"tkdry_l1",minval(p_tkdry(:,1,1)),maxval(p_tkdry(:,1,1))
      print*,"tkdry_l2",minval(p_tkdry(:,1,2)),maxval(p_tkdry(:,1,2))
      print*,"tksatf_l1",minval(p_tksatf(:,1,1)),maxval(p_tksatf(:,1,1))
      print*,"tksatf_l2",minval(p_tksatf(:,1,2)),maxval(p_tksatf(:,1,2))
      print*,"tksatu_l1",minval(p_tksatu(:,1,1)),maxval(p_tksatu(:,1,1))
      print*,"tksatu_l2",minval(p_tksatu(:,1,2)),maxval(p_tksatu(:,1,2))
      print*,"Standard Deviation:"
      print*,"slope",minval(p_slope(:,2)),maxval(p_slope(:,2))
      print*,"lai",minval(p_lai(:,2)),maxval(p_lai(:,2))
      print*,"k_s_l1",minval(p_k_s(:,2,1)),maxval(p_k_s(:,2,1))
      print*,"k_s_l2",minval(p_k_s(:,2,2)),maxval(p_k_s(:,2,2))
      print*,"k_solids_l1",minval(p_k_sl(:,2,1)),maxval(p_k_sl(:,2,1))
      print*,"k_solids_l2",minval(p_k_sl(:,2,2)),maxval(p_k_sl(:,2,2))
      print*,"tkdry_l1",minval(p_tkdry(:,2,1)),maxval(p_tkdry(:,2,1))
      print*,"tkdry_l2",minval(p_tkdry(:,2,2)),maxval(p_tkdry(:,2,2))
      print*,"tksatf_l1",minval(p_tksatf(:,2,1)),maxval(p_tksatf(:,2,1))
      print*,"tksatf_l2",minval(p_tksatf(:,2,2)),maxval(p_tksatf(:,2,2))
      print*,"tksatu_l1",minval(p_tksatu(:,2,1)),maxval(p_tksatu(:,2,1))
      print*,"tksatu_l2",minval(p_tksatu(:,2,2)),maxval(p_tksatu(:,2,2))


!------------------------------------
! 5.6 Calculate the number of meshes to be refined based on the threshold
!------------------------------------
      ref = 0
      IsInRfArea = .false.
   
      do i = nmp(iter-1)+1,nmp(iter),1

         if((fraction_mainarea(i,1) == 0.).or.(fraction_mainarea(i,1) == maxlc))then
            cycle
         end if

if (refine_num_landtypes == .True. .and. n_landtypes(i)>th_num_landtypes) then

   ref(i) = 1
   ref_tr(i,1) = 1
end if

if (refine_area_mainland == .True. .and. fraction_mainarea(i, 2) < th_area_mainland) then
   ref(i) = 1
   ref_tr(i,2) = 1
end if

if (refine_lai_m == .True. .and. (p_lai(i, 1) > th_lai_m)) then
   ref(i) = 1
   ref_tr(i,3) = 1
end if

if (refine_lai_s == .True. .and. (p_lai(i, 2) > th_lai_s)) then
   ref(i) = 1
   ref_tr(i,4) = 1
end if

if (refine_slope_m == .True. .and. (p_slope(i, 1) > th_slope_m)) then

   ref(i) = 1
   ref_tr(i,5) = 1
end if

if (refine_slope_s == .True. .and. (p_slope(i, 2) > th_slope_s)) then

   ref(i) = 1
   ref_tr(i,6) = 1
end if

if (refine_k_s_m == .True. .and. ((p_k_s(i, 1, 1) > th_k_s_m).or.(p_k_s(i, 1, 2) > th_k_s_m))) then
   ref(i) = 1
   ref_tr(i,7) = 1
end if

if (refine_k_s_s == .True. .and. ((p_k_s(i, 2, 1) > th_k_s_s).or.(p_k_s(i, 2, 2) > th_k_s_s))) then
   ref(i) = 1
   ref_tr(i,8) = 1
end if

if (refine_k_solids_m == .True. .and. ((p_k_sl(i, 1, 1) > th_k_solids_m).or.(p_k_sl(i, 1, 2) > th_k_solids_m))) then
   ref(i) = 1
   ref_tr(i,9) = 1
end if

if (refine_k_solids_s == .True. .and. ((p_k_sl(i, 2, 1) > th_k_solids_s).or.(p_k_sl(i, 2, 2) > th_k_solids_s))) then
   ref(i) = 1
   ref_tr(i,10) = 1
end if

if (refine_tkdry_m == .True. .and. ((p_tkdry(i, 1, 1) > th_tkdry_m).or.(p_tkdry(i, 1, 2) > th_tkdry_m))) then
   ref(i) = 1
   ref_tr(i,11) = 1
end if

if (refine_tkdry_s == .True. .and. ((p_tkdry(i, 2, 1) > th_tkdry_s).or.(p_tkdry(i, 2, 2) > th_tkdry_s))) then
   ref(i) = 1
   ref_tr(i,12) = 1
end if

if (refine_tksatf_m == .True. .and. ((p_tksatf(i, 1, 1) > th_tksatf_m).or.(p_tksatf(i, 1, 2) > th_tksatf_m))) then
   ref(i) = 1
   ref_tr(i,13) = 1
end if

if (refine_tksatf_s == .True. .and. ((p_tksatf(i, 2, 1) > th_tksatf_s).or.(p_tksatf(i, 2, 2) > th_tksatf_s))) then
   ref(i) = 1
   ref_tr(i,14) = 1
end if

if (refine_tksatu_m == .True. .and. ((p_tksatu(i, 1, 1) > th_tksatu_m).or.(p_tksatu(i, 1, 2) > th_tksatu_m))) then
   ref(i) = 1
   ref_tr(i,15) = 1
end if

if (refine_tksatu_s == .True. .and. ((p_tksatu(i, 2, 1) > th_tksatu_s).or.(p_tksatu(i, 2, 2) > th_tksatu_s))) then

   ref(i) = 1
   ref_tr(i,16) = 1
end if

      end do

      num_ref = INT(sum(ref))

!------------------------------------
! 5.7 Set flag based on the number of grids
!------------------------------------

      if(num_ref == 0)then
         isover = .true.
         print*,"There are no more grids to refine"
      end if

      refed = 0

      if(iter >= max_iter)then
         isover = .true.
         print*,"Reached the preset maximum number of iterations"
      end if

   end do

   print*,"All iterations are completed, and the total number of iterations is",iter

!--------------------------------------------------------------------------
! 6. Record the end result
!--------------------------------------------------------------------------
      write(it,'(i2.2)')iter
      !lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/result/gridfile.nc4"
      lndname = trim(base_dir) // trim(EXPNME) // "/makegrid/result/gridfile_NXP" // trim(nxpc) // "_sjx.nc4"
      print*,lndname
      CALL CHECK(NF90_CREATE(trim(lndname), ior(nf90_clobber,nf90_netcdf4), ncID))
      CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", nmp(iter)  , spDimID))
      CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", nwp(iter)  , lpDimID))
      CALL CHECK(NF90_DEF_DIM(ncID, "dim_a"     , 3          , thDimID))
      CALL CHECK(NF90_DEF_DIM(ncID, "dim_d"     , 16         , sxDimID))
      CALL CHECK(NF90_DEF_DIM(ncID, "dim_c"     , 7          , seDimID))
      CALL CHECK(NF90_DEF_VAR(ncID, "GLONM"     , NF90_FLOAT, (/ spDimID /), varid(1)))
      CALL CHECK(NF90_DEF_VAR(ncID, "GLATM"     , NF90_FLOAT, (/ spDimID /), varid(2)))
      CALL CHECK(NF90_DEF_VAR(ncID, "GLONW"     , NF90_FLOAT, (/ lpDimID /), varid(3)))
      CALL CHECK(NF90_DEF_VAR(ncID, "GLATW"     , NF90_FLOAT, (/ lpDimID /), varid(4)))
      CALL CHECK(NF90_DEF_VAR(ncID, "itab_m%iw" , NF90_INT  , (/ thDimID, spDimID /), varid(5)))
      CALL CHECK(NF90_DEF_VAR(ncID, "itab_w%im" , NF90_INT  , (/ seDimID, lpDimID /), varid(6)))
      CALL CHECK(NF90_DEF_VAR(ncID, "ref_tr"    , NF90_INT  , (/ spDimID, sxDimID /), varid(7)))
      CALL CHECK(NF90_ENDDEF(ncID))

      CALL CHECK(NF90_PUT_VAR(ncID, varid(1), mp_new(1:nmp(iter),1)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(2), mp_new(1:nmp(iter),2)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(3), wp_new(1:nwp(iter),1)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(4), wp_new(1:nwp(iter),2)))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(5), ngrmw_new(1:3,1:nmp(iter))))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(6), ngrwm_new(1:7,1:nwp(iter))))
      CALL CHECK(NF90_PUT_VAR(ncID, varid(7), ref_tr(1:nmp(iter),1:16)))
      CALL CHECK(NF90_CLOSE(ncID))


   deallocate(mp)
   deallocate(mp_new)
   deallocate(wp)
   deallocate(wp_new)
   deallocate(mp_id)
   deallocate(mp_id_new)
   deallocate(mp_ii)
   deallocate(mp_ii_new)
   deallocate(ngrmw)
   deallocate(ngrmw_new)


   end subroutine refine_sjx


   subroutine CHECK(STATUS)
      INTEGER, intent (in) :: STATUS
      if  (STATUS .NE. NF90_NOERR) then 
         print *, NF90_STRERROR(STATUS)
         stop 'stopped'
      endif
   end subroutine CHECK


   INTEGER FUNCTION IsCrossLine(x1, x2)

      implicit none

      real(r8),intent(in) :: x1,x2
      IsCrossLine = 0

      if(abs(x1-x2) > 180.)then
         IsCrossLine = 1
      end if

   END FUNCTION IsCrossLine


   SUBROUTINE CheckLon(x)

      implicit none

      real(r8),intent(out) :: x

      if(x > 180.)then
         x = x - 360.
      else if(x < -180.)then
         x = x + 360.
      end if

   END SUBROUTINE CheckLon


   SUBROUTINE GetRefineJW(a,b,c,ref_jw)

      implicit none

      real(r8) :: a(2),b(2),c(2),dx,dy
      real(r8) :: lon(nlons_source),lat(nlats_source)
      real(r8) :: minlon,minlat,maxlon,maxlat
      integer,intent(out) :: ref_jw(nlons_source,nlats_source)
      integer :: i,j, icl(3)

      dy = 180. / nlats_source
      dx = 360. / nlons_source

      do i = 1,nlons_source,1
         lon(i) = -180. + (2 * i - 1) * dx / 2.
      end do

      do i = 1,nlats_source,1
         lat(i) = 90. - (2 * i - 1) * dy / 2.
      end do

      icl = 0
      icl(1) = IsCrossLine(a(1),b(1))
      icl(2) = IsCrossLine(a(1),c(1))
      icl(3) = IsCrossLine(b(1),c(1))

      if(INT(sum(icl)) > 0)then
         if(a(1) < 0.)then
            a(1) = a(1) + 360.
         end if
         if(b(1) < 0.)then
            b(1) = b(1) + 360.
         end if
         if(b(1) < 0.)then
            b(1) = b(1) + 360.
         end if
         do i = 1,nlons_source,1
            if(lon(i) < 0.)then
               lon(i) = lon(i) + 360.
            end if
         end do
      end if

      minlon = min(a(1),b(1),c(1))
      minlat = min(a(2),b(2),c(2))
      maxlon = max(a(1),b(1),c(1))
      maxlat = max(a(2),b(2),c(2))

      do i = 1,nlons_source,1
         if((lon(i)>minlon).and.(lon(i)<maxlon))then
            do j = 1,nlats_source,1
               if((lat(j)>minlat).and.(lat(j)<maxlat))then
                  ref_jw(i,j) = 1
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

      real(r8),intent(out) :: area(nlons_source,nlats_source)
      integer :: i,j
      real(r8) :: re,pi,deg2rad,global,dx,dy,error
      real(r8) :: lats(nlats_source),latn(nlats_source)
      real(r8) :: lonw(nlons_source),lone(nlons_source)

      re = 6.37122e6 * 0.001                    ! kilometer
      pi = 4.*atan(1.)
      deg2rad = pi/180.
      global = 0.

      dx = 360. / nlons_source
      dy = 180. / nlats_source

      do i = 1,nlons_source,1
         lone(i) = -180. + i * dx
         lonw(i) = -180. + (i - 1) * dx
      end do

      do i = 1,nlats_source,1
         latn(i) = 90. - (i - 1) * dy
         lats(i) = 90. - i * dy 
      end do

!$OMP PARALLEL DO NUM_THREADS(96) SCHEDULE(DYNAMIC,1)&
!$OMP PRIVATE(i,j,dx,dy)
      do j = 1,nlats_source,1
         do i = 1,nlons_source,1
            if(lone(i)<lonw(i))then   ! west edge is more western than data line
               dx = (lone(i)-lonw(i)+360.0)*deg2rad
            else
               dx = (lone(i)-lonw(i))*deg2rad
            endif
            if(latn(j)>lats(j)) then          ! north to south grid
               dy = sin(latn(j)*deg2rad) - sin(lats(j)*deg2rad)
            else                              ! south to north grid
               dy = sin(lats(j)*deg2rad) - sin(latn(j)*deg2rad)
            end if
            area(i,j) = dx*dy*re*re

         end do
      end do
!$OMP END PARALLEL DO

      global = sum(area(:,:))

      dx = (180. - (-180.)) * deg2rad
      dy = sin(90.*deg2rad) - sin(-90.*deg2rad)
      error = dx*dy*re*re
      if(abs(global-error)/error > 1.0e-7) then
         print*, 'CELLAREA error: correct area is ',error, &
              ' but summed area of grid cells is ', global
      end if

      return

   END SUBROUTINE CellArea


   REAL FUNCTION IsInUstrGrid(ustr,lon,lat)

      implicit none

      integer :: inc(4),i,j,iscross_l(2),stat,num_inter
      integer,allocatable :: iscross_g(:)                   ! Determine whether the unstructured grid line segment crosses the latitude and longitude grid
      real(r8),intent(in) :: lon,lat
      real(r8) :: minlat,maxlat,dx,dy,minlon,maxlon
      real(r8),intent(in) :: ustr(3,2)     ! Vertex of unstructured grid element
      real(r8) :: ustr_move(3,2)                ! Non-structural grid longitude
      real(r8) :: point(4,2)                      ! Vertex of the latitude and longitude grid
      real(r8) :: center_point(2)                 ! Unstructured grid cell center point
      real(r8) :: interarea_points(20,2)          ! Vertex of the area where two grids intersect
      real(r8),allocatable :: inter_points(:,:,:)           ! The intersection of two grids
      real(r8),allocatable :: area(:)   ! The area of a triangle consisting of two adjacent points of an unstructured grid and the vertices of any latitude and longitude grid
      real(r8),dimension(5) :: area_i   ! The area of a triangle composed of two adjacent points of the latitude and longitude grid and the vertices of any unstructured grid
      real(r8) :: tmpa(2),tmpb(2),tmpc(2),tmpd(3,2)

      maxlat = 0.
      minlat = 0.
      maxlon = 0.
      minlon = 0.

      iscross_l = 0
      iscross_g = 0
      inc = 0                 ! Determine the number of latitude and longitude grid vertices in an unstructured grid
      IsInUstrGrid = 0        ! The position relationship between latitude and longitude grid and unstructured grid is determined
      center_point = 0        ! Center point of latitude and longitude grid
      num_inter = 0           ! Number of polygon vertices of overlapping area of two meshes
      interarea_points = 0    ! Two kinds of mesh overlap polygon vertices

      dx = 360. / nlons_source
      dy = 180. / nlats_source

      ! Calculate the vertex coordinates of the latitude and longitude grid
      point(1,1) = lon + dx/2.
      point(1,2) = lat + dy/2.
      point(2,1) = lon - dx/2.
      point(2,2) = lat + dy/2.
      point(3,1) = lon - dx/2.
      point(3,2) = lat - dy/2.
      point(4,1) = lon + dx/2.
      point(4,2) = lat - dy/2.

      ! Ensure that the absolute latitude of the latitude grid does not exceed 90°
      do i = 1,4,1
         if(point(i,2) > 90.)then
            point(i,2) = 90.
         else if(point(i,2) < -90.)then
            point(i,2) = -90.
         end if
      end do

!------------------------------------------------------------------
! Preliminary screening based on latitude
!------------------------------------------------------------------

      maxlat = maxval(ustr(1:3,2))
      minlat = minval(ustr(1:3,2))

      if((maxlat > 85.).or.(minlat < -85.))then
         IsInUstrGrid = -1
         return
      end if

      if((point(4,2) > maxlat).or.(point(1,2) < minlat))then
         IsInUstrGrid = -1
         return
      end if

      ustr_move = ustr

      allocate(inter_points(3,3,2))
      allocate(area(4))
      allocate(iscross_g(3))
      inter_points = 0
      iscross_g = 0
      area = 0

!-------------------------------------------------------------------------------------
! Determine whether the two grids cross ±180° longitude
!-------------------------------------------------------------------------------------
      if(point(1,1) > 180.)then
         iscross_l(1) = 1
      else if(point(2,1) < -180.)then
         iscross_l(1) = -1
      end if

      iscross_l(2) = IsCrossLine2(ustr_move(:,1),3)

!-------------------------------------------------------------------------------------
! Move the grid point longitude according to the above judgment
!-------------------------------------------------------------------------------------
      if((iscross_l(1) /= 0).or.(iscross_l(2) == 1))then
         if(iscross_l(1) == -1)then
            stat = MoveLons(point(:,1),4,2)
            iscross_l(1) = 1
         else
            stat = MoveLons(point(:,1),4,1)
         end if
         if(iscross_l(2) == 1)then
            stat = MoveLons(ustr_move(:,1),3,1)
         else
            stat = MoveLons(ustr_move(:,1),3,2)
         end if
      end if

!-------------------------------------------------------------------------------------
! Filter the grid by longitude
!-------------------------------------------------------------------------------------
      minlon = minval(ustr_move(1:3,1))
      maxlon = maxval(ustr_move(1:3,1))

      if((point(2,1) > maxlon).or.(point(1,1) < minlon))then
         IsInUstrGrid = -1
         return
      end if

!-------------------------------------------------------------------------------------
! Start to determine the position relationship between the two grids and record key points
!-------------------------------------------------------------------------------------
      tmpa = ustr_move(1, 1:2)
      tmpb = ustr_move(2, 1:2)
      tmpd = inter_points(1, 1:3, 1:2)
      iscross_g(1) = IsCrossGrid(point, tmpa, tmpb, tmpd)
      !iscross_g(1) = IsCrossGrid(point,ustr_move(1,:),ustr_move(2,:),inter_points(1,:,:))
      tmpa = ustr_move(2, 1:2)
      tmpb = ustr_move(3, 1:2)
      tmpd = inter_points(2, 1:3, 1:2)
      iscross_g(2) = IsCrossGrid(point, tmpa, tmpb, tmpd)
      !iscross_g(2) = IsCrossGrid(point,ustr_move(2,:),ustr_move(3,:),inter_points(2,:,:))
      tmpa = ustr_move(1, 1:2)
      tmpb = ustr_move(3, 1:2)
      tmpd = inter_points(3, 1:3, 1:2)
      iscross_g(3) = IsCrossGrid(point, tmpa, tmpb, tmpd)
      !iscross_g(3) = IsCrossGrid(point,ustr_move(1,:),ustr_move(3,:),inter_points(3,:,:))      

      ! Calculate the number of vertices in the unstructured grid
      ! If it is 4, it contains, otherwise it intersects
      do i = 1,4,1
         area = 0.
         tmpa = ustr_move(1, 1:2)
         tmpb = ustr_move(2, 1:2)
         tmpc = point(i, 1:2)
         area(1) = GetTriangleArea(tmpa, tmpb, tmpc)
         !area(1) = GetTriangleArea(ustr_move(1,:),ustr_move(2,:),point(i,:))
         tmpa = ustr_move(2, 1:2)
         tmpb = ustr_move(3, 1:2)
         tmpc = point(i, 1:2)
         area(2) = GetTriangleArea(tmpa, tmpb, tmpc)
         !area(2) = GetTriangleArea(ustr_move(2,:),ustr_move(3,:),point(i,:))
         tmpa = ustr_move(1, 1:2)
         tmpb = ustr_move(3, 1:2)
         tmpc = point(i, 1:2)
         area(3) = GetTriangleArea(tmpa, tmpb, tmpc)
         !area(3) = GetTriangleArea(ustr_move(1,:),ustr_move(3,:),point(i,:))
         tmpa = ustr_move(1, 1:2)
         tmpb = ustr_move(2, 1:2)
         tmpc = ustr_move(3, 1:2)
         area(4) = GetTriangleArea(tmpa, tmpb, tmpc)
         !area(4) = GetTriangleArea(ustr_move(1,:),ustr_move(2,:),ustr_move(3,:))

         if(abs(area(1) + area(2) + area(3) - area(4)) < 0.0000001)then
            inc(i) = 1
         end if
      end do
                
      !print*,"inc",inc

!----------------------------------------------------------------------------------------------------                
! Calculate the overlapping area between the unstructured grid and the latitude and longitude grid
! 1. Find points in the latitude and longitude grid of the unstructured grid
! 2. Find points in the unstructured grid where the latitude and longitude grid is located
! 3. Find the intersection point between the unstructured grid and the latitude and longitude grid
! 4. Sort them and combine them into convex polygons (or triangles)
! 5. Calculate the area of the figure
!---------------------------------------------------------------------------------------------------

      if(sum(inc) == 4)then   ! If contain
         IsInUstrGrid = 1
         return
      else if((sum(inc) == 0).and.(sum(iscross_g) == 0))then  ! It does not intersect if it is not included
         IsInUstrGrid = -1
         return
      else
         do i = 1,4,1
            if(inc(i) == 1)then
               num_inter = num_inter + 1
               interarea_points(num_inter,:) = point(i,:)
            end if
         end do

         center_point = 0.
         do i = 1,4,1
            center_point = center_point + point(i,:)
         end do
         center_point = center_point / 4

         do i = 1,3,1
            area_i = 0.
            do j = 1,3,1
               tmpa = ustr_move(i, 1:2)
               tmpb = point(j, 1:2)
               tmpc = point(j + 1, 1:2)
               area_i(j) = GetTriangleArea(tmpa, tmpb, tmpc)
               !area_i(j) = GetTriangleArea(ustr_move(i,:),point(j,:),point(j+1,:))
               tmpb = point(j, 1:2)
               tmpc = point(j + 1, 1:2)
               area_i(5) = area_i(5) + GetTriangleArea(center_point, tmpb, tmpc)
               !area_i(5) = area_i(5) + GetTriangleArea(center_point,point(j,:),point(j+1,:))
            end do

            tmpa = ustr_move(i, 1:2)
            tmpb = point(1, 1:2)
            tmpc = point(4, 1:2)
            area_i(4) = GetTriangleArea(tmpa, tmpb, tmpc)
            !area_i(4) = GetTriangleArea(ustr_move(i,:),point(1,:),point(4,:))
            tmpb = point(1, 1:2)
            tmpc = point(4, 1:2)
            area_i(5) = area_i(5) + GetTriangleArea(center_point, tmpb, tmpc)
            !area_i(5) = area_i(5) + GetTriangleArea(center_point,point(1,:),point(4,:))

            if(abs(sum(area_i(1:4)) - area_i(5)) < 0.00006)then
               num_inter = num_inter + 1
               interarea_points(num_inter,:) = ustr_move(i,:)
            end if
         end do

         do i = 1,3,1
            if(inter_points(i,3,1) /= 0)then
               num_inter = num_inter + 1
               interarea_points(num_inter:num_inter + inter_points(i, 3, 1) - 1, :) = &
                  inter_points(i, 1:inter_points(i, 3, 1), :)
               num_inter = num_inter + inter_points(i,3,1) - 1
            end if
         end do

         if(num_inter == 0)then
            IsInUstrGrid = 0.
            return
         end if

         CALL SortPoints(interarea_points, num_inter)

         IsInUstrGrid = GetAreaPercent(interarea_points, num_inter, point)

      end if

   END FUNCTION IsInUstrGrid


   ! Determine whether the grid crosses the 189° and -180° longitude lines
   INTEGER FUNCTION IsCrossLine2(lons,num)

      implicit none

      integer :: num,i,j
      real(r8),dimension(num) :: lons

      IsCrossLine2 = 0

      do i = 1,num - 1,1
         do j = i + 1,num,1
            if(abs(lons(j) - lons(i)) > 180.)then
               IsCrossLine2 = 1
               return
            end if
         end do
      end do

   END FUNCTION IsCrossLine2


   ! Adjust the longitude of the grid points
   INTEGER FUNCTION MoveLons(lons,num,lor)         ! lor = left or right

      implicit none

      integer :: i
      integer,intent(in) :: num,lor
      real(r8),dimension(num) :: lons

      if(lor == 1)then         
         do i = 1,num,1
            if(lons(i) < 0)then
               lons(i) = lons(i) + 360.
            end if
         end do
      else if((lor == 2).and.(lons(1) < 0))then    
         do i = 1,num,1
            lons(i) = lons(i) + 360.
         end do
      end if

      MoveLons = 1

   END FUNCTION MoveLons


   ! Determine whether the unstructured grid line segment crosses the latitude and longitude grid
   INTEGER FUNCTION IsCrossGrid(point,a,b,inter_point)

      implicit none

      real(r8),dimension(2) :: a,b   
      real(r8),dimension(2) :: x,y
      real(r8),dimension(4,2) :: point
      real(r8),dimension(3,2) :: inter_point
      real(r8) :: x1,x2,y1,y2,m,n,num

      IsCrossGrid = 0
      inter_point = 0
      num = 0

      x(1) = max(a(1),b(1))
      x(2) = min(a(1),b(1))
      y(1) = max(a(2),b(2))
      y(2) = min(a(2),b(2))

      if(a(1) == b(1))then
         if(a(1)>point(2,1).and.a(1)<point(1,1))then
            if((y(1)>point(1,2)).and.((y(2)>point(4,2))).and.(y(2)<point(1,2)))then
               inter_point(1,1) = a(1)
               inter_point(1,2) = point(1,2)
               inter_point(3,1) = 1
               IsCrossGrid = 1
               return
            else if((y(1)>point(1,2)).and.((y(2)<point(4,2))))then
               inter_point(1,1) = a(1)
               inter_point(1,2) = point(1,2)
               inter_point(2,1) = a(1)
               inter_point(2,2) = point(4,2)
               inter_point(3,1) = 2
               IsCrossGrid = 2
               return
            else if((y(1)>point(4,2)).and.(y(1)<point(1,2)).and.(y(2)<point(4,2)))then
               inter_point(1,1) = a(1)
               inter_point(1,2) = point(4,2)
               inter_point(3,1) = 1
               IsCrossGrid = 1
               return
            end if
         end if
         IsCrossGrid = 0
         return
      else if(a(2) == b(2))then
         if(a(2)>point(4,2).and.a(2)<point(1,2))then
            if((x(1)>point(1,1)).and.((x(2)>point(2,1))).and.(x(2)<point(1,1)))then
               inter_point(1,1) = point(1,1)
               inter_point(1,2) = a(2)
               inter_point(3,1) = 1
               IsCrossGrid = 1
               return
            else if((x(1)>point(1,1)).and.((x(2)<point(2,1))))then
               inter_point(1,1) = point(1,1)
               inter_point(1,2) = a(2)
               inter_point(2,1) = point(2,1)
               inter_point(2,2) = a(2)
               inter_point(3,1) = 2
               IsCrossGrid = 2
               return
            else if((x(1)>point(2,1)).and.(x(1)<point(1,1)).and.(x(2)<point(2,1)))then
               inter_point(1,1) = point(2,1)
               inter_point(1,2) = a(2)
               inter_point(3,1) = 1
               IsCrossGrid = 1
               return
            end if
         end if
         IsCrossGrid = 0
         return
      else
         m = (a(2) - b(2)) / (a(1) - b(1))
         n = a(2) - m * a(1)
         y1 = m * point(1,1) + n
         y2 = m * point(2,1) + n
         x1 = (point(1,2) - n) / m
         x2 = (point(4,2) - n) / m
         if((y1>y(2)).and.(y1<y(1)).and.(y1>point(4,2)).and.(y1<point(1,2)))then
            IsCrossGrid = IsCrossGrid + 1
            inter_point(IsCrossGrid,1) = point(1,1)
            inter_point(IsCrossGrid,2) = y1
         end if
         if((y2>y(2)).and.(y2<y(1)).and.(y2>point(4,2)).and.(y2<point(1,2)))then
            IsCrossGrid = IsCrossGrid + 1
            inter_point(IsCrossGrid,1) = point(2,1)
            inter_point(IsCrossGrid,2) = y2
         end if
         if((x1>x(2)).and.(x1<x(1)).and.(x1>point(2,1)).and.(x1<point(1,1)))then
            IsCrossGrid = IsCrossGrid + 1
            inter_point(IsCrossGrid,1) = x1
            inter_point(IsCrossGrid,2) = point(1,2)
         end if
         if((x2>x(2)).and.(x2<x(1)).and.(x2>point(2,1)).and.(x2<point(1,1)))then
            IsCrossGrid = IsCrossGrid + 1
            inter_point(IsCrossGrid,1) = x2
            inter_point(IsCrossGrid,2) = point(4,2)
         end if

         if(IsCrossGrid > 2)then
            IsCrossGrid = 0
         end if

         inter_point(3,1) = IsCrossGrid
      end if

   END FUNCTION IsCrossGrid


   ! Calculate the area of the triangle (by latitude and longitude, not the actual area)
   REAL FUNCTION GetTriangleArea(a,b,c)

      implicit none

      real(r8),dimension(2),intent(in) :: a,b,c
      real(r8) :: aa,bb,cc,p

      GetTriangleArea = 0

      aa = sqrt((c(1)-b(1))*(c(1)-b(1))+(c(2)-b(2))*(c(2)-b(2)))
      bb = sqrt((c(1)-a(1))*(c(1)-a(1))+(c(2)-a(2))*(c(2)-a(2)))
      cc = sqrt((a(1)-b(1))*(a(1)-b(1))+(a(2)-b(2))*(a(2)-b(2)))

      p = (aa + bb + cc) / 2

      GetTriangleArea = sqrt(p*(p-aa)*(p-bb)*(p-cc))

   END FUNCTION GetTriangleArea


   ! Gets the proportion of unstructured grids that contain latitude and longitude grids
   REAL FUNCTION GetAreaPercent(inter_point,num,point)

      implicit none

      integer,intent(in) :: num
      real(r8),dimension(20,2),intent(in) :: inter_point
      real(r8),dimension(4,2),intent(in) :: point
      real(r8),dimension(2) :: center_point
      integer :: i
      real(r8) :: inter_area, tmpa(2), tmpb(2)

      GetAreaPercent = 0
      inter_area = 0
      center_point = 0

      do i = 1,num,1
         center_point = center_point + inter_point(i,:)
      end do
      center_point = center_point / num

      do i = 1,num-1,1
         tmpa = inter_point(i, 1:2)
         tmpb = inter_point(i + 1, 1:2)
         inter_area = inter_area + GetTriangleArea(center_point, tmpa, tmpb)
         !inter_area = inter_area + GetTriangleArea(center_point,inter_point(i,:),inter_point(i+1,:))
      end do

      tmpa = inter_point(1, 1:2)
      tmpb = inter_point(num, 1:2)
      inter_area = inter_area + GetTriangleArea(center_point, tmpa, tmpb)
      !inter_area = inter_area + GetTriangleArea(center_point,inter_point(1,:),inter_point(num,:))

      GetAreaPercent = inter_area / (abs((point(1,1) - point(2,1))*(point(1,2) - point(4,2))))
      ! The area ratio is the overlapping triangle divided by the area of the latitude and longitude grid

      !if(GetAreaPercent >= 1)then
      !        GetAreaPercent = 0
      !end if

   END FUNCTION GetAreaPercent


   ! Sort points into polygons
   SUBROUTINE SortPoints(points,num)

      implicit none

      integer :: i,j,x
      real(r8) :: angle_x,pi
      integer,intent(in) :: num
      integer,dimension(num) :: sort_i
      real(r8), dimension(20, 2),intent(out) :: points
      real(r8), dimension(num, 2) :: points_i
      real(r8),dimension(2) :: center_point
      real(r8),dimension(num) :: angle

      center_point = 0

      pi = 3.1415926535

      do i = 1,num,1
         center_point = center_point + points(i,:)
         !print*,points(i,:)
         sort_i(i) = i
      end do
      center_point = center_point / num

      !print*,"center_point",center_point

      do i = 1,num,1
         points_i(i,:) = points(i,:) - center_point
         if(points_i(i,2) >= 0)then
            if(points_i(i,1) == 0)then
               angle(i) = pi / 2
            else
               angle(i) = atan(points_i(i,2) / points_i(i,1))
               if(points_i(i,1) < 0)then
                  angle(i) = angle(i) + pi
               end if
            end if
         else
            if(points_i(i,1) == 0)then
               angle(i) = 1.5 * pi
            else if(points_i(i,1) < 0)then
               angle(i) = atan(points_i(i,2) / points_i(i,1)) + pi
            else
               angle(i) = atan(points_i(i,2) / points_i(i,1)) + 2 * pi
            end if
         end if
      end do

      do i = 1,num - 1,1
         do j = i + 1,num,1
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

      do i = 1,num,1
         points(i,:) = points_i(sort_i(i),:) + center_point
      end do

   END SUBROUTINE SortPoints


   SUBROUTINE IsInRefineArea(IsInRfArea,sjx_points,lbx_points,ngrmw,wp,nmp)

      implicit none

      integer :: i,j,n
      integer,intent(in) :: sjx_points,lbx_points,nmp
      integer,dimension(3,sjx_points*5),intent(in) :: ngrmw
      real(r8),dimension(lbx_points*5,2),intent(in) :: wp
      real(r8) :: sjx(3,2)

      logical,dimension(sjx_points),intent(out) :: IsInRfArea

      sjx = 0.
      IsInRfArea = .false.

      do i = 1,nmp,1
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

END module  MOD_refine_sjx
