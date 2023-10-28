! 1. Create the nlons_source*nlats_source slope and elevation threshold file
! 2. Create a triangle mesh or polygon mesh threshold file

module MOD_GetThreshold

   use consts_coms
   use netcdf

   implicit none

contains
  SUBROUTINE GetThreshold()

   integer :: sjx_points               ! Number of triangular grids
   integer :: maxid(1)                 ! Record the largest land type
   integer :: num_all                  ! Each unstructured grid contains the total number of latitude and longitude grids
   integer :: row,col,L
   integer :: i,j
   !integer :: maxlc                    ! Land type maximum number,17 or 24
   integer :: nmDimID,loDimID,laDimID,spDimID,upDimID
   integer :: iunit,ncid,ncvarid(10),varid,dimID_sjx,thDimID,dimID_lbx

   real(r8),allocatable :: mp(:,:)     ! Center point and area of triangular grid
   real(r8),allocatable :: mp_ii(:,:)  ! The triangular grid contains the latitude and longitude grid
   real(r8),allocatable :: area(:)     ! Area of each land type in unstructured grid 

   ! Threshold array
   real(r8),allocatable :: fraction_mainarea(:,:),dem(:,:)
   real(r8),allocatable :: temp(:,:),slope_max(:,:)
   real(r8),allocatable :: slope_avg(:,:),lai(:,:)
   real(r8),allocatable :: p_slope(:,:),p_lai(:,:),p_k_s(:,:,:)
   real(r8),allocatable :: p_k_sl(:,:,:),p_tkdry(:,:,:),p_tksatu(:,:,:)
   real(r8),allocatable :: p_tksatf(:,:,:),k_s(:,:,:),k_sl(:,:,:)
   real(r8),allocatable :: tkdry(:,:,:),tksatf(:,:,:),tksatu(:,:,:)
   real(r8),allocatable :: landtypes(:,:)
   integer,allocatable :: nlaa(:),n_landtypes(:)
   
   integer,allocatable :: mp_id(:,:)   ! The triangular mesh contains the number of polygonal meshes and the starting position in mp_ii
   character(LEN=256) :: lndname,inputfile,outputfile,nxpc

   write(nxpc, '(I3.3)') NXP
   outputfile = trim(base_dir) // trim(EXPNME) // "/makegrid/threshold/threshold_" // trim(nxpc) // ".nc4"
   inputfile = trim(base_dir) // trim(EXPNME) // '/makegrid/gridfile/gridfile_NXP' // trim(nxpc) // '.nc4'
   print*,inputfile

   CALL CHECK(NF90_OPEN(trim(inputfile),nf90_nowrite,ncid))

   CALL CHECK(NF90_INQ_DIMID(ncid,"sjx_points",dimID_sjx))
   CALL CHECK(NF90_INQUIRE_DIMENSION(ncid,dimID_sjx,len=sjx_points))
   print*,"sjx_points = ",sjx_points
   lndname =  trim(base_dir) // trim(EXPNME) // "/makegrid/contain/initial/mp"
   
   CALL CHECK(NF90_CLOSE(ncid))

   allocate(mp(sjx_points,3))
   allocate(mp_id(sjx_points,2))
   allocate(n_landtypes(sjx_points))
   mp = 0.
   mp_id = 0
   n_landtypes = 0

   iunit = 100
   OPEN(iunit,file=trim(lndname)//'.bin',form='unformatted',status='old')
   READ(iunit) mp
   close(iunit)

   iunit = 100
   OPEN(iunit,file=trim(lndname)//'_id.bin',form='unformatted',status='old')
   READ(iunit) mp_id
   close(iunit)

   num_all = INT(sum(mp_id(:,1)))
   allocate(mp_ii(num_all,4))

   iunit = 100
   OPEN(iunit,file=trim(lndname)//'_ii.bin',form='unformatted',status='old')
   READ(iunit) mp_ii
   close(iunit)

   allocate(landtypes(nlons_source,nlats_source))
   landtypes = 0.

   if(lcs == "igbp")then
      !maxlc = 17
      lndname = trim(source_dir) // 'landtype_igbp_update.nc'
   else
      !maxlc = 24
      lndname = trim(source_dir) // 'landtype_usgs_update.nc'
   end if
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"landtype",varid))
   CALL CHECK(NF90_GET_VAR(ncid,varid,landtypes))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"landtypes",minval(landtypes),maxval(landtypes)

   allocate(lai(nlons_source,nlats_source))
   allocate(k_s(nlons_source,nlats_source,2))
   allocate(k_sl(nlons_source,nlats_source,2))
   allocate(tkdry(nlons_source,nlats_source,2))
   allocate(tksatf(nlons_source,nlats_source,2))
   allocate(tksatu(nlons_source,nlats_source,2))
   lai = 0.
   k_s = 0.
   k_sl = 0.
   tkdry = 0.
   tksatf = 0.
   tksatu = 0.

   lndname = trim(source_dir) // 'lai.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"lai",varid))
   CALL CHECK(NF90_GET_VAR(ncid,varid,lai))
   CALL CHECK(NF90_CLOSE(ncid))

   lndname = trim(source_dir) // 'k_s.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"k_s_l1",ncvarid(1)))
   CALL CHECK(NF90_INQ_VARID(ncid,"k_s_l2",ncvarid(2)))
   CALL CHECK(NF90_GET_VAR(ncid,ncvarid(1),k_s(:,:,1)))
   CALL CHECK(NF90_GET_VAR(ncid,ncvarid(2),k_s(:,:,2)))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"k_s_l1",minval(k_s(:,:,1)),maxval(k_s(:,:,1))
   print*,"k_s_l2",minval(k_s(:,:,2)),maxval(k_s(:,:,2))

   lndname = trim(source_dir) // 'k_solids.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"k_solids_l1",ncvarid(1)))
   CALL CHECK(NF90_INQ_VARID(ncid,"k_solids_l2",ncvarid(2)))
   CALL CHECK(NF90_GET_VAR(ncid,ncvarid(1),k_sl(:,:,1)))
   CALL CHECK(NF90_GET_VAR(ncid,ncvarid(2),k_sl(:,:,2)))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"k_solids_l1",minval(k_sl(:,:,1)),maxval(k_sl(:,:,1))
   print*,"k_solids_l2",minval(k_sl(:,:,2)),maxval(k_sl(:,:,2))

   lndname = trim(source_dir) // 'tkdry.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"tkdry_l1",ncvarid(1)))
   CALL CHECK(NF90_INQ_VARID(ncid,"tkdry_l2",ncvarid(2)))
   CALL CHECK(NF90_GET_VAR(ncid,ncvarid(1),tkdry(:,:,1)))
   CALL CHECK(NF90_GET_VAR(ncid,ncvarid(2),tkdry(:,:,2)))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"tkdry_l1",minval(tkdry(:,:,1)),maxval(tkdry(:,:,1))
   print*,"tkdry_l2",minval(tkdry(:,:,2)),maxval(tkdry(:,:,2))

   lndname = trim(source_dir) // 'tksatf.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"tksatf_l1",ncvarid(1)))
   CALL CHECK(NF90_INQ_VARID(ncid,"tksatf_l2",ncvarid(2)))
   CALL CHECK(NF90_GET_VAR(ncid,ncvarid(1),tksatf(:,:,1)))
   CALL CHECK(NF90_GET_VAR(ncid,ncvarid(2),tksatf(:,:,2)))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"tksatf_l1",minval(tksatf(:,:,1)),maxval(tksatf(:,:,1))
   print*,"tksatf_l2",minval(tksatf(:,:,2)),maxval(tksatf(:,:,2))

   lndname = trim(source_dir) // 'tksatu.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"tksatu_l1",ncvarid(1)))
   CALL CHECK(NF90_INQ_VARID(ncid,"tksatu_l2",ncvarid(2)))
   CALL CHECK(NF90_GET_VAR(ncid,ncvarid(1),tksatu(:,:,1)))
   CALL CHECK(NF90_GET_VAR(ncid,ncvarid(2),tksatu(:,:,2)))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,"tksatu_l1",minval(tksatu(:,:,1)),maxval(tksatu(:,:,1))
   print*,"tksatu_l2",minval(tksatu(:,:,2)),maxval(tksatu(:,:,2))

   allocate(dem(nlons_source,nlats_source))             ! Longitude and latitude grid surface elevation
   dem = 0.

   lndname = trim(source_dir) // 'dem.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"topo",varid))
   CALL CHECK(NF90_GET_VAR(ncid,varid,dem))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,minval(dem),maxval(dem)


   allocate(slope_max(nlons_source,nlats_source))          ! The maximum slope value in the eight neighborhood of the latitude and longitude grid
   allocate(slope_avg(nlons_source,nlats_source))          ! The mean slope value in the eight neighborhood of the latitude and longitude grid
   slope_max=0.
   slope_avg=0.
   lndname = trim(source_dir) // 'slope_max.nc'
   print*,lndname
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"slope_max",varid))
   CALL CHECK(NF90_GET_VAR(ncid,varid,slope_max))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,minval(slope_max),maxval(slope_max)

   lndname = trim(source_dir) // 'slope_avg.nc'
   CALL CHECK(NF90_OPEN(lndname,nf90_nowrite,ncid))
   CALL CHECK(NF90_INQ_VARID(ncid,"slope_avg",varid))
   CALL CHECK(NF90_GET_VAR(ncid,varid,slope_avg))
   CALL CHECK(NF90_CLOSE(ncid))
   print*,minval(slope_avg),maxval(slope_avg)

   allocate(nlaa(0:maxlc))                           ! Whether the unstructured grid contains an indication of the land type
   allocate(area(0:maxlc))                         ! Area of each land type in unstructured grid
   allocate(fraction_mainarea(sjx_points,2))   ! Number and proportion of main land types in unstructured grid
   allocate(p_slope(sjx_points,3))           ! Maximum slope of latitude and longitude grid in unstructured grid
   allocate(p_lai(sjx_points,3))               ! Mean value, maximum value and number of lai in unstructured grids
   allocate(p_k_s(sjx_points,3,2))
   allocate(p_k_sl(sjx_points,3,2))
   allocate(p_tkdry(sjx_points,3,2))
   allocate(p_tksatf(sjx_points,3,2))
   allocate(p_tksatu(sjx_points,3,2))
   fraction_mainarea = 0.
   p_slope = 0.
   p_lai = 0.
   p_k_s = 0.
   p_k_sl = 0.
   p_tkdry = 0.
   p_tksatf = 0.
   p_tksatu = 0.

!$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,row,col,nlaa,L,area,maxid)
   do i = 1,sjx_points,1
      nlaa = 0
      area = 0.
      maxid = 0
      if(mp_id(i,1) == 0)then
         cycle
      end if
      do j = 0,mp_id(i,1)-1,1
         col = mp_ii(mp_id(i,2)+j,1)
         row = mp_ii(mp_id(i,2)+j,2)
         if(col == 0. .or. row == 0.)then
            cycle
         end if
         L = int(landtypes(col,row))
         nlaa(L) = 1
         area(L) = area(L) + mp_ii(mp_id(i,2)+j,4)

         if((L /= 0).and.(L /= maxlc))then

            nlaa(L) = 1
            area(L) = area(L) + mp_ii(mp_id(i,2)+j,4)

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
      n_landtypes(i) = INT(sum(nlaa))
      maxid = maxloc(area) - 1
      fraction_mainarea(i,1) = maxid(1)
      fraction_mainarea(i,2) = area(maxid(1))/mp(i,3)

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
!$OMP PRIVATE(i,j,row,col)
   do i = 1,sjx_points,1
      if(mp_id(i,1) == 0)then
         cycle
      end if
      do j = 0,mp_id(i,1)-1,1
         col = mp_ii(mp_id(i,2)+j,1)
         row = mp_ii(mp_id(i,2)+j,2)
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

   print*,outputfile
   CALL CHECK(NF90_CREATE(trim(outputfile), ior(nf90_clobber,nf90_netcdf4), ncID))
   CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points"      , sjx_points, upDimID))
   CALL CHECK(NF90_DEF_DIM(ncID, "dima"             , 2          , nmDimID))
   CALL CHECK(NF90_DEF_VAR(ncID, "num_landtypes"    , NF90_INT   , (/ upDimID /), ncVarID(1)))
   CALL CHECK(NF90_DEF_VAR(ncID, "fraction_mainarea", NF90_float , (/ upDimID, nmDimID /), ncVarID(2)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_slope"          , NF90_float , (/ upDimID, nmDimID /), ncVarID(3)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_lai"            , NF90_float , (/ upDimID, nmDimID /), ncVarID(4)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_k_s"            , NF90_float , (/ upDimID, nmDimID, nmDimID /), ncVarID(5)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_k_sl"           , NF90_float , (/ upDimID, nmDimID, nmDimID /), ncVarID(6)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_tkdry"          , NF90_float , (/ upDimID, nmDimID, nmDimID /), ncVarID(7)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_tksatf"         , NF90_float , (/ upDimID, nmDimID, nmDimID /), ncVarID(8)))
   CALL CHECK(NF90_DEF_VAR(ncID, "p_tksatu"         , NF90_float , (/ upDimID, nmDimID, nmDimID /), ncVarID(9)))
   CALL CHECK(NF90_ENDDEF(ncID))
   CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1), n_landtypes))
   CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(2), fraction_mainarea))
   CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(3), p_slope(:,1:2)))
   CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(4), p_lai(:,1:2)))
   CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5), p_k_s(:,1:2,:)))
   CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6), p_k_sl(:,1:2,:)))
   CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(7), p_tkdry(:,1:2,:)))
   CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(8), p_tksatf(:,1:2,:)))
   CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(9), p_tksatu(:,1:2,:)))
   CALL CHECK(NF90_CLOSE(ncID))
   
   print*,"n_landtypes",minval(n_landtypes),maxval(n_landtypes)
   print*,"mainarea_id",minval(fraction_mainarea(:,1)),maxval(fraction_mainarea(:,1))
   print*,"fraction_mainarea",minval(fraction_mainarea(:,2)),maxval(fraction_mainarea(:,2))
   print*,"Mean:"
   print*,"p_slope_s",minval(p_slope(:,1)),maxval(p_slope(:,1))
   print*,"p_lai",minval(p_lai(:,1)),maxval(p_lai(:,1))
   print*,"p_k_s_l1",minval(p_k_s(:,1,1)),maxval(p_k_s(:,1,1))
   print*,"p_k_s_l2",minval(p_k_s(:,1,2)),maxval(p_k_s(:,1,2))
   print*,"p_k_solids_l1",minval(p_k_sl(:,1,1)),maxval(p_k_sl(:,1,1))
   print*,"p_k_solids_l2",minval(p_k_sl(:,1,2)),maxval(p_k_sl(:,1,2))
   print*,"p_tkdry_l1",minval(p_tkdry(:,1,1)),maxval(p_tkdry(:,1,1))
   print*,"p_tkdry_l2",minval(p_tkdry(:,1,2)),maxval(p_tkdry(:,1,2))
   print*,"p_tksatf_l1",minval(p_tksatf(:,1,1)),maxval(p_tksatf(:,1,1))
   print*,"p_tksatf_l2",minval(p_tksatf(:,1,2)),maxval(p_tksatf(:,1,2))
   print*,"p_tksatu_l1",minval(p_tksatu(:,1,1)),maxval(p_tksatu(:,1,1))
   print*,"p_tksatu_l2",minval(p_tksatu(:,1,2)),maxval(p_tksatu(:,1,2))
   print*,"Standard Deviation:"
   print*,"p_slope_s",minval(p_slope(:,2)),maxval(p_slope(:,2))
   print*,"p_lai",minval(p_lai(:,2)),maxval(p_lai(:,2))
   print*,"p_k_s_l1",minval(p_k_s(:,2,1)),maxval(p_k_s(:,2,1))
   print*,"p_k_s_l2",minval(p_k_s(:,2,2)),maxval(p_k_s(:,2,2))
   print*,"p_k_solids_l1",minval(p_k_sl(:,2,1)),maxval(p_k_sl(:,2,1))
   print*,"p_k_solids_l2",minval(p_k_sl(:,2,2)),maxval(p_k_sl(:,2,2))
   print*,"p_tkdry_l1",minval(p_tkdry(:,2,1)),maxval(p_tkdry(:,2,1))
   print*,"p_tkdry_l2",minval(p_tkdry(:,2,2)),maxval(p_tkdry(:,2,2))
   print*,"p_tksatf_l1",minval(p_tksatf(:,2,1)),maxval(p_tksatf(:,2,1))
   print*,"p_tksatf_l2",minval(p_tksatf(:,2,2)),maxval(p_tksatf(:,2,2))
   print*,"p_tksatu_l1",minval(p_tksatu(:,2,1)),maxval(p_tksatu(:,2,1))
   print*,"p_tksatu_l2",minval(p_tksatu(:,2,2)),maxval(p_tksatu(:,2,2))

end subroutine GetThreshold

   subroutine CHECK(STATUS)
      INTEGER, intent (in) :: STATUS
      if  (STATUS .NE. NF90_NOERR) then 
         print *, NF90_STRERROR(STATUS)
         stop 'stopped'
      endif
   end subroutine CHECK
 
   
end module MOD_GetThreshold
