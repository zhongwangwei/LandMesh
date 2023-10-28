! 1. Calculate the inclusion relationship between unstructured grid and structured grid
! 2. Calculate the patchtypes file required for the mpi version
Module MOD_Get_Contain_Patch

    use consts_coms
    USE NETCDF

    implicit none

contains
    subroutine Get_Contain_PatchId()

        implicit none

        real(r8), allocatable :: lon_i(:), lat_i(:)               ! Center point of latitude and longitude grid
        real(r8), allocatable :: area_fine_gridcell(:, :)         ! Area of each unit of the latitude and longitude grid
        real(r8), allocatable :: mp(:, :), wp(:, :)               ! Center point and area of triangular mesh and polygon mesh
        real(r8), allocatable :: ustr_ii(:, :)                    ! Latitude and longitude meshes contained in triangular and polygonal meshes (before processing)
        real(r8), allocatable :: ustr_ii_new(:, :)                ! Latitude and longitude meshes contained in triangular and polygonal meshes (after processing)
        real(r8), allocatable :: patches_fraction(:, :)           ! 
        real(r8), allocatable :: lon_e(:), lon_w(:), lat_n(:), lat_s(:)    
        real(r8) :: dbx(7, 2), sjx(7, 2)                          ! Output the latitude and longitude of the polygon with the triangular mesh
        real(r8) :: maxlat, minlat                                 
        real(r8) :: isinply                                       
        real(r8) :: dx, dy                                         
        integer, allocatable :: ustr_id(:, :)                     ! The unstructured grid contains the number of latitude and longitude grids 
								  !and the starting position in ustr_ii (before processing).
        integer, allocatable :: ustr_id_new(:, :)                 ! The unstructured grid contains the number of latitude and longitude grids 
								  !and the starting position in ustr_ii (after processing).
        integer, allocatable :: ngrmw(:, :), ngrwm(:, :)          ! Neighborhood array of center points of triangular and polygon meshes
        real(r8), allocatable :: landtypes(:, :)                  ! Land type
        integer, allocatable :: seaorland(:, :)                   ! Determine whether the latitude and longitude grid is land or sea
        integer, allocatable :: patchtypes(:, :)                  ! Record the sequence number of the unstructured grid where the latitude and longitude grids reside
        integer, allocatable :: map(:, :, :)                      ! Mapping of source grid and dest grid
        integer :: sjx_points, lbx_points, i, j, k, l, ncid, varid(6), sum_land, sum_sea
        integer :: num_i, id, row, col, x, y
        integer :: ustr_points                                    ! Number of grids
        integer :: numpatch                                       ! Each unstructured grid contains the total number of latitude and longitude grids
        integer :: num_null                                       ! The number of null values that need to be processed
        integer :: lodimid, ladimid, dimID_sjx, dimID_lbx
        character(LEN = 256) :: outputfile, flnm, nxpc
        character(LEN = 20) :: p_name(6) = (/"GLONW", "GLATW", "GLONM", "GLATM", "itab_w%im", "itab_m%iw"/)
        logical :: ispart                                         ! Indicates whether the latitude and longitude grid is fully contained
        logical,allocatable :: IsInDmArea(:,:)                    ! Indicates whether the latitude and longitude grid is in the refinement region

        isinply = 0.

        if(mode == 3)then    ! Triangular grid
            outputfile = trim(base_dir) // trim(EXPNME) // '/makegrid/patchtype/patchtype_sjx.nc4'      ! 最终输出文件
        else if(mode == 6)then     ! Polygonal mesh
            outputfile = trim(base_dir) // trim(EXPNME) // '/makegrid/patchtype/patchtype_lbx.nc4'
        end if

        print*, "Start reading unstructured grid data......"
        print*, ""

        write(nxpc, '(I3.3)')NXP

        if(refine == .false.)then
           flnm = trim(base_dir) // trim(EXPNME) // '/makegrid/gridfile/gridfile_NXP' // trim(nxpc) // '.nc4'
        else if(mode == 3)then
           flnm = trim(base_dir) // trim(EXPNME) // '/makegrid/result/gridfile_NXP' // trim(nxpc) // '_sjx.nc4'
        else if(mode == 6)then
           flnm = trim(base_dir) // trim(EXPNME) // '/makegrid/result/gridfile_NXP' // trim(nxpc) // '_lbx.nc4'
        end if

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

        if(mode == 3)then
          ustr_points = sjx_points
        else if(mode == 6)then
          ustr_points = lbx_points
        end if

        allocate(wp(lbx_points, 3))          ! Longitude, latitude, area
        allocate(mp(sjx_points, 3))
        allocate(ustr_id(ustr_points, 2))    ! The unstructured grid contains the number of latitude and longitude grids
        allocate(ngrwm(8, lbx_points))       ! Neighborhood array, the extra dimension represents adjacent points
        allocate(ngrmw(4, sjx_points))
        wp = 0.
        mp = 0.
        ustr_id = 0
        ngrmw = 0
        ngrwm = 0

        CALL CHECK(NF90_GET_VAR(ncid, varid(1), wp(:, 1)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(2), wp(:, 2)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(3), mp(:, 1)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(4), mp(:, 2)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(5), ngrwm(1:7, :)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(6), ngrmw(1:3, :)))
        CALL CHECK(NF90_CLOSE(ncid))

        ! Gets the number of adjacent w(m) points for each m(w) point
        call GetNgrNum(sjx_points, lbx_points, ngrmw, ngrwm)

        allocate(landtypes(nlons_source, nlats_source))
        landtypes = 0.

        print*,nlons_source, nlats_source
        if(lcs == "igbp")then
           flnm = trim(source_dir) // 'landtype_igbp_update.nc'
        else
           flnm = trim(source_dir) // 'landtype_usgs_update.nc'
        end if
        print*,flnm
        CALL CHECK(NF90_OPEN(flnm, nf90_nowrite, ncid))
        CALL CHECK(NF90_INQ_VARID(ncid, "landtype", varid(1)))
        CALL CHECK(NF90_GET_VAR(ncid, varid(1), landtypes))
        CALL CHECK(NF90_CLOSE(ncid))
        print*,"landtypes",minval(landtypes),maxval(landtypes)

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
        deallocate(landtypes)

        allocate(area_fine_gridcell(nlons_source, nlats_source))
        area_fine_gridcell(:, :) = 0.
        print*, "Start calculating the area of the latitude and longitude grid......"
        call cellarea(area_fine_gridcell, nlons_source, nlats_source)
        print*, "The latitude and longitude grid area is calculated"

        maxlat = maxval(mp(:, 2))
        minlat = minval(mp(:, 2))

        num_i = 0

        allocate(IsInDmArea(nlons_source,nlats_source))
        IsInDmArea = .false.
        CALL IsInDomainArea(IsInDmArea,lon_i,lat_i,dx,dy)

        !stop

        ! --------------------------------------------------------
        ! Obtain the number and proportion of latitude and longitude grids in unstructured grids
        ! 1. Calculate the array size
        ! 2. Allocate memory
        ! 3. Calculate inclusion relationships
        ! --------------------------------------------------------
        ! First calculate the array size
        if(mode == 3)then
          print*, "Start to calculate the triangle grid and latitude and longitude grid contain relation array size......"
!$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
!$OMP PRIVATE(i,j,k,l,sjx,isinply,ispart)
          do i = 1, nlons_source, 1
              num_i = num_i + 1
              do j = 1, nlats_source, 1                    ! Loop through the initial grid cell
                  if(IsInDmArea(i,j) == .false.)then
                     cycle
                  end if
                  ispart = .false.
                  if(seaorland(i, j) == 0)then
                      cycle
                  end if                           ! Only traverse the ocean grid
                  do k = 1, ustr_points, 1           ! Loop through the triangular grid
                      if(ispart == .true.)then
                          exit
                      end if
                      if(ngrmw(4, k) == 3)then       ! Determine if there are three recorded w points in the neighborhood of point m
                          sjx = 0.
                          do l = 1, 3, 1
                              sjx(l, 1:2) = wp(ngrmw(l, k), 1:2)  ! Record the vertex information of the triangle mesh
                          end do

                          isinply = IsInUstrGrid(sjx, lon_i(i), lat_i(j), 3, maxlat, minlat)

                          if(isinply == 1.)then      ! Complete inclusion
                              ispart = .true.
                          end if

                          if(isinply > 0.)then       ! Record the number of latitude and longitude grids contained in each unstructured grid
                              ustr_id(k, 1) = ustr_id(k, 1) + 1
                          end if
                      end if
                  end do
              end do
          end do
!$OMP END PARALLEL DO

        else if(mode == 6)then
          print*, "Start to calculate the size of the hexagon grid and latitude and longitude grid contain relation array......"
!$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1)&
!$OMP PRIVATE(i,j,k,l,dbx,isinply,ispart)
          do i = 1, nlons_source, 1
              num_i = num_i + 1
              !print*, num_i
              do j = 1, nlats_source, 1                    ! Loop through the initial grid cell
                  if(IsInDmArea(i,j) == .false.)then
                     cycle
                  end if
                  ispart = .false.
                  if(seaorland(i, j) == 0)then
                      cycle
                  end if                           ! Only traverse the ocean grid
                  do k = 1, ustr_points, 1           ! Loop through the polygonal mesh
                      if(ispart == .true.)then
                          exit
                      end if
                      if(ngrwm(8, k) > 4)then       ! Determine whether there are more than four recorded m points in the neighborhood of point w
                          dbx = 0.
                          do l = 1, ngrwm(8, k), 1
                              dbx(l, 1:2) = mp(ngrwm(l, k), 1:2)  ! Record the vertex information of polygon mesh
                          end do

                          isinply = IsInUstrGrid(dbx, lon_i(i), lat_i(j), ngrwm(8, k), maxlat, minlat)

                          if(isinply == 1.)then
                              ispart = .true.
                          end if

                          if(isinply > 0.)then            ! Record the number of latitude and longitude grids contained in each unstructured grid
                              ustr_id(k, 1) = ustr_id(k, 1) + 1
                          end if
                      end if
                  end do
              end do
          end do
!$OMP END PARALLEL DO
        end if
        print*, "The size of the containing relational array is calculated......"

        numpatch = INT(sum(ustr_id(:, 1))) * 2
        allocate(ustr_ii(numpatch, 4))
        ustr_ii = 0.

        ! ustr_id(:,1) records the position of the mp or wp array in the ustr_ii array
        ustr_id(1, 2) = 1

        do i = 2, ustr_points, 1
          ustr_id(i, 2) = ustr_id(i - 1, 2) + ustr_id(i - 1, 1) * 2
        end do

        ustr_id(:, 1) = 0
        num_i = 0

        if(mode == 3)then
          print*, "The inclusion relationship between the triangular mesh and the latitude and longitude mesh is calculated......"
          !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
          !$OMP PRIVATE(i,j,k,l,ispart,sjx,isinply,id)
          do i = 1, nlons_source, 1
              num_i = num_i + 1
              !print*, num_i
              do j = 1, nlats_source, 1
                  if(IsInDmArea(i,j) == .false.)then
                     cycle
                  end if
                  ispart = .false.
                  if(seaorland(i, j) == 0)then
                      cycle
                  end if
                  do k = 1, sjx_points, 1
                      if(ispart == .true.)then
                          exit
                      end if
                      if(ngrmw(4, k) == 3)then
                          sjx = 0.
                          do l = 1, 3, 1
                              sjx(l, 1:2) = wp(ngrmw(l, k), 1:2)
                          end do

                          isinply = IsInUstrGrid(sjx, lon_i(i), lat_i(j), 3, maxlat, minlat)

                          if(isinply == 1.)then
                              ustr_id(k, 1) = ustr_id(k, 1) + 1
                              mp(k, 3) = mp(k, 3) + area_fine_gridcell(i, j)     ! Unstructured grid area

                              id = ustr_id(k, 2) + ustr_id(k, 1) - 1
                              ustr_ii(id, 1) = i
                              ustr_ii(id, 2) = j
                              ustr_ii(id, 3) = isinply
                              ustr_ii(id, 4) = area_fine_gridcell(i, j)

                              ispart = .true.
                          else if(isinply > 0.)then
                              ustr_id(k, 1) = ustr_id(k, 1) + 1
                              mp(k, 3) = mp(k, 3) + area_fine_gridcell(i, j) * isinply

                              id = ustr_id(k, 2) + ustr_id(k, 1) - 1
                              ustr_ii(id, 1) = i
                              ustr_ii(id, 2) = j
                              ustr_ii(id, 3) = isinply
                              ustr_ii(id, 4) = area_fine_gridcell(i, j) * isinply

                          end if
                      end if
                  end do
              end do
          end do
          !$OMP END PARALLEL DO
        else if(mode == 6)then
          print*, "The inclusion relationship between the polygon mesh and the latitude and longitude mesh is calculated......"
          !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1)&
          !$OMP PRIVATE(i,j,k,l,ispart,dbx,isinply,id)
          do i = 1, nlons_source, 1
              num_i = num_i + 1
              !print*, num_i
              do j = 1, nlats_source, 1
                  if(IsInDmArea(i,j) == .false.)then
                     cycle
                  end if
                  ispart = .false.
                  if(seaorland(i, j) == 0)then
                      cycle
                  end if
                  do k = 1, lbx_points, 1
                      if(ispart == .true.)then
                          exit
                      end if
                      if(ngrwm(8, k) > 4)then
                          dbx = 0.
                          do l = 1, ngrwm(8, k), 1
                              dbx(l, 1:2) = mp(ngrwm(l, k), 1:2)
                          end do

                          isinply = IsInUstrGrid(dbx, lon_i(i), lat_i(j), ngrwm(8, k), maxlat, minlat)

                          if(isinply == 1.)then
                              ustr_id(k, 1) = ustr_id(k, 1) + 1
                              wp(k, 3) = wp(k, 3) + area_fine_gridcell(i, j)

                              id = ustr_id(k, 2) + ustr_id(k, 1) - 1
                              ustr_ii(id, 1) = i
                              ustr_ii(id, 2) = j
                              ustr_ii(id, 3) = isinply
                              ustr_ii(id, 4) = area_fine_gridcell(i, j)

                              ispart = .true.

                          else if(isinply > 0.)then
                              ustr_id(k, 1) = ustr_id(k, 1) + 1
                              wp(k, 3) = wp(k, 3) + area_fine_gridcell(i, j) * isinply

                              id = ustr_id(k, 2) + ustr_id(k, 1) - 1
                              ustr_ii(id, 1) = i
                              ustr_ii(id, 2) = j
                              ustr_ii(id, 3) = isinply
                              ustr_ii(id, 4) = area_fine_gridcell(i, j) * isinply

                          end if
                      end if
                  end do
              end do
          end do
          !$OMP END PARALLEL DO
        end if

        do i = 1, ustr_points, 1
          num_null = 0
          do j = ustr_id(i, 2), ustr_id(i, 2) + ustr_id(i, 1) - 1, 1
              if(ustr_ii(j, 3) == 0.)then
                  num_null = num_null + 1
              end if
          end do
          ustr_id(i, 1) = ustr_id(i, 1) - num_null
        end do

        allocate(ustr_id_new(ustr_points, 2))
        ustr_id_new = 0
        ustr_id_new(:, 1) = ustr_id(:, 1)
        ustr_id_new(1, 2) = 1

        do i = 2, ustr_points, 1
          ustr_id_new(i, 2) = ustr_id_new(i - 1, 2) + ustr_id_new(i - 1, 1)
        end do

        numpatch = INT(sum(ustr_id_new(:, 1)))
        allocate(ustr_ii_new(numpatch, 4))
        ustr_ii_new = 0

        do i = 1, ustr_points, 1
          do j = 1, 4, 1
              if(ustr_id_new(i, 1) /= 0)then
                  ustr_ii_new(ustr_id_new(i, 2):ustr_id_new(i, 2) + ustr_id_new(i, 1) - 1, j) = &
                          ustr_ii(ustr_id(i, 2):ustr_id(i, 2) + ustr_id_new(i, 1) - 1, j)
              end if
          end do
        end do

        if(mode == 3)then

          print*, "The triangular grid contains the array. The calculation is complete and the storage begins......"
          print*, ""
          flnm = trim(base_dir) // trim(EXPNME) // '/makegrid/contain'
          call SaveFile(mp, ustr_id_new, ustr_ii_new, ustr_points, numpatch, flnm)
          print*, "Unstructured grid storage complete"
          print*, ""

        else if(mode == 6)then

          print*, "The polygonal mesh contains the array. The calculation is complete and the storage begins......"
          print*, ""
          flnm = trim(base_dir) // trim(EXPNME) // '/makegrid/contain'
          call SaveFile(wp, ustr_id_new, ustr_ii_new, ustr_points, numpatch, flnm)
          print*, "Unstructured grid storage complete"
          print*, ""

        end if
        !end if

        !-----------------------------------------------------------
        ! Start calculating the patchID file required for mpi mode
        !-----------------------------------------------------------

        dx = nlons_source / nlons_dest
        dy = nlats_source / nlats_dest

        allocate(lon_e(nlons_dest))
        allocate(lon_w(nlons_dest))
        allocate(lat_n(nlats_dest))
        allocate(lat_s(nlats_dest))

        do i = 1, nlons_dest, 1
            lon_e(i) = -180. + 360. * i / nlons_dest
            lon_w(i) = -180. + 360. * (i - 1) / nlons_dest
        end do

        do i = 1, nlats_dest, 1
            lat_s(i) = 90. - 180. * i / nlats_dest
            lat_n(i) = 90. - 180. * (i - 1) / nlats_dest
        end do
      
        allocate(map(nlons_source,nlats_source,4))     ! Record the four borders, north, south, east
        map(1, :, 3) = 1
        map(nlons_source, :, 1) = nlons_dest
        map(:, 1, 2) = 1
        map(:, nlats_source, 4) = nlats_dest

        row = 1
        col = 1

        ! Calculate the boundaries of the target latitude and longitude grid
        do i = 1, nlons_dest, 1
            if(i * dx <= row)then
                map(row, :, 1) = i
            else
                map(row + 1, :, 3) = i
                map(row + 1, :, 1) = i
                row = row + 1
            end if
        end do

        do j = 1, nlats_dest, 1
            if(j * dy <= col)then
                map(:, col, 4) = j
            else
                map(:, col + 1, 4) = j
                map(:, col + 1, 2) = j
                col = col + 1
            end if
        end do

        allocate(patches_fraction(nlons_dest, nlats_dest))   ! Area ratio
        allocate(patchtypes(nlons_dest, nlats_dest))         ! Output file
        patches_fraction = 0.
        patchtypes = 0

        do i = 1, ustr_points, 1
            if(ustr_id(i, 1) == 0.)then
                cycle
            end if
            do j = 0, ustr_id(i, 1) - 1, 1
                row = int(ustr_ii(ustr_id(i, 2) + j, 1))
                col = int(ustr_ii(ustr_id(i, 2) + j, 2))
                if((row == 0) .or. (col == 0))then
                    cycle
                end if

                ! The ownership is determined by the size of the contain proportion
                if(patches_fraction(row, col) < ustr_ii(ustr_id(i, 2) + j, 3))then
                    patches_fraction(row, col) = ustr_ii(ustr_id(i, 2) + j, 3)

                    ! Assign values according to boundaries
                    do x = map(row, col, 3), map(row, col, 1), 1
                        do y = map(row, col, 2), map(row, col, 4), 1
                            patchtypes(x, y) = i
                        end do
                    end do

                end if
            end do
        end do

        varid = 0

        print*, outputfile
        CALL CHECK(NF90_CREATE(trim(outputfile), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, "nlon", nlons_dest, loDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "nlat", nlats_dest, laDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "elmindex", NF90_INT, (/ loDimID, laDimID /), varid(1)))
        CALL CHECK(NF90_DEF_VAR(ncID, "lon_e", NF90_DOUBLE, (/ loDimID /), varid(2)))
        CALL CHECK(NF90_DEF_VAR(ncID, "lon_w", NF90_DOUBLE, (/ loDimID /), varid(3)))
        CALL CHECK(NF90_DEF_VAR(ncID, "lat_n", NF90_DOUBLE, (/ laDimID /), varid(4)))
        CALL CHECK(NF90_DEF_VAR(ncID, "lat_s", NF90_DOUBLE, (/ laDimID /), varid(5)))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, varid(1), patchtypes(:, :)))
        CALL CHECK(NF90_PUT_VAR(ncID, varid(2), lon_e))
        CALL CHECK(NF90_PUT_VAR(ncID, varid(3), lon_w))
        CALL CHECK(NF90_PUT_VAR(ncID, varid(4), lat_n))
        CALL CHECK(NF90_PUT_VAR(ncID, varid(5), lat_s))
        CALL CHECK(NF90_CLOSE(ncID))

    END SUBROUTINE Get_Contain_PatchId

    SUBROUTINE CHECK(STATUS)
        INTEGER, intent (in) :: STATUS
        if  (STATUS .NE. NF90_NOERR) then
            print *, NF90_STRERROR(STATUS)
            stop 'stopped'
        endif
    END SUBROUTINE CHECK


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
                if((ngrmw(j, i) /= 1).and.(ngrmw(j, i) /= 0))then
                    flag = flag + 1
                end if
            end do  ! Point 1 m and point w are zero
            ngrmw(4, i) = flag
        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC,1) &
        !$OMP PRIVATE(i,j,flag)
        do i = 1, lbx_points, 1
            flag = 0
            do j = 1, 7, 1
                if((ngrwm(j, i) /= 1).and.(ngrwm(j, i) /= 0))then
                    flag = flag + 1
                end if
            end do
            ngrwm(8, i) = flag
        end do
        !$OMP END PARALLEL DO

    END SUBROUTINE GetNgrNum


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

        !$OMP PARALLEL DO NUM_THREADS(96) SCHEDULE(DYNAMIC,1)&
        !$OMP PRIVATE(i,j,dx,dy)
        do j = 1, nlats, 1
            do i = 1, nlons, 1
                if(lone(i)<lonw(i))then   ! west edge is more western than data line
                    ! The western edge is west of the dateline
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
                ! The arc length formula solves the area
            end do
        end do
        !$OMP END PARALLEL DO

        global = sum(area(:, :))

        ! Ensure that the total area of the grid cells is the same 
	! as the area of the grid defined by their edges
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

        integer :: ncID, idDimID, infoDimID, ncVarID, iunit
        integer, intent(in) :: num, num_ii
        character(len = 256) :: outputfile(6), na, flnm
        real(r8), dimension(num, 3), intent(in) :: ustr
        real(r8), dimension(num_ii, 4), intent(in) :: ustr_ii
        integer, dimension(num, 2), intent(in) :: ustr_id

      if(mode == 6)then
         outputfile(1) = trim(flnm)//"/wp.nc4"
         outputfile(2) = trim(flnm)//"/wp_ii.nc4"
         outputfile(3) = trim(flnm)//"/wp_id.nc4"
         outputfile(4) = trim(flnm)//"/wp.bin"
         outputfile(5) = trim(flnm)//"/wp_ii.bin"
         outputfile(6) = trim(flnm)//"/wp_id.bin"
         na = "lbx_points"
      else if(mode == 3)then
         outputfile(1) = trim(flnm)//"/mp.nc4"
         outputfile(2) = trim(flnm)//"/mp_ii.nc4"
         outputfile(3) = trim(flnm)//"/mp_id.nc4"
         outputfile(4) = trim(flnm)//"/mp.bin"
         outputfile(5) = trim(flnm)//"/mp_ii.bin"
         outputfile(6) = trim(flnm)//"/mp_id.bin"
         na = "sjx_points"
      end if

        print*, outputfile(1)
        CALL CHECK(NF90_CREATE(trim(outputfile(1)), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, na, num, idDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dima", 3, infoDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "info", NF90_FLOAT, (/ idDimID, infoDimID /), ncVarID))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID, ustr))
        CALL CHECK(NF90_CLOSE(ncID))

        print*, outputfile(2)
        CALL CHECK(NF90_CREATE(trim(outputfile(2)), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, na, num_ii, idDimID))
        CALL CHECK(NF90_DEF_DIM(ncID, "dima", 4, infoDimID))
        CALL CHECK(NF90_DEF_VAR(ncID, "info", NF90_FLOAT, (/ idDimID, infoDimID /), ncVarID))
        CALL CHECK(NF90_ENDDEF(ncID))
        CALL CHECK(NF90_PUT_VAR(ncID, ncVarID, ustr_ii))
        CALL CHECK(NF90_CLOSE(ncID))

        print*, outputfile(3)
        CALL CHECK(NF90_CREATE(trim(outputfile(3)), ior(nf90_clobber, nf90_netcdf4), ncID))
        CALL CHECK(NF90_DEF_DIM(ncID, na, num, idDimID))
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

        if(refine == .false.)then
           if(mode == 3)then
              call execute_command_line('cp '//trim(outputfile(1))//' '//trim(flnm)//'/initial/mp.nc4')
              call execute_command_line('cp '//trim(outputfile(2))//' '//trim(flnm)//'/initial/mp_ii.nc4')
              call execute_command_line('cp '//trim(outputfile(3))//' '//trim(flnm)//'/initial/mp_id.nc4')
              call execute_command_line('cp '//trim(outputfile(4))//' '//trim(flnm)//'/initial/mp.bin')
              call execute_command_line('cp '//trim(outputfile(5))//' '//trim(flnm)//'/initial/mp_ii.bin')
              call execute_command_line('cp '//trim(outputfile(6))//' '//trim(flnm)//'/initial/mp_id.bin')
           else if(mode == 6)then
              call execute_command_line('cp '//trim(outputfile(1))//' '//trim(flnm)//'/initial/wp.nc4')
              call execute_command_line('cp '//trim(outputfile(2))//' '//trim(flnm)//'/initial/wp_ii.nc4')
              call execute_command_line('cp '//trim(outputfile(3))//' '//trim(flnm)//'/initial/wp_id.nc4')
              call execute_command_line('cp '//trim(outputfile(4))//' '//trim(flnm)//'/initial/wp.bin')
              call execute_command_line('cp '//trim(outputfile(5))//' '//trim(flnm)//'/initial/wp_ii.bin')
              call execute_command_line('cp '//trim(outputfile(6))//' '//trim(flnm)//'/initial/wp_id.bin')
           end if
        end if

    END SUBROUTINE SaveFile


    REAL FUNCTION IsInUstrGrid(ustr, lon, lat, num_points, maxlat_m, minlat_m)

        implicit none

        integer :: inc(4), i, j, iscross_l(2), ispole
        integer, intent(in) :: num_points
        integer :: num_points_i, sjxorlbx_i, num_inter
        real(r8), intent(in) :: lon, lat, maxlat_m, minlat_m
        real(r8), dimension(7, 2), intent(in) :: ustr          ! Vertex of unstructured grid element
        real(r8), allocatable :: ustr_move(:, :)               ! Non-structural grid longitude
        real(r8), dimension(4, 2) :: point                     ! Vertex of the latitude and longitude grid
        real(r8), dimension(2) :: center_point                 ! Unstructured grid cell center point
        integer, allocatable :: iscross_g(:)                   ! Determine whether the unstructured grid line segment crosses the latitude and longitude grid
        real(r8), dimension(20, 2) :: interarea_points         ! Vertex of the area where two grids intersect
        real(r8), allocatable :: inter_points(:, :, :)         ! The intersection of two grids
        real(r8), allocatable :: area(:)   ! The area of a triangle consisting of two adjacent points of an 
					   !unstructured grid and the vertices of any latitude and longitude grid
        real(r8), dimension(5) :: area_i   ! The area of a triangle composed of two adjacent points of the latitude 
					   !and longitude grid and the vertices of any unstructured grid
        real(r8) :: dx, dy, maxlat, minlat, maxlon, minlon
        real(r8) :: tmpa(2),tmpb(2),tmpc(2),tmpd(3,2)          ! Intermediate variable

        inc = 0                 ! Determine the number of latitude and longitude grid vertices in an unstructured grid
        IsInUstrGrid = 0        ! The position relationship between latitude and longitude grid and unstructured grid is determined
        center_point = 0.       ! Center point of latitude and longitude grid
        ispole = 0              ! Determine whether an unstructured grid contains poles
        interarea_points = 0.   ! Two kinds of mesh overlap polygon vertices
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

        ! Ensure that the absolute latitude of the latitude grid does not exceed 90°
        do i = 1, 4, 1
            if(point(i, 2) > 90.)then
                point(i, 2) = 90.
            else if(point(i, 2) < -90.)then
                point(i, 2) = -90.
            end if
        end do

        ! Filter the grid by latitude
        if((point(3, 2) > maxlat).or.(point(1, 2) < minlat))then
            IsInUstrGrid = -1
            return
        end if

        !print*,ustr
        !print*,maxlat,minlat

        !------------------------------------------------------------------
        ! Handle the top and bottom mesh of the triangular grid
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
        ! Handle the top and bottom of the polygonal mesh
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

        !print*,""
        !print*,"point",point
        !print*,"ustr",ustr
        !print*,"ustr_move",ustr_move

        allocate(inter_points(num_points_i, 3, 2))
        allocate(area(num_points_i + 1))
        allocate(iscross_g(num_points_i))
        inter_points = 0
        iscross_g = 0
        area = 0

        !-------------------------------------------------------------------------------------
        ! Determine whether the two grids cross ±180° longitude
        !-------------------------------------------------------------------------------------
        if((point(1, 1) > 180.).or.(point(2, 1) < -180.))then
            iscross_l(1) = 1
        end if

        iscross_l(2) = IsCrossLine(ustr_move(:, 1), num_points_i)

        if(iscross_l(2) == 1)then
            iscross_l(1) = 1
        end if
        !print*,iscross_l

        !-------------------------------------------------------------------------------------
        ! Move the grid point longitude according to the above judgment
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

        !print*,ustr

        !-------------------------------------------------------------------------------------
        ! Filter the grid by longitude
        !-------------------------------------------------------------------------------------
        minlon = minval(ustr_move(1:num_points_i, 1))
        maxlon = maxval(ustr_move(1:num_points_i, 1))

        if((point(2, 1) > maxlon).or.(point(1, 1) < minlon))then
            IsInUstrGrid = -1
            return
        end if
         
        !-------------------------------------------------------------------------------------
        ! Start to determine the position relationship between the two grids and record key points
        !-------------------------------------------------------------------------------------
        ! First determine that two grids intersect
        do i = 1, num_points_i - 1, 1
            tmpa = ustr_move(i, 1:2)
            tmpb = ustr_move(i + 1, 1:2)
            tmpd = inter_points(i, 1:3, 1:2)

            iscross_g(i) = IsCrossGrid(point, tmpa, tmpb, tmpd)
        end do

        tmpa = ustr_move(1, 1:2)
        tmpb = ustr_move(num_points_i, 1:2)
        tmpd = inter_points(num_points_i, 1:3, 1:2)

        iscross_g(num_points_i) = IsCrossGrid(point, tmpa, tmpb, tmpd)

        ! Calculate the number of vertices in the unstructured grid
        ! If it is 4, it contains, otherwise it intersects

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
                end do

                tmpa = ustr_move(1, 1:2)
                tmpb = ustr_move(num_points_i, 1:2)
                tmpc = point(i, 1:2)
                area(num_points_i) = GetTriangleArea(tmpa, tmpb, tmpc)

                tmpa = ustr_move(1, 1:2)
                tmpb = ustr_move(num_points_i, 1:2)
                area(num_points_i + 1) = area(num_points_i + 1) + GetTriangleArea(tmpa, tmpb, center_point)

                if(abs(sum(area(1:num_points_i)) - area(num_points_i + 1)) < 0.00005)then
                    inc(i) = 1
                end if
            end if
        end do

      if  (no_caculate_fraction) then
         if((sum(inc) > 0).or.(sum(iscross_g) > 0))then
            IsInUstrGrid = 1.
            return
         else
            IsInUstrGrid = 0
            return
         end if
      endif

        if(sum(inc) == 4)then   ! If contain
            IsInUstrGrid = 1.
            return
        else if((sum(inc) == 0).and.(sum(iscross_g) == 0))then  ! It does not intersect if it is not included
            IsInUstrGrid = 0.
            return
        else
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
                end do

                tmpa = ustr_move(i, 1:2)
                tmpb = point(1, 1:2)
                tmpc = point(4, 1:2)
                area_i(4) = GetTriangleArea(tmpa, tmpb, tmpc)

                tmpb = point(1, 1:2)
                tmpc = point(4, 1:2)
                area_i(5) = area_i(5) + GetTriangleArea(center_point, tmpb, tmpc)
                if(abs(sum(area_i(1:4)) - area_i(5)) < 0.05)then
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

            if(num_inter == 0)then
               IsInUstrGrid = 0.
               return
            end if

            call SortPoints(interarea_points, num_inter)

            IsInUstrGrid = GetAreaPercent(interarea_points, num_inter, point)

        end if

    END FUNCTION IsInUstrGrid


    ! Determine whether the grid crosses the 189° and -180° longitude lines
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


    ! Adjust the longitude of the grid points
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


    ! Determine whether the unstructured grid line segment crosses the latitude and longitude grid
    integer function IsCrossGrid(point, a, b, inter_point)

        implicit none

        real(r8), dimension(2), intent(in) :: a, b     
        real(r8), dimension(2) :: x, y
        real(r8), dimension(4, 2), intent(in) :: point
        real(r8), dimension(3, 2), intent(out) :: inter_point
        real(r8) :: x1, x2, y1, y2, m, n, num

        !print*,""
        !print*,"point",point
        !print*,"ab",a,b

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
        !print*,m,n,y1,y2,x1,x2
        !print*,IsCrossGrid

    end function IsCrossGrid


    ! Calculate the area of the triangle (by latitude and longitude, not the actual area)
    REAL FUNCTION GetTriangleArea(a, b, c)

        implicit none

        real(r8), dimension(2), intent(in) :: a, b, c
        real(r8) :: aa, bb, cc, p

        GetTriangleArea = 0.

        aa = sqrt((c(1) - b(1)) * (c(1) - b(1)) + (c(2) - b(2)) * (c(2) - b(2)))
        bb = sqrt((c(1) - a(1)) * (c(1) - a(1)) + (c(2) - a(2)) * (c(2) - a(2)))
        cc = sqrt((a(1) - b(1)) * (a(1) - b(1)) + (a(2) - b(2)) * (a(2) - b(2)))

        p = (aa + bb + cc) / 2

        GetTriangleArea = sqrt(p * (p - aa) * (p - bb) * (p - cc))

    END FUNCTION GetTriangleArea


    ! Gets the proportion of unstructured grids that contain latitude and longitude grids
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
        end do

        tmpa = inter_point(1, 1:2)
        tmpb = inter_point(num, 1:2)
        inter_area = inter_area + GetTriangleArea(center_point, tmpa, tmpb)

        GetAreaPercent = inter_area / (abs((point(1, 1) - point(2, 1)) * (point(1, 2) - point(4, 2))))
        ! The area ratio is the overlapping triangle divided by the area of the latitude and longitude grid

        if(GetAreaPercent >= 1)then
            GetAreaPercent = 1.
        end if

    END FUNCTION GetAreaPercent


    ! Sort points into polygons
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
            !print*,points(i,:)
            sort_i(i) = i
        end do

        center_point = center_point / num

        !print*,"center_point",center_point

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

        !print*,"angle2",angle

        do i = 1, num, 1
            points(i, :) = points_i(sort_i(i), :) + center_point
        end do

    END SUBROUTINE SortPoints


    SUBROUTINE IsInDomainArea(IsInDmArea,lon_i,lat_i,dx,dy)

      implicit none

      integer :: i,j,n
      real(r8),intent(in) :: dx,dy
      !integer,intent(in) :: nlons,nlats
      real(r8),dimension(nlons_source),intent(in) :: lon_i
      real(r8),dimension(nlats_source),intent(in) :: lat_i
      !real(r8),intent(in) :: edge_e,edge_w,edges,edgen
      real(r8),dimension(nlons_source) :: lone,lonw
      real(r8),dimension(nlats_source) :: latn,lats
      logical,dimension(nlons_source,nlats_source),intent(out) :: IsInDmArea


      do i = 1,nlons_source,1
         lone(i) = lon_i(i) + dx / 2.
         lonw(i) = lon_i(i) - dx / 2.
      end do

      do j = 1,nlats_source,1
         latn(j) = lat_i(j) + dy / 2.
         lats(j) = lat_i(j) - dy / 2.
      end do

      IsInDmArea = .false.

      do n = 1,ndm_domain,1
         do i = 1,nlons_source,1
            do j = 1,nlats_source,1
               if((lone(i) < edgee(n)).and.(lonw(i) > edgew(n)).and.&
                        (latn(j) < edgen(n)).and.(lats(j) > edges(n)))then
                  IsInDmArea(i,j) = .true.
               end if
            end do
         end do
      end do

   END SUBROUTINE IsInDomainArea

END Module MOD_Get_Contain_Patch


