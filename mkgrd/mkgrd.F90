!DESCRIPTION
!===========
! This program is the unstructure mesh generation tool for land surface models (e.g.,CoLM).

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
!* :SUBROUTINE:"***" :  ***

!Initial Author
!--------------
! Hanwen Fan and Zhongwang Wei @ SYSU

!REVISION HISTORY
!----------------
! 2023.10.28  Hanwen Fan and Zhongwang Wei @ SYSU
! 2023.02.21  Zhongwang Wei @ SYSU
! 2021.12.02  Zhongwang Wei @ SYSU 
! 2020.10.01  Zhongwang Wei @ SYSU


!===============================================================================

! Portions of this software are copied or derived from the OLAM software
! package.  The following copyright notice pertains to OLAM and its derivatives,
! including this code:  

!----------------------------------------------------------------------------
! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University;
! Colorado State University Research Foundation ; ATMET, LLC

! This software is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option)
! any later version.

! This software is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.

! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA
! (http://www.gnu.org/licenses/gpl.html)
!----------------------------------------------------------------------------

!===============================================================================
program main
    use consts_coms
    use refine_vars
    use MOD_Get_Contain_Patch, only : Get_Contain_PatchId
    use MOD_GetContain  , only : Get_Contain
    use MOD_GetThreshold, only : GetThreshold
    use MOD_Refine_LBX  , only : refine_lbx
    use MOD_Refine_LBX_Step2, only : refine_lbx_step2
    use MOD_Refine_SJX  , only : refine_sjx
    use MOD_domain_judge ,only : domain_judge
    implicit none
    character(pathlen) :: nlfile = 'mkgrd.mnl'
    character(LEN=256) :: finfolist

    ! If run is sequential, default choice is to set io6 to standard output unit 6.

    io6 = 6

    ! Initialize HDF5 library

    !call h5open_f(hdferr)
    CALL getarg(1, nlfile)

    ! Read Fortran namelist
    call read_nl(nlfile)

    call execute_command_line('mkdir -p '//trim(base_dir)//'/'//trim(expnme))
    call execute_command_line('mkdir -p '//trim(base_dir)//'/'//trim(expnme)//"/makegrid"//"/")
    call execute_command_line('mkdir -p '//trim(base_dir)//'/'//trim(expnme)//"/makegrid/contain"//"/")
    call execute_command_line('mkdir -p '//trim(base_dir)//'/'//trim(expnme)//"/makegrid/contain/initial"//"/")
    call execute_command_line('mkdir -p '//trim(base_dir)//'/'//trim(expnme)//"/makegrid/gridfile"//"/")
    call execute_command_line('mkdir -p '//trim(base_dir)//'/'//trim(expnme)//"/makegrid/patchtype"//"/")
    call execute_command_line('mkdir -p '//trim(base_dir)//'/'//trim(expnme)//"/makegrid/result"//"/")
    call execute_command_line('mkdir -p '//trim(base_dir)//'/'//trim(expnme)//"/makegrid/threshold"//"/")
    call execute_command_line('mkdir -p '//trim(base_dir)//'/'//trim(expnme)//"/makegrid/tmp"//"/")

    ! Initialize, execute, and end olam run

    call init_consts()

    call gridinit()

    finfolist =trim(base_dir)//'/'//trim(expnme)//'/makegrid/result/'//'clmkg.infolist'
      OPEN(100,file=trim(finfolist),form='formatted')
      write(100,*) '&makegrid'
      write(100,*) 'EXPNAME  = ', expnme 
      write(100,*) 'NXP = ', nxp
      write(100,*) 'mode = ',mode
      write(100,*) 'nlons_source = ', nlons_source
      write(100,*) 'nlats_source = ', nlats_source
      write(100,*) 'nlons_dest = ', nlons_dest
      write(100,*) 'nlats_dest = ', nlats_dest
      write(100,*) 'refine = ', refine
      write(100,*) 'no_caculate_fraction = ', no_caculate_fraction
      write(100,*) 'ndm_domain = ', ndm_domain
      write(100,*) 'edgee = ', edgee(1:ndm_domain)
      write(100,*) 'edgew = ', edgew(1:ndm_domain)
      write(100,*) 'edges = ', edges(1:ndm_domain)
      write(100,*) 'edgen = ', edgen(1:ndm_domain)
      write(100,*) 'lcs = ', lcs
      write(100,*) 'maxlc = ', maxlc

      if (refine) then
         write(100,*) ''
         write(100,*) '&refine'
         write(100,*) 'ndm_refine = ', ndm_refine
         write(100,*) 'edgee_rf = ', edgee_rf(1:ndm_refine)
         write(100,*) 'edgew_rf = ', edgew_rf(1:ndm_refine)
         write(100,*) 'edges_rf = ', edges_rf(1:ndm_refine)
         write(100,*) 'edgen_rf = ', edgen_rf(1:ndm_refine)
         write(100,*) 'max_iter = ', max_iter
         write(100,*) 'max_sa_iter = ', max_sa_iter
         write(100,*) 'refine_num_landtypes = ', refine_num_landtypes
         write(100,*) 'refine_area_mainland = ', refine_area_mainland
         write(100,*) 'refine_lai_m = ', refine_lai_m
         write(100,*) 'refine_lai_s = ', refine_lai_s
         write(100,*) 'refine_slope_m = ', refine_slope_m
         write(100,*) 'refine_slope_s = ', refine_slope_s
         write(100,*) 'refine_k_s_m = ', refine_k_s_m
         write(100,*) 'refine_k_s_s = ', refine_k_s_s
         write(100,*) 'refine_k_solids_m = ', refine_k_solids_m
         write(100,*) 'refine_k_solids_s = ', refine_k_solids_s
         write(100,*) 'refine_tkdry_m = ', refine_tkdry_m
         write(100,*) 'refine_tkdry_s = ', refine_tkdry_s
         write(100,*) 'refine_tksatf_m = ', refine_tksatf_m
         write(100,*) 'refine_tksatf_s = ', refine_tksatf_s
         write(100,*) 'refine_tksatu_m = ', refine_tksatu_m
         write(100,*) 'refine_tksatu_s = ', refine_tksatu_s
         write(100,*) 'th_num_landtypes = ', th_num_landtypes
         write(100,*) 'th_area_mainland = ', th_area_mainland
         write(100,*) 'th_lai_m = ', th_lai_m
         write(100,*) 'th_lai_s = ', th_lai_s
         write(100,*) 'th_slope_m = ', th_slope_m
         write(100,*) 'th_slope_s = ', th_slope_s
         write(100,*) 'th_k_s_m = ', th_k_s_m
         write(100,*) 'th_k_s_s = ', th_k_s_s
         write(100,*) 'th_k_solids_m = ', th_k_solids_m
         write(100,*) 'th_k_solids_s = ', th_k_solids_s
         write(100,*) 'th_tkdry_m = ', th_tkdry_m
         write(100,*) 'th_tkdry_s = ', th_tkdry_s
         write(100,*) 'th_tksatf_m = ', th_tksatf_m
         write(100,*) 'th_tksatf_s = ', th_tksatf_s
         write(100,*) 'th_tksatu_m = ', th_tksatu_m
         write(100,*) 'th_tksatu_s = ', th_tksatu_s
      end if
      write(100,*) '/'
      CLOSE(100)

    write(io6, *) 'make grid run complete'

    if (refine) then

        write(io6, *) 'Determine whether there is intersection between the refined region and domain'
        call domain_judge()

        write(io6, *) 'make grid with refine mesh'
        call Get_Contain()
        call GetThreshold()

        if(mode == 6)then
            call refine_lbx()
            if(max_iter > 1)then
                call refine_lbx_step2()
            endif
        else if(mode == 3)then
            call refine_sjx()
        endif

        call Get_Contain_PatchId()

    else
        print*, 'make grid with basic mesh'
        call Get_Contain_PatchId()
    end if

    write(io6, '(A)')
    write(io6, '(A)') "make basic grid end"

end program main


subroutine read_nl(file)
    use consts_coms
    use refine_vars
    implicit none

    character(*), intent(in) :: file

    logical :: fexists
    namelist /mkgrd/ nl
    namelist /mkrefine/ rl
    ! OPEN THE NAMELIST FILE
    inquire(file = file, exist = fexists)
    print*, file
    if (.not. fexists) then
        write(*, *) "The namelist file " // trim(file) // " is missing."
        stop "Stopping model run."
    endif
    open(10, status = 'OLD', file = file)
    ! READ GRID POINT, MODEL OPTIONS, AND PLOTTING INFORMATION FROM THE NAMELIST
    REWIND(10)
    !print*,nl
    read(10, nml = mkgrd)
    close(10)
    write(*, nml = mkgrd)

    expnme = nl%expnme
    !----------------------------------------------------------

    ! Variables in the following section either must not be changed on a history
    ! restart or changing them would be irrelevant.  Thus, they are only copied
    ! from the namelist if a history file is not being read.
    nxp = nl%nxp
    openmp = nl%openmp
    refine = nl%refine
    nlons_source = nl%nlons_source
    nlats_source = nl%nlats_source
    nlons_dest = nl%nlons_dest
    nlats_dest = nl%nlats_dest
    no_caculate_fraction = nl%no_caculate_fraction
    mode = nl%mode
    base_dir  = NL%base_dir
    source_dir = NL%source_dir
    ndm_domain = NL%ndm_domain
    edgee(1:ndm_domain) = nl%edgee(1:ndm_domain)
    edgew(1:ndm_domain) = nl%edgew(1:ndm_domain)
    edges(1:ndm_domain) = nl%edges(1:ndm_domain)
    edgen(1:ndm_domain) = nl%edgen(1:ndm_domain)
    lcs = nl%lcs
    maxlc = nl%maxlc
    if (refine) then
        open(10, status = 'OLD', file = file)
        ! READ GRID POINT, MODEL OPTIONS, AND PLOTTING INFORMATION FROM THE NAMELIST
        REWIND(10)
        !print*,nl
        read(10, nml = mkrefine)
        close(10)
        write(*, nml = mkrefine)

        ndm_refine            = RL%ndm_refine
        edgee_rf(1:ndm_refine)   = RL%edgee_rf(1:ndm_refine)
        edgew_rf(1:ndm_refine)   = RL%edgew_rf(1:ndm_refine)
        edges_rf(1:ndm_refine)   = RL%edges_rf(1:ndm_refine)
        edgen_rf(1:ndm_refine)   = RL%edgen_rf(1:ndm_refine)

        max_iter              = rl%max_iter
        max_sa_iter           = rl%max_sa_iter
        th_num_landtypes      = rl%th_num_landtypes

        refine_num_landtypes  =  rl%refine_num_landtypes
        refine_area_mainland  =  rl%refine_area_mainland
        refine_lai_m          =  rl%refine_lai_m
        refine_lai_s          =  rl%refine_lai_s
        refine_slope_m        =  rl%refine_slope_m
        refine_slope_s        =  rl%refine_slope_s
        refine_k_s_m          =  rl%refine_k_s_m
        refine_k_s_s          =  rl%refine_k_s_s
        refine_k_solids_m     =  rl%refine_k_solids_m
        refine_k_solids_s     =  rl%refine_k_solids_s
        refine_tkdry_m        =  rl%refine_tkdry_m
        refine_tkdry_s        =  rl%refine_tkdry_s
        refine_tksatf_m       =  rl%refine_tksatf_m
        refine_tksatf_s       =  rl%refine_tksatf_s
        refine_tksatu_m       =  rl%refine_tksatu_m
        refine_tksatu_s       =  rl%refine_tksatu_s

        th_area_mainland      = rl%th_area_mainland
        th_lai_m              = rl%th_lai_m
        th_lai_s              = rl%th_lai_s
        th_slope_m            = rl%th_slope_m
        th_slope_s            = rl%th_slope_s
        th_k_s_m              = rl%th_k_s_m
        th_k_s_s              = rl%th_k_s_s
        th_k_solids_m         = rl%th_k_solids_m
        th_k_solids_s         = rl%th_k_solids_s
        th_tkdry_m            = rl%th_tkdry_m
        th_tkdry_s            = rl%th_tkdry_s
        th_tksatf_m           = rl%th_tksatf_m
        th_tksatf_s           = rl%th_tksatf_s
        th_tksatu_m           = rl%th_tksatu_m
        th_tksatu_s           = rl%th_tksatu_s

    end if

end subroutine read_nl


subroutine init_consts()
    use consts_coms
    implicit none


    ! Standard (Earth) values

    erad = 6371.22e3          ! Earth radius [m]
    omega = 7.29212e-5         ! Earth rotational angular velocity
    xscale = 1.0

    ! Secondary values

    erad2 = erad * 2.
    erad4 = erad * 4.
    eradsq = erad * erad
    erador5 = erad / sqrt(5.)
    eradi = 1.0 / erad
    erad2sq = erad2 * erad2
    dlat = erad * pio180
    omega2 = omega * 2.

end subroutine init_consts


!===============================================================================
subroutine gridinit()

    use consts_coms, only : io6, nxp

    use mem_ijtabs, only : fill_jtabs

    use mem_delaunay, only : copy_tri_grid, copyback_tri_grid, nmd, nud, nwd

    use mem_grid, only : nma, nua, nva, nwa, alloc_grid1, alloc_grid2

    implicit none


    ! Horizontal grid setup


    ! Now generate global grid

    write(io6, '(/,a)') 'gridinit calling icosahedron'
    call icosahedron(nxp)  ! global spherical domain; calls 2 allocs
    write(io6, '(/,a)') 'gridinit after icosahedron'
    write(io6, '(a,i0)')    ' nmd = ', nmd
    write(io6, '(a,i0)')    ' nud = ', nud
    write(io6, '(a,i0)')    ' nwd = ', nwd


    ! Store a temporaty copy of the full Delaunay mesh
    ! to be used later to construct the surface grid
    call copy_tri_grid()

    call voronoi()
    call pcvt()
    write(io6, '(/,a)') 'gridinit after voronoi'
    write(io6, '(a,i8)')   ' nma = ', nma
    write(io6, '(a,i8)')   ' nua = ', nua
    write(io6, '(a,i8)')   ' nwa = ', nwa

    ! Allocate remaining GRID FOOTPRINT arrays for full domain

    write(io6, '(/,a)') 'gridinit calling alloc_grid1 for full domain'

    call alloc_grid1(nma, nva, nwa)

    ! Initialize dtlm, dtsm, ndtrat, and nacoust,
    ! and compute the timestep schedule for all grid operations.

    !write(io6,'(/,a)') 'gridinit calling modsched'

    !call modsched()

    write(io6, '(/,a)') 'gridinit calling fill_jtabs'
    print*, nma, nva, nwa

    call fill_jtabs(nma, nva, nwa, 0)

    ! Fill remaining GRID FOOTPRINT geometry for full domain

    write(io6, '(/,a)') 'gridinit calling grid_geometry'

    call grid_geometry_hex()

    call copyback_tri_grid()



    !write(io6,'(/,a)') 'gridinit calling makesfc3'

    !call makesfc3()

    ! Allocate remaining unstructured grid geometry arrays

    write(io6, '(/,a)') 'gridinit calling alloc_grid2'

    call alloc_grid2(nma, nva, nwa)

    ! Set up control volumes in atmospheric grid

    !write(io6,'(/,a)') 'gridinit calling ctrlvols'

    !call ctrlvols_hex()



    ! Write GRIDFILE and SFCGRILE

    write(io6, '(/,a)') 'gridinit calling gridfile_write'
    call gridfile_write()

    !write(io6,'(/,a)') 'calling sfcgile_write'
    !call sfcgfile_write()

    write(io6, '(/,a)') 'gridinit completed'

end subroutine gridinit


!===============================================================================

subroutine gridfile_write()
    use netcdf
    use MOD_Get_distance
    use consts_coms, only : pathlen
    use consts_coms, only : io6,base_dir,EXPNME,NXP,mode,refine
    use mem_ijtabs, only : mloops, itab_m, itab_w
    use mem_grid, only : nma, nwa, glatw, glonw, glatm, glonm

    implicit none

    ! This routine writes the grid variables to the gridfile.
    integer :: im, iw

    ! Scratch arrays for copying output

    real(r8), allocatable :: mp(:,:),wp(:,:)
    integer, allocatable :: ngrmw(:,:),ngrwm(:,:)
    real(r8),allocatable :: dismm(:,:),disww(:,:)

    character(pathlen) :: flnm,nxpc
    integer :: lbx_points, sjx_points, ncvarid(8)
    integer :: spdimid, lpdimid, sedimid, thdimid, ncid

    allocate (ngrmw(3, nma)) ; ngrmw = 0
    allocate (ngrwm(7, nwa)) ; ngrwm = 0
    do im = 1, nma
       ngrmw(1:3, im) = itab_m(im)%iw(1:3)
    enddo
    do iw = 1, nwa
       ngrwm(1:7, iw) = itab_w(iw)%im(1:7)
    enddo

    allocate(mp(nma,2))
    allocate(wp(nwa,2))
    mp(:,1) = GLONM
    mp(:,2) = GLATM
    wp(:,1) = GLONW
    wp(:,2) = GLATW
    
    allocate(dismm(nwa,7)); dismm = 0.
    allocate(disww(nma,3)); disww = 0.
    Call Get_Dis(mp,wp,ngrmw,ngrwm,dismm,disww,nma,nwa)

    deallocate(mp)
    deallocate(wp)
    ! Execute the command
    call execute_command_line('mkdir -p '//trim(base_dir)//'/'//trim(expnme))
    write(nxpc, '(I3.3)')NXP
    ! Open gridfile
    flnm = trim(base_dir) // trim(EXPNME) // '/makegrid/gridfile/gridfile_NXP' // trim(nxpc) // '.nc4'

    write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(io6, *) 'grid_write: opening file:', flnm
    write(io6, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    sjx_points = nma
    lbx_points = nwa
    CALL CHECK(NF90_CREATE(trim(flnm), ior(nf90_clobber, nf90_netcdf4), ncID))
    CALL CHECK(NF90_DEF_DIM(ncID, "sjx_points", sjx_points, spDimID))
    CALL CHECK(NF90_DEF_DIM(ncID, "lbx_points", lbx_points, lpDimID))
    CALL CHECK(NF90_DEF_DIM(ncID, "dimb", 3, thDimID))
    CALL CHECK(NF90_DEF_DIM(ncID, "dimc", 7, seDimID))
    CALL CHECK(NF90_DEF_VAR(ncID, "GLONW", NF90_FLOAT, (/ lpDimID /), ncVarID(1)))
    CALL CHECK(NF90_DEF_VAR(ncID, "GLATW", NF90_FLOAT, (/ lpDimID /), ncVarID(2)))
    CALL CHECK(NF90_DEF_VAR(ncID, "GLONM", NF90_FLOAT, (/ spDimID /), ncVarID(3)))
    CALL CHECK(NF90_DEF_VAR(ncID, "GLATM", NF90_FLOAT, (/ spDimID /), ncVarID(4)))
    CALL CHECK(NF90_DEF_VAR(ncID, "itab_m%iw", NF90_INT, (/ thDimID, spDimID /), ncVarID(5)))
    CALL CHECK(NF90_DEF_VAR(ncID, "itab_w%im", NF90_INT, (/ seDimID, lpDimID /), ncVarID(6)))
    CALL CHECK(NF90_DEF_VAR(ncID, "dis_w%iw", NF90_FLOAT, (/ spDimID, thDimID /), ncVarID(7)))
    CALL CHECK(NF90_DEF_VAR(ncID, "dis_m%im", NF90_FLOAT, (/ lpDimID, seDimID /), ncVarID(8)))
    CALL CHECK(NF90_ENDDEF(ncID))

    CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(1), GLONW))
    CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(2), GLATW))
    CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(3), GLONM))
    CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(4), GLATM))

    !allocate (iscr2(3, nma)) ; iscr2 = 0
    !do im = 1, nma
    !    iscr2(1:3, im) = itab_m(im)%iw(1:3)
    !enddo
    CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(5), ngrmw))

    !allocate (iscr2(7, nwa))

    !do iw = 1, nwa
    !    iscr2(1:7, iw) = itab_w(iw)%im(1:7)
    !enddo
    CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(6), ngrwm))
    CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(7), disww))
    CALL CHECK(NF90_PUT_VAR(ncID, ncVarID(8), dismm))
    CALL CHECK(NF90_CLOSE(ncID))
    deallocate(ngrmw)
    deallocate(ngrwm)
    deallocate(dismm)
    deallocate(disww)

    if(refine == .false.)then
       if(mode == 3)then
          call execute_command_line('cp '//trim(flnm)//' '//trim(base_dir)//trim(EXPNME)//'/makegrid/result/gridfile_NXP'//trim(nxpc)//'_sjx.nc4')
       else if(mode == 6)then
          call execute_command_line('cp '//trim(flnm)//' '//trim(base_dir)//trim(EXPNME)//'/makegrid/result/gridfile_NXP'//trim(nxpc)//'_lbx.nc4')
       end if
    end if

end subroutine gridfile_write
subroutine CHECK(STATUS)
    use netcdf
    INTEGER, intent (in) :: STATUS
    if  (STATUS .NE. NF90_NOERR) then 
        print *, NF90_STRERROR(STATUS)
        stop 'stopped'
    endif
end subroutine CHECK
subroutine voronoi()

    use mem_ijtabs, only : mloops, itab_m, itab_v, itab_w, alloc_itabs

    use mem_delaunay, only : itab_md, itab_ud, itab_wd, &
            xemd, yemd, zemd, nmd, nud, nwd

    use mem_grid, only : nma, nua, nva, nwa, mma, mua, mva, mwa, &
            xem, yem, zem, xew, yew, zew, &
            alloc_xyzem, alloc_xyzew

    use consts_coms, only : pi2, erad

    implicit none

    integer :: im1, im2
    integer :: iw1, iw2, iw3, im, iw
    integer :: iwd, iv, iud, iud1, iud2, imd, npoly, j, j1
    real :: expansion

    ! Interchange grid dimensions

    nma = nwd
    nua = nud
    nva = nud
    nwa = nmd

    mma = nma
    mua = nua
    mva = nva
    mwa = nwa

    ! Allocate Voronoi set of itabs

    call alloc_itabs(nma, nva, nwa, 0)

    ! Allocate XEW,YEW,ZEW arrays, and fill their values from XEMD,YEMD,ZEMD, which
    ! still have the OLD nmad dimension which is the NEW nwa dimension

    call move_alloc(xemd, xew)
    call move_alloc(yemd, yew)
    call move_alloc(zemd, zew)

    ! Allocate XEM,YEM,ZEM to NEW nma dimension

    call alloc_xyzem(nma)

    ! Since XEM,YEM,ZEM have just been re-allocated, initialize their values to be
    ! barycenters of Delaunay triangles whose vertices are at XEW,YEW,ZEW

    do iwd = 2, nwd
        im = iwd

        ! Indices of 3 M points surrounding WD point

        if (any(itab_wd(iwd)%im(1:3) < 2)) cycle

        iw1 = itab_wd(iwd)%im(1)
        iw2 = itab_wd(iwd)%im(2)
        iw3 = itab_wd(iwd)%im(3)

        xem(im) = (xew(iw1) + xew(iw2) + xew(iw3)) / 3.
        yem(im) = (yew(iw1) + yew(iw2) + yew(iw3)) / 3.
        zem(im) = (zew(iw1) + zew(iw2) + zew(iw3)) / 3.

        ! push M point coordinates out to earth radius
        expansion = erad / sqrt(xem(im) ** 2  &
                + yem(im) ** 2  &
                + zem(im) ** 2)

        xem(im) = xem(im) * expansion
        yem(im) = yem(im) * expansion
        zem(im) = zem(im) * expansion

    enddo

    ! Loop over V points

    do iv = 2, nva
        iud = iv

        itab_v(iv)%loop(1:mloops) = itab_ud(iud)%loop(1:mloops)

        itab_v(iv)%ivp = itab_ud(iud)%iup
        itab_v(iv)%ivglobe = iv
        itab_v(iv)%mrlv = itab_ud(iud)%mrlu

        itab_v(iv)%im(1:6) = itab_ud(iud)%iw(1:6)
        itab_v(iv)%iw(1:2) = itab_ud(iud)%im(1:2)

        itab_v(iv)%iv(1:4) = itab_ud(iud)%iu(1:4)
        ! itab_v(iv)%iv(1:12) = itab_ud(iud)%iu(1:12)

        ! For periodic Cartesian hex domain, compute coordinates for outer M points

        im1 = itab_v(iv)%im(1)
        im2 = itab_v(iv)%im(2)
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        if (itab_wd(im1)%npoly < 3) then ! itab_m(im1)%npoly not filled yet
            xem(im1) = xew(iw1) + xew(iw2) - xem(im2)
            yem(im1) = yew(iw1) + yew(iw2) - yem(im2)
            zem(im1) = 0.
        elseif (itab_wd(im2)%npoly < 3) then ! itab_m(im2)%npoly not filled yet
            xem(im2) = xew(iw1) + xew(iw2) - xem(im1)
            yem(im2) = yew(iw1) + yew(iw2) - yem(im1)
            zem(im2) = 0.
        endif

        ! Extract information from IMD1 neighbor

        imd = itab_ud(iud)%im(1)
        npoly = itab_md(imd)%npoly

        do j = 1, npoly
            j1 = j + 1
            if (j == npoly) j1 = 1

            iud1 = itab_md(imd)%iu(j)
            iud2 = itab_md(imd)%iu(j1)

            ! IW(3) and IW(4) neighbors of IV

            if (iud2 == iv) then
                iw1 = itab_ud(iud1)%im(1)
                iw2 = itab_ud(iud1)%im(2)

                if (iw1 == imd) then
                    itab_v(iv)%iw(3) = iw2
                else
                    itab_v(iv)%iw(3) = iw1
                endif
            endif

            if (iud1 == iv) then
                iw1 = itab_ud(iud2)%im(1)
                iw2 = itab_ud(iud2)%im(2)

                if (iw1 == imd) then
                    itab_v(iv)%iw(4) = iw2
                else
                    itab_v(iv)%iw(4) = iw1
                endif
            endif

        enddo

    enddo

    ! Loop over W points

    do iw = 2, nwa
        imd = iw

        itab_w(iw)%loop(1:mloops) = itab_md(imd)%loop(1:mloops)

        itab_w(iw)%iwp = iw

        itab_w(iw)%npoly = itab_md(imd)%npoly
        itab_w(iw)%iwglobe = iw

        itab_w(iw)%mrlw = itab_md(imd)%mrlm
        itab_w(iw)%mrlw_orig = itab_md(imd)%mrlm_orig
        itab_w(iw)%ngr = itab_md(imd)%ngr

        npoly = itab_w(iw)%npoly

        ! Loop over IM/IV neighbors of IW

        do j = 1, itab_w(iw)%npoly
            im = itab_md(imd)%iw(j)
            iwd = im
            iv = itab_md(imd)%iu(j)

            iw1 = itab_v(iv)%iw(1)
            iw2 = itab_v(iv)%iw(2)

            itab_w(iw)%im(j) = im
            itab_w(iw)%iv(j) = iv

            if (iw1 == iw) then
                itab_w(iw)%iw(j) = iw2
                itab_w(iw)%dirv(j) = -1.
            else
                itab_w(iw)%iw(j) = iw1
                itab_w(iw)%dirv(j) = 1.
            endif

        enddo

    enddo

    ! Loop over M points

    do im = 2, nma
        iwd = im

        itab_m(im)%loop(1:mloops) = itab_wd(iwd)%loop(1:mloops)

        itab_m(im)%imp = im

        itab_m(im)%npoly = itab_wd(iwd)%npoly
        itab_m(im)%imglobe = im

        itab_m(im)%mrlm = itab_wd(iwd)%mrlw
        itab_m(im)%ngr = itab_wd(iwd)%ngr

        itab_m(im)%mrlm_orig = itab_wd(iwd)%mrlw_orig
        itab_m(im)%mrow = itab_wd(iwd)%mrow

        itab_m(im)%iv(1:3) = itab_wd(iwd)%iu(1:3)
        itab_m(im)%iw(1:3) = itab_wd(iwd)%im(1:3)
    enddo

    deallocate(itab_md, itab_ud, itab_wd)

end subroutine voronoi

!===============================================================================

subroutine pcvt()

    ! Iterative procedure for defining centroidal voronoi cells

    use mem_ijtabs, only : itab_m
    use mem_grid, only : nma, xem, yem, zem, xew, yew, zew
    use consts_coms, only : erad, eradi

    implicit none

    integer :: im, iw1, iw2, iw3

    real :: raxis, raxisi, expansion
    real :: sinwlat, coswlat, sinwlon, coswlon
    real :: dxe, dye, dze
    real :: xebc, yebc, zebc
    real :: x1, x2, x3, y1, y2, y3
    real :: dx12, dx13, dx23
    real :: s1, s2, s3
    real :: xcc, ycc

    ! Compute XEM,YEM,ZEM location as circumcentric coordinates of 3 W points.
    ! This establishes W cell as voronoi.

    ! Loop over all M points

    !$omp parallel
    !$omp do private(iw1,iw2,iw3,xebc,yebc,zebc,raxis,raxisi, &
    !$omp            sinwlat,coswlat,sinwlon,coswlon,dxe,dye,dze,x1,y1, &
    !$omp            x2,y2,x3,y3,dx12,dx13,dx23,s1,s2,s3,ycc,xcc)
    do im = 2, nma

        ! Indices of 3 W points surrounding M point

        if (any(itab_m(im)%iw(1:3) < 2)) cycle

        iw1 = itab_m(im)%iw(1)
        iw2 = itab_m(im)%iw(2)
        iw3 = itab_m(im)%iw(3)

        ! These were initialized to be the barycenter of each triangle

        xebc = xem(im)
        yebc = yem(im)
        zebc = zem(im)


        ! For global domain, transform from sphere to PS plane

        raxis = sqrt(xebc ** 2 + yebc ** 2)
        raxisi = 1.0 / raxis

        sinwlat = zebc * eradi
        coswlat = raxis * eradi

        sinwlon = yebc * raxisi
        coswlon = xebc * raxisi

        ! Transform 3 W points to PS coordinates

        dxe = xew(iw1) - xebc
        dye = yew(iw1) - yebc
        dze = zew(iw1) - zebc
        call de_ps(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, x1, y1)

        dxe = xew(iw2) - xebc
        dye = yew(iw2) - yebc
        dze = zew(iw2) - zebc
        call de_ps(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, x2, y2)

        dxe = xew(iw3) - xebc
        dye = yew(iw3) - yebc
        dze = zew(iw3) - zebc
        call de_ps(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, x3, y3)

        ! Compute intermediate quanties

        dx12 = x2 - x1
        dx13 = x3 - x1
        dx23 = x3 - x2

        s1 = x1**2 + y1**2
        s2 = x2**2 + y2**2
        s3 = x3**2 + y3**2

        ! Algebraic solution for circumcenter Y coordinate

        ycc = .5 * (dx13 * s2 - dx12 * s3 - dx23 * s1) &
                / (dx13 * y2 - dx12 * y3 - dx23 * y1)

        ! Algebraic solution for circumcenter X coordinate

        if (abs(dx12) > abs(dx13)) then
            xcc = (s2 - s1 - ycc * 2. * (y2 - y1)) / (2. * dx12)
        else
            xcc = (s3 - s1 - ycc * 2. * (y3 - y1)) / (2. * dx13)
        endif

        ! For global domain, transform circumcenter from PS to earth coordinates

        call ps_de(dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, xcc, ycc)

        xem(im) = dxe + xebc
        yem(im) = dye + yebc
        zem(im) = dze + zebc

    enddo
    !$omp end do nowait

    ! Adjust each M point to the Earth's radius for global domain



    !$omp do private(expansion)
    do im = 2, nma

        expansion = erad / sqrt(xem(im) ** 2 &
                + yem(im) ** 2 &
                + zem(im) ** 2)

        xem(im) = xem(im) * expansion
        yem(im) = yem(im) * expansion
        zem(im) = zem(im) * expansion

    enddo
    !$omp end do nowait


    !$omp end parallel

end subroutine pcvt

!===============================================================================

subroutine grid_geometry_hex()

    use mem_ijtabs, only : itab_m, itab_v, itab_w
    use mem_grid, only : nma, nva, nwa, xev, yev, zev, xem, yem, zem, &
            xew, yew, zew, unx, uny, unz, wnx, wny, wnz, &
            vnx, vny, vnz, glonw, glatw, dnu, dniu, dnv, dniv, arw0, &
            arm0, glonm, glatm, glatv, glonv
    use consts_coms, only : erad, piu180, pio2
    use consts_coms, only : r8

    implicit none

    integer :: im, iv, iw
    integer :: im1, im2
    integer :: iw1, iw2
    integer :: j, npoly
    real :: expansion
    real :: raxis
    real :: dvm1, dvm2
    integer :: j1, j2
    integer :: npoly1, npoly2, np
    real :: xv, yv, frac, alpha
    real :: xw1, xw2, yw1, yw2
    integer :: iwp, ivp, imp
    logical :: dops
    real :: quarter_kite(2, nva)
    integer :: lwork
    integer :: info

    real(r8) :: b(7), fo(7), vnx_ps(7), vny_ps(7), vnz_ps(7), vrot_x(7), vrot_y(7)
    real(r8), allocatable :: work(:)
    real(r8), allocatable :: a(:, :)
    real(r8) :: wsize(1), vdotw, vmag, fact

    ! Loop over all M points

    !$omp parallel

    !$omp do private(raxis)
    do im = 2, nma

        ! Latitude and longitude at M points
        raxis = sqrt(xem(im) ** 2 + yem(im) ** 2)  ! dist from earth axis
        glatm(im) = atan2(zem(im), raxis) * piu180
        glonm(im) = atan2(yem(im), xem(im)) * piu180

        ! Fill global index (replaced later if this run is parallel)
        itab_m(im)%imglobe = im

    enddo
    !$omp end do

    ! Loop over all V points

    !$omp do private(im1,im2,iw1,iw2,expansion,raxis,dvm1,dvm2,frac)
    do iv = 2, nva

        ! Fill global index (replaced later if this run is parallel)

        itab_v(iv)%ivglobe = iv

        ! M-point indices of two end points of V segment

        im1 = itab_v(iv)%im(1)
        im2 = itab_v(iv)%im(2)

        ! W-point indices on either side of V segment

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        ! V point is midway between W points of Voronoi cells

        xev(iv) = .5 * (xew(iw1) + xew(iw2))
        yev(iv) = .5 * (yew(iw1) + yew(iw2))
        zev(iv) = .5 * (zew(iw1) + zew(iw2))

        !  push V point coordinates out to earth radius
        expansion = erad / sqrt(xev(iv) ** 2 &
                + yev(iv) ** 2 &
                + zev(iv) ** 2)

        xev(iv) = xev(iv) * expansion
        yev(iv) = yev(iv) * expansion
        zev(iv) = zev(iv) * expansion


        ! Latitude and longitude at V point

        raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis
        glatv(iv) = atan2(zev(iv), raxis) * piu180
        glonv(iv) = atan2(yev(iv), xev(iv)) * piu180


        ! Normal distance across U face
        ! Unit vector components of U face
        !x Convert normal distance across U face to geodesic arc length

        dnu(iv) = sqrt((xem(im1) - xem(im2))**2 &
                + (yem(im1) - yem(im2))**2 &
                + (zem(im1) - zem(im2))**2)

        unx(iv) = (xem(im2) - xem(im1)) / dnu(iv)
        uny(iv) = (yem(im2) - yem(im1)) / dnu(iv)
        unz(iv) = (zem(im2) - zem(im1)) / dnu(iv)

        !x   dnu(iv) = erad2 * asin(dnu(iv) / erad2)
        dniu(iv) = 1. / dnu(iv)

        ! Normal distance across V face
        ! Unit vector components of V face
        !x Convert normal distance across V face to geodesic arc length

        dnv(iv) = sqrt((xew(iw1) - xew(iw2))**2 &
                + (yew(iw1) - yew(iw2))**2 &
                + (zew(iw1) - zew(iw2))**2)

        vnx(iv) = (xew(iw2) - xew(iw1)) / dnv(iv)
        vny(iv) = (yew(iw2) - yew(iw1)) / dnv(iv)
        vnz(iv) = (zew(iw2) - zew(iw1)) / dnv(iv)

        !x   dnv(iv) = erad2 * asin(dnv(iv) / erad2)
        dniv(iv) = 1. / dnv(iv)

        ! Skip this U point if iw1 < 2 or iw2 < 2

        if (iw1 < 2 .or. iw2 < 2) cycle

        itab_v(iv)%mrlv = max(itab_w(iw1)%mrlw, itab_w(iw2)%mrlw)

        ! Compute IM1 and IM2 values of quarter kite area,
        ! and add to ARM0 and ARW0 arrays

        dvm1 = sqrt((xev(iv) - xem(im1))**2 &
                + (yev(iv) - yem(im1))**2 &
                + (zev(iv) - zem(im1))**2)

        dvm2 = sqrt((xev(iv) - xem(im2))**2 &
                + (yev(iv) - yem(im2))**2 &
                + (zev(iv) - zem(im2))**2)

        ! Fractional distance along V edge where intersection with U edge is located

        frac = dvm1 * dniu(iv)

        if (im1 > 1 .and. im2 > 1 .and. (frac < .0001 .or. frac > .9999)) then
            print*, 'Non-intersecting U-V edges detected in grid geometry'
            print*, 'FRAC  = ', frac
            print*, 'IW1 = ', iw1, ' IW2 = ', iw2
            print*, 'IV    = ', iv
            print*, 'GLATV = ', glatv(iv)
            print*, 'GLONV = ', glonv(iv)

            print*, 'dnu(iv),dniu(iv) ', dnu(iv), dniu(iv)

            stop 'STOP U-V edges'
        endif

        quarter_kite(1, iv) = .25 * dvm1 * dnv(iv)
        quarter_kite(2, iv) = .25 * dvm2 * dnv(iv)

    enddo
    !$omp end do
    !$omp end parallel

    !dir$ novector
    do iv = 2, nva
        im1 = itab_v(iv)%im(1)
        im2 = itab_v(iv)%im(2)

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        arm0(im1) = arm0(im1) + 2. * quarter_kite(1, iv)
        arm0(im2) = arm0(im2) + 2. * quarter_kite(2, iv)

        arw0(iw1) = arw0(iw1) + quarter_kite(1, iv) + quarter_kite(2, iv)
        arw0(iw2) = arw0(iw2) + quarter_kite(1, iv) + quarter_kite(2, iv)
    enddo

    ! Lateral boundary copy of arw0

    do iw = 2, nwa
        iwp = itab_w(iw)%iwp
        if (iw /= iwp) arw0(iw) = arw0(iwp)
    enddo

    ! Lateral boundary copy of arm0

    do im = 2, nma
        imp = itab_m(im)%imp
        if (im /= imp) arm0(im) = arm0(imp)
    enddo

    !$omp parallel
    !$omp do private(raxis)
    do iw = 2, nwa

        ! Fill global index (replaced later if this run is parallel)

        itab_w(iw)%iwglobe = iw

        ! Fill outward unit vector components and latitude and longitude of W point

        raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)

        glatw(iw) = atan2(zew(iw), raxis) * piu180
        glonw(iw) = atan2(yew(iw), xew(iw)) * piu180

        wnx(iw) = xew(iw) / erad
        wny(iw) = yew(iw) / erad
        wnz(iw) = zew(iw) / erad
    enddo
    !$omp end do

    !$omp single
    iw = minloc(arw0(2:), 1) + 1
    write(*, *)
    write(*, '(A,f0.4,A)')       " Minimum atmos grid spacing is ", .001 * sqrt(arw0(iw)), " km"
    write(*, '(A,I0,2(A,f0.4))') " at iw = ", iw, ", Lat = ", glatw(iw), ", Lon = ", glonw(iw)

    iw = maxloc(arw0(2:), 1) + 1
    write(*, *)
    write(*, '(A,f0.4,A)')       " Maximum atmos grid spacing is ", .001 * sqrt(arw0(iw)), " km"
    write(*, '(A,I0,2(A,f0.4))') " at iw = ", iw, ", Lat = ", glatw(iw), ", Lon = ", glonw(iw)
    write(*, *)
    !$omp end single

    !$omp do private(npoly,j2,j1,iv,iw1,iw2,dops,npoly1,npoly2,np,&
    !$omp            im1,xw1,xw2,yw1,yw2,xv,yv,alpha)
    do iw = 2, nwa

        ! Number of polygon edges/vertices

        npoly = itab_w(iw)%npoly

        ! Loop over all polygon edges

        do j2 = 1, npoly
            j1 = j2 - 1
            if (j2 == 1) j1 = npoly

            iv = itab_w(iw)%iv(j2)
            iw2 = itab_w(iw)%iw(j2)
            iw1 = itab_w(iw)%iw(j1)

            ! Fractional area of arw0(iw) that is occupied by M and V sectors

            if (itab_v(iv)%iw(1) == iw1) then
                itab_w(iw)%farm(j1) = itab_w(iw)%farm(j1) + quarter_kite(1, iv) / arw0(iw)
                itab_w(iw)%farm(j2) = itab_w(iw)%farm(j2) + quarter_kite(2, iv) / arw0(iw)
            else
                itab_w(iw)%farm(j1) = itab_w(iw)%farm(j1) + quarter_kite(2, iv) / arw0(iw)
                itab_w(iw)%farm(j2) = itab_w(iw)%farm(j2) + quarter_kite(1, iv) / arw0(iw)
            endif
            itab_w(iw)%farv(j2) = (quarter_kite(1, iv) + quarter_kite(2, iv)) / arw0(iw)

            !----------------------------------------
            ! NEW SECTION JULY 2011
            !----------------------------------------

            ! Special - skip gradient calculation if we are at the periodic
            ! domain border and iw1 and iw2 do not share a common vertex

            if (iw == itab_w(iw)%iwp) then
                dops = .true.
            else
                dops = .false.
                npoly1 = itab_w(iw1)%npoly
                npoly2 = itab_w(iw2)%npoly

                do np = 1, npoly1
                    im1 = itab_w(iw1)%im(np)
                    if (im1 == 1) cycle
                    if (any(itab_w(iw2)%im(1:npoly2) == itab_w(iw1)%im(np))) then
                        dops = .true.
                        exit
                    endif
                enddo
            endif

            if (dops) then

                ! Evaluate x,y coordinates of IW1 and IW2 points on polar stereographic plane
                ! tangent at IW

                call e_ps(xew(iw1), yew(iw1), zew(iw1), glatw(iw), glonw(iw), xw1, yw1)
                call e_ps(xew(iw2), yew(iw2), zew(iw2), glatw(iw), glonw(iw), xw2, yw2)
                call e_ps(xev(iv), yev(iv), zev(iv), glatw(iw), glonw(iw), xv, yv)

                ! Coefficients for eastward and northward components of gradient

                itab_w(iw)%gxps1(j1) = yw2 / (xw1 * yw2 - xw2 * yw1)
                itab_w(iw)%gxps2(j1) = -yw1 / (xw1 * yw2 - xw2 * yw1)

                itab_w(iw)%gyps1(j1) = -xw2 / (xw1 * yw2 - xw2 * yw1)
                itab_w(iw)%gyps2(j1) = xw1 / (xw1 * yw2 - xw2 * yw1)

                !----------------------------------------

                if (itab_w(iw)%dirv(j2) < 0.) then
                    alpha = atan2(yw2, xw2)   ! VC(iv) direction counterclockwise from east

                    itab_v(iv)%cosv(1) = cos(alpha)
                    itab_v(iv)%sinv(1) = sin(alpha)

                    itab_v(iv)%dxps(1) = xv
                    itab_v(iv)%dyps(1) = yv
                else
                    alpha = atan2(-yw2, -xw2) ! VC(iv) direction counterclockwise from east

                    itab_v(iv)%cosv(2) = cos(alpha)
                    itab_v(iv)%sinv(2) = sin(alpha)

                    itab_v(iv)%dxps(2) = xv
                    itab_v(iv)%dyps(2) = yv
                endif

            endif

            ! Earth-grid components of rotated polar stereographic easterly
            ! (or cartesian positive x) horizontal unit vector

            itab_w(iw)%unx_w = -sin(glonw(iw))
            itab_w(iw)%uny_w = cos(glonw(iw))



            ! Earth-grid components of rotated polar stereographic northerly
            ! (or cartesian positive y) horizontal unit vector

            itab_w(iw)%vnx_w = -sin(glatw(iw)) * cos(glonw(iw))
            itab_w(iw)%vny_w = -sin(glatw(iw)) * sin(glonw(iw))
            itab_w(iw)%vnz_w = cos(glatw(iw))

            !----------------------------------------
            ! END NEW SECTION JULY 2011
            !----------------------------------------

        enddo

    enddo
    !$omp end do

    ! Loop over all V points

    !$omp do private(ivp,iw1,iw2)
    do iv = 2, nva

        ! Let's not do this section on the boundary cells

        ivp = itab_v(iv)%ivp
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        ! FARW(1) and FARW(2) interpolation coefficients for ARW and VOLT
        ! (taking V control volume to be full DNU(IV) * DNV(IV) rectangle)

        itab_v(iv)%farw(1) = 2. * (quarter_kite(1, iv) + quarter_kite(2, iv)) / arw0(iw1)
        itab_v(iv)%farw(2) = 2. * (quarter_kite(1, iv) + quarter_kite(2, iv)) / arw0(iw2)

        itab_v(iv)%mrlv = max(itab_w(iw1)%mrlw, itab_w(iw2)%mrlw)

    enddo  ! IV
    !$omp end do



    ! Scale eastward and northward gradient components by farm

    !$omp do private(j)
    do iw = 2, nwa

        ! The gradient components are not computed at the lateral boundaries
        if (iw /= itab_w(iw)%iwp) cycle

        do j = 1, itab_w(iw)%npoly

            itab_w(iw)%gxps1(j) = itab_w(iw)%gxps1(j) * itab_w(iw)%farm(j)
            itab_w(iw)%gyps1(j) = itab_w(iw)%gyps1(j) * itab_w(iw)%farm(j)

            itab_w(iw)%gxps2(j) = itab_w(iw)%gxps2(j) * itab_w(iw)%farm(j)
            itab_w(iw)%gyps2(j) = itab_w(iw)%gyps2(j) * itab_w(iw)%farm(j)

        enddo
    enddo
    !$omp end do

    ! Coefficients for converting earth-cartesian velocity to V and W


    !$omp do private(npoly, fo, a, b, work, info, j, iv, vdotw, vmag, fact, &
    !$omp            vnx_ps, vny_ps, vnz_ps, vrot_x, vrot_y, wsize, lwork)
    do iw = 2, nwa

        npoly = itab_w(iw)%npoly

        ! Default coefficients from Perot
        fo(1:npoly) = 2.0_r8 * itab_w(iw)%farv(1:npoly)

        if (allocated(a)) then
            if (size(a, 2) /= npoly) deallocate(a)
        endif

        if (.not. allocated(a)) allocate(a(3, npoly))

        do j = 1, npoly
            iv = itab_w(iw)%iv(j)

            ! Compute the components of the V unit normals perpendicular to W

            vdotw = vnx(iv) * wnx(iw) + vny(iv) * wny(iw) + vnz(iv) * wnz(iw)

            vnx_ps(j) = vnx(iv) - vdotw * wnx(iw)
            vny_ps(j) = vny(iv) - vdotw * wny(iw)
            vnz_ps(j) = vnz(iv) - vdotw * wnz(iw)

            ! Normalize these new vectors to unit length

            vmag = sqrt(vnx_ps(j)**2 + vny_ps(j)**2 + vnz_ps(j)**2)

            vnx_ps(j) = vnx_ps(j) / vmag
            vny_ps(j) = vny_ps(j) / vmag
            vnz_ps(j) = vnz_ps(j) / vmag

            ! Rotate these new unit normals to a coordinate system with Z aligned with W

            if (wnz(iw) >= 0.0) then

                fact = (wny(iw) * vnx_ps(j) - wnx(iw) * vny_ps(j)) / (1.0 + wnz(iw))

                vrot_x(j) = vnx_ps(j) * wnz(iw) - vnz_ps(j) * wnx(iw) + wny(iw) * fact
                vrot_y(j) = vny_ps(j) * wnz(iw) - vnz_ps(j) * wny(iw) - wnx(iw) * fact

            else

                fact = (wny(iw) * vnx_ps(j) - wnx(iw) * vny_ps(j)) / (1._r8 - wnz(iw))

                vrot_x(j) = -vnx_ps(j) * wnz(iw) + vnz_ps(j) * wnx(iw) + wny(iw) * fact
                vrot_y(j) = -vny_ps(j) * wnz(iw) + vnz_ps(j) * wny(iw) - wnx(iw) * fact

            endif

        enddo

        a(1, 1:npoly) = vrot_x(1:npoly) * vrot_x(1:npoly)
        a(2, 1:npoly) = vrot_y(1:npoly) * vrot_y(1:npoly)
        a(3, 1:npoly) = vrot_x(1:npoly) * vrot_y(1:npoly)

        b(1) = 1._r8 - sum(fo(1:npoly) * a(1, :))
        b(2) = 1._r8 - sum(fo(1:npoly) * a(2, :))
        b(3) = - sum(fo(1:npoly) * a(3, :))

        call dgels('N', 3, npoly, 1, a, 3, b, 7, wsize, -1, info)
        lwork = nint(wsize(1)) + 1

        if (allocated(work)) then
            if (size(work) < lwork) deallocate(work)
        endif
        if (.not. allocated(work)) allocate(work(lwork))

        call dgels('N', 3, npoly, 1, a, 3, b, 7, work, size(work), info)

        ! Vector b is now the correction to the coefficients fo
        b(1:npoly) = b(1:npoly) + fo(1:npoly)

        if (info == 0 .and. all(b(1:npoly) > 0.05_r8) .and. all(b(1:npoly) < 0.7_r8)) then

            ! itab_w(iw)%ecvec_vx(1:npoly) = b(1:npoly) * vnx_ps(1:npoly)
            ! itab_w(iw)%ecvec_vy(1:npoly) = b(1:npoly) * vny_ps(1:npoly)
            ! itab_w(iw)%ecvec_vz(1:npoly) = b(1:npoly) * vnz_ps(1:npoly)

            fact = sum(b(1:npoly) * vnx_ps(1:npoly) * vnx(itab_w(iw)%iv(1:npoly)))
            fact = (1._r8 - wnx(iw)**2) / max(fact, 1.e-30_r8)

            itab_w(iw)%ecvec_vx(1:npoly) = b(1:npoly) * vnx_ps(1:npoly) * fact

            fact = sum(b(1:npoly) * vny_ps(1:npoly) * vny(itab_w(iw)%iv(1:npoly)))
            fact = (1._r8 - wny(iw)**2) / max(fact, 1.e-30_r8)

            itab_w(iw)%ecvec_vy(1:npoly) = b(1:npoly) * vny_ps(1:npoly) * fact

            fact = sum(b(1:npoly) * vnz_ps(1:npoly) * vnz(itab_w(iw)%iv(1:npoly)))
            fact = (1._r8 - wnz(iw)**2) / max(fact, 1.e-30_r8)

            itab_w(iw)%ecvec_vz(1:npoly) = b(1:npoly) * vnz_ps(1:npoly) * fact

        else

            write(*, *) "Problem optimizing vector coefficients for iw = ", iw
            write(*, *) glatw(iw), glonw(iw)
            write(*, *) info
            write(*, *) real(b (1:npoly))
            write(*, *) real(fo(1:npoly))
            write(*, *) "Using default coefficients."

            itab_w(iw)%ecvec_vx(1:npoly) = fo(1:npoly) * vnx_ps(1:npoly)
            itab_w(iw)%ecvec_vy(1:npoly) = fo(1:npoly) * vny_ps(1:npoly)
            itab_w(iw)%ecvec_vz(1:npoly) = fo(1:npoly) * vnz_ps(1:npoly)

        endif

    enddo
    !omp end do

    if (allocated(a))    deallocate(a)
    if (allocated(work)) deallocate(work)
    !$omp end parallel
end subroutine grid_geometry_hex

