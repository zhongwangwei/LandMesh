module MOD_Area_judge
    USE consts_coms, only: io6, r8, refine, edgee, edgew, edges, edgen, ndm_domain, nlons_source, nlats_source, openmp, lcs
    USE refine_vars, only: edgee_rf, edgew_rf, edges_rf, edgen_rf, ndm_refine
    USE MOD_data_preprocess, only: landtypes, lon_vertex, lat_vertex
    implicit none
    ! 判断是否是海洋>网格，判断是否在domain, 判断是否在refine
    integer, allocatable, public :: seaorland(:, :)
    logical, allocatable, public :: IsInDmArea_grid(:, :), IsInRfArea_grid(:, :)
    integer, public :: minlon_DmArea, maxlon_DmArea, maxlat_DmArea, minlat_DmArea
    integer, public :: minlon_RfArea, maxlon_RfArea, maxlat_RfArea, minlat_RfArea
    contains

    SUBROUTINE Area_judge()

        IMPLICIT NONE
        integer :: i, j, sum_land, sum_sea
        integer :: numpatch_Dm, numpatch_Rf
        write(io6, *) 'IsInArea_grid_Calculation start'
        allocate(IsInDmArea_grid(nlons_source, nlats_source)); IsInDmArea_grid = .false.
        CALL IsInArea_grid_Calculation(edgee, edgew, edges, edgen, ndm_domain, IsInDmArea_grid, minlon_DmArea, maxlon_DmArea, maxlat_DmArea, minlat_DmArea, numpatch_Dm) 
        write(io6, *) 'IsInArea_grid_Calculation complete'
        
        write(io6, *) 'sea or land judge start'
        allocate(seaorland(nlons_source, nlats_source)); seaorland = 0
        !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
        !$OMP PRIVATE(i, j)
        do j = maxlat_DmArea, minlat_DmArea, 1 ! 影响先内层循环
            do i = minlon_DmArea, maxlon_DmArea, 1
                if (IsInDmArea_grid(i, j)) then
                    ! 包含关系计算区域内海陆网格判断, 只对0处理就好了，maxlc是冰川或者水体可以先不弄
                    if (landtypes(i, j) /= 0) seaorland(i, j) = 1 ! 1 表示为陆地网格
                end if
            end do
        end do
        !$OMP END PARALLEL DO
        sum_land = sum(seaorland)
        ! sum_sea  = numpatch_Dm - sum_land
        print*,"包含关系计算区域内, 陆地网格个数为",sum_land ! 海洋网格个数为",sum_sea,

        write(io6, *) 'make grid with refine mesh in the MOD_Area_judge.F90'
        if (refine) then! 计算并记录位于计算细化区域的经纬度网格
            allocate(IsInRfArea_grid(nlons_source, nlats_source)); IsInRfArea_grid = .false. 
            ! CALL IsInArea_grid_Calculation(edgee_rf, edgew_rf, edges_rf, edgen_rf, ndm_refine, IsInRfArea_grid)
            CALL IsInArea_grid_Calculation(edgee_rf, edgew_rf, edges_rf, edgen_rf, ndm_refine, IsInRfArea_grid, minlon_RfArea, maxlon_RfArea, maxlat_RfArea, minlat_RfArea, numpatch_rf) 
            !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
            !$OMP PRIVATE(i, j)
            do j = maxlat_RfArea, minlat_RfArea, 1 ! 影响先内层循环
                do i = minlon_RfArea, maxlon_RfArea, 1
                    if (IsInRfArea_grid(i, j)) then
                        if (IsInDmArea_grid(i, j) .eqv. .false.) then ! 保证细化区域网格位于包含关系计算区域
                            print*,"ERROR!!! the refine area exceed the domain area!!!"
                            stop ! 错误！！！细化区域和包含关系计算区域之间没有交集！！！”
                        end if
                    end if
                end do
            end do
            !$OMP END PARALLEL DO
            ! “细化区域完全位于包含关系计算区域内……”
            print*,"the refine area Completely locate at the domain area!!!"
        end if
    END SUBROUTINE Area_judge
    
    SUBROUTINE IsInArea_grid_Calculation(edgee_temp, edgew_temp, edges_temp, edgen_temp, ndm, IsInArea_grid, minlon_source, maxlon_source, maxlat_source, minlat_source, numpatch)
        ! 判断经纬度网格是否位于细化区域/包含关系计算区域
        implicit none
        integer, intent(in) :: ndm
        real(r8), dimension(:), intent(in) :: edgee_temp, edgew_temp, edges_temp, edgen_temp
        integer :: n, temp1, temp2, temp3, temp4
        logical, intent(inout) :: IsInArea_grid(:, :)
        integer, intent(out) :: minlon_source, maxlon_source, maxlat_source, minlat_source, numpatch

        minlon_source = 1; maxlon_source = nlons_source
        maxlat_source = 1; minlat_source = nlats_source
        numpatch = 0
        do n = 1, ndm, 1
            ! lon from -180 to 180   lat from 90 to -90
            CALL Source_Find(edgew_temp(n), lon_vertex, 'lon', temp1)! minlon_source
            CALL Source_Find(edgee_temp(n), lon_vertex, 'lon', temp2)! maxlon_source
            CALL Source_Find(edgen_temp(n), lat_vertex, 'lat', temp3)! maxlat_source
            CALL Source_Find(edges_temp(n), lat_vertex, 'lat', temp4)! minlat_source 纬度大值反而是小的索引值

            ! Adjust the boundaries ! different from MOD_GetContain.F90 Line256-257
            temp2 = temp2 - 2
            temp4 = temp4 - 2
            if (temp2 == nlons_source-1) temp2 = temp2 + 1
            if (temp4 == nlats_source-1) temp4 = temp4 + 1
            ! calculate the size of ustrgrid need and adjust IsInArea_ustr
            IsInArea_grid(temp1:temp2,temp3:temp4) = .true.
            numpatch = numpatch + (temp2 - temp1 + 1) * (temp4 - temp3 + 1)
            ! Update boundaries if necessary
            if (temp1 > minlon_source) minlon_source = temp1
            if (temp2 < maxlon_source) maxlon_source = temp2
            if (temp3 > maxlat_source) maxlat_source = temp3
            if (temp4 < minlat_source) minlat_source = temp4
           
        end do

    END SUBROUTINE IsInArea_grid_Calculation


    SUBROUTINE Source_Find(temp, seq_lonlat, str1, source)
        implicit none

        real(r8), intent(in) :: temp
        real(r8), dimension(:), intent(in) :: seq_lonlat
        character(LEN = 3), intent(in) :: str1
        integer :: i, gridnum_perdegree, minsource, maxsource
        integer, intent(out) :: source

        gridnum_perdegree = 120
        if (lcs == 'igbp') gridnum_perdegree = 240
        if (trim(str1) == 'lon') then ! 针对lon序列! 经度度是从-180到180 !注意出现跨越180经线的>问题
            minsource = ( floor(temp)   - (-180)) * gridnum_perdegree
            maxsource = ( ceiling(temp) - (-180)) * gridnum_perdegree
            minsource = max(1,              minsource-10)
            maxsource = min(1+nlons_source, maxsource+10)
            do i = minsource, maxsource, 1 ! because
                if (temp <= seq_lonlat(i)) then ! 满足条件就退出，等号是针对极点出现的情况
                    source = i ! 找下/左边界,上/右边界都是一样的
                    return ! exit
                end if
            end do
            print*, "temp <= seq_lonlat(i) is not exist!!!!!!!!!!!!!!!"
            print*, "temp = ", temp
            print*, "minsource = ", minsource, "maxsource = ", maxsource
            print*, "seq_lonlat(minsource) = ", seq_lonlat(minsource)
            print*, "seq_lonlat(maxsource) = ", seq_lonlat(maxsource)
            source = i
        else
            minsource = ( 90 - ceiling(temp)) * gridnum_perdegree
            maxsource = ( 90 -   floor(temp)) * gridnum_perdegree
            minsource = max(1,              minsource-10)
            maxsource = min(1+nlats_source, maxsource+10)
            do i = minsource, maxsource, 1
                if (temp >= seq_lonlat(i)) then ! 满足条件就退出，等号是针对极点出现的情况
                    source = i ! 找下/左边界,上/右边界都是一样的
                    return ! exit
                end if
            end do
            print*, "temp >= seq_lonlat(i) is not exist!!!!!!!!!!!!!!!"
            print*, "temp = ", temp
            print*, "minsource = ", minsource, "maxsource = ", maxsource
            print*, "seq_lonlat(minsource) = ", seq_lonlat(minsource)
            print*, "seq_lonlat(maxsource) = ", seq_lonlat(maxsource)
            source = i
        end if

    END SUBROUTINE Source_Find

END Module MOD_Area_judge
