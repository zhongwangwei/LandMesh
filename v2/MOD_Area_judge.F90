module MOD_Area_judge
    USE consts_coms, only: io6, r8, refine, edgee, edgew, edges, edgen, ndm_domain, nlons_source, nlats_source, openmp
    USE refine_vars, only: edgee_rf, edgew_rf, edges_rf, edgen_rf, ndm_refine
    USE MOD_data_preprocess, only: landtypes, lon_vertex, lat_vertex
    implicit none
    ! 判断是否是海洋>网格，判断是否在domain, 判断是否在refine
    integer, allocatable, public :: seaorland(:, :)
    logical, allocatable, public :: IsInDmArea_grid(:, :), IsInRfArea_grid(:, :)
    contains

    SUBROUTINE Area_judge()

       IMPLICIT NONE
       integer :: i, j, sum_land, sum_sea

       write(io6, *) 'IsInArea_grid_Calculation start'
       allocate(IsInDmArea_grid(nlons_source, nlats_source)); IsInDmArea_grid = .false. 
       CALL IsInArea_grid_Calculation(edgee, edgew, edges, edgen, ndm_domain, IsInDmArea_grid) 
       write(io6, *) 'IsInArea_grid_Calculation complete'
      
       write(io6, *) 'sea or land judge start'
       allocate(seaorland(nlons_source, nlats_source)); seaorland = 0
       sum_sea = 0
       sum_land = 0
       do j = 1, nlats_source, 1 ! 影响先内层循环
           do i = 1, nlons_source, 1
               if (IsInDmArea_grid(i, j) .eqv. .true.) then
                   if (landtypes(i, j) /= 0) then ! 包含关系计算区域内海陆网格判断，只对0处理就好了，maxlc是冰川或者水体可以先不弄
                       seaorland(i, j) = 1 ! 1 表示为陆地网格
                       sum_land = sum_land + 1
                   else
                       sum_sea = sum_sea + 1
                   end if
               end if
           end do
       end do

       print*,"包含关系计算区域内，海洋网格个数为",sum_sea,"，陆地网格个数为",sum_land
       if (sum_land + sum_sea == 0) stop "stop for sum_land + sum_sea == 0 !ERROR!"
       ! print*, "sum_land + sum_sea = ", sum_land + sum_sea
       ! stop "stop for seaorland test" 
       write(io6, *) 'make grid with refine mesh in the MOD_Area_judge.F90'
       if (refine .eqv. .True.) then! 计算并记录位于计算细化区域的经纬度网格
           allocate(IsInRfArea_grid(nlons_source, nlats_source)); IsInRfArea_grid = .false. 
           CALL IsInArea_grid_Calculation(edgee_rf, edgew_rf, edges_rf, edgen_rf, ndm_refine, IsInRfArea_grid) 
           !$OMP PARALLEL DO NUM_THREADS(openmp) SCHEDULE(DYNAMIC)&
           !$OMP PRIVATE(i, j)
           do j = 1, nlats_source, 1 ! 影响先内层循环
               do i = 1, nlons_source, 1
                   if (IsInDmArea_grid(i, j) .eqv. .False.) then
                      if (IsInRfArea_grid(i, j) .eqv. .true.) then ! 保证细化区域网格位于包含关系计算区域
                          print*,"ERROR!!! the refinement area exceed the inclusion relation calculation area!!!"
                          stop ! 错误！！！细化区域和包含关系计算区域之间没有交集！！！”
                      end if
                   end if
               end do
           end do
           !$OMP END PARALLEL DO
           ! “细化区域完全位于包含关系计算区域内……”
           print*,"the refinement area Completely locate the inclusion relation calculation area!!!"
       end if


    END SUBROUTINE Area_judge
    
    SUBROUTINE IsInArea_grid_Calculation(edgee_temp, edgew_temp, edges_temp, edgen_temp, ndm, IsInArea_grid)
      ! 判断经纬度网格是否位于细化区域/包含关系计算区域
      implicit none
      integer, intent(in) :: ndm 
      real(r8), dimension(:), intent(in) :: edgee_temp, edgew_temp, edges_temp, edgen_temp
      integer :: n, i, j
      logical, intent(inout) :: IsInArea_grid(:, :)
      ! old scheme from FHW code
      do n = 1, ndm, 1
         do i = 1, nlons_source, 1    
            do j = 1, nlats_source, 1
               if((lon_vertex(i+1) <= edgee_temp(n)).and.(lon_vertex(i)   >= edgew_temp(n)).and.&
                  (lat_vertex(j)   <= edgen_temp(n)).and.(lat_vertex(j+1) >= edges_temp(n)))then
                  IsInArea_grid(i, j) = .true.
               end if
            end do
         end do
      end do

      ! New scheme for Rui Zhang but have a minor bug    
      ! IsInArea_grid = .false. ! 全数组赋值
      ! gridnum_perdegree = 120
      ! do n = 1, ndm, 1
         ! 经度度是从-180到180
      !    temp1 = (edgew_temp(n) - (-180)) * gridnum_perdegree ! 获取从-180经度到左边界edgew_rf的格点数
      !    temp2 = (edgee_temp(n) - (-180)) * gridnum_perdegree ! 获取左边界到右左边界的格数
         ! 纬度是从90到-90
      !    temp3 = (90 - edgen_temp(n)) * gridnum_perdegree ! 获取从90纬度到上边界edgen_rf的格点数
      !    temp4 = (90 - edges_temp(n)) * gridnum_perdegree ! 获取上边界到下边界的格数
         ! close in the left and open in the right
      !    IsInArea_grid(temp1+1:temp2,temp3+1:temp4) = .true.
      !    print*,'temp1=',temp1
      !    print*,'temp2=',temp2
      !    print*,'temp3=',temp3 
      !    print*,'temp4=',temp4
      ! end do

    END SUBROUTINE IsInArea_grid_Calculation

END Module MOD_Area_judge
