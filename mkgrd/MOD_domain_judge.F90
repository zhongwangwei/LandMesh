! 判断细化区域是否与计算包含关系区域有交界
! 若没有，则中断并提示
module MOD_domain_judge

   !use netcdf
   USE consts_coms
   USE refine_vars

   implicit none

   contains
   SUBROUTINE domain_judge ()

      integer :: i,j
      real(r8) :: dx,dy
      real(r8),allocatable :: lon_i(:),lat_i(:)

      logical :: intersection
      logical,allocatable :: IsInRfArea(:,:)
      logical,allocatable :: IsInDmArea(:,:)

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

      allocate(IsInRfArea(nlons_source,nlats_source))
      allocate(IsInDmArea(nlons_source,nlats_source))
      IsInRfArea = .false.
      IsInDmArea = .false.
      ! 计算并记录位于细化区域/计算包含关系区域内的经纬度网格
      CALL IsInRefineArea(IsInRfArea,lon_i,lat_i,dx,dy)
      CALL IsInDomainArea(IsInDmArea,lon_i,lat_i,dx,dy)

      intersection = .false.

      do i = 1,nlons_source,1
         do j = 1,nlats_source,1
            if((IsInRfArea(i,j) .eqv. .true.).and.(IsInDmArea(i,j) .eqv. .true.))then
               intersection = .true.
               exit
            end if
         end do
      end do

      if(intersection .eqv. .false.)then
         print*,"ERROR!!! There is no intersection between the refinement area and the inclusion relation calculation area!!!"
         stop
      else
         print*,"There is intersection between the refinement area and the inclusion relation calculation area..."
      end if

      deallocate(lon_i)
      deallocate(lat_i)
      deallocate(IsInRfArea)
      deallocate(IsInDmArea)

   END SUBROUTINE domain_judge

   ! Determine whether the latitude and longitude grid is located in the refinement region
    SUBROUTINE IsInRefineArea(IsInRfArea,lon_i,lat_i,dx,dy)

      use refine_vars

      implicit none

      integer :: i,j,n
      real(r8),intent(in) :: dx,dy
      real(r8),dimension(nlons_source),intent(in) :: lon_i
      real(r8),dimension(nlats_source),intent(in) :: lat_i
      real(r8),dimension(nlons_source) :: lone,lonw
      real(r8),dimension(nlats_source) :: latn,lats
      logical,dimension(nlons_source,nlats_source),intent(out) :: IsInRfArea


      do i = 1,nlons_source,1
         lone(i) = lon_i(i) + dx / 2.
         lonw(i) = lon_i(i) - dx / 2.
      end do

      do j = 1,nlats_source,1
         latn(j) = lat_i(j) + dy / 2.
         lats(j) = lat_i(j) - dy / 2.
      end do

      IsInRfArea = .false.

      do n = 1,ndm_refine,1
         do i = 1,nlons_source,1
            do j = 1,nlats_source,1
               if((lone(i) < edgee_rf(n)).and.(lonw(i) > edgew_rf(n)).and.&
                        (latn(j) < edgen_rf(n)).and.(lats(j) > edges_rf(n)))then
                  IsInRfArea(i,j) = .true.
               end if
            end do
         end do
      end do

   END SUBROUTINE IsInRefineArea


   ! 判断经纬度网格是否位于计算包含关系区域
   SUBROUTINE IsInDomainArea(IsInDmArea,lon_i,lat_i,dx,dy)

      implicit none

      integer :: i,j,n
      real(r8),intent(in) :: dx,dy
      real(r8),dimension(nlons_source),intent(in) :: lon_i
      real(r8),dimension(nlats_source),intent(in) :: lat_i
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

END Module MOD_domain_judge
