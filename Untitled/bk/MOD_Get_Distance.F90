module MOD_Get_Distance

   use consts_coms, only : r8, erad, pio180 ! add erad and pio180

   implicit none

   contains
   
   SUBROUTINE Get_Dis(mp, wp, ngrmw, ngrwm, dismm, disww, sjx_points, lbx_points)
      
      use refine_vars

      implicit none

      integer,intent(in) :: sjx_points,lbx_points
      integer,intent(in) :: ngrmw(3,sjx_points), ngrwm(7,lbx_points)
      !ngrmw(sjx_points,3), ngrwm(lbx_points,7) ! in the further after Array transpose!
      real(r8),intent(in) :: mp(sjx_points,2), wp(lbx_points,2) ! 三角形中心的的经纬度! 多边形中心的经纬度
      real(r8),intent(out) :: dismm(lbx_points,7), disww(sjx_points,3)
      integer :: i, j, edges! edges to make sure edge num of polygon
      real(r8) :: x1, x2, y1, y2, h(8)! array only use in the MOD_Get_Distance
      h = 0
      edges = 3 
      disww = 0.
      dismm = 0.
      ! tri start from 2
      do i = 2, sjx_points, 1 
         h(1:3) = ngrmw(:,i)!ngrmw(i, 1:3)
         h(4)   = ngrmw(1,i)!ngrmw(i, 1)
         do j = 1, 3, 1
            x1 = wp( h(j)   ,1)
            y1 = wp( h(j)   ,2)
            x2 = wp( h(j+1) ,1)
            y2 = wp( h(j+1) ,2)
            disww(i,j) = Distance(x1, x2, y1, y2) ! 存三角形顶点之间的位置距离
         end do
      end do  
      ! polygon start from 2
      do i = 2, lbx_points, 1
      if(ngrwm(6, i) == 1)then! 1 mean this polygon do not connect other node 
            edges = 5
         else if(ngrwm(7, i) == 1)then
            edges = 6
         else
            edges= 7
         end if
         h(1:edges) = ngrwm(1:edges, i)!ngrwm(i, 1:edges)
         h(edges+1) = ngrwm(1, i)!ngrwm(i, 1)
         do j = 1, edges, 1
            x1 = mp( h(j)   ,1)
            y1 = mp( h(j)   ,2)
            x2 = mp( h(j+1) ,1)
            y2 = mp( h(j+1) ,2)
            dismm(i,j) = Distance(x1, x2, y1, y2) ! 存多边形顶点之间的位置距离
         end do
      end do

   END SUBROUTINE Get_dis

   Function Distance(x1, x2, y1, y2)

      implicit none
      
      real(r8),intent(in) :: x1, x2, y1, y2
      real(r8) :: px1, py1, px2, py2, v, Distance
      px1 = x1 * pio180 ! lon1
      py1 = y1 * pio180 ! lat1
      px2 = x2 * pio180 ! lon2
      py2 = y2 * pio180 ! lat2
      v = sin(py1) * sin(py2) + cos(py1) * cos(py2) * cos(px1 - px2)
      Distance = erad * acos(v) !erad = 6371.22e3 Earth radius [km]
      return

   END Function Distance
      
END Module MOD_Get_Distance
