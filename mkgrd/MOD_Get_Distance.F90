module MOD_Get_Distance

   use consts_coms

   implicit none

   contains
   SUBROUTINE Get_Dis(mp,wp,ngrmw,ngrwm,dismm,disww,sjx_points,lbx_points)
      
      use refine_vars

      implicit none

      integer :: i,j
      real(r8) :: dis,x1,x2,y1,y2

      integer,intent(in) :: sjx_points,lbx_points
      integer,intent(in) :: ngrmw(3,sjx_points)
      integer,intent(in) :: ngrwm(7,lbx_points)
      real(r8),intent(in) :: mp(sjx_points,2)
      real(r8),intent(in) :: wp(lbx_points,2)
      real(r8),intent(out) :: dismm(lbx_points,7)
      real(r8),intent(out) :: disww(sjx_points,3)

      integer :: n_ngrmw(sjx_points)
      integer :: n_ngrwm(lbx_points)

      n_ngrmw = 3
      n_ngrwm = 7

      do i = 1,sjx_points,1
         do j = 1,3,1
            if(ngrmw(j,i) == 1)then
               n_ngrmw(i) = n_ngrmw(i) - 1 
            end if
         end do
      end do

      do i = 1,lbx_points,1
         do j = 1,7,1
            if(ngrwm(j,i) == 1)then
               n_ngrwm(i) = n_ngrwm(i) - 1
            end if
         end do
      end do

      do i = 1,sjx_points,1
         if(n_ngrmw(i) /= 3)then
            cycle
         end if
         do j = 1,3,1
            x1 = wp(ngrmw(j,i),1)
            y1 = wp(ngrmw(j,i),2)
            if(j == 3)then
               x2 = wp(ngrmw(1,i),1)
               y2 = wp(ngrmw(1,i),2)
            else
               x2 = wp(ngrmw(j+1,i),1)
               y2 = wp(ngrmw(j+1,i),2)
            end if
         
            disww(i,j) = Distance(x1,x2,y1,y2)
         end do
      end do

      do i = 1,lbx_points,1
         if(n_ngrwm(i) < 5)then
            cycle
         end if
         do j = 1,n_ngrwm(i),1
            x1 = mp(ngrwm(j,i),1)
            y1 = mp(ngrwm(j,i),2)
            if(j == n_ngrwm(i))then
               x2 = mp(ngrwm(1,i),1)
               y2 = mp(ngrwm(1,i),2)
            else
               x2 = mp(ngrwm(j+1,i),1)
               y2 = mp(ngrwm(j+1,i),2)
            end if

            dismm(i,j) = Distance(x1,x2,y1,y2)
         end do
      end do

   END SUBROUTINE Get_dis


   Function HaverSin(theta)

      real(r8),intent(in) :: theta
      real(r8) :: HaverSin

      HaverSin = sin(theta/2.)*sin(theta/2.)
      
      return

   END FUNCTION HaverSin


   Function Distance(x1,x2,y1,y2)

      real(r8),intent(in) :: x1,x2,y1,y2
      real(r8) :: EARTH_RADIUS,PI
      real(r8) :: rx1,ry1,rx2,ry2,vx,vy
      real(r8) :: h,Distance

      EARTH_RADIUS = 6371.0
      PI = 3.14

      rx1 = x1*PI/180
      ry1 = y1*PI/180
      rx2 = x2*PI/180
      ry2 = y2*PI/180

      vx = abs(rx1-rx2)
      vy = abs(ry1-ry2)

      h = HaverSin(vy)+cos(ry1)*cos(ry2)*HaverSin(vx)
      Distance = 2 * EARTH_RADIUS * asin(sqrt(h))
      
      return

   END Function Distance

      
END Module MOD_Get_Distance
