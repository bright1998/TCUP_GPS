c    -----------------
      subroutine cflc
c    -----------------
      use common
      implicit double precision(a-h,o-z)

      if(ncflconst .eq. 0) then
         dtmin = 1.e10
         do j=2,j0-1
         do i=2,i0-1
            vx1 = vx(i,j)
            vy1 = vy(i,j)
            vz1 = vz(i,j)
            pr1 = pr(i,j)
            dt = dmin1(dxm(i-1),dym(j-1))/
     &           (Cs(i,j) + dsqrt(vx1**2 + vy1**2 + vz1**2))
            if(dt .lt. dtmin) dtmin = dt
         enddo
         enddo
         dt = cn*dtmin
      else
         dt = dtcon
      endif

      return
      end subroutine cflc
