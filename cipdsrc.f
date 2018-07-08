c    --------------------
      subroutine cipdsrc
c    --------------------
      use common
      implicit double precision(a-h,o-z)

c Vx
      do j=3,j0-2
      do i=3,i0-3
         vxdxmn(i,j) = vxdxmh(i,j) - dt*
     &(  vxdxmh(i,j)*(-vxmh(i-1,j) + vxmh(i+1,j))
     & + vxdymh(i,j)*(-vymx(i-1,j) + vymx(i+1,j)))/(dxm(i-1) + dxm(i))
      enddo
      enddo

      do j=4,j0-3
      do i=2,i0-2
         vxdymn(i,j) = vxdymh(i,j) - dt*
     &(  vxdxmh(i,j)*(-vxmh(i,j-1) + vxmh(i,j+1))
     & + vxdymh(i,j)*(-vymx(i,j-1) + vymx(i,j+1)))/(dy(j-1) + dy(j))
      enddo
      enddo

c Vy
      do j=2,j0-2
      do i=4,i0-3
         vydxmn(i,j) = vydxmh(i,j) - dt*
     &(  vydxmh(i,j)*(-vxmy(i-1,j) + vxmy(i+1,j))
     & + vydymh(i,j)*(-vymh(i-1,j) + vymh(i+1,j)))/(dx(i-1) + dx(i))
      enddo
      enddo

      do j=3,j0-3
      do i=3,i0-2
         vydymn(i,j) = vydymh(i,j) - dt*
     &(  vydxmh(i,j)*(-vxmy(i,j-1) + vxmy(i,j+1))
     & + vydymh(i,j)*(-vymh(i,j-1) + vymh(i,j+1)))/(dym(j-1) + dym(j))
      enddo
      enddo

c Vz, Pressure, Temperature
      do j=4,j0-3
      do i=4,i0-3
         vzdxn(i,j) = vzdxh(i,j) - dt*
     &(  vzdxh(i,j)*(-vx(i-1,j) + vx(i+1,j))
     & + vzdyh(i,j)*(-vy(i-1,j) + vy(i+1,j)))/(dx(i-1) + dx(i))

         vzdyn(i,j) = vzdyh(i,j) - dt*
     &(  vzdxh(i,j)*(-vx(i,j-1) + vx(i,j+1))
     & + vzdyh(i,j)*(-vy(i,j-1) + vy(i,j+1)))/(dy(j-1) + dy(j))

         prdxn(i,j) = prdxh(i,j) - dt*
     &(  prdxh(i,j)*(-vx(i-1,j) + vx(i+1,j))
     & + prdyh(i,j)*(-vy(i-1,j) + vy(i+1,j)))/(dx(i-1) + dx(i))

         prdyn(i,j) = prdyh(i,j) - dt*
     &(  prdxh(i,j)*(-vx(i,j-1) + vx(i,j+1))
     & + prdyh(i,j)*(-vy(i,j-1) + vy(i,j+1)))/(dy(j-1) + dy(j))

         tedxn(i,j) = tedxh(i,j) - dt*
     &(  tedxh(i,j)*(-vx(i-1,j) + vx(i+1,j))
     & + tedyh(i,j)*(-vy(i-1,j) + vy(i+1,j)))/(dx(i-1) + dx(i))

         tedyn(i,j) = tedyh(i,j) - dt*
     &(  tedxh(i,j)*(-vx(i,j-1) + vx(i,j+1))
     & + tedyh(i,j)*(-vy(i,j-1) + vy(i,j+1)))/(dy(j-1) + dy(j))
      enddo
      enddo

      return
      end subroutine cipdsrc
