c    -------------------
      subroutine prbase
c    -------------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine spheat(pre,tem,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: pre,tem
            integer :: is,ie,js,je
         end subroutine spheat

         subroutine eos(den,pre,tem,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: pre,tem
            double precision,dimension(:,:),intent(out) :: den
            integer :: is,ie,js,je
         end subroutine eos

         subroutine bdfrdx(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdfrdx

         subroutine bdfrex(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdfrex

         subroutine bdperx(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdperx

         subroutine bdsymx(da,margin,mbnd,men,coff)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: coff
            integer :: margin,mbnd,men
         end subroutine bdsymx

         subroutine bdconx(da,margin,mbnd,men,const)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: const
            integer :: margin,mbnd,men
         end subroutine bdconx

         subroutine bdinix(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdinix

         subroutine bdfrdy(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdfrdy

         subroutine bdfrey(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdfrey

         subroutine bdpery(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdpery

         subroutine bdsymy(da,margin,mbnd,men,coff)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: coff
            integer :: margin,mbnd,men
         end subroutine bdsymy

         subroutine bdcony(da,margin,mbnd,men,const)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: const
            integer :: margin,mbnd,men
         end subroutine bdcony

         subroutine bdiniy(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdiniy
      end interface

      dimension ::  s(1:i0,1:j0),c1(1:i0,1:j0),c2(1:i0,1:j0),
     &             d1(1:i0,1:j0),d2(1:i0,1:j0), e(1:i0,1:j0)
      dimension :: r(1:i0,1:j0),work(1:i0,1:j0)

      dimension :: DIV(1:i0,1:j0),
     &             averoxi(1:i0,1:j0),averoyi(1:i0,1:j0)
      dimension :: forcex(1:i0,1:j0),forcey(1:i0,1:j0),forcez(1:i0,1:j0)

      dti = 1.d0/dt
      dti2 = dti**2

      call spheat(prh2,teh2,5,i0-4,5,j0-4)
      call eos(roh,prh2,teh2,5,i0-4,5,j0-4)

      do j=5,j0-4
      do i=5,i0-4
         delvx =-vxmh2(i-1,j) + vxmh2(i,j)
         delvy =-vymh2(i,j-1) + vymh2(i,j)
         DIV(i,j) = delvx*dxmi(i-1) + delvy*dymi(j-1)
      enddo
      enddo

      do j=5,j0-4
      do i=5,i0-5
         averoxi(i,j) = 2.d0/(roh(i,j) + roh(i+1,j))
      enddo
      enddo

      do j=5,j0-5
      do i=5,i0-4
         averoyi(i,j) = 2.d0/(roh(i,j) + roh(i,j+1))
      enddo
      enddo

C Force except pressure gradient force
      do j=1,j0
      do i=1,i0-1
         forcex(i,j) = Fexxm(i,j)
      enddo
      enddo
      do j=1,j0-1
      do i=1,i0
         forcey(i,j) = Fexym(i,j)
      enddo
      enddo
      do j=1,j0
      do i=1,i0
         forcez(i,j) = Fexz(i,j)
      enddo
      enddo

      do j=6,j0-5
      do i=6,i0-5
         dxdi = dxmi(i-1)*dxi(i)
         dxci = dxmi(i-1)*dxi(i-1)
         dydi = dymi(j-1)*dyi(j)
         dyci = dymi(j-1)*dyi(j-1)

          s(i,j) = prh2(i,j)*dti2/(roh(i,j)*Cs(i,j)**2) - DIV(i,j)*dti
     &           + (forcex(i-1,j) - forcex(i,j))*dxmi(i-1)
     &           + (forcey(i,j-1) - forcey(i,j))*dymi(j-1)
         c1(i,j) =-dxci*averoxi(i-1,j)
         c2(i,j) =-dyci*averoyi(i,j-1)
         d1(i,j) =-dxdi*averoxi(i,j)
         d2(i,j) =-dydi*averoyi(i,j)
          e(i,j) = dti2/(roh(i,j)*Cs(i,j)**2)
     &           + dxci*averoxi(i-1,j) + dxdi*averoxi(i,j)
     &           + dyci*averoyi(i,j-1) + dydi*averoyi(i,j)
      enddo
      enddo

      sigmar0 = 0.d0
      do j=6,j0-5
      do i=6,i0-5
         r(i,j) = c2(i,j)*prh2(i,j-1) + c1(i,j)*prh2(i-1,j)
     &          +  e(i,j)*prh2(i,j)   + d1(i,j)*prh2(i+1,j)
     &          + d2(i,j)*prh2(i,j+1) -  s(i,j)
         sigmar0 = sigmar0 + abs(r(i,j))
      enddo
      enddo

      do j=5,j0-4
      do i=5,i0-4
         work(i,j) = prh2(i,j)
      enddo
      enddo

      do n=1,itemax
C Red
         do j=6,j0-5
            na = int(mod(j,2))
            do i=6+na,i0-5,2
               r(i,j) = c2(i,j)*work(i,j-1) + c1(i,j)*work(i-1,j)
     &                +  e(i,j)*work(i,j)   + d1(i,j)*work(i+1,j)
     &                + d2(i,j)*work(i,j+1) -  s(i,j)
            enddo
            do i=6+na,i0-5,2
               workn = work(i,j) - r(i,j)/e(i,j)
               work(i,j) = (1.d0 - omega)*work(i,j) + omega*workn
            enddo
         enddo
C Black
         do j=6,j0-5
            na = int(mod(j,2))
            do i=7-na,i0-5,2
               r(i,j) = c2(i,j)*work(i,j-1) + c1(i,j)*work(i-1,j)
     &                +  e(i,j)*work(i,j)   + d1(i,j)*work(i+1,j)
     &                + d2(i,j)*work(i,j+1) -  s(i,j)
            enddo
            do i=7-na,i0-5,2
               workn = work(i,j) - r(i,j)/e(i,j)
               work(i,j) = (1.d0 - omega)*work(i,j) + omega*workn
            enddo
         enddo

C Boundary condition
         if(nbnd1 .eq. 1) then
            call bdfrdx(work,pr0,5,0,0)
         elseif(nbnd1 .eq. 2) then
            call bdfrex(work,5,0,0)
         elseif(nbnd1 .eq. 3) then
            call bdperx(work,5,0,0)
         elseif(nbnd1 .eq. 4) then
            call bdsymx(work,5,0,0,1.d0)
         elseif(nbnd1 .eq. 5) then
            call bdfrex(work,5,0,0)
         elseif(nbnd1 .eq. 6) then
            call bdfrex(work,5,0,0)
         endif

         if(nbnd2 .eq. 1) then
            call bdfrdx(work,pr0,5,1,0)
         elseif(nbnd2 .eq. 2) then
            call bdfrex(work,5,1,0)
         elseif(nbnd2 .eq. 3) then
            call bdperx(work,5,1,0)
         elseif(nbnd2 .eq. 4) then
            call bdsymx(work,5,1,0,1.d0)
         elseif(nbnd2 .eq. 5) then
            call bdfrex(work,5,1,0)
         elseif(nbnd2 .eq. 6) then
            call bdfrex(work,5,1,0)
         endif

         if(nbnd3 .eq. 1) then
            call bdfrdy(work,pr0,5,0,0)
         elseif(nbnd3 .eq. 2) then
            call bdfrey(work,5,0,0)
         elseif(nbnd3 .eq. 3) then
            call bdpery(work,5,0,0)
         elseif(nbnd3 .eq. 4) then
            call bdsymy(work,5,0,0,1.d0)
         elseif(nbnd3 .eq. 5) then
            call bdfrey(work,5,0,0)
         elseif(nbnd3 .eq. 6) then
            call bdfrey(work,5,0,0)
         endif

         if(nbnd4 .eq. 1) then
            call bdfrdy(work,pr0,5,1,0)
         elseif(nbnd4 .eq. 2) then
            call bdfrey(work,5,1,0)
         elseif(nbnd4 .eq. 3) then
            call bdpery(work,5,1,0)
         elseif(nbnd4 .eq. 4) then
            call bdsymy(work,5,1,0,1.d0)
         elseif(nbnd4 .eq. 5) then
            call bdfrey(work,5,1,0)
         elseif(nbnd4 .eq. 6) then
            call bdfrey(work,5,1,0)
         endif

         sigmar = 0.d0
         do j=6,j0-5
         do i=6,i0-5
            sigmar = sigmar + abs(r(i,j))
         enddo
         enddo

         if(sigmar/sigmar0 .lt. epsit) goto 100
      enddo

 100  continue
      do j=6,j0-5
      do i=6,i0-5
         prn(i,j) = work(i,j)
      enddo
      enddo

      do j=6,j0-5
      do i=6,i0-6
         vxmn(i,j) = vxmh2(i,j)
     &             + dt*(-averoxi(i,j)*(-prn(i,j) + prn(i+1,j))*dxi(i)
     &                                 + forcex(i,j))
      enddo
      enddo

      do j=6,j0-6
      do i=6,i0-5
         vymn(i,j) = vymh2(i,j)
     &             + dt*(-averoyi(i,j)*(-prn(i,j) + prn(i,j+1))*dyi(j)
     &                                 + forcey(i,j))
      enddo
      enddo

      do j=5,j0-4
      do i=5,i0-4
         vzn(i,j) = vzh2(i,j) + dt*forcez(i,j)
      enddo
      enddo

      do j=6,j0-5
      do i=6,i0-5
         ten(i,j) = teh2(i,j)
     &            + (1.d0/(roh(i,j)*Cp(i,j)) + rmuj(i,j))
     &             *(prn(i,j) - prh2(i,j))
      enddo
      enddo

      return
      end subroutine prbase
