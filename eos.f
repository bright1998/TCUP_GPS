c    -----------------------------------------
      subroutine eos(den,pre,tem,is,ie,js,je)
c    -----------------------------------------
C Equation of State
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(in) :: pre,tem
      double precision,dimension(:,:),intent(out) :: den
      integer :: is,ie,js,je

c den : Density
c pre : Pressure
c tem : Temperature
c Cs  : Sound Velocity
c rmuj: Joule Thomson Coefficient
c gm  : Ratio of Specific Heat
c rgc : Gas Constant

c Start of Editable Block -------------------------------------------
c Example: Ideal Gas Case
      do j=js,je
      do i=is,ie
         den(i,j) = rmw*pre(i,j)/(tem(i,j)*rgc)
      enddo
      enddo

      do j=js,je
      do i=is,ie
         Cs(i,j) = dsqrt(gm(i,j)*pre(i,j)/den(i,j))
      enddo
      enddo

      do j=js,je
      do i=is,ie
         rmuj(i,j) = 0.d0
      enddo
      enddo
c End of Editable Block ---------------------------------------------

      return
      end subroutine eos
