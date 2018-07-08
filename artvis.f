c    ----------------------------------------
      subroutine artvis(den,DIV,is,ie,js,je)
c    ----------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(in) :: den,DIV
      integer :: is,ie,js,je

C Initialization
      do j=js,je
      do i=is,ie
         rmua(i,j)     = 0.d0
         rlambdaa(i,j) = 0.d0
      enddo
      enddo

c      do j=js,je
c      do i=is,ie
c         vis = den(i,j)*(0.125d0*(gm(i,j) + 1.d0)*lenga)**2
c     &                 *dmin1(DIV(i,j),0.d0)
c         rmua(i,j)     =-0.75d0*vis
c         rlambdaa(i,j) =  0.5d0*vis
c      enddo
c      enddo

      return
      end subroutine artvis
