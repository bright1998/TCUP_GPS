c    ----------------------------------------------------------
      subroutine gradient(ain,adx,ady,aout1,aout2,is,ie,js,je)
c    ----------------------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(in) :: ain
      double precision,dimension(:),intent(in) :: adx,ady
      double precision,dimension(:,:),intent(out) :: aout1,aout2
      integer :: is,ie,js,je

      do j=js,je
      do i=is+1,ie-1
         aout1(i,j) = (-ain(i-1,j) + ain(i+1,j))/(adx(i-1) + adx(i))
      enddo
      enddo

      do j=js+1,je-1
      do i=is,ie
         aout2(i,j) = (-ain(i,j-1) + ain(i,j+1))/(ady(j-1) + ady(j))
      enddo
      enddo

      return
      end subroutine gradient
