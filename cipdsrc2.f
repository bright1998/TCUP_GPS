c    --------------------------------------------------------
      subroutine cipdsrc2x(da,dan,dadx,dadxn,dl,is,ie,js,je)
c    --------------------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(in) :: da,dan,dadx
      double precision,dimension(:,:),intent(out) :: dadxn
      double precision,dimension(:),intent(in) :: dl
      integer :: is,ie,js,je

      do j=js,je
      do i=is,ie
         dadxn(i,j) = dadx(i,j) + 
     &(-(dan(i-1,j) - da(i-1,j)) + (dan(i+1,j) - da(i+1,j)))
     &/(dl(i-1) + dl(i))
      enddo
      enddo

      return
      end subroutine cipdsrc2x

c    --------------------------------------------------------
      subroutine cipdsrc2y(da,dan,dady,dadyn,dl,is,ie,js,je)
c    --------------------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(in) :: da,dan,dady
      double precision,dimension(:,:),intent(out) :: dadyn
      double precision,dimension(:),intent(in) :: dl
      integer :: is,ie,js,je

      do j=js,je
      do i=is,ie
         dadyn(i,j) = dady(i,j) + 
     &(-(dan(i,j-1) - da(i,j-1)) + (dan(i,j+1) - da(i,j+1)))
     &/(dl(j-1) + dl(j))
      enddo
      enddo

      return
      end subroutine cipdsrc2y
