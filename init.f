c    -----------------
      subroutine init
c    -----------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine chck(da,damin,is,ie,js,je)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: damin
            integer :: is,ie,js,je
         end subroutine chck

         subroutine convertx(da,dam,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: da
            double precision,dimension(:,:),intent(inout) :: dam
            integer :: is,ie,js,je
         end subroutine convertx

         subroutine converty(da,dam,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: da
            double precision,dimension(:,:),intent(inout) :: dam
            integer :: is,ie,js,je
         end subroutine converty

         subroutine spheat(pre,tem,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: pre,tem
            integer :: is,ie,js,je
         end subroutine spheat

         subroutine eos(den,pre,tem,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: pre,tem
            double precision,dimension(:,:),intent(out) :: den
            integer :: is,ie,js,je
         end subroutine eos
      end interface

      do j=1,j0
      do i=1,i0
         pr(i,j) = 1.d0*1.01325e+5
         te(i,j) = 283.15d0
         vx(i,j) = 0.d0
         vy(i,j) = 0.d0
         vz(i,j) = 0.d0
      enddo
      enddo

      call chck(pr,pmin,1,i0,1,j0)
      call chck(te,tmin,1,i0,1,j0)
      call convertx(vx,vxm,1,i0-1,1,j0)
      call converty(vy,vym,1,i0,1,j0-1)

      do j=1,j0
      do i=2,i0-1
         prdx(i,j) = (-pr(i-1,j) + pr(i+1,j))/(dx(i-1) + dx(i))
         tedx(i,j) = (-te(i-1,j) + te(i+1,j))/(dx(i-1) + dx(i))
         vzdx(i,j) = (-vz(i-1,j) + vz(i+1,j))/(dx(i-1) + dx(i))
      enddo
      enddo

      do j=1,j0
      do i=2,i0-2
         vxdxm(i,j) = (-vxm(i-1,j) + vxm(i+1,j))/(dxm(i-1) + dxm(i))
      enddo
      enddo

      do j=1,j0-1
      do i=2,i0-1
         vydxm(i,j) = (-vym(i-1,j) + vym(i+1,j))/(dx(i-1) + dx(i))
      enddo
      enddo

      do j=2,j0-1
      do i=1,i0
         prdy(i,j) = (-pr(i,j-1) + pr(i,j+1))/(dy(j-1) + dy(j))
         tedy(i,j) = (-te(i,j-1) + te(i,j+1))/(dy(j-1) + dy(j))
         vzdy(i,j) = (-vz(i,j-1) + vz(i,j+1))/(dy(j-1) + dy(j))
      enddo
      enddo

      do j=2,j0-1
      do i=1,i0-1
         vxdym(i,j) = (-vxm(i,j-1) + vxm(i,j+1))/(dy(j-1) + dy(j))
      enddo
      enddo

      do j=2,j0-2
      do i=1,i0
         vydym(i,j) = (-vym(i,j-1) + vym(i,j+1))/(dym(j-1) + dym(j))
      enddo
      enddo

      call bndc

      call spheat(pr,te,1,i0,1,j0)
      call eos(ro,pr,te,1,i0,1,j0)

C Preserve the initial value
      do j=1,j0
      do i=1,i0
         pr0(i,j) = pr(i,j)
         te0(i,j) = te(i,j)
         vx0(i,j) = vx(i,j)
         vy0(i,j) = vy(i,j)
         vz0(i,j) = vz(i,j)
         ro0(i,j) = ro(i,j)
      enddo
      enddo

      do j=1,j0
      do i=1,i0-1
         vxm0(i,j) = vxm(i,j)
      enddo
      enddo

      do j=1,j0-1
      do i=1,i0
         vym0(i,j) = vym(i,j)
      enddo
      enddo

      do j=1,j0
      do i=1,i0
         prdx0(i,j) = prdx(i,j)
         tedx0(i,j) = tedx(i,j)
         vzdx0(i,j) = vzdx(i,j)
         prdy0(i,j) = prdy(i,j)
         tedy0(i,j) = tedy(i,j)
         vzdy0(i,j) = vzdy(i,j)
      enddo
      enddo

      do j=1,j0
      do i=1,i0-1
         vxdxm0(i,j) = vxdxm(i,j)
         vxdym0(i,j) = vxdym(i,j)
      enddo
      enddo

      do j=1,j0-1
      do i=1,i0
         vydxm0(i,j) = vydxm(i,j)
         vydym0(i,j) = vydym(i,j)
      enddo
      enddo

      return
      end subroutine init
