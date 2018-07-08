c    -----------------
      subroutine read
c    -----------------
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
      end interface

C Read the data for Restart Calculation
      open(22,file='out/restart.data',form='unformatted')
      read(22)tm,ns,nf,nbcast
      read(22)
     &((  vx(i,j),i=1,i0),j=1,j0),
     &((  vy(i,j),i=1,i0),j=1,j0),
     &((  vz(i,j),i=1,i0),j=1,j0),
     &((vzdx(i,j),i=1,i0),j=1,j0),
     &((vzdy(i,j),i=1,i0),j=1,j0),
     &((  pr(i,j),i=1,i0),j=1,j0),
     &((prdx(i,j),i=1,i0),j=1,j0),
     &((prdy(i,j),i=1,i0),j=1,j0),
     &((  te(i,j),i=1,i0),j=1,j0),
     &((tedx(i,j),i=1,i0),j=1,j0),
     &((tedy(i,j),i=1,i0),j=1,j0)
      read(22)
     &((  vxm(i,j),i=1,i0-1),j=1,j0),
     &((vxdxm(i,j),i=1,i0-1),j=1,j0),
     &((vxdym(i,j),i=1,i0-1),j=1,j0)
      read(22)
     &((  vym(i,j),i=1,i0),j=1,j0-1),
     &((vydxm(i,j),i=1,i0),j=1,j0-1),
     &((vydym(i,j),i=1,i0),j=1,j0-1)
      close(22)

      write(6,*)' '
      write(6,*)
     &'read data for restart calculation'

      call spheat(pr,te,1,i0,1,j0)
      call eos(ro,pr,te,1,i0,1,j0)

      return
      end subroutine read
