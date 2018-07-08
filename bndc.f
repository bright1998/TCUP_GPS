c    -----------------
      subroutine bndc
c    -----------------
      use common
      implicit double precision(a-h,o-z)

      interface
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

      if(nbnd1 .eq. 1) then
         call bdfrdx(vxm,vxm0,6,0,1)
         call bdinix(vxdxm,vxdxm0,6,0,1)
         call bdfrdx(vxdym,vxdym0,6,0,1)

         call bdfrdx(vym,vym0,6,0,0)
         call bdinix(vydxm,vydxm0,6,0,0)
         call bdfrdx(vydym,vydym0,6,0,0)

         call bdfrdx(vz,vz0,6,0,0)
         call bdinix(vzdx,vzdx0,6,0,0)
         call bdfrdx(vzdy,vzdy0,6,0,0)

         call bdfrdx(pr,pr0,6,0,0)
         call bdinix(prdx,prdx0,6,0,0)
         call bdfrdx(prdy,prdy0,6,0,0)

         call bdfrdx(te,te0,6,0,0)
         call bdinix(tedx,tedx0,6,0,0)
         call bdfrdx(tedy,tedy0,6,0,0)
      elseif(nbnd1 .eq. 2) then
         call bdfrex(vxm,6,0,1)
         call bdconx(vxdxm,6,0,1,0.d0)
         call bdfrex(vxdym,6,0,1)

         call bdfrex(vym,6,0,0)
         call bdconx(vydxm,6,0,0,0.d0)
         call bdfrex(vydym,6,0,0)

         call bdfrex(vz,6,0,0)
         call bdconx(vzdx,6,0,0,0.d0)
         call bdfrex(vzdy,6,0,0)

         call bdfrex(pr,6,0,0)
         call bdconx(prdx,6,0,0,0.d0)
         call bdfrex(prdy,6,0,0)

         call bdfrex(te,6,0,0)
         call bdconx(tedx,6,0,0,0.d0)
         call bdfrex(tedy,6,0,0)
      elseif(nbnd1 .eq. 3) then
         call bdperx(vxm,6,0,1)
         call bdperx(vxdxm,6,0,1)
         call bdperx(vxdym,6,0,1)

         call bdperx(vym,6,0,0)
         call bdperx(vydxm,6,0,0)
         call bdperx(vydym,6,0,0)

         call bdperx(vz,6,0,0)
         call bdperx(vzdx,6,0,0)
         call bdperx(vzdy,6,0,0)

         call bdperx(pr,6,0,0)
         call bdperx(prdx,6,0,0)
         call bdperx(prdy,6,0,0)

         call bdperx(te,6,0,0)
         call bdperx(tedx,6,0,0)
         call bdperx(tedy,6,0,0)
      elseif(nbnd1 .eq. 4) then
         call bdsymx(vxm,6,0,1,-1.d0)
         call bdsymx(vxdxm,6,0,1, 1.d0)
         call bdsymx(vxdym,6,0,1,-1.d0)

         call bdsymx(vym,6,0,0, 1.d0)
         call bdsymx(vydxm,6,0,0,-1.d0)
         call bdsymx(vydym,6,0,0, 1.d0)

         call bdsymx(vz,6,0,0, 1.d0)
         call bdsymx(vzdx,6,0,0,-1.d0)
         call bdsymx(vzdy,6,0,0, 1.d0)

         call bdsymx(pr,6,0,0, 1.d0)
         call bdsymx(prdx,6,0,0,-1.d0)
         call bdsymx(prdy,6,0,0, 1.d0)

         call bdsymx(te,6,0,0, 1.d0)
         call bdsymx(tedx,6,0,0,-1.d0)
         call bdsymx(tedy,6,0,0, 1.d0)
      elseif(nbnd1 .eq. 5) then
         call bdconx(vxm,6,0,1,0.d0)
         call bdconx(vxdxm,6,0,1,0.d0)
         call bdconx(vxdym,6,0,1,0.d0)

         call bdfrex(vym,6,0,0)
         call bdconx(vydxm,6,0,0,0.d0)
         call bdfrex(vydym,6,0,0)

         call bdfrex(vz,6,0,0)
         call bdconx(vzdx,6,0,0,0.d0)
         call bdfrex(vzdy,6,0,0)

         call bdfrex(pr,6,0,0)
         call bdconx(prdx,6,0,0,0.d0)
         call bdfrex(prdy,6,0,0)

         if(nT_w1 .eq. 0) then
            call bdconx(te,6,0,0,T_wall1)
            call bdconx(tedx,6,0,0,0.d0)
            call bdconx(tedy,6,0,0,0.d0)
         elseif(nT_w1 .eq. 1) then
            call bdfrex(te,6,0,0)
            call bdconx(tedx,6,0,0,0.d0)
            call bdfrex(tedy,6,0,0)
         endif
      elseif(nbnd1 .eq. 6) then
         call bdconx(vxm,6,0,1,0.d0)
         call bdconx(vxdxm,6,0,1,0.d0)
         call bdconx(vxdym,6,0,1,0.d0)

         call bdconx(vym,6,0,0,0.d0)
         call bdconx(vydxm,6,0,0,0.d0)
         call bdconx(vydym,6,0,0,0.d0)

         call bdconx(vz,6,0,0,0.d0)
         call bdconx(vzdx,6,0,0,0.d0)
         call bdconx(vzdy,6,0,0,0.d0)

         call bdfrex(pr,6,0,0)
         call bdconx(prdx,6,0,0,0.d0)
         call bdfrex(prdy,6,0,0)

         if(nT_w1 .eq. 0) then
            call bdconx(te,6,0,0,T_wall1)
            call bdconx(tedx,6,0,0,0.d0)
            call bdconx(tedy,6,0,0,0.d0)
         elseif(nT_w1 .eq. 1) then
            call bdfrex(te,6,0,0)
            call bdconx(tedx,6,0,0,0.d0)
            call bdfrex(tedy,6,0,0)
         endif
      endif

      if(nbnd2 .eq. 1) then
         call bdfrdx(vxm,vxm0,6,1,1)
         call bdinix(vxdxm,vxdxm0,6,1,1)
         call bdfrdx(vxdym,vxdym0,6,1,1)

         call bdfrdx(vym,vym0,6,1,0)
         call bdinix(vydxm,vydxm0,6,1,0)
         call bdfrdx(vydym,vydym0,6,1,0)

         call bdfrdx(vz,vz0,6,1,0)
         call bdinix(vzdx,vzdx0,6,1,0)
         call bdfrdx(vzdy,vzdy0,6,1,0)

         call bdfrdx(pr,pr0,6,1,0)
         call bdinix(prdx,prdx0,6,1,0)
         call bdfrdx(prdy,prdy0,6,1,0)

         call bdfrdx(te,te0,6,1,0)
         call bdinix(tedx,tedx0,6,1,0)
         call bdfrdx(tedy,tedy0,6,1,0)
      elseif(nbnd2 .eq. 2) then
         call bdfrex(vxm,6,1,1)
         call bdconx(vxdxm,6,1,1,0.d0)
         call bdfrex(vxdym,6,1,1)

         call bdfrex(vym,6,1,0)
         call bdconx(vydxm,6,1,0,0.d0)
         call bdfrex(vydym,6,1,0)

         call bdfrex(vz,6,1,0)
         call bdconx(vzdx,6,1,0,0.d0)
         call bdfrex(vzdy,6,1,0)

         call bdfrex(pr,6,1,0)
         call bdconx(prdx,6,1,0,0.d0)
         call bdfrex(prdy,6,1,0)

         call bdfrex(te,6,1,0)
         call bdconx(tedx,6,1,0,0.d0)
         call bdfrex(tedy,6,1,0)
      elseif(nbnd2 .eq. 3) then
         call bdperx(vxm,6,1,1)
         call bdperx(vxdxm,6,1,1)
         call bdperx(vxdym,6,1,1)

         call bdperx(vym,6,1,0)
         call bdperx(vydxm,6,1,0)
         call bdperx(vydym,6,1,0)

         call bdperx(vz,6,1,0)
         call bdperx(vzdx,6,1,0)
         call bdperx(vzdy,6,1,0)

         call bdperx(pr,6,1,0)
         call bdperx(prdx,6,1,0)
         call bdperx(prdy,6,1,0)

         call bdperx(te,6,1,0)
         call bdperx(tedx,6,1,0)
         call bdperx(tedy,6,1,0)
      elseif(nbnd2 .eq. 4) then
         call bdsymx(vxm,6,1,1,-1.d0)
         call bdsymx(vxdxm,6,1,1, 1.d0)
         call bdsymx(vxdym,6,1,1,-1.d0)

         call bdsymx(vym,6,1,0, 1.d0)
         call bdsymx(vydxm,6,1,0,-1.d0)
         call bdsymx(vydym,6,1,0, 1.d0)

         call bdsymx(vz,6,1,0, 1.d0)
         call bdsymx(vzdx,6,1,0,-1.d0)
         call bdsymx(vzdy,6,1,0, 1.d0)

         call bdsymx(pr,6,1,0, 1.d0)
         call bdsymx(prdx,6,1,0,-1.d0)
         call bdsymx(prdy,6,1,0, 1.d0)

         call bdsymx(te,6,1,0, 1.d0)
         call bdsymx(tedx,6,1,0,-1.d0)
         call bdsymx(tedy,6,1,0, 1.d0)
      elseif(nbnd2 .eq. 5) then
         call bdconx(vxm,6,1,1,0.d0)
         call bdconx(vxdxm,6,1,1,0.d0)
         call bdconx(vxdym,6,1,1,0.d0)

         call bdfrex(vym,6,1,0)
         call bdconx(vydxm,6,1,0,0.d0)
         call bdfrex(vydym,6,1,0)

         call bdfrex(vz,6,1,0)
         call bdconx(vzdx,6,1,0,0.d0)
         call bdfrex(vzdy,6,1,0)

         call bdfrex(pr,6,1,0)
         call bdconx(prdx,6,1,0,0.d0)
         call bdfrex(prdy,6,1,0)

         if(nT_w2 .eq. 0) then
            call bdconx(te,6,1,0,T_wall2)
            call bdconx(tedx,6,1,0,0.d0)
            call bdconx(tedy,6,1,0,0.d0)
         elseif(nT_w2 .eq. 1) then
            call bdfrex(te,6,1,0)
            call bdconx(tedx,6,1,0,0.d0)
            call bdfrex(tedy,6,1,0)
         endif
      elseif(nbnd2 .eq. 6) then
         call bdconx(vxm,6,1,1,0.d0)
         call bdconx(vxdxm,6,1,1,0.d0)
         call bdconx(vxdym,6,1,1,0.d0)

         call bdconx(vym,6,1,0,0.d0)
         call bdconx(vydxm,6,1,0,0.d0)
         call bdconx(vydym,6,1,0,0.d0)

         call bdconx(vz,6,1,0,0.d0)
         call bdconx(vzdx,6,1,0,0.d0)
         call bdconx(vzdy,6,1,0,0.d0)

         call bdfrex(pr,6,1,0)
         call bdconx(prdx,6,1,0,0.d0)
         call bdfrex(prdy,6,1,0)

         if(nT_w2 .eq. 0) then
            call bdconx(te,6,1,0,T_wall2)
            call bdconx(tedx,6,1,0,0.d0)
            call bdconx(tedy,6,1,0,0.d0)
         elseif(nT_w2 .eq. 1) then
            call bdfrex(te,6,1,0)
            call bdconx(tedx,6,1,0,0.d0)
            call bdfrex(tedy,6,1,0)
         endif
      endif

      if(nbnd3 .eq. 1) then
         call bdfrdy(vxm,vxm0,6,0,0)
         call bdfrdy(vxdxm,vxdxm0,6,0,0)
         call bdiniy(vxdym,vxdym0,6,0,0)

         call bdfrdy(vym,vym0,6,0,1)
         call bdfrdy(vydxm,vydxm0,6,0,1)
         call bdiniy(vydym,vydym0,6,0,1)

         call bdfrdy(vz,vz0,6,0,0)
         call bdfrdy(vzdx,vzdx0,6,0,0)
         call bdiniy(vzdy,vzdy0,6,0,0)

         call bdfrdy(pr,pr0,6,0,0)
         call bdfrdy(prdx,prdx0,6,0,0)
         call bdiniy(prdy,prdy0,6,0,0)

         call bdfrdy(te,te0,6,0,0)
         call bdfrdy(tedx,tedx0,6,0,0)
         call bdiniy(tedy,tedy0,6,0,0)
      elseif(nbnd3 .eq. 2) then
         call bdfrey(vxm,6,0,0)
         call bdfrey(vxdxm,6,0,0)
         call bdcony(vxdym,6,0,0,0.d0)

         call bdfrey(vym,6,0,1)
         call bdfrey(vydxm,6,0,1)
         call bdcony(vydym,6,0,1,0.d0)

         call bdfrey(vz,6,0,0)
         call bdfrey(vzdx,6,0,0)
         call bdcony(vzdy,6,0,0,0.d0)

         call bdfrey(pr,6,0,0)
         call bdfrey(prdx,6,0,0)
         call bdcony(prdy,6,0,0,0.d0)

         call bdfrey(te,6,0,0)
         call bdfrey(tedx,6,0,0)
         call bdcony(tedy,6,0,0,0.d0)
      elseif(nbnd3 .eq. 3) then
         call bdpery(vxm,6,0,0)
         call bdpery(vxdxm,6,0,0)
         call bdpery(vxdym,6,0,0)

         call bdpery(vym,6,0,1)
         call bdpery(vydxm,6,0,1)
         call bdpery(vydym,6,0,1)

         call bdpery(vz,6,0,0)
         call bdpery(vzdx,6,0,0)
         call bdpery(vzdy,6,0,0)

         call bdpery(pr,6,0,0)
         call bdpery(prdx,6,0,0)
         call bdpery(prdy,6,0,0)

         call bdpery(te,6,0,0)
         call bdpery(tedx,6,0,0)
         call bdpery(tedy,6,0,0)
      elseif(nbnd3 .eq. 4) then
         call bdsymy(vxm,6,0,0, 1.d0)
         call bdsymy(vxdxm,6,0,0, 1.d0)
         call bdsymy(vxdym,6,0,0,-1.d0)

         call bdsymy(vym,6,0,1,-1.d0)
         call bdsymy(vydxm,6,0,1,-1.d0)
         call bdsymy(vydym,6,0,1, 1.d0)

         call bdsymy(vz,6,0,0, 1.d0)
         call bdsymy(vzdx,6,0,0, 1.d0)
         call bdsymy(vzdy,6,0,0,-1.d0)

         call bdsymy(pr,6,0,0, 1.d0)
         call bdsymy(prdx,6,0,0, 1.d0)
         call bdsymy(prdy,6,0,0,-1.d0)

         call bdsymy(te,6,0,0, 1.d0)
         call bdsymy(tedx,6,0,0, 1.d0)
         call bdsymy(tedy,6,0,0,-1.d0)
      elseif(nbnd3 .eq. 5) then
         call bdfrey(vxm,6,0,0)
         call bdfrey(vxdxm,6,0,0)
         call bdcony(vxdym,6,0,0,0.d0)

         call bdcony(vym,6,0,1,0.d0)
         call bdcony(vydxm,6,0,1,0.d0)
         call bdcony(vydym,6,0,1,0.d0)

         call bdfrey(vz,6,0,0)
         call bdfrey(vzdx,6,0,0)
         call bdcony(vzdy,6,0,0,0.d0)

         call bdfrey(pr,6,0,0)
         call bdfrey(prdx,6,0,0)
         call bdcony(prdy,6,0,0,0.d0)

         if(nT_w3 .eq. 0) then
            call bdcony(te,6,0,0,T_wall3)
            call bdcony(tedx,6,0,0,0.d0)
            call bdcony(tedy,6,0,0,0.d0)
         elseif(nT_w3 .eq. 1) then
            call bdfrey(te,6,0,0)
            call bdfrey(tedx,6,0,0)
            call bdcony(tedy,6,0,0,0.d0)
         endif
      elseif(nbnd3 .eq. 6) then
         call bdcony(vxm,6,0,0,0.d0)
         call bdcony(vxdxm,6,0,0,0.d0)
         call bdcony(vxdym,6,0,0,0.d0)

         call bdcony(vym,6,0,1,0.d0)
         call bdcony(vydxm,6,0,1,0.d0)
         call bdcony(vydym,6,0,1,0.d0)

         call bdcony(vz,6,0,0,0.d0)
         call bdcony(vzdx,6,0,0,0.d0)
         call bdcony(vzdy,6,0,0,0.d0)

         call bdfrey(pr,6,0,0)
         call bdfrey(prdx,6,0,0)
         call bdcony(prdy,6,0,0,0.d0)

         if(nT_w3 .eq. 0) then
            call bdcony(te,6,0,0,T_wall3)
            call bdcony(tedx,6,0,0,0.d0)
            call bdcony(tedy,6,0,0,0.d0)
         elseif(nT_w3 .eq. 1) then
            call bdfrey(te,6,0,0)
            call bdfrey(tedx,6,0,0)
            call bdcony(tedy,6,0,0,0.d0)
         endif
      endif

      if(nbnd4 .eq. 1) then
         call bdfrdy(vxm,vxm0,6,1,0)
         call bdfrdy(vxdxm,vxdxm0,6,1,0)
         call bdiniy(vxdym,vxdym0,6,1,0)

         call bdfrdy(vym,vym0,6,1,1)
         call bdfrdy(vydxm,vydxm0,6,1,1)
         call bdiniy(vydym,vydym0,6,1,1)

         call bdfrdy(vz,vz0,6,1,0)
         call bdfrdy(vzdx,vzdx0,6,1,0)
         call bdiniy(vzdy,vzdy0,6,1,0)

         call bdfrdy(pr,pr0,6,1,0)
         call bdfrdy(prdx,prdx0,6,1,0)
         call bdiniy(prdy,prdy0,6,1,0)

         call bdfrdy(te,te0,6,1,0)
         call bdfrdy(tedx,tedx0,6,1,0)
         call bdiniy(tedy,tedy0,6,1,0)
      elseif(nbnd4 .eq. 2) then
         call bdfrey(vxm,6,1,0)
         call bdfrey(vxdxm,6,1,0)
         call bdcony(vxdym,6,1,0,0.d0)

         call bdfrey(vym,6,1,1)
         call bdfrey(vydxm,6,1,1)
         call bdcony(vydym,6,1,1,0.d0)

         call bdfrey(vz,6,1,0)
         call bdfrey(vzdx,6,1,0)
         call bdcony(vzdy,6,1,0,0.d0)

         call bdfrey(pr,6,1,0)
         call bdfrey(prdx,6,1,0)
         call bdcony(prdy,6,1,0,0.d0)

         call bdfrey(te,6,1,0)
         call bdfrey(tedx,6,1,0)
         call bdcony(tedy,6,1,0,0.d0)
      elseif(nbnd4 .eq. 3) then
         call bdpery(vxm,6,1,0)
         call bdpery(vxdxm,6,1,0)
         call bdpery(vxdym,6,1,0)

         call bdpery(vym,6,1,1)
         call bdpery(vydxm,6,1,1)
         call bdpery(vydym,6,1,1)

         call bdpery(vz,6,1,0)
         call bdpery(vzdx,6,1,0)
         call bdpery(vzdy,6,1,0)

         call bdpery(pr,6,1,0)
         call bdpery(prdx,6,1,0)
         call bdpery(prdy,6,1,0)

         call bdpery(te,6,1,0)
         call bdpery(tedx,6,1,0)
         call bdpery(tedy,6,1,0)
      elseif(nbnd4 .eq. 4) then
         call bdsymy(vxm,6,1,0, 1.d0)
         call bdsymy(vxdxm,6,1,0, 1.d0)
         call bdsymy(vxdym,6,1,0,-1.d0)

         call bdsymy(vym,6,1,1,-1.d0)
         call bdsymy(vydxm,6,1,1,-1.d0)
         call bdsymy(vydym,6,1,1, 1.d0)

         call bdsymy(vz,6,1,0, 1.d0)
         call bdsymy(vzdx,6,1,0, 1.d0)
         call bdsymy(vzdy,6,1,0,-1.d0)

         call bdsymy(pr,6,1,0, 1.d0)
         call bdsymy(prdx,6,1,0, 1.d0)
         call bdsymy(prdy,6,1,0,-1.d0)

         call bdsymy(te,6,1,0, 1.d0)
         call bdsymy(tedx,6,1,0, 1.d0)
         call bdsymy(tedy,6,1,0,-1.d0)
      elseif(nbnd4 .eq. 5) then
         call bdfrey(vxm,6,1,0)
         call bdfrey(vxdxm,6,1,0)
         call bdcony(vxdym,6,1,0,0.d0)

         call bdcony(vym,6,1,1,0.d0)
         call bdcony(vydxm,6,1,1,0.d0)
         call bdcony(vydym,6,1,1,0.d0)

         call bdfrey(vz,6,1,0)
         call bdfrey(vzdx,6,1,0)
         call bdcony(vzdy,6,1,0,0.d0)

         call bdfrey(pr,6,1,0)
         call bdfrey(prdx,6,1,0)
         call bdcony(prdy,6,1,0,0.d0)

         if(nT_w4 .eq. 0) then
            call bdcony(te,6,1,0,T_wall4)
            call bdcony(tedx,6,1,0,0.d0)
            call bdcony(tedy,6,1,0,0.d0)
         elseif(nT_w4 .eq. 1) then
            call bdfrey(te,6,1,0)
            call bdfrey(tedx,6,1,0)
            call bdcony(tedy,6,1,0,0.d0)
         endif
      elseif(nbnd4 .eq. 6) then
         call bdcony(vxm,6,1,0,0.d0)
         call bdcony(vxdxm,6,1,0,0.d0)
         call bdcony(vxdym,6,1,0,0.d0)

         call bdcony(vym,6,1,1,0.d0)
         call bdcony(vydxm,6,1,1,0.d0)
         call bdcony(vydym,6,1,1,0.d0)

         call bdcony(vz,6,1,0,0.d0)
         call bdcony(vzdx,6,1,0,0.d0)
         call bdcony(vzdy,6,1,0,0.d0)

         call bdfrey(pr,6,1,0)
         call bdfrey(prdx,6,1,0)
         call bdcony(prdy,6,1,0,0.d0)

         if(nT_w4 .eq. 0) then
            call bdcony(te,6,1,0,T_wall4)
            call bdcony(tedx,6,1,0,0.d0)
            call bdcony(tedy,6,1,0,0.d0)
         elseif(nT_w4 .eq. 1) then
            call bdfrey(te,6,1,0)
            call bdfrey(tedx,6,1,0)
            call bdcony(tedy,6,1,0,0.d0)
         endif
      endif

      return
      end subroutine bndc

c    -------------------------------------------
      subroutine bdfrdx(da,da0,margin,mbnd,men)
c    -------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision,dimension(:,:),intent(in) :: da0
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         ibnd = 1+margin
         do j=1,j0
         do i=1,margin
            da(ibnd-i,j) = (da0(ibnd-i,j) - da0(ibnd,j)) + da(ibnd,j)
         enddo
         enddo
      else
         ibnd = i0-margin-men
         do j=1,j0
         do i=1,margin
            da(ibnd+i,j) = (da0(ibnd+i,j) - da0(ibnd,j)) + da(ibnd,j)
         enddo
         enddo
      endif

      return
      end subroutine bdfrdx

c    ---------------------------------------
      subroutine bdfrex(da,margin,mbnd,men)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         do j=1,j0
         do i=1,margin
            da(i,j) = da(margin+1,j)
         enddo
         enddo
      else
         ibnd1 = i0+1-men
         do j=1,j0
         do i=1,margin
            da(ibnd1-i,j) = da(ibnd1-margin-1,j)
         enddo
         enddo
      endif

      return
      end subroutine bdfrex

c    ---------------------------------------
      subroutine bdperx(da,margin,mbnd,men)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         ibnd = 1+margin
         do j=1,j0
         do i=1,margin
            da(ibnd-i,j) = da(i0-margin-i,j)
         enddo
         enddo
      else
         ibnd = i0-margin-men
         do j=1,j0
         do i=1,margin
            if(men .eq. 0) then
               da(ibnd+i,j) = da(margin+1+i,j)
            else
               da(ibnd+i,j) = da(margin+i,j)
            endif
         enddo
         enddo
      endif

      return
      end subroutine bdperx

c    --------------------------------------------
      subroutine bdsymx(da,margin,mbnd,men,coff)
c    --------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision :: coff
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         ibnd1 = 1+margin !(4 for men=0,1)
         ibnd2 = 1+margin-men !(4 for men=0, 3 for men=1)
         do j=1,j0
         do i=1,margin
            da(ibnd1-i,j) = coff*da(ibnd2+i,j)
         enddo
         enddo
      else
         ibnd1 = i0-margin !(i0-3 for men=0,1)
         ibnd2 = i0-margin+men !(i0-3 for men=0, i0-2 for men=1)
         do j=1,j0
         do i=1,margin
            da(ibnd1+i,j) = coff*da(ibnd2-i,j)
         enddo
         enddo
      endif

      end subroutine bdsymx

c    ---------------------------------------------
      subroutine bdconx(da,margin,mbnd,men,const)
c    ---------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision :: const
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         do j=1,j0
         do i=1,margin
            da(i,j) = const
         enddo
         enddo
      else
         ibnd1 = i0+1-men
         do j=1,j0
         do i=1,margin
            da(ibnd1-i,j) = const
         enddo
         enddo
      endif

      return
      end subroutine bdconx

c    -------------------------------------------
      subroutine bdinix(da,da0,margin,mbnd,men)
c    -------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision,dimension(:,:),intent(in) :: da0
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         ibnd = 1+margin
         do j=1,j0
         do i=1,margin
            da(ibnd-i,j) = da0(ibnd-i,j)
         enddo
         enddo
      else
         ibnd = i0-margin-men
         do j=1,j0
         do i=1,margin
            da(ibnd+i,j) = da0(ibnd+i,j)
         enddo
         enddo
      endif

      return
      end subroutine bdinix

c    -------------------------------------------
      subroutine bdfrdy(da,da0,margin,mbnd,men)
c    -------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision,dimension(:,:),intent(in) :: da0
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         jbnd = 1+margin
         do j=1,margin
         do i=1,i0
            da(i,jbnd-j) = (da0(i,jbnd-j) - da0(i,jbnd)) + da(i,jbnd)
         enddo
         enddo
      else
         jbnd = j0-margin-men
         do j=1,margin
         do i=1,i0
            da(i,jbnd+j) = (da0(i,jbnd+j) - da0(i,jbnd)) + da(i,jbnd)
         enddo
         enddo
      endif

      return
      end subroutine bdfrdy

c    ---------------------------------------
      subroutine bdfrey(da,margin,mbnd,men)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         do j=1,margin
         do i=1,i0
            da(i,j) = da(i,margin+1)
         enddo
         enddo
      else
         jbnd1 = j0+1-men
         do j=1,margin
         do i=1,i0
            da(i,jbnd1-j) = da(i,jbnd1-margin-1)
         enddo
         enddo
      endif

      return
      end subroutine bdfrey

c    ---------------------------------------
      subroutine bdpery(da,margin,mbnd,men)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         jbnd = 1+margin
         do j=1,margin
            do i=1,i0
               da(i,jbnd-j) = da(i,j0-margin-j)
            enddo
         enddo
      else
         jbnd = j0-margin-men
         do j=1,margin
            do i=1,i0
               if(men .eq. 0) then
                  da(i,jbnd+j) = da(i,margin+1+j)
               else
                  da(i,jbnd+j) = da(i,margin+j)
               endif
            enddo
         enddo
      endif

      return
      end subroutine bdpery

c    --------------------------------------------
      subroutine bdsymy(da,margin,mbnd,men,coff)
c    --------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision :: coff
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         jbnd1 = 1+margin
         jbnd2 = 1+margin-men
         do j=1,margin
            do i=1,i0
               da(i,jbnd1-j) = coff*da(i,jbnd2+j)
            enddo
         enddo
      else
         jbnd1 = j0-margin
         jbnd2 = j0-margin+men
         do j=1,margin
            do i=1,i0
               da(i,jbnd1+j) = coff*da(i,jbnd2-j)
            enddo
         enddo
      endif

      end subroutine bdsymy

c    ---------------------------------------------
      subroutine bdcony(da,margin,mbnd,men,const)
c    ---------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision :: const
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         do j=1,margin
         do i=1,i0
            da(i,j) = const
         enddo
         enddo
      else
         jbnd1 = j0+1-men
         do j=1,margin
         do i=1,i0
            da(i,jbnd1-j) = const
         enddo
         enddo
      endif

      return
      end subroutine bdcony

c    -------------------------------------------
      subroutine bdiniy(da,da0,margin,mbnd,men)
c    -------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision,dimension(:,:),intent(in) :: da0
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         jbnd = 1+margin
         do j=1,margin
         do i=1,i0
            da(i,jbnd-j) = da0(i,jbnd-j)
         enddo
         enddo
      else
         jbnd = j0-margin-men
         do j=1,margin
         do i=1,i0
            da(i,jbnd+j) = da0(i,jbnd+j)
         enddo
         enddo
      endif

      return
      end subroutine bdiniy
