c    ----------------
      subroutine cip
c    ----------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine cipadv(cadx,cady,da,dan,dadx,dady,dadxn,dadyn,
     &u,v,is,ie,js,je)
            double precision,dimension(:),intent(in) :: cadx,cady
            double precision,dimension(:,:),intent(in) :: da,dadx,dady,
     &u,v
            double precision,dimension(:,:),intent(inout) :: dan,dadxn,
     &dadyn
            integer :: is,ie,js,je
         end subroutine cipadv

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

         subroutine convertxm(dam,da,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: dam
            double precision,dimension(:,:),intent(inout) :: da
            integer :: is,ie,js,je
         end subroutine convertxm

         subroutine converty(da,dam,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: da
            double precision,dimension(:,:),intent(inout) :: dam
            integer :: is,ie,js,je
         end subroutine converty

         subroutine convertym(dam,da,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: dam
            double precision,dimension(:,:),intent(inout) :: da
            integer :: is,ie,js,je
         end subroutine convertym

         subroutine cipdsrc2x(da,dan,dadx,dadxn,dl,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: da,dan,dadx
            double precision,dimension(:,:),intent(out) :: dadxn
            double precision,dimension(:),intent(in) :: dl
            integer :: is,ie,js,je
         end subroutine cipdsrc2x

         subroutine cipdsrc2y(da,dan,dady,dadyn,dl,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: da,dan,dady
            double precision,dimension(:,:),intent(out) :: dadyn
            double precision,dimension(:),intent(in) :: dl
            integer :: is,ie,js,je
         end subroutine cipdsrc2y

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

         subroutine spheat(pre,tem,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: pre,tem
            integer :: is,ie,js,je
         end subroutine spheat

         subroutine eos(den,pre,tem,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: pre,tem
            double precision,dimension(:,:),intent(out) :: den
            integer :: is,ie,js,je
         end subroutine eos

         subroutine viscos(pre,tem,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: pre,tem
            integer :: is,ie,js,je
         end subroutine viscos

         subroutine artvis(den,DIV,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: den,DIV
            integer :: is,ie,js,je
         end subroutine artvis

         subroutine tcond(pre,tem,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: pre,tem
            integer :: is,ie,js,je
         end subroutine tcond
      end interface

      dimension :: DIV(1:i0,1:j0),DIVh(1:i0,1:j0)
      dimension :: averoxi(1:i0,1:j0),averoyi(1:i0,1:j0)

C Advection Phase

C Pressure
C Log convert
      do j=1,j0
      do i=1,i0
         prdx(i,j) = prdx(i,j)/pr(i,j)
         prdy(i,j) = prdy(i,j)/pr(i,j)
         pr(i,j) = dlog(pr(i,j))
      enddo
      enddo

      call cipadv(dx,dy,pr,prh,prdx,prdy,prdxh,prdyh,vx,vy,
     &2,i0-1,2,j0-1)

C Exp convert
      do j=2,j0-1
      do i=2,i0-1
         prh(i,j) = dexp(prh(i,j))
         prdxh(i,j) = prh(i,j)*prdxh(i,j)
         prdyh(i,j) = prh(i,j)*prdyh(i,j)
      enddo
      enddo

C Temperature
C Log convert
      do j=1,j0
      do i=1,i0
         tedx(i,j) = tedx(i,j)/te(i,j)
         tedy(i,j) = tedy(i,j)/te(i,j)
         te(i,j) = dlog(te(i,j))
      enddo
      enddo

      call cipadv(dx,dy,te,teh,tedx,tedy,tedxh,tedyh,vx,vy,
     &2,i0-1,2,j0-1)

C Exp convert
      do j=2,j0-1
      do i=2,i0-1
         teh(i,j) = dexp(teh(i,j))
         tedxh(i,j) = teh(i,j)*tedxh(i,j)
         tedyh(i,j) = teh(i,j)*tedyh(i,j)
      enddo
      enddo

      call convertx(vy,vymx,1,i0-1,1,j0)
      call converty(vx,vxmy,1,i0,1,j0-1)

C Vx
      call cipadv(dxm,dy,vxm,vxmh,vxdxm,vxdym,vxdxmh,vxdymh,vxm,vymx,
     &2,i0-2,2,j0-1)
C Vy
      call cipadv(dx,dym,vym,vymh,vydxm,vydym,vydxmh,vydymh,vxmy,vym,
     &2,i0-1,2,j0-2)
C Vz
      call cipadv(dx,dy,vz,vzh,vzdx,vzdy,vzdxh,vzdyh,vx,vy,
     &2,i0-1,2,j0-1)

      call chck(prh,pmin,2,i0-1,2,j0-1)
      call chck(teh,tmin,2,i0-1,2,j0-1)

      call convertxm(vxmh,vxh,3,i0-2,2,j0-1)
      call convertym(vymh,vyh,2,i0-1,3,j0-2)

      call convertx(vyh,vymx,2,i0-2,3,j0-2)
      call converty(vxh,vxmy,3,i0-2,2,j0-2)

      call cipdsrc

C Diffusion Phase
      call spheat(prh,teh,2,i0-1,2,j0-1)
      call eos(roh,prh,teh,2,i0-1,2,j0-1)

      do j=3,j0-2
      do i=3,i0-2
         delvx =-vxmh(i-1,j) + vxmh(i,j)
         delvy =-vymh(i,j-1) + vymh(i,j)
         DIV(i,j) = delvx*dxmi(i-1) + delvy*dymi(j-1)
      enddo
      enddo

C Viscosity Stress Tensor
      call viscos(prh,teh,2,i0-1,2,j0-1)
      call artvis(roh,DIV,3,i0-2,3,j0-2)

C Thermal Conductivity
      call tcond(prh,teh,2,i0-1,2,j0-1)

      do j=2,j0-1
      do i=2,i0-2
         qxm(i,j) =-kappaxm(i,j)*(-teh(i,j) + teh(i+1,j))*dxi(i)
      enddo
      enddo

      do j=2,j0-2
      do i=2,i0-1
         qym(i,j) =-kappaym(i,j)*(-teh(i,j) + teh(i,j+1))*dyi(j)
      enddo
      enddo

C (i,j)
      do j=3,j0-2
      do i=3,i0-2
         delvx =-vxmh(i-1,j) + vxmh(i,j)
         vst(i,j,1,1) = (rlambda(i,j) + rlambdaa(i,j))*DIV(i,j) 
     &                + 2.d0*(rmu(i,j) + rmua(i,j))*delvx*dxmi(i-1)
      enddo
      enddo

C (i,j)
      do j=3,j0-2
      do i=3,i0-2
         delvy =-vymh(i,j-1) + vymh(i,j)
         vst(i,j,2,2) = (rlambda(i,j) + rlambdaa(i,j))*DIV(i,j)
     &                + 2.d0*(rmu(i,j) + rmua(i,j))*delvy*dymi(j-1)
      enddo
      enddo

C (i+1/2,j+1/2)
      do j=3,j0-3
      do i=3,i0-3
         avemu = 0.25d0*
     &(( rmu(i,j) +  rmu(i+1,j) +  rmu(i,j+1) +  rmu(i+1,j+1)) + 
     & (rmua(i,j) + rmua(i+1,j) + rmua(i,j+1) + rmua(i+1,j+1)))
         delvx =-vxmh(i,j) + vxmh(i,j+1)
         delvy =-vymh(i,j) + vymh(i+1,j)

         vst(i,j,1,2) = avemu*(delvx*dyi(j) + delvy*dxi(i))
         vst(i,j,2,1) = vst(i,j,1,2)
      enddo
      enddo

C (i+1/2,j)
      do j=3,j0-2
      do i=3,i0-3
         avemu = 0.5d0*(( rmu(i,j) +  rmu(i+1,j)) + 
     &                  (rmua(i,j) + rmua(i+1,j)))
         delvz =-vzh(i,j) + vzh(i+1,j)

         vst(i,j,3,1) = avemu*delvz*dxi(i)
      enddo
      enddo

C (i,j+1/2)
      do j=3,j0-3
      do i=3,i0-2
         avemu = 0.5d0*(( rmu(i,j) +  rmu(i,j+1)) + 
     &                  (rmua(i,j) + rmua(i,j+1)))
         delvz =-vzh(i,j) + vzh(i,j+1)

         vst(i,j,3,2) = avemu*delvz*dyi(j)
      enddo
      enddo

      do j=2,j0-1
      do i=2,i0-2
         averoxi(i,j) = 2.d0/(roh(i,j) + roh(i+1,j))
      enddo
      enddo

      do j=2,j0-2
      do i=2,i0-1
         averoyi(i,j) = 2.d0/(roh(i,j) + roh(i,j+1))
      enddo
      enddo

      do j=4,j0-3
      do i=3,i0-3
         vxmh2(i,j) = vxmh(i,j) + dt*averoxi(i,j)*
     &((-vst(i,j,1,1) + vst(i+1,j,1,1))*dxi(i) + 
     & (-vst(i,j-1,1,2) + vst(i,j,1,2))*dymi(j-1))
      enddo
      enddo

      do j=3,j0-3
      do i=4,i0-3
         vymh2(i,j) = vymh(i,j) + dt*averoyi(i,j)*
     &((-vst(i-1,j,2,1) + vst(i,j,2,1))*dxmi(i-1) + 
     & (-vst(i,j,2,2) + vst(i,j+1,2,2))*dyi(j))
      enddo
      enddo

      do j=4,j0-3
      do i=4,i0-3
         vzh2(i,j) = vzh(i,j) + dt/roh(i,j)*
     &((-vst(i-1,j,3,1) + vst(i,j,3,1))*dxmi(i-1) + 
     & (-vst(i,j-1,3,2) + vst(i,j,3,2))*dymi(j-1))
      enddo
      enddo

      call convertxm(vxmh2,vxh2,4,i0-3,4,j0-3)
      call convertym(vymh2,vyh2,4,i0-3,4,j0-3)

      do j=5,j0-4
      do i=5,i0-4
C (i-1/2,j)
         vst11xm1 = 0.5d0*(vst(i-1,j,1,1) + vst(i,j,1,1))
C (i+1/2,j)
         vst11xm2 = 0.5d0*(vst(i,j,1,1) + vst(i+1,j,1,1))
C (i-1/2,j)
         vst12xm1 = 0.5d0*(vst(i-1,j-1,1,2) + vst(i-1,j,1,2))
C (i+1/2,j)
         vst12xm2 = 0.5d0*(vst(i,j-1,1,2) + vst(i,j,1,2))

C (i,j-1/2)
         vst22ym1 = 0.5d0*(vst(i,j-1,2,2) + vst(i,j,2,2))
C (i,j+1/2)
         vst22ym2 = 0.5d0*(vst(i,j,2,2) + vst(i,j+1,2,2))
C (i,j-1/2)
         vst21ym1 = 0.5d0*(vst(i-1,j-1,2,1) + vst(i,j-1,2,1))
C (i,j+1/2)
         vst21ym2 = 0.5d0*(vst(i-1,j,2,1) + vst(i,j,2,1))

C (i-1/2,j)
         vyxm1 = 0.5d0*(vyh2(i-1,j) + vyh2(i,j))
C (i+1/2,j)
         vyxm2 = 0.5d0*(vyh2(i,j) + vyh2(i+1,j))
C (i-1/2,j)
         vzxm1 = 0.5d0*(vzh2(i-1,j) + vzh2(i,j))
C (i+1/2,j)
         vzxm2 = 0.5d0*(vzh2(i,j) + vzh2(i+1,j))

C (i,j-1/2)
         vxym1 = 0.5d0*(vxh2(i,j-1) + vxh2(i,j))
C (i,j+1/2)
         vxym2 = 0.5d0*(vxh2(i,j) + vxh2(i,j+1))
C (i,j-1/2)
         vzym1 = 0.5d0*(vzh2(i,j-1) + vzh2(i,j))
C (i,j+1/2)
         vzym2 = 0.5d0*(vzh2(i,j) + vzh2(i,j+1))

         v21 =  vxh(i,j)**2 +  vyh(i,j)**2 +  vzh(i,j)**2
         v22 = vxh2(i,j)**2 + vyh2(i,j)**2 + vzh2(i,j)**2

         htheta = (-vst11xm1*vxmh2(i-1,j) + vst11xm2*vxmh2(i,j)
     &            - vst12xm1*vyxm1 + vst12xm2*vyxm2
     &            - vst(i-1,j,3,1)*vzxm1 + vst(i,j,3,1)*vzxm2
     &            + qxm(i-1,j) - qxm(i,j))*dxmi(i-1)
     &          + (-vst21ym1*vxym1 + vst21ym2*vxym2 
     &            - vst22ym1*vymh2(i,j-1) + vst22ym2*vymh2(i,j)
     &            - vst(i,j-1,3,2)*vzym1 + vst(i,j,3,2)*vzym2
     &            + qym(i,j-1) - qym(i,j))*dymi(j-1)
     &          - 0.5d0*dti*roh(i,j)*(v22 - v21)

         teh2(i,j) = teh(i,j) + dt*htheta/(roh(i,j)*Cv(i,j))
         prh2(i,j) = prh(i,j) + dt*(gm(i,j) - 1.d0)
     &             /(1.d0 + roh(i,j)*Cp(i,j)*rmuj(i,j))*htheta
      enddo
      enddo

      call cipdsrc2x(vxmh,vxmh2,vxdxmn,vxdxmh2,dxm,4,i0-4,4,j0-3)
      call cipdsrc2x(vymh,vymh2,vydxmn,vydxmh2,dx,5,i0-4,3,j0-3)
      call cipdsrc2x(vzh,vzh2,vzdxn,vzdxh2,dx,5,i0-4,4,j0-3)
      call cipdsrc2x(prh,prh2,prdxn,prdxh2,dx,6,i0-5,5,j0-4)
      call cipdsrc2x(teh,teh2,tedxn,tedxh2,dx,6,i0-5,5,j0-4)

      call cipdsrc2y(vxmh,vxmh2,vxdymn,vxdymh2,dy,3,i0-3,5,j0-4)
      call cipdsrc2y(vymh,vymh2,vydymn,vydymh2,dym,4,i0-3,4,j0-4)
      call cipdsrc2y(vzh,vzh2,vzdyn,vzdyh2,dy,4,i0-3,5,j0-4)
      call cipdsrc2y(prh,prh2,prdyn,prdyh2,dy,5,i0-4,6,j0-5)
      call cipdsrc2y(teh,teh2,tedyn,tedyh2,dy,5,i0-4,6,j0-5)

C Acoustic Phase
      call prbase

      call cipdsrc2x(vxmh2,vxmn,vxdxmh2,vxdxmn,dxm,7,i0-7,6,j0-5)
      call cipdsrc2x(vymh2,vymn,vydxmh2,vydxmn,dx,7,i0-6,6,j0-6)
      call cipdsrc2x(vzh2,vzn,vzdxh2,vzdxn,dx,6,i0-5,5,j0-4)
      call cipdsrc2x(prh2,prn,prdxh2,prdxn,dx,7,i0-6,6,j0-5)
      call cipdsrc2x(teh2,ten,tedxh2,tedxn,dx,7,i0-6,6,j0-5)

      call cipdsrc2y(vxmh2,vxmn,vxdymh2,vxdymn,dy,6,i0-6,7,j0-6)
      call cipdsrc2y(vymh2,vymn,vydymh2,vydymn,dym,6,i0-5,7,j0-7)
      call cipdsrc2y(vzh2,vzn,vzdyh2,vzdyn,dy,5,i0-4,6,j0-5)
      call cipdsrc2y(prh2,prn,prdyh2,prdyn,dy,6,i0-5,7,j0-6)
      call cipdsrc2y(teh2,ten,tedyh2,tedyn,dy,6,i0-5,7,j0-6)

      call chck(prn,pmin,6,i0-5,6,j0-5)
      call chck(ten,tmin,6,i0-5,6,j0-5)

C update
      do j=7,j0-6
      do i=7,i0-6
         vz(i,j) = vzn(i,j)
         pr(i,j) = prn(i,j)
         te(i,j) = ten(i,j)
         vzdx(i,j) = vzdxn(i,j)
         prdx(i,j) = prdxn(i,j)
         tedx(i,j) = tedxn(i,j)
         vzdy(i,j) = vzdyn(i,j)
         prdy(i,j) = prdyn(i,j)
         tedy(i,j) = tedyn(i,j)
      enddo
      enddo

      do j=7,j0-6
      do i=7,i0-7
         vxm(i,j) = vxmn(i,j)
         vxdxm(i,j) = vxdxmn(i,j)
         vxdym(i,j) = vxdymn(i,j)
      enddo
      enddo

      do j=7,j0-7
      do i=7,i0-6
         vym(i,j) = vymn(i,j)
         vydxm(i,j) = vydxmn(i,j)
         vydym(i,j) = vydymn(i,j)
      enddo
      enddo

      call bndc

      call convertxm(vxm,vx,2,i0-1,1,j0)
      call convertym(vym,vy,1,i0,2,j0-1)

      if(nbnd1 .eq. 1) then
         call bdfrdx(vx,vx0,1,0,0)
      elseif(nbnd1 .eq. 2) then
         call bdfrex(vx,1,0,0)
      elseif(nbnd1 .eq. 3) then
         call bdperx(vx,1,0,0)
      elseif(nbnd1 .eq. 4) then
         call bdsymx(vx,1,0,0,-1.d0)
      elseif(nbnd1 .eq. 5) then
         call bdconx(vx,1,0,0,0.d0)
      elseif(nbnd1 .eq. 6) then
         call bdconx(vx,1,0,0,0.d0)
      endif

      if(nbnd2 .eq. 1) then
         call bdfrdx(vx,vx0,1,1,0)
      elseif(nbnd2 .eq. 2) then
         call bdfrex(vx,1,1,0)
      elseif(nbnd2 .eq. 3) then
         call bdperx(vx,1,1,0)
      elseif(nbnd2 .eq. 4) then
         call bdsymx(vx,1,1,0,-1.d0)
      elseif(nbnd2 .eq. 5) then
         call bdconx(vx,1,1,0,0.d0)
      elseif(nbnd2 .eq. 6) then
         call bdconx(vx,1,1,0,0.d0)
      endif

      if(nbnd3 .eq. 1) then
         call bdfrdy(vy,vy0,1,0,0)
      elseif(nbnd3 .eq. 2) then
         call bdfrey(vy,1,0,0)
      elseif(nbnd3 .eq. 3) then
         call bdpery(vy,1,0,0)
      elseif(nbnd3 .eq. 4) then
         call bdsymy(vy,1,0,0,-1.d0)
      elseif(nbnd3 .eq. 5) then
         call bdcony(vy,1,0,0,0.d0)
      elseif(nbnd3 .eq. 6) then
         call bdcony(vy,1,0,0,0.d0)
      endif

      if(nbnd4 .eq. 1) then
         call bdfrdy(vy,vy0,1,1,0)
      elseif(nbnd4 .eq. 2) then
         call bdfrey(vy,1,1,0)
      elseif(nbnd4 .eq. 3) then
         call bdpery(vy,1,1,0)
      elseif(nbnd4 .eq. 4) then
         call bdsymy(vy,1,1,0,-1.d0)
      elseif(nbnd4 .eq. 5) then
         call bdcony(vy,1,1,0,0.d0)
      elseif(nbnd4 .eq. 6) then
         call bdcony(vy,1,1,0,0.d0)
      endif

      call spheat(pr,te,1,i0,1,j0)
      call eos(ro,pr,te,1,i0,1,j0)

      return
      end subroutine cip
