c    ----------------------------------------
      subroutine spheat(pre,tem,is,ie,js,je)
c    ----------------------------------------
C Calculate Specific Heat
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(in) :: pre,tem
      integer :: is,ie,js,je

C Specific Heat at Constant Pressure [J/kg/K]
      if(ncp .eq. 0) then
         do j=js,je
         do i=is,ie
            Cp(i,j) = concp
         enddo
         enddo
      elseif(ncp .eq. 1) then
         do j=js,je
         do i=is,ie
            Cp(i,j) = concp
            do k=1,5
               Cp(i,j) = Cp(i,j) + coCp(k)*tem(i,j)**k
            enddo
         enddo
         enddo
      elseif(ncp .eq. 2) then
         do j=js,je
         do i=is,ie
c Start of Editable Block -------------------------------------------
c            Cp(i,j) = concp
c End of Editable Block ---------------------------------------------
         enddo
         enddo
      endif

C Specific Heat at Constant Volume [J/kg/K]
      do j=js,je
      do i=is,ie
C Mayer's Relation
         Cv(i,j) = Cp(i,j) - rgc/rmw
      enddo
      enddo

C Ratio of Specific Heat
      do j=js,je
      do i=is,ie
         gm(i,j) = Cp(i,j)/Cv(i,j)
      enddo
      enddo

      return
      end subroutine spheat

c    ----------------------------------------
      subroutine viscos(pre,tem,is,ie,js,je)
c    ----------------------------------------
C Calculate Viscosity
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(in) :: pre,tem
      integer :: is,ie,js,je

      if(nmu .eq. 0) then
         do j=js,je
         do i=is,ie
            rmu(i,j) = conmu
         enddo
         enddo
      elseif(nmu .eq. 1) then
         do j=js,je
         do i=is,ie
            rmu(i,j) = conmu
            do k=1,5
               rmu(i,j) = rmu(i,j) + comu(k)*tem(i,j)**k
            enddo
         enddo
         enddo
      elseif(nmu .eq. 2) then
C Air ---Sutherland's Law ---
C rmu0: [Pa s] !Viscosity at T=tref, P=101.325[kPa]
C tref: [K] !Reference Temperature
C stem: [K] Sutherland Temperature
         do j=js,je
         do i=is,ie
            rmu(i,j) = (tem(i,j)/tref)**1.5*
     &                 (tref + stem)/(tem(i,j) + stem)*rmu0
         enddo
         enddo
      elseif(nmu .eq. 3) then
         do j=js,je
         do i=is,ie
c Start of Editable Block -------------------------------------------
C Water
c            rmu(i,j) = 1.e-3*dexp(-1.5668d0 + 230.298/
c     &                                        (tem(i,j) - 146.797))
c End of Editable Block ---------------------------------------------
         enddo
         enddo
      endif

c Stokes relation
      do j=js,je
      do i=is,ie
         rlambda(i,j) = -2.d0/3.d0*rmu(i,j)
      enddo
      enddo

      return
      end subroutine viscos

c    ---------------------------------------
      subroutine tcond(pre,tem,is,ie,js,je)
c    ---------------------------------------
C Calculate Thermal Conductivity
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(in) :: pre,tem
      integer :: is,ie,js,je

c Temperature [K] at (i+1/2,j) & (i,j+1/2)
      dimension :: texm(1:i0,1:j0),teym(1:i0,1:j0)

c (i+1/2,j)
      do j=js,je
      do i=is,ie-1
         texm(i,j) = 0.5d0*(tem(i,j) + tem(i+1,j))
      enddo
      enddo

c (i,j+1/2)
      do j=js,je-1
      do i=is,ie
         teym(i,j) = 0.5d0*(tem(i,j) + tem(i,j+1))
      enddo
      enddo

C kappaxm: Thermal Conductivity [W/m/K] at (i+1/2,j)
C kappaym: Thermal Conductivity [W/m/K] at (i,j+1/2)
      if(nkap .eq. 0) then
         do j=js,je
         do i=is,ie-1
            kappaxm(i,j) = conkap
         enddo
         enddo
         do j=js,je-1
         do i=is,ie
            kappaym(i,j) = conkap
         enddo
         enddo
      elseif(nkap .eq. 1) then
         do j=js,je
         do i=is,ie-1
            kappaxm(i,j) = conkap
            do k=1,5
               kappaxm(i,j) = kappaxm(i,j) + cokap(k)*texm(i,j)**k
            enddo
         enddo
         enddo
         do j=js,je-1
         do i=is,ie
            kappaym(i,j) = conkap
            do k=1,5
               kappaym(i,j) = kappaym(i,j) + cokap(k)*teym(i,j)**k
            enddo
         enddo
         enddo
      elseif(nkap .eq. 2) then
         do j=js,je
         do i=is,ie-1
c Start of Editable Block -------------------------------------------
c            kappaxm(i,j) = conkap
c End of Editable Block ---------------------------------------------
         enddo
         enddo
         do j=js,je-1
         do i=is,ie
c Start of Editable Block -------------------------------------------
c            kappaym(i,j) = conkap
c End of Editable Block ---------------------------------------------
         enddo
         enddo
      endif

      return
      end subroutine tcond

