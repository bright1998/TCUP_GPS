c    -----------------
      subroutine prms
c    -----------------
      use common
      implicit double precision(a-h,o-z)

      allocate(coCp(1:5),comu(1:5),cokap(1:5))

      open(2,file='analysis.param',form='formatted')
      read(2,*)nsmax
      read(2,*)tmmax
      read(2,*)ftm
      read(2,*)ncflconst
c cn & dtlim are used if ncflconst = 0
      read(2,*)cn
      read(2,*)dtlim
c dtcon is used if ncflconst != 0
      read(2,*)dtcon
      read(2,*)nrestart

      open(3,file='grid.param',form='formatted')
      read(3,*)ic,jc
      read(3,*)xcmin,xcmax
      read(3,*)ycmin,ycmax
      read(3,*)i0,rax
      read(3,*)j0,ray

      open(5,file='thermo_prop.param',form='formatted')
C Molecular Weight [kg/kmol]
      read(5,*)rmw
C Specific Heat at Constant Pressure [J/kg/K]
C ncp = 0: Constant
C       1: Polynomial
C       2: User Definition
      read(5,*)ncp
      read(5,*)concp
      read(5,*)(coCp(i),i=1,5)
C Viscosity
C nmu = 0: Constant
C       1: Polynomial
C       2: Sutherland's Law
C       3: User Definition
      read(5,*)nmu
      read(5,*)conmu
      read(5,*)(comu(i),i=1,5)
      read(5,*)rmu0,tref,stem
C Thermal Conductivity [W/m/K]
C nkap = 0: Constant
C        1: Polynomial
C        2: User Definition
      read(5,*)nkap
      read(5,*)conkap
      read(5,*)(cokap(i),i=1,5)

      open(6,file='other.param',form='formatted')
      read(6,*)epsit
      read(6,*)itemax
C Over-relaxation Coefficient
      read(6,*)omega
C Lower Limit of Pressure and Temperature
      read(6,*)pmin,tmin
C Thickness Parameter of Capturing Shock Wave
      read(6,*)lenga

      open(7,file='bndc.param',form='formatted')
      read(7,*)nbnd1,nbnd2
      read(7,*)nbnd3,nbnd4

      if((nbnd1 .eq. 5) .or. (nbnd1 .eq. 6)) then
         open(8,file='wall1.param',form='formatted')
         read(8,*)nT_w1
         read(8,*)T_wall1
         close(8)
      endif
      if((nbnd2 .eq. 5) .or. (nbnd2 .eq. 6)) then
         open(8,file='wall2.param',form='formatted')
         read(8,*)nT_w2
         read(8,*)T_wall2
         close(8)
      endif
      if((nbnd3 .eq. 5) .or. (nbnd3 .eq. 6)) then
         open(8,file='wall3.param',form='formatted')
         read(8,*)nT_w3
         read(8,*)T_wall3
         close(8)
      endif
      if((nbnd4 .eq. 5) .or. (nbnd4 .eq. 6)) then
         open(8,file='wall4.param',form='formatted')
         read(8,*)nT_w4
         read(8,*)T_wall4
         close(8)
      endif

      call allocate

      close(2)
      close(3)
      close(5)
      close(6)
      close(7)

      return
      end subroutine prms
