c    ---------------
      module common
c    ---------------
      implicit double precision(a-h,o-z)

c main.f
      double precision :: tm,pi,pi4,pi4i,rgc
      integer :: ns,nf,nbcast

c prms.f
      integer :: nsmax,ic,i0,jc,j0,ncflconst,itemax,nrestart
      integer :: nbnd1,nbnd2,nbnd3,nbnd4
      integer :: ncp,nmu,nkap
      double precision :: xcmin,xcmax,rax
      double precision :: ycmin,ycmax,ray
      double precision :: tmmax,ftm,dtlim,cn,dtcon
      double precision :: epsit
      double precision :: pmin,tmin,lenga
      double precision :: nT_w1,nT_w2,nT_w3,nT_w4
      double precision :: T_wall1,T_wall2,T_wall3,T_wall4
      double precision :: rmw,concp,conmu,rmu0,tref,stem,conkap

c cflc.f
      double precision :: dt

c prbase.f
      double precision :: omega

C Arraies on Normal Grids
      double precision,allocatable :: x(:),dx(:),dxi(:)
      double precision,allocatable :: y(:),dy(:),dyi(:)

      double precision,allocatable ::  ro(:,:), ro0(:,:),roh(:,:)
      double precision,allocatable ::  vx(:,:), vx0(:,:),
     &                                vxh(:,:),vxh2(:,:)
      double precision,allocatable ::  vy(:,:), vy0(:,:),
     &                                vyh(:,:),vyh2(:,:)
      double precision,allocatable ::  vz(:,:), vz0(:,:),
     &                                vzh(:,:),vzh2(:,:),vzn(:,:)
      double precision,allocatable ::  pr(:,:), pr0(:,:),
     &                                prh(:,:),prh2(:,:),prn(:,:)
      double precision,allocatable ::  te(:,:), te0(:,:),
     &                                teh(:,:),teh2(:,:),ten(:,:)

C Arraies on Staggered Grids
      double precision,allocatable :: xm(:),dxm(:),dxmi(:)
      double precision,allocatable :: ym(:),dym(:),dymi(:)

      double precision,allocatable ::  vxm(:,:), vxm0(:,:),
     &                                vxmh(:,:),vxmh2(:,:),vxmn(:,:)
      double precision,allocatable ::  vym(:,:), vym0(:,:),
     &                                vymh(:,:),vymh2(:,:),vymn(:,:)

C For CIP Method
      double precision,allocatable ::  vxdxm(:,:), vxdxm0(:,:),
     &                                vxdxmh(:,:),vxdxmh2(:,:),
     &                                vxdxmn(:,:)
      double precision,allocatable ::  vydxm(:,:), vydxm0(:,:),
     &                                vydxmh(:,:),vydxmh2(:,:),
     &                                vydxmn(:,:)
      double precision,allocatable ::   vzdx(:,:),  vzdx0(:,:),
     &                                 vzdxh(:,:), vzdxh2(:,:),
     &                                 vzdxn(:,:)
      double precision,allocatable ::   prdx(:,:),  prdx0(:,:),
     &                                 prdxh(:,:), prdxh2(:,:),
     &                                 prdxn(:,:)
      double precision,allocatable ::   tedx(:,:),  tedx0(:,:),
     &                                 tedxh(:,:), tedxh2(:,:),
     &                                 tedxn(:,:)

      double precision,allocatable ::  vxdym(:,:), vxdym0(:,:),
     &                                vxdymh(:,:),vxdymh2(:,:),
     &                                vxdymn(:,:)
      double precision,allocatable ::  vydym(:,:), vydym0(:,:),
     &                                vydymh(:,:),vydymh2(:,:),
     &                                vydymn(:,:)
      double precision,allocatable ::   vzdy(:,:),  vzdy0(:,:),
     &                                 vzdyh(:,:), vzdyh2(:,:),
     &                                 vzdyn(:,:)
      double precision,allocatable ::   prdy(:,:),  prdy0(:,:),
     &                                 prdyh(:,:), prdyh2(:,:),
     &                                 prdyn(:,:)
      double precision,allocatable ::   tedy(:,:),  tedy0(:,:),
     &                                 tedyh(:,:), tedyh2(:,:),
     &                                 tedyn(:,:)

      double precision,allocatable :: vxmy(:,:),vymx(:,:)

      double precision,allocatable :: Cp(:,:),Cv(:,:),gm(:,:)
C Coefficients of Polynomial for Specific Heat,Viscosity,Thermal Conductivity
      double precision,allocatable :: coCp(:),comu(:),cokap(:)
C Sound Velocity
      double precision,allocatable :: Cs(:,:)
C Viscosity
      double precision,allocatable :: rmu(:,:),rlambda(:,:)
C Artificial Viscosity
      double precision,allocatable :: rmua(:,:),rlambdaa(:,:)
C Joule Thomson Coefficient
      double precision,allocatable :: rmuj(:,:)
C Viscosity Stress Tensor
      double precision,allocatable :: vst(:,:,:,:)

C Source Term
      double precision,allocatable :: Fexx(:,:),Fexy(:,:),Fexz(:,:)
      double precision,allocatable :: Fexxm(:,:),Fexym(:,:)

C Thermal Conduction
      double precision,allocatable :: kappaxm(:,:),kappaym(:,:)
      double precision,allocatable :: qxm(:,:),qym(:,:)

      end module common
