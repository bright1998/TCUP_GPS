c    --------------
      program main
c    --------------
      use common
      implicit double precision(a-h,o-z)

      write(6,*)'*****************************
     &***************************************'
      write(6,*)'*** Start Calculation by GPS 
     &(Graphically Performable Simulator) ***'
      write(6,*)'***               Solving by 
     &Thermo CIP-CUP Scheme               ***'
      write(6,*)'*****************************
     &***************************************'
      write(6,*)' '

      call prms

      call grid

      tm = 0.d0
      ns = 0
      nf = 0
      nbcast = 0

      pi = 4.d0*DATAN(1.d0)
      pi4 = 4.d0*pi
      pi4i = 1.d0/pi4
C Gas Constant [J/K/kmol]
      rgc = 8.31447e+3
      write(6,*)'Gas Constant= ',rgc,'[J/K/kmol]'

      call init
C Source Term
      call source

      call pgwrite
      call write
      if(nrestart .eq. 1) then
         call read
      endif

      write(6,*)' '
      write(6,*)'Maximum step number'
      write(6,*)'nsmax= ',nsmax
      write(6,*)' '
      write(6,*)'============ Start Calculation ============'

C Main Loop
 50   continue

      call cflc

      write(6,*)'ns= ',ns,' tm= ',tm,' dt= ',dt
      ns = ns + 1
      tm = tm + dt

      if((ncflconst .eq. 0) .and. (dt .le. dtlim)) goto 500

      call cip

      if(tm .gt. (ftm*dfloat(nf))) then
         call pgwrite
         call write
      endif
      if(ns .ge. nsmax) goto 100
      if(tm .ge. tmmax) goto 600

      goto 50

 100  write(6,*)'$$ stopped due to come to maximum step number',nsmax
      goto 600

 500  write(6,*)' '
      write(6,*)'!! abnormal stopped due to small dt, less than ',dtlim
      goto 200

 600  call restart

 200  continue

      stop
      end program
