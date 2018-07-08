c    ------------------
      subroutine write
c    ------------------
      use common
      implicit double precision(a-h,o-z)

      character(len=3) cnm
c      character(len=12) fn
      character(len=13) fn

c      fn='out/'//cnm(nf)//'.data'
      fn='out/x'//cnm(nf)//'.data'
c      fn='out/y'//cnm(nf)//'.data'
      open(99,file=fn,form='formatted')
      j = (j0 + 1)/2
      do i=1,i0
         write(99,990,advance="NO")x(i),ro(i,j),vx(i,j),vy(i,j),vz(i,j),
     &pr(i,j),te(i,j)
         write(99,*)
      enddo
c      i = (i0 + 1)/2
c      do j=1,j0
c         write(99,990,advance="NO")y(j),ro(i,j),vx(i,j),vy(i,j),vz(i,j),
c     &pr(i,j),te(i,j)
c         write(99,*)
c      enddo
      close(99)

 990  FORMAT(7F24.15)

      write(6,*)' '
      write(6,*)
     &'wrote data:( nf=',nf,' ns=',ns,' tm=',tm,
     &') in "','*'//cnm(nf)//'.data','"'

      nf = nf + 1

      return
      end subroutine write
