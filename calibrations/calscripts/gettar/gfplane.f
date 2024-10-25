
      parameter(nmax=1000000)
      integer ima(nmax)
      character file1*100,cifu*3,cspec*3,cnum*3,cexp*5
      character cfa(nmax)*3,csa(nmax)*3

      open(unit=1,file='listall',status='old')

      nt=0
      do i=1,1000
         read(1,*,end=666) file1
         read(file1(40:47),'(i8)') imth
         open(unit=2,file=file1,status='old')
c         print *,file1
         do j=1,6
            read(2,*)
         enddo
         do j=1,100
            read(2,*,end=667) cifu,x2,x3,cspec
            nt=nt+1
c            print *,nt,imth,cifu," ",cspec
            ima(nt)=imth
            cfa(nt)=cifu
            csa(nt)=cspec
         enddo
 667     continue
         close(2)
      enddo
 666  continue
      close(1)

c rback 20191201 011 exp01 026 318 201912 1
c 20191101 010 exp01 013

      open(unit=1,file="j1",status='old')
      open(unit=11,file="out",status='unknown')
      do i=1,1000000
         read(1,*,end=668) i1,cnum,cexp,cifu
c         print *,i1
         do j=1,nt
            if(ima(j).eq.i1.and.cfa(j).eq.cifu) then
               write(file1,"(i8)") i1
               write(11,1101) i1,cnum,cexp,cifu,csa(j),file1(1:6)
c               write(*,1101) i1,cnum,cexp,cifu,csa(j),file1(1:6)
            endif
         enddo
      enddo
 668  continue
      close(1)

 1101 format("rback ",i8,1x,a3,1x,a5,1x,a3,1x,a3,1x,a6,1x,"1")

      end
