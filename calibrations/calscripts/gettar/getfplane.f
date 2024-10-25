
      real xifu(100),yifu(100)
      character aifu(100)*3,as*34,ass*3,aid*3,a1*30,ai*3

      open(unit=1,file='fplane.ifu0',status='old')
      nfp=0
      do i=1,100
         read(1,*,end=666) a1,x2,x3
         nfp=nfp+1
         aifu(nfp)=a1
         xifu(nfp)=x2
         yifu(nfp)=x3
      enddo
 666  continue
      close(1)

      open(unit=1,file='in',status='old')
      open(unit=11,file='out',status='unknown')
      do i=1,4
         write(11,"(a7)") "# Empty"
      enddo
      write(11,"(a56)") 
     $     "# IFUSLOT X_FP Y_FP SPECID SPECSLOT IFUID IFUROT PLATESC"
      write(11,1101) "000",0.,0.,"900","900","900",0.,1.
      do i=1,100
         read(1,*,end=667) a1,i2,i3,i4,i5
         ai="000"
         if(i2.lt.100) write(ai(2:3),1002) i2
         if(i2.ge.100) write(ai(1:3),1003) i2
         as="000"
         if(i3.lt.10) write(as(3:3),1001) i3
         if(i3.ge.10.and.i3.lt.100) write(as(2:3),1002) i3
         if(i3.ge.100) write(as(1:3),1003) i3
         ass="000"
         if(i4.lt.10) write(ass(3:3),1001) i4
         if(i4.ge.10.and.i4.lt.100) write(ass(2:3),1002) i4
         if(i4.ge.100) write(ass(1:3),1003) i4
         aid="000"
         if(i5.lt.10) write(aid(3:3),1001) i5
         if(i5.ge.10.and.i5.lt.100) write(aid(2:3),1002) i5
         if(i5.ge.100) write(aid(1:3),1003) i5
         do j=1,nfp
            if(aifu(j).eq.ai) then
               write(11,1101) ai,xifu(j),yifu(j),as,ass,aid,0.,1.
               goto 888
            endif
         enddo
         print *,"Not in fplane.ifu0: ",ai
 888     continue
      enddo
 667  continue
      close(1)
      close(11)

 1001 format(i1)
 1002 format(i2)
 1003 format(i3)
 1101 format(a3,2(1x,f8.3),3(2x,a3),2(1x,f8.3))

      end
