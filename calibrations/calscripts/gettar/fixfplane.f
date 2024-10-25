
      real xifu(100),yifu(100)
      character aifu(100)*3,a1*3,a4*3,a5*3,a6*3

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
      do i=1,6
         read(1,*)
      enddo
      open(unit=11,file='out',status='unknown')
      do i=1,4
         write(11,"(a7)") "# Empty"
      enddo
      write(11,"(a56)") 
     $     "# IFUSLOT X_FP Y_FP SPECID SPECSLOT IFUID IFUROT PLATESC"
      write(11,1101) "000",0.,0.,"900","900","900",0.,1.
      do i=1,100
         read(1,*,end=667) a1,x2,x3,a4,a5,a6
         do j=1,nfp
            if(aifu(j).eq.a1) then
               write(11,1101) a1,xifu(j),yifu(j),a4,a5,a6,0.,1.
               goto 888
            endif
         enddo
         print *,"Not in fplane.ifu0: ",a1
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
