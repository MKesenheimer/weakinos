      subroutine jirocx(ric,foc,gc,vmass,vwidth , jiroc)
c
c This subroutine computes an off-shell vector current from an anti-
c Rarita Schwinger fermion and an antifermion. The vector boson propagator is given in Feynman gauge
c for a massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex foc(6)         : flow-out antifermion               <foc|
c       complex ric(18)        : flow-in antiRSfermion              |ric>
c       complex gc(2)          : coupling constants                  gfrv
c       real    vmass          : mass  of OUTPUT vector v
c       real    vwidth         : width of OUTPUT vector v
c
c output:
c       complex jiroc(6)         : vector current          j^mu(<foc|v|ric>)
c   
c- by Yoshitaro Takaesu - 2011/03/26
c  
      implicit none
      double complex foc(6),ric(18),gc(2),jiroc(6),c0,c1,c2,c3,cs,d
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex gcc(2)
      double complex cm2        ! mass**2- I Gamma mass (Fabio)


      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

      double complex sq1(2,2),sq2(2,2),sqg1(0:3,2,2),sqg2(0:3,2,2)
      double complex gsq1(0:3,2,2),gsq2(0:3,2,2)
      double complex rc(0:3,4),sqgr(0:3,0:3,4),gsqgr(0:3,4)
      integer i,j

#ifdef HELAS_CHECK
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      if ( abs(foc(1))+abs(foc(2))+abs(foc(3))+abs(foc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : foc in jirocx is zero spinor'
      endif
      if ( abs(foc(5))+abs(foc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : foc in jirocx has zero momentum'
      endif
      if (abs(ric(1))+abs(ric(2))+abs(ric(3))+abs(ric(4))+abs(ric(5))
     & +abs(ric(6))+abs(ric(7))+abs(ric(8))abs(ric(9))+abs(ric(10))
     & +abs(ric(11))+abs(ric(12))+abs(ric(13))+abs(ric(14))+abs(ric(15))
     & +abs(ric(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ric in jirocx is zero spinor'
      endif
      if ( abs(ric(17))+abs(ric(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ric in jirocx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in jirocx is zero coupling'
      endif
      if ( vmass.lt.rZero ) then
         write(stdo,*) ' helas-error : vmass in jirocx is negative'
         write(stdo,*) '             : vmass = ',vmass
      endif
      if ( vwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : vwidth in jirocx is negative'
         write(stdo,*) '             : vwidth = ',vwidth
      endif
#endif

      jiroc(5) = foc(5)-ric(17)
      jiroc(6) = foc(6)-ric(18)

      q(0) = dble( jiroc(5))
      q(1) = dble( jiroc(6))
      q(2) = dimag(jiroc(6))
      q(3) = dimag(jiroc(5))
      
      
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

#ifdef HELAS_CHECK
      if ( abs(jiroc(5))+abs(jiroc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : jiroc in jirocx has zero momentum'
      endif
      if ( vwidth.eq.rZero .and. q2.eq.vm2 ) then
         write(stdo,*)
     &        ' helas-error : jiroc in jirocx is on vmass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         jiroc(1) = cZero
         jiroc(2) = cZero
         jiroc(3) = cZero
         jiroc(4) = cZero
         return
      endif
#endif

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1))     

      rc(0,1) = ric(1)
      rc(0,2) = ric(2)
      rc(0,3) = ric(3)
      rc(0,4) = ric(4)
      rc(1,1) = ric(5)
      rc(1,2) = ric(6)
      rc(1,3) = ric(7)
      rc(1,4) = ric(8)
      rc(2,1) = ric(9)
      rc(2,2) = ric(10)
      rc(2,3) = ric(11)
      rc(2,4) = ric(12)
      rc(3,1) = ric(13)
      rc(3,2) = ric(14)
      rc(3,3) = ric(15)
      rc(3,4) = ric(16)

      
      sq1(1,1) = q(0)-q(3)
      sq1(1,2) = -(q(1)-cImag*q(2))
      sq1(2,1) = -(q(1)+cImag*q(2))
      sq1(2,2) = q(0)+q(3)

      sq2(1,1) = q(0)+q(3)
      sq2(1,2) = q(1)-cImag*q(2)
      sq2(2,1) = q(1)+cImag*q(2)
      sq2(2,2) = q(0)-q(3)

      sqg1(0,1,1) = sq1(1,1)
      sqg1(0,1,2) = sq1(1,2)
      sqg1(0,2,1) = sq1(2,1)
      sqg1(0,2,2) = sq1(2,2)
      sqg1(1,1,1) = -sq1(1,2)
      sqg1(1,1,2) = -sq1(1,1)
      sqg1(1,2,1) = -sq1(2,2)
      sqg1(1,2,2) = -sq1(2,1)
      sqg1(2,1,1) = -cImag*sq1(1,2)
      sqg1(2,1,2) = cImag*sq1(1,1)
      sqg1(2,2,1) = -cImag*sq1(2,2)
      sqg1(2,2,2) = cImag*sq1(2,1)
      sqg1(3,1,1) = -sq1(1,1)
      sqg1(3,1,2) = sq1(1,2)
      sqg1(3,2,1) = -sq1(2,1)
      sqg1(3,2,2) = sq1(2,2)

      sqg2(0,1,1) = sq2(1,1)
      sqg2(0,1,2) = sq2(1,2)
      sqg2(0,2,1) = sq2(2,1)
      sqg2(0,2,2) = sq2(2,2)
      sqg2(1,1,1) = sq2(1,2)
      sqg2(1,1,2) = sq2(1,1)
      sqg2(1,2,1) = sq2(2,2)
      sqg2(1,2,2) = sq2(2,1)
      sqg2(2,1,1) = cImag*sq2(1,2)
      sqg2(2,1,2) = -cImag*sq2(1,1)
      sqg2(2,2,1) = cImag*sq2(2,2)
      sqg2(2,2,2) = -cImag*sq2(2,1)
      sqg2(3,1,1) = sq2(1,1)
      sqg2(3,1,2) = -sq2(1,2)
      sqg2(3,2,1) = sq2(2,1)
      sqg2(3,2,2) = -sq2(2,2)

      gsq1(0,1,1) = sq2(1,1)
      gsq1(0,1,2) = sq2(1,2)
      gsq1(0,2,1) = sq2(2,1)
      gsq1(0,2,2) = sq2(2,2)
      gsq1(1,1,1) = sq2(2,1)
      gsq1(1,1,2) = sq2(2,2)
      gsq1(1,2,1) = sq2(1,1)
      gsq1(1,2,2) = sq2(1,2)
      gsq1(2,1,1) = -cImag*sq2(2,1)
      gsq1(2,1,2) = -cImag*sq2(2,2)
      gsq1(2,2,1) = cImag*sq2(1,1)
      gsq1(2,2,2) = cImag*sq2(1,2)
      gsq1(3,1,1) = sq2(1,1)
      gsq1(3,1,2) = sq2(1,2)
      gsq1(3,2,1) = -sq2(2,1)
      gsq1(3,2,2) = -sq2(2,2)

      gsq2(0,1,1) = sq1(1,1)
      gsq2(0,1,2) = sq1(1,2)
      gsq2(0,2,1) = sq1(2,1)
      gsq2(0,2,2) = sq1(2,2)
      gsq2(1,1,1) = -sq1(2,1)
      gsq2(1,1,2) = -sq1(2,2)
      gsq2(1,2,1) = -sq1(1,1)
      gsq2(1,2,2) = -sq1(1,2)
      gsq2(2,1,1) = cImag*sq1(2,1)
      gsq2(2,1,2) = cImag*sq1(2,2)
      gsq2(2,2,1) = -cImag*sq1(1,1)
      gsq2(2,2,2) = -cImag*sq1(1,2)
      gsq2(3,1,1) = -sq1(1,1)
      gsq2(3,1,2) = -sq1(1,2)
      gsq2(3,2,1) = sq1(2,1)
      gsq2(3,2,2) = sq1(2,2)


      do i=0,3
       do j=0,3    
        sqgr(i,j,1) = rc(j,1)*(sqg1(i,1,1)-gsq1(i,1,1))+rc(j,2)
     &*(sqg1(i,1,2)-gsq1(i,1,2))
        sqgr(i,j,2) = rc(j,1)*(sqg1(i,2,1)-gsq1(i,2,1))+rc(j,2)
     &*(sqg1(i,2,2)-gsq1(i,2,2))   
        sqgr(i,j,3) = rc(j,3)*(sqg2(i,1,1)-gsq2(i,1,1))+rc(j,4)
     &*(sqg2(i,1,2)-gsq2(i,1,2))
        sqgr(i,j,4) = rc(j,3)*(sqg2(i,2,1)-gsq2(i,2,1))+rc(j,4)
     &*(sqg2(i,2,2)-gsq2(i,2,2))
       enddo
      enddo

      do i=0,3
       gsqgr(i,1) = sqgr(i,0,3)-sqgr(i,1,4)+cImag*sqgr(i,2,4)
     &              -sqgr(i,3,3) 
       gsqgr(i,2) = sqgr(i,0,4)-sqgr(i,1,3)-cImag*sqgr(i,2,3)
     &              +sqgr(i,3,4) 
       gsqgr(i,3) = sqgr(i,0,1)+sqgr(i,1,2)-cImag*sqgr(i,2,2)
     &              +sqgr(i,3,1)
       gsqgr(i,4) = sqgr(i,0,2)+sqgr(i,1,1)+cImag*sqgr(i,2,1)
     &              -sqgr(i,3,2) 
      enddo

      jiroc(1) = cZero
      jiroc(2) = cZero
      jiroc(3) = cZero
      jiroc(4) = cZero

      if ( vmass.ne.rZero ) then

       d = -rOne/dcmplx( q2-vm2, vmass*vwidth )

       if ( gcc(2).ne.cZero ) then
        do i=1,4
         jiroc(i) = dconjg(gcc(2))*d*(foc(1)*gsqgr(i-1,1)+foc(2)
     &             *gsqgr(i-1,2))
        enddo
       endif
      
       if ( gcc(1).ne.cZero ) then
        do i=1,4
         jiroc(i) = jiroc(i)+dconjg(gcc(1))*d*(foc(3)*gsqgr(i-1,3)
     &             +foc(4)*gsqgr(i-1,4))
        enddo
       endif

      else

       d = -dcmplx( rOne/q2, rZero )

      if ( gcc(2).ne.cZero ) then
        do i=1,4
         jiroc(i) = dconjg(gcc(2))*d*(foc(1)*gsqgr(i-1,1)+foc(2)
     &             *gsqgr(i-1,2))
        enddo
       endif
    
       if ( gcc(1).ne.cZero ) then
        do i=1,4
         jiroc(i) = jiroc(i)+dconjg(gcc(1))*d*(foc(3)*gsqgr(i-1,3)
     &            +foc(4)*gsqgr(i-1,4))
        enddo
       endif
         

      end if
c
      return
      end
