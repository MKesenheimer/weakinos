      subroutine jiroxx(ri,fo,gc,vmass,vwidth , jiro)
c
c This subroutine computes an off-shell vector current from a Rarita-
c Schwinger fermion and a fermion. The vector boson propagator is given
c in Feynman gauge for a massless vector and in unitary gauge for a 
c massive vector.
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex ri(18)         : flow-in  RS fermion               |ri>v
c       complex gc(2)          : coupling constants                 GFRV
c       real    vmass          : mass  of OUTPUT vector v'(q)
c       real    vwidth         : width of OUTPUT vector v'(q)
c
c output:
c       complex jiro(6)        : vector current    j^mu(<fo|v[v',q]|ri>)
c   
c- by Y.Takaesu - 2010/12/16
c  
      implicit none
      double complex fo(6),ri(18),gc(2),jiro(6),c0,c1,c2,c3,cs,d
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex cm2 ! mass**2- I Gamma mass (Fabio)


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
      if ( abs(fo(1))+abs(fo(2))+abs(fo(3))+abs(fo(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fo in jiroxx is zero spinor'
      endif
      if ( abs(fo(5))+abs(fo(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in jiroxx has zero momentum'
      endif
      if (abs(ri(1))+abs(ri(2))+abs(ri(3))+abs(ri(4))+abs(ri(5))
     & +abs(ri(6))+abs(ri(7))+abs(ri(8))abs(ri(9))+abs(ri(10))
     & +abs(ri(11))+abs(ri(12))+abs(ri(13))+abs(ri(14))+abs(ri(15))
     & +abs(ri(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ri in jiroxx is zero spinor'
      endif
      if ( abs(ri(17))+abs(ri(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ri in jiroxx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in jiroxx is zero coupling'
      endif
      if ( vmass.lt.rZero ) then
         write(stdo,*) ' helas-error : vmass in jiroxx is negative'
         write(stdo,*) '             : vmass = ',vmass
      endif
      if ( vwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : vwidth in jiroxx is negative'
         write(stdo,*) '             : vwidth = ',vwidth
      endif
#endif

      jiro(5) = fo(5)-ri(17)
      jiro(6) = fo(6)-ri(18)

      q(0) = dble( jiro(5))
      q(1) = dble( jiro(6))
      q(2) = dimag(jiro(6))
      q(3) = dimag(jiro(5))
      
      
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

#ifdef HELAS_CHECK
      if ( abs(jiro(5))+abs(jiro(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : jiro in jiroxx has zero momentum'
      endif
      if ( vwidth.eq.rZero .and. q2.eq.vm2 ) then
         write(stdo,*)
     &        ' helas-error : jiro in jiroxx is on vmass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         jiro(1) = cZero
         jiro(2) = cZero
         jiro(3) = cZero
         jiro(4) = cZero
         return
      endif
#endif

      rc(0,1) = ri(1)
      rc(0,2) = ri(2)
      rc(0,3) = ri(3)
      rc(0,4) = ri(4)
      rc(1,1) = ri(5)
      rc(1,2) = ri(6)
      rc(1,3) = ri(7)
      rc(1,4) = ri(8)
      rc(2,1) = ri(9)
      rc(2,2) = ri(10)
      rc(2,3) = ri(11)
      rc(2,4) = ri(12)
      rc(3,1) = ri(13)
      rc(3,2) = ri(14)
      rc(3,3) = ri(15)
      rc(3,4) = ri(16)

      
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

      jiro(1) = cZero
      jiro(2) = cZero
      jiro(3) = cZero
      jiro(4) = cZero

      if ( vmass.ne.rZero ) then

       d = -rOne/dcmplx( q2-vm2, vmass*vwidth )

       if ( gc(2).ne.cZero ) then
        do i=1,4
         jiro(i) = dconjg(gc(2))*d*(fo(1)*gsqgr(i-1,1)+fo(2)
     &             *gsqgr(i-1,2))
        enddo
       endif
      
       if ( gc(1).ne.cZero ) then
        do i=1,4
         jiro(i) = jiro(i)+dconjg(gc(1))*d*(fo(3)*gsqgr(i-1,3)+fo(4)
     &             *gsqgr(i-1,4))
        enddo
       endif

      else

       d = -dcmplx( rOne/q2, rZero )

      if ( gc(2).ne.cZero ) then
        do i=1,4
         jiro(i) = dconjg(gc(2))*d*(fo(1)*gsqgr(i-1,1)+fo(2)
     &             *gsqgr(i-1,2))
        enddo
       endif
    
       if ( gc(1).ne.cZero ) then
        do i=1,4
         jiro(i) = jiro(i)+dconjg(gc(1))*d*(fo(3)*gsqgr(i-1,3)+fo(4)
     &             *gsqgr(i-1,4))
        enddo
       endif
         

      end if
c
      return
      end
