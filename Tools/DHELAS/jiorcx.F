      subroutine jiorcx(fic,roc,gc,vmass,vwidth , jiorc)
c
c This subroutine computes an off-shell vector current from an antiRarita-
c Schwinger fermion and an antifermion. The vector boson propagator is given
c in Feynman gauge for a massless vector and in unitary gauge for a 
c massive vector.
c
c input:
c       complex fic(6)         : flow-in  antifermion              |fic>
c       complex roc(18)        : flow-out antiRSfermion           v<roc|
c       complex gc(2)          : coupling constants                 GFRV
c       real    vmass          : mass  of OUTPUT vector v'(q)
c       real    vwidth         : width of OUTPUT vector v'(q)
c
c output:
c       complex jiorc(6)        : vector current    j^mu(<roc|[q,v']v|fic>)
c   
c- by Yoshitaro Takaesu - 2011/03/27
c  
      implicit none
      double complex fic(6),roc(18),gc(2),jiorc(6),c0,c1,c2,c3,cs,d
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex gcc(2)
      double complex cm2        ! mass**2- I Gamma mass (Fabio)

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

      double complex sql(2,2),sqr(2,2),sqgl(0:3,2,2),sqgr(0:3,2,2)
      double complex gsql(0:3,2,2),gsqr(0:3,2,2)
      double complex rc(0:3,4),rspg(0:3,0:3,4),rspgg(0:3,4)
      integer i,j

#ifdef HELAS_CHECK
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      if ( abs(fic(1))+abs(fic(2))+abs(fic(3))+abs(fic(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fic in jiorcx is zero spinor'
      endif
      if ( abs(fic(5))+abs(fic(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fic in jiorcx has zero momentum'
      endif
      if (abs(roc(1))+abs(roc(2))+abs(roc(3))+abs(roc(4))+abs(roc(5))
     & +abs(roc(6))+abs(roc(7))+abs(roc(8))abs(roc(9))+abs(roc(10))
     & +abs(roc(11))+abs(roc(12))+abs(roc(13))+abs(roc(14))+abs(roc(15))
     & +abs(roc(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : roc in jiorcx is zero spinor'
      endif
      if ( abs(roc(17))+abs(roc(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in jiorcx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in jiorcx is zero coupling'
      endif
      if ( vmass.lt.rZero ) then
         write(stdo,*) ' helas-error : vmass in jiorcx is negative'
         write(stdo,*) '             : vmass = ',vmass
      endif
      if ( vwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : vwidth in jiorcx is negative'
         write(stdo,*) '             : vwidth = ',vwidth
      endif
#endif

      jiorc(5) = roc(17)-fic(5)
      jiorc(6) = roc(18)-fic(6)

      q(0) = dble( jiorc(5))
      q(1) = dble( jiorc(6))
      q(2) = dimag(jiorc(6))
      q(3) = dimag(jiorc(5))
            
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

#ifdef HELAS_CHECK
      if ( abs(jiorc(5))+abs(jiorc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : jiorc in jiorcx has zero momentum'
      endif
      if ( vwidth.eq.rZero .and. q2.eq.vm2 ) then
         write(stdo,*)
     &        ' helas-error : jiorc in jiorcx is on vmass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         jiorc(1) = cZero
         jiorc(2) = cZero
         jiorc(3) = cZero
         jiorc(4) = cZero
         return
      endif
#endif

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1))

      rc(0,1) = roc(1)
      rc(0,2) = roc(2)
      rc(0,3) = roc(3)
      rc(0,4) = roc(4)
      rc(1,1) = roc(5)
      rc(1,2) = roc(6)
      rc(1,3) = roc(7)
      rc(1,4) = roc(8)
      rc(2,1) = roc(9)
      rc(2,2) = roc(10)
      rc(2,3) = roc(11)
      rc(2,4) = roc(12)
      rc(3,1) = roc(13)
      rc(3,2) = roc(14)
      rc(3,3) = roc(15)
      rc(3,4) = roc(16)
      
      sql(1,1) = q(0)-q(3)
      sql(1,2) = -(q(1)-cImag*q(2))
      sql(2,1) = -(q(1)+cImag*q(2))
      sql(2,2) = q(0)+q(3)

      sqr(1,1) = q(0)+q(3)
      sqr(1,2) = q(1)-cImag*q(2)
      sqr(2,1) = q(1)+cImag*q(2)
      sqr(2,2) = q(0)-q(3)

      sqgl(0,1,1) = sql(1,1)
      sqgl(0,1,2) = sql(1,2)
      sqgl(0,2,1) = sql(2,1)
      sqgl(0,2,2) = sql(2,2)
      sqgl(1,1,1) = -sql(1,2)
      sqgl(1,1,2) = -sql(1,1)
      sqgl(1,2,1) = -sql(2,2)
      sqgl(1,2,2) = -sql(2,1)
      sqgl(2,1,1) = -cImag*sql(1,2)
      sqgl(2,1,2) = cImag*sql(1,1)
      sqgl(2,2,1) = -cImag*sql(2,2)
      sqgl(2,2,2) = cImag*sql(2,1)
      sqgl(3,1,1) = -sql(1,1)
      sqgl(3,1,2) = sql(1,2)
      sqgl(3,2,1) = -sql(2,1)
      sqgl(3,2,2) = sql(2,2)

      sqgr(0,1,1) = sqr(1,1)
      sqgr(0,1,2) = sqr(1,2)
      sqgr(0,2,1) = sqr(2,1)
      sqgr(0,2,2) = sqr(2,2)
      sqgr(1,1,1) = sqr(1,2)
      sqgr(1,1,2) = sqr(1,1)
      sqgr(1,2,1) = sqr(2,2)
      sqgr(1,2,2) = sqr(2,1)
      sqgr(2,1,1) = cImag*sqr(1,2)
      sqgr(2,1,2) = -cImag*sqr(1,1)
      sqgr(2,2,1) = cImag*sqr(2,2)
      sqgr(2,2,2) = -cImag*sqr(2,1)
      sqgr(3,1,1) = sqr(1,1)
      sqgr(3,1,2) = -sqr(1,2)
      sqgr(3,2,1) = sqr(2,1)
      sqgr(3,2,2) = -sqr(2,2)

      gsql(0,1,1) = sqr(1,1)
      gsql(0,1,2) = sqr(1,2)
      gsql(0,2,1) = sqr(2,1)
      gsql(0,2,2) = sqr(2,2)
      gsql(1,1,1) = sqr(2,1)
      gsql(1,1,2) = sqr(2,2)
      gsql(1,2,1) = sqr(1,1)
      gsql(1,2,2) = sqr(1,2)
      gsql(2,1,1) = -cImag*sqr(2,1)
      gsql(2,1,2) = -cImag*sqr(2,2)
      gsql(2,2,1) = cImag*sqr(1,1)
      gsql(2,2,2) = cImag*sqr(1,2)
      gsql(3,1,1) = sqr(1,1)
      gsql(3,1,2) = sqr(1,2)
      gsql(3,2,1) = -sqr(2,1)
      gsql(3,2,2) = -sqr(2,2)

      gsqr(0,1,1) = sql(1,1)
      gsqr(0,1,2) = sql(1,2)
      gsqr(0,2,1) = sql(2,1)
      gsqr(0,2,2) = sql(2,2)
      gsqr(1,1,1) = -sql(2,1)
      gsqr(1,1,2) = -sql(2,2)
      gsqr(1,2,1) = -sql(1,1)
      gsqr(1,2,2) = -sql(1,2)
      gsqr(2,1,1) = cImag*sql(2,1)
      gsqr(2,1,2) = cImag*sql(2,2)
      gsqr(2,2,1) = -cImag*sql(1,1)
      gsqr(2,2,2) = -cImag*sql(1,2)
      gsqr(3,1,1) = -sql(1,1)
      gsqr(3,1,2) = -sql(1,2)
      gsqr(3,2,1) = sql(2,1)
      gsqr(3,2,2) = sql(2,2)

      do i=0,3
       do j=0,3    
        rspg(i,j,1) = rc(j,1)*(sqgl(i,1,1)-gsql(i,1,1))+rc(j,2)
     &*(sqgl(i,2,1)-gsql(i,2,1))
        rspg(i,j,2) = rc(j,1)*(sqgl(i,1,2)-gsql(i,1,2))+rc(j,2)
     &*(sqgl(i,2,2)-gsql(i,2,2))   
        rspg(i,j,3) = rc(j,3)*(sqgr(i,1,1)-gsqr(i,1,1))+rc(j,4)
     &*(sqgr(i,2,1)-gsqr(i,2,1))
        rspg(i,j,4) = rc(j,3)*(sqgr(i,1,2)-gsqr(i,1,2))+rc(j,4)
     &*(sqgr(i,2,2)-gsqr(i,2,2))
       enddo
      enddo

      do i=0,3
       rspgg(i,1) = rspg(i,0,3)+rspg(i,1,4)+cImag*rspg(i,2,4)
     &              +rspg(i,3,3) 
       rspgg(i,2) = rspg(i,0,4)+rspg(i,1,3)-cImag*rspg(i,2,3)
     &              -rspg(i,3,4) 
       rspgg(i,3) = rspg(i,0,1)-rspg(i,1,2)-cImag*rspg(i,2,2)
     &              -rspg(i,3,1)
       rspgg(i,4) = rspg(i,0,2)-rspg(i,1,1)+cImag*rspg(i,2,1)
     &              +rspg(i,3,2) 
      enddo

      jiorc(1) = cZero
      jiorc(2) = cZero
      jiorc(3) = cZero
      jiorc(4) = cZero

      if ( vmass.ne.rZero ) then

       d = -1d0/dcmplx( q2-vm2, vmass*vwidth )

       if ( gcc(2).ne.cZero ) then
        do i=1,4
         jiorc(i) = jiorc(i)+gcc(2)*d*(fic(3)*rspgg(i-1,3)+fic(4)
     &             *rspgg(i-1,4))
        enddo
       endif
      
       if ( gcc(1).ne.cZero ) then
        do i=1,4
         jiorc(i) = jiorc(i)+gcc(1)*d*(fic(1)*rspgg(i-1,1)+fic(2)
     &             *rspgg(i-1,2))
        enddo
       endif

      else

       d = dcmplx( -1d0/q2, rZero )

      if ( gcc(2).ne.cZero ) then
        do i=1,4
         jiorc(i) = jiorc(i)+gcc(2)*d*(fic(3)*rspgg(i-1,3)+fic(4)
     &             *rspgg(i-1,4))
        enddo
       endif
    
       if ( gcc(1).ne.cZero ) then
        do i=1,4
         jiorc(i) = jiorc(i)+gcc(1)*d*(fic(1)*rspgg(i-1,1)+fic(2)
     &             *rspgg(i-1,2))
        enddo
       endif
         
      endif

      return
      end
