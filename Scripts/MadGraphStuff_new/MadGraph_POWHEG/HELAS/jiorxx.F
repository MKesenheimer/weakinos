      subroutine jiorxx(fi,ro,gc,vmass,vwidth , jior)
c
c This subroutine computes an off-shell vector current from a Rarita-
c Schwinger fermion and a fermion. The vector boson propagator is given
c in Feynman gauge for a massless vector and in unitary gauge for a 
c massive vector.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex ro(18)         : flow-out RS fermion               v<ro|
c       complex gc(2)          : coupling constants                 GFRV
c       real    vmass          : mass  of OUTPUT vector v'(q)
c       real    vwidth         : width of OUTPUT vector v'(q)
c
c output:
c       complex jior(6)        : vector current    j^mu(<ro|[q,v']v|fi>)
c   
c- by Y.Takaesu - 2009/11/27
c  
      implicit none
      double complex fi(6),ro(18),gc(2),jior(6),c0,c1,c2,c3,cs,d
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex cm2 ! mass**2- I Gamma mass (Fabio)

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
      if ( abs(fi(1))+abs(fi(2))+abs(fi(3))+abs(fi(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fi in jiorxx is zero spinor'
      endif
      if ( abs(fi(5))+abs(fi(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fi in jiorxx has zero momentum'
      endif
      if (abs(ro(1))+abs(ro(2))+abs(ro(3))+abs(ro(4))+abs(ro(5))
     & +abs(ro(6))+abs(ro(7))+abs(ro(8))abs(ro(9))+abs(ro(10))
     & +abs(ro(11))+abs(ro(12))+abs(ro(13))+abs(ro(14))+abs(ro(15))
     & +abs(ro(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ro in jiorxx is zero spinor'
      endif
      if ( abs(ro(17))+abs(ro(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in jiorxx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in jiorxx is zero coupling'
      endif
      if ( vmass.lt.rZero ) then
         write(stdo,*) ' helas-error : vmass in jiorxx is negative'
         write(stdo,*) '             : vmass = ',vmass
      endif
      if ( vwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : vwidth in jiorxx is negative'
         write(stdo,*) '             : vwidth = ',vwidth
      endif
#endif

      jior(5) = ro(17)-fi(5)
      jior(6) = ro(18)-fi(6)

      q(0) = dble( jior(5))
      q(1) = dble( jior(6))
      q(2) = dimag(jior(6))
      q(3) = dimag(jior(5))
            
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

#ifdef HELAS_CHECK
      if ( abs(jior(5))+abs(jior(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : jior in jiorxx has zero momentum'
      endif
      if ( vwidth.eq.rZero .and. q2.eq.vm2 ) then
         write(stdo,*)
     &        ' helas-error : jior in jiorxx is on vmass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         jior(1) = cZero
         jior(2) = cZero
         jior(3) = cZero
         jior(4) = cZero
         return
      endif
#endif

      rc(0,1) = ro(1)
      rc(0,2) = ro(2)
      rc(0,3) = ro(3)
      rc(0,4) = ro(4)
      rc(1,1) = ro(5)
      rc(1,2) = ro(6)
      rc(1,3) = ro(7)
      rc(1,4) = ro(8)
      rc(2,1) = ro(9)
      rc(2,2) = ro(10)
      rc(2,3) = ro(11)
      rc(2,4) = ro(12)
      rc(3,1) = ro(13)
      rc(3,2) = ro(14)
      rc(3,3) = ro(15)
      rc(3,4) = ro(16)
      
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

      jior(1) = cZero
      jior(2) = cZero
      jior(3) = cZero
      jior(4) = cZero

      if ( vmass.ne.rZero ) then

       d = -1d0/dcmplx( q2-vm2, vmass*vwidth )

       if ( gc(2).ne.cZero ) then
        do i=1,4
         jior(i) = jior(i)+gc(2)*d*(fi(3)*rspgg(i-1,3)+fi(4)
     &             *rspgg(i-1,4))
        enddo
       endif
      
       if ( gc(1).ne.cZero ) then
        do i=1,4
         jior(i) = jior(i)+gc(1)*d*(fi(1)*rspgg(i-1,1)+fi(2)
     &             *rspgg(i-1,2))
        enddo
       endif

      else

       d = dcmplx( -1d0/q2, rZero )

      if ( gc(2).ne.cZero ) then
        do i=1,4
         jior(i) = jior(i)+gc(2)*d*(fi(3)*rspgg(i-1,3)+fi(4)
     &             *rspgg(i-1,4))
        enddo
       endif
    
       if ( gc(1).ne.cZero ) then
        do i=1,4
         jior(i) = jior(i)+gc(1)*d*(fi(1)*rspgg(i-1,1)+fi(2)
     &             *rspgg(i-1,2))
        enddo
       endif
         
      endif

      return
      end
