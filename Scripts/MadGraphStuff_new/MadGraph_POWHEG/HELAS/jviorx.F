      subroutine jviorx(fi,ro,vc,gc,vmass,vwidth , jvior)
c
c This subroutine computes an off-shell vector current from an external
c Rarita-Schwinger fermion, an internal fermion, and a vector boson.
c The vector boson propagator is given in Feynman gauge for a massless
c vector and in unitary gauge for a massive vector.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex ro(18)         : flow-out Rarita Schwinger fermion  <fo|
c       complex vc(6)          : vector boson                         v
c       complex gc(2)          : coupling constants                gfrvv
c       real    vmass          : mass  of OUTPUT vector v
c       real    vwidth         : width of OUTPUT vector v
c
c output:
c       complex jvior(6)       : vector current          j^mu(<ro|vv'|fi>)
c   
c- by Y.Takaesu - 2009/11/12
c  
      implicit none
      double complex fi(6),ro(18),gc(2),jvior(6),c0,c1,c2,c3,cs,d
      double complex vc(6)
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex cm2 ! mass**2- I Gamma mass (Fabio)


      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

      double complex sql(2,2),sqr(2,2),sqgl(0:3,2,2),sqgr(0:3,2,2)
      double complex gsql(0:3,2,2),gsqr(0:3,2,2)
      double complex rc(0:3,4),rspg(0:3,0:3,4),rspgg(0:3,4)
      double complex pv(4),spvl(2,2),svcl(2,2),sssl(2,2),vertex1
      double complex spvr(2,2),svcr(2,2),sssr(2,2),vertex2
      integer i,j

#ifdef HELAS_CHECK
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      if ( abs(fi(1))+abs(fi(2))+abs(fi(3))+abs(fi(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fi in jviorx is zero spinor'
      endif
      if ( abs(fi(5))+abs(fi(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fi in jviorx has zero momentum'
      endif
      if (abs(ro(1))+abs(ro(2))+abs(ro(3))+abs(ro(4))+abs(ro(5))
     & +abs(ro(6))+abs(ro(7))+abs(ro(8))abs(ro(9))+abs(ro(10))
     & +abs(ro(11))+abs(ro(12))+abs(ro(13))+abs(ro(14))+abs(ro(15))
     & +abs(ro(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ro in jviorx is zero spinor'
      endif
      if ( abs(ro(17))+abs(ro(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ro in jviorx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in jviorx is zero coupling'
      endif
      if ( vmass.lt.rZero ) then
         write(stdo,*) ' helas-error : vmass in jviorx is negative'
         write(stdo,*) '             : vmass = ',vmass
      endif
      if ( vwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : vwidth in jviorx is negative'
         write(stdo,*) '             : vwidth = ',vwidth
      endif
#endif

      jvior(5) = ro(17)-fi(5)+vc(5)
      jvior(6) = ro(18)-fi(6)+vc(6)

      q(0) = dble(jvior(5))
      q(1) = dble(jvior(6))
      q(2) = dimag(jvior(6))
      q(3) = dimag(jvior(5))
      
      
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

#ifdef HELAS_CHECK
      if ( abs(jvior(5))+abs(jvior(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : jvior in jviorx has zero momentum'
      endif
      if ( vwidth.eq.rZero .and. q2.eq.vm2 ) then
         write(stdo,*)
     &        ' helas-error : jvior in jviorx is on vmass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         jvior(1) = cZero
         jvior(2) = cZero
         jvior(3) = cZero
         jvior(4) = cZero
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

      sql(1,1) = vc(1)-vc(4)
      sql(1,2) = -(vc(2)-cImag*vc(3))
      sql(2,1) = -(vc(2)+cImag*vc(3))
      sql(2,2) = vc(1)+vc(4)

      sqr(1,1) = vc(1)+vc(4)
      sqr(1,2) = vc(2)-cImag*vc(3)
      sqr(2,1) = vc(2)+cImag*vc(3)
      sqr(2,2) = vc(1)-vc(4)

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

      vertex2 = cZero
      vertex1 = cZero

      if ( gc(2).ne.cZero ) then

        spvl(1,1) = q(0)-q(3)
        spvl(1,2) = -(q(1)-cImag*q(2))
        spvl(2,1) = -(q(1)+cImag*q(2))
        spvl(2,2) = q(0)+q(3)
        
        svcl(1,1) = vc(1)+vc(4)
        svcl(1,2) = vc(2)-cImag*vc(3)
        svcl(2,1) = vc(2)+cImag*vc(3)
        svcl(2,2) = vc(1)-vc(4)

        sssl(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
        sssl(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
        sssl(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
        sssl(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

      vertex2 = gc(2) 
     &          *( rc(0,1) *( sssl(1,1)*fi(3)+sssl(1,2)*fi(4) )
     &            +rc(0,2) *( sssl(2,1)*fi(3)+sssl(2,2)*fi(4) )           
     &            -rc(1,1) *( sssl(1,2)*fi(3)+sssl(1,1)*fi(4) )
     &            -rc(1,2) *( sssl(2,2)*fi(3)+sssl(2,1)*fi(4) )
     &            -rc(2,1) *( sssl(1,2)*fi(3)-sssl(1,1)*fi(4) )*cImag
     &            -rc(2,2) *( sssl(2,2)*fi(3)-sssl(2,1)*fi(4) )*cImag
     &            -rc(3,1) *( sssl(1,1)*fi(3)-sssl(1,2)*fi(4) )
     &            -rc(3,2) *( sssl(2,1)*fi(3)-sssl(2,2)*fi(4) ) 
     &           )

        spvl(1,1) = vc(1)-vc(4)
        spvl(1,2) = -(vc(2)-cImag*vc(3))
        spvl(2,1) = -(vc(2)+cImag*vc(3))
        spvl(2,2) = vc(1)+vc(4)
        
        svcl(1,1) = q(0)+q(3)
        svcl(1,2) = q(1)-cImag*q(2)
        svcl(2,1) = q(1)+cImag*q(2)
        svcl(2,2) = q(0)-q(3)

        sssl(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
        sssl(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
        sssl(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
        sssl(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

      vertex2 = vertex2
     &        -gc(2) 
     &          *( rc(0,1) *( sssl(1,1)*fi(3)+sssl(1,2)*fi(4) )
     &            +rc(0,2) *( sssl(2,1)*fi(3)+sssl(2,2)*fi(4) )           
     &            -rc(1,1) *( sssl(1,2)*fi(3)+sssl(1,1)*fi(4) )
     &            -rc(1,2) *( sssl(2,2)*fi(3)+sssl(2,1)*fi(4) )
     &            -rc(2,1) *( sssl(1,2)*fi(3)-sssl(1,1)*fi(4) )*cImag
     &            -rc(2,2) *( sssl(2,2)*fi(3)-sssl(2,1)*fi(4) )*cImag
     &            -rc(3,1) *( sssl(1,1)*fi(3)-sssl(1,2)*fi(4) )
     &            -rc(3,2) *( sssl(2,1)*fi(3)-sssl(2,2)*fi(4) ) 
     &           )
      
      endif

      if ( gc(1).ne.cZero ) then

        spvr(1,1) = q(0)+q(3)
        spvr(1,2) = q(1)-cImag*q(2)
        spvr(2,1) = q(1)+cImag*q(2)
        spvr(2,2) = q(0)-q(3)
        
        svcr(1,1) = vc(1)-vc(4)
        svcr(1,2) = -(vc(2)-cImag*vc(3))
        svcr(2,1) = -(vc(2)+cImag*vc(3))
        svcr(2,2) = vc(1)+vc(4)

        sssr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
        sssr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
        sssr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
        sssr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

      vertex1 = gc(1) 
     &          *( rc(0,3) *( sssr(1,1)*fi(1)+sssr(1,2)*fi(2) )
     &            +rc(0,4) *( sssr(2,1)*fi(1)+sssr(2,2)*fi(2) )           
     &            +rc(1,3) *( sssr(1,2)*fi(1)+sssr(1,1)*fi(2) )
     &            +rc(1,4) *( sssr(2,2)*fi(1)+sssr(2,1)*fi(2) )
     &            +rc(2,3) *( sssr(1,2)*fi(1)-sssr(1,1)*fi(2) )*cImag
     &            +rc(2,4) *( sssr(2,2)*fi(1)-sssr(2,1)*fi(2) )*cImag
     &            +rc(3,3) *( sssr(1,1)*fi(1)-sssr(1,2)*fi(2) )
     &            +rc(3,4) *( sssr(2,1)*fi(1)-sssr(2,2)*fi(2) ) 
     &           )

        spvr(1,1) = vc(1)+vc(4)
        spvr(1,2) = vc(2)-cImag*vc(3)
        spvr(2,1) = vc(2)+cImag*vc(3)
        spvr(2,2) = vc(1)-vc(4)
        
        svcr(1,1) = q(0)-q(3)
        svcr(1,2) = -(q(1)-cImag*q(2))
        svcr(2,1) = -(q(1)+cImag*q(2))
        svcr(2,2) = q(0)+q(3)

        sssr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
        sssr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
        sssr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
        sssr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

      vertex1 = vertex1 
     &        -gc(1) 
     &          *( rc(0,3) *( sssr(1,1)*fi(1)+sssr(1,2)*fi(2) )
     &            +rc(0,4) *( sssr(2,1)*fi(1)+sssr(2,2)*fi(2) )           
     &            +rc(1,3) *( sssr(1,2)*fi(1)+sssr(1,1)*fi(2) )
     &            +rc(1,4) *( sssr(2,2)*fi(1)+sssr(2,1)*fi(2) )
     &            +rc(2,3) *( sssr(1,2)*fi(1)-sssr(1,1)*fi(2) )*cImag
     &            +rc(2,4) *( sssr(2,2)*fi(1)-sssr(2,1)*fi(2) )*cImag
     &            +rc(3,3) *( sssr(1,1)*fi(1)-sssr(1,2)*fi(2) )
     &            +rc(3,4) *( sssr(2,1)*fi(1)-sssr(2,2)*fi(2) ) 
     &           )

      endif
   
      jvior(1) = cZero
      jvior(2) = cZero
      jvior(3) = cZero
      jvior(4) = cZero

      if ( vmass.ne.rZero ) then

       d = -rOne/dcmplx( q2-vm2, vmass*vwidth )

       if ( gc(2).ne.cZero ) then
        do i=1,4
         jvior(i) = jvior(i)+gc(2)*d*(fi(3)*rspgg(i-1,3)+fi(4)
     &             *rspgg(i-1,4))+d*q(i-1)/vmass**2*vertex2
        enddo
       endif
      
       if ( gc(1).ne.cZero ) then
        do i=1,4
         jvior(i) = jvior(i)+gc(1)*d*(fi(1)*rspgg(i-1,1)+fi(2)
     &             *rspgg(i-1,2))+d*q(i-1)/vmass**2*vertex1
        enddo
       endif

      else

       d = dcmplx( -rOne/q2, rZero )

      if ( gc(2).ne.cZero ) then
        do i=1,4
         jvior(i) = jvior(i)+gc(2)*d*(fi(3)*rspgg(i-1,3)+fi(4)
     &             *rspgg(i-1,4))
        enddo
       endif
    
       if ( gc(1).ne.cZero ) then
        do i=1,4
         jvior(i) = jvior(i)+gc(1)*d*(fi(1)*rspgg(i-1,1)+fi(2)
     &             *rspgg(i-1,2))
        enddo
       endif
         
      end if

      return
      end
