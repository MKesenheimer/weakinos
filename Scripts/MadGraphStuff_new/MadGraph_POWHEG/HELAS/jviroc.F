      subroutine jviroc(ric,foc,vc,gc,vmass,vwidth , jvirc)
c
c This subroutine computes an off-shell vector current from a flow-IN anti-
c Rarita-Schwinger fermion, a flow-OUT antifermion, and a vector boson. 
c The vector boson propagator is given in Feynman gauge c for a massless vector
c  and in unitary gauge for a massive vector.
c
c input:
c       complex ric(6)         : flow-in  antiRSfermion            |ric>
c       complex foc(18)        : flow-out antifermion              <foc|
c       complex vc(6)          : vector boson                         v
c       complex gc(2)          : coupling constants                gfrvv
c       real    vmass          : mass  of OUTPUT vector v
c       real    vwidth         : width of OUTPUT vector v
c
c output:
c       complex jvirc(6)       : vector current          j^mu(<foc|vv'|ric>)
c   
c- by Yoshitaro Takaesu - 2011/03/27
c  
      implicit none
      double complex ric(18),foc(6),gc(2),jvirc(6),c0,c1,c2,c3,cs,d
      double complex vc(6),gcc(2)
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex cm2 ! mass**2- I Gamma mass (Fabio)


      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

      double complex sq1(2,2),sq2(2,2),sqg1(0:3,2,2),sqg2(0:3,2,2)
      double complex gsq1(0:3,2,2),gsq2(0:3,2,2)
      double complex rc(0:3,4),sqgr(0:3,0:3,4),gsqgr(0:3,4)
      double complex spvl(2,2),svcl(2,2),sssl(2,2)
      double complex spvr(2,2),svcr(2,2),sssr(2,2),vertex1,vertex2
      integer i,j

#ifdef HELAS_CHECK
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      if ( abs(foc(1))+abs(foc(2))+abs(foc(3))+abs(foc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : foc in jviroc is zero spinor'
      endif
      if ( abs(foc(5))+abs(foc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : foc in jviroc has zero momentum'
      endif
      if (abs(ric(1))+abs(ric(2))+abs(ric(3))+abs(ric(4))+abs(ric(5))
     & +abs(ric(6))+abs(ric(7))+abs(ric(8))abs(ric(9))+abs(ric(10))
     & +abs(ric(11))+abs(ric(12))+abs(ric(13))+abs(ric(14))+abs(ric(15))
     & +abs(ric(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ric in jviroc is zero spinor'
      endif
      if ( abs(ric(17))+abs(ric(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ric in jviroc has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in jviroc is zero coupling'
      endif
      if ( vmass.lt.rZero ) then
         write(stdo,*) ' helas-error : vmass in jviroc is negative'
         write(stdo,*) '             : vmass = ',vmass
      endif
      if ( vwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : vwidth in jviroc is negative'
         write(stdo,*) '             : vwidth = ',vwidth
      endif
#endif

      jvirc(5) = -ric(17)+foc(5)+vc(5)
      jvirc(6) = -ric(18)+foc(6)+vc(6)

      q(0) = dble(jvirc(5))
      q(1) = dble(jvirc(6))
      q(2) = dimag(jvirc(6))
      q(3) = dimag(jvirc(5))
      
      
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

#ifdef HELAS_CHECK
      if ( abs(jvirc(5))+abs(jvirc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : jvirc in jviroc has zero momentum'
      endif
      if ( vwidth.eq.rZero .and. q2.eq.vm2 ) then
         write(stdo,*)
     &        ' helas-error : jvirc in jviroc is on vmass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         jvirc(1) = cZero
         jvirc(2) = cZero
         jvirc(3) = cZero
         jvirc(4) = cZero
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

      
      sq1(1,1) = vc(1)-vc(4)
      sq1(1,2) = -(vc(2)-cImag*vc(3))
      sq1(2,1) = -(vc(2)+cImag*vc(3))
      sq1(2,2) = vc(1)+vc(4)

      sq2(1,1) = vc(1)+vc(4)
      sq2(1,2) = vc(2)-cImag*vc(3)
      sq2(2,1) = vc(2)+cImag*vc(3)
      sq2(2,2) = vc(1)-vc(4)

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

       vertex1 = cZero
       vertex2 = cZero

      if ( gcc(2).ne.cZero ) then

        spvl(1,1) = q(0)+q(3)
        spvl(1,2) = q(1)-cImag*q(2)
        spvl(2,1) = q(1)+cImag*q(2)
        spvl(2,2) = q(0)-q(3)
        
        svcl(1,1) = vc(1)-vc(4)
        svcl(1,2) = -(vc(2)-cImag*vc(3))
        svcl(2,1) = -(vc(2)+cImag*vc(3))
        svcl(2,2) = vc(1)+vc(4)

        sssl(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
        sssl(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
        sssl(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
        sssl(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)     

      vertex2 = dconjg(gcc(2)) 
     &          *( foc(1) *( sssl(1,1)*rc(0,3)+sssl(1,2)*rc(0,4) )
     &            +foc(2) *( sssl(2,1)*rc(0,3)+sssl(2,2)*rc(0,4) )
     &            -foc(1) *( sssl(2,1)*rc(1,3)+sssl(2,2)*rc(1,4) )
     &            -foc(2) *( sssl(1,1)*rc(1,3)+sssl(1,2)*rc(1,4) )
     &            -foc(1) *(-sssl(2,1)*rc(2,3)-sssl(2,2)*rc(2,4) )*cImag
     &            -foc(2) *( sssl(1,1)*rc(2,3)+sssl(1,2)*rc(2,4) )*cImag
     &            -foc(1) *( sssl(1,1)*rc(3,3)+sssl(1,2)*rc(3,4) )
     &            -foc(2) *(-sssl(2,1)*rc(3,3)-sssl(2,2)*rc(3,4) )
     &           )


        spvl(1,1) = vc(1)+vc(4)
        spvl(1,2) = vc(2)-cImag*vc(3)
        spvl(2,1) = vc(2)+cImag*vc(3)
        spvl(2,2) = vc(1)-vc(4)
        
        svcl(1,1) = q(0)-q(3)
        svcl(1,2) = -(q(1)-cImag*q(2))
        svcl(2,1) = -(q(1)+cImag*q(2))
        svcl(2,2) = q(0)+q(3)

        sssl(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
        sssl(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
        sssl(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
        sssl(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

      vertex2 = vertex2
     &        -dconjg(gcc(2)) 
     &          *( foc(1) *( sssl(1,1)*rc(0,3)+sssl(1,2)*rc(0,4) )
     &            +foc(2) *( sssl(2,1)*rc(0,3)+sssl(2,2)*rc(0,4) )
     &            -foc(1) *( sssl(2,1)*rc(1,3)+sssl(2,2)*rc(1,4) )
     &            -foc(2) *( sssl(1,1)*rc(1,3)+sssl(1,2)*rc(1,4) )
     &            -foc(1) *(-sssl(2,1)*rc(2,3)-sssl(2,2)*rc(2,4) )*cImag
     &            -foc(2) *( sssl(1,1)*rc(2,3)+sssl(1,2)*rc(2,4) )*cImag
     &            -foc(1) *( sssl(1,1)*rc(3,3)+sssl(1,2)*rc(3,4) )
     &            -foc(2) *(-sssl(2,1)*rc(3,3)-sssl(2,2)*rc(3,4) )
     &           )
      
      endif

      if ( gcc(1).ne.cZero ) then

        spvr(1,1) = q(0)-q(3)
        spvr(1,2) = -(q(1)-cImag*q(2))
        spvr(2,1) = -(q(1)+cImag*q(2))
        spvr(2,2) = q(0)+q(3)
        
        svcr(1,1) = vc(1)+vc(4)
        svcr(1,2) = vc(2)-cImag*vc(3)
        svcr(2,1) = vc(2)+cImag*vc(3)
        svcr(2,2) = vc(1)-vc(4)

        sssr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
        sssr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
        sssr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
        sssr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

      vertex1 =
     &        +dconjg(gcc(1)) 
     &          *( foc(3) *( sssr(1,1)*rc(0,1)+sssr(1,2)*rc(0,2) )
     &            +foc(4) *( sssr(2,1)*rc(0,1)+sssr(2,2)*rc(0,2) )
     &            +foc(3) *( sssr(2,1)*rc(1,1)+sssr(2,2)*rc(1,2) )
     &            +foc(4) *( sssr(1,1)*rc(1,1)+sssr(1,2)*rc(1,2) )
     &            +foc(3) *(-sssr(2,1)*rc(2,1)-sssr(2,2)*rc(2,2) )*cImag
     &            +foc(4) *( sssr(1,1)*rc(2,1)+sssr(1,2)*rc(2,2) )*cImag
     &            +foc(3) *( sssr(1,1)*rc(3,1)+sssr(1,2)*rc(3,2) )
     &            +foc(4) *(-sssr(2,1)*rc(3,1)-sssr(2,2)*rc(3,2) )
     &           )

        spvr(1,1) = vc(1)-vc(4)
        spvr(1,2) = -(vc(2)-cImag*vc(3))
        spvr(2,1) = -(vc(2)+cImag*vc(3))
        spvr(2,2) = vc(1)+vc(4)
        
        svcr(1,1) = q(0)+q(3)
        svcr(1,2) = q(1)-cImag*q(2)
        svcr(2,1) = q(1)+cImag*q(2)
        svcr(2,2) = q(0)-q(3)

        sssr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
        sssr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
        sssr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
        sssr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

      vertex1 = vertex1 
     &        -dconjg(gcc(1)) 
     &          *( foc(3) *( sssr(1,1)*rc(0,1)+sssr(1,2)*rc(0,2) )
     &            +foc(4) *( sssr(2,1)*rc(0,1)+sssr(2,2)*rc(0,2) )
     &            +foc(3) *( sssr(2,1)*rc(1,1)+sssr(2,2)*rc(1,2) )
     &            +foc(4) *( sssr(1,1)*rc(1,1)+sssr(1,2)*rc(1,2) )
     &            +foc(3) *(-sssr(2,1)*rc(2,1)-sssr(2,2)*rc(2,2) )*cImag
     &            +foc(4) *( sssr(1,1)*rc(2,1)+sssr(1,2)*rc(2,2) )*cImag
     &            +foc(3) *( sssr(1,1)*rc(3,1)+sssr(1,2)*rc(3,2) )
     &            +foc(4) *(-sssr(2,1)*rc(3,1)-sssr(2,2)*rc(3,2) )
     &           )

      endif


      jvirc(1) = cZero
      jvirc(2) = cZero
      jvirc(3) = cZero
      jvirc(4) = cZero

      if ( vmass.ne.rZero ) then

       d = rOne/dcmplx( q2-vm2, vmass*vwidth )
       d = -d

       if ( gcc(2).ne.cZero ) then
        do i=1,4
         jvirc(i) = jvirc(i)+dconjg(gcc(2))*d*(foc(1)*gsqgr(i-1,1)
     &             +foc(2)*gsqgr(i-1,2))-d*q(i-1)/vm2*vertex2
        enddo
       endif
      
       if ( gcc(1).ne.cZero ) then
        do i=1,4
         jvirc(i) = jvirc(i)+dconjg(gcc(1))*d*(foc(3)*gsqgr(i-1,3)
     &            +foc(4)*gsqgr(i-1,4))-d*q(i-1)/vm2*vertex1
        enddo
       endif

      else

       d = dcmplx( rOne/q2, rZero )
       d = -d

      if ( gcc(2).ne.cZero ) then
        do i=1,4
         jvirc(i) = jvirc(i)+dconjg(gcc(2))*d*(foc(1)*gsqgr(i-1,1)
     &             +foc(2)*gsqgr(i-1,2))
        enddo
       endif
    
       if ( gcc(1).ne.cZero ) then
        do i=1,4
         jvirc(i) = jvirc(i)+dconjg(gcc(1))*d*(foc(3)*gsqgr(i-1,3)
     &             +foc(4)*gsqgr(i-1,4))
        enddo
       endif
         

      end if
c
      return
      end
