      subroutine jsiorc(fic,roc,sc,gc,vmass,vwidth , jsioc)
c
c This subroutine computes an off-shell vector current from an anti-
c Rarita-Schwinger fermion, antifermion, and a scalar. The vector
c boson propagator is given in Feynman gauge for a massless vector and
c in unitary gauge for a massive vector.
c
c input:
c       complex fic(6)         : flow-in  antifermion               |fic>
c       complex roc(18)        : flow-out antiRSfermion             <roc|
c       complex gc(2)          : coupling constants                 gfrvs
c       real    vmass          : mass  of OUTPUT vector v
c       real    vwidth         : width of OUTPUT vector v
c
c output:
c       complex jsioc(6)        : vector current          j^mu(<roc|s|fic>)
c   
c- by Yoshitaro Takaesu - 2011/03/27
c  
      implicit none
      double complex fic(6),roc(18),gc(2),sc(3),jsioc(6)  
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex gcc(2),c0,c1,c2,c3,cs,d
      double complex cm2 ! mass**2- I Gamma mass (Fabio)

      double complex rsq(0:3,4),rsqg(4),sq1(2,2),sq2(2,2),rsqgf,rf(4)
      double complex rc(0:3,4)
      integer i,j

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

#ifdef HELAS_CHECK
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      if ( abs(fic(1))+abs(fic(2))+abs(fic(3))+abs(fic(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fic in jsiorc is zero spinor'
      endif
      if ( abs(fic(5))+abs(fic(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fic in jsiorc has zero momentum'
      endif
      if (abs(roc(1))+abs(roc(2))+abs(roc(3))+abs(roc(4))+abs(roc(5))
     & +abs(roc(6))+abs(roc(7))+abs(roc(8))abs(roc(9))+abs(roc(10))
     & +abs(roc(11))+abs(roc(12))+abs(roc(13))+abs(roc(14))+abs(roc(15))
     & +abs(roc(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : roc in jsiorc is zero spinor'
      endif
      if ( abs(roc(17))+abs(roc(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in jsiorc has zero momentum'
      endif
	if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in jsiorc is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in jsiorc has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in jsiorc is zero coupling'
      endif
      if ( vmass.lt.rZero ) then
         write(stdo,*) ' helas-error : vmass in jsiorc is negative'
         write(stdo,*) '             : vmass = ',vmass
      endif
      if ( vwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : vwidth in jsiorcx is negative'
         write(stdo,*) '             : vwidth = ',vwidth
      endif
#endif

      jsioc(5) = roc(17)-fic(5)+sc(2)
      jsioc(6) = roc(18)-fic(6)+sc(3)

      q(0) = dble( jsioc(5))
      q(1) = dble( jsioc(6))
      q(2) = dimag(jsioc(6))
      q(3) = dimag(jsioc(5))
      
      
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

#ifdef HELAS_CHECK
      if ( abs(jsioc(5))+abs(jsioc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : jsioc in jsiorc has zero momentum'
      endif
      if ( vwidth.eq.rZero .and. q2.eq.vm2 ) then
         write(stdo,*)
     &        ' helas-error : jsioc in jsiorc is on vmass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         jsioc(1) = cZero
         jsioc(2) = cZero
         jsioc(3) = cZero
         jsioc(4) = cZero
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

      rf(1) = gcc(1)*(rc(0,1)*fic(1)+rc(0,2)*fic(2))
     &       +gcc(2)*(rc(0,3)*fic(3)+rc(0,4)*fic(4))
      rf(2) = gcc(1)*(rc(1,1)*fic(1)+rc(1,2)*fic(2))
     &       +gcc(2)*(rc(1,3)*fic(3)+rc(1,4)*fic(4))
      rf(3) = gcc(1)*(rc(2,1)*fic(1)+rc(2,2)*fic(2))
     &       +gcc(2)*(rc(2,3)*fic(3)+rc(2,4)*fic(4))
      rf(4) = gcc(1)*(rc(3,1)*fic(1)+rc(3,2)*fic(2))
     &       +gcc(2)*(rc(3,3)*fic(3)+rc(3,4)*fic(4))

      sq1(1,1) = q(0)-q(3)
      sq1(1,2) = -(q(1)-cImag*q(2))
      sq1(2,1) = -(q(1)+cImag*q(2))
      sq1(2,2) = q(0)+q(3)

      sq2(1,1) = q(0)+q(3)
      sq2(1,2) = q(1)-cImag*q(2)
      sq2(2,1) = q(1)+cImag*q(2)
      sq2(2,2) = q(0)-q(3)

      do i=0,3
        rsq(i,1) = rc(i,3)*sq2(1,1)+rc(i,4)*sq2(2,1)
        rsq(i,2) = rc(i,3)*sq2(1,2)+rc(i,4)*sq2(2,2)
        rsq(i,3) = rc(i,1)*sq1(1,1)+rc(i,2)*sq1(2,1)
        rsq(i,4) = rc(i,1)*sq1(1,2)+rc(i,2)*sq1(2,2)
      enddo

      rsqg(1)=rsq(0,3)+rsq(1,4)+cImag*rsq(2,4)+rsq(3,3)
      rsqg(2)=rsq(0,4)+rsq(1,3)-cImag*rsq(2,3)-rsq(3,4)
      rsqg(3)=rsq(0,1)-rsq(1,2)-cImag*rsq(2,2)-rsq(3,1)
      rsqg(4)=rsq(0,2)-rsq(1,1)+cImag*rsq(2,1)+rsq(3,2)

      rsqgf = sc(1)*(gcc(1)*(rsqg(1)*fic(1)+rsqg(2)*fic(2))
     &              +gcc(2)*(rsqg(3)*fic(3)+rsqg(4)*fic(4)))
      
      jsioc(1) = cZero
      jsioc(2) = cZero
      jsioc(3) = cZero
      jsioc(4) = cZero

      if ( vmass.ne.rZero ) then

       d = rOne/dcmplx( q2-vm2, vmass*vwidth )

       do i=1,4
	  jsioc(i) = d*(2*rf(i)*sc(1)-q(i-1)/vmass**2*rsqgf)
	 enddo

      else

       d = dcmplx( rOne/q2, rZero )
	 do i=1,4
	  jsioc(i) = d*2*rf(i)*sc(1)
	 enddo         

      end if

      return
      end
 
