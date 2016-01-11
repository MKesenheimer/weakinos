      subroutine hirocx(ric,foc,gc,smass,swidth , hiroc)
c
c This subroutine computes an off-shell scalar current from a flowing-IN
c antiRarita-Schwinger fermion and a flowing-OUT antifermion.
c
c input:
c       complex ric(18)        : flow-in  antiRSfermion            |ric>v
c       complex foc(6)         : flow-out antifermion              <foc|
c       complex gc(2)          : coupling constants                 GFRS 
c       real    smass          : mass  of OUTPUT scalar s(q)
c       real    swidth         : width of OUTPUT scalar s(q)
c
c output:
c       complex hiroc(3)         : scalar current          j(<foc|svq|ric>)
c  
c- by Yoshitaro Takaesu - 2011/03/26
c 

      implicit none
      double complex ric(18),foc(6),hiroc(3),gc(2),dn,rc(0:3,4),sq1(2,2)
      double complex sq2(2,2),cImag,sqr(0:3,4),gsqr(4),gcc(2)
      double precision q(0:3),smass,swidth,q2
      integer i

      parameter (cImag = (0d0,1d0))

#ifdef HELAS_CHECK
      double precision rZero, cZero
      parameter( rZero = 0.0d0 )
      double complex cZero
      parameter( cZero = ( 0.0d0, 0.0d0 ) )
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      if ( abs(foc(1))+abs(foc(2))+abs(foc(3))+abs(foc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : foc in hirocx is zero spinor'
      endif
      if ( abs(foc(5))+abs(foc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : foc in hirocx has zero momentum'
      endif
       if (abs(ric(1))+abs(ric(2))+abs(ric(3))+abs(ric(4))+abs(ric(5))
     & +abs(ric(6))+abs(ric(7))+abs(ric(8))abs(ric(9))+abs(ric(10))
     & +abs(ric(11))+abs(ric(12))+abs(ric(13))+abs(ric(14))+abs(ric(15))
     & +abs(ric(16)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-warn  : ric in hirocx is zero spinor'
      endif
      if ( abs(ric(17))+abs(ric(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ric in hirocx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in hirocx is zero coupling'
      endif
      if ( smass.lt.rZero ) then
         write(stdo,*) ' helas-error : smass in hirocx is negative'
         write(stdo,*) '             : smass = ',smass
      endif
      if ( swidth.lt.rZero ) then
         write(stdo,*) ' helas-error : swidth in hirocx is negative'
         write(stdo,*) '             : swidth = ',swidth
      endif
#endif

      hiroc(2) = foc(5)-ric(17)
      hiroc(3) = foc(6)-ric(18)

      q(0) = dble( hiroc(2))
      q(1) = dble( hiroc(3))
      q(2) = dimag(hiroc(3))
      q(3) = dimag(hiroc(2))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)

#ifdef HELAS_CHECK
      if ( abs(hiroc(2))+abs(hiroc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : hiroc in hirocx has zero momentum'
      endif
      if ( swidth.eq.rZero .and. q2.eq.smass**2 ) then
         write(stdo,*)
     &        ' helas-error : hiroc in hirocx is on smass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         hiroc(1) = cZero
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

      do i=0,3
        sqr(i,1) = rc(i,3)*sq1(1,1)+rc(i,4)*sq1(1,2)
        sqr(i,2) = rc(i,3)*sq1(2,1)+rc(i,4)*sq1(2,2)
	sqr(i,3) = rc(i,1)*sq2(1,1)+rc(i,2)*sq2(1,2)
	sqr(i,4) = rc(i,1)*sq2(2,1)+rc(i,2)*sq2(2,2)
      enddo

      gsqr(1)=sqr(0,3)-sqr(1,4)+cImag*sqr(2,4)-sqr(3,3)
      gsqr(2)=sqr(0,4)-sqr(1,3)-cImag*sqr(2,3)+sqr(3,4)
      gsqr(3)=sqr(0,1)+sqr(1,2)-cImag*sqr(2,2)+sqr(3,1)
      gsqr(4)=sqr(0,2)+sqr(1,1)+cImag*sqr(2,1)-sqr(3,2)

      dn = -1d0/dcmplx( q2-smass**2, smass*swidth )
      dn = -dn

      hiroc(1) = dn*( dconjg(gcc(1))*(gsqr(3)*foc(3)+gsqr(4)*foc(4))
     &              +dconjg(gcc(2))*(gsqr(1)*foc(1)+gsqr(2)*foc(2)))
     &          
      
      return
      end
