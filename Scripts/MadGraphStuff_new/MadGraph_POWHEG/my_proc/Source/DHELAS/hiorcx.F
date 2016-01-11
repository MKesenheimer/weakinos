      subroutine hiorcx(fic,roc,gc,smass,swidth , hiorc)
c
c This subroutine computes an off-shell scalar current from a flowing-IN
c antifermion and a flowing-OUT antiRarita-Schwinger fermion.
c
c input:
c       complex fic(6)         : flow-in  antifermion              |fic>
c       complex roc(18)        : flow-out antiRSfermion           v<roc|
c       complex gc(2)          : coupling constants                 GFRS 
c       real    smass          : mass  of OUTPUT scalar s(q)
c       real    swidth         : width of OUTPUT scalar s(q)
c
c output:
c       complex hio(3)         : scalar current           j(<roc|sqv|fic>)
c
c- by Yoshitaro Takaesu - 2011/03/26
c
      implicit none
      double complex fic(6),roc(18),hiorc(3),gc(2),dn,rc(0:3,4),sq1(2,2) 
      double complex sq2(2,2),cImag,rsq(0:3,4),rsqg(4),gcc(2)
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
      if ( abs(fic(1))+abs(fic(2))+abs(fic(3))+abs(fic(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fic in hiorcx is zero spinor'
      endif
      if ( abs(fic(5))+abs(fic(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fic in hiorcx has zero momentum'
      endif
       if (abs(roc(1))+abs(roc(2))+abs(roc(3))+abs(roc(4))+abs(roc(5))
     & +abs(roc(6))+abs(roc(7))+abs(roc(8))abs(roc(9))+abs(roc(10))
     & +abs(roc(11))+abs(roc(12))+abs(roc(13))+abs(roc(14))+abs(roc(15))
     & +abs(roc(16)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-warn  : roc in hiorcx is zero spinor'
      endif
      if ( abs(roc(17))+abs(roc(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : roc in hiorcx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in hiorcx is zero coupling'
      endif
      if ( smass.lt.rZero ) then
         write(stdo,*) ' helas-error : smass in hiorcx is negative'
         write(stdo,*) '             : smass = ',smass
      endif
      if ( swidth.lt.rZero ) then
         write(stdo,*) ' helas-error : swidth in hiorcx is negative'
         write(stdo,*) '             : swidth = ',swidth
      endif
#endif

      hiorc(2) = roc(17)-fic(5)
      hiorc(3) = roc(18)-fic(6)

      q(0) = dble( hiorc(2))
      q(1) = dble( hiorc(3))
      q(2) = dimag(hiorc(3))
      q(3) = dimag(hiorc(2))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)

#ifdef HELAS_CHECK
      if ( abs(hiorc(2))+abs(hiorc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : hiorc in hiorcx has zero momentum'
      endif
      if ( swidth.eq.rZero .and. q2.eq.smass**2 ) then
         write(stdo,*)
     &        ' helas-error : hiorc in hiorcx is on smass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         hiorc(1) = cZero
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

      dn = 1d0/dcmplx( q2-smass**2, smass*swidth )
      dn = -dn

      hiorc(1) = dn*( gcc(1)*(rsqg(1)*fic(1)+rsqg(2)*fic(2))
     &              +gcc(2)*(rsqg(3)*fic(3)+rsqg(4)*fic(4)))
      
      return
      end
