      subroutine fvorcx(roc,vc,gc,fmass,fwidth , fvorc)
c
c This subroutine computes an off-shell antifermion wavefunction from a
c flowing-OUT antiRarita-Schwinger fermion and a vector boson.
c
c
c input:
c       complex roc(18)        : flow-out antiRS fermion           v<roc|
c       complex vc(6)          : input    vector                   v'(q)
c       complex gc(2)          : coupling constants                 GFRV
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fvorc(6)        : off-shell fermion       <roc,[q,v']v,f'|
c
c- by Yoshitaro Takaesu - 2011/03/25
c
      implicit none
      double complex roc(18),vc(6),gc(2),gcc(2),fvorc(6),d
      double precision fmass,fwidth

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

      integer i
      double precision pf(0:3),pf2,sabsri,pv(0:3)
      double complex svcl(2,2),svcr(2,2),spvl(2,2),spvr(2,2),rc(0:3,4)
      double complex eprc(0:3,4),geprc(4),sssl(2,2),sssr(2,2)
      double complex spfl(2,2),spfr(2,2)

#ifdef HELAS_CHECK
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      if ( abs(roc(1))+abs(roc(2))+abs(roc(3))+abs(roc(4))+abs(roc(5))
     & +abs(roc(6))+abs(roc(7))+abs(roc(8))abs(roc(9))+abs(roc(10))
     & +abs(roc(11))+abs(roc(12))+abs(roc(13))+abs(roc(14))+abs(roc(15))
     & +abs(roc(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : roc in fvorcx is zero spinor'
      endif
      if ( abs(roc(17))+abs(roc(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : roc in fvorcx has zero momentum'
      endif
      if ( abs(vc(1))+abs(vc(2)+abs(vc(3))+abs(vc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc in fvorcx is zero vector'
      endif
      if ( abs(vc(5))+abs(vc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc in fvorcx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in fvorcx is zero coupling'
      endif
      if ( fmass.lt.rZero ) then
         write(stdo,*) ' helas-error : fmass in fvorcx is negative'
         write(stdo,*) '             : fmass = ',fmass
      endif
      if ( fwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : fwidth in fvorcx is negative'
         write(stdo,*) '             : fwidth = ',fwidth
      endif
#endif

      fvorc(5) = roc(17)+vc(5)
      fvorc(6) = roc(18)+vc(6)

      pf(0) = dble( fvorc(5))
      pf(1) = dble( fvorc(6))
      pf(2) = dimag(fvorc(6))
      pf(3) = dimag(fvorc(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

#ifdef HELAS_CHECK
      if ( abs(fvorc(5))+abs(fvorc(6)).eq.rZero ) then
          write(stdo,*)
     &        ' helas-error : fvorc in fvorcx has zero momentum'
       endif
       if ( fwidth.eq.rZero .and. pf2.eq.fmass**2 ) then
          write(stdo,*)
     &         ' helas-error : fvorc in fvorcx is on fmass pole'
          write(stdo,*)
     &         '             : p     = ',pf(0),pf(1),pf(2),pf(3)
          write(stdo,*)
     &         '             : abs(p)= ',sqrt(abs(pf2))
         fvorc(1) = cZero
         fvorc(2) = cZero
         fvorc(3) = cZero
         fvorc(4) = cZero
         return
      endif
#endif

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1)) 

      pv(0) = dble( vc(5))
      pv(1) = dble( vc(6))
      pv(2) = dimag(vc(6))
      pv(3) = dimag(vc(5))

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

      fvorc(1) = cZero
      fvorc(2) = cZero
      fvorc(3) = cZero
      fvorc(4) = cZero

      svcl(1,1) = vc(1)+vc(4)
      svcl(1,2) = vc(2)-cImag*vc(3)
      svcl(2,1) = vc(2)+cImag*vc(3)
      svcl(2,2) = vc(1)-vc(4)

      svcr(1,1) = vc(1)-vc(4)
      svcr(1,2) = -(vc(2)-cImag*vc(3))
      svcr(2,1) = -(vc(2)+cImag*vc(3))
      svcr(2,2) = vc(1)+vc(4)
        
      spvl(1,1) = pv(0)-pv(3)
      spvl(1,2) = -(pv(1)-cImag*pv(2))
      spvl(2,1) = -(pv(1)+cImag*pv(2))
      spvl(2,2) = pv(0)+pv(3)

      spvr(1,1) = pv(0)+pv(3)
      spvr(1,2) = pv(1)-cImag*pv(2)
      spvr(2,1) = pv(1)+cImag*pv(2)
      spvr(2,2) = pv(0)-pv(3)

      spfl(1,1) = pf(0)+pf(3)
      spfl(1,2) = pf(1)-cImag*pf(2)
      spfl(2,1) = pf(1)+cImag*pf(2)
      spfl(2,2) = pf(0)-pf(3)

      spfr(1,1) = pf(0)-pf(3)
      spfr(1,2) = -(pf(1)-cImag*pf(2))
      spfr(2,1) = -(pf(1)+cImag*pf(2))
      spfr(2,2) = pf(0)+pf(3)

      d = -rOne/dcmplx(pf2-fmass**2,fmass*fwidth)

      if ( gcc(1).ne.cZero ) then
        
       sssl(1,1) = svcl(1,1)*spvl(1,1)+svcl(1,2)*spvl(2,1)
       sssl(1,2) = svcl(1,1)*spvl(1,2)+svcl(1,2)*spvl(2,2) 
       sssl(2,1) = svcl(2,1)*spvl(1,1)+svcl(2,2)*spvl(2,1) 
       sssl(2,2) = svcl(2,1)*spvl(1,2)+svcl(2,2)*spvl(2,2) 
        
       do i=0,3
        eprc(i,1) = cZero
        eprc(i,2) = cZero
        eprc(i,3) = sssl(1,1)*rc(i,3)+sssl(2,1)*rc(i,4)
        eprc(i,4) = sssl(1,2)*rc(i,3)+sssl(2,2)*rc(i,4)
       enddo

       geprc(1) = eprc(0,3)+eprc(1,4)+cImag*eprc(2,4)+eprc(3,3)
       geprc(2) = eprc(0,4)+eprc(1,3)-cImag*eprc(2,3)-eprc(3,4)
       geprc(3) = cZero
       geprc(4) = cZero

       fvorc(1) = fvorc(1)-d*gcc(1)*fmass*geprc(1)
       fvorc(2) = fvorc(2)-d*gcc(1)*fmass*geprc(2)
       fvorc(3) = fvorc(3)-d*gcc(1)*(spfr(1,1)*geprc(1)
     &            +spfr(2,1)*geprc(2))
       fvorc(4) = fvorc(4)-d*gcc(1)*(spfr(1,2)*geprc(1)
     &            +spfr(2,2)*geprc(2))

       sssl(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
       sssl(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
       sssl(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
       sssl(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2) 
        
       do i=0,3
        eprc(i,1) = cZero
        eprc(i,2) = cZero
        eprc(i,3) = sssl(1,1)*rc(i,3)+sssl(2,1)*rc(i,4)
        eprc(i,4) = sssl(1,2)*rc(i,3)+sssl(2,2)*rc(i,4)
       enddo

       geprc(1) = eprc(0,3)+eprc(1,4)+cImag*eprc(2,4)+eprc(3,3)
       geprc(2) = eprc(0,4)+eprc(1,3)-cImag*eprc(2,3)-eprc(3,4)
       geprc(3) = cZero
       geprc(4) = cZero

       fvorc(1) = fvorc(1)+d*gcc(1)*fmass*geprc(1)
       fvorc(2) = fvorc(2)+d*gcc(1)*fmass*geprc(2)
       fvorc(3) = fvorc(3)+d*gcc(1)*(spfr(1,1)*geprc(1)
     &            +spfr(2,1)*geprc(2))
       fvorc(4) = fvorc(4)+d*gcc(1)*(spfr(1,2)*geprc(1)
     &            +spfr(2,2)*geprc(2))
     
      endif

      if ( gcc(2).ne.cZero ) then

       sssr(1,1) = svcr(1,1)*spvr(1,1)+svcr(1,2)*spvr(2,1)
       sssr(1,2) = svcr(1,1)*spvr(1,2)+svcr(1,2)*spvr(2,2) 
       sssr(2,1) = svcr(2,1)*spvr(1,1)+svcr(2,2)*spvr(2,1) 
       sssr(2,2) = svcr(2,1)*spvr(1,2)+svcr(2,2)*spvr(2,2)

       do i=0,3
        eprc(i,1) = sssr(1,1)*rc(i,1)+sssr(2,1)*rc(i,2)
        eprc(i,2) = sssr(1,2)*rc(i,1)+sssr(2,2)*rc(i,2)
        eprc(i,3) = cZero
        eprc(i,4) = cZero
       enddo

       geprc(1) = cZero
       geprc(2) = cZero
       geprc(3) = eprc(0,1)-eprc(1,2)-cImag*eprc(2,2)-eprc(3,1)
       geprc(4) = eprc(0,2)-eprc(1,1)+cImag*eprc(2,1)+eprc(3,2)  

       fvorc(1) = fvorc(1)-d*gcc(2)*(spfl(1,1)*geprc(3)
     &            +spfl(2,1)*geprc(4))
       fvorc(2) = fvorc(2)-d*gcc(2)*(spfl(1,2)*geprc(3)
     &            +spfl(2,2)*geprc(4))
       fvorc(3) = fvorc(3)-d*gcc(2)*fmass*geprc(3)
       fvorc(4) = fvorc(4)-d*gcc(2)*fmass*geprc(4)

       sssr(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
       sssr(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
       sssr(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
       sssr(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

       do i=0,3
        eprc(i,1) = sssr(1,1)*rc(i,1)+sssr(2,1)*rc(i,2)
        eprc(i,2) = sssr(1,2)*rc(i,1)+sssr(2,2)*rc(i,2)
        eprc(i,3) = cZero
        eprc(i,4) = cZero
       enddo

       geprc(1) = cZero
       geprc(2) = cZero
       geprc(3) = eprc(0,1)-eprc(1,2)-cImag*eprc(2,2)-eprc(3,1)
       geprc(4) = eprc(0,2)-eprc(1,1)+cImag*eprc(2,1)+eprc(3,2)  

       fvorc(1) = fvorc(1)+d*gcc(2)*(spfl(1,1)*geprc(3)
     &            +spfl(2,1)*geprc(4))
       fvorc(2) = fvorc(2)+d*gcc(2)*(spfl(1,2)*geprc(3)
     &            +spfl(2,2)*geprc(4))
       fvorc(3) = fvorc(3)+d*gcc(2)*fmass*geprc(3)
       fvorc(4) = fvorc(4)+d*gcc(2)*fmass*geprc(4)

      endif   

      return
      end
