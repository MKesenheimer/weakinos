      subroutine fvvirc(ric,vc1,vc2,gc,fmass,fwidth , fvvic)
c
c This subroutine computes an off-shell antifermion wavefunction from a
c flowing-IN antiRarita-Schwinger fermion and a vector boson pair.
c
c input:
c       complex ric(18)        : flow-in  antiRSfermion             |ric>
c       complex vc1(6)         : input    vector                     v1
c       complex vc2(6)         : input    vector                     v2
c       complex gc(2)          : coupling constants                  gfrvv
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fvvic(6)        : off-shell antifermion          |f',v1,v2,ric>
c
c- by Yoshitaro Takaesu - 2011/03/27
c
      implicit none
      double complex ric(18),vc1(6),vc2(6),gc(2),fvvic(6),d,gcc(2)
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
      sabsri=0
      do i=1,16
       sabsri=abs(ric(i))+sabsri
      enddo
      if ( sabsri.eq.rZero ) then
         write(stdo,*)
     & ' helas-warn  : ric in fvvirc is zero vector-spinor'
      endif
      if ( abs(ric(17))+abs(ric(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ric in fvvirc has zero momentum'
      endif
      if ( abs(vc1(1))+abs(vc1(2))+abs(vc1(3))+abs(vc1(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc1 in fvvirc is zero vector'
      endif
      if ( abs(vc1(5))+abs(vc1(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc1 in fvvirc has zero momentum'
      endif
      if ( abs(vc2(1))+abs(vc2(2))+abs(vc2(3))+abs(vc2(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc2 in fvvirc is zero vector'
      endif
      if ( abs(vc2(5))+abs(vc2(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc2 in fvvirc has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in fvvirc is zero coupling'
      endif
      if ( fmass.lt.rZero ) then
         write(stdo,*) ' helas-error : fmass in fvvirc is negative'
         write(stdo,*) '             : fmass = ',fmass
      endif
      if ( fwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : fwidth in fvvirc is negative'
         write(stdo,*) '             : fwidth = ',fwidth
      endif
#endif

      fvvic(5) = -ric(17)+vc1(5)+vc2(5)
      fvvic(6) = -ric(18)+vc1(6)+vc2(6)

      pf(0) = dble( fvvic(5))
      pf(1) = dble( fvvic(6))
      pf(2) = dimag(fvvic(6))
      pf(3) = dimag(fvvic(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

#ifdef HELAS_CHECK
      if ( abs(fvvic(5))+abs(fvvic(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fvvic in fvvirc has zero momentum'
      endif
      if ( fwidth.eq.rZero .and. pf2.eq.fmass**2 ) then
         write(stdo,*)
     &        ' helas-error : fvvic in fvvirc is on fmass pole'
         write(stdo,*)
     &        '             : p     = ',pf(0),pf(1),pf(2),pf(3)
         write(stdo,*)
     &        '             : abs(p)= ',dsqrt(dabs(pf2))
         fvvic(1) = cZero
         fvvic(2) = cZero
         fvvic(3) = cZero
         fvvic(4) = cZero
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

      fvvic(1) = cZero
      fvvic(2) = cZero
      fvvic(3) = cZero
      fvvic(4) = cZero

      svcl(1,1) = vc2(1)+vc2(4)
      svcl(1,2) = vc2(2)-cImag*vc2(3)
      svcl(2,1) = vc2(2)+cImag*vc2(3)
      svcl(2,2) = vc2(1)-vc2(4)

      svcr(1,1) = vc2(1)-vc2(4)
      svcr(1,2) = -(vc2(2)-cImag*vc2(3))
      svcr(2,1) = -(vc2(2)+cImag*vc2(3))
      svcr(2,2) = vc2(1)+vc2(4)
        
      spvl(1,1) = vc1(1)-vc1(4)
      spvl(1,2) = -(vc1(2)-cImag*vc1(3))
      spvl(2,1) = -(vc1(2)+cImag*vc1(3))
      spvl(2,2) = vc1(1)+vc1(4)

      spvr(1,1) = vc1(1)+vc1(4)
      spvr(1,2) = vc1(2)-cImag*vc1(3)
      spvr(2,1) = vc1(2)+cImag*vc1(3)
      spvr(2,2) = vc1(1)-vc1(4)

      spfl(1,1) = pf(0)+pf(3)
      spfl(1,2) = pf(1)-cImag*pf(2)
      spfl(2,1) = pf(1)+cImag*pf(2)
      spfl(2,2) = pf(0)-pf(3)

      spfr(1,1) = pf(0)-pf(3)
      spfr(1,2) = -(pf(1)-cImag*pf(2))
      spfr(2,1) = -(pf(1)+cImag*pf(2))
      spfr(2,2) = pf(0)+pf(3)

      d = 1d0/dcmplx(pf2-fmass**2,fmass*fwidth)
      d = -d

      if ( gcc(2).ne.cZero ) then
        
       sssl(1,1) = svcl(1,1)*spvl(1,1)+svcl(1,2)*spvl(2,1)
       sssl(1,2) = svcl(1,1)*spvl(1,2)+svcl(1,2)*spvl(2,2) 
       sssl(2,1) = svcl(2,1)*spvl(1,1)+svcl(2,2)*spvl(2,1) 
       sssl(2,2) = svcl(2,1)*spvl(1,2)+svcl(2,2)*spvl(2,2) 
        
       do i=0,3
        eprc(i,1) = cZero
        eprc(i,2) = cZero
        eprc(i,3) = sssl(1,1)*rc(i,3)+sssl(1,2)*rc(i,4)
        eprc(i,4) = sssl(2,1)*rc(i,3)+sssl(2,2)*rc(i,4)
       enddo

       geprc(1) = eprc(0,3)-eprc(1,4)+cImag*eprc(2,4)-eprc(3,3)
       geprc(2) = eprc(0,4)-eprc(1,3)-cImag*eprc(2,3)+eprc(3,4)
       geprc(3) = cZero
       geprc(4) = cZero

       fvvic(1) = fvvic(1)-d*dconjg(gcc(2))*fmass*geprc(1)
       fvvic(2) = fvvic(2)-d*dconjg(gcc(2))*fmass*geprc(2)
       fvvic(3) = fvvic(3)-d*dconjg(gcc(2))*(spfl(1,1)*geprc(1)
     &           +spfl(1,2)*geprc(2))
       fvvic(4) = fvvic(4)-d*dconjg(gcc(2))*(spfl(2,1)*geprc(1)
     &           +spfl(2,2)*geprc(2))

       sssl(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
       sssl(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
       sssl(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
       sssl(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

       do i=0,3
        eprc(i,1) = cZero
        eprc(i,2) = cZero
        eprc(i,3) = sssl(1,1)*rc(i,3)+sssl(1,2)*rc(i,4)
        eprc(i,4) = sssl(2,1)*rc(i,3)+sssl(2,2)*rc(i,4)
       enddo

       geprc(1) = eprc(0,3)-eprc(1,4)+cImag*eprc(2,4)-eprc(3,3)
       geprc(2) = eprc(0,4)-eprc(1,3)-cImag*eprc(2,3)+eprc(3,4)
       geprc(3) = cZero
       geprc(4) = cZero  

       fvvic(1) = fvvic(1)+d*dconjg(gcc(2))*fmass*geprc(1)
       fvvic(2) = fvvic(2)+d*dconjg(gcc(2))*fmass*geprc(2)
       fvvic(3) = fvvic(3)+d*dconjg(gcc(2))*(spfl(1,1)*geprc(1)
     &          +spfl(1,2)*geprc(2))
       fvvic(4) = fvvic(4)+d*dconjg(gcc(2))*(spfl(2,1)*geprc(1)
     &          +spfl(2,2)*geprc(2))
     
      endif

      if ( gcc(1).ne.cZero ) then

       sssr(1,1) = svcr(1,1)*spvr(1,1)+svcr(1,2)*spvr(2,1)
       sssr(1,2) = svcr(1,1)*spvr(1,2)+svcr(1,2)*spvr(2,2) 
       sssr(2,1) = svcr(2,1)*spvr(1,1)+svcr(2,2)*spvr(2,1) 
       sssr(2,2) = svcr(2,1)*spvr(1,2)+svcr(2,2)*spvr(2,2)

       do i=0,3
        eprc(i,1) = sssr(1,1)*rc(i,1)+sssr(1,2)*rc(i,2)
        eprc(i,2) = sssr(2,1)*rc(i,1)+sssr(2,2)*rc(i,2)
        eprc(i,3) = cZero
        eprc(i,4) = cZero
       enddo

       geprc(1) = cZero
       geprc(2) = cZero
       geprc(3) = eprc(0,1)+eprc(1,2)-cImag*eprc(2,2)+eprc(3,1)
       geprc(4) = eprc(0,2)+eprc(1,1)+cImag*eprc(2,1)-eprc(3,2)  

       fvvic(1) = fvvic(1)-d*dconjg(gcc(1))*(spfr(1,1)*geprc(3)
     &          +spfr(1,2) *geprc(4))
       fvvic(2) = fvvic(2)-d*dconjg(gcc(1))*(spfr(2,1)*geprc(3)
     &           +spfr(2,2)*geprc(4))
       fvvic(3) = fvvic(3)-d*dconjg(gcc(1))*fmass*geprc(3)
       fvvic(4) = fvvic(4)-d*dconjg(gcc(1))*fmass*geprc(4)

       sssr(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
       sssr(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
       sssr(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
       sssr(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

       do i=0,3
        eprc(i,1) = sssr(1,1)*rc(i,1)+sssr(1,2)*rc(i,2)
        eprc(i,2) = sssr(2,1)*rc(i,1)+sssr(2,2)*rc(i,2)
        eprc(i,3) = cZero
        eprc(i,4) = cZero
       enddo

       geprc(1) = cZero
       geprc(2) = cZero
       geprc(3) = eprc(0,1)+eprc(1,2)-cImag*eprc(2,2)+eprc(3,1)
       geprc(4) = eprc(0,2)+eprc(1,1)+cImag*eprc(2,1)-eprc(3,2)  

       fvvic(1) = fvvic(1)+d*dconjg(gcc(1))*(spfr(1,1)*geprc(3)
     &           +spfr(1,2)*geprc(4))
       fvvic(2) = fvvic(2)+d*dconjg(gcc(1))*(spfr(2,1)*geprc(3)
     &           +spfr(2,2)*geprc(4))
       fvvic(3) = fvvic(3)+d*dconjg(gcc(1))*fmass*geprc(3)
       fvvic(4) = fvvic(4)+d*dconjg(gcc(1))*fmass*geprc(4)

      endif
     
      return
      end
