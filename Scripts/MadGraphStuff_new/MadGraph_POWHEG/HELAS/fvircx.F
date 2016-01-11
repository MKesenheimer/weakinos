      subroutine fvircx(ric,vc,gc,fmass,fwidth , fvirc)
c
c This subroutine computes an off-shell antifermion wavefunction from a
c flowing-IN antiRatita-Schwinger fermion and a vector boson.
c
c input:
c       complex ric(18)        : flow-in  antiRSfermion             |ric>
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fvirc(6)        : off-shell fermion             |f',v,ric>
c
c- by Yoshitaro Takaesu - 2011/03/25
c
      implicit none
      double complex ric(18),vc(6),gc(2),gcc(2),fvirc(6),d
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
     & ' helas-warn  : ric in fvircx is zero vector-spinor'
      endif
      if ( abs(ric(17))+abs(ric(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ric in fvircx has zero momentum'
      endif
      if ( abs(vc(1))+abs(vc(2))+abs(vc(3))+abs(vc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc in fvircx is zero vector'
      endif
      if ( abs(vc(5))+abs(vc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc in fvircx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in fvircx is zero coupling'
      endif
      if ( fmass.lt.rZero ) then
         write(stdo,*) ' helas-error : fmass in fvircx is negative'
         write(stdo,*) '             : fmass = ',fmass
      endif
      if ( fwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : fwidth in fvircx is negative'
         write(stdo,*) '             : fwidth = ',fwidth
      endif
#endif

      fvirc(5) = ric(17)-vc(5)
      fvirc(6) = ric(18)-vc(6)

      pf(0) = dble( fvirc(5))
      pf(1) = dble( fvirc(6))
      pf(2) = dimag(fvirc(6))
      pf(3) = dimag(fvirc(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

#ifdef HELAS_CHECK
      if ( abs(fvirc(5))+abs(fvirc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fvirc in fvircx has zero momentum'
      endif
      if ( fwidth.eq.rZero .and. pf2.eq.fmass**2 ) then
         write(stdo,*)
     &        ' helas-error : fvirc in fvircx is on fmass pole'
         write(stdo,*)
     &        '             : p     = ',pf(0),pf(1),pf(2),pf(3)
         write(stdo,*)
     &        '             : abs(p)= ',dsqrt(dabs(pf2))
         fvirc(1) = cZero
         fvirc(2) = cZero
         fvirc(3) = cZero
         fvirc(4) = cZero
         return
      endif
#endif

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1)) 

      pv(0) = dble( vc(5))
      pv(1) = dble( vc(6))
      pv(2) = dimag(vc(6))
      pv(3) = dimag(vc(5))

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

      fvirc(1) = cZero
      fvirc(2) = cZero
      fvirc(3) = cZero
      fvirc(4) = cZero

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

      d = 1d0/dcmplx(pf2-fmass**2,fmass*fwidth)

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

       fvirc(1) = fvirc(1)+d*dconjg(gcc(2))*fmass*geprc(1)
       fvirc(2) = fvirc(2)+d*dconjg(gcc(2))*fmass*geprc(2)
       fvirc(3) = fvirc(3)+d*dconjg(gcc(2))*(spfl(1,1)*geprc(1)
     &           +spfl(1,2)*geprc(2))
       fvirc(4) = fvirc(4)+d*dconjg(gcc(2))*(spfl(2,1)*geprc(1)
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

       fvirc(1) = fvirc(1)-d*dconjg(gcc(2))*fmass*geprc(1)
       fvirc(2) = fvirc(2)-d*dconjg(gcc(2))*fmass*geprc(2)
       fvirc(3) = fvirc(3)-d*dconjg(gcc(2))*(spfl(1,1)*geprc(1)
     &          +spfl(1,2)*geprc(2))
       fvirc(4) = fvirc(4)-d*dconjg(gcc(2))*(spfl(2,1)*geprc(1)
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

       fvirc(1) = fvirc(1)+d*dconjg(gcc(1))*(spfr(1,1)*geprc(3)
     &          +spfr(1,2) *geprc(4))
       fvirc(2) = fvirc(2)+d*dconjg(gcc(1))*(spfr(2,1)*geprc(3)
     &           +spfr(2,2)*geprc(4))
       fvirc(3) = fvirc(3)+d*dconjg(gcc(1))*fmass*geprc(3)
       fvirc(4) = fvirc(4)+d*dconjg(gcc(1))*fmass*geprc(4)

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

       fvirc(1) = fvirc(1)-d*dconjg(gcc(1))*(spfr(1,1)*geprc(3)
     &           +spfr(1,2)*geprc(4))
       fvirc(2) = fvirc(2)-d*dconjg(gcc(1))*(spfr(2,1)*geprc(3)
     &           +spfr(2,2)*geprc(4))
       fvirc(3) = fvirc(3)-d*dconjg(gcc(1))*fmass*geprc(3)
       fvirc(4) = fvirc(4)-d*dconjg(gcc(1))*fmass*geprc(4)

      endif
     
      return
      end
