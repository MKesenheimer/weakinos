c
c ----------------------------------------------------------------------
c
      subroutine fvorgx(fo,vc,gc,fmass,fwidth, fvorg)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-out external fermion and a vector boson, for the NLSP-boson-
c Goldstino vertex.
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fvorg(6)         : off-shell fermion             <fo,v,f'|
c
c- by Yoshitaro Takaesu - 2010/11/15
c

      implicit none
      double complex  fo(6), vc(6), gc(2), fvorg(6), d, svcl(2,2)
      double complex  svcr(2,2), spvl(2,2), spvr(2,2), spfl(2,2)
      double complex  spfr(2,2), sssl(2,2), sssr(2,2), ssml(2,2)
      double complex ssmr(2,2), qvpl(2,2), qvpr(2,2), vqpl(2,2)
      double complex vqpr(2,2)
      double precision  fmass, fwidth
      double precision  pf(0:3), pv(0:3), pf2, pdotpG

      double precision rOne
      parameter( rOne = 1.0d0 )
      double complex ci, cZero
      parameter( ci = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      fvorg(1) = cZero
      fvorg(2) = cZero
      fvorg(3) = cZero
      fvorg(4) = cZero
      fvorg(5) = fo(5) + vc(5)
      fvorg(6) = fo(6) + vc(6)

      pv(0) = dble( vc(5))
      pv(1) = dble( vc(6))
      pv(2) = dimag(vc(6))
      pv(3) = dimag(vc(5))

      pf(0) = dble( fvorg(5))
      pf(1) = dble( fvorg(6))
      pf(2) = dimag(fvorg(6))
      pf(3) = dimag(fvorg(5))
      pf2 = pf(0)**2 - pf(1)**2 - pf(2)**2 - pf(3)**2

      svcl(1,1) = vc(1)-vc(4)
      svcl(1,2) = -(vc(2)-ci*vc(3))
      svcl(2,1) = -(vc(2)+ci*vc(3))
      svcl(2,2) = vc(1)+vc(4)

      svcr(1,1) = vc(1)+vc(4)
      svcr(1,2) = vc(2)-ci*vc(3)
      svcr(2,1) = vc(2)+ci*vc(3)
      svcr(2,2) = vc(1)-vc(4)
        
      spvl(1,1) = pv(0)-pv(3)
      spvl(1,2) = -(pv(1)-ci*pv(2))
      spvl(2,1) = -(pv(1)+ci*pv(2))
      spvl(2,2) = pv(0)+pv(3)

      spvr(1,1) = pv(0)+pv(3)
      spvr(1,2) = pv(1)-ci*pv(2)
      spvr(2,1) = pv(1)+ci*pv(2)
      spvr(2,2) = pv(0)-pv(3)

      spfl(1,1) = pf(0)-pf(3)
      spfl(1,2) = -(pf(1)-ci*pf(2))
      spfl(2,1) = -(pf(1)+ci*pf(2))
      spfl(2,2) = pf(0)+pf(3)

      spfr(1,1) = pf(0)+pf(3)
      spfr(1,2) = pf(1)-ci*pf(2)
      spfr(2,1) = pf(1)+ci*pf(2)
      spfr(2,2) = pf(0)-pf(3)

      d = -rOne/dcmplx(pf2-fmass**2,fmass*fwidth)
      
      sssl(1,1) = spvl(1,1)*svcr(1,1)+spvl(1,2)*svcr(2,1)
      sssl(1,2) = spvl(1,1)*svcr(1,2)+spvl(1,2)*svcr(2,2) 
      sssl(2,1) = spvl(2,1)*svcr(1,1)+spvl(2,2)*svcr(2,1) 
      sssl(2,2) = spvl(2,1)*svcr(1,2)+spvl(2,2)*svcr(2,2)

      sssr(1,1) = spvr(1,1)*svcl(1,1)+spvr(1,2)*svcl(2,1)
      sssr(1,2) = spvr(1,1)*svcl(1,2)+spvr(1,2)*svcl(2,2) 
      sssr(2,1) = spvr(2,1)*svcl(1,1)+spvr(2,2)*svcl(2,1) 
      sssr(2,2) = spvr(2,1)*svcl(1,2)+spvr(2,2)*svcl(2,2)

      ssml(1,1) = svcl(1,1)*spvr(1,1)+svcl(1,2)*spvr(2,1)
      ssml(1,2) = svcl(1,1)*spvr(1,2)+svcl(1,2)*spvr(2,2) 
      ssml(2,1) = svcl(2,1)*spvr(1,1)+svcl(2,2)*spvr(2,1) 
      ssml(2,2) = svcl(2,1)*spvr(1,2)+svcl(2,2)*spvr(2,2)
      
      ssmr(1,1) = svcr(1,1)*spvl(1,1)+svcr(1,2)*spvl(2,1)
      ssmr(1,2) = svcr(1,1)*spvl(1,2)+svcr(1,2)*spvl(2,2) 
      ssmr(2,1) = svcr(2,1)*spvl(1,1)+svcr(2,2)*spvl(2,1) 
      ssmr(2,2) = svcr(2,1)*spvl(1,2)+svcr(2,2)*spvl(2,2)

      qvpl(1,1) = sssl(1,1)*spfl(1,1)+sssl(1,2)*spfl(2,1)
      qvpl(1,2) = sssl(1,1)*spfl(1,2)+sssl(1,2)*spfl(2,2) 
      qvpl(2,1) = sssl(2,1)*spfl(1,1)+sssl(2,2)*spfl(2,1) 
      qvpl(2,2) = sssl(2,1)*spfl(1,2)+sssl(2,2)*spfl(2,2)

      qvpr(1,1) = sssr(1,1)*spfr(1,1)+sssr(1,2)*spfr(2,1)
      qvpr(1,2) = sssr(1,1)*spfr(1,2)+sssr(1,2)*spfr(2,2) 
      qvpr(2,1) = sssr(2,1)*spfr(1,1)+sssr(2,2)*spfr(2,1) 
      qvpr(2,2) = sssr(2,1)*spfr(1,2)+sssr(2,2)*spfr(2,2)

      vqpl(1,1) = ssml(1,1)*spfl(1,1)+ssml(1,2)*spfl(2,1)
      vqpl(1,2) = ssml(1,1)*spfl(1,2)+ssml(1,2)*spfl(2,2) 
      vqpl(2,1) = ssml(2,1)*spfl(1,1)+ssml(2,2)*spfl(2,1) 
      vqpl(2,2) = ssml(2,1)*spfl(1,2)+ssml(2,2)*spfl(2,2)

      vqpr(1,1) = ssmr(1,1)*spfr(1,1)+ssmr(1,2)*spfr(2,1)
      vqpr(1,2) = ssmr(1,1)*spfr(1,2)+ssmr(1,2)*spfr(2,2) 
      vqpr(2,1) = ssmr(2,1)*spfr(1,1)+ssmr(2,2)*spfr(2,1) 
      vqpr(2,2) = ssmr(2,1)*spfr(1,2)+ssmr(2,2)*spfr(2,2)

      fvorg(1) = d*gc(1)*( fo(1)*fmass*(sssl(1,1)-ssml(1,1))
     & +fo(2)*fmass*(sssl(2,1)-ssml(2,1)) )
      fvorg(2) = d*gc(1)*( fo(1)*fmass*(sssl(1,2)-ssml(1,2))
     & +fo(2)*fmass*(sssl(2,2)-ssml(2,2)) )
      fvorg(3) = d*gc(1)*( fo(1)*(qvpl(1,1)-vqpl(1,1))
     & +fo(2)*(qvpl(2,1)-vqpl(2,1)) )
      fvorg(4) = d*gc(1)*( fo(1)*(qvpl(1,2)-vqpl(1,2))
     & +fo(2)*(qvpl(2,2)-vqpl(2,2)) )

      if (gc(2) .ne. 0d0) then
         fvorg(1) = fvorg(1)+d*gc(2)*( fo(3)*(qvpr(1,1)-vqpr(1,1))
     &        +fo(4)*(qvpr(2,1)-vqpr(2,1)) )
         fvorg(2) = fvorg(2)+d*gc(2)*( fo(3)*(qvpr(1,2)-vqpr(1,2))
     &        +fo(4)*(qvpr(2,2)-vqpr(2,2)) )
         fvorg(3) = fvorg(3)+d*gc(2)*( fo(3)*fmass*(sssr(1,1)-ssmr(1,1))
     &        +fo(4)*fmass*(sssr(2,1)-ssmr(2,1)) )
         fvorg(4) = fvorg(4)+d*gc(2)*( fo(3)*fmass*(sssr(1,2)-ssmr(1,2))
     &        +fo(4)*fmass*(sssr(2,2)-ssmr(2,2)) )
      endif
      
      return          
      end
