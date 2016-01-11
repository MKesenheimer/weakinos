ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c FINITE PART OF THE INTEGRATE DIPOLES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Returns fi(8):
c fi(1), regular piece at z!=1
c fi(2), delta piece at z=1
c fi(3), plus distribution at z!=1
c fi(4), plus distribution at z=1
c fi(5), mass correction to plus distribution at z!=1
c fi(6), mass correction to plus distribution at z=1
c fi(7), delta piece at z=1-4m^2/s
c fi(8), plus distribution at z!=1-4m^2/s
c fi(9), plus distribution at z=1-4m^2/s
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      SUBROUTINE finiteiiqed(sik,sikzone,x,id1,id2,qint1,qint2,leg,pdg,fi)
c  calculates the finite terms when both emitter
c  and spectator are in the initial state.
      implicit none
      include "coupl.inc"
      include "dipole.inc"

c Arguments
      REAL*8 sik,sikzone,x,fi(9)
      INTEGER id1,id2,qint1,qint2,pdg,leg
c Global
      REAL*8 ddilog
      external ddilog
c Local
      INTEGER i,qint(2),nc
      REAl*8 gammaE,pi,gsq,mass,Pff,Pfg,Pgf,Cff,Cfg,Cgf
      REAL*8 L,L_one,rs,qem(2)
      PARAMETER (gammaE=0.57721566d0)
      PARAMETER (pi=3.1415926535897932385d0)
      LOGICAL pdf

      gsq = Dreal(gal(1))**2
c  extract the electromagnetic charge out of integers
      qint(1)=qint1
      qint(2)=qint2
      do i=1,2
       if(qint(i).eq.23) then
         qem(i)=2d0/3d0
       elseif(qint(i).eq.-23) then
         qem(i)=-2d0/3d0
       elseif(qint(i).eq.13) then
         qem(i)=1d0/3d0
       elseif(qint(i).eq.-13) then
         qem(i)=-1d0/3d0
       elseif(qint(i).eq.1) then
         qem(i)=1d0
       elseif(qint(i).eq.-1) then
         qem(i)=-1d0
       endif
      enddo

c determine the color factor for final-final a -> f fbar splitting
      if(abs(qem(1)).eq.1d0.or.id2.eq.1) then
        nc=1d0
      elseif(abs(qem(1)).ne.1.and.id1.eq.0) then
        nc=3d0
      endif

      gsq=nc*qem(1)*qem(2)*gsq

c determine whether pdf is present
      if(abs(pdg).lt.10.or.pdg.eq.21) then
        pdf=.true.
      else
        pdf=.false.
      endif

  
      mass= massint(leg)


      L=dlog(mu**2/sik)
      L_one=dlog(mu**2/sikzone)

      do i=1,9
       fi(i)=0d0
      enddo

      Pff=(1d0+x**2)/(1d0-x)
      Pfg=x**2+(1d0-x)**2
      Pgf=(1d0+(1d0-x)**2)/x
      Cff=Pff*(dlog((1d0-x)/x)-3d0/4d0)+(9d0+5*x)/4d0
      Cfg=Pfg*dlog((1d0-x)/x)-8d0*x**2+8d0*x-1d0
      Cgf=-Cff

c fermion-fermion splitting (i:fermion, ij~: fermion)
      if(id1 .eq.1.and.id2.eq.1) then
        if(scheme .eq. 'HV') then
          rs=0d0
        elseif(scheme .eq. 'DR') then
          rs=-1d0/2d0
        endif

        if(rscheme.eq.'MREG') then
         if(pdf) then
          fi(2)=gsq/(8.*pi**2)*(-pi**2/3.+2d0)+gsq/(8.*pi**2)*((1d0 + Dlog(mass**2/sikzone))*DLog(mass_a**2/sikzone)
     &       -0.5d0*DLog(mass**2/sikzone)**2+0.5d0*Dlog(mass**2/sikzone))
          fi(3)=gsq/(8.*pi**2)*( Pff*(DLOG(sikzone/muf**2)-1d0
     &         +2d0*DLOG(1d0-x)+1d0 -Cff))+gsq/(8.*pi**2)*(1d0-x)
          fi(4)=fi(3)
         else
          fi(2)=gsq/(8.*pi**2)*(-pi**2/3.+2d0)+gsq/(8.*pi**2)*((1d0 + Dlog(mass**2/sikzone))*DLog(mass_a**2/sikzone)
     &       -0.5d0*DLog(mass**2/sikzone)**2+0.5d0*Dlog(mass**2/sikzone))
          fi(3)=gsq/(8.*pi**2)*((1d0+x**2)/(1d0-x)*(DLOG(sikzone/mass**2)-1d0))+gsq/(8.*pi**2)*(1d0-x)
          fi(4)=fi(3)
c          fi(2)=fi(2)+gsq/(8.*pi**2)*((-7 + 8*alpha_ii - alpha_ii**2 - 6*Log(alpha_ii) - 4*Log(alpha_ii)**2)/4.)



         endif
c     adding the alpha dependent term
          if((1d0-x).gt.alpha_ii) then
             fi(1)=(gsq*(1+x**2)*DLog(alpha_ii/(1-x)))/
     &           (8.*pi**2*(1-x))
          endif
      
        else
          fi(2)=(gsq*(2*L_one**2-pi**2+4*rs+6*DLog(mu**2/muf**2)))/
     &        (32.*pi**2)
          fi(1)=-(gsq*(-1 + x + (-1 - x)*(L -dLog(mu**2/muf**2)) 
     &        + 2*(1 + x)*dLog(1 - x)))/(8.*pi**2)
          fi(4)=-(gsq*((2*(L_one - dLog(mu**2/muf**2)))/(1 - x)
     &        + (4*dLog(1 - x))/(-1 + x)))/(8.*pi**2)
          fi(3)=-(gsq*((2*(L - dLog(mu**2/muf**2)))/(1 - x)
     &        + (4*dLog(1 - x))/(-1 + x)))/(8.*pi**2)

c     adding the alpha dependent term
          if((1d0-x).gt.alpha_ii) then
             fi(1)=fi(1)+(gsq*(1+x**2)*DLog(alpha_ii/(1-x)))/
     &           (8.*pi**2*(1-x))
          endif
        endif



c fermion-photon splitting (i:fermion, ij~:photon)
      elseif(id1.eq.1.and.id2.eq.0) then

       if(rscheme.eq.'MREG') then

        if(pdf) then
         fi(1)=gsq/(8.*pi**2)*((DLog(sikzone*(1d0-x)**2/(x**2*muf**2))
     &          +2d0*DLog(x)+1d0)*Pgf -Cgf-2d0*(1d0-x)/x)
        else
         fi(1)=gsq/(8.*pi**2)*(DLog(sikzone*(1d0-x)**2/(x**2*mass**2))
     &         *Pgf-2d0*(1d0-x)/x)
        endif
c     adding the alpha dependent term
         if((1d0-x).gt.alpha_ii) then
            fi(1)=fi(1)+(gsq*(1+(1-x)**2)*DLog(alpha_ii/(1-x)))/
     &           (8.*pi**2*x)
         endif

       else

         fi(1)=(gsq*(x**2 - L*(2 + (-2 + x)*x)
     &        +(2 + (-2 + x)*x)*dLog(mu**2/muf**2) + 
     &        2*(2 + (-2 + x)*x)*DLog(1 - x)))/(4.*pi**2*x)
c     adding the alpha dependent term
         if((1d0-x).gt.alpha_ii) then
            fi(1)=fi(1)+(gsq*(1+(1-x)**2)*DLog(alpha_ii/(1-x)))/
     &           (8.*pi**2*x)
         endif
       endif


c photon-fermion splitting (i:photon, ij~:fermion)
      elseif(id1.eq.0.and.id2.eq.1) then
       if(rscheme.eq.'MREG') then
        if(photonpdf) then
         fi(1)=gsq/(8.*pi**2)*( Pfg*Dlog(sikzone*(1d0-x)**2/muf**2)
     &         +2d0*x*(1d0-x) -Cfg)
        else
         fi(1)=gsq/(8.*pi**2)*( Pfg*Dlog(sikzone*(1d0-x)**2/mass**2)
     &         +2d0*x*(1d0-x) )
        endif
       else
         fi(1)=(gsq*(-L + 2*(1 + L)*x - 2*(1 + L)*x**2 
     &        + (1 + 2*(-1 + x)*x)*DLog(mu**2/muf**2) + 
     &        (2 + 4*(-1 + x)*x)*DLog(1 - x)))/(16.*pi**2)
       endif
c adding the alpha dependent term
         if((1d0-x).gt.alpha_ii) then
            fi(1)=fi(1)+(gsq*((1-x)**2+x**2)*DLog(alpha_ii/(1-x)))/
     &           (8.*pi**2)
         endif

      endif

      end


 
      SUBROUTINE finiteifqed(mk,sik,sikzone,x,id1,id2,qint1,qint2,leg,pdg,sab,sabt,fi)
c  calculates the finite terms when the emitter is
c  in the initial state and the spectator is in the
c  final state.
      implicit none
      include "coupl.inc"
      include "dipole.inc"

c Arguments
      REAL*8 mk,sik,sikzone,x,fi(9)
      INTEGER id1,id2,qint1,qint2,leg,pdg
c Global
      REAL*8 ddilog,lambda_tr
      external ddilog,lambda_tr
c Local
      INTEGER i,qint(2),nc
      REAl*8 gammaE,pi,gsq,qem(2),L_one,mass
      REAL*8 zp,musq_k,L,rs,muk,musq_k_one,muk_one
      REAL*8 massb,Piasq,Piasqbar,lambda_ia,gamma,Pfg,Pgf,Pff,Cff,Cfg,Cgf
      REAL*8 rhoia,Ria,z1,z2,c0,c1,c2,c3,c4,c5,b0,x0,sab,sabt
      REAL*8 Piasq_one, Piasqbar_one, lambda_ia_one, gamma_one
      REAL*8 Ria_one, z1_one,z2_one
      REAL*8 A1,y1,y2,y3,y4,y5,y6,yp,A, y7,y8,B1,c6,c61,c62
      PARAMETER (gammaE=0.57721566d0)
      PARAMETER (pi=3.1415926535897932385d0)
      LOGICAL pdf

      gsq = Dreal(gal(1))**2
c  extract the electromagnetic charge out of integers
      qint(1)=qint1
      qint(2)=qint2
      do i=1,2
       if(qint(i).eq.23) then
         qem(i)=2d0/3d0
       elseif(qint(i).eq.-23) then
         qem(i)=-2d0/3d0
       elseif(qint(i).eq.13) then
         qem(i)=1d0/3d0
       elseif(qint(i).eq.-13) then
         qem(i)=-1d0/3d0
       elseif(qint(i).eq.1) then
         qem(i)=1d0
       elseif(qint(i).eq.-1) then
         qem(i)=-1d0
       endif
      enddo

c determine the color factor for final-final a -> f fbar splitting
      if(abs(qem(1)).eq.1d0) then
        nc=1d0
      elseif(abs(qem(1)).ne.1.and.id1.eq.0) then
        nc=3d0
      endif

      gsq=nc*qem(1)*qem(2)*gsq
      musq_k=mk**2/sik
      musq_k_one=mk**2/sikzone
      muk=Sqrt(musq_k)
      muk_one=Sqrt(musq_k)
      zp=(1d0-x)/(1d0-x+musq_k)
      L=DLog(mu**2/sik)
      L_one=DLog(mu**2/sikzone)
      x0=0d0

c determine whether pdf is present
      if(abs(pdg).lt.10.or.pdg.eq.21) then
        pdf=.true.
      else
        pdf=.false.
      endif

      mass=massint(leg)

      do i=1,9
       fi(i)=0d0
      enddo

      Pff=(1d0+x**2)/(1d0-x)
      Pfg=x**2+(1d0-x)**2
      Pgf=(1d0+(1d0-x)**2)/x
      Cff=Pff*(dlog((1d0-x)/x)-3d0/4d0)+(9d0+5*x)/4d0
      Cfg=Pfg*dlog((1d0-x)/x)-8d0*x**2+8d0*x-1d0
      Cgf=-Cff

c fermion-fermion splitting (i:fermion, ij~: fermion)
      if(id1 .eq.1.and.id2.eq.1) then
        if(scheme .eq. 'HV') then
          rs=0d0
        elseif(scheme .eq. 'DR') then
          rs=-1d0/2d0
        endif

        Piasqbar=-sik
        Piasqbar_one=-sikzone


        Piasq=Piasqbar+mk**2+mass**2
        Piasq_one=Piasqbar_one+mk**2+mass**2

        lambda_ia=lambda_tr(Piasq,mass**2,mk**2)
        lambda_ia_one=lambda_tr(Piasq_one,mass**2,mk**2)

        Ria=sqrt((piasqbar+2*mass**2*x)**2-4d0*mass**2*piasq*x**2)/sqrt(lambda_ia)
        Ria_one=sqrt((piasqbar_one+2*mass**2*x)**2-4d0*mass**2*piasq_one*x**2)/sqrt(lambda_ia_one)
        rhoia=sqrt(lambda_tr(sabt,massint(1)**2,massint(2)**2)/
     &        lambda_tr(sab,massint(1)**2,massint(2)**2))
        z1=(piasqbar*(piasqbar-x*(piasqbar+2*mk**2)) - sqrt(piasqbar**2*(1d0-x)**2
     &     -4d0*mk**2*mass_a**2*x**2)*sqrt(lambda_ia)*Ria)/(2d0*piasqbar*(piasqbar
     &     -x*(piasq-mass**2)))
        z1_one=(piasqbar_one*(piasqbar_one-x*(piasqbar_one+2*mk**2)) - sqrt(piasqbar_one**2*(1d0-x)**2
     &     -4d0*mk**2*mass_a**2*x**2)*sqrt(lambda_ia_one)*Ria_one)/(2d0*piasqbar_one*(piasqbar_one
     &     -x*(piasq_one-mass**2)))
        z2=(piasqbar*(piasqbar-x*(piasqbar+2*mk**2)) + sqrt(piasqbar**2*(1d0-x)**2
     &     -4d0*mk**2*mass_a**2*x**2)*sqrt(lambda_ia)*Ria)/(2d0*piasqbar*(piasqbar
     &     -x*(piasq-mass**2)))
        z2_one=(piasqbar_one*(piasqbar_one-x*(piasqbar_one+2*mk**2)) + sqrt(piasqbar_one**2*(1d0-x)**2
     &     -4d0*mk**2*mass_a**2*x**2)*sqrt(lambda_ia_one)*Ria_one)/(2d0*piasqbar_one*(piasqbar_one
     &     -x*(piasq_one-mass**2)))
        if(z2_one.ge.1d0) then
         z2_one=1d0-((mass_a**2 + mass**2*(-1 + x)**2)*x)/(piasq_one*(-1 + x))
         endif


        c0=(piasqbar_one+sqrt(lambda_ia_one))/(piasqbar_one-sqrt(lambda_ia_one))
c        if(piasqbar_one+sqrt(lambda_ia_one).le.0d0.or.c0.le.0d0) then
c         c0=(-sikzone*(2.*musq_k_one*mass**2/sikzone))/(piasqbar_one-sqrt(lambda_ia_one))
c        endif
        c1=(2d0*mk**2-piasqbar_one-sqrt(lambda_ia_one))/(2d0*mk**2-piasqbar_one+sqrt(lambda_ia_one))
        b0=-4d0*mass**2*piasqbar_one/(lambda_ia_one*(1d0+abs(piasqbar_one)/sqrt(lambda_ia_one))**2
     &     +4d0*mass**2*(piasqbar_one+mk**2))
        c2=-(2*mass**2+piasqbar_one-sqrt(lambda_ia_one))/(2d0*mass**2)*b0
        c3=-(2*mass**2+piasqbar_one+sqrt(lambda_ia_one))/(2d0*mass**2)*b0
        c4= 2d0*(mk**2-mass**2-sqrt(lambda_ia_one))/(2d0*mass**2-piasqbar_one-sqrt(lambda_ia_one))*b0
        c5= 2d0*(mk**2-mass**2+sqrt(lambda_ia_one))/(2d0*mass**2-piasqbar_one+sqrt(lambda_ia_one))*b0
        gamma=mass/(sqrt(-piasqbar_one-mk**2))

        A=Sqrt(-piasqbar_one-mk**2)
        A1=2d0*piasqbar_one/sqrt(4d0*A**2*mass**2+4d0*mass**2*piasqbar_one+piasqbar_one**2)
        y1=(sqrt(4d0*A**2*mass**2+4d0*mass**2*piasqbar_one+piasqbar_one**2)+2d0*A*mass)/(-piasqbar_one)
        y2=(sqrt(4d0*A**2*mass**2+4d0*mass**2*piasqbar_one+piasqbar_one**2)-2d0*A*mass)/(piasq_one)
        y3=(2d0*mass-A)/A
        y4=(3d0*mass+A)/(A-mass)
        y5=(2d0*mass+A)/A
        y6=(3d0*mass-A)/(A+mass)
        y7=-((2*A*mass + Sqrt(4*A**2*mass**2 + 4*mass**2*Piasqbar_one + Piasqbar_one**2))/Piasqbar_one)
        y8=(-2*A*mass + Sqrt(4*A**2*mass**2 + Piasqbar_one*(4*mass**2 + Piasqbar_one)))/Piasqbar_one
        B1=(2*Piasqbar_one)/Sqrt(4*A**2*mass**2 + Piasqbar_one*(4*mass**2 + Piasqbar_one))
        c6=dSqrt(Piasqbar_one*(Piasqbar_one - 4*mass**2*(-1 + alpha_if)) + 4*A**2*mass**2*(-1 + alpha_if)**2)
        c62=c6 + Piasqbar_one + 2*mass**2*(-1 + alpha_if)
        if(c6+Piasqbar_one.le.0d0) then
          c62=0d0
        endif
        c61=1d0 - (2*c6)/(c6 + Piasqbar_one)
        if(c61.lt.0d0) then
          c61=-mass**2/Piasqbar_one*(1d0-alpha_if)
        endif
        yp=-2d0*mass*A/piasqbar_one*(1d0-alpha_if)-sqrt((piasqbar_one
     &      +2*mass**2*(1d0-alpha_if))**2-4d0*mass**2*piasq_one*(1d0-alpha_if)**2)/piasqbar_one

        if(x.gt.1d0-2.*mk*mass_a/abs(Piasqbar)) then
          return
        endif

       if(rscheme.eq.'MREG') then
c massive case: mk>0
        if(mk.gt.0d0) then


          fi(2)=gsq/(8.*pi**2)*(2d0*dlog(mass_a*mk/(abs(piasqbar_one)))+piasqbar_one**2*(1d0-
     &          abs(piasqbar_one)/sqrt(lambda_ia_one))/(4d0*mass**2*(piasqbar_one+mk**2))+2d0 
     &          +piasqbar_one**2*(3d0*piasqbar_one+2d0*mk**2)/(2d0*sqrt(lambda_ia_one)*(piasqbar_one+mk**2)**2)
     &          *(1d0/gamma*dlog((piasqbar_one*gamma**2+2d0*mass**2+gamma*sqrt(lambda_ia_one))/
     &          (piasqbar_one*gamma**2+gamma*abs(piasqbar_one)))+dlog(-piasqbar_one/mk**2)
     &          +dlog((-sikzone*(2.*musq_k_one*mass**2/sikzone)-2d0*piasqbar_one*gamma**2-2d0*mass**2)/(piasqbar_one*
     &          (1d0-2d0*gamma**2)+abs(piasqbar_one))))
     &          +piasqbar_one/sqrt(lambda_ia_one)*(2d0*dlog(mass_a*mass/b0/sqrt(lambda_ia_one))*
     &          dlog(c1/c0)-dlog((mass**2+mk**2-piasqbar_one)/mass**2)*dlog(c1)+0.5d0*dlog(c0*c1)*dlog(c1/c0)
     &          -0.5d0*dlog(c0)-2d0*(ddilog(c0)-ddilog(c1)+ddilog(c2)-ddilog(c3)+ddilog(c4)-ddilog(c5))))
          fi(1)=0d0
          fi(3)=rhoia/x*gsq/(8.*pi**2)*(-piasqbar/sqrt(lambda_ia)/Ria*(2d0/(1d0-x)*
     &          dlog((1d0-z1)*(2d0-x-z2)/((1d0-z2)*(2d0-x-z1)))+Ria*(1d0+x)*dlog((1d0-z2)/(1d0-z1))
     &          +2d0*mass**2*x**2/piasqbar*(1d0/(1d0-z2)-1d0/(1d0-z1))))
          fi(4)=gsq/(8.*pi**2)*(-piasqbar_one/sqrt(lambda_ia_one)/Ria_one*(2d0/(1d0-x)*
     &          dlog((1d0-z1_one)*(2d0-x-z2_one)/((1d0-z2_one)*(2d0-x-z1_one)))+Ria_one*(1d0+x)*dlog((1d0-z2_one)/(1d0-z1_one))
     &          +2d0*mass**2*x**2/piasqbar_one*(1d0/(1d0-z2_one)-1d0/(1d0-z1_one))))

         if(pdf) then
          fi(3)=fi(3)-gsq/(8.*pi**2)*( Pff*(dlog(muf**2/mass**2)+2d0*dlog(1d0-x)+1d0) -Cff)
          fi(4)=fi(4)-gsq/(8.*pi**2)*( Pff*(dlog(muf**2/mass**2)+2d0*dlog(1d0-x)+1d0) -Cff)
         endif

c adding the alpha dependent term
         if(alpha_if.ne.1d0) then
          fi(2)=fi(2)-gsq/(8.*pi**2)*((-1d0)*((A1*(-2*DLog(-(((-1 + y1)*(y2 - y4))/((-1 + y2)
     &          *(y1 - y4))))*DLog(-1 + y4) + DLog((-1 + y4)/(-y1 + y4))**2 
     &          + 2*DLog(-(((-1 + y1)*(y2 - y5))/((-1 + y2)*(y1 - y5))))*DLog(-1 + y5)
     &          - DLog((-1 + y5)/(-y1 + y5))**2 + 
     &          2*DLog(-1 + (2*A)/(A + mass))*DLog(((-1 + y2)*(y1 - yp))
     &          /((-1 + y1)*(y2 - yp))) + 2*DLog(((-y2 + y4)*(y1 - yp))
     &          /((y1 - y4)*(y2 - yp)))*DLog(y4 - yp) - 2*DLog(((-y2 + y5)
     &          *(y1 - yp))/((y1 - y5)*(y2 - yp)))*DLog(y5 - yp) - 
     &          DLog((-y4 + yp)/(y1 - y4))**2 + DLog((-y5 + yp)/(y1 - y5))**2 
     &          + 2*ddilog((-y1 + y4)/(-1 + y4)) + 2*ddilog((-1 + y4)/(-y2 + y4))
     &          - 2*ddilog((-y1 + y5)/(-1 + y5)) - 2*ddilog((-1 + y5)/(-y2 + y5)) - 
     &          2*ddilog((-y1 + y4)/(y4 - yp)) + 2*ddilog((-y1 + y5)/(y5 - yp)) 
     &          - 2*(DLog(1 - y3)*DLog(((-1 + y1)*(y2 - y3))/((-1 + y2)*(y1 - y3))) 
     &          + DLog(((y1 - y3)*(y2 - yp))/((y2 - y3)*(y1 - yp)))*DLog(-y3 + yp) + 
     &          ddilog((-1 + y3)/(-y1 + y3)) - ddilog((-1 + y3)/(-y2 + y3)) 
     &          - ddilog((-y3 + yp)/(y1 - y3)) + ddilog((-y3 + yp)/(y2 - y3))) 
     &          - 2*ddilog((-y4 + yp)/(y2 - y4)) + 2*ddilog((-y5 + yp)/(y2 - y5)) + 
     &          2*(DLog(1 - y6)*DLog(((-1 + y1)*(y2 - y6))/((-1 + y2)*(y1 - y6))) 
     &          + DLog(((y1 - y6)*(y2 - yp))/((y2 - y6)*(y1 - yp)))*DLog(-y6 + yp)
     &          + ddilog((-1 + y6)/(-y1 + y6)) - ddilog((-1 + y6)/(-y2 + y6)) - 
     &          ddilog((-y6 + yp)/(y1 - y6)) + ddilog((-y6 + yp)/(y2 - y6)))))/2.)
     &          +
     &          (Piasqbar_one*(-(A**2*Piasqbar_one*(c6 + Piasqbar_one)) + 2*mass*(A*(2*A**2
     &          - Piasqbar_one)*Piasqbar_one*dLog(2*mass*(-A + mass)*Piasqbar_one) + A**4*mass*(-3 
     &          + alpha_if)*(-1 + alpha_if)*dLog(c61) + 
     &          Piasqbar_one*(-2*A**2 + Piasqbar_one)*(mass*dLog(-8*A**2*mass**2*Piasqbar_one 
     &          + 4*A**4*(c62)) 
     &          - mass*dLog(8*A**2*mass**2*(-Piasqbar_one + A**2*(-1 + alpha_if))) + 
     &          A*dLog(2*mass*(A*c6 + mass*Piasqbar_one - 2*A**2*mass*(-1 + alpha_if)))))))
     &          /(4.*A**4*mass**2*Sqrt(lambda_tr(Piasq_one,mk**2,mass**2)))
     &          +
     &          (B1*(-4*DLog((1d0+y8)/(1d0-y8))/2d0*dLog(2d0) 
     &           - 2*dLog(-2 + (2*A)/mass)*dLog(-1 + y7) 
     &         + dLog(4d0)*dLog(1 + y7) - dLog((2*mass)/(A + 2*mass - A*y7))**2 + 
     &         4*DLOG((2d0-(2*mass)/A)/((2*mass)/A))/2d0*(dLog(-1 + y7) - dLog(1 - y8)) + 2*dLog(-1 + A/mass)*dLog(1 - y8) + 
     &         2*dLog((2*mass)/A)*(dLog(-((A*(-1 + y7))/(-2*mass + A*(-1 + y7))))
     &         - dLog((A*(-1 + y8))/(-2*mass + A*(-1 + y8)))) - dLog(2/(1 + y8))**2 + 
     &         2*dLog(2 - (2*mass)/A)*(-dLog((A*(-1 + y7))/(A - 2*mass + A*y7)) 
     &         + dLog((A*(-1 + y8))/(A - 2*mass + A*y8))) + 
     &         2*dLog(1 + (2*mass)/A - yp)*(-dLog((A*(y7 - yp))/(A + 2*mass - A*y7))
     &         + dLog((A*(y8 - yp))/(-2*mass + A*(-1 + y8)))) + Log((1 + yp)/(1 + y8))**2 + 
     &         2*(dLog((A*(y7 - yp))/(A - 2*mass + A*y7)) - Log((A*(y8 - yp))/(A - 2*mass + A*y8)))
     &         *dLog(1 - (2*mass)/A + yp) + dLog((A + 2*mass - A*yp)/(A + 2*mass - A*y7))**2 + 
     &         2*(dLog(y7 - yp)*dLog(-1 + yp) - dLog((y7 - yp)/(-1 + y7))*dLog(-1 + yp) 
     &         + dLog((y8 - yp)/(-1 + y8))*dLog(-1 + yp) - Log(y7 - yp)*Log(1 + yp) + 
     &         dLog((y7 - yp)/(1 + y7))*Log(1 + yp) + DLOG((1d0+(2*mass)/A - yp)/(1d0-(2*mass)/A + yp))*
     &         (dLog(y7 - yp) - dLog(-y8 + yp)) - dLog(-1 + yp)*dLog(-y8 + yp) + 
     &         dLog(1 + yp)*dLog(-y8 + yp) - dLog(1 + yp)*dLog((-y8 + yp)/(1 + y8)) + 
     &         (dLog(y7 - yp) - dLog(-y8 + yp))*dLog(-(((1 + yp)*(A - 2*mass + A*yp))
     &         /((-2*mass + A*(-1 + yp))*(-1 + yp))))) - 2*ddilog(2/(1 + y7)) - 
     &         2*ddilog((A + 2*mass - A*y7)/(2.*mass)) 
     &         - 2*ddilog((2*(A - mass))/(A - 2*mass + A*y7)) - 2*ddilog((1 + y8)/2.) - 
     &         2*ddilog((2*mass)/(A + 2*mass - A*y8)) + 2*ddilog((2*(A - mass))
     &         /(A - 2*mass + A*y8)) - 2*ddilog((-1 + yp)/(-1 + y7)) + 
     &         2*ddilog((-1 + yp)/(-1 + y8)) + 2*ddilog((1 + y8)/(1 + yp)) +
     &         2*ddilog((1 + yp)/(1 + y7)) + 2*ddilog((A + 2*mass - A*y7)/(A + 2*mass - A*yp)) + 
     &         2*ddilog((A + 2*mass - A*yp)/(A + 2*mass - A*y8)) + 2*ddilog((A - 2*mass + A*yp)
     &         /(A - 2*mass + A*y7)) - 2*ddilog((A - 2*mass + A*yp)/(A - 2*mass + A*y8))))/2.
     &          +
     &          2d0*(-1 + alpha_if - Log(alpha_if))
     &            )


         endif 

        elseif(mk.eq.0d0) then
          Pff=(1d0+x**2)/(1d0-x)
          fi(2)=gsq/(8.*pi**2)*(pi**2/6d0-1d0)
          fi(3)=gsq/(8.*pi**2)*(Pff*(dlog(sik/mass**2)-1d0)-2d0/(1d0-x)*
     &          dlog(2d0-x)+(1d0+x)*dlog(1d0-x)+1d0-x)
          fi(4)=gsq/(8.*pi**2)*(Pff*(dlog(sikzone/mass**2)-1d0)-2d0/(1d0-x)*
     &          dlog(2d0-x)+(1d0+x)*dlog(1d0-x)+1d0-x)
         if(pdf) then
          fi(3)=fi(3)-gsq/(8.*pi**2)*( Pff*(dlog(muf**2/mass**2)+2d0*dlog(1d0-x)+1d0) -Cff)
          fi(4)=fi(4)-gsq/(8.*pi**2)*( Pff*(dlog(muf**2/mass**2)+2d0*dlog(1d0-x)+1d0) -Cff)
         endif
c adding the alpha dependent term
         if(zp.gt.alpha_if) then
            fi(2)=fi(2)-gsq/(8.*pi**2)*(-Log(zp)**2 + Log(1 + zp)**2 + Log(alpha_if)**2
     &            + (3*Log(alpha_if/zp))/2. - Log(1 + alpha_if)**2 + 
     &            2*ddilog(zp/(1 + zp)) - 2*ddilog(alpha_if/(1 + alpha_if)) )
            fi(3)=fi(3)-(gsq*(-((1+x)*DLog(zp/alpha_if)) + 
     &           (2*DLog(((1+alpha_if-x)*zp)/(alpha_if*(1-x+zp))))/
     &           (1 - x)))/(8.*pi**2)
            fi(4)=fi(3)
         endif

        endif



       else
         fi(2)=(gsq*(2*L_one**2-pi**2+4*rs+6*DLog(mu**2/muf**2) + 
     &        2*DLog(1 + musq_k_one)*(2*L_one + DLog(1 + musq_k_one)) + 
     &        8*ddilog(1/(1 + musq_k_one))))/(32.*pi**2)
         fi(1)=-(gsq*(-((-1 + x**2)*(L - DLog(mu**2/muf**2)))+ 
     &        (-1+x**2)*DLog(1 - x) - 2*DLog(2 - x) + 
     &        (-1+x)*(-1 + x +(1+x)*DLog((-1+x)/(-1-musq_k+x)))))/
     &        (8.*pi**2*(-1 + x))
         fi(3)=(gsq*(L-DLog(mu**2/muf**2)-2*DLog(1 - x)))/
     &        (4.*pi**2*(-1 + x))
         fi(4)=(gsq*(L_one-DLog(mu**2/muf**2)-2*DLog(1 - x)))/
     &        (4.*pi**2*(-1 + x))
         if (mk.ne.0d0) then
            fi(5)=gsq/(4d0*pi**2)*Dlog((2d0-x)/(2d0-x+musq_k))/(1d0-x)
            fi(6)=gsq/(4d0*pi**2)*Dlog(1d0/(1d0+musq_k_one))/(1d0-x)
         endif
c adding the alpha dependent term
         if(zp.gt.alpha_if) then
            fi(1)=fi(1)-(gsq*(-((1+x)*DLog(zp/alpha_if)) + 
     &           (2*DLog(((1+alpha_if-x)*zp)/(alpha_if*(1-x+zp))))/
     &           (1 - x)))/(8.*pi**2)
         endif

        endif

c fermion-photon splitting (i:fermion, ij~: photon)
      elseif(id1.eq.1.and.id2.eq.0 ) then
       if(rscheme.eq.'MREG') then
        if(pdf) then
         fi(1)=gsq/(8.*pi**2)*(dlog(sikzone*(1d0-x)/muf**2/x**2)*Pgf-2d0*(1d0-x)/x
     &          +Pgf*(2d0*dlog(x)+1d0) -Cgf+Pgf*dlog((1d0-x)/(1d0-x+musq_k)))
        else
         fi(1)=gsq/(8.*pi**2)*(dlog(sikzone*(1d0-x)/mass**2/x**2)*Pgf-2d0*(1d0-x)/x
     &           +Pgf*dlog((1d0-x)/(1d0-x+musq_k)))
        endif

       else    
        fi(1)=(gsq*(x**2 - L*(2 + (-2 + x)*x) + 
     &        (2 + (-2 + x)*x)*DLog(mu**2/muf**2) + 
     &        (2 + (-2 + x)*x)*DLog(1 - x) + 
     &        (2 + (-2 + x)*x)*DLog((-1 + x)/(-1 - musq_k + x))))/
     &        (4.*pi**2*x)
       endif
c adding the alpha dependent term
       if(zp.gt.alpha_if) then
        if(musq_k.gt.0d0) then
         fi(1)=fi(1)-(gsq*((2*musq_k*DLog((1 - zp)/(1 - alpha_if)))/x 
     &             + ((1 + (1 - x)**2)*DLog(zp/alpha_if))/x))/(8.*pi**2)
        elseif(musq_k.eq.0d0) then
        fi(1)=fi(1)-(gsq*(2-2*x+x**2)*DLog(zp/alpha_if))/(8.*pi**2*x)
        endif
       endif

c photon-fermion splitting (i:photon, ij~: fermion)
      elseif(id1.eq.0.and.id2.eq.1) then

       if(rscheme.eq.'MREG') then
        if(photonpdf) then
         fi(1)=gsq/(8.*pi**2)*(-Pfg*dlog(muf**2/sikzone/(1d0-x)*(1d0+musq_k/(1d0-x)))
     &         +2d0*x*(1d0-x) -Cfg)
        else
         fi(1)=gsq/(8.*pi**2)*(-Pfg*dlog(mass**2/sikzone/(1d0-x)*(1d0+musq_k/(1d0-x)))
     &         +2d0*x*(1d0-x) )
        endif

       else
        fi(1)=(gsq*(-L + 2*(1 + L)*x - 2*(1 + L)*x**2 + 
     &       (1 + 2*(-1 + x)*x)*DLog(mu**2/muf**2) + 
     &       (1 + 2*(-1 + x)*x)*DLog(1 - x) + 
     &       (1 + 2*(-1 + x)*x)*DLog((-1 + x)/(-1 - musq_k + x))))
     &        /(16.*pi**2)
       endif
c adding the alpha dependent term
       if(zp.gt.alpha_if) then
        fi(1)=fi(1)-(gsq*((1-x)**2+x**2)*DLog(zp/alpha_if))/(8.*pi**2)
       endif

      endif

      end



      SUBROUTINE finitefiqed(mi,sik,sikzone,x,id,id1,qint1,qint2,leg,sab,sabt,part,fi)
c  calculates the finite terms when emitter is in
c  final state and spectator is in initial state
      implicit none

      include "coupl.inc"
      include 'dipole.inc'
c Arguments
      REAL*8 mi,qsq,sik,x,fi(9),sikzone,zi
      INTEGER id,id1,qint1,qint2,leg,part
c Global
      REAL*8 ddilog,lambda_tr
      external ddilog,lambda_tr
c Local
      INTEGER i,qint(2),nc
      REAl*8 rhoi,rhok,rho,musq_i,gsq,gammaE,pi,qem(2)
      REAL*8 L,L_one,rs,musq_i_one,mass
      PARAMETER (gammaE=0.57721566d0)
      PARAMETER (pi=3.1415926535897932385d0)
      REAL*8 massb,Piasq,Piasqbar,lambda_ia,gamma,Pfg,Pgf,Riax0,Pff
      REAL*8 rhoia,Ria,z1,z2,b1,b2,b3,b4,b5,b6,b0,x0,sab,sabt
      REAL*8 Piasq_one, Piasqbar_one, lambda_ia_one,Ria_one,z1_one,z2_one
      REAL*8 A1,y1,y2,y3,y4,y5,y6,yp,A

      musq_i=mi**2/sik
      musq_i_one=mi**2/sikzone
      L_one=DLog(mu**2/sikzone)
      gsq = Dreal(gal(1))**2
c  extract the electromagnetic charge out of integers
      qint(1)=qint1
      qint(2)=qint2
      do i=1,2
       if(qint(i).eq.23) then
         qem(i)=2d0/3d0
       elseif(qint(i).eq.-23) then
         qem(i)=-2d0/3d0
       elseif(qint(i).eq.13) then
         qem(i)=1d0/3d0
       elseif(qint(i).eq.-13) then
         qem(i)=-1d0/3d0
       elseif(qint(i).eq.1) then
         qem(i)=1d0
       elseif(qint(i).eq.-1) then
         qem(i)=-1d0
       endif
      enddo

c determine the color factor for final-final a -> f fbar splitting
      if(abs(qem(1)).eq.1d0.or.id.eq.1) then
        nc=1d0
      else
        nc=3d0
      endif

      mass=massint(leg)

      gsq=nc*qem(1)*qem(2)*gsq

      L=DLog(mu**2/sik)

      do i=1,9
       fi(i)=0d0
      enddo


C  massive fermion (ij~: fermion)
      if(mi .gt. 0d0 .and. id.eq.1) then

      if(rscheme.eq.'MREG') then
c        x0=-piasqbar/(2*mass*(mass-sqrt(piasqbar)))+1d-3
        Piasqbar=-sik
        Piasqbar_one=-sikzone
        Piasq=Piasqbar+mi**2+mass**2
        Piasq_one=Piasqbar_one+mi**2+mass**2
        lambda_ia=lambda_tr(Piasq,mass**2,mi**2)
        lambda_ia_one=lambda_tr(Piasq_one,mass**2,mi**2)
c         print*, 'hier', lambda_ia_one, Piasqbar_one**2-4.*mi**2*mass**2
        Ria=sqrt((piasqbar+2*mass**2*x)**2-4d0*mass**2*piasq*x**2)/sqrt(lambda_ia)
        Ria_one=sqrt((piasqbar_one+2*mass**2*x)**2-4d0*mass**2*piasq_one*x**2)/sqrt(lambda_ia_one)
        Riax0=sqrt((piasqbar+2*mass**2*x0)**2-4d0*mass**2*piasqbar*x0**2)/sqrt(lambda_ia)
        rhoia=sqrt(lambda_tr(sabt,massint(1)**2,massint(2)**2)/
     &        lambda_tr(sab,massint(1)**2,massint(2)**2))
        z1=(piasqbar*(piasqbar-x*(piasqbar+2*mi**2)) - sqrt(piasqbar**2*(1d0-x)**2
     &     -4d0*mi**2*mass_a**2*x**2)*sqrt(lambda_ia)*Ria)/(2d0*piasqbar*(piasqbar
     &     -x*(piasq-mass**2)))
        z1_one=(piasqbar_one*(piasqbar_one-x*(piasqbar_one+2*mi**2)) - sqrt(piasqbar_one**2*(1d0-x)**2
     &     -4d0*mi**2*mass_a**2*x**2)*sqrt(lambda_ia_one)*Ria_one)/(2d0*piasqbar_one*(piasqbar_one
     &     -x*(piasq_one-mass**2)))
        z2=(piasqbar*(piasqbar-x*(piasqbar+2*mi**2)) + sqrt(piasqbar**2*(1d0-x)**2
     &     -4d0*mi**2*mass_a**2*x**2)*sqrt(lambda_ia)*Ria)/(2d0*piasqbar*(piasqbar
     &     -x*(piasq-mass**2)))
        z2_one=(piasqbar_one*(piasqbar_one-x*(piasqbar_one+2*mi**2)) + sqrt(piasqbar_one**2*(1d0-x)**2
     &     -4d0*mi**2*mass_a**2*x**2)*sqrt(lambda_ia_one)*Ria_one)/(2d0*piasqbar_one*(piasqbar_one
     &     -x*(piasq_one-mass**2)))
        b0=-4d0*mass**2*piasqbar_one/(lambda_ia_one*(1d0+abs(piasqbar_one)/sqrt(lambda_ia_one))**2
     &     +4d0*mass**2*(piasqbar_one+mi**2))
        b1=(2d0*mi**2-piasqbar_one-sqrt(lambda_ia_one))/(2d0*mi**2-piasqbar_one+sqrt(lambda_ia_one))
        b2=(2d0*mi**2+piasqbar_one+sqrt(lambda_ia_one))/(-piasqbar_one+sqrt(lambda_ia_one))*b0
        if(-piasqbar_one-sqrt(lambda_ia_one).ne.0d0) then
         b3=(2d0*mi**2+piasqbar_one-sqrt(lambda_ia_one))/(-piasqbar_one-sqrt(lambda_ia_one))*b0
        else
         b3=(2d0*mi**2+piasqbar_one-sqrt(lambda_ia_one))/(sikzone*(2.*musq_i_one*mass**2/sikzone))*b0
        endif
        b4= 2d0*(mi**2-mass**2-sqrt(lambda_ia_one))/(2d0*mass**2-piasqbar_one-sqrt(lambda_ia_one))*b0
c        b4= 2d0*(mi**2-mass**2-sqrt(lambda_ia_one))/(2d0*mass**2+sikzone*(2.*musq_i_one*mass**2/sikzone))*b0
        b5= 2d0*(mi**2-mass**2+sqrt(lambda_ia_one))/(2d0*mass**2-piasqbar_one+sqrt(lambda_ia_one))*b0
        A=Sqrt(-piasqbar_one-mi**2)
        A1=2d0*piasqbar_one/sqrt(4d0*A**2*mass**2+4d0*mass**2*piasqbar_one+piasqbar_one**2)
        y1=(sqrt(4d0*A**2*mass**2+4d0*mass**2*piasqbar_one+piasqbar_one**2)+2d0*A*mass)/(-piasqbar_one)
        y2=(sqrt(4d0*A**2*mass**2+4d0*mass**2*piasqbar_one+piasqbar_one**2)-2d0*A*mass)/(piasq_one)
        y3=(2d0*mass-A)/A
        y4=(3d0*mass+A)/(A-mass)
        y5=(2d0*mass+A)/A
        y6=(3d0*mass-A)/(A+mass)
        yp=-2d0*mass*A/piasqbar_one*(1d0-alpha_fi)-sqrt((piasqbar_one
     &      +2*mass**2*(1d0-alpha_fi))**2-4d0*mass**2*piasq_one*(1d0-alpha_fi)**2)/piasqbar_one

        if(x.gt.1d0-2.*mi*mass_a/abs(Piasqbar)) then
          return
        endif
       
         fi(2)=gsq/(8.*pi**2)*(2d0*dlog(mass_a*mi/abs(piasqbar_one))
     &         +piasqbar_one**2/(2d0*(piasqbar_one+mi**2)**2)*dlog(-piasqbar_one
     &         /mi**2)- piasqbar_one**2/(2d0*(piasqbar_one+mi**2)*
     &         (piasqbar_one))+2d0+piasqbar_one/sqrt(lambda_ia_one)*
     &         (2*dlog(b1)*dlog(b0*sqrt(lambda_ia_one*(mass**2+mi**2-piasqbar_one))
     &         /(mass_a*mass**2))-0.5d0*dlog(b1)**2+pi**2/3.+2d0*(-ddilog(b1)
     &         +ddilog(b2)-ddilog(b3)+ddilog(b4)-ddilog(b5))))
         fi(1)=0d0
         fi(3)=rhoia/x*gsq/(8.*pi**2)*(-piasqbar/sqrt(lambda_ia)/Ria/(1d0-x)*
     &         (2d0*dlog((2d0-x-z1)/(2d0-x-z2))+(z1-z2)*(1d0+(z1+z2)/2
     &         -2d0*mi**2*x/(piasqbar*(1d0-x)))))
         fi(4)=gsq/(8.*pi**2)*(-piasqbar_one/sqrt(lambda_ia_one)/Ria_one/(1d0-x)*
     &         (2d0*dlog((2d0-x-z1_one)/(2d0-x-z2_one))+(z1_one-z2_one)*(1d0+(z1_one+z2_one)/2
     &         -2d0*mi**2*x/(piasqbar_one*(1d0-x)))))
c         print*, -piasqbar_one-sqrt(lambda_ia_one), sikzone*(1.-sqrt((1.-4.*musq_i_one*mass**2/sikzone))),sikzone*(2.*musq_i_one*mass**2/sikzone)

c adding the alpha dependent term
          fi(2)=fi(2)-gsq/(8.*pi**2)*((A1*(-2*DLog(-(((-1 + y1)*(y2 - y4))/((-1 + y2)
     &          *(y1 - y4))))*DLog(-1 + y4) + DLog((-1 + y4)/(-y1 + y4))**2 
     &          + 2*DLog(-(((-1 + y1)*(y2 - y5))/((-1 + y2)*(y1 - y5))))*DLog(-1 + y5)
     &          - DLog((-1 + y5)/(-y1 + y5))**2 + 
     &          2*DLog(-1 + (2*A)/(A + mass))*DLog(((-1 + y2)*(y1 - yp))
     &          /((-1 + y1)*(y2 - yp))) + 2*DLog(((-y2 + y4)*(y1 - yp))
     &          /((y1 - y4)*(y2 - yp)))*DLog(y4 - yp) - 2*DLog(((-y2 + y5)
     &          *(y1 - yp))/((y1 - y5)*(y2 - yp)))*DLog(y5 - yp) - 
     &          DLog((-y4 + yp)/(y1 - y4))**2 + DLog((-y5 + yp)/(y1 - y5))**2 
     &          + 2*ddilog((-y1 + y4)/(-1 + y4)) + 2*ddilog((-1 + y4)/(-y2 + y4))
     &          - 2*ddilog((-y1 + y5)/(-1 + y5)) - 2*ddilog((-1 + y5)/(-y2 + y5)) - 
     &          2*ddilog((-y1 + y4)/(y4 - yp)) + 2*ddilog((-y1 + y5)/(y5 - yp)) 
     &          - 2*(DLog(1 - y3)*DLog(((-1 + y1)*(y2 - y3))/((-1 + y2)*(y1 - y3))) 
     &          + DLog(((y1 - y3)*(y2 - yp))/((y2 - y3)*(y1 - yp)))*DLog(-y3 + yp) + 
     &          ddilog((-1 + y3)/(-y1 + y3)) - ddilog((-1 + y3)/(-y2 + y3)) 
     &          - ddilog((-y3 + yp)/(y1 - y3)) + ddilog((-y3 + yp)/(y2 - y3))) 
     &          - 2*ddilog((-y4 + yp)/(y2 - y4)) + 2*ddilog((-y5 + yp)/(y2 - y5)) + 
     &          2*(DLog(1 - y6)*DLog(((-1 + y1)*(y2 - y6))/((-1 + y2)*(y1 - y6))) 
     &          + DLog(((y1 - y6)*(y2 - yp))/((y2 - y6)*(y1 - yp)))*DLog(-y6 + yp)
     &          + ddilog((-1 + y6)/(-y1 + y6)) - ddilog((-1 + y6)/(-y2 + y6)) - 
     &          ddilog((-y6 + yp)/(y1 - y6)) + ddilog((-y6 + yp)/(y2 - y6)))))/2.
     &          -2*DLog(A**2) + (A**2*piasqbar_one*(A**2 + piasqbar_one)*(-1d0 + alpha_fi)
     &          + (-piasqbar_one + A**2*(-1 + alpha_fi))*(piasqbar_one**2*(DLog(-piasqbar_one)
     &          - DLog(-piasqbar_one + A**2*(-1 + alpha_fi))) + 4*A**4*Log(A**2*alpha_fi)))
     &          /(-2*A**4*piasqbar_one + 2*A**6*(-1 + alpha_fi)) )


      else

       rs=0d0
         fi(2)=(gsq*((1d0+DLOG(musq_i_one/(musq_i_one+1d0)))*L_one
     &        -2d0*ddilog(-musq_i_one)-pi**2/3d0+2d0+
     &        DLOG(musq_i_one)**2/2d0+DLOG(1+musq_i_one)**2/2d0-
     &        2d0*DLOG(musq_i_one)*DLOG(1+musq_i_one)+DLOG(musq_i_one)))
     &        /(8.*pi**2)
c adding the alpha dependent term
            fi(2)=fi(2)+gsq/(8.*pi**2)*2.*DLog(alpha_fi)*
     &           (DLog((1d0+musq_i_one)/(musq_i_one))-1d0)
       if(x.gt.(1d0-alpha_fi)) then
             fi(1)=gsq/(8.*pi**2)*((1.-x)/(2.*(1-x+musq_i)**2)+2./
     &            (1-x)*DLog(((2.-x+musq_i)*musq_i_one)
     &            /((1.+musq_i_one)*(1-x+musq_i))))
             fi(3)=gsq/(8.*pi**2)*2./(1-x)*
     &            (DLOG((1+musq_i_one)/(musq_i_one))-1d0)
             fi(4)=fi(3)
       endif
      endif

C  massless fermion (ij~: fermion)
      elseif(mi.eq.0d0 .and. id.eq.1) then
       if(rscheme.eq.'MREG') then
         fi(2)=gsq/(8.*pi**2)*(-pi**2/2d0+3d0/2d0)
c adding the alpha dependent term
         if (alpha_fi.ne.1d0) then
            fi(2)=fi(2)+gsq/(8.*pi**2)*(-3./2.*Dlog(alpha_fi)-
     &           Dlog(alpha_fi)**2-2.*ddilog(-alpha_fi))
         endif
         if(x.gt.(1d0-alpha_fi)) then
            fi(3)=gsq/(8.*pi**2)/(1.-x)*(2d0*dlog((2d0-x)/(1d0-x))-3./2.)
            fi(4)=fi(3)
         endif

         if(NCS)then
          Pff=(1d0+zi**2)/(1d0-zi)
          fi(5)=gsq/(8.*pi**2)*(Pff*(dlog(-Piasq*zi/mi**2)-1d0) 
     &          -2d0*dlog(2d0-zi)/(1d0-zi)+(1d0+zi)*dlog(1d0-zi) ) 
          fi(6)=gsq/(8.*pi**2)*(Pff*(dlog(-Piasq_one*zi/mi**2)-1d0) 
     &          -2d0*dlog(2d0-zi)/(1d0-zi)+(1d0+zi)*dlog(1d0-zi))
          fi(7)=gsq/(8.*pi**2)*(1d0/(1d0-x)*(2d0/(2d0-x-zi)-1d0-zi))
          fi(8)=fi(7)
          if(photonfrag.and.part.eq.1) then
           fi(5)=fi(5)+gsq/(8.*pi**2)*Pff*(dlog(mi**2/muf**2*(1d0-zi)**2)+1d0)
           fi(6)=fi(5)
          endif
         endif

       else
        if(scheme .eq. 'HV') then
          rs=0d0
        elseif(scheme .eq. 'DR') then
          rs=-1d0/2d0
        endif
         fi(2)=(gsq*(6*(7+L_one*(3+L_one))-7*pi**2+12*rs))/(96.*pi**2)
c adding the alpha dependent term
         if (alpha_fi.ne.1d0) then
            fi(2)=fi(2)+gsq/(8.*pi**2)*(-3./2.*Dlog(alpha_fi)-
     &           Dlog(alpha_fi)**2)
         endif
         if(x.gt.(1d0-alpha_fi)) then
            fi(1)=gsq/(4d0*pi**2*(1-x))*DLog(2d0-x)
            fi(3)=-(gsq*(3/(1-x)+(4*DLog(1-x))/(1-x)))/(16.*pi**2)
            fi(4)=fi(3)
         endif

         if(NCS)then
          Pff=(1d0+zi**2)/(1d0-zi)
          fi(5)=gsq/(8.*pi**2)*(Pff*dlog(zi)+(zi**2+3d0)/(1d0-zi)*dlog(1d0-zi)
     &          -2d0*dlog(2d0-zi)/(1d0-zi)+1d0-zi)
          fi(6)=fi(5)
          fi(7)=gsq/(8.*pi**2)*(1d0/(1d0-x)*(2d0/(2d0-x-zi)-1d0-zi))
          fi(8)=fi(7)
           if(photonfrag.and.part.eq.1) then
            fi(5)=fi(5)+gsq/(8.*pi**2)*dlog(mu/muf)*Pff
            fi(6)=fi(5)
           endif
         endif

        endif



C  photon (ij~: photon)
      elseif(id.eq.0) then
        if(scheme .eq. 'HV') then
          rs=0d0
        elseif(scheme .eq. 'DR') then
          rs=-1d0/6d0
        endif
C  1. a-> QQ (i:massive fermion )
        if(id1.eq.1.and.mi.gt.0d0) then
         fi(7)=-(gsq*(5*Sqrt(1-4*musq_i_one)+4*musq_i_one*Sqrt(1-4*musq_i_one)+DLog(64d0) + 
     &            3*DLog(musq_i_one) - 6*DLog(1+Sqrt(1-4*musq_i_one))))/(36.*pi**2)
c adding the alpha dependent term
         fi(7)=fi(7)+(gsq*((5 - 16*musq_i_one**2 - 
     &            5*Sqrt(((1 - 4*musq_i_one)*(alpha_fi - 4*musq_i_one))/alpha_fi) - 
     &            4*musq_i_one*(4 + Sqrt(((1 - 4*musq_i_one)*(alpha_fi - 4*musq_i_one))/
     &           alpha_fi**3)))/Sqrt(1 - 4*musq_i_one)-6*DLog(1+Sqrt(1-4*musq_i_one)) + 
     &           6*DLog(Sqrt(alpha_fi)+Sqrt(alpha_fi-4*musq_i_one))))/(36.*pi**2)
       if(x.lt.1d0-4d0*musq_i)then
        if(x.gt.(1d0-alpha_fi)) then
         fi(8)=(gsq*Sqrt(1 + (4*musq_i)/(-1 + x))*(1 + 2*musq_i - x))
     &             /(12.*pi**2*(-1 + x)**2)
         fi(9)=fi(8)
        endif
       endif

c  2. a-> qq (i:massless fermion)
        elseif(id1.eq.1.and.mi.eq.0d0) then
          if(rscheme.eq.'MREG') then
            fi(2)=gsq/(8.*pi**2)*(-10d0/9d0-2d0/3d0*dlog(mass**2/sikzone))
          else
            fi(2)=-(gsq*(10 + 6*L_one))/(72.*pi**2)
          endif
c adding the alpha dependent term
            if(alpha_fi.ne.1d0) then
               fi(2)=fi(2)+gsq/(4d0*pi**2)*DLog(alpha_fi)/3.
            endif
            if(x.gt.(1d0-alpha_fi)) then
               fi(3)=gsq/(12.*pi**2*(1 - x))
               fi(4)=fi(3)
            endif
          endif


       endif

      end



      SUBROUTINE finiteffqed(mi,mk,sik,zi,id,id1,qint1,qint2,leg,part,fi)
c  calculates the finite terms when both emitter
c  and spectator are in the final state
      implicit none

      include "coupl.inc"
      include 'dipole.inc'
c Arguments
      REAL*8 mi,mk,sik,fi(9),zi
      INTEGER id,id1,qint1,qint2,leg,part
c Global
      REAL*8 ddilog,lambda_tr
      external ddilog,lambda_tr
c Local
      REAl*8 rhoi,rhok,rho,musq_i,musq_k,gsq,gammaE,pi,qem(2)
      REAL*8 xp,yp,vik,muk,Qik,Qik2,L,yl,rs,mui
      REAL*8 a,b,c,d,xa,c1,c2,v,lambda_ik,rhoj
      REAL*8 a1,a2,a3,pijsqbar,lambda_ij,mass
      REAL*8 y1,y2,y3,y4,y5,xl,ck,ci,Pff,eta,sigma,xi
      REAL*8 y2_ncs, c_ncs, c1_ncs,c2_ncs,c3_ncs
      integer i,qint(2),nc
      PARAMETER (gammaE=0.57721566d0)
      PARAMETER (pi=3.1415926535897932385d0)
      complex*16 test

      musq_i=mi**2/sik
c      print*, musq_i
      mui=Sqrt(musq_i)
      musq_k=mk**2/sik
      muk=Sqrt(musq_k)
      Qik2=sik-mi**2-mk**2
      Qik=Sqrt(Qik2)
      L=DLog(mu**2/sik)
      gsq = Dreal(gal(1))**2
      mass=massint(leg)
c  extract the electromagnetic charge out of integers
      qint(1)=qint1
      qint(2)=qint2
      do i=1,2
       if(qint(i).eq.23) then
         qem(i)=2d0/3d0
       elseif(qint(i).eq.-23) then
         qem(i)=-2d0/3d0
       elseif(qint(i).eq.13) then
         qem(i)=1d0/3d0
       elseif(qint(i).eq.-13) then
         qem(i)=-1d0/3d0
       elseif(qint(i).eq.1) then
         qem(i)=1d0
       elseif(qint(i).eq.-1) then
         qem(i)=-1d0
       endif
      enddo

c determine the color factor for final-final a -> f fbar splitting
      if(abs(qem(1)).eq.1d0.or.id.eq.1) then
        nc=1d0
      else
        nc=3d0
      endif

      gsq=nc*qem(1)*qem(2)*gsq


      vik=Sqrt(lambda_tr(1d0,musq_i,musq_k))/(1d0-musq_i-musq_k)
      lambda_ik=lambda_tr(1d0,musq_i,musq_k)
      rho=Sqrt((1d0-vik)/(1d0+vik))
      xp=1d0-2d0*muk*(1d0-muk)/(1d0-musq_i-musq_k)

C  ij~ massive and massive spectator
      if(mi.gt.0d0 .and. mk .gt. 0d0.and.id.eq.1) then
      rs=0d0
      rhoj=Sqrt((1d0-vik+2d0*musq_i/(1d0-musq_i-musq_k))/
     &     (1d0+vik+2d0*musq_i/(1d0-musq_i-musq_k)))
      rhoi=Sqrt((1d0-vik+2d0*musq_k/(1d0-musq_i-musq_k))/
     &     (1d0+vik+2d0*musq_k/(1d0-musq_i-musq_k)))
       c=2.*musq_i/(1-musq_i-musq_k)
       c1=2.*musq_k/(1-musq_i-musq_k)
       c2=4.*musq_k/(1-musq_i-musq_k)**2
       v=Sqrt(lambda_tr(1d0,musq_i,musq_k))/(1d0-musq_i-musq_k)
       xp=1d0-(2d0*muk*(1d0-muk))*(1d0-musq_i-musq_k)
      if(rscheme.eq.'MREG') then
       pijsqbar=sik-mi**2-mk**2-mass_a**2
       lambda_ij=lambda_tr(sik,mi**2,mk**2)
       a1=(pijsqbar+2.*mi**2-sqrt(lambda_ij))/(pijsqbar+2.*mi**2+sqrt(lambda_ij))
       a2=(pijsqbar-sqrt(lambda_ij))/(pijsqbar+sqrt(lambda_ij))
       a3= mi/(sqrt(sik)-mk)
       fi(2)=gsq/(8.*pi**2)*(Dlog(mass_a**2*a3**3/mi**2)-2.*Dlog(1-a3**2)
     &       +a3**2/2d0+3./2.+pijsqbar/sqrt(lambda_ij)*(Dlog(a1)
     &       *Dlog(mass_a**2*mk**2/lambda_ij/a2)+2d0*ddilog(a1)
     &       +4d0*ddilog(-sqrt(a2/a1))-4d0*ddilog(-sqrt(a1*a2))+
     &       0.5d0*Dlog(a1)**2-pi**2/3d0))

      else
       fi(2)=gsq/(8.*pi**2)*(2d0 *(-((Log(rhoj)*Log((2*xp**2*Sqrt(lambda_ik))
     &       /(c*rhoj*v)))/v) + ddilog((-2*v)/(rhoj**2*(1 + c + v)))/v + 
     &       (-Pi**2 + 6*Log((-c2 + (2 + c + c1)*(1 + c1 - v))/(1 + c - v))
     &       *Log(2*v) + 6*Log((2*v)/rhoj)*Log(1 + c + v) + 
     &       6*(Log(1 + c - v)*Log(1/(2.*v)) + Log(1 + c1 - v)
     &       *Log(1/(2.*(1 + c - v)*v)) + Log(1 + c - v)*Log((1 + c1 + v)/rhoj)) + 
     &       6*Log(rhoi*rhoj)*Log(-c2 + (2 + c + c1)*(1 + c1 - xp)) 
     &       + 6*Log(1/rhoi)*Log(1 + c1 - xp) - 3*Log((1 + c1 - xp)/(1 + c1 - v))**2 + 
     &       6*Log(rhoj)*Log(1 + c + xp) + 3*Log((1 + c + xp)/(1 + c - v))**2 - 
     &       3*Log(1/rhoj)*Log(((v - xp)**2*(v + xp)**2)
     &       /(4.*rhoj*v**2*(1 + c1 - xp)**2)) + 
     &       3*Log((c2 + (2 + c + c1)*(-1 - c1 + xp))/((1 + c - v)
     &       *(-1 - c1 + v)))**2 + 6*ddilog(rhoi) - 6*ddilog(rhoj) - 
     &       6*ddilog(rhoi*rhoj) + 6*ddilog(((1 + c - v)*(1 + c1 - v))
     &       /(-c2 + (2 + c + c1)*(1 + c1 - xp))) + 
     &       6*ddilog((-c2 + (2 + c + c1)*(1 + c1 - xp))/((1 + c + v)
     &       *(1 + c1 + v))) - 6*ddilog((1 + c1 - v)/(1 + c1 - xp)) - 
     &        6*ddilog((1 + c1 - xp)/(1 + c1 + v)) + 6*ddilog((1 + c - v)
     &       /(1 + c + xp)) + 6*ddilog((1 + c + xp)/(1 + c + v)))/(6.*v) )
     &       +Log(mui - mui*muk) + (3 + mui**2/(-1 + muk)**2 - 4*Log(1 - mui**2 
     &       - 2*muk + muk**2))/2. +((v + 2*Log(rhoj))*L)/v )
       endif

C adding alpha dependent terms
      if(alpha_ff.ne.1d0) then
       yp=1d0-2*muk*(1d0-muk)/(1d0-mui**2-muk**2)
       ci=mui**2/(1d0-mui**2-muk**2)
       ck=muk/(1d0-mui**2-muk**2) !Note, different definition of ck, compared to ci.
       xl= 1 - (2*(1 - muk)*muk)/(1 - mui**2 - muk**2) 
     &     - (1 - (2*(1 - muk)*muk)/(1 - mui**2 - muk**2))*alpha_ff + 
     &     Sqrt((1 - (2*(1 - muk)*muk)/(1 - mui**2 - muk**2) 
     &     - (1 - (2*(1 - muk)*muk)/(1 - mui**2 - muk**2))*alpha_ff)*(1 + 4*ck*muk 
     &     + (2*(1 - muk)*muk)/(1 - mui**2 - muk**2)
     &     - (1 - (2*(1 - muk)*muk)/(1 - mui**2 - muk**2))*alpha_ff))
       y1=1 + 2*ck*(-1 + muk) - Sqrt(1 + 4*ck*(muk + ck*(-1 + muk**2)))
       y2=1 + 2*ck*(-1 + muk) + Sqrt(1 + 4*ck*(muk + ck*(-1 + muk**2)))
       y3=-2*ck*(1 + ck*(-1 + mui**2 + muk**2))
       y4=-2*ck
       y5=-2*ck - 2/(-1 + mui**2 + muk**2)
       A1=-2d0/Sqrt(1 + 4*ck*(muk + ck*(-1 + muk**2)))
       A2=2d0/Sqrt(1 + 4*ck*(muk + ck*(-1 + muk**2)))
       A3=2d0/(1 - mui**2 - muk**2)
       fi(2)=fi(2)-gsq/(8.*pi**2)*((2*DLog(A3)*(A1*DLog(-y1) + A2*DLog(y2)) 
     &       + 2*A1*DLog(xl - y1)*DLog(xl - y3) + 2*A2*DLog(-xl + y2)*DLog(xl - y3)
     &       - 2*A1*DLog(xl - y3)*DLog((xl - y1)/(y1 - y3)) 
     &       + A1*DLog((xl - y3)/(y1 - y3))**2 - 2*A2*DLog(xl - y3)
     &       *DLog((-xl + y2)/(y2 - y3)) + 2*A2*DLog(y2/(y2 - y3))*DLog(-y3)
     &       - A1*DLog(-(y3/(y1 - y3)))**2 + 2*A1*DLog(-y3)*Log(y1/(-y1 + y3))
     &       - 2*A1*DLog(xl - y1)*DLog(xl - y4) - 2*A2*DLog(-xl + y2)*DLog(xl - y4) + 
     &       2*A1*DLog(xl - y4)*DLog((xl - y1)/(y1 - y4)) - A1*DLog((xl - y4)/(y1 - y4))**2
     &       + 2*A2*DLog(xl - y4)*DLog((-xl + y2)/(y2 - y4))
     &       - 2*A2*DLog(y2/(y2 - y4))*DLog(-y4) + A1*DLog(-(y4/(y1 - y4)))**2 - 
     &       2*A1*DLog(-y4)*DLog(y1/(-y1 + y4)) - 2*A1*DLog(xl - y1)
     &       *DLog((A3*(-xl + y3))/((xl - y4)*(xl - y5))) - 2*A2*DLog(-xl + y2)
     &       *DLog((A3*(-xl + y3))/((xl - y4)*(xl - y5))) - A2*DLog((xl - y5)/(y2 - y5))**2 - 
     &       2*A1*DLog(y1/(y1 - y5))*DLog(y5) + A2*DLog(-(y5/(y2 - y5)))**2 
     &       - 2*A2*DLog(-xl + y2)*DLog(-xl + y5) + 2*A2*DLog((xl - y2)/(y2 - y5))
     &       *DLog(-xl + y5) - 2*A1*DLog(-xl + y5)*DLog(-y1 + y5) 
     &        - 2*A2*DLog(y5)*DLog(y2/(-y2 + y5)))/2.
     &       + A1*(-ddilog(1 - y1/y3) + ddilog((-y1 + y3)/(-xl + y3)) + ddilog(1 - y1/y4) 
     &       - ddilog((-y1 + y4)/(-xl + y4)) + ddilog((xl - y5)/(y1 - y5))
     &       - ddilog(-(y5/(y1 - y5)))) +  A2*(-ddilog((xl - y3)/(y2 - y3))
     &       + ddilog(-(y3/(y2 - y3))) + ddilog((xl - y4)/(y2 - y4)) - ddilog(-(y4/(y2 - y4)))
     &       + ddilog(1 - y2/y5) - ddilog((-y2 + y5)/(-xl + y5)))
     &       -
     &       ((-(ci/(ci + yp)) + ci/(ci + yp*alpha_ff) - Log(ci + yp) 
     &       - 4*Log(alpha_ff) + Log(ci + yp*alpha_ff))/2.)
     &       )


      endif


C ij~ massive and massless spectator
      elseif(mi.gt.0d0.and.mk.eq.0d0.and.id.eq.1) then
        rs=0d0
      if(rscheme.eq.'MREG') then
         fi(2)=gsq/(8.*pi**2)*(-pi**2/3d0+3d0/2d0)+gsq/(8.*pi**2)*((1d0 + Dlog(mi**2/sik))*DLog(mass_a**2/sik)
     &       -0.5d0*DLog(mi**2/sik)**2+0.5d0*Dlog(mi**2/sik))
      else

        fi(2)=gsq/(8.*pi**2)*(2 *( Pi**2/6.+3*dLog(mui)**2
     &      -4*DLog(mui)*Log(1 - mui**2) + 2*ddilog(1 - mui**(-2)) - ddilog(mui**2) ) + 
     &       (3 + mui**2 + 2*DLog(mui) - 4*DLog(1 - mui**2))/2.+L + 2*Log(mui)*L)
      endif 

C adding alpha dependent terms
       if(alpha_ff.ne.1d0) then
        fi(2)=fi(2)+gsq/(8.*pi**2)*(2d0*((-Pi**2/6. - 2*Log(mui)*Log(alpha_ff) 
     &        + DLog((-1 + mui**2)*(-1 + alpha_ff))*Log(mui**2 
     &        + alpha_ff - mui**2*alpha_ff) - 
     &        ddilog(1 - mui**(-2)) + ddilog(alpha_ff - alpha_ff/mui**2)
     &        + ddilog(mui**2 + alpha_ff - mui**2*alpha_ff)))
     &        +(((-1 + alpha_ff)*mui**2*(-1 + mui**2))/(-alpha_ff + (-1 + alpha_ff)*mui**2) 
     &        + 4*dLog(alpha_ff) - dLog(alpha_ff + mui**2 - alpha_ff*mui**2))/2.
     &        )
        endif


C ij~ massless fermion and massive spectator
      elseif(mi.eq.0d0 .and. mk.gt.0d0 .and. id.eq.1) then
        if(scheme .eq. 'HV') then
          rs=0d0
        elseif(scheme .eq. 'DR') then
          rs=-1d0/2d0
        endif

      if(rscheme.eq.'MREG') then
       pijsqbar=sik-mass**2-mk**2-mass_a**2
       lambda_ij=lambda_tr(sik,mass**2,mk**2)
       a1=(pijsqbar+2.*mass**2-sqrt(lambda_ij))/(pijsqbar+2.*mass**2+sqrt(lambda_ij))
       a2=(pijsqbar-sqrt(lambda_ij))/(pijsqbar+sqrt(lambda_ij))
       a3= mass/(sqrt(pijsqbar)-mk)

       fi(2)=gsq/(8.*pi**2)*(Dlog(mass_a**2*a3**3/mass**2)-2.*Dlog(1-a3**2)
     &       +a3**2/2d0+3./2.+pijsqbar/sqrt(lambda_ij)*(Dlog(a1)
     &       *Dlog(mass_a**2*mk**2/lambda_ij/a2)+2d0*ddilog(a1)
     &       +4d0*ddilog(-sqrt(a2/a1))-4d0*ddilog(-sqrt(a1*a2))+
     &       0.5d0*Dlog(a1)**2-pi**2/3d0))

         if(NCS) then
          sigma=dsqrt(1d0+4d0*mk**2/pijsqbar*zi*(1d0-zi))
          xi=mk**2/(2d0*pijsqbar*zi*(1d0-zi))
          if(zi.lt.0.5d0) then
           eta=(1d0-(xi+1d0+dsqrt(xi*(xi+2)))**(-1d0))*zi
          elseif(zi.gt.0.5d0) then
           eta=(1d0-(xi+1d0+dsqrt(xi*(xi+2)))**(-1d0))*(1d0-zi)
          endif
          Pff=(1d0+zi**2)/(1d0-zi)
          fi(3)=fi(3)+gsq/(8.*pi**2)*(-Pff*dlog(mi**2/zi/pijsqbar*(1d0-eta))
     &          +(1d0+zi)*dlog(1d0-zi)-2d0*zi/(1d0-zi)+(1d0+zi)*
     &          dlog(1d0+mk/pijsqbar/eta)-2d0/(1d0-zi)/sigma*(dlog(1d0+
     &          pijsqbar*eta*(1d0-zi*eta)/mk**2/(1d0-zi))-2d0*dlog(1d0-2*zi*eta/sigma)
     &          +sigma*dlog(mk**2/pijsqbar/eta*(1d0-zi))))
          if(photonfrag.and.part.eq.1) then
           fi(3)=fi(3)+gsq/(8.*pi**2)*Pff*(dlog(mi**2/muf**2*(1d0-zi)**2)+1d0)
          endif
         endif

      else
        fi(2)=gsq/(8.*pi**2)*(2d0*( (-7*Pi**2)/12.+DLog(-1 + 1/muk)*Log(1-muk)
     &        + dLog(1 + muk)*(4*dLog(1 - muk) - 3*dLog(muk) + 2*dLog(1 + muk)) + 
     &         2*ddilog(1/(1 + muk)) ) + 3.5d0 - 3*DLog(1 - muk)
     &        +(-Pi**2 + 6*L*(3 - 4*dLog(1 - muk**2) + L))/12.+rs)

         if(NCS) then
           xi=mk**2/(2d0*pijsqbar*zi*(1d0-zi))
           y2_ncs=1d0/(xi+1+dsqrt(xi*(xi+2d0)))
           c_ncs=1d0-musq_k
           c1_ncs=1d0-musq_k**2
           c2_ncs=dsqrt(c_ncs**2*(y2_ncs**2+1d0)-2d0*c1_ncs*y2_ncs)
           c3_ncs=c1_ncs-c_ncs*(c_ncs*y2_ncs+c2_ncs)
           Pff=(1d0+zi**2)/(1d0-zi)
c           fi(3)=fi(3)+gsq/(8.*pi**2)*(1 - zi - ((1 + zi**2)
c     &                *dLog(-(c_ncs*y2_ncs*(-1 + zi)*zi)))/(-1 + zi) + 
c     &                (-(zi*(-1 + zi**2)*(c1_ncs*(-dLog(-c_ncs**2 + c1_ncs) 
c     &                 + dLog(c3_ncs)) + c_ncs*(c2_ncs + c_ncs*(-1 + Log(4d0))
c     &                 + 2*c_ncs*dLog(c_ncs*(-c_ncs**2 + c1_ncs)) 
c     &                 - c_ncs*(2*Log(c3_ncs) + dLog(2*c_ncs*(c_ncs + c2_ncs)
c     &                 - 2*c1_ncs*y2_ncs))))) + 2*c_ncs**2
c     &                 *dLog(1 + (y2_ncs*zi)/(1 - zi)))/(c_ncs**2*(-1 + zi)*zi))
           fi(3)=fi(3)+gsq/(8.*pi**2)*(1 - zi + Pff*dLog(-(c_ncs*y2_ncs*(-1 + zi)*zi)) 
     &          + (dLog(2d0) - 2*(-1 + zi**2)*dLog(c_ncs) + dLog(-c_ncs**2 + c1_ncs) 
     &          - dLog(c3_ncs) + zi**2*(-dLog(-8*c_ncs**2 + 8*c1_ncs) + dLog(c3_ncs) 
     &          + dLog(4*c_ncs*(c_ncs + c2_ncs) - 4*c1_ncs*y2_ncs))
     &          - dLog(c_ncs*(c_ncs + c2_ncs) - c1_ncs*y2_ncs) - 2*dLog(1 - zi) + 
     &          2*dLog(1 + (-1 + y2_ncs)*zi))/(-1 + zi))
           if(photonfrag.and.part.eq.1) then
            fi(3)=fi(3)+gsq/(8.*pi**2)*dlog(mu/muf)*Pff
           endif
         endif

      endif

C adding alpha dependent terms
         if (alpha_ff.ne.1d0) then
            xa=(alpha_ff - muk*alpha_ff + Sqrt((-1 + alpha_ff)*(-(1 + muk)**2 
     &          + (-1 + muk)**2*alpha_ff)))/(1 + muk)
            fi(2)=fi(2)+gsq/(8.*pi**2)*(2d0*((DLog(2 - 2*muk)*DLog(2/(1 + muk)) 
     &            - DLog((2*muk)/(1 + muk))*DLog(-1 + 2/(1 + muk)) - 
     &             DLog(1 - xa)*DLog(-4/((-1 + muk**2)*(1 + xa)**2))
     &             - DLog(1 + xa)*DLog(-((-1 + muk**2)*(1 + xa))) - 
     &             2*ddilog(1/(1 + muk)) + 2*ddilog((1 + xa)/2.)))
     &             +(3d0*Log(alpha_ff))/2.
     &             )
         endif


C ij~ massless fermion and massless spectator
      elseif(mi.eq.0d0 .and. mk.eq.0d0 .and. id.eq.1) then
        if(scheme .eq. 'HV') then
          rs=0d0
        elseif(scheme .eq. 'DR') then
          rs=-1d0/2d0
        endif
          Pff=(1d0+zi**2)/(1d0-zi)
        if(rscheme.eq.'MREG') then
         fi(2)=gsq/(8.*pi**2)*(-pi**2/3d0+3d0/2d0)+gsq/(8.*pi**2)*((1d0 + Dlog(mass**2/sik))*DLog(mass_a**2/sik)
     &       -0.5d0*DLog(mass**2/sik)**2+0.5d0*Dlog(mass**2/sik))
         if(NCS) then
          fi(3)=fi(3)+gsq/(8.*pi**2)*(Pff*(dlog(sik*zi/mass**2)-1d0)
     &            +(1d0+zi)*dlog(1d0-zi))
          if(photonfrag.and.part.eq.1) then
           fi(3)=fi(3)+gsq/(8.*pi**2)*Pff*(dlog(mi**2/muf**2*(1d0-zi)**2)+1d0)
          endif
         endif
        else
         fi(2)=gsq/(8.*pi**2)*(7./2.-pi**2/6.+3./2.*L+1./2.*L**2-pi**2/12.+rs)
         if(NCS) then
          fi(3)=fi(3)+gsq/(8.*pi**2)*(Pff*dlog(zi)+(zi**2+3d0)/(1d0-zi)*dlog(1d0-zi)
     &                +(1d0-zi))
          if(photonfrag.and.part.eq.1) then
           fi(3)=fi(3)+gsq/(8.*pi**2)*dlog(mu/muf)*Pff
          endif
         endif
        endif

C adding alpha dependent terms
         if (alpha_ff.ne.1d0) then
            fi(2)=fi(2)+gsq/(8.*pi**2)*(-3./2.*Dlog(alpha_ff)-Dlog(alpha_ff)**2
     &            -2.*ddilog(1.-alpha_ff))
         endif





C ij~ photon and massive spectator
      elseif(mk.gt.0d0.and. id.eq.0) then
        if(scheme .eq. 'HV') then
          rs=0d0
        elseif(scheme .eq. 'DR') then
          rs=-1d0/6d0
        endif
C  1. a->FF splitting
       if(id1.eq.1.and.mi.gt.0d0) then
          fi(2)=(gsq*((2*(1 - (4*musq_i)/(-1 + muk)**2)**1.5*muk)/
     &         (1 + muk) + (8*Sqrt(1 - (4*musq_i)/(-1 + muk)**2)*
     &         (1 + mui - muk)*(-1 + mui + muk))/(3.*(-1 + muk)**2)
     &         - 2*DLog(mui) + 2*DLog((1 +
     &         Sqrt(1 - (4*musq_i)/(-1 + muk)**2))/2.)+2*DLog(1 - muk)- 
     &         (2*musq_k*((8*musq_i*Sqrt(1 - (4*musq_i)/(-1 + muk)**2))/
     &         (-1 + musq_k)-DLog((1-Sqrt(1-(4*musq_i)/(-1 + muk)**2))/
     &         (1 + Sqrt(1 - (4*musq_i)/(-1 + muk)**2))) + 
     &         ((-1 + 4*musq_i + musq_k)/(-1 + musq_k))**1.5*
     &         DLog((-Sqrt(1 - (4*musq_i)/(-1 + muk)**2) +
     &         Sqrt((-1 + 4*musq_i + musq_k)/(-1 + musq_k)))/
     &         (Sqrt(1 - (4*musq_i)/(-1 + muk)**2) +
     &         Sqrt((-1 + 4*musq_i + musq_k)/(-1 + musq_k))))))
     &         /(1 - musq_k)))/(12.*Pi**2)


C adding alpha dependent terms
          if(alpha_ff.ne.1d0) then
           a=1.-musq_k
           c=-1.+2*musq_i+musq_k
           yp=1.-2*muk*(1.-muk)/(1.-2*musq_i-musq_k)
           b=sqrt(c**2*yp**2-4.*musq_i**2)
           d=sqrt(alpha_ff*c**2*yp**2-4.*musq_i**2)
           fi(2)=fi(2)-gsq*((c*Sqrt(-c + 2*musq_i)*(4*(b - d)*musq_i**2 
     &          + c**2*yp*(-2*alpha_ff*b-d*(-2+yp) + alpha_ff**2*b*yp)+ 
     &          4*c*musq_i*(d - d*yp + b*(-1 + alpha_ff*yp))) - 
     &          2*b*d*Sqrt(-c + 2*musq_i)*
     &          (c**2 - 2*(1 + c)*musq_i + 4*musq_i**2)
     &          *DATAN((2*musq_i)/Sqrt(-4*musq_i**2 + c**2*yp**2)) + 
     &          2*b*d*Sqrt(-c + 2*musq_i)*
     &          (c**2 - 2*(1 + c)*musq_i + 4*musq_i**2)
     &          *DATAN((2*musq_i)/
     &          Sqrt(-4*musq_i**2 + alpha_ff**2*c**2*yp**2)) + 
     &          b*d*(2*a*(c + c**2 + 2*musq_i - 4*musq_i**2)*
     &          DLog(-((-c - 2*musq_i)**2.5
     &          *(1 + c - 2*musq_i)*(-1 + yp))) - 
     &          2*a*(c + c**2 + 2*musq_i - 4*musq_i**2)*
     &          DLog((-c - 2*musq_i)**2.5
     &          *(1 + c - 2*musq_i)*(1 - alpha_ff*yp)) + 
     &          2*c*Sqrt(-c + 2*musq_i)*DLog(-2*(b + c*yp)) 
     &          + 3*c**2*Sqrt(-c + 2*musq_i)*DLog(-2*(b + c*yp)) - 
     &          4*c*musq_i*Sqrt(-c + 2*musq_i)*DLog(-2*(b + c*yp)) 
     &          - 2*c*Sqrt(-c + 2*musq_i)*Log(-2*(d + alpha_ff*c*yp)) - 
     &          3*c**2*Sqrt(-c + 2*musq_i)*DLog(-2*(d + alpha_ff*c*yp)) 
     &          + 4*c*musq_i*Sqrt(-c + 2*musq_i)*
     &          DLog(-2*(d + alpha_ff*c*yp)) - 
     &          2*a*c*DLog(-4*musq_i**2 - b*Sqrt(c**2 - 4*musq_i**2) +
     &          c**2*yp) - 2*a*c**2*DLog(-4*musq_i**2 -
     &          b*Sqrt(c**2 - 4*musq_i**2) + c**2*yp) - 
     &          4*a*musq_i*DLog(-4*musq_i**2 - b*Sqrt(c**2 - 4*musq_i**2)
     &          + c**2*yp) + 8*a*musq_i**2*DLog(-4*musq_i**2 -
     &          b*Sqrt(c**2 - 4*musq_i**2) + c**2*yp) + 
     &          2*a*c*DLog(-4*musq_i**2 - d*Sqrt(c**2 - 4*musq_i**2) +
     &          alpha_ff*c**2*yp) + 2*a*c**2*Log(-4*musq_i**2 -
     &          d*Sqrt(c**2 - 4*musq_i**2) + alpha_ff*c**2*yp) + 
     &          4*a*musq_i*DLog(-4*musq_i**2 - d*Sqrt(c**2 - 4*musq_i**2)
     &          + alpha_ff*c**2*yp) - 8*a*musq_i**2*DLog(-4*musq_i**2 -
     &          d*Sqrt(c**2 - 4*musq_i**2) + alpha_ff*c**2*yp)))/
     &          (3.*c*(-c+2*musq_i)**1.5*Sqrt(-4*musq_i**2 + c**2*yp**2)
     &          *Sqrt(-4*musq_i**2+alpha_ff**2*c**2*yp**2)))/
     &          (4.*Pi**2)
             endif

C  2. a->ff splitting
       elseif(id1.eq.1.and. mi.eq.0d0) then
          fi(2)=(gsq*((-1 + muk)*(-6*L*(1 + muk) - 4*(4 + muk)
     &         ) + 12*(-1 + musq_k)*DLog(1 - muk) + 
     &         12*musq_k*DLog((2*muk)/(1 + muk))))/
     &         (72.*(-1 + musq_k)*pi**2)
       if(rscheme.eq.'MREG') then
          fi(2)=(gsq*((-1 + muk)*(- 4*(4 + muk)
     &         ) + 12*(-1 + musq_k)*DLog(1 - muk) + 
     &         12*musq_k*DLog((2*muk)/(1 + muk))))/
     &         (72.*(-1 + musq_k)*pi**2)
     &         +gsq/(8.*pi**2)*(-2d0/3d0*dlog(mass**2/sik))
        endif
C     adding alpha dependent terms
          if(alpha_ff.ne.1d0) then
             fi(2)=fi(2)-(gsq*(-((-1 + alpha_ff)*(-1 + muk)**2) + 
     &            DLog(alpha_ff) + musq_k*(DLog(4d0) + DLog(1/alpha_ff)
     &            + 2*DLog(muk) - 2*DLog(1 + muk) - 
     &            2*DLog(1+alpha_ff-(2*alpha_ff)/(1+muk))))/
     &            (3.*(-1+musq_k)))/(4.*pi**2)
          endif


       endif

C ij~ photon and massless spectator
      elseif(mk.eq.0d0.and.id.eq.0) then
        if(scheme .eq. 'HV') then
          rs=0d0
        elseif(scheme .eq. 'DR') then
          rs=-1d0/6d0
        endif
C  1. a->FF splitting
       if(id1.eq.1.and.mi.gt.0d0) then
            fi(2)=(gsq*(-16*Sqrt(1-4*musq_i)+16*musq_i*Sqrt(1-4*musq_i)+
     &           12*DLog(Sqrt(musq_i))+ 
     &           12*DLog((1 + Sqrt(1 - 4*musq_i))/2.)))/(72.*pi**2)
C adding alpha dependent terms
            if (alpha_ff.ne.1d0) then
               fi(2)=fi(2)-(gsq*(Sqrt(1-4*musq_i)+
     &              Sqrt(-4*musq_i**2+alpha_ff**2*(1-2*musq_i)**2) + 
     &              (2*Sqrt(-4*musq_i**2+alpha_ff**2*(1-2*musq_i)**2))/
     &              (-alpha_ff+2*(-1+alpha_ff)*musq_i) + 
     &              (-1 + 2*musq_i)*(-2*DATAN((2*musq_i)/
     &              Sqrt(1 - 4*musq_i)) + 
     &              2*DATan((2*musq_i)/
     &              Sqrt(-4*musq_i**2 + alpha_ff**2*(1 - 2*musq_i)**2))+ 
     &              DLog(-2*(-1 + 2*musq_i + Sqrt(1 - 4*musq_i))) - 
     &              DLog(-2*(alpha_ff*(-1+2*musq_i)+
     &              Sqrt(-4*musq_i**2+alpha_ff**2*(1-2*musq_i)**2))))))/
     &              (12.*pi**2)
            endif

C  2. a->ff splitting
       elseif(id1.eq.1.and. mi.eq.0d0) then
        if(rscheme.eq.'MREG') then
          fi(2)=-(gsq*(16d0/9d0+2d0/3d0*dlog(mass**2/sik)  ))/(8.*pi**2)
        else
          fi(2)=-(gsq*(8 + 3*L))/(72.*pi**2)
        endif
          
C adding alpha dependent terms
            if (alpha_ff.ne.1d0) then
               fi(2)=fi(2)-(gsq*(-1 + alpha_ff - DLog(alpha_ff)))/
     &              (12.*pi**2)
            endif


       endif
      endif

      end



      FUNCTION DDILOG(X)
*
* imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp
*
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(0:19)
      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12)
      DATA C( 0) / 0.42996 69356 08136 97D0/
      DATA C( 1) / 0.40975 98753 30771 05D0/
      DATA C( 2) /-0.01858 84366 50145 92D0/
      DATA C( 3) / 0.00145 75108 40622 68D0/
      DATA C( 4) /-0.00014 30418 44423 40D0/
      DATA C( 5) / 0.00001 58841 55418 80D0/
      DATA C( 6) /-0.00000 19078 49593 87D0/
      DATA C( 7) / 0.00000 02419 51808 54D0/
      DATA C( 8) /-0.00000 00319 33412 74D0/
      DATA C( 9) / 0.00000 00043 45450 63D0/
      DATA C(10) /-0.00000 00006 05784 80D0/
      DATA C(11) / 0.00000 00000 86120 98D0/
      DATA C(12) /-0.00000 00000 12443 32D0/
      DATA C(13) / 0.00000 00000 01822 56D0/
      DATA C(14) /-0.00000 00000 00270 07D0/
      DATA C(15) / 0.00000 00000 00040 42D0/
      DATA C(16) /-0.00000 00000 00006 10D0/
      DATA C(17) / 0.00000 00000 00000 93D0/
      DATA C(18) /-0.00000 00000 00000 14D0/
      DATA C(19) /+0.00000 00000 00000 02D0/
      IF(X .EQ. 1) THEN
       H=PI6
      ELSEIF(X .EQ. -1) THEN
       H=-PI12
      ELSE
       T=-X
       IF(T .LE. -2) THEN
        Y=-1/(1+T)
        S=1
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2)
       ELSEIF(T .LT. -1) THEN
        Y=-1-T
        S=-1
        A=LOG(-T)
        A=-PI6+A*(A+LOG(1+1/T))
       ELSE IF(T .LE. -HF) THEN
        Y=-(1+T)/T
        S=1
        A=LOG(-T)
        A=-PI6+A*(-HF*A+LOG(1+T))
       ELSE IF(T .LT. 0) THEN
        Y=-T/(1+T)
        S=-1
        A=HF*LOG(1+T)**2
       ELSE IF(T .LE. 1) THEN
        Y=T
        S=1
        A=0
       ELSE
        Y=1/T
        S=-1
        A=PI6+HF*LOG(T)**2
       ENDIF
       H=Y+Y-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 19,0,-1
       B0=C(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=-(S*(B0-H*B2)+A)
      ENDIF
      DDILOG=H
      RETURN
      END


