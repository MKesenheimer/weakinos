      Subroutine epsiiqed(id1,id2,qint1,qint2,leg,sik,epssq,eps)
C  Calculates the 1/eps^2 and 1/eps term 
C  for initial state emitter and initial state spectator
      implicit none

      include 'coupl.inc'
      include 'dipole.inc'
C  Arguments
      integer id1,id2,qint1,qint2,leg
      real*8 sik,epssq,eps

c  Local
      REAL*8 pi,L,gsq,qem(2), mass
      PARAMETER (pi=3.1415926535897932385d0)
      INTEGER i,qint(2),nc

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

      L=dlog(mu**2/sik)
      epssq=0d0
      eps=0d0


C  fermion-fermion splitting (i:fermion, ij~: fermion)
      if(id1.eq.1 .and. id2.eq.1) then
        if(rscheme.eq.'MREG') then
c         eps=gsq/(8.*pi**2)*((1d0 + Dlog(mass**2/sik))*DLog(mass_a**2/sik)
c     &       -0.5d0*DLog(mass**2/sik)**2+0.5d0*Dlog(mass**2/sik)
         return
        else
         epssq=gsq/(8d0*pi**2)
         eps=gsq*(3d0 + 2d0*L)/(16.*pi**2)
        endif

C  fermion-photon splitting (i:fermion, ij~: photon)
      elseif(id1.eq.1.and.id2.eq.0) then
       return

C  photon-fermion splitting (i:photon, ij~: fermion)
      elseif(id1.eq.0.and.id2.eq.1) then
       return

      endif
      end


      Subroutine epsifqed(id1,id2,qint1,qint2,leg,sik,mk,epssq,eps)
C  Calculates the 1/eps^2 and 1/eps term 
C  for initial state emitter and final state spectator
      implicit none

      include 'coupl.inc'
      include 'dipole.inc'
C  Arguments
      integer id1,id2,qint1,qint2,leg
      real*8 sik,mk,epssq,eps

c  Local
      REAL*8 pi,L,gsq,musq_k,qem(2)
      PARAMETER (pi=3.1415926535897932385d0)
      INTEGER i,qint(2),nc

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

      L=dlog(mu**2/sik)
      musq_k=mk**2/sik
      epssq=0d0
      eps=0d0

C  fermion-fermion splitting (i:fermion, ij~: fermion)
      if(id1.eq.1.and.id2.eq.1) then
       if(rscheme.eq.'MREG') then
c        eps=gsq/(8.*pi**2)*((1d0 + Dlog(mass**2/sik))*DLog(mass_a**2/sik)
c     &       -0.5d0*DLog(mass**2/sik)**2+0.5d0*Dlog(mass**2/sik)
        return
       else
        epssq=gsq/(8d0*pi**2)
        eps=(gsq*(3+2*L+2*DLog(1+musq_k)))/(16.*pi**2)
       endif


C  fermion-photon splitting (i:fermion, ij~: photon)
      elseif(id1.eq.1.and.id2.eq.0) then
       return

C  photon-fermion splitting (i:photon, ij~: fermion)
      elseif(id1.eq.0.and.id2.eq.1) then
       return

      endif
      end


      Subroutine epsfiqed(id1,id2,qint1,qint2,sik,mi,zi,part,epssq,eps)
C  Calculates the 1/eps^2 and 1/eps term 
C  for final state emitter and initia; state spectator
      implicit none

      include 'coupl.inc'
      include 'dipole.inc'
C  Arguments
      integer id1,id2,qint1,qint2,part
      real*8 sik,mi,epssq,eps,zi

c  Local
      REAL*8 pi,L,gsq,musq_i,qem(2)
      PARAMETER (pi=3.1415926535897932385d0)
      INTEGER i,qint(2),nc

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
      if(abs(qem(1)).eq.1d0.or.id1.eq.1) then
        nc=1d0
      else
        nc=3d0
      endif

      gsq=nc*qem(1)*qem(2)*gsq

      L=dlog(mu**2/sik)
      musq_i=mi**2/sik
      epssq=0d0
      eps=0d0


      if(rscheme.eq.'MREG') then
       return
      endif

C  massive fermion (ij~:fermion)
      if(id1.eq.1.and.mi.gt.0d0) then
        eps=gsq*(1+DLog(musq_i/(1+musq_i)))/(8.*pi**2)


C  massless fermion( ij~:fermion)
      elseif(id1.eq.1.and.mi.eq.0d0) then
        epssq= gsq/(8d0*pi**2)
        eps=gsq*(3 + 2*L)/(16.*pi**2)
        if(ncs) then
         eps=eps -gsq/(8d0*pi**2)*(1.+zi**2)/(1.-zi)
         if(photonfrag.and.part.eq.1) then
         eps=eps +gsq/(8d0*pi**2)*(1.+zi**2)/(1.-zi)
         endif
        endif


C  photon splitting (ij~: photon)
C  1. a-> QQ (i:massive fermion )
      elseif(id1.eq.0.and. mi.gt.0d0 ) then
        return
c  2. a-> qq (i:massless fermion)
      elseif(id1.eq.0.and.id2.eq.1.and.mi.eq.0d0) then
        eps=-gsq/(24.*pi**2)

      endif
      end


      Subroutine epsffqed(mi,mk,sik,zi,id,id1,qint1,qint2,part,epssq,eps)
C  Calculates the 1/eps^2 and 1/eps term 
C  for final state emitter and initia; state spectator
      implicit none

      include 'coupl.inc'
      include 'dipole.inc'
c Arguments
      REAL*8 mi,mk,sik,epssq,eps,zi
      INTEGER id,id1,qint1,qint2,part
c Global
      REAL*8 ddilog,lambda_tr
      external ddilog,lambda_tr
c  Local
      REAL*8 pi,L,gsq,vik,rho,musq_i,musq_k,qem(2),Qik,Qik2
      PARAMETER (pi=3.1415926535897932385d0)
      INTEGER i,qint(2),NC

      gsq = DReal(gal(1))**2
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

      Qik2=sik-mi**2-mk**2
      Qik=Sqrt(Qik2)
      L=dlog(mu**2/sik)
      musq_i=mi**2/sik
      musq_k=mk**2/sik
      vik=Sqrt(lambda_tr(1d0,musq_i,musq_k))/(1d0-musq_i-musq_k)
      rho=Sqrt((1d0-vik)/(1d0+vik))

      epssq=0d0
      eps=0d0

      if(rscheme.eq.'MREG') then
       return
      endif

C  ij~ massive and massive spectator
      if(mi.gt.0d0.and.mk.gt.0d0) then
        eps=(gsq*(1d0 + DLog(rho)/vik))/(8.*pi**2)

C  ij~ massive and massless spectator
      elseif(mi.gt.0d0.and.mk.eq.0d0) then
        eps=gsq/(8.*pi**2)*(dlog(musq_i)+1d0)

C  ij~ massless fermion and massive spectator
      elseif(id.eq.1.and.mi.eq.0d0.and.mk.gt.0d0) then
        epssq=gsq/(8d0*pi**2)
        eps=gsq/(8.*pi**2)*(3d0/2d0 - 2*Log(1 - musq_k) + L)
        if(ncs) then
         eps=eps -gsq/(8d0*pi**2)*(1.+zi**2)/(1.-zi)
         if(photonfrag.and.part.eq.1) then
         eps=eps +gsq/(8d0*pi**2)*(1.+zi**2)/(1.-zi)
         endif
        endif

C  ij~ massless fermion and massless spectator
      elseif(id.eq.1.and.mi.eq.0d0.and.mk.eq.0d0) then
        epssq=gsq/(8d0*pi**2)
        eps=(gsq*(3 + 2*L))/(16.*pi**2)
        if(ncs) then
         eps=eps -gsq/(8d0*pi**2)*(1.+zi**2)/(1.-zi)
         if(photonfrag.and.part.eq.1) then
         eps=eps +gsq/(8d0*pi**2)*(1.+zi**2)/(1.-zi)
         endif
        endif


C  ij~ photon and massive spectator
      elseif(id.eq.0.and.mk.gt.0d0) then
C  1. a->FF splitting
       if(id1.eq.1.and.mi.gt.0d0) then
         return
       endif

C  2. a->ff aplitting
       elseif(id1.eq.1.and.mi.eq.0d0) then
         eps=-gsq/(24.*pi**2)

       endif

C  ij~ photon and massless spectator
      if(id.eq.0.and.mk.eq.0d0) then
C  1. a->FF cplitting
       if(id1.eq.1.and.mi.gt.0d0) then
         return
       endif

C  2. a->ff aplitting
       elseif(id1.eq.1.and.mi.eq.0d0) then
         eps=-gsq/(24.*pi**2)

      endif
      end
