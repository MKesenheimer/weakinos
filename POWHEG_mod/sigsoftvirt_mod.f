c MK: copied and modified version of sigsoftvirt.f, revision 3154
c modified on the basis of disquark/sigsoftvirt.f
c changes marked with "! CH, MK:" or "! MK:"
c bigger changes over a whole section are marked with !===...

      subroutine btildevirt(resvirt)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_br.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      include 'pwhg_pdf.h'
      integer pdgfac ! MK: added
      parameter(pdgfac=1000000)
      real * 8 resvirt(flst_nborn)
      real * 8 c(-6:6),gamma(-6:6),gammap(-6:6),gammatilde(-6:6)
      integer j,jb,fl1,fl2,fl,leg,legi,legj,jres,kres
      real * 8 Q,I,s,etot,e,eij,virt_arr(flst_nborn),tot,kij,arglog,
     1     loglog,pborn(0:3,nlegborn)
      real * 8 pdfb1(-pdf_nparton:pdf_nparton),
     1     pdfb2(-pdf_nparton:pdf_nparton)      
      logical ini
      data ini/.true./
      real * 8 ddilog,dotp
      external ddilog,dotp
      save c,gamma,gammap,gammatilde,ini
      real * 8 tiny
      parameter (tiny=1d-14)
      real * 8 ll,llew,rescfac
      real * 8 Intm_ep,Int0m_0,Int0m_ep,Intmm_0,Intmm_ep
      external Intm_ep,Int0m_0,Int0m_ep,Intmm_0,Intmm_ep
      real * 8 tmp,tmpew,Qem,Iem,ch_i,ch_j,ch1,ch2
      logical pwhg_isfinite
      real * 8 chargeofparticle
      external pwhg_isfinite,chargeofparticle
      integer sigma
      external sigma
      logical is_charged, is_coloured
      external is_charged, is_coloured
c     if virtual.f, will not do smartsig when this flag is active
      if (flg_novirtual) then
         resvirt = 0d0
         return
      endif
c from 2.100 of FNO2007
      if(ini) then
         do j=-6,6
            if(j.eq.0) then
               c(j)=ca
               gamma(j)=(11*ca-4*tf*st_nlight)/6
               gammap(j)=(67d0/9-2*pi**2/3)*ca-23d0/9*tf*st_nlight
               gammatilde(j)=ca/6d0
            else
               c(j)=cf
               gamma(j)=3d0/2*cf
               gammap(j)=(13d0/2-2*pi**2/3)*cf
               gammatilde(j)=cf/2
            endif
         enddo
         ini=.false.
      endif
      if(.not.flg_minlo) then
c get pdfs
         call pdfcall(1,kn_xb1,pdfb1)
         call pdfcall(2,kn_xb2,pdfb2)         
      endif

      tot=0

      call sigvirtual(virt_arr)

c Strong Q contribution
      Q = 0
c EW Q contribution
      Qem = 0

      do jres=1,flst_nreson
         kres=flst_reslist(jres)
         if(kres.ne.0) then
            call boost2reson(kn_cmpborn(:,kres),nlegborn,
     1           kn_cmpborn,pborn)
            etot=pborn(0,kres)
         else
            pborn=kn_cmpborn
            etot=2*kn_cmpborn(0,1)
         endif
         s=etot**2
         if(.not.flg_minlo) then
            ll = log(par_csicut**2*s/st_muren2)
            if(flg_with_em) llew = log(par_csicut**2*s/em_muren2)
            rescfac = 1
         endif
         do jb=1,flst_nborn
c Strong Q contribution
            Q = 0
c EW Q contribution
            Qem = 0

            if(flg_minlo) then
               call setlocalscales(jb,2,rescfac)
               ll = log(par_csicut**2*s/st_muren2)
               if(flg_with_em) llew = log(par_csicut**2*s/em_muren2)
c get pdfs
               call pdfcall(1,kn_xb1,pdfb1)
               call pdfcall(2,kn_xb2,pdfb2)      
            endif

            fl1=flst_born(1,jb)
            fl2=flst_born(2,jb)
c     initial-state parton contribution
            if(kres.eq.0) then
               if(is_coloured(fl1)) then
                  Q=Q-log(st_mufact2/st_muren2)*(gamma(fl1)+2*c(fl1)*
     1                 log(par_csicut))
                  if(flg_drscheme) then
                     Q=Q-gammatilde(fl1)
                  endif
               endif
               if(is_coloured(fl2)) then
                  Q=Q-log(st_mufact2/st_muren2)*
     1                 (gamma(fl2)+2*c(fl2)*log(par_csicut))
                  if(flg_drscheme) then
                     Q=Q-gammatilde(fl2)
                  endif
               endif
c EW photonic contribution
               ch1=chargeofparticle(fl1)
               ch2=chargeofparticle(fl2)
               if (flg_with_em.and.ch1.ne.0) then
                  Qem=Qem-log(st_mufact2/em_muren2)*( 
     1               (1.5d0+2*log(par_csicut))*ch1**2)
               endif
               if (flg_with_em.and.ch2.ne.0) then
                  Qem=Qem-log(st_mufact2/em_muren2)*( 
     3             + (1.5d0+2*log(par_csicut))*ch2**2)
               endif
            endif
c     loop on final-state massless partons
            do leg=3,nlegborn
               if(leg.ne.kres.and.flst_bornres(leg,jb).ne.kres) cycle
               fl=flst_born(leg,jb)
               if(is_coloured(fl)) then
                  if(kn_masses(leg).eq.0) then
                     if(flg_drscheme) then
                        Q=Q-gammatilde(fl)
                     endif
                     e=pborn(0,leg)
                     Q=Q+gammap(fl)
     1                    -log(s/st_muren2)*(gamma(fl)-2*c(fl)
     3                    *log(2*e/(par_csicut*etot)))
     3                    +2*c(fl)*(log(2*e/etot)**2-log(par_csicut)**2)
     4                    -2*gamma(fl)*log(2*e/etot)
                  else
                     ! CH: to include Sparticles here: simply extract the
                     ! flavour-index the flag is_coloured(fl) restricts 
                     ! this automatically to sq and go 
                     Q = Q - c(mod(flst_born(leg,jb),pdgfac))*
     1                    (ll-0.5*Intm_ep(pborn(0,leg)))
#ifdef DEBUG
                     print*,"in sigsoftvirt_mod.f:165"
                     print*,"uncomment to continue"
                     print*,"Q",Q
                     stop
#endif
                  endif
               endif
c WEW
               if (flg_with_em.and.chargeofparticle(fl).ne.0) then
                  if(kn_masses(leg).eq.0) then
                     write(*,*) ' not prepared to handle this!'
                     call pwhg_exit(-1)
                  else
                     Qem = Qem - 
     1                    (llew-0.5*Intm_ep(kn_cmpborn(0,leg))) 
     2                    * chargeofparticle(fl)**2
                  endif
               endif
            enddo
            Q=Q*br_born(jb)
            Qem=Qem*br_born(jb)
            Iem=0
            I=0
            do legi=1,nlegborn
               if(legi.ne.kres.and.flst_bornres(legi,jb).ne.kres)
     1              cycle
               do legj=legi+1,nlegborn
                  if(legj.ne.kres.and.flst_bornres(legj,jb).ne.kres)
     1                 cycle
c     both particles are colored
                  if(is_coloured(flst_born(legi,jb)).and.
     1                 is_coloured(flst_born(legj,jb))     ) then
c     massless-massless case
                     if (kn_masses(legi).eq.0.and.kn_masses(legj).eq.0)
     1                    then
                        kij=dotp(pborn(0,legi),pborn(0,legj))
                        eij=pborn(0,legi)*pborn(0,legj)
                        arglog = abs(1-kij/(2*eij))
                        if (arglog.lt.tiny) then
                           loglog = 0d0
                        else
                           loglog = log(arglog)*log(kij/(2*eij))
                        endif                     
                        I=I+(1d0/2*ll**2
     1                       +ll*log(kij/(2*eij))
     2                       -ddilog(kij/(2*eij))
     3                       +1d0/2*log(kij/(2*eij))**2
     4                       -loglog)*
     5                       br_bornjk(legi,legj,jb)
                     endif
c     massless-massive case
                     if (kn_masses(legi).eq.0.and.kn_masses(legj).gt.0)
     1                    then
                        I = I + (0.5*(0.5*ll**2 - pi**2/6)
     1                     +0.5*Int0m_0(pborn(0,legi),pborn(0,legj))*ll
     2                     -0.5*Int0m_ep(pborn(0,legi),pborn(0,legj)))*
     3                       br_bornjk(legi,legj,jb)
                     endif
c     massive-massless case
                     if (kn_masses(legi).gt.0.and.kn_masses(legj).eq.0)
     1                    then
                        I = I + (0.5*(0.5*ll**2 - pi**2/6)
     1                     +0.5*Int0m_0(pborn(0,legj),pborn(0,legi))*ll
     2                    -0.5*Int0m_ep(pborn(0,legj),pborn(0,legi)))*
     3                    br_bornjk(legj,legi,jb)
                     endif
c     massive-massive case
                     if (kn_masses(legi).gt.0.and.kn_masses(legj).gt.0)
     1                    then
                        I = I + (0.5*ll*
     1                    Intmm_0(pborn(0,legi),pborn(0,legj))
     2                    -0.5*Intmm_ep(pborn(0,legi),pborn(0,legj)))
     3                    *br_bornjk(legi,legj,jb)
                     
                     endif
                  endif
                  ch_i=chargeofparticle(flst_born(legi,jb))
                  ch_j=chargeofparticle(flst_born(legj,jb))

c     both particles are charged
                  if (flg_with_em.and.ch_i.ne.0.and.ch_j.ne.0) then
c     massless-massless case
                     if (kn_masses(legi).eq.0.and.kn_masses(legj).eq.0)
     1                    then
                        kij=dotp(kn_cmpborn(0,legi),kn_cmpborn(0,legj))
                        eij=kn_cmpborn(0,legi)*kn_cmpborn(0,legj)
                        arglog = abs(1-kij/(2*eij))
                        if (arglog.lt.tiny) then
                           loglog = 0d0
                        else
                           loglog = log(arglog)*log(kij/(2*eij))
                        endif                     
                        Iem=Iem+(1d0/2*llew**2
     1              +llew*log(kij/(2*eij))
     2              -ddilog(kij/(2*eij))+1d0/2*log(kij/(2*eij))**2
     3              -loglog)*
     4              br_born(jb)*ch_i*ch_j*(-1)*
     5              sigma(flst_born(legi,jb),legi)*
     6              sigma(flst_born(legj,jb),legj)
                     endif
c     massless-massive case
                     if (kn_masses(legi).eq.0.and.kn_masses(legj).gt.0)
     1                    then
                        Iem = Iem + (0.5*(0.5*llew**2 - pi**2/6)
     1          +0.5*Int0m_0(kn_cmpborn(0,legi),kn_cmpborn(0,legj))*llew
     2          -0.5*Int0m_ep(kn_cmpborn(0,legi),kn_cmpborn(0,legj)))*
     3              br_born(jb)*ch_i*ch_j*(-1)*
     4              sigma(flst_born(legi,jb),legi)*
     5              sigma(flst_born(legj,jb),legj)
                     endif
c     massive-massless case
                     if (kn_masses(legi).gt.0.and.kn_masses(legj).eq.0)
     1                    then
                        Iem = Iem + (0.5*(0.5*llew**2 - pi**2/6)
     1          +0.5*Int0m_0(kn_cmpborn(0,legj),kn_cmpborn(0,legi))*llew
     2          -0.5*Int0m_ep(kn_cmpborn(0,legj),kn_cmpborn(0,legi)))*
     3              br_born(jb)*ch_i*ch_j*(-1)*
     4              sigma(flst_born(legi,jb),legi)*
     5              sigma(flst_born(legj,jb),legj)
                     endif
c     massive-massive case
                     if (kn_masses(legi).gt.0.and.kn_masses(legj).gt.0)
     1                    then
                        Iem = Iem + (0.5*llew*
     1            Intmm_0(kn_cmpborn(0,legi),kn_cmpborn(0,legj))
     2            -0.5*Intmm_ep(kn_cmpborn(0,legi),kn_cmpborn(0,legj)))*
     3              br_born(jb)*ch_i*ch_j*(-1)*
     4              sigma(flst_born(legi,jb),legi)*
     5              sigma(flst_born(legj,jb),legj)
                     endif
                  endif
               enddo
            enddo
            I=I*2
c WEW
            Iem=Iem*2
            if(jres.eq.1) then
c     first time
               resvirt(jb)=( (Q+I+virt_arr(jb))*st_alpha
     1              + (Qem+Iem)*em_alpha )/(2*pi)
     2              *pdfb1(fl1)*pdfb2(fl2)*kn_jacborn * rescfac
            else
c     others
               resvirt(jb)=resvirt(jb)+(   (Q+I)*st_alpha
     1              + (Qem+Iem)*em_alpha )/(2*pi)
     2              *pdfb1(fl1)*pdfb2(fl2)*kn_jacborn * rescfac
            endif
            tot=tot+resvirt(jb)
         enddo
      enddo
      if (.not.pwhg_isfinite(tot)) then
         do jb=1,flst_nborn
            resvirt(jb)=0d0
         enddo
      endif
      end



c I_ep(k,m)
      function Int0m_ep(k,m)
c               /                                           [                  ]
c               |           d phi                           [   k.m       k0   ]
c Int0m_ep = -2 | d cos th  ----- log(sin th sin phi)  l0^2 [ ------- - ------ ]
c               |             pi                            [ k.l m.l   k.l l0 ]
c               /                                           [                  ]
c the range in phi is 0<phi<pi.
      implicit none
      real * 8 Int0m_ep, k(0:3),m(0:3)
      real * 8 kh(0:3),mh(0:3),khmh,b
      integer j
      real * 8 ddilog,dotp
      external ddilog,dotp
      do j=0,3
         kh(j)=k(j)/k(0)
         mh(j)=m(j)/m(0)
      enddo
      b=sqrt(1-dotp(mh,mh))
      khmh=dotp(kh,mh)
      Int0m_ep=-2*(log((1-b)/(1+b))**2/4+log(khmh/(1+b))*log(khmh/(1-b))
     # +ddilog(1-khmh/(1+b))+ddilog(1-khmh/(1-b)))
      end



c I_0(k,m)
c            /                       [                  ]
c            |           d phi       [   k.m       k0   ]
c  I_0     = | d cos th  -----  l0^2 [ ------- - ------ ]
c            |             pi        [ k.l m.l   k.l l0 ]
c            /                       [                  ]
c the range in phi is 0<phi<pi.
      function Int0m_0(k,m)
      implicit none
      real * 8 Int0m_0, k(0:3),m(0:3)
      real * 8 kh(0:3),mh(0:3)
      integer j
      real * 8 dotp
      external ddilog,dotp
      do j=0,3
         kh(j)=k(j)/k(0)
         mh(j)=m(j)/m(0)
      enddo
      Int0m_0=log( dotp(kh,mh)**2/dotp(mh,mh) )
      end



      function Intm_ep(p)
c               /                                              2
c               |           d phi                             p  
c  Intm_ep = -2 | d cos th  ----- log(sin th sin phi) l0^2  -------
c               |             pi                            (p.l)^2
c               /
c the range in phi is 0<phi<pi.
      implicit none
      real * 8 p(0:3),Intm_ep
      real * 8 beta2,beta
      integer j
      beta2=0
      do j=1,3
         beta2=beta2+(p(j)/p(0))**2
      enddo
      beta=sqrt(beta2)
      if(beta.lt.1d-4) then
         Intm_ep=2*(2+2*beta**2/3+2*beta**4/5+2*beta**6/7)
      else
         Intm_ep=2*log((1+beta)/(1-beta))/beta
      endif
      end

      function Intmm_ep(p1,p2)
c               / 
c               |           d phi                             p1.p2  
c Intmm_ep= -2  | d cos th  ----- log(sin th sin phi) k0^2 ----------
c               |             pi                           p1.k  p2.k
c               /
c p1^2>0, p2^2>0.
c The range in phi is 0<phi<pi.
      implicit none
      real * 8 p1(0:3),p2(0:3), intmm_ep
      real * 8 i,z,zm,zmob,zp,z1,z1ob,z2,a,b,x1,x2,bv1(3),bv2(3),
     1   bb1,bb2,b1b2,b2
      integer j
      real * 8 ddilog,tmp
      external ddilog
      i(z)=-1d0/2*log((z-zm)*(zp-z)/((zp+z)*(z+zm)))**2
     # -2*ddilog(2*zm/(zp-zm)*(zp-z)/(zm+z))
     # -2*ddilog(-2*zp/(zp-zm)*(zm+z)/(zp-z))
      do j=1,3
         bv1(j)=p1(j)/p1(0)
         bv2(j)=p2(j)/p2(0)
      enddo
      bb1=0
      bb2=0
      b1b2=0
      do j=1,3
         bb1=bb1+bv1(j)**2
         bb2=bb2+bv2(j)**2
         b1b2=b1b2+bv1(j)*bv2(j)
      enddo
      if(bb1.gt.bb2) then
         tmp=bb1
         bb1=bb2
         bb2=tmp
      endif
      a=bb1+bb2-2*b1b2
      x1=bb1-b1b2
      x2=bb2-b1b2
      b=bb1*bb2-b1b2**2
      zp=sqrt(a)+sqrt(a-b)
      zm=sqrt(a)-sqrt(a-b)
c zm over b (introduce to handle back-to-back case without generating 0/0)
      zmob=1/zp
      z1=sqrt(x1**2+b)-x1
      z2=sqrt(x2**2+b)+x2
      if(bb1.eq.0) then
c handle limiting case, b1=0
         b2=sqrt(bb2)
         Intmm_ep=(-0.5d0*log((1-b2)/(1+b2))**2
     1      +2*ddilog(0d0)-2*ddilog(-2*b2/(1-b2)))
     1        * (1-b1b2)/sqrt(a-b)
      else
c     z1 over b (same as before
         z1ob=1/(sqrt(x1**2+b)+x1)
         Intmm_ep=(i(z2)
     1 +1d0/2*log((z1ob-zmob)*(zp-z1)/((zp+z1)*(z1ob+zmob)))**2
     2 +2*ddilog(2*zmob/(zp-zm)*(zp-z1)/(zmob+z1ob))
     3 +2*ddilog(-2*zp/(zp-zm)*(zm+z1)/(zp-z1)) )
     4  *(1-b1b2)/sqrt(a-b)
      endif
      end


      function Intmm_0(p1,p2)
c            / 
c            |           d phi        p1.p2  
c  Intmm_0 = | d cos th  ----- k0^2 ----------
c            |             pi       p1.k  p2.k
c            /
c the range in phi is 0<phi<pi.
      implicit none
      real * 8 p1(0:3),p2(0:3), Intmm_0
      real * 8 kh(3),gh(3),k2,g2,kg,beta,ombeta2
      integer j
      do j=1,3
         kh(j)=p1(j)/p1(0)
         gh(j)=p2(j)/p2(0)
      enddo
      k2=0
      g2=0
      kg=0
      do j=1,3
         k2=k2+kh(j)**2
         g2=g2+gh(j)**2
         kg=kg+kh(j)*gh(j)
      enddo
      ombeta2=(1-k2)*(1-g2)/(1-kg)**2
      beta=sqrt(1-ombeta2)
c      Intmm_0=log((1+beta)/(1-beta))/beta
c Equivalent:
      Intmm_0=log((1+beta)**2/ombeta2)/beta
      end
