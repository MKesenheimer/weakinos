c This file has been generated automatically by FeynRules

c Version: 1.5.9   Date: 8. 7. 2009,    11:55

c Internal parameters definition :

c Parameters that should not be recomputed event by event.
      if(readlha) then
      aEW = 1/aEWM1
      sw2 = 1 - MW**2/MZ**2
      ee = 3.5449077018110318*sqrt(aEW)
      cw = sqrt(1. - 1.*sw2)
      sw = sqrt(sw2)
      gw = ee/sw
      g1 = ee/cw
      G = 3.5449077018110318*sqrt(aS)
      v = (2*MW*sw)/ee
      CKM11 = (0,1)*CKMI11 + CKMR11
      CKM12 = (0,1)*CKMI12 + CKMR12
      CKM21 = (0,1)*CKMI21 + CKMR21
      CKM22 = (0,1)*CKMI22 + CKMR22
      l6 = (0,1)*lI6 + lR6
      l7 = (0,1)*lI7 + lR7
      GU11 = (0,1)*GUI11 + GUR11
      GU12 = (0,1)*GUI12 + GUR12
      GU13 = (0,1)*GUI13 + GUR13
      GU21 = (0,1)*GUI21 + GUR21
      GU22 = (0,1)*GUI22 + GUR22
      GU23 = (0,1)*GUI23 + GUR23
      GU31 = (0,1)*GUI31 + GUR31
      GU32 = (0,1)*GUI32 + GUR32
      GU33 = (0,1)*GUI33 + GUR33
      GD11 = (0,1)*GDI11 + GDR11
      GD12 = (0,1)*GDI12 + GDR12
      GD13 = (0,1)*GDI13 + GDR13
      GD21 = (0,1)*GDI21 + GDR21
      GD22 = (0,1)*GDI22 + GDR22
      GD23 = (0,1)*GDI23 + GDR23
      GD31 = (0,1)*GDI31 + GDR31
      GD32 = (0,1)*GDI32 + GDR32
      GD33 = (0,1)*GDI33 + GDR33
      GL11 = (0,1)*GLI11 + GLR11
      GL12 = (0,1)*GLI12 + GLR12
      GL13 = (0,1)*GLI13 + GLR13
      GL21 = (0,1)*GLI21 + GLR21
      GL22 = (0,1)*GLI22 + GLR22
      GL23 = (0,1)*GLI23 + GLR23
      GL31 = (0,1)*GLI31 + GLR31
      GL32 = (0,1)*GLI32 + GLR32
      GL33 = (0,1)*GLI33 + GLR33
      mu1 = -(l1*v**2)
      mu2 = mhc**2 - (l3*v**2)/2.
      mu3 = -(l6*v**2)/2.
      Sqrt2 = sqrt(2.)
      SqrtPi = sqrt(3.141592653589793)
      Pgw2 = gw**2
      Psw2 = sw**2
      Pcw2 = cw**2
      CONJCKM21 = conjg(CKM21)
      CONJCKM11 = conjg(CKM11)
      CONJCKM22 = conjg(CKM22)
      CONJCKM12 = conjg(CKM12)
      CONJl7 = conjg(l7)
      PTH112 = TH11**2
      PTH212 = TH21**2
      PTH312 = TH31**2
      PTH113 = TH11**3
      CONJl6 = conjg(l6)
      PTH213 = TH21**3
      PTH313 = TH31**3
      PTH114 = TH11**4
      PTH214 = TH21**4
      PTH314 = TH31**4
      PTH122 = TH12**2
      PTH222 = TH22**2
      PTH322 = TH32**2
      PTH123 = TH12**3
      PTH223 = TH22**3
      PTH323 = TH32**3
      PTH124 = TH12**4
      PTH224 = TH22**4
      PTH324 = TH32**4
      PTH132 = TH13**2
      PTH232 = TH23**2
      PTH332 = TH33**2
      PTH133 = TH13**3
      PTH233 = TH23**3
      PTH333 = TH33**3
      PTH134 = TH13**4
      PTH234 = TH23**4
      PTH334 = TH33**4
      CONJGD23 = conjg(GD23)
      CONJGD33 = conjg(GD33)
      CONJGD13 = conjg(GD13)
      CONJGD21 = conjg(GD21)
      CONJGD31 = conjg(GD31)
      CONJGD11 = conjg(GD11)
      CONJGD22 = conjg(GD22)
      CONJGD32 = conjg(GD32)
      CONJGD12 = conjg(GD12)
      CONJGL11 = conjg(GL11)
      CONJGL21 = conjg(GL21)
      CONJGL31 = conjg(GL31)
      CONJGL12 = conjg(GL12)
      CONJGL22 = conjg(GL22)
      CONJGL32 = conjg(GL32)
      CONJGL13 = conjg(GL13)
      CONJGL23 = conjg(GL23)
      CONJGL33 = conjg(GL33)
      CONJGU32 = conjg(GU32)
      CONJGU12 = conjg(GU12)
      CONJGU22 = conjg(GU22)
      CONJGU33 = conjg(GU33)
      CONJGU13 = conjg(GU13)
      CONJGU23 = conjg(GU23)
      CONJGU31 = conjg(GU31)
      CONJGU11 = conjg(GU11)
      CONJGU21 = conjg(GU21)
      endif
c Parameters that should be recomputed at an event by even basis.
      GGT1 = sqrt(G**2)


c Definition of the EW coupling used in the write out of aqed
      gal(1) = ee
      gal(2) = ee


c Definition of DUM symbols
      DUM0 = 0
      DUM1 = 1


