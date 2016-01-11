c This file has been generated automatically by FeynRules

c Version: 1.4.0   Date: 15. 7. 2009,    16:23

c Internal parameters definition :

c Parameters that should not be recomputed event by event.
      if(readlha) then
      aEW = 1/aEWM1
      sw2 = sw**2
      cw = sqrt(1. - 1.*sw2)
      ee = 3.5449077018110318*sqrt(aEW)
      gw = ee/sw
      g1 = ee/cw
      G = 3.5449077018110318*sqrt(aS)
      v = (2*MW*sw)/ee
      lam = MH**2/(2.*v**2)
      muH = sqrt(lam*v**2)
      ytau = (sqrt(2.)*ymtau)/v
      yc = (sqrt(2.)*ymc)/v
      yt = (sqrt(2.)*ymt)/v
      yb = (sqrt(2.)*ymb)/v
      c12 = sqrt(1. - 1.*s12**2)
      c23 = sqrt(1. - 1.*s23**2)
      c13 = sqrt(1. - 1.*s13**2)
      CKM11 = c12*c13
      CKM12 = c13*s12
      CKM13 = s13
      CKM21 = -(c23*s12) - c12*s13*s23
      CKM22 = c12*c23 - s12*s13*s23
      CKM23 = c13*s23
      CKM31 = -(c12*c23*s13) + s12*s23
      CKM32 = -(c23*s12*s13) - c12*s23
      CKM33 = c13*c23
      YukD11 = 0
      YukD12 = 0
      YukD13 = 0
      YukD21 = 0
      YukD22 = 0
      YukD23 = 0
      YukD31 = 0
      YukD32 = 0
      YukD33 = (sqrt(2.)*ymb)/v
      gmass = 3.5449077018110318*sqrt(amass)
      mn = (0.006332573977646111*Log(LR**2))/R
      mQ2 = (g1**2/16. + 3*gmass**2 + (27*gw**2)/16.)*mn
      mu2 = (g1**2 + 3*gmass**2)*mn
      md2 = (g1**2/4. + 3*gmass**2)*mn
      mL2 = ((9*g1**2)/16. + (27*gw**2)/16.)*mn
      me2 = (9*g1**2*mn)/4.
      mQt = (-3*mn*ymt**2)/(2.*v**2)
      mt0 = (-3*mn*ymt**2)/v**2
      cH2 = (mn*((3*g1**2)/4. + (3*gw**2)/2. - MH**2/v**2))/R
      MM11 = R**(-2) - (0.015039773079174007*g1**2)/R**2 + (g1**2*v
     +**2)/4. - (0.0010554289962743518*g1**2*Log(LR**2))/R**2
      MM12 = (g1*gw*v**2)/4.
      MM22 = R**(-2) - (0.0019281760357915392*gw**2)/R**2 + (gw**2*
     +v**2)/4. + (0.04749430483234583*gw**2*Log(LR**2))/R**2
      Delta = sqrt(4.*MM12**2 + (MM11 - 1.*MM22)**2)
      Sqrt2 = sqrt(2.)
      SqrtPi = sqrt(3.141592653589793)
      Pgw2 = gw**2
      Psw2 = sw**2
      Pcw2 = cw**2
      CONJCKM23 = conjg(CKM23)
      CONJCKM33 = conjg(CKM33)
      CONJCKM13 = conjg(CKM13)
      CONJCKM21 = conjg(CKM21)
      CONJCKM31 = conjg(CKM31)
      CONJCKM11 = conjg(CKM11)
      CONJCKM22 = conjg(CKM22)
      CONJCKM32 = conjg(CKM32)
      CONJCKM12 = conjg(CKM12)
      endif
c Parameters that should be recomputed at an event by even basis.
      GGT1 = sqrt(G**2)
      G1G1T2 = sqrt(1.5)*sqrt(G**2)


c Definition of the EW coupling used in the write out of aqed
      gal(1) = ee
      gal(2) = ee


c Definition of DUM symbols
      DUM0 = 0
      DUM1 = 1


