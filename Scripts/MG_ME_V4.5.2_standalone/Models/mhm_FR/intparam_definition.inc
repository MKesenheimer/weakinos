c This file has been generated automatically by FeynRules

c Version: 1.5.8   Date: 7. 7. 2009,    13:57

c Internal parameters definition :

c Parameters that should not be recomputed event by event.
      if(readlha) then
      G = 3.5449077018110318*sqrt(aS)
      aEW = 1/aEWM1
      EE = 3.5449077018110318*sqrt(aEW)
      cM = MW/MZ
      sM = sqrt(1. - 1.*cM**2)
      RM = MW/MWP
      x = (1 + RM**2 - sqrt(1. - 6.*RM**2 + RM**4))/(2.*RM)
      t = sqrt((-1.*(-1. + cM**2)*(2.*cM**4*(x**2 + x**4) + x**2*(2
     +. + x**2 - 1.*sqrt(4. + x**4)) + cM**2*(-1.*x**6 + 2.*x**2*(
     +-4. + sqrt(4. + x**4)) + x**4*(-4. + sqrt(4. + x**4)) + 2.*(
     +-2. + sqrt(4. + x**4)))))/(cM**2*x**2*(x**2 + cM**4*(1. + x*
     +*2)**2 - 1.*cM**2*(2. + 3.*x**2 + x**4))))/sqrt(2.)
      EpL = (sqrt(2.)*x)/sqrt(2. - 1.*x**2 + sqrt(4. + x**4))
      g0 = EE*sqrt(1. + t**(-2) + x**2)
      gp = EE*sqrt(1. + t**2 + t**2*x**2)
      gt = EE*sqrt(1. + x**(-2) + 1/(t**2*x**2))
      fpi = (2*sqrt(2.)*MW)/(gt*sqrt(2. + x**2 - 1.*sqrt(4. + x**4)
     +))
      EpRnu = 0
      EpRel = 0
      EpRup = 0
      EpRdo = 0
      EpRmu = sqrt((mmu**2*(-1.*MF**2 - 1.*EpL**2*MF**2 + mmu**2))/
     +(MF**2*(-1.*EpL**2*MF**2 + mmu**2)))
      EpRta = sqrt((mta**2*(-1.*MF**2 - 1.*EpL**2*MF**2 + mta**2))/
     +(MF**2*(-1.*EpL**2*MF**2 + mta**2)))
      EpRch = sqrt((mch**2*(mch**2 - 1.*MF**2 - 1.*EpL**2*MF**2))/(
     +MF**2*(mch**2 - 1.*EpL**2*MF**2)))
      EpRst = sqrt((mst**2*(-1.*MF**2 - 1.*EpL**2*MF**2 + mst**2))/
     +(MF**2*(-1.*EpL**2*MF**2 + mst**2)))
      EpRto = sqrt((mto**2*(-1.*MF**2 - 1.*EpL**2*MF**2 + mto**2))/
     +(MF**2*(-1.*EpL**2*MF**2 + mto**2)))
      EpRbo = sqrt((mbo**2*(mbo**2 - 1.*MF**2 - 1.*EpL**2*MF**2))/(
     +MF**2*(mbo**2 - 1.*EpL**2*MF**2)))
      MGZ = MZ
      MGW = MW
      MGZP = MZP
      MGWP = MWP
      MhNu = (sqrt(1. + EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**4))*MF)
     +/sqrt(2.)
      v0g = t/sqrt(1. + t**2*(1. + x**2))
      v1g = (t*x)/sqrt(1. + t**2*(1. + x**2))
      v2g = 1/Sqrt(1 + t**2*(1 + x**2))
      v0Z = ((-1 + t**2)*x**2 + sqrt(4. + (-1. + t**2)**2*x**4))/(t
     +*sqrt(4. + (-2. + (-1. + t**2)*x**2 + sqrt(4. + (-1. + t**2)
     +**2*x**4))**2/(t**2*x**2) + ((-1. + t**2)*x**2 + sqrt(4. + (
     +-1. + t**2)**2*x**4))**2/t**2))
      v1Z = -((-2 + (-1 + t**2)*x**2 + sqrt(4. + (-1. + t**2)**2*x*
     +*4))/(sqrt(2.)*sqrt(4. + (-1. + t**2)**2*x**6 - 2.*sqrt(4. +
     + (-1. + t**2)**2*x**4) + (-1. + t**2)*x**4*(-1. + t**2 + sqr
     +t(4. + (-1. + t**2)**2*x**4)) + x**2*(4. - 1.*sqrt(4. + (-1.
     + + t**2)**2*x**4) + t**2*sqrt(4. + (-1. + t**2)**2*x**4)))))
      v2Z = -2/sqrt(4. + (-2. + (-1. + t**2)*x**2 + sqrt(4. + (-1. 
     ++ t**2)**2*x**4))**2/(t**2*x**2) + ((-1. + t**2)*x**2 + sqrt
     +(4. + (-1. + t**2)**2*x**4))**2/t**2)
      v0ZP = ((-1 + t**2)*x**2 - sqrt(4. + (-1. + t**2)**2*x**4))/(
     +t*sqrt(4. + (-1.*(-1. + t**2)*x**2 + sqrt(4. + (-1. + t**2)*
     +*2*x**4))**2/t**2 + (2. - 1.*(-1. + t**2)*x**2 + sqrt(4. + (
     +-1. + t**2)**2*x**4))**2/(t**2*x**2)))
      v1ZP = (2 - (-1 + t**2)*x**2 + sqrt(4. + (-1. + t**2)**2*x**4
     +))/(t*x*sqrt(4. + (-1.*(-1. + t**2)*x**2 + sqrt(4. + (-1. + 
     +t**2)**2*x**4))**2/t**2 + (2. - 1.*(-1. + t**2)*x**2 + sqrt(
     +4. + (-1. + t**2)**2*x**4))**2/(t**2*x**2)))
      v2ZP = -2/sqrt(4. + (-1.*(-1. + t**2)*x**2 + sqrt(4. + (-1. +
     + t**2)**2*x**4))**2/t**2 + (2. - 1.*(-1. + t**2)*x**2 + sqrt
     +(4. + (-1. + t**2)**2*x**4))**2/(t**2*x**2))
      v0W = (2 - x**2 + sqrt(4. + x**4))/(sqrt(2.)*sqrt(4. + x**4 +
     + 2.*sqrt(4. + x**4) - 1.*x**2*sqrt(4. + x**4)))
      v1W = 1/Sqrt(1 + (2 - x**2 + sqrt(4. + x**4))**2/(4.*x**2))
      v0WP = -(-2 + x**2 + sqrt(4. + x**4))/(2.*x*sqrt(1. + (0.25*(
     +-2. + x**2 + sqrt(4. + x**4))**2)/x**2))
      v1WP = 1/Sqrt(1 + (-2 + x**2 + sqrt(4. + x**4))**2/(4.*x**2))
      Npi = sqrt(v1W**2 + (-1.*v1W + v0W*x)**2)
      v0pi = (-v1W + v0W*x)/Npi
      v1pi = v1W/Npi
      v0piP = -v1pi
      v1piP = v0pi
      Npi0 = sqrt((-1.*v1Z + v0Z*x)**2 + (v1Z - 1.*t*v2Z*x)**2)
      v0pi0 = (-v1Z + v0Z*x)/Npi0
      v1pi0 = (v1Z - t*v2Z*x)/Npi0
      v0pi0P = -v1pi0
      v1pi0P = v0pi0
      vL0nu = -(1 - EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**4))/(2.*EpL
     +*sqrt(1. + (0.25*(1. - 1.*EpL**2 + sqrt(1. + 2.*EpL**2 + EpL
     +**4))**2)/EpL**2))
      vL0el = -(1 - EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**4))/(2.*EpL
     +*sqrt(1. + (0.25*(1. - 1.*EpL**2 + sqrt(1. + 2.*EpL**2 + EpL
     +**4))**2)/EpL**2))
      vL0mu = -(1 - EpL**2 + EpRmu**2 + sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRmu**2) + (1. + EpRmu**2)**2))/(2.*EpL*sqrt(1. + (0.25
     +*(1. - 1.*EpL**2 + EpRmu**2 + sqrt(EpL**4 - 2.*EpL**2*(-1. +
     + EpRmu**2) + (1. + EpRmu**2)**2))**2)/EpL**2))
      vL0ta = -(1 - EpL**2 + EpRta**2 + sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRta**2) + (1. + EpRta**2)**2))/(2.*EpL*sqrt(1. + (0.25
     +*(1. - 1.*EpL**2 + EpRta**2 + sqrt(EpL**4 - 2.*EpL**2*(-1. +
     + EpRta**2) + (1. + EpRta**2)**2))**2)/EpL**2))
      vL0up = -(1 - EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**4))/(2.*EpL
     +*sqrt(1. + (0.25*(1. - 1.*EpL**2 + sqrt(1. + 2.*EpL**2 + EpL
     +**4))**2)/EpL**2))
      vL0do = -(1 - EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**4))/(2.*EpL
     +*sqrt(1. + (0.25*(1. - 1.*EpL**2 + sqrt(1. + 2.*EpL**2 + EpL
     +**4))**2)/EpL**2))
      vL0ch = -(1 - EpL**2 + EpRch**2 + sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRch**2) + (1. + EpRch**2)**2))/(2.*EpL*sqrt(1. + (0.25
     +*(1. - 1.*EpL**2 + EpRch**2 + sqrt(EpL**4 - 2.*EpL**2*(-1. +
     + EpRch**2) + (1. + EpRch**2)**2))**2)/EpL**2))
      vL0st = -(1 - EpL**2 + EpRst**2 + sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRst**2) + (1. + EpRst**2)**2))/(2.*EpL*sqrt(1. + (0.25
     +*(1. - 1.*EpL**2 + EpRst**2 + sqrt(EpL**4 - 2.*EpL**2*(-1. +
     + EpRst**2) + (1. + EpRst**2)**2))**2)/EpL**2))
      vL0to = -(1 - EpL**2 + EpRto**2 + sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRto**2) + (1. + EpRto**2)**2))/(2.*EpL*sqrt(1. + (0.25
     +*(1. - 1.*EpL**2 + EpRto**2 + sqrt(EpL**4 - 2.*EpL**2*(-1. +
     + EpRto**2) + (1. + EpRto**2)**2))**2)/EpL**2))
      vL0bo = -(1 - EpL**2 + EpRbo**2 + sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRbo**2) + (1. + EpRbo**2)**2))/(2.*EpL*sqrt(1. + (0.25
     +*(1. - 1.*EpL**2 + EpRbo**2 + sqrt(EpL**4 - 2.*EpL**2*(-1. +
     + EpRbo**2) + (1. + EpRbo**2)**2))**2)/EpL**2))
      vL1nu = 1/Sqrt(1 + (1 - EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**4
     +))**2/(4.*EpL**2))
      vL1el = 1/Sqrt(1 + (1 - EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**4
     +))**2/(4.*EpL**2))
      vL1mu = 1/Sqrt(1 + (1 - EpL**2 + EpRmu**2 + sqrt(EpL**4 - 2.*
     +EpL**2*(-1. + EpRmu**2) + (1. + EpRmu**2)**2))**2/(4.*EpL**2
     +))
      vL1ta = 1/Sqrt(1 + (1 - EpL**2 + EpRta**2 + sqrt(EpL**4 - 2.*
     +EpL**2*(-1. + EpRta**2) + (1. + EpRta**2)**2))**2/(4.*EpL**2
     +))
      vL1up = 1/Sqrt(1 + (1 - EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**4
     +))**2/(4.*EpL**2))
      vL1do = 1/Sqrt(1 + (1 - EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**4
     +))**2/(4.*EpL**2))
      vL1ch = 1/Sqrt(1 + (1 - EpL**2 + EpRch**2 + sqrt(EpL**4 - 2.*
     +EpL**2*(-1. + EpRch**2) + (1. + EpRch**2)**2))**2/(4.*EpL**2
     +))
      vL1st = 1/Sqrt(1 + (1 - EpL**2 + EpRst**2 + sqrt(EpL**4 - 2.*
     +EpL**2*(-1. + EpRst**2) + (1. + EpRst**2)**2))**2/(4.*EpL**2
     +))
      vL1to = 1/Sqrt(1 + (1 - EpL**2 + EpRto**2 + sqrt(EpL**4 - 2.*
     +EpL**2*(-1. + EpRto**2) + (1. + EpRto**2)**2))**2/(4.*EpL**2
     +))
      vL1bo = 1/Sqrt(1 + (1 - EpL**2 + EpRbo**2 + sqrt(EpL**4 - 2.*
     +EpL**2*(-1. + EpRbo**2) + (1. + EpRbo**2)**2))**2/(4.*EpL**2
     +))
      vR1nu = 0
      vR1el = 0
      vR1up = 0
      vR1do = 0
      vR1mu = -(-1 - EpL**2 + EpRmu**2 + sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRmu**2) + (1. + EpRmu**2)**2))/(2.*EpRmu*sqrt(1. + (0
     +.25*(-1. - 1.*EpL**2 + EpRmu**2 + sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRmu**2) + (1. + EpRmu**2)**2))**2)/EpRmu**2))
      vR1ta = -(-1 - EpL**2 + EpRta**2 + sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRta**2) + (1. + EpRta**2)**2))/(2.*EpRta*sqrt(1. + (0
     +.25*(-1. - 1.*EpL**2 + EpRta**2 + sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRta**2) + (1. + EpRta**2)**2))**2)/EpRta**2))
      vR1ch = -(-1 - EpL**2 + EpRch**2 + sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRch**2) + (1. + EpRch**2)**2))/(2.*EpRch*sqrt(1. + (0
     +.25*(-1. - 1.*EpL**2 + EpRch**2 + sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRch**2) + (1. + EpRch**2)**2))**2)/EpRch**2))
      vR1st = -(-1 - EpL**2 + EpRst**2 + sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRst**2) + (1. + EpRst**2)**2))/(2.*EpRst*sqrt(1. + (0
     +.25*(-1. - 1.*EpL**2 + EpRst**2 + sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRst**2) + (1. + EpRst**2)**2))**2)/EpRst**2))
      vR1to = -(-1 - EpL**2 + EpRto**2 + sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRto**2) + (1. + EpRto**2)**2))/(2.*EpRto*sqrt(1. + (0
     +.25*(-1. - 1.*EpL**2 + EpRto**2 + sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRto**2) + (1. + EpRto**2)**2))**2)/EpRto**2))
      vR1bo = -(-1 - EpL**2 + EpRbo**2 + sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRbo**2) + (1. + EpRbo**2)**2))/(2.*EpRbo*sqrt(1. + (0
     +.25*(-1. - 1.*EpL**2 + EpRbo**2 + sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRbo**2) + (1. + EpRbo**2)**2))**2)/EpRbo**2))
      vR2nu = 1
      vR2el = 1
      vR2up = 1
      vR2do = 1
      vR2mu = 1/Sqrt(1 + (-1 - EpL**2 + EpRmu**2 + sqrt(EpL**4 - 2.
     +*EpL**2*(-1. + EpRmu**2) + (1. + EpRmu**2)**2))**2/(4.*EpRmu
     +**2))
      vR2ta = 1/Sqrt(1 + (-1 - EpL**2 + EpRta**2 + sqrt(EpL**4 - 2.
     +*EpL**2*(-1. + EpRta**2) + (1. + EpRta**2)**2))**2/(4.*EpRta
     +**2))
      vR2ch = 1/Sqrt(1 + (-1 - EpL**2 + EpRch**2 + sqrt(EpL**4 - 2.
     +*EpL**2*(-1. + EpRch**2) + (1. + EpRch**2)**2))**2/(4.*EpRch
     +**2))
      vR2st = 1/Sqrt(1 + (-1 - EpL**2 + EpRst**2 + sqrt(EpL**4 - 2.
     +*EpL**2*(-1. + EpRst**2) + (1. + EpRst**2)**2))**2/(4.*EpRst
     +**2))
      vR2to = 1/Sqrt(1 + (-1 - EpL**2 + EpRto**2 + sqrt(EpL**4 - 2.
     +*EpL**2*(-1. + EpRto**2) + (1. + EpRto**2)**2))**2/(4.*EpRto
     +**2))
      vR2bo = 1/Sqrt(1 + (-1 - EpL**2 + EpRbo**2 + sqrt(EpL**4 - 2.
     +*EpL**2*(-1. + EpRbo**2) + (1. + EpRbo**2)**2))**2/(4.*EpRbo
     +**2))
      vL0hNu = (1 - EpL**2 - sqrt(1. + 2.*EpL**2 + EpL**4))/(2.*EpL
     +*sqrt(1. + (0.25*(-1. + EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**
     +4))**2)/EpL**2))
      vL0hEl = (1 - EpL**2 - sqrt(1. + 2.*EpL**2 + EpL**4))/(2.*EpL
     +*sqrt(1. + (0.25*(-1. + EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**
     +4))**2)/EpL**2))
      vL0hMu = (1 - EpL**2 + EpRmu**2 - sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRmu**2) + (1. + EpRmu**2)**2))/(2.*EpL*sqrt(1. + (0.25
     +*(-1. + EpL**2 - 1.*EpRmu**2 + sqrt(EpL**4 - 2.*EpL**2*(-1. 
     ++ EpRmu**2) + (1. + EpRmu**2)**2))**2)/EpL**2))
      vL0hTa = (1 - EpL**2 + EpRta**2 - sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRta**2) + (1. + EpRta**2)**2))/(2.*EpL*sqrt(1. + (0.25
     +*(-1. + EpL**2 - 1.*EpRta**2 + sqrt(EpL**4 - 2.*EpL**2*(-1. 
     ++ EpRta**2) + (1. + EpRta**2)**2))**2)/EpL**2))
      vL0hUp = (1 - EpL**2 - sqrt(1. + 2.*EpL**2 + EpL**4))/(2.*EpL
     +*sqrt(1. + (0.25*(-1. + EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**
     +4))**2)/EpL**2))
      vL0hDo = (1 - EpL**2 - sqrt(1. + 2.*EpL**2 + EpL**4))/(2.*EpL
     +*sqrt(1. + (0.25*(-1. + EpL**2 + sqrt(1. + 2.*EpL**2 + EpL**
     +4))**2)/EpL**2))
      vL0hCh = (1 - EpL**2 + EpRch**2 - sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRch**2) + (1. + EpRch**2)**2))/(2.*EpL*sqrt(1. + (0.25
     +*(-1. + EpL**2 - 1.*EpRch**2 + sqrt(EpL**4 - 2.*EpL**2*(-1. 
     ++ EpRch**2) + (1. + EpRch**2)**2))**2)/EpL**2))
      vL0hSt = (1 - EpL**2 + EpRst**2 - sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRst**2) + (1. + EpRst**2)**2))/(2.*EpL*sqrt(1. + (0.25
     +*(-1. + EpL**2 - 1.*EpRst**2 + sqrt(EpL**4 - 2.*EpL**2*(-1. 
     ++ EpRst**2) + (1. + EpRst**2)**2))**2)/EpL**2))
      vL0hTo = (1 - EpL**2 + EpRto**2 - sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRto**2) + (1. + EpRto**2)**2))/(2.*EpL*sqrt(1. + (0.25
     +*(-1. + EpL**2 - 1.*EpRto**2 + sqrt(EpL**4 - 2.*EpL**2*(-1. 
     ++ EpRto**2) + (1. + EpRto**2)**2))**2)/EpL**2))
      vL0hBo = (1 - EpL**2 + EpRbo**2 - sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRbo**2) + (1. + EpRbo**2)**2))/(2.*EpL*sqrt(1. + (0.25
     +*(-1. + EpL**2 - 1.*EpRbo**2 + sqrt(EpL**4 - 2.*EpL**2*(-1. 
     ++ EpRbo**2) + (1. + EpRbo**2)**2))**2)/EpL**2))
      vL1hNu = -(1/Sqrt(1 + (-1 + EpL**2 + sqrt(1. + 2.*EpL**2 + Ep
     +L**4))**2/(4.*EpL**2)))
      vL1hEl = -(1/Sqrt(1 + (-1 + EpL**2 + sqrt(1. + 2.*EpL**2 + Ep
     +L**4))**2/(4.*EpL**2)))
      vL1hMu = -(1/Sqrt(1 + (-1 + EpL**2 - EpRmu**2 + sqrt(EpL**4 -
     + 2.*EpL**2*(-1. + EpRmu**2) + (1. + EpRmu**2)**2))**2/(4.*Ep
     +L**2)))
      vL1hTa = -(1/Sqrt(1 + (-1 + EpL**2 - EpRta**2 + sqrt(EpL**4 -
     + 2.*EpL**2*(-1. + EpRta**2) + (1. + EpRta**2)**2))**2/(4.*Ep
     +L**2)))
      vL1hUp = -(1/Sqrt(1 + (-1 + EpL**2 + sqrt(1. + 2.*EpL**2 + Ep
     +L**4))**2/(4.*EpL**2)))
      vL1hDo = -(1/Sqrt(1 + (-1 + EpL**2 + sqrt(1. + 2.*EpL**2 + Ep
     +L**4))**2/(4.*EpL**2)))
      vL1hCh = -(1/Sqrt(1 + (-1 + EpL**2 - EpRch**2 + sqrt(EpL**4 -
     + 2.*EpL**2*(-1. + EpRch**2) + (1. + EpRch**2)**2))**2/(4.*Ep
     +L**2)))
      vL1hSt = -(1/Sqrt(1 + (-1 + EpL**2 - EpRst**2 + sqrt(EpL**4 -
     + 2.*EpL**2*(-1. + EpRst**2) + (1. + EpRst**2)**2))**2/(4.*Ep
     +L**2)))
      vL1hTo = -(1/Sqrt(1 + (-1 + EpL**2 - EpRto**2 + sqrt(EpL**4 -
     + 2.*EpL**2*(-1. + EpRto**2) + (1. + EpRto**2)**2))**2/(4.*Ep
     +L**2)))
      vL1hBo = -(1/Sqrt(1 + (-1 + EpL**2 - EpRbo**2 + sqrt(EpL**4 -
     + 2.*EpL**2*(-1. + EpRbo**2) + (1. + EpRbo**2)**2))**2/(4.*Ep
     +L**2)))
      vR1hNu = -1
      vR1hEl = -1
      vR1hUp = -1
      vR1hDo = -1
      vR1hMu = (-1 - EpL**2 + EpRmu**2 - sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRmu**2) + (1. + EpRmu**2)**2))/(2.*EpRmu*sqrt(1. + (0
     +.25*(1. + EpL**2 - 1.*EpRmu**2 + sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRmu**2) + (1. + EpRmu**2)**2))**2)/EpRmu**2))
      vR1hTa = (-1 - EpL**2 + EpRta**2 - sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRta**2) + (1. + EpRta**2)**2))/(2.*EpRta*sqrt(1. + (0
     +.25*(1. + EpL**2 - 1.*EpRta**2 + sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRta**2) + (1. + EpRta**2)**2))**2)/EpRta**2))
      vR1hCh = (-1 - EpL**2 + EpRch**2 - sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRch**2) + (1. + EpRch**2)**2))/(2.*EpRch*sqrt(1. + (0
     +.25*(1. + EpL**2 - 1.*EpRch**2 + sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRch**2) + (1. + EpRch**2)**2))**2)/EpRch**2))
      vR1hSt = (-1 - EpL**2 + EpRst**2 - sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRst**2) + (1. + EpRst**2)**2))/(2.*EpRst*sqrt(1. + (0
     +.25*(1. + EpL**2 - 1.*EpRst**2 + sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRst**2) + (1. + EpRst**2)**2))**2)/EpRst**2))
      vR1hTo = (-1 - EpL**2 + EpRto**2 - sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRto**2) + (1. + EpRto**2)**2))/(2.*EpRto*sqrt(1. + (0
     +.25*(1. + EpL**2 - 1.*EpRto**2 + sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRto**2) + (1. + EpRto**2)**2))**2)/EpRto**2))
      vR1hBo = (-1 - EpL**2 + EpRbo**2 - sqrt(EpL**4 - 2.*EpL**2*(-
     +1. + EpRbo**2) + (1. + EpRbo**2)**2))/(2.*EpRbo*sqrt(1. + (0
     +.25*(1. + EpL**2 - 1.*EpRbo**2 + sqrt(EpL**4 - 2.*EpL**2*(-1
     +. + EpRbo**2) + (1. + EpRbo**2)**2))**2)/EpRbo**2))
      vR2hNu = 0
      vR2hEl = 0
      vR2hUp = 0
      vR2hDo = 0
      vR2hMu = -(1/Sqrt(1 + (1 + EpL**2 - EpRmu**2 + sqrt(EpL**4 - 
     +2.*EpL**2*(-1. + EpRmu**2) + (1. + EpRmu**2)**2))**2/(4.*EpR
     +mu**2)))
      vR2hTa = -(1/Sqrt(1 + (1 + EpL**2 - EpRta**2 + sqrt(EpL**4 - 
     +2.*EpL**2*(-1. + EpRta**2) + (1. + EpRta**2)**2))**2/(4.*EpR
     +ta**2)))
      vR2hCh = -(1/Sqrt(1 + (1 + EpL**2 - EpRch**2 + sqrt(EpL**4 - 
     +2.*EpL**2*(-1. + EpRch**2) + (1. + EpRch**2)**2))**2/(4.*EpR
     +ch**2)))
      vR2hSt = -(1/Sqrt(1 + (1 + EpL**2 - EpRst**2 + sqrt(EpL**4 - 
     +2.*EpL**2*(-1. + EpRst**2) + (1. + EpRst**2)**2))**2/(4.*EpR
     +st**2)))
      vR2hTo = -(1/Sqrt(1 + (1 + EpL**2 - EpRto**2 + sqrt(EpL**4 - 
     +2.*EpL**2*(-1. + EpRto**2) + (1. + EpRto**2)**2))**2/(4.*EpR
     +to**2)))
      vR2hBo = -(1/Sqrt(1 + (1 + EpL**2 - EpRbo**2 + sqrt(EpL**4 - 
     +2.*EpL**2*(-1. + EpRbo**2) + (1. + EpRbo**2)**2))**2/(4.*EpR
     +bo**2)))
      Sqrt2 = sqrt(2.)
      SqrtPi = sqrt(3.141592653589793)
      Pv0W2 = v0W**2
      Pv1W2 = v1W**2
      PEE2 = EE**2
      Pg02 = g0**2
      Pv0W4 = v0W**4
      Pgt2 = gt**2
      Pv1W4 = v1W**4
      Pv0WP2 = v0WP**2
      Pv1WP2 = v1WP**2
      Pv0W3 = v0W**3
      Pv1W3 = v1W**3
      Pv0WP3 = v0WP**3
      Pv1WP3 = v1WP**3
      Pv0WP4 = v0WP**4
      Pv1WP4 = v1WP**4
      Pv0Z2 = v0Z**2
      Pv1Z2 = v1Z**2
      Pv0ZP2 = v0ZP**2
      Pv1ZP2 = v1ZP**2
      PvL0do2 = vL0do**2
      PvL1do2 = vL1do**2
      PvL0st2 = vL0st**2
      PvL1st2 = vL1st**2
      PvR1st2 = vR1st**2
      PvR2st2 = vR2st**2
      PvL0bo2 = vL0bo**2
      PvL1bo2 = vL1bo**2
      PvR1bo2 = vR1bo**2
      PvR2bo2 = vR2bo**2
      PvL0hDo2 = vL0hDo**2
      PvL1hDo2 = vL1hDo**2
      PvL0hSt2 = vL0hSt**2
      PvL1hSt2 = vL1hSt**2
      PvR1hSt2 = vR1hSt**2
      PvR2hSt2 = vR2hSt**2
      PvL0hBo2 = vL0hBo**2
      PvL1hBo2 = vL1hBo**2
      PvR1hBo2 = vR1hBo**2
      PvR2hBo2 = vR2hBo**2
      PvL0el2 = vL0el**2
      PvL1el2 = vL1el**2
      PvL0mu2 = vL0mu**2
      PvL1mu2 = vL1mu**2
      PvR1mu2 = vR1mu**2
      PvR2mu2 = vR2mu**2
      PvL0ta2 = vL0ta**2
      PvL1ta2 = vL1ta**2
      PvR1ta2 = vR1ta**2
      PvR2ta2 = vR2ta**2
      PvL0hEl2 = vL0hEl**2
      PvL1hEl2 = vL1hEl**2
      PvL0hMu2 = vL0hMu**2
      PvL1hMu2 = vL1hMu**2
      PvR1hMu2 = vR1hMu**2
      PvR2hMu2 = vR2hMu**2
      PvL0hTa2 = vL0hTa**2
      PvL1hTa2 = vL1hTa**2
      PvR1hTa2 = vR1hTa**2
      PvR2hTa2 = vR2hTa**2
      PvL0hUp2 = vL0hUp**2
      PvL1hUp2 = vL1hUp**2
      PvL0hCh2 = vL0hCh**2
      PvL1hCh2 = vL1hCh**2
      PvR1hCh2 = vR1hCh**2
      PvR2hCh2 = vR2hCh**2
      PvL0hTo2 = vL0hTo**2
      PvL1hTo2 = vL1hTo**2
      PvR1hTo2 = vR1hTo**2
      PvR2hTo2 = vR2hTo**2
      PvL0up2 = vL0up**2
      PvL1up2 = vL1up**2
      PvL0ch2 = vL0ch**2
      PvL1ch2 = vL1ch**2
      PvR1ch2 = vR1ch**2
      PvR2ch2 = vR2ch**2
      PvL0to2 = vL0to**2
      PvL1to2 = vL1to**2
      PvR1to2 = vR1to**2
      PvR2to2 = vR2to**2
      PvL0hNu2 = vL0hNu**2
      PvL1hNu2 = vL1hNu**2
      PvL0nu2 = vL0nu**2
      PvL1nu2 = vL1nu**2
      endif
c Parameters that should be recomputed at an event by even basis.
      GGT1 = sqrt(G**2)


c Definition of the EW coupling used in the write out of aqed
      gal(1) = ee
      gal(2) = ee


c Definition of DUM symbols
      DUM0 = 0
      DUM1 = 1


