      subroutine coup7(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 7
      if(readlha) then

      MGVX2 = -6*lam*v
      MGVX1 = -6*lam
      MGVX16 = (ee**2*v)/(2.*Psw2)
      MGVX30 = ee**2*v + (ee**2*Pcw2*v)/(2.*Psw2) + (ee**2*Psw2*v)/
     +(2.*Pcw2)
      MGVX15 = ee**2/(2.*Psw2)
      MGVX29 = ee**2 + (ee**2*Pcw2)/(2.*Psw2) + (ee**2*Psw2)/(2.*Pc
     +w2)

      MGVX193(1) = -((CONJCKM32*gw)/Sqrt2)
      MGVX193(2) = -((CONJCKM32*gw)/Sqrt2)
      MGVX194(1) = -((CONJCKM12*gw)/Sqrt2)
      MGVX194(2) = -((CONJCKM12*gw)/Sqrt2)
      MGVX195(1) = (cw*gw)/2. + (g1*sw*y3)/2.
      MGVX195(2) = (g1*sw*y5)/2.
      MGVX198(1) = -(cw*gw)/2. + (g1*sw*y3)/2.
      MGVX198(2) = (g1*sw*y4)/2.
      MGVX204(1) = -(cw*gw)/2. + (g1*sw*y3)/2.
      MGVX204(2) = -(cw*gw)/2. + (g1*sw*y3)/2.
      MGVX201(1) = (cw*gw)/2. + (g1*sw*y3)/2.
      MGVX201(2) = (cw*gw)/2. + (g1*sw*y3)/2.
      MGVX66(1) = (g1*sw*y5)/2.
      MGVX66(2) = (g1*sw*y5)/2.
      MGVX57(1) = (g1*sw*y4)/2.
      MGVX57(2) = (g1*sw*y4)/2.
      endif

      return
      end
