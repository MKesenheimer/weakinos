      subroutine coup2(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 2
      if(readlha) then
      MGVX3 = -(gw*sw)
      MGVX5 = cw*gw
      MGVX4 = Pgw2*Psw2
      MGVX7 = cw*Pgw2*sw
      MGVX6 = -Pgw2
      MGVX8 = Pcw2*Pgw2

      MGVX53 = -6*lam*v
      MGVX52 = -6*lam
      MGVX55 = (ee**2*v)/(2.*Psw2)
      MGVX57 = ee**2*v + (ee**2*Pcw2*v)/(2.*Psw2) + (ee**2*Psw2*v)/
     +(2.*Pcw2)
      MGVX54 = ee**2/(2.*Psw2)
      MGVX56 = ee**2 + (ee**2*Pcw2)/(2.*Psw2) + (ee**2*Psw2)/(2.*Pc
     +w2)

      MGVX18(1) = ee
      MGVX18(2) = ee
      MGVX24(1) = -(ee/(Sqrt2*sw))
      MGVX24(2) = 0
      MGVX43(1) = (cw*ee)/(2.*sw) - (ee*sw)/(2.*cw)
      MGVX43(2) = -((ee*sw)/cw)
      MGVX49(1) = -(cw*ee)/(2.*sw) - (ee*sw)/(2.*cw)
      MGVX49(2) = 0
      MGVX59(1) = -(ytau/Sqrt2)
      MGVX59(2) = -(ytau/Sqrt2)
      MGVX15(1) = ee/3.
      MGVX15(2) = ee/3.
      MGVX21(1) = (-2*ee)/3.
      MGVX21(2) = (-2*ee)/3.
      MGVX58(1) = -(yb/Sqrt2)
      MGVX58(2) = -(yb/Sqrt2)
      MGVX60(1) = -(yc/Sqrt2)
      MGVX60(2) = -(yc/Sqrt2)
      MGVX61(1) = -(yt/Sqrt2)
      MGVX61(2) = -(yt/Sqrt2)
      MGVX27(1) = -((CKM21*ee)/(Sqrt2*sw))
      MGVX27(2) = 0
      MGVX28(1) = -((CKM22*ee)/(Sqrt2*sw))
      MGVX28(2) = 0
      MGVX30(1) = -((CKM11*ee)/(Sqrt2*sw))
      MGVX30(2) = 0
      MGVX31(1) = -((CKM12*ee)/(Sqrt2*sw))
      MGVX31(2) = 0
      MGVX36(1) = -((CONJCKM21*ee)/(Sqrt2*sw))
      MGVX36(2) = 0
      MGVX37(1) = -((CONJCKM11*ee)/(Sqrt2*sw))
      MGVX37(2) = 0
      MGVX38(1) = -((CONJCKM22*ee)/(Sqrt2*sw))
      MGVX38(2) = 0
      MGVX39(1) = -((CONJCKM12*ee)/(Sqrt2*sw))
      MGVX39(2) = 0
      MGVX40(1) = (cw*ee)/(2.*sw) + (ee*sw)/(6.*cw)
      MGVX40(2) = -(ee*sw)/(3.*cw)
      MGVX46(1) = -(cw*ee)/(2.*sw) + (ee*sw)/(6.*cw)
      MGVX46(2) = (2*ee*sw)/(3.*cw)
      endif

      return
      end
