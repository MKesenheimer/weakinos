      subroutine coup4(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 4
      if(readlha) then


      MGVX78(1) = -(gw*sw)/2. - (cw*g1*y3)/2.
      MGVX78(2) = -(cw*g1*y4)/2.
      MGVX90(1) = -(gw*sw)/2. - (cw*g1*y3)/2.
      MGVX90(2) = -(gw*sw)/2. - (cw*g1*y3)/2.
      MGVX87(1) = (gw*sw)/2. - (cw*g1*y3)/2.
      MGVX87(2) = (gw*sw)/2. - (cw*g1*y3)/2.
      MGVX60(1) = -(cw*g1*y5)/2.
      MGVX60(2) = -(cw*g1*y5)/2.
      MGVX51(1) = -(cw*g1*y4)/2.
      MGVX51(2) = -(cw*g1*y4)/2.
      MGVX81(1) = -(g1*y3)/2.
      MGVX81(2) = 0
      MGVX291(1) = 0
      MGVX291(2) = -(g1*y5)/2.
      MGVX282(1) = 0
      MGVX282(2) = -(g1*y4)/2.
      MGVX23(1) = -(yb/Sqrt2)
      MGVX23(2) = -(yb/Sqrt2)
      MGVX25(1) = -(yc/Sqrt2)
      MGVX25(2) = -(yc/Sqrt2)
      MGVX26(1) = -(yt/Sqrt2)
      MGVX26(2) = -(yt/Sqrt2)
      MGVX123(1) = -((CKM23*gw)/Sqrt2)
      MGVX123(2) = 0
      MGVX124(1) = -((CKM21*gw)/Sqrt2)
      MGVX124(2) = 0
      MGVX125(1) = -((CKM22*gw)/Sqrt2)
      MGVX125(2) = 0
      MGVX132(1) = -((CKM23*gw)/Sqrt2)
      MGVX132(2) = -((CKM23*gw)/Sqrt2)
      endif

      return
      end
