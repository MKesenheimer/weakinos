      subroutine coup3(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 3
      if(readlha) then


      MGVX69(1) = -(cw*g1*y2)/2.
      MGVX69(2) = -(cw*g1*y2)/2.
      MGVX72(1) = (g1*sw*y2)/2.
      MGVX72(2) = (g1*sw*y2)/2.
      MGVX219(1) = (gw*sw)/2. - (cw*g1*y1)/2.
      MGVX219(2) = -(cw*g1*y2)/2.
      MGVX222(1) = -(gw*sw)/2. - (cw*g1*y1)/2.
      MGVX222(2) = 0
      MGVX225(1) = -(g1*y1)/2.
      MGVX225(2) = 0
      MGVX231(1) = (gw*sw)/2. - (cw*g1*y1)/2.
      MGVX231(2) = (gw*sw)/2. - (cw*g1*y1)/2.
      MGVX234(1) = -(gw*sw)/2. - (cw*g1*y1)/2.
      MGVX234(2) = -(gw*sw)/2. - (cw*g1*y1)/2.
      MGVX243(1) = -(gw/Sqrt2)
      MGVX243(2) = 0
      MGVX246(1) = -(gw/Sqrt2)
      MGVX246(2) = -(gw/Sqrt2)
      MGVX267(1) = (cw*gw)/2. + (g1*sw*y1)/2.
      MGVX267(2) = (g1*sw*y2)/2.
      MGVX270(1) = -(cw*gw)/2. + (g1*sw*y1)/2.
      MGVX270(2) = 0
      MGVX273(1) = (cw*gw)/2. + (g1*sw*y1)/2.
      MGVX273(2) = (cw*gw)/2. + (g1*sw*y1)/2.
      MGVX276(1) = -(cw*gw)/2. + (g1*sw*y1)/2.
      MGVX276(2) = -(cw*gw)/2. + (g1*sw*y1)/2.
      MGVX303(1) = 0
      MGVX303(2) = -(g1*y2)/2.
      MGVX75(1) = (gw*sw)/2. - (cw*g1*y3)/2.
      MGVX75(2) = -(cw*g1*y5)/2.
      endif

      return
      end
