      subroutine coup2(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 2
      if(readlha) then
      MGVX9 = -(gw*sw)
      MGVX18 = cw*gw
      MGVX10 = gw
      MGVX17 = Pgw2*Psw2
      MGVX33 = Pgw2*sw
      MGVX28 = cw*Pgw2*sw
      MGVX13 = -Pgw2
      MGVX35 = cw*Pgw2
      MGVX32 = Pcw2*Pgw2
      MGVX38 = Pgw2
      MGVX14 = (-3*Pgw2)/2.
      MGVX37 = (3*Pgw2)/2.


      MGVX24(1) = -(ytau/Sqrt2)
      MGVX24(2) = -(ytau/Sqrt2)
      MGVX39(1) = gw/2.
      MGVX39(2) = 0
      MGVX45(1) = -gw/2.
      MGVX45(2) = 0
      endif

      return
      end
