      subroutine coup3(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 3
      if(readlha) then
      MGVX28 = EE*g0*Pv0W2*v0Z + EE*gt*Pv1W2*v1Z
      MGVX36 = EE*g0*Pv0W2*v0ZP + EE*gt*Pv1W2*v1ZP
      MGVX30 = EE*g0*v0W*v0WP*v0Z + EE*gt*v1W*v1WP*v1Z
      MGVX38 = EE*g0*v0W*v0WP*v0ZP + EE*gt*v1W*v1WP*v1ZP
      MGVX31 = EE*g0*Pv0WP2*v0Z + EE*gt*Pv1WP2*v1Z
      MGVX39 = EE*g0*Pv0WP2*v0ZP + EE*gt*Pv1WP2*v1ZP
      MGVX9 = -(Pg02*Pv0W4) - Pgt2*Pv1W4
      MGVX19 = -(Pg02*Pv0W3*v0WP) - Pgt2*Pv1W3*v1WP
      MGVX25 = -(Pg02*Pv0W2*Pv0WP2) - Pgt2*Pv1W2*Pv1WP2
      MGVX32 = Pg02*Pv0W2*Pv0Z2 + Pgt2*Pv1W2*Pv1Z2
      MGVX40 = Pg02*Pv0W2*v0Z*v0ZP + Pgt2*Pv1W2*v1Z*v1ZP
      MGVX44 = Pg02*Pv0W2*Pv0ZP2 + Pgt2*Pv1W2*Pv1ZP2
      MGVX26 = -(Pg02*Pv0WP3*v0W) - Pgt2*Pv1WP3*v1W
      MGVX34 = Pg02*Pv0Z2*v0W*v0WP + Pgt2*Pv1Z2*v1W*v1WP
      MGVX42 = Pg02*v0W*v0WP*v0Z*v0ZP + Pgt2*v1W*v1WP*v1Z*v1ZP


      endif

      return
      end
