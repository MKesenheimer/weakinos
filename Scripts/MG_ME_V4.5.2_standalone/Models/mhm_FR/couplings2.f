      subroutine coup2(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 2
      if(readlha) then
      MGVX3 = EE*Pv0W2 + EE*Pv1W2
      MGVX5 = -(EE*v0W*v0WP) - EE*v1W*v1WP
      MGVX4 = EE*v0W*v0WP + EE*v1W*v1WP
      MGVX10 = EE*Pv0WP2 + EE*Pv1WP2
      MGVX7 = -(g0*Pv0W2*v0Z) - gt*Pv1W2*v1Z
      MGVX17 = -(g0*v0W*v0WP*v0Z) - gt*v1W*v1WP*v1Z
      MGVX12 = g0*v0W*v0WP*v0Z + gt*v1W*v1WP*v1Z
      MGVX22 = -(g0*Pv0WP2*v0Z) - gt*Pv1WP2*v1Z
      MGVX8 = -(g0*Pv0W2*v0ZP) - gt*Pv1W2*v1ZP
      MGVX18 = -(g0*v0W*v0WP*v0ZP) - gt*v1W*v1WP*v1ZP
      MGVX13 = g0*v0W*v0WP*v0ZP + gt*v1W*v1WP*v1ZP
      MGVX23 = -(g0*Pv0WP2*v0ZP) - gt*Pv1WP2*v1ZP
      MGVX6 = PEE2*Pv0W2 + PEE2*Pv1W2
      MGVX16 = PEE2*v0W*v0WP + PEE2*v1W*v1WP
      MGVX20 = PEE2*Pv0WP2 + PEE2*Pv1WP2


      endif

      return
      end
