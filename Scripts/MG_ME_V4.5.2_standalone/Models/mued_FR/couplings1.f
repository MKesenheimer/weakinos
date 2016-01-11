      subroutine coup1(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 1
      MGVX3 = G
      MGVX4 = GGT1
      MGVX6 = G**2
      MGVX7 = G1G1T2


      GG(1) = -G
      GG(2) = -G
      MGVX111(1) = -G
      MGVX111(2) = 0
      MGVX297(1) = 0
      MGVX297(2) = -G

      return
      end
