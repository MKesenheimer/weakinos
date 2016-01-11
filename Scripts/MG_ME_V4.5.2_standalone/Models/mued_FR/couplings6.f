      subroutine coup6(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 6
      if(readlha) then


      MGVX160(1) = -((CONJCKM33*gw)/Sqrt2)
      MGVX160(2) = 0
      MGVX161(1) = -((CONJCKM13*gw)/Sqrt2)
      MGVX161(2) = 0
      MGVX162(1) = -((CONJCKM21*gw)/Sqrt2)
      MGVX162(2) = 0
      MGVX163(1) = -((CONJCKM31*gw)/Sqrt2)
      MGVX163(2) = 0
      MGVX164(1) = -((CONJCKM11*gw)/Sqrt2)
      MGVX164(2) = 0
      MGVX174(1) = -((CONJCKM22*gw)/Sqrt2)
      MGVX174(2) = 0
      MGVX175(1) = -((CONJCKM32*gw)/Sqrt2)
      MGVX175(2) = 0
      MGVX176(1) = -((CONJCKM12*gw)/Sqrt2)
      MGVX176(2) = 0
      MGVX186(1) = -((CONJCKM23*gw)/Sqrt2)
      MGVX186(2) = -((CONJCKM23*gw)/Sqrt2)
      MGVX187(1) = -((CONJCKM33*gw)/Sqrt2)
      MGVX187(2) = -((CONJCKM33*gw)/Sqrt2)
      MGVX188(1) = -((CONJCKM13*gw)/Sqrt2)
      MGVX188(2) = -((CONJCKM13*gw)/Sqrt2)
      MGVX189(1) = -((CONJCKM21*gw)/Sqrt2)
      MGVX189(2) = -((CONJCKM21*gw)/Sqrt2)
      MGVX190(1) = -((CONJCKM31*gw)/Sqrt2)
      MGVX190(2) = -((CONJCKM31*gw)/Sqrt2)
      MGVX191(1) = -((CONJCKM11*gw)/Sqrt2)
      MGVX191(2) = -((CONJCKM11*gw)/Sqrt2)
      MGVX192(1) = -((CONJCKM22*gw)/Sqrt2)
      MGVX192(2) = -((CONJCKM22*gw)/Sqrt2)
      endif

      return
      end
