      subroutine coup5(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 5
      if(readlha) then


      MGVX133(1) = -((CKM21*gw)/Sqrt2)
      MGVX133(2) = -((CKM21*gw)/Sqrt2)
      MGVX134(1) = -((CKM22*gw)/Sqrt2)
      MGVX134(2) = -((CKM22*gw)/Sqrt2)
      MGVX126(1) = -((CKM33*gw)/Sqrt2)
      MGVX126(2) = 0
      MGVX127(1) = -((CKM31*gw)/Sqrt2)
      MGVX127(2) = 0
      MGVX128(1) = -((CKM32*gw)/Sqrt2)
      MGVX128(2) = 0
      MGVX135(1) = -((CKM33*gw)/Sqrt2)
      MGVX135(2) = -((CKM33*gw)/Sqrt2)
      MGVX136(1) = -((CKM31*gw)/Sqrt2)
      MGVX136(2) = -((CKM31*gw)/Sqrt2)
      MGVX137(1) = -((CKM32*gw)/Sqrt2)
      MGVX137(2) = -((CKM32*gw)/Sqrt2)
      MGVX129(1) = -((CKM13*gw)/Sqrt2)
      MGVX129(2) = 0
      MGVX130(1) = -((CKM11*gw)/Sqrt2)
      MGVX130(2) = 0
      MGVX131(1) = -((CKM12*gw)/Sqrt2)
      MGVX131(2) = 0
      MGVX138(1) = -((CKM13*gw)/Sqrt2)
      MGVX138(2) = -((CKM13*gw)/Sqrt2)
      MGVX139(1) = -((CKM11*gw)/Sqrt2)
      MGVX139(2) = -((CKM11*gw)/Sqrt2)
      MGVX140(1) = -((CKM12*gw)/Sqrt2)
      MGVX140(2) = -((CKM12*gw)/Sqrt2)
      MGVX159(1) = -((CONJCKM23*gw)/Sqrt2)
      MGVX159(2) = 0
      endif

      return
      end
