      subroutine coup18(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 18
      if(readlha) then


      MGVX416(1) = -(CONJNN11*CONJRd44*ee*Sqrt2)/(3.*cw)
      MGVX416(2) = 0
      MGVX417(1) = -(CONJNN21*CONJRd44*ee*Sqrt2)/(3.*cw)
      MGVX417(2) = 0
      MGVX418(1) = -(CONJNN31*CONJRd44*ee*Sqrt2)/(3.*cw)
      MGVX418(2) = 0
      MGVX419(1) = -(CONJNN41*CONJRd44*ee*Sqrt2)/(3.*cw)
      MGVX419(2) = 0
      MGVX472(1) = 0
      MGVX472(2) = -(ee*NN11*Rd44*Sqrt2)/(3.*cw)
      MGVX473(1) = 0
      MGVX473(2) = -(ee*NN21*Rd44*Sqrt2)/(3.*cw)
      MGVX474(1) = 0
      MGVX474(2) = -(ee*NN31*Rd44*Sqrt2)/(3.*cw)
      MGVX475(1) = 0
      MGVX475(2) = -(ee*NN41*Rd44*Sqrt2)/(3.*cw)
      MGVX433(1) = 0
      MGVX433(2) = -(CONJRd51*ee*NN11)/(3.*cw*Sqrt2) + (CONJRd51*ee
     +*NN12)/(Sqrt2*sw)
      MGVX440(1) = 0
      MGVX440(2) = -(CONJRd51*ee*NN21)/(3.*cw*Sqrt2) + (CONJRd51*ee
     +*NN22)/(Sqrt2*sw)
      MGVX447(1) = 0
      MGVX447(2) = -(CONJRd51*ee*NN31)/(3.*cw*Sqrt2) + (CONJRd51*ee
     +*NN32)/(Sqrt2*sw)
      MGVX454(1) = 0
      MGVX454(2) = -(CONJRd51*ee*NN41)/(3.*cw*Sqrt2) + (CONJRd51*ee
     +*NN42)/(Sqrt2*sw)
      MGVX632(1) = 0
      MGVX632(2) = -((CONJRd51*ee*UU11)/sw)
      MGVX640(1) = 0
      MGVX640(2) = -((CONJRd51*ee*UU21)/sw)
      MGVX476(1) = -((CONJUU11*ee*Rd51)/sw)
      MGVX476(2) = 0
      endif

      return
      end
