      subroutine coup19(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 19
      if(readlha) then


      MGVX477(1) = -((CONJUU21*ee*Rd51)/sw)
      MGVX477(2) = 0
      MGVX478(1) = -(CONJNN11*ee*Rd51)/(3.*cw*Sqrt2) + (CONJNN12*ee
     +*Rd51)/(Sqrt2*sw)
      MGVX478(2) = 0
      MGVX479(1) = -(CONJNN21*ee*Rd51)/(3.*cw*Sqrt2) + (CONJNN22*ee
     +*Rd51)/(Sqrt2*sw)
      MGVX479(2) = 0
      MGVX480(1) = -(CONJNN31*ee*Rd51)/(3.*cw*Sqrt2) + (CONJNN32*ee
     +*Rd51)/(Sqrt2*sw)
      MGVX480(2) = 0
      MGVX481(1) = -(CONJNN41*ee*Rd51)/(3.*cw*Sqrt2) + (CONJNN42*ee
     +*Rd51)/(Sqrt2*sw)
      MGVX481(2) = 0
      MGVX633(1) = 0
      MGVX633(2) = -((CONJRd62*ee*UU11)/sw)
      MGVX641(1) = 0
      MGVX641(2) = -((CONJRd62*ee*UU21)/sw)
      MGVX434(1) = 0
      MGVX434(2) = -(CONJRd62*ee*NN11)/(3.*cw*Sqrt2) + (CONJRd62*ee
     +*NN12)/(Sqrt2*sw)
      MGVX441(1) = 0
      MGVX441(2) = -(CONJRd62*ee*NN21)/(3.*cw*Sqrt2) + (CONJRd62*ee
     +*NN22)/(Sqrt2*sw)
      MGVX448(1) = 0
      MGVX448(2) = -(CONJRd62*ee*NN31)/(3.*cw*Sqrt2) + (CONJRd62*ee
     +*NN32)/(Sqrt2*sw)
      MGVX455(1) = 0
      MGVX455(2) = -(CONJRd62*ee*NN41)/(3.*cw*Sqrt2) + (CONJRd62*ee
     +*NN42)/(Sqrt2*sw)
      MGVX482(1) = -((CONJUU11*ee*Rd62)/sw)
      MGVX482(2) = 0
      MGVX483(1) = -((CONJUU21*ee*Rd62)/sw)
      MGVX483(2) = 0
      MGVX484(1) = -(CONJNN11*ee*Rd62)/(3.*cw*Sqrt2) + (CONJNN12*ee
     +*Rd62)/(Sqrt2*sw)
      MGVX484(2) = 0
      MGVX485(1) = -(CONJNN21*ee*Rd62)/(3.*cw*Sqrt2) + (CONJNN22*ee
     +*Rd62)/(Sqrt2*sw)
      MGVX485(2) = 0
      endif

      return
      end
