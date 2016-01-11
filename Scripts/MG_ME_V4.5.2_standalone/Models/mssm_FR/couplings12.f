      subroutine coup12(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 12
      if(readlha) then


      MGVX541(1) = (CONJNN41*CONJUU12*ee*Sinbet)/(cw*Sqrt2) - (CONJ
     +NN43*CONJUU11*ee*Sinbet)/sw + (CONJNN42*CONJUU12*ee*Sinbet)/
     +(Sqrt2*sw)
      MGVX541(2) = -((Cosbet*ee*NN44*VV11)/sw) - (Cosbet*ee*NN41*VV
     +12)/(cw*Sqrt2) - (Cosbet*ee*NN42*VV12)/(Sqrt2*sw)
      MGVX542(1) = (CONJNN11*CONJUU22*ee*Sinbet)/(cw*Sqrt2) - (CONJ
     +NN13*CONJUU21*ee*Sinbet)/sw + (CONJNN12*CONJUU22*ee*Sinbet)/
     +(Sqrt2*sw)
      MGVX542(2) = -((Cosbet*ee*NN14*VV21)/sw) - (Cosbet*ee*NN11*VV
     +22)/(cw*Sqrt2) - (Cosbet*ee*NN12*VV22)/(Sqrt2*sw)
      MGVX543(1) = (CONJNN21*CONJUU22*ee*Sinbet)/(cw*Sqrt2) - (CONJ
     +NN23*CONJUU21*ee*Sinbet)/sw + (CONJNN22*CONJUU22*ee*Sinbet)/
     +(Sqrt2*sw)
      MGVX543(2) = -((Cosbet*ee*NN24*VV21)/sw) - (Cosbet*ee*NN21*VV
     +22)/(cw*Sqrt2) - (Cosbet*ee*NN22*VV22)/(Sqrt2*sw)
      MGVX544(1) = (CONJNN31*CONJUU22*ee*Sinbet)/(cw*Sqrt2) - (CONJ
     +NN33*CONJUU21*ee*Sinbet)/sw + (CONJNN32*CONJUU22*ee*Sinbet)/
     +(Sqrt2*sw)
      MGVX544(2) = -((Cosbet*ee*NN34*VV21)/sw) - (Cosbet*ee*NN31*VV
     +22)/(cw*Sqrt2) - (Cosbet*ee*NN32*VV22)/(Sqrt2*sw)
      MGVX545(1) = (CONJNN41*CONJUU22*ee*Sinbet)/(cw*Sqrt2) - (CONJ
     +NN43*CONJUU21*ee*Sinbet)/sw + (CONJNN42*CONJUU22*ee*Sinbet)/
     +(Sqrt2*sw)
      MGVX545(2) = -((Cosbet*ee*NN44*VV21)/sw) - (Cosbet*ee*NN41*VV
     +22)/(cw*Sqrt2) - (Cosbet*ee*NN42*VV22)/(Sqrt2*sw)
      MGVX626(1) = 0
      MGVX626(2) = -((CONJRl51*ee*UU11)/sw)
      MGVX627(1) = 0
      MGVX627(2) = -((CONJRl42*ee*UU11)/sw)
      MGVX628(1) = 0
      MGVX628(2) = -((CONJRl13*ee*UU11)/sw) + CONJRl16*UU12*yukl3
      MGVX629(1) = 0
      MGVX629(2) = -((CONJRl63*ee*UU11)/sw) + CONJRl66*UU12*yukl3
      MGVX634(1) = 0
      MGVX634(2) = -((CONJRl51*ee*UU21)/sw)
      MGVX635(1) = 0
      MGVX635(2) = -((CONJRl42*ee*UU21)/sw)
      MGVX636(1) = 0
      MGVX636(2) = -((CONJRl13*ee*UU21)/sw) + CONJRl16*UU22*yukl3
      MGVX637(1) = 0
      MGVX637(2) = -((CONJRl63*ee*UU21)/sw) + CONJRl66*UU22*yukl3
      MGVX642(1) = 0
      MGVX642(2) = -((CONJRLn31*ee*VV11)/sw)
      MGVX643(1) = 0
      MGVX643(2) = -((CONJRLn22*ee*VV11)/sw)
      endif

      return
      end
