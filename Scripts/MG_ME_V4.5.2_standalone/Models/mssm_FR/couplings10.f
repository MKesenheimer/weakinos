      subroutine coup10(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 10
      if(readlha) then


      MGVX511(1) = (CONJNN21*ee*Rl51)/(cw*Sqrt2) + (CONJNN22*ee*Rl5
     +1)/(Sqrt2*sw)
      MGVX511(2) = 0
      MGVX512(1) = (CONJNN31*ee*Rl51)/(cw*Sqrt2) + (CONJNN32*ee*Rl5
     +1)/(Sqrt2*sw)
      MGVX512(2) = 0
      MGVX513(1) = (CONJNN41*ee*Rl51)/(cw*Sqrt2) + (CONJNN42*ee*Rl5
     +1)/(Sqrt2*sw)
      MGVX513(2) = 0
      MGVX514(1) = -((CONJUU11*ee*Rl63)/sw) + CONJUU12*Rl66*yukl3
      MGVX514(2) = 0
      MGVX515(1) = -((CONJUU21*ee*Rl63)/sw) + CONJUU22*Rl66*yukl3
      MGVX515(2) = 0
      MGVX516(1) = (CONJNN11*ee*Rl63)/(cw*Sqrt2) + (CONJNN12*ee*Rl6
     +3)/(Sqrt2*sw) - CONJNN13*Rl66*yukl3
      MGVX516(2) = -((ee*NN11*Rl66*Sqrt2)/cw) - NN13*Rl63*yukl3
      MGVX517(1) = (CONJNN21*ee*Rl63)/(cw*Sqrt2) + (CONJNN22*ee*Rl6
     +3)/(Sqrt2*sw) - CONJNN23*Rl66*yukl3
      MGVX517(2) = -((ee*NN21*Rl66*Sqrt2)/cw) - NN23*Rl63*yukl3
      MGVX518(1) = (CONJNN31*ee*Rl63)/(cw*Sqrt2) + (CONJNN32*ee*Rl6
     +3)/(Sqrt2*sw) - CONJNN33*Rl66*yukl3
      MGVX518(2) = -((ee*NN31*Rl66*Sqrt2)/cw) - NN33*Rl63*yukl3
      MGVX519(1) = (CONJNN41*ee*Rl63)/(cw*Sqrt2) + (CONJNN42*ee*Rl6
     +3)/(Sqrt2*sw) - CONJNN43*Rl66*yukl3
      MGVX519(2) = -((ee*NN41*Rl66*Sqrt2)/cw) - NN43*Rl63*yukl3
      MGVX520(1) = (CONJNN11*ee*RLn13)/(cw*Sqrt2) - (CONJNN12*ee*RL
     +n13)/(Sqrt2*sw)
      MGVX520(2) = 0
      MGVX521(1) = (CONJNN21*ee*RLn13)/(cw*Sqrt2) - (CONJNN22*ee*RL
     +n13)/(Sqrt2*sw)
      MGVX521(2) = 0
      MGVX522(1) = (CONJNN31*ee*RLn13)/(cw*Sqrt2) - (CONJNN32*ee*RL
     +n13)/(Sqrt2*sw)
      MGVX522(2) = 0
      MGVX523(1) = (CONJNN41*ee*RLn13)/(cw*Sqrt2) - (CONJNN42*ee*RL
     +n13)/(Sqrt2*sw)
      MGVX523(2) = 0
      MGVX524(1) = -((CONJVV11*ee*RLn13)/sw)
      MGVX524(2) = RLn13*UU12*yukl3
      MGVX525(1) = -((CONJVV21*ee*RLn13)/sw)
      MGVX525(2) = RLn13*UU22*yukl3
      endif

      return
      end
