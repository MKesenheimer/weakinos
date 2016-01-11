      subroutine coup11(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 11
      if(readlha) then


      MGVX526(1) = (CONJNN11*ee*RLn22)/(cw*Sqrt2) - (CONJNN12*ee*RL
     +n22)/(Sqrt2*sw)
      MGVX526(2) = 0
      MGVX527(1) = (CONJNN21*ee*RLn22)/(cw*Sqrt2) - (CONJNN22*ee*RL
     +n22)/(Sqrt2*sw)
      MGVX527(2) = 0
      MGVX528(1) = (CONJNN31*ee*RLn22)/(cw*Sqrt2) - (CONJNN32*ee*RL
     +n22)/(Sqrt2*sw)
      MGVX528(2) = 0
      MGVX529(1) = (CONJNN41*ee*RLn22)/(cw*Sqrt2) - (CONJNN42*ee*RL
     +n22)/(Sqrt2*sw)
      MGVX529(2) = 0
      MGVX530(1) = -((CONJVV11*ee*RLn22)/sw)
      MGVX530(2) = 0
      MGVX531(1) = -((CONJVV21*ee*RLn22)/sw)
      MGVX531(2) = 0
      MGVX532(1) = (CONJNN11*ee*RLn31)/(cw*Sqrt2) - (CONJNN12*ee*RL
     +n31)/(Sqrt2*sw)
      MGVX532(2) = 0
      MGVX533(1) = (CONJNN21*ee*RLn31)/(cw*Sqrt2) - (CONJNN22*ee*RL
     +n31)/(Sqrt2*sw)
      MGVX533(2) = 0
      MGVX534(1) = (CONJNN31*ee*RLn31)/(cw*Sqrt2) - (CONJNN32*ee*RL
     +n31)/(Sqrt2*sw)
      MGVX534(2) = 0
      MGVX535(1) = (CONJNN41*ee*RLn31)/(cw*Sqrt2) - (CONJNN42*ee*RL
     +n31)/(Sqrt2*sw)
      MGVX535(2) = 0
      MGVX536(1) = -((CONJVV11*ee*RLn31)/sw)
      MGVX536(2) = 0
      MGVX537(1) = -((CONJVV21*ee*RLn31)/sw)
      MGVX537(2) = 0
      MGVX538(1) = (CONJNN11*CONJUU12*ee*Sinbet)/(cw*Sqrt2) - (CONJ
     +NN13*CONJUU11*ee*Sinbet)/sw + (CONJNN12*CONJUU12*ee*Sinbet)/
     +(Sqrt2*sw)
      MGVX538(2) = -((Cosbet*ee*NN14*VV11)/sw) - (Cosbet*ee*NN11*VV
     +12)/(cw*Sqrt2) - (Cosbet*ee*NN12*VV12)/(Sqrt2*sw)
      MGVX539(1) = (CONJNN21*CONJUU12*ee*Sinbet)/(cw*Sqrt2) - (CONJ
     +NN23*CONJUU11*ee*Sinbet)/sw + (CONJNN22*CONJUU12*ee*Sinbet)/
     +(Sqrt2*sw)
      MGVX539(2) = -((Cosbet*ee*NN24*VV11)/sw) - (Cosbet*ee*NN21*VV
     +12)/(cw*Sqrt2) - (Cosbet*ee*NN22*VV12)/(Sqrt2*sw)
      MGVX540(1) = (CONJNN31*CONJUU12*ee*Sinbet)/(cw*Sqrt2) - (CONJ
     +NN33*CONJUU11*ee*Sinbet)/sw + (CONJNN32*CONJUU12*ee*Sinbet)/
     +(Sqrt2*sw)
      MGVX540(2) = -((Cosbet*ee*NN34*VV11)/sw) - (Cosbet*ee*NN31*VV
     +12)/(cw*Sqrt2) - (Cosbet*ee*NN32*VV12)/(Sqrt2*sw)
      endif

      return
      end
