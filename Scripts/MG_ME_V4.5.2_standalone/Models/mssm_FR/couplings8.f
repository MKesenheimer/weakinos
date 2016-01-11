      subroutine coup8(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 8
      if(readlha) then


      MGVX445(1) = 0
      MGVX445(2) = (CONJRLn22*ee*NN31)/(cw*Sqrt2) - (CONJRLn22*ee*N
     +N32)/(Sqrt2*sw)
      MGVX446(1) = 0
      MGVX446(2) = (CONJRLn13*ee*NN31)/(cw*Sqrt2) - (CONJRLn13*ee*N
     +N32)/(Sqrt2*sw)
      MGVX449(1) = 0
      MGVX449(2) = (CONJRl51*ee*NN41)/(cw*Sqrt2) + (CONJRl51*ee*NN4
     +2)/(Sqrt2*sw)
      MGVX450(1) = 0
      MGVX450(2) = (CONJRl42*ee*NN41)/(cw*Sqrt2) + (CONJRl42*ee*NN4
     +2)/(Sqrt2*sw)
      MGVX451(1) = 0
      MGVX451(2) = (CONJRLn31*ee*NN41)/(cw*Sqrt2) - (CONJRLn31*ee*N
     +N42)/(Sqrt2*sw)
      MGVX452(1) = 0
      MGVX452(2) = (CONJRLn22*ee*NN41)/(cw*Sqrt2) - (CONJRLn22*ee*N
     +N42)/(Sqrt2*sw)
      MGVX453(1) = 0
      MGVX453(2) = (CONJRLn13*ee*NN41)/(cw*Sqrt2) - (CONJRLn13*ee*N
     +N42)/(Sqrt2*sw)
      MGVX488(1) = -((CONJUU11*ee*Rl13)/sw) + CONJUU12*Rl16*yukl3
      MGVX488(2) = 0
      MGVX489(1) = -((CONJUU21*ee*Rl13)/sw) + CONJUU22*Rl16*yukl3
      MGVX489(2) = 0
      MGVX490(1) = (CONJNN11*ee*Rl13)/(cw*Sqrt2) + (CONJNN12*ee*Rl1
     +3)/(Sqrt2*sw) - CONJNN13*Rl16*yukl3
      MGVX490(2) = -((ee*NN11*Rl16*Sqrt2)/cw) - NN13*Rl13*yukl3
      MGVX491(1) = (CONJNN21*ee*Rl13)/(cw*Sqrt2) + (CONJNN22*ee*Rl1
     +3)/(Sqrt2*sw) - CONJNN23*Rl16*yukl3
      MGVX491(2) = -((ee*NN21*Rl16*Sqrt2)/cw) - NN23*Rl13*yukl3
      MGVX492(1) = (CONJNN31*ee*Rl13)/(cw*Sqrt2) + (CONJNN32*ee*Rl1
     +3)/(Sqrt2*sw) - CONJNN33*Rl16*yukl3
      MGVX492(2) = -((ee*NN31*Rl16*Sqrt2)/cw) - NN33*Rl13*yukl3
      MGVX493(1) = (CONJNN41*ee*Rl13)/(cw*Sqrt2) + (CONJNN42*ee*Rl1
     +3)/(Sqrt2*sw) - CONJNN43*Rl16*yukl3
      MGVX493(2) = -((ee*NN41*Rl16*Sqrt2)/cw) - NN43*Rl13*yukl3
      MGVX494(1) = 0
      MGVX494(2) = -((ee*NN11*Rl24*Sqrt2)/cw)
      MGVX495(1) = 0
      MGVX495(2) = -((ee*NN21*Rl24*Sqrt2)/cw)
      endif

      return
      end
