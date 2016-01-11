      subroutine coup9(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 9
      if(readlha) then


      MGVX496(1) = 0
      MGVX496(2) = -((ee*NN31*Rl24*Sqrt2)/cw)
      MGVX497(1) = 0
      MGVX497(2) = -((ee*NN41*Rl24*Sqrt2)/cw)
      MGVX498(1) = 0
      MGVX498(2) = -((ee*NN11*Rl35*Sqrt2)/cw)
      MGVX499(1) = 0
      MGVX499(2) = -((ee*NN21*Rl35*Sqrt2)/cw)
      MGVX500(1) = 0
      MGVX500(2) = -((ee*NN31*Rl35*Sqrt2)/cw)
      MGVX501(1) = 0
      MGVX501(2) = -((ee*NN41*Rl35*Sqrt2)/cw)
      MGVX502(1) = -((CONJUU11*ee*Rl42)/sw)
      MGVX502(2) = 0
      MGVX503(1) = -((CONJUU21*ee*Rl42)/sw)
      MGVX503(2) = 0
      MGVX504(1) = (CONJNN11*ee*Rl42)/(cw*Sqrt2) + (CONJNN12*ee*Rl4
     +2)/(Sqrt2*sw)
      MGVX504(2) = 0
      MGVX505(1) = (CONJNN21*ee*Rl42)/(cw*Sqrt2) + (CONJNN22*ee*Rl4
     +2)/(Sqrt2*sw)
      MGVX505(2) = 0
      MGVX506(1) = (CONJNN31*ee*Rl42)/(cw*Sqrt2) + (CONJNN32*ee*Rl4
     +2)/(Sqrt2*sw)
      MGVX506(2) = 0
      MGVX507(1) = (CONJNN41*ee*Rl42)/(cw*Sqrt2) + (CONJNN42*ee*Rl4
     +2)/(Sqrt2*sw)
      MGVX507(2) = 0
      MGVX508(1) = -((CONJUU11*ee*Rl51)/sw)
      MGVX508(2) = 0
      MGVX509(1) = -((CONJUU21*ee*Rl51)/sw)
      MGVX509(2) = 0
      MGVX510(1) = (CONJNN11*ee*Rl51)/(cw*Sqrt2) + (CONJNN12*ee*Rl5
     +1)/(Sqrt2*sw)
      MGVX510(2) = 0
      endif

      return
      end
