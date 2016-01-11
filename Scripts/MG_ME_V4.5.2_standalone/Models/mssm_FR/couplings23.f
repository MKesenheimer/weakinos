      subroutine coup23(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 23
      if(readlha) then


      MGVX598(1) = -(CONJNN41*ee*Ru42)/(3.*cw*Sqrt2) - (CONJNN42*ee
     +*Ru42)/(Sqrt2*sw)
      MGVX598(2) = 0
      MGVX599(1) = -((CONJVV11*ee*Ru42)/sw)
      MGVX599(2) = 0
      MGVX600(1) = -((CONJVV21*ee*Ru42)/sw)
      MGVX600(2) = 0
      MGVX647(1) = 0
      MGVX647(2) = -((CONJRu51*ee*VV11)/sw)
      MGVX654(1) = 0
      MGVX654(2) = -((CONJRu51*ee*VV21)/sw)
      MGVX602(1) = 0
      MGVX602(2) = -(CONJRu51*ee*NN11)/(3.*cw*Sqrt2) - (CONJRu51*ee
     +*NN12)/(Sqrt2*sw)
      MGVX603(1) = 0
      MGVX603(2) = -(CONJRu51*ee*NN21)/(3.*cw*Sqrt2) - (CONJRu51*ee
     +*NN22)/(Sqrt2*sw)
      MGVX604(1) = 0
      MGVX604(2) = -(CONJRu51*ee*NN31)/(3.*cw*Sqrt2) - (CONJRu51*ee
     +*NN32)/(Sqrt2*sw)
      MGVX605(1) = 0
      MGVX605(2) = -(CONJRu51*ee*NN41)/(3.*cw*Sqrt2) - (CONJRu51*ee
     +*NN42)/(Sqrt2*sw)
      MGVX607(1) = -(CONJNN11*ee*Ru51)/(3.*cw*Sqrt2) - (CONJNN12*ee
     +*Ru51)/(Sqrt2*sw)
      MGVX607(2) = 0
      MGVX608(1) = -(CONJNN21*ee*Ru51)/(3.*cw*Sqrt2) - (CONJNN22*ee
     +*Ru51)/(Sqrt2*sw)
      MGVX608(2) = 0
      MGVX609(1) = -(CONJNN31*ee*Ru51)/(3.*cw*Sqrt2) - (CONJNN32*ee
     +*Ru51)/(Sqrt2*sw)
      MGVX609(2) = 0
      MGVX610(1) = -(CONJNN41*ee*Ru51)/(3.*cw*Sqrt2) - (CONJNN42*ee
     +*Ru51)/(Sqrt2*sw)
      MGVX610(2) = 0
      MGVX611(1) = -((CONJVV11*ee*Ru51)/sw)
      MGVX611(2) = 0
      MGVX612(1) = -((CONJVV21*ee*Ru51)/sw)
      MGVX612(2) = 0
      endif

      return
      end
