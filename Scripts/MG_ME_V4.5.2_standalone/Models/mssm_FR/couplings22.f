      subroutine coup22(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 22
      if(readlha) then


      MGVX582(1) = (2*CONJNN31*CONJRu35*ee*Sqrt2)/(3.*cw)
      MGVX582(2) = 0
      MGVX583(1) = (2*CONJNN41*CONJRu35*ee*Sqrt2)/(3.*cw)
      MGVX583(2) = 0
      MGVX585(1) = 0
      MGVX585(2) = (2*ee*NN11*Ru35*Sqrt2)/(3.*cw)
      MGVX586(1) = 0
      MGVX586(2) = (2*ee*NN21*Ru35*Sqrt2)/(3.*cw)
      MGVX587(1) = 0
      MGVX587(2) = (2*ee*NN31*Ru35*Sqrt2)/(3.*cw)
      MGVX588(1) = 0
      MGVX588(2) = (2*ee*NN41*Ru35*Sqrt2)/(3.*cw)
      MGVX590(1) = 0
      MGVX590(2) = -(CONJRu42*ee*NN11)/(3.*cw*Sqrt2) - (CONJRu42*ee
     +*NN12)/(Sqrt2*sw)
      MGVX591(1) = 0
      MGVX591(2) = -(CONJRu42*ee*NN21)/(3.*cw*Sqrt2) - (CONJRu42*ee
     +*NN22)/(Sqrt2*sw)
      MGVX592(1) = 0
      MGVX592(2) = -(CONJRu42*ee*NN31)/(3.*cw*Sqrt2) - (CONJRu42*ee
     +*NN32)/(Sqrt2*sw)
      MGVX593(1) = 0
      MGVX593(2) = -(CONJRu42*ee*NN41)/(3.*cw*Sqrt2) - (CONJRu42*ee
     +*NN42)/(Sqrt2*sw)
      MGVX646(1) = 0
      MGVX646(2) = -((CONJRu42*ee*VV11)/sw)
      MGVX653(1) = 0
      MGVX653(2) = -((CONJRu42*ee*VV21)/sw)
      MGVX595(1) = -(CONJNN11*ee*Ru42)/(3.*cw*Sqrt2) - (CONJNN12*ee
     +*Ru42)/(Sqrt2*sw)
      MGVX595(2) = 0
      MGVX596(1) = -(CONJNN21*ee*Ru42)/(3.*cw*Sqrt2) - (CONJNN22*ee
     +*Ru42)/(Sqrt2*sw)
      MGVX596(2) = 0
      MGVX597(1) = -(CONJNN31*ee*Ru42)/(3.*cw*Sqrt2) - (CONJNN32*ee
     +*Ru42)/(Sqrt2*sw)
      MGVX597(2) = 0
      endif

      return
      end
