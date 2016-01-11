      subroutine coup16(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 16
      if(readlha) then

      MGVX54 = ee
      MGVX90 = (ee*TH21)/(2.*sw) + ((0,0.5)*ee*TH31)/sw
      MGVX92 = (ee*TH22)/(2.*sw) + ((0,0.5)*ee*TH32)/sw
      MGVX94 = (ee*TH23)/(2.*sw) + ((0,0.5)*ee*TH33)/sw
      MGVX112 = (cw*ee)/(2.*sw) - (ee*sw)/(2.*cw)

      endif

      return
      end
