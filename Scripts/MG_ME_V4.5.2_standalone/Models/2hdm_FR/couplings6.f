      subroutine coup6(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 6
      if(readlha) then

      MGVX99 = (ee**2*TH23)/(2.*sw) - ((0,0.5)*ee**2*TH33)/sw
      MGVX121 = -(ee**2*TH23)/(2.*cw) + ((0,0.5)*ee**2*TH33)/cw
      MGVX93 = (ee**2*TH23)/(2.*sw) + ((0,0.5)*ee**2*TH33)/sw
      MGVX118 = -(ee**2*TH23)/(2.*cw) - ((0,0.5)*ee**2*TH33)/cw
      MGVX53 = 2*ee**2
      MGVX111 = (cw*ee**2)/sw - (ee**2*sw)/cw
      MGVX101 = ee**2/(2.*Psw2)
      MGVX122 = -ee**2 + (ee**2*Pcw2)/(2.*Psw2) + (ee**2*Psw2)/(2.*
     +Pcw2)
      MGVX113 = ((0,-0.5)*cw*ee*TH22*TH31)/sw - ((0,0.5)*ee*sw*TH22
     +*TH31)/cw + ((0,0.5)*cw*ee*TH21*TH32)/sw + ((0,0.5)*ee*sw*TH
     +21*TH32)/cw
      MGVX114 = ((0,-0.5)*cw*ee*TH23*TH31)/sw - ((0,0.5)*ee*sw*TH23
     +*TH31)/cw + ((0,0.5)*cw*ee*TH21*TH33)/sw + ((0,0.5)*ee*sw*TH
     +21*TH33)/cw
      MGVX115 = ((0,-0.5)*cw*ee*TH23*TH32)/sw - ((0,0.5)*ee*sw*TH23
     +*TH32)/cw + ((0,0.5)*cw*ee*TH22*TH33)/sw + ((0,0.5)*ee*sw*TH
     +22*TH33)/cw
      MGVX96 = -(ee*TH21)/(2.*sw) + ((0,0.5)*ee*TH31)/sw
      MGVX98 = -(ee*TH22)/(2.*sw) + ((0,0.5)*ee*TH32)/sw
      MGVX100 = -(ee*TH23)/(2.*sw) + ((0,0.5)*ee*TH33)/sw
      MGVX54 = ee
      MGVX90 = (ee*TH21)/(2.*sw) + ((0,0.5)*ee*TH31)/sw
      MGVX92 = (ee*TH22)/(2.*sw) + ((0,0.5)*ee*TH32)/sw
      MGVX94 = (ee*TH23)/(2.*sw) + ((0,0.5)*ee*TH33)/sw
      MGVX112 = (cw*ee)/(2.*sw) - (ee*sw)/(2.*cw)

      endif

      return
      end
