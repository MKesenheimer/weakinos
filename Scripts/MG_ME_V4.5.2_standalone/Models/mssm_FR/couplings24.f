      subroutine coup24(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 24
      if(readlha) then

      MGVX768 = -(Cosalp*Cosbet*ee**2*Sinbet*vd)/(4.*Pcw2*Psw2) - (
     +ee**2*Sinalp*vd*Cos(bet)**2)/(4.*Pcw2*Psw2) + (ee**2*Sinalp*
     +vd*Sin(bet)**2)/(4.*Pcw2*Psw2) + (Cosalp*ee**2*Tanbet*vd*Sin
     +(bet)**2)/(4.*Pcw2*Psw2)

      MGVX648(1) = CONJRu63*CONJUU12*yukd3
      MGVX648(2) = -((CONJRu63*ee*VV11)/sw) + CONJRu66*VV12*yuku3
      MGVX655(1) = CONJRu63*CONJUU22*yukd3
      MGVX655(2) = -((CONJRu63*ee*VV21)/sw) + CONJRu66*VV22*yuku3
      MGVX614(1) = (2*CONJNN11*CONJRu66*ee*Sqrt2)/(3.*cw) - CONJNN1
     +4*CONJRu63*yuku3
      MGVX614(2) = -(CONJRu63*ee*NN11)/(3.*cw*Sqrt2) - (CONJRu63*ee
     +*NN12)/(Sqrt2*sw) - CONJRu66*NN14*yuku3
      MGVX615(1) = (2*CONJNN21*CONJRu66*ee*Sqrt2)/(3.*cw) - CONJNN2
     +4*CONJRu63*yuku3
      MGVX615(2) = -(CONJRu63*ee*NN21)/(3.*cw*Sqrt2) - (CONJRu63*ee
     +*NN22)/(Sqrt2*sw) - CONJRu66*NN24*yuku3
      MGVX616(1) = (2*CONJNN31*CONJRu66*ee*Sqrt2)/(3.*cw) - CONJNN3
     +4*CONJRu63*yuku3
      MGVX616(2) = -(CONJRu63*ee*NN31)/(3.*cw*Sqrt2) - (CONJRu63*ee
     +*NN32)/(Sqrt2*sw) - CONJRu66*NN34*yuku3
      MGVX617(1) = (2*CONJNN41*CONJRu66*ee*Sqrt2)/(3.*cw) - CONJNN4
     +4*CONJRu63*yuku3
      MGVX617(2) = -(CONJRu63*ee*NN41)/(3.*cw*Sqrt2) - (CONJRu63*ee
     +*NN42)/(Sqrt2*sw) - CONJRu66*NN44*yuku3
      MGVX621(1) = -(CONJNN11*ee*Ru63)/(3.*cw*Sqrt2) - (CONJNN12*ee
     +*Ru63)/(Sqrt2*sw) - CONJNN14*Ru66*yuku3
      MGVX621(2) = (2*ee*NN11*Ru66*Sqrt2)/(3.*cw) - NN14*Ru63*yuku3
      MGVX622(1) = -(CONJNN21*ee*Ru63)/(3.*cw*Sqrt2) - (CONJNN22*ee
     +*Ru63)/(Sqrt2*sw) - CONJNN24*Ru66*yuku3
      MGVX622(2) = (2*ee*NN21*Ru66*Sqrt2)/(3.*cw) - NN24*Ru63*yuku3
      MGVX623(1) = -(CONJNN31*ee*Ru63)/(3.*cw*Sqrt2) - (CONJNN32*ee
     +*Ru63)/(Sqrt2*sw) - CONJNN34*Ru66*yuku3
      MGVX623(2) = (2*ee*NN31*Ru66*Sqrt2)/(3.*cw) - NN34*Ru63*yuku3
      MGVX624(1) = -(CONJNN41*ee*Ru63)/(3.*cw*Sqrt2) - (CONJNN42*ee
     +*Ru63)/(Sqrt2*sw) - CONJNN44*Ru66*yuku3
      MGVX624(2) = (2*ee*NN41*Ru66*Sqrt2)/(3.*cw) - NN44*Ru63*yuku3
      MGVX619(1) = -((CONJVV11*ee*Ru63)/sw) + CONJVV12*Ru66*yuku3
      MGVX619(2) = Ru63*UU12*yukd3
      MGVX620(1) = -((CONJVV21*ee*Ru63)/sw) + CONJVV22*Ru66*yuku3
      MGVX620(2) = Ru63*UU22*yukd3
      MGVX36(1) = ee/(2.*cw*sw) - (ee*sw)/(3.*cw)
      MGVX36(2) = -(ee*sw)/(3.*cw)
      MGVX42(1) = -ee/(2.*cw*sw) + (2*ee*sw)/(3.*cw)
      MGVX42(2) = (2*ee*sw)/(3.*cw)
      endif

      return
      end
