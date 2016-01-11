      subroutine coup3(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 3
      if(readlha) then


      MGVX350(1) = -((CONJNN11*CONJRl24*ee*Sqrt2)/cw)
      MGVX350(2) = 0
      MGVX351(1) = -((CONJNN21*CONJRl24*ee*Sqrt2)/cw)
      MGVX351(2) = 0
      MGVX352(1) = -((CONJNN31*CONJRl24*ee*Sqrt2)/cw)
      MGVX352(2) = 0
      MGVX353(1) = -((CONJNN41*CONJRl24*ee*Sqrt2)/cw)
      MGVX353(2) = 0
      MGVX354(1) = -((CONJNN11*CONJRl35*ee*Sqrt2)/cw)
      MGVX354(2) = 0
      MGVX355(1) = -((CONJNN21*CONJRl35*ee*Sqrt2)/cw)
      MGVX355(2) = 0
      MGVX356(1) = -((CONJNN31*CONJRl35*ee*Sqrt2)/cw)
      MGVX356(2) = 0
      MGVX357(1) = -((CONJNN41*CONJRl35*ee*Sqrt2)/cw)
      MGVX357(2) = 0
      MGVX358(1) = -((CONJNN11*CONJVV12*Cosbet*ee)/(cw*Sqrt2)) - (C
     +ONJNN14*CONJVV11*Cosbet*ee)/sw - (CONJNN12*CONJVV12*Cosbet*e
     +e)/(Sqrt2*sw)
      MGVX358(2) = -((ee*NN13*Sinbet*UU11)/sw) + (ee*NN11*Sinbet*UU
     +12)/(cw*Sqrt2) + (ee*NN12*Sinbet*UU12)/(Sqrt2*sw)
      MGVX359(1) = -((CONJNN11*CONJVV22*Cosbet*ee)/(cw*Sqrt2)) - (C
     +ONJNN14*CONJVV21*Cosbet*ee)/sw - (CONJNN12*CONJVV22*Cosbet*e
     +e)/(Sqrt2*sw)
      MGVX359(2) = -((ee*NN13*Sinbet*UU21)/sw) + (ee*NN11*Sinbet*UU
     +22)/(cw*Sqrt2) + (ee*NN12*Sinbet*UU22)/(Sqrt2*sw)
      MGVX360(1) = (CONJNN11*CONJNN13*Cosalp*ee)/cw - (CONJNN11*CON
     +JNN14*ee*Sinalp)/cw - (CONJNN12*CONJNN13*Cosalp*ee)/sw + (CO
     +NJNN12*CONJNN14*ee*Sinalp)/sw
      MGVX360(2) = (Cosalp*ee*NN11*NN13)/cw - (ee*NN11*NN14*Sinalp)
     +/cw - (Cosalp*ee*NN12*NN13)/sw + (ee*NN12*NN14*Sinalp)/sw
      MGVX361(1) = -((CONJNN11*CONJNN14*Cosalp*ee)/cw) - (CONJNN11*
     +CONJNN13*ee*Sinalp)/cw + (CONJNN12*CONJNN14*Cosalp*ee)/sw + 
     +(CONJNN12*CONJNN13*ee*Sinalp)/sw
      MGVX361(2) = -((Cosalp*ee*NN11*NN14)/cw) - (ee*NN11*NN13*Sina
     +lp)/cw + (Cosalp*ee*NN12*NN14)/sw + (ee*NN12*NN13*Sinalp)/sw
      MGVX362(1) = ((0,1)*CONJNN11*CONJNN14*Cosbet*ee)/cw - ((0,1)*
     +CONJNN11*CONJNN13*ee*Sinbet)/cw - ((0,1)*CONJNN12*CONJNN14*C
     +osbet*ee)/sw + ((0,1)*CONJNN12*CONJNN13*ee*Sinbet)/sw
      MGVX362(2) = ((0,-1)*Cosbet*ee*NN11*NN14)/cw + ((0,1)*ee*NN11
     +*NN13*Sinbet)/cw + ((0,1)*Cosbet*ee*NN12*NN14)/sw - ((0,1)*e
     +e*NN12*NN13*Sinbet)/sw
      MGVX363(1) = (CONJNN13*CONJNN21*Cosalp*ee)/(2.*cw) + (CONJNN1
     +1*CONJNN23*Cosalp*ee)/(2.*cw) - (CONJNN14*CONJNN21*ee*Sinalp
     +)/(2.*cw) - (CONJNN11*CONJNN24*ee*Sinalp)/(2.*cw) - (CONJNN1
     +3*CONJNN22*Cosalp*ee)/(2.*sw) - (CONJNN12*CONJNN23*Cosalp*ee
     +)/(2.*sw) + (CONJNN14*CONJNN22*ee*Sinalp)/(2.*sw) + (CONJNN1
     +2*CONJNN24*ee*Sinalp)/(2.*sw)
      MGVX363(2) = (Cosalp*ee*NN13*NN21)/(2.*cw) + (Cosalp*ee*NN11*
     +NN23)/(2.*cw) - (ee*NN14*NN21*Sinalp)/(2.*cw) - (ee*NN11*NN2
     +4*Sinalp)/(2.*cw) - (Cosalp*ee*NN13*NN22)/(2.*sw) - (Cosalp*
     +ee*NN12*NN23)/(2.*sw) + (ee*NN14*NN22*Sinalp)/(2.*sw) + (ee*
     +NN12*NN24*Sinalp)/(2.*sw)
      MGVX364(1) = -(CONJNN14*CONJNN21*Cosalp*ee)/(2.*cw) - (CONJNN
     +11*CONJNN24*Cosalp*ee)/(2.*cw) - (CONJNN13*CONJNN21*ee*Sinal
     +p)/(2.*cw) - (CONJNN11*CONJNN23*ee*Sinalp)/(2.*cw) + (CONJNN
     +14*CONJNN22*Cosalp*ee)/(2.*sw) + (CONJNN12*CONJNN24*Cosalp*e
     +e)/(2.*sw) + (CONJNN13*CONJNN22*ee*Sinalp)/(2.*sw) + (CONJNN
     +12*CONJNN23*ee*Sinalp)/(2.*sw)
      MGVX364(2) = -(Cosalp*ee*NN14*NN21)/(2.*cw) - (Cosalp*ee*NN11
     +*NN24)/(2.*cw) - (ee*NN13*NN21*Sinalp)/(2.*cw) - (ee*NN11*NN
     +23*Sinalp)/(2.*cw) + (Cosalp*ee*NN14*NN22)/(2.*sw) + (Cosalp
     +*ee*NN12*NN24)/(2.*sw) + (ee*NN13*NN22*Sinalp)/(2.*sw) + (ee
     +*NN12*NN23*Sinalp)/(2.*sw)
      endif

      return
      end
