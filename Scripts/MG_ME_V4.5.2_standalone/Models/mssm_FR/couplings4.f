      subroutine coup4(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 4
      if(readlha) then


      MGVX365(1) = ((0,0.5)*CONJNN14*CONJNN21*Cosbet*ee)/cw + ((0,0
     +.5)*CONJNN11*CONJNN24*Cosbet*ee)/cw - ((0,0.5)*CONJNN13*CONJ
     +NN21*ee*Sinbet)/cw - ((0,0.5)*CONJNN11*CONJNN23*ee*Sinbet)/c
     +w - ((0,0.5)*CONJNN14*CONJNN22*Cosbet*ee)/sw - ((0,0.5)*CONJ
     +NN12*CONJNN24*Cosbet*ee)/sw + ((0,0.5)*CONJNN13*CONJNN22*ee*
     +Sinbet)/sw + ((0,0.5)*CONJNN12*CONJNN23*ee*Sinbet)/sw
      MGVX365(2) = ((0,-0.5)*Cosbet*ee*NN14*NN21)/cw - ((0,0.5)*Cos
     +bet*ee*NN11*NN24)/cw + ((0,0.5)*ee*NN13*NN21*Sinbet)/cw + ((
     +0,0.5)*ee*NN11*NN23*Sinbet)/cw + ((0,0.5)*Cosbet*ee*NN14*NN2
     +2)/sw + ((0,0.5)*Cosbet*ee*NN12*NN24)/sw - ((0,0.5)*ee*NN13*
     +NN22*Sinbet)/sw - ((0,0.5)*ee*NN12*NN23*Sinbet)/sw
      MGVX366(1) = (CONJNN13*CONJNN31*Cosalp*ee)/(2.*cw) + (CONJNN1
     +1*CONJNN33*Cosalp*ee)/(2.*cw) - (CONJNN14*CONJNN31*ee*Sinalp
     +)/(2.*cw) - (CONJNN11*CONJNN34*ee*Sinalp)/(2.*cw) - (CONJNN1
     +3*CONJNN32*Cosalp*ee)/(2.*sw) - (CONJNN12*CONJNN33*Cosalp*ee
     +)/(2.*sw) + (CONJNN14*CONJNN32*ee*Sinalp)/(2.*sw) + (CONJNN1
     +2*CONJNN34*ee*Sinalp)/(2.*sw)
      MGVX366(2) = (Cosalp*ee*NN13*NN31)/(2.*cw) + (Cosalp*ee*NN11*
     +NN33)/(2.*cw) - (ee*NN14*NN31*Sinalp)/(2.*cw) - (ee*NN11*NN3
     +4*Sinalp)/(2.*cw) - (Cosalp*ee*NN13*NN32)/(2.*sw) - (Cosalp*
     +ee*NN12*NN33)/(2.*sw) + (ee*NN14*NN32*Sinalp)/(2.*sw) + (ee*
     +NN12*NN34*Sinalp)/(2.*sw)
      MGVX367(1) = -(CONJNN14*CONJNN31*Cosalp*ee)/(2.*cw) - (CONJNN
     +11*CONJNN34*Cosalp*ee)/(2.*cw) - (CONJNN13*CONJNN31*ee*Sinal
     +p)/(2.*cw) - (CONJNN11*CONJNN33*ee*Sinalp)/(2.*cw) + (CONJNN
     +14*CONJNN32*Cosalp*ee)/(2.*sw) + (CONJNN12*CONJNN34*Cosalp*e
     +e)/(2.*sw) + (CONJNN13*CONJNN32*ee*Sinalp)/(2.*sw) + (CONJNN
     +12*CONJNN33*ee*Sinalp)/(2.*sw)
      MGVX367(2) = -(Cosalp*ee*NN14*NN31)/(2.*cw) - (Cosalp*ee*NN11
     +*NN34)/(2.*cw) - (ee*NN13*NN31*Sinalp)/(2.*cw) - (ee*NN11*NN
     +33*Sinalp)/(2.*cw) + (Cosalp*ee*NN14*NN32)/(2.*sw) + (Cosalp
     +*ee*NN12*NN34)/(2.*sw) + (ee*NN13*NN32*Sinalp)/(2.*sw) + (ee
     +*NN12*NN33*Sinalp)/(2.*sw)
      MGVX368(1) = ((0,0.5)*CONJNN14*CONJNN31*Cosbet*ee)/cw + ((0,0
     +.5)*CONJNN11*CONJNN34*Cosbet*ee)/cw - ((0,0.5)*CONJNN13*CONJ
     +NN31*ee*Sinbet)/cw - ((0,0.5)*CONJNN11*CONJNN33*ee*Sinbet)/c
     +w - ((0,0.5)*CONJNN14*CONJNN32*Cosbet*ee)/sw - ((0,0.5)*CONJ
     +NN12*CONJNN34*Cosbet*ee)/sw + ((0,0.5)*CONJNN13*CONJNN32*ee*
     +Sinbet)/sw + ((0,0.5)*CONJNN12*CONJNN33*ee*Sinbet)/sw
      MGVX368(2) = ((0,-0.5)*Cosbet*ee*NN14*NN31)/cw - ((0,0.5)*Cos
     +bet*ee*NN11*NN34)/cw + ((0,0.5)*ee*NN13*NN31*Sinbet)/cw + ((
     +0,0.5)*ee*NN11*NN33*Sinbet)/cw + ((0,0.5)*Cosbet*ee*NN14*NN3
     +2)/sw + ((0,0.5)*Cosbet*ee*NN12*NN34)/sw - ((0,0.5)*ee*NN13*
     +NN32*Sinbet)/sw - ((0,0.5)*ee*NN12*NN33*Sinbet)/sw
      MGVX369(1) = (CONJNN13*CONJNN41*Cosalp*ee)/(2.*cw) + (CONJNN1
     +1*CONJNN43*Cosalp*ee)/(2.*cw) - (CONJNN14*CONJNN41*ee*Sinalp
     +)/(2.*cw) - (CONJNN11*CONJNN44*ee*Sinalp)/(2.*cw) - (CONJNN1
     +3*CONJNN42*Cosalp*ee)/(2.*sw) - (CONJNN12*CONJNN43*Cosalp*ee
     +)/(2.*sw) + (CONJNN14*CONJNN42*ee*Sinalp)/(2.*sw) + (CONJNN1
     +2*CONJNN44*ee*Sinalp)/(2.*sw)
      MGVX369(2) = (Cosalp*ee*NN13*NN41)/(2.*cw) + (Cosalp*ee*NN11*
     +NN43)/(2.*cw) - (ee*NN14*NN41*Sinalp)/(2.*cw) - (ee*NN11*NN4
     +4*Sinalp)/(2.*cw) - (Cosalp*ee*NN13*NN42)/(2.*sw) - (Cosalp*
     +ee*NN12*NN43)/(2.*sw) + (ee*NN14*NN42*Sinalp)/(2.*sw) + (ee*
     +NN12*NN44*Sinalp)/(2.*sw)
      MGVX370(1) = -(CONJNN14*CONJNN41*Cosalp*ee)/(2.*cw) - (CONJNN
     +11*CONJNN44*Cosalp*ee)/(2.*cw) - (CONJNN13*CONJNN41*ee*Sinal
     +p)/(2.*cw) - (CONJNN11*CONJNN43*ee*Sinalp)/(2.*cw) + (CONJNN
     +14*CONJNN42*Cosalp*ee)/(2.*sw) + (CONJNN12*CONJNN44*Cosalp*e
     +e)/(2.*sw) + (CONJNN13*CONJNN42*ee*Sinalp)/(2.*sw) + (CONJNN
     +12*CONJNN43*ee*Sinalp)/(2.*sw)
      MGVX370(2) = -(Cosalp*ee*NN14*NN41)/(2.*cw) - (Cosalp*ee*NN11
     +*NN44)/(2.*cw) - (ee*NN13*NN41*Sinalp)/(2.*cw) - (ee*NN11*NN
     +43*Sinalp)/(2.*cw) + (Cosalp*ee*NN14*NN42)/(2.*sw) + (Cosalp
     +*ee*NN12*NN44)/(2.*sw) + (ee*NN13*NN42*Sinalp)/(2.*sw) + (ee
     +*NN12*NN43*Sinalp)/(2.*sw)
      MGVX371(1) = ((0,0.5)*CONJNN14*CONJNN41*Cosbet*ee)/cw + ((0,0
     +.5)*CONJNN11*CONJNN44*Cosbet*ee)/cw - ((0,0.5)*CONJNN13*CONJ
     +NN41*ee*Sinbet)/cw - ((0,0.5)*CONJNN11*CONJNN43*ee*Sinbet)/c
     +w - ((0,0.5)*CONJNN14*CONJNN42*Cosbet*ee)/sw - ((0,0.5)*CONJ
     +NN12*CONJNN44*Cosbet*ee)/sw + ((0,0.5)*CONJNN13*CONJNN42*ee*
     +Sinbet)/sw + ((0,0.5)*CONJNN12*CONJNN43*ee*Sinbet)/sw
      MGVX371(2) = ((0,-0.5)*Cosbet*ee*NN14*NN41)/cw - ((0,0.5)*Cos
     +bet*ee*NN11*NN44)/cw + ((0,0.5)*ee*NN13*NN41*Sinbet)/cw + ((
     +0,0.5)*ee*NN11*NN43*Sinbet)/cw + ((0,0.5)*Cosbet*ee*NN14*NN4
     +2)/sw + ((0,0.5)*Cosbet*ee*NN12*NN44)/sw - ((0,0.5)*ee*NN13*
     +NN42*Sinbet)/sw - ((0,0.5)*ee*NN12*NN43*Sinbet)/sw
      MGVX372(1) = -((CONJNN21*CONJVV12*Cosbet*ee)/(cw*Sqrt2)) - (C
     +ONJNN24*CONJVV11*Cosbet*ee)/sw - (CONJNN22*CONJVV12*Cosbet*e
     +e)/(Sqrt2*sw)
      MGVX372(2) = -((ee*NN23*Sinbet*UU11)/sw) + (ee*NN21*Sinbet*UU
     +12)/(cw*Sqrt2) + (ee*NN22*Sinbet*UU12)/(Sqrt2*sw)
      MGVX373(1) = -((CONJNN21*CONJVV22*Cosbet*ee)/(cw*Sqrt2)) - (C
     +ONJNN24*CONJVV21*Cosbet*ee)/sw - (CONJNN22*CONJVV22*Cosbet*e
     +e)/(Sqrt2*sw)
      MGVX373(2) = -((ee*NN23*Sinbet*UU21)/sw) + (ee*NN21*Sinbet*UU
     +22)/(cw*Sqrt2) + (ee*NN22*Sinbet*UU22)/(Sqrt2*sw)
      MGVX374(1) = (CONJNN21*CONJNN23*Cosalp*ee)/cw - (CONJNN21*CON
     +JNN24*ee*Sinalp)/cw - (CONJNN22*CONJNN23*Cosalp*ee)/sw + (CO
     +NJNN22*CONJNN24*ee*Sinalp)/sw
      MGVX374(2) = (Cosalp*ee*NN21*NN23)/cw - (ee*NN21*NN24*Sinalp)
     +/cw - (Cosalp*ee*NN22*NN23)/sw + (ee*NN22*NN24*Sinalp)/sw
      MGVX375(1) = -((CONJNN21*CONJNN24*Cosalp*ee)/cw) - (CONJNN21*
     +CONJNN23*ee*Sinalp)/cw + (CONJNN22*CONJNN24*Cosalp*ee)/sw + 
     +(CONJNN22*CONJNN23*ee*Sinalp)/sw
      MGVX375(2) = -((Cosalp*ee*NN21*NN24)/cw) - (ee*NN21*NN23*Sina
     +lp)/cw + (Cosalp*ee*NN22*NN24)/sw + (ee*NN22*NN23*Sinalp)/sw
      MGVX376(1) = ((0,1)*CONJNN21*CONJNN24*Cosbet*ee)/cw - ((0,1)*
     +CONJNN21*CONJNN23*ee*Sinbet)/cw - ((0,1)*CONJNN22*CONJNN24*C
     +osbet*ee)/sw + ((0,1)*CONJNN22*CONJNN23*ee*Sinbet)/sw
      MGVX376(2) = ((0,-1)*Cosbet*ee*NN21*NN24)/cw + ((0,1)*ee*NN21
     +*NN23*Sinbet)/cw + ((0,1)*Cosbet*ee*NN22*NN24)/sw - ((0,1)*e
     +e*NN22*NN23*Sinbet)/sw
      MGVX377(1) = (CONJNN23*CONJNN31*Cosalp*ee)/(2.*cw) + (CONJNN2
     +1*CONJNN33*Cosalp*ee)/(2.*cw) - (CONJNN24*CONJNN31*ee*Sinalp
     +)/(2.*cw) - (CONJNN21*CONJNN34*ee*Sinalp)/(2.*cw) - (CONJNN2
     +3*CONJNN32*Cosalp*ee)/(2.*sw) - (CONJNN22*CONJNN33*Cosalp*ee
     +)/(2.*sw) + (CONJNN24*CONJNN32*ee*Sinalp)/(2.*sw) + (CONJNN2
     +2*CONJNN34*ee*Sinalp)/(2.*sw)
      MGVX377(2) = (Cosalp*ee*NN23*NN31)/(2.*cw) + (Cosalp*ee*NN21*
     +NN33)/(2.*cw) - (ee*NN24*NN31*Sinalp)/(2.*cw) - (ee*NN21*NN3
     +4*Sinalp)/(2.*cw) - (Cosalp*ee*NN23*NN32)/(2.*sw) - (Cosalp*
     +ee*NN22*NN33)/(2.*sw) + (ee*NN24*NN32*Sinalp)/(2.*sw) + (ee*
     +NN22*NN34*Sinalp)/(2.*sw)
      MGVX378(1) = -(CONJNN24*CONJNN31*Cosalp*ee)/(2.*cw) - (CONJNN
     +21*CONJNN34*Cosalp*ee)/(2.*cw) - (CONJNN23*CONJNN31*ee*Sinal
     +p)/(2.*cw) - (CONJNN21*CONJNN33*ee*Sinalp)/(2.*cw) + (CONJNN
     +24*CONJNN32*Cosalp*ee)/(2.*sw) + (CONJNN22*CONJNN34*Cosalp*e
     +e)/(2.*sw) + (CONJNN23*CONJNN32*ee*Sinalp)/(2.*sw) + (CONJNN
     +22*CONJNN33*ee*Sinalp)/(2.*sw)
      MGVX378(2) = -(Cosalp*ee*NN24*NN31)/(2.*cw) - (Cosalp*ee*NN21
     +*NN34)/(2.*cw) - (ee*NN23*NN31*Sinalp)/(2.*cw) - (ee*NN21*NN
     +33*Sinalp)/(2.*cw) + (Cosalp*ee*NN24*NN32)/(2.*sw) + (Cosalp
     +*ee*NN22*NN34)/(2.*sw) + (ee*NN23*NN32*Sinalp)/(2.*sw) + (ee
     +*NN22*NN33*Sinalp)/(2.*sw)
      MGVX379(1) = ((0,0.5)*CONJNN24*CONJNN31*Cosbet*ee)/cw + ((0,0
     +.5)*CONJNN21*CONJNN34*Cosbet*ee)/cw - ((0,0.5)*CONJNN23*CONJ
     +NN31*ee*Sinbet)/cw - ((0,0.5)*CONJNN21*CONJNN33*ee*Sinbet)/c
     +w - ((0,0.5)*CONJNN24*CONJNN32*Cosbet*ee)/sw - ((0,0.5)*CONJ
     +NN22*CONJNN34*Cosbet*ee)/sw + ((0,0.5)*CONJNN23*CONJNN32*ee*
     +Sinbet)/sw + ((0,0.5)*CONJNN22*CONJNN33*ee*Sinbet)/sw
      MGVX379(2) = ((0,-0.5)*Cosbet*ee*NN24*NN31)/cw - ((0,0.5)*Cos
     +bet*ee*NN21*NN34)/cw + ((0,0.5)*ee*NN23*NN31*Sinbet)/cw + ((
     +0,0.5)*ee*NN21*NN33*Sinbet)/cw + ((0,0.5)*Cosbet*ee*NN24*NN3
     +2)/sw + ((0,0.5)*Cosbet*ee*NN22*NN34)/sw - ((0,0.5)*ee*NN23*
     +NN32*Sinbet)/sw - ((0,0.5)*ee*NN22*NN33*Sinbet)/sw
      endif

      return
      end
