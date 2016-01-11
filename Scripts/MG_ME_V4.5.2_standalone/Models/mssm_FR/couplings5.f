      subroutine coup5(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 5
      if(readlha) then


      MGVX380(1) = (CONJNN23*CONJNN41*Cosalp*ee)/(2.*cw) + (CONJNN2
     +1*CONJNN43*Cosalp*ee)/(2.*cw) - (CONJNN24*CONJNN41*ee*Sinalp
     +)/(2.*cw) - (CONJNN21*CONJNN44*ee*Sinalp)/(2.*cw) - (CONJNN2
     +3*CONJNN42*Cosalp*ee)/(2.*sw) - (CONJNN22*CONJNN43*Cosalp*ee
     +)/(2.*sw) + (CONJNN24*CONJNN42*ee*Sinalp)/(2.*sw) + (CONJNN2
     +2*CONJNN44*ee*Sinalp)/(2.*sw)
      MGVX380(2) = (Cosalp*ee*NN23*NN41)/(2.*cw) + (Cosalp*ee*NN21*
     +NN43)/(2.*cw) - (ee*NN24*NN41*Sinalp)/(2.*cw) - (ee*NN21*NN4
     +4*Sinalp)/(2.*cw) - (Cosalp*ee*NN23*NN42)/(2.*sw) - (Cosalp*
     +ee*NN22*NN43)/(2.*sw) + (ee*NN24*NN42*Sinalp)/(2.*sw) + (ee*
     +NN22*NN44*Sinalp)/(2.*sw)
      MGVX381(1) = -(CONJNN24*CONJNN41*Cosalp*ee)/(2.*cw) - (CONJNN
     +21*CONJNN44*Cosalp*ee)/(2.*cw) - (CONJNN23*CONJNN41*ee*Sinal
     +p)/(2.*cw) - (CONJNN21*CONJNN43*ee*Sinalp)/(2.*cw) + (CONJNN
     +24*CONJNN42*Cosalp*ee)/(2.*sw) + (CONJNN22*CONJNN44*Cosalp*e
     +e)/(2.*sw) + (CONJNN23*CONJNN42*ee*Sinalp)/(2.*sw) + (CONJNN
     +22*CONJNN43*ee*Sinalp)/(2.*sw)
      MGVX381(2) = -(Cosalp*ee*NN24*NN41)/(2.*cw) - (Cosalp*ee*NN21
     +*NN44)/(2.*cw) - (ee*NN23*NN41*Sinalp)/(2.*cw) - (ee*NN21*NN
     +43*Sinalp)/(2.*cw) + (Cosalp*ee*NN24*NN42)/(2.*sw) + (Cosalp
     +*ee*NN22*NN44)/(2.*sw) + (ee*NN23*NN42*Sinalp)/(2.*sw) + (ee
     +*NN22*NN43*Sinalp)/(2.*sw)
      MGVX382(1) = ((0,0.5)*CONJNN24*CONJNN41*Cosbet*ee)/cw + ((0,0
     +.5)*CONJNN21*CONJNN44*Cosbet*ee)/cw - ((0,0.5)*CONJNN23*CONJ
     +NN41*ee*Sinbet)/cw - ((0,0.5)*CONJNN21*CONJNN43*ee*Sinbet)/c
     +w - ((0,0.5)*CONJNN24*CONJNN42*Cosbet*ee)/sw - ((0,0.5)*CONJ
     +NN22*CONJNN44*Cosbet*ee)/sw + ((0,0.5)*CONJNN23*CONJNN42*ee*
     +Sinbet)/sw + ((0,0.5)*CONJNN22*CONJNN43*ee*Sinbet)/sw
      MGVX382(2) = ((0,-0.5)*Cosbet*ee*NN24*NN41)/cw - ((0,0.5)*Cos
     +bet*ee*NN21*NN44)/cw + ((0,0.5)*ee*NN23*NN41*Sinbet)/cw + ((
     +0,0.5)*ee*NN21*NN43*Sinbet)/cw + ((0,0.5)*Cosbet*ee*NN24*NN4
     +2)/sw + ((0,0.5)*Cosbet*ee*NN22*NN44)/sw - ((0,0.5)*ee*NN23*
     +NN42*Sinbet)/sw - ((0,0.5)*ee*NN22*NN43*Sinbet)/sw
      MGVX383(1) = -((CONJNN31*CONJVV12*Cosbet*ee)/(cw*Sqrt2)) - (C
     +ONJNN34*CONJVV11*Cosbet*ee)/sw - (CONJNN32*CONJVV12*Cosbet*e
     +e)/(Sqrt2*sw)
      MGVX383(2) = -((ee*NN33*Sinbet*UU11)/sw) + (ee*NN31*Sinbet*UU
     +12)/(cw*Sqrt2) + (ee*NN32*Sinbet*UU12)/(Sqrt2*sw)
      MGVX384(1) = -((CONJNN31*CONJVV22*Cosbet*ee)/(cw*Sqrt2)) - (C
     +ONJNN34*CONJVV21*Cosbet*ee)/sw - (CONJNN32*CONJVV22*Cosbet*e
     +e)/(Sqrt2*sw)
      MGVX384(2) = -((ee*NN33*Sinbet*UU21)/sw) + (ee*NN31*Sinbet*UU
     +22)/(cw*Sqrt2) + (ee*NN32*Sinbet*UU22)/(Sqrt2*sw)
      MGVX385(1) = (CONJNN31*CONJNN33*Cosalp*ee)/cw - (CONJNN31*CON
     +JNN34*ee*Sinalp)/cw - (CONJNN32*CONJNN33*Cosalp*ee)/sw + (CO
     +NJNN32*CONJNN34*ee*Sinalp)/sw
      MGVX385(2) = (Cosalp*ee*NN31*NN33)/cw - (ee*NN31*NN34*Sinalp)
     +/cw - (Cosalp*ee*NN32*NN33)/sw + (ee*NN32*NN34*Sinalp)/sw
      MGVX386(1) = -((CONJNN31*CONJNN34*Cosalp*ee)/cw) - (CONJNN31*
     +CONJNN33*ee*Sinalp)/cw + (CONJNN32*CONJNN34*Cosalp*ee)/sw + 
     +(CONJNN32*CONJNN33*ee*Sinalp)/sw
      MGVX386(2) = -((Cosalp*ee*NN31*NN34)/cw) - (ee*NN31*NN33*Sina
     +lp)/cw + (Cosalp*ee*NN32*NN34)/sw + (ee*NN32*NN33*Sinalp)/sw
      MGVX387(1) = ((0,1)*CONJNN31*CONJNN34*Cosbet*ee)/cw - ((0,1)*
     +CONJNN31*CONJNN33*ee*Sinbet)/cw - ((0,1)*CONJNN32*CONJNN34*C
     +osbet*ee)/sw + ((0,1)*CONJNN32*CONJNN33*ee*Sinbet)/sw
      MGVX387(2) = ((0,-1)*Cosbet*ee*NN31*NN34)/cw + ((0,1)*ee*NN31
     +*NN33*Sinbet)/cw + ((0,1)*Cosbet*ee*NN32*NN34)/sw - ((0,1)*e
     +e*NN32*NN33*Sinbet)/sw
      MGVX388(1) = (CONJNN33*CONJNN41*Cosalp*ee)/(2.*cw) + (CONJNN3
     +1*CONJNN43*Cosalp*ee)/(2.*cw) - (CONJNN34*CONJNN41*ee*Sinalp
     +)/(2.*cw) - (CONJNN31*CONJNN44*ee*Sinalp)/(2.*cw) - (CONJNN3
     +3*CONJNN42*Cosalp*ee)/(2.*sw) - (CONJNN32*CONJNN43*Cosalp*ee
     +)/(2.*sw) + (CONJNN34*CONJNN42*ee*Sinalp)/(2.*sw) + (CONJNN3
     +2*CONJNN44*ee*Sinalp)/(2.*sw)
      MGVX388(2) = (Cosalp*ee*NN33*NN41)/(2.*cw) + (Cosalp*ee*NN31*
     +NN43)/(2.*cw) - (ee*NN34*NN41*Sinalp)/(2.*cw) - (ee*NN31*NN4
     +4*Sinalp)/(2.*cw) - (Cosalp*ee*NN33*NN42)/(2.*sw) - (Cosalp*
     +ee*NN32*NN43)/(2.*sw) + (ee*NN34*NN42*Sinalp)/(2.*sw) + (ee*
     +NN32*NN44*Sinalp)/(2.*sw)
      MGVX389(1) = -(CONJNN34*CONJNN41*Cosalp*ee)/(2.*cw) - (CONJNN
     +31*CONJNN44*Cosalp*ee)/(2.*cw) - (CONJNN33*CONJNN41*ee*Sinal
     +p)/(2.*cw) - (CONJNN31*CONJNN43*ee*Sinalp)/(2.*cw) + (CONJNN
     +34*CONJNN42*Cosalp*ee)/(2.*sw) + (CONJNN32*CONJNN44*Cosalp*e
     +e)/(2.*sw) + (CONJNN33*CONJNN42*ee*Sinalp)/(2.*sw) + (CONJNN
     +32*CONJNN43*ee*Sinalp)/(2.*sw)
      MGVX389(2) = -(Cosalp*ee*NN34*NN41)/(2.*cw) - (Cosalp*ee*NN31
     +*NN44)/(2.*cw) - (ee*NN33*NN41*Sinalp)/(2.*cw) - (ee*NN31*NN
     +43*Sinalp)/(2.*cw) + (Cosalp*ee*NN34*NN42)/(2.*sw) + (Cosalp
     +*ee*NN32*NN44)/(2.*sw) + (ee*NN33*NN42*Sinalp)/(2.*sw) + (ee
     +*NN32*NN43*Sinalp)/(2.*sw)
      MGVX390(1) = ((0,0.5)*CONJNN34*CONJNN41*Cosbet*ee)/cw + ((0,0
     +.5)*CONJNN31*CONJNN44*Cosbet*ee)/cw - ((0,0.5)*CONJNN33*CONJ
     +NN41*ee*Sinbet)/cw - ((0,0.5)*CONJNN31*CONJNN43*ee*Sinbet)/c
     +w - ((0,0.5)*CONJNN34*CONJNN42*Cosbet*ee)/sw - ((0,0.5)*CONJ
     +NN32*CONJNN44*Cosbet*ee)/sw + ((0,0.5)*CONJNN33*CONJNN42*ee*
     +Sinbet)/sw + ((0,0.5)*CONJNN32*CONJNN43*ee*Sinbet)/sw
      MGVX390(2) = ((0,-0.5)*Cosbet*ee*NN34*NN41)/cw - ((0,0.5)*Cos
     +bet*ee*NN31*NN44)/cw + ((0,0.5)*ee*NN33*NN41*Sinbet)/cw + ((
     +0,0.5)*ee*NN31*NN43*Sinbet)/cw + ((0,0.5)*Cosbet*ee*NN34*NN4
     +2)/sw + ((0,0.5)*Cosbet*ee*NN32*NN44)/sw - ((0,0.5)*ee*NN33*
     +NN42*Sinbet)/sw - ((0,0.5)*ee*NN32*NN43*Sinbet)/sw
      MGVX391(1) = -((CONJNN41*CONJVV12*Cosbet*ee)/(cw*Sqrt2)) - (C
     +ONJNN44*CONJVV11*Cosbet*ee)/sw - (CONJNN42*CONJVV12*Cosbet*e
     +e)/(Sqrt2*sw)
      MGVX391(2) = -((ee*NN43*Sinbet*UU11)/sw) + (ee*NN41*Sinbet*UU
     +12)/(cw*Sqrt2) + (ee*NN42*Sinbet*UU12)/(Sqrt2*sw)
      MGVX392(1) = -((CONJNN41*CONJVV22*Cosbet*ee)/(cw*Sqrt2)) - (C
     +ONJNN44*CONJVV21*Cosbet*ee)/sw - (CONJNN42*CONJVV22*Cosbet*e
     +e)/(Sqrt2*sw)
      MGVX392(2) = -((ee*NN43*Sinbet*UU21)/sw) + (ee*NN41*Sinbet*UU
     +22)/(cw*Sqrt2) + (ee*NN42*Sinbet*UU22)/(Sqrt2*sw)
      MGVX393(1) = (CONJNN41*CONJNN43*Cosalp*ee)/cw - (CONJNN41*CON
     +JNN44*ee*Sinalp)/cw - (CONJNN42*CONJNN43*Cosalp*ee)/sw + (CO
     +NJNN42*CONJNN44*ee*Sinalp)/sw
      MGVX393(2) = (Cosalp*ee*NN41*NN43)/cw - (ee*NN41*NN44*Sinalp)
     +/cw - (Cosalp*ee*NN42*NN43)/sw + (ee*NN42*NN44*Sinalp)/sw
      MGVX394(1) = -((CONJNN41*CONJNN44*Cosalp*ee)/cw) - (CONJNN41*
     +CONJNN43*ee*Sinalp)/cw + (CONJNN42*CONJNN44*Cosalp*ee)/sw + 
     +(CONJNN42*CONJNN43*ee*Sinalp)/sw
      MGVX394(2) = -((Cosalp*ee*NN41*NN44)/cw) - (ee*NN41*NN43*Sina
     +lp)/cw + (Cosalp*ee*NN42*NN44)/sw + (ee*NN42*NN43*Sinalp)/sw
      endif

      return
      end
