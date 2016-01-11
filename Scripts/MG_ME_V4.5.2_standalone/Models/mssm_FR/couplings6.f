      subroutine coup6(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 6
      if(readlha) then


      MGVX395(1) = ((0,1)*CONJNN41*CONJNN44*Cosbet*ee)/cw - ((0,1)*
     +CONJNN41*CONJNN43*ee*Sinbet)/cw - ((0,1)*CONJNN42*CONJNN44*C
     +osbet*ee)/sw + ((0,1)*CONJNN42*CONJNN43*ee*Sinbet)/sw
      MGVX395(2) = ((0,-1)*Cosbet*ee*NN41*NN44)/cw + ((0,1)*ee*NN41
     +*NN43*Sinbet)/cw + ((0,1)*Cosbet*ee*NN42*NN44)/sw - ((0,1)*e
     +e*NN42*NN43*Sinbet)/sw
      MGVX396(1) = -((CONJNN11*CONJRl16*ee*Sqrt2)/cw) - CONJNN13*CO
     +NJRl13*yukl3
      MGVX396(2) = (CONJRl13*ee*NN11)/(cw*Sqrt2) + (CONJRl13*ee*NN1
     +2)/(Sqrt2*sw) - CONJRl16*NN13*yukl3
      MGVX397(1) = -((CONJNN11*CONJRl66*ee*Sqrt2)/cw) - CONJNN13*CO
     +NJRl63*yukl3
      MGVX397(2) = (CONJRl63*ee*NN11)/(cw*Sqrt2) + (CONJRl63*ee*NN1
     +2)/(Sqrt2*sw) - CONJRl66*NN13*yukl3
      MGVX398(1) = -((CONJNN21*CONJRl16*ee*Sqrt2)/cw) - CONJNN23*CO
     +NJRl13*yukl3
      MGVX398(2) = (CONJRl13*ee*NN21)/(cw*Sqrt2) + (CONJRl13*ee*NN2
     +2)/(Sqrt2*sw) - CONJRl16*NN23*yukl3
      MGVX399(1) = -((CONJNN21*CONJRl66*ee*Sqrt2)/cw) - CONJNN23*CO
     +NJRl63*yukl3
      MGVX399(2) = (CONJRl63*ee*NN21)/(cw*Sqrt2) + (CONJRl63*ee*NN2
     +2)/(Sqrt2*sw) - CONJRl66*NN23*yukl3
      MGVX400(1) = -((CONJNN31*CONJRl16*ee*Sqrt2)/cw) - CONJNN33*CO
     +NJRl13*yukl3
      MGVX400(2) = (CONJRl13*ee*NN31)/(cw*Sqrt2) + (CONJRl13*ee*NN3
     +2)/(Sqrt2*sw) - CONJRl16*NN33*yukl3
      MGVX401(1) = -((CONJNN31*CONJRl66*ee*Sqrt2)/cw) - CONJNN33*CO
     +NJRl63*yukl3
      MGVX401(2) = (CONJRl63*ee*NN31)/(cw*Sqrt2) + (CONJRl63*ee*NN3
     +2)/(Sqrt2*sw) - CONJRl66*NN33*yukl3
      MGVX402(1) = -((CONJNN41*CONJRl16*ee*Sqrt2)/cw) - CONJNN43*CO
     +NJRl13*yukl3
      MGVX402(2) = (CONJRl13*ee*NN41)/(cw*Sqrt2) + (CONJRl13*ee*NN4
     +2)/(Sqrt2*sw) - CONJRl16*NN43*yukl3
      MGVX403(1) = -((CONJNN41*CONJRl66*ee*Sqrt2)/cw) - CONJNN43*CO
     +NJRl63*yukl3
      MGVX403(2) = (CONJRl63*ee*NN41)/(cw*Sqrt2) + (CONJRl63*ee*NN4
     +2)/(Sqrt2*sw) - CONJRl66*NN43*yukl3
      MGVX420(1) = -((CONJUU11*CONJVV12*Cosalp*ee)/(Sqrt2*sw)) + (C
     +ONJUU12*CONJVV11*ee*Sinalp)/(Sqrt2*sw)
      MGVX420(2) = (ee*Sinalp*UU12*VV11)/(Sqrt2*sw) - (Cosalp*ee*UU
     +11*VV12)/(Sqrt2*sw)
      MGVX421(1) = ((0,1)*CONJUU11*CONJVV12*Cosbet*ee)/(Sqrt2*sw) +
     + ((0,1)*CONJUU12*CONJVV11*ee*Sinbet)/(Sqrt2*sw)
      MGVX421(2) = ((0,-1)*ee*Sinbet*UU12*VV11)/(Sqrt2*sw) - ((0,1)
     +*Cosbet*ee*UU11*VV12)/(Sqrt2*sw)
      MGVX422(1) = -((CONJUU11*CONJVV22*Cosalp*ee)/(Sqrt2*sw)) + (C
     +ONJUU12*CONJVV21*ee*Sinalp)/(Sqrt2*sw)
      MGVX422(2) = (ee*Sinalp*UU22*VV11)/(Sqrt2*sw) - (Cosalp*ee*UU
     +21*VV12)/(Sqrt2*sw)
      MGVX423(1) = ((0,1)*CONJUU11*CONJVV22*Cosbet*ee)/(Sqrt2*sw) +
     + ((0,1)*CONJUU12*CONJVV21*ee*Sinbet)/(Sqrt2*sw)
      MGVX423(2) = ((0,-1)*ee*Sinbet*UU22*VV11)/(Sqrt2*sw) - ((0,1)
     +*Cosbet*ee*UU21*VV12)/(Sqrt2*sw)
      MGVX424(1) = -((CONJUU21*CONJVV12*Cosalp*ee)/(Sqrt2*sw)) + (C
     +ONJUU22*CONJVV11*ee*Sinalp)/(Sqrt2*sw)
      MGVX424(2) = (ee*Sinalp*UU12*VV21)/(Sqrt2*sw) - (Cosalp*ee*UU
     +11*VV22)/(Sqrt2*sw)
      MGVX425(1) = ((0,1)*CONJUU21*CONJVV12*Cosbet*ee)/(Sqrt2*sw) +
     + ((0,1)*CONJUU22*CONJVV11*ee*Sinbet)/(Sqrt2*sw)
      MGVX425(2) = ((0,-1)*ee*Sinbet*UU12*VV21)/(Sqrt2*sw) - ((0,1)
     +*Cosbet*ee*UU11*VV22)/(Sqrt2*sw)
      endif

      return
      end
