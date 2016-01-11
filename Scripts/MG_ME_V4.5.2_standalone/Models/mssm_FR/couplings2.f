      subroutine coup2(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 2
      if(readlha) then
      MGVX3 = -ee
      MGVX5 = (cw*ee)/sw
      MGVX4 = ee**2
      MGVX7 = (cw*ee**2)/sw
      MGVX6 = -(ee**2/Psw2)
      MGVX8 = (ee**2*Pcw2)/Psw2


      MGVX12(1) = ee
      MGVX12(2) = ee
      MGVX27(1) = -(ee/(Sqrt2*sw))
      MGVX27(2) = 0
      MGVX39(1) = ee/(2.*cw*sw) - (ee*sw)/cw
      MGVX39(2) = -((ee*sw)/cw)
      MGVX45(1) = -ee/(2.*cw*sw)
      MGVX45(2) = 0
      MGVX343(1) = -ee
      MGVX343(2) = -ee
      MGVX346(1) = -((CONJUU12*CONJVV11*Cosalp*ee)/(Sqrt2*sw)) - (C
     +ONJUU11*CONJVV12*ee*Sinalp)/(Sqrt2*sw)
      MGVX346(2) = -((Cosalp*ee*UU12*VV11)/(Sqrt2*sw)) - (ee*Sinalp
     +*UU11*VV12)/(Sqrt2*sw)
      MGVX347(1) = -((CONJUU12*CONJVV21*Cosalp*ee)/(Sqrt2*sw)) - (C
     +ONJUU11*CONJVV22*ee*Sinalp)/(Sqrt2*sw)
      MGVX347(2) = -((Cosalp*ee*UU22*VV11)/(Sqrt2*sw)) - (ee*Sinalp
     +*UU21*VV12)/(Sqrt2*sw)
      MGVX348(1) = -((CONJUU22*CONJVV11*Cosalp*ee)/(Sqrt2*sw)) - (C
     +ONJUU21*CONJVV12*ee*Sinalp)/(Sqrt2*sw)
      MGVX348(2) = -((Cosalp*ee*UU12*VV21)/(Sqrt2*sw)) - (ee*Sinalp
     +*UU11*VV22)/(Sqrt2*sw)
      MGVX349(1) = -((CONJUU22*CONJVV21*Cosalp*ee)/(Sqrt2*sw)) - (C
     +ONJUU21*CONJVV22*ee*Sinalp)/(Sqrt2*sw)
      MGVX349(2) = -((Cosalp*ee*UU22*VV21)/(Sqrt2*sw)) - (ee*Sinalp
     +*UU21*VV22)/(Sqrt2*sw)
      endif

      return
      end
