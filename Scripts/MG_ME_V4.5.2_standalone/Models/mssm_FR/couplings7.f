      subroutine coup7(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 7
      if(readlha) then


      MGVX426(1) = -((CONJUU21*CONJVV22*Cosalp*ee)/(Sqrt2*sw)) + (C
     +ONJUU22*CONJVV21*ee*Sinalp)/(Sqrt2*sw)
      MGVX426(2) = (ee*Sinalp*UU22*VV21)/(Sqrt2*sw) - (Cosalp*ee*UU
     +21*VV22)/(Sqrt2*sw)
      MGVX427(1) = ((0,1)*CONJUU21*CONJVV22*Cosbet*ee)/(Sqrt2*sw) +
     + ((0,1)*CONJUU22*CONJVV21*ee*Sinbet)/(Sqrt2*sw)
      MGVX427(2) = ((0,-1)*ee*Sinbet*UU22*VV21)/(Sqrt2*sw) - ((0,1)
     +*Cosbet*ee*UU21*VV22)/(Sqrt2*sw)
      MGVX428(1) = 0
      MGVX428(2) = (CONJRl51*ee*NN11)/(cw*Sqrt2) + (CONJRl51*ee*NN1
     +2)/(Sqrt2*sw)
      MGVX429(1) = 0
      MGVX429(2) = (CONJRl42*ee*NN11)/(cw*Sqrt2) + (CONJRl42*ee*NN1
     +2)/(Sqrt2*sw)
      MGVX430(1) = 0
      MGVX430(2) = (CONJRLn31*ee*NN11)/(cw*Sqrt2) - (CONJRLn31*ee*N
     +N12)/(Sqrt2*sw)
      MGVX431(1) = 0
      MGVX431(2) = (CONJRLn22*ee*NN11)/(cw*Sqrt2) - (CONJRLn22*ee*N
     +N12)/(Sqrt2*sw)
      MGVX432(1) = 0
      MGVX432(2) = (CONJRLn13*ee*NN11)/(cw*Sqrt2) - (CONJRLn13*ee*N
     +N12)/(Sqrt2*sw)
      MGVX435(1) = 0
      MGVX435(2) = (CONJRl51*ee*NN21)/(cw*Sqrt2) + (CONJRl51*ee*NN2
     +2)/(Sqrt2*sw)
      MGVX436(1) = 0
      MGVX436(2) = (CONJRl42*ee*NN21)/(cw*Sqrt2) + (CONJRl42*ee*NN2
     +2)/(Sqrt2*sw)
      MGVX437(1) = 0
      MGVX437(2) = (CONJRLn31*ee*NN21)/(cw*Sqrt2) - (CONJRLn31*ee*N
     +N22)/(Sqrt2*sw)
      MGVX438(1) = 0
      MGVX438(2) = (CONJRLn22*ee*NN21)/(cw*Sqrt2) - (CONJRLn22*ee*N
     +N22)/(Sqrt2*sw)
      MGVX439(1) = 0
      MGVX439(2) = (CONJRLn13*ee*NN21)/(cw*Sqrt2) - (CONJRLn13*ee*N
     +N22)/(Sqrt2*sw)
      MGVX442(1) = 0
      MGVX442(2) = (CONJRl51*ee*NN31)/(cw*Sqrt2) + (CONJRl51*ee*NN3
     +2)/(Sqrt2*sw)
      MGVX443(1) = 0
      MGVX443(2) = (CONJRl42*ee*NN31)/(cw*Sqrt2) + (CONJRl42*ee*NN3
     +2)/(Sqrt2*sw)
      MGVX444(1) = 0
      MGVX444(2) = (CONJRLn31*ee*NN31)/(cw*Sqrt2) - (CONJRLn31*ee*N
     +N32)/(Sqrt2*sw)
      endif

      return
      end
