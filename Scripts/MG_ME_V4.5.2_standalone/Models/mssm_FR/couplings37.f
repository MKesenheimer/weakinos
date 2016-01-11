      subroutine coup37(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 37
      if(readlha) then

      MGVX288 = (CONJRLn22*ee**2*Rl42)/(cw*Sqrt2)
      MGVX205 = (CONJRl51*ee**2*Rl51)/(cw*sw) - (2*ee**2*sw)/cw
      MGVX191 = (CONJRl51*ee**2*Rl51)/(2.*Psw2)
      MGVX305 = (2*ee**2*Psw2)/Pcw2 - (2*CONJRl51*ee**2*Rl51)/Pcw2 
     ++ (CONJRl51*ee**2*Rl51)/(2.*Pcw2*Psw2)
      MGVX120 = -((CONJRl51*ee**2*RLn31)/(Sqrt2*sw))
      MGVX277 = (CONJRl51*ee**2*RLn31)/(cw*Sqrt2)
      MGVX156 = -((CONJRLn31*ee**2*Rl51)/(Sqrt2*sw))
      MGVX289 = (CONJRLn31*ee**2*Rl51)/(cw*Sqrt2)
      MGVX206 = (CONJRl63*ee**2*Rl63)/(cw*sw) - (2*ee**2*sw)/cw
      MGVX193 = (CONJRl63*ee**2*Rl63)/(2.*Psw2)
      MGVX306 = (2*ee**2*Psw2)/Pcw2 - (2*CONJRl63*ee**2*Rl63)/Pcw2 
     ++ (CONJRl63*ee**2*Rl63)/(2.*Pcw2*Psw2)
      MGVX115 = -((CONJRl63*ee**2*RLn13)/(Sqrt2*sw))
      MGVX275 = (CONJRl63*ee**2*RLn13)/(cw*Sqrt2)
      MGVX158 = -((CONJRLn13*ee**2*Rl63)/(Sqrt2*sw))
      MGVX290 = (CONJRLn13*ee**2*Rl63)/(cw*Sqrt2)

      endif

      return
      end
