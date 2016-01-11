      subroutine coup36(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 36
      if(readlha) then

      MGVX114 = -((CONJRl13*ee**2*RLn13)/(Sqrt2*sw))
      MGVX274 = (CONJRl13*ee**2*RLn13)/(cw*Sqrt2)
      MGVX235 = (CONJRl63*ee**2*Rl13)/(cw*sw)
      MGVX189 = (CONJRl63*ee**2*Rl13)/(2.*Psw2)
      MGVX320 = (-2*CONJRl63*ee**2*Rl13)/Pcw2 + (CONJRl63*ee**2*Rl1
     +3)/(2.*Pcw2*Psw2)
      MGVX152 = -((CONJRLn13*ee**2*Rl13)/(Sqrt2*sw))
      MGVX287 = (CONJRLn13*ee**2*Rl13)/(cw*Sqrt2)
      MGVX202 = (-2*ee**2*sw)/cw
      MGVX302 = (2*ee**2*Psw2)/Pcw2
      MGVX204 = (CONJRl42*ee**2*Rl42)/(cw*sw) - (2*ee**2*sw)/cw
      MGVX190 = (CONJRl42*ee**2*Rl42)/(2.*Psw2)
      MGVX304 = (2*ee**2*Psw2)/Pcw2 - (2*CONJRl42*ee**2*Rl42)/Pcw2 
     ++ (CONJRl42*ee**2*Rl42)/(2.*Pcw2*Psw2)
      MGVX118 = -((CONJRl42*ee**2*RLn22)/(Sqrt2*sw))
      MGVX276 = (CONJRl42*ee**2*RLn22)/(cw*Sqrt2)
      MGVX154 = -((CONJRLn22*ee**2*Rl42)/(Sqrt2*sw))

      endif

      return
      end
