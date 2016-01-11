      subroutine coup35(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 35
      if(readlha) then

      MGVX281 = -(CONJRd62*ee**2*Ru42)/(3.*cw*Sqrt2)
      MGVX294 = -(CONJRu42*ee**2*Rd62)/(3.*cw*Sqrt2)
      MGVX142 = (Cosalp*Cosbet*ee**2)/(2.*sw) + (ee**2*Sinalp*Sinbe
     +t)/(2.*sw)
      MGVX285 = -(Cosalp*Cosbet*ee**2)/(2.*cw) - (ee**2*Sinalp*Sinb
     +et)/(2.*cw)
      MGVX160 = (Cosbet*ee**2*Sinalp)/(2.*sw) - (Cosalp*ee**2*Sinbe
     +t)/(2.*sw)
      MGVX291 = -(Cosbet*ee**2*Sinalp)/(2.*cw) + (Cosalp*ee**2*Sinb
     +et)/(2.*cw)
      MGVX54 = 2*ee**2
      MGVX207 = ee**2/(cw*sw) - (2*ee**2*sw)/cw
      MGVX310 = (-2*ee**2)/Pcw2 + ee**2/(2.*Pcw2*Psw2) + (2*ee**2*P
     +sw2)/Pcw2
      MGVX201 = (CONJRl13*ee**2*Rl13)/(cw*sw) - (2*ee**2*sw)/cw
      MGVX188 = (CONJRl13*ee**2*Rl13)/(2.*Psw2)
      MGVX301 = (2*ee**2*Psw2)/Pcw2 - (2*CONJRl13*ee**2*Rl13)/Pcw2 
     ++ (CONJRl13*ee**2*Rl13)/(2.*Pcw2*Psw2)
      MGVX237 = (CONJRl13*ee**2*Rl63)/(cw*sw)
      MGVX192 = (CONJRl13*ee**2*Rl63)/(2.*Psw2)
      MGVX321 = (-2*CONJRl13*ee**2*Rl63)/Pcw2 + (CONJRl13*ee**2*Rl6
     +3)/(2.*Pcw2*Psw2)

      endif

      return
      end
