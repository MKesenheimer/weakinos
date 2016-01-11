      subroutine coup31(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 31
      if(readlha) then

      MGVX322 = -(ee**2*Sinalp*vd)/(2.*Pcw2*Psw2) + (Cosalp*ee**2*T
     +anbet*vd)/(2.*Pcw2*Psw2)
      MGVX177 = (Cosalp*ee**2*vd)/(2.*Psw2) + (ee**2*Sinalp*Tanbet*
     +vd)/(2.*Psw2)
      MGVX307 = (Cosalp*ee**2*vd)/(2.*Pcw2*Psw2) + (ee**2*Sinalp*Ta
     +nbet*vd)/(2.*Pcw2*Psw2)
      MGVX63 = (2*ee**2)/9.
      MGVX88 = (8*ee**2)/9.
      MGVX220 = (CONJRd13*ee**2*Rd13)/(3.*cw*sw) - (2*ee**2*sw)/(9.
     +*cw)
      MGVX234 = (CONJRd13*ee**2*Rd23)/(3.*cw*sw)
      MGVX130 = (CONJRd13*ee**2*Ru13)/(3.*Sqrt2*sw)
      MGVX138 = (CONJRd13*ee**2*Ru63)/(3.*Sqrt2*sw)
      MGVX232 = (CONJRd23*ee**2*Rd13)/(3.*cw*sw)
      MGVX162 = (CONJRu13*ee**2*Rd13)/(3.*Sqrt2*sw)
      MGVX170 = (CONJRu63*ee**2*Rd13)/(3.*Sqrt2*sw)
      MGVX222 = (CONJRd23*ee**2*Rd23)/(3.*cw*sw) - (2*ee**2*sw)/(9.
     +*cw)
      MGVX131 = (CONJRd23*ee**2*Ru13)/(3.*Sqrt2*sw)
      MGVX139 = (CONJRd23*ee**2*Ru63)/(3.*Sqrt2*sw)

      endif

      return
      end
