      subroutine coup33(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 33
      if(readlha) then

      MGVX270 = (2*CONJRu63*ee**2*Ru63)/(3.*cw*sw) - (8*ee**2*sw)/(
     +9.*cw)
      MGVX180 = ee**2/(2.*Psw2)
      MGVX311 = ee**2/(2.*Pcw2*Psw2)
      MGVX143 = ((0,-0.5)*ee**2)/sw
      MGVX286 = ((0,0.5)*ee**2)/cw
      MGVX111 = ((0,0.5)*ee**2)/sw
      MGVX273 = ((0,-0.5)*ee**2)/cw
      MGVX182 = (CONJRd13*ee**2*Rd13)/(2.*Psw2)
      MGVX312 = (2*ee**2*Psw2)/(9.*Pcw2) - (2*CONJRd13*ee**2*Rd13)/
     +(3.*Pcw2) + (CONJRd13*ee**2*Rd13)/(2.*Pcw2*Psw2)
      MGVX184 = (CONJRd13*ee**2*Rd23)/(2.*Psw2)
      MGVX319 = (-2*CONJRd13*ee**2*Rd23)/(3.*Pcw2) + (CONJRd13*ee**
     +2*Rd23)/(2.*Pcw2*Psw2)
      MGVX279 = -(CONJRd13*ee**2*Ru13)/(3.*cw*Sqrt2)
      MGVX283 = -(CONJRd13*ee**2*Ru63)/(3.*cw*Sqrt2)
      MGVX183 = (CONJRd23*ee**2*Rd13)/(2.*Psw2)
      MGVX318 = (-2*CONJRd23*ee**2*Rd13)/(3.*Pcw2) + (CONJRd23*ee**
     +2*Rd13)/(2.*Pcw2*Psw2)

      endif

      return
      end
