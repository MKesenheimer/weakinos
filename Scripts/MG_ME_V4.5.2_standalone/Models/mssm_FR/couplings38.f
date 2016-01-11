      subroutine coup38(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 38
      if(readlha) then

      MGVX195 = (CONJRu13*ee**2*Ru13)/(2.*Psw2)
      MGVX323 = (8*ee**2*Psw2)/(9.*Pcw2) - (4*CONJRu13*ee**2*Ru13)/
     +(3.*Pcw2) + (CONJRu13*ee**2*Ru13)/(2.*Pcw2*Psw2)
      MGVX199 = (CONJRu13*ee**2*Ru63)/(2.*Psw2)
      MGVX329 = (-4*CONJRu13*ee**2*Ru63)/(3.*Pcw2) + (CONJRu13*ee**
     +2*Ru63)/(2.*Pcw2*Psw2)
      MGVX198 = (CONJRu63*ee**2*Ru13)/(2.*Psw2)
      MGVX328 = (-4*CONJRu63*ee**2*Ru13)/(3.*Pcw2) + (CONJRu63*ee**
     +2*Ru13)/(2.*Pcw2*Psw2)
      MGVX324 = (8*ee**2*Psw2)/(9.*Pcw2)
      MGVX196 = (CONJRu42*ee**2*Ru42)/(2.*Psw2)
      MGVX326 = (8*ee**2*Psw2)/(9.*Pcw2) - (4*CONJRu42*ee**2*Ru42)/
     +(3.*Pcw2) + (CONJRu42*ee**2*Ru42)/(2.*Pcw2*Psw2)
      MGVX197 = (CONJRu51*ee**2*Ru51)/(2.*Psw2)
      MGVX327 = (8*ee**2*Psw2)/(9.*Pcw2) - (4*CONJRu51*ee**2*Ru51)/
     +(3.*Pcw2) + (CONJRu51*ee**2*Ru51)/(2.*Pcw2*Psw2)
      MGVX200 = (CONJRu63*ee**2*Ru63)/(2.*Psw2)
      MGVX330 = (8*ee**2*Psw2)/(9.*Pcw2) - (4*CONJRu63*ee**2*Ru63)/
     +(3.*Pcw2) + (CONJRu63*ee**2*Ru63)/(2.*Pcw2*Psw2)
      MGVX62 = -ee/3.
      MGVX86 = (2*ee)/3.

      endif

      return
      end
