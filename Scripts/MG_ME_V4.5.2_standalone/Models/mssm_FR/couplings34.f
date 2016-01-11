      subroutine coup34(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 34
      if(readlha) then

      MGVX292 = -(CONJRu13*ee**2*Rd13)/(3.*cw*Sqrt2)
      MGVX296 = -(CONJRu63*ee**2*Rd13)/(3.*cw*Sqrt2)
      MGVX185 = (CONJRd23*ee**2*Rd23)/(2.*Psw2)
      MGVX313 = (2*ee**2*Psw2)/(9.*Pcw2) - (2*CONJRd23*ee**2*Rd23)/
     +(3.*Pcw2) + (CONJRd23*ee**2*Rd23)/(2.*Pcw2*Psw2)
      MGVX280 = -(CONJRd23*ee**2*Ru13)/(3.*cw*Sqrt2)
      MGVX284 = -(CONJRd23*ee**2*Ru63)/(3.*cw*Sqrt2)
      MGVX293 = -(CONJRu13*ee**2*Rd23)/(3.*cw*Sqrt2)
      MGVX297 = -(CONJRu63*ee**2*Rd23)/(3.*cw*Sqrt2)
      MGVX314 = (2*ee**2*Psw2)/(9.*Pcw2)
      MGVX186 = (CONJRd51*ee**2*Rd51)/(2.*Psw2)
      MGVX316 = (2*ee**2*Psw2)/(9.*Pcw2) - (2*CONJRd51*ee**2*Rd51)/
     +(3.*Pcw2) + (CONJRd51*ee**2*Rd51)/(2.*Pcw2*Psw2)
      MGVX282 = -(CONJRd51*ee**2*Ru51)/(3.*cw*Sqrt2)
      MGVX295 = -(CONJRu51*ee**2*Rd51)/(3.*cw*Sqrt2)
      MGVX187 = (CONJRd62*ee**2*Rd62)/(2.*Psw2)
      MGVX317 = (2*ee**2*Psw2)/(9.*Pcw2) - (2*CONJRd62*ee**2*Rd62)/
     +(3.*Pcw2) + (CONJRd62*ee**2*Rd62)/(2.*Pcw2*Psw2)

      endif

      return
      end
