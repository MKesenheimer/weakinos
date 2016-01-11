      subroutine coup26(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 26
      if(readlha) then

      MGVX706 = CONJRd26*CONJtd33*Ru13*Sinbet + CONJRd23*Cosbet*Ru1
     +6*tu33 - (CONJRd23*ee**2*Ru13*Sinbet*vd)/(Psw2*Sqrt2) + (CON
     +JRd23*Pyukd32*Ru13*Sinbet*vd)/Sqrt2 + (CONJRd23*Pyuku32*Ru13
     +*Sinbet*vd)/Sqrt2 + CONJRd26*Cosbet*MMMU*Ru13*yukd3 + CONJMM
     +MU*CONJRd23*Ru16*Sinbet*yuku3 + (CONJRd26*Cosbet*Ru16*vd*yuk
     +d3*yuku3)/Sqrt2 + (CONJRd26*Ru16*Sinbet*Tanbet*vd*yukd3*yuku
     +3)/Sqrt2
      MGVX710 = CONJRd26*CONJtd33*Ru63*Sinbet + CONJRd23*Cosbet*Ru6
     +6*tu33 - (CONJRd23*ee**2*Ru63*Sinbet*vd)/(Psw2*Sqrt2) + (CON
     +JRd23*Pyukd32*Ru63*Sinbet*vd)/Sqrt2 + (CONJRd23*Pyuku32*Ru63
     +*Sinbet*vd)/Sqrt2 + CONJRd26*Cosbet*MMMU*Ru63*yukd3 + CONJMM
     +MU*CONJRd23*Ru66*Sinbet*yuku3 + (CONJRd26*Cosbet*Ru66*vd*yuk
     +d3*yuku3)/Sqrt2 + (CONJRd26*Ru66*Sinbet*Tanbet*vd*yukd3*yuku
     +3)/Sqrt2
      MGVX792 = -((CONJRd51*ee**2*Ru51*Sinbet*vd)/(Psw2*Sqrt2))
      MGVX790 = -((CONJRd62*ee**2*Ru42*Sinbet*vd)/(Psw2*Sqrt2))
      MGVX767 = -(Cosalp*Cosbet*ee**2*Sinbet*vd)/(2.*Pcw2) - (Cosal
     +p*Cosbet*ee**2*Sinbet*vd)/Psw2 + (Cosalp*Cosbet*ee**2*Sinbet
     +*vd)/(4.*Pcw2*Psw2) - (ee**2*Sinalp*vd*Cos(bet)**2)/(2.*Pcw2
     +) + (ee**2*Sinalp*vd*Cos(bet)**2)/(4.*Pcw2*Psw2) + (ee**2*Si
     +nalp*vd*Sin(bet)**2)/(2.*Pcw2) + (ee**2*Sinalp*vd*Sin(bet)**
     +2)/Psw2 - (ee**2*Sinalp*vd*Sin(bet)**2)/(4.*Pcw2*Psw2) + (Co
     +salp*ee**2*Tanbet*vd*Sin(bet)**2)/(2.*Pcw2) - (Cosalp*ee**2*
     +Tanbet*vd*Sin(bet)**2)/(4.*Pcw2*Psw2)
      MGVX752 = -(Cosbet*ee**2*Sinalp*Sinbet*vd)/(2.*Pcw2) - (Cosbe
     +t*ee**2*Sinalp*Sinbet*vd)/Psw2 + (Cosbet*ee**2*Sinalp*Sinbet
     +*vd)/(4.*Pcw2*Psw2) + (Cosalp*ee**2*vd*Cos(bet)**2)/(2.*Pcw2
     +) - (Cosalp*ee**2*vd*Cos(bet)**2)/(4.*Pcw2*Psw2) - (Cosalp*e
     +e**2*vd*Sin(bet)**2)/(2.*Pcw2) - (Cosalp*ee**2*vd*Sin(bet)**
     +2)/Psw2 + (Cosalp*ee**2*vd*Sin(bet)**2)/(4.*Pcw2*Psw2) + (ee
     +**2*Sinalp*Tanbet*vd*Sin(bet)**2)/(2.*Pcw2) - (ee**2*Sinalp*
     +Tanbet*vd*Sin(bet)**2)/(4.*Pcw2*Psw2)
      MGVX722 = CONJRl16*CONJte33*RLn13*Sinbet - (CONJRl13*ee**2*RL
     +n13*Sinbet*vd)/(Psw2*Sqrt2) + (CONJRl13*Pyukl32*RLn13*Sinbet
     +*vd)/Sqrt2 + CONJRl16*Cosbet*MMMU*RLn13*yukl3
      MGVX787 = -((CONJRl42*ee**2*RLn22*Sinbet*vd)/(Psw2*Sqrt2))
      MGVX788 = -((CONJRl51*ee**2*RLn31*Sinbet*vd)/(Psw2*Sqrt2))
      MGVX723 = CONJRl66*CONJte33*RLn13*Sinbet - (CONJRl63*ee**2*RL
     +n13*Sinbet*vd)/(Psw2*Sqrt2) + (CONJRl63*Pyukl32*RLn13*Sinbet
     +*vd)/Sqrt2 + CONJRl66*Cosbet*MMMU*RLn13*yukl3
      MGVX687 = (CONJRd16*CONJtd33*Rd13*Sinalp)/Sqrt2 + (CONJRd13*R
     +d16*Sinalp*td33)/Sqrt2 + (CONJRd13*ee**2*Rd13*Sinalp*vd)/(6.
     +*Pcw2) - (CONJRd13*ee**2*Rd13*Sinalp*vd)/(4.*Pcw2*Psw2) + CO
     +NJRd13*Pyukd32*Rd13*Sinalp*vd - (CONJRd16*ee**2*Rd16*Sinalp*
     +vd)/(6.*Pcw2) + CONJRd16*Pyukd32*Rd16*Sinalp*vd + (CONJRd13*
     +Cosalp*ee**2*Rd13*Tanbet*vd)/(6.*Pcw2) - (CONJRd13*Cosalp*ee
     +**2*Rd13*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRd16*Cosalp*ee**2*
     +Rd16*Tanbet*vd)/(6.*Pcw2) + (CONJRd16*Cosalp*MMMU*Rd13*yukd3
     +)/Sqrt2 + (CONJMMMU*CONJRd13*Cosalp*Rd16*yukd3)/Sqrt2
      MGVX691 = (CONJRd16*CONJtd33*Rd23*Sinalp)/Sqrt2 + (CONJRd13*R
     +d26*Sinalp*td33)/Sqrt2 + (CONJRd13*ee**2*Rd23*Sinalp*vd)/(6.
     +*Pcw2) - (CONJRd13*ee**2*Rd23*Sinalp*vd)/(4.*Pcw2*Psw2) + CO
     +NJRd13*Pyukd32*Rd23*Sinalp*vd - (CONJRd16*ee**2*Rd26*Sinalp*
     +vd)/(6.*Pcw2) + CONJRd16*Pyukd32*Rd26*Sinalp*vd + (CONJRd13*
     +Cosalp*ee**2*Rd23*Tanbet*vd)/(6.*Pcw2) - (CONJRd13*Cosalp*ee
     +**2*Rd23*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRd16*Cosalp*ee**2*
     +Rd26*Tanbet*vd)/(6.*Pcw2) + (CONJRd16*Cosalp*MMMU*Rd23*yukd3
     +)/Sqrt2 + (CONJMMMU*CONJRd13*Cosalp*Rd26*yukd3)/Sqrt2
      MGVX689 = (CONJRd26*CONJtd33*Rd13*Sinalp)/Sqrt2 + (CONJRd23*R
     +d16*Sinalp*td33)/Sqrt2 + (CONJRd23*ee**2*Rd13*Sinalp*vd)/(6.
     +*Pcw2) - (CONJRd23*ee**2*Rd13*Sinalp*vd)/(4.*Pcw2*Psw2) + CO
     +NJRd23*Pyukd32*Rd13*Sinalp*vd - (CONJRd26*ee**2*Rd16*Sinalp*
     +vd)/(6.*Pcw2) + CONJRd26*Pyukd32*Rd16*Sinalp*vd + (CONJRd23*
     +Cosalp*ee**2*Rd13*Tanbet*vd)/(6.*Pcw2) - (CONJRd23*Cosalp*ee
     +**2*Rd13*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRd26*Cosalp*ee**2*
     +Rd16*Tanbet*vd)/(6.*Pcw2) + (CONJRd26*Cosalp*MMMU*Rd13*yukd3
     +)/Sqrt2 + (CONJMMMU*CONJRd23*Cosalp*Rd16*yukd3)/Sqrt2
      MGVX693 = (CONJRd26*CONJtd33*Rd23*Sinalp)/Sqrt2 + (CONJRd23*R
     +d26*Sinalp*td33)/Sqrt2 + (CONJRd23*ee**2*Rd23*Sinalp*vd)/(6.
     +*Pcw2) - (CONJRd23*ee**2*Rd23*Sinalp*vd)/(4.*Pcw2*Psw2) + CO
     +NJRd23*Pyukd32*Rd23*Sinalp*vd - (CONJRd26*ee**2*Rd26*Sinalp*
     +vd)/(6.*Pcw2) + CONJRd26*Pyukd32*Rd26*Sinalp*vd + (CONJRd23*
     +Cosalp*ee**2*Rd23*Tanbet*vd)/(6.*Pcw2) - (CONJRd23*Cosalp*ee
     +**2*Rd23*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRd26*Cosalp*ee**2*
     +Rd26*Tanbet*vd)/(6.*Pcw2) + (CONJRd26*Cosalp*MMMU*Rd23*yukd3
     +)/Sqrt2 + (CONJMMMU*CONJRd23*Cosalp*Rd26*yukd3)/Sqrt2
      MGVX769 = -(CONJRd35*ee**2*Rd35*Sinalp*vd)/(6.*Pcw2) - (CONJR
     +d35*Cosalp*ee**2*Rd35*Tanbet*vd)/(6.*Pcw2)

      endif

      return
      end
