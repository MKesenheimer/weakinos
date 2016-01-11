      subroutine coup28(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 28
      if(readlha) then

      MGVX736 = -((CONJRu16*CONJtu33*Cosalp*Ru13)/Sqrt2) - (CONJRu1
     +3*Cosalp*Ru16*tu33)/Sqrt2 - (CONJRu13*ee**2*Ru13*Sinalp*vd)/
     +(3.*Pcw2) + (CONJRu13*ee**2*Ru13*Sinalp*vd)/(4.*Pcw2*Psw2) +
     + (CONJRu16*ee**2*Ru16*Sinalp*vd)/(3.*Pcw2) - (CONJRu13*Cosal
     +p*ee**2*Ru13*Tanbet*vd)/(3.*Pcw2) + (CONJRu13*Cosalp*ee**2*R
     +u13*Tanbet*vd)/(4.*Pcw2*Psw2) - CONJRu13*Cosalp*Pyuku32*Ru13
     +*Tanbet*vd + (CONJRu16*Cosalp*ee**2*Ru16*Tanbet*vd)/(3.*Pcw2
     +) - CONJRu16*Cosalp*Pyuku32*Ru16*Tanbet*vd - (CONJRu16*MMMU*
     +Ru13*Sinalp*yuku3)/Sqrt2 - (CONJMMMU*CONJRu13*Ru16*Sinalp*yu
     +ku3)/Sqrt2
      MGVX742 = -((CONJRu16*CONJtu33*Cosalp*Ru63)/Sqrt2) - (CONJRu1
     +3*Cosalp*Ru66*tu33)/Sqrt2 - (CONJRu13*ee**2*Ru63*Sinalp*vd)/
     +(3.*Pcw2) + (CONJRu13*ee**2*Ru63*Sinalp*vd)/(4.*Pcw2*Psw2) +
     + (CONJRu16*ee**2*Ru66*Sinalp*vd)/(3.*Pcw2) - (CONJRu13*Cosal
     +p*ee**2*Ru63*Tanbet*vd)/(3.*Pcw2) + (CONJRu13*Cosalp*ee**2*R
     +u63*Tanbet*vd)/(4.*Pcw2*Psw2) - CONJRu13*Cosalp*Pyuku32*Ru63
     +*Tanbet*vd + (CONJRu16*Cosalp*ee**2*Ru66*Tanbet*vd)/(3.*Pcw2
     +) - CONJRu16*Cosalp*Pyuku32*Ru66*Tanbet*vd - (CONJRu16*MMMU*
     +Ru63*Sinalp*yuku3)/Sqrt2 - (CONJMMMU*CONJRu13*Ru66*Sinalp*yu
     +ku3)/Sqrt2
      MGVX739 = -((CONJRu66*CONJtu33*Cosalp*Ru13)/Sqrt2) - (CONJRu6
     +3*Cosalp*Ru16*tu33)/Sqrt2 - (CONJRu63*ee**2*Ru13*Sinalp*vd)/
     +(3.*Pcw2) + (CONJRu63*ee**2*Ru13*Sinalp*vd)/(4.*Pcw2*Psw2) +
     + (CONJRu66*ee**2*Ru16*Sinalp*vd)/(3.*Pcw2) - (CONJRu63*Cosal
     +p*ee**2*Ru13*Tanbet*vd)/(3.*Pcw2) + (CONJRu63*Cosalp*ee**2*R
     +u13*Tanbet*vd)/(4.*Pcw2*Psw2) - CONJRu63*Cosalp*Pyuku32*Ru13
     +*Tanbet*vd + (CONJRu66*Cosalp*ee**2*Ru16*Tanbet*vd)/(3.*Pcw2
     +) - CONJRu66*Cosalp*Pyuku32*Ru16*Tanbet*vd - (CONJRu66*MMMU*
     +Ru13*Sinalp*yuku3)/Sqrt2 - (CONJMMMU*CONJRu63*Ru16*Sinalp*yu
     +ku3)/Sqrt2
      MGVX778 = (CONJRu24*ee**2*Ru24*Sinalp*vd)/(3.*Pcw2) + (CONJRu
     +24*Cosalp*ee**2*Ru24*Tanbet*vd)/(3.*Pcw2)
      MGVX780 = (CONJRu35*ee**2*Ru35*Sinalp*vd)/(3.*Pcw2) + (CONJRu
     +35*Cosalp*ee**2*Ru35*Tanbet*vd)/(3.*Pcw2)
      MGVX782 = -(CONJRu42*ee**2*Ru42*Sinalp*vd)/(3.*Pcw2) + (CONJR
     +u42*ee**2*Ru42*Sinalp*vd)/(4.*Pcw2*Psw2) - (CONJRu42*Cosalp*
     +ee**2*Ru42*Tanbet*vd)/(3.*Pcw2) + (CONJRu42*Cosalp*ee**2*Ru4
     +2*Tanbet*vd)/(4.*Pcw2*Psw2)
      MGVX784 = -(CONJRu51*ee**2*Ru51*Sinalp*vd)/(3.*Pcw2) + (CONJR
     +u51*ee**2*Ru51*Sinalp*vd)/(4.*Pcw2*Psw2) - (CONJRu51*Cosalp*
     +ee**2*Ru51*Tanbet*vd)/(3.*Pcw2) + (CONJRu51*Cosalp*ee**2*Ru5
     +1*Tanbet*vd)/(4.*Pcw2*Psw2)
      MGVX745 = -((CONJRu66*CONJtu33*Cosalp*Ru63)/Sqrt2) - (CONJRu6
     +3*Cosalp*Ru66*tu33)/Sqrt2 - (CONJRu63*ee**2*Ru63*Sinalp*vd)/
     +(3.*Pcw2) + (CONJRu63*ee**2*Ru63*Sinalp*vd)/(4.*Pcw2*Psw2) +
     + (CONJRu66*ee**2*Ru66*Sinalp*vd)/(3.*Pcw2) - (CONJRu63*Cosal
     +p*ee**2*Ru63*Tanbet*vd)/(3.*Pcw2) + (CONJRu63*Cosalp*ee**2*R
     +u63*Tanbet*vd)/(4.*Pcw2*Psw2) - CONJRu63*Cosalp*Pyuku32*Ru63
     +*Tanbet*vd + (CONJRu66*Cosalp*ee**2*Ru66*Tanbet*vd)/(3.*Pcw2
     +) - CONJRu66*Cosalp*Pyuku32*Ru66*Tanbet*vd - (CONJRu66*MMMU*
     +Ru63*Sinalp*yuku3)/Sqrt2 - (CONJMMMU*CONJRu63*Ru66*Sinalp*yu
     +ku3)/Sqrt2
      MGVX696 = -((CONJRd16*CONJtd33*Cosalp*Rd13)/Sqrt2) - (CONJRd1
     +3*Cosalp*Rd16*td33)/Sqrt2 - (CONJRd13*Cosalp*ee**2*Rd13*vd)/
     +(6.*Pcw2) + (CONJRd13*Cosalp*ee**2*Rd13*vd)/(4.*Pcw2*Psw2) -
     + CONJRd13*Cosalp*Pyukd32*Rd13*vd + (CONJRd16*Cosalp*ee**2*Rd
     +16*vd)/(6.*Pcw2) - CONJRd16*Cosalp*Pyukd32*Rd16*vd + (CONJRd
     +13*ee**2*Rd13*Sinalp*Tanbet*vd)/(6.*Pcw2) - (CONJRd13*ee**2*
     +Rd13*Sinalp*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRd16*ee**2*Rd16
     +*Sinalp*Tanbet*vd)/(6.*Pcw2) + (CONJRd16*MMMU*Rd13*Sinalp*yu
     +kd3)/Sqrt2 + (CONJMMMU*CONJRd13*Rd16*Sinalp*yukd3)/Sqrt2
      MGVX698 = -((CONJRd16*CONJtd33*Cosalp*Rd23)/Sqrt2) - (CONJRd1
     +3*Cosalp*Rd26*td33)/Sqrt2 - (CONJRd13*Cosalp*ee**2*Rd23*vd)/
     +(6.*Pcw2) + (CONJRd13*Cosalp*ee**2*Rd23*vd)/(4.*Pcw2*Psw2) -
     + CONJRd13*Cosalp*Pyukd32*Rd23*vd + (CONJRd16*Cosalp*ee**2*Rd
     +26*vd)/(6.*Pcw2) - CONJRd16*Cosalp*Pyukd32*Rd26*vd + (CONJRd
     +13*ee**2*Rd23*Sinalp*Tanbet*vd)/(6.*Pcw2) - (CONJRd13*ee**2*
     +Rd23*Sinalp*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRd16*ee**2*Rd26
     +*Sinalp*Tanbet*vd)/(6.*Pcw2) + (CONJRd16*MMMU*Rd23*Sinalp*yu
     +kd3)/Sqrt2 + (CONJMMMU*CONJRd13*Rd26*Sinalp*yukd3)/Sqrt2
      MGVX697 = -((CONJRd26*CONJtd33*Cosalp*Rd13)/Sqrt2) - (CONJRd2
     +3*Cosalp*Rd16*td33)/Sqrt2 - (CONJRd23*Cosalp*ee**2*Rd13*vd)/
     +(6.*Pcw2) + (CONJRd23*Cosalp*ee**2*Rd13*vd)/(4.*Pcw2*Psw2) -
     + CONJRd23*Cosalp*Pyukd32*Rd13*vd + (CONJRd26*Cosalp*ee**2*Rd
     +16*vd)/(6.*Pcw2) - CONJRd26*Cosalp*Pyukd32*Rd16*vd + (CONJRd
     +23*ee**2*Rd13*Sinalp*Tanbet*vd)/(6.*Pcw2) - (CONJRd23*ee**2*
     +Rd13*Sinalp*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRd26*ee**2*Rd16
     +*Sinalp*Tanbet*vd)/(6.*Pcw2) + (CONJRd26*MMMU*Rd13*Sinalp*yu
     +kd3)/Sqrt2 + (CONJMMMU*CONJRd23*Rd16*Sinalp*yukd3)/Sqrt2
      MGVX699 = -((CONJRd26*CONJtd33*Cosalp*Rd23)/Sqrt2) - (CONJRd2
     +3*Cosalp*Rd26*td33)/Sqrt2 - (CONJRd23*Cosalp*ee**2*Rd23*vd)/
     +(6.*Pcw2) + (CONJRd23*Cosalp*ee**2*Rd23*vd)/(4.*Pcw2*Psw2) -
     + CONJRd23*Cosalp*Pyukd32*Rd23*vd + (CONJRd26*Cosalp*ee**2*Rd
     +26*vd)/(6.*Pcw2) - CONJRd26*Cosalp*Pyukd32*Rd26*vd + (CONJRd
     +23*ee**2*Rd23*Sinalp*Tanbet*vd)/(6.*Pcw2) - (CONJRd23*ee**2*
     +Rd23*Sinalp*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRd26*ee**2*Rd26
     +*Sinalp*Tanbet*vd)/(6.*Pcw2) + (CONJRd26*MMMU*Rd23*Sinalp*yu
     +kd3)/Sqrt2 + (CONJMMMU*CONJRd23*Rd26*Sinalp*yukd3)/Sqrt2
      MGVX754 = (CONJRd35*Cosalp*ee**2*Rd35*vd)/(6.*Pcw2) - (CONJRd
     +35*ee**2*Rd35*Sinalp*Tanbet*vd)/(6.*Pcw2)
      MGVX755 = (CONJRd44*Cosalp*ee**2*Rd44*vd)/(6.*Pcw2) - (CONJRd
     +44*ee**2*Rd44*Sinalp*Tanbet*vd)/(6.*Pcw2)
      MGVX756 = -(CONJRd51*Cosalp*ee**2*Rd51*vd)/(6.*Pcw2) + (CONJR
     +d51*Cosalp*ee**2*Rd51*vd)/(4.*Pcw2*Psw2) + (CONJRd51*ee**2*R
     +d51*Sinalp*Tanbet*vd)/(6.*Pcw2) - (CONJRd51*ee**2*Rd51*Sinal
     +p*Tanbet*vd)/(4.*Pcw2*Psw2)

      endif

      return
      end
