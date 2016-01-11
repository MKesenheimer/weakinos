      subroutine coup29(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 29
      if(readlha) then

      MGVX757 = -(CONJRd62*Cosalp*ee**2*Rd62*vd)/(6.*Pcw2) + (CONJR
     +d62*Cosalp*ee**2*Rd62*vd)/(4.*Pcw2*Psw2) + (CONJRd62*ee**2*R
     +d62*Sinalp*Tanbet*vd)/(6.*Pcw2) - (CONJRd62*ee**2*Rd62*Sinal
     +p*Tanbet*vd)/(4.*Pcw2*Psw2)
      MGVX751 = (3*ee**2*Sinalp*Tanbet*vd*Cos(alp)**2)/(4.*Pcw2*Psw
     +2) - (3*ee**2*vd*Cos(alp)**3)/(4.*Pcw2*Psw2) + (3*Cosalp*ee*
     +*2*vd*Sin(alp)**2)/(4.*Pcw2*Psw2) - (3*ee**2*Tanbet*vd*Sin(a
     +lp)**3)/(4.*Pcw2*Psw2)
      MGVX725 = -((CONJRl16*CONJte33*Cosalp*Rl13)/Sqrt2) - (CONJRl1
     +3*Cosalp*Rl16*te33)/Sqrt2 - (CONJRl13*Cosalp*ee**2*Rl13*vd)/
     +(2.*Pcw2) + (CONJRl13*Cosalp*ee**2*Rl13*vd)/(4.*Pcw2*Psw2) -
     + CONJRl13*Cosalp*Pyukl32*Rl13*vd + (CONJRl16*Cosalp*ee**2*Rl
     +16*vd)/(2.*Pcw2) - CONJRl16*Cosalp*Pyukl32*Rl16*vd + (CONJRl
     +13*ee**2*Rl13*Sinalp*Tanbet*vd)/(2.*Pcw2) - (CONJRl13*ee**2*
     +Rl13*Sinalp*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRl16*ee**2*Rl16
     +*Sinalp*Tanbet*vd)/(2.*Pcw2) + (CONJRl16*MMMU*Rl13*Sinalp*yu
     +kl3)/Sqrt2 + (CONJMMMU*CONJRl13*Rl16*Sinalp*yukl3)/Sqrt2
      MGVX727 = -((CONJRl16*CONJte33*Cosalp*Rl63)/Sqrt2) - (CONJRl1
     +3*Cosalp*Rl66*te33)/Sqrt2 - (CONJRl13*Cosalp*ee**2*Rl63*vd)/
     +(2.*Pcw2) + (CONJRl13*Cosalp*ee**2*Rl63*vd)/(4.*Pcw2*Psw2) -
     + CONJRl13*Cosalp*Pyukl32*Rl63*vd + (CONJRl16*Cosalp*ee**2*Rl
     +66*vd)/(2.*Pcw2) - CONJRl16*Cosalp*Pyukl32*Rl66*vd + (CONJRl
     +13*ee**2*Rl63*Sinalp*Tanbet*vd)/(2.*Pcw2) - (CONJRl13*ee**2*
     +Rl63*Sinalp*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRl16*ee**2*Rl66
     +*Sinalp*Tanbet*vd)/(2.*Pcw2) + (CONJRl16*MMMU*Rl63*Sinalp*yu
     +kl3)/Sqrt2 + (CONJMMMU*CONJRl13*Rl66*Sinalp*yukl3)/Sqrt2
      MGVX726 = -((CONJRl66*CONJte33*Cosalp*Rl13)/Sqrt2) - (CONJRl6
     +3*Cosalp*Rl16*te33)/Sqrt2 - (CONJRl63*Cosalp*ee**2*Rl13*vd)/
     +(2.*Pcw2) + (CONJRl63*Cosalp*ee**2*Rl13*vd)/(4.*Pcw2*Psw2) -
     + CONJRl63*Cosalp*Pyukl32*Rl13*vd + (CONJRl66*Cosalp*ee**2*Rl
     +16*vd)/(2.*Pcw2) - CONJRl66*Cosalp*Pyukl32*Rl16*vd + (CONJRl
     +63*ee**2*Rl13*Sinalp*Tanbet*vd)/(2.*Pcw2) - (CONJRl63*ee**2*
     +Rl13*Sinalp*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRl66*ee**2*Rl16
     +*Sinalp*Tanbet*vd)/(2.*Pcw2) + (CONJRl66*MMMU*Rl13*Sinalp*yu
     +kl3)/Sqrt2 + (CONJMMMU*CONJRl63*Rl16*Sinalp*yukl3)/Sqrt2
      MGVX758 = (CONJRl24*Cosalp*ee**2*Rl24*vd)/(2.*Pcw2) - (CONJRl
     +24*ee**2*Rl24*Sinalp*Tanbet*vd)/(2.*Pcw2)
      MGVX759 = (CONJRl35*Cosalp*ee**2*Rl35*vd)/(2.*Pcw2) - (CONJRl
     +35*ee**2*Rl35*Sinalp*Tanbet*vd)/(2.*Pcw2)
      MGVX760 = -(CONJRl42*Cosalp*ee**2*Rl42*vd)/(2.*Pcw2) + (CONJR
     +l42*Cosalp*ee**2*Rl42*vd)/(4.*Pcw2*Psw2) + (CONJRl42*ee**2*R
     +l42*Sinalp*Tanbet*vd)/(2.*Pcw2) - (CONJRl42*ee**2*Rl42*Sinal
     +p*Tanbet*vd)/(4.*Pcw2*Psw2)
      MGVX761 = -(CONJRl51*Cosalp*ee**2*Rl51*vd)/(2.*Pcw2) + (CONJR
     +l51*Cosalp*ee**2*Rl51*vd)/(4.*Pcw2*Psw2) + (CONJRl51*ee**2*R
     +l51*Sinalp*Tanbet*vd)/(2.*Pcw2) - (CONJRl51*ee**2*Rl51*Sinal
     +p*Tanbet*vd)/(4.*Pcw2*Psw2)
      MGVX728 = -((CONJRl66*CONJte33*Cosalp*Rl63)/Sqrt2) - (CONJRl6
     +3*Cosalp*Rl66*te33)/Sqrt2 - (CONJRl63*Cosalp*ee**2*Rl63*vd)/
     +(2.*Pcw2) + (CONJRl63*Cosalp*ee**2*Rl63*vd)/(4.*Pcw2*Psw2) -
     + CONJRl63*Cosalp*Pyukl32*Rl63*vd + (CONJRl66*Cosalp*ee**2*Rl
     +66*vd)/(2.*Pcw2) - CONJRl66*Cosalp*Pyukl32*Rl66*vd + (CONJRl
     +63*ee**2*Rl63*Sinalp*Tanbet*vd)/(2.*Pcw2) - (CONJRl63*ee**2*
     +Rl63*Sinalp*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRl66*ee**2*Rl66
     +*Sinalp*Tanbet*vd)/(2.*Pcw2) + (CONJRl66*MMMU*Rl63*Sinalp*yu
     +kl3)/Sqrt2 + (CONJMMMU*CONJRl63*Rl66*Sinalp*yukl3)/Sqrt2
      MGVX747 = -(Cosalp*ee**2*vd)/(4.*Pcw2*Psw2) + (ee**2*Sinalp*T
     +anbet*vd)/(4.*Pcw2*Psw2)
      MGVX735 = -((CONJRu16*CONJtu33*Ru13*Sinalp)/Sqrt2) - (CONJRu1
     +3*Ru16*Sinalp*tu33)/Sqrt2 + (CONJRu13*Cosalp*ee**2*Ru13*vd)/
     +(3.*Pcw2) - (CONJRu13*Cosalp*ee**2*Ru13*vd)/(4.*Pcw2*Psw2) -
     + (CONJRu16*Cosalp*ee**2*Ru16*vd)/(3.*Pcw2) - (CONJRu13*ee**2
     +*Ru13*Sinalp*Tanbet*vd)/(3.*Pcw2) + (CONJRu13*ee**2*Ru13*Sin
     +alp*Tanbet*vd)/(4.*Pcw2*Psw2) - CONJRu13*Pyuku32*Ru13*Sinalp
     +*Tanbet*vd + (CONJRu16*ee**2*Ru16*Sinalp*Tanbet*vd)/(3.*Pcw2
     +) - CONJRu16*Pyuku32*Ru16*Sinalp*Tanbet*vd + (CONJRu16*Cosal
     +p*MMMU*Ru13*yuku3)/Sqrt2 + (CONJMMMU*CONJRu13*Cosalp*Ru16*yu
     +ku3)/Sqrt2
      MGVX741 = -((CONJRu16*CONJtu33*Ru63*Sinalp)/Sqrt2) - (CONJRu1
     +3*Ru66*Sinalp*tu33)/Sqrt2 + (CONJRu13*Cosalp*ee**2*Ru63*vd)/
     +(3.*Pcw2) - (CONJRu13*Cosalp*ee**2*Ru63*vd)/(4.*Pcw2*Psw2) -
     + (CONJRu16*Cosalp*ee**2*Ru66*vd)/(3.*Pcw2) - (CONJRu13*ee**2
     +*Ru63*Sinalp*Tanbet*vd)/(3.*Pcw2) + (CONJRu13*ee**2*Ru63*Sin
     +alp*Tanbet*vd)/(4.*Pcw2*Psw2) - CONJRu13*Pyuku32*Ru63*Sinalp
     +*Tanbet*vd + (CONJRu16*ee**2*Ru66*Sinalp*Tanbet*vd)/(3.*Pcw2
     +) - CONJRu16*Pyuku32*Ru66*Sinalp*Tanbet*vd + (CONJRu16*Cosal
     +p*MMMU*Ru63*yuku3)/Sqrt2 + (CONJMMMU*CONJRu13*Cosalp*Ru66*yu
     +ku3)/Sqrt2
      MGVX738 = -((CONJRu66*CONJtu33*Ru13*Sinalp)/Sqrt2) - (CONJRu6
     +3*Ru16*Sinalp*tu33)/Sqrt2 + (CONJRu63*Cosalp*ee**2*Ru13*vd)/
     +(3.*Pcw2) - (CONJRu63*Cosalp*ee**2*Ru13*vd)/(4.*Pcw2*Psw2) -
     + (CONJRu66*Cosalp*ee**2*Ru16*vd)/(3.*Pcw2) - (CONJRu63*ee**2
     +*Ru13*Sinalp*Tanbet*vd)/(3.*Pcw2) + (CONJRu63*ee**2*Ru13*Sin
     +alp*Tanbet*vd)/(4.*Pcw2*Psw2) - CONJRu63*Pyuku32*Ru13*Sinalp
     +*Tanbet*vd + (CONJRu66*ee**2*Ru16*Sinalp*Tanbet*vd)/(3.*Pcw2
     +) - CONJRu66*Pyuku32*Ru16*Sinalp*Tanbet*vd + (CONJRu66*Cosal
     +p*MMMU*Ru13*yuku3)/Sqrt2 + (CONJMMMU*CONJRu63*Cosalp*Ru16*yu
     +ku3)/Sqrt2
      MGVX777 = -(CONJRu24*Cosalp*ee**2*Ru24*vd)/(3.*Pcw2) + (CONJR
     +u24*ee**2*Ru24*Sinalp*Tanbet*vd)/(3.*Pcw2)

      endif

      return
      end
