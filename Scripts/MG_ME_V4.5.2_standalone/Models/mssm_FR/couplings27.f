      subroutine coup27(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 27
      if(readlha) then

      MGVX770 = -(CONJRd44*ee**2*Rd44*Sinalp*vd)/(6.*Pcw2) - (CONJR
     +d44*Cosalp*ee**2*Rd44*Tanbet*vd)/(6.*Pcw2)
      MGVX771 = (CONJRd51*ee**2*Rd51*Sinalp*vd)/(6.*Pcw2) - (CONJRd
     +51*ee**2*Rd51*Sinalp*vd)/(4.*Pcw2*Psw2) + (CONJRd51*Cosalp*e
     +e**2*Rd51*Tanbet*vd)/(6.*Pcw2) - (CONJRd51*Cosalp*ee**2*Rd51
     +*Tanbet*vd)/(4.*Pcw2*Psw2)
      MGVX772 = (CONJRd62*ee**2*Rd62*Sinalp*vd)/(6.*Pcw2) - (CONJRd
     +62*ee**2*Rd62*Sinalp*vd)/(4.*Pcw2*Psw2) + (CONJRd62*Cosalp*e
     +e**2*Rd62*Tanbet*vd)/(6.*Pcw2) - (CONJRd62*Cosalp*ee**2*Rd62
     +*Tanbet*vd)/(4.*Pcw2*Psw2)
      MGVX765 = (-3*ee**2*Sinalp*vd*Cos(alp)**2)/(4.*Pcw2*Psw2) - (
     +3*ee**2*Tanbet*vd*Cos(alp)**3)/(4.*Pcw2*Psw2) + (3*Cosalp*ee
     +**2*Tanbet*vd*Sin(alp)**2)/(4.*Pcw2*Psw2) + (3*ee**2*vd*Sin(
     +alp)**3)/(4.*Pcw2*Psw2)
      MGVX750 = (-5*ee**2*Sinalp*Tanbet*vd*Cos(alp)**2)/(4.*Pcw2*Ps
     +w2) + (ee**2*vd*Cos(alp)**3)/(4.*Pcw2*Psw2) - (5*Cosalp*ee**
     +2*vd*Sin(alp)**2)/(4.*Pcw2*Psw2) + (ee**2*Tanbet*vd*Sin(alp)
     +**3)/(4.*Pcw2*Psw2)
      MGVX766 = (5*ee**2*Sinalp*vd*Cos(alp)**2)/(4.*Pcw2*Psw2) + (e
     +e**2*Tanbet*vd*Cos(alp)**3)/(4.*Pcw2*Psw2) - (5*Cosalp*ee**2
     +*Tanbet*vd*Sin(alp)**2)/(4.*Pcw2*Psw2) - (ee**2*vd*Sin(alp)*
     +*3)/(4.*Pcw2*Psw2)
      MGVX712 = (CONJRl16*CONJte33*Rl13*Sinalp)/Sqrt2 + (CONJRl13*R
     +l16*Sinalp*te33)/Sqrt2 + (CONJRl13*ee**2*Rl13*Sinalp*vd)/(2.
     +*Pcw2) - (CONJRl13*ee**2*Rl13*Sinalp*vd)/(4.*Pcw2*Psw2) + CO
     +NJRl13*Pyukl32*Rl13*Sinalp*vd - (CONJRl16*ee**2*Rl16*Sinalp*
     +vd)/(2.*Pcw2) + CONJRl16*Pyukl32*Rl16*Sinalp*vd + (CONJRl13*
     +Cosalp*ee**2*Rl13*Tanbet*vd)/(2.*Pcw2) - (CONJRl13*Cosalp*ee
     +**2*Rl13*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRl16*Cosalp*ee**2*
     +Rl16*Tanbet*vd)/(2.*Pcw2) + (CONJRl16*Cosalp*MMMU*Rl13*yukl3
     +)/Sqrt2 + (CONJMMMU*CONJRl13*Cosalp*Rl16*yukl3)/Sqrt2
      MGVX717 = (CONJRl16*CONJte33*Rl63*Sinalp)/Sqrt2 + (CONJRl13*R
     +l66*Sinalp*te33)/Sqrt2 + (CONJRl13*ee**2*Rl63*Sinalp*vd)/(2.
     +*Pcw2) - (CONJRl13*ee**2*Rl63*Sinalp*vd)/(4.*Pcw2*Psw2) + CO
     +NJRl13*Pyukl32*Rl63*Sinalp*vd - (CONJRl16*ee**2*Rl66*Sinalp*
     +vd)/(2.*Pcw2) + CONJRl16*Pyukl32*Rl66*Sinalp*vd + (CONJRl13*
     +Cosalp*ee**2*Rl63*Tanbet*vd)/(2.*Pcw2) - (CONJRl13*Cosalp*ee
     +**2*Rl63*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRl16*Cosalp*ee**2*
     +Rl66*Tanbet*vd)/(2.*Pcw2) + (CONJRl16*Cosalp*MMMU*Rl63*yukl3
     +)/Sqrt2 + (CONJMMMU*CONJRl13*Cosalp*Rl66*yukl3)/Sqrt2
      MGVX713 = (CONJRl66*CONJte33*Rl13*Sinalp)/Sqrt2 + (CONJRl63*R
     +l16*Sinalp*te33)/Sqrt2 + (CONJRl63*ee**2*Rl13*Sinalp*vd)/(2.
     +*Pcw2) - (CONJRl63*ee**2*Rl13*Sinalp*vd)/(4.*Pcw2*Psw2) + CO
     +NJRl63*Pyukl32*Rl13*Sinalp*vd - (CONJRl66*ee**2*Rl16*Sinalp*
     +vd)/(2.*Pcw2) + CONJRl66*Pyukl32*Rl16*Sinalp*vd + (CONJRl63*
     +Cosalp*ee**2*Rl13*Tanbet*vd)/(2.*Pcw2) - (CONJRl63*Cosalp*ee
     +**2*Rl13*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRl66*Cosalp*ee**2*
     +Rl16*Tanbet*vd)/(2.*Pcw2) + (CONJRl66*Cosalp*MMMU*Rl13*yukl3
     +)/Sqrt2 + (CONJMMMU*CONJRl63*Cosalp*Rl16*yukl3)/Sqrt2
      MGVX773 = -(CONJRl24*ee**2*Rl24*Sinalp*vd)/(2.*Pcw2) - (CONJR
     +l24*Cosalp*ee**2*Rl24*Tanbet*vd)/(2.*Pcw2)
      MGVX774 = -(CONJRl35*ee**2*Rl35*Sinalp*vd)/(2.*Pcw2) - (CONJR
     +l35*Cosalp*ee**2*Rl35*Tanbet*vd)/(2.*Pcw2)
      MGVX775 = (CONJRl42*ee**2*Rl42*Sinalp*vd)/(2.*Pcw2) - (CONJRl
     +42*ee**2*Rl42*Sinalp*vd)/(4.*Pcw2*Psw2) + (CONJRl42*Cosalp*e
     +e**2*Rl42*Tanbet*vd)/(2.*Pcw2) - (CONJRl42*Cosalp*ee**2*Rl42
     +*Tanbet*vd)/(4.*Pcw2*Psw2)
      MGVX776 = (CONJRl51*ee**2*Rl51*Sinalp*vd)/(2.*Pcw2) - (CONJRl
     +51*ee**2*Rl51*Sinalp*vd)/(4.*Pcw2*Psw2) + (CONJRl51*Cosalp*e
     +e**2*Rl51*Tanbet*vd)/(2.*Pcw2) - (CONJRl51*Cosalp*ee**2*Rl51
     +*Tanbet*vd)/(4.*Pcw2*Psw2)
      MGVX718 = (CONJRl66*CONJte33*Rl63*Sinalp)/Sqrt2 + (CONJRl63*R
     +l66*Sinalp*te33)/Sqrt2 + (CONJRl63*ee**2*Rl63*Sinalp*vd)/(2.
     +*Pcw2) - (CONJRl63*ee**2*Rl63*Sinalp*vd)/(4.*Pcw2*Psw2) + CO
     +NJRl63*Pyukl32*Rl63*Sinalp*vd - (CONJRl66*ee**2*Rl66*Sinalp*
     +vd)/(2.*Pcw2) + CONJRl66*Pyukl32*Rl66*Sinalp*vd + (CONJRl63*
     +Cosalp*ee**2*Rl63*Tanbet*vd)/(2.*Pcw2) - (CONJRl63*Cosalp*ee
     +**2*Rl63*Tanbet*vd)/(4.*Pcw2*Psw2) - (CONJRl66*Cosalp*ee**2*
     +Rl66*Tanbet*vd)/(2.*Pcw2) + (CONJRl66*Cosalp*MMMU*Rl63*yukl3
     +)/Sqrt2 + (CONJMMMU*CONJRl63*Cosalp*Rl66*yukl3)/Sqrt2
      MGVX762 = (ee**2*Sinalp*vd)/(4.*Pcw2*Psw2) + (Cosalp*ee**2*Ta
     +nbet*vd)/(4.*Pcw2*Psw2)

      endif

      return
      end
