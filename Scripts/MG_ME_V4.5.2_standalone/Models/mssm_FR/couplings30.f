      subroutine coup30(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 30
      if(readlha) then

      MGVX779 = -(CONJRu35*Cosalp*ee**2*Ru35*vd)/(3.*Pcw2) + (CONJR
     +u35*ee**2*Ru35*Sinalp*Tanbet*vd)/(3.*Pcw2)
      MGVX781 = (CONJRu42*Cosalp*ee**2*Ru42*vd)/(3.*Pcw2) - (CONJRu
     +42*Cosalp*ee**2*Ru42*vd)/(4.*Pcw2*Psw2) - (CONJRu42*ee**2*Ru
     +42*Sinalp*Tanbet*vd)/(3.*Pcw2) + (CONJRu42*ee**2*Ru42*Sinalp
     +*Tanbet*vd)/(4.*Pcw2*Psw2)
      MGVX783 = (CONJRu51*Cosalp*ee**2*Ru51*vd)/(3.*Pcw2) - (CONJRu
     +51*Cosalp*ee**2*Ru51*vd)/(4.*Pcw2*Psw2) - (CONJRu51*ee**2*Ru
     +51*Sinalp*Tanbet*vd)/(3.*Pcw2) + (CONJRu51*ee**2*Ru51*Sinalp
     +*Tanbet*vd)/(4.*Pcw2*Psw2)
      MGVX744 = -((CONJRu66*CONJtu33*Ru63*Sinalp)/Sqrt2) - (CONJRu6
     +3*Ru66*Sinalp*tu33)/Sqrt2 + (CONJRu63*Cosalp*ee**2*Ru63*vd)/
     +(3.*Pcw2) - (CONJRu63*Cosalp*ee**2*Ru63*vd)/(4.*Pcw2*Psw2) -
     + (CONJRu66*Cosalp*ee**2*Ru66*vd)/(3.*Pcw2) - (CONJRu63*ee**2
     +*Ru63*Sinalp*Tanbet*vd)/(3.*Pcw2) + (CONJRu63*ee**2*Ru63*Sin
     +alp*Tanbet*vd)/(4.*Pcw2*Psw2) - CONJRu63*Pyuku32*Ru63*Sinalp
     +*Tanbet*vd + (CONJRu66*ee**2*Ru66*Sinalp*Tanbet*vd)/(3.*Pcw2
     +) - CONJRu66*Pyuku32*Ru66*Sinalp*Tanbet*vd + (CONJRu66*Cosal
     +p*MMMU*Ru63*yuku3)/Sqrt2 + (CONJMMMU*CONJRu63*Cosalp*Ru66*yu
     +ku3)/Sqrt2
      MGVX703 = CONJRu16*CONJtu33*Cosbet*Rd13 + CONJRu13*Rd16*Sinbe
     +t*td33 - (CONJRu13*ee**2*Rd13*Sinbet*vd)/(Psw2*Sqrt2) + (CON
     +JRu13*Pyukd32*Rd13*Sinbet*vd)/Sqrt2 + (CONJRu13*Pyuku32*Rd13
     +*Sinbet*vd)/Sqrt2 + CONJMMMU*CONJRu13*Cosbet*Rd16*yukd3 + CO
     +NJRu16*MMMU*Rd13*Sinbet*yuku3 + (CONJRu16*Cosbet*Rd16*vd*yuk
     +d3*yuku3)/Sqrt2 + (CONJRu16*Rd16*Sinbet*Tanbet*vd*yukd3*yuku
     +3)/Sqrt2
      MGVX707 = CONJRu66*CONJtu33*Cosbet*Rd13 + CONJRu63*Rd16*Sinbe
     +t*td33 - (CONJRu63*ee**2*Rd13*Sinbet*vd)/(Psw2*Sqrt2) + (CON
     +JRu63*Pyukd32*Rd13*Sinbet*vd)/Sqrt2 + (CONJRu63*Pyuku32*Rd13
     +*Sinbet*vd)/Sqrt2 + CONJMMMU*CONJRu63*Cosbet*Rd16*yukd3 + CO
     +NJRu66*MMMU*Rd13*Sinbet*yuku3 + (CONJRu66*Cosbet*Rd16*vd*yuk
     +d3*yuku3)/Sqrt2 + (CONJRu66*Rd16*Sinbet*Tanbet*vd*yukd3*yuku
     +3)/Sqrt2
      MGVX704 = CONJRu16*CONJtu33*Cosbet*Rd23 + CONJRu13*Rd26*Sinbe
     +t*td33 - (CONJRu13*ee**2*Rd23*Sinbet*vd)/(Psw2*Sqrt2) + (CON
     +JRu13*Pyukd32*Rd23*Sinbet*vd)/Sqrt2 + (CONJRu13*Pyuku32*Rd23
     +*Sinbet*vd)/Sqrt2 + CONJMMMU*CONJRu13*Cosbet*Rd26*yukd3 + CO
     +NJRu16*MMMU*Rd23*Sinbet*yuku3 + (CONJRu16*Cosbet*Rd26*vd*yuk
     +d3*yuku3)/Sqrt2 + (CONJRu16*Rd26*Sinbet*Tanbet*vd*yukd3*yuku
     +3)/Sqrt2
      MGVX708 = CONJRu66*CONJtu33*Cosbet*Rd23 + CONJRu63*Rd26*Sinbe
     +t*td33 - (CONJRu63*ee**2*Rd23*Sinbet*vd)/(Psw2*Sqrt2) + (CON
     +JRu63*Pyukd32*Rd23*Sinbet*vd)/Sqrt2 + (CONJRu63*Pyuku32*Rd23
     +*Sinbet*vd)/Sqrt2 + CONJMMMU*CONJRu63*Cosbet*Rd26*yukd3 + CO
     +NJRu66*MMMU*Rd23*Sinbet*yuku3 + (CONJRu66*Cosbet*Rd26*vd*yuk
     +d3*yuku3)/Sqrt2 + (CONJRu66*Rd26*Sinbet*Tanbet*vd*yukd3*yuku
     +3)/Sqrt2
      MGVX791 = -((CONJRu51*ee**2*Rd51*Sinbet*vd)/(Psw2*Sqrt2))
      MGVX789 = -((CONJRu42*ee**2*Rd62*Sinbet*vd)/(Psw2*Sqrt2))
      MGVX716 = CONJRLn13*Rl16*Sinbet*te33 - (CONJRLn13*ee**2*Rl13*
     +Sinbet*vd)/(Psw2*Sqrt2) + (CONJRLn13*Pyukl32*Rl13*Sinbet*vd)
     +/Sqrt2 + CONJMMMU*CONJRLn13*Cosbet*Rl16*yukl3
      MGVX785 = -((CONJRLn22*ee**2*Rl42*Sinbet*vd)/(Psw2*Sqrt2))
      MGVX786 = -((CONJRLn31*ee**2*Rl51*Sinbet*vd)/(Psw2*Sqrt2))
      MGVX721 = CONJRLn13*Rl66*Sinbet*te33 - (CONJRLn13*ee**2*Rl63*
     +Sinbet*vd)/(Psw2*Sqrt2) + (CONJRLn13*Pyukl32*Rl63*Sinbet*vd)
     +/Sqrt2 + CONJMMMU*CONJRLn13*Cosbet*Rl66*yukl3
      MGVX194 = -(ee**2*Sinalp*vd)/(2.*Psw2) + (Cosalp*ee**2*Tanbet
     +*vd)/(2.*Psw2)

      endif

      return
      end
