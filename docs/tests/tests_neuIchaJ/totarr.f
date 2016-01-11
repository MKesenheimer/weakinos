        equivalence (rad_osres_totarr(1,1,1), rad_totosres(1))
        equivalence (rad_osres_totarr(2,1,1), rad_etotosres(1))
        equivalence (rad_osres_totarr(1,2,1), rad_totabsosres(1))
        equivalence (rad_osres_totarr(2,2,1), rad_etotabsosres(1))
        equivalence (rad_osres_totarr(1,3,1), rad_totpososres(1))
        equivalence (rad_osres_totarr(2,3,1), rad_etotpososres(1))
        equivalence (rad_osres_totarr(1,4,1), rad_totnegosres(1))
        equivalence (rad_osres_totarr(2,4,1), rad_etotnegosres(1))
        equivalence (rad_osres_totarr(1,5,1), rad_totosresgen(1))
        equivalence (rad_osres_totarr(1,5,1), rad_etotosresgen(1))

        equivalence (rad_osres_totarr(1,1,2), rad_totosres(2))
        equivalence (rad_osres_totarr(2,1,2), rad_etotosres(2))
        equivalence (rad_osres_totarr(1,2,2), rad_totabsosres(2))
        equivalence (rad_osres_totarr(2,2,2), rad_etotabsosres(2))
        equivalence (rad_osres_totarr(1,3,2), rad_totpososres(2))
        equivalence (rad_osres_totarr(2,3,2), rad_etotpososres(2))
        equivalence (rad_osres_totarr(1,4,2), rad_totnegosres(2))
        equivalence (rad_osres_totarr(2,4,2), rad_etotnegosres(2))
        equivalence (rad_osres_totarr(1,5,2), rad_totosresgen(2))
        equivalence (rad_osres_totarr(1,5,2), rad_etotosresgen(2))

        equivalence (rad_osres_totarr(1,1,3), rad_totosres(3))
        equivalence (rad_osres_totarr(2,1,3), rad_etotosres(3))
        equivalence (rad_osres_totarr(1,2,3), rad_totabsosres(3))
        equivalence (rad_osres_totarr(2,2,3), rad_etotabsosres(3))
        equivalence (rad_osres_totarr(1,3,3), rad_totpososres(3))
        equivalence (rad_osres_totarr(2,3,3), rad_etotpososres(3))
        equivalence (rad_osres_totarr(1,4,3), rad_totnegosres(3))
        equivalence (rad_osres_totarr(2,4,3), rad_etotnegosres(3))
        equivalence (rad_osres_totarr(1,5,3), rad_totosresgen(3))
        equivalence (rad_osres_totarr(1,5,3), rad_etotosresgen(3))

        equivalence (rad_osres_totarr(1,1,4), rad_totosres(4))
        equivalence (rad_osres_totarr(2,1,4), rad_etotosres(4))
        equivalence (rad_osres_totarr(1,2,4), rad_totabsosres(4))
        equivalence (rad_osres_totarr(2,2,4), rad_etotabsosres(4))
        equivalence (rad_osres_totarr(1,3,4), rad_totpososres(4))
        equivalence (rad_osres_totarr(2,3,4), rad_etotpososres(4))
        equivalence (rad_osres_totarr(1,4,4), rad_totnegosres(4))
        equivalence (rad_osres_totarr(2,4,4), rad_etotnegosres(4))
        equivalence (rad_osres_totarr(1,5,4), rad_totosresgen(4))
        equivalence (rad_osres_totarr(1,5,4), rad_etotosresgen(4))



oder


        ! TEST
        ! equivalence statements
        ! btilde contributions 1:4
        equivalence (rad_totarr(1,1), rad_totbtl)
        equivalence (rad_totarr(2,1), rad_etotbtl)
        equivalence (rad_totarr(1,2), rad_totabsbtl)
        equivalence (rad_totarr(2,2), rad_etotabsbtl)
        equivalence (rad_totarr(1,3), rad_totposbtl)
        equivalence (rad_totarr(2,3), rad_etotposbtl)
        equivalence (rad_totarr(1,4), rad_totnegbtl)
        equivalence (rad_totarr(2,4), rad_etotnegbtl)

        ! regular contributions 5:8
        equivalence (rad_totarr(1,5), rad_totreg)
        equivalence (rad_totarr(2,5), rad_etotreg)
        equivalence (rad_totarr(1,6), rad_totabsreg)
        equivalence (rad_totarr(2,6), rad_etotabsreg)
        equivalence (rad_totarr(1,7), rad_totposreg)
        equivalence (rad_totarr(2,7), rad_etotposreg)
        equivalence (rad_totarr(1,8), rad_totnegreg)
        equivalence (rad_totarr(2,8), rad_etotnegreg)

        ! remnant contributions 9:15
        equivalence (rad_totarr(1,9) , rad_totrem)
        equivalence (rad_totarr(2,9) , rad_etotrem)
        equivalence (rad_totarr(1,10), rad_totabsrem)
        equivalence (rad_totarr(2,10), rad_etotabsrem)
        equivalence (rad_totarr(1,11), rad_totposrem)
        equivalence (rad_totarr(2,11), rad_etotposrem)
        equivalence (rad_totarr(1,12), rad_totnegrem)
        equivalence (rad_totarr(2,12), rad_etotnegrem)

        equivalence (rad_totarr(1,13), rad_totbtlgen)
        equivalence (rad_totarr(2,13), rad_etotbtlgen)
        equivalence (rad_totarr(1,14), rad_totreggen)
        equivalence (rad_totarr(2,14), rad_etotreggen)
        equivalence (rad_totarr(1,15), rad_totremgen)
        equivalence (rad_totarr(2,15), rad_etotremgen)

        ! totals 16:17
        equivalence (rad_totarr(1,16), rad_tot)
        equivalence (rad_totarr(2,16), rad_etot)
        equivalence (rad_totarr(1,17), rad_totgen)
        equivalence (rad_totarr(2,17), rad_etotgen)

        ! on-shell contributions
        equivalence (rad_osres_totarr(1,1,1:cnosres), rad_totosres(1:cnosres))
        equivalence (rad_osres_totarr(2,1,1:cnosres), rad_etotosres(1:cnosres))
        equivalence (rad_osres_totarr(1,2,1:cnosres), rad_totabsosres(1:cnosres))
        equivalence (rad_osres_totarr(2,2,1:cnosres), rad_etotabsosres(1:cnosres))
        equivalence (rad_osres_totarr(1,3,1:cnosres), rad_totpososres(1:cnosres))
        equivalence (rad_osres_totarr(2,3,1:cnosres), rad_etotpososres(1:cnosres))
        equivalence (rad_osres_totarr(1,4,1:cnosres), rad_totnegosres(1:cnosres))
        equivalence (rad_osres_totarr(2,4,1:cnosres), rad_etotnegosres(1:cnosres))
        equivalence (rad_osres_totarr(1,5,1:cnosres), rad_totosresgen(1:cnosres))
        equivalence (rad_osres_totarr(1,5,1:cnosres), rad_etotosresgen(1:cnosres))
