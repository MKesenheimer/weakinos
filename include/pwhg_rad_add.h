c############### pwhg_add_rad.h ########################################
c last modified by MK, date 05.12.2015
c additional flavor list for on-shell resonant diagrams

c definitions
        ! results of on-shell resonant contributions for radiation phase
        double precision rad_osres_arr(maxprocreal,cnosres)
        ! we want to work with positive contributions, so save the sign
        ! in an additional array
        integer rad_osres_sign(maxprocreal,cnosres)
        integer rad_reg_sign(maxprocreal)
        integer rad_damp_rem_sign(maxalr)
        integer rad_realosres

        ! introduce this array to save the rnd-numbers
        ! for the generic phase-spaces:
        double precision rad_xradosres(ndiminteg)

        ! on-shell contributions to radiation (with errors)
        double precision rad_totosres(cnosres), rad_etotosres(cnosres)
        double precision rad_totabsosres(cnosres)
        double precision rad_etotabsosres(cnosres)
        double precision rad_totpososres(cnosres)
        double precision rad_etotpososres(cnosres)
        double precision rad_totnegosres(cnosres)
        double precision rad_etotnegosres(cnosres)
        double precision rad_totosresgen(cnosres)
        double precision rad_etotosresgen(cnosres)

        ! meaning is the same as for btilde
        double precision rad_totrem, rad_etotrem
        double precision rad_totposrem, rad_etotposrem
        double precision rad_totnegrem, rad_etotnegrem
        double precision rad_totabsrem, rad_etotabsrem
        double precision rad_totreg, rad_etotreg
        double precision rad_totposreg, rad_etotposreg
        double precision rad_totnegreg, rad_etotnegreg
        double precision rad_totabsreg, rad_etotabsreg

        ! similar to btlgen, the following variables are set either
        ! to totabs...(if neg. weights are kept) or to totreg....
        double precision rad_totreggen,rad_etotreggen
        double precision rad_totremgen,rad_etotremgen

        ! save all contributions into two arrays which then can be
        ! accessed simpler and faster
        ! Note: using equivalence statements is deprecated and should
        ! not be used.
        integer ntot
        parameter (ntot = 17)
        double precision rad_totarr(2,ntot)

        integer ntot_osres
        parameter (ntot_osres = 5)
        double precision rad_osres_totarr(2,ntot_osres,cnosres)

        ! common blocks
        ! doubles
        common/pwhg_add_rad/rad_osres_arr
        common/pwhg_add_rad/rad_xradosres
        common/pwhg_add_rad/rad_totosres, rad_etotosres
        common/pwhg_add_rad/rad_totabsosres, rad_etotabsosres
        common/pwhg_add_rad/rad_totpososres, rad_etotpososres
        common/pwhg_add_rad/rad_totnegosres, rad_etotnegosres
        common/pwhg_add_rad/rad_totosresgen, rad_etotosresgen

        common/pwhg_add_rad/rad_totrem, rad_etotrem
        common/pwhg_add_rad/rad_totposrem, rad_etotposrem
        common/pwhg_add_rad/rad_totnegrem, rad_etotnegrem
        common/pwhg_add_rad/rad_totabsrem, rad_etotabsrem
        common/pwhg_add_rad/rad_totreg, rad_etotreg
        common/pwhg_add_rad/rad_totposreg, rad_etotposreg
        common/pwhg_add_rad/rad_totnegreg, rad_etotnegreg
        common/pwhg_add_rad/rad_totabsreg, rad_etotabsreg

        common/pwhg_add_rad/rad_totreggen, rad_etotreggen
        common/pwhg_add_rad/rad_totremgen, rad_etotremgen

        common/pwhg_add_rad/rad_totarr, rad_osres_totarr
        
        ! integers
        common/pwhg_add_rad/rad_osres_sign
        common/pwhg_add_rad/rad_reg_sign
        common/pwhg_add_rad/rad_damp_rem_sign
        common/pwhg_add_rad/rad_realosres

c############### end pwhg_add_rad.h ####################################