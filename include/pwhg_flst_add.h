c############### pwhg_add_flst.h #######################################
c last modified by MK, date 05.12.2015
c additional flavor list for on-shell resonant diagrams

c definitions
        ! flavor list of on-shell resonant diagrams
        integer flst_osres(nlegreal,maxprocreal)
        ! number of on-shell resonant diagrams
        integer flst_nosres

        common/pwhg_add_flst/ flst_osres,flst_nosres


c############### end pwhg_add_flst.h ###################################