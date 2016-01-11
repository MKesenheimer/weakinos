c############### cgengrids_add.h #######################################
c last modified by MK, 08.12.2015

c definitions for integrating the on-shell resonant contributions
        double precision xgridosres(0:nintervals,ndiminteg)
        double precision xgrid0osres(0:nintervals,ndiminteg)
        double precision ymaxosres(nintervals,ndiminteg)
        double precision ymaxratosres(nintervals,ndiminteg)
        double precision xmmmosres(0:nintervals,ndiminteg)
        double precision xintosres,xaccosres(0:nintervals,ndiminteg)
        integer ifoldosres(ndiminteg)
        integer nhitsosres(1:nintervals,ndiminteg)
        common/cgengrids_add/ xgridosres,xgrid0osres,ymaxosres
        common/cgengrids_add/ ymaxratosres,xmmmosres,xintosres
        common/cgengrids_add/ xaccosres
        common/cgengrids_add/ ifoldosres,nhitsosres
        save /cgengrids_add/


c############### end cgengrids_add.h ###################################