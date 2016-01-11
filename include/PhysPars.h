c############### PhysPars.h ############################################
c last modified by MK, 18.12.2014

c read from slha files:
#include "SLHA.h"
c formcalc model definitions
#include "model_mssm.h"
#include "model_sm.h"
c madgraph definitions
#include "coupl.inc"
#include "sm_read_values.inc"

c common definitions
c SM parameters needed for neutralino pair production
        double precision  par_alpha, par_gf
        double precision  par_Zmass,par_Wmass,par_Zwidth,par_Wwidth
        double precision  par_Zmass2,par_Wmass2
        double precision  par_MU,par_MD,par_MS,par_MC,par_MB,par_MT
        
        common/par_common/ 
     &        par_alpha,par_gf,                ! fine structure and fermi constant
     &        par_Zmass,par_Zwidth,            ! mass of Z boson
     &        par_Wmass,par_Wwidth,            ! mass of W boson
     &        par_Zmass2,par_Wmass2,           ! squared masses
     &        par_MU,par_MD,par_MS,
     &        par_MC,par_MB,par_MT

c MSSM parameters needed for neutralino pair production
c parameters will be read from SLHA file in init_couplings.F
        double precision par_Fin1mass,par_Fin2mass
        double precision par_MSf(2,4,3),par_MNeu(4),par_MCha(2)
        double precision par_MGl
        double precision par_Sfwidth(2,4,3)
        double precision par_Glwidth
        common /par_process/  
     &          par_Fin1mass,par_Fin2mass,
     &          par_MSf,par_MNeu,par_MCha,
     &          par_MGl,par_Sfwidth,par_Glwidth

c MSSM parameters in SLHALib
        double complex slhadata(nslhadata)
        ! decay widths
        integer partid(0:7,1024)
        double precision width(1024)
        integer nparticles, nchannels
        common /par_slha/ slhadata, partid, width, nparticles, nchannels

c complex unit
        double complex ii
        parameter (ii = (0D0,1D0))

c additional parameters


c############### end PhysPars.h ########################################