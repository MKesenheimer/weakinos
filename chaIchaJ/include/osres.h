c############### osres.h ###############################################
c last modified by MK, date 05.12.2015
c global definitions for treating the on-shell resonant diagrams

c for the case of weakino pair-production, there are four resonances
c but if you have less resonances you can simply initialize the variables
c defined here in a different way.
c                        mi
c                        / 
c   ____      !         /
c  |    |  sij=mij^2   /
c  | s  |--------------
c  |____|\             \
c         \             \
c          \             \
c          mk            mj

c definitions
        ! store the number of resonances here
        integer nosres, cnosres ! constant
        ! dummy variable to initialize the array
        parameter (cnosres=4)
        ! this is the human readable array which stores the id of the 
        ! on-shell resonance, e.g. "ul35", "dl35".
        ! this variable is set in init_processes
        character*4 osresID(cnosres)
        !data osresID/"dl35","dl45","ul35","ul45"/

        ! the mass array of the on-shell resonant particles
        double precision osresM(cnosres)
        ! the width of the on-shell resonant particle (regulator)
        ! this must be set in init_couplings.f
        double precision wreg

        ! store in an additional array which legs get resonant
        ! e.g. osreslegs(3,1) = 3 and osreslegs(3,2) = 5 means that for 
        ! resonance with ID osresID = 3 (dl35) the resonant particle 
        ! (left handed down-type squark) decays into the final state 
        ! particles with numbers 3 & 5
        ! osreslegs(chan,1): first resonant particle (index is i used in our code)
        ! osreslegs(chan,2): second resonant particle (index is j used in our code)
        ! osreslegs(chan,3): spectator (index is k used in our code) 
        integer osreslegs(cnosres,1:3)

        common/c_onshell/ osresM, osreslegs, wreg
        common/c_onshell/ osresID
        common/c_onshell/ nosres


#if defined(DSUB_I) || defined(DSUB_II) || defined(DSUB_II_TEST)
#define WDLR wreg
#define WDRR wreg
#define WULR wreg
#define WURR wreg

#else
#define WDLR WDL
#define WDRR WDR
#define WULR WUL
#define WURR WUR
#endif

c############### end osres.h ###########################################
