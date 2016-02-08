c############### dsubtraction.h ########################################
c Variables to apply the diagram subtractions scheme

#ifndef DSUBTRACTION_H
#define DSUBTRACTION_H

      character*4 CHAN                       ! distuinguish between u-, t-channels, left and right
      complex*16 JAMPR(1)                    ! sum of resonant diagrams
      double precision MATRIX_RESONANT       ! sum of resonant diagrams squared
      
#if defined(DSUB_I)
      double precision RATIO35L              ! the ratio of the Breit-Wigner functions if the left handed u-channel diagrams get resonant
      double precision RATIO35R              ! "-" right handed u-channel diagrams get resonant
      double precision RATIO45L              ! "-" left handed t-channel diagrams get resonant
      double precision RATIO45R              ! "-" right handed t-channel diagrams get resonant
      double precision THETA35L              ! the theta functions
      double precision THETA35R
      double precision THETA45L
      double precision THETA45R
      double precision COUNTER35L            ! counter term for left handed diagram with resonance at s35 = msq^2
      double precision COUNTER35R            ! counter term for right handed diagram with resonance at s35 = msq^2
      double precision COUNTER45L            ! counter term for left handed diagram with resonance at s45 = msq^2
      double precision COUNTER45R            ! counter term for right handed diagram with resonance at s45 = msq^2
      double precision P_OS(0:3,NEXTERNAL)   ! remapped momenta to on-shell kinematic
      double precision S,S35,S45             ! invariants
      double precision  momsum2sq            ! external function to calculate (pi+pj)^2
      external momsum2sq
#endif

#if defined(DSUB_I) || defined(DSUB_II) || defined(DSUB_II_TEST)
#ifndef _EXTFUNCT_               
      double precision MATRIX_DG_NINJD_RES   ! gives back matrix element squared 
      external MATRIX_DG_NINJD_RES           ! from resonant diagrams in  
      double precision MATRIX_DXG_NINJDX_RES ! dependence on CHAN:
      external MATRIX_DXG_NINJDX_RES         ! u channel diagrams (s35 = msq^2), CHAN=0: left, CHAN=1: right
      double precision MATRIX_GD_NINJD_RES   ! t channel diagrams (s45 = msq^2), CHAN=2: left, CHAN=3: right
      external MATRIX_GD_NINJD_RES
      double precision MATRIX_GDX_NINJDX_RES 
      external MATRIX_GDX_NINJDX_RES
      double precision MATRIX_GU_NINJU_RES
      external MATRIX_GU_NINJU_RES
      double precision MATRIX_GUX_NINJUX_RES
      external MATRIX_GUX_NINJUX_RES
      double precision MATRIX_UG_NINJU_RES
      external MATRIX_UG_NINJU_RES
      double precision MATRIX_UXG_NINJUX_RES
      external MATRIX_UXG_NINJUX_RES
#endif
#endif


c WREG is now defined in osres.h
#include "osres.h"

#endif

c############### end dsubtraction.h ####################################
