c############### init_couplings.f ######################################
c last modified by MK, 22.06.2015
c weakino pair production

c############### init_couplings subroutine #############################
c setup SM and MSSM parameters in PhysPars.h
c requires encoded process ID in lprup(1)
      subroutine init_couplings
        implicit none
#include "PhysPars.h"
#include "LesHouches.h"
#include "pwhg_st.h"
#include "pwhg_math.h"
#include "nlegborn.h"
#include "osres.h"
        double precision powheginput
        external powheginput
        ! final state IDs
        integer id_fin1, id_fin2
        logical verbose
        parameter (verbose=.true.)
        
        ! extract particle types from process id
        call decode_pair(lprup(1),id_fin1,id_fin2)
#ifdef DEBUG
        print*, "extracted from process ID: ", lprup(1)
        print*, "the final state pair: ", id_fin1, ", ", id_fin2
#endif

        ! set PhysPars parameters
        call set_phys_params
        ! set madgraph parameters.
        ! Note: calling this subroutine for the first time reads from 
        ! the provided slha file.
        call set_madgraph_params
        
        ! override SLHA parameters by explicit parameters from powheg.input
        ! this can be easily extended
        print*, "checking for parameters in powheg.input"
        print*, "(if given, they override SLHA SM parameters)"
        if(powheginput("#inv_alpha").gt.0) then
          print*,"Note: using inv_alpha from powheg.input"
          par_alpha = 1d0/powheginput("inv_alpha")
          alpha = par_alpha
          call set_madgraph_params ! overwrite the current madgraph alpha with the new one
        endif
        if(powheginput("#Zmass").gt.0) then
          print*,"Note: using Zmass from powheg.input"
          par_Zmass = powheginput("Zmass")
          zmass = par_Zmass
          call set_madgraph_params
        endif
        if(powheginput("#Wmass").gt.0) then
          print*,"Note: using Wmass from powheg.input"
          par_Wmass = powheginput("Wmass")
          wmass = par_Wmass
          call set_madgraph_params
        endif
        if(powheginput("#Zwidth").gt.0) then
          print*,"Note: using Zwidth from powheg.input"
          par_Zwidth = powheginput("Zwidth")
          zwidth = par_Zwidth ! madgraph
          call set_madgraph_params
        endif
        if(powheginput("#Wwidth").gt.0) then
          print*,"Note: using Wwidth from powheg.input"
          par_Wwidth = powheginput("Wwidth")
          wwidth = par_Wwidth ! madgraph
          call set_madgraph_params
        endif

        ! set remaining FormCalc parameters
        call set_formcalc_params

        ! hack to overwrite the madgraph alpha 
        ! (par_alpha is the one we want to use)
        if(alfa .ne. alpha) then
          alpha = par_alpha
          call set_madgraph_params
        endif

        ! set general MadGraph couplings
        call set_general_couplings(id_fin1,id_fin2)
        ! set the parameters for treating the on-shell resonances
        call set_osres_params


        ! Print out all parameters if you want
        if(verbose) then
          print*
          print*, "===================================================="
          print*, "             set POWHEG parameters                  "
          print*, "===================================================="
          print*
          print*, "Z mass          = ",par_Zmass
          print*, "W mass          = ",par_Wmass
          print*, "Z width         = ",par_Zwidth
          print*, "W width         = ",par_Wwidth
          print*, "1/alpha         = ",1d0/par_alpha
          print*
          print*, "final state masses = ", par_Fin1mass
          print*, "                     ", par_Fin2mass
          print*
          print*, "OS regulator WREG = ", wreg
          print*        
#ifdef DEBUGQ
          print*
          print*, "===================================================="
          print*, "             set FormCalc parameters                "
          print*, "===================================================="
          print*
          call printFormcalcModelParameters
          print*
          print*
          print*, "===================================================="
          print*, "             set MadGraph parameters                "
          print*, "===================================================="
          print*
          call printMadgraphModelParameters
          print*
#endif
        endif
        
        
        ! check if the MadGraph parameters are the same as the ones
        ! used by FormCalc
        if( (MZ                 .ne. zmass        ) .or.
     &      (MW                 .ne. wmass        ) .or.
     &      (WZ                 .ne. zwidth       ) .or.
     &      (WW                 .ne. wwidth       ) .or.
     &      (alfa               .ne. alpha        ) .or.
     &      (EL                 .ne. dreal(gal(1))) .or.
     &      (dabs(MNeu(1))      .ne. dabs(mn1)    ) .or.
     &      (dabs(MNeu(2))      .ne. dabs(mn2)    ) .or.
     &      (dabs(MNeu(3))      .ne. dabs(mn3)    ) .or.
     &      (dabs(MNeu(4))      .ne. dabs(mn4)    ) .or.
     &      (MCha(1)            .ne. mx1          ) .or.
     &      (MCha(2)            .ne. mx2          ) .or.
     &      (MSf(1,3,1)         .ne. mul          ) .or.
     &      (MSf(2,3,1)         .ne. mur          ) .or.
     &      (MSf(1,4,1)         .ne. mdl          ) .or.
     &      (MSf(2,4,1)         .ne. mdr          ) .or.
     &      (MSf(1,3,2)         .ne. mcl          ) .or.
     &      (MSf(2,3,2)         .ne. mcr          ) .or.
     &      (MSf(1,4,2)         .ne. msl          ) .or.
     &      (MSf(2,4,2)         .ne. msr          ) .or.
     &      (MSf(1,3,3)         .ne. mtl          ) .or.
     &      (MSf(2,3,3)         .ne. mtr          ) .or.
     &      (MSf(1,4,3)         .ne. mbl          ) .or.
     &      (MSf(2,4,3)         .ne. mbr          ) .or.
     &      (MGl                .ne. mgo          ) .or.
     &      (WGl                .ne. wgo          ) .or.
     &      (WSf(1,4,1)         .ne. wdl          ) .or.
     &      (WSf(2,4,1)         .ne. wdr          ) .or.
     &      (WSf(1,3,1)         .ne. wul          ) .or.
     &      (WSf(2,3,1)         .ne. wur          ) .or.
     &      (WSf(1,4,2)         .ne. wsl          ) .or.
     &      (WSf(2,4,2)         .ne. wsr          ) .or.
     &      (WSf(1,3,2)         .ne. wcl          ) .or.
     &      (WSf(2,3,2)         .ne. wcr          ) .or.
     &      (WSf(1,4,3)         .ne. wb1          ) .or.
     &      (WSf(2,4,3)         .ne. wb2          ) .or.
     &      (WSf(1,3,3)         .ne. wt1          ) .or.
     &      (WSf(2,3,3)         .ne. wt2          ) .or.
     &      (cdabs(ZNeu(1,1))   .ne. dabs(bwmix(1,1)) ) .or. ! check only if the absolute
     &      (cdabs(ZNeu(1,2))   .ne. dabs(bwmix(1,2)) ) .or. ! value is equal, since the
     &      (cdabs(ZNeu(1,3))   .ne. dabs(bwmix(1,3)) ) .or. ! FormCalc-mixing matrices
     &      (cdabs(ZNeu(1,4))   .ne. dabs(bwmix(1,4)) ) .or. ! can be complex
     &      (cdabs(ZNeu(2,1))   .ne. dabs(bwmix(2,1)) ) .or.
     &      (cdabs(ZNeu(2,2))   .ne. dabs(bwmix(2,2)) ) .or.
     &      (cdabs(ZNeu(2,3))   .ne. dabs(bwmix(2,3)) ) .or.
     &      (cdabs(ZNeu(2,4))   .ne. dabs(bwmix(2,4)) ) .or.
     &      (cdabs(ZNeu(3,1))   .ne. dabs(bwmix(3,1)) ) .or.
     &      (cdabs(ZNeu(3,2))   .ne. dabs(bwmix(3,2)) ) .or.
     &      (cdabs(ZNeu(3,3))   .ne. dabs(bwmix(3,3)) ) .or.
     &      (cdabs(ZNeu(3,4))   .ne. dabs(bwmix(3,4)) ) .or.
     &      (cdabs(ZNeu(4,1))   .ne. dabs(bwmix(4,1)) ) .or.
     &      (cdabs(ZNeu(4,2))   .ne. dabs(bwmix(4,2)) ) .or.
     &      (cdabs(ZNeu(4,3))   .ne. dabs(bwmix(4,3)) ) .or.
     &      (cdabs(ZNeu(4,4))   .ne. dabs(bwmix(4,4)) ) .or.
     &      (UCha(1,1)          .ne. uumix(1,1)   ) .or.
     &      (UCha(1,2)          .ne. uumix(1,2)   ) .or.
     &      (UCha(2,1)          .ne. uumix(2,1)   ) .or.
     &      (UCha(2,2)          .ne. uumix(2,2)   ) .or.
     &      (VCha(1,1)          .ne. vvmix(1,1)   ) .or.
     &      (VCha(1,2)          .ne. vvmix(1,2)   ) .or.
     &      (VCha(2,1)          .ne. vvmix(2,1)   ) .or.
     &      (VCha(2,2)          .ne. vvmix(2,2)   ) ) then
          print*, "Error: Madgraph uses parameters different to "
     &          //"PhysPars.h or FormCalc paramaters"
          stop
        endif

        ! check if all paramaters are properly set
        if ( (par_alpha    .eq. invalid) .or.
     &       (par_Zmass    .eq. invalid) .or.
     &       (par_Wmass    .eq. invalid) .or.
     &       (par_Zwidth   .eq. invalid) .or.
     &       (par_Wwidth   .eq. invalid) .or.
     &       (par_Zmass2   .eq. invalid) .or.
     &       (par_Wmass2   .eq. invalid) .or.
     &       (par_Fin1mass .eq. invalid) .or.
     &       (par_Fin2mass .eq. invalid) ) then
          print*,"WARNING: One or more parameters are not set:"
          print*,"par_alpha: ",par_alpha
          print*,"par_Zmass: ",par_Zmass
          print*,"par_Wmass: ",par_Wmass
          print*,"par_Zwidth: ",par_Zwidth
          print*,"par_Wwidth: ",par_Wwidth
          print*,"par_Zmass2: ",par_Zmass2
          print*,"par_Wmass2: ",par_Wmass2
          print*,"par_Fin1mass: ",par_Fin1mass
          print*,"par_Fin2mass: ",par_Fin2mass
          print*,"If possible set the invalid parameters in "//
     &           "powheg.input"
          stop
        endif

        ! init LoopTools
        print*, "===================================================="
        print*, "             init LoopTools                         "
        call ltini
        !call setversionkey(0)
        print*
      end
c############### end subroutine init_couplings #########################

c############### subroutine init_madgraph_params #######################
c sets the Madgraph parameters
      subroutine set_madgraph_params
        implicit none
#include "PhysPars.h"
        character*100 slhafilename
        logical called, readslha
        data called/.false./
        save called
        ! reads from slha file only when called this subroutine
        ! for the first time. When this subroutine is called
        ! more than once, the SM Input parameters can be overwritten
        ! with parameters from powheg.input.
        if(called) then
          print*, "Note: Overwrite madgraph paramaters."
          readslha = .false.
        else
          call powheginputstring("SLHA",slhafilename)
          called   = .true.
          readslha = .true.
        endif
        call setpara(slhafilename,readslha)
      end
c############### end subroutine init_madgraph_params ###################

c############### subroutine set_ebe_couplings ##########################
c set the madgraph and formcalc couplings event by event
      subroutine set_ebe_couplings
        implicit none
#include "PhysPars.h"
#include "pwhg_st.h"
#include "pwhg_math.h"
        ! MadGraph
        ! QCD coupling constant
        alphaS = st_alpha
        G      = dsqrt(alphaS*4d0*pi)
        GG(1)  = -G
        GG(2)  = -G
        ! e.m. coupling
        alpha  = par_alpha
        Gal(1) = dcmplx(dsqrt(alpha*4d0*pi))
        Gal(2) = dcmplx(dsqrt(alpha*4d0*pi))
        ! FormCalc
        ! QCD coupling
        AlfaS  = alphaS
        AlfaS2 = AlfaS**2
        GS = dsqrt(AlfaS*4d0*pi)
        ! e.m. coupling
        Alfa = par_alpha
        Alfa2   = Alfa**2
        EL = dsqrt(Alfa*4d0*pi)
        ! initialize the counter term wich cures the susy breaking
        ! introduced by DR
        call init_susy_restore
#ifdef DEBUGQ
        print*,"AlfaS: ", AlfaS
        print*,"Alfa:  ", Alfa
        stop
#endif
        return
      end
c############### end subroutine set_ebe_couplings ######################

c############### set_phys_params subroutine ############################
c extracts mass and mixing matrices from slhadata
c assumes slhadata to be filled already
c sets the PhysPars.h parameters
      subroutine set_phys_params
        implicit none
#include "PhysPars.h"
#include "pwhg_math.h"
#include "LesHouches.h"
        double precision powheginput
        external powheginput
        integer i,j,k   ! loops
        ! decay widths
        integer particle, channel
        ! for reading parameters from SHLA file
        integer error
        character*100 slhafilename
        ! final state IDs
        integer id_fin1, id_fin2
        ! weinberg angle
        double precision sinW
        
        par_alpha    = invalid
        par_Zmass    = invalid
        par_Wmass    = invalid
        par_Zmass2   = invalid
        par_Wmass2   = invalid
        par_Zwidth   = invalid
        par_Wwidth   = invalid
        par_Fin1mass = invalid
        par_Fin2mass = invalid

        ! get SLHA filename
        call powheginputstring("SLHA",slhafilename)
        ! read in parameters from SLHA file
        call SLHARead(error, slhadata, slhafilename, 1)
        if (error.ne.0)  stop "error reading parameters from SLHA file"

        ! read all particle widths from SLHA-file
        nchannels = SLHADecayTable(slhadata, 0, width, partid, 7, 1024)

#ifdef DEBUGQ
        do channel = 1, nchannels
          nparticles = partid(0,channel)
          if(nparticles .eq. 1) then
            print*,"==> total width(",partid(1,channel),") = "
     &            ,width(channel)
          else
            print*, "--> partial width(",partid(1,channel),"to"
            do particle = 2, nparticles
              print*,partid(particle,channel)
            enddo
            print*, "            ) = ", width(channel)
          endif
          print*
        enddo
        stop
#endif

        ! SM parameters from SLHA-file
        print*, "note: ignoring alpha_s value given in SLHA file"
        par_Zmass   = dreal(SMInputs_MZ)
        par_Wmass   = dreal(Mass_MW)
        par_gf      = dreal(SMInputs_GF)
#define SLHA_ALPHA
#ifdef SLHA_ALPHA
        par_alpha   = 1d0/dreal(SMInputs_invAlfaMZ)
#else
        sinW        = dsqrt(1d0-par_Wmass**2/par_Zmass**2)
        par_alpha   = dsqrt(2d0)*par_Wmass**2*sinW**2*par_gf/pi
#endif
        !print*,1/par_alpha
        !stop
        par_Zmass2  = par_Zmass**2
        par_Wmass2  = par_Wmass**2

        ! particle widths from SLHA-file
        do channel = 1, nchannels
          nparticles = partid(0,channel)
          if(nparticles .eq. 1) then
            if(partid(1,channel).eq.23) then
              par_Zwidth = width(channel)
            else if(partid(1,channel).eq.24) then
              par_Wwidth = width(channel)
            endif
          endif
        enddo

        if( par_Wwidth .eq. invalid) then
          print*, "WARNING: W-width not found in SLHA file."
          print*, "         Use zero width."
          par_Wwidth = 0D0
        endif

        if( par_Zwidth .eq. invalid) then
          print*, "WARNING: Z-width not found in SLHA file."
          print*, "         Use zero width."
          par_Zwidth = 0D0
        endif
        
#ifdef DEBUGQ
        print*, "par_Zwidth", par_Zwidth
        print*, "par_Wwidth", par_Wwidth
#endif
        
        ! extract particle types from process id
        call decode_pair(lprup(1),id_fin1,id_fin2)
#ifdef DEBUG
        print*, "extracted from process ID: ", lprup(1)
        print*, "the final state pair: ", id_fin1, ", ", id_fin2
#endif
        
        select case(abs(id_fin1))
          case (1000022) ! neutralino1
            par_Fin1mass = dreal(Mass_MNeu(1))
          case (1000023) ! neutralino2
            par_Fin1mass = dreal(Mass_MNeu(2))
          case (1000025) ! neutralino3
            par_Fin1mass = dreal(Mass_MNeu(3))
          case (1000035) ! neutralino4
            par_Fin1mass = dreal(Mass_MNeu(4))
          case (1000024) ! chargino1+
            par_Fin1mass = dreal(Mass_MCha(1))
          case (1000037) ! chargino2+
            par_Fin1mass = dreal(Mass_MCha(2))
          case (2000011) ! slepton+
            par_Fin1mass = dreal(Mass_MSf(2,2,1))  
          case default
            write(*,*) "encountered unknown particle ID ", id_fin1
            stop
        end select
        select case(abs(id_fin2))
          case (1000022) ! neutralino1
            par_Fin2mass = dreal(Mass_MNeu(1))
          case (1000023) ! neutralino2
            par_Fin2mass = dreal(Mass_MNeu(2))
          case (1000025) ! neutralino3
            par_Fin2mass = dreal(Mass_MNeu(3))
          case (1000035) ! neutralino4
            par_Fin2mass = dreal(Mass_MNeu(4))
          case (1000024) ! chargino1-
            par_Fin2mass = dreal(Mass_MCha(1))
          case (1000037) ! chargino2-
            par_Fin2mass = dreal(Mass_MCha(2))
          case (2000011) ! slepton-
            par_Fin2mass = dreal(Mass_MSf(2,2,1))
          case default
            write(*,*) "encountered unknown particle ID ", id_fin2
            stop
        end select
        
        if( (dabs(par_Fin1mass) + dabs(par_Fin2mass)) .lt. par_Zmass )
     &  then
          print*, "Error: Final state masses low, W/Z decays might "//
     &            "be more suitable"
          stop
        endif
        
        ! Sfermion masses
        do i=1,2     ! left, right
          do j=1,4   ! type (nu,e,u,d)
            do k=1,3 ! generation
              par_MSf(i,j,k) = dreal(Mass_MSf(i,j,k))
            enddo
          enddo
        enddo
        
        ! Neutralino masses
        do i=1,4
          par_MNeu(i)  = dreal(Mass_MNeu(i))
        enddo
        
        ! Chargino masses
        do i=1,2
          par_MCha(i)  = dreal(Mass_MCha(i))
        enddo
        
        ! gluino mass
        par_MGl = dreal(Mass_MGl)

c small mass to regularize LoopTools integrals
#define EPSMASS 1D-30
        
        print*,"Note (PhysPars): neglecting u,d,c,s-quark masses"
        !epsMass = powheginput("epsMass")
        par_MU = EPSMASS ! Note: setting this exactly to zero can cause 
        par_MC = EPSMASS ! problems in calculating propagators, since the
        par_MD = EPSMASS ! formcalc code already neglects all light 
        par_MS = EPSMASS ! masses analytically but leaves some masses in propagators

c b quark mass
#ifdef NEGLECTBMASS
        print*,"Note (PhysPars): neglecting b-quark mass"
        par_MB = EPSMASS
#else
        ! read mass paramater from slha file
        par_MB  = dreal(Mass_Mf(4,3)) 
#endif

        ! top mass
        !par_MT  = dreal(Mass_Mf(3,3))
        par_MT = dreal(SMInputs_Mt)

c particle widths from SLHA-file
c NOTE: FormCalc does not use squark widths in the current software 
c version.
        do channel = 1, nchannels
          nparticles = partid(0,channel)
          if(nparticles .eq. 1) then
            select case(partid(1,channel))
              case (1000021) ! gluino
                par_Glwidth = width(channel)
#ifndef NEGLECT_SWIDTH
              print*, "WARNING: FormCalc does not use squark-widths."  
              case (1000001) ! sdown_L
                par_Sfwidth(1,4,1) = width(channel)
              case (2000001) ! sdown_R
                par_Sfwidth(2,4,1) = width(channel)
              case (1000002) ! sup_L
                par_Sfwidth(1,3,1) = width(channel)
              case (2000002) ! sup_R
                par_Sfwidth(2,3,1) = width(channel)
              case (1000003) ! sstrange_L
                par_Sfwidth(1,4,2) = width(channel)
              case (2000003) ! sstrange_R
                par_Sfwidth(2,4,2) = width(channel)
              case (1000004) ! scharm_L
                par_Sfwidth(1,3,2) = width(channel)
              case (2000004) ! scharm_R
                par_Sfwidth(2,3,2) = width(channel)
              case (1000005) ! sbottom_1
                par_Sfwidth(1,4,3) = width(channel)
              case (2000005) ! sbottom_2
                par_Sfwidth(2,4,3) = width(channel)
              case (1000006) ! stop_1
                par_Sfwidth(1,3,3) = width(channel)
              case (2000006) ! stop_2
                par_Sfwidth(2,3,3) = width(channel)
#endif
            end select
          endif
        enddo
#ifdef NEGLECT_SWIDTH
        do i=1,2 ! sfermion index s
          do j=3,4 ! type (u or d)
            do k=1,3 ! generation
              WSf(i,j,k) = 0d0
            enddo
          enddo
        enddo
#endif
        
#ifdef DEBUGQ
        print*, "sd_L width =", par_Sfwidth(1,4,1)
        print*, "sd_R width =", par_Sfwidth(2,4,1)
        print*, "go width   =", par_Glwidth
        stop
#endif
      end
c############### end subroutine set_phys_params ########################

c############### set_params subroutine #################################
c extracts mass and mixing matrices from slhadata
c assumes slhadata to be filled already
c sets the FormCalc parameters
      subroutine set_formcalc_params
        implicit none
#include "PhysPars.h"
#include "pwhg_math.h"
        double precision powheginput
        external powheginput
        integer i,j,k,l ! loops
        real ident(4,4) ! identity matrix
        real eps        ! auxilliary variable for checking unitarity of ZNeu
        real cP1        ! (checkLimitPlus 1) upper limit for checking unitarity
        real cM1        ! (checkLimitMinus 1) lower limit -"-
        real cP0        ! (checkLimitPlus 0) upper limit for checking unitarity
        real cM0        ! (checkLimitMinus 0) lower limit -"-
        ! Avoid multiple calls to this subroutine.
        logical called
        data called/.false./
        save called
        if(called) then
          return
        else
          called=.true.
        endif

c CKM-Matrix
c if you change this, you have to change something in
c MadGraphs couplings.f as well.
#define NO_CKM_MIXING
#ifdef NO_CKM_MIXING
        print*, "Note (FormCalc): no CKM mixing"
        do i=1,3
          do j=1,3
            if(i.eq.j)  CKM(i,j) = 1D0
            if(i.ne.j)  CKM(i,j) = 0D0
          enddo
        enddo
#else
        CKM(1,1) = 0.97427d0
        CKM(1,2) = 0.22534d0
        CKM(1,3) = 0.00351d0
        CKM(2,1) = 0.22520d0
        CKM(2,2) = 0.97344d0
        CKM(2,3) = 0.04120d0
        CKM(3,1) = 0.00867d0
        CKM(3,2) = 0.04040d0
        CKM(3,3) = 0.99915d0
#endif

        ! Squark masses and widths
        ! Note: The squark widths are not needed in the current 
        ! implementation of FormCalc. They will be initialized here
        ! anyway.
        do i=1,2 ! sfermion index s
          do j=1,4 ! type (nu,e,u,d)
            do k=1,3 ! generation
              MSf(i,j,k)  = par_MSf(i,j,k)
              MSf2(i,j,k) = MSf(i,j,k)**2
              WSf(i,j,k)  = par_Sfwidth(i,j,k)
            enddo
          enddo
        enddo
       
        ! Neutralino masses
        do i=1,4 !neutralino index s
          MNeu(i)  = par_MNeu(i)
          MNeu2(i) = MNeu(i)**2
        enddo
        
        ! Chargino masses
        do i=1,2 !chargino index c
          MCha(i)  = par_MCha(i)
          MCha2(i) = MCha(i)**2
        enddo
        
        ! Z/W mass
        MZ  = par_Zmass
        MZ2 = MZ**2
        MW  = par_Wmass
        MW2 = MW**2

        ! Quark masses
        MU   = par_MU
        MU2  = MU**2
        MC   = par_MC
        MC2  = MC**2
        MD   = par_MD
        MD2  = MD**2
        MS   = par_MS
        MS2  = MS**2
        MB   = par_MB
        MB2  = MB**2
        MT   = par_MT
        MT2  = MT**2
        
        Mf(1,1) = 0D0
        Mf(2,1) = 0D0
        Mf(3,1) = MU
        Mf(4,1) = MD
        Mf(1,2) = 0D0
        Mf(2,2) = 0D0
        Mf(3,2) = MC
        Mf(4,2) = MS
        Mf(1,3) = 0D0
        Mf(2,3) = 0D0
        Mf(3,3) = MT
        Mf(4,3) = MB

        Mf2(1,1) = 0D0
        Mf2(2,1) = 0D0
        Mf2(3,1) = MU2
        Mf2(4,1) = MD2
        Mf2(1,2) = 0D0
        Mf2(2,2) = 0D0
        Mf2(3,2) = MC2
        Mf2(4,2) = MS2
        Mf2(1,3) = 0D0
        Mf2(2,3) = 0D0
        Mf2(3,3) = MT2
        Mf2(4,3) = MB2
        
        ! Gluino mass and width
        MGl  = par_Mgl
        MGl2 = MGl**2
        WGl  = par_Glwidth
        
        ! Z/W Widths (added by hand to formcalc)
        WZ  = par_Zwidth
        WW  = par_Wwidth

        ! Neutralino mixing matrix
        ! SLHA name: NMix_ZNeu(n1,n2)
        do i=1,4 ! neutralino index n1
          do j=1,4 ! neutralino index n2
            if(MNeu(i).lt.0D0) then ! check if we need a complex mixing matrix
              ZNeu(i,j) = ii*dcmplx(NMix_ZNeu(i,j))
            else
              ZNeu(i,j) = dcmplx(NMix_ZNeu(i,j))
            endif  
          enddo
        enddo
        
        ! Chargino mixing matrix
        do i=1,2 ! chargino index c1
          do j=1,2 ! chargino index c2
            VCha(i,j) =  dcmplx(VMix_VCha(i,j))
            UCha(i,j) =  dcmplx(UMix_UCha(i,j))
          enddo
        enddo

        ! Squark mixing matrices (u, d, s, c - sfermions)
        do i=1,2 ! sfermion index s1 (left, right)
          do j=1,2 ! sfermion index s2
            do k=2,4 ! type index t (u- or d-type)
              do l=1,2 ! generation index g (note: generation 3 is in line 448)
                if(i.eq.j)  USf(i,j,k,l) = (1D0,0D0)
                if(i.ne.j)  USf(i,j,k,l) = (0D0,0D0)
              enddo
            enddo
          enddo
        enddo
        do i=1,2 ! sfermion index s1
          do j=1,2 ! sfermion index s2
            ! Mixing between t-sfermions (mixing)
            USf(i,j,3,3) = dcmplx(StopMix_USf(i,j))
            ! Mixing between b-sfermions (mixing)
            USf(i,j,4,3) = dcmplx(SbotMix_USf(i,j))
          enddo
        enddo
        
        ! mixing angles
        CW  = MW/MZ
        CW2 = CW**2
        SW2 = 1 - CW2
        SW  = dsqrt(SW2)
        
        ! sin(beta), cos(beta), tan(beta)
        ! Note: the matrix elements do not depend on TB
        !TB = HMix_TB
        TB = MinPar_TB
        CB = dabs(1.D0/dsqrt(1+TB**2))
        SB = dabs(TB*CB)
        SB2 = SB**2
        CB2 = CB**2
        TB2 = TB**2
        S2B = 2*TB*CB2
        C2B = (CB - SB)*(CB + SB)
        
        ! Gaugino mass parameter
        !Mino3 = dcmplx(Mgl)
        Mino3 = dcmplx(MinPar_M12)
        
        ! trilinear couplings
        ! Note: the matrix elements do not depend on Af
        do i=2,4 ! type 
          do j=1,3 ! generation
            do k=1,3 ! generation
              Af(i,j,k) = dcmplx(Af_Af(j,k,i))
            enddo
          enddo
        enddo
        
        ! mass term mixing the 2 Higgs doublets
        ! Note: the matrix elements do not depend on MUE
        MUE  = dcmplx(HMix_MUE)

        ! couplings
        ! set Formcalc couplings equal to the ones provided by POWHEG
        Alfa    = par_alpha
        Alfa2   = Alfa**2
        AlfaS   = SMInputs_AlfasMZ ! dummy value, gets overwritten by set_ebe_couplings
        AlfaS2  = AlfaS**2
        AlfaMZ  = 1D0/1.27934000D2
        AlfaSMZ = SMInputs_AlfasMZ
        EL = dsqrt(4*pi*Alfa)
        GS = dsqrt(4*pi*AlfaS)
      end
c############### end subroutine set_params #############################

c############### init_susy_restore subroutine ##########################
c initialization of susy restoring counter terms and yukawa couplings
c quark-squark-gluino (not needed here):
c     gsy -> gsy + dZgs1y
c                  dZgs1y = gs*Alfas/(3*Pi)
c and SUSY-restoring term for Yukawa coupling
c quark-squark-neutralino:
c     gy -> gy + dZe1y
c                dZe1y = - e*Alfas/(6*Pi)
c the counter terms at amplitude squared niveau is implemented in
c ddbar_SquaredME.F and uubar_SquaredME.F.
      subroutine init_susy_restore
        implicit none
#include "PhysPars.h"
#include "pwhg_st.h"
#include "pwhg_math.h"
#include "susy_restore.h"
        ! the counter term for quark-squark-gluino is not needed here
        dZgs1y = 0 ! gs*Alfas/(3D0*Pi)
        ! counter term for quark-squark-neutralino yukawa coupling
        dZe1y  = -el*Alfas/(6D0*Pi)
        
#ifdef CHECKUV
        dZgs1y = 0
        dZe1y  = 0
#endif
        ! not needed in the current implementation, but can be 
        ! useful for debugging purpose
        ! the yukawa couplings
        gsy  = gs+dZgs1y
        ely  = el+dZe1y
        ! the squared yukawa-couplings. Double check the code for
        ! consistency! In the matrix elements must only occur terms
        ! prop to Alfa**2 + Alfa**2*AlfaS (for weakino pair production)
        ! or      AlfaS**2 + AlfaS**3 (for disquark)
        gsy2 = gs**2+2D0*gs*dZgs1y
        ely2 = el**2+2D0*el*dZe1y
        
#ifdef DEBUGQ
        print*,"dZgs1y",dZgs1y
        print*,"dZe1y",dZe1y
        print*,"GS",GS
        print*,"GSY",GSY
        print*,"EL",EL
        print*,"ELY",ELY
        print*,"st_alpha", st_alpha
        print*,"AlfaS", AlfaS
        print*,"Alfa",Alfa
        print*,"EL2/4Pi",EL**2/(4*Pi)
        !stop
#endif
        return
      end
c############### end init_susy_restore subroutine ######################

c############### subroutine set_osres_params ###########################
c set the paramaters for treating the on-shell resonances
      subroutine set_osres_params
        implicit none
#include "PhysPars.h"
#include "osres.h"
        integer i
        double precision m_avrg, ewi
        double precision powheginput
        external powheginput

        osresM(1) = mdl
        osresM(2) = mdl
        osresM(3) = mul
        osresM(4) = mul
#ifdef NINJ
        osresM(5) = mdr
        osresM(6) = mdr
        osresM(7) = mur
        osresM(8) = mur
#endif
        
#if defined(DSUB_I) || defined(DSUB_II) || defined(DSUB_II_TEST)
        m_avrg = 0D0
        do i=1,nosres
          m_avrg = m_avrg + osresM(i)
        enddo
        m_avrg = m_avrg/nosres

        ! set the width of the immediate particle (regulator wreg)
        ewi  = 1D-2 ! (default)
        if(powheginput("#ewi").gt.0) then
          ewi = powheginput("ewi")
        endif
        wreg = ewi*m_avrg ! Prospino like definition
#else
        wreg = 0D0
#endif

        
#ifdef DEBUGQ
        print*,"wreg : ", wreg
        do i=1,nosres
          print*,"osresID: ", osresID(i)
          print*,"osresM : ", osresM(i)
        enddo 
        stop
#endif
      end
c############### end subroutine set_osres_params #######################
