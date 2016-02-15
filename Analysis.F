c############### subroutine init_hist ##################################
c book all histograms, will be filled later
      subroutine init_hist
        implicit none
#include "LesHouches.h"
#include "pwhg_math.h"
        ! tell this analysis file which program is running it
        character *6 WHCPRG
        common/cWHCPRG/WHCPRG

        call inihists
        
        ! no cuts
        call bookupeqbins("sigma_total_nocuts",1D0,-0.5D0,0.5D0)
        call bookupeqbins("inv_mass_pair_nocuts",10D0,200D0,1000D0)
        call bookupeqbins("pt_finA_nocuts",6D0,0D0,500D0)
        call bookupeqbins("pt_finB_nocuts",6D0,0D0,500D0)
        call bookupeqbins("eta_finA_nocuts",0.2D0,-4D0,4D0)
        call bookupeqbins("number_of_visible_jets_nocuts",
     &                    1D0,-0.5D0,2.5D0)
        call bookupeqbins("pt_jet_nocuts",2D0,0D0,100D0)
        call bookupeqbins("eta_jet_nocuts",0.2D0,-4D0,4D0)
        
        ! jet cuts
        call bookupeqbins("sigma_total_jcuts",1D0,-0.5D0,0.5D0)
        call bookupeqbins("inv_mass_pair_jcuts",10D0,200D0,1000D0)
        call bookupeqbins("pt_finA_jcuts",6D0,0D0,500D0)
        call bookupeqbins("pt_finB_jcuts",6D0,0D0,500D0)
        call bookupeqbins("eta_finA_jcuts",0.2D0,-4D0,4D0)
        call bookupeqbins("number_of_visible_jets_jcuts",
     &                    1D0,-0.5D0,2.5D0)
        call bookupeqbins("pt_jet_jcuts",50D0,0D0,1000D0)
        call bookupeqbins("eta_jet_jcuts",0.2D0,-4D0,4D0)
        
        if(WHCPRG.eq.'PYTHIA') then
          ! lepton cuts
          call bookupeqbins('sigma_total_lcuts',1d0,-0.5d0,0.5d0)
          call bookupeqbins('azimuth_diff_lepton_pair_lcuts',
     &                       pi/20d0,-pi,pi)
          call bookupeqbins('pt_lepton_A_lcuts',6d0,0d0,600d0)
          call bookupeqbins('eta_lepton_A_lcuts',0.2d0,-4d0,4d0)
          call bookupeqbins('pt_lepton_B_lcuts',6d0,0d0,600d0)
          call bookupeqbins('eta_lepton_B_lcuts',0.2d0,-4d0,4d0)
        endif  
        
      end
c############### end subroutine init_hist ##############################

c############### subroutine analysis ###################################
c extract all data required for the histograms
c calculate quantities that shall be plotted
c fill histograms
      subroutine analysis(dsig)
        implicit none
#include "pwhg_math.h"
#include "LesHouches.h"
        double precision dsig ! total cross section
        ! tells if subroutine is called for the first time
        logical ini
        data ini/.true./
        save ini
        ! tell this analysis file which program is running it
        character *6 WHCPRG
        common/cWHCPRG/WHCPRG
        data WHCPRG/'NLO'/ ! default, should be set when analysis is called
        ! external functions
        double precision powheginput
        external powheginput
        ! kinematical quantities
        double precision p_finA(0:3), p_finB(0:3), p_finAfinB(0:3)
        double precision pt_finA, pt_finB
        double precision eta_finA
        double precision y_finA
        double precision m_finAfinB
        double precision pt_finAfinB, y_finAfinB
        double precision p_decA(0:3), p_decB(0:3)
        double precision pt_decA, pt_decB
        double precision eta_decA, eta_decB
        double precision y_decA
        double precision phi_decA, phi_decB, phi_decAdecB
        ! jet parameters
        integer j, mjets, num_vis_jets, maxjet
        parameter (maxjet=2048)
        double precision kt(maxjet), eta(maxjet), rap(maxjet)
        double precision phi(maxjet), pj(4,maxjet), ptrel(maxjet)
        ! ptJmin is a technical cut
        double precision ptJmin, etaJmax, R
        common /jetcuts/ptJmin, etaJmax, R
        ! pt_min for lepton cuts (physical cut)
        double precision ptLmin, etaLmax
        ! IDs
        integer id_finA, id_finB
        ! IDs of leptonic decay products
        integer id_decA, id_decB
        save id_decA, id_decB
        ! indices
        integer mu, ij
        ! cuts
        logical JcutsPassed, Lcutspassed
        ! decay
        integer number_lepA,number_lepB
        logical decayed, isvaliddecayspectrum
        save decayed
        decayed=.true.
        
        ! at first call
        if (ini) then            
          ! check program that runs analysis
          print*,
          print*, "===================================================="
          if(WHCPRG.eq."NLO") then
            print*, "fixed-order NLO analysis"
          elseif(WHCPRG.eq."LHE") then
            print*, "LHE analysis"
            call init_couplings    
          elseif(WHCPRG.eq."HERWIG") then
            print*, "HERWIG ANALYSIS"
            print*, "not implemented analysis"
            print*, "no plots will be present at the end of the run"
            call init_couplings
          elseif(WHCPRG.eq."PYTHIA") then
            print*, "PYTHIA ANALYSIS"
            call init_couplings
          endif
          print*, "===================================================="
 
          ! determine the decay products from powheg.input
          id_decA = powheginput("dec1")
          id_decB = powheginput("dec2")
          
          ! extract final state types from process id
          call decode_pair(lprup(1),id_finA,id_finB)

          ! check final state
          if (.not.((id_finA.eq.1000022 .or. id_finA.eq. 1000023  .or.
     &               id_finA.eq.1000025 .or. id_finA.eq. 1000035  .or.
     &               id_finA.eq.1000024 .or. id_finA.eq.-1000024  .or.
     &               id_finA.eq.1000037 .or. id_finA.eq.-1000037) .and.
     &              (id_finB.eq.1000022 .or. id_finB.eq. 1000023  .or.
     &               id_finB.eq.1000025 .or. id_finB.eq. 1000035  .or.
     &               id_finB.eq.1000024 .or. id_finB.eq.-1000024  .or.
     &               id_finB.eq.1000037 .or. id_finB.eq.-1000037))) then
            print*, "analysis: final state not supported"
            stop
          endif
          
          ! check decay products
          ! only leptons are allowed
          if (.not.((iabs(id_decA).eq.11)   .or.
     &              (iabs(id_decA).eq.13)   .or. 
     &              (iabs(id_decA).eq.15))  .or. 
     &        .not.((iabs(id_decB).eq.11)   .or.
     &              (iabs(id_decB).eq.13)   .or. 
     &              (iabs(id_decB).eq.15)) ) then
            print*,"unsupported decay channel."
            print*,"id_decA: ", id_decA
            print*,"id_decB: ", id_decB
            print*,"In the current version are only charged leptons "//
     &             "allowed."
            stop
          endif
          ! the two leptons have to be different in the current implementation
          !if (id_decA.eq.id_decB) then
          !  print*,"the leptons have to be different"
          !  print*,"change type or sign of one lepton in powheg.input"
          !  stop
          !endif  
          
          ! if everything is ok print status
          print*, "extracted from process ID: ", lprup(1)
          print*, "the final state pair: ", id_finA, ", ", id_finB
          print*, "the decay products:   ", id_decA, ", ", id_decB
          
          ! find out if allowed decay took place
          ! determine the number of the specified leptons in decay products
          call getnumber(id_decA,number_lepA)
          call getnumber(id_decB,number_lepB)
          ! allowed decay requires the specified leptons in final state
          if ((number_lepA.eq.0).and.(number_lepB.eq.0)) then
            decayed=.false.
            print*, "particles with id",id_decA,id_decB,"not found"
            print*, "-> no allowed decay took place"
          endif
          ! only if we already know that a decay took place
          if (decayed) then
            ! check mass spectrum if decays are possible
            if (.not.(isvaliddecayspectrum(id_finA,id_finB,
     &                                     id_decA,id_decB))) then
              decayed=.false.
              print*, "histograms for decay products are empty"
            endif
          endif
        endif
        
        ! jets
        ! initialize
        do j = 1,maxjet
          kt(j)    = 0D0
          eta(j)   = 0D0
          rap(j)   = 0D0
          phi(j)   = 0D0
          pj(:,j)  = 0D0
          ptrel(j) = 0D0
        enddo
        
        ! init cuts (change here if needed)
        ! jet cuts
        ptJmin  = 1D0
        etaJmax = 4.5D0
        R       = 0.4D0
        mjets   = 40
        ! lepton cuts
        ptLmin  = 20D0 ! physical cut
        etaLmax = 2.5D0
        
        if (ini) then
          ! write out informations
          print*, "===================================================="
          ! write out all cuts
          print*, "cuts in analysis routines"
          print*
          print*, "cuts (jet):"
          print*, "hardest jet   p_t >",ptJmin," GeV"
          print*, "hardest jet |eta| <",etaJmax
          print*, "radius parameter R = ",R
          print*
          print*, "cuts (leptons:)"
          print*, "hardest lepton   p_t >",ptLmin," GeV"
          print*, "hardest lepton |eta| <",etaLmax
          print*
          print*, "===================================================="
        endif
        
        ini=.false.
        
        ! from pico to femto
        dsig = dsig*1000D0

        ! find final state 4-momenta
        call get4momentum_fin_pair(id_finA,p_finA,id_finB,p_finB)
        
        ! find lepton 4-momenta of decay products
        if(decayed) then
          ! undecayed lepton with largest pt
          call get4momentum_largestpt(id_decA,p_decA)
          ! undecayed lepton with largest pt
          call get4momentum_largestpt(id_decB,p_decB)
        endif

        ! calculate final state kinematical quantities:
        ! final state A, B transverse momentum
        pt_finA = dsqrt(p_finA(1)**2+p_finA(2)**2)
        pt_finB = dsqrt(p_finB(1)**2+p_finB(2)**2)
        ! final state A pseudorapidity
        call get_pseudorap(p_finA,eta_finA)
        ! final state A rapidity
        call getrapidity(p_finA,y_finA)
        ! invariant mass of the final state system
        do mu=0,3
          p_finAfinB(mu) = p_finB(mu)+p_finA(mu)
        enddo
        call getinvmass(p_finAfinB,m_finAfinB)
        ! transverse momentum of the final state system
        pt_finAfinB = dsqrt((p_finB(1)+p_finA(1))**2
     &                     +(p_finB(2)+p_finA(2))**2)
        ! rapidity of the final state system
        call getrapidity(p_finAfinB,y_finAfinB) 
        
        ! build jets
        call buildjets(mjets,kt,eta,rap,phi,pj,ptrel)
        ! determine number of visible jets
        num_vis_jets = 0
        do ij = 1,mjets
          if ((kt(ij).gt.ptJmin).and.(dabs(eta(ij)).lt.etaJmax)) then
            num_vis_jets=num_vis_jets+1
          endif
        enddo
        
        if(decayed) then
          ! calculate lepton kinematical quantities:
          ! lepton A, B transverse momentum
          pt_decA = sqrt(p_decA(1)**2 + p_decA(2)**2)
          pt_decB = sqrt(p_decB(1)**2 + p_decB(2)**2)
          ! rapidity of lepton A, B
          call getrapidity(p_decA,y_decA)
          ! pseudorapidity of lepton A, B
          call get_pseudorap(p_decA,eta_decA)
          call get_pseudorap(p_decB,eta_decB)
          ! azimuth of lepton A, B
          call getazimuth(p_decA,phi_decA)
          call getazimuth(p_decB,phi_decB)
          ! azimuthal difference of the leptons
          phi_decAdecB = phi_decB-phi_decA
          if(phi_decAdecB.lt.(-pi)) phi_decAdecB = phi_decAdecB+2d0*pi
          if(phi_decAdecB.gt.pi)    phi_decAdecB = phi_decAdecB-2d0*pi
        endif
        
        ! no cuts
        call filld("sigma_total_nocuts",0D0,dsig)
        call filld("inv_mass_pair_nocuts",m_finAfinB,dsig)
        call filld("pt_finA_nocuts",pt_finA,dsig)
        call filld("pt_finB_nocuts",pt_finB,dsig)
        call filld("eta_finA_nocuts",eta_finA,dsig)
        call filld("number_of_visible_jets_nocuts",
     &             dble(num_vis_jets),dsig)
        call filld("pt_jet_nocuts",kt(1),dsig)
        call filld("eta_jet_nocuts",eta(1),dsig)

        ! determine if jet cuts have passed
        JcutsPassed = .true.
        if(mjets.lt.1) JcutsPassed = .false.
        if(kt(1).le.ptJmin) JcutsPassed = .false.
        if(dabs(eta(1)).ge.etaJmax) JcutsPassed = .false.
        
        ! fill histograms if jet cuts have passed
        if(JcutsPassed) then
          call filld("sigma_total_jcuts",0D0,dsig)
          call filld("inv_mass_pair_jcuts",m_finAfinB,dsig)
          call filld("pt_finA_jcuts",pt_finA,dsig)
          call filld("pt_finB_jcuts",pt_finB,dsig)
          call filld("eta_finA_jcuts",eta_finA,dsig)
          call filld("number_of_visible_jets_jcuts",
     &               dble(num_vis_jets),dsig)
          call filld("pt_jet_jcuts",kt(1),dsig)
          call filld("eta_jet_jcuts",eta(1),dsig)
        endif
        
        ! determine if lepton cuts have passed
        Lcutspassed = .true.
        if(decayed.eqv..false.) Lcutspassed = .false.
        if(pt_decA.le.ptLmin) Lcutspassed = .false.
        if(pt_decB.le.ptLmin) Lcutspassed = .false.
        if(dabs(eta_decA).ge.etaLmax) Lcutspassed = .false.
        if(dabs(eta_decB).ge.etaLmax) Lcutspassed = .false.
        
        ! cutsL
        if(Lcutspassed) then
          call filld('sigma_total_lcuts',0d0,dsig)
          call filld('azimuth_diff_lepton_pair_lcuts',phi_decAdecB,dsig)
          call filld('pt_lepton_A_lcuts',pt_decA,dsig)
          call filld('eta_lepton_A_lcuts',eta_decA,dsig)
          call filld('pt_lepton_B_lcuts',pt_decB,dsig)
          call filld('eta_lepton_B_lcuts',eta_decB,dsig)
        endif
      end
c############### end subroutine analysis ###############################

c############### subroutine getnumber ##################################
c takes PDG ID
c returns number of undecayed particles with this ID in n
      subroutine getnumber(PDG_ID,n)
        implicit none
#include "hepevt.h"
        integer PDG_ID, n, ihep

        n=0
        do ihep=1,nhep
          if(isthep(ihep).eq.1) then
            if(idhep(ihep).eq.PDG_ID) then
              n=n+1
            endif
          endif
        enddo
      
#ifdef DEBUGQ
        print*,"number of particles with ID",PDG_ID,": ",n
#endif
      end
c############### end subroutine getnumber ##############################

c############### subroutine get4momentum_largestpt #####################
c takes PDG ID
c returns four-momentum of undecayed particle with this ID and largest pt
      subroutine get4momentum_largestpt(PDG_ID,p)
        implicit none
        integer PDG_ID, n, nmax, mu
        parameter (nmax=200)
        integer list(nmax)
        double precision p(0:3)

        do mu=0,3
          p(mu) = 0d0
        enddo

        call getptsortedlist(PDG_ID,n,list)
        if (list(1).ne.0) then  ! only if particle with PDG_ID is found
          ! take four-momentum of first particle in list (largest pt)
          call get4momentum_ihep(list(1),p)
        endif
      end
c############### end subroutine get4momentum_largestpt #################

c############### subroutine getptsortedlist ############################
c takes PDG ID
c returns number of undecayed particles with this ID in n
c returns list with ihep numbers of these particles sorted by pt in list
      subroutine getptsortedlist(PDG_ID,n,list)
        implicit none
#include "hepevt.h"
        integer PDG_ID, n, nmax, ihep
        parameter (nmax=200)
        integer list(nmax)

        do n=1,nmax
          list(n) = 0
        enddo

        ! determine number of undecayed particles with PDG_ID in event
        n = 0
        do ihep=1,nhep
          if(isthep(ihep).eq.1) then
            if(idhep(ihep).eq.PDG_ID) then
              n = n+1
            endif
          endif
        enddo
        ! check that maximal number of allowed particles is not exceeded
        if(n.gt.nmax) then
          write(*,*) 'more than', nmax, 'particles with ID', PDG_ID
          write(*,*) 'could not fill pt sorted list'
        else
          ! fill list with found particles
          n = 0
          do ihep=1,nhep
            if(isthep(ihep).eq.1) then
              if(idhep(ihep).eq.PDG_ID) then
                n = n+1
                list(n)=ihep
              endif
            endif
          enddo
          ! sort by pt
          call sortbypt(n,list)
        endif
      end
c############### end subroutine getptsortedlist ########################

c############### subroutine sortbypt ###################################
c sorts array iarr of ihep numbers by pt (decreasing)
      subroutine sortbypt(n,iarr)
        implicit none
#include "hepevt.h"
        integer n,iarr(n)
        integer j,k
        double precision tmp,pt(nmxhep)
        logical touched
        
        do j=1,n
          pt(j) = sqrt(phep(1,iarr(j))**2+phep(2,iarr(j))**2)
        enddo
        ! bubble sort
        touched = .true.
        do while(touched)
          touched = .false.
          do j=1,n-1
            if(pt(j).lt.pt(j+1)) then
              k = iarr(j)
              iarr(j) = iarr(j+1)
              iarr(j+1) = k
              tmp = pt(j)
              pt(j) = pt(j+1)
              pt(j+1) = tmp
              touched = .true.
            endif
          enddo
        enddo
      end
c############### end subroutine sortbypt ###############################

c############### subroutine get4momentum_ihep ##########################
c takes ihep number of a particle in the HEPEVT common block
c puts its 4-momentum into p
c the phep momenta is in the laboratory frame
      subroutine get4momentum_ihep(ihep,p)
        implicit none
#include "hepevt.h"
#include "nlegborn.h"
#include "pwhg_kn.h"
        integer ihep, mu
        double precision p(0:3)

        do mu=0,3
          p(mu)=0D0
        enddo
      
        if(ihep.lt.0 .or. ihep.gt.nhep) then
          print*, "no particle with ihep number", ihep,
     &              "in current event"
          print*, "could not determine four-momentum"
          stop
        else
          ! p is defined in the laboratory frame (~kn_pborn(i))
          p(0) = phep(4,ihep)
          do mu=1,3
            p(mu) = phep(mu,ihep)
          enddo
        endif
      
#ifdef DEBUGQ
        print*,"p(",ihep,"): ",p(:)
        print*,"phep(1): ", phep(4,1),phep(1,1),phep(2,1),phep(3,1)
        print*,"phep(2): ", phep(4,2),phep(1,2),phep(2,2),phep(3,2)
        print*,"phep(3): ", phep(4,3),phep(1,3),phep(2,3),phep(3,3)
        print*,"phep(4): ", phep(4,4),phep(1,4),phep(2,4),phep(3,4)
        print*,"phep(5): ", phep(4,5),phep(1,5),phep(2,5),phep(3,5)
        print*,"kn_preal(1): ",kn_preal(:,1)
        print*,"kn_preal(2): ",kn_preal(:,2)
        print*,"kn_preal(3): ",kn_preal(:,3)
        print*,"kn_preal(4): ",kn_preal(:,4)
        print*,"kn_preal(5): ",kn_preal(:,5)
        print*,"kn_pborn(1): ",kn_pborn(:,1)
        print*,"kn_pborn(2): ",kn_pborn(:,2)
        print*,"kn_pborn(3): ",kn_pborn(:,3)
        print*,"kn_pborn(4): ",kn_pborn(:,4)
#endif
      end
c############### end subroutine get4momentum_ihep ######################

c############### subroutine get4momentum_fin_pair ######################
c takes particle data group identity code of the final state pair
c puts 4-momentum of first and second particle
c in HEPEVT common block with this ID into p1 and p2
c this routine can be used even if the IDs of particle 1 and 2 are the same
c (this was not possible in dislepton_jet)
      subroutine get4momentum_fin_pair(PDG_ID1,p1,PDG_ID2,p2)
        implicit none
#include "hepevt.h"
#include "PhysPars.h"
        integer PDG_ID1,PDG_ID2
        ! output
        double precision p1(0:3), p2(0:3)
        double precision momsq
        external momsq
        ! indices
        integer ihep, imu
        ! store id of found particle?
        integer found1, found2
        ! check on-shell condition
        double precision relerror

        do imu=0,3
          p1(imu) = 0D0
          p2(imu) = 0D0
        enddo
        found1 = -1
        found2 = -1

        do ihep=1,nhep
          if(.not.(found2 .eq. found1)) then
            if(PDG_ID2.eq.idhep(ihep)) then
              call get4momentum_ihep(ihep,p2)
              found2=ihep
            endif
          endif
          if(found1 .eq. -1) then
            if(PDG_ID1.eq.idhep(ihep)) then
              call get4momentum_ihep(ihep,p1)
              found1=ihep
            endif
          endif
        enddo
        
#ifdef DEBUGQ
        print*,"PDG_ID1",PDG_ID1
        print*,"p1",p1
        print*,"PDG_ID2",PDG_ID2
        print*,"p2",p2
        ! check on-shell condition
        print*,"m1",dsqrt(momsq(p1)),par_Fin1mass
        print*,"m2",dsqrt(momsq(p2)),par_Fin2mass
        !stop
#endif

        ! check if the extracted particles are indeed the final state particles
        relerror = (dsqrt(momsq(p1))-par_Fin1mass)
     &             /(dsqrt(momsq(p1))+par_Fin1mass)
        if( dabs(relerror) .gt. 1D-6 ) then
          print*,"error: wrong particle extracted in "
          print*,"subroutine get4momentum_fin_pair"
          print*,"p1**2 = ",momsq(p1)
          print*,"m1**2 = ",par_Fin1mass**2
          stop
        endif
        relerror = (dsqrt(momsq(p2))-par_Fin2mass)
     &             /(dsqrt(momsq(p2))+par_Fin2mass)
        if( dabs(relerror) .gt. 1D-6 ) then
          print*,"error: wrong particle extracted in "
          print*,"subroutine get4momentum_fin_pair"
          print*,"p2**2 = ",momsq(p2)
          print*,"m2**2 = ",par_Fin2mass**2
          stop
        endif

        if(found1 .eq. found2) then
          print*, "error in get4momentum_fin_pair."
          print*, "idhep is the same for particle 1 and particle 2,"
          print*, "but should be different."
          print*, "idhep1", found1
          print*, "idhep2", found2
          stop
        endif
        
        if(found1 .eq. -1) then
          print*, "particle with PDG ID", PDG_ID1, "not found"
          print*, "could not determine four-momentum"
          stop
        endif
        
        if(found2 .eq. -1) then
          print*, "particle with PDG ID", PDG_ID2, "not found"
          print*, "could not determine four-momentum"
          stop
        endif
      end
c############### end subroutine get4momentum_fin_pair ##################

c############### subroutine getazimuth #################################
c calculates azimuth phi in [0,2pi] from x-axis
      subroutine getazimuth(p,phi)
        implicit none
#include "pwhg_math.h"
        double precision p(0:3),phi
        phi=0D0
        if(p(1).GT.0D0) then
          if(p(2).GE.0D0) phi = datan(p(2)/p(1))
          if(p(2).LT.0D0) phi = datan(p(2)/p(1))+2D0*pi
        elseif(p(1).LT.0D0) then
          phi = datan(p(2)/p(1))+pi
        elseif(p(1).EQ.0D0) then
          if(p(2).GT.0D0) phi = pi/2D0
          if(p(2).LT.0D0) phi = -pi/2D0
        endif
      end
c############### end subroutine getazimuth #############################

c############### subroutine getrapidity ################################
      subroutine getrapidity(p,y)
        implicit none
        double precision p(0:3),y
        y = 0.5D0*dlog((p(0)+p(3))/(p(0)-p(3)))
      end
c############### end subroutine getrapidity ############################

c############### subroutine getinvmass #################################
      subroutine getinvmass(p,m)
        implicit none
        double precision p(0:3),m
        m = dsqrt(dabs((p(0)**2-p(3)**2)-p(1)**2-p(2)**2))
      end
c############### end subroutine getinvmass #############################

c############### subroutine get_pseudorap ##############################
      subroutine get_pseudorap(p,eta)
        implicit none
        double precision p(0:3),eta,pt,th
        double precision eps
        parameter (eps=1.d-5)
        pt = dsqrt(p(1)**2+p(2)**2)
        if(pt.lt.eps.and.dabs(p(3)).lt.eps)then
          eta = sign(1.D0,p(3))*1.d8
        elseif(pt.lt.eps) then   ! MK: added this elseif
          eta = sign(1.D0,p(3))*1.d8
        else
          th  = datan2(pt,p(3))
          eta = -dlog(dtan(th/2.D0))
        endif
      end
c############### end subroutine get_pseudorap ##########################

c############### subroutine isvaliddecayspectrum #######################
c checks mass spectrum to ensure that the extracted decay for the 
c final states is Xi -> f + ~fj
c assumes MADGRAPH parameters to be filled already
      logical function isvaliddecayspectrum(id_finA,id_finB,
     &                                      id_decA,id_decB)
        implicit none
#include "PhysPars.h"
        ! IDs
        integer id_finA, id_finB, id_decA, id_decB
        ! masses
        double precision m_finA, m_finB
         ! sum of masses of left and right handed decay products
        double precision m_dec(2)

        isvaliddecayspectrum = .true.
        ! invalid default values
        m_finA = invalid
        m_finB = invalid

        ! Note: neutralino masses can be negative
        ! set final state masses
        m_finA = dabs(par_Fin1mass)
        m_finB = dabs(par_Fin2mass)
        
        ! masses of decay products for decay product A
        select case(abs(id_decA))
        case (11) ! electron/positron
          m_dec(1) = 0D0 ! neglected
          m_dec(2) = 0D0
        case (13) ! muon
          m_dec(1) = 0D0 ! neglected
          m_dec(2) = 0D0
        case (15) ! tau
          m_dec(1) = 1.77682d0  ! tau mass (PDG 2012)
          m_dec(2) = 1.77682d0  ! tau mass (PDG 2012)
        case default
          write(*,*) 'encountered unsupported decay ID ', id_decA
          stop
        end select
        
        ! masses of decay products for decay product A
        select case(abs(id_decB))
        case (11) ! electron/positron
          m_dec(1) = 0D0 ! neglected
          m_dec(2) = 0D0
        case (13) ! muon
          m_dec(1) = 0D0 ! neglected
          m_dec(2) = 0D0
        case (15) ! tau
          m_dec(1) = 1.77682d0  ! tau mass (PDG 2012)
          m_dec(2) = 1.77682d0  ! tau mass (PDG 2012)
        case default
          write(*,*) 'encountered unsupported decay ID ', id_decB
          stop
        end select
        
        
        ! masses of decay products for final state A
        select case(abs(id_finA))
        case (1000022)                          ! n1
          m_dec(:) = par_MSf(:,2,1) + m_dec(:)  ! selectron + lepton mass
        case (1000023)                          ! n2
          m_dec(:) = par_MSf(:,2,1) + m_dec(:)  ! selectron + lepton mass
        case (1000025)                          ! n3
          m_dec(:) = par_MSf(:,2,1) + m_dec(:)  ! selectron + lepton mass
        case (1000035)                          ! n4
          m_dec(:) = par_MSf(:,2,1) + m_dec(:)  ! selectron + lepton mass
        case (1000024)                          ! x1
          m_dec(:) = par_MSf(:,1,1) + m_dec(:)  ! sneutrino + lepton mass
        case (1000037)                          ! x2
          m_dec(:) = par_MSf(:,1,1) + m_dec(:)  ! sneutrino + lepton mass
        case default
          write(*,*) 'encountered unknown final state ID ', id_finA
          stop
        end select
      
        ! masses of decay products for final state B
        select case(abs(id_finB))
        case (1000022)                          ! n1
          m_dec(:) = par_MSf(:,2,1) + m_dec(:)  ! selectron + lepton mass
        case (1000023)                          ! n2
          m_dec(:) = par_MSf(:,2,1) + m_dec(:)  ! selectron + lepton mass
        case (1000025)                          ! n3
          m_dec(:) = par_MSf(:,2,1) + m_dec(:)  ! selectron + lepton mass
        case (1000035)                          ! n4
          m_dec(:) = par_MSf(:,2,1) + m_dec(:)  ! selectron + lepton mass
        case (1000024)                          ! x1
          m_dec(:) = par_MSf(:,1,1) + m_dec(:)  ! sneutrino + lepton mass
        case (1000037)                          ! x2
          m_dec(:) = par_MSf(:,1,1) + m_dec(:)  ! sneutrino + lepton mass
        case default
          write(*,*) 'encountered unknown final state ID ', id_finB
          stop
        end select
        
        ! check if final state particle A can decay in lepton
        ! if not, the mass spectrum for this analysis is invalid
        ! and the lepton-lepton histograms are empty
        if ( (m_finA.lt.m_dec(1)) .and.
     &       (m_finA.lt.m_dec(2)) ) then
          print*, "decay for final state A into lepton not possible"
          isvaliddecayspectrum = .false.
        endif  
         
        ! same for final state B
        if ( (m_finB.lt.m_dec(1)) .and.
     &       (m_finB.lt.m_dec(2)) ) then
          print*, "decay for final state B into lepton not possible"
          isvaliddecayspectrum = .false.
        endif 

        ! if spectrum failed check write message
        if (.not.(isvaliddecayspectrum)) then
          print*, "spectrum not suitable for allowed decay"
          print*, "actual decays in events generated by PYTHIA"
        endif

      end
c############### end subroutine isvaliddecayspectrum ###################
      

c############### subroutine get_pseudorap ##############################
c jets are orderd by decreasing pt
      subroutine buildjets(mjets,kt,eta,rap,phi,pjet,ptrel)
        implicit none
#include   "hepevt.h"
        ! input
        integer mjets
        ! output
        double precision kt(mjets),eta(mjets),rap(mjets),phi(mjets)
        double precision pjet(4,mjets),ptrel(mjets)
        integer   maxtrack,maxjet
        parameter (maxtrack=2048,maxjet=2048)
        double precision  ptrack(4,maxtrack),pj(4,maxjet)
        integer   jetvec(maxtrack),itrackhep(maxtrack)
        integer   ntracks,njets
        integer   j,k,mu
        double precision palg,pp,tmp
        double precision ptJmin,etaJmax,R
        common /jetcuts/ptJmin,etaJmax,R
        logical isjetparticle

        ! Initialize arrays and counters for output jets
        do j=1,maxtrack
          do mu=1,4
            ptrack(mu,j) = 0D0
          enddo
          jetvec(j) = 0
        enddo
        ntracks = 0
        do j=1,mjets
          do mu=1,4
            pjet(mu,j) = 0D0
            pj(mu,j) = 0D0
          enddo
        enddo
        ! Extract final state particles to feed to jet finder
        do j=1,nhep
#ifdef DEBUGQ
#define DEBUG1
          print*,"idhep = ", idhep(j)
#endif
          if (isthep(j).eq.1 .and. isjetparticle(idhep(j))) then
            if(ntracks.eq.maxtrack) then
              print*, "analyze: need to increase maxtrack!"
              print*, "ntracks: ",ntracks
              stop
            endif
            ntracks = ntracks+1
            itrackhep(ntracks) = j
            do mu=1,4
              ptrack(mu,ntracks) = phep(mu,j)
            enddo
          endif
        enddo
#ifdef DEBUGQ
#define DEBUG1
        print*,"ntracks = ", ntracks
#endif
        if (ntracks.eq.0) then
          mjets=0
          return
        endif

        ! note: ptJmin and R are in common block "jetcuts"
        ! anti-kT (c.f. ZZ code):
        palg=-1
        call fastjetppgenkt(ptrack,ntracks,R,palg,ptJmin,pjet,njets,
     &                      jetvec)

#ifdef DEBUG1
        print*,"phep(1): ", phep(1:4,1)
        print*,"phep(2): ", phep(1:4,2)
        print*,"phep(3): ", phep(1:4,3)
        print*,"phep(4): ", phep(1:4,4)
        print*,"phep(5): ", phep(1:4,5)
        print*,"ptrack: ", ptrack(:,1)
        print*,"ntracks: ", ntracks
        print*,"R: ", R
        print*,"palg: ", palg
        print*,"ptJmin: ", ptJmin
        print*,"pjet: ", pjet(:,1)
        print*,"jetvec: ", jetvec(1)
        print*,"mjets: ", mjets
        print*,"njets: ", njets
        stop
#endif
     
        mjets = min(mjets,njets)
        if(njets.eq.0) return
        ! check consistency
        do k=1,ntracks
          if(jetvec(k).gt.0) then
            do mu=1,4
              pj(mu,jetvec(k)) = pj(mu,jetvec(k))+ptrack(mu,k)
            enddo
          endif
        enddo
        tmp = 0
        do j=1,mjets
          do mu=1,4
            tmp = tmp+dabs(pj(mu,j)-pjet(mu,j))
          enddo
        enddo
        if(tmp.gt.1d-4) then
          print*, "error in buildjets:894"
        endif
        ! end check consistency

        ! computing arrays of useful kinematics quantities for hardest jets
        do j=1,mjets
          kt(j) = dsqrt(pjet(1,j)**2+pjet(2,j)**2)
          pp = dsqrt(kt(j)**2+pjet(3,j)**2)
          eta(j) = 0.5D0*dlog((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
          rap(j) = 0.5D0*dlog((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
          phi(j) = datan2(pjet(2,j),pjet(1,j))
        enddo

        call computeptrel(ptrack,ntracks,rap,kt,phi,mjets,jetvec,ptrel)

      end
c############### end subroutine get_pseudorap ##########################

c############### subroutine isjetparticle ##############################
c determines if particle contributes to a jet
      logical function isjetparticle(PDG_ID)
        implicit none
        integer PDG_ID
        logical injetwhitelist, injetblacklist
        isjetparticle = .false.
        ! find out if particle is in
        ! jetwhitelist: particles that contribute to a jet
        ! jetblacklist: particles that do not contribute to a jet
        if (injetwhitelist(PDG_ID)) then
          isjetparticle = .true.
        else if (injetblacklist(PDG_ID)) then
          isjetparticle = .false.
        else
          print*, "could not determine if particle with ID", PDG_ID
          print*, "contributes to a jet"
          stop
        endif
      end
c############### end subroutine isjetparticle ##########################

c############### subroutine injetwhitelist #############################
c determines if particle is in whitelist
c of particles that contribute to a jet
      logical function injetwhitelist(PDG_ID)
        implicit none
        integer PDG_ID, id
        injetwhitelist=.false.
        id=iabs(PDG_ID)
        ! the following particles contribute to a jet
        ! quarks
        if (1.le.id .and. id.le.6) injetwhitelist = .true.
        ! photon is not in whitelist
        if (id.eq.22) injetwhitelist = .false.
        ! gluon
        if (id.eq.21) injetwhitelist = .true.
        if (id.eq.0) injetwhitelist = .true.  ! POWHEG numbering
        ! diquarks
        if ((1000.le.id .and. id.le.6000)
     &      .and. (mod(id,100).eq.1 .or. mod(id,100).eq.3)) then
          injetwhitelist = .true.
        endif
        ! mesons
        if ((111.le.id .and. id.le.600)
     &      .or. (10000.le.id .and. id.le.999999)
     &      .or. (9000000.le.id .and. id.le.9100000)) then
          injetwhitelist = .true.
        endif
        ! baryons
        if ((1000.le.id .and. id.le.6000)
     &      .and. (mod(id,10).eq.2 .or. mod(id,10).eq.4)) then
          injetwhitelist = .true.
        endif
      end
c############### end subroutine injetwhitelist #########################

c############### subroutine injetblacklist #############################
c     determines if particle is in blacklist
c     of particles that do not contribute to a jet
      logical function injetblacklist(PDG_ID)
        implicit none
        integer PDG_ID, id
        injetblacklist=.false.
        id=iabs(PDG_ID)
        ! the following particles may occur
        ! but do not contribute to a jet
        ! leptons
        if (id.eq.11) injetblacklist=.true.  ! e
        if (id.eq.12) injetblacklist=.true.  ! nue
        if (id.eq.13) injetblacklist=.true.  ! mu
        if (id.eq.14) injetblacklist=.true.  ! numu
        if (id.eq.15) injetblacklist=.true.  ! tau
        if (id.eq.16) injetblacklist=.true.  ! nu_tau
        ! sleptons
        if (id.eq.1000011) injetblacklist=.true.  ! ~e_L
        if (id.eq.1000012) injetblacklist=.true.  ! ~snue_L
        if (id.eq.1000013) injetblacklist=.true.  ! ~mu_L
        if (id.eq.1000014) injetblacklist=.true.  ! ~snumu_L
        if (id.eq.1000015) injetblacklist=.true.  ! ~stau_1
        if (id.eq.1000016) injetblacklist=.true.  ! ~snu_tau_L
        if (id.eq.2000011) injetblacklist=.true.  ! ~e_R
        if (id.eq.2000013) injetblacklist=.true.  ! ~mu_R
        if (id.eq.2000015) injetblacklist=.true.  ! ~stau_2
        ! fintralinos
        if (id.eq.1000022) injetblacklist=.true.  ! ~n(1)
        if (id.eq.1000023) injetblacklist=.true.  ! ~n(2)
        if (id.eq.1000025) injetblacklist=.true.  ! ~n(3)
        if (id.eq.1000035) injetblacklist=.true.  ! ~n(4)
        ! charginos
        if (id.eq.1000024) injetblacklist=.true.  ! ~x(1)
        if (id.eq.1000037) injetblacklist=.true.  ! ~x(2)
        ! photon
        if (id.eq.22) injetblacklist=.true.
      end
c############### end subroutine injetblacklist #########################

c############### subroutine computeptrel ###############################
c pt rel is sum of the ptrack momentum projection ortogonal to the jet
c momentum in the frame where the jet has zero rapidity
      subroutine computeptrel(ptracks,ntracks,rapjets,ktjets,phijets,
     &     njets,jetvec,ptrel)
        implicit none
        integer ntracks,njets,jetvec(ntracks)
        double precision ptracks(4,ntracks),rapjets(njets)
        double precision ktjets(njets),phijets(njets),ptrel(njets)
        integer j,i
        double precision yj,kj1,kj2,y,pt(3)
        do j=1,njets
          ptrel(j)=0
        enddo
        do i=1,ntracks
          j = jetvec(i)
          if(j.gt.0.and.j.le.njets) then
            ! Track i belongs to jet j
            yj = rapjets(j)
            kj1 = ktjets(j)*dcos(phijets(j))
            kj2 = ktjets(j)*dsin(phijets(j))
            ! rapidity of track i
            y = 0.5D0*dlog((ptracks(4,i)+ptracks(3,i))
     &                 /(ptracks(4,i)-ptracks(3,i)))
            ! rapidity of track i in frame where the jet has zero rapidity
            y = y-yj
            ! find momentum of track i in frame where the jet has zero rapidity
            pt(1) = ptracks(1,i)
            pt(2) = ptracks(2,i)
            pt(3) = dsqrt(pt(1)**2+pt(2)**2)*sinh(y)
            ! pt rel is sum of the ptrack momentum projection ortogonal to the jet
            ! momentum in the frame where the jet has zero rapidity
            ptrel(j) = dsqrt(((pt(1)*kj2-pt(2)*kj1)**2+
     &                     (         -pt(3)*kj2)**2+
     &                     (pt(3)*kj1          )**2)/
     &                     (kj1**2+kj2**2)) + ptrel(j)
          endif
        enddo
      end
c############### end subroutine computeptrel ###########################