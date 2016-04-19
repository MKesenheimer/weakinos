c############### subroutine init_hist ##################################
c book all histograms, will be filled later
      subroutine init_hist
        implicit none
#include "LesHouches.h"
#include "pwhg_math.h"
#include "PhysPars.h"
        ! tell this analysis file which program is running it
        character*6 WHCPRG
        common/cWHCPRG/WHCPRG
        integer ndecay,i
        parameter (ndecay = 4)
        ! external functions
        double precision powheginput
        external powheginput
        double precision mm, mp, ptm, ptp, ptjm, ptjp, em, ep, ym, yp
        double precision mdm, mdp, ptdm, ptdp, edm, edp
        double precision phim, phip
        integer n1, n2, n3
        character*4 cstr(ndecay)
        data cstr/"dec1","dec2","dec3","dec4"/
        ! correlations
        integer corI, corJ
        common /icorr/ corI,corJ

        call inihists

        ! check program that runs analysis
        print*,
        print*, "===================================================="
        if(whcprg.eq."NLO") then
          print*, "fixed-order NLO analysis"
        elseif(whcprg.eq."LHE") then
          print*, "LHE analysis"
          call init_couplings
        elseif(whcprg.eq."HERWIG") then
          print*, "HERWIG analysis"
          print*, "not implemented analysis"
          print*, "no plots will be present at the end of the run"
          call init_couplings
        elseif(whcprg.eq."PYTHIA") then
          print*, "PYTHIA analysis"
          call init_couplings
        endif
        print*, "===================================================="

        ! for which particles do we want to calculate correlations?
        ! corI and corJ is saved in the common block icorr
        if(powheginput("#corI").gt.0) then
          corI = powheginput("corI")
        endif
        if(powheginput("#corJ").gt.0) then
          corJ = powheginput("corJ")
        endif
          
        ! number of bins
        n1 = 80
        ! mass(A+B), binwidth = 10GeV
        mm = dabs(par_Fin1mass) + dabs(par_Fin2mass)
        mp = mm + 800D0
        ! Pt(A), Pt(B), Pt(miss), binwidth = 10GeV
        ptm = 0D0
        ptp = ptm + 800D0
        ! Pt(J), binwidth = 2.5GeV
        ptjm = 0D0
        ptjp = ptjm + 200D0
        ! Pt(decI), binwidth = 10GeV
        ptdm = 0D0
        ptdp = ptdm + 800D0
        
        ! number of bins
        n2 = 40
        n3 = 48
        ! eta(A), eta(B), eta(J), binwidth = 0.25
        em = -5D0
        ep = em + 10D0
        ! eta(decI), binwidth = 0.125
        edm = -3D0
        edp = edm + 6D0
        ! y(A), binwidth = 0.2
        ym = -4D0
        yp = ym + 8D0
        ! mass(decI+decJ), binwidth = 2.5GeV
        mdm = 0D0
        mdp = mdm + 200D0
        ! dPhi(decI+decJ), binwidth = 2pi/40 ~ 0.15708
        phim = -pi
        phip = pi

        ! book all possible histograms although we don't need
        ! all of them for every process

        ! no cuts
        call bookupeqbins("xsec",1D0,-0.5D0,0.5D0)
        call bookupeqbins("mass(A+B)",(mp-mm)/n1,mm,mp)
        call bookupeqbins("Pt(A)",(ptp-ptm)/n1,ptm,ptp)
        call bookupeqbins("Pt(B)",(ptp-ptm)/n1,ptm,ptp)
        call bookupeqbins("eta(A)",(ep-em)/n2,em,ep)
        call bookupeqbins("eta(B)",(ep-em)/n2,em,ep)
        call bookupeqbins("y(A)",(yp-ym)/n2,ym,yp)
        call bookupeqbins("y(B)",(yp-ym)/n2,ym,yp)
        call bookupeqbins("N(J)",1D0,-0.5D0,2.5D0)
        call bookupeqbins("Pt(J)",(ptjp-ptjm)/n1,ptjm,ptjp)
        call bookupeqbins("eta(J)",(ep-em)/n2,em,ep)

        ! jet cuts
        call bookupeqbins("xsec(Jcuts)",1D0,-0.5D0,0.5D0)
        call bookupeqbins("N(J,Jcuts)",1D0,-0.5D0,2.5D0)
        call bookupeqbins("Pt(J,Jcuts)",(ptjp-ptjm)/n1,ptjm,ptjp)
        call bookupeqbins("eta(J,Jcuts)",(ep-em)/n2,em,ep)
     
        ! decays and missing Pt
        if(WHCPRG.eq."PYTHIA") then
          ! no cuts
          call bookupeqbins("dPhi("//cstr(corI)//"+"//cstr(corJ)//")",
     &                     (phip-phim)/n2,phim,phip)
          call bookupeqbins("mass("//cstr(corI)//"+"//cstr(corJ)//")",
     &                     (mdp-mdm)/n1,mdm,mdp)
          do i=1,ndecay
            call bookupeqbins("Pt("//cstr(i)//")",
     &                       (ptdp-ptdm)/n1,ptdm,ptdp)
            call bookupeqbins("eta("//cstr(i)//")",
     &                       (edp-edm)/n3,edm,edp)
          enddo
          ! decay cuts
          call bookupeqbins("xsec(Dcuts)",1d0,-0.5d0,0.5d0)
          call bookupeqbins("dPhi("//cstr(corI)//"+"//cstr(corJ)//
     &                     ",Dcuts)",(phip-phim)/n2,phim,phip)
          call bookupeqbins("mass("//cstr(corI)//"+"//cstr(corJ)//
     &                     ",Dcuts)",(mdp-mdm)/n1,mdm,mdp)
          call bookupeqbins("Pt(miss)",(ptp-ptm)/n1,ptm,ptp)
          do i=1,ndecay
            call bookupeqbins("Pt("//cstr(i)//",Dcuts)",
     &                       (ptdp-ptdm)/n1,ptdm,ptdp)
            call bookupeqbins("eta("//cstr(i)//",Dcuts)",
     &                       (edp-edm)/n3,edm,edp)
          enddo
        endif
      end
c############### end subroutine init_hist ##############################

c############### subroutine analysis ###################################
c extract all data required for the histograms
c calculate quantities that shall be plotted and fill histograms
      subroutine analysis(dsig0)
        implicit none
#include "pwhg_math.h"
#include "LesHouches.h"
#include "hepevt.h"
#include "pwhg_weights.h"
        double precision dsig(7), dsig0 ! total cross section
        ! tells if subroutine is called for the first time
        logical ini
        data ini/.true./
        save ini
        ! tell this analysis file which program is running it
        character *6 WHCPRG
        common/cWHCPRG/WHCPRG
        data WHCPRG/'NLO'/ ! default, should be set when analysis is called
        ! external functions
        double precision powheginput, pt, rsepn
        logical isdetectable
        external powheginput, pt, rsepn, isdetectable
        ! kinematical quantities
        double precision p_finA(0:3), p_finB(0:3), p_finAfinB(0:3)
        double precision pt_finA, pt_finB
        double precision eta_finA,eta_finB
        double precision y_finA,y_finB
        double precision m_finAfinB
        double precision pt_finAfinB, y_finAfinB
        ! jet parameters
        integer mjets, num_vis_jets, maxjet
        parameter (maxjet=2048)
        double precision kt(maxjet), eta(maxjet), rap(maxjet)
        double precision phi(maxjet), pj(4,maxjet), ptrel(maxjet)
        ! PDG and HEP IDs
        integer id_finA, id_finB
        integer ihep_finA, ihep_finB
        ! decay list of finA and finB
        integer ihep_decchainA(nhep,0:2), ihep_decchainB(nhep,0:2)
        integer num_decchainA, num_decchainB
        ! the considered decay chains
        character*100 str_user_decchainA, str_user_decchainB
        integer user_decchainA(100,0:2), user_decchainB(100,0:2)
        save user_decchainA, user_decchainB
        integer nmothersA, nmothersB
        save nmothersA, nmothersB
        integer nmatchA, nmatchB
        ! decays
        ! number of decay products to consider
        integer ndecay
        parameter (ndecay = 4)
        ! IDs of decay products
        integer id_dec(ndecay)
        save id_dec
        ! momenta, pt, eta, y of decay products
        double precision p_dec(0:3,ndecay), pt_dec(ndecay)
        double precision eta_dec(ndecay), y_dec(ndecay)
        double precision phi_dec(ndecay), dphi_decIdecJ
        double precision p_decIdecJ(0:3), m_decIdecJ
        ! variables to store which correlation (I,J) we want to consider
        integer corI, corJ
        common /icorr/ corI,corJ
        ! handling the decays
        character*4 cstr(ndecay)
        data cstr/"dec1","dec2","dec3","dec4"/
        ! counter to keep track how many decays and cuts have took place
        ! is used to determine the branching ratios of the decays
        integer cnt_decays_tot, cnt_decays(ndecay), cnt_decays_proc
        save cnt_decays_tot, cnt_decays, cnt_decays_proc
        ! acces the branching corrected ratios of PYTHIA (this is set in setup-PYTHIA-lhef.f)
        double precision bratio
        common/pybratio/bratio
        ! variables for missing energy
        double precision p_miss(0:3), pt_miss, pt_miss1
        ! cuts
        ! ptJmin is a technical cut
        integer mJmin
        double precision ptJmin, etaJmax, R
        common /jetcuts/ptJmin, etaJmax, R
        ! physical Jet Cuts
        double precision ptJminPhys
        ! decay cuts
        double precision RDsepmin
        double precision RDmin(ndecay), RDsep(ndecay,ndecay)
        double precision ptDmin, etaDmax
        double precision ptmissmin
        save ptDmin, etaDmax
        logical JcutsPassed, DcutsPassed
        ! flags for decay
        logical fdecA, fdecB, decayed
        save fdecA, fdecB, decayed
        ! number of events processed
        integer nevents
        ! indices
        integer mu, i, j        
        ! debug
        logical isvetoed
        isvetoed = .true.

        ! at first call
        if (ini) then
          ! reset number of events processed
          nevents = 0
          ! check if multiple weights are present
          if(weights_num.eq.0) then
            call setupmulti(1)
          else
            call setupmulti(weights_num)
#ifdef DEBUG
            print*,"[DEBUG]:195 in Analysis.f multiple weights present, check and uncomment!"
            stop
#endif
          endif
          ! determine the decay products from powheg.input
          ! remember that not every number is needed for every process
          ! it is assumed that dec1 and dec2 are daughters of fin1, and
          ! dec3 and dec4 are daughters of fin2
          do i=1,ndecay
            id_dec(i) = int(powheginput(cstr(i))) ! if zero: don't observe this decay
          enddo
          ! the considered decay chain
          call powheginputstring("chainA",str_user_decchainA)
          call store_decchain(str_user_decchainA,user_decchainA)
          call powheginputstring("chainB",str_user_decchainB)
          call store_decchain(str_user_decchainB,user_decchainB)
          ! count the number of mothers in the user_decchainA
          nmothersA = 0
          nmothersB = 0
          do i=1,nhep
            if(user_decchainA(i,0) .ne. 0) then
              nmothersA = nmothersA + 1
            else
              exit
            endif  
          enddo
          do i=1,nhep
            if(user_decchainB(i,0) .ne. 0) then
              nmothersB = nmothersB + 1
            else
              exit
            endif  
          enddo
          ! init the counters
          cnt_decays_tot = 1 ! to not divide by zero
          cnt_decays_proc = 0
          do i=1,ndecay
            cnt_decays(i) = 0
          enddo
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
          ! if everything is ok print status
          print*, "extracted from process ID: ", lprup(1)
          print*, "the final state pair: ", id_finA, ", ", id_finB
          print*, "the decay products:   ", id_dec(1), ", ", id_dec(2)
          print*, "                      ", id_dec(3), ", ", id_dec(4)
        endif

        ! init cuts (change here if needed)
        ! cuts are mainly taken from arXiv:1501.07110v2
        ! should not be in if(ini)... since we want to make sure
        ! that the variables are always filled
        ! jet cuts
        mJmin    = 1
        ptJmin   = 1D0
        etaJmax  = 4.5D0
        R        = 0.4D0
        mjets    = 40
        ! physical pt_Jet cut (for pt_miss)
        ptJminPhys = 20D0
        ! lepton cuts
        ptDmin   = 10D0 ! physical cut
        etaDmax  = 2.5D0
        RDsepmin = 0.05D0
        ! the missing pt threshold
        ptmissmin = 100D0

        if(ini) then
          ! write out all cuts
          print*, "===================================================="
          print*, "cuts in analysis routines"
          print*
          print*, "cuts (jet):"
          print*, "hardest jet p_t (selection) >",ptJmin," GeV"
          print*, "hardest jet p_t (physical)  >",ptJminPhys," GeV"
          print*, "hardest jet |eta| <",etaJmax
          print*, "radius parameter R = ",R
          print*
          print*, "cuts (leptons:)"
          print*, "hardest lepton   p_t >",ptDmin," GeV"
          print*, "hardest lepton |eta| <",etaDmax
          print*, "R-separation   R_min >",RDsepmin
          print*,
          print*, "p_t-threshold        >",ptmissmin," GeV"
          print*
          print*, "===================================================="
          ini=.false.
        endif

        ! from pico to femto
        dsig0 = dsig0*1000D0

        ! multiple weights analysis
        dsig(:) = 0D0
        if(weights_num.eq.0) then
          dsig(1) = dsig0
        else
          dsig(1:weights_num) = weights_val(1:weights_num)
        endif
        ! return if no entries in dsig
        if(sum(abs(dsig)).eq.0) return

        ! find final state 4-momenta
        call get4momentum_fin_pair(id_finA,p_finA,id_finB,p_finB)

        ! calculate final state kinematical quantities:
        ! final state A, B transverse momentum
        pt_finA = pt(p_finA)
        pt_finB = pt(p_finB)
        ! final states pseudorapidity
        call getpseudorap(p_finA,eta_finA)
        call getpseudorap(p_finB,eta_finB)
        ! final states rapidity
        call getrapidity(p_finA,y_finA)
        call getrapidity(p_finB,y_finB)
        ! invariant mass of the final state system
        do mu=0,3
          p_finAfinB(mu) = p_finB(mu)+p_finA(mu)
        enddo
        call getinvmass(p_finAfinB,m_finAfinB)
        ! transverse momentum of the final state system
        pt_finAfinB = pt(p_finAfinB)
        ! rapidity of the final state system
        call getrapidity(p_finAfinB,y_finAfinB) 
        
        ! build jets
        ! initialize
        do j = 1,maxjet
          kt(j)    = 0D0
          eta(j)   = 0D0
          rap(j)   = 0D0
          phi(j)   = 0D0
          pj(:,j)  = 0D0
          ptrel(j) = 0D0
        enddo
        call buildjets(mjets,kt,eta,rap,phi,pj,ptrel)
        ! determine number of visible jets
        num_vis_jets = 0
        do i = 1,mjets
          if ((kt(i).gt.ptJmin).and.(dabs(eta(i)).lt.etaJmax)) then
            num_vis_jets=num_vis_jets+1
          endif
        enddo
        
        ! find out if allowed decay took place
        call isdecayed(id_finA,fdecA)
        call isdecayed(id_finB,fdecB)
        if (fdecA.and.fdecB) then
          decayed = .true.
          cnt_decays_tot = cnt_decays_tot + 1
        else
          decayed = .false.
        endif
        
        ! handle the decay products
        if(decayed) then
          ! store the whole decay chain of the final states in an array 
          ! and search for the decay products we want
          call getihep(id_finA,2,.true.,ihep_finA)
          call getihep(id_finB,2,.false.,ihep_finB)
          ! first decay
          ihep_decchainA(1,0) = ihep_finA
          ihep_decchainA(1,1) = jdahep(1,ihep_decchainA(1,0))
          ihep_decchainA(1,2) = jdahep(2,ihep_decchainA(1,0))
          ihep_decchainB(1,0) = ihep_finB
          ihep_decchainB(1,1) = jdahep(1,ihep_decchainB(1,0))
          ihep_decchainB(1,2) = jdahep(2,ihep_decchainB(1,0))
          num_decchainA = 0
          num_decchainB = 0
          ! further decays
          do i=2,nhep
            ihep_decchainA(i,0) = ihep_decchainA(i/2,1+mod(i,2))
            ihep_decchainA(i,1) = jdahep(1,ihep_decchainA(i,0))
            ihep_decchainA(i,2) = jdahep(2,ihep_decchainA(i,0))
            ! count how many "zero" decays took place
            if(ihep_decchainA(i,0).eq.0) then
              num_decchainA = num_decchainA + 1
            else
              num_decchainA = 0
            endif
          enddo
          do i=2,nhep
            ihep_decchainB(i,0) = ihep_decchainB(i/2,1+mod(i,2))
            ihep_decchainB(i,1) = jdahep(1,ihep_decchainB(i,0))
            ihep_decchainB(i,2) = jdahep(2,ihep_decchainB(i,0))
            ! count how many "zero" decays took place
            if(ihep_decchainB(i,0).eq.0) then
              num_decchainB = num_decchainB + 1
            else
              num_decchainB = 0
            endif  
          enddo
          
          ! reset
          p_dec(:,:) = 0D0
          nmatchA = 0
          nmatchB = 0
          
          do i=1,nhep-num_decchainA
          do j=1,nmothersA
            ! if the mothers are the same, check the daughters
            if(idhep(ihep_decchainA(i,0)).eq.user_decchainA(j,0)) then
              if((idhep(ihep_decchainA(i,1)).eq.user_decchainA(j,1).and.
     &            idhep(ihep_decchainA(i,2)).eq.user_decchainA(j,2)).or.
     &           (idhep(ihep_decchainA(i,1)).eq.user_decchainA(j,2).and.
     &            idhep(ihep_decchainA(i,2)).eq.user_decchainA(j,1)))
     &           then
                   nmatchA = nmatchA + 1
              endif  
            endif  
          enddo
          enddo

          do i=1,nhep-num_decchainB
          do j=1,nmothersB
            ! if the mothers are the same, check the daughters
            if(idhep(ihep_decchainB(i,0)).eq.user_decchainB(j,0)) then
              if((idhep(ihep_decchainB(i,1)).eq.user_decchainB(j,1).and.
     &            idhep(ihep_decchainB(i,2)).eq.user_decchainB(j,2)).or.
     &           (idhep(ihep_decchainB(i,1)).eq.user_decchainB(j,2).and.
     &            idhep(ihep_decchainB(i,2)).eq.user_decchainB(j,1)))
     &           then
                   nmatchB = nmatchB + 1
              endif  
            endif  
          enddo
          enddo
          
          ! now check if the requested particles are in the decay chain.
          ! remember: dec1 and dec2 are assumed to be daughters of finA,
          ! dec3 and dec4 are assumed to be daughters of finB
          do i=1,nhep-num_decchainA
            ! if the momentum of the first decay particle was not already
            ! assigned and the ihep pdg number is equal to the 
            ! requested decay number id_dec
            if(p_dec(0,1).eq.0D0 .and.
     &         idhep(ihep_decchainA(i,0)).eq.id_dec(1)) then
              call get4momentum_ihep(ihep_decchainA(i,0),p_dec(:,1))
              isvetoed = .false.
            endif
            ! same for second decay product
            if(p_dec(0,2).eq.0D0 .and.
     &         idhep(ihep_decchainA(i,0)).eq.id_dec(2)) then
              call get4momentum_ihep(ihep_decchainA(i,0),p_dec(:,2))
              isvetoed = .false.
            endif
          enddo
          ! decay chain of final state particle B
          do i=1,nhep-num_decchainB
            ! if the momentum of the first decay particle was not already
            ! assigned and the ihep pdg number is equal to the 
            ! requested decay number id_dec
            if(p_dec(0,3).eq.0D0 .and.
     &         idhep(ihep_decchainB(i,0)).eq.id_dec(3)) then
              call get4momentum_ihep(ihep_decchainB(i,0),p_dec(:,3))
              isvetoed = .false.
            endif
            ! same for second decay product
            if(p_dec(0,4).eq.0D0 .and.
     &         idhep(ihep_decchainB(i,0)).eq.id_dec(4)) then
              call get4momentum_ihep(ihep_decchainB(i,0),p_dec(:,4))
              isvetoed = .false.
            endif
          enddo
          ! check if the decay chain matches
          if(nmothersA.ne.nmatchA) then
            p_dec(:,1:2) = 0D0
            isvetoed = .true.
          endif
          if(nmothersB.ne.nmatchB) then
            p_dec(:,3:4) = 0D0
            isvetoed = .true.
          endif 
          ! increment the counter of how many decays took place
          ! this is used to calculate the branching ratio
          ! do this before cuts are applied
          do i=1,ndecay
            ! number of decays of individual particles
            if(p_dec(0,i).ne.0D0) then
              cnt_decays(i) = cnt_decays(i) + 1
            endif
          enddo
          ! number of decays with correct final state particles
          if(p_dec(0,1).ne.0D0 .and. p_dec(0,2).ne.0D0 .and.
     &       p_dec(0,3).ne.0D0 .and. p_dec(0,4).ne.0D0) then
            cnt_decays_proc = cnt_decays_proc + 1
          endif
          ! final test: check if we have the correct decay channel
          ! p_dec is only filled if the requested particle was found
          ! in the decay chain
          ! (all requested momenta must be filled)
          ! if the decay ID is not specified in powheg.input, don't
          ! take its momentum into account
          do i=1,ndecay
            if(p_dec(0,i).eq.0D0 .and. id_dec(i).ne.0) then
              p_dec(:,:) = 0D0
              isvetoed = .true. 
            endif
          enddo

#ifdef DEBUGQ
          if(.not.isvetoed) then
          print*,"number of good decays: ",cnt_decays_proc
          print*,"number of total decays: ",cnt_decays_tot
          print*,"BR of total process: ",real(cnt_decays_proc)/real(cnt_decays_tot)
          print*,"Decay chain of A", id_finA
          !print*,id_finA,p_finA
          do i=1,2
          !print*,id_dec(i),p_dec(:,i)
          enddo
          do j=1,nhep-num_decchainA
          if(idhep(ihep_decchainA(j,0)) .ne. 0) then
            print*,idhep(ihep_decchainA(j,0)),"->",
     &       idhep(ihep_decchainA(j,1)),idhep(ihep_decchainA(j,2))
          endif
          enddo

          print*,"Decay chain of B", id_finB
          !print*,id_finB,p_finB
          do i=3,4
          !print*,id_dec(i),p_dec(:,i)
          enddo
          do j=1,nhep-num_decchainB
          if(idhep(ihep_decchainB(j,0)) .ne. 0) then
            print*,idhep(ihep_decchainB(j,0)),"->",
     &       idhep(ihep_decchainB(j,1)),idhep(ihep_decchainB(j,2))
          endif
          enddo
          print*
          endif
#endif

          ! assign observables for the decay products
          do i=1,ndecay
            pt_dec(i) = pt(p_dec(:,i))
            if(pt_dec(i).ne.0D0) then
              call getrapidity(p_dec(:,i),y_dec(i))
              call getpseudorap(p_dec(:,i),eta_dec(i))
              call getazimuth(p_dec(:,i),phi_dec(i))
            endif  
          enddo
          ! compute min R separations for the decay products
          do i=1,ndecay
            RDmin(i) = 1D50 
            do j=1,ndecay
              RDsep(i,j) = 1D50
              if(i.ne.j .and. p_dec(0,i).ne.0D0 .and.
     &           p_dec(0,j).ne.0D0) then
                RDsep(i,j) = rsepn(p_dec(:,i),p_dec(:,j))
                RDmin(i) = min(RDmin(i),RDsep(i,j))
              endif  
            enddo
          enddo
          ! calculate only for particle corI & corJ the difference in phi
          ! -> the user can simply change the correlation indices corI and corJ
          if((pt_dec(corI).ne.0D0) .and. (pt_dec(corJ).ne.0D0)) then
            dphi_decIdecJ = phi_dec(corI)-phi_dec(corJ)
            if(dphi_decIdecJ.lt.(-pi)) then
              dphi_decIdecJ=dphi_decIdecJ+2d0*pi
            endif
            if(dphi_decIdecJ.gt.pi) then
              dphi_decIdecJ=dphi_decIdecJ-2d0*pi
            endif
            do mu=0,3
              p_decIdecJ(mu) = p_dec(mu,corI)+p_dec(mu,corJ)
            enddo
            call getinvmass(p_decIdecJ,m_decIdecJ)
          endif
          
          ! calculate missing Pt if event is not vetoed 
          ! (not the correct final state or whatever)
          if(.not.isvetoed) then
            p_miss(:) = 0D0
            pt_miss   = 0D0
            ! a experimental-like approach: negative momentum sum
            ! over all preselected and detectable particles.
            ! this is missing p since the incoming partons have 
            ! momentum sum 0
            do i=1,ndecay
              if(isdetectable(id_dec(i))) then
                do mu=0,3
                  p_miss(mu) = p_miss(mu) - p_dec(mu,i)
                enddo
              endif
            enddo
            pt_miss1 = pt(p_miss(:)) ! store for debugging
            ! add the detectable jets
            if(mjets.ge.mJmin) then
              do i=1,mjets
                if(kt(i).gt.ptJminPhys .and. eta(i).lt.etaJmax) then
                  p_miss(0) = p_miss(0) - pj(4,i)
                  do mu=1,3
                    p_miss(mu) = p_miss(mu) - pj(mu,i)
                  enddo
                endif
              enddo
            endif
            pt_miss = pt(p_miss(:))
#ifdef DEBUGQ
              print*,pt_miss1,pt_miss
#endif
          endif
        endif
        
        ! fill histograms
        ! no cuts
        call filld("xsec",0D0,dsig)
        call filld("mass(A+B)",m_finAfinB,dsig)
        call filld("Pt(A)",pt_finA,dsig)
        call filld("Pt(B)",pt_finB,dsig)
        call filld("eta(A)",eta_finA,dsig)
        call filld("eta(B)",eta_finB,dsig)
        call filld("y(A)",y_finA,dsig)
        call filld("y(B)",y_finB,dsig)
        call filld("N(J)",dble(num_vis_jets),dsig)
        call filld("Pt(J)",kt(1),dsig)
        call filld("eta(J)",eta(1),dsig)

        ! determine if jet cuts have passed
        JcutsPassed = .true.
        if(mjets.lt.mJmin) JcutsPassed = .false.
        if(kt(1).le.ptJmin) JcutsPassed = .false.
        if(dabs(eta(1)).ge.etaJmax) JcutsPassed = .false.

        ! fill histograms if jet cuts have passed
        if(JcutsPassed) then
          call filld("xsec(Jcuts)",0D0,dsig)
          call filld("N(J,Jcuts)",dble(num_vis_jets),dsig)
          call filld("Pt(J,Jcuts)",kt(1),dsig)
          call filld("eta(J,Jcuts)",eta(1),dsig)
        endif
        
        ! if the branching ratios for the decays are modified in PYTHIA 
        ! by hand correct the weight here:
        ! NOTE: If a modified decay is requested, the
        ! variable bratio should be set in setup-PYTHIA-lhef.f.
        if(bratio.eq.0D0) then
          bratio = 1D0
        endif  
        dsig(:) = dsig(:)*bratio
        
        ! determine if decay cuts have passed
        DcutsPassed = .true.
        if(isvetoed) DcutsPassed = .false.
        if(.not.decayed) DcutsPassed = .false.
        if(pt_miss.le.ptmissmin) DcutsPassed = .false.
        do i=1,ndecay
          ! apply cut only if decay ID is not empty and if particle
          ! is detectable (don't take cuts of neutrinos.. into account)
          if(isdetectable(id_dec(i)) .and. id_dec(i).ne.0) then
            if(pt_dec(i).le.ptDmin) DcutsPassed = .false.
            if(dabs(eta_dec(i)).ge.etaDmax) DcutsPassed = .false.
            if(RDmin(i).le.RDsepmin) DcutsPassed = .false.
          endif  
        enddo
        
        ! histogramms for decay products, no cuts
        if(decayed) then
          call filld("dPhi("//cstr(corI)//"+"//cstr(corJ)//")",
     &               dphi_decIdecJ,dsig)
          call filld("mass("//cstr(corI)//"+"//cstr(corJ)//")",
     &               m_decIdecJ,dsig)
          do i=1,ndecay
            call filld("Pt("//cstr(i)//")",pt_dec(i),dsig)
            call filld("eta("//cstr(i)//")",eta_dec(i),dsig)
          enddo
        endif
        
        ! histogramms for decay products, with Dcuts
        if(DcutsPassed) then
          call filld("dPhi("//cstr(corI)//"+"//cstr(corJ)//",Dcuts)",
     &                 dphi_decIdecJ,dsig)
          call filld("mass("//cstr(corI)//"+"//cstr(corJ)//",Dcuts)",
     &                 m_decIdecJ,dsig)
          do i=1,ndecay
            call filld("Pt("//cstr(i)//",Dcuts)",pt_dec(i),dsig)
            call filld("eta("//cstr(i)//",Dcuts)",eta_dec(i),dsig)
          enddo
          call filld("xsec(Dcuts)",0d0,dsig)
          call filld("Pt(miss)",pt_miss,dsig)
          !print*,"cuts passed"
        endif
        
        !if(.not.isvetoed .and. .not.DcutsPassed) then
        !  print*,"cuts not passed"
        !endif  

        ! print the ratio of the number of requested decays and the total
        ! number of decays (this is the branching ratio).
        ! this could be used to compare with the branching ratio of 
        ! the SLHA file.
        ! do not overload the screen
        if(mod(nevents,20000).eq.1) then
          if(WHCPRG.eq."PYTHIA") then
            print*
            print*,"BR(finA/B->dec(i)+X) = ",
     &             (real(cnt_decays(1:ndecay)))/cnt_decays_tot
            print*,"BR(finA->dec(1)+dec(2)+X) x "//
     &             "BR(finB->dec(3)+dec(4)+X) = ",
     &             (((real(cnt_decays(1)))/cnt_decays_tot)*
     &              ((real(cnt_decays(3)))/cnt_decays_tot))
            print*, "BR(proc,Dcuts) x bratio = ",
     &              (real(cnt_decays_proc))/cnt_decays_tot*bratio
            print*
          endif
        endif
        
        ! increment number of events processed (for printing informations)
        nevents = nevents + 1
      end
c############### end subroutine analysis ###############################

c############### subroutine getihep ####################################
c takes PDG ID and returns the index of this particle in the HEPEVT-list
c call this routine multiple times to get for one particle type
c all occurent indices in the HEPEVT-list
c st is the hepevent status:
c 1: particle undecayed in final state
c 2: particle decayed and not in final state
c use the flag reset to reset the search and to begin at the beginning
      subroutine getihep(PDG_ID,st,reset,HEP_i)
        implicit none
#include "hepevt.h"
        integer PDG_ID, HEP_i, ihep, st
        ! store here the index of the already found particle id's
        ! (to keep track which particle was already found,
        !  -> it is possible to call this subroutine several times with 
        !  the same PDG_ID and get the next HEPEVT-id)
        integer foundlist(nmxhep)
        logical reset
        save foundlist
        if(reset) then
          do ihep=1,nhep
            foundlist = 0
          enddo
        endif
        HEP_i = 0
        do ihep=1,nhep
          if(isthep(ihep).eq.st) then
            if( (foundlist(ihep).ne.idhep(ihep)) .and.
     &          (idhep(ihep).eq.PDG_ID) ) then
              HEP_i = ihep
              foundlist(ihep) = idhep(ihep)
              return
            endif  
          endif
        enddo
      end
c############### end subroutine getihep ################################

c############### subroutine getdaughters ###############################
c takes PDG ID of mother and returns the PDG ID and the momenta of two
c daughters.
c calling this routine for the first time reset must be set to .true.
c calling this routine a second time with reset equals to .false. will
c return the daughters of the next mother.
c example:
c Mother 1 with PDG = 1000023 decays into 1000022 and -14
c Mother 2 with PDG = 1000023 decays into 1000022 and -16
c Mother 3 with PDG = 1000024 decays into 1000022 and -15
c
c call getdaughters(1000023,.true.,D1PDG_ID,D2PDG_ID)
c   -> D1PDG_ID = 1000022 and D2PDG_ID = -14
c call getdaughters(1000023,.false.,D1PDG_ID,D2PDG_ID)
c   -> D1PDG_ID = 1000022 and D2PDG_ID = -16
c call getdaughters(1000024,.false.,D1PDG_ID,D2PDG_ID)
c   -> D1PDG_ID = 1000022 and D2PDG_ID = -15
      subroutine getdaughters(MPDG_ID,reset,D1PDG_ID,D2PDG_ID,pD1,pD2)
        implicit none
#include "hepevt.h"
        integer MPDG_ID, D1PDG_ID, D2PDG_ID, ihep
        double precision pD1(0:3), pD2(0:3)
        logical reset
        call getihep(MPDG_ID,2,reset,ihep)
        D1PDG_ID = idhep(jdahep(1,ihep))
        D2PDG_ID = idhep(jdahep(2,ihep))
        call get4momentum_ihep(jdahep(1,ihep),pD1)
        call get4momentum_ihep(jdahep(2,ihep),pD2)
      end
c############### end subroutine getdaughters ###########################

c############### subroutine isdecayed ##################################
c takes PDG ID and returns if the particle decayed
      subroutine isdecayed(PDG_ID,fdec)
        implicit none
#include "hepevt.h"
        integer PDG_ID, ihep
        logical fdec
        fdec = .false.
        do ihep=1,nhep
          if(isthep(ihep).eq.2) then ! hepevent status 2 -> particle decayed
            if(idhep(ihep).eq.PDG_ID) then
              fdec = .true.
            endif
          endif
        enddo
      end
c############### end subroutine isdecayed ##############################

c############### subroutine get4momentum_sortedpt ######################
c takes PDG ID
c returns four-momentum of undecayed particle with this ID sorted by
c pt. Select with iorder the pt occurence:
c 1: highest pt
c 2: second highest pt
c 3: third highest pt ...
      subroutine get4momentum_sortedpt(PDG_ID,iorder,p)
        implicit none
        integer PDG_ID, n, nmax, mu
        parameter (nmax=200)
        integer list(nmax),iorder
        double precision p(0:3)
        do mu=0,3
          p(mu) = 0d0
        enddo
        call getptsortedlist(PDG_ID,n,list)
        if (list(iorder).ne.0) then  ! only if particle with PDG_ID is found
          ! take four-momentum of first particle in list (largest pt)
          call get4momentum_ihep(list(iorder),p)
        endif
      end
c############### end subroutine get4momentum_sortedpt ##################

c############### subroutine get4momentum_largestpt #####################
c takes PDG ID
c returns four-momentum of undecayed particle with this ID and largest pt
      subroutine get4momentum_largestpt(PDG_ID,p)
        implicit none
        integer PDG_ID
        double precision p(0:3)
        call get4momentum_sortedpt(PDG_ID,1,p)
      end
c############### end subroutine get4momentum_largestpt #################

c############### subroutine getptsortedlist ############################
c takes PDG ID
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
          print*,"more than", nmax, "particles with ID", PDG_ID
          print*,"could not fill pt sorted list"
        else
          ! fill list with found particles
          n = 0
          do ihep=1,nhep
            if(isthep(ihep).eq.1) then
              if(idhep(ihep).eq.PDG_ID) then
                n = n+1
                list(n) = ihep
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
        double precision tmp,pt_m(nmxhep)
        logical touched
        double precision pt
        external pt
        do j=1,n
          pt_m(j) = dsqrt(phep(1,iarr(j))**2+phep(2,iarr(j))**2)
        enddo
        ! bubble sort
        touched = .true.
        do while(touched)
          touched = .false.
          do j=1,n-1
            if(pt_m(j).lt.pt_m(j+1)) then
              k = iarr(j)
              iarr(j) = iarr(j+1)
              iarr(j+1) = k
              tmp = pt_m(j)
              pt_m(j) = pt_m(j+1)
              pt_m(j+1) = tmp
              touched = .true.
            endif
          enddo
        enddo
      end
c############### end subroutine sortbypt ###############################

c############### subroutine getPt ######################################
c calculate Pt, it is assumed p(0) = E, p(1:3) = vec{p}
      double precision function pt(p)
        implicit none
        double precision p(0:3)
        pt = dsqrt(p(1)**2+p(2)**2)
      end
c############### end subroutine getPt ##################################

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
          print*, "Warning: no particle with ihep number", ihep,
     &              "in current event"
          print*, "could not determine four-momentum"
          return
          !stop
        else
          ! p is defined in the laboratory frame (~kn_pborn(i))
          call convert_hep_momenta(phep(:,ihep),p(:))
        endif
      end
c############### end subroutine get4momentum_ihep ######################

c############### subroutine convert_hep_momenta ########################
c converts the hep momentum standard to the conventional momentum
c standard: hep (p_i, E, m) -> conventional (E, p_i)
      subroutine convert_hep_momenta(p_hep,p_con)
        implicit none
        integer mu
        double precision p_hep(1:5)
        double precision p_con(0:3)
        p_con(0) = p_hep(4)
        do mu=1,3
          p_con(mu) = p_hep(mu)
        enddo
      end
c############### end subroutine convert_hep_momenta ####################

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
        print*,"m1",dsqrt(momsq(p1)),dabs(par_Fin1mass)
        print*,"m2",dsqrt(momsq(p2)),dabs(par_Fin2mass)
        !stop
#endif
        ! check if the extracted particles are indeed the final state particles
        relerror = (dsqrt(momsq(p1))-dabs(par_Fin1mass))
     &             /(dsqrt(momsq(p1))+dabs(par_Fin1mass))
        if( dabs(relerror) .gt. 1D-6 ) then
          print*,"error: wrong particle extracted in "
          print*,"subroutine get4momentum_fin_pair"
          print*,"p1**2 = ",momsq(p1)
          print*,"m1**2 = ",par_Fin1mass**2
          stop
        endif
        relerror = (dsqrt(momsq(p2))-dabs(par_Fin2mass))
     &             /(dsqrt(momsq(p2))+dabs(par_Fin2mass))
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

c############### function rsepn ########################################
c calculate the separation in the lego plot between the two momenta
c p1 and p2
      function rsepn(p1,p2)
        implicit none
        include 'pwhg_math.h'  
        double precision rsepn,p1(0:3),p2(0:3)
        double precision y1,phi1,y2,phi2
        double precision delphi
        y1 = 0D0
        y2 = 0D0
        delphi = 0D0
        if(p1(0).eq.0D0 .or. p2(0).eq.0D0) then
          goto 111
        endif

        call getazimuth(p1(:),phi1)
        call getazimuth(p2(:),phi2)
        call getrapidity(p1(:),y1)
        call getrapidity(p2(:),y2)
      
        delphi = dabs(phi1-phi2)
        if (delphi.gt.pi) then
          delphi = 2*pi-delphi
        endif
        if (delphi.lt.0 .or. delphi.gt.pi) then
          print*,' problem in rsepn. delphi = ',delphi
        endif
#ifdef DEBUGQ
        print*,p1
        print*,p2
        print*,y1,y2
        print*,phi1,phi2
        print*,"rsepn",dsqrt((y1-y2)**2 + delphi**2)
#endif
 111    continue
        rsepn = dsqrt((y1-y2)**2 + delphi**2)
      end
c############### function rsepn ########################################

c############### subroutine getazimuth #################################
c calculates azimuth phi in [0,2pi] from x-axis
      subroutine getazimuth(p,phi)
        implicit none
#include "pwhg_math.h"
        double precision p(0:3),phi
        phi=0D0
        if(p(1).gt.0D0) then
          if(p(2).ge.0D0) phi = datan(p(2)/p(1))
          if(p(2).lt.0D0) phi = datan(p(2)/p(1))+2D0*pi
        elseif(p(1).lt.0D0) then
          phi = datan(p(2)/p(1))+pi
        elseif(p(1).EQ.0D0) then
          if(p(2).gt.0D0) phi = pi/2D0
          if(p(2).lt.0D0) phi = -pi/2D0
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

c############### subroutine getpseudorap ###############################
      subroutine getpseudorap(p,eta)
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
c############### end subroutine getpseudorap ###########################      

c############### subroutine buildjets ##################################
c jets are orderd by decreasing pt
      subroutine buildjets(mjets,kt,eta,rap,phi,pjet,ptrel)
        implicit none
#include "hepevt.h"
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
          if(isthep(j).eq.1 .and. isjetparticle(idhep(j))) then
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
        if (ntracks.eq.0) then
          mjets=0
          return
        endif

        ! note: ptJmin and R are in common block "jetcuts"
        ! palg=-1 -> anti-kT (c.f. ZZ code):
        palg=-1
        call fastjetppgenkt(ptrack,ntracks,R,palg,ptJmin,pjet,njets,
     &                      jetvec)
        mjets = min(mjets,njets)
#ifdef DEBUGQ
        print*,"ntracks: ", ntracks
        !print*,"id,ptrack(1): ",idhep(itrackhep(1)),ptrack(4,1),ptrack(1:3,1)
        print*,jetvec(1),ptrack(:,jetvec(1))
        print*,"pjet(1): ", pjet(4,1),pjet(1:3,1)
        !print*,"phep(9): ",phep(1:4,9)
        call nstop(1)
#endif
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
        ! computing arrays of useful kinematics quantities for hardest jets
        do j=1,mjets
          kt(j) = dsqrt(pjet(1,j)**2+pjet(2,j)**2)
          pp = dsqrt(kt(j)**2+pjet(3,j)**2)
          eta(j) = 0.5D0*dlog((pjet(4,j)+pjet(3,j))/
     &                        (pjet(4,j)-pjet(3,j)))
          rap(j) = 0.5D0*dlog((pjet(4,j)+pjet(3,j))/
     &                        (pjet(4,j)-pjet(3,j)))
          phi(j) = datan2(pjet(2,j),pjet(1,j))
        enddo
        call computeptrel(ptrack,ntracks,rap,kt,phi,mjets,jetvec,ptrel)
      end
c############### end subroutine buildjets ##############################

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
#ifdef DEBUGQ
          print*,"is jetparticle: ",PDG_ID
#endif
        else if (injetblacklist(PDG_ID)) then
          isjetparticle = .false.
        else
          print*, "could not determine if particle with ID", PDG_ID
          print*, "contributes to a jet"
          stop
        endif
      end
c############### end subroutine isjetparticle ##########################

c############### function injetwhitelist ###############################
c determines if particle is in whitelist of particles that contribute
c to a jet
      logical function injetwhitelist(PDG_ID)
        implicit none
        integer PDG_ID, id
        injetwhitelist = .false.
        id = iabs(PDG_ID)
        ! the following particles contribute to a jet
        ! quarks
        if(1.le.id .and. id.le.6) injetwhitelist = .true.
        ! photon is not in whitelist
        if(id.eq.22) injetwhitelist = .false.
        ! gluon
        if(id.eq.21) injetwhitelist = .true.
        if(id.eq.0) injetwhitelist = .true.  ! POWHEG numbering
        ! diquarks
        if((1000.le.id .and. id.le.6000)
     &      .and. (mod(id,100).eq.1 .or. mod(id,100).eq.3)) then
          injetwhitelist = .true.
        endif
        ! mesons
        if((111.le.id .and. id.le.600)
     &      .or. (10000.le.id .and. id.le.999999)
     &      .or. (9000000.le.id .and. id.le.9100000)) then
          injetwhitelist = .true.
        endif
        ! baryons
        if((1000.le.id .and. id.le.6000)
     &      .and. (mod(id,10).eq.2 .or. mod(id,10).eq.4)) then
          injetwhitelist = .true.
        endif
      end
c############### end function injetwhitelist ###########################

c############### function injetblacklist ###############################
c determines if particle is in blacklist of particles that do not
c contribute to a jet
      logical function injetblacklist(PDG_ID)
        implicit none
        integer PDG_ID, id
        injetblacklist = .false.
        id = iabs(PDG_ID)
        ! the following particles may occur
        ! but do not contribute to a jet
        ! leptons
        if(id.eq.11) injetblacklist=.true.  ! e
        if(id.eq.12) injetblacklist=.true.  ! nue
        if(id.eq.13) injetblacklist=.true.  ! mu
        if(id.eq.14) injetblacklist=.true.  ! numu
        if(id.eq.15) injetblacklist=.true.  ! tau
        if(id.eq.16) injetblacklist=.true.  ! nu_tau
        ! sleptons
        if(id.eq.1000011) injetblacklist=.true.  ! ~e_L
        if(id.eq.1000012) injetblacklist=.true.  ! ~snue_L
        if(id.eq.1000013) injetblacklist=.true.  ! ~mu_L
        if(id.eq.1000014) injetblacklist=.true.  ! ~snumu_L
        if(id.eq.1000015) injetblacklist=.true.  ! ~stau_1
        if(id.eq.1000016) injetblacklist=.true.  ! ~snu_tau_L
        if(id.eq.2000011) injetblacklist=.true.  ! ~e_R
        if(id.eq.2000013) injetblacklist=.true.  ! ~mu_R
        if(id.eq.2000015) injetblacklist=.true.  ! ~stau_2
        ! neutralinos
        if(id.eq.1000022) injetblacklist=.true.  ! ~n(1)
        if(id.eq.1000023) injetblacklist=.true.  ! ~n(2)
        if(id.eq.1000025) injetblacklist=.true.  ! ~n(3)
        if(id.eq.1000035) injetblacklist=.true.  ! ~n(4)
        ! charginos
        if(id.eq.1000024) injetblacklist=.true.  ! ~x(1)
        if(id.eq.1000037) injetblacklist=.true.  ! ~x(2)
        ! photon
        if(id.eq.22) injetblacklist=.true.
        ! W,Z,Higgs
        if(id.eq.23) injetblacklist=.true.
        if(id.eq.24) injetblacklist=.true.
        if(id.eq.25) injetblacklist=.true.
      end
c############### end function injetblacklist ###########################

c############### subroutine injetblacklist #############################
c determines if particle is detectable or not (for PT_miss)
      logical function isdetectable(PDG_ID)
        implicit none
        integer PDG_ID, id
        isdetectable = .true.
        id = iabs(PDG_ID)
        ! leptons
        if (id.eq.12) isdetectable=.false.  ! nue
        if (id.eq.14) isdetectable=.false.  ! numu
        if (id.eq.16) isdetectable=.false.  ! nu_tau
        ! sleptons
        if (id.eq.1000012) isdetectable=.false.  ! ~snue_L
        if (id.eq.1000014) isdetectable=.false.  ! ~snumu_L
        if (id.eq.1000016) isdetectable=.false.  ! ~snu_tau_L
        ! neutralinos
        if (id.eq.1000022) isdetectable=.false.  ! ~n(1)
        if (id.eq.1000023) isdetectable=.false.  ! ~n(2)
        if (id.eq.1000025) isdetectable=.false.  ! ~n(3)
        if (id.eq.1000035) isdetectable=.false.  ! ~n(4)
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
     &                   /(ptracks(4,i)-ptracks(3,i)))
            ! rapidity of track i in frame where the jet has zero rapidity
            y = y-yj
            ! find momentum of track i in frame where the jet has zero rapidity
            pt(1) = ptracks(1,i)
            pt(2) = ptracks(2,i)
            pt(3) = dsqrt(pt(1)**2+pt(2)**2)*sinh(y)
            ! pt rel is sum of the ptrack momentum projection ortogonal to the jet
            ! momentum in the frame where the jet has zero rapidity
            ptrel(j) = dsqrt(((pt(1)*kj2-pt(2)*kj1)**2+
     &                       (         -pt(3)*kj2)**2+
     &                       (pt(3)*kj1          )**2)/
     &                       (kj1**2+kj2**2)) + ptrel(j)
          endif
        enddo
      end
c############### end subroutine computeptrel ###########################


c############### subroutine store_decchain #############################
      subroutine store_decchain(str_user_decchain,user_decchain)
        ! the considered decay chains
        character*100 str_user_decchain
        integer user_decchain(100,0:2)
        ! operations on strings of decay chains
        character*10 substr
        integer npos,len1
        
        user_decchain(:,:) = 0
        
        len1 = len(trim(str_user_decchain))
        
        do i=1,len(trim(str_user_decchain))
          ! get the mother (string in front of "->")
          npos = index(str_user_decchain,"->")
          substr = str_user_decchain(1:npos-1)
          ! convert string to int
          read(substr(:),*) user_decchain(i,0)
          str_user_decchain = str_user_decchain(npos+2:len1)
          ! get the first daughter
          npos = index(str_user_decchain,",")
          substr = str_user_decchain(1:npos-1)
          read(substr(:),*) user_decchain(i,1)
          str_user_decchain = str_user_decchain(npos+1:len1)
          ! get the second daughter
          npos = index(str_user_decchain,";")
          substr = str_user_decchain(1:npos-1)
          read(substr(:),*) user_decchain(i,2)
          str_user_decchain = str_user_decchain(npos+1:len1)
          print*,user_decchain(i,0),"->",user_decchain(i,1),user_decchain(i,2)
          if(len(trim(str_user_decchain)).eq.0) then
            exit
          endif
        enddo
      end
c############### end subroutine store_decchain #########################
