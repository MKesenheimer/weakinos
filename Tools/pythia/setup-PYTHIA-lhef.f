c     note on SUSY extension:
c     to make Pythia read the SLHA spectrum, we rely on this assumption:
c     the routines in this file are called exactly once by
c     main_pythia in ../main-PYTHIA.f in a specific order:
c     call setup_PYTHIA_tune
c     call PYINIT
c     call setup_PYTHIA_parameters


c     print some Pythia configuration
      SUBROUTINE PYTHIA_PRINT_SOME_PARAMS
      IMPLICIT real *8(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
      WRITE(*,*) "***************************************************"
      WRITE(*,*) " PYTHIA SM parameters"
      WRITE(*,*) "***************************************************"
      WRITE(*,*) "1/ALPHEM(0)   = ", (1d0/PYALEM(0d0))
      WRITE(*,*) "1/ALPHEM(MZ2) = ", (1d0/PYALEM(PMAS(23,1)**2))
      WRITE(*,*) "MW = ", PMAS(24,1)
      WRITE(*,*) "MZ = ", PMAS(23,1)
      WRITE(*,*) "S2W = ", PARU(102)
      WRITE(*,*) "gfermi = ",PARU(105)
      WRITE(*,*) "***************************************************"
      END



c     tells PYTHIA to use SLHA data and where to find them
      subroutine setup_PYTHIA_tune
      implicit none
      ! PYTHIA MSSM and subprocess common blocks
      real *8 RMSS,CKIN
      integer IMSS,MSEL,MSELPD,MSUB,KFIN
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      integer LUNSPC, tuneid
      COMMON/PYSLHAIO/LUNSPC  ! SLHA spectrum file descriptor
      character*100 slhafilename
      call powheginputstring("SLHA",slhafilename)
      
      ! 100       A : Rick Field's CDF Tune A                     (Oct 2002)
      ! 103      DW : Rick Field's CDF Tune DW                    (Apr 2006)
      ! 320 Perugia 0 : "Perugia" update of S0-Pro                (Feb 2009)

c      tuneid = 320
      tuneid = -1 ! no tune setting
      if (tuneid.ge.0) then
        write(*,*) "setting Pythia tune ", tuneid
        call PYTUNE(tuneid)
      else
        write(*,*) "no Pythia tune is set"
      endif

      ! set SUSY parameters using SLHA input
      ! switch on SUSY MSSM input from an SLHA file
      IMSS(1)=11
      ! open the SLHA file and tell Pythia on which LUN to find it:
      LUNSPC=22
      OPEN(LUNSPC,file=slhafilename,status='old')
      IMSS(21)=LUNSPC
      ! If the file also contains an SLHA decay table that you want to use,
      ! uncomment the next line:
      IMSS(22)=LUNSPC

      end



c     sets several options in PYTHIA
c     hands over SM parameters from init_couplings to PYTHIA
      subroutine setup_PYTHIA_parameters
      implicit none
      include 'sm_read_values.inc'
      include 'coupl.inc'
      include 'hepevt.h'
      include 'LesHouches.h'
      real *8 parp,pari
      integer mstp,msti
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
      integer MSTU,MSTJ
      real *8 PARU,PARJ
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer KCHG
      real *8 PMAS, PARF, VCKM
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      integer MDCY,MDME,KFDP
      real *8 brat
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer pycomp
      external pycomp
      ! multiple interactions (mostly untested)
      logical mult_inter
      parameter (mult_inter=.false.)
      integer maxev
      common/mcmaxev/maxev
      ! SLHA file
      integer LUNSPC
      COMMON/PYSLHAIO/LUNSPC  ! SLHA spectrum file descriptor
      logical update_pythia_parameters
      parameter (update_pythia_parameters=.true.)
      ! further MADGRAPH common block
      double precision  sw, sw2, cw, tw, s2w, c2w, e, gx, gz,
     &                  qe, qu, qd, t3e, t3v, t3u, t3d, rt2, ppi,
     &                  r_lc(2,2), r_bc(2,2), r_tc(2,2)
      common /ewparam/  sw, sw2, cw, tw, s2w, c2w, e, gx, gz,
     &                  qe, qu, qd, t3e, t3v, t3u, t3d, rt2, ppi,
     &                  r_lc, r_bc, r_tc

      ! close Pythia specific handle to SLHA input file
      close(LUNSPC)
      ! update also SM parameters consistently
      if (update_pythia_parameters) then
        call init_couplings
        MSTU(101) = 0         ! running alpha_em off
        PARU(101) = alpha     ! fixed value for alpha_em
        PARU(102) = sw2
        PARU(105) = gfermi
        ! special parameters (might not be used at all)
        PARJ(123)=PMAS(23,1)  ! set Z mass for e+e- routines
        PARJ(124)=PMAS(23,2)  ! set Z witdh for e+e- routines
        call PYTHIA_PRINT_SOME_PARAMS
      endif

      ! multiple interactions
      ! (MI can increase a lot the execution time)
      if(.not.mult_inter) then
         write(*,*) "deselecting MPI in Pythia"
         mstp(81)=20   ! No Multiple interactions. Force a call to PYEVNW 
      else
         write(*,*) "selecting MPI in Pythia"
         mstp(81)=21   ! MPI on in the PYEVNW MPI scenario
      endif

      write(*,*) "other Pythia settings:"

c      mstj(41)=12 ! photon radiation off quarks and leptons
      mstj(41)=11 ! no photon radiation off quarks and leptons
      write(*,*) "control of photon radiation mstj(41): ", mstj(41)

c      mstp(61)=0                ! No IS shower
      mstp(61)=1                ! IS shower just QCD in hadronic events
      write(*,*) "IS shower control mstp(61):           ", mstp(61)

c      mstp(71)=0                ! No FS shower
      write(*,*) "FS shower control mstp(71):           ", mstp(71)

c      mstp(91)=0                ! No Primordial kt
      write(*,*) "primordial kt control mstp(91):       ", mstp(91)

      mstp(131)=0               ! No Pile Up
      write(*,*) "pile-up control mstp(131):            ", mstp(131)

c      mstp(111)=0               ! No hadronization
      write(*,*) "hadronisation control mstp(111):      ", mstp(111)

      ! decays without hadronization need a patched PYTHIA code
      mstp(41)=1               ! force all resonance decays
c      mstp(41)=0               ! prevent all resonance decays
      write(*,*) "resonance decays:                     ", mstp(41)
      
#ifdef DEBUGQ
        print*, "[DEBUG] Stop"
        stop
#endif

c      mstp(64) =3   ! use Lambda_MC for IS shower > 6.4.19
c      mstp(64) =1   ! use Lambda_MSbar (default)

      ! number of warnings printed on the shell
      mstu(26)=20

c      call PYLIST(12)  ! to see the PYTHIA decay table

      ! tolerate 2% of killed events
      mstu(22)=maxev/50

      end



c     opens input file and counts number of events, setting maxev
      subroutine getmaxev(maxev)
      implicit none
      integer maxev
      call opencount(maxev)
      end



      subroutine UPINIT
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      real *8 parp,pari
      integer mstp,msti
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
      integer MSTU,MSTJ
      real *8 PARU,PARJ
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MDCY,MDME,KFDP
      real *8 brat
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer pycomp
      external pycomp
      integer maxev
      common/mcmaxev/maxev

      nevhep=0
      ! read the header first, so lprup is set
      call lhefreadhdr(97)
      ! Make PI0 stable as in herwig default
      mdcy(pycomp(111),1)=0
      if (lprup(1).eq.10015)  mdcy(pycomp(15),1)=0

      end



      subroutine UPEVNT
      implicit none
      call lhefreadev(97)
      end



c     PYTHIA routine to abort event
      subroutine upveto
      implicit none
      end



      subroutine pyabeg
      implicit none
      call init_hist
      end



c     writes PYTHIA analysis output into .top file
      subroutine pyaend
      implicit none
      character *20 pwgprefix
      character *100 filename
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      filename=pwgprefix(1:lprefix)//'POWHEG+PYTHIA-output'
      call pwhgsetout
      call pwhgtopout(filename)
      close(99)
      end



c     jumps to next event and calls analysis
      subroutine pyanal
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      integer mint
      real *8 vint
      COMMON/PYINT1/MINT(400),VINT(400)
      ! check parameters
      logical verbose
      parameter (verbose=.false.)

      if(mint(51).ne.0) then
         if(verbose) then
            write(*,*) 'Killed event'
            write(*,*) 'Scalup= ',scalup
            call pylist(7)      !hepeup
            call pylist(2)      !all the event
         endif
         return
      endif
      nevhep=nevhep+1
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup
      end
