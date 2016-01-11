c#####################Print Model Parameters############################
c prints all model paramters from madgraph and formcalc

c############### subroutine printFormcalcModelParameters ###############

      subroutine printFormcalcModelParameters

#include "types.h"
#include "PhysPars.h"

        print*,"== SM Parameters == "
        print*
        print*,"CKM matrix"
        print*,"CKM(1,1) = ",CKM(1,1)
        print*,"CKM(1,2) = ",CKM(1,2)
        print*,"CKM(1,3) = ",CKM(1,3)
        print*,"CKM(2,1) = ",CKM(2,1)
        print*,"CKM(2,2) = ",CKM(2,2)
        print*,"CKM(2,3) = ",CKM(2,3)
        print*,"CKM(3,1) = ",CKM(3,1)
        print*,"CKM(3,2) = ",CKM(3,2)
        print*,"CKM(3,3) = ",CKM(3,3)
        print*
        print*,"Strong and e.m. couplings"
        print*,"Alfa^-1    = ", 1/Alfa
        print*,"AlfaS^-1   = ", 1/AlfaS
        print*,"AlfaMZ^-1  = ", 1/AlfaMZ
        print*,"AlfaSMZ^-1 = ", 1/AlfaSMZ
        print*
        print*,"EL = ", EL
        print*,"GS = ", GS
        print*
        print*,"Mixing angles"
        print*,"CW  = ", CW
        print*,"CW2 = ", CW2
        print*,"SW  = ", SW
        print*,"SW2 = ", SW2
        print*,"CB  = ", CB
        print*,"CB2 = ", CB2
        print*
        print*,"Quark masses"
        print*,"MU = ", MU
        print*,"MC = ", MC
        print*,"MT = ", MT
        print*,"MD = ", MD
        print*,"MS = ", MS
        print*,"MB = ", MB
        print*
        print*,"Z/W mass"
        print*,"MZ = ", MZ
        print*,"MW = ", MW
        print*
        print*

        print*,"== MSSM Parameters == "
        print*
        print*,"Neutralino mixing matrix:"
        print*,"ZNeu(1,1) = ", ZNeu(1,1)
        print*,"ZNeu(1,2) = ", ZNeu(1,2)
        print*,"ZNeu(1,3) = ", ZNeu(1,3)
        print*,"ZNeu(1,4) = ", ZNeu(1,4)
        print*,"ZNeu(2,1) = ", ZNeu(2,1)
        print*,"ZNeu(2,2) = ", ZNeu(2,2)
        print*,"ZNeu(2,3) = ", ZNeu(2,3)
        print*,"ZNeu(2,4) = ", ZNeu(2,4)
        print*,"ZNeu(3,1) = ", ZNeu(3,1)
        print*,"ZNeu(3,2) = ", ZNeu(3,2)
        print*,"ZNeu(3,3) = ", ZNeu(3,3)
        print*,"ZNeu(3,4) = ", ZNeu(3,4)
        print*,"ZNeu(4,1) = ", ZNeu(4,1)
        print*,"ZNeu(4,2) = ", ZNeu(4,2)
        print*,"ZNeu(4,3) = ", ZNeu(4,3)
        print*,"ZNeu(4,4) = ", ZNeu(4,4)
        print*
        print*,"Chargino mixing matrix:"
        print*,"UCha(1,1) = ", UCha(1,1)
        print*,"UCha(1,2) = ", UCha(1,2)
        print*,"UCha(2,1) = ", UCha(2,1)
        print*,"UCha(2,2) = ", UCha(2,2)
        print*
        print*,"VCha(1,1) = ", VCha(1,1)
        print*,"VCha(1,2) = ", VCha(1,2)
        print*,"VCha(2,1) = ", VCha(2,1)
        print*,"VCha(2,2) = ", VCha(2,2)
        print*
        print*,"Neutralino masses"
        print*,"MNeu(1) = ",MNeu(1)
        print*,"MNeu(2) = ",MNeu(2)
        print*,"MNeu(3) = ",MNeu(3)
        print*,"MNeu(4) = ",MNeu(4)
        print*
        print*,"Sfermion masses"
        print*,"MSf(1,2,1) = ",MSf(1,2,1) !el
        print*,"MSf(2,2,1) = ",MSf(2,2,1) !er
        print*,"MSf(1,3,1) = ",MSf(1,3,1) !ul
        print*,"MSf(2,3,1) = ",MSf(2,3,1) !ur
        print*,"MSf(1,4,1) = ",MSf(1,4,1) !dl
        print*,"MSf(2,4,1) = ",MSf(2,4,1) !dr
        print*,"MSf(1,3,2) = ",MSf(1,3,2) !cl
        print*,"MSf(2,3,2) = ",MSf(2,3,2) !cr
        print*,"MSf(1,4,2) = ",MSf(1,4,2) !sl
        print*,"MSf(2,4,2) = ",MSf(2,4,2) !sr
        print*,"MSf(1,3,3) = ",MSf(1,3,3) !tl
        print*,"MSf(2,3,3) = ",MSf(2,3,3) !tr
        print*,"MSf(1,4,3) = ",MSf(1,4,3) !bl
        print*,"MSf(2,4,3) = ",MSf(2,4,3) !br
        print*
        print*,"Squark widths"
        print*,"WSf(1,3,1) = ",WSf(1,3,1)
        print*,"WSf(2,3,1) = ",WSf(2,3,1)
        print*,"WSf(1,4,1) = ",WSf(1,4,1)
        print*,"WSf(2,4,1) = ",WSf(2,4,1)
        print*,"WSf(1,3,2) = ",WSf(1,3,2)
        print*,"WSf(2,3,2) = ",WSf(2,3,2)
        print*,"WSf(1,4,2) = ",WSf(1,4,2)
        print*,"WSf(2,4,2) = ",WSf(2,4,2)
        print*,"WSf(1,3,3) = ",WSf(1,3,3)
        print*,"WSf(2,3,3) = ",WSf(2,3,3)
        print*,"WSf(1,4,3) = ",WSf(1,4,3)
        print*,"WSf(2,4,3) = ",WSf(2,4,3)
        print*
        print*,"Gluino Mass"
        print*,"MGl = ", MGl
        print*
        print*,"Sfermion mixing matrix USf(s,s',t,g):"
        print*,"su:"
        print*,"USf(1,1,3,1) = ",USf(1,1,3,1)
        print*,"USf(1,2,3,1) = ",USf(1,2,3,1)
        print*,"USf(2,1,3,1) = ",USf(2,1,3,1)
        print*,"USf(2,2,3,1) = ",USf(2,2,3,1)
        print*,"sd:"
        print*,"USf(1,1,4,1) = ",USf(1,1,4,1)
        print*,"USf(1,2,4,1) = ",USf(1,2,4,1)
        print*,"USf(2,1,4,1) = ",USf(2,1,4,1)
        print*,"USf(2,2,4,1) = ",USf(2,2,4,1)
        print*,"sc:"
        print*,"USf(1,1,3,2) = ",USf(1,1,3,2)
        print*,"USf(1,2,3,2) = ",USf(1,2,3,2)
        print*,"USf(2,1,3,2) = ",USf(2,1,3,2)
        print*,"USf(2,2,3,2) = ",USf(2,2,3,2)
        print*,"ss:"
        print*,"USf(1,1,4,2) = ",USf(1,1,4,2)
        print*,"USf(1,2,4,2) = ",USf(1,2,4,2)
        print*,"USf(2,1,4,2) = ",USf(2,1,4,2)
        print*,"USf(2,2,4,2) = ",USf(2,2,4,2)
        print*,"st:"
        print*,"USf(1,1,3,3) = ",USf(1,1,3,3)
        print*,"USf(1,2,3,3) = ",USf(1,2,3,3)
        print*,"USf(2,1,3,3) = ",USf(2,1,3,3)
        print*,"USf(2,2,3,3) = ",USf(2,2,3,3)
        print*,"sb:"
        print*,"USf(1,1,4,3) = ",USf(1,1,4,3)
        print*,"USf(1,2,4,3) = ",USf(1,2,4,3)
        print*,"USf(2,1,4,3) = ",USf(2,1,4,3)
        print*,"USf(2,2,4,3) = ",USf(2,2,4,3)
        print*
        print*
!         print*,"MSSM mixing angles:" ! Note: the matrix elements do not depent on beta
!         print*,"sin(beta) [SB] = ",SB
!         print*,"cos(beta) [CB] = ",CB
!         print*,"tan(beta) [TB] = ",TB
!         print*
!         print*
        
      end

c############### end subroutine printFormcalcModelParameters ###########

c############### subroutine printMadgraphModelParameters ###############

      subroutine printMadgraphModelParameters

#include "types.h"
#include "PhysPars.h"

        print*,"== SM Parameters == "
        print*
        print*,"CKM matrix"
        print*,"Note: intrinsic diagonal CKM matrix is used"
        print*,"      by MadGraph for model MSSM"
        print*
        print*,"Strong and e.m. couplings"
        print*,"Alpha^-1    = ", 1/Alpha
        print*,"AlphaS^-1   = ", 1/AlphaS
        print*
        print*,"EL = ", dreal(gal(1))
        print*,"GS = ", dreal(-gg(1))
        print*
        print*,"Mixing angles"
        print*,"CW  = ", gwwz/gw
        print*,"CW2 = ", (gwwz/gw)**2
        print*,"SW  = ", gwwa/gw
        print*,"SW2 = ", (gwwa/gw)**2
        print*
        print*,"Quark masses"
        print*,"MU = ", 0d0 ! Note: umass, dmass and smass do not occur in madgraph
        print*,"MC = ", cmass
        print*,"MT = ", tmass
        print*,"MD = ", 0d0
        print*,"MS = ", 0d0
        print*,"MB = ", bmass
        print*
        print*,"Z/W mass"
        print*,"MZ = ", zmass
        print*,"MW = ", wmass
        print*
        print*

        print*,"== MSSM Parameters == "
        print*
        print*,"Neutralino mixing matrix:"
        print*,"ZNeu(1,1) = ", bwmix(1,1)
        print*,"ZNeu(1,2) = ", bwmix(1,2)
        print*,"ZNeu(1,3) = ", bwmix(1,3)
        print*,"ZNeu(1,4) = ", bwmix(1,4)
        print*,"ZNeu(2,1) = ", bwmix(2,1)
        print*,"ZNeu(2,2) = ", bwmix(2,2)
        print*,"ZNeu(2,3) = ", bwmix(2,3)
        print*,"ZNeu(2,4) = ", bwmix(2,4)
        print*,"ZNeu(3,1) = ", bwmix(3,1)
        print*,"ZNeu(3,2) = ", bwmix(3,2)
        print*,"ZNeu(3,3) = ", bwmix(3,3)
        print*,"ZNeu(3,4) = ", bwmix(3,4)
        print*,"ZNeu(4,1) = ", bwmix(4,1)
        print*,"ZNeu(4,2) = ", bwmix(4,2)
        print*,"ZNeu(4,3) = ", bwmix(4,3)
        print*,"ZNeu(4,4) = ", bwmix(4,4)
        print*
        print*,"Chargino mixing matrix:"
        print*,"UCha(1,1) = ", uumix(1,1)
        print*,"UCha(1,2) = ", uumix(1,2)
        print*,"UCha(2,1) = ", uumix(2,1)
        print*,"UCha(2,2) = ", uumix(2,2)
        print*
        print*,"VCha(1,1) = ", vvmix(1,1)
        print*,"VCha(1,2) = ", vvmix(1,2)
        print*,"VCha(2,1) = ", vvmix(2,1)
        print*,"VCha(2,2) = ", vvmix(2,2)
        print*
        print*,"Neutralino masses"
        print*,"MNeu(1) = ",mn1
        print*,"MNeu(2) = ",mn2
        print*,"MNeu(3) = ",mn3
        print*,"MNeu(4) = ",mn4
        print*
        print*,"Sfermion masses"
        print*,"MSf(1,3,1) = ",mul
        print*,"MSf(2,3,1) = ",mur
        print*,"MSf(1,4,1) = ",mdl
        print*,"MSf(2,4,1) = ",mdr
        print*,"MSf(1,3,2) = ",mcl
        print*,"MSf(2,3,2) = ",mcr
        print*,"MSf(1,4,2) = ",msl
        print*,"MSf(2,4,2) = ",msr
        print*,"MSf(1,3,3) = ",mtl
        print*,"MSf(2,3,3) = ",mtr
        print*,"MSf(1,4,3) = ",mbl
        print*,"MSf(2,4,3) = ",mbr
        print*
        print*,"Squark widths"
        print*,"WSf(1,3,1) = ",wul
        print*,"WSf(2,3,1) = ",wur
        print*,"WSf(1,4,1) = ",wdl
        print*,"WSf(2,4,1) = ",wdr
        print*,"WSf(1,3,2) = ",wcl
        print*,"WSf(2,3,2) = ",wcr
        print*,"WSf(1,4,2) = ",wsl
        print*,"WSf(2,4,2) = ",wsr
        print*,"WSf(1,3,3) = ",wt1
        print*,"WSf(2,3,3) = ",wt2
        print*,"WSf(1,4,3) = ",wb1
        print*,"WSf(2,4,3) = ",wb2
        print*
        print*,"Gluino Mass"
        print*,"MGl = ", mgo
        print*
        print*,"Sfermion mixing matrix USf(s,s',t,g):"
        print*,"Note: Sfermion matrices for u, d, c and s are not"
        print*,"      present in MadGraph"
!         print*,"su:"
!         print*,"USf(1,1,3,1) = ",1d0
!         print*,"USf(1,2,3,1) = ",0d0
!         print*,"USf(2,1,3,1) = ",0d0
!         print*,"USf(2,2,3,1) = ",1d0
!         print*,"sd:"
!         print*,"USf(1,1,4,1) = ",0d0
!         print*,"USf(1,2,4,1) = ",1d0
!         print*,"USf(2,1,4,1) = ",1d0
!         print*,"USf(2,2,4,1) = ",0d0
!         print*,"sc:"
!         print*,"USf(1,1,3,2) = ",1d0
!         print*,"USf(1,2,3,2) = ",0d0
!         print*,"USf(2,1,3,2) = ",0d0
!         print*,"USf(2,2,3,2) = ",1d0
!         print*,"ss:"
!         print*,"USf(1,1,4,2) = ",0d0
!         print*,"USf(1,2,4,2) = ",1d0
!         print*,"USf(2,1,4,2) = ",1d0
!         print*,"USf(2,2,4,2) = ",0d0
        print*,"st:"
        print*,"USf(1,1,3,3) = ",mix_t(1,1)
        print*,"USf(1,2,3,3) = ",mix_t(1,2)
        print*,"USf(2,1,3,3) = ",mix_t(2,1)
        print*,"USf(2,2,3,3) = ",mix_t(2,2)
        print*,"sb:"
        print*,"USf(1,1,4,3) = ",mix_b(1,1)
        print*,"USf(1,2,4,3) = ",mix_b(1,2)
        print*,"USf(2,1,4,3) = ",mix_b(2,1)
        print*,"USf(2,2,4,3) = ",mix_b(2,2)
        print*
        print*
!         print*,"MSSM mixing angles:" ! Note: the matrix elements do not depent on beta
!         print*,"sin(beta) [SB] = ",abs(tanb/sqrt(1+tanb**2))
!         print*,"cos(beta) [CB] = ",abs(1.D0/sqrt(1+tanb**2))
!         print*,"tan(beta) [TB] = ",tanb
        
      end

c############### end subroutine printMadgraphModelParameters ###########
