c############### model_sm.h ############################################
* model_sm.h
* common blocks for the model parameters
* this file is part of FormCalc
* last modified 23 Dec 14 th by MK

        !Z mass, W mass, cos(theta_w), sin(theta_w)
	RealType MZ, MZ2, MW, MW2, CW, CW2, SW, SW2
	RealType WZ, WW
        ! Fermi constant, fine structure, fine structure constant at MZ, strong coupling
	RealType GF, Alfa, Alfa2, AlfaMZ, AlfasMZ, Alfas, Alfas2
        ! Quark masses
	RealType MU, MU2, MC, MC2, MT, MT2
	RealType MD, MD2, MS, MS2, MB, MB2
	RealType Mf(4,3),Mf2(4,3)
        ! CKM matrix (here: identity)
	RealType CKM(3,3)
        ! strong and e.m. coupling constants
	RealType EL, GS

        common /smpara/ MZ, MZ2, MW, MW2, CW, CW2, SW, SW2
        common /smpara/ WZ, WW
        common /smpara/ GF, Alfa, Alfa2, AlfaMZ, AlfasMZ, Alfas, Alfas2
        common /smpara/ MU, MU2, MC, MC2, MT, MT2
        common /smpara/ MD, MD2, MS, MS2, MB, MB2
        common /smpara/ Mf,Mf2
	common /smpara/ CKM, EL, GS

c############### end model_sm.h ########################################
