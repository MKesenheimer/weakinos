c############### susy_restore.h ########################################
c last modified by MK, 19 May 2015

#ifndef SUSY_RESTORE_H
#define SUSY_RESTORE_H

        ! Added SUSY-restoring term for Yukawa coupling
        ! quark-squark-gluino (not needed here):
        !     gsy -> gsy + dZgs1y
        !                  dZgs1y = gs*Alfas/(3*Pi)
        ! and SUSY-restoring term for Yukawa coupling
        ! quark-squark-neutralino:
        !     gy -> gy + dZe1y
        !                dZe1y = - e*Alfas/(6*Pi)
        RealType gsy,gsy2,dZgs1y,ely,ely2,dZe1y
        common/SUSYrestore/gsy,gsy2,dZgs1y,ely,ely2,dZe1y
        
#endif

c############### end susy_restore.h ####################################