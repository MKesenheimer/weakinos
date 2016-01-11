c variables needed by madgraph to read from slha file
      double precision  alpha,gfermi,alphas,yt,yb,yl
      common/values/    alpha,gfermi,alphas,yt,yb,yl
      
c neutralino and chargino mixing matrices (used to check if all parameters
c from madgraph are equal to the ones in formcalc)
c stop and sbottom mixing matrices, mixing angle tan(beta)
      double precision bwmix(4,4),uumix(2,2),vvmix(2,2),mix_t(2,2),mix_b(2,2)
      double precision tanb
      double precision twid
      common/printmad/ bwmix,uumix,vvmix,mix_t,mix_b,tanb,twid