      integer    maxparticles
      parameter (maxparticles=35)
      integer    maxclasses
      parameter (maxclasses=15)
      integer eventinfo(0:maxclasses,maxparticles)
      common/eventinfo/eventinfo
      integer    max1
      parameter (max1=35)        !max1=maxparticles
      integer    max2
      parameter (max2=595)       !max2=35*34/2
      integer    max3
      parameter (max3=6545)      !max3=35*34*33/6
      integer    max4
      parameter (max4=52360)     !max4=35*34*33*32/24
      integer s1(max1),s2(2,max2),s3(3,max3),s4(4,max4)
      common/pointers/s1,s2,s3,s4
      double precision etmin(max1),etmax(max1)
      double precision ptmin(Max1),ptmax(Max1)
      double precision emin(max1),emax(max1)
      double precision mommin(Max1),mommax(Max1)
      double precision etamin(Max1),etamax(Max1)
      double precision ymin(Max1),ymax(Max1)
      double precision costhmin(Max1),costhmax(Max1)
      double precision phimin(Max1),phimax(Max1)
      double precision X1min(Max1),X1max(Max1)
      double precision X2min(Max1),X2max(Max1)
      double precision X3min(Max1),X3max(Max1)
      double precision dRijmin(Max2),dRijmax(Max2)
      double precision mijmin(Max2),mijmax(Max2)
      double precision delta_phimin(Max2),delta_phimax(Max2)
      double precision delta_etamin(Max2),delta_etamax(Max2)
      double precision kTijmin(Max2),kTijmax(Max2)
      double precision XY1min(Max2),XY1max(Max2)
      double precision XY2min(Max2),XY2max(Max2)
      double precision XY3min(Max2),XY3max(Max2)
      double precision XYZ1min(Max3),XYZ1max(Max3)
      double precision XYZ2min(Max3),XYZ2max(Max3)
      double precision XYZ3min(Max3),XYZ3max(Max3)
      double precision XYZA1min(Max4),XYZA1max(Max4)
      double precision XYZA2min(Max4),XYZA2max(Max4)
      double precision XYZA3min(Max4),XYZA3max(Max4)
      double precision etmissmin,etmissmax
      common/to_cuts/  etmin,etmax,etmissmin,etmissmax,etamax,etamin
     &                ,dRijmin,dRijmax,mijmin,mijmax,ptmin,ptmax,delta_phimin,delta_phimax
     &                ,emin,emax,mommin,mommax,ymax,ymin
     &                ,kTijmin,kTijmax,costhmin,costhmax,phimin,phimax,delta_etamin,delta_etamax
     &                ,X1min,X1max,X2min,X2max,X3min,X3max
     &                ,XY1min,XY1max,XY2min,XY2max,XY3min,XY3max
     &                ,XYZ1min,XYZ1max,XYZ2min,XYZ2max,XYZ3min,XYZ3max
     &                ,XYZA1min,XYZA1max,XYZA2min,XYZA2max,XYZA3min,XYZA3max
      character*8 name(0:maxclasses,maxparticles)
      integer lname(0:maxclasses,maxparticles)
      common/names/name,lname
      logical etmin_active(max1),etmax_active(max1)
      logical ptmin_active(max1),ptmax_active(max1)
      logical emin_active(max1),emax_active(max1)
      logical mommin_active(max1),mommax_active(max1)
      logical etamin_active(max1),etamax_active(max1)
      logical ymin_active(max1),ymax_active(max1)
      logical costhmin_active(max1),costhmax_active(max1)
      logical phimin_active(max1),phimax_active(max1)
      logical X1min_active(max1),X1max_active(max1)
      logical X2min_active(max1),X2max_active(max1)
      logical X3min_active(max1),X3max_active(max1)
      logical dRijmin_active(max2),dRijmax_active(max2)
      logical mijmin_active(max2),mijmax_active(max2)
      logical delta_phimin_active(max2),delta_phimax_active(max2)
      logical delta_etamin_active(max2),delta_etamax_active(max2)
      logical kTijmin_active(max2),kTijmax_active(max2)
      logical XY1min_active(max2),XY1max_active(max2)
      logical XY2min_active(max2),XY2max_active(max2)
      logical XY3min_active(max2),XY3max_active(max2)
      logical XYZ1min_active(max3),XYZ1max_active(max3)
      logical XYZ2min_active(max3),XYZ2max_active(max3)
      logical XYZ3min_active(max3),XYZ3max_active(max3)
      logical XYZA1min_active(max4),XYZA1max_active(max4)
      logical XYZA2min_active(max4),XYZA2max_active(max4)
      logical XYZA3min_active(max4),XYZA3max_active(max4)
      common/active/   etmin_active,etmax_active,ptmin_active,ptmax_active,etamin_active,etamax_active
     &                ,X1min_active,X1max_active,X2min_active,X2max_active,X3min_active,X3max_active
     &                ,emin_active,emax_active,mommin_active,mommax_active,ymin_active,ymax_active
     &                ,costhmin_active,costhmax_active,phimin_active,phimax_active,kTijmin_active
     &                ,dRijmin_active,dRijmax_active,mijmin_active,mijmax_active,kTijmax_active
     &                ,delta_phimin_active,delta_phimax_active,delta_etamin_active,delta_etamax_active
     &                ,XY1min_active,XY1max_active,XY2min_active,XY2max_active,XY3min_active,XY3max_active
     &                ,XYZ1min_active,XYZ1max_active,XYZ2min_active,XYZ2max_active,XYZ3min_active,XYZ3max_active
     &                ,XYZA1min_active,XYZA1max_active,XYZA2min_active,XYZA2max_active,XYZA3min_active,XYZA3max_active
      integer ptplot(maxclasses),etaplot(maxclasses),etplot(maxclasses)
      integer mijplot(maxclasses),dRijplot(maxclasses),delta_phiplot(maxclasses)
      integer eplot(maxclasses),yplot(maxclasses),momplot(maxclasses),kTijplot(maxclasses)
      integer costhplot(maxclasses),phiplot(maxclasses),delta_etaplot(maxclasses)
      integer X1plot(maxclasses),X2plot(maxclasses),X3plot(maxclasses)
      integer XY1plot(maxclasses),XY2plot(maxclasses),XY3plot(maxclasses)
      integer XYZ1plot(maxclasses),XYZ2plot(maxclasses),XYZ3plot(maxclasses)
      integer XYZA1plot(maxclasses),XYZA2plot(maxclasses),XYZA3plot(maxclasses)
      logical gnu
      common/plotoptions/gnu,X1plot,X2plot,X3plot,XY1plot,XY2plot,XY3plot,
     &                   XYZ1plot,XYZ2plot,XYZ3plot,XYZA1plot,XYZA2plot,XYZA3plot,
     &                   ptplot,etplot,etaplot,mijplot,dRijplot,delta_phiplot,
     &                   eplot,momplot,yplot,costhplot,phiplot,delta_etaplot,kTijplot
      double precision ptrange(3), etarange(3), etrange(3), mijrange(3),delta_etarange(3)
      double precision erange(3), yrange(3), momrange(3), costhrange(3), phirange(3)
      double precision dRijrange(3), delta_phirange(3), X1range(3), X2range(3),kTijrange(3)
      double precision X3range(3), XY1range(3), XY2range(3), XY3range(3)
      double precision XYZ1range(3), XYZ2range(3), XYZ3range(3), etmissrange(3)
      double precision XYZA1range(3), XYZA2range(3), XYZA3range(3)
      common/ranges/ptrange,etarange,etrange,mijrange,dRijrange,delta_phirange,
     &              X1range,X2range,X3range,XY1range,XY2range,XY3range,etmissrange,
     &              XYZ1range,XYZ2range,XYZ3range,XYZA1range,XYZA2range,XYZA3range,
     &              erange,momrange,yrange,costhrange,phirange,delta_etarange,kTijrange

