c This file has been generated automatically by FeynRules
c Version: 1.5.8   Date: 7. 7. 2009,    13:57

c External Parameters :
      double precision aEWM1,Gf,aS,ZM,MF
      common/external_params/ aEWM1,Gf,aS,ZM,MF

c Internal Parameters :
      double precision aEW,EE,cM,sM,RM,x,t,EpL,g0,gp,gt,fpi,EpRnu,E
     +pRel,EpRup,EpRdo,EpRmu,EpRta,EpRch,EpRst,EpRto,EpRbo,MGZ,MGW
     +,MGZP,MGWP,MhNu,v0g,v1g,v2g,v0Z,v1Z,v2Z,v0ZP,v1ZP,v2ZP,v0W,v
     +1W,v0WP,v1WP,Npi,v0pi,v1pi,v0piP,v1piP,Npi0,v0pi0,v1pi0,v0pi
     +0P,v1pi0P,vL0nu,vL0el,vL0mu,vL0ta,vL0up,vL0do,vL0ch,vL0st,vL
     +0to,vL0bo,vL1nu,vL1el,vL1mu,vL1ta,vL1up,vL1do,vL1ch,vL1st,vL
     +1to,vL1bo,vR1nu,vR1el,vR1up,vR1do,vR1mu,vR1ta,vR1ch,vR1st,vR
     +1to,vR1bo,vR2nu,vR2el,vR2up,vR2do,vR2mu,vR2ta,vR2ch,vR2st,vR
     +2to,vR2bo,vL0hNu,vL0hEl,vL0hMu,vL0hTa,vL0hUp,vL0hDo,vL0hCh,v
     +L0hSt,vL0hTo,vL0hBo,vL1hNu,vL1hEl,vL1hMu,vL1hTa,vL1hUp,vL1hD
     +o,vL1hCh,vL1hSt,vL1hTo,vL1hBo,vR1hNu,vR1hEl,vR1hUp,vR1hDo,vR
     +1hMu,vR1hTa,vR1hCh,vR1hSt,vR1hTo,vR1hBo,vR2hNu,vR2hEl,vR2hUp
     +,vR2hDo,vR2hMu,vR2hTa,vR2hCh,vR2hSt,vR2hTo,vR2hBo,Sqrt2,Sqrt
     +Pi,Pv0W2,Pv1W2,PEE2,Pg02,Pv0W4,Pgt2,Pv1W4,Pv0WP2,Pv1WP2,Pv0W
     +3,Pv1W3,Pv0WP3,Pv1WP3,Pv0WP4,Pv1WP4,Pv0Z2,Pv1Z2,Pv0ZP2,Pv1ZP
     +2,PvL0do2,PvL1do2,PvL0st2,PvL1st2,PvR1st2,PvR2st2,PvL0bo2,Pv
     +L1bo2,PvR1bo2,PvR2bo2,PvL0hDo2,PvL1hDo2,PvL0hSt2,PvL1hSt2,Pv
     +R1hSt2,PvR2hSt2,PvL0hBo2,PvL1hBo2,PvR1hBo2,PvR2hBo2,PvL0el2,
     +PvL1el2,PvL0mu2,PvL1mu2,PvR1mu2,PvR2mu2,PvL0ta2,PvL1ta2,PvR1
     +ta2,PvR2ta2,PvL0hEl2,PvL1hEl2,PvL0hMu2,PvL1hMu2,PvR1hMu2,PvR
     +2hMu2,PvL0hTa2,PvL1hTa2,PvR1hTa2,PvR2hTa2,PvL0hUp2,PvL1hUp2,
     +PvL0hCh2,PvL1hCh2,PvR1hCh2,PvR2hCh2,PvL0hTo2,PvL1hTo2,PvR1hT
     +o2,PvR2hTo2,PvL0up2,PvL1up2,PvL0ch2,PvL1ch2,PvR1ch2,PvR2ch2,
     +PvL0to2,PvL1to2,PvR1to2,PvR2to2,PvL0hNu2,PvL1hNu2,PvL0nu2,Pv
     +L1nu2,GGT1
      common/internal_params_R/ aEW,EE,cM,sM,RM,x,t,EpL,g0,gp,gt,fp
     +i,EpRnu,EpRel,EpRup,EpRdo,EpRmu,EpRta,EpRch,EpRst,EpRto,EpRb
     +o,MGZ,MGW,MGZP,MGWP,MhNu,v0g,v1g,v2g,v0Z,v1Z,v2Z,v0ZP,v1ZP,v
     +2ZP,v0W,v1W,v0WP,v1WP,Npi,v0pi,v1pi,v0piP,v1piP,Npi0,v0pi0,v
     +1pi0,v0pi0P,v1pi0P,vL0nu,vL0el,vL0mu,vL0ta,vL0up,vL0do,vL0ch
     +,vL0st,vL0to,vL0bo,vL1nu,vL1el,vL1mu,vL1ta,vL1up,vL1do,vL1ch
     +,vL1st,vL1to,vL1bo,vR1nu,vR1el,vR1up,vR1do,vR1mu,vR1ta,vR1ch
     +,vR1st,vR1to,vR1bo,vR2nu,vR2el,vR2up,vR2do,vR2mu,vR2ta,vR2ch
     +,vR2st,vR2to,vR2bo,vL0hNu,vL0hEl,vL0hMu,vL0hTa,vL0hUp,vL0hDo
     +,vL0hCh,vL0hSt,vL0hTo,vL0hBo,vL1hNu,vL1hEl,vL1hMu,vL1hTa,vL1
     +hUp,vL1hDo,vL1hCh,vL1hSt,vL1hTo,vL1hBo,vR1hNu,vR1hEl,vR1hUp,
     +vR1hDo,vR1hMu,vR1hTa,vR1hCh,vR1hSt,vR1hTo,vR1hBo,vR2hNu,vR2h
     +El,vR2hUp,vR2hDo,vR2hMu,vR2hTa,vR2hCh,vR2hSt,vR2hTo,vR2hBo,S
     +qrt2,SqrtPi,Pv0W2,Pv1W2,PEE2,Pg02,Pv0W4,Pgt2,Pv1W4,Pv0WP2,Pv
     +1WP2,Pv0W3,Pv1W3,Pv0WP3,Pv1WP3,Pv0WP4,Pv1WP4,Pv0Z2,Pv1Z2,Pv0
     +ZP2,Pv1ZP2,PvL0do2,PvL1do2,PvL0st2,PvL1st2,PvR1st2,PvR2st2,P
     +vL0bo2,PvL1bo2,PvR1bo2,PvR2bo2,PvL0hDo2,PvL1hDo2,PvL0hSt2,Pv
     +L1hSt2,PvR1hSt2,PvR2hSt2,PvL0hBo2,PvL1hBo2,PvR1hBo2,PvR2hBo2
     +,PvL0el2,PvL1el2,PvL0mu2,PvL1mu2,PvR1mu2,PvR2mu2,PvL0ta2,PvL
     +1ta2,PvR1ta2,PvR2ta2,PvL0hEl2,PvL1hEl2,PvL0hMu2,PvL1hMu2,PvR
     +1hMu2,PvR2hMu2,PvL0hTa2,PvL1hTa2,PvR1hTa2,PvR2hTa2,PvL0hUp2,
     +PvL1hUp2,PvL0hCh2,PvL1hCh2,PvR1hCh2,PvR2hCh2,PvL0hTo2,PvL1hT
     +o2,PvR1hTo2,PvR2hTo2,PvL0up2,PvL1up2,PvL0ch2,PvL1ch2,PvR1ch2
     +,PvR2ch2,PvL0to2,PvL1to2,PvR1to2,PvR2to2,PvL0hNu2,PvL1hNu2,P
     +vL0nu2,PvL1nu2,GGT1
