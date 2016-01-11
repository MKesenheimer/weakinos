c############### set_general_couplings.f ##################################
c last modified by MK, date

c############### subroutine set_general_couplings #########################
      subroutine set_general_couplings(fin1,fin2)
        implicit none
        
#include "finalstate.h"

        include 'coupl.inc'
        integer fin1,fin2

        final1 = fin1
        final2 = fin2
        
        if(fin1.eq.1000022) then
          MNI = mn1
          
          GDLNIP = GDLN1P
          GDRNIP = GDRN1P
          GDLNIM = GDLN1M
          GDRNIM = GDRN1M
          
          GULNIP = GULN1P
          GURNIP = GURN1P
          GULNIM = GULN1M
          GURNIM = GURN1M
          
          GB1NIP = GB1N1P
          GB2NIP = GB2N1P
          GB1NIM = GB1N1M
          GB2NIM = GB2N1M
          
          GT1NIP = GT1N1P
          GT2NIP = GT2N1P
          GT1NIM = GT1N1M
          GT2NIM = GT2N1M
          
        else if(fin1.eq.1000023) then
          MNI = mn2
          
          GDLNIP = GDLN2P
          GDRNIP = GDRN2P
          GDLNIM = GDLN2M
          GDRNIM = GDRN2M
          
          GULNIP = GULN2P
          GURNIP = GURN2P
          GULNIM = GULN2M
          GURNIM = GURN2M
          
          GB1NIP = GB1N2P
          GB2NIP = GB2N2P
          GB1NIM = GB1N2M
          GB2NIM = GB2N2M
          
          GT1NIP = GT1N2P
          GT2NIP = GT2N2P
          GT1NIM = GT1N2M
          GT2NIM = GT2N2M
          
        else if(fin1.eq.1000025) then
          MNI = mn3
          
          GDLNIP = GDLN3P
          GDRNIP = GDRN3P
          GDLNIM = GDLN3M
          GDRNIM = GDRN3M
          
          GULNIP = GULN3P
          GURNIP = GURN3P
          GULNIM = GULN3M
          GURNIM = GURN3M
          
          GB1NIP = GB1N3P
          GB2NIP = GB2N3P
          GB1NIM = GB1N3M
          GB2NIM = GB2N3M
          
          GT1NIP = GT1N3P
          GT2NIP = GT2N3P
          GT1NIM = GT1N3M
          GT2NIM = GT2N3M
          
        else if(fin1.eq.1000035) then
          MNI = mn4
          
          GDLNIP = GDLN4P
          GDRNIP = GDRN4P
          GDLNIM = GDLN4M
          GDRNIM = GDRN4M
          
          GULNIP = GULN4P
          GURNIP = GURN4P
          GULNIM = GULN4M
          GURNIM = GURN4M
          
          GB1NIP = GB1N4P
          GB2NIP = GB2N4P
          GB1NIM = GB1N4M
          GB2NIM = GB2N4M
          
          GT1NIP = GT1N4P
          GT2NIP = GT2N4P
          GT1NIM = GT1N4M
          GT2NIM = GT2N4M
          
        else if(abs(fin1).eq.1000024) then
          MXI = mx1
          
          GDLXIP = GDLX1P
          GDLXIM = GDLX1M
          
          GULXIP = GULX1P
          GULXIM = GULX1M
          
          GB1XIP = GB1X1P
          GB2XIP = GB2X1P
          GB1XIM = GB1X1M
          GB2XIM = GB2X1M
          
          GT1XIP = GT1X1P
          GT2XIP = GT2X1P
          GT1XIM = GT1X1M
          GT2XIM = GT2X1M
          
        else if(abs(fin1).eq.1000037) then
          MXI = mx2
          
          GDLXIP = GDLX2P
          GDLXIM = GDLX2M
          
          GULXIP = GULX2P
          GULXIM = GULX2M
          
          GB1XIP = GB1X2P
          GB2XIP = GB2X2P
          GB1XIM = GB1X2M
          GB2XIM = GB2X2M
          
          GT1XIP = GT1X2P
          GT2XIP = GT2X2P
          GT1XIM = GT1X2M
          GT2XIM = GT2X2M
          
        else
          print*,"error in set_general_couplings."
          print*,"fin1", fin1
          stop
        endif
        
        if(fin2.eq.1000022) then
          MNJ = mn1
          
          GDLNJP = GDLN1P
          GDRNJP = GDRN1P
          GDLNJM = GDLN1M
          GDRNJM = GDRN1M
          
          GULNJP = GULN1P
          GURNJP = GURN1P
          GULNJM = GULN1M
          GURNJM = GURN1M
          
          GB1NJP = GB1N1P
          GB2NJP = GB2N1P
          GB1NJM = GB1N1M
          GB2NJM = GB2N1M
          
          GT1NJP = GT1N1P
          GT2NJP = GT2N1P
          GT1NJM = GT1N1M
          GT2NJM = GT2N1M
          
        else if(fin2.eq.1000023) then
          MNJ = mn2
          
          GDLNJP = GDLN2P
          GDRNJP = GDRN2P
          GDLNJM = GDLN2M
          GDRNJM = GDRN2M
          
          GULNJP = GULN2P
          GURNJP = GURN2P
          GULNJM = GULN2M
          GURNJM = GURN2M
          
          GB1NJP = GB1N2P
          GB2NJP = GB2N2P
          GB1NJM = GB1N2M
          GB2NJM = GB2N2M
          
          GT1NJP = GT1N2P
          GT2NJP = GT2N2P
          GT1NJM = GT1N2M
          GT2NJM = GT2N2M
          
        else if(fin2.eq.1000025) then
          MNJ = mn3
          
          GDLNJP = GDLN3P
          GDRNJP = GDRN3P
          GDLNJM = GDLN3M
          GDRNJM = GDRN3M
          
          GULNJP = GULN3P
          GURNJP = GURN3P
          GULNJM = GULN3M
          GURNJM = GURN3M
          
          GB1NJP = GB1N3P
          GB2NJP = GB2N3P
          GB1NJM = GB1N3M
          GB2NJM = GB2N3M
          
          GT1NJP = GT1N3P
          GT2NJP = GT2N3P
          GT1NJM = GT1N3M
          GT2NJM = GT2N3M
          
        else if(fin2.eq.1000035) then
          MNJ = mn4
          
          GDLNJP = GDLN4P
          GDRNJP = GDRN4P
          GDLNJM = GDLN4M
          GDRNJM = GDRN4M
          
          GULNJP = GULN4P
          GURNJP = GURN4P
          GULNJM = GULN4M
          GURNJM = GURN4M
          
          GB1NJP = GB1N4P
          GB2NJP = GB2N4P
          GB1NJM = GB1N4M
          GB2NJM = GB2N4M
          
          GT1NJP = GT1N4P
          GT2NJP = GT2N4P
          GT1NJM = GT1N4M
          GT2NJM = GT2N4M
          
        else if(abs(fin2).eq.1000024) then
          MXJ = mx1
          
          GDLXJP = GDLX1P
          GDLXJM = GDLX1M
          
          GULXJP = GULX1P
          GULXJM = GULX1M
          
          GB1XJP = GB1X1P
          GB2XJP = GB2X1P
          GB1XJM = GB1X1M
          GB2XJM = GB2X1M
          
          GT1XJP = GT1X1P
          GT2XJP = GT2X1P
          GT1XJM = GT1X1M
          GT2XJM = GT2X1M
          
        else if(abs(fin2).eq.1000037) then
          MXJ = mx2
          
          GDLXJP = GDLX2P
          GDLXJM = GDLX2M
          
          GULXJP = GULX2P
          GULXJM = GULX2M
          
          GB1XJP = GB1X2P
          GB2XJP = GB2X2P
          GB1XJM = GB1X2M
          GB2XJM = GB2X2M
          
          GT1XJP = GT1X2P
          GT2XJP = GT2X2P
          GT1XJM = GT1X2M
          GT2XJM = GT2X2M
          
        else
          print*,"error in set_general_couplings."
          print*,"fin2", fin2
          stop
        endif
        
        if((fin1.eq.1000022.and.fin2.eq.1000022) .or.
     &     (fin1.eq.1000022.and.fin2.eq.1000022) ) then
          GH1NIJ = GH1N11
          GH2NIJ = GH2N11
          GH3NIJ = GH3N11
          GZNIJ  = GZN11
        else if((fin1.eq.1000022.and.fin2.eq.1000023) .or.
     &          (fin1.eq.1000023.and.fin2.eq.1000022) ) then
          GH1NIJ = GH1N12
          GH2NIJ = GH2N12
          GH3NIJ = GH3N12
          GZNIJ  = GZN12
        else if((fin1.eq.1000022.and.fin2.eq.1000025) .or.
     &          (fin1.eq.1000025.and.fin2.eq.1000022) ) then
          GH1NIJ = GH1N13
          GH2NIJ = GH2N13
          GH3NIJ = GH3N13
          GZNIJ  = GZN13
        else if((fin1.eq.1000022.and.fin2.eq.1000035) .or.
     &          (fin1.eq.1000035.and.fin2.eq.1000022)) then
          GH1NIJ = GH1N14
          GH2NIJ = GH2N14
          GH3NIJ = GH3N14
          GZNIJ  = GZN14
        
        else if((fin1.eq.1000023.and.fin2.eq.1000023) .or.
     &          (fin1.eq.1000023.and.fin2.eq.1000023)) then
          GH1NIJ = GH1N22
          GH2NIJ = GH2N22
          GH3NIJ = GH3N22
          GZNIJ  = GZN22
        else if((fin1.eq.1000023.and.fin2.eq.1000025) .or.
     &          (fin1.eq.1000025.and.fin2.eq.1000023)) then
          GH1NIJ = GH1N23
          GH2NIJ = GH2N23
          GH3NIJ = GH3N23
          GZNIJ  = GZN23
        else if((fin1.eq.1000023.and.fin2.eq.1000035) .or.
     &          (fin1.eq.1000035.and.fin2.eq.1000023)) then
          GH1NIJ = GH1N24
          GH2NIJ = GH2N24
          GH3NIJ = GH3N24
          GZNIJ  = GZN24
          
        else if((fin1.eq.1000025.and.fin2.eq.1000025) .or.
     &          (fin1.eq.1000025.and.fin2.eq.1000025)) then
          GH1NIJ = GH1N33
          GH2NIJ = GH2N33
          GH3NIJ = GH3N33
          GZNIJ  = GZN33
        else if((fin1.eq.1000025.and.fin2.eq.1000035) .or.
     &          (fin1.eq.1000035.and.fin2.eq.1000025)) then
          GH1NIJ = GH1N34
          GH2NIJ = GH2N34
          GH3NIJ = GH3N34
          GZNIJ  = GZN34

        else if((fin1.eq.1000035.and.fin2.eq.1000035) .or.
     &          (fin1.eq.1000035.and.fin2.eq.1000035)) then
          GH1NIJ = GH1N44
          GH2NIJ = GH2N44
          GH3NIJ = GH3N44
          GZNIJ  = GZN44
          
        else if((abs(fin1).eq.1000024.and.abs(fin2).eq.1000024) .or.
     &          (abs(fin1).eq.1000024.and.abs(fin2).eq.1000024)) then
          GH1XIJ = GH1X11
          GH2XIJ = GH2X11
          GH3XIJ = GH3X11
          GZXIJ  = GZX11
          GH1XII = GH1X11
          GH2XII = GH2X11
          GH3XII = GH3X11
          GZXII  = GZX11
        else if((abs(fin1).eq.1000024.and.abs(fin2).eq.1000037) .or.
     &          (abs(fin1).eq.1000037.and.abs(fin2).eq.1000024)) then
          GH1XIJ = GH1X12
          GH2XIJ = GH2X12
          GH3XIJ = GH3X12
          GZXIJ  = GZX12
        else if((abs(fin1).eq.1000037.and.abs(fin2).eq.1000037) .or.
     &          (abs(fin1).eq.1000037.and.abs(fin2).eq.1000037)) then
          GH1XIJ = GH1X22
          GH2XIJ = GH2X22
          GH3XIJ = GH3X22
          GZXIJ  = GZX22
          GH1XII = GH1X22
          GH2XII = GH2X22
          GH3XII = GH3X22
          GZXII  = GZX22
          
        else if((fin1.eq.1000022.and.abs(fin2).eq.1000024) .or.
     &          (abs(fin1).eq.1000024.and.fin2.eq.1000022)) then
          GWNIXJ = GWN1X1
          GWXJNI = GWX1N1
          GHNIXJ = GHN1X1
          GHXJNI = GHX1N1
        else if((fin1.eq.1000022.and.abs(fin2).eq.1000037) .or.
     &          (abs(fin1).eq.1000037.and.fin2.eq.1000022)) then
          GWNIXJ = GWN1X2
          GWXJNI = GWX2N1
          GHNIXJ = GHN1X2
          GHXJNI = GHX2N1
     
        else if((fin1.eq.1000023.and.abs(fin2).eq.1000024) .or.
     &          (abs(fin1).eq.1000024.and.fin2.eq.1000023)) then
          GWNIXJ = GWN2X1
          GWXJNI = GWX1N2
          GHNIXJ = GHN1X1
          GHXJNI = GHX1N2
        else if((fin1.eq.1000023.and.abs(fin2).eq.1000037) .or.
     &          (abs(fin1).eq.1000037.and.fin2.eq.1000023)) then
          GWNIXJ = GWN2X2
          GWXJNI = GWX2N2
          GHNIXJ = GHN2X2
          GHXJNI = GHX2N2
          
        else if((fin1.eq.1000025.and.abs(fin2).eq.1000024) .or.
     &          (abs(fin1).eq.1000024.and.fin2.eq.1000025)) then
          GWNIXJ = GWN3X1
          GWXJNI = GWX1N3
          GHNIXJ = GHN3X1
          GHXJNI = GHX1N3
        else if((fin1.eq.1000025.and.abs(fin2).eq.1000037) .or.
     &          (abs(fin1).eq.1000037.and.fin2.eq.1000025)) then
          GWNIXJ = GWN3X2
          GWXJNI = GWX2N3
          GHNIXJ = GHN3X2
          GHXJNI = GHX2N3
          
        else if((fin1.eq.1000035.and.abs(fin2).eq.1000024) .or.
     &          (abs(fin1).eq.1000024.and.fin2.eq.1000035)) then
          GWNIXJ = GWN4X1
          GWXJNI = GWX1N4
          GHNIXJ = GHN4X1
          GHXJNI = GHX1N4
        else if((fin1.eq.1000035.and.abs(fin2).eq.1000037) .or.
     &          (abs(fin1).eq.1000037.and.fin2.eq.1000035)) then
          GWNIXJ = GWN4X2
          GWXJNI = GWX2N4
          GHNIXJ = GHN4X2
          GHXJNI = GHX1N4
          
        else
          print*,"error in set_general_couplings."
          print*,"fin1", fin1
          print*,"fin2", fin2
          stop
        endif
        
      end

c############### end subroutine set_general_couplings #####################
