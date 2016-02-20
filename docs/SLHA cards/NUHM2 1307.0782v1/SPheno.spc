# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.3.8  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 20.02.2016,  11:51
Block SPINFO         # Program information
     1   SPheno      # spectrum calculator
     2   v3.3.8      # version number
#
Block SPhenoINFO     # SPheno specific information
    1      2         # using 2-loop RGEs
    2      1         # using running masses for boundary conditions at mZ
Block MODSEL  # Model selection
    1    1    # mSUGRA model
Block MINPAR  # Input parameters
    1    1.00000000E+04  # m0      
    2    5.00000000E+02  # m12     
    3    7.24717518E+00  # tanb at m_Z   
    4    1.00000000E+00  # cos(phase_mu)
    5   -1.60000000E+04  # A0
#
Block gauge Q=  1.00000000E+04  # (GUT scale)
   1    4.74187592E-01  # g'(M_GUT)^DRbar
   2    6.35907945E-01  # g(M_GUT)^DRbar
   3    9.96371026E-01  # g3(M_GUT)^DRbar
Block Yu Q=  1.00000000E+04  # (GUT scale)
  1  1     6.50884004E-06   # Y_u(M_GUT)^DRbar
  2  2     3.30649155E-03   # Y_c(M_GUT)^DRbar
  3  3     8.04356324E-01   # Y_t(M_GUT)^DRbar
Block Yd Q=  1.00000000E+04  # (GUT scale)
  1  1     8.89673850E-05   # Y_d(M_GUT)^DRbar
  2  2     1.69037991E-03   # Y_s(M_GUT)^DRbar
  3  3     8.37724160E-02   # Y_b(M_GUT)^DRbar
Block Ye Q=  1.00000000E+04  # (GUT scale)
  1  1     2.02283110E-05   # Y_e(M_GUT)^DRbar
  2  2     4.18257329E-03   # Y_mu(M_GUT)^DRbar
  3  3     7.03377492E-02   # Y_tau(M_GUT)^DRbar
Block EXTPAR  # non-universal input parameters
   23    6.00000000E+03  # mu
   25    7.00000000E+00  # tan(beta)
   26    9.50000000E+02  # m_A, pole mass
Block SMINPUTS  # SM parameters
         1     1.27931205E+02  # alpha_em^-1(MZ)^MSbar
         2     1.16637900E-05  # G_mu [GeV^-2]
         3     1.18400000E-01  # alpha_s(MZ)^MSbar
         4     9.11876000E+01  # m_Z(pole)
         5     4.00000000E+00  # m_b(m_b), MSbar
         6     1.73300000E+02  # m_t(pole)
         7     1.77682000E+00  # m_tau(pole)
         8     0.00000000E+00  # m_nu_3
        11     5.10998930E-04  # m_e(pole)
        12     0.00000000E+00  # m_nu_1
        13     1.05658372E-01  # m_muon(pole)
        14     0.00000000E+00  # m_nu_2
        21     5.00000000E-03  # m_d(2 GeV), MSbar
        22     2.50000000E-03  # m_u(2 GeV), MSbar
        23     9.50000000E-02  # m_s(2 GeV), MSbar
        24     1.27000000E+00  # m_c(m_c), MSbar
Block gauge Q=  1.00000000E+03  # (SUSY scale)
   1    3.59471792E-01  # g'(Q)^DRbar
   2    6.31464134E-01  # g(Q)^DRbar
   3    1.04102081E+00  # g3(Q)^DRbar
Block Yu Q=  1.00000000E+03  # (SUSY scale)
  1  1     6.98416204E-06   # Y_u(Q)^DRbar
  2  2     3.54795321E-03   # Y_c(Q)^DRbar
  3  3     8.38246236E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.00000000E+03  # (SUSY scale)
  1  1     9.82256395E-05   # Y_d(Q)^DRbar
  2  2     1.86628613E-03   # Y_s(Q)^DRbar
  3  3     9.15726388E-02   # Y_b(Q)^DRbar
Block Ye Q=  1.00000000E+03  # (SUSY scale)
  1  1     2.06934464E-05   # Y_e(Q)^DRbar
  2  2     4.27874518E-03   # Y_mu(Q)^DRbar
  3  3     7.19391095E-02   # Y_tau(Q)^DRbar
Block Au Q=  1.00000000E+03  # (SUSY scale)
  1  1    -1.52009687E+04   # A_u(Q)^DRbar
  2  2    -1.52009515E+04   # A_c(Q)^DRbar
  3  3    -1.43320558E+04   # A_t(Q)^DRbar
Block Ad Q=  1.00000000E+03  # (SUSY scale)
  1  1    -1.60898566E+04   # A_d(Q)^DRbar
  2  2    -1.60898469E+04   # A_s(Q)^DRbar
  3  3    -1.57936942E+04   # A_b(Q)^DRbar
Block Ae Q=  1.00000000E+03  # (SUSY scale)
  1  1    -1.60094534E+04   # A_e(Q)^DRbar
  2  2    -1.60094282E+04   # A_mu(Q)^DRbar
  3  3    -1.60023344E+04   # A_tau(Q)^DRbar
Block MSOFT Q=  1.00000000E+03  # soft SUSY breaking masses at Q
   1    4.80582506E+02  # M_1
   2    4.96269115E+02  # M_2
   3    5.52155868E+02  # M_3
  21   -3.76495351E+07  # M^2_(H,d)
  22   -5.65366404E+07  # M^2_(H,u)
  31    9.99820617E+03  # M_(L,11)
  32    9.99819521E+03  # M_(L,22)
  33    9.99510965E+03  # M_(L,33)
  34    9.99887653E+03  # M_(E,11)
  35    9.99885450E+03  # M_(E,22)
  36    9.99264935E+03  # M_(E,33)
  41    9.96872802E+03  # M_(Q,11)
  42    9.96871958E+03  # M_(Q,22)
  43    9.61728914E+03  # M_(Q,33)
  44    9.97120573E+03  # M_(U,11)
  45    9.97119263E+03  # M_(U,22)
  46    9.24925766E+03  # M_(U,33)
  47    9.97066114E+03  # M_(D,11)
  48    9.97065725E+03  # M_(D,22)
  49    9.96147594E+03  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73300000E+02  # m_t(pole)
        23     9.11876000E+01  # m_Z(pole)
        24     8.04061001E+01  # W+
        15     1.77682000E+00  # m_tau(pole)
        25     3.20741621E+01  # h0
        35     9.34202300E+02  # H0
        36     9.50000000E+02  # A0
        37     9.32816130E+02  # H+
   1000001     1.01006434E+04  # ~d_L
   2000001     1.00780665E+04  # ~d_R
   1000002     1.01003052E+04  # ~u_L
   2000002     1.00797865E+04  # ~u_R
   1000003     1.01006418E+04  # ~s_L
   2000003     1.00780653E+04  # ~s_R
   1000004     1.01003133E+04  # ~c_L
   2000004     1.00797721E+04  # ~c_R
   1000005     9.95408257E+03  # ~b_1
   2000005     1.00764988E+04  # ~b_2
   1000006     9.74011321E+03  # ~t_1
   2000006     1.00139300E+04  # ~t_2
   1000011     1.00179474E+04  # ~e_L-
   2000011     1.00091734E+04  # ~e_R-
   1000012     1.00172137E+04  # ~nu_eL
   1000013     1.00179550E+04  # ~mu_L-
   2000013     1.00091568E+04  # ~mu_R-
   1000014     1.00172106E+04  # ~nu_muL
   1000015     1.00052571E+04  # ~tau_1-
   2000015     1.00193322E+04  # ~tau_2-
   1000016     1.00163315E+04  # ~nu_tauL
   1000021     7.49439380E+02  # ~g
   1000022     4.96681760E+02  # ~chi_10
   1000023     5.50327765E+02  # ~chi_20
   1000025    -6.03287678E+03  # ~chi_30
   1000035     6.03327094E+03  # ~chi_40
   1000024     5.50331913E+02  # ~chi_1+
   1000037     6.03344860E+03  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -2.03176871E-01   # alpha
Block Hmix Q=  1.00000000E+03  # Higgs mixing parameters
   1    6.00000000E+03  # mu
   2    7.00000000E+00  # tan[beta](Q)
   3    2.46107279E+02  # v(Q)
   4    3.96575161E+06  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     4.79479449E-01   # Re[R_st(1,1)]
   1  2     8.77553108E-01   # Re[R_st(1,2)]
   2  1    -8.77553108E-01   # Re[R_st(2,1)]
   2  2     4.79479449E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     9.98937795E-01   # Re[R_sb(1,1)]
   1  2     4.60790729E-02   # Re[R_sb(1,2)]
   2  1    -4.60790729E-02   # Re[R_sb(2,1)]
   2  2     9.98937795E-01   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     3.97863105E-01   # Re[R_sta(1,1)]
   1  2     9.17444794E-01   # Re[R_sta(1,2)]
   2  1    -9.17444794E-01   # Re[R_sta(2,1)]
   2  2     3.97863105E-01   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1     9.99961515E-01   # Re[N(1,1)]
   1  2    -4.21219540E-03   # Re[N(1,2)]
   1  3     7.49918914E-03   # Re[N(1,3)]
   1  4    -1.72873256E-03   # Re[N(1,4)]
   2  1     4.31323722E-03   # Re[N(2,1)]
   2  2     9.99904095E-01   # Re[N(2,2)]
   2  3    -1.28106081E-02   # Re[N(2,3)]
   2  4     3.01412688E-03   # Re[N(2,4)]
   3  1    -4.05076694E-03   # Re[N(3,1)]
   3  2     6.94449877E-03   # Re[N(3,2)]
   3  3     7.07047682E-01   # Re[N(3,3)]
   3  4     7.07120174E-01   # Re[N(3,4)]
   4  1    -6.47736182E-03   # Re[N(4,1)]
   4  2     1.12174865E-02   # Re[N(4,2)]
   4  3     7.07010061E-01   # Re[N(4,3)]
   4  4    -7.07084850E-01   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -9.99835279E-01   # Re[U(1,1)]
   1  2     1.81497813E-02   # Re[U(1,2)]
   2  1     1.81497813E-02   # Re[U(2,1)]
   2  2     9.99835279E-01   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -9.99990856E-01   # Re[V(1,1)]
   1  2     4.27634141E-03   # Re[V(1,2)]
   2  1     4.27634141E-03   # Re[V(2,1)]
   2  2     9.99990856E-01   # Re[V(2,2)]
Block SPhenoLowEnergy  # low energy observables
    1    3.45375823E-04   # BR(b -> s gamma)
    2    1.59030866E-06   # BR(b -> s mu+ mu-)
    3    3.53314221E-05   # BR(b -> s nu nu)
    4    2.06331432E-15   # BR(Bd -> e+ e-)
    5    8.81424405E-11   # BR(Bd -> mu+ mu-)
    6    1.84493259E-08   # BR(Bd -> tau+ tau-)
    7    6.95759096E-14   # BR(Bs -> e+ e-)
    8    2.97228048E-09   # BR(Bs -> mu+ mu-)
    9    6.30360626E-07   # BR(Bs -> tau+ tau-)
   10    9.65814734E-05   # BR(B_u -> tau nu)
   11    9.97799619E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
   12    5.41822560E-01   # |Delta(M_Bd)| [ps^-1] 
   13    1.93595147E+01   # |Delta(M_Bs)| [ps^-1] 
   16    2.15925651E-03   # epsilon_K
   17    2.28167954E-15   # Delta(M_K)
   18    2.48546099E-11   # BR(K^0 -> pi^0 nu nu)
   19    8.30345937E-11   # BR(K^+ -> pi^+ nu nu)
   20    1.84754322E-17   # Delta(g-2)_electron/2
   21    7.89883549E-13   # Delta(g-2)_muon/2
   22    2.23468428E-10   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39   -2.23195517E-04   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.88854100E-01   # C7
     0305 4422   00   2    -2.14215578E-01   # C7
     0305 4322   00   2    -5.20652058E-04   # C7'
     0305 6421   00   0    -9.52389567E-02   # C8
     0305 6421   00   2    -1.24888321E-01   # C8
     0305 6321   00   2    -6.08948399E-04   # C8'
 03051111 4133   00   0     1.17430609E+00   # C9 e+e-
 03051111 4133   00   2     1.17418399E+00   # C9 e+e-
 03051111 4233   00   2    -1.17106599E-06   # C9' e+e-
 03051111 4137   00   0    -3.99876551E+00   # C10 e+e-
 03051111 4137   00   2    -4.00214849E+00   # C10 e+e-
 03051111 4237   00   2     3.15382412E-05   # C10' e+e-
 03051313 4133   00   0     1.17430609E+00   # C9 mu+mu-
 03051313 4133   00   2     1.17418399E+00   # C9 mu+mu-
 03051313 4233   00   2    -1.17108171E-06   # C9' mu+mu-
 03051313 4137   00   0    -3.99876551E+00   # C10 mu+mu-
 03051313 4137   00   2    -4.00214849E+00   # C10 mu+mu-
 03051313 4237   00   2     3.15382569E-05   # C10' mu+mu-
 03051212 4137   00   0     1.50585488E+00   # C11 nu_1 nu_1
 03051212 4137   00   2     1.50666978E+00   # C11 nu_1 nu_1
 03051212 4237   00   2    -7.59176692E-06   # C11' nu_1 nu_1
 03051414 4137   00   0     1.50585488E+00   # C11 nu_2 nu_2
 03051414 4137   00   2     1.50666978E+00   # C11 nu_2 nu_2
 03051414 4237   00   2    -7.59176313E-06   # C11' nu_2 nu_2
 03051616 4137   00   0     1.50585488E+00   # C11 nu_3 nu_3
 03051616 4137   00   2     1.50666978E+00   # C11 nu_3 nu_3
 03051616 4237   00   2    -7.59069707E-06   # C11' nu_3 nu_3
Block IMFWCOEF Q=  1.60000000E+02  # Im(Wilson coefficients) at scale Q
#    id        order  M        value         comment
     0305 4422   00   0     3.86063341E-07   # C7
     0305 4422   00   2     5.27764156E-07   # C7
     0305 4322   00   2     7.44832997E-08   # C7'
     0305 6421   00   0     3.30688670E-07   # C8
     0305 6421   00   2    -9.05214596E-08   # C8
     0305 6321   00   2     3.07018313E-08   # C8'
 03051111 4133   00   2     9.34223902E-08   # C9 e+e-
 03051111 4233   00   2     1.17180260E-08   # C9' e+e-
 03051111 4137   00   2     9.22068477E-07   # C10 e+e-
 03051111 4237   00   2    -3.15500038E-07   # C10' e+e-
 03051313 4133   00   2     9.34223035E-08   # C9 mu+mu-
 03051313 4233   00   2     1.17180258E-08   # C9' mu+mu-
 03051313 4137   00   2     9.22068568E-07   # C10 mu+mu-
 03051313 4237   00   2    -3.15500038E-07   # C10' mu+mu-
 03051212 4137   00   2    -2.11867583E-07   # C11 nu_1 nu_1
 03051212 4237   00   2     7.59459844E-08   # C11' nu_1 nu_1
 03051414 4137   00   2    -2.11867577E-07   # C11 nu_2 nu_2
 03051414 4237   00   2     7.59459844E-08   # C11' nu_2 nu_2
 03051616 4137   00   2    -2.11866121E-07   # C11 nu_3 nu_3
 03051616 4237   00   2     7.59459844E-08   # C11' nu_3 nu_3
