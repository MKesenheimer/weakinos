# SUSY Les Houches Accord - MSSM spectrum + Decays
# SPheno2.1
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.unizh.ch
Block MODSEL  # Select model
    1    1
Block MINPAR  # Input parameters
    1   0.100000000E+03  # m0      
    2   0.250000000E+03  # m12     
    3   0.100000000E+02  # tanb    
    4   0.100000000E+01  # Sign(mu)
    5  -0.100000000E+03  # A0
#
Block SMINPUTS  # Mass spectrum
          4   0.911870000E+02  # MZ
          5   0.425000000E+01  # mb(mb)
          6   0.175000000E+03  # t
#
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         25   0.112265330E+03  # h0
         35   0.400481169E+03  # H0
         36   0.399769788E+03  # A0
         37   0.409091863E+03  # H+
    1000001   0.569498044E+03  # ~d_L
    2000001   0.547601268E+03  # ~d_R
    1000002   0.564892619E+03  # ~u_L
    2000002   0.547790210E+03  # ~u_R
    1000003   0.569498964E+03  # ~s_L
    2000003   0.547594947E+03  # ~s_R
    1000004   0.564902784E+03  # ~c_L
    2000004   0.547775859E+03  # ~c_R
    1000005   0.516225720E+03  # ~b_1
    2000005   0.547471349E+03  # ~b_2
    1000006   0.398749215E+03  # ~t_1
    2000006   0.589079372E+03  # ~t_2
    1000011   0.206630723E+03  # ~e_L-
    2000011   0.143872558E+03  # ~e_R-
    1000012   0.190463420E+03  # ~nu_eL
    1000013   0.206645846E+03  # ~mu_L-
    2000013   0.143838140E+03  # ~mu_R-
    1000014   0.190460263E+03  # ~nu_muL
    1000015   0.134514453E+03  # ~tau_L-
    2000015   0.210401949E+03  # ~tau_R-
    1000016   0.189568253E+03  # ~nu_tauL
    1000021   0.594297332E+03  # ~g
    1000022   0.977190851E+02  # ~chi_10
    1000023   0.183160462E+03  # ~chi_20
    1000025  -0.364508626E+03  # ~chi_30
    1000035   0.382044178E+03  # ~chi_40
    1000024   0.181597212E+03  # ~chi_1+
    1000037   0.382233627E+03  # ~chi_2+
# Higgs mixing
Block alpha # Higgs mixing parameters
     -0.113991975E+00  # alpha
Block hmix Q= 0.484786694E+03  # Higgs mixing parameters
   1   0.357788560E+03  # mu
Block stopmix  # stop mixing matrix
  1  1   0.552785989E+00  # R_st(1,1)
  1  2   0.833323257E+00  # R_st(1,2)
  2  1  -0.833323257E+00  # R_st(2,1)
  2  2   0.552785989E+00  # R_st(2,2)
Block sbotmix  # sbottom mixing matrix
  1  1   0.951803418E+00  # R_sb(1,1)
  1  2   0.306708745E+00  # R_sb(1,2)
  2  1  -0.306708745E+00  # R_sb(2,1)
  2  2   0.951803418E+00  # R_sb(2,2)
Block staumix  # stau mixing matrix
  1  1   0.270045675E+00  # R_sta(1,1)
  1  2   0.962847513E+00  # R_sta(1,2)
  2  1  -0.962847513E+00  # R_sta(2,1)
  2  2   0.270045675E+00  # R_sta(2,2)
Block nmix  # neutralino mixing matrix
  1  1   0.985490276E+00  # N(1,1)  
  1  2  -0.550659970E-01  # N(1,2) 
  1  3   0.150665297E+00  # N(1,3)  
  1  4  -0.554672864E-01  # N(1,4)  
  2  1  -0.104134398E+00  # N(2,1) SIGN 
  2  2  -0.941334568E+00  # N(2,2) SIGN 
  2  3   0.278005195E+00  # N(2,3) SIGN
  2  4  -0.160494144E+00  # N(2,4) SIGN
  3  1  -0.611177267E-01  # N(3,1) 
  3  2   0.906968322E-01  # N(3,2)  
  3  3   0.694659828E+00  # N(3,3)  
  3  4   0.710975690E+00  # N(3,4)  
  4  1  -0.119287746E+00  # N(4,1)  
  4  2   0.320360815E+00  # N(4,2)  
  4  3   0.646112068E+00  # N(4,3)  
  4  4  -0.682406461E+00  # N(4,4)  
Block Umix  # chargino U mixing matrix
  1  1  -0.917012658E+00  # U(1,1)
  1  2   0.398858101E+00  # U(1,2)
  2  1   0.398858101E+00  # U(2,1)
  2  2   0.917012658E+00  # U(2,2)
Block Vmix  # chargino V mixing matrix
  1  1  -0.970316957E+00  # V(1,1)
  1  2   0.241836727E+00  # V(1,2)
  2  1   0.241836727E+00  # V(2,1)
  2  2   0.970316957E+00  # V(2,2)
Block gauge Q= 0.484786694E+03  # (SUSY scale)
   1   0.361493012E+00  # g'(Q)^DRbar
   2   0.644983706E+00  # g(Q)^DRbar
   3   0.110334304E+01  # g3(Q)^DRbar
Block au Q= 0.484786694E+03     # (SUSY scale)
  3  3   0.1000000000+02  # At
Block ad Q= 0.484786694E+03     # (SUSY scale)
  3  3   0.1000000000+02  # Ab
Block ae Q= 0.484786694E+03     # (SUSY scale)
  3  3   0.1000000000+02  # Atau
Block yatop Q= 0.484786694E+03   
   1   0.893480549E+00  # Y_t(Q)^DRbar
   2  -0.501095259E+03  # A_t(Q)^DRbar
Block yabot Q= 0.484786694E+03   
   1   0.125450458E+00  # Y_b(Q)^DRbar
   2  -0.851503980E+03  # A_b(Q)^DRbar
Block yatau Q= 0.484786694E+03   
   1   0.100790510E+00  # Y_tau(Q)^DRbar
   2  -0.250194146E+03  # A_tau(Q)^DRbar
Block SPINFO         # Program information
    1    SPheno      # spectrum calculator
    2    2.1         # version number
DECAY        25 1.017523300E-11  # h0 decays
#          BR          NDA      ID1     ID2    ID3
    6.000000000E-01    3          5        -5      16  # b b 
    5.000000000E-02    3          2        -1      15  # c c
    5.000000000E-02    3         -2         1     -15  # s s
    5.500000000E-02    3          2        -1      13  # d d
DECAY   1000022 1.017523300E-11  # chi^0_1 decays
#          BR          NDA      ID1     ID2    ID3
    6.000000000E-01    3          5        -5      16  # b b nu_tau
    5.000000000E-02    3          2        -1      15  # u d tau-
    5.000000000E-02    3         -2         1     -15  # u d tau+
    5.500000000E-02    3          2        -1      13  # u d mu-
    5.500000000E-02    3         -2         1     -13  # u d mu+
    5.000000000E-03    3          2        -1      11  # u d e-
    5.000000000E-03    3         -2         1     -11  # u d e+
    5.000000000E-02    3         11       -15      16  # e- tau+ nu_tau
    5.000000000E-02    3        -11        15      16  # e+ tau- nu_tau
    4.000000000E-02    3         13       -15      16  # mu- tau+ nu_tau
    4.000000000E-02    3        -13        15      16  # mu+ tau- nu_tau
    
