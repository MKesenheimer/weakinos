# SUSY Les Houches Accord 1.0 - MSSM spectrum + Decays
# SPheno 2.2.3
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.unizh.ch
# Created: 25.09.2008,  15:48
#
Block SPINFO         # Program information
     1   SPheno      # spectrum calculator
     2   2.2.3       # version number
#
#Block SPhenoINFO     # SPheno specific information
#    1      2         # using 2-loop RGEs
#    2      1         # using running masses for boundary conditions at mZ
Block MODSEL  # Model selection
    1    1    # mSUGRA model
Block MINPAR  # Input parameters
    1    7.00000000E+01  # m0             29
    2    2.50000000E+02  # m12            30
    3    1.04054857E+01  # tanb at m_Z    31
    4    1.00000000E+00  # Sign(mu)       32
    5   -3.00000000E+02  # A0             33
#
Block SMINPUTS  # SM parameters
         1     1.27931500E+02  # alpha_em^-1(MZ)^MSbar      13
         2     1.16637000E-05  # G_mu [GeV^-2]              14
         3     1.17600000E-01  # alpha_s(MZ)^MSbar          15
         4     9.11876000E+01  # m_Z(pole)                  16
         5     4.20000000E+00  # m_b(m_b), MSbar            28
         6     1.71200000E+02  # m_t(pole)                  27
         7     1.77684000E+00  # m_tau(pole)                26

Block EXTPAR    #  **
         0     1.00000000E+00  #          35    **
        33     3.30000000E+00  #          50    **
        49     4.90000000E+01  #          62    **
        53     5.30000000E+01  # N5(3)    65    ** OK
        67     6.70000000E+01  #          72    **
        68     6.80000000E+01  #          73    **
        69     6.90000000E+01  #          74    ** OK
        70     7.00000000E+01  #              ** WRONG - missing

Block QEXTPAR    #  **
         1     1.00000000E+00  #              **
         2     2.00000000E+00  #              **
         3     3.00000000E+00  #              **
        11     1.10000000E+01  #              **
        12     1.20000000E+01  #              **
        13     1.30000000E+01  #              **
        47     4.70000000E+01  #              **

Block MASS  # Mass spectrum
#   PDG code      mass          particle
        24     8.03848627E+01  # W+
        25     1.08040876E+02  # h0
        35     3.89354899E+02  # H0
        36     3.88844776E+02  # A0
        37     3.97520808E+02  # H+
   1000001     5.68838451E+02  # ~d_L
   2000001     5.45242169E+02  # ~d_R
   1000002     5.63323368E+02  # ~u_L
   2000002     5.45425859E+02  # ~u_R
   1000003     5.68838522E+02  # ~s_L
   2000003     5.45236726E+02  # ~s_R
   1000004     5.63328315E+02  # ~c_L
   2000004     5.45416876E+02  # ~c_R
   1000005     5.16901339E+02  # ~b_1
   2000005     5.44935162E+02  # ~b_2
   1000006     4.42518564E+02  # ~t_1
   2000006     5.65855921E+02  # ~t_2
   1000011     1.90660786E+02  # ~e_L-
   2000011     1.26177277E+02  # ~e_R-
   1000012     1.73179672E+02  # ~nu_eL
   1000013     1.90678032E+02  # ~mu_L-
   2000013     1.26122625E+02  # ~mu_R-
   1000014     1.73172885E+02  # ~nu_muL
   1000015     1.11140270E+02  # ~tau_1-
   2000015     1.94666072E+02  # ~tau_2-
   1000016     1.71262280E+02  # ~nu_tauL
   1000021     6.06482068E+02  # ~g
   1000022     9.79251254E+01  # ~chi_10
   1000023     1.82396027E+02  # ~chi_20
   1000025    -3.62598042E+02  # ~chi_30
   1000035     3.80736776E+02  # ~chi_40
   1000024     1.81973013E+02  # ~chi_1+
   1000037     3.81460640E+02  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -1.11966411E-01   # alpha
Block Hmix Q=  1.00000000E+03  # Higgs mixing parameters
   1    3.57577731E+02  # mu
   2    1.00000000E+01  # tan[beta](Q)
   3    2.43103272E+02  # v(Q)
   4    1.39746126E+05  # m^2_A(Q)
Block stopmix  # stop mixing matrix
   1  1     4.57945887E-01   # R_st(1,1)        160
   1  2     8.88980070E-01   # R_st(1,2)        162
   2  1    -8.88980070E-01   # R_st(2,1)        161
   2  2     4.57945887E-01   # R_st(2,2)        163
Block sbotmix  # sbottom mixing matrix
   1  1     8.96399268E-01   # R_sb(1,1)        164
   1  2     4.43247507E-01   # R_sb(1,2)        166
   2  1    -4.43247507E-01   # R_sb(2,1)        165
   2  2     8.96399268E-01   # R_sb(2,2)        167
Block staumix  # stau mixing matrix
   1  1     2.98848252E-01   # R_sta(1,1)       156
   1  2     9.54300646E-01   # R_sta(1,2)       158
   2  1    -9.54300646E-01   # R_sta(2,1)       157
   2  2     2.98848252E-01   # R_sta(2,2)       159
Block Nmix  # neutralino mixing matrix
   1  1    -9.85377596E-01   # N(1,1)
   1  2     5.89966467E-02   # N(1,2)
   1  3    -1.49834044E-01   # N(1,3)
   1  4     5.56789690E-02   # N(1,4)
   2  1    -1.06662944E-01   # N(2,1)
   2  2    -9.43572571E-01   # N(2,2)
   2  3     2.71991494E-01   # N(2,3)
   2  4    -1.55930905E-01   # N(2,4)
   3  1     6.02396509E-02   # N(3,1)
   3  2    -8.99159886E-02   # N(3,2)
   3  3    -6.95521064E-01   # N(3,3)
   3  4    -7.10307503E-01   # N(3,4)
   4  1     1.18428010E-01   # N(4,1)
   4  2    -3.13217679E-01   # N(4,2)
   4  3    -6.47935827E-01   # N(4,3)
   4  4     6.84140816E-01   # N(4,4)
Block Umix  # chargino U mixing matrix
   1  1    -9.16840379E-01   # U(1,1)
   1  2     3.99253953E-01   # U(1,2)
   2  1     3.99253953E-01   # U(2,1)     150
   2  2     9.16840379E-01   # U(2,2)     151
Block Vmix  # chargino V mixing matrix
   1  1    -9.73205395E-01   # V(1,1)     152
   1  2     2.29937513E-01   # V(1,2)     153
   2  1     2.29937513E-01   # V(2,1)     154
   2  2     9.73205395E-01   # V(2,2)     155
Block gauge Q=  1.00000000E+03  # (SUSY scale)
   1    3.63627762E-01  # g'(Q)^DRbar
   2    6.48022350E-01  # g(Q)^DRbar
   3    1.07924824E+00  # g3(Q)^DRbar
Block Yu Q=  1.00000000E+03  # (SUSY scale)
  1  1     9.00514746E-06   # Y_u(Q)^DRbar
  2  2     3.60205999E-03   # Y_c(Q)^DRbar
  3  3     9.06594294E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.73838792E-04   # Y_d(Q)^DRbar
  2  2     3.47677618E-03   # Y_s(Q)^DRbar
  3  3     1.45209040E-01   # Y_b(Q)^DRbar
Block Ye Q=  1.00000000E+03  # (SUSY scale)
  1  1     3.02792953E-05   # Y_e(Q)^DRbar
  2  2     6.26080127E-03   # Y_mu(Q)^DRbar
  3  3     1.05301319E-01   # Y_tau(Q)^DRbar
Block Au Q=  1.00000000E+03  # (SUSY scale)
  1  1    -5.00058028E+02   # A_u(Q)^DRbar
  2  2    -5.00052695E+02   # A_c(Q)^DRbar
  3  3    -2.46832508E+02   # A_t(Q)^DRbar
Block Ad Q=  1.00000000E+03  # (SUSY scale)
  1  1    -6.70532174E+02   # A_d(Q)^DRbar
  2  2    -6.70526894E+02   # A_s(Q)^DRbar
  3  3    -5.82087115E+02   # A_b(Q)^DRbar
Block Ae Q=  1.00000000E+03  # (SUSY scale)
  1  1    -3.45520645E+02   # A_e(Q)^DRbar
  2  2    -3.45507110E+02   # A_mu(Q)^DRbar
  3  3    -3.41692275E+02   # A_tau(Q)^DRbar
Block MSOFT Q=  1.00000000E+03  # soft SUSY breaking masses at Q
   1    1.04021792E+02  # M_1
   2    1.94743157E+02  # M_2
   3    5.71234546E+02  # M_3
  21    2.56793857E+04  # M^2_(H,d)
  22   -1.17802275E+05  # M^2_(H,u)
  31    1.81805165E+02  # M_(L,11)
  32    1.81799135E+02  # M_(L,22)
  33    1.80105439E+02  # M_(L,33)
  34    1.16333332E+02  # M_(E,11)
  35    1.16314204E+02  # M_(E,22)
  36    1.10835762E+02  # M_(E,33)
  41    5.24813521E+02  # M_(Q,11)
  42    5.24811294E+02  # M_(Q,22)
  43    4.80636219E+02  # M_(Q,33)
  44    5.05857573E+02  # M_(U,11)
  45    5.05855500E+02  # M_(U,22)
  46    4.08745269E+02  # M_(U,33)
  47    5.03724964E+02  # M_(D,11)
  48    5.03722342E+02  # M_(D,22)
  49    4.99511050E+02  # M_(D,33)
#Block SPhenoLowEnergy  # low energy observables
#    1    4.41190700E+00   # BR(b -> s gamma)
#    2    7.40053852E-09   # (g-2)_muon
#    3    2.87132998E-04   # Delta(rho)
