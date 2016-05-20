########################################################################
#                          -*- Makefile -*-                            #
########################################################################
# General Conventions for Makefiles
SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .c .f .F .cc .cpp .h .hh .inc .o .a
.DEFAULT_GOAL := all

########################################################################
## Flags

## Compiler and additional compiler Flags
FC  =  gfortran
CXX =  g++
CC  =  gcc
FCFLAGS  = -g
CXXFLAGS = -g
CFLAGS   = -g
LDFLAGS  = -ff2c -g

# recommended compiler flags
REC_FCFLAGS   = -fno-automatic -fno-range-check
REC_FCFLAGS  += -ffixed-line-length-none -lgfortran -DU77EXT=0 -DQuad=0
REC_FCFLAGS  += -ff2c -fno-second-underscore
REC_FCFLAGS  += $(FCFLAGS)
REC_CXXFLAGS  = -fomit-frame-pointer -ffast-math -Wall -m64
REC_CXXFLAGS += $(CXXFLAGS)
REC_CFLAGS    = -fomit-frame-pointer -ffast-math -Wall -m64
REC_CFLAGS   += -DNOUNDERSCORE=0 -DBIGENDIAN=0
REC_CFLAGS   += $(CFLAGS)

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
  #Mac OSX
  REC_CFLAGS   += -stdlib=libstdc++ -mmacosx-version-min=10.6 -Qunused-arguments
  REC_CXXFLAGS += -stdlib=libstdc++ -mmacosx-version-min=10.6 -Qunused-arguments
endif

## PDF
## choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
## wheter to link precompiled LHAPDF library statically or compile it 
## from source (static, none)
PDF = lhapdf
PDFSTATIC = none

## path to static lhapdf library
#STATICLIBSLHAPDF = /opt/lib/libLHAPDF.a

## choose analysis: none, default
ANALYSIS = default

## path to LHAPDF config executable
LHAPDF_CONFIG = lhapdf-config

## path to fastjet config executable
FASTJET_CONFIG = fastjet-config

## Pythia8
# if you want to use Pythia8 instead of Pythia6 uncomment the following
#USE_PYTHIA8 = yes
# path to pythia8-config executable
# (uncomment only if you want to use a locally installed Pythia8)
#PYTHIA8_CONFIG = pythia8-config

# (not necessary anymore)
## path to fortran and c++ libraries
# (use $locate libstdc++.a to find the library)
#ifeq ($(UNAME), Darwin)
#  #Mac OSX
#  LIBGFORTRANPATH = /usr/local/lib
#  LIBSTDCPPPATH   = /usr/local/lib
#else
#  #Linux 
#  LIBGFORTRANPATH = /usr/lib/gcc/x86_64-linux-gnu/4.8
#  LIBSTDCPPPATH   = /usr/lib/gcc/x86_64-linux-gnu/4.8
#endif

## warning for type-conversions -> basically useless, as those occur in
## too many places
#WARN  = -Wconversion -Wall -Wtabs -Wall -Wimplicit-interface
## -fbounds-check sometimes causes a weird error due to non-lazy
## evaluation of boolean in gfortran.
#WARN += -fbounds-check
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT = -O2

### generate directory build, if not yet existing
$(shell mkdir -p build)

### version control
GV = $(shell git rev-parse --short HEAD)
$(shell echo "      subroutine print_version" > ../version.f)
$(shell echo "        print*,'        git version $(GV)'" >> ../version.f)
$(shell echo "      end" >> ../version.f)

########################################################################
## Runtime flags

## Preprocessor
# it might be advisable to use the -ffree-line-length-none 
# or -ffixed-line-length-none options
CPP = -cpp

## For debugging uncomment the following
#DEBUG = -ggdb -pg -D DEBUG
## more verbose debugging
DEBUG = -ggdb -pg -D DEBUG -D DEBUGV
## If you wish to show all possible debug output
#DEBUG = -ggdb -pg -D DEBUG -D DEBUGV -D DEBUGQ -D CHECK

## More debugging flags
# Check UV finiteness in FormCalc Virtuals: -D CHECKUV
# Check if MadGraph and FormCalc born amplitudes are the same during 
# calculation and abort if they are not: -D CHECK_FORM_MAD
# Note: Use this only for testing born amplitudes, not for the final 
# program and not while checking real amplitudes!
#USRFLAGS += -D CHECK_FORM_MAD

## Additional and mandatory flags
# Use MadGraph Borns (default): -D MAD_BORN 
# Use FormCalc Borns: -D FORM_BORN
USRFLAGS += -D MAD_BORN

## Treatment of quark masses and squark widths
# The current version uses vanishing masses for the incoming quarks and
# zero widths for the internal squarks.
# For later software versions the possibility exists to make the
# calculations with finite b-quark mass or finite squark widths.
USRFLAGS += -D NEGLECTBMASS
USRFLAGS += -D NEGLECT_SWIDTH
# If you want to change these flags, you probably have to change a lot
# for the FormCalc amplitudes. For example, the b-quark amplitudes are 
# currently calculated by crossing relations from d-quark amplitudes.
# Since FormCalc does not use squark widths, you have to introduce
# squark widths by hand.

## Subtraction schemes
# There are divergences originating from on-shell produced squark lines
# in the real processes. The resonant contribution can easily exceed 
# the full NLO corrections. Keeping it as part of the real corrections
# would spoil the predictivity of the NLO calculation. Therefore, these
# contributions have to be removed in a consistent way.
# There exist several methods to cope with them, the most relevant 
# being the following:
# - Diagram Removal type I, implemented with preprocessor flag -D DR_I
# - Diagram Removal type II, implemented with preprocessor flag -D DR_II
# - Diagram Subtraction, implemented with preprocessor flag -D DSUB_I,
#   no correction of the jacobian, which is in principal not ok, but 
#   could be used for testing the convergence.
# - Diagram Subtraction with corrected jacobian, flag -D DSUB_II.
#   With this flag defined, the real amplitude is split into two parts.
#   In Reals.f the matrix element |M_NR|**2 + 2 M_R* M_NR is calculated,
#   in Reals_osres.f the remaining terms and the on-shell counter terms
#   are added: |M_R|**2 - BW*corrfac*|M_R^OS|**2. The latter gets 
#   integrated seperately in sigosres.f and is included at a later point
#   to the main result.
#USRFLAGS += -D DR_I
#USRFLAGS += -D DR_II
#USRFLAGS += -D DSUB_I
#USRFLAGS += -D DSUB_II_TEST
USRFLAGS += -D DSUB_II

########################################################################
#           -*- no editing is required below this line -*-             #
########################################################################

########################################################################
## Runtime flags which should not be changed (don't touch this unless
## you know what you do!) 
CHANNEL += -D NINJ

########################################################################
## Paths

WORKINGDIR = $(shell pwd)/..
PWHGDIR    = $(WORKINGDIR)/..
SUBWORKINGDIR = $(shell pwd)

# directories
TOOLS   = $(WORKINGDIR)/Tools
REAL    = $(SUBWORKINGDIR)/MadGraph_Reals
VIRT    = $(SUBWORKINGDIR)/FormCalc_Virtuals
VIRTSUB = $(VIRT)/uubar_n1n1_squaredme $(VIRT)/ddbar_n1n1_squaredme	\
	$(VIRT)/uubar_n1n2_squaredme $(VIRT)/ddbar_n1n2_squaredme	\
	$(VIRT)/uubar_n2n2_squaredme $(VIRT)/ddbar_n2n2_squaredme	\
	$(VIRT)/uubar_n1n3_squaredme $(VIRT)/ddbar_n1n3_squaredme	\
	$(VIRT)/uubar_n1n4_squaredme $(VIRT)/ddbar_n1n4_squaredme	\
	$(VIRT)/uubar_n2n3_squaredme $(VIRT)/ddbar_n2n3_squaredme	\
	$(VIRT)/uubar_n2n4_squaredme $(VIRT)/ddbar_n2n4_squaredme	\
	$(VIRT)/uubar_n3n3_squaredme $(VIRT)/ddbar_n3n3_squaredme	\
	$(VIRT)/uubar_n3n4_squaredme $(VIRT)/ddbar_n3n4_squaredme	\
	$(VIRT)/uubar_n4n4_squaredme $(VIRT)/ddbar_n4n4_squaredme
REALSUB  = $(REAL)/SubProcesses_nInJ
MODIFIED = $(WORKINGDIR)/POWHEG_mod
SCRIPTS  = $(WORKINGDIR)/Scripts

FJ     = $(TOOLS)/fastjet
PY     = $(TOOLS)/pythia6425
PY8    = $(TOOLS)/pythia8215
PYF    = $(TOOLS)/pythiafortran
HER    = $(TOOLS)/herwig
LT     = $(TOOLS)/LoopTools-2.12
SLHA   = $(TOOLS)/SLHALib-2.2
DHELAS = $(TOOLS)/DHELAS
MADMODEL  = $(REAL)/MODEL

# includes
UINCLUDE    = $(WORKINGDIR)/include
SINCLUDE    = $(SUBWORKINGDIR)/include
PINCLUDE    = $(PWHGDIR)/include
VINCLUDE    = $(VIRT)/include
MDINCLUDE   = $(REAL)/include
SLHAINCLUDE = $(SLHA)/src
LTINCLUDE   = $(LT)/src/include
PY8INCLUDE  = $(PY8)/include

########################################################################
## search for the files and set paths

vpath %.f $(PWHGDIR) $(WORKINGDIR) $(TOOLS)/functions $(DHELAS) $(PY) $(PYF) $(HER)
vpath %.F $(WORKINGDIR) $(TOOLS)/FormCalcTools $(DHELAS)
vpath %.F $(VIRT)/RenConsts $(VIRTSUB)
vpath %.F $(MADMODEL) $(REALSUB)
vpath %.f $(MADMODEL) $(REALSUB)
vpath %.f $(MODIFIED)
vpath %.cc $(FJ) $(PY8) $(PYF)
vpath %.o $(SUBWORKINGDIR)/build

########################################################################
## Source files

### PDFs ###
ifeq ("$(PDF)","lhapdf")
  PDFPACK       = lhapdfif.o
  LIBSLHAPDF    = -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)		\
                  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF
ifeq  ("$(PDFSTATIC)","static")
  LIBSLHAPDF += $(STATICLIBSLHAPDF)
endif
  LIBS += $(LIBSLHAPDF) 
else
  PDFPACK = mlmpdfif.o hvqpdfpho.o
endif

### ANALYSIS ###
ifeq ("$(ANALYSIS)","default")
  USE_FASTJET  = yes
  PWHGANALYSIS = Analysis.o
endif
ifeq ("$(ANALYSIS)","none")
  PWHGANALYSIS = Analysis-dummy.o
endif

ifdef USE_FASTJET
  ifeq ($(FASTJET_CONFIG),)
    $(error cannot find fastjet-config (did you install fastjet?))
  endif
  LIBS += $(shell $(FASTJET_CONFIG) --libs --plugins )
  FJCXXFLAGS  += $(shell $(FASTJET_CONFIG) --cxxflags)
  PWHGANALYSIS+= fastjetfortran.o
endif

# the user choosed to use Pythia8
ifeq ("$(USE_PYTHIA8)","yes")
# the user choosed to use a locally installed Pythia8
ifdef PYTHIA8_CONFIG
  LIBS += $(shell $(PYTHIA8_CONFIG) --libs --plugins )
  PYCXXFLAGS  += $(shell $(PYTHIA8_CONFIG) --cxxflags)
else
  LIBS += $(TOOLS)/libpythia8.a
  PYCXXFLAGS += -I$(PY8INCLUDE)
endif
  PWHGANALYSIS+= pythiafortran.o
endif

### USER Files ###
USER = Born.o Born_phsp.o Born_lh.o Reals.o Reals_lh.o Virtuals.o 	\
       init_processes.o init_couplings.o Reals_osres.o 			\
       Real_osres_phsp.o phi1_2.o sigosres.o version.o

### Functions in Tools/functions/
FUNC = funcprocess.o Print.o funcbasic.o off_to_on.o

### MadGraph mandatory files
MMADOBJ = couplings.o hardstop.o printout.o read_slha.o 		\
	set_general_couplings.o 

### Formcalc Files in Tools/FormcalcTools/###
FORMCALC = VecSet.o

### PWHG FILES in ../ ###
PWHG = pwhg_main.o pwhg_init.o btilde.o utils.o gen_Born_phsp.o		\
	test_Sudakov.o pt2maxreg.o maxrat.o Bornzerodamp.o		\
	random.o boostrot.o bra_ket_subroutines.o			\
	init_phys.o powheginput.o sigcollremn.o				\
	pwhg_bookhist-multi.o pwhg_analysis_driver.o checkmomzero.o	\
	setstrongcoupl.o integrator.o mwarn.o				\
	sigcollsoft.o validflav.o decidemerging.o          		\
	loadstorenlo.o setlocalscales.o ubprojections.o      		\
	mint_upb.o reshufflemoms.o opencount.o

# modified POWHEG-files
PWHG += cernroutines_mod.o find_regions_mod.o gen_real_phsp_mod.o	\
	sigborn_mod.o sigreal_mod.o sigsoftvirt_mod.o			\
	sigvirtual_mod.o bbinit_mod.o LesHouches_mod.o 			\
	LesHouchesreg_mod.o gen_index_mod.o lhefwrite_mod.o		\
	gen_radiation_mod.o sigremnants_mod.o pwhgreweight_mod.o 	\
	pdfcalls.o newunit_mod.o

### LHEF ###
# analyze LHEF output
LHEF = lhef_analysis.o boostrot.o random.o cernroutines_mod.o		\
     opencount.o powheginput.o pwhg_bookhist-multi.o 			\
     lhefread.o newunit_mod.o pwhg_analysis_driver.o 			\
     funcprocess.o init_couplings.o Print.o				\
     bra_ket_subroutines.o off_to_on.o funcbasic.o

### PYTHIA ###
# target to read event file, shower events with PYTHIA + analysis
# NOTE: pythia-6.4.25-decay.f is based on POWHEG-BOX/pythia-6.4.25.f,
# to which we applied a simple patch to enable slepton decays
# (for clarity a separate patch file is included which is not needed otherwise)
PYTHIA = boostrot.o powheginput.o init_couplings.o funcprocess.o 	\
	pwhg_bookhist-multi.o lhefread.o newunit_mod.o pdfdummies.o	\
	pwhg_analysis_driver.o random.o cernroutines_mod.o opencount.o	\
	bra_ket_subroutines.o off_to_on.o funcbasic.o Print.o
ifeq ("$(USE_PYTHIA8)","yes")
  # TODO:
  PYTHIA += main-PYTHIA_mod.o
else
  PYTHIA += main-PYTHIA.o setup-PYTHIA-lhef.o
  PYTHIA += pythia-6.4.25-decay.o
endif

### HERWIG ###
# this is work in progress
HERWIG = main-HERWIG.o setup-HERWIG-lhef.o herwig6510.o boostrot.o 	\
     opencount.o powheginput.o pwhg_bookhist-multi.o random.o		\
     lhefread.o newunit_mod.o pwhg_analysis_driver.o 			\
     funcprocess.o init_couplings.o cernroutines_mod.o			\
     bra_ket_subroutines.o off_to_on.o funcbasic.o

# MadGraph born amplitudes
MADOBJ += sborn_proc.o
MADOBJ += b_cl_nInJ_001_001.o  b_cl_nInJ_006_001.o  b_sf_nInJ_001_003.o  b_sf_nInJ_006_003.o	\
	b_cl_nInJ_001_002.o  b_cl_nInJ_006_002.o  b_sf_nInJ_001_004.o  b_sf_nInJ_006_004.o	\
	b_cl_nInJ_002_001.o  b_cl_nInJ_007_001.o  b_sf_nInJ_002_003.o  b_sf_nInJ_007_003.o	\
	b_cl_nInJ_002_002.o  b_cl_nInJ_007_002.o  b_sf_nInJ_002_004.o  b_sf_nInJ_007_004.o	\
	b_cl_nInJ_003_001.o  b_cl_nInJ_008_001.o  b_sf_nInJ_003_003.o  b_sf_nInJ_008_003.o	\
	b_cl_nInJ_003_002.o  b_cl_nInJ_008_002.o  b_sf_nInJ_003_004.o  b_sf_nInJ_008_004.o	\
	b_cl_nInJ_004_001.o  b_cl_nInJ_009_001.o  b_sf_nInJ_004_003.o  b_sf_nInJ_009_003.o	\
	b_cl_nInJ_004_002.o  b_cl_nInJ_009_002.o  b_sf_nInJ_004_004.o  b_sf_nInJ_009_004.o	\
	b_cl_nInJ_005_001.o  b_cl_nInJ_010_001.o  b_sf_nInJ_005_003.o  b_sf_nInJ_010_003.o	\
	b_cl_nInJ_005_002.o  b_cl_nInJ_010_002.o  b_sf_nInJ_005_004.o  b_sf_nInJ_010_004.o
MADOBJ += sborn_cl_nInJ_001.o  sborn_cl_nInJ_007.o  sborn_sf_nInJ_002.o  sborn_sf_nInJ_008.o	\
	sborn_cl_nInJ_002.o  sborn_cl_nInJ_008.o  sborn_sf_nInJ_003.o  sborn_sf_nInJ_009.o	\
	sborn_cl_nInJ_003.o  sborn_cl_nInJ_009.o  sborn_sf_nInJ_004.o  sborn_sf_nInJ_010.o	\
	sborn_cl_nInJ_004.o  sborn_cl_nInJ_010.o  sborn_sf_nInJ_005.o				\
	sborn_cl_nInJ_005.o  sborn_sf_nInJ_006.o						\
	sborn_cl_nInJ_006.o  sborn_sf_nInJ_001.o  sborn_sf_nInJ_007.o

# real Amplitudes
MADOBJ += sreal_proc.o 
MADOBJ += r_ddx_nInJg.o r_dg_nInJd.o r_dxd_nInJg.o r_dxg_nInJdx.o	\
	  r_gd_nInJd.o r_gdx_nInJdx.o r_gu_nInJu.o   r_gux_nInJux.o	\
	  r_ug_nInJu.o r_uux_nInJg.o  r_uxg_nInJux.o r_uxu_nInJg.o

# on-shell subtraction scheme
MADOBJ += sreal_proc_res.o
MADOBJ += r_dg_nInJd_res.o r_dxg_nInJdx_res.o r_gd_nInJd_res.o		\
	   r_gdx_nInJdx_res.o r_gu_nInJu_res.o r_gux_nInJux_res.o 	\
	   r_ug_nInJu_res.o r_uxg_nInJux_res.o

### Virtual Amplitudes in Formcalc_Virtuals/ ###
VIRTUALS = CalcRenConst.o uubar_n1n1_abbr0h.o uubar_n1n1_abbr0s.o	\
          uubar_n1n1_abbr1a.o uubar_n1n1_abbr1h.o uubar_n1n1_abbr1s.o	\
          uubar_n1n1_born.o uubar_n1n1_born_Sfe5.o			\
          uubar_n1n1_box_Sfe5.o uubar_n1n1_box_Sfe5_Sfe6.o		\
          uubar_n1n1_self_Sfe5.o uubar_n1n1_self_Sfe5_Sfe6.o		\
          uubar_n1n1_self_Sfe5_Sfe6_Sfe7.o				\
          uubar_n1n1_SquaredME.o uubar_n1n1_vert_Sfe5_Sfe6.o		\
          uubar_n1n1_vert.o uubar_n1n1_vert_Sfe5.o			\
          ddbar_n1n1_abbr0h.o ddbar_n1n1_abbr0s.o            		\
          ddbar_n1n1_abbr1a.o ddbar_n1n1_abbr1h.o ddbar_n1n1_abbr1s.o	\
          ddbar_n1n1_born.o ddbar_n1n1_born_Sfe5.o			\
          ddbar_n1n1_box_Sfe5.o ddbar_n1n1_box_Sfe5_Sfe6.o		\
          ddbar_n1n1_self_Sfe5.o ddbar_n1n1_self_Sfe5_Sfe6.o		\
          ddbar_n1n1_self_Sfe5_Sfe6_Sfe7.o				\
          ddbar_n1n1_SquaredME.o ddbar_n1n1_vert_Sfe5_Sfe6.o		\
          ddbar_n1n1_vert.o ddbar_n1n1_vert_Sfe5.o
VIRTUALS += uubar_n1n2_abbr0h.o uubar_n1n2_abbr0s.o			\
          uubar_n1n2_abbr1a.o uubar_n1n2_abbr1h.o uubar_n1n2_abbr1s.o	\
          uubar_n1n2_born.o uubar_n1n2_born_Sfe5.o			\
          uubar_n1n2_box_Sfe5.o uubar_n1n2_box_Sfe5_Sfe6.o		\
          uubar_n1n2_self_Sfe5.o uubar_n1n2_self_Sfe5_Sfe6.o		\
          uubar_n1n2_self_Sfe5_Sfe6_Sfe7.o				\
          uubar_n1n2_SquaredME.o uubar_n1n2_vert_Sfe5_Sfe6.o		\
          uubar_n1n2_vert.o uubar_n1n2_vert_Sfe5.o			\
          ddbar_n1n2_abbr0h.o ddbar_n1n2_abbr0s.o            		\
          ddbar_n1n2_abbr1a.o ddbar_n1n2_abbr1h.o ddbar_n1n2_abbr1s.o	\
          ddbar_n1n2_born.o ddbar_n1n2_born_Sfe5.o			\
          ddbar_n1n2_box_Sfe5.o ddbar_n1n2_box_Sfe5_Sfe6.o		\
          ddbar_n1n2_self_Sfe5.o ddbar_n1n2_self_Sfe5_Sfe6.o		\
          ddbar_n1n2_self_Sfe5_Sfe6_Sfe7.o				\
          ddbar_n1n2_SquaredME.o ddbar_n1n2_vert_Sfe5_Sfe6.o		\
          ddbar_n1n2_vert.o ddbar_n1n2_vert_Sfe5.o
VIRTUALS += uubar_n2n2_abbr0h.o uubar_n2n2_abbr0s.o			\
          uubar_n2n2_abbr1a.o uubar_n2n2_abbr1h.o uubar_n2n2_abbr1s.o	\
          uubar_n2n2_born.o uubar_n2n2_born_Sfe5.o			\
          uubar_n2n2_box_Sfe5.o uubar_n2n2_box_Sfe5_Sfe6.o		\
          uubar_n2n2_self_Sfe5.o uubar_n2n2_self_Sfe5_Sfe6.o		\
          uubar_n2n2_self_Sfe5_Sfe6_Sfe7.o				\
          uubar_n2n2_SquaredME.o uubar_n2n2_vert_Sfe5_Sfe6.o		\
          uubar_n2n2_vert.o uubar_n2n2_vert_Sfe5.o			\
          ddbar_n2n2_abbr0h.o ddbar_n2n2_abbr0s.o            		\
          ddbar_n2n2_abbr1a.o ddbar_n2n2_abbr1h.o ddbar_n2n2_abbr1s.o	\
          ddbar_n2n2_born.o ddbar_n2n2_born_Sfe5.o			\
          ddbar_n2n2_box_Sfe5.o ddbar_n2n2_box_Sfe5_Sfe6.o		\
          ddbar_n2n2_self_Sfe5.o ddbar_n2n2_self_Sfe5_Sfe6.o		\
          ddbar_n2n2_self_Sfe5_Sfe6_Sfe7.o				\
          ddbar_n2n2_SquaredME.o ddbar_n2n2_vert_Sfe5_Sfe6.o		\
          ddbar_n2n2_vert.o ddbar_n2n2_vert_Sfe5.o
VIRTUALS += uubar_n1n3_abbr0h.o uubar_n1n3_abbr0s.o                     \
          uubar_n1n3_abbr1a.o uubar_n1n3_abbr1h.o uubar_n1n3_abbr1s.o   \
          uubar_n1n3_born.o uubar_n1n3_born_Sfe5.o                      \
          uubar_n1n3_box_Sfe5.o uubar_n1n3_box_Sfe5_Sfe6.o              \
          uubar_n1n3_self_Sfe5.o uubar_n1n3_self_Sfe5_Sfe6.o            \
          uubar_n1n3_self_Sfe5_Sfe6_Sfe7.o                              \
          uubar_n1n3_SquaredME.o uubar_n1n3_vert_Sfe5_Sfe6.o            \
          uubar_n1n3_vert.o uubar_n1n3_vert_Sfe5.o                      \
          ddbar_n1n3_abbr0h.o ddbar_n1n3_abbr0s.o                       \
          ddbar_n1n3_abbr1a.o ddbar_n1n3_abbr1h.o ddbar_n1n3_abbr1s.o   \
          ddbar_n1n3_born.o ddbar_n1n3_born_Sfe5.o                      \
          ddbar_n1n3_box_Sfe5.o ddbar_n1n3_box_Sfe5_Sfe6.o              \
          ddbar_n1n3_self_Sfe5.o ddbar_n1n3_self_Sfe5_Sfe6.o            \
          ddbar_n1n3_self_Sfe5_Sfe6_Sfe7.o                              \
          ddbar_n1n3_SquaredME.o ddbar_n1n3_vert_Sfe5_Sfe6.o            \
          ddbar_n1n3_vert.o ddbar_n1n3_vert_Sfe5.o
VIRTUALS += uubar_n1n4_abbr0h.o uubar_n1n4_abbr0s.o                     \
          uubar_n1n4_abbr1a.o uubar_n1n4_abbr1h.o uubar_n1n4_abbr1s.o   \
          uubar_n1n4_born.o uubar_n1n4_born_Sfe5.o                      \
          uubar_n1n4_box_Sfe5.o uubar_n1n4_box_Sfe5_Sfe6.o              \
          uubar_n1n4_self_Sfe5.o uubar_n1n4_self_Sfe5_Sfe6.o            \
          uubar_n1n4_self_Sfe5_Sfe6_Sfe7.o                              \
          uubar_n1n4_SquaredME.o uubar_n1n4_vert_Sfe5_Sfe6.o            \
          uubar_n1n4_vert.o uubar_n1n4_vert_Sfe5.o                      \
          ddbar_n1n4_abbr0h.o ddbar_n1n4_abbr0s.o                       \
          ddbar_n1n4_abbr1a.o ddbar_n1n4_abbr1h.o ddbar_n1n4_abbr1s.o   \
          ddbar_n1n4_born.o ddbar_n1n4_born_Sfe5.o                      \
          ddbar_n1n4_box_Sfe5.o ddbar_n1n4_box_Sfe5_Sfe6.o              \
          ddbar_n1n4_self_Sfe5.o ddbar_n1n4_self_Sfe5_Sfe6.o            \
          ddbar_n1n4_self_Sfe5_Sfe6_Sfe7.o                              \
          ddbar_n1n4_SquaredME.o ddbar_n1n4_vert_Sfe5_Sfe6.o            \
          ddbar_n1n4_vert.o ddbar_n1n4_vert_Sfe5.o
VIRTUALS += uubar_n2n3_abbr0h.o uubar_n2n3_abbr0s.o                     \
          uubar_n2n3_abbr1a.o uubar_n2n3_abbr1h.o uubar_n2n3_abbr1s.o   \
          uubar_n2n3_born.o uubar_n2n3_born_Sfe5.o                      \
          uubar_n2n3_box_Sfe5.o uubar_n2n3_box_Sfe5_Sfe6.o              \
          uubar_n2n3_self_Sfe5.o uubar_n2n3_self_Sfe5_Sfe6.o            \
          uubar_n2n3_self_Sfe5_Sfe6_Sfe7.o                              \
          uubar_n2n3_SquaredME.o uubar_n2n3_vert_Sfe5_Sfe6.o            \
          uubar_n2n3_vert.o uubar_n2n3_vert_Sfe5.o                      \
          ddbar_n2n3_abbr0h.o ddbar_n2n3_abbr0s.o                       \
          ddbar_n2n3_abbr1a.o ddbar_n2n3_abbr1h.o ddbar_n2n3_abbr1s.o   \
          ddbar_n2n3_born.o ddbar_n2n3_born_Sfe5.o                      \
          ddbar_n2n3_box_Sfe5.o ddbar_n2n3_box_Sfe5_Sfe6.o              \
          ddbar_n2n3_self_Sfe5.o ddbar_n2n3_self_Sfe5_Sfe6.o            \
          ddbar_n2n3_self_Sfe5_Sfe6_Sfe7.o                              \
          ddbar_n2n3_SquaredME.o ddbar_n2n3_vert_Sfe5_Sfe6.o            \
          ddbar_n2n3_vert.o ddbar_n2n3_vert_Sfe5.o
VIRTUALS += uubar_n2n4_abbr0h.o uubar_n2n4_abbr0s.o                     \
          uubar_n2n4_abbr1a.o uubar_n2n4_abbr1h.o uubar_n2n4_abbr1s.o   \
          uubar_n2n4_born.o uubar_n2n4_born_Sfe5.o                      \
          uubar_n2n4_box_Sfe5.o uubar_n2n4_box_Sfe5_Sfe6.o              \
          uubar_n2n4_self_Sfe5.o uubar_n2n4_self_Sfe5_Sfe6.o            \
          uubar_n2n4_self_Sfe5_Sfe6_Sfe7.o                              \
          uubar_n2n4_SquaredME.o uubar_n2n4_vert_Sfe5_Sfe6.o            \
          uubar_n2n4_vert.o uubar_n2n4_vert_Sfe5.o                      \
          ddbar_n2n4_abbr0h.o ddbar_n2n4_abbr0s.o                       \
          ddbar_n2n4_abbr1a.o ddbar_n2n4_abbr1h.o ddbar_n2n4_abbr1s.o   \
          ddbar_n2n4_born.o ddbar_n2n4_born_Sfe5.o                      \
          ddbar_n2n4_box_Sfe5.o ddbar_n2n4_box_Sfe5_Sfe6.o              \
          ddbar_n2n4_self_Sfe5.o ddbar_n2n4_self_Sfe5_Sfe6.o            \
          ddbar_n2n4_self_Sfe5_Sfe6_Sfe7.o                              \
          ddbar_n2n4_SquaredME.o ddbar_n2n4_vert_Sfe5_Sfe6.o            \
          ddbar_n2n4_vert.o ddbar_n2n4_vert_Sfe5.o
VIRTUALS += uubar_n3n3_abbr0h.o uubar_n3n3_abbr0s.o                     \
          uubar_n3n3_abbr1a.o uubar_n3n3_abbr1h.o uubar_n3n3_abbr1s.o   \
          uubar_n3n3_born.o uubar_n3n3_born_Sfe5.o                      \
          uubar_n3n3_box_Sfe5.o uubar_n3n3_box_Sfe5_Sfe6.o              \
          uubar_n3n3_self_Sfe5.o uubar_n3n3_self_Sfe5_Sfe6.o            \
          uubar_n3n3_self_Sfe5_Sfe6_Sfe7.o                              \
          uubar_n3n3_SquaredME.o uubar_n3n3_vert_Sfe5_Sfe6.o            \
          uubar_n3n3_vert.o uubar_n3n3_vert_Sfe5.o                      \
          ddbar_n3n3_abbr0h.o ddbar_n3n3_abbr0s.o                       \
          ddbar_n3n3_abbr1a.o ddbar_n3n3_abbr1h.o ddbar_n3n3_abbr1s.o   \
          ddbar_n3n3_born.o ddbar_n3n3_born_Sfe5.o                      \
          ddbar_n3n3_box_Sfe5.o ddbar_n3n3_box_Sfe5_Sfe6.o              \
          ddbar_n3n3_self_Sfe5.o ddbar_n3n3_self_Sfe5_Sfe6.o            \
          ddbar_n3n3_self_Sfe5_Sfe6_Sfe7.o                              \
          ddbar_n3n3_SquaredME.o ddbar_n3n3_vert_Sfe5_Sfe6.o            \
          ddbar_n3n3_vert.o ddbar_n3n3_vert_Sfe5.o
VIRTUALS += uubar_n3n4_abbr0h.o uubar_n3n4_abbr0s.o                     \
          uubar_n3n4_abbr1a.o uubar_n3n4_abbr1h.o uubar_n3n4_abbr1s.o   \
          uubar_n3n4_born.o uubar_n3n4_born_Sfe5.o                      \
          uubar_n3n4_box_Sfe5.o uubar_n3n4_box_Sfe5_Sfe6.o              \
          uubar_n3n4_self_Sfe5.o uubar_n3n4_self_Sfe5_Sfe6.o            \
          uubar_n3n4_self_Sfe5_Sfe6_Sfe7.o                              \
          uubar_n3n4_SquaredME.o uubar_n3n4_vert_Sfe5_Sfe6.o            \
          uubar_n3n4_vert.o uubar_n3n4_vert_Sfe5.o                      \
          ddbar_n3n4_abbr0h.o ddbar_n3n4_abbr0s.o                       \
          ddbar_n3n4_abbr1a.o ddbar_n3n4_abbr1h.o ddbar_n3n4_abbr1s.o   \
          ddbar_n3n4_born.o ddbar_n3n4_born_Sfe5.o                      \
          ddbar_n3n4_box_Sfe5.o ddbar_n3n4_box_Sfe5_Sfe6.o              \
          ddbar_n3n4_self_Sfe5.o ddbar_n3n4_self_Sfe5_Sfe6.o            \
          ddbar_n3n4_self_Sfe5_Sfe6_Sfe7.o                              \
          ddbar_n3n4_SquaredME.o ddbar_n3n4_vert_Sfe5_Sfe6.o            \
          ddbar_n3n4_vert.o ddbar_n3n4_vert_Sfe5.o
VIRTUALS += uubar_n4n4_abbr0h.o uubar_n4n4_abbr0s.o                     \
          uubar_n4n4_abbr1a.o uubar_n4n4_abbr1h.o uubar_n4n4_abbr1s.o   \
          uubar_n4n4_born.o uubar_n4n4_born_Sfe5.o                      \
          uubar_n4n4_box_Sfe5.o uubar_n4n4_box_Sfe5_Sfe6.o              \
          uubar_n4n4_self_Sfe5.o uubar_n4n4_self_Sfe5_Sfe6.o            \
          uubar_n4n4_self_Sfe5_Sfe6_Sfe7.o                              \
          uubar_n4n4_SquaredME.o uubar_n4n4_vert_Sfe5_Sfe6.o            \
          uubar_n4n4_vert.o uubar_n4n4_vert_Sfe5.o                      \
          ddbar_n4n4_abbr0h.o ddbar_n4n4_abbr0s.o                       \
          ddbar_n4n4_abbr1a.o ddbar_n4n4_abbr1h.o ddbar_n4n4_abbr1s.o   \
          ddbar_n4n4_born.o ddbar_n4n4_born_Sfe5.o                      \
          ddbar_n4n4_box_Sfe5.o ddbar_n4n4_box_Sfe5_Sfe6.o              \
          ddbar_n4n4_self_Sfe5.o ddbar_n4n4_self_Sfe5_Sfe6.o            \
          ddbar_n4n4_self_Sfe5_Sfe6_Sfe7.o                              \
          ddbar_n4n4_SquaredME.o ddbar_n4n4_vert_Sfe5_Sfe6.o            \
          ddbar_n4n4_vert.o ddbar_n4n4_vert_Sfe5.o

SOURCESMAIN = $(PWHG) $(PDFPACK) $(USER) $(PWHGANALYSIS) $(FUNC) $(VIRTUALS)  \
	$(FORMCALC) $(MMADOBJ) $(MADOBJ)

SOURCESLHEF = $(LHEF) $(PWHGANALYSIS) $(MMADOBJ)

SOURCESPYTHIA = $(PYTHIA) $(PWHGANALYSIS) $(MMADOBJ)

SOURCESHERWIG = $(HERWIG) $(PWHGANALYSIS) $(MMADOBJ)
          
# Additional FormCalc Files (maybe custom ones)
#VIRTUALS += BornEps_n1x1.o

########################################################################
## Libraries
### LibSLHA in Tools/SLHALib-2.2 ###
LIBS += $(TOOLS)/libSLHA.a

### LoopTools in Tools/Looptools-2.12 ###
LIBS += $(TOOLS)/libooptools.a

### DHELAS in Tools/DHELAS ###
LIBS += $(TOOLS)/libdhelas3.a

#not necessary anymore
#LIBS += $(LIBGFORTRANPATH)/libgfortranbegin.a
#LIBS += $(LIBGFORTRANPATH)/libgfortran.a
#LIBS += $(LIBSTDCPPPATH)/libstdc++.a

### Libraries for linking c++ code
LIBS += -ldl -lstdc++

########################################################################
## combine all flags, libraries and includes

ALL_FCFLAGS   = $(REC_FCFLAGS) $(OPT) $(WARN) $(CPP) $(DEBUG) $(USRFLAGS) $(CHANNEL)
ALL_FCFLAGS  += -I$(UINCLUDE) -I$(PINCLUDE) -I$(LTINCLUDE) -I$(MDINCLUDE)
ALL_FCFLAGS  += -I$(WORKINGDIR) -I$(VINCLUDE) -I$(SLHAINCLUDE)
ALL_FCFLAGS  += -I$(SINCLUDE)

ALL_CFLAGS    = $(REC_CFLAGS) $(OPT)
ALL_CFLAGS   += -I$(LTINCLUDE)

ALL_CXXFLAGS  = $(REC_CXXFLAGS) $(OPT)

LINKER = $(CPPFLAGS) $(LIBS) $(LDFLAGS)

HEADERS  = $(wildcard *.h $(UINCLUDE)/*.h  $(PINCLUDE)/*.h $(VINCLUDE)/*.h)
HEADERS += $(wildcard *.h $(LTINCLUDE)/*.h $(SLHAINCLUDE)/*.h)
HEADERS += $(wildcard *.h $(MDINCLUDE)/*.h $(SINCLUDE)/*.h)
HEADERS += $(wildcard *.inc $(MDINCLUDE)/*.inc $(UINCLUDE)/*.inc)
     
########################################################################
## Rules, generate objects

%.o: %.f $(HEADERS)
	@echo "Compiling:" $<
	@$(FC) $(ALL_FCFLAGS) -c -o build/$@ $<

%.o: %.F $(HEADERS)
	@echo "Compiling:" $<
	@$(FC) $(ALL_FCFLAGS) -c -o build/$@ $<

pythiafortran.o: pythiafortran.cc
	@echo "Compiling:" $<
	@$(CXX) $(ALL_CXXFLAGS) $(PYCXXFLAGS) -c -o build/$@ $^

fastjetfortran.o: fastjetfortran.cc
	@echo "Compiling:" $<
	@$(CXX) $(ALL_CXXFLAGS) $(FJCXXFLAGS) -c -o build/$@ $^

%.o: %.c $(HEADERS)
	@echo "Compiling:" $<
	@$(CC) $(ALL_CFLAGS) -c -o build/$@ $<

########################################################################
## Rules, link
## type make -j4 [rule] to speed up the compilation

pwhg_main_ninj: $(SOURCESMAIN)
	$(FC) $(ALL_FCFLAGS) $(patsubst %,build/%,$(SOURCESMAIN)) $(LINKER) -o $@

lhef_analysis:$(SOURCESLHEF)
	$(FC) $(ALL_FCFLAGS) $(patsubst %,build/%,$(SOURCESLHEF)) $(LINKER) -o $@ 

main-PYTHIA-lhef: $(SOURCESPYTHIA)
	$(FC) $(ALL_FCFLAGS) $(patsubst %,build/%,$(SOURCESPYTHIA)) $(LINKER) -o $@

# not implemented yet
main-HERWIG-lhef: $(SOURCESHERWIG)
	$(FC) $(ALL_FCFLAGS) $(patsubst %,build/%,$(SOURCESHERWIG)) $(LINKER) -o $@

merge-pwg-stat: 
	$(FC) $(TOOLS)/programs/merge-pwg-stat.f -o $@

merge-data: 
	$(FC) $(TOOLS)/programs/merge-data.f -o $@

copy-scripts:
	cp $(SCRIPTS)/runparallel.sh  $(SUBWORKINGDIR) && \
	cp $(SCRIPTS)/clean_rundir.sh $(SUBWORKINGDIR) && \
	cp $(SCRIPTS)/reweight.sh     $(SUBWORKINGDIR) && \
	cp $(SCRIPTS)/merge.sh        $(SUBWORKINGDIR)

all: pwhg_main_ninj copy-scripts merge-pwg-stat merge-data lhef_analysis main-PYTHIA-lhef

clean:
	rm -f build/*.o pwhg_main_nixj lhef_analysis main-PYTHIA-lhef merge-pwg-stat merge-data runparallel.sh clean_rundir.sh reweight.sh merge.sh

clean-results:
	find ./testrun ! \( -name '*.slha' -o -name '*.input' -o -name 'pwgseeds.dat' \) \
	-type f -exec rm -f {} +
	
clean-all: clean clean-results
	rm -f pwhg_main_ninj lhef_analysis main-PYTHIA-lhef

do-pwhg do: pwhg_main_ninj copy-scripts merge-pwg-stat merge-data
	cd testrun && ../pwhg_main_ninj

do-lhef: lhef_analysis
	cd testrun && ../lhef_analysis
	
do-pythia: main-PYTHIA-lhef
	cd testrun && ../main-PYTHIA-lhef
	
do-all: do-pwhg do-lhef do-pythia
	
########################################################################
#                       -*- End of Makefile -*-                        #
########################################################################
