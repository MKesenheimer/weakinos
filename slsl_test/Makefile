########################################################################
#                          -*- Makefile -*-                            #
########################################################################

########################################################################
## Flags

## Compiler and Compiler Flags
# TODO: use this convention:
#  CC          C compiler command
#  CFLAGS      C compiler flags
#  LDFLAGS     linker flags, e.g. -L<lib dir> if you have libraries in a
#              nonstandard directory <lib dir>
#  LIBS        libraries to pass to the linker, e.g. -l<library>
#  CPPFLAGS    C/C++/Objective C preprocessor flags, e.g. -I<include dir> if
#              you have headers in a nonstandard directory <include dir>
#  CPP         C preprocessor
#  CXX         C++ compiler command
#  CXXFLAGS    C++ compiler flags
#  CXXCPP      C++ preprocessor
#  FC          Fortran compiler command
#  FCFLAGS     Fortran compiler flags

FCOMP    =  gfortran
FFLAGS   = -fno-automatic -fno-range-check
FFLAGS  += -ffixed-line-length-none -lgfortran -ff2c -DU77EXT=0 -DQuad=0
CCOMP    =  gcc
CFLAGS   = -g -fomit-frame-pointer -ffast-math -Wall -m64
CFLAGS  += -stdlib=libstdc++ -mmacosx-version-min=10.6 -Qunused-arguments
CFLAGS  += -DNOUNDERSCORE=0 -DBIGENDIAN=0
CXXCOMP  =  g++
CXXFLAGS = -g -fomit-frame-pointer -ffast-math -Wall -m64
CXXFLAGS = -stdlib=libstdc++ -mmacosx-version-min=10.6 -Qunused-arguments

## PDF
## choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
## wheter to link precompiled LHAPDF library statically or compile it 
## from source (static, none)
PDF = lhapdf
PDFSTATIC = none

## choose analysis: none, default
ANALYSIS = default

## warning for type-conversions -> basically useless, as those occur in
## too many places
#WARN  = -Wconversion -Wall -Wtabs -Wall -Wimplicit-interface
## -fbounds-check sometimes causes a weird error due to non-lazy
## evaluation of boolean in gfortran.
#WARN += -fbounds-check
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT = -O2


### generate directory build and run, if not yet existing
$(shell mkdir -p build)
$(shell mkdir -p testrun)


########################################################################
## Runtime flags

## Preprocessor
# it might be advisable to use the -ffree-line-length-none 
# or -ffixed-line-length-none options
PREPRO = -cpp

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
#FLAGS += -D CHECK_FORM_MAD

## Additional and mandatory flags
# Use MadGraph Borns (default): -D MAD_BORN 
# Use FormCalc Borns: -D FORM_BORN
# Note: Use always MadGraph Borns to check for IR-finiteness
FLAGS += -D MAD_BORN

## Treatment of quark masses and squark widths
# The current version uses vanishing masses for the incoming quarks and
# zero widths for the internal squarks.
# For later software versions the possibility exists to make the
# calculations with finite b-quark mass or finite squark widths.
FLAGS += -D NEGLECTBMASS
FLAGS += -D NEGLECT_SWIDTH
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
# - Diagram Subtraction, implemented with preprocessor flag -D DSUB
FLAGS += -D DSUB

########################################################################
## Paths

WORKINGDIR = $(shell pwd)/..
PWHGDIR    = $(WORKINGDIR)/..
SUBWORKINGDIR = $(shell pwd)
#Linux
#LIBGFORTRANPATH = /usr/lib/gcc/x86_64-linux-gnu/4.8
#LIBSTDCPPPATH   = /usr/lib/gcc/x86_64-linux-gnu/4.8
#Mac OSX
LIBGFORTRANPATH = /usr/local/lib
LIBSTDCPPPATH   = /usr/local/lib

# directories
TOOLS   = $(WORKINGDIR)/Tools
REAL    = $(SUBWORKINGDIR)/MadGraph_Reals
VIRT    = $(SUBWORKINGDIR)/FormCalc_Virtuals
VIRTSUB = $(VIRT)/uubar_squaredme $(VIRT)/ddbar_squaredme
REALSUB = $(REAL)/SubProcesses


FJ     = $(TOOLS)/fastjet
PY     = $(TOOLS)/pythia
HER    = $(TOOLS)/herwig
LT     = $(TOOLS)/LoopTools-2.12
SLHA   = $(TOOLS)/SLHALib-2.2
DHELAS = $(TOOLS)/DHELAS
MADMODEL  = $(REAL)/MODEL

# includes
UINCLUDE  = $(WORKINGDIR)/include
PINCLUDE  = $(PWHGDIR)/include
VINCLUDE  = $(VIRT)/include
MDINCLUDE = $(REAL)/include
SLHAINCLUDE = $(SLHA)/src
LTINCLUDE = $(LT)/src/include

########################################################################
## search for the files and set paths

vpath %.f $(PWHGDIR) $(WORKINGDIR) $(TOOLS)/functions $(DHELAS) $(PY) $(HER)
vpath %.F $(WORKINGDIR) $(TOOLS)/FormCalcTools $(DHELAS)
vpath %.F $(VIRT)/RenConsts $(VIRTSUB)
vpath %.F $(MADMODEL) $(REALSUB)
vpath %.f $(MADMODEL) $(REALSUB)
vpath %.cc $(FJ)
vpath %.o $(SUBWORKINGDIR)/build

########################################################################
## Includes

FF  = $(FCOMP) $(FFLAGS) $(OPT) $(WARN) $(DEBUG) $(FLAGS) $(PREPRO) -I$(UINCLUDE) 
FF += -I$(PINCLUDE) -I$(LTINCLUDE) -I$(WORKINGDIR) -I$(VINCLUDE) -I$(SLHAINCLUDE)
FF += -I$(MDINCLUDE)
CC  = $(CCOMP) $(CFLAGS) $(OPT) -I$(LTINCLUDE)
CXX = $(CXXCOMP) $(CXXFLAGS) $(OPT)
INCLUDE  = $(wildcard *.h $(UINCLUDE)/*.h  $(PINCLUDE)/*.h $(VINCLUDE)/*.h)
INCLUDE += $(wildcard *.h $(UINCLUDE)/*.inc $(LTINCLUDE) $(SLHAINCLUDE)/*.h)
INCLUDE += $(wildcard *.h $(MDINCLUDE)/*.h)
INCLUDE += $(wildcard *.inc $(MDINCLUDE)/*.inc)

########################################################################
## Libraries
#LIBS  += $(LIBGFORTRANPATH)/libgfortranbegin.a
LIBS  += $(LIBGFORTRANPATH)/libgfortran.a
LIBS  += $(LIBSTDCPPPATH)/libstdc++.a

### LibSLHA in Tools/SLHALib-2.2 ###
LIBS   += $(TOOLS)/libSLHA.a


### LoopTools in Tools/Looptools-2.12 ###
LIBS  += $(TOOLS)/libooptools.a


### DHELAS in Tools/DHELAS ###
LIBS  += $(TOOLS)/libdhelas3.a

########################################################################
## Linker files

### USER Files ###
USER = Born.o Born_phsp.o Born_lh.o Reals.o Reals_lh.o Virtuals.o	\
       init_processes_slsl.o init_couplings.o

       
### Functions in Tools/functions/
FUNC = hep2pwg.o funcprocess.o Print.o funcbasic.o off_to_on.o


### MadGraph Amplitudes
MADOBJ = couplings.o hardstop.o printout.o read_slha.o set_general_couplings.o

# MadGraph born amplitudes
MADOBJ += b_cl_001_001.o b_sf_003_003.o realmtrx_008.o sborn_cl_002.o \
	b_cl_001_002.o b_sf_003_004.o realmtrx_009.o sborn_cl_003.o \
	b_cl_002_001.o b_sf_004_003.o realmtrx_010.o sborn_cl_004.o \
	b_cl_002_002.o b_sf_004_004.o realmtrx_011.o sborn_cl_005.o \
	b_cl_003_001.o b_sf_005_003.o realmtrx_012.o sborn_cl_006.o \
	b_cl_003_002.o b_sf_005_004.o realmtrx_013.o sborn_cl_007.o \
	b_cl_004_001.o b_sf_006_003.o realmtrx_014.o sborn_cl_008.o \
	b_cl_004_002.o b_sf_006_004.o realmtrx_015.o sborn_cl_009.o \
	b_cl_005_001.o b_sf_007_003.o realmtrx_016.o sborn_cl_010.o \
	b_cl_005_002.o b_sf_007_004.o realmtrx_017.o sborn_proc.o   \
	b_cl_006_001.o b_sf_008_003.o realmtrx_018.o sborn_sf_001.o \
	b_cl_006_002.o b_sf_008_004.o realmtrx_019.o sborn_sf_002.o \
	b_cl_007_001.o b_sf_009_003.o realmtrx_020.o sborn_sf_003.o \
	b_cl_007_002.o b_sf_009_004.o realmtrx_021.o sborn_sf_004.o \
	b_cl_008_001.o b_sf_010_003.o realmtrx_022.o sborn_sf_005.o \
	b_cl_008_002.o b_sf_010_004.o realmtrx_023.o sborn_sf_006.o \
	b_cl_009_001.o realmtrx_024.o sborn_sf_007.o \
	b_cl_009_002.o realmtrx_001.o realmtrx_025.o sborn_sf_008.o \
	b_cl_010_001.o realmtrx_002.o realmtrx_026.o sborn_sf_009.o \
	b_cl_010_002.o realmtrx_003.o realmtrx_027.o sborn_sf_010.o \
	b_sf_001_003.o realmtrx_004.o realmtrx_028.o sreal_proc.o  \
	b_sf_001_004.o realmtrx_005.o realmtrx_029.o \
	b_sf_002_003.o realmtrx_006.o realmtrx_030.o \
	b_sf_002_004.o realmtrx_007.o sborn_cl_001.o

### PDFs ###
ifeq ("$(PDF)","lhapdf")
  LHAPDF_CONFIG = lhapdf-config
  PDFPACK       = lhapdfif.o
  LIBSLHAPDF    = -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)		\
                  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF
ifeq  ("$(PDFSTATIC)","static")
  ## if you want to provide a static lhapdf library, specify it here:
  LIBSLHAPDF += /opt/lib/libLHAPDF.a
endif
  LIBS += $(LIBSLHAPDF) 
else
  PDFPACK = mlmpdfif.o hvqpdfpho.o
endif


### ANALYSIS ###
ifeq ("$(ANALYSIS)","default")
  USE_FASTJET=yes
  PWHGANALYSIS = Analysis.o
endif
ifeq ("$(ANALYSIS)","none")
  PWHGANALYSIS = Analysis-dummy.o
endif

ifdef USE_FASTJET
  FASTJET_CONFIG=$(shell which fastjet-config)
  ifeq ($(FASTJET_CONFIG),)
    $(error cannot find fastjet-config (did you install fastjet?))
  endif
  LIBS += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
  FJCXXFLAGS  += $(shell $(FASTJET_CONFIG) --cxxflags)
  PWHGANALYSIS+= fastjetfortran.o
  # fastjetsisconewrap.o fastjetktwrap.o
endif

### Formcalc Files in Tools/FormcalcTools/###
FORMCALC = VecSet.o


### Virtual Amplitudes in Formcalc_Virtuals/ ###
VIRTUALS = CalcRenConst.o uubar_abbr0h.o uubar_abbr0s.o     	        \
          uubar_abbr1a.o uubar_abbr1h.o uubar_abbr1s.o                  \
          uubar_born.o uubar_box_Sfe5_Sfe6.o                            \
          uubar_SquaredME.o                                             \
          uubar_vert.o uubar_vert_Sfe5.o uubar_vert_Sfe5_Sfe6.o		\
          ddbar_abbr0h.o ddbar_abbr0s.o            			\
          ddbar_abbr1a.o ddbar_abbr1h.o ddbar_abbr1s.o                  \
          ddbar_born.o ddbar_box_Sfe5_Sfe6.o                            \
          ddbar_SquaredME.o                                             \
          ddbar_vert.o ddbar_vert_Sfe5.o ddbar_vert_Sfe5_Sfe6.o	

### PWHG FILES in ../ ###
PWHG = pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o utils.o    \
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o			\
	test_Sudakov.o pt2maxreg.o maxrat.o				\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o bra_ket_subroutines.o			\
	init_phys.o powheginput.o pdfcalls.o sigcollremn.o		\
	pwhg_bookhist-multi.o pwhg_analysis_driver.o checkmomzero.o	\
	setstrongcoupl.o integrator.o newunit.o mwarn.o			\
	sigcollsoft.o validflav.o decidemerging.o          		\
        loadstorenlo.o setlocalscales.o ubprojections.o      		\
        mint_upb.o pwhgreweight.o reshufflemoms.o opencount.o		
PWHG += cernroutines.o find_regions.o gen_real_phsp.o sigborn.o		\
	sigreal.o sigsoftvirt.o sigvirtual.o


### LHEF ###
# analyze LHEF output
LHEF = lhef_analysis.o boostrot.o random.o cernroutines.o		\
     opencount.o powheginput.o pwhg_bookhist-multi.o 			\
     lhefread.o newunit.o pwhg_analysis_driver.o 			\
     funcprocess.o init_couplings.o Print.o				\
     bra_ket_subroutines.o off_to_on.o funcbasic.o


### PYTHIA ###
# target to read event file, shower events with PYTHIA + analysis
# NOTE: pythia-6.4.25-decay.f is based on POWHEG-BOX/pythia-6.4.25.f,
# to which we applied a simple patch to enable slepton decays
# (for clarity a separate patch file is included which is not needed otherwise)
PYTHIA = main-PYTHIA.o setup-PYTHIA-lhef.o pythia6428.o	\
	boostrot.o powheginput.o init_couplings.o funcprocess.o 	\
	pwhg_bookhist-multi.o lhefread.o newunit.o pdfdummies.o		\
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o	\
	bra_ket_subroutines.o off_to_on.o funcbasic.o Print.o

### HERWIG ###
# this is work in progress
HERWIG = main-HERWIG.o setup-HERWIG-lhef.o herwig6510.o boostrot.o 	\
     opencount.o powheginput.o pwhg_bookhist-multi.o random.o		\
     lhefread.o newunit.o pwhg_analysis_driver.o 			\
     funcprocess.o init_couplings.o cernroutines.o			\
     bra_ket_subroutines.o off_to_on.o funcbasic.o

########################################################################
## Generate objects
%.o: %.f $(INCLUDE)
	$(FF) -c -o build/$@ $<

%.o: %.F $(INCLUDE)
	$(FF) -c -o build/$@ $<

%.o: %.cc
	$(CXX) -c -o build/$@ $^ $(FJCXXFLAGS)

%.o: %.c $(INCLUDE)
	$(CC) -c -o build/$@ $<
	
########################################################################
## Rules
## type make -j4 [rule] to speed up the compilation
BUILDMAIN = $(PWHG) $(PDFPACK) $(USER) $(PWHGANALYSIS) $(FUNC) $(VIRTUALS)  \
	$(FORMCALC) $(LTOBJ) $(DHELASOBJ) $(MADOBJ)

BUILDLHEF = $(LHEF) $(PWHGANALYSIS) $(MADOBJ)

BUILDPYTHIA = $(PYTHIA) $(PWHGANALYSIS)  $(MADOBJ)

BUILDHERWIG = $(HERWIG) $(PWHGANALYSIS)  $(MADOBJ)

pwhg_main: $(BUILDMAIN)
	$(FF) $(patsubst %,build/%,$(BUILDMAIN)) $(LDFLAGS) $(LIBS) -o $@ 

lhef_analysis:$(BUILDLHEF)
	$(FF) $(patsubst %,build/%,$(BUILDLHEF)) $(LIBS) -o $@ 

main-PYTHIA-lhef: $(BUILDPYTHIA)
	$(FF) $(patsubst %,build/%,$(BUILDPYTHIA)) $(LIBS) -o $@

# not implemented yet
main-HERWIG-lhef: $(BUILDHERWIG)
	$(FF) $(patsubst %,build/%,$(BUILDHERWIG)) $(LIBS) -o $@

libdhelas3.a:
	cd $(DHELAS) && make FC="$(FF)" F77="$(FF)"
	
liblooptools.a:
	cd $(LT) && make FC="$(FF)" F77="$(FF)"

libSLHA.a:
	cd $(SLHA) && make FC="$(FF)" F77="$(FF)"

libs: libdhelas3.a liblooptools.a libSLHA.a

all: pwhg_main lhef_analysis main-PYTHIA-lhef

clean:
	rm -f build/*.o

clean-results:
	cd testrun && rm -f *.out *.lhe *.top *.gp *.eps *.dat FlavRegList	\
	bornequiv realequiv realequivregions pwhg_checklimits   	\
	virtequiv

clean-libs:
	cd $(TOOLS) && rm -f *.a
	cd $(LT) && make clean
	cd $(DHELAS) && make clean
	cd $(SLHA) && make clean
	
clean-all: clean clean-results
	rm -f pwhg_main lhef_analysis main-PYTHIA-lhef

do-pwhg do: pwhg_main
	cd testrun && ../pwhg_main

do-lhef: lhef_analysis
	cd testrun && ../lhef_analysis
	
do-pythia: main-PYTHIA-lhef
	cd testrun && ../main-PYTHIA-lhef
	
do-all: do-pwhg do-lhef do-pythia
	
########################################################################
#                       -*- End of Makefile -*-                        #
########################################################################
