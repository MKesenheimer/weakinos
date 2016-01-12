# --- variables defined by configure ---

SRC = ./src
PREFIX = x86_64-Darwin
LIBDIRSUFFIX = 

EXE = 
DEF = -D
NOUNDERSCORE = 0

FC = gfortran
FFLAGS =  -O1 -g -ffixed-line-length-none     -fno-range-check -ff2c  $(DEF)QUAD=1 $(DEF)U77EXT=0

CC = clang
CFLAGS = -O3 -g -fomit-frame-pointer -ffast-math -Wall -m64   -DNOUNDERSCORE=$(NOUNDERSCORE) -DBIGENDIAN=0

CXX = clang++
CXXFLAGS = -O3 -g -fomit-frame-pointer -ffast-math -Wall -m64   

ML = 1
MCC = mcc
MCFLAGS = -st -b64
MCLIBS = -lpthread

LDFLAGS = -L/usr/local/lib/gcc/x86_64-apple-darwin14.4.0/5.1.0 -L/usr/local/lib/gcc/x86_64-apple-darwin14.4.0/5.1.0/../../.. -lgfortran -lSystem -lgcc_ext.10.5 -lgcc -lquadmath -lm -lgcc_ext.10.5 -lgcc -lSystem -L/usr/local/lib/gcc/x86_64-apple-darwin14.4.0/5.1.0 -L/usr/local/lib/gcc/x86_64-apple-darwin14.4.0/5.1.0/../../.. -lgfortran -lSystem -lgcc_ext.10.5 -lgcc -lquadmath -lm -lgcc_ext.10.5 -lgcc -lSystem -m64

AR = ar
RANLIB = ranlib
NM = nm
DLLTOOL = 

PARALLEL = -j 2

# --- end defs by configure ---


BLD = build$(QUADSUFFIX)

LIBDIR = $(PREFIX)/lib$(LIBDIRSUFFIX)
INCLUDEDIR = $(PREFIX)/include
BINDIR = $(PREFIX)/bin

LIB = libooptools$(QUADSUFFIX).a
FE = lt$(QUADSUFFIX)$(EXE)
MFE = LoopTools$(QUADSUFFIX)$(EXE)
INCLUDE = $(BLD)/looptools.h $(BLD)/clooptools.h

ARGS = $(PARALLEL) \
  LIB="$(LIB)" \
  FE="$(FE)" \
  MFE="$(MFE)" \
  EXE="$(EXE)" \
  DEF="$(DEF)" \
  NOUNDERSCORE="$(NOUNDERSCORE)" \
  XFC="$(FC$(QUADSUFFIX)) $(FFLAGS) $(FFLAGS-quad) -I." \
  F90="$(F90)" \
  CC="$(CC)" \
  CFLAGS="$(CFLAGS) $(CFLAGS-quad)" \
  CXX="$(CXX)" \
  CXXFLAGS="$(CXXFLAGS)" \
  ML="$(ML)" \
  MCC="$(MCC)" \
  MCFLAGS="$(MCFLAGS)" \
  MCLIBS="$(MCLIBS)" \
  AR="$(AR)" \
  RANLIB="$(RANLIB)" \
  NM="$(NM)" \
  DLLTOOL="$(DLLTOOL)" \
  LDFLAGS="$(LDFLAGS)" \
  LIBPATH="$(LIBPATH)" \
  OBJS-quad="$(OBJS-quad)"


default all lib frontend mma: force
	cd $(BLD) && $(MAKE) $(ARGS) $@
	cp -p $(BLD)/$(LIB) ../

install: lib frontend
	-mkdir $(PREFIX)
	-mkdir $(LIBDIR) $(BINDIR) $(INCLUDEDIR)
	cp -p $(BLD)/$(LIB) $(LIBDIR)
	cp -p $(BLD)/$(LIB) ../
	cp -p $(INCLUDE$(QUADSUFFIX)) $(INCLUDEDIR)
	strip $(BLD)/$(FE)
	cp -p $(BLD)/fcc $(BLD)/$(FE) $(BINDIR)
	rm -f $(BINDIR)/f++
	ln -s fcc $(BINDIR)/f++
	test ! -f $(BLD)/$(MFE) || { strip $(BLD)/$(MFE) ; cp -p $(BLD)/$(MFE) $(BINDIR); }

force: $(BLD)/timestamp

$(BLD)/timestamp:
	-mkdir $(BLD)
	find "`cd $(SRC) ; pwd`" -follow -exec ln -sf {} $(BLD) \;
	touch $(BLD)/timestamp

clean:
	rm -fr $(BLD)
	rm -fr $(QUADSUFFIX)

