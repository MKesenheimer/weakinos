# ----------------------------------------------------------------------------
#
# Makefile for DHELAS 3.0 library
# Feb. 28, 2001
#
# ----------------------------------------------------------------------------
#
# Use a TAB to precede shell commands (e.g., f90, ar, etc).
# Note: this Makefile uses features that *may not*
# be supported by make utilities other than GNU make.
#
# ----------------------------------------------------------------------------
# coupsm.F and coupsm.inc are not to be included if HELAS is used with 
# the 2003 version of MadEvent
# ----------------------------------------------------------------------------

#FC            = gfortran

FFLAGS        = $(XFFLAGS) -I.

LD	      = ld

DEST	      = ../

EXTHDRS	      =

HDRS	      =

INSTALL	      = /etc/install

LIBRARY	      = ../$(OBJ)/libdhelas3.a

MAKEFILE      = Makefile

OBJS	      =	httsxx.o hstlxx.o hvvshx.o \
		vvvtlx.o jvvtlx.o uvvvlx.o sstlxx.o usslxx.o \
		hvvhxx.o hvvvxx.o httaxx.o \
		vvsshx.o jvsshx.o ttssxx.o utssxx.o \
		ttsaxx.o jvvsxx.o vvvsxx.o jvshxx.o utsaxx.o \
		vvshxx.o jvtaxx.o vvtaxx.o uvvaxx.o \
		momntx.o mom2cx.o boostx.o rotxxx.o \
		ixxxxx.o oxxxxx.o vxxxxx.o sxxxxx.o \
		iovxxx.o fvixxx.o fvoxxx.o jioxxx.o j3xxxx.o \
		iosxxx.o fsixxx.o fsoxxx.o hioxxx.o \
		vvvxxx.o jvvxxx.o gggxxx.o jggxxx.o \
		vvsxxx.o jvsxxx.o hvvxxx.o \
		vssxxx.o jssxxx.o hvsxxx.o \
		sssxxx.o hssxxx.o \
		wwwwxx.o jwwwxx.o w3w3xx.o jw3wxx.o \
		ggggxx.o jgggxx.o \
		vvssxx.o jvssxx.o hvvsxx.o \
		ssssxx.o hsssxx.o \
		eaixxx.o eaoxxx.o jeexxx.o \
		ioscxx.o fsicxx.o fsocxx.o hiocxx.o \
		iovcxx.o fvicxx.o fvocxx.o jiocxx.o \
		iovgox.o fvigox.o fvogox.o jiogox.o \
		iovdmx.o fvidmx.o fvodmx.o jiodmx.o \
		txxxxx.o txxxx2.o \
		iotxkk.o iovtkk.o vvtxkk.o ftixkk.o ftoxkk.o \
                iovkxx.o vvvkxx.o jiokxx.o jvvkxx.o \
                sstxxx.o hstxxx.o ussxxx.o \
                iotxxx.o ftixxx.o ftoxxx.o uioxxx.o \
                vvtxxx.o uvvxxx.o jvtxxx.o \
                iovtxx.o fvtixx.o fvtoxx.o jiotxx.o uiovxx.o \
                vvvtxx.o jvvtxx.o uvvvxx.o \
                ggggtx.o jgggtx.o uggggx.o \
                w3w3nx.o jw3wnx.o jwwwnx.o wwwwnx.o \
		irxxxx.o orxxxx.o \
		iorsxx.o fsorxx.o hiorxx.o \
		irosxx.o fsirxx.o hiroxx.o \
		iorvxx.o fvorxx.o jiorxx.o \
		irovxx.o fvirxx.o jiroxx.o \
		iorvsx.o fvsorx.o jsiorx.o hviorx.o \
		irovsx.o fvsirx.o jsirox.o hvirox.o \
		iorvvx.o fvvorx.o jviorx.o \
		irovvx.o fvvirx.o jvirox.o \
		iorsgx.o irosgx.o fsorgx.o fsirgx.o hiorgx.o \
		hirogx.o \
		iorvgx.o fvorgx.o jiorgx.o \
                irovgx.o fvirgx.o jirogx.o \
                iorvvg.o fvvorg.o jviorg.o \
                irovvg.o fvvirg.o jvirog.o \
                iorscx.o fsorcx.o hiorcx.o \
                iroscx.o fsircx.o hirocx.o \
                iorvcx.o fvorcx.o jiorcx.o \
                irovcx.o fvircx.o jirocx.o \
                iorvsc.o fvsorc.o jsiorc.o hviorc.o \
                irovsc.o fvsirc.o jsiroc.o hviroc.o \
                iorvvc.o fvvorc.o jviorc.o \
                irovvc.o fvvirc.o jviroc.o 

PRINT	      = pr

SHELL	      = /bin/sh

SYSHDRS	      =

MFLAGS        = -e

#.F.o:
%.o : %.F
	$(FC) $(FFLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

all:		$(LIBRARY)

helas:
		-rm -f *.o
		export FFLAGS='-O +cpp' && \
		export LIBRARY=libdhelas3.a && \
		$(MAKE) $(MFLAGS)

helas_check:
		-rm -f *.o
		export FFLAGS='-O +cpp -DHELAS_CHECK' && \
		export LIBRARY=libdhelas3_check.a && \
		$(MAKE) $(MFLAGS)

install-helas:
		export LIBRARY=libdhelas3.a && \
		$(MAKE) $(MFLAGS) install	

install-helas_check:
		export LIBRARY=libdhelas3.a && \
		$(MAKE) $(MFLAGS) install	

$(LIBRARY):	$(OBJS)
		@echo  "Loading $(LIBRARY) ... "
		@ar cru $(LIBRARY) $(OBJS)
		@echo "done"

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(LIBRARY) core tags

install:	$(LIBRARY)
	        @echo Installing $(LIBRARY) in $(DEST)
	        @if [ $(DEST) != . ]; then \
	        (rm -f $(DEST)/$(LIBRARY); $(INSTALL) -f $(DEST) $(LIBRARY)); fi

$(DEST)/$(LIBRARY): $(SRCS) $(HDRS) $(EXTHDRS)
	        @$(MAKE) -f $(MAKEFILE) ROOT=$(ROOT) DEST=$(DEST) install
