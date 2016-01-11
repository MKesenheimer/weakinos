FF = gfortran
FFLAGS = -O2 -ffixed-line-length-132 
OBJECTS= particles.o sort.o readproc.o tops.o fields.o create.o\
colors.o tsums.o fsums.o ftsums.o writeamp_powheg.o gencolor.o eigen_lapack.o\
lapack.o fermi4.o drawfeyn.o post.o permute.o vert3.o vert4.o \
qed.o qcd.o qfd.o readmodel.o nvert3.o nvert4.o writedsig.o crossing.o \
write_configs.o leshouche.o mainBorns.o

OBJECTS1= particles.o sort.o readproc.o tops.o fields.o create.o\
colors.o tsums.o fsums.o ftsums.o writeamp_powheg2.o gencolor.o eigen_lapack.o\
lapack.o fermi4.o drawfeyn.o post.o permute.o vert3.o vert4.o \
qed.o qcd.o qfd.o readmodel.o nvert3.o nvert4.o writedsig.o crossing.o \
write_configs.o leshouche.o mainReal.o

.f.o:; $(FF) $(FFLAGS) -c $*.f

all: mg2_clB mg2_real

mg2_clB:$(OBJECTS) makefile
	$(FF) $(FFLAGS) -o mg2_clB $(OBJECTS)

mg2_real:$(OBJECTS1) makefile
	$(FF) $(FFLAGS) -o mg2_real $(OBJECTS1)
