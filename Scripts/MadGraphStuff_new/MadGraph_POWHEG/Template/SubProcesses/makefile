
#F77           = gfortran
FFLAGS        = $(XFFLAGS)
LIBDIR        = ../$(OBJ)
LIBRARY	      = $(LIBDIR)/libmadgraph.a
CODE          =  sborn_proc.o $(patsubst %.f,%.o,$(wildcard sborn_sf_???.f)) $(patsubst %.f,%.o,$(wildcard sborn_cl_???.f)) $(patsubst %.f,%.o,$(wildcard b_sf_???_???.f)) $(patsubst %.f,%.o,$(wildcard b_cl_???_???.f)) $(patsubst %.f,%.o,$(wildcard realmtrx_???.f)) sreal_proc.o nexper.o

.f.o: ; $(F77) $(FFLAGS) -c $*.f

all: $(LIBDIR)libmadgraph.a

$(LIBDIR)libmadgraph.a: $(CODE)
	ar cru libmadgraph.a $(CODE)
	ranlib libmadgraph.a
	mv libmadgraph.a  $(LIBDIR)

