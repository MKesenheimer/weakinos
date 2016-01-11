F77   = f77
FFLAGS= -O -ffixed-line-length-132
PLOT    =  plot_events.o dbook.o rw_events.o kin_func.o set_cuts.o

.f.o: ; $(F77) $(FFLAGS) -c $*.f

plot_events: $(PLOT)
	$(F77) $(FFLAGS) -o plot_events $(PLOT)

restore:
	./restorelinks.sh
