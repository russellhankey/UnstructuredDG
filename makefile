FCOMP    = gfortran
OPTS     = -c -fbounds-check
LINKOPTS = -O3 -o
OBJS = main.o boundcond.o  cfl.o compflux.o connect.o detector.o  \
getrusanovflux.o read_data.o init.o preprocess.o  setup.o tec.o  calright.o caljacob.o tecpost.o tecplotter.o teccccc.o


rundg:$(OBJS)
	$(FCOMP) $(LINKOPTS) ./rundir/dg_fortran.exe $(OBJS)
	
setup.o:setup.f90
	$(FCOMP) $(OPTS) setup.f90

boundcond.o:setup.o boundcond.f90
	$(FCOMP) $(OPTS) boundcond.f90
	
cfl.o:setup.o cfl.f90
	$(FCOMP) $(OPTS) cfl.f90

compflux.o:setup.o compflux.f90
	$(FCOMP) $(OPTS) compflux.f90

getrusanovflux.o:setup.o getrusanovflux.f90
	$(FCOMP) $(OPTS) getrusanovflux.f90

read_data.o:read_data.f90
	$(FCOMP) $(OPTS) read_data.f90
	
init.o:setup.o init.f90
	$(FCOMP) $(OPTS) init.f90
	
main.o:setup.o main.f90
	$(FCOMP) $(OPTS) main.f90
	
detector.o:setup.o detector.f90
	$(FCOMP) $(OPTS) detector.f90

connect.o:connect.f90
	$(FCOMP) $(OPTS) connect.f90

preprocess.o:preprocess.f90
	$(FCOMP) $(OPTS) preprocess.f90

tec.o:setup.o tec.f90
	$(FCOMP) $(OPTS) tec.f90
	
teccccc.o:setup.o tec.f90
	$(FCOMP) $(OPTS) teccccc.f90
	
tecplotter.o:setup.o tecplotter.f90
	$(FCOMP) $(OPTS) tecplotter.f90
	
calright.o:setup.o calright.f90
	$(FCOMP) $(OPTS) calright.f90
	
caljacob.o:setup.o caljacob.f90 
	$(FCOMP) $(OPTS) caljacob.f90
	
tecpost.o:tecpost.f90
	$(FCOMP) $(OPTS) tecpost.f90


clean: 
	rm *.o  *.mod ./rundir/dg_fortran.exe

pltclean:
	rm ./rundir/*.out ./rundir/testec*
