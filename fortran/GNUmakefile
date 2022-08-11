.SUFFIXES :
.SUFFIXES : .f90 .o

FC = /opt/local/bin/gfortran
FFLAGS = -O2
MGOBJS = fd_module.o mg_module.o
QGOBJS = $(MGOBJS) ode_module.o qg_module.o
DENKFOBJS = la_module.o denkf_module.o
DENKFLDFLAGS = -L/opt/local/lib
DENKFLIBS = -lopenblas -llapack

MG = test_mg
QG = run_qg
DENKF = run_denkf
EXES = $(MG) $(QG) $(DENKF)

all : $(QG) $(DENKF)

$(MG) : $(MGOBJS) $(MG).o
	$(FC) $^ -o $@

$(QG) : $(QGOBJS) $(QG).o
	$(FC) $^ -o $@

$(DENKF) : $(DENKFOBJS) $(DENKF).o
	$(FC) $(DENKFLDFLAGS) $^ -o $@ $(DENKFLIBS)

clean :
	rm -f *.o *.mod $(EXES)

.f90.o :
	$(FC) $(FFLAGS) $< -c
