.SUFFIXES :
.SUFFIXES : .f90 .o

FC = /opt/local/bin/gfortran
FFLAGS = -O2
MGOBJS = fd_module.o mg_module.o
QGOBJS = $(MGOBJS) math_module.o ode_module.o qg_module.o
INVOBJS = la_module.o
GENOBSOBJS = math_module.o genobs.o
DENKFOBJS = math_module.o la_module.o denkf_module.o
LDFLAGS = -L/opt/local/lib
BLASLAPACK = -lopenblas -llapack

MG = test_mg
QG = run_qg
INV = test_inv
GENLOC = genloc
GENOBS = genobs
DENKF = run_denkf
EXES = $(QG) $(GENLOC) $(GENOBS) $(DENKF)

all : $(EXES)

$(MG) : $(MGOBJS) $(MG).o
	$(FC) $^ -o $@

$(QG) : $(QGOBJS) $(QG).o
	$(FC) $^ -o $@

$(INV) : $(INVOBJS) $(INV).o
	$(FC) $(LDFLAGS) $^ -o $@ $(BLASLAPACK)

$(GENLOC) : $(GENLOC).o
	$(FC) $< -o $@

$(GENOBS) : $(GENOBSOBJS) $(GENOBS).o
	$(FC) $(LDFLAGS) $^ -o $@

$(DENKF) : $(DENKFOBJS) $(DENKF).o
	$(FC) $(LDFLAGS) $^ -o $@ $(BLASLAPACK)

clean :
	rm -f *.o *.mod $(EXES)

.f90.o :
	$(FC) $(FFLAGS) $< -c
