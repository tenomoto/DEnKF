.SUFFIXES :
.SUFFIXES : .f90 .o

FC = /opt/local/bin/gfortran-mp-devel
FCFLAGS = -O0
OBJS = fd_module.o mg_module.o ode_module.o qg_module.o
SRCS = $(OBJS:.o=.f90)

#TARGET = test_mg
TARGET = run_qg

all : $(TARGET)

$(TARGET) : $(OBJS) $(TARGET).f90
	$(FC) -o $@ $^

.f90.o :
	$(FC) $(FFLAGS) $< -c
