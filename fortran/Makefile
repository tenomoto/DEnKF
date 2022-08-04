.SUFFIXES :
.SUFFIXES : .f90 .o

FC = /usr/local/bin/gfortran10
FFLAGS = -O2
OBJS = fd_module.o mg_module.o ode_module.o qg_module.o
SRCS = $(OBJS:.o=.f90)

#TARGET = test_mg
TARGET = run_qg

all : $(TARGET)

$(TARGET) : $(OBJS) $(TARGET).o
	$(FC) $> -o $@

clean :
	rm -f *.o *.mod $(TARGET)

.f90.o :
	$(FC) $(FFLAGS) $< -c
