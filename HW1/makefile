F_EXT = f95

FC = gfortran

FCFLAGS = -O3
EXEC = a
# Object Files for build
a.out: \
set_precision.o \
set_constants.o \
set_inputs.o \
subroutines.o \
exact_q1d_type.o \
test_func.o
		$(FC) $(FCFLAGS) $^ -o $@

set_constants.o : set_precision.o
set_inputs.o : set_precision.o set_constants.o
subroutines.o : set_precision.o set_constants.o set_inputs.o
exact_q1d_type.o : set_precision.o set_constants.o set_inputs.o subroutines.o
test_func.o: set_precision.o set_constants.o set_inputs.o subroutines.o exact_q1d_type.o

%.o: %.$(F_EXT)
		$(FC) $(FCFLAGS) -o $@ -c $<

debug:
		$(eval FCFLAGS = -ggdb3 -00 -fbounds-check -Wall -Wextra\
		-fimplicit-none -pedantic -ffpe-trap=invalid,zero,overflow\
		-Wconversion-extra)

clean:
		rm -rf *.o *.mod
