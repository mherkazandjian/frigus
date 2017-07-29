#FFLAGS_gnu = -Og

FFLAGS_gnu = -Wall -fcheck=all -g -fbacktrace -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow
#,inexact,denormal

FFLAGS_intel = -extend_source 132 
FC     = gfortran
OBJ    = types.o sort.o energy_levels.o radiation.o collisions.o  testing_data.o read_data.o linear_algebra.o matrix_construction.o level_population.o
LIBS   = -llapack -L/usr/lib/lapack -lblas -L/usr/lib/libblas

target_h2 = cooling_H2.f03
target_hd = cooling_HD.f03

h2: $(OBJ)
	$(FC) $(FFLAGS_gnu) $(target_h2) $(OBJ) $(LIBS)
	
hd: $(OBJ)
	$(FC) $(FFLAGS_gnu) $(target_hd) $(OBJ) $(LIBS)
	
%.o: %.f03
	$(FC) $(FFLAGS_gnu) -c $^ -o $@
	
clean:
	rm -fvr a.out *.o *.mod
