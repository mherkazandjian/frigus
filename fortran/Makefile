FFLAGS_gnu = -Wall -fcheck=all -g -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow 
FFLAGS_intel = -extend_source 132 
FC     = gfortran

target = level_population.f03

default:
	#$(FC) $(target) -llapack -L/usr/lib/lapack -lblas -L/usr/lib/libblas
	$(FC) $(FFLAGS_gnu) -c types.f03
	$(FC) $(FFLAGS_gnu) -c sort.f03
	$(FC) $(FFLAGS_gnu) -c energy_levels.f03
	$(FC) $(FFLAGS_gnu) -c radiation.f03
	$(FC) $(FFLAGS_gnu) -c collisions.f03
	$(FC) $(FFLAGS_gnu) -c matrix_construction.f03	
	$(FC) $(FFLAGS_gnu) -c linear_algebra.f03		
	$(FC) $(FFLAGS_gnu) $(target) types.o sort.o energy_levels.o radiation.o collisions.o matrix_construction.o linear_algebra.o -llapack -L/usr/lib/lapack -lblas -L/usr/lib/libblas
object:
	 $(FC) $(FFLAGS_gnu) -c types.f03
	 $(FC) $(FFLAGS_gnu) -c sort.f03
	 $(FC) $(FFLAGS_gnu) -c energy_levels.f03
	 $(FC) $(FFLAGS_gnu) -c radiation.f03
	 $(FC) $(FFLAGS_gnu) -c collisions.f03
	 $(FC) $(FFLAGS_gnu) -c matrix_construction.f03
	 $(FC) $(FFLAGS_gnu) -c linear_algebra.f03		


clean:
	rm -fvr a.out *.o *.mod

