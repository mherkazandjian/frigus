FFLAGS = -Wall -fcheck=all -g -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow
FC     = gfortran

target = level_population.f03

default:
	#$(FC) $(target) -llapack -L/usr/lib/lapack -lblas -L/usr/lib/libblas
	$(FC) $(FFLAGS) -c types.f03
	$(FC) $(FFLAGS) -c sort.f03
	$(FC) $(FFLAGS) -c energy_levels.f03
	$(FC) $(FFLAGS) -c radiation.f03
	$(FC) $(FFLAGS) -c collisions.f03
	$(FC) $(FFLAGS) -c matrix_construction.f03	
	$(FC) $(FFLAGS) $(target) types.o sort.o energy_levels.o radiation.o collisions.o matrix_construction.o #printauxil.f gamd.f90  
object:
	 $(FC) $(FFLAGS) -c types.f03
	 $(FC) $(FFLAGS) -c sort.f03
	 $(FC) $(FFLAGS) -c energy_levels.f03
	 $(FC) $(FFLAGS) -c radiation.f03
	 $(FC) $(FFLAGS) -c collisions.f03
	 $(FC) $(FFLAGS) -c matrix_construction.f03



clean:
	rm -fvr a.out *.o *.mod

