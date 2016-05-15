#FFLAGS = -Wall
FC     = gfortran

target = level_population.f03

default:
	#$(FC) $(target) -llapack -L/usr/lib/lapack -lblas -L/usr/lib/libblas
	$(FC) -O -Wall -fcheck=all -g -fbacktrace $(target) types.o sort.o energy_levels.o radiation.o collisions.o matrix_construction.o printauxil.f gamd.f90  
object:
	 $(FC) -c types.f03
	 $(FC) -c sort.f03
	 $(FC) -c energy_levels.f03
	 $(FC) -c radiation.f03
	 $(FC) -c collisions.f03
	 $(FC) -c matrix_construction.f03



clean:
	rm -fvr a.out *.o *.mod

