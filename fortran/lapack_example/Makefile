#FFLAGS = -Wall
FC     = gfortran

target = lapack_prb.f90

default:
	$(FC) $(target) -llapack -L/usr/lib/lapack -lblas -L/usr/lib/libblas


clean:
	rm -fvr a.out

