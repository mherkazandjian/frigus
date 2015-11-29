#FFLAGS = -Wall
FC     = gfortran

target = cf.f90

default:
	$(FC) $(target) -llapack -L/usr/lib/lapack -lblas -L/usr/lib/libblas


clean:
	rm -fvr a.out

