# This is an commentary line in a makefile
# Start of the makefile
fcomp = gfortran
#gfortran
filename = cf

all: 
	${fcomp} -o ${filename}.out ${filename}.f03
	./${filename}.out
# global.mod: global.o global.f90
# 	ifort -c global.f90
# global.o: global.f90
# 	ifort -c global.f90
# main.o: global.mod main.f90
# 	ifort -c main.f90
# function1.o: global.mod function1.f90
# 	ifort -c function1.f90
# subroutine1.o: subroutine1.f90
# 	ifort -c subroutine1.f90
clean:
	rm  ${filename}.out #main.o function1.o subroutine1.o execname
# End of the makefile