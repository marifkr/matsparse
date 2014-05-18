#
# Makefile
#
# The compiler options
FC = gfortran
CC =
#
FFLAGS = #-g
#
LIBPATH =
#
INCLUDEPATH =
#
LIB =
#
PROG = kjor
# The object files
OBJECTS = MatSparse.o main.o
#
kjor : $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)
	./kjor A-ascii.txt b-vector-ascii.txt
MatSparse.o : MatSparse.f90
	$(FC) $(FFLAGS) -c $?

main.o : main.f90
	$(FC) $(FFLAGS) -c $?

clean :
	rm -f *.o
	rm -r *.mod
	rm -f kjor