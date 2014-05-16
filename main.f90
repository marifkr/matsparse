PROGRAM main
	USE MatSparse
	IMPLICIT NONE
	TYPE(matrix)		:: X1	
		CHARACTER(LEN=80)		::	Avector
		CHARACTER(LEN=80)		::	bvector
!		CHARACTER(LEN=80)		::	xvector
		INTEGER					::	nargs !number of command line arguments
		INTEGER					::	MAXIT, MAXITER
		REAL					::	TOL, omega
	
	!// Count command line arguments
	nargs = COMMAND_ARGUMENT_COUNT()
	!// Check if enough filenames are provided
	IF (nargs/=2) THEN
		!// Print error message and stop the program.
		PRINT*, 'Please provide filenames for Avector and bvector'
		STOP
	END IF
	
	!// Get the command line arguments and store filenames to Avector, bvector and xvector
	CALL GET_COMMAND_ARGUMENT(1,Avector)
	CALL GET_COMMAND_ARGUMENT(2,bvector)
	!CALL GET_COMMAND_ARGUMENT(3,xvector)
	
	!// Text displayed in terminal for user to set MAXIT, TOL and omega
	!PRINT*, 'Set an integer value for the maximum number of iterations, MAXIT='
	!READ(*,*) MAXIT
	!PRINT*, 'Give a floating point number for the tolerance convergence limit, TOL='
	!READ(*,*) TOL
	!PRINT*, 'Give a floating point number for the constant omega, omega='
	!READ(*,*) omega
	
	!X1 = matrix(AAvector,bbvector) !,xvector)
	CALL X1%scan(Avector, bvector)

END PROGRAM main