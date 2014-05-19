PROGRAM main
	USE MatSparse
	IMPLICIT NONE
	TYPE(matrix)		:: X1	
	CHARACTER(LEN=80)		::	Avector
	CHARACTER(LEN=80)		::	bvector
	CHARACTER(LEN=80)		::	xvector
	INTEGER					::	nargs !number of command line arguments
	INTEGER					::	maxit, maxiter
	REAL					::	tol, omega
	
	!// Count command line arguments
	nargs = COMMAND_ARGUMENT_COUNT()
	!// Check if enough filenames are provided
	IF (nargs/=3) THEN
		!// Print error message and stop the program.
		PRINT*, 'Please provide filenames for A-vector, b-vector and x-vector'
		STOP
	END IF
	
	!// Get the command line arguments and store filenames to A-vector, b-vector and x-vector
	CALL GET_COMMAND_ARGUMENT(1,Avector)
	CALL GET_COMMAND_ARGUMENT(2,bvector)
	CALL GET_COMMAND_ARGUMENT(3,xvector)
	
	!// Text displayed in terminal for user to set maximum number of Jacobi iterations
	!// maxit, and successive over relaxation iterations, maxiter, the tolerance
	!// convergence limit, tol and the floating point number omega.
	PRINT*, 'Set an integer value for the maximum number of Jacobi iterations, maxit='
	READ(*,*) maxit
	PRINT*, 'Set an integer value for the maximum number of SOR iterations, maxiter='
	READ(*,*) maxiter
	PRINT*, 'Give a floating point number for the tolerance convergence limit, tol='
	READ(*,*) tol
	PRINT*, 'Give a floating point number for the constant omega:'
	READ(*,*) omega
	
	CALL X1%scan(Avector, bvector)
		
	CALL X1%Jacobi(tol, maxit)
	
	CALL X1%SOR(tol, omega, maxiter)
	
	CALL X1%dump(xvector, tol, omega, maxit, maxiter)
	
	PRINT*, 'Whoop! hele gjennom!'
	STOP
END PROGRAM main