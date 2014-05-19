!/////////////////////////////////////////////////////////////
!//
!// Module MatSparse containing subroutines for 
!// - reading number of rows, columns (n) and non-zero 
!//   values (nnz) from two files containing a matrix A,
!//   and a vector b, and reading A and b.
!//	- performing Jacobi and SOR iterations to find
!//   the x-vector in Ax=b.
!// - dumping the x-vector to a file for both Jacobi and
!//   SOR iterations.
!//
!////////////////////////////////////////////////////////////
MODULE MatSparse
	IMPLICIT NONE
		TYPE	::	matrix
			REAL(KIND=8), POINTER	::	A(:) => null(), b(:) => null()
			REAL(KIND=8), POINTER	::	xoldJac(:) => null(), xnewJac(:) => null()
			REAL(KIND=8), POINTER	::	xoldSOR(:) => null(), xnewSOR(:) => null()
			REAL					::	tol, omega
			INTEGER, POINTER		::	irow(:) => null(), jcol(:) => null()
			INTEGER					::	itsJac = 0, itsSOR = 0, n = 0, nnz = 0
			INTEGER					::	maxit, maxiter
		CONTAINS
			PROCEDURE				:: scan => scanfromfile
			PROCEDURE				:: dump => dumptofile
			PROCEDURE				:: Jacobi1_it => Jacobi1iteration
			PROCEDURE				:: SOR1_it => SOR1iteration
			PROCEDURE				:: Jacobi => Jacobi_it
			PROCEDURE				:: SOR => SOR_it
		END TYPE matrix
	CONTAINS
		SUBROUTINE scanfromfile(mat1,Avector,bvector)
			!/////////////////////////////////////////////////////////////
			!//
			!// scanfromfile reads files assigned in the main program
			!// to find n, nnz, A and b. The routine also provides irow
			!// and jcol to find positions of non-zero values in A.
			!// Space for xvectors needed in Jacobi and SOR iterations
			!// is allocated.
			!//
			!/////////////////////////////////////////////////////////////
			IMPLICIT NONE
			CLASS(matrix)					::	mat1
			CHARACTER(LEN=80), INTENT(IN)	::	Avector, bvector !filenames
			CHARACTER(LEN=80)				::	line
			INTEGER, PARAMETER				::	lun = 10, lyn = 11
			INTEGER							::	i, j, k, l, m, n, nnz, res 
			INTEGER							::	rowval = 0, rowvalold = 0
			!-----------------------------------------------------------------------------
			!// Confirm entry to module
			PRINT*, 'Got to module!'
		!---------------------------------------------------------------------------
			!// First b
			
			!// Open file containing the b-vector
			OPEN(UNIT=lun, FILE=bvector, FORM='FORMATTED', IOSTAT=res)
			!// Test if the file was opened correctly
			IF (res/= 0) THEN
				!// Print an error message
				PRINT*, 'Error in opening the b-vector file, status:', res
				!// Stop the program
				STOP
			END IF
			!// Read n from the first line in the file
			READ(UNIT=lun, FMT='(A)'), line
			IF(res/=0) THEN
				!// Print an error message
				PRINT*, 'Error in reading b-file, status:', res
				!// Close the file and stop the program
				CLOSE(UNIT=lun)
				STOP
			END IF
			!// Find position of :
			j=INDEX(line,':')
			!// Read everything after :
			READ(line(j+1:),*), mat1%n

			ALLOCATE(mat1%b(mat1%n))
			IF (res/=0) THEN
				!//Print error message
				PRINT*, 'Error in allocating b, status:', res
				!// Close file and stop program
				CLOSE(UNIT=lun)
				STOP
			END IF
			
			!// Read the data in to the b-vector
			DO i = 1,mat1%n
				READ(UNIT=lun, FMT='(A)',IOSTAT=res), line
				IF(res/=0) THEN
					PRINT*, 'Error in reading b-file, status:', res
					!// Close the file and stop the program
					CLOSE(UNIT=lun)
					STOP
				END IF
				!// Find position of =
				j=INDEX(line,'=')
				!// Read everything after =, and place in b-vector
				READ(line(j+1:),*), mat1%b(i)
			END DO
			!// Done reading file containing b-vector - close the file
			CLOSE(UNIT=lun)
			
		! -----------------------------------------------------------------------------
			!// Now A

			!// Open file containing the sparse A-matrix
			OPEN(UNIT=lyn, FILE=Avector, FORM='FORMATTED', IOSTAT=res)
			!// Test if the file was opened correctly
			IF (res/= 0) THEN
				!// Print an error message
				PRINT*, 'Error in opening A-file, status:', res
				!// Stop the program
				STOP
			END IF

			!// Read n and nnz from the first two lines in the file
			READ(UNIT=lyn, FMT='(A)'), line
			IF(res/=0) THEN
				!// Print an error message
				PRINT*, 'Error in reading A-file, status:', res
				!// Close the file and stop the program
				CLOSE(UNIT=lyn)
				STOP
			END IF
			!// Find position of :
			j=INDEX(line,':')
			!// Read everything after :
			READ(line(j+1:),*), mat1%n

			READ(UNIT=lyn, FMT='(A)'), line
			IF(res/=0) THEN
				!// Print an error message
				PRINT*, 'Error in reading A-file, status:', res
				!// Close the file and stop the program
				CLOSE(UNIT=lyn)
				STOP
			END IF
			!// Find position of :
			j=INDEX(line,':')
			!// Read everything after :
			READ(line(j+1:),*), mat1%nnz

			!// Allocating space for A, irow and jcol
			!// If allocation did not go well, close file and stop program
			ALLOCATE(mat1%A(mat1%nnz), STAT=res)
			IF (res/=0) THEN
				!//Print error message
				PRINT*, 'Error in allocating A, status:', res
				!// Close file and stop program
				CLOSE(UNIT=lyn)
				STOP
			END IF
			ALLOCATE(mat1%irow(mat1%n+1), STAT=res)
			IF (res/=0) THEN
				!//Print error message
				PRINT*, 'Error in allocating irow, status:', res
				!// Close file and stop program
				CLOSE(UNIT=lyn)
				STOP
			END IF
			ALLOCATE(mat1%jcol(mat1%nnz), STAT=res)
			IF (res/=0) THEN
				!//Print error message
				PRINT*, 'Error in allocating jcol, status:', res
				!// Close file and stop program
				CLOSE(UNIT=lyn)
				STOP
			END IF

			!// Read the data in to the A-vector
			DO i = 1,mat1%nnz
				READ(UNIT=lyn, FMT='(A)',IOSTAT=res), line
				IF(res/=0) THEN
					PRINT*, 'Error in reading A-file, status:', res
					!// Close the file and stop the program
					CLOSE(UNIT=lyn)
					STOP
				END IF
				!// Find position of =
				j=INDEX(line,'=')
				!// Read everything after =, and place in b-vector
				READ(line(j+1:),*), mat1%A(i)
	
				!// Now make jcol and irow with help from (row,column)
				
				!// Find positions of '(' and ',' and')'
				k=INDEX(line,',')
				l=INDEX(line,'(')
				m=INDEX(line,')')
	
				!// Read row and then column
				READ(line(l+1:k-1),*), rowval
				READ(line(k+1:m-1),*), mat1%jcol(i)
	
				!// Place number of value on new row in irow
				IF (rowval>rowvalold) THEN
					mat1%irow(rowval) = i
				END IF
				rowvalold = rowval	
			END DO
			
			!// Set the last value of irow
			mat1%irow(mat1%n+1) = mat1%nnz+1

			!// Done reading file containing A-vector - close the file
			CLOSE(UNIT=lyn)
				
			!// Print confirmation that subroutine ran successfully
			PRINT*, 'scanfromfile completed'
			
		!---------------------------------------------------------------------------------	
			!// Allocate space for x-vectors and check that it went well
			!// If not -> stop program
			ALLOCATE(mat1%xoldJac(mat1%n), STAT=res);mat1%xoldJac = 0.0
			IF (res/=0) THEN
				!//Print error message and stop program file
				PRINT*, 'Error in allocating xoldJac, status:'
				STOP
			END IF
			ALLOCATE(mat1%xnewJac(mat1%n), STAT=res);mat1%xnewJac = 0.0
			IF (res/=0) THEN
				!//Print error message and stop program file
				PRINT*, 'Error in allocating xnewJac, status:'
				STOP
			END IF
			ALLOCATE(mat1%xoldSOR(mat1%n), STAT=res);mat1%xoldSOR = 0.0
			IF (res/=0) THEN
				!//Print error message and stop program file
				PRINT*, 'Error in allocating xoldSOR, status:'
				STOP
			END IF
			ALLOCATE(mat1%xnewSOR(mat1%n), STAT=res);mat1%xnewSOR = 0.0
			IF (res/=0) THEN
				!//Print error message and stop program file
				PRINT*, 'Error in allocating xnewSOR, status:'
				STOP
			END IF
			
		END SUBROUTINE scanfromfile
! --------------------------- NEW SUBROUTINE----------------------------------------------
		SUBROUTINE Jacobi1iteration(mat1)
			!//////////////////////////////////////////////////////////
			!//
			!// Jacobi1iteration performs a single iteration with
			!// the Jacobi method.
			!// To be used for more iterations by Jacobi_it subroutine.
			!// Producing the resulting x-vector, xoldJac.
			!//
			!//////////////////////////////////////////////////////////
			IMPLICIT NONE
			CLASS(matrix)			::	mat1
			REAL					::	sum, diagA = 0
			INTEGER					::	i, j, k, l
			!-----------------------------------------------------------------------------
			DO i = 1,mat1%n
			!// Initialize sum to zero
			sum = 0
				DO k = mat1%irow(i),mat1%irow(i+1)-1
					j = mat1%jcol(k)
					!// Avoiding the diagonal element in the sum
					IF (i/=j) THEN
						sum = sum + mat1%A(k)*mat1%xoldJac(j)
					END IF
					!// Making a(i,i)
					IF (i==j) THEN
						diagA = mat1%A(k)
					END IF
				END DO
				mat1%xnewJac(i) = (1/diagA)*(mat1%b(i) - sum)
			END DO
		END SUBROUTINE Jacobi1iteration
! --------------------------- NEW SUBROUTINE ---------------------------------------------
		SUBROUTINE SOR1iteration(mat1, omega)
			!/////////////////////////////////////////////////////////////////
			!//
			!// SOR1iteration performs a single iteration with the successive
			!// over relaxation method, containing the Gauss-Seidel method.
			!// To be used for more iterations by SOR_it and produces the
			!// resulting x-vector, xoldSOR.
			!//
			!/////////////////////////////////////////////////////////////////
			IMPLICIT NONE
			CLASS(matrix)			::	mat1
			REAL, INTENT(IN)		::	omega
			REAL					::	diagA = 0, sum1, sum2
			INTEGER					::	i, j, k, l
			!-----------------------------------------------------------------------------
			DO i = 1,mat1%n
			!// Initialize sums to zero
			sum1 = 0
			sum2 = 0
				DO k = mat1%irow(i),mat1%irow(i+1)-1
					j = mat1%jcol(k)
					!// Avoiding the diagonal element in the sums
					IF (i>j) THEN
						sum1 = sum1 + mat1%A(k)*mat1%xnewSOR(j)
					END IF
					IF (i<j) THEN
						sum2 = sum2 + mat1%A(k)*mat1%xoldSOR(j)
					END IF
					!// Making a(i,i)
					IF (i==j) THEN
						diagA = mat1%A(k)
					END IF
				END DO
				mat1%xnewSOR(i) = (1/diagA)*(mat1%b(i) - sum1 - sum2)
			END DO
			!// Make the transformation from Gauss-Seidel iteration to SOR iteration
			DO l = 1,mat1%n
				mat1%xnewSOR(l) = omega*mat1%xnewSOR(l) + (1-omega)*mat1%xoldSOR(l)
			END DO
		END SUBROUTINE SOR1iteration
! --------------------------- NEW SUBROUTINE ---------------------------------------------
		SUBROUTINE Jacobi_it(mat1, tol, maxit)
			!///////////////////////////////////////////////////////////////////
			!//
			!// Jacobi_it calls the Jacobi1iteration subroutine
			!// as many times as it needs to do iterations.
			!// It calculates the difference between new iteration
			!// and old one to see if value is smaller than tol,
			!// in which case the loop is exited and complete.
			!//
			!///////////////////////////////////////////////////////////////////
			
			CLASS(matrix)					::	mat1
			REAL, INTENT(IN)				::	tol
			INTEGER, INTENT(IN)				::	maxit
			INTEGER							::	i, j, k, itsJac
			INTEGER							::	count = 1
			REAL							:: 	diff
			!-----------------------------------------------------------------------------
			DO i = 1, maxit
				mat1%itsJac = i
				diff = 0
				CALL mat1%Jacobi1_it()
				DO j = 1,mat1%n
					diff = SQRT(diff + (mat1%xnewJac(j) - mat1%xoldJac(j))**2)
				END DO
				!// Make the new x-vector the old one in case of a new iteration
				DO k = 1,mat1%n
					mat1%xoldJac(k) = mat1%xnewJac(k)
				END DO
				!// If the difference between new and old iteration is less
				!// than tol, exit the loop.
				IF (diff<tol) THEN
					EXIT
				END IF
				IF (i==count) THEN
					PRINT*,'i=', i, 'diffJac=', diff
					count = count + 20
				END IF
				IF (i>640) THEN
					PRINT*,'i=', i, 'diffJac=', diff
				END IF
			END DO
			
			PRINT*, 'diffJac=', diff
			PRINT*, 'xJac=', mat1%xoldJac(1:10)
			PRINT*, 'Jacobi iterations completed:', mat1%itsJac
		END SUBROUTINE Jacobi_it
! --------------------------- NEW SUBROUTINE ---------------------------------------------
		SUBROUTINE SOR_it(mat1, tol, omega, maxiter)
			!////////////////////////////////////////////////////////////////
			!//
			!// SOR_it calls the SOR1iteration subroutine as 
			!// many times as it needs to do iterations.
			!// It calculates the difference between new iteration
			!// and old one to see if the value is smaller than tol,
			!// in which case the loop is exited and complete.
			!//
			!////////////////////////////////////////////////////////////////
			
			CLASS(matrix)					::	mat1
			REAL, INTENT(IN)				::	tol, omega
			INTEGER, INTENT(IN)				::	maxiter
			INTEGER							::	i, j, k, itsSOR
			INTEGER							::	count = 1
			REAL							:: 	diff
			!-----------------------------------------------------------------------------
			DO i = 1, maxiter
				mat1%itsSOR= i
				diff = 0
				CALL mat1%SOR1_it(omega)
				DO j = 1,mat1%n
					diff = diff + (mat1%xnewSOR(j) - mat1%xoldSOR(j))**2
				END DO
				diff = SQRT(diff)
				!// Make the new x-vector the old one in case of a new iteration
				DO k = 1,mat1%n
					mat1%xoldSOR(k) = mat1%xnewSOR(k)
				END DO
				!// If the difference between new and old iteration is less
				!// than tol, exit the loop.
				IF (diff<tol) THEN
					EXIT
				END IF
			END DO
			
			PRINT*, 'diffSOR=', diff
			PRINT*, 'xSOR=', mat1%xoldSOR(1:10)
			PRINT*, 'SOR iterations completed:', mat1%itsSOR			
		END SUBROUTINE SOR_it
! --------------------------- NEW SUBROUTINE ---------------------------------------------
		SUBROUTINE dumptofile(mat1, xvector, tol, omega, maxit, maxiter)
			!//////////////////////////////////////////////////////////////////////
			!//
			!// dumptofile writes the xvectors to a file with appropriate headers.
			!//
			!//////////////////////////////////////////////////////////////////////
			CLASS(matrix)					::	mat1
			CHARACTER(LEN=80), INTENT(IN)	::	xvector	!filename
			CHARACTER(LEN=80), DIMENSION(8)	::	text	!array with text for file
			REAL, INTENT(IN)				::	tol, omega
			REAL(KIND=8), DIMENSION(8)		::	invars	!array with values for file
			INTEGER, INTENT(IN)				::	maxit, maxiter
			INTEGER							::	i, j
			INTEGER							::	olun = 12, res
			!-----------------------------------------------------------------------------
			OPEN(UNIT=olun, FILE=xvector, FORM='FORMATTED', IOSTAT=res)
			!// Check if the file was opened correctly
			IF (res/=0) THEN
				!// Print an error message and stop program
				PRINT*, 'Error in opening outfile for x-vector, status:', res
				STOP
			END IF
			
			!// Make arrays for shorter writing code
			!// Array containing strings for value identification
			text(1) = 'tolerance convergence limit:'
			text(2) = 'value used for the constant omega:'
			text(3) = 'maximum number of Jacobi iterations:'
			text(4) = 'maximum number of SOR iterations:'
			text(5) = 'number of Jacobi iterations, itsJac:'
			text(6) = 'number of SOR iterations, itsSOR:'
			text(7) = 'number of rows, length of b, n:'
			text(8) = 'number of non-zero values, length of A, nnz:'
			
			!// Array containing tol, omega, maxit, maxiter, itsJac, itsSOR, n and nnz
			invars(1) = tol
			invars(2) = omega
			invars(3) = maxit
			invars(4) = maxiter
			invars(5) = mat1%itsJac
			invars(6) = mat1%itsSOR
			invars(7) = mat1%n
			invars(8) = mat1%nnz
			
			!// Write tol, omega and number of maximum iterations for Jacobi and SOR
			DO i = 1, size(invars)
				WRITE(UNIT=olun, FMT='(A45, 4X, F10.5)', IOSTAT=res) &
					TRIM(text(i)), invars(i)
				!// Check if writing went well
				IF (res/=0) THEN
					!// Print error message
					PRINT*, 'Error in writing file, status:', res
					!// Close the file and stop the program
					CLOSE(UNIT=olun)
					STOP
				END IF
			END DO
			
			!// Write headers for x from Jacobi and SOR iterations respectively
			WRITE(UNIT=olun, FMT='(A)', IOSTAT=res) &
				'x from Jacobi      x from SOR'
			!// Check if writing went well
			IF (res/=0) THEN
				!// Print error message
				PRINT*, 'Error in writing x, status:', res
				!// Close the file and stop the program
				CLOSE(UNIT=olun)
				STOP
			END IF
			
			DO j = 1, mat1%n
				WRITE(UNIT=olun, FMT='(A2,I3,A2,F9.4,3X,A2,I3,A2,F9.4)', IOSTAT=res) &
					'x(' , j, ')=', mat1%xoldJac(j), 'x(', j, ')=', mat1%xoldSOR(j)
				!// Check if writing went well
				IF (res/=0) THEN
					!// Print error message
					PRINT*, 'Error in writing x to file, status:', res
					!// Close the file and stop the program
					CLOSE(UNIT=olun)
					STOP
				END IF
			END DO
			
			!// Done writing x to file, close the file
			CLOSE(UNIT=olun)
			PRINT*, 'Successfully dumped x'
			PRINT*, 'filename=', TRIM(xvector)
		END SUBROUTINE dumptofile
END MODULE MatSparse

