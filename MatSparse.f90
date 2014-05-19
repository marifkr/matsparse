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
!			CHARACTER(LEN=80)		::	Avector, bvector, xvector !filenames
			REAL(KIND=8), POINTER	::	A(:) => null(), b(:) => null()
			REAL(KIND=8), POINTER	::	xoldJac(:) => null(), xnewJac(:) => null()
			REAL(KIND=8), POINTER	::	xoldSOR(:) => null(), xnewSOR(:) => null()
			REAL					::	tol, omega
			INTEGER, POINTER		::	irow(:) => null(), jcol(:) => null()
			INTEGER					::	n = 0, nnz = 0
			INTEGER					::	maxit, maxiter
			!// More variables here for x-vectors and other I might need
		CONTAINS
			PROCEDURE				:: scan => scanfromfile
			PROCEDURE				:: dump => dumptofile
			PROCEDURE				:: Jacobi1_it => Jacobi1iteration
			PROCEDURE				:: SOR1_it => SOR1iteration
!			PROCEDURE				:: Jacobi => Jacobi_it
!			PROCEDURE				:: SOR => SOR_it
		END TYPE matrix
		!// More global variables here
	CONTAINS
		SUBROUTINE scanfromfile(mat1,Avector,bvector)
			!/////////////////////////////////////////////////////////////
			!//
			!// scanfromfile reads files assigned in the main program
			!// to find n, nnz, A and b. The routine also provides irow
			!// and jcol to find positions of non-zero values in A.
			!//
			!/////////////////////////////////////////////////////////////
			IMPLICIT NONE
			CLASS(matrix)					::	mat1
			CHARACTER(LEN=80), INTENT(IN)	::	Avector, bvector !filenames
			CHARACTER(LEN=80)				::	line
			INTEGER, PARAMETER				::	lun = 10, lyn = 11
			INTEGER							::	i, j, k, l, m, n, nnz, res, rowval, rowvalold
			
			!// Confirm entry to module
			PRINT*, 'Got to module!'
		!---------------------------------------------------------------------------
			!// First b
			
			!// Open file containing the b-vector
			OPEN(UNIT=lun, FILE=bvector, FORM='FORMATTED', IOSTAT=res)
			!// Test if the file was opened correctly
			IF (res/= 0) THEN
				!// Print an error message
				PRINT *, 'Error in opening the b-vector file, status:', res
				!// Stop the program
				STOP
			END IF
			!// Read n from the first line in the file
			READ(UNIT=lun, FMT='(A)'), line
			IF(res/=0) THEN
				!// Print an error message
				PRINT *, 'Error in reading file, status:', res
				!// Close the file and stop the program
				CLOSE(UNIT=lun)
				STOP
			END IF
			!// Find position of :
			j=INDEX(line,':')
			!// Read everything after :
			READ(line(j+1:),*), mat1%n
			!// Print n
		!	PRINT*, 'n=', mat1%n

			ALLOCATE(mat1%b(mat1%n))
			! SHOULD check if it worked!!!

			!// Read the data in to the b-vector
			DO i = 1,mat1%n
				READ(UNIT=lun, FMT='(A)',IOSTAT=res), line
				IF(res/=0) THEN
					PRINT*, 'Error in reading file, status:', res
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
			!// Print parts of b to check that reading has been done correctly
		!	PRINT*, 'b=', mat1%b(1:10)

		! -----------------------------------------------------------------------------
			!// Now A

			!// Open file containing the sparse A-matrix
			OPEN(UNIT=lyn, FILE=Avector, FORM='FORMATTED', IOSTAT=res)
			!// Test if the file was opened correctly
			IF (res/= 0) THEN
				!// Print an error message
				PRINT *, 'Error in opening the file, status:', res
				!// Stop the program
				STOP
			END IF

			!// Read n and nnz from the first two lines in the file
			READ(UNIT=lyn, FMT='(A)'), line
			IF(res/=0) THEN
				!// Print an error message
				PRINT *, 'Error in reading A-file, status:', res
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
				PRINT *, 'Error in reading A-file, status:', res
				!// Close the file and stop the program
				CLOSE(UNIT=lyn)
				STOP
			END IF
			!// Find position of :
			j=INDEX(line,':')
			!// Read everything after :
			READ(line(j+1:),*), mat1%nnz
			!// Print number of rows, n and number of non-zero values in matrix A, nnz
		!	PRINT*, 'n=', ',', mat1%n, 'nnz=', mat1%nnz

			!// Allocating space for A, irow and jcol
			ALLOCATE(mat1%A(mat1%nnz))
			! SHOULD check if it worked!!!
			ALLOCATE(mat1%irow(mat1%n+1))
			ALLOCATE(mat1%jcol(mat1%nnz))

			rowvalold = 0

			!// Read the data in to the A-vector
			DO i = 1,mat1%nnz
				READ(UNIT=lyn, FMT='(A)',IOSTAT=res), line
				IF(res/=0) THEN
					PRINT*, 'Error in reading A-file, status:', res
					!// Close the file and stop the program
					CLOSE(UNIT=lun)
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
			!// Print parts of A, irow and jcol to check that reading has been done correctly
		!	PRINT*, 'A=', mat1%A(1:10)
		!	PRINT*, 'irow=', mat1%irow(1:10)
		!	PRINT*, 'jcol=', mat1%jcol(1:10)	
			!// Print confirmation that subroutine ran successfully
			PRINT*, 'scanfromfile completed'
				
			!// Allocate space for xold and xnew for both Jacobi and SOR iterations
			ALLOCATE(mat1%xoldJac(mat1%n));mat1%xoldJac = 0.0
			ALLOCATE(mat1%xnewJac(mat1%n));mat1%xnewJac = 0.0
			ALLOCATE(mat1%xoldSOR(mat1%n));mat1%xoldSOR = 0.0
			ALLOCATE(mat1%xnewSOR(mat1%n));mat1%xnewSOR = 0.0
			
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
			
			!// Make the new x-vector the old one in case of a new iteration
			DO l = 1,mat1%n
				mat1%xoldJac(l) = mat1%xnewJac(l)
			END DO
		
			PRINT*, 'xJac=', mat1%xoldJac(1:10)
			PRINT*, 'Jacobi1iteration completed'
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

			!// Make the new x-vector the old one in case of a new iteration
			!// and make it go from being Gauss-Seidel to SOR
			DO l = 1,mat1%n
				mat1%xoldSOR(l) = omega*mat1%xnewSOR(l) + (1-omega)*mat1%xoldSOR(l)
			END DO
			
			PRINT*, 'xSOR=', mat1%xoldSOR(1:10)
			PRINT*, 'SOR1iteration completed'
		END SUBROUTINE SOR1iteration
! --------------------------- NEW SUBROUTINE ---------------------------------------------
!		SUBROUTINE Jacobi_it(mat1, maxit, tol)
			!/////////////////
			!//
			

!		END SUBROUTINE Jacobi_it
! --------------------------- NEW SUBROUTINE ---------------------------------------------
!		SUBROUTINE SOR_it(mat1, maxiter)
			!////////////////////////////////////////////////////////////////
			!//
			!//
			
!		END SUBROUTINE SOR_it
! --------------------------- NEW SUBROUTINE ---------------------------------------------
		SUBROUTINE dumptofile(mat1, xvector, tol, omega, maxit, maxiter)
			!//////////////////////////////////////////////////////////////////////
			!//
			!// dumptofile writes the xvectors to a file with appropriate headers.
			!//
			!//////////////////////////////////////////////////////////////////////
			CLASS(matrix)					::	mat1
			CHARACTER(LEN=80), INTENT(IN)	::	xvector	!filename
			!CHARACTER(LEN=80)				::	toltext, omegatext, maxittext, maxitertext
			CHARACTER(LEN=80), DIMENSION(6)	::	text
			REAL, INTENT(IN)				::	tol, omega
			REAL(KIND=8), DIMENSION(6)		::	invars
			INTEGER, INTENT(IN)				::	maxit, maxiter
			INTEGER							::	i, j, k, l
			INTEGER							::	olun = 12, res
			
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
			text(5) = 'n:'
			text(6) = 'nnz:'
			
			!// Array containing tol, omega, maxit and maxiter
			PRINT*, tol, omega, maxit, maxiter
			invars(1) = tol
			invars(2) = omega
			invars(3) = maxit
			invars(4) = maxiter
			invars(5) = mat1%n
			invars(6) = mat1%nnz
			
			!// Write tol, omega and number of maximum iterations for Jacobi and SOR
			DO i = 1, size(invars)
				WRITE(UNIT=olun, FMT='(A37, 4X, F10.5)', IOSTAT=res) &
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

