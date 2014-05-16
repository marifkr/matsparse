!/////////////////////////////////////////////////////////////
!//
!// Module MatSparse containing subroutines for 
!// - reading number of rows, columns (n) and non-zero 
!//   values (nnz) from two files containing a matrix A,
!//   and a vector b, and reading A and b.
!//	- performing Jacobi and SOR iterations to find
!//   the x-vector in Ax=b.
!// - dumping the x-vector to a file.
!//
!////////////////////////////////////////////////////////////
MODULE MatSparse
	IMPLICIT NONE
		TYPE	::	matrix
			CHARACTER(LEN=80)		:: Avector, bvector!, xvector !filenames
			REAL(KIND=8), POINTER	:: A(:), b(:), x(:), xnew(:), xold(:)
			INTEGER, POINTER		:: irow(:)
			INTEGER, POINTER		:: jcol(:)
			INTEGER					:: n
			INTEGER					:: nnz
			!// More variables here for x-vectors and other I might need
		CONTAINS
			PROCEDURE				:: scan => scanfromfile
!			PROCEDURE				:: dump => dumptofile
!			PROCEDURE				:: Jacobi1_it => Jacobi1iteration
!			PROCEDURE				:: SOR1_it => SOR1iteration
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
			CLASS(matrix)						::	mat1
				CHARACTER(LEN=80), INTENT(IN)	::	Avector, bvector!, xvector !filenames
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
					PRINT *, 'Error in opening the file, status:', res
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
	
				ALLOCATE(mat1%b(n))
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
					PRINT *, 'Error in reading file, status:', res
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
					PRINT *, 'Error in reading file, status:', res
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
	
				ALLOCATE(mat1%A(mat1%nnz))
				! SHOULD check if it worked!!!
				ALLOCATE(mat1%irow(mat1%n+1))
				ALLOCATE(mat1%jcol(mat1%nnz))

	
				rowvalold = 0
	
				!// Read the data in to the A-vector
				DO i = 1,mat1%nnz
					READ(UNIT=lyn, FMT='(A)',IOSTAT=res), line
					IF(res/=0) THEN
						PRINT*, 'Error in reading file, status:', res
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
	
				mat1%irow(mat1%n+1) = mat1%nnz+1
	
				!// Done reading file containing A-vector - close the file
				CLOSE(UNIT=lyn)
				!// Print parts of A, irow and jcol to check that reading has been done correctly
			!	PRINT*, 'A=', mat1%A(1:10)
			!	PRINT*, 'irow=', mat1%irow(1:10)
			!	PRINT*, 'jcol=', mat1%jcol(1:10)	
				!// Print confirmation that program ran successfully
				PRINT*, 'scanfromfile completed'

		END SUBROUTINE scanfromfile
		!// More subroutine declarations here

! --------------------------- NEW SUBROUTINE----------------------------------------------
!		SUBROUTINE Jacobi1iteration
			!//////////////////////////////////////////////////////////
			!//
			!// Jacobi1iteration performs a single iteration with
			!// the Jacobi method.
			!//
			!//////////////////////////////////////////////////////////
!			IMPLICIT NONE
!			CLASS(matrix)	::	mat1
!				INTEGER						::	i, j, k
				
!				mat1%xold() = 0
!				DO i = mat1%n
					


!				PRINT*, 'Jacobi1iteration completed'
!		END SUBROUTINE Jacobi1iteration
! --------------------------- NEW SUBROUTINE ---------------------------------------------
END MODULE MatSparse