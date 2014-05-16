!PROGRAM scanfromfile
SUBROUTINE scanfromfile
	IMPLICIT NONE
	CHARACTER(LEN=80)				::	line
!	REAL(KIND=8),ALLOCATABLE		::	A(:)
!	REAL(KIND=8),ALLOCATABLE		::	b(:)
!	REAL(KIND=8),ALLOCATABLE		::	irow(:)
!	REAL(KIND=8),ALLOCATABLE		::	jcol(:)
	INTEGER, PARAMETER				::	lun = 10, lyn = 11
	INTEGER							::	i, j, k, l, m, n, nnz, res, rowval, rowvalold
	
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
	READ(line(j+1:),*), n
	!// Print n
!	PRINT*, 'n=', n
	
	ALLOCATE(b(n))
	! SHOULD check if it worked!!!
	
	!// Read the data in to the b-vector
	DO i = 1,n
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
		READ(line(j+1:),*), b(i)
	END DO
	!// Done reading file containing b-vector - close the file
	CLOSE(UNIT=lun)
	!// Print parts of b to check that reading has been done correctly
!	PRINT*, 'b=', b(1:10)
	
! -----------------------------------------------------------------------------

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
	READ(line(j+1:),*), n
	
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
	READ(line(j+1:),*), nnz
	!// Print n and nnz
!	PRINT*, 'n=', ',', n, 'nnz=', nnz
	
	ALLOCATE(A(nnz))
	! SHOULD check if it worked!!!
	ALLOCATE(irow(n+1))
	ALLOCATE(jcol(nnz))

	
	rowvalold = 0
	
	!// Read the data in to the A-vector
	DO i = 1,nnz
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
		READ(line(j+1:),*), A(i)
		
		!// Now make jcol and tempirow with help from (row,column)
		k=INDEX(line,',')
		l=INDEX(line,'(')
		m=INDEX(line,')')
		
		!// Read row from between '(' and ',' and then column from between ',' and ')'
		READ(line(l+1:k-1),*), rowval
		READ(line(k+1:m-1),*), jcol(i)
		
		!// Place number of value on new row in irow
		IF (rowval>rowvalold) THEN
			irow(rowval) = i
		END IF
		rowvalold = rowval	
	END DO
	
	irow(n+1) = nnz+1
	
	!// Done reading file containing A-vector - close the file
	CLOSE(UNIT=lyn)
	!// Print parts of A to check that reading has been done correctly
!	PRINT*, 'A=', A(1:10)
	PRINT*, 'irow=', irow(:)
!	PRINT*, 'jcol=', jcol(1:10)	
!// Print confirmation that program ran successfully
PRINT*, 'scanfromfile ferdig'
END SUBROUTINE scanfromfile
!END PROGRAM scanfromfile