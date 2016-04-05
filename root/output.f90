! print value of root for each iteration
SUBROUTINE PRINT_OUTPUT(root, error, numIte, info)
	USE PARAMS,          ,ONLY: MAX_ITERATION, TOLERANCE 
	IMPLICIT NONE
	REAL, INTENT(IN)         :: root, error
	INTEGER, INTENT(IN)      :: numIte
	INTEGER, INTENT(OUT)     :: info

	WRITE (6,1000) numIte, root, error
	IF ( numIte == MAX_ITERATION .AND. ABS(error) > TOLERANCE ) THEN
		info = 1
		WRITE (6,'(A)') "Maximum iteration reaches!"
		WRITE (6,2000) root, ABS(error)
	ELSE IF ( ABS(error) <= TOLERANCE ) THEN
		info = 2
		WRITE (6,3000) root, error, numIte
	END IF

	1000 FORMAT ("Iteration =", I3, " x =", ES15.8, " error =", ES15.8)
	2000 FORMAT ("x =", ES15.8, " error =", ES15.8,/)
	3000 FORMAT ("Root =", ES15.8, " error =", ES15.8, " after", I3, " iterations",/)
END SUBROUTINE PRINT_OUTPUT
