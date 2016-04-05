! count the number of root within the interval
SUBROUTINE EXAUSITVE_SEARCH( F, lBracket, rBracket, nroot )
	USE PARAMS          ,ONLY : NUM_GRID
	IMPLICIT NONE
	REAL                     :: F
	REAL, INTENT(IN)         :: lBracket, rBracket
	INTEGER, INTENT(OUT)     :: nroot

	INTEGER                  :: i
	REAL                     :: dstep
	REAL                     :: grid(NUM_GRID), fgrid(NUM_GRID)

	nroot = 0
	dstep = (rBracket - lBracket)/(NUM_GRID - 1)
	grid = (/ (lBracket + dstep * i, i = 0, NUM_GRID - 1) /)

	DO i = 1, NUM_GRID
		fgrid(i) = F(grid(i))
	END DO

	DO i = 1, NUM_GRID - 1
		IF ( fgrid(i) * fgrid(i+1) < 0 ) nroot = nroot + 1
	END DO
END SUBROUTINE EXAUSITVE_SEARCH
