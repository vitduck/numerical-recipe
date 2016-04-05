PROGRAM ROOT
	USE solve
	IMPLICIT NONE 
	INTEGER            :: numGrid, numRoot
	REAL               :: lBracket, rBracket
	REAL               :: guess, first_guess, second_guess, third_guess
	REAL               :: root
	
	lBracket = -9.5; rBracket = -3.5
	guess = -6.5
	first_guess = -8.5;	second_guess = -7.0; third_guess = -5.5
	
	WRITE (6,'(A)') "Using exhautive search to count root:"
	CALL exSearch( f, lBracket, rBracket, numGrid, numRoot )
	WRITE (6,'("There are", I3, " Root(s)",/)') numRoot
	
	WRITE (6,'(A)') "Bisection method:"
	CALL bisection( f, lBracket, rBracket, root )

	WRITE (6,'(A)') "Newton-Raphson method:"
	CALL newton_raphson( f, df, guess, root )

	WRITE (6,'(A)') "Secant method:"
	CALL secant( f, first_guess, second_guess, root )
	
	WRITE (6,'(A)') "False Position method:"
	CALL falsePosition( f, lBracket, rBracket, root )
	
	WRITE (6,'(A)') "Muller method:"
	CALL mueller( f, first_guess, second_guess, third_guess, root )
	
	WRITE (6,'(A)') "Bisection + Newton-Raphson method:"
	CALL hybrid_newton_bisect( f, df, lBracket, rBracket, root )
CONTAINS
REAL FUNCTION F(x) 
	IMPLICIT NONE 
	REAL :: x
	F = COS(x) - x*SIN(x)
END FUNCTION f
REAL FUNCTION DF(x) 
	IMPLICIT NONE 
	REAL :: x
	DF = - x*COS(x) - 2*SIN(x)
END FUNCTION DF
END PROGRAM ROOT
