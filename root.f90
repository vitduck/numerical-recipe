MODULE Root
    IMPLICIT NONE 

    INTEGER, PARAMETER       :: IMAX      = 30
    INTEGER, PARAMETER       :: NGRID     = 100
    REAL, PARAMETER          :: TOLERANCE = 1.0E-5
    REAL, PARAMETER          :: PI        = 3.141592653589793

CONTAINS

SUBROUTINE Bracket( f, a, b, nroot, debug )
    IMPLICIT NONE

    REAL                     :: f
    REAL, INTENT(IN)         :: a, b
    INTEGER, INTENT(OUT)     :: nroot  
    INTEGER, INTENT(IN)      :: debug

    REAL                     :: h
    REAL                     :: grid(NGRID), fgrid(NGRID)
    INTEGER                  :: i

    OPTIONAL                 :: debug

    ! initialization 
    nroot = 0 
    h     = (b - a)/(NGRID - 1)
    grid  = (/ (a + h * i, i = 0, NGRID - 1) /)

    ! debug header 
    IF ( PRESENT(debug) .AND. debug == 1 ) &
        PRINT 1000, 'n', 'x', 'f(x)'

    ! value of f at each grid point 
    DO i = 1, NGRID
        fgrid(i) = f(grid(i))
        ! debug message 
        IF ( PRESENT(debug) .AND. debug == 1 ) &
            PRINT 2000, i, grid(i), fgrid(i) 
        END DO

    ! count the number of times f changes sign, 
    ! i.e. number of root  
    DO i = 1, NGRID - 1
        IF ( fgrid(i) * fgrid(i+1) < 0 ) nroot = nroot + 1
    END DO
    
    1000 FORMAT (A3, 2(A10, 3X))
    2000 FORMAT (I3, 2(F13.6))
END SUBROUTINE Bracket

SUBROUTINE Bisection( f, a, b, x_n, debug )
    IMPLICIT NONE

    REAL                     :: f
    REAL, INTENT(IN)         :: a, b
    REAL, INTENT(OUT)        :: x_n
    INTEGER, INTENT(IN)      :: debug

   REAL                      :: a_n, b_n, x_n_1, error
   INTEGER                   :: nstep = 0

    OPTIONAL debug 

    ! assumption: [a, b]
    IF ( a > b ) THEN
        PRINT '(A)', "WRONG BRACKET ORDER"
        RETURN 
    ! intermeidate value theorem 
    ELSE IF ( f(a)*f(b) > 0 ) THEN
        PRINT '(A)', "ROOT IS NOT BRACKETED"
        RETURN
    END IF

    ! initialization 
    a_n   = a
    b_n   = b 
    x_n_1 = a_n

    ! debug header 
    IF ( PRESENT(debug) .AND. debug == 1 ) & 
        PRINT 1000, 'n', 'a', 'b', 'x', 'f(x)', 'error'

    DO
        nstep = nstep + 1

        ! next guess 
        x_n = 0.5 * ( a_n + b_n )

        ! relative error 
        error = ABS((x_n - x_n_1)/x_n)

        ! debug 
        IF ( PRESENT(debug) .AND. debug == 1 ) & 
            PRINT 2000, nstep, a_n, b_n, x_n, f(x_n), error 

        ! termination
        x_n_1 = x_n
        IF ( error < TOLERANCE ) RETURN 

        ! reset the bracket 
        IF ( f(x_n)*f(a_n) < 0.0 ) THEN
            b_n = x_n
        ELSE
            a_n = x_n
        END IF
    END DO

    1000 FORMAT (A3, 3(A10, 3X), A13, 4X, A13 )
    2000 FORMAT (I3, 3F13.6, 2(4X, ES13.6))
END SUBROUTINE Bisection

SUBROUTINE False_Position( f, a, b, x_n, debug )
    IMPLICIT NONE
    
    REAL                     :: f
    REAL, INTENT(IN)         :: a, b
    REAL, INTENT(OUT)        :: x_n
    INTEGER, INTENT(IN)      :: debug

    REAL                     :: a_n, b_n, x_n_1, error
    INTEGER                  :: nstep = 0

    OPTIONAL debug   

    ! assumption: [a, b]
    IF ( a > b ) THEN
        PRINT '(A)', "WRONG BRACKET ORDER"
        RETURN 
    ! intermeidate value theorem 
    ELSE IF ( f(a)*f(b) > 0 ) THEN
        PRINT '(A)', "ROOT IS NOT BRACKETED"
        RETURN
    END IF
    
    ! initializationn
    a_n   = a 
    b_n   = b
    x_n_1 = a_n

    ! debug header 
    IF ( PRESENT(debug) .AND. debug == 1 ) & 
        PRINT 1000, 'n', 'a', 'b', 'x', 'f(x)', 'error'

    DO
        nstep = nstep + 1

        ! next guess 
        x_n = (a_n*f(b_n) - b_n*f(a_n))/(f(b_n) - f(a_n))
        
        ! relative error 
        error = ABS((x_n - x_n_1)/x_n)

        ! debug 
        IF ( PRESENT(debug) .AND. debug == 1 ) & 
            PRINT 2000, nstep, a_n, b_n, x_n, f(x_n), error 

        ! termination
        x_n_1 = x_n
        IF ( error < TOLERANCE ) RETURN 
        
        ! reset the bracket 
        IF ( f(x_n)*f(a_n) < 0.0 ) THEN
            b_n = x_n
        ELSE
            a_n = x_n
        END IF
    END DO

    1000 FORMAT (A3, 3(A10, 3X), A13, 4X, A13 )
    2000 FORMAT (I3, 3F13.6, 2(4X, ES13.6))
END SUBROUTINE False_Position

SUBROUTINE Newton_Raphson( f, df, x_0, x_n, debug )
    IMPLICIT NONE

    REAL                     :: f, df
    REAL, INTENT(IN)         :: x_0
    REAL, INTENT(OUT)        :: x_n
    INTEGER, INTENT(IN)      :: debug 
    
    REAL                     :: x_n_1, error
    INTEGER                  :: nstep = 0

    OPTIONAL debug 

    ! initializationn
    x_n_1 = x_0

    ! debug header 
    IF ( PRESENT(debug) .AND. debug == 1 ) & 
        PRINT 1000, 'n', 'x', 'f(x)', 'error'

    ! iteration 
    DO
        nstep = nstep + 1

        ! next guess 
        x_n = x_n_1 - f(x_n_1)/df(x_n_1) 

        ! relative error 
        error = ABS((x_n - x_n_1)/x_n_1)

        ! debug 
        IF ( PRESENT(debug) ) & 
            PRINT 2000, nstep, x_n, f(x_n), error

        ! termination
        x_n_1 = x_n 
        IF ( error < TOLERANCE ) RETURN
    END DO

    1000 FORMAT (A3, (A10, 3X), A13, 4X, A13 )
    2000 FORMAT (I3, F13.6, 2(4X, ES13.6))
END SUBROUTINE Newton_Raphson

SUBROUTINE Secant( f, x_0, x_1, x_n, debug )
    IMPLICIT NONE

    REAL                     :: f
    REAL, INTENT(IN)         :: x_0, x_1
    REAL, INTENT(OUT)        :: x_n
    INTEGER, INTENT(IN)      :: debug 

    REAL                     :: x_n_2, x_n_1, error  
    INTEGER                  :: nstep = 0 

    OPTIONAL debug 

    !initialization
    x_n_2  = x_0
    x_n_1  = x_1

    ! debug header 
    IF ( PRESENT(debug) .AND. debug == 1 ) & 
        PRINT 1000, 'n', 'x', 'f(x)', 'error'

    DO
        nstep = nstep + 1
        
        ! next guess 
        x_n = x_n - f(x_n_1) * (x_n_1 - x_n_2)/(f(x_n_1) - f(x_n_2))

        ! relative error 
        error  = ABS((x_n - x_n_1)/x_n_1)
    
        ! debug 
        IF ( PRESENT(debug) .AND. debug == 1 ) & 
            PRINT 2000, nstep, x_n, f(x_n), error 

        ! termination
        x_n_2 = x_n_1
        x_n_1 = x_n 
        IF ( error < TOLERANCE ) RETURN 
    END DO
    
    1000 FORMAT (A3, (A10, 3X), A13, 4X, A13 )
    2000 FORMAT (I3, F13.6, 2(4X, ES13.6))
END SUBROUTINE Secant 

SUBROUTINE Mueller( f, x_0, x_1, x_2, x_n, debug )
    IMPLICIT NONE

    REAL                     :: f
    REAL, INTENT(IN)         :: x_0, x_1, x_2
    REAL, INTENT(OUT)        :: x_n
    INTEGER, INTENT(IN)      :: debug

    REAL                     :: x_n_3, x_n_2, x_n_1, error
    REAL                     :: a, b, c, delta
    INTEGER                  :: nstep = 0

    OPTIONAL debug 

    ! initialiazation
    x_n_3 = x_0
    x_n_2 = x_1
    x_n_1 = x_2

    ! debug header 
    IF ( PRESENT(debug) .AND. debug == 1 ) & 
        PRINT 1000, 'n', 'x', 'f(x)', 'error'

    DO
        nstep = nstep + 1

        a = ((f(x_n_3) - f(x_n_1))*(x_n_2 - x_n_1) - (f(x_n_2) - f(x_n_1))*(x_n_3 - x_n_1)) / &
            ((x_n_3 - x_n_1)*(x_n_2 - x_n_1)*(x_n_3 - x_n_2))
        b = ((f(x_n_2) - f(x_n_1))*(x_n_3 - x_n_1)**2 - (f(x_n_3) - f(x_n_1))*(x_n_2 - x_n_1)**2) / &
            ((x_n_3 - x_n_1)*(x_n_2 - x_n_1)*(x_n_3 - x_n_2))
        c = f(x_n_1)
        
        ! special Delta 
        delta = b**2 - 4*a*c

        ! next guess 
        IF ( b >= 0 ) error = -2*c/(b + SQRT(delta))
        IF ( b  < 0 ) error = -2*c/(b - SQRT(delta))
        x_n = x_n_1 + error 
       
        ! debug 
        IF ( PRESENT(debug) .AND. debug == 1 ) & 
            PRINT 2000, nstep, x_n, f(x_n), error 

        ! termination
        x_n_3 = x_n_2
        x_n_2 = x_n_1
        x_n_1 = x_n
        error = ABS(error) 
        IF ( error < TOLERANCE ) RETURN 
    END DO
    
    1000 FORMAT (A3, (A10, 3X), A13, 4X, A13 )
    2000 FORMAT (I3, F13.6, 2(4X, ES13.6))
END SUBROUTINE Mueller
END MODULE Root
