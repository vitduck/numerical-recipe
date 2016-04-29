PROGRAM BENCHMARK
    USE ROOT

    REAL                     :: a = 0, b = 1
    REAL                     :: x, x_0 = 0.25, x_1 = 0.5, x_2 = 0.75
    REAL, EXTERNAL           :: f, df

    INTEGER                  :: nroot 
    INTEGER                  :: debug = 1

    CALL BRACKET( f, a, b, nroot )
    PRINT 1000, nroot
    PRINT *

    PRINT 2000, a, b
    CALL BISECTION( f, a, b, x, debug )
    PRINT *

    PRINT 3000, a, b 
    CALL FALSE_POSITION( f, a, b, x, debug )
    PRINT *
    
    PRINT 4000, x_0
    CALL NEWTON_RAPHSON( f, df, x_0, x, debug )
    PRINT *

    PRINT 5000, x_0, x_1
    CALL SECANT( f, x_0, x_1, x, debug )
    PRINT *

    PRINT 6000, x_0, x_1, x_2
    CALL MUELLER( f, x_0, x_1, x_2, x, debug )

    ! format
    1000 FORMAT ('There is ', I2, ' root(s)')
    2000 FORMAT ('Bisection method: a = ', F7.3, ', b = ', F7.3)
    3000 FORMAT ('False position method: a = ', F7.3, ', b = ', F7.3)
    4000 FORMAT ('Newton-Raphson method: x_0 = ', F7.3)
    5000 FORMAT ('Secant method: x_0 = ', F7.3, ', x_1 = ', F7.3)
    6000 FORMAT ('Mueller method: x_0 = ', F7.3, ', x_1 = ', F7.3, ', x_2 = ', F7.3 )
END PROGRAM BENCHMARK 

REAL FUNCTION f(x) 
    IMPLICIT NONE 

    REAL, INTENT(IN)         :: x 

    f = COS(x) - x*SIN(x) 
END FUNCTION f

REAL FUNCTION df(x) 
    IMPLICIT NONE 

    REAL, INTENT(IN)         :: x 

    df = -2*SIN(x) - x*COS(x) 
END FUNCTION df
