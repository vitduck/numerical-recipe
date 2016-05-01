MODULE INTERPOLATION 

CONTAINS 

REAL FUNCTION LAGRANGE( x, x_fit, f_fit, order )
    IMPLICIT NONE 

    REAL, INTENT(IN)         :: x
    REAL, INTENT(IN)         :: x_fit(:), f_fit(:)
    INTEGER, INTENT(IN)      :: order

    REAL                     :: coefficient 
    INTEGER                  :: i, j 

    ! initialization 
    LAGRANGE = 0 

    ! sanity check 
    IF ( UBOUND( x_fit, 1 ) /= UBOUND( f_fit, 1 ) ) &
        STOP 'Incompatible number of data between x and f(x)'

    IF ( (UBOUND( x_fit, 1 ) - LBOUND( x_fit, 1 )) < order ) THEN 
        PRINT 1000, order
        1000 FORMAT ('Not enough data points to fit to ', I2, 'th order polynomial')
        STOP
    END IF 

    ! loop through n data point 
    DO i = 1, UBOUND( x_fit, 1 ) 
        coefficient = 1.0     

        ! kronecker delta
        DO j = 1, UBOUND( x_fit, 1 ) 
            IF ( j /= i ) & 
                coefficient = coefficient * (x - x_fit(j))/(x_fit(i) - x_fit(j))
        END DO 

        LAGRANGE = LAGRANGE + coefficient * f_fit(i)
    END DO   
END FUNCTION LAGRANGE

END MODULE INTERPOLATION 
