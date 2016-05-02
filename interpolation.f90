MODULE INTERPOLATION 
    INTERFACE LAGRANGE 
        MODULE PROCEDURE LAGRANGE_SCALAR
        MODULE PROCEDURE LAGRANGE_VECTOR
    END INTERFACE

CONTAINS 

FUNCTION LAGRANGE_SCALAR(x, x_fit, f_fit, order)
    IMPLICIT NONE 

    REAL, INTENT(IN)         :: x
    REAL, INTENT(IN)         :: x_fit(:), f_fit(:)
    INTEGER, INTENT(IN)      :: order
    REAL                     :: LAGRANGE_SCALAR

    REAL                     :: coefficient 
    INTEGER                  :: j, k

    ! initialization 
    LAGRANGE_SCALAR = 0 

    ! sanity check
    CALL LAGRANGE_CHECK(x_fit, f_fit, order) 
   
    ! loop through n data point 
    DO j = 1, UBOUND(x_fit, 1) 
        coefficient = 1.0     

        ! kronecker delta
        DO k = 1, UBOUND(x_fit, 1) 
            IF ( j /= k ) & 
                coefficient = coefficient * (x - x_fit(k))/(x_fit(j) - x_fit(k))
        END DO 

        LAGRANGE_SCALAR = LAGRANGE_SCALAR + coefficient * f_fit(j)
    END DO   
END FUNCTION LAGRANGE_SCALAR

FUNCTION LAGRANGE_VECTOR(x, x_fit, f_fit, order) 
    IMPLICIT NONE 

    REAL, INTENT(IN)         :: x(:)
    REAL, INTENT(IN)         :: x_fit(:), f_fit(:)
    INTEGER, INTENT(IN)      :: order
    REAL                     :: LAGRANGE_VECTOR(LBOUND(x, 1):UBOUND(x, 1))
    
    INTEGER                  :: i  

    ! initialization 
    LAGRANGE_VECTOR = 0 

    ! sanity check
    CALL LAGRANGE_CHECK( x_fit, f_fit, order ) 

    DO i = 1, UBOUND(x, 1) 
        LAGRANGE_VECTOR(i) = LAGRANGE_SCALAR(x(i), x_fit, f_fit, order)
    END DO
END FUNCTION LAGRANGE_VECTOR

SUBROUTINE LAGRANGE_CHECK(x_fit, f_fit, order) 
    IMPLICIT NONE 
    
    REAL, INTENT(IN)         :: x_fit(:), f_fit(:)
    INTEGER, INTENT(IN)      :: order

    IF ( UBOUND(x_fit, 1) /= UBOUND(f_fit, 1) ) &
        STOP 'Incompatible number of data between x and f(x)'

    IF ( (UBOUND(x_fit, 1) - LBOUND(x_fit, 1)) < order ) THEN 
        PRINT 1000, order
        1000 FORMAT ('Not enough data points to fit to ', I2, 'th order polynomial')
        STOP
    END IF 
END SUBROUTINE LAGRANGE_CHECK
END MODULE INTERPOLATION 
