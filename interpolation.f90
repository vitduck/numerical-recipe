MODULE INTERPOLATION 
    INTERFACE LAGRANGE 
        MODULE PROCEDURE LAGRANGE_SCALAR
        MODULE PROCEDURE LAGRANGE_VECTOR
    END INTERFACE

CONTAINS 

FUNCTION LAGRANGE_SCALAR(x, t, y, order)
    IMPLICIT NONE 

    REAL, INTENT(IN)         :: x
    REAL, INTENT(IN)         :: t(:), y(:)
    INTEGER, INTENT(IN)      :: order
    REAL                     :: LAGRANGE_SCALAR

    REAL                     :: coefficient 
    INTEGER                  :: j, k

    ! initialization 
    LAGRANGE_SCALAR = 0 

    ! sanity check
    CALL LAGRANGE_CHECK(t, y, order) 
   
    ! loop through n data point 
    DO j = 1, UBOUND(t, 1) 
        coefficient = 1.0     

        ! kronecker delta
        DO k = 1, UBOUND(t, 1) 
            IF ( j /= k ) & 
                coefficient = coefficient * (x - t(k))/(t(j) - t(k))
        END DO 

        LAGRANGE_SCALAR = LAGRANGE_SCALAR + coefficient * y(j)
    END DO   
END FUNCTION LAGRANGE_SCALAR

FUNCTION LAGRANGE_VECTOR(x, t, y, order) 
    IMPLICIT NONE 

    REAL, INTENT(IN)         :: x(:)
    REAL, INTENT(IN)         :: t(:), y(:)
    INTEGER, INTENT(IN)      :: order
    REAL                     :: LAGRANGE_VECTOR(LBOUND(x, 1):UBOUND(x, 1))
    
    INTEGER                  :: i  

    ! initialization 
    LAGRANGE_VECTOR(:) = 0 

    ! sanity check
    CALL LAGRANGE_CHECK(t, y, order) 

    DO i = 1, UBOUND(x, 1) 
        LAGRANGE_VECTOR(i) = LAGRANGE_SCALAR(x(i), t, y, order)
    END DO
END FUNCTION LAGRANGE_VECTOR

SUBROUTINE LAGRANGE_CHECK(t, y, order) 
    IMPLICIT NONE 
    
    REAL, INTENT(IN)         :: t(:), y(:)
    INTEGER, INTENT(IN)      :: order

    IF ( UBOUND(t, 1) /= UBOUND(y, 1) ) &
        STOP 'Incompatible number of data between x and f(x)'

    IF ( (UBOUND(t, 1) - LBOUND(t, 1)) < order ) THEN 
        PRINT 1000, order
        1000 FORMAT ('Not enough data points to fit to ', I2, 'th order polynomial')
        STOP
    END IF 
END SUBROUTINE LAGRANGE_CHECK

END MODULE INTERPOLATION 
