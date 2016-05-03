MODULE INTERPOLATION 
    INTERFACE LAGRANGE_INTERPOLATION
        MODULE PROCEDURE LAGRANGE_SCALAR
        MODULE PROCEDURE LAGRANGE_VECTOR
    END INTERFACE

CONTAINS 

FUNCTION LAGRANGE_SCALAR(x, t, y) 
    IMPLICIT NONE 
    
    REAL, INTENT(IN)         :: x
    REAL, INTENT(IN)         :: t(:), y(:)
    REAL                     :: LAGRANGE_SCALAR

    ! sanity check
    IF ( SIZE(t(:)) /= SIZE(y(:)) ) &
        STOP 'Incompatible number of data between x and f(x)'

    LAGRANGE_SCALAR = LAGRANGE_POLYNOMIAL(x, t, y)
END FUNCTION LAGRANGE_SCALAR

FUNCTION LAGRANGE_VECTOR(x, t, y)
    IMPLICIT NONE 

    REAL, INTENT(IN)         :: x(:)
    REAL, INTENT(IN)         :: t(:), y(:)
    REAL                     :: LAGRANGE_VECTOR(SIZE(x(:)))
    
    INTEGER                  :: i  

    ! initialization 
    LAGRANGE_VECTOR(:) = 0 
    
    ! sanity check
    IF ( SIZE(t(:)) /= SIZE(y(:)) ) &
        STOP 'Incompatible number of data between x and f(x)'

    DO i = 1, SIZE(x(:))
        LAGRANGE_VECTOR(i) = LAGRANGE_POLYNOMIAL(x(i), t(:), y(:))
    END DO
END FUNCTION LAGRANGE_VECTOR

FUNCTION LAGRANGE_POLYNOMIAL(x, t, y)
    IMPLICIT NONE 

    REAL, INTENT(IN)         :: x
    REAL, INTENT(IN)         :: t(:), y(:)
    REAL                     :: LAGRANGE_POLYNOMIAL

    REAL                     :: coefficient 
    INTEGER                  :: j, k

    ! initialization 
    LAGRANGE_POLYNOMIAL = 0 
      
    ! loop through n data point 
    DO j = 1, SIZE(t(:))
        coefficient = 1.0     

        ! kronecker delta
        DO k = 1, SIZE(t(:)) 
            IF ( j /= k ) & 
                coefficient = coefficient * (x - t(k))/(t(j) - t(k))
        END DO 

        LAGRANGE_POLYNOMIAL = LAGRANGE_POLYNOMIAL + coefficient * y(j)
    END DO   
END FUNCTION LAGRANGE_POLYNOMIAL

!FUNCTION CUBIC_SPLINE(x, t, y)
    !IMPLICIT NONE 

    !REAL, INTENT(IN)         :: x
    !REAL, INTENT(IN)         :: t(:), y(:)
    !REAL                     :: CUBIC_SPLINE

    !REAL                     :: S2(LBOUND(t(:),1):UBOUND(t(:),1))
    !REAL                     :: h 
    !REAL                     :: A, B, C
    !INTEGER                  :: i

    !! initialization of second derivative 
    !CALL CUBIC_SPLINE_INIT(t(:), y(:), S2(:))

    !! bracket x value 
    !i = CUBIC_SPLINE_INDEX(x, t(:))  

    !! horner polynomial 
    !! S_i(x) = y_i + (x-t_i)[C + (x-t_i)[B + (x-t_i)A]]
    !h = t(i+1) - t(i) 
    !A = (S2(i+1) - S2(i))/(6*h)
    !B = S2(i)/2 
    !C = -h*S2(i+1)/6 - h*S2(i)/3 + (y(i+1) - y(i))/h

    !CUBIC_SPLINE = y(i) + (x-t(i))*(C + (x-t(i))*(B + (x-t(i)*A)))
!END FUNCTION CUBIC_SPLINE

!-----------
! AUXILARY !
!-----------
!SUBROUTINE CUBIC_SPLINE_INIT(t, y, s2)
    !IMPLICIT NONE 

    !REAL, INTENT(IN)         :: t(:), y(:) 
    !REAL, INTENT(OUT)        :: s2(:) 

    !REAL                     :: beta(LBOUND(t(:),1):UBOUND(t(:),1))
    !REAL                     :: a, b, c, r
    !INTEGER                  :: i  

    !! natural spline
    !! this gives rise to a N-2 linear system 
    !s2(LBOUND(t(:),1)) = 0.0 
    !s2(UBOUND(t(:),1)) = 0.0 

    !! initialization 
    !! note: i = 2, N-1
    !! beta(1) and beta(N) is undefined because of boundary condition 
    !beta(2) = 2*(t(3) - t(1))
    !s2(2)   = 6*((y(3) - y(2))/(t(3)-t(2)) - (y(2) - y(1))/(t(2) - t(1)))  

    !! forward elimination 
    !DO i = 3, UBOUND(t(:),1) - 1 
        !! off-diagonal term ( symmetric matrix A(i) = C(i-1) )
        !a = t(i) -t(i-1) 
        !c = a 

        !! diagonal term   
        !b = 2.0*(t(i+1) - t(i-1))
        
        !! right-hand side 
        !r = 6.0*((y(i+1)-y(i))/(t(i+1)-t(i)) - (y(i)-y(i-1))/(t(i)-t(i-1)))

        !! eliminiation 
        !beta(i) = b - a*c/beta(i-1) 
        !s2(i)   = r - a*s2(i-1)/beta(i-1)
    !END DO 

    !! backward substitution 
    !s2(UBOUND(t(:),1)-1) = s2(UBOUND(t(:),1)-1)/beta(UBOUND(t(:),1)-1)
    !DO i = UBOUND(t(:),1) - 2, 2, -1
        !! recover the off-diagonal term 
        !c     = s2(i+1) - s2(i) 
        !s2(i) = ( s2(i) - c*s2(i+1))/beta(i) 
    !END DO 
!END SUBROUTINE CUBIC_SPLINE_INIT

FUNCTION CUBIC_SPLINE_INDEX(x, t)
    IMPLICIT NONE 

    REAL, INTENT(IN)         :: x
    REAL, INTENT(IN)         :: t(:)
    INTEGER                  :: CUBIC_SPLINE_INDEX

    INTEGER                  :: a, b, m 

    ! initialization 
    ! preseving the index of t(:)
    a = 1
    b = SIZE(t(:))

    ! bisection 
    DO WHILE ( (b - a) > 1 )
    !   ! mid point 
        m = (a + b)/2 
        IF ( x > t(m) ) THEN 
            !set lower bound 
            a = m  
        ELSE 
            ! set upper bound 
            b = m 
        END IF
    END DO 

    ! points at the boundaries 
    IF ( x == t(a) ) THEN 
        CUBIC_SPLINE_INDEX = 1 
    ELSE IF ( x == t(b) ) THEN 
        CUBIC_SPLINE_INDEX = b-1
    ! somewhere in the midle 
    ELSE 
        CUBIC_SPLINE_INDEX = a 
    END IF 
END FUNCTION CUBIC_SPLINE_INDEX

END MODULE INTERPOLATION 
