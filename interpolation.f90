MODULE INTERPOLATION 
    INTERFACE LAGRANGE_INTERPOLATION
        MODULE PROCEDURE LAGRANGE
        MODULE PROCEDURE LAGRANGE_ARRAY
    END INTERFACE

    INTERFACE CUBIC_SPLINE_INTERPOLATION 
        MODULE PROCEDURE CUBIC_SPLINE
        MODULE PROCEDURE CUBIC_SPLINE_ARRAY
    END INTERFACE CUBIC_SPLINE_INTERPOLATION 

CONTAINS 

    !------------------------! 
    ! Lagrange Interpolation ! 
    !------------------------! 
    FUNCTION LAGRANGE(x, t, y) 
        IMPLICIT NONE 

        REAL, INTENT(IN)               :: x
        REAL, DIMENSION(:), INTENT(IN) :: t, y
        REAL                           :: LAGRANGE

        REAL                           :: coefficient
        INTEGER                        :: j, k

        ! sanity check
        IF ( SIZE(t) /= SIZE(y) ) STOP 'Incompatible number of data between x and f(x)'
        
        ! initialization 
        LAGRANGE = 0 
        
        ! loop through n data point 
        DO j = 1, SIZE(t)
            coefficient = 1.0     

            ! kronecker delta
            DO k = 1, SIZE(t) 
                IF ( k /= j ) & 
                    coefficient = coefficient * (x - t(k)) / (t(j) - t(k))
            END DO 

            LAGRANGE = LAGRANGE + coefficient * y(j)
        END DO
    END FUNCTION LAGRANGE

    FUNCTION LAGRANGE_ARRAY(x, t, y)
        IMPLICIT NONE 

        REAL, DIMENSION(:), INTENT(IN) :: x, t, y
        REAL, DIMENSION(SIZE(x))       :: LAGRANGE_ARRAY

        INTEGER                        :: i

        ! initialization 
        LAGRANGE_ARRAY = 0 

        ! sanity check
        IF ( SIZE(t) /= SIZE(y) ) STOP 'Incompatible number of data between x and f(x)'

        DO i = 1, SIZE(x)
            LAGRANGE_ARRAY(i) = LAGRANGE(x(i), t, y)
        END DO
    END FUNCTION LAGRANGE_ARRAY

    !----------------------------! 
    ! Cubic Spline Interpolation ! 
    !----------------------------! 
    FUNCTION CUBIC_SPLINE(x, t, y) 
        IMPLICIT NONE 

        REAL, INTENT(IN)               :: x 
        REAL, DIMENSION(:), INTENT(IN) :: t, y
        REAL                           :: CUBIC_SPLINE

        REAL, DIMENSION(SIZE(t))       :: S2

        ! sanity check
        IF ( SIZE(t) /= SIZE(y) ) STOP 'Incompatible number of data between x and f(x)'

        ! solve for the set of second derivatives 
        CALL CUBIC_SPLINE_INIT(t, y, S2)

        CUBIC_SPLINE = CUBIC_SPLINE_POLYNOMIAL(x, t, y, S2)
    END FUNCTION CUBIC_SPLINE

    FUNCTION CUBIC_SPLINE_ARRAY(x, t, y) 
        IMPLICIT NONE 

        REAL, DIMENSION(:), INTENT(IN) :: x, t, y
        REAL, DIMENSION(SIZE(x))       :: CUBIC_SPLINE_ARRAY

        REAL, DIMENSION(SIZE(t))       :: S2 
        INTEGER                        :: i

        ! sanity check
        IF ( SIZE(t) /= SIZE(y) ) STOP 'Incompatible number of data between x and f(x)'

        ! solve for the set of second derivatives 
        CALL CUBIC_SPLINE_INIT(t, y, S2)

        DO i = 1, SIZE(x)  
            CUBIC_SPLINE_ARRAY(i) = CUBIC_SPLINE_POLYNOMIAL(x(i), t, y, S2)
        END DO 
    END FUNCTION CUBIC_SPLINE_ARRAY

    FUNCTION CUBIC_SPLINE_POLYNOMIAL(x, t, y, S2)
        IMPLICIT NONE 

        REAL, INTENT(IN)               :: x
        REAL, DIMENSION(:), INTENT(IN) :: t, y
        REAL, DIMENSION(:), INTENT(IN) :: S2
        REAL                           :: CUBIC_SPLINE_POLYNOMIAL

        REAL                           :: h, A, B, C, D
        INTEGER                        :: i

        ! bracket x value 
        i = CUBIC_SPLINE_INDEX(x, t)  

        h = t(i+1)  - t(i) 
        A = S2(i+1) / (6 * h)
        B = S2(i)   / (6 * h)
        C = y(i+1)  / h       - S2(i+1) * h / 6
        D = y(i)    / h       - S2(i)   * h / 6

        CUBIC_SPLINE_POLYNOMIAL = & 
            A * (x - t(i))**3 - B * (x - t(i+1))**3 + C * (x - t(i)) - D * (x - t(i+1))
    END FUNCTION CUBIC_SPLINE_POLYNOMIAL

    SUBROUTINE CUBIC_SPLINE_INIT(t, y, S2)
        IMPLICIT NONE 

        REAL, DIMENSION(:), INTENT(IN)  :: t, y
        REAL, DIMENSION(:), INTENT(OUT) :: S2 

        REAL, DIMENSION(SIZE(t))        :: beta
        REAL                            :: a_i, b_i, c_i, r_i
        INTEGER                         :: i 

        ! natural spline
        ! this gives rise to a N-2 linear system of equations  
        S2 = 0.0

        ! note: i = 2, N-1
        ! beta(1) and beta(N) is undefined because of boundary condition 
        beta(2) = 2 * (t(3) - t(1))
        S2(2)   = 6 * ((y(3) - y(2)) / (t(3) - t(2)) - (y(2) - y(1)) / (t(2) - t(1)))  

        ! forward elimination 
        DO i = 3, SIZE(t) - 1 
            ! off-diagonal term ( symmetric matrix A(i) = C(i-1) )
            a_i = t(i) - t(i-1) 
            c_i = a_i

            ! diagonal term   
            b_i = 2.0 * (t(i+1) - t(i-1))

            ! right-hand side 
            r_i = 6.0 * ((y(i+1) - y(i)) / (t(i+1) - t(i)) - (y(i) - y(i-1)) / (t(i) - t(i-1)))

            ! eliminiation 
            beta(i) = b_i - a_i * c_i     / beta(i-1) 
            S2(i)   = r_i - a_i * S2(i-1) / beta(i-1)
        END DO 

        ! backward substitution 
        S2(SIZE(t) - 1) = S2(SIZE(t) - 1) / beta(SIZE(t) - 1)
        DO i = SIZE(t) - 2, 2, -1
            c_i   = t(i+1) - t(i) 
            S2(i) = (S2(i) - c_i * S2(i+1)) / beta(i) 
        END DO 
    END SUBROUTINE CUBIC_SPLINE_INIT

    FUNCTION CUBIC_SPLINE_INDEX(x, t)
        IMPLICIT NONE 

        REAL, INTENT(IN)               :: x
        REAL, DIMENSION(:), INTENT(IN) :: t
        INTEGER                        :: CUBIC_SPLINE_INDEX

        INTEGER                        :: a, b, m 

        ! initialization 
        ! preseving the index of t(:)
        a = 1
        b = SIZE(t)

        ! bisection 
        DO WHILE ( (b - a) > 1 )
        !   ! mid point 
            m = (a + b) / 2 
            IF ( x > t(m) ) THEN 
                a = m  
            ELSE 
                b = m 
            END IF
        END DO 

        ! points at the boundaries 
        IF ( x == t(a) ) THEN 
            CUBIC_SPLINE_INDEX = 1 
        ELSE IF ( x == t(b) ) THEN 
            CUBIC_SPLINE_INDEX = b-1
        ! otherwise lower bound 
        ELSE 
            CUBIC_SPLINE_INDEX = a 
        END IF 
    END FUNCTION CUBIC_SPLINE_INDEX

END MODULE INTERPOLATION 
