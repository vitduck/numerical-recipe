PROGRAM BESSEL
    USE INTERPOLATION 

    IMPLICIT NONE

    INTEGER, PARAMETER       :: N = 10, NGRID = 90

    REAL                     :: t(N), J0(N), J1(N), J2(N)
    REAL                     :: x(NGRID+1), Lx(NGRID+1), S3(NGRID+1)  

    INTEGER                  :: i
    REAL                     :: dgrid 

    ! (t, y) for interpolation 
     t(:) = [ ( REAL(i), i = 1, N ) ]
    J0(:) = BESSEL_J0(t(:)) 
    J1(:) = BESSEL_J1(t(:)) 
    J2(:) = BESSEL_JN(2, t(:)) 

    ! 1d grid 
    dgrid = (t(N) - t(1))/NGRID
    x(:) = [ ( t(1) + dgrid*i , i = 0, NGRID ) ]

    ! Lagrange interpolation (vector call)
    ! N poly requires N+1 (x, f(x)) pair
    Lx(:) = LAGRANGE_INTERPOLATION(x(:), t(3:8), J1(3:8))

    ! cubic spline interpolation (vector call) 
    S3(:) = CUBIC_SPLINE_INTERPOLATION(x(:), t(:), J1(:))

    ! debug   
    DO i = 1, NGRID+1
        PRINT 1000, x(i), BESSEL_J1(x(i)), Lx(i), S3(i)
        1000 FORMAT (4F13.6)
    END DO 

CONTAINS 

SUBROUTINE PRINT_TABLE
    IMPLICIT NONE 

    INTEGER                  :: i 
    
    DO i = 1, N
        PRINT 1000, i, J0(i), J1(i), J2(i)
    END DO 

    1000 FORMAT (I3, (3F13.6))
END SUBROUTINE Print_Table

END PROGRAM BESSEL
