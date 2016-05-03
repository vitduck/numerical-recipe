PROGRAM BESSEL
    USE INTERPOLATION 

    IMPLICIT NONE

    INTEGER, PARAMETER       :: N = 10

    REAL                     :: J0(N), J1(N), J2(N)
    REAL                     :: x(N), Lx(N), S3(N)  

    INTEGER                  :: i

    ! initialization 
     x(:) = [( REAL(i), i = 1, N )]
    J0(:) = BESSEL_J0(x(:)) 
    J1(:) = BESSEL_J1(x(:)) 
    J2(:) = BESSEL_JN(2, x(:)) 

    ! Lagrange interpolation
    ! N poly requires N+1 (x, f(x)) pair  
    ! vector call 
    Lx(:) = LAGRANGE_INTERPOLATION(x(:), x(4:6), J1(4:6))

    ! cubic spline interpolation  
    !DO i = 1, N 
        !S3(i) = CUBIC_SPLINE(temp(i), x(:), J1(:))
    !END DO 

    !! debug   
    DO i = 1, N 
        PRINT 1000, i, J1(i), Lx(i)
        1000 FORMAT (I3, (2F13.6))
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
