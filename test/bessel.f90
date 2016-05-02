PROGRAM BESSEL
    USE INTERPOLATION 

    IMPLICIT NONE

    INTEGER, PARAMETER       :: N = 10
    INTEGER, PARAMETER       :: ORDER = 3

    REAL, ALLOCATABLE        :: J0(:), J1(:), J2(:)
    REAL, ALLOCATABLE        :: x(:), x_fit(:), f_fit(:) 
    REAL, ALLOCATABLE        :: Lx(:)  

    INTEGER                  :: i

    ALLOCATE(x(N), J0(N), J1(N), J2(N)) 
    ALLOCATE(Lx(N))

    ! N poly requires N+1 (x, f(x)) pair  
    ALLOCATE(x_fit(ORDER+1))
    ALLOCATE(f_fit(ORDER+1))

    ! initialization 
    x  = [ ( REAL(i), i = 1, N ) ]
    J0 = BESSEL_J0( x ) 
    J1 = BESSEL_J1( x ) 
    J2 = BESSEL_JN( 2,x ) 

    ! data for interpolation 
    x_fit =  x(4:6)
    f_fit = J1(4:6)

    ! vector call 
    Lx(:) = LAGRANGE(x(:), x_fit(:), f_fit(:), ORDER) 

    ! scalar call 
    !DO i = 1, N
        !Lx(i) = LAGRANGE(x(i), x_fit(:), f_fit(:), ORDER) 
    !END DO 

    DO i = 1, N 
        PRINT 1000, i, J1(i), Lx(i) 
        1000 FORMAT (I3, (2F13.6))
    END DO 
    
    ! free memory 
    DEALLOCATE(x, J0, J1, J2)
    DEALLOCATE(Lx) 
    DEALLOCATE(x_fit, f_fit)

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
