PROGRAM BESSEL
    USE INTERPOLATION 

    IMPLICIT NONE

    INTEGER, PARAMETER       :: N = 10
    INTEGER, PARAMETER       :: ORDER = 3

    REAL, ALLOCATABLE        :: J0(:), J1(:), J2(:)
    REAL, ALLOCATABLE        :: x(:), t(:), y(:) 
    REAL, ALLOCATABLE        :: Lx(:)  

    INTEGER                  :: i

    ALLOCATE(x(N), J0(N), J1(N), J2(N)) 
    ALLOCATE(Lx(N))

    ! N poly requires N+1 (x, f(x)) pair  
    ALLOCATE(t(ORDER+1))
    ALLOCATE(y(ORDER+1))

    ! initialization 
    x  = [ ( REAL(i), i = 1, N ) ]
    J0 = BESSEL_J0( x ) 
    J1 = BESSEL_J1( x ) 
    J2 = BESSEL_JN( 2,x ) 

    ! data for interpolation 
    t(:) =  x(4:6)
    y(:) = J1(4:6)

    ! vector call 
    Lx(:) = LAGRANGE(x(:), t(:), y(:), ORDER) 

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
    DEALLOCATE(t, y)

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
