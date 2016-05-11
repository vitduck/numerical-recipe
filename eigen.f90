MODULE EIGEN 

CONTAINS 

    SUBROUTINE TRIDIAGONAL(A, B, C, X, R) 
        IMPLICIT NONE 

        REAL, DIMENSION(:), INTENT(IN)  :: A, B, C, R
        REAL, DIMENSION(:), INTENT(OUT) :: X 

        REAL, DIMENSION(:), ALLOCATABLE :: beta
        INTEGER                         :: i, NSIZE

        ! initialization 
        NSIZE = SIZE(B)
        ALLOCATE(beta(NSIZE))

        beta(1) = B(1) 
        X(1)    = R(1) 

        ! Forward elimination of x_i  
        DO i = 2, NSIZE 
            beta(i) = B(i) - A(i) * C(i-1) / beta(i-1)
            X(i)    = R(i) - A(i) * X(i-1) / beta(i-1)
        END DO 

        ! Backward substitution 
        X(NSIZE) = X(NSIZE) / beta(NSIZE) 
        DO i = NSIZE - 1, 1, -1 
            X(i) = (X(i) - C(i) * X(i+1)) / beta(i)
        END DO 

        DEALLOCATE(beta)
    END SUBROUTINE TRIDIAGONAL 

END MODULE EIGEN 
