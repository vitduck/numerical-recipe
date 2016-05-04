MODULE EIGEN 

CONTAINS 

SUBROUTINE TRIDIAGONAL(A, B, C, X, R) 
    IMPLICIT NONE 

    REAL, INTENT(IN)         :: A(:), B(:), C(:), R(:)
    REAL, INTENT(OUT)        :: X(:) 

    REAL                     :: beta(SIZE(B(:))) 
    INTEGER                  :: i

    ! initialization 
    beta(1) = B(1) 
    X(1)  = R(1) 

    ! Forward elimination of x_i  
    DO i = 2, n 
        beta(i) = B(i) - A(i)*C(i-1)/beta(i-1)
        X(i)    = R(i) - A(i)*X(i-1)/beta(i-1)
    END DO 

    ! Backward substitution 
    X(n) = X(n)/beta(n) 
    DO i = n-1, 1, -1 
        X(i) = (X(i) - C(i)*X(i+1))/beta(i)
    END DO 
END SUBROUTINE TRIDIAGONAL 

END MODULE EIGEN 
