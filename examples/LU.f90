PROGRAM solve_linear
    USE PRECISION, ONLY: wp
    
    IMPLICIT NONE 

    REAL(wp), ALLOCATABLE :: a(:,:), b(:,:), y(:)
    REAL(wp), ALLOCATABLE :: ipvt(:)
    REAL(wp)              :: relerr, dnrm2
    REAL(wp)              :: t_start, t_end 
    REAL(wp), PARAMETER   :: one = 1.0_wp, zero = 0.0_wp

    INTEGER               :: i, j
    INTEGER               :: info 
    INTEGER, PARAMETER    :: n = 4000

    CHARACTER(80)         :: msg

    ! allocate array 
    ALLOCATE(a(n,n), b(n,1), y(n), ipvt(n), STAT=info, ERRMSG=msg)

    IF ( info /= 0 ) THEN 
        WRITE(*, '("Error in array allocation: ", A)') msg
        STOP
    END IF

    ! fill matrices with uniformed random numbers
    CALL RANDOM_NUMBER(a)
    CALL RANDOM_NUMBER(b)
    CALL RANDOM_NUMBER(y)

    ! B = A*Y (BLAS2)
    CALL DGEMV('N', n, n, one, a, n, y, 1, zero, b, 1)

    ! start measurement 
    CALL CPU_TIME(t_start)

    ! LU factorization 
    CALL DGETRF(n, n, a, n ,ipvt, info)
    
    ! DGETRF return value
    IF ( info < 0 ) THEN 
        WRITE(*, '("Argument", i3, "has an illegal value")') info
    ELSE IF ( info > 0 ) THEN 
        WRITE(*, '("Zero diagonal element at", i3, "Division by zero will occur")') info
    END IF

    ! solve linear system
    CALL DGETRS('N', n, 1, a, n, ipvt, b, n, info)

    ! DGETRF return value
    IF ( info < 0 ) THEN 
        WRITE(*, '("Argument", i3, "has an illegal value")') info
    END IF 

    CALL CPU_TIME(t_end)

    ! difference between lapack solution and initial generator 
    b(:,1) = b(:,1) - y
    relerr = DNRM2(n, b, 1)/DNRM2(n, y, 1)
    
    WRITE(*,'("Relative error: ", F15.8)') relerr
    WRITE(*,'("Elapsed time", F15.8)') t_end - t_start

    ! free memeory 
    DEALLOCATE(a, b, y, ipvt, STAT=info )

    IF (info /= 0 ) THEN 
        WRITE(*, '("Deallocation failed")')
    END IF 
END PROGRAM
