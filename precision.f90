MODULE PRECISION
    IMPLICIT NONE 

    INTRINSIC KIND 
    
    ! standard IEEE arithmetic 
    !INTEGER, PARAMETER :: skind = SELECT_REAL_KIND(p=6, r=37)
    !INTEGER, PARAMETER :: dkind = SELECT_REAL_KIND(p=15, r=307)

    ! standard precision set by compiler 
    INTEGER, PARAMETER :: skind = KIND(0.0E0)
    INTEGER, PARAMETER :: dkind = KIND(0.0D0)
    
    ! Default precision of module
    INTEGER, PARAMETER :: wp = dkind 
END MODULE
