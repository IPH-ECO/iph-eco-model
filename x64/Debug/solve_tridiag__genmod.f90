        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:41:57 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_TRIDIAG__genmod
          INTERFACE 
            SUBROUTINE SOLVE_TRIDIAG(A,B,C,V,X,N)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: A(N)
              REAL(KIND=8), INTENT(IN) :: B(N)
              REAL(KIND=8), INTENT(IN) :: C(N)
              REAL(KIND=8), INTENT(IN) :: V(N)
              REAL(KIND=8), INTENT(OUT) :: X(N)
            END SUBROUTINE SOLVE_TRIDIAG
          END INTERFACE 
        END MODULE SOLVE_TRIDIAG__genmod
