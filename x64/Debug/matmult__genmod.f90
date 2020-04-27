        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:41:55 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MATMULT__genmod
          INTERFACE 
            SUBROUTINE MATMULT(A,B,C,N,M)
              INTEGER(KIND=4), INTENT(IN) :: M
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: A(N,N)
              REAL(KIND=8), INTENT(IN) :: B(N,M)
              REAL(KIND=8), INTENT(OUT) :: C(N,M)
            END SUBROUTINE MATMULT
          END INTERFACE 
        END MODULE MATMULT__genmod
