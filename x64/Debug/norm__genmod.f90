        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:41:58 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NORM__genmod
          INTERFACE 
            SUBROUTINE NORM(A,DIM,RES)
              INTEGER(KIND=4), INTENT(IN) :: DIM
              REAL(KIND=8), INTENT(IN) :: A(DIM)
              REAL(KIND=8), INTENT(INOUT) :: RES
            END SUBROUTINE NORM
          END INTERFACE 
        END MODULE NORM__genmod
