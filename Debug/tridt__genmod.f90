        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:49:53 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRIDT__genmod
          INTERFACE 
            SUBROUTINE TRIDT(DIM,MAT,D,SOL)
              INTEGER(KIND=4), INTENT(IN) :: DIM
              REAL(KIND=4), INTENT(IN) :: MAT(3,DIM)
              REAL(KIND=4), INTENT(IN) :: D(DIM)
              REAL(KIND=4), INTENT(OUT) :: SOL(DIM)
            END SUBROUTINE TRIDT
          END INTERFACE 
        END MODULE TRIDT__genmod
