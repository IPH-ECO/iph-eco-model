        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:41:57 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DGTSV__genmod
          INTERFACE 
            SUBROUTINE DGTSV(N,NRHS,DL,D,DU,B,LDB,INFO)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NRHS
              REAL(KIND=8) :: DL(*)
              REAL(KIND=8) :: D(*)
              REAL(KIND=8) :: DU(*)
              REAL(KIND=8) :: B(LDB,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGTSV
          END INTERFACE 
        END MODULE DGTSV__genmod
