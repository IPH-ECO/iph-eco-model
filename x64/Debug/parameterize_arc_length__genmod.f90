        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:41:56 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PARAMETERIZE_ARC_LENGTH__genmod
          INTERFACE 
            SUBROUTINE PARAMETERIZE_ARC_LENGTH(M,DATA_NUM,P_DATA,T_DATA)
              INTEGER(KIND=4) :: DATA_NUM
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: P_DATA(M,DATA_NUM)
              REAL(KIND=8) :: T_DATA(DATA_NUM)
            END SUBROUTINE PARAMETERIZE_ARC_LENGTH
          END INTERFACE 
        END MODULE PARAMETERIZE_ARC_LENGTH__genmod
