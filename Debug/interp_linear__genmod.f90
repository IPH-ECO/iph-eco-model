        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:49:48 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INTERP_LINEAR__genmod
          INTERFACE 
            SUBROUTINE INTERP_LINEAR(M,DATA_NUM,T_DATA,P_DATA,INTERP_NUM&
     &,T_INTERP,P_INTERP)
              INTEGER(KIND=4) :: INTERP_NUM
              INTEGER(KIND=4) :: DATA_NUM
              INTEGER(KIND=4) :: M
              REAL(KIND=4) :: T_DATA(DATA_NUM)
              REAL(KIND=4) :: P_DATA(M,DATA_NUM)
              REAL(KIND=4) :: T_INTERP(INTERP_NUM)
              REAL(KIND=4) :: P_INTERP(M,INTERP_NUM)
            END SUBROUTINE INTERP_LINEAR
          END INTERFACE 
        END MODULE INTERP_LINEAR__genmod
