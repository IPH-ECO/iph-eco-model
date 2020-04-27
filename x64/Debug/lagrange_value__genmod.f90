        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:41:56 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LAGRANGE_VALUE__genmod
          INTERFACE 
            SUBROUTINE LAGRANGE_VALUE(DATA_NUM,T_DATA,INTERP_NUM,       &
     &T_INTERP,L_INTERP)
              INTEGER(KIND=4) :: INTERP_NUM
              INTEGER(KIND=4) :: DATA_NUM
              REAL(KIND=8) :: T_DATA(DATA_NUM)
              REAL(KIND=8) :: T_INTERP(INTERP_NUM)
              REAL(KIND=8) :: L_INTERP(DATA_NUM,INTERP_NUM)
            END SUBROUTINE LAGRANGE_VALUE
          END INTERFACE 
        END MODULE LAGRANGE_VALUE__genmod
