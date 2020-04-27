        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:13 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZOOBENTHOS__genmod
          INTERFACE 
            SUBROUTINE ZOOBENTHOS(HYDROPARAM,MESHPARAM,LIMNOPARAM,DT,   &
     &DTDAY)
              USE LIMNOLOGYVARS
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: DTDAY
            END SUBROUTINE ZOOBENTHOS
          END INTERFACE 
        END MODULE ZOOBENTHOS__genmod
