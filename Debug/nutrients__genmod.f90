        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:51:18 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NUTRIENTS__genmod
          INTERFACE 
            SUBROUTINE NUTRIENTS(HYDROPARAM,MESHPARAM,LIMNOPARAM,DT,    &
     &DTDAY)
              USE LIMNOLOGYVARS
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              REAL(KIND=4) :: DT
              REAL(KIND=4) :: DTDAY
            END SUBROUTINE NUTRIENTS
          END INTERFACE 
        END MODULE NUTRIENTS__genmod
