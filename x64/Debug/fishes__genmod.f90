        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:24 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FISHES__genmod
          INTERFACE 
            SUBROUTINE FISHES(HYDROPARAM,MESHPARAM,LIMNOPARAM,DT,DTDAY, &
     &JULDAY)
              USE LIMNOLOGYVARS
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: DTDAY
              INTEGER(KIND=4) :: JULDAY
            END SUBROUTINE FISHES
          END INTERFACE 
        END MODULE FISHES__genmod
