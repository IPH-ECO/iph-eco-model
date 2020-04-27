        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:51:25 2020
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
              REAL(KIND=4) :: DT
              REAL(KIND=4) :: DTDAY
              INTEGER(KIND=4) :: JULDAY
            END SUBROUTINE FISHES
          END INTERFACE 
        END MODULE FISHES__genmod
