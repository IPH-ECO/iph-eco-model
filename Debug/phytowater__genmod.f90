        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:50:56 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHYTOWATER__genmod
          INTERFACE 
            SUBROUTINE PHYTOWATER(HYDROPARAM,MESHPARAM,METEOPARAM,      &
     &LIMNOPARAM,DT,DTDAY)
              USE LIMNOLOGYVARS
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              REAL(KIND=4) :: DT
              REAL(KIND=4) :: DTDAY
            END SUBROUTINE PHYTOWATER
          END INTERFACE 
        END MODULE PHYTOWATER__genmod
