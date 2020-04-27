        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:05 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SEDIMENT__genmod
          INTERFACE 
            SUBROUTINE SEDIMENT(HYDROPARAM,MESHPARAM,METEOPARAM,        &
     &LIMNOPARAM,DT,DTDAY)
              USE LIMNOLOGYVARS
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: DTDAY
            END SUBROUTINE SEDIMENT
          END INTERFACE 
        END MODULE SEDIMENT__genmod
