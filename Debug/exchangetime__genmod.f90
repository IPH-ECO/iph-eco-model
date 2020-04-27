        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:50:47 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EXCHANGETIME__genmod
          INTERFACE 
            SUBROUTINE EXCHANGETIME(HYDROPARAM,MESHPARAM,METEOPARAM,DT)
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              REAL(KIND=4) :: DT
            END SUBROUTINE EXCHANGETIME
          END INTERFACE 
        END MODULE EXCHANGETIME__genmod
