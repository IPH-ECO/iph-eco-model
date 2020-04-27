        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:23 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TURBULENCE__genmod
          INTERFACE 
            SUBROUTINE TURBULENCE(HYDROPARAM,MESHPARAM,METEOPARAM,DT)
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              REAL(KIND=8) :: DT
            END SUBROUTINE TURBULENCE
          END INTERFACE 
        END MODULE TURBULENCE__genmod
