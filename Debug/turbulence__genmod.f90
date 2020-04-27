        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:52:14 2020
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
              REAL(KIND=4) :: DT
            END SUBROUTINE TURBULENCE
          END INTERFACE 
        END MODULE TURBULENCE__genmod
