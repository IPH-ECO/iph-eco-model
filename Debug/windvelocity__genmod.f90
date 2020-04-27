        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:52:36 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WINDVELOCITY__genmod
          INTERFACE 
            SUBROUTINE WINDVELOCITY(HYDROPARAM,MESHPARAM,METEOPARAM)
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
            END SUBROUTINE WINDVELOCITY
          END INTERFACE 
        END MODULE WINDVELOCITY__genmod
