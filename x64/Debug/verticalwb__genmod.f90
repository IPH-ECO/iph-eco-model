        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:55 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VERTICALWB__genmod
          INTERFACE 
            SUBROUTINE VERTICALWB(HYDROPARAM,MESHPARAM,METEOPARAM,DT)
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              REAL(KIND=8) :: DT
            END SUBROUTINE VERTICALWB
          END INTERFACE 
        END MODULE VERTICALWB__genmod
