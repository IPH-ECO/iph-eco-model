        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:23 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TWOEQ_KC__genmod
          INTERFACE 
            SUBROUTINE TWOEQ_KC(HYDROPARAM,MESHPARAM,METEOPARAM,DT)
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              REAL(KIND=8) :: DT
            END SUBROUTINE TWOEQ_KC
          END INTERFACE 
        END MODULE TWOEQ_KC__genmod
