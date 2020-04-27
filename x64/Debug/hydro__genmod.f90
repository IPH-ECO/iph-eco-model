        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:49:54 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HYDRO__genmod
          INTERFACE 
            SUBROUTINE HYDRO(HYDROPARAM,MESHPARAM,METEOPARAM,DT,TIME,   &
     &SIMTIME)
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: TIME
              REAL(KIND=8) :: SIMTIME
            END SUBROUTINE HYDRO
          END INTERFACE 
        END MODULE HYDRO__genmod
