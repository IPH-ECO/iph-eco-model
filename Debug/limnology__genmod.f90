        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:52:08 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LIMNOLOGY__genmod
          INTERFACE 
            SUBROUTINE LIMNOLOGY(HYDROPARAM,MESHPARAM,METEOPARAM,       &
     &LIMNOPARAM,SIMPARAM)
              USE SIMULATIONMODEL
              USE LIMNOLOGYVARS
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              TYPE (SIMULATIONPARAM) :: SIMPARAM
            END SUBROUTINE LIMNOLOGY
          END INTERFACE 
        END MODULE LIMNOLOGY__genmod
