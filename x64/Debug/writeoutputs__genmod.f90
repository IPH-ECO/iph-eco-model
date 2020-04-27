        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:18 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITEOUTPUTS__genmod
          INTERFACE 
            SUBROUTINE WRITEOUTPUTS(SIMPARAM,HYDROPARAM,MESHPARAM,      &
     &LIMNOPARAM,METEOPARAM)
              USE METEOROLOGICAL
              USE LIMNOLOGYVARS
              USE MESHVARS
              USE HYDRODYNAMIC
              USE SIMULATIONMODEL
              TYPE (SIMULATIONPARAM) :: SIMPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
            END SUBROUTINE WRITEOUTPUTS
          END INTERFACE 
        END MODULE WRITEOUTPUTS__genmod
