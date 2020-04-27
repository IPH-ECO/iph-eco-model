        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:56 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATELIMNVARS__genmod
          INTERFACE 
            SUBROUTINE UPDATELIMNVARS(HYDROPARAM,MESHPARAM,LIMNOPARAM,  &
     &SIMPARAM)
              USE SIMULATIONMODEL
              USE LIMNOLOGYVARS
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              TYPE (SIMULATIONPARAM) :: SIMPARAM
            END SUBROUTINE UPDATELIMNVARS
          END INTERFACE 
        END MODULE UPDATELIMNVARS__genmod
