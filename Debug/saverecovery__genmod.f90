        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:52:09 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SAVERECOVERY__genmod
          INTERFACE 
            SUBROUTINE SAVERECOVERY(SIM,SIMPARAM,MESHPARAM,HYDROPARAM,  &
     &LIMNOPARAM)
              USE LIMNOLOGYVARS
              USE HYDRODYNAMIC
              USE MESHVARS
              USE SIMULATIONMODEL
              USE DOMAIN_TYPES
              TYPE (SIMULATION) :: SIM
              TYPE (SIMULATIONPARAM) :: SIMPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
            END SUBROUTINE SAVERECOVERY
          END INTERFACE 
        END MODULE SAVERECOVERY__genmod
