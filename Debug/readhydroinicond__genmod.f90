        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:52:12 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READHYDROINICOND__genmod
          INTERFACE 
            SUBROUTINE READHYDROINICOND(HYDROPARAM,HYDROCONFIGURATION,  &
     &SIMPARAM,MESHPARAM)
              USE MESHVARS
              USE SIMULATIONMODEL
              USE DOMAIN_TYPES
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (HYDRODYNAMICCONFIGURATION) :: HYDROCONFIGURATION
              TYPE (SIMULATIONPARAM) :: SIMPARAM
            END SUBROUTINE READHYDROINICOND
          END INTERFACE 
        END MODULE READHYDROINICOND__genmod
