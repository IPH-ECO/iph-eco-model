        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:48 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE STARTSIMULATION__genmod
          INTERFACE 
            SUBROUTINE STARTSIMULATION(SIM) BIND(C, NAME = '            &
     &startSimulation')
              USE DOMAIN_TYPES
              TYPE (SIMULATION) :: SIM
            END SUBROUTINE STARTSIMULATION
          END INTERFACE 
        END MODULE STARTSIMULATION__genmod
