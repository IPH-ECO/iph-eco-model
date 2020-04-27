        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:51:21 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READOUTPUTS__genmod
          INTERFACE 
            SUBROUTINE READOUTPUTS(SIM,SIMPARAM,KMAX,NEDGE,NELEM)
              USE SIMULATIONMODEL
              USE DOMAIN_TYPES
              TYPE (SIMULATION) :: SIM
              TYPE (SIMULATIONPARAM) :: SIMPARAM
              INTEGER(KIND=4) :: KMAX
              INTEGER(KIND=4) :: NEDGE
              INTEGER(KIND=4) :: NELEM
            END SUBROUTINE READOUTPUTS
          END INTERFACE 
        END MODULE READOUTPUTS__genmod
