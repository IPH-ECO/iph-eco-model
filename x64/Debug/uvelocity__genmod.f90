        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:12 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UVELOCITY__genmod
          INTERFACE 
            SUBROUTINE UVELOCITY(HYDROPARAM,MESHPARAM,DT)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              REAL(KIND=8) :: DT
            END SUBROUTINE UVELOCITY
          END INTERFACE 
        END MODULE UVELOCITY__genmod
