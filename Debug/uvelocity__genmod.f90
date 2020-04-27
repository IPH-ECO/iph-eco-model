        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:52:37 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UVELOCITY__genmod
          INTERFACE 
            SUBROUTINE UVELOCITY(HYDROPARAM,MESHPARAM,DT)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              REAL(KIND=4) :: DT
            END SUBROUTINE UVELOCITY
          END INTERFACE 
        END MODULE UVELOCITY__genmod
