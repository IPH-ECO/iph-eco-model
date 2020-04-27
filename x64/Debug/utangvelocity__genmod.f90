        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:05 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UTANGVELOCITY__genmod
          INTERFACE 
            SUBROUTINE UTANGVELOCITY(HYDROPARAM,MESHPARAM,DT)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              REAL(KIND=8) :: DT
            END SUBROUTINE UTANGVELOCITY
          END INTERFACE 
        END MODULE UTANGVELOCITY__genmod
