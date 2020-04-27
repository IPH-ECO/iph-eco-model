        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:34 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WVELOCITY__genmod
          INTERFACE 
            SUBROUTINE WVELOCITY(HYDROPARAM,MESHPARAM,DT)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              REAL(KIND=8) :: DT
            END SUBROUTINE WVELOCITY
          END INTERFACE 
        END MODULE WVELOCITY__genmod
