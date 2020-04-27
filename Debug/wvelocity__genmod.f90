        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:51:17 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WVELOCITY__genmod
          INTERFACE 
            SUBROUTINE WVELOCITY(HYDROPARAM,MESHPARAM,DT)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              REAL(KIND=4) :: DT
            END SUBROUTINE WVELOCITY
          END INTERFACE 
        END MODULE WVELOCITY__genmod
