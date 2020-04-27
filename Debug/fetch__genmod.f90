        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:50:46 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FETCH__genmod
          INTERFACE 
            SUBROUTINE FETCH(MESHPARAM,HYDROPARAM)
              USE HYDRODYNAMIC
              USE MESHVARS
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
            END SUBROUTINE FETCH
          END INTERFACE 
        END MODULE FETCH__genmod
