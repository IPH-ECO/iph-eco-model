        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:51:23 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ALLOCATEHYDROVARS__genmod
          INTERFACE 
            SUBROUTINE ALLOCATEHYDROVARS(HYDROPARAM,MESHPARAM)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
            END SUBROUTINE ALLOCATEHYDROVARS
          END INTERFACE 
        END MODULE ALLOCATEHYDROVARS__genmod
