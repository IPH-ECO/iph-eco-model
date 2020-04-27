        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:51:39 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETHYDROBOUNDARYCONDITIONS__genmod
          INTERFACE 
            SUBROUTINE GETHYDROBOUNDARYCONDITIONS(HYDROPARAM,MESHPARAM, &
     &TIME)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              REAL(KIND=4) :: TIME
            END SUBROUTINE GETHYDROBOUNDARYCONDITIONS
          END INTERFACE 
        END MODULE GETHYDROBOUNDARYCONDITIONS__genmod
