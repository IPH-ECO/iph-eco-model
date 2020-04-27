        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:15 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EXPLICITTERMS__genmod
          INTERFACE 
            SUBROUTINE EXPLICITTERMS(HYDROPARAM,MESHPARAM,DT)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              REAL(KIND=8) :: DT
            END SUBROUTINE EXPLICITTERMS
          END INTERFACE 
        END MODULE EXPLICITTERMS__genmod
