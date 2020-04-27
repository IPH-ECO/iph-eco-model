        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:51:10 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EXPLICITTERMS__genmod
          INTERFACE 
            SUBROUTINE EXPLICITTERMS(HYDROPARAM,MESHPARAM,DT)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              REAL(KIND=4) :: DT
            END SUBROUTINE EXPLICITTERMS
          END INTERFACE 
        END MODULE EXPLICITTERMS__genmod
