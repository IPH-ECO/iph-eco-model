        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:51:16 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RATIOS__genmod
          INTERFACE 
            SUBROUTINE RATIOS(HYDROPARAM,MESHPARAM,LIMNOPARAM)
              USE LIMNOLOGYVARS
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
            END SUBROUTINE RATIOS
          END INTERFACE 
        END MODULE RATIOS__genmod
