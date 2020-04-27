        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:40 2020
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
