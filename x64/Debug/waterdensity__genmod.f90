        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:36 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WATERDENSITY__genmod
          INTERFACE 
            SUBROUTINE WATERDENSITY(HYDROPARAM,MESHPARAM,LIMNOPARAM)
              USE LIMNOLOGYVARS
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
            END SUBROUTINE WATERDENSITY
          END INTERFACE 
        END MODULE WATERDENSITY__genmod
