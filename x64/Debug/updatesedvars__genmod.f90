        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:18 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATESEDVARS__genmod
          INTERFACE 
            SUBROUTINE UPDATESEDVARS(HYDROPARAM,MESHPARAM,LIMNOPARAM)
              USE LIMNOLOGYVARS
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
            END SUBROUTINE UPDATESEDVARS
          END INTERFACE 
        END MODULE UPDATESEDVARS__genmod
