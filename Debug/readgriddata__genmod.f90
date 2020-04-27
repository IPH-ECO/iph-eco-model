        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:51:09 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READGRIDDATA__genmod
          INTERFACE 
            SUBROUTINE READGRIDDATA(GRIDDATACONFIG,MESHPARAM,HYDROPARAM)
              USE HYDRODYNAMIC
              USE MESHVARS
              USE DOMAIN_TYPES
              TYPE (GRIDDATACONFIGURATION) :: GRIDDATACONFIG
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
            END SUBROUTINE READGRIDDATA
          END INTERFACE 
        END MODULE READGRIDDATA__genmod
