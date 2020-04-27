        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:03 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ALLOCATEMETEOVARS__genmod
          INTERFACE 
            SUBROUTINE ALLOCATEMETEOVARS(MESHPARAM,METEOPARAM)
              USE METEOROLOGICAL
              USE MESHVARS
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
            END SUBROUTINE ALLOCATEMETEOVARS
          END INTERFACE 
        END MODULE ALLOCATEMETEOVARS__genmod
