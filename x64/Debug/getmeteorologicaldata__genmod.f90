        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:02 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETMETEOROLOGICALDATA__genmod
          INTERFACE 
            SUBROUTINE GETMETEOROLOGICALDATA(MESHPARAM,METEOPARAM,TIME, &
     &AIRTEMPREF)
              USE METEOROLOGICAL
              USE MESHVARS
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              REAL(KIND=8) :: TIME
              REAL(KIND=8) :: AIRTEMPREF
            END SUBROUTINE GETMETEOROLOGICALDATA
          END INTERFACE 
        END MODULE GETMETEOROLOGICALDATA__genmod
