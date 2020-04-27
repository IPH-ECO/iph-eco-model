        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:17 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SALINITY__genmod
          INTERFACE 
            SUBROUTINE SALINITY(HYDROPARAM,MESHPARAM,METEOPARAM,        &
     &LIMNOPARAM,DT,DTDAY)
              USE LIMNOLOGYVARS
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: DTDAY
            END SUBROUTINE SALINITY
          END INTERFACE 
        END MODULE SALINITY__genmod
