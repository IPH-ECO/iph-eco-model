        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:42 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WATERTEMP__genmod
          INTERFACE 
            SUBROUTINE WATERTEMP(HYDROPARAM,MESHPARAM,METEOPARAM,       &
     &LIMNOPARAM,DT,DTDAY)
              USE LIMNOLOGYVARS
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: DTDAY
            END SUBROUTINE WATERTEMP
          END INTERFACE 
        END MODULE WATERTEMP__genmod
