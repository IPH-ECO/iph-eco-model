        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:52:17 2020
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
              REAL(KIND=4) :: DT
              REAL(KIND=4) :: DTDAY
            END SUBROUTINE WATERTEMP
          END INTERFACE 
        END MODULE WATERTEMP__genmod
