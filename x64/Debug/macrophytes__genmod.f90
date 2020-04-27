        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:16 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MACROPHYTES__genmod
          INTERFACE 
            SUBROUTINE MACROPHYTES(HYDROPARAM,MESHPARAM,METEOPARAM,     &
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
            END SUBROUTINE MACROPHYTES
          END INTERFACE 
        END MODULE MACROPHYTES__genmod
