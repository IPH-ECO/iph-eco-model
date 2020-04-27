        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:42 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WTEMP__genmod
          INTERFACE 
            SUBROUTINE WTEMP(HYDROPARAM,MESHPARAM,METEOPARAM,LIMNOPARAM,&
     &SOURCE)
              USE LIMNOLOGYVARS
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              REAL(KIND=8) :: SOURCE(MESHPARAM%KMAX,MESHPARAM%NELEM)
            END SUBROUTINE WTEMP
          END INTERFACE 
        END MODULE WTEMP__genmod
