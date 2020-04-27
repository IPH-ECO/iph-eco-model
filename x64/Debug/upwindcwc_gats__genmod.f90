        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:20 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPWINDCWC_GATS__genmod
          INTERFACE 
            SUBROUTINE UPWINDCWC_GATS(INDEX,ULOADVAR,DVAR,VAREST,VARESTP&
     &,DT,DTDAY,HYDROPARAM,MESHPARAM)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              INTEGER(KIND=4) :: INDEX
              REAL(KIND=8) :: ULOADVAR(MESHPARAM%KMAX,MESHPARAM%NELEM)
              REAL(KIND=8) :: DVAR(MESHPARAM%KMAX,MESHPARAM%NELEM)
              REAL(KIND=8) :: VAREST(MESHPARAM%KMAX,MESHPARAM%NELEM)
              REAL(KIND=8) :: VARESTP(MESHPARAM%KMAX,MESHPARAM%NELEM)
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: DTDAY
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
            END SUBROUTINE UPWINDCWC_GATS
          END INTERFACE 
        END MODULE UPWINDCWC_GATS__genmod
