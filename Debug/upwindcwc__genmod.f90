        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:52:34 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPWINDCWC__genmod
          INTERFACE 
            SUBROUTINE UPWINDCWC(INDEX,ULOADVAR,DVAR,VAREST,VARESTP,DT, &
     &DTDAY,HYDROPARAM,MESHPARAM)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              INTEGER(KIND=4) :: INDEX
              REAL(KIND=4) :: ULOADVAR(MESHPARAM%KMAX,MESHPARAM%NELEM)
              REAL(KIND=4) :: DVAR(MESHPARAM%KMAX,MESHPARAM%NELEM)
              REAL(KIND=4) :: VAREST(MESHPARAM%KMAX,MESHPARAM%NELEM)
              REAL(KIND=4) :: VARESTP(MESHPARAM%KMAX,MESHPARAM%NELEM)
              REAL(KIND=4) :: DT
              REAL(KIND=4) :: DTDAY
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
            END SUBROUTINE UPWINDCWC
          END INTERFACE 
        END MODULE UPWINDCWC__genmod
