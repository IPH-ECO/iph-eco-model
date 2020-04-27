        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:13 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPWINDCWC_HR__genmod
          INTERFACE 
            SUBROUTINE UPWINDCWC_HR(INDEX,ULOADVAR,DVAR,VAREST,VARESTP, &
     &DT,DTDAY,HYDROPARAM,MESHPARAM,PSI_FLAG)
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
              INTEGER(KIND=4) :: PSI_FLAG
            END SUBROUTINE UPWINDCWC_HR
          END INTERFACE 
        END MODULE UPWINDCWC_HR__genmod
