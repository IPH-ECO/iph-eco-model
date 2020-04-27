        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:51:33 2020
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
              REAL(KIND=4) :: ULOADVAR(MESHPARAM%KMAX,MESHPARAM%NELEM)
              REAL(KIND=4) :: DVAR(MESHPARAM%KMAX,MESHPARAM%NELEM)
              REAL(KIND=4) :: VAREST(MESHPARAM%KMAX,MESHPARAM%NELEM)
              REAL(KIND=4) :: VARESTP(MESHPARAM%KMAX,MESHPARAM%NELEM)
              REAL(KIND=4) :: DT
              REAL(KIND=4) :: DTDAY
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              INTEGER(KIND=4) :: PSI_FLAG
            END SUBROUTINE UPWINDCWC_HR
          END INTERFACE 
        END MODULE UPWINDCWC_HR__genmod
