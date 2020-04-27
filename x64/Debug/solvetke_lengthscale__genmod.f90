        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:22 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVETKE_LENGTHSCALE__genmod
          INTERFACE 
            SUBROUTINE SOLVETKE_LENGTHSCALE(HYDROPARAM,MESHPARAM,       &
     &METEOPARAM,DT,Q2PP,Q2LPP,Q2P,Q2LP,Q2,Q2L)
              USE METEOROLOGICAL
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (METEOROLOGICALPARAM) :: METEOPARAM
              REAL(KIND=8) :: DT
              REAL(KIND=8), INTENT(IN) :: Q2PP(MESHPARAM%KMAX+1,        &
     &MESHPARAM%NELEM)
              REAL(KIND=8), INTENT(IN) :: Q2LPP(MESHPARAM%KMAX+1,       &
     &MESHPARAM%NELEM)
              REAL(KIND=8), INTENT(IN) :: Q2P(MESHPARAM%KMAX+1,MESHPARAM&
     &%NELEM)
              REAL(KIND=8), INTENT(IN) :: Q2LP(MESHPARAM%KMAX+1,        &
     &MESHPARAM%NELEM)
              REAL(KIND=8), INTENT(OUT) :: Q2(MESHPARAM%KMAX+1,MESHPARAM&
     &%NELEM)
              REAL(KIND=8), INTENT(OUT) :: Q2L(MESHPARAM%KMAX+1,        &
     &MESHPARAM%NELEM)
            END SUBROUTINE SOLVETKE_LENGTHSCALE
          END INTERFACE 
        END MODULE SOLVETKE_LENGTHSCALE__genmod
