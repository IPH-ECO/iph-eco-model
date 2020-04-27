        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:24 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CGOPNH__genmod
          INTERFACE 
            SUBROUTINE CGOPNH(A,B,DT,HSUBLAYER,HYDROPARAM,MESHPARAM)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              REAL(KIND=8), INTENT(IN) :: A(MESHPARAM%KMAX+1,MESHPARAM% &
     &NELEM)
              REAL(KIND=8), INTENT(OUT) :: B(MESHPARAM%KMAX+1,MESHPARAM%&
     &NELEM)
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: HSUBLAYER
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
            END SUBROUTINE CGOPNH
          END INTERFACE 
        END MODULE CGOPNH__genmod
