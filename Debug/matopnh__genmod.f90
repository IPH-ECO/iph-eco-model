        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:50:52 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MATOPNH__genmod
          INTERFACE 
            SUBROUTINE MATOPNH(A,B,DT,HSUBLAYER,HYDROPARAM,MESHPARAM)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              REAL(KIND=4), INTENT(IN) :: A(MESHPARAM%KMAX+1,MESHPARAM% &
     &NELEM)
              REAL(KIND=4), INTENT(OUT) :: B(MESHPARAM%KMAX+1,MESHPARAM%&
     &NELEM)
              REAL(KIND=4) :: DT
              REAL(KIND=4) :: HSUBLAYER
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
            END SUBROUTINE MATOPNH
          END INTERFACE 
        END MODULE MATOPNH__genmod
