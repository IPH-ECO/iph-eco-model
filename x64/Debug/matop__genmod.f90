        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:30 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MATOP__genmod
          INTERFACE 
            SUBROUTINE MATOP(A,B,DT,HYDROPARAM,MESHPARAM)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              REAL(KIND=8), INTENT(IN) :: A(MESHPARAM%NELEM)
              REAL(KIND=8), INTENT(OUT) :: B(MESHPARAM%NELEM)
              REAL(KIND=8) :: DT
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
            END SUBROUTINE MATOP
          END INTERFACE 
        END MODULE MATOP__genmod
