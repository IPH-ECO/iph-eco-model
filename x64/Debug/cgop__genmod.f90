        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:55 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CGOP__genmod
          INTERFACE 
            SUBROUTINE CGOP(A,B,DT,HYDROPARAM,MESHPARAM)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              REAL(KIND=8), INTENT(IN) :: A(MESHPARAM%NELEM)
              REAL(KIND=8), INTENT(OUT) :: B(MESHPARAM%NELEM)
              REAL(KIND=8) :: DT
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
            END SUBROUTINE CGOP
          END INTERFACE 
        END MODULE CGOP__genmod
