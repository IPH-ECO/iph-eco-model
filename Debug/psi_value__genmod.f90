        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:49:48 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PSI_VALUE__genmod
          INTERFACE 
            SUBROUTINE PSI_VALUE(PSI_FLAG,R_FACE,COURANT,FI_SMALL,      &
     &PSI_FACE)
              INTEGER(KIND=4), INTENT(IN) :: PSI_FLAG
              REAL(KIND=4), INTENT(IN) :: R_FACE
              REAL(KIND=4), INTENT(IN) :: COURANT
              REAL(KIND=4), INTENT(IN) :: FI_SMALL
              REAL(KIND=4), INTENT(OUT) :: PSI_FACE
            END SUBROUTINE PSI_VALUE
          END INTERFACE 
        END MODULE PSI_VALUE__genmod
