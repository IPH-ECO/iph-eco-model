        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:41:56 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PSI_VALUE__genmod
          INTERFACE 
            SUBROUTINE PSI_VALUE(PSI_FLAG,R_FACE,COURANT,FI_SMALL,      &
     &PSI_FACE)
              INTEGER(KIND=4), INTENT(IN) :: PSI_FLAG
              REAL(KIND=8), INTENT(IN) :: R_FACE
              REAL(KIND=8), INTENT(IN) :: COURANT
              REAL(KIND=8), INTENT(IN) :: FI_SMALL
              REAL(KIND=8), INTENT(OUT) :: PSI_FACE
            END SUBROUTINE PSI_VALUE
          END INTERFACE 
        END MODULE PSI_VALUE__genmod
