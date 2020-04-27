        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:42 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BULKTRANSFER__genmod
          INTERFACE 
            SUBROUTINE BULKTRANSFER(H,PI,L,CDN,CEN,CD,CE)
              REAL(KIND=8), INTENT(IN) :: H
              REAL(KIND=8) :: PI
              REAL(KIND=8) :: L
              REAL(KIND=8), INTENT(IN) :: CDN
              REAL(KIND=8), INTENT(IN) :: CEN
              REAL(KIND=8), INTENT(OUT) :: CD
              REAL(KIND=8), INTENT(OUT) :: CE
            END SUBROUTINE BULKTRANSFER
          END INTERFACE 
        END MODULE BULKTRANSFER__genmod
