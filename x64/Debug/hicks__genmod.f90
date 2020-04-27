        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:42 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HICKS__genmod
          INTERFACE 
            SUBROUTINE HICKS(CDN,CHEN,SFRT,SFRRHO,SPHAIR,SPHWTR,AIRDENS,&
     &LATHEAT,TV,WS,AIRT,G,PI,CD,CHE,LMO)
              REAL(KIND=8), INTENT(IN) :: CDN
              REAL(KIND=8), INTENT(IN) :: CHEN
              REAL(KIND=8), INTENT(IN) :: SFRT
              REAL(KIND=8), INTENT(IN) :: SFRRHO
              REAL(KIND=8), INTENT(IN) :: SPHAIR
              REAL(KIND=8), INTENT(IN) :: SPHWTR
              REAL(KIND=8), INTENT(IN) :: AIRDENS
              REAL(KIND=8), INTENT(IN) :: LATHEAT
              REAL(KIND=8), INTENT(IN) :: TV
              REAL(KIND=8), INTENT(IN) :: WS
              REAL(KIND=8), INTENT(IN) :: AIRT
              REAL(KIND=8) :: G
              REAL(KIND=8) :: PI
              REAL(KIND=8), INTENT(OUT) :: CD
              REAL(KIND=8), INTENT(OUT) :: CHE
              REAL(KIND=8), INTENT(OUT) :: LMO
            END SUBROUTINE HICKS
          END INTERFACE 
        END MODULE HICKS__genmod
