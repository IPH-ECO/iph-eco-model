        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:41:58 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TENSION__genmod
          INTERFACE 
            SUBROUTINE TENSION(SURFTENSIONFLAG,BOTTOMTENSIONFLAG,IEDGE, &
     &KSURF,KBOTTOM,G,UB,VB,US,VS,CHEZY,RHOAIR,RHO0,CV,WINDEDGE,        &
     &WINDEDGEXY,GAMMAT,GAMMAB)
              INTEGER(KIND=4) :: SURFTENSIONFLAG
              INTEGER(KIND=4) :: BOTTOMTENSIONFLAG
              INTEGER(KIND=4) :: IEDGE
              INTEGER(KIND=4) :: KSURF
              INTEGER(KIND=4) :: KBOTTOM
              REAL(KIND=8) :: G
              REAL(KIND=8) :: UB
              REAL(KIND=8) :: VB
              REAL(KIND=8) :: US
              REAL(KIND=8) :: VS
              REAL(KIND=8) :: CHEZY
              REAL(KIND=8) :: RHOAIR
              REAL(KIND=8) :: RHO0
              REAL(KIND=8) :: CV
              REAL(KIND=8) :: WINDEDGE(2)
              REAL(KIND=8) :: WINDEDGEXY(2)
              REAL(KIND=8) :: GAMMAT
              REAL(KIND=8) :: GAMMAB
            END SUBROUTINE TENSION
          END INTERFACE 
        END MODULE TENSION__genmod
