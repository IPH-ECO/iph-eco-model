        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:49:50 2020
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
              REAL(KIND=4) :: G
              REAL(KIND=4) :: UB
              REAL(KIND=4) :: VB
              REAL(KIND=4) :: US
              REAL(KIND=4) :: VS
              REAL(KIND=4) :: CHEZY
              REAL(KIND=4) :: RHOAIR
              REAL(KIND=4) :: RHO0
              REAL(KIND=4) :: CV
              REAL(KIND=4) :: WINDEDGE(2)
              REAL(KIND=4) :: WINDEDGEXY(2)
              REAL(KIND=4) :: GAMMAT
              REAL(KIND=4) :: GAMMAB
            END SUBROUTINE TENSION
          END INTERFACE 
        END MODULE TENSION__genmod
