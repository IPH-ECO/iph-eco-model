        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:50:49 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETWQBOUNDARYCONDITIONS__genmod
          INTERFACE 
            SUBROUTINE GETWQBOUNDARYCONDITIONS(HYDROPARAM,LIMNOPARAM,   &
     &TIME)
              USE LIMNOLOGYVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              REAL(KIND=4) :: TIME
            END SUBROUTINE GETWQBOUNDARYCONDITIONS
          END INTERFACE 
        END MODULE GETWQBOUNDARYCONDITIONS__genmod
