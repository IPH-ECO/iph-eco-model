        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:46 2020
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
              REAL(KIND=8) :: TIME
            END SUBROUTINE GETWQBOUNDARYCONDITIONS
          END INTERFACE 
        END MODULE GETWQBOUNDARYCONDITIONS__genmod
