        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:57 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READWQINICOND__genmod
          INTERFACE 
            SUBROUTINE READWQINICOND(HYDROPARAM,MESHPARAM,LIMNOPARAM,   &
     &SIMPARAM)
              USE SIMULATIONMODEL
              USE LIMNOLOGYVARS
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              TYPE (SIMULATIONPARAM) :: SIMPARAM
            END SUBROUTINE READWQINICOND
          END INTERFACE 
        END MODULE READWQINICOND__genmod
