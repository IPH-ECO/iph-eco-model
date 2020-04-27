        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:51:56 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VELOCITIES__genmod
          INTERFACE 
            SUBROUTINE VELOCITIES(HYDROPARAM,MESHPARAM)
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
            END SUBROUTINE VELOCITIES
          END INTERFACE 
        END MODULE VELOCITIES__genmod
