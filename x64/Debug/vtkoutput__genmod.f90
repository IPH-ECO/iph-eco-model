        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:42:54 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VTKOUTPUT__genmod
          INTERFACE 
            SUBROUTINE VTKOUTPUT(SIMPARAM,HYDROPARAM,MESHPARAM,         &
     &LIMNOPARAM)
              USE LIMNOLOGYVARS
              USE MESHVARS
              USE HYDRODYNAMIC
              USE SIMULATIONMODEL
              TYPE (SIMULATIONPARAM) :: SIMPARAM
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
            END SUBROUTINE VTKOUTPUT
          END INTERFACE 
        END MODULE VTKOUTPUT__genmod
