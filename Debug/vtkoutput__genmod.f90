        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:50:24 2020
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
