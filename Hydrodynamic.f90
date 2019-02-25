  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>\brief Compute the hydrodynamic simulation
  !>\details Compute the hydrodynamic simulation. It is called from main
  !>\author Carlos Ruberto Fragoso
  !>\version 1.3
  !>\attention Main modification list:\n
  !>\-# Implemetation Carlos Ruberto Fragoso
  !>\remark All variables are passed through the following modules
  !> -#  Hydrodynamic
  !> -#  MeshVars
 !>\section sec_flow_simhydro Flowchart
  !>\image html SIMUL_ONLY_HYDRO_flowchart.jpg "Main flowchart version 1.3 " width = 6cm
Subroutine HydrodynamicModel
    
    Use Hydrodynamic
    Use MeshVars
    
    Implicit none

   
    
    !2. Use Casulli's Semi-Implicit Solution (TRIM) to solve the Depth-Averaged Navier-Stokes Equation    
    Call Hydro  
     
    ! 3. Compute average and tangential velocicities
    Call Velocities
 !               
	!IF (Pref_Sal==1) THEN
	!	! Call Salinity           ! Salinity model
 !   ENDIF
 !       !WRITE(*,*) IT,'Compute waterTEMP'
 !   Call WaterTemp              ! Water temperature model
    ! WRITE(*,*) IT,'WRITE TIME SERIES'
      
    
End Subroutine HydrodynamicModel
    