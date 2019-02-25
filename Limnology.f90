Subroutine Limnology(HydroParam,MeshParam,MeteoParam,LimnoParam,simParam)
	
    Use Hydrodynamic
    Use MeshVars
    Use Meteorological
	Use LimnologyVars
    Use SimulationModel
    
	Implicit None
	
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam
    type(LimnologyParam) :: LimnoParam
    type(SimulationParam) :: simParam
    
    Call Density(HydroParam,MeshParam,MeteoParam,LimnoParam,simParam%dt,simParam%dtday)
    
    !Call Salinity(HydroParam,MeshParam,MeteoParam,LimnoParam,simParam%dt,simParam%dtday)
    !
    Call WaterTemp(HydroParam,MeshParam,MeteoParam,LimnoParam,simParam%dt,simParam%dtday) 
    !
    !If (LimnoParam%iLimno == 1) Then 
    !
    !    If (LimnoParam%iSed==1) Call Sediment(HydroParam,MeshParam,MeteoParam,LimnoParam,simParam%dt,simParam%dtday)  
    !    
    !    !If (LimnoParam%InclFish==1) Call Fishes(HydroParam,MeshParam,LimnoParam,simParam%dt,simParam%dtday,simParam%Julday)
    !    
    !    If (LimnoParam%InclZoo==1) Call Zooplankton(HydroParam,MeshParam,LimnoParam,simParam%dt,simParam%dtday)
    !    
    !    If (LimnoParam%InclPhyt == 1) Call PhytoWater(HydroParam,MeshParam,MeteoParam,LimnoParam,simParam%dt,simParam%dtday)
    !    
    !    If (LimnoParam%InclMacr == 1) Call Macrophytes(HydroParam,MeshParam,MeteoParam,LimnoParam,simParam%dt,simParam%dtday)
    !    
    !    If (LimnoParam%InclBac == 1) Call Bacterioplankton(HydroParam,MeshParam,LimnoParam,simParam%dt,simParam%dtday)
    !    
    !    If (LimnoParam%InclMatInorg == 1) Call InorgMatter(HydroParam,MeshParam,LimnoParam,simParam%dt,simParam%dtday)
    !    
    !    If (LimnoParam%InclMatOrg == 1) Call OrgMatter(HydroParam,MeshParam,LimnoParam,simParam%dt,simParam%dtday)
    !    
    !    Call Nutrients(HydroParam,MeshParam,LimnoParam,simParam%dt,simParam%dtday)
    !    
    !    Call Gases(HydroParam,MeshParam,LimnoParam,simParam%dt,simParam%dtday)
    ! 
    !    !Call Ions
    ! 
	   ! Call UpdateLimnVars(HydroParam,MeshParam,LimnoParam,simParam)
    !
    !EndIf
	
	Return
End