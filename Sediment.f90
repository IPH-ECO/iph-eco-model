SUBROUTINE Sediment(HydroParam,MeshParam,MeteoParam,LimnoParam,dt,dtday)
    
    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars
    Use Meteorological

    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    type(MeteorologicalParam) :: MeteoParam
    Real:: dt,dtday
    
    
    Call SedimentFluxes(HydroParam,MeshParam,LimnoParam,dt,dtday)
    
    !-------------------------------------------------------------------------------------!
    !                           Abiotic processes in the sediment					      !
    !-------------------------------------------------------------------------------------!
    Call AbioticSediment(HydroParam,MeshParam,LimnoParam,dt,dtday)
    !-------------------------------------------------------------------------------------!
    !                           Foodweb in the sediment                                   !
    !-------------------------------------------------------------------------------------!
    If (LimnoParam%InclBent== 1) Call Zoobenthos(HydroParam,MeshParam,LimnoParam,dt,dtday)
    
    If (LimnoParam%InclPhyt == 1) Call PhytoSediment(HydroParam,MeshParam,LimnoParam,dt,dtday)
    
    If (LimnoParam%InclBac == 1) Call BacSediment(HydroParam,MeshParam,LimnoParam,dt,dtday)
    
    Call UpdateSedVars(HydroParam,MeshParam,LimnoParam)
      

Return
End
