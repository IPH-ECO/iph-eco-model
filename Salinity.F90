Subroutine Salinity(HydroParam,MeshParam,MeteoParam,LimnoParam,dt,dtday)
    
    ! This routine calculates the Salt Concentration in Water
    ! Called in routine 0-MAIN
    ! Call the following routine: UpWindCWC or LTS
    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars
    Use Meteorological

    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam
    type(LimnologyParam) :: LimnoParam
    Integer:: index,iElem,iLayer
    Real:: dt,dtday
        
    
    If (LimnoParam%iSal==1) Then 
        index = 2           ! Boundary Condition Index for Salinity (If no Boundary condition index = 0)
        LimnoParam%sDSalP = LimnoParam%sDSal      ! Save the previous Salt Concentration Field
        
        ! 1. Salt-Load Boundary Condition
        HydroParam%uLoadVarEst = LimnoParam%uDLoadSal
        
        ! 2. Source/Sink Term for Salinity
        HydroParam%dVarEst = 0.
    
        ! 3. Solver Transport equation
        If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDSal,LimnoParam%sDSalP,dt,dtday,HydroParam,MeshParam) 
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDSal,LimnoParam%sDSalP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
        ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDSal,LimnoParam%sDSalP,dt,dtday,HydroParam,MeshParam) 
            !Call LTS(iter,ContCell,sDSal,sDSalP)      ! Local Time Stepping Scheme
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDSal,LimnoParam%sDSalP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
        EndIf
    
    EndIf

	Return
End Subroutine Salinity
