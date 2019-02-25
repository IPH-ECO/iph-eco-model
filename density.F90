Subroutine Density(HydroParam,MeshParam,MeteoParam,LimnoParam,dt,dtday)
    
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
        LimnoParam%sDSalP = HydroParam%sDRhoW      ! Save the previous Salt Concentration Field
    
        Do iElem = 1,MeshParam%nElem
            ! 1. Salt-Load Boundary Condition
            !Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                !If (LimnoParam%IndexWQ(iElem,index)>0) then
                !    HydroParam%uLoadVarEst(iLayer,iElem) = LimnoParam%uDLoadSal(iLayer,LimnoParam%IndexWQ(iElem,index))
                !Else
                !    HydroParam%uLoadVarEst(iLayer,iElem) = LimnoParam%sDSal(iLayer,iElem)
                !EndIf 
            !EndDo
            ! 2. Source/Sink Term for Salinity
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                ! Obs.:  There is no Source/Sink Term for Salinity
                HydroParam%dVarEst(iLayer,iElem,1) = 0.
	        EndDo
        EndDo
    
        ! 3. Solver Transport equation
        If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),HydroParam%sDRhoW,LimnoParam%sDSalP,dt,dtday,HydroParam,MeshParam) 
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),HydroParam%sDRhoW,LimnoParam%sDSalP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
        ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),HydroParam%sDRhoW,LimnoParam%sDSalP,dt,dtday,HydroParam,MeshParam) 
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),HydroParam%sDRhoW,LimnoParam%sDSalP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
        EndIf
    
    EndIf

	Return
End Subroutine Density
