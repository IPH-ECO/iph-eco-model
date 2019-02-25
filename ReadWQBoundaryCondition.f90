!> This subroutine reads the simulation parameters. 
Subroutine ReadWQBoundaryCondition(MeshParam,HydroParam,LimnoParam,hydroConfiguration,wqConfiguration,IniTime,FinalTime)
    
    Use domain_types
    Use MeshVars
    Use LimnologyVars
    Use Hydrodynamic
    
    
    Integer i,j,k,iElem,iEdge,iNode1,iNode2,Face,iLayer,iRange
    Integer:: idate(6)
    character(len=20):: str
    Integer:: NObjTemp,NObjSal,NObjIM,NObjDom,NObjPom,NObjPO4,NObjPAIM,NObjNH4,NObjNO3,NObjSiO2,NObjO2,NObjDBO,NObjDic
    Integer:: NObjTempRange,NObjSalRange,NObjIMRange,NObjDomRange,NObjPomRange,NObjPO4Range,NObjPAIMRange,NObjNH4Range,NObjNO3Range,NObjSiO2Range,NObjO2Range,NObjDBORange,NObjDicRange
    Integer(c_int):: numberOfBoundaryConditions
    type(HydrodynamicConfiguration) :: hydroConfiguration
    type(WaterQualityConfiguration) :: wqConfiguration
    type(BoundaryCondition), dimension(:), pointer :: wqboundaryConditions
    type(NonVerticallyIntegratedRange), pointer :: nonVerticallyIntegratedCell(:)
    !integer, pointer :: objectIds(:)
    type(HydrodynamicBoundaryConditionCell), pointer :: boundaryConditionCells(:)
    type(TimeSeries), pointer :: hydrotimeSeries(:)
    Integer:: IniTime !< Initial time of simulation (unix time)
    Integer:: FinalTime !< Final time of simulation (unix time)
    type(MeshGridParam) :: MeshParam
    type(LimnologyParam) :: LimnoParam
    type(HydrodynamicParam) :: HydroParam
    
    LimnoParam%NuDLoadTemp = 0 !< Number of boundary conditions of Water temperature
    NObjTemp = 0 !pegar o número de linhas referente ao maior arquivo
    NObjTempRange = 0 !pegar o número máximo de faixas na vertical
            
    LimnoParam%NuDLoadSal = 0 !< Number of boundary conditions of Water temperature
    NObjSal = 0 !pegar o número de linhas referente ao maior arquivo
    NObjSalRange = 0 !pegar o número máximo de faixas na vertical
    
    LimnoParam%NuDLoadPom = 0 !< Number of boundary conditions of Water temperature
    NObjPom = 0 !pegar o número de linhas referente ao maior arquivo
    NObjPomRange = 0 !pegar o número máximo de faixas na vertical
                
    LimnoParam%NuDLoadDom = 0 !< Number of boundary conditions of Water temperature
    NObjDom = 0 !pegar o número de linhas referente ao maior arquivo
    NObjDomRange = 0 !pegar o número máximo de faixas na vertical
            
    LimnoParam%NuDLoadIM = 0 !< Number of boundary conditions of Water temperature
    NObjIM = 0 !pegar o número de linhas referente ao maior arquivo
    NObjIMRange = 0 !pegar o número máximo de faixas na vertical
            
    LimnoParam%NuPLoadPO4 = 0 !< Number of boundary conditions of Water temperature
    NObjPO4 = 0 !pegar o número de linhas referente ao maior arquivo
    NObjPO4Range = 0 !pegar o número máximo de faixas na vertical
            
    LimnoParam%NuDLoadIM = 0 !< Number of boundary conditions of Water temperature
    NObjIM = 0 !pegar o número de linhas referente ao maior arquivo
    NObjIMRange = 0 !pegar o número máximo de faixas na vertical
            
    LimnoParam%NuNLoadNH4 = 0 !< Number of boundary conditions of Water temperature
    NObjNH4 = 0 !pegar o número de linhas referente ao maior arquivo
    NObjNH4Range = 0 !pegar o número máximo de faixas na vertical
            
    LimnoParam%NuNLoadNO3 = 0 !< Number of boundary conditions of Water temperature
    NObjNO3 = 0 !pegar o número de linhas referente ao maior arquivo
    NObjNO3Range = 0 !pegar o número máximo de faixas na vertical
            
    LimnoParam%NuSiLoadSiO2 = 0 !< Number of boundary conditions of Water temperature
    NObjSiO2 = 0 !pegar o número de linhas referente ao maior arquivo
    NObjSiO2Range = 0 !pegar o número máximo de faixas na vertical
            
    LimnoParam%NuO2LoadO2 = 0 !< Number of boundary conditions of Water temperature
    NObjO2 = 0 !pegar o número de linhas referente ao maior arquivo
    NObjO2Range = 0 !pegar o número máximo de faixas na vertical
            
    LimnoParam%NuDBOLoadDBO = 0 !< Number of boundary conditions of Water temperature
    NObjDBO = 0 !pegar o número de linhas referente ao maior arquivo
    NObjDBORange = 0 !pegar o número máximo de faixas na vertical
            
    LimnoParam%NuDicLoadDic = 0 !< Number of boundary conditions of Water temperature
    NObjDic = 0 !pegar o número de linhas referente ao maior arquivo
    NObjDicRange = 0 !pegar o número máximo de faixas na vertical
    
    call c_f_pointer(wqConfiguration%boundaryConditions, wqboundaryConditions, [wqConfiguration%numberOfBoundaryConditions])
    
    Do i = 1,wqConfiguration%numberOfBoundaryConditions
        ! 1. Water Temperature (Index == 1)
        ! 1.1 Mapping cells/layers with boundary condition of Water Temperature
        If (trim(wqboundaryConditions(i)%conditionType) == "Water Temperature") Then ! Mapping cells/layers with boundary condition of Water Temperature
            
            Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                If (wqboundaryConditions(i)%verticallyIntegrated) Then
                    LimnoParam%NuDLoadTemp = LimnoParam%NuDLoadTemp + 1
                    NObjTempRange = Max(NObjTempRange,1)
                    NObjTemp = Max(NObjTemp, wqboundaryConditions(i)%timeSeriesListSize)
                Else
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        LimnoParam%NuDLoadTemp = LimnoParam%NuDLoadTemp + 1
                        NObjTempRange = Max(NObjTempRange,wqboundaryConditions(i)%rangesSize)
                        If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                            NObjTemp = Max(NObjTemp,1)
                        Else
                            NObjTemp = Max(NObjTemp, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                        EndIf
                    EndDo
                EndIf
            EndDo
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Salinity") Then !Salinity
            
            Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                If (wqboundaryConditions(i)%verticallyIntegrated) Then
                    LimnoParam%NuDLoadSal = LimnoParam%NuDLoadSal + 1
                    NObjSalRange = Max(NObjSalRange,1)
                    NObjSal = Max(NObjSal, wqboundaryConditions(i)%timeSeriesListSize)
                Else
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        LimnoParam%NuDLoadSal = LimnoParam%NuDLoadSal + 1
                        NObjSalRange = Max(NObjSalRange,wqboundaryConditions(i)%rangesSize)
                        If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                            NObjSal = Max(NObjSal,1)
                        Else
                            NObjSal = Max(NObjSal, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                        EndIf
                    EndDo
                EndIf
            EndDo    
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Organic Matter") Then !Organic Matter
            
            If (LimnoParam%InclMatOrgSplit == 1) Then
            
                Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                    If (wqboundaryConditions(i)%verticallyIntegrated) Then
                        LimnoParam%NuDLoadPom = LimnoParam%NuDLoadPom + 1
                        NObjPomRange = Max(NObjPomRange,1)
                        NObjPom = Max(NObjPom, wqboundaryConditions(i)%timeSeriesListSize)
                    Else
                        call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                        Do iRange = 1,wqboundaryConditions(i)%rangesSize
                            LimnoParam%NuDLoadPom = LimnoParam%NuDLoadPom + 1
                            NObjPomRange = Max(NObjPomRange,wqboundaryConditions(i)%rangesSize)
                            If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                                NObjPom = Max(NObjPom,1)
                            Else
                                NObjPom = Max(NObjPom, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                            EndIf
                        EndDo
                    EndIf
                EndDo 
            
                Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                    If (wqboundaryConditions(i)%verticallyIntegrated) Then
                        LimnoParam%NuDLoadDom = LimnoParam%NuDLoadDom + 1
                        NObjDomRange = Max(NObjDomRange,1)
                        NObjDom = Max(NObjDom, wqboundaryConditions(i)%timeSeriesListSize)
                    Else
                        call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                        Do iRange = 1,wqboundaryConditions(i)%rangesSize
                            LimnoParam%NuDLoadDom = LimnoParam%NuDLoadDom + 1
                            NObjDomRange = Max(NObjDomRange,wqboundaryConditions(i)%rangesSize)
                            If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                                NObjDom = Max(NObjDom,1)
                            Else
                                NObjDom = Max(NObjDom, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                            EndIf
                        EndDo
                    EndIf
                EndDo  
                
            Else
            
                Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                    If (wqboundaryConditions(i)%verticallyIntegrated) Then
                        LimnoParam%NuDLoadPom = LimnoParam%NuDLoadPom + 1
                        NObjPomRange = Max(NObjPomRange,1)
                        NObjPom = Max(NObjPom, wqboundaryConditions(i)%timeSeriesListSize)
                    Else
                        call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                        Do iRange = 1,wqboundaryConditions(i)%rangesSize
                            LimnoParam%NuDLoadPom = LimnoParam%NuDLoadPom + 1
                            NObjPomRange = Max(NObjPomRange,wqboundaryConditions(i)%rangesSize)
                            If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                                NObjPom = Max(NObjPom,1)
                            Else
                                NObjPom = Max(NObjPom, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                            EndIf
                        EndDo
                    EndIf
                EndDo 
                
            EndIf
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Inorganic Matter") Then !Inorganic Matter
            
            Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                If (wqboundaryConditions(i)%verticallyIntegrated) Then
                    LimnoParam%NuDLoadIM = LimnoParam%NuDLoadIM + 1
                    NObjIMRange = Max(NObjIMRange,1)
                    NObjIM = Max(NObjIM, wqboundaryConditions(i)%timeSeriesListSize)
                Else
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        LimnoParam%NuDLoadIM = LimnoParam%NuDLoadIM + 1
                        NObjIMRange = Max(NObjIMRange,wqboundaryConditions(i)%rangesSize)
                        If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                            NObjIM = Max(NObjIM,1)
                        Else
                            NObjIM = Max(NObjIM, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                        EndIf
                    EndDo
                EndIf
            EndDo    
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Orthophosphate (PO4)") Then !Orthophosphate (PO4)
            
            Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                If (wqboundaryConditions(i)%verticallyIntegrated) Then
                    LimnoParam%NuPLoadPO4 = LimnoParam%NuPLoadPO4 + 1
                    NObjPO4Range = Max(NObjPO4Range,1)
                    NObjPO4 = Max(NObjPO4, wqboundaryConditions(i)%timeSeriesListSize)
                Else
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        LimnoParam%NuPLoadPO4 = LimnoParam%NuPLoadPO4 + 1
                        NObjPO4Range = Max(NObjPO4Range,wqboundaryConditions(i)%rangesSize)
                        If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                            NObjPO4 = Max(NObjPO4,1)
                        Else
                            NObjPO4 = Max(NObjPO4, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                        EndIf
                    EndDo
                EndIf
            EndDo 
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Particulate Absorbed Inorganic Phosphorus (PAIM)") Then !Particulate Absorbed Inorganic Phosphorus (PAIM)
            
            Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                If (wqboundaryConditions(i)%verticallyIntegrated) Then
                    LimnoParam%NuDLoadIM = LimnoParam%NuDLoadIM + 1
                    NObjIMRange = Max(NObjIMRange,1)
                    NObjIM = Max(NObjIM, wqboundaryConditions(i)%timeSeriesListSize)
                Else
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        LimnoParam%NuDLoadIM = LimnoParam%NuDLoadIM + 1
                        NObjIMRange = Max(NObjIMRange,wqboundaryConditions(i)%rangesSize)
                        If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                            NObjIM = Max(NObjIM,1)
                        Else
                            NObjIM = Max(NObjIM, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                        EndIf
                    EndDo
                EndIf
            EndDo    
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Ammonium (NH4)") Then !Ammonium (NH4)
            
            Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                If (wqboundaryConditions(i)%verticallyIntegrated) Then
                    LimnoParam%NuNLoadNH4 = LimnoParam%NuNLoadNH4 + 1
                    NObjNH4Range = Max(NObjNH4Range,1)
                    NObjNH4 = Max(NObjNH4, wqboundaryConditions(i)%timeSeriesListSize)
                Else
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        LimnoParam%NuNLoadNH4 = LimnoParam%NuNLoadNH4 + 1
                        NObjNH4Range = Max(NObjNH4Range,wqboundaryConditions(i)%rangesSize)
                        If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                            NObjNH4 = Max(NObjNH4,1)
                        Else
                            NObjNH4 = Max(NObjNH4, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                        EndIf
                    EndDo
                EndIf
            EndDo    
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Nitrate (NO3)") Then !Nitrate (NO3)
            
            Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                If (wqboundaryConditions(i)%verticallyIntegrated) Then
                    LimnoParam%NuNLoadNO3 = LimnoParam%NuNLoadNO3 + 1
                    NObjNO3Range = Max(NObjNO3Range,1)
                    NObjNO3 = Max(NObjNO3, wqboundaryConditions(i)%timeSeriesListSize)
                Else
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        LimnoParam%NuNLoadNO3 = LimnoParam%NuNLoadNO3 + 1
                        NObjNO3Range = Max(NObjNO3Range,wqboundaryConditions(i)%rangesSize)
                        If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                            NObjNO3 = Max(NObjNO3,1)
                        Else
                            NObjNO3 = Max(NObjNO3, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                        EndIf
                    EndDo
                EndIf
            EndDo    
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Silicate (SiO2)") Then !Silicate (SiO2)
            
            Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                If (wqboundaryConditions(i)%verticallyIntegrated) Then
                    LimnoParam%NuSiLoadSiO2 = LimnoParam%NuSiLoadSiO2 + 1
                    NObjSiO2Range = Max(NObjSiO2Range,1)
                    NObjSiO2 = Max(NObjSiO2, wqboundaryConditions(i)%timeSeriesListSize)
                Else
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        LimnoParam%NuSiLoadSiO2 = LimnoParam%NuSiLoadSiO2 + 1
                        NObjSiO2Range = Max(NObjSiO2Range,wqboundaryConditions(i)%rangesSize)
                        If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                            NObjSiO2 = Max(NObjSiO2,1)
                        Else
                            NObjSiO2 = Max(NObjSiO2, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                        EndIf
                    EndDo
                EndIf
            EndDo    
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Dissolved Oxygen (O2)") Then !Silicate (SiO2)    
            
            Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                If (wqboundaryConditions(i)%verticallyIntegrated) Then
                    LimnoParam%NuO2LoadO2 = LimnoParam%NuO2LoadO2 + 1
                    NObjO2Range = Max(NObjO2Range,1)
                    NObjO2 = Max(NObjO2, wqboundaryConditions(i)%timeSeriesListSize)
                Else
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        LimnoParam%NuO2LoadO2 = LimnoParam%NuO2LoadO2 + 1
                        NObjO2Range = Max(NObjO2Range,wqboundaryConditions(i)%rangesSize)
                        If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                            NObjO2 = Max(NObjO2,1)
                        Else
                            NObjO2 = Max(NObjO2, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                        EndIf
                    EndDo
                EndIf
            EndDo    
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Biological Oxygen Demand (BOD)") Then !Biological Oxygen Demand (BOD)
            
            Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                If (wqboundaryConditions(i)%verticallyIntegrated) Then
                    LimnoParam%NuDBOLoadDBO = LimnoParam%NuDBOLoadDBO + 1
                    NObjDBORange = Max(NObjDBORange,1)
                    NObjDBO = Max(NObjDBO, wqboundaryConditions(i)%timeSeriesListSize)
                Else
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        LimnoParam%NuDBOLoadDBO = LimnoParam%NuDBOLoadDBO + 1
                        NObjDBORange = Max(NObjDBORange,wqboundaryConditions(i)%rangesSize)
                        If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                            NObjDBO = Max(NObjDBO,1)
                        Else
                            NObjDBO = Max(NObjDBO, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                        EndIf
                    EndDo
                EndIf
            EndDo    
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Dissolved Inorganic Carbon") Then !Dissolved Inorganic Carbon
            
            Do j = 1,wqboundaryConditions(i)%cellsLength
                    
                If (wqboundaryConditions(i)%verticallyIntegrated) Then
                    LimnoParam%NuDicLoadDic = LimnoParam%NuDicLoadDic + 1
                    NObjDicRange = Max(NObjDicRange,1)
                    NObjDic = Max(NObjDic, wqboundaryConditions(i)%timeSeriesListSize)
                Else
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        LimnoParam%NuDicLoadDic = LimnoParam%NuDicLoadDic + 1
                        NObjDicRange = Max(NObjDicRange,wqboundaryConditions(i)%rangesSize)
                        If (nonVerticallyIntegratedCell(iRange)%function==1) Then
                            NObjDic = Max(NObjDic,1)
                        Else
                            NObjDic = Max(NObjDic, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                        EndIf
                    EndDo
                EndIf
            EndDo    
            
        EndIf
        
    EndDo

    ! 1.2 Allocating Boundary conditions
    Allocate(LimnoParam%uDLoadTemp(MeshParam%KMax,MeshParam%nElem))
    If (NObjTemp<2) NObjTemp = 2
    Allocate (LimnoParam%TempnTime(LimnoParam%NuDLoadTemp))    
    Allocate (LimnoParam%TempValue(LimnoParam%NuDLoadTemp,NObjTemp))
    Allocate (LimnoParam%TempTime(LimnoParam%NuDLoadTemp,NObjTemp))
    Allocate (LimnoParam%TempSmallm(LimnoParam%NuDLoadTemp))
    Allocate (LimnoParam%TempCapitalM(LimnoParam%NuDLoadTemp))
    Allocate (LimnoParam%NTempIndex(LimnoParam%NuDLoadTemp,2)) 
    
    Allocate(LimnoParam%uDLoadSal(MeshParam%KMax,MeshParam%nElem))
    If (NObjSal<2) NObjSal = 2
    Allocate (LimnoParam%SalnTime(LimnoParam%NuDLoadSal))    
    Allocate (LimnoParam%SalValue(LimnoParam%NuDLoadSal,NObjSal))
    Allocate (LimnoParam%SalTime(LimnoParam%NuDLoadSal,NObjSal))
    Allocate (LimnoParam%SalSmallm(LimnoParam%NuDLoadSal))
    Allocate (LimnoParam%SalCapitalM(LimnoParam%NuDLoadSal))
    Allocate (LimnoParam%NSalIndex(LimnoParam%NuDLoadSal,2)) 
    
    Allocate(LimnoParam%uDLoadPom(MeshParam%KMax,MeshParam%nElem))
    If (NObjPom<2) NObjPom = 2
    Allocate (LimnoParam%PomnTime(LimnoParam%NuDLoadPom))    
    Allocate (LimnoParam%PomValue(LimnoParam%NuDLoadPom,NObjPom))
    Allocate (LimnoParam%PomTime(LimnoParam%NuDLoadPom,NObjPom))
    Allocate (LimnoParam%PomSmallm(LimnoParam%NuDLoadPom))
    Allocate (LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom))
    Allocate (LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,2)) 
    
    Allocate(LimnoParam%uDLoadDom(MeshParam%KMax,MeshParam%nElem))
    If (NObjDom<2) NObjDom = 2
    Allocate (LimnoParam%DomnTime(LimnoParam%NuDLoadDom))    
    Allocate (LimnoParam%DomValue(LimnoParam%NuDLoadDom,NObjPom))
    Allocate (LimnoParam%DomTime(LimnoParam%NuDLoadDom,NObjPom))
    Allocate (LimnoParam%DomSmallm(LimnoParam%NuDLoadDom))
    Allocate (LimnoParam%DomCapitalM(LimnoParam%NuDLoadDom))
    Allocate (LimnoParam%NDomIndex(LimnoParam%NuDLoadDom,2)) 
    
    Allocate(LimnoParam%uDLoadIM(MeshParam%KMax,MeshParam%nElem))
    If (NObjIM<2) NObjIM = 2
    Allocate (LimnoParam%IMnTime(LimnoParam%NuDLoadIM))    
    Allocate (LimnoParam%IMValue(LimnoParam%NuDLoadIM,NObjIM))
    Allocate (LimnoParam%IMTime(LimnoParam%NuDLoadIM,NObjIM))
    Allocate (LimnoParam%IMSmallm(LimnoParam%NuDLoadIM))
    Allocate (LimnoParam%IMCapitalM(LimnoParam%NuDLoadIM))
    Allocate (LimnoParam%NIMIndex(LimnoParam%NuDLoadIM,2))   

    Allocate(LimnoParam%uPLoadPO4(MeshParam%KMax,MeshParam%nElem))
    If (NObjPO4<2) NObjPO4 = 2
    Allocate (LimnoParam%PO4nTime(LimnoParam%NuPLoadPO4))    
    Allocate (LimnoParam%PO4Value(LimnoParam%NuPLoadPO4,NObjPO4))
    Allocate (LimnoParam%PO4Time(LimnoParam%NuPLoadPO4,NObjPO4))
    Allocate (LimnoParam%PO4Smallm(LimnoParam%NuPLoadPO4))
    Allocate (LimnoParam%PO4CapitalM(LimnoParam%NuPLoadPO4))
    Allocate (LimnoParam%NPO4Index(LimnoParam%NuPLoadPO4,2)) 
    
    Allocate(LimnoParam%uPLoadPAIM(MeshParam%KMax,MeshParam%nElem))
    If (NObjPAIM<2) NObjPAIM = 2
    Allocate (LimnoParam%PAIMnTime(LimnoParam%NuPLoadPAIM))    
    Allocate (LimnoParam%PAIMValue(LimnoParam%NuPLoadPAIM,NObjPAIM))
    Allocate (LimnoParam%PAIMTime(LimnoParam%NuPLoadPAIM,NObjPAIM))
    Allocate (LimnoParam%PAIMSmallm(LimnoParam%NuPLoadPAIM))
    Allocate (LimnoParam%PAIMCapitalM(LimnoParam%NuPLoadPAIM))
    Allocate (LimnoParam%NPAIMIndex(LimnoParam%NuPLoadPAIM,2))   
    
    Allocate(LimnoParam%uNLoadNH4(MeshParam%KMax,MeshParam%nElem))
    If (NObjNH4<2) NObjNH4 = 2
    Allocate (LimnoParam%NH4nTime(LimnoParam%NuNLoadNH4))    
    Allocate (LimnoParam%NH4Value(LimnoParam%NuNLoadNH4,NObjNH4))
    Allocate (LimnoParam%NH4Time(LimnoParam%NuNLoadNH4,NObjNH4))
    Allocate (LimnoParam%NH4Smallm(LimnoParam%NuNLoadNH4))
    Allocate (LimnoParam%NH4CapitalM(LimnoParam%NuNLoadNH4))
    Allocate (LimnoParam%NNH4Index(LimnoParam%NuNLoadNH4,2))   
    
    Allocate(LimnoParam%uNLoadNO3(MeshParam%KMax,MeshParam%nElem))
    If (NObjNO3<2) NObjNO3 = 2
    Allocate (LimnoParam%NO3nTime(LimnoParam%NuNLoadNO3))    
    Allocate (LimnoParam%NO3Value(LimnoParam%NuNLoadNO3,NObjNO3))
    Allocate (LimnoParam%NO3Time(LimnoParam%NuNLoadNO3,NObjNO3))
    Allocate (LimnoParam%NO3Smallm(LimnoParam%NuNLoadNO3))
    Allocate (LimnoParam%NO3CapitalM(LimnoParam%NuNLoadNO3))
    Allocate (LimnoParam%NNO3Index(LimnoParam%NuNLoadNO3,2))                 
            
    Allocate(LimnoParam%uSiLoadSiO2(MeshParam%KMax,MeshParam%nElem))
    If (NObjSiO2<2) NObjSiO2 = 2
    Allocate (LimnoParam%SiO2nTime(LimnoParam%NuSiLoadSiO2))    
    Allocate (LimnoParam%SiO2Value(LimnoParam%NuSiLoadSiO2,NObjSiO2))
    Allocate (LimnoParam%SiO2Time(LimnoParam%NuSiLoadSiO2,NObjSiO2))
    Allocate (LimnoParam%SiO2Smallm(LimnoParam%NuSiLoadSiO2))
    Allocate (LimnoParam%SiO2CapitalM(LimnoParam%NuSiLoadSiO2))
    Allocate (LimnoParam%NSiO2Index(LimnoParam%NuSiLoadSiO2,2))    
                
    Allocate(LimnoParam%uO2LoadO2(MeshParam%KMax,MeshParam%nElem))
    If (NObjO2<2) NObjO2 = 2
    Allocate (LimnoParam%O2nTime(LimnoParam%NuO2LoadO2))    
    Allocate (LimnoParam%O2Value(LimnoParam%NuO2LoadO2,NObjO2))
    Allocate (LimnoParam%O2Time(LimnoParam%NuO2LoadO2,NObjO2))
    Allocate (LimnoParam%O2Smallm(LimnoParam%NuO2LoadO2))
    Allocate (LimnoParam%O2CapitalM(LimnoParam%NuO2LoadO2))
    Allocate (LimnoParam%NO2Index(LimnoParam%NuO2LoadO2,2))       
    
    Allocate(LimnoParam%uDLoadDBO(MeshParam%KMax,MeshParam%nElem))
    If (NObjDBO<2) NObjDBO = 2
    Allocate (LimnoParam%DBOnTime(LimnoParam%NuDBOLoadDBO))    
    Allocate (LimnoParam%DBOValue(LimnoParam%NuDBOLoadDBO,NObjDBO))
    Allocate (LimnoParam%DBOTime(LimnoParam%NuDBOLoadDBO,NObjDBO))
    Allocate (LimnoParam%DBOSmallm(LimnoParam%NuDBOLoadDBO))
    Allocate (LimnoParam%DBOCapitalM(LimnoParam%NuDBOLoadDBO))
    Allocate (LimnoParam%NDBOIndex(LimnoParam%NuDBOLoadDBO,2)) 
    
    Allocate(LimnoParam%uDLoadDic(MeshParam%KMax,MeshParam%nElem))
    If (NObjDic<2) NObjDic = 2
    Allocate (LimnoParam%DicnTime(LimnoParam%NuDicLoadDic))    
    Allocate (LimnoParam%DicValue(LimnoParam%NuDicLoadDic,NObjDic))
    Allocate (LimnoParam%DicTime(LimnoParam%NuDicLoadDic,NObjDic))
    Allocate (LimnoParam%DicSmallm(LimnoParam%NuDicLoadDic))
    Allocate (LimnoParam%DicCapitalM(LimnoParam%NuDicLoadDic))
    Allocate (LimnoParam%NDicIndex(LimnoParam%NuDicLoadDic,2))             
            
                
    LimnoParam%NuDLoadTemp = 0
    LimnoParam%NuDLoadSal = 0
    LimnoParam%NuDLoadPom = 0
    LimnoParam%NuDLoadDom = 0
    LimnoParam%NuDLoadIM = 0
    LimnoParam%NuPLoadPO4 = 0
    LimnoParam%NuPLoadPAIM = 0
    LimnoParam%NuNLoadNH4 = 0
    LimnoParam%NuNLoadNO3 = 0
    LimnoParam%NuSiLoadSiO2 = 0
    LimnoParam%NuO2LoadO2 = 0
    LimnoParam%NuDBOLoadDBO = 0
    LimnoParam%NuDicLoadDic = 0
    
    Do i = 1,wqConfiguration%numberOfBoundaryConditions
            !
        call c_f_pointer(wqboundaryConditions(i)%cells, boundaryConditionCells, [wqboundaryConditions(i)%cellsLength])
            
        If (trim(wqboundaryConditions(i)%conditionType) == "Water Temperature") Then !Water Temperature
            
            If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                iRange = 1
                If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuDLoadTemp = LimnoParam%NuDLoadTemp + 1
                        LimnoParam%NTempIndex(LimnoParam%NuDLoadTemp,1) = 1
                        LimnoParam%NTempIndex(LimnoParam%NuDLoadTemp,2) = boundaryConditionCells(k)%cellId + 1
                        LimnoParam%TempnTime(LimnoParam%NuDLoadTemp) = 2
                        LimnoParam%TempSmallm(LimnoParam%NuDLoadTemp) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%TempCapitalM(LimnoParam%NuDLoadTemp) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%TempTime(LimnoParam%NuDLoadTemp,1) = IniTime
                        LimnoParam%TempTime(LimnoParam%NuDLoadTemp,2) = FinalTime
                        LimnoParam%TempValue(LimnoParam%NuDLoadTemp,:) = wqboundaryConditions(i)%constantValue
                    EndDo
                        
                ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuDLoadTemp = LimnoParam%NuDLoadTemp + 1
                        LimnoParam%NTempIndex(LimnoParam%NuDLoadTemp,1) = 1
                        LimnoParam%NTempIndex(LimnoParam%NuDLoadTemp,2) = boundaryConditionCells(k)%cellId + 1
                        If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                            LimnoParam%TempnTime(LimnoParam%NuDLoadTemp) = 2
                            LimnoParam%TempSmallm(LimnoParam%NuDLoadTemp) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%TempCapitalM(LimnoParam%NuDLoadTemp) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%TempTime(LimnoParam%NuDLoadTemp,1) = IniTime
                            LimnoParam%TempTime(LimnoParam%NuDLoadTemp,2) = FinalTime
                            LimnoParam%TempValue(LimnoParam%NuDLoadTemp,:) = wqboundaryConditions(i)%constantValue
                        Else
                            LimnoParam%TempnTime(LimnoParam%NuDLoadTemp) = wqboundaryConditions(i)%timeSeriesListSize
                            LimnoParam%TempSmallm(LimnoParam%NuDLoadTemp) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%TempCapitalM(LimnoParam%NuDLoadTemp) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                            Do j = 1, LimnoParam%TempnTime(LimnoParam%NuDLoadTemp)
                                LimnoParam%TempTime(LimnoParam%NuDLoadTemp,j) = hydrotimeSeries(j)%timeStamp
                                LimnoParam%TempValue(LimnoParam%NuDLoadTemp,j) = hydrotimeSeries(j)%value1 
                            EndDo                                   
                        EndIf
                    EndDo
                EndIf
                    
                    
            Else ! Non Vertically Integrated
                call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                Do iRange = 1,wqboundaryConditions(i)%rangesSize
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDLoadTemp = LimnoParam%NuDLoadTemp + 1
                            LimnoParam%NTempIndex(LimnoParam%NuDLoadTemp,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NTempIndex(LimnoParam%NuDLoadTemp,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%TempnTime(LimnoParam%NuDLoadTemp) = 2
                            iElem = boundaryConditionCells(k)%cellId + 1
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                    LimnoParam%TempSmallm(LimnoParam%NuDLoadTemp) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%TempSmallm(LimnoParam%NuDLoadTemp) = HydroParam%ElSmallm(iElem)
                                EndIf
                            EndDo
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                    LimnoParam%TempCapitalM(LimnoParam%NuDLoadTemp) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%TempCapitalM(LimnoParam%NuDLoadTemp) = HydroParam%ElCapitalM(iElem)
                                EndIf
                            EndDo                                   
                            LimnoParam%TempTime(LimnoParam%NuDLoadTemp,1) = IniTime
                            LimnoParam%TempTime(LimnoParam%NuDLoadTemp,2) = FinalTime
                            LimnoParam%TempValue(LimnoParam%NuDLoadTemp,:) = nonVerticallyIntegratedCell(iRange)%Value
                        EndDo
                            
                    ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDLoadTemp = LimnoParam%NuDLoadTemp + 1
                            LimnoParam%NTempIndex(LimnoParam%NuDLoadTemp,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NTempIndex(LimnoParam%NuDLoadTemp,2) = boundaryConditionCells(k)%cellId + 1
                            If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                LimnoParam%TempnTime(LimnoParam%NuDLoadTemp) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%TempSmallm(LimnoParam%NuDLoadTemp) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%TempSmallm(LimnoParam%NuDLoadTemp) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%TempCapitalM(LimnoParam%NuDLoadTemp) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%TempCapitalM(LimnoParam%NuDLoadTemp) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%TempTime(LimnoParam%NuDLoadTemp,1) = IniTime
                                LimnoParam%TempTime(LimnoParam%NuDLoadTemp,2) = FinalTime
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                LimnoParam%TempValue(LimnoParam%NuDLoadTemp,:) = hydrotimeSeries(1)%value1     
                            Else
                                LimnoParam%TempnTime(LimnoParam%NuDLoadTemp) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%TempSmallm(LimnoParam%NuDLoadTemp) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%TempSmallm(LimnoParam%NuDLoadTemp) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%TempCapitalM(LimnoParam%NuDLoadTemp) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%TempCapitalM(LimnoParam%NuDLoadTemp) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                Do j = 1, LimnoParam%TempnTime(LimnoParam%NuDLoadTemp)
                                    LimnoParam%TempTime(LimnoParam%NuDLoadTemp,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%TempValue(LimnoParam%NuDLoadTemp,j) = hydrotimeSeries(j)%value1 
                                EndDo
                            EndIf
                        EndDo
                    EndIf
                EndDo
            EndIf
                
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Salinity") Then !Salinity
            
            If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                iRange = 1
                If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuDLoadSal = LimnoParam%NuDLoadSal + 1
                        LimnoParam%NSalIndex(LimnoParam%NuDLoadSal,1) = 1
                        LimnoParam%NSalIndex(LimnoParam%NuDLoadSal,2) = boundaryConditionCells(k)%cellId + 1
                        LimnoParam%SalnTime(LimnoParam%NuDLoadSal) = 2
                        LimnoParam%SalSmallm(LimnoParam%NuDLoadSal) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%SalCapitalM(LimnoParam%NuDLoadSal) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%SalTime(LimnoParam%NuDLoadSal,1) = IniTime
                        LimnoParam%SalTime(LimnoParam%NuDLoadSal,2) = FinalTime
                        LimnoParam%SalValue(LimnoParam%NuDLoadSal,:) = wqboundaryConditions(i)%constantValue
                    EndDo
                        
                ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuDLoadSal = LimnoParam%NuDLoadSal + 1
                        LimnoParam%NSalIndex(LimnoParam%NuDLoadSal,1) = 1
                        LimnoParam%NSalIndex(LimnoParam%NuDLoadSal,2) = boundaryConditionCells(k)%cellId + 1
                        If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                            LimnoParam%SalnTime(LimnoParam%NuDLoadSal) = 2
                            LimnoParam%SalSmallm(LimnoParam%NuDLoadSal) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%SalCapitalM(LimnoParam%NuDLoadSal) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%SalTime(LimnoParam%NuDLoadSal,1) = IniTime
                            LimnoParam%SalTime(LimnoParam%NuDLoadSal,2) = FinalTime
                            LimnoParam%SalValue(LimnoParam%NuDLoadSal,:) = wqboundaryConditions(i)%constantValue
                        Else
                            LimnoParam%SalnTime(LimnoParam%NuDLoadSal) = wqboundaryConditions(i)%timeSeriesListSize
                            LimnoParam%SalSmallm(LimnoParam%NuDLoadSal) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%SalCapitalM(LimnoParam%NuDLoadSal) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                            Do j = 1, LimnoParam%SalnTime(LimnoParam%NuDLoadSal)
                                LimnoParam%SalTime(LimnoParam%NuDLoadSal,j) = hydrotimeSeries(j)%timeStamp
                                LimnoParam%SalValue(LimnoParam%NuDLoadSal,j) = hydrotimeSeries(j)%value1 
                            EndDo                                   
                        EndIf
                    EndDo
                EndIf
                    
                    
            Else ! Non Vertically Integrated
                call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                Do iRange = 1,wqboundaryConditions(i)%rangesSize
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDLoadSal = LimnoParam%NuDLoadSal + 1
                            LimnoParam%NSalIndex(LimnoParam%NuDLoadSal,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NSalIndex(LimnoParam%NuDLoadSal,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%SalnTime(LimnoParam%NuDLoadSal) = 2
                            iElem = boundaryConditionCells(k)%cellId + 1
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                    LimnoParam%SalSmallm(LimnoParam%NuDLoadSal) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%SalSmallm(LimnoParam%NuDLoadSal) = HydroParam%ElSmallm(iElem)
                                EndIf
                            EndDo
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                    LimnoParam%SalCapitalM(LimnoParam%NuDLoadSal) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%SalCapitalM(LimnoParam%NuDLoadSal) = HydroParam%ElCapitalM(iElem)
                                EndIf
                            EndDo                                   
                            LimnoParam%SalTime(LimnoParam%NuDLoadSal,1) = IniTime
                            LimnoParam%SalTime(LimnoParam%NuDLoadSal,2) = FinalTime
                            LimnoParam%SalValue(LimnoParam%NuDLoadSal,:) = nonVerticallyIntegratedCell(iRange)%Value
                        EndDo
                            
                    ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDLoadSal = LimnoParam%NuDLoadSal + 1
                            LimnoParam%NSalIndex(LimnoParam%NuDLoadSal,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NSalIndex(LimnoParam%NuDLoadSal,2) = boundaryConditionCells(k)%cellId + 1
                            If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                LimnoParam%SalnTime(LimnoParam%NuDLoadSal) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%SalSmallm(LimnoParam%NuDLoadSal) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%SalSmallm(LimnoParam%NuDLoadSal) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%SalCapitalM(LimnoParam%NuDLoadSal) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%SalCapitalM(LimnoParam%NuDLoadSal) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%SalTime(LimnoParam%NuDLoadSal,1) = IniTime
                                LimnoParam%SalTime(LimnoParam%NuDLoadSal,2) = FinalTime
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                LimnoParam%SalValue(LimnoParam%NuDLoadSal,:) = hydrotimeSeries(1)%value1     
                            Else
                                LimnoParam%SalnTime(LimnoParam%NuDLoadSal) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%SalSmallm(LimnoParam%NuDLoadSal) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%SalSmallm(LimnoParam%NuDLoadSal) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%SalCapitalM(LimnoParam%NuDLoadSal) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%SalCapitalM(LimnoParam%NuDLoadSal) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                Do j = 1, LimnoParam%SalnTime(LimnoParam%NuDLoadSal)
                                    LimnoParam%SalTime(LimnoParam%NuDLoadSal,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%SalValue(LimnoParam%NuDLoadSal,j) = hydrotimeSeries(j)%value1 
                                EndDo
                            EndIf
                        EndDo
                    EndIf
                EndDo
            EndIf
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Organic Matter") Then !Organic Matter
            
            If (LimnoParam%InclMatOrgSplit == 1) Then
                
                If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                    iRange = 1
                    If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDLoadPom = LimnoParam%NuDLoadPom + 1
                            LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,1) = 1
                            LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%PomnTime(LimnoParam%NuDLoadPom) = 2
                            LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%PomTime(LimnoParam%NuDLoadPom,1) = IniTime
                            LimnoParam%PomTime(LimnoParam%NuDLoadPom,2) = FinalTime
                            LimnoParam%PomValue(LimnoParam%NuDLoadPom,:) = wqboundaryConditions(i)%constantValue
                        EndDo
                        
                    ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDLoadPom = LimnoParam%NuDLoadPom + 1
                            LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,1) = 1
                            LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,2) = boundaryConditionCells(k)%cellId + 1
                            If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                                LimnoParam%PomnTime(LimnoParam%NuDLoadPom) = 2
                                LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                                LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                                LimnoParam%PomTime(LimnoParam%NuDLoadPom,1) = IniTime
                                LimnoParam%PomTime(LimnoParam%NuDLoadPom,2) = FinalTime
                                LimnoParam%PomValue(LimnoParam%NuDLoadPom,:) = wqboundaryConditions(i)%constantValue
                            Else
                                LimnoParam%PomnTime(LimnoParam%NuDLoadPom) = wqboundaryConditions(i)%timeSeriesListSize
                                LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                                LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                                call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                                Do j = 1, LimnoParam%PomnTime(LimnoParam%NuDLoadPom)
                                    LimnoParam%PomTime(LimnoParam%NuDLoadPom,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%PomValue(LimnoParam%NuDLoadPom,j) = hydrotimeSeries(j)%value1 
                                EndDo                                   
                            EndIf
                        EndDo
                    EndIf
                    
                    
                Else ! Non Vertically Integrated
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                            Do k = 1,wqboundaryConditions(i)%cellsLength
                                LimnoParam%NuDLoadPom = LimnoParam%NuDLoadPom + 1
                                LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,1) = wqboundaryConditions(i)%rangesSize
                                LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,2) = boundaryConditionCells(k)%cellId + 1
                                LimnoParam%PomnTime(LimnoParam%NuDLoadPom) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%PomTime(LimnoParam%NuDLoadPom,1) = IniTime
                                LimnoParam%PomTime(LimnoParam%NuDLoadPom,2) = FinalTime
                                LimnoParam%PomValue(LimnoParam%NuDLoadPom,:) = nonVerticallyIntegratedCell(iRange)%Value
                            EndDo
                            
                        ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                            Do k = 1,wqboundaryConditions(i)%cellsLength
                                LimnoParam%NuDLoadPom = LimnoParam%NuDLoadPom + 1
                                LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,1) = wqboundaryConditions(i)%rangesSize
                                LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,2) = boundaryConditionCells(k)%cellId + 1
                                If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                    LimnoParam%PomnTime(LimnoParam%NuDLoadPom) = 2
                                    iElem = boundaryConditionCells(k)%cellId + 1
                                    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                        If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                            LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = HydroParam%ElSmallm(iElem)
                                        EndIf
                                    EndDo
                                    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                        If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                            LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = HydroParam%ElCapitalM(iElem)
                                        EndIf
                                    EndDo                                   
                                    LimnoParam%PomTime(LimnoParam%NuDLoadPom,1) = IniTime
                                    LimnoParam%PomTime(LimnoParam%NuDLoadPom,2) = FinalTime
                                    call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                    LimnoParam%PomValue(LimnoParam%NuDLoadPom,:) = hydrotimeSeries(1)%value1     
                                Else
                                    LimnoParam%PomnTime(LimnoParam%NuDLoadPom) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                    iElem = boundaryConditionCells(k)%cellId + 1
                                    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                        If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                            LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = HydroParam%ElSmallm(iElem)
                                        EndIf
                                    EndDo
                                    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                        If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                            LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = HydroParam%ElCapitalM(iElem)
                                        EndIf
                                    EndDo                                   
                                    call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                    Do j = 1, LimnoParam%PomnTime(LimnoParam%NuDLoadPom)
                                        LimnoParam%PomTime(LimnoParam%NuDLoadPom,j) = hydrotimeSeries(j)%timeStamp
                                        LimnoParam%PomValue(LimnoParam%NuDLoadPom,j) = hydrotimeSeries(j)%value1 
                                    EndDo
                                EndIf
                            EndDo
                        EndIf
                    EndDo
                EndIf     
                
                If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                    iRange = 1
                    If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDLoadDom = LimnoParam%NuDLoadDom + 1
                            LimnoParam%NDomIndex(LimnoParam%NuDLoadDom,1) = 1
                            LimnoParam%NDomIndex(LimnoParam%NuDLoadDom,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%DomnTime(LimnoParam%NuDLoadDom) = 2
                            LimnoParam%DomSmallm(LimnoParam%NuDLoadDom) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%DomCapitalM(LimnoParam%NuDLoadDom) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%DomTime(LimnoParam%NuDLoadDom,1) = IniTime
                            LimnoParam%DomTime(LimnoParam%NuDLoadDom,2) = FinalTime
                            LimnoParam%DomValue(LimnoParam%NuDLoadDom,:) = wqboundaryConditions(i)%constantValue
                        EndDo
                        
                    ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDLoadDom = LimnoParam%NuDLoadDom + 1
                            LimnoParam%NDomIndex(LimnoParam%NuDLoadDom,1) = 1
                            LimnoParam%NDomIndex(LimnoParam%NuDLoadDom,2) = boundaryConditionCells(k)%cellId + 1
                            If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                                LimnoParam%DomnTime(LimnoParam%NuDLoadDom) = 2
                                LimnoParam%DomSmallm(LimnoParam%NuDLoadDom) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                                LimnoParam%DomCapitalM(LimnoParam%NuDLoadDom) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                                LimnoParam%DomTime(LimnoParam%NuDLoadDom,1) = IniTime
                                LimnoParam%DomTime(LimnoParam%NuDLoadDom,2) = FinalTime
                                LimnoParam%DomValue(LimnoParam%NuDLoadDom,:) = wqboundaryConditions(i)%constantValue
                            Else
                                LimnoParam%DomnTime(LimnoParam%NuDLoadDom) = wqboundaryConditions(i)%timeSeriesListSize
                                LimnoParam%DomSmallm(LimnoParam%NuDLoadDom) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                                LimnoParam%DomCapitalM(LimnoParam%NuDLoadDom) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                                call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                                Do j = 1, LimnoParam%DomnTime(LimnoParam%NuDLoadDom)
                                    LimnoParam%DomTime(LimnoParam%NuDLoadDom,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%DomValue(LimnoParam%NuDLoadDom,j) = hydrotimeSeries(j)%value1 
                                EndDo                                   
                            EndIf
                        EndDo
                    EndIf
                    
                    
                Else ! Non Vertically Integrated
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                            Do k = 1,wqboundaryConditions(i)%cellsLength
                                LimnoParam%NuDLoadDom = LimnoParam%NuDLoadDom + 1
                                LimnoParam%NDomIndex(LimnoParam%NuDLoadDom,1) = wqboundaryConditions(i)%rangesSize
                                LimnoParam%NDomIndex(LimnoParam%NuDLoadDom,2) = boundaryConditionCells(k)%cellId + 1
                                LimnoParam%DomnTime(LimnoParam%NuDLoadDom) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%DomSmallm(LimnoParam%NuDLoadDom) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%DomSmallm(LimnoParam%NuDLoadDom) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%DomCapitalM(LimnoParam%NuDLoadDom) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%DomCapitalM(LimnoParam%NuDLoadDom) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%DomTime(LimnoParam%NuDLoadDom,1) = IniTime
                                LimnoParam%DomTime(LimnoParam%NuDLoadDom,2) = FinalTime
                                LimnoParam%DomValue(LimnoParam%NuDLoadDom,:) = nonVerticallyIntegratedCell(iRange)%Value
                            EndDo
                            
                        ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                            Do k = 1,wqboundaryConditions(i)%cellsLength
                                LimnoParam%NuDLoadDom = LimnoParam%NuDLoadDom + 1
                                LimnoParam%NDomIndex(LimnoParam%NuDLoadDom,1) = wqboundaryConditions(i)%rangesSize
                                LimnoParam%NDomIndex(LimnoParam%NuDLoadDom,2) = boundaryConditionCells(k)%cellId + 1
                                If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                    LimnoParam%DomnTime(LimnoParam%NuDLoadDom) = 2
                                    iElem = boundaryConditionCells(k)%cellId + 1
                                    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                        If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                            LimnoParam%DomSmallm(LimnoParam%NuDLoadDom) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            LimnoParam%DomSmallm(LimnoParam%NuDLoadDom) = HydroParam%ElSmallm(iElem)
                                        EndIf
                                    EndDo
                                    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                        If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                            LimnoParam%DomCapitalM(LimnoParam%NuDLoadDom) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            LimnoParam%DomCapitalM(LimnoParam%NuDLoadDom) = HydroParam%ElCapitalM(iElem)
                                        EndIf
                                    EndDo                                   
                                    LimnoParam%DomTime(LimnoParam%NuDLoadDom,1) = IniTime
                                    LimnoParam%DomTime(LimnoParam%NuDLoadDom,2) = FinalTime
                                    call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                    LimnoParam%DomValue(LimnoParam%NuDLoadDom,:) = hydrotimeSeries(1)%value1     
                                Else
                                    LimnoParam%DomnTime(LimnoParam%NuDLoadDom) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                    iElem = boundaryConditionCells(k)%cellId + 1
                                    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                        If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                            LimnoParam%DomSmallm(LimnoParam%NuDLoadDom) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            LimnoParam%DomSmallm(LimnoParam%NuDLoadDom) = HydroParam%ElSmallm(iElem)
                                        EndIf
                                    EndDo
                                    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                        If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                            LimnoParam%DomCapitalM(LimnoParam%NuDLoadDom) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            LimnoParam%DomCapitalM(LimnoParam%NuDLoadDom) = HydroParam%ElCapitalM(iElem)
                                        EndIf
                                    EndDo                                   
                                    call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                    Do j = 1, LimnoParam%DomnTime(LimnoParam%NuDLoadDom)
                                        LimnoParam%DomTime(LimnoParam%NuDLoadDom,j) = hydrotimeSeries(j)%timeStamp
                                        LimnoParam%DomValue(LimnoParam%NuDLoadDom,j) = hydrotimeSeries(j)%value1 
                                    EndDo
                                EndIf
                            EndDo
                        EndIf
                    EndDo
                EndIf                     
                
            Else

                If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                    iRange = 1
                    If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDLoadPom = LimnoParam%NuDLoadPom + 1
                            LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,1) = 1
                            LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%PomnTime(LimnoParam%NuDLoadPom) = 2
                            LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%PomTime(LimnoParam%NuDLoadPom,1) = IniTime
                            LimnoParam%PomTime(LimnoParam%NuDLoadPom,2) = FinalTime
                            LimnoParam%PomValue(LimnoParam%NuDLoadPom,:) = wqboundaryConditions(i)%constantValue
                        EndDo
                        
                    ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDLoadPom = LimnoParam%NuDLoadPom + 1
                            LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,1) = 1
                            LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,2) = boundaryConditionCells(k)%cellId + 1
                            If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                                LimnoParam%PomnTime(LimnoParam%NuDLoadPom) = 2
                                LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                                LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                                LimnoParam%PomTime(LimnoParam%NuDLoadPom,1) = IniTime
                                LimnoParam%PomTime(LimnoParam%NuDLoadPom,2) = FinalTime
                                LimnoParam%PomValue(LimnoParam%NuDLoadPom,:) = wqboundaryConditions(i)%constantValue
                            Else
                                LimnoParam%PomnTime(LimnoParam%NuDLoadPom) = wqboundaryConditions(i)%timeSeriesListSize
                                LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                                LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                                call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                                Do j = 1, LimnoParam%PomnTime(LimnoParam%NuDLoadPom)
                                    LimnoParam%PomTime(LimnoParam%NuDLoadPom,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%PomValue(LimnoParam%NuDLoadPom,j) = hydrotimeSeries(j)%value1 
                                EndDo                                   
                            EndIf
                        EndDo
                    EndIf
                    
                    
                Else ! Non Vertically Integrated
                    call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                    Do iRange = 1,wqboundaryConditions(i)%rangesSize
                        If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                            Do k = 1,wqboundaryConditions(i)%cellsLength
                                LimnoParam%NuDLoadPom = LimnoParam%NuDLoadPom + 1
                                LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,1) = wqboundaryConditions(i)%rangesSize
                                LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,2) = boundaryConditionCells(k)%cellId + 1
                                LimnoParam%PomnTime(LimnoParam%NuDLoadPom) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%PomTime(LimnoParam%NuDLoadPom,1) = IniTime
                                LimnoParam%PomTime(LimnoParam%NuDLoadPom,2) = FinalTime
                                LimnoParam%PomValue(LimnoParam%NuDLoadPom,:) = nonVerticallyIntegratedCell(iRange)%Value
                            EndDo
                            
                        ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                            Do k = 1,wqboundaryConditions(i)%cellsLength
                                LimnoParam%NuDLoadPom = LimnoParam%NuDLoadPom + 1
                                LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,1) = wqboundaryConditions(i)%rangesSize
                                LimnoParam%NPomIndex(LimnoParam%NuDLoadPom,2) = boundaryConditionCells(k)%cellId + 1
                                If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                    LimnoParam%PomnTime(LimnoParam%NuDLoadPom) = 2
                                    iElem = boundaryConditionCells(k)%cellId + 1
                                    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                        If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                            LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = HydroParam%ElSmallm(iElem)
                                        EndIf
                                    EndDo
                                    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                        If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                            LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = HydroParam%ElCapitalM(iElem)
                                        EndIf
                                    EndDo                                   
                                    LimnoParam%PomTime(LimnoParam%NuDLoadPom,1) = IniTime
                                    LimnoParam%PomTime(LimnoParam%NuDLoadPom,2) = FinalTime
                                    call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                    LimnoParam%PomValue(LimnoParam%NuDLoadPom,:) = hydrotimeSeries(1)%value1     
                                Else
                                    LimnoParam%PomnTime(LimnoParam%NuDLoadPom) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                    iElem = boundaryConditionCells(k)%cellId + 1
                                    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                        If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                            LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            LimnoParam%PomSmallm(LimnoParam%NuDLoadPom) = HydroParam%ElSmallm(iElem)
                                        EndIf
                                    EndDo
                                    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                        If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                            LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            LimnoParam%PomCapitalM(LimnoParam%NuDLoadPom) = HydroParam%ElCapitalM(iElem)
                                        EndIf
                                    EndDo                                   
                                    call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                    Do j = 1, LimnoParam%PomnTime(LimnoParam%NuDLoadPom)
                                        LimnoParam%PomTime(LimnoParam%NuDLoadPom,j) = hydrotimeSeries(j)%timeStamp
                                        LimnoParam%PomValue(LimnoParam%NuDLoadPom,j) = hydrotimeSeries(j)%value1 
                                    EndDo
                                EndIf
                            EndDo
                        EndIf
                    EndDo
                EndIf                     
                
            EndIf
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Inorganic Matter") Then !Inorganic Matter
            
            If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                iRange = 1
                If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuDLoadIM = LimnoParam%NuDLoadIM + 1
                        LimnoParam%NIMIndex(LimnoParam%NuDLoadIM,1) = 1
                        LimnoParam%NIMIndex(LimnoParam%NuDLoadIM,2) = boundaryConditionCells(k)%cellId + 1
                        LimnoParam%IMnTime(LimnoParam%NuDLoadIM) = 2
                        LimnoParam%IMSmallm(LimnoParam%NuDLoadIM) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%IMCapitalM(LimnoParam%NuDLoadIM) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%IMTime(LimnoParam%NuDLoadIM,1) = IniTime
                        LimnoParam%IMTime(LimnoParam%NuDLoadIM,2) = FinalTime
                        LimnoParam%IMValue(LimnoParam%NuDLoadIM,:) = wqboundaryConditions(i)%constantValue
                    EndDo
                        
                ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuDLoadIM = LimnoParam%NuDLoadIM + 1
                        LimnoParam%NIMIndex(LimnoParam%NuDLoadIM,1) = 1
                        LimnoParam%NIMIndex(LimnoParam%NuDLoadIM,2) = boundaryConditionCells(k)%cellId + 1
                        If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                            LimnoParam%IMnTime(LimnoParam%NuDLoadIM) = 2
                            LimnoParam%IMSmallm(LimnoParam%NuDLoadIM) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%IMCapitalM(LimnoParam%NuDLoadIM) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%IMTime(LimnoParam%NuDLoadIM,1) = IniTime
                            LimnoParam%IMTime(LimnoParam%NuDLoadIM,2) = FinalTime
                            LimnoParam%IMValue(LimnoParam%NuDLoadIM,:) = wqboundaryConditions(i)%constantValue
                        Else
                            LimnoParam%IMnTime(LimnoParam%NuDLoadIM) = wqboundaryConditions(i)%timeSeriesListSize
                            LimnoParam%IMSmallm(LimnoParam%NuDLoadIM) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%IMCapitalM(LimnoParam%NuDLoadIM) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                            Do j = 1, LimnoParam%IMnTime(LimnoParam%NuDLoadIM)
                                LimnoParam%IMTime(LimnoParam%NuDLoadIM,j) = hydrotimeSeries(j)%timeStamp
                                LimnoParam%IMValue(LimnoParam%NuDLoadIM,j) = hydrotimeSeries(j)%value1 
                            EndDo                                   
                        EndIf
                    EndDo
                EndIf
                    
                    
            Else ! Non Vertically Integrated
                call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                Do iRange = 1,wqboundaryConditions(i)%rangesSize
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDLoadIM = LimnoParam%NuDLoadIM + 1
                            LimnoParam%NIMIndex(LimnoParam%NuDLoadIM,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NIMIndex(LimnoParam%NuDLoadIM,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%IMnTime(LimnoParam%NuDLoadIM) = 2
                            iElem = boundaryConditionCells(k)%cellId + 1
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                    LimnoParam%IMSmallm(LimnoParam%NuDLoadIM) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%IMSmallm(LimnoParam%NuDLoadIM) = HydroParam%ElSmallm(iElem)
                                EndIf
                            EndDo
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                    LimnoParam%IMCapitalM(LimnoParam%NuDLoadIM) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%IMCapitalM(LimnoParam%NuDLoadIM) = HydroParam%ElCapitalM(iElem)
                                EndIf
                            EndDo                                   
                            LimnoParam%IMTime(LimnoParam%NuDLoadIM,1) = IniTime
                            LimnoParam%IMTime(LimnoParam%NuDLoadIM,2) = FinalTime
                            LimnoParam%IMValue(LimnoParam%NuDLoadIM,:) = nonVerticallyIntegratedCell(iRange)%Value
                        EndDo
                            
                    ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDLoadIM = LimnoParam%NuDLoadIM + 1
                            LimnoParam%NIMIndex(LimnoParam%NuDLoadIM,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NIMIndex(LimnoParam%NuDLoadIM,2) = boundaryConditionCells(k)%cellId + 1
                            If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                LimnoParam%IMnTime(LimnoParam%NuDLoadIM) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%IMSmallm(LimnoParam%NuDLoadIM) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%IMSmallm(LimnoParam%NuDLoadIM) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%IMCapitalM(LimnoParam%NuDLoadIM) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%IMCapitalM(LimnoParam%NuDLoadIM) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%IMTime(LimnoParam%NuDLoadIM,1) = IniTime
                                LimnoParam%IMTime(LimnoParam%NuDLoadIM,2) = FinalTime
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                LimnoParam%IMValue(LimnoParam%NuDLoadIM,:) = hydrotimeSeries(1)%value1     
                            Else
                                LimnoParam%IMnTime(LimnoParam%NuDLoadIM) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%IMSmallm(LimnoParam%NuDLoadIM) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%IMSmallm(LimnoParam%NuDLoadIM) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%IMCapitalM(LimnoParam%NuDLoadIM) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%IMCapitalM(LimnoParam%NuDLoadIM) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                Do j = 1, LimnoParam%IMnTime(LimnoParam%NuDLoadIM)
                                    LimnoParam%IMTime(LimnoParam%NuDLoadIM,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%IMValue(LimnoParam%NuDLoadIM,j) = hydrotimeSeries(j)%value1 
                                EndDo
                            EndIf
                        EndDo
                    EndIf
                EndDo
            EndIf            
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Orthophosphate (PO4)") Then !Orthophosphate (PO4)
            
            If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                iRange = 1
                If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuPLoadPO4 = LimnoParam%NuPLoadPO4 + 1
                        LimnoParam%NPO4Index(LimnoParam%NuPLoadPO4,1) = 1
                        LimnoParam%NPO4Index(LimnoParam%NuPLoadPO4,2) = boundaryConditionCells(k)%cellId + 1
                        LimnoParam%PO4nTime(LimnoParam%NuPLoadPO4) = 2
                        LimnoParam%PO4Smallm(LimnoParam%NuPLoadPO4) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%PO4CapitalM(LimnoParam%NuPLoadPO4) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%PO4Time(LimnoParam%NuPLoadPO4,1) = IniTime
                        LimnoParam%PO4Time(LimnoParam%NuPLoadPO4,2) = FinalTime
                        LimnoParam%PO4Value(LimnoParam%NuPLoadPO4,:) = wqboundaryConditions(i)%constantValue
                    EndDo
                        
                ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuPLoadPO4 = LimnoParam%NuPLoadPO4 + 1
                        LimnoParam%NPO4Index(LimnoParam%NuPLoadPO4,1) = 1
                        LimnoParam%NPO4Index(LimnoParam%NuPLoadPO4,2) = boundaryConditionCells(k)%cellId + 1
                        If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                            LimnoParam%PO4nTime(LimnoParam%NuPLoadPO4) = 2
                            LimnoParam%PO4Smallm(LimnoParam%NuPLoadPO4) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%PO4CapitalM(LimnoParam%NuPLoadPO4) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%PO4Time(LimnoParam%NuPLoadPO4,1) = IniTime
                            LimnoParam%PO4Time(LimnoParam%NuPLoadPO4,2) = FinalTime
                            LimnoParam%PO4Value(LimnoParam%NuPLoadPO4,:) = wqboundaryConditions(i)%constantValue
                        Else
                            LimnoParam%PO4nTime(LimnoParam%NuPLoadPO4) = wqboundaryConditions(i)%timeSeriesListSize
                            LimnoParam%PO4Smallm(LimnoParam%NuPLoadPO4) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%PO4CapitalM(LimnoParam%NuPLoadPO4) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                            Do j = 1, LimnoParam%PO4nTime(LimnoParam%NuPLoadPO4)
                                LimnoParam%PO4Time(LimnoParam%NuPLoadPO4,j) = hydrotimeSeries(j)%timeStamp
                                LimnoParam%PO4Value(LimnoParam%NuPLoadPO4,j) = hydrotimeSeries(j)%value1 
                            EndDo                                   
                        EndIf
                    EndDo
                EndIf
                    
                    
            Else ! Non Vertically Integrated
                call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                Do iRange = 1,wqboundaryConditions(i)%rangesSize
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuPLoadPO4 = LimnoParam%NuPLoadPO4 + 1
                            LimnoParam%NPO4Index(LimnoParam%NuPLoadPO4,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NPO4Index(LimnoParam%NuPLoadPO4,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%PO4nTime(LimnoParam%NuPLoadPO4) = 2
                            iElem = boundaryConditionCells(k)%cellId + 1
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                    LimnoParam%PO4Smallm(LimnoParam%NuPLoadPO4) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%PO4Smallm(LimnoParam%NuPLoadPO4) = HydroParam%ElSmallm(iElem)
                                EndIf
                            EndDo
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                    LimnoParam%PO4CapitalM(LimnoParam%NuPLoadPO4) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%PO4CapitalM(LimnoParam%NuPLoadPO4) = HydroParam%ElCapitalM(iElem)
                                EndIf
                            EndDo                                   
                            LimnoParam%PO4Time(LimnoParam%NuPLoadPO4,1) = IniTime
                            LimnoParam%PO4Time(LimnoParam%NuPLoadPO4,2) = FinalTime
                            LimnoParam%PO4Value(LimnoParam%NuPLoadPO4,:) = nonVerticallyIntegratedCell(iRange)%Value
                        EndDo
                            
                    ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuPLoadPO4 = LimnoParam%NuPLoadPO4 + 1
                            LimnoParam%NPO4Index(LimnoParam%NuPLoadPO4,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NPO4Index(LimnoParam%NuPLoadPO4,2) = boundaryConditionCells(k)%cellId + 1
                            If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                LimnoParam%PO4nTime(LimnoParam%NuPLoadPO4) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%PO4Smallm(LimnoParam%NuPLoadPO4) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%PO4Smallm(LimnoParam%NuPLoadPO4) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%PO4CapitalM(LimnoParam%NuPLoadPO4) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%PO4CapitalM(LimnoParam%NuPLoadPO4) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%PO4Time(LimnoParam%NuPLoadPO4,1) = IniTime
                                LimnoParam%PO4Time(LimnoParam%NuPLoadPO4,2) = FinalTime
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                LimnoParam%PO4Value(LimnoParam%NuPLoadPO4,:) = hydrotimeSeries(1)%value1     
                            Else
                                LimnoParam%PO4nTime(LimnoParam%NuPLoadPO4) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%PO4Smallm(LimnoParam%NuPLoadPO4) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%PO4Smallm(LimnoParam%NuPLoadPO4) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%PO4CapitalM(LimnoParam%NuPLoadPO4) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%PO4CapitalM(LimnoParam%NuPLoadPO4) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                Do j = 1, LimnoParam%PO4nTime(LimnoParam%NuPLoadPO4)
                                    LimnoParam%PO4Time(LimnoParam%NuPLoadPO4,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%PO4Value(LimnoParam%NuPLoadPO4,j) = hydrotimeSeries(j)%value1 
                                EndDo
                            EndIf
                        EndDo
                    EndIf
                EndDo
            EndIf            
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Particulate Absorbed Inorganic Phosphorus (PAIM)") Then !Particulate Absorbed Inorganic Phosphorus (PAIM)
            
            If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                iRange = 1
                If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuPLoadPAIM = LimnoParam%NuPLoadPAIM + 1
                        LimnoParam%NPAIMIndex(LimnoParam%NuPLoadPAIM,1) = 1
                        LimnoParam%NPAIMIndex(LimnoParam%NuPLoadPAIM,2) = boundaryConditionCells(k)%cellId + 1
                        LimnoParam%PAIMnTime(LimnoParam%NuPLoadPAIM) = 2
                        LimnoParam%PAIMSmallm(LimnoParam%NuPLoadPAIM) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%PAIMCapitalM(LimnoParam%NuPLoadPAIM) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%PAIMTime(LimnoParam%NuPLoadPAIM,1) = IniTime
                        LimnoParam%PAIMTime(LimnoParam%NuPLoadPAIM,2) = FinalTime
                        LimnoParam%PAIMValue(LimnoParam%NuPLoadPAIM,:) = wqboundaryConditions(i)%constantValue
                    EndDo
                        
                ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuPLoadPAIM = LimnoParam%NuPLoadPAIM + 1
                        LimnoParam%NPAIMIndex(LimnoParam%NuPLoadPAIM,1) = 1
                        LimnoParam%NPAIMIndex(LimnoParam%NuPLoadPAIM,2) = boundaryConditionCells(k)%cellId + 1
                        If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                            LimnoParam%PAIMnTime(LimnoParam%NuPLoadPAIM) = 2
                            LimnoParam%PAIMSmallm(LimnoParam%NuPLoadPAIM) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%PAIMCapitalM(LimnoParam%NuPLoadPAIM) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%PAIMTime(LimnoParam%NuPLoadPAIM,1) = IniTime
                            LimnoParam%PAIMTime(LimnoParam%NuPLoadPAIM,2) = FinalTime
                            LimnoParam%PAIMValue(LimnoParam%NuPLoadPAIM,:) = wqboundaryConditions(i)%constantValue
                        Else
                            LimnoParam%PAIMnTime(LimnoParam%NuPLoadPAIM) = wqboundaryConditions(i)%timeSeriesListSize
                            LimnoParam%PAIMSmallm(LimnoParam%NuPLoadPAIM) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%PAIMCapitalM(LimnoParam%NuPLoadPAIM) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                            Do j = 1, LimnoParam%PAIMnTime(LimnoParam%NuPLoadPAIM)
                                LimnoParam%PAIMTime(LimnoParam%NuPLoadPAIM,j) = hydrotimeSeries(j)%timeStamp
                                LimnoParam%PAIMValue(LimnoParam%NuPLoadPAIM,j) = hydrotimeSeries(j)%value1 
                            EndDo                                   
                        EndIf
                    EndDo
                EndIf
                    
                    
            Else ! Non Vertically Integrated
                call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                Do iRange = 1,wqboundaryConditions(i)%rangesSize
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuPLoadPAIM = LimnoParam%NuPLoadPAIM + 1
                            LimnoParam%NPAIMIndex(LimnoParam%NuPLoadPAIM,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NPAIMIndex(LimnoParam%NuPLoadPAIM,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%PAIMnTime(LimnoParam%NuPLoadPAIM) = 2
                            iElem = boundaryConditionCells(k)%cellId + 1
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                    LimnoParam%PAIMSmallm(LimnoParam%NuPLoadPAIM) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%PAIMSmallm(LimnoParam%NuPLoadPAIM) = HydroParam%ElSmallm(iElem)
                                EndIf
                            EndDo
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                    LimnoParam%PAIMCapitalM(LimnoParam%NuPLoadPAIM) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%PAIMCapitalM(LimnoParam%NuPLoadPAIM) = HydroParam%ElCapitalM(iElem)
                                EndIf
                            EndDo                                   
                            LimnoParam%PAIMTime(LimnoParam%NuPLoadPAIM,1) = IniTime
                            LimnoParam%PAIMTime(LimnoParam%NuPLoadPAIM,2) = FinalTime
                            LimnoParam%PAIMValue(LimnoParam%NuPLoadPAIM,:) = nonVerticallyIntegratedCell(iRange)%Value
                        EndDo
                            
                    ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuPLoadPAIM = LimnoParam%NuPLoadPAIM + 1
                            LimnoParam%NPAIMIndex(LimnoParam%NuPLoadPAIM,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NPAIMIndex(LimnoParam%NuPLoadPAIM,2) = boundaryConditionCells(k)%cellId + 1
                            If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                LimnoParam%PAIMnTime(LimnoParam%NuPLoadPAIM) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%PAIMSmallm(LimnoParam%NuPLoadPAIM) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%PAIMSmallm(LimnoParam%NuPLoadPAIM) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%PAIMCapitalM(LimnoParam%NuPLoadPAIM) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%PAIMCapitalM(LimnoParam%NuPLoadPAIM) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%PAIMTime(LimnoParam%NuPLoadPAIM,1) = IniTime
                                LimnoParam%PAIMTime(LimnoParam%NuPLoadPAIM,2) = FinalTime
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                LimnoParam%PAIMValue(LimnoParam%NuPLoadPAIM,:) = hydrotimeSeries(1)%value1     
                            Else
                                LimnoParam%PAIMnTime(LimnoParam%NuPLoadPAIM) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%PAIMSmallm(LimnoParam%NuPLoadPAIM) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%PAIMSmallm(LimnoParam%NuPLoadPAIM) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%PAIMCapitalM(LimnoParam%NuPLoadPAIM) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%PAIMCapitalM(LimnoParam%NuPLoadPAIM) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                Do j = 1, LimnoParam%PAIMnTime(LimnoParam%NuPLoadPAIM)
                                    LimnoParam%PAIMTime(LimnoParam%NuPLoadPAIM,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%PAIMValue(LimnoParam%NuPLoadPAIM,j) = hydrotimeSeries(j)%value1 
                                EndDo
                            EndIf
                        EndDo
                    EndIf
                EndDo
            EndIf                   
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Ammonium (NH4)") Then !Ammonium (NH4)
            
            If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                iRange = 1
                If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuNLoadNH4 = LimnoParam%NuNLoadNH4 + 1
                        LimnoParam%NNH4Index(LimnoParam%NuNLoadNH4,1) = 1
                        LimnoParam%NNH4Index(LimnoParam%NuNLoadNH4,2) = boundaryConditionCells(k)%cellId + 1
                        LimnoParam%NH4nTime(LimnoParam%NuNLoadNH4) = 2
                        LimnoParam%NH4Smallm(LimnoParam%NuNLoadNH4) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%NH4CapitalM(LimnoParam%NuNLoadNH4) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%NH4Time(LimnoParam%NuNLoadNH4,1) = IniTime
                        LimnoParam%NH4Time(LimnoParam%NuNLoadNH4,2) = FinalTime
                        LimnoParam%NH4Value(LimnoParam%NuNLoadNH4,:) = wqboundaryConditions(i)%constantValue
                    EndDo
                        
                ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuNLoadNH4 = LimnoParam%NuNLoadNH4 + 1
                        LimnoParam%NNH4Index(LimnoParam%NuNLoadNH4,1) = 1
                        LimnoParam%NNH4Index(LimnoParam%NuNLoadNH4,2) = boundaryConditionCells(k)%cellId + 1
                        If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                            LimnoParam%NH4nTime(LimnoParam%NuNLoadNH4) = 2
                            LimnoParam%NH4Smallm(LimnoParam%NuNLoadNH4) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%NH4CapitalM(LimnoParam%NuNLoadNH4) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%NH4Time(LimnoParam%NuNLoadNH4,1) = IniTime
                            LimnoParam%NH4Time(LimnoParam%NuNLoadNH4,2) = FinalTime
                            LimnoParam%NH4Value(LimnoParam%NuNLoadNH4,:) = wqboundaryConditions(i)%constantValue
                        Else
                            LimnoParam%NH4nTime(LimnoParam%NuNLoadNH4) = wqboundaryConditions(i)%timeSeriesListSize
                            LimnoParam%NH4Smallm(LimnoParam%NuNLoadNH4) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%NH4CapitalM(LimnoParam%NuNLoadNH4) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                            Do j = 1, LimnoParam%NH4nTime(LimnoParam%NuNLoadNH4)
                                LimnoParam%NH4Time(LimnoParam%NuNLoadNH4,j) = hydrotimeSeries(j)%timeStamp
                                LimnoParam%NH4Value(LimnoParam%NuNLoadNH4,j) = hydrotimeSeries(j)%value1 
                            EndDo                                   
                        EndIf
                    EndDo
                EndIf
                    
                    
            Else ! Non Vertically Integrated
                call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                Do iRange = 1,wqboundaryConditions(i)%rangesSize
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuNLoadNH4 = LimnoParam%NuNLoadNH4 + 1
                            LimnoParam%NNH4Index(LimnoParam%NuNLoadNH4,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NNH4Index(LimnoParam%NuNLoadNH4,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%NH4nTime(LimnoParam%NuNLoadNH4) = 2
                            iElem = boundaryConditionCells(k)%cellId + 1
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                    LimnoParam%NH4Smallm(LimnoParam%NuNLoadNH4) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%NH4Smallm(LimnoParam%NuNLoadNH4) = HydroParam%ElSmallm(iElem)
                                EndIf
                            EndDo
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                    LimnoParam%NH4CapitalM(LimnoParam%NuNLoadNH4) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%NH4CapitalM(LimnoParam%NuNLoadNH4) = HydroParam%ElCapitalM(iElem)
                                EndIf
                            EndDo                                   
                            LimnoParam%NH4Time(LimnoParam%NuNLoadNH4,1) = IniTime
                            LimnoParam%NH4Time(LimnoParam%NuNLoadNH4,2) = FinalTime
                            LimnoParam%NH4Value(LimnoParam%NuNLoadNH4,:) = nonVerticallyIntegratedCell(iRange)%Value
                        EndDo
                            
                    ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuNLoadNH4 = LimnoParam%NuNLoadNH4 + 1
                            LimnoParam%NNH4Index(LimnoParam%NuNLoadNH4,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NNH4Index(LimnoParam%NuNLoadNH4,2) = boundaryConditionCells(k)%cellId + 1
                            If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                LimnoParam%NH4nTime(LimnoParam%NuNLoadNH4) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%NH4Smallm(LimnoParam%NuNLoadNH4) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%NH4Smallm(LimnoParam%NuNLoadNH4) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%NH4CapitalM(LimnoParam%NuNLoadNH4) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%NH4CapitalM(LimnoParam%NuNLoadNH4) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%NH4Time(LimnoParam%NuNLoadNH4,1) = IniTime
                                LimnoParam%NH4Time(LimnoParam%NuNLoadNH4,2) = FinalTime
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                LimnoParam%NH4Value(LimnoParam%NuNLoadNH4,:) = hydrotimeSeries(1)%value1     
                            Else
                                LimnoParam%NH4nTime(LimnoParam%NuNLoadNH4) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%NH4Smallm(LimnoParam%NuNLoadNH4) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%NH4Smallm(LimnoParam%NuNLoadNH4) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%NH4CapitalM(LimnoParam%NuNLoadNH4) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%NH4CapitalM(LimnoParam%NuNLoadNH4) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                Do j = 1, LimnoParam%NH4nTime(LimnoParam%NuNLoadNH4)
                                    LimnoParam%NH4Time(LimnoParam%NuNLoadNH4,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%NH4Value(LimnoParam%NuNLoadNH4,j) = hydrotimeSeries(j)%value1 
                                EndDo
                            EndIf
                        EndDo
                    EndIf
                EndDo
            EndIf                     
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Nitrate (NO3)") Then !Nitrate (NO3)
            
            If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                iRange = 1
                If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuNLoadNO3 = LimnoParam%NuNLoadNO3 + 1
                        LimnoParam%NNO3Index(LimnoParam%NuNLoadNO3,1) = 1
                        LimnoParam%NNO3Index(LimnoParam%NuNLoadNO3,2) = boundaryConditionCells(k)%cellId + 1
                        LimnoParam%NO3nTime(LimnoParam%NuNLoadNO3) = 2
                        LimnoParam%NO3Smallm(LimnoParam%NuNLoadNO3) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%NO3CapitalM(LimnoParam%NuNLoadNO3) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%NO3Time(LimnoParam%NuNLoadNO3,1) = IniTime
                        LimnoParam%NO3Time(LimnoParam%NuNLoadNO3,2) = FinalTime
                        LimnoParam%NO3Value(LimnoParam%NuNLoadNO3,:) = wqboundaryConditions(i)%constantValue
                    EndDo
                        
                ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuNLoadNO3 = LimnoParam%NuNLoadNO3 + 1
                        LimnoParam%NNO3Index(LimnoParam%NuNLoadNO3,1) = 1
                        LimnoParam%NNO3Index(LimnoParam%NuNLoadNO3,2) = boundaryConditionCells(k)%cellId + 1
                        If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                            LimnoParam%NO3nTime(LimnoParam%NuNLoadNO3) = 2
                            LimnoParam%NO3Smallm(LimnoParam%NuNLoadNO3) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%NO3CapitalM(LimnoParam%NuNLoadNO3) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%NO3Time(LimnoParam%NuNLoadNO3,1) = IniTime
                            LimnoParam%NO3Time(LimnoParam%NuNLoadNO3,2) = FinalTime
                            LimnoParam%NO3Value(LimnoParam%NuNLoadNO3,:) = wqboundaryConditions(i)%constantValue
                        Else
                            LimnoParam%NO3nTime(LimnoParam%NuNLoadNO3) = wqboundaryConditions(i)%timeSeriesListSize
                            LimnoParam%NO3Smallm(LimnoParam%NuNLoadNO3) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%NO3CapitalM(LimnoParam%NuNLoadNO3) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                            Do j = 1, LimnoParam%NO3nTime(LimnoParam%NuNLoadNO3)
                                LimnoParam%NO3Time(LimnoParam%NuNLoadNO3,j) = hydrotimeSeries(j)%timeStamp
                                LimnoParam%NO3Value(LimnoParam%NuNLoadNO3,j) = hydrotimeSeries(j)%value1 
                            EndDo                                   
                        EndIf
                    EndDo
                EndIf
                    
                    
            Else ! Non Vertically Integrated
                call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                Do iRange = 1,wqboundaryConditions(i)%rangesSize
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuNLoadNO3 = LimnoParam%NuNLoadNO3 + 1
                            LimnoParam%NNO3Index(LimnoParam%NuNLoadNO3,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NNO3Index(LimnoParam%NuNLoadNO3,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%NO3nTime(LimnoParam%NuNLoadNO3) = 2
                            iElem = boundaryConditionCells(k)%cellId + 1
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                    LimnoParam%NO3Smallm(LimnoParam%NuNLoadNO3) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%NO3Smallm(LimnoParam%NuNLoadNO3) = HydroParam%ElSmallm(iElem)
                                EndIf
                            EndDo
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                    LimnoParam%NO3CapitalM(LimnoParam%NuNLoadNO3) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%NO3CapitalM(LimnoParam%NuNLoadNO3) = HydroParam%ElCapitalM(iElem)
                                EndIf
                            EndDo                                   
                            LimnoParam%NO3Time(LimnoParam%NuNLoadNO3,1) = IniTime
                            LimnoParam%NO3Time(LimnoParam%NuNLoadNO3,2) = FinalTime
                            LimnoParam%NO3Value(LimnoParam%NuNLoadNO3,:) = nonVerticallyIntegratedCell(iRange)%Value
                        EndDo
                            
                    ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuNLoadNO3 = LimnoParam%NuNLoadNO3 + 1
                            LimnoParam%NNO3Index(LimnoParam%NuNLoadNO3,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NNO3Index(LimnoParam%NuNLoadNO3,2) = boundaryConditionCells(k)%cellId + 1
                            If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                LimnoParam%NO3nTime(LimnoParam%NuNLoadNO3) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%NO3Smallm(LimnoParam%NuNLoadNO3) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%NO3Smallm(LimnoParam%NuNLoadNO3) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%NO3CapitalM(LimnoParam%NuNLoadNO3) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%NO3CapitalM(LimnoParam%NuNLoadNO3) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%NO3Time(LimnoParam%NuNLoadNO3,1) = IniTime
                                LimnoParam%NO3Time(LimnoParam%NuNLoadNO3,2) = FinalTime
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                LimnoParam%NO3Value(LimnoParam%NuNLoadNO3,:) = hydrotimeSeries(1)%value1     
                            Else
                                LimnoParam%NO3nTime(LimnoParam%NuNLoadNO3) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%NO3Smallm(LimnoParam%NuNLoadNO3) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%NO3Smallm(LimnoParam%NuNLoadNO3) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%NO3CapitalM(LimnoParam%NuNLoadNO3) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%NO3CapitalM(LimnoParam%NuNLoadNO3) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                Do j = 1, LimnoParam%NO3nTime(LimnoParam%NuNLoadNO3)
                                    LimnoParam%NO3Time(LimnoParam%NuNLoadNO3,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%NO3Value(LimnoParam%NuNLoadNO3,j) = hydrotimeSeries(j)%value1 
                                EndDo
                            EndIf
                        EndDo
                    EndIf
                EndDo
            EndIf                       
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Silicate (SiO2)") Then !Silicate (SiO2)
            
            If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                iRange = 1
                If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuSiLoadSiO2 = LimnoParam%NuSiLoadSiO2 + 1
                        LimnoParam%NSiO2Index(LimnoParam%NuSiLoadSiO2,1) = 1
                        LimnoParam%NSiO2Index(LimnoParam%NuSiLoadSiO2,2) = boundaryConditionCells(k)%cellId + 1
                        LimnoParam%SiO2nTime(LimnoParam%NuSiLoadSiO2) = 2
                        LimnoParam%SiO2Smallm(LimnoParam%NuSiLoadSiO2) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%SiO2CapitalM(LimnoParam%NuSiLoadSiO2) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%SiO2Time(LimnoParam%NuSiLoadSiO2,1) = IniTime
                        LimnoParam%SiO2Time(LimnoParam%NuSiLoadSiO2,2) = FinalTime
                        LimnoParam%SiO2Value(LimnoParam%NuSiLoadSiO2,:) = wqboundaryConditions(i)%constantValue
                    EndDo
                        
                ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuSiLoadSiO2 = LimnoParam%NuSiLoadSiO2 + 1
                        LimnoParam%NSiO2Index(LimnoParam%NuSiLoadSiO2,1) = 1
                        LimnoParam%NSiO2Index(LimnoParam%NuSiLoadSiO2,2) = boundaryConditionCells(k)%cellId + 1
                        If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                            LimnoParam%SiO2nTime(LimnoParam%NuSiLoadSiO2) = 2
                            LimnoParam%SiO2Smallm(LimnoParam%NuSiLoadSiO2) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%SiO2CapitalM(LimnoParam%NuSiLoadSiO2) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%SiO2Time(LimnoParam%NuSiLoadSiO2,1) = IniTime
                            LimnoParam%SiO2Time(LimnoParam%NuSiLoadSiO2,2) = FinalTime
                            LimnoParam%SiO2Value(LimnoParam%NuSiLoadSiO2,:) = wqboundaryConditions(i)%constantValue
                        Else
                            LimnoParam%SiO2nTime(LimnoParam%NuSiLoadSiO2) = wqboundaryConditions(i)%timeSeriesListSize
                            LimnoParam%SiO2Smallm(LimnoParam%NuSiLoadSiO2) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%SiO2CapitalM(LimnoParam%NuSiLoadSiO2) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                            Do j = 1, LimnoParam%SiO2nTime(LimnoParam%NuSiLoadSiO2)
                                LimnoParam%SiO2Time(LimnoParam%NuSiLoadSiO2,j) = hydrotimeSeries(j)%timeStamp
                                LimnoParam%SiO2Value(LimnoParam%NuSiLoadSiO2,j) = hydrotimeSeries(j)%value1 
                            EndDo                                   
                        EndIf
                    EndDo
                EndIf
                    
                    
            Else ! Non Vertically Integrated
                call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                Do iRange = 1,wqboundaryConditions(i)%rangesSize
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuSiLoadSiO2 = LimnoParam%NuSiLoadSiO2 + 1
                            LimnoParam%NSiO2Index(LimnoParam%NuSiLoadSiO2,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NSiO2Index(LimnoParam%NuSiLoadSiO2,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%SiO2nTime(LimnoParam%NuSiLoadSiO2) = 2
                            iElem = boundaryConditionCells(k)%cellId + 1
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                    LimnoParam%SiO2Smallm(LimnoParam%NuSiLoadSiO2) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%SiO2Smallm(LimnoParam%NuSiLoadSiO2) = HydroParam%ElSmallm(iElem)
                                EndIf
                            EndDo
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                    LimnoParam%SiO2CapitalM(LimnoParam%NuSiLoadSiO2) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%SiO2CapitalM(LimnoParam%NuSiLoadSiO2) = HydroParam%ElCapitalM(iElem)
                                EndIf
                            EndDo                                   
                            LimnoParam%SiO2Time(LimnoParam%NuSiLoadSiO2,1) = IniTime
                            LimnoParam%SiO2Time(LimnoParam%NuSiLoadSiO2,2) = FinalTime
                            LimnoParam%SiO2Value(LimnoParam%NuSiLoadSiO2,:) = nonVerticallyIntegratedCell(iRange)%Value
                        EndDo
                            
                    ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuSiLoadSiO2 = LimnoParam%NuSiLoadSiO2 + 1
                            LimnoParam%NSiO2Index(LimnoParam%NuSiLoadSiO2,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NSiO2Index(LimnoParam%NuSiLoadSiO2,2) = boundaryConditionCells(k)%cellId + 1
                            If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                LimnoParam%SiO2nTime(LimnoParam%NuSiLoadSiO2) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%SiO2Smallm(LimnoParam%NuSiLoadSiO2) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%SiO2Smallm(LimnoParam%NuSiLoadSiO2) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%SiO2CapitalM(LimnoParam%NuSiLoadSiO2) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%SiO2CapitalM(LimnoParam%NuSiLoadSiO2) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%SiO2Time(LimnoParam%NuSiLoadSiO2,1) = IniTime
                                LimnoParam%SiO2Time(LimnoParam%NuSiLoadSiO2,2) = FinalTime
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                LimnoParam%SiO2Value(LimnoParam%NuSiLoadSiO2,:) = hydrotimeSeries(1)%value1     
                            Else
                                LimnoParam%SiO2nTime(LimnoParam%NuSiLoadSiO2) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%SiO2Smallm(LimnoParam%NuSiLoadSiO2) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%SiO2Smallm(LimnoParam%NuSiLoadSiO2) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%SiO2CapitalM(LimnoParam%NuSiLoadSiO2) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%SiO2CapitalM(LimnoParam%NuSiLoadSiO2) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                Do j = 1, LimnoParam%SiO2nTime(LimnoParam%NuSiLoadSiO2)
                                    LimnoParam%SiO2Time(LimnoParam%NuSiLoadSiO2,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%SiO2Value(LimnoParam%NuSiLoadSiO2,j) = hydrotimeSeries(j)%value1 
                                EndDo
                            EndIf
                        EndDo
                    EndIf
                EndDo
            EndIf                         
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Dissolved Oxygen (O2)") Then !Silicate (SiO2)    
            
            If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                iRange = 1
                If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuO2LoadO2 = LimnoParam%NuO2LoadO2 + 1
                        LimnoParam%NO2Index(LimnoParam%NuO2LoadO2,1) = 1
                        LimnoParam%NO2Index(LimnoParam%NuO2LoadO2,2) = boundaryConditionCells(k)%cellId + 1
                        LimnoParam%O2nTime(LimnoParam%NuO2LoadO2) = 2
                        LimnoParam%O2Smallm(LimnoParam%NuO2LoadO2) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%O2CapitalM(LimnoParam%NuO2LoadO2) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%O2Time(LimnoParam%NuO2LoadO2,1) = IniTime
                        LimnoParam%O2Time(LimnoParam%NuO2LoadO2,2) = FinalTime
                        LimnoParam%O2Value(LimnoParam%NuO2LoadO2,:) = wqboundaryConditions(i)%constantValue
                    EndDo
                        
                ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuO2LoadO2 = LimnoParam%NuO2LoadO2 + 1
                        LimnoParam%NO2Index(LimnoParam%NuO2LoadO2,1) = 1
                        LimnoParam%NO2Index(LimnoParam%NuO2LoadO2,2) = boundaryConditionCells(k)%cellId + 1
                        If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                            LimnoParam%O2nTime(LimnoParam%NuO2LoadO2) = 2
                            LimnoParam%O2Smallm(LimnoParam%NuO2LoadO2) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%O2CapitalM(LimnoParam%NuO2LoadO2) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%O2Time(LimnoParam%NuO2LoadO2,1) = IniTime
                            LimnoParam%O2Time(LimnoParam%NuO2LoadO2,2) = FinalTime
                            LimnoParam%O2Value(LimnoParam%NuO2LoadO2,:) = wqboundaryConditions(i)%constantValue
                        Else
                            LimnoParam%O2nTime(LimnoParam%NuO2LoadO2) = wqboundaryConditions(i)%timeSeriesListSize
                            LimnoParam%O2Smallm(LimnoParam%NuO2LoadO2) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%O2CapitalM(LimnoParam%NuO2LoadO2) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                            Do j = 1, LimnoParam%O2nTime(LimnoParam%NuO2LoadO2)
                                LimnoParam%O2Time(LimnoParam%NuO2LoadO2,j) = hydrotimeSeries(j)%timeStamp
                                LimnoParam%O2Value(LimnoParam%NuO2LoadO2,j) = hydrotimeSeries(j)%value1 
                            EndDo                                   
                        EndIf
                    EndDo
                EndIf
                    
                    
            Else ! Non Vertically Integrated
                call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                Do iRange = 1,wqboundaryConditions(i)%rangesSize
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuO2LoadO2 = LimnoParam%NuO2LoadO2 + 1
                            LimnoParam%NO2Index(LimnoParam%NuO2LoadO2,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NO2Index(LimnoParam%NuO2LoadO2,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%O2nTime(LimnoParam%NuO2LoadO2) = 2
                            iElem = boundaryConditionCells(k)%cellId + 1
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                    LimnoParam%O2Smallm(LimnoParam%NuO2LoadO2) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%O2Smallm(LimnoParam%NuO2LoadO2) = HydroParam%ElSmallm(iElem)
                                EndIf
                            EndDo
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                    LimnoParam%O2CapitalM(LimnoParam%NuO2LoadO2) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%O2CapitalM(LimnoParam%NuO2LoadO2) = HydroParam%ElCapitalM(iElem)
                                EndIf
                            EndDo                                   
                            LimnoParam%O2Time(LimnoParam%NuO2LoadO2,1) = IniTime
                            LimnoParam%O2Time(LimnoParam%NuO2LoadO2,2) = FinalTime
                            LimnoParam%O2Value(LimnoParam%NuO2LoadO2,:) = nonVerticallyIntegratedCell(iRange)%Value
                        EndDo
                            
                    ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuO2LoadO2 = LimnoParam%NuO2LoadO2 + 1
                            LimnoParam%NO2Index(LimnoParam%NuO2LoadO2,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NO2Index(LimnoParam%NuO2LoadO2,2) = boundaryConditionCells(k)%cellId + 1
                            If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                LimnoParam%O2nTime(LimnoParam%NuO2LoadO2) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%O2Smallm(LimnoParam%NuO2LoadO2) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%O2Smallm(LimnoParam%NuO2LoadO2) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%O2CapitalM(LimnoParam%NuO2LoadO2) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%O2CapitalM(LimnoParam%NuO2LoadO2) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%O2Time(LimnoParam%NuO2LoadO2,1) = IniTime
                                LimnoParam%O2Time(LimnoParam%NuO2LoadO2,2) = FinalTime
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                LimnoParam%O2Value(LimnoParam%NuO2LoadO2,:) = hydrotimeSeries(1)%value1     
                            Else
                                LimnoParam%O2nTime(LimnoParam%NuO2LoadO2) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%O2Smallm(LimnoParam%NuO2LoadO2) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%O2Smallm(LimnoParam%NuO2LoadO2) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%O2CapitalM(LimnoParam%NuO2LoadO2) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%O2CapitalM(LimnoParam%NuO2LoadO2) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                Do j = 1, LimnoParam%O2nTime(LimnoParam%NuO2LoadO2)
                                    LimnoParam%O2Time(LimnoParam%NuO2LoadO2,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%O2Value(LimnoParam%NuO2LoadO2,j) = hydrotimeSeries(j)%value1 
                                EndDo
                            EndIf
                        EndDo
                    EndIf
                EndDo
            EndIf                 
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Biological Oxygen Demand (BOD)") Then !Biological Oxygen Demand (BOD)
            
            If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                iRange = 1
                If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuDBOLoadDBO = LimnoParam%NuDBOLoadDBO + 1
                        LimnoParam%NDBOIndex(LimnoParam%NuDBOLoadDBO,1) = 1
                        LimnoParam%NDBOIndex(LimnoParam%NuDBOLoadDBO,2) = boundaryConditionCells(k)%cellId + 1
                        LimnoParam%DBOnTime(LimnoParam%NuDBOLoadDBO) = 2
                        LimnoParam%DBOSmallm(LimnoParam%NuDBOLoadDBO) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%DBOCapitalM(LimnoParam%NuDBOLoadDBO) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%DBOTime(LimnoParam%NuDBOLoadDBO,1) = IniTime
                        LimnoParam%DBOTime(LimnoParam%NuDBOLoadDBO,2) = FinalTime
                        LimnoParam%DBOValue(LimnoParam%NuDBOLoadDBO,:) = wqboundaryConditions(i)%constantValue
                    EndDo
                        
                ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuDBOLoadDBO = LimnoParam%NuDBOLoadDBO + 1
                        LimnoParam%NDBOIndex(LimnoParam%NuDBOLoadDBO,1) = 1
                        LimnoParam%NDBOIndex(LimnoParam%NuDBOLoadDBO,2) = boundaryConditionCells(k)%cellId + 1
                        If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                            LimnoParam%DBOnTime(LimnoParam%NuDBOLoadDBO) = 2
                            LimnoParam%DBOSmallm(LimnoParam%NuDBOLoadDBO) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%DBOCapitalM(LimnoParam%NuDBOLoadDBO) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%DBOTime(LimnoParam%NuDBOLoadDBO,1) = IniTime
                            LimnoParam%DBOTime(LimnoParam%NuDBOLoadDBO,2) = FinalTime
                            LimnoParam%DBOValue(LimnoParam%NuDBOLoadDBO,:) = wqboundaryConditions(i)%constantValue
                        Else
                            LimnoParam%DBOnTime(LimnoParam%NuDBOLoadDBO) = wqboundaryConditions(i)%timeSeriesListSize
                            LimnoParam%DBOSmallm(LimnoParam%NuDBOLoadDBO) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%DBOCapitalM(LimnoParam%NuDBOLoadDBO) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                            Do j = 1, LimnoParam%DBOnTime(LimnoParam%NuDBOLoadDBO)
                                LimnoParam%DBOTime(LimnoParam%NuDBOLoadDBO,j) = hydrotimeSeries(j)%timeStamp
                                LimnoParam%DBOValue(LimnoParam%NuDBOLoadDBO,j) = hydrotimeSeries(j)%value1 
                            EndDo                                   
                        EndIf
                    EndDo
                EndIf
                    
                    
            Else ! Non Vertically Integrated
                call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                Do iRange = 1,wqboundaryConditions(i)%rangesSize
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDBOLoadDBO = LimnoParam%NuDBOLoadDBO + 1
                            LimnoParam%NDBOIndex(LimnoParam%NuDBOLoadDBO,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NDBOIndex(LimnoParam%NuDBOLoadDBO,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%DBOnTime(LimnoParam%NuDBOLoadDBO) = 2
                            iElem = boundaryConditionCells(k)%cellId + 1
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                    LimnoParam%DBOSmallm(LimnoParam%NuDBOLoadDBO) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%DBOSmallm(LimnoParam%NuDBOLoadDBO) = HydroParam%ElSmallm(iElem)
                                EndIf
                            EndDo
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                    LimnoParam%DBOCapitalM(LimnoParam%NuDBOLoadDBO) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%DBOCapitalM(LimnoParam%NuDBOLoadDBO) = HydroParam%ElCapitalM(iElem)
                                EndIf
                            EndDo                                   
                            LimnoParam%DBOTime(LimnoParam%NuDBOLoadDBO,1) = IniTime
                            LimnoParam%DBOTime(LimnoParam%NuDBOLoadDBO,2) = FinalTime
                            LimnoParam%DBOValue(LimnoParam%NuDBOLoadDBO,:) = nonVerticallyIntegratedCell(iRange)%Value
                        EndDo
                            
                    ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDBOLoadDBO = LimnoParam%NuDBOLoadDBO + 1
                            LimnoParam%NDBOIndex(LimnoParam%NuDBOLoadDBO,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NDBOIndex(LimnoParam%NuDBOLoadDBO,2) = boundaryConditionCells(k)%cellId + 1
                            If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                LimnoParam%DBOnTime(LimnoParam%NuDBOLoadDBO) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%DBOSmallm(LimnoParam%NuDBOLoadDBO) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%DBOSmallm(LimnoParam%NuDBOLoadDBO) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%DBOCapitalM(LimnoParam%NuDBOLoadDBO) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%DBOCapitalM(LimnoParam%NuDBOLoadDBO) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%DBOTime(LimnoParam%NuDBOLoadDBO,1) = IniTime
                                LimnoParam%DBOTime(LimnoParam%NuDBOLoadDBO,2) = FinalTime
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                LimnoParam%DBOValue(LimnoParam%NuDBOLoadDBO,:) = hydrotimeSeries(1)%value1     
                            Else
                                LimnoParam%DBOnTime(LimnoParam%NuDBOLoadDBO) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%DBOSmallm(LimnoParam%NuDBOLoadDBO) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%DBOSmallm(LimnoParam%NuDBOLoadDBO) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%DBOCapitalM(LimnoParam%NuDBOLoadDBO) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%DBOCapitalM(LimnoParam%NuDBOLoadDBO) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                Do j = 1, LimnoParam%DBOnTime(LimnoParam%NuDBOLoadDBO)
                                    LimnoParam%DBOTime(LimnoParam%NuDBOLoadDBO,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%DBOValue(LimnoParam%NuDBOLoadDBO,j) = hydrotimeSeries(j)%value1 
                                EndDo
                            EndIf
                        EndDo
                    EndIf
                EndDo
            EndIf            
            
        ElseIf (trim(wqboundaryConditions(i)%conditionType) == "Dissolved Inorganic Carbon") Then !Dissolved Inorganic Carbon
            
            If (wqboundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                iRange = 1
                If (wqboundaryConditions(i)%conditionFunction == 1) Then !Constant
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuDicLoadDic = LimnoParam%NuDicLoadDic + 1
                        LimnoParam%NDicIndex(LimnoParam%NuDicLoadDic,1) = 1
                        LimnoParam%NDicIndex(LimnoParam%NuDicLoadDic,2) = boundaryConditionCells(k)%cellId + 1
                        LimnoParam%DicnTime(LimnoParam%NuDicLoadDic) = 2
                        LimnoParam%DicSmallm(LimnoParam%NuDicLoadDic) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%DicCapitalM(LimnoParam%NuDicLoadDic) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                        LimnoParam%DicTime(LimnoParam%NuDicLoadDic,1) = IniTime
                        LimnoParam%DicTime(LimnoParam%NuDicLoadDic,2) = FinalTime
                        LimnoParam%DicValue(LimnoParam%NuDicLoadDic,:) = wqboundaryConditions(i)%constantValue
                    EndDo
                        
                ElseIf (wqboundaryConditions(i)%conditionFunction == 2) Then !Time-series
                    Do k = 1,wqboundaryConditions(i)%cellsLength
                        LimnoParam%NuDicLoadDic = LimnoParam%NuDicLoadDic + 1
                        LimnoParam%NDicIndex(LimnoParam%NuDicLoadDic,1) = 1
                        LimnoParam%NDicIndex(LimnoParam%NuDicLoadDic,2) = boundaryConditionCells(k)%cellId + 1
                        If (wqboundaryConditions(i)%timeSeriesListSize<2) Then
                            LimnoParam%DicnTime(LimnoParam%NuDicLoadDic) = 2
                            LimnoParam%DicSmallm(LimnoParam%NuDicLoadDic) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%DicCapitalM(LimnoParam%NuDicLoadDic) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%DicTime(LimnoParam%NuDicLoadDic,1) = IniTime
                            LimnoParam%DicTime(LimnoParam%NuDicLoadDic,2) = FinalTime
                            LimnoParam%DicValue(LimnoParam%NuDicLoadDic,:) = wqboundaryConditions(i)%constantValue
                        Else
                            LimnoParam%DicnTime(LimnoParam%NuDicLoadDic) = wqboundaryConditions(i)%timeSeriesListSize
                            LimnoParam%DicSmallm(LimnoParam%NuDicLoadDic) = HydroParam%ElSmallm(boundaryConditionCells(k)%cellId + 1)
                            LimnoParam%DicCapitalM(LimnoParam%NuDicLoadDic) = HydroParam%ElCapitalM(boundaryConditionCells(k)%cellId + 1)
                            call c_f_pointer(wqboundaryConditions(i)%timeSeriesList, hydrotimeSeries, [wqboundaryConditions(i)%timeSeriesListSize])
                            Do j = 1, LimnoParam%DicnTime(LimnoParam%NuDicLoadDic)
                                LimnoParam%DicTime(LimnoParam%NuDicLoadDic,j) = hydrotimeSeries(j)%timeStamp
                                LimnoParam%DicValue(LimnoParam%NuDicLoadDic,j) = hydrotimeSeries(j)%value1 
                            EndDo                                   
                        EndIf
                    EndDo
                EndIf
                    
                    
            Else ! Non Vertically Integrated
                call c_f_pointer(wqboundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [wqboundaryConditions(i)%rangesSize])
                Do iRange = 1,wqboundaryConditions(i)%rangesSize
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDicLoadDic = LimnoParam%NuDicLoadDic + 1
                            LimnoParam%NDicIndex(LimnoParam%NuDicLoadDic,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NDicIndex(LimnoParam%NuDicLoadDic,2) = boundaryConditionCells(k)%cellId + 1
                            LimnoParam%DicnTime(LimnoParam%NuDicLoadDic) = 2
                            iElem = boundaryConditionCells(k)%cellId + 1
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                    LimnoParam%DicSmallm(LimnoParam%NuDicLoadDic) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%DicSmallm(LimnoParam%NuDicLoadDic) = HydroParam%ElSmallm(iElem)
                                EndIf
                            EndDo
                            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                    LimnoParam%DicCapitalM(LimnoParam%NuDicLoadDic) = iLayer !Vertical Integrated (review)
                                    Exit
                                Else
                                    LimnoParam%DicCapitalM(LimnoParam%NuDicLoadDic) = HydroParam%ElCapitalM(iElem)
                                EndIf
                            EndDo                                   
                            LimnoParam%DicTime(LimnoParam%NuDicLoadDic,1) = IniTime
                            LimnoParam%DicTime(LimnoParam%NuDicLoadDic,2) = FinalTime
                            LimnoParam%DicValue(LimnoParam%NuDicLoadDic,:) = nonVerticallyIntegratedCell(iRange)%Value
                        EndDo
                            
                    ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                        Do k = 1,wqboundaryConditions(i)%cellsLength
                            LimnoParam%NuDicLoadDic = LimnoParam%NuDicLoadDic + 1
                            LimnoParam%NDicIndex(LimnoParam%NuDicLoadDic,1) = wqboundaryConditions(i)%rangesSize
                            LimnoParam%NDicIndex(LimnoParam%NuDicLoadDic,2) = boundaryConditionCells(k)%cellId + 1
                            If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                LimnoParam%DicnTime(LimnoParam%NuDicLoadDic) = 2
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%DicSmallm(LimnoParam%NuDicLoadDic) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%DicSmallm(LimnoParam%NuDicLoadDic) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%DicCapitalM(LimnoParam%NuDicLoadDic) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%DicCapitalM(LimnoParam%NuDicLoadDic) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                LimnoParam%DicTime(LimnoParam%NuDicLoadDic,1) = IniTime
                                LimnoParam%DicTime(LimnoParam%NuDicLoadDic,2) = FinalTime
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                LimnoParam%DicValue(LimnoParam%NuDicLoadDic,:) = hydrotimeSeries(1)%value1     
                            Else
                                LimnoParam%DicnTime(LimnoParam%NuDicLoadDic) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                iElem = boundaryConditionCells(k)%cellId + 1
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
                                    If (HydroParam%Ze(iLayer,iElem)>=nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                        LimnoParam%DicSmallm(LimnoParam%NuDicLoadDic) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%DicSmallm(LimnoParam%NuDicLoadDic) = HydroParam%ElSmallm(iElem)
                                    EndIf
                                EndDo
                                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                                    If (HydroParam%Ze(iLayer+1,iElem)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                        LimnoParam%DicCapitalM(LimnoParam%NuDicLoadDic) = iLayer !Vertical Integrated (review)
                                        Exit
                                    Else
                                        LimnoParam%DicCapitalM(LimnoParam%NuDicLoadDic) = HydroParam%ElCapitalM(iElem)
                                    EndIf
                                EndDo                                   
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                Do j = 1, LimnoParam%DicnTime(LimnoParam%NuDicLoadDic)
                                    LimnoParam%DicTime(LimnoParam%NuDicLoadDic,j) = hydrotimeSeries(j)%timeStamp
                                    LimnoParam%DicValue(LimnoParam%NuDicLoadDic,j) = hydrotimeSeries(j)%value1 
                                EndDo
                            EndIf
                        EndDo
                    EndIf
                EndDo
            EndIf                      
                
        EndIf
        
    EndDo
    
End Subroutine ReadWQBoundaryCondition
    
    