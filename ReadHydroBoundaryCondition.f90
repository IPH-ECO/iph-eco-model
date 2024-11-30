!> This subroutine reads the simulation parameters. 
Subroutine ReadHydroBoundaryCondition(HydroParam,hydroConfiguration,IniTime,FinalTime,MeshParam)
    
    Use domain_types
    Use Hydrodynamic
    Use MeshVars
    
    Implicit none
    Integer i,j,k,iEdge,iNode1,iNode2,Face,iLayer,iRange,nWaterLvl
    Integer:: idate(6)
    character(len=20):: str
    Integer NObjInflow, NObjWaterLevel, NObjRange
    Integer(c_int):: numberOfBoundaryConditions
    type(HydrodynamicConfiguration) :: hydroConfiguration
    type(BoundaryCondition), dimension(:), pointer :: boundaryConditions
    type(NonVerticallyIntegratedRange), pointer :: nonVerticallyIntegratedCell(:)
    !integer, pointer :: objectIds(:)
    type(HydrodynamicBoundaryConditionCell), pointer :: boundaryConditionCells(:)
    type(TimeSeries), pointer :: hydrotimeSeries(:)
    Integer:: IniTime !< Initial time of simulation (unix time)
    Integer:: FinalTime !< Final time of simulation (unix time)
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    
    call c_f_pointer(hydroConfiguration%boundaryConditions, boundaryConditions, [hydroConfiguration%numberOfBoundaryConditions])
    
    
    !Allocate(HydroParam%IndexInflow(HydroParam%NInflow,2))
    Allocate(HydroParam%IndexInflow(MeshParam%nEdge,2))
    HydroParam%IndexInflow = 0
    
    
    HydroParam%IndexWaterLevel = 0.d0
    
    HydroParam%NInflow = 0.d0
    HydroParam%NWaterLevel = 0.d0
    NObjWaterLevel = 0.d0
    NObjInflow = 0.d0
    NObjRange = 0.d0
    
    ! 1. Mapping cells (Water level) or faces (Inflow/outflow) with boundary condition
    Do i = 1,hydroConfiguration%numberOfBoundaryConditions
        call c_f_pointer(boundaryConditions(i)%cells, boundaryConditionCells, [boundaryConditions(i)%cellsLength])
        
        If (trim(boundaryConditions(i)%conditionType) == "waterLevel") Then !Water level
            Do j = 1,boundaryConditions(i)%cellsLength
                HydroParam%NWaterLevel = HydroParam%NWaterLevel + 1
                NObjWaterLevel = Max(NObjWaterLevel, boundaryConditions(i)%timeSeriesListSize)
                
            EndDo
        ElseIf (trim(boundaryConditions(i)%conditionType) == "waterFlow") Then !Inflow/Outflow
            Do j = 1,boundaryConditions(i)%cellsLength
                If (boundaryConditions(i)%verticallyIntegrated) Then
                    HydroParam%NInflow = HydroParam%NInflow + 1
                    NObjRange = Max(NObjRange,1)
                    NObjInflow = Max(NObjInflow, boundaryConditions(i)%timeSeriesListSize)
                    
                Else
                    NObjRange = Max(NObjRange,boundaryConditions(i)%rangesSize)
                    call c_f_pointer(boundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [boundaryConditions(i)%rangesSize])
 

                    Do iRange = 1,boundaryConditions(i)%rangesSize
                        HydroParam%NInflow = HydroParam%NInflow + 1
                            
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then ! Constant
                            NObjInflow = Max(NObjInflow, 1)
                        ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then ! Time-series
                            NObjInflow = Max(NObjInflow, nonVerticallyIntegratedCell(iRange)%timeSeriesListSize)
                        EndIf
                    EndDo
                EndIf
                    
            EndDo
            
            
            
        ElseIf (trim(boundaryConditions(i)%conditionType) == "normalDepth") Then !Normal Depth
            Do j = 1,boundaryConditions(i)%cellsLength
                HydroParam%NInflow = HydroParam%NInflow + 1
                !NObjInflow = Max(NObjInflow, boundaryConditions(i)%timeSeriesListSize)
                Do iEdge = 1,4
                    Face = MeshParam%Edge(iEdge,boundaryConditionCells(j)%cellId + 1)
                    iNode1 = MeshParam%Quadri(MeshParam%EdgeDef(1,iEdge),boundaryConditionCells(j)%cellId + 1)         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                    iNode2 = MeshParam%Quadri(MeshParam%EdgeDef(2,iEdge),boundaryConditionCells(j)%cellId + 1)         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                    If (iNode1==boundaryConditionCells(j)%verticeIds(1).and.iNode2==boundaryConditionCells(j)%verticeIds(2)) Then
                        HydroParam%IndexInflow(Face,1) = HydroParam%NInflow
                        HydroParam%IndexInflow(Face,2) = boundaryConditionCells(j)%cellId + 1
                        Exit
                    EndIf
                EndDo
            EndDo
        EndIf
        
    EndDo
    
    ! 2. Reading Hydrodynamic Boundary conditions
    If (NObjWaterLevel<2) NObjWaterLevel = 2
    Allocate(HydroParam%WaterLevelValue(HydroParam%NWaterLevel,NObjWaterLevel))
    Allocate (HydroParam%WaterLevelTime(HydroParam%NWaterLevel,NObjWaterLevel))
    Allocate (HydroParam%WaterLevelnTime(HydroParam%NWaterLevel))
    Allocate(HydroParam%IndexWaterLevel(HydroParam%NWaterLevel,2))
    Allocate (HydroParam%WaterLevel(HydroParam%NWaterLevel))
    If (NObjInflow<2) NObjInflow = 2
    Allocate (HydroParam%InFlowValue(HydroParam%NInflow,NObjInflow))
    Allocate (HydroParam%InFlowTime(HydroParam%NInflow,NObjInflow))
    Allocate (HydroParam%InFlownTime(HydroParam%NInflow))
    Allocate (HydroParam%InFlowSmallm(HydroParam%NInflow))
    Allocate (HydroParam%InFlowCapitalM(HydroParam%NInflow))
    Allocate (HydroParam%NRange(HydroParam%NInflow)) 
    
    HydroParam%IndexWaterLevel = 0
    HydroParam%IndexWaterLevelEdge = 0
    
    HydroParam%IndexInflowEdge = 0
    HydroParam%NInflow = 0
    HydroParam%NWaterLevel = 0
    !NObjRange = 0
    nWaterLvl = 0
    
    Do i = 1,hydroConfiguration%numberOfBoundaryConditions
        call c_f_pointer(boundaryConditions(i)%cells, boundaryConditionCells, [boundaryConditions(i)%cellsLength])
            
        
        If (trim(boundaryConditions(i)%conditionType) == "waterLevel") Then !Water level
            If (boundaryConditions(i)%conditionFunction == 1) Then !Constant
                Do k = 1,boundaryConditions(i)%cellsLength
                    HydroParam%NWaterLevel = HydroParam%NWaterLevel + 1
                    HydroParam%WaterLevelnTime(HydroParam%NWaterLevel) = 2
                    Do iEdge = 1,4
                        Face = MeshParam%Edge(iEdge,boundaryConditionCells(k)%cellId + 1)
                        iNode1 = MeshParam%Quadri(MeshParam%EdgeDef(1,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                        iNode2 = MeshParam%Quadri(MeshParam%EdgeDef(2,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                        If (iNode1==boundaryConditionCells(k)%verticeIds(1).and.iNode2==boundaryConditionCells(k)%verticeIds(2)) Then
                            HydroParam%IndexWaterLevel(HydroParam%NWaterLevel,1) = Face
                            HydroParam%IndexWaterLevel(HydroParam%NWaterLevel,2) = boundaryConditionCells(k)%cellId + 1
                            nWaterLvl = nWaterLvl + 1
                            HydroParam%IndexWaterLevelEdge(Face) = nWaterLvl
                            !HydroParam%IndexWaterLevelEdge(Face) = 1
                        EndIf
                    EndDo
                    
                    HydroParam%WaterLevelTime(HydroParam%NWaterLevel,1) = IniTime
                    HydroParam%WaterLevelTime(HydroParam%NWaterLevel,2) = FinalTime
                    HydroParam%WaterLevelValue(HydroParam%NWaterLevel,:) = boundaryConditions(i)%constantValue
                    
                EndDo
            ElseIf (boundaryConditions(i)%conditionFunction == 2) Then !Time-series
                Do k = 1,boundaryConditions(i)%cellsLength
                HydroParam%NWaterLevel = HydroParam%NWaterLevel + 1
                If (boundaryConditions(i)%timeSeriesListSize<2) Then
                    HydroParam%WaterLevelnTime(HydroParam%NWaterLevel) = 2
                    
                    Do iEdge = 1,4
                        Face = MeshParam%Edge(iEdge,boundaryConditionCells(k)%cellId + 1)
                        iNode1 = MeshParam%Quadri(MeshParam%EdgeDef(1,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                        iNode2 = MeshParam%Quadri(MeshParam%EdgeDef(2,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                        If (iNode1==boundaryConditionCells(k)%verticeIds(1).and.iNode2==boundaryConditionCells(k)%verticeIds(2)) Then
                            HydroParam%IndexWaterLevel(HydroParam%NWaterLevel,1) = Face
                            HydroParam%IndexWaterLevel(HydroParam%NWaterLevel,2) = boundaryConditionCells(k)%cellId + 1
                            nWaterLvl = nWaterLvl + 1
                            HydroParam%IndexWaterLevelEdge(Face) = nWaterLvl
                            !HydroParam%IndexWaterLevelEdge(Face) = 1
                        EndIf
                    EndDo 
                    
                    HydroParam%WaterLevelTime(HydroParam%NWaterLevel,1) = IniTime
                    HydroParam%WaterLevelTime(HydroParam%NWaterLevel,2) = FinalTime
                    call c_f_pointer(boundaryConditions(i)%timeSeriesList, hydrotimeSeries, [boundaryConditions(i)%timeSeriesListSize])
                    HydroParam%WaterLevelValue(HydroParam%NWaterLevel,:) = hydrotimeSeries(1)%value1
                Else
                    HydroParam%WaterLevelnTime(HydroParam%NWaterLevel) = boundaryConditions(i)%timeSeriesListSize
                    
                    Do iEdge = 1,4
                        Face = MeshParam%Edge(iEdge,boundaryConditionCells(k)%cellId + 1)
                        iNode1 = MeshParam%Quadri(MeshParam%EdgeDef(1,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                        iNode2 = MeshParam%Quadri(MeshParam%EdgeDef(2,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                        If (iNode1==boundaryConditionCells(k)%verticeIds(1).and.iNode2==boundaryConditionCells(k)%verticeIds(2)) Then
                            HydroParam%IndexWaterLevel(HydroParam%NWaterLevel,1) = Face
                            HydroParam%IndexWaterLevel(HydroParam%NWaterLevel,2) = boundaryConditionCells(k)%cellId + 1
                            nWaterLvl = nWaterLvl + 1
                            HydroParam%IndexWaterLevelEdge(Face) = nWaterLvl
                            !HydroParam%IndexWaterLevelEdge(Face) = 1
                        EndIf
                    EndDo

                    call c_f_pointer(boundaryConditions(i)%timeSeriesList, hydrotimeSeries, [boundaryConditions(i)%timeSeriesListSize])
                        Do j = 1, HydroParam%WaterLevelnTime(HydroParam%NWaterLevel)
                            ! Verifing consistency of time series
                            !If (j == 1) Then !First date of time series
                            !    If (IniTime<hydrotimeSeries(j)%timeStampSize) Then
                            !        Call unix2c(hydrotimeSeries(j)%timeStampSize, idate)
                            !        Print*, 'First date of time series ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//')'
                            !        Call unix2c(IniTime, idate)
                            !        Print*, 'is greater than initial time of simulation ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//').' 
                            !        Print*, 'Please check Boundary condition number ('//trim(str(i))//').'
                            !        Pause
                            !        Stop 
                            !    EndIf
                            !ElseIf (j == WaterLevelnTime(NWaterLevel)) Then
                            !    If (FinalTime>hydrotimeSeries(j)%timeStampSize) Then
                            !        Call unix2c(hydrotimeSeries(j)%timeStampSize, idate)
                            !        Print*, 'Last date of time series ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//')'
                            !        Call unix2c(FinalTime, idate)
                            !        Print*, 'is lesser than initial time of simulation ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//').' 
                            !        Print*, 'Please check Boundary condition number ('//trim(str(i))//').'
                            !        Pause
                            !        Stop 
                            !    EndIf
                            !EndIf
                            ! Storing stamp time and values of time series 
                            HydroParam%WaterLevelTime(HydroParam%NWaterLevel,j) = hydrotimeSeries(j)%timeStamp
                            HydroParam%WaterLevelValue(HydroParam%NWaterLevel,j) = hydrotimeSeries(j)%value1
                        EndDo
                    
                    EndIf
                EndDo
            EndIf
        ElseIf (trim(boundaryConditions(i)%conditionType) == "waterFlow") Then !Inflow/Outflow
            
            If (boundaryConditions(i)%verticallyIntegrated) Then ! Vertically Integrated
                !HydroParam%NRange(i) = 1
                !iRange = 1
                If (boundaryConditions(i)%conditionFunction == 1) Then !Constant
                    Do k = 1,boundaryConditions(i)%cellsLength
                        HydroParam%NInflow = HydroParam%NInflow + 1
                        HydroParam%InFlownTime(HydroParam%NInflow) = 2
                        Do iEdge = 1,4
                            Face = MeshParam%Edge(iEdge,boundaryConditionCells(k)%cellId + 1)
                            iNode1 = MeshParam%Quadri(MeshParam%EdgeDef(1,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                            iNode2 = MeshParam%Quadri(MeshParam%EdgeDef(2,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                            If (iNode1==boundaryConditionCells(k)%verticeIds(1).and.iNode2==boundaryConditionCells(k)%verticeIds(2)) Then
                                HydroParam%InFlowSmallm(HydroParam%NInflow) = HydroParam%Smallm(Face)
                                HydroParam%InFlowCapitalM(HydroParam%NInflow) = HydroParam%CapitalM(Face)
                                HydroParam%IndexInflow(HydroParam%NInflow,1) = Face
                                HydroParam%IndexInflow(HydroParam%NInflow,2) = boundaryConditionCells(k)%cellId + 1
                                HydroParam%IndexInflowEdge(Face) = 1
                            EndIf
                        EndDo
                        HydroParam%InFlowTime(HydroParam%NInflow,1) = IniTime
                        HydroParam%InFlowTime(HydroParam%NInflow,2) = FinalTime
                        HydroParam%InFlowValue(HydroParam%NInflow,:) = boundaryConditions(i)%constantValue
                    EndDo
                ElseIf (boundaryConditions(i)%conditionFunction == 2) Then !Time-series
                    Do k = 1,boundaryConditions(i)%cellsLength
                        HydroParam%NInflow = HydroParam%NInflow + 1
                        If (boundaryConditions(i)%timeSeriesListSize<2) Then
                            HydroParam%InFlownTime(HydroParam%NInflow) = 2
                            Do iEdge = 1,4
                                Face = MeshParam%Edge(iEdge,boundaryConditionCells(k)%cellId + 1)
                                iNode1 = MeshParam%Quadri(MeshParam%EdgeDef(1,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                                iNode2 = MeshParam%Quadri(MeshParam%EdgeDef(2,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                                If (iNode1==boundaryConditionCells(k)%verticeIds(1).and.iNode2==boundaryConditionCells(k)%verticeIds(2)) Then
                                    HydroParam%InFlowSmallm(HydroParam%NInflow) = HydroParam%Smallm(Face)
                                    HydroParam%InFlowCapitalM(HydroParam%NInflow) = HydroParam%CapitalM(Face)
                                    HydroParam%IndexInflow(HydroParam%NInflow,1) = Face
                                    HydroParam%IndexInflow(HydroParam%NInflow,2) = boundaryConditionCells(k)%cellId + 1
                                    HydroParam%IndexInflowEdge(Face) = 1
                                EndIf
                            EndDo 
                            HydroParam%InFlowTime(HydroParam%NInflow,1) = IniTime
                            HydroParam%InFlowTime(HydroParam%NInflow,2) = FinalTime
                            call c_f_pointer(boundaryConditions(i)%timeSeriesList, hydrotimeSeries, [boundaryConditions(i)%timeSeriesListSize])
                            HydroParam%InFlowValue(HydroParam%NInflow,:) = hydrotimeSeries(1)%value1
                        Else
                            HydroParam%InFlownTime(HydroParam%NInflow) = boundaryConditions(i)%timeSeriesListSize
                            Do iEdge = 1,4
                                Face = MeshParam%Edge(iEdge,boundaryConditionCells(k)%cellId + 1)
                                iNode1 = MeshParam%Quadri(MeshParam%EdgeDef(1,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                                iNode2 = MeshParam%Quadri(MeshParam%EdgeDef(2,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                                If (iNode1==boundaryConditionCells(k)%verticeIds(1).and.iNode2==boundaryConditionCells(k)%verticeIds(2)) Then
                                    HydroParam%InFlowSmallm(HydroParam%NInflow) = HydroParam%Smallm(Face)
                                    HydroParam%InFlowCapitalM(HydroParam%NInflow) = HydroParam%CapitalM(Face)
                                    HydroParam%IndexInflow(HydroParam%NInflow,1) = Face
                                    HydroParam%IndexInflow(HydroParam%NInflow,2) = boundaryConditionCells(k)%cellId + 1
                                    HydroParam%IndexInflowEdge(Face) = 1
                                EndIf
                            EndDo
                            call c_f_pointer(boundaryConditions(i)%timeSeriesList, hydrotimeSeries, [boundaryConditions(i)%timeSeriesListSize])
                            Do j = 1, HydroParam%InFlownTime(HydroParam%NInflow)
                                ! Verifing consistency of time series
                                !If (j == 1) Then !First date of time series
                                !    If (IniTime<hydrotimeSeries(j)%timeStampSize) Then
                                !        Call unix2c(hydrotimeSeries(j)%timeStampSize, idate)
                                !        Print*, 'First date of time series ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//')'
                                !        Call unix2c(IniTime, idate)
                                !        Print*, 'is greater than initial time of simulation ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//').' 
                                !        Print*, 'Please check Boundary condition number ('//trim(str(i))//').'
                                !        Pause
                                !        Stop 
                                !    EndIf
                                !ElseIf (j == WaterLevelnTime(NWaterLevel)) Then
                                !    If (FinalTime>hydrotimeSeries(j)%timeStampSize) Then
                                !        Call unix2c(hydrotimeSeries(j)%timeStampSize, idate)
                                !        Print*, 'Last date of time series ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//')'
                                !        Call unix2c(FinalTime, idate)
                                !        Print*, 'is lesser than initial time of simulation ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//').' 
                                !        Print*, 'Please check Boundary condition number ('//trim(str(i))//').'
                                !        Pause
                                !        Stop 
                                !    EndIf
                                !EndIf
                                ! Storing stamp time and values of time series 
                                HydroParam%InFlowTime(HydroParam%NInflow,j) = hydrotimeSeries(j)%timeStamp
                                HydroParam%InFlowValue(HydroParam%NInflow,j) = hydrotimeSeries(j)%value1 
                            EndDo                            
                        EndIf
                    EndDo
                EndIf
            Else ! Non Vertically Integrated
                call c_f_pointer(boundaryConditions(i)%ranges, nonVerticallyIntegratedCell, [boundaryConditions(i)%rangesSize])
                !HydroParam%NRange(i) = boundaryConditions(i)%rangesSize
                Do iRange = 1,boundaryConditions(i)%rangesSize
                    If (nonVerticallyIntegratedCell(iRange)%function == 1) Then !Constant
                        Do k = 1,boundaryConditions(i)%cellsLength
                            HydroParam%NInflow = HydroParam%NInflow + 1
                            HydroParam%InFlownTime(HydroParam%NInflow) = 2
                            Do iEdge = 1,4
                                Face = MeshParam%Edge(iEdge,boundaryConditionCells(k)%cellId + 1)
                                iNode1 = MeshParam%Quadri(MeshParam%EdgeDef(1,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                                iNode2 = MeshParam%Quadri(MeshParam%EdgeDef(2,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                                If (iNode1==boundaryConditionCells(k)%verticeIds(1).and.iNode2==boundaryConditionCells(k)%verticeIds(2)) Then
                                    HydroParam%IndexInflow(HydroParam%NInflow,1) = Face
                                    HydroParam%IndexInflow(HydroParam%NInflow,2) = boundaryConditionCells(k)%cellId + 1
                                    HydroParam%IndexInflowEdge(Face) = 1
                                    Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face) 
                                        If (HydroParam%Z(iLayer+1,Face)>nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                            HydroParam%InFlowSmallm(HydroParam%NInflow) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            HydroParam%InFlowSmallm(HydroParam%NInflow) = HydroParam%Smallm(Face)
                                        EndIf
                                    EndDo
                                    Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face)
                                        If (HydroParam%Z(iLayer+1,Face)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                            HydroParam%InFlowCapitalM(HydroParam%NInflow) = iLayer !Vertical Integrated (review)
                                            Exit
                                        Else
                                            HydroParam%InFlowCapitalM(HydroParam%NInflow) = HydroParam%CapitalM(Face)
                                        EndIf
                                    EndDo     

                                    Exit
                                EndIf
                            EndDo  
                        EndDo
                        HydroParam%InFlowTime(HydroParam%NInflow,1) = IniTime
                        HydroParam%InFlowTime(HydroParam%NInflow,2) = FinalTime
                        HydroParam%InFlowValue(HydroParam%NInflow,:) = nonVerticallyIntegratedCell(iRange)%value
                    ElseIf (nonVerticallyIntegratedCell(iRange)%function == 2) Then !Time-series
                        Do k = 1,boundaryConditions(i)%cellsLength
                            HydroParam%NInflow = HydroParam%NInflow + 1
                            If (nonVerticallyIntegratedCell(iRange)%timeSeriesListSize<2) Then
                                HydroParam%InFlownTime(HydroParam%NInflow) = 2
                                Do iEdge = 1,4
                                    Face = MeshParam%Edge(iEdge,boundaryConditionCells(k)%cellId + 1)
                                    iNode1 = MeshParam%Quadri(MeshParam%EdgeDef(1,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                                    iNode2 = MeshParam%Quadri(MeshParam%EdgeDef(2,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                                    If (iNode1==boundaryConditionCells(k)%verticeIds(1).and.iNode2==boundaryConditionCells(k)%verticeIds(2)) Then
                                        HydroParam%IndexInflow(HydroParam%NInflow,1) = Face
                                        HydroParam%IndexInflow(HydroParam%NInflow,2) = boundaryConditionCells(k)%cellId + 1
                                        HydroParam%IndexInflowEdge(Face) = 1
                                        Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face) 
                                            If (HydroParam%Z(iLayer+1,Face)>nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                                HydroParam%InFlowSmallm(HydroParam%NInflow) = iLayer !Vertical Integrated (review)
                                                Exit
                                            Else
                                                HydroParam%InFlowSmallm(HydroParam%NInflow) = HydroParam%Smallm(Face)
                                            EndIf
                                        EndDo
                                        Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face)
                                            If (HydroParam%Z(iLayer+1,Face)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                                HydroParam%InFlowCapitalM(HydroParam%NInflow) = iLayer !Vertical Integrated (review)
                                                Exit
                                            Else
                                                HydroParam%InFlowCapitalM(HydroParam%NInflow) = HydroParam%CapitalM(Face)
                                            EndIf
                                        EndDo                                    
                                        Exit
                                    EndIf
                                EndDo  
                                HydroParam%InFlowTime(HydroParam%NInflow,1) = IniTime
                                HydroParam%InFlowTime(HydroParam%NInflow,2) = FinalTime
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                HydroParam%InFlowValue(HydroParam%NInflow,:) = hydrotimeSeries(1)%value1 
                            Else
                                HydroParam%InFlownTime(HydroParam%NInflow) = nonVerticallyIntegratedCell(iRange)%timeSeriesListSize
                                Do iEdge = 1,4
                                    Face = MeshParam%Edge(iEdge,boundaryConditionCells(k)%cellId + 1)
                                    iNode1 = MeshParam%Quadri(MeshParam%EdgeDef(1,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                                    iNode2 = MeshParam%Quadri(MeshParam%EdgeDef(2,iEdge),boundaryConditionCells(k)%cellId + 1)         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                                    If (iNode1==boundaryConditionCells(k)%verticeIds(1).and.iNode2==boundaryConditionCells(k)%verticeIds(2)) Then
                                        HydroParam%IndexInflow(HydroParam%NInflow,1) = Face
                                        HydroParam%IndexInflow(HydroParam%NInflow,2) = boundaryConditionCells(k)%cellId + 1
                                        HydroParam%IndexInflowEdge(Face) = 1
                                        Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face) 
                                            If (HydroParam%Z(iLayer+1,Face)>nonVerticallyIntegratedCell(iRange)%minimumElevation) Then
                                                HydroParam%InFlowSmallm(HydroParam%NInflow) = iLayer !Vertical Integrated (review)
                                                Exit
                                            Else
                                                HydroParam%InFlowSmallm(HydroParam%NInflow) = HydroParam%Smallm(Face)
                                            EndIf
                                        EndDo
                                        Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face)
                                            If (HydroParam%Z(iLayer+1,Face)>=nonVerticallyIntegratedCell(iRange)%maximumElevation) Then
                                                HydroParam%InFlowCapitalM(HydroParam%NInflow) = iLayer !Vertical Integrated (review)
                                                Exit
                                            Else
                                                HydroParam%InFlowCapitalM(HydroParam%NInflow) = HydroParam%CapitalM(Face)
                                            EndIf
                                        EndDo                                    
                                        Exit
                                    EndIf
                                EndDo  
                                call c_f_pointer(nonVerticallyIntegratedCell(iRange)%timeSeriesList, hydrotimeSeries, [nonVerticallyIntegratedCell(iRange)%timeSeriesListSize])
                                Do j = 1, HydroParam%InFlownTime(HydroParam%NInflow)
                                    ! Verifing consistency of time series
                                    !If (j == 1) Then !First date of time series
                                    !    If (IniTime<hydrotimeSeries(j)%timeStampSize) Then
                                    !        Call unix2c(hydrotimeSeries(j)%timeStampSize, idate)
                                    !        Print*, 'First date of time series ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//')'
                                    !        Call unix2c(IniTime, idate)
                                    !        Print*, 'is greater than initial time of simulation ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//').' 
                                    !        Print*, 'Please check Boundary condition number ('//trim(str(i))//').'
                                    !        Pause
                                    !        Stop 
                                    !    EndIf
                                    !ElseIf (j == WaterLevelnTime(NWaterLevel)) Then
                                    !    If (FinalTime>hydrotimeSeries(j)%timeStampSize) Then
                                    !        Call unix2c(hydrotimeSeries(j)%timeStampSize, idate)
                                    !        Print*, 'Last date of time series ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//')'
                                    !        Call unix2c(FinalTime, idate)
                                    !        Print*, 'is lesser than initial time of simulation ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//').' 
                                    !        Print*, 'Please check Boundary condition number ('//trim(str(i))//').'
                                    !        Pause
                                    !        Stop 
                                    !    EndIf
                                    !EndIf
                                    ! Storing stamp time and values of time series 
                                    HydroParam%InFlowTime(HydroParam%NInflow,j) = hydrotimeSeries(j)%timeStamp
                                    HydroParam%InFlowValue(HydroParam%NInflow,j) = hydrotimeSeries(j)%value1 
                                EndDo
                                
                                
                            EndIf
                        EndDo
                    EndIf
                EndDo
            EndIf
 
        ElseIf (trim(boundaryConditions(i)%conditionType) == "normalDepth") Then !Normal Depth
            If (boundaryConditions(i)%conditionFunction == 1) Then !Constant
                Do j = 1,boundaryConditions(i)%cellsLength
                    HydroParam%NInflow = HydroParam%NInflow + 1
                    HydroParam%InFlownTime(HydroParam%NInflow) = -999
                    Do iEdge = 1,4
                        Face = MeshParam%Edge(iEdge,boundaryConditionCells(j)%cellId + 1)
                        iNode1 = MeshParam%Quadri(MeshParam%EdgeDef(1,iEdge),boundaryConditionCells(j)%cellId + 1)         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                        iNode2 = MeshParam%Quadri(MeshParam%EdgeDef(2,iEdge),boundaryConditionCells(j)%cellId + 1)         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                        If (iNode1==boundaryConditionCells(j)%verticeIds(1).and.iNode2==boundaryConditionCells(j)%verticeIds(2)) Then
                            HydroParam%IndexInflow(HydroParam%NInflow,1) = Face
                            HydroParam%IndexInflow(HydroParam%NInflow,2) = boundaryConditionCells(j)%cellId + 1
                            HydroParam%IndexInflowEdge(Face) = 1
                            HydroParam%InFlowSmallm(HydroParam%NInflow) = 1 !Vertical Integrated (review)
                            HydroParam%InFlowCapitalM(HydroParam%NInflow) = MeshParam%KMax !Vertical Integrated (review)
                            Exit
                        EndIf
                    EndDo
                    HydroParam%InFlowTime(HydroParam%NInflow,1) = -999
                    HydroParam%InFlowTime(HydroParam%NInflow,2) = -999
                    HydroParam%InFlowValue(HydroParam%NInflow,:) = boundaryConditions(i)%constantValue
                    
                EndDo
            EndIf
            
        EndIf
    EndDo    
    
    
    
    !Open(96,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='OLD',ACTION='READ')
    !Do i=1,HydroParam%InflownTime(1)
    !    READ(96,*) HydroParam%irrgMirim(i), HydroParam%irrgMangueira(i)    
    !EndDo
        
   
End Subroutine ReadHydroBoundaryCondition
    
    