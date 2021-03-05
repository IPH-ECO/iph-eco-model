Subroutine VTKOutput(simParam,HydroParam,MeshParam,LimnoParam)
    
    ! Write Results in VisIt Format (.vtk)
    ! Based on:
    ! [1] Childs abd Others. VisIt: An End-User Tool for Visualizing and Analysing Very Large Data.
    !     In: High Performance Visualization -- Enabling Extreme-Scale Scientific Insight. Oct/2012. Pages: 357-372
    ! [2] VTK File Format. Taken from The VTK User's Guide. Kitware: www.kitware.com
    
    ! Input:
    ! it     -> Time Iteration
    ! time -> Time Step
    ! Output:
    ! ".vtk" File with Mesh and Variables Information
    
    ! List of Modifications: 
    !   -> 07.11.2014: Routine Implementation (Rafael Cavalcanti)
    !   -> 03.12.2014: 3D Extension           (Rafael Cavalcanti)
    ! Programmer: Rafael Cavalcanti
    Use SimulationModel
    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars
    Use LIB_VTK_IO
    
    Implicit None
    Integer:: iNode, iElem, iLayer
    Integer:: k, Sum, bb, gg, mm, zz, ben
    Integer:: TotNumberOfPoints, nElem3D, satLayers
    Real:: Vel(2), lw(2),DIR_VENTO,FF
    Real:: V, NearZero
    Real:: idate(6)
    Character(200):: FileName
    Character(200):: VarName
    Character(200):: OutputTime
    type(SimulationParam) :: simParam
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    
  
    integer(I4P):: E_IO ! Input/Output inquiring flag: 0 if IO is done, > 0 if IO is not done
                        !        Call unix2c(hydrotimeSeries(j)%timeStampSize, idate)
                        !        Print*, 'First date of time series ('//trim(str(idate(1)))//'/'//trim(str(idate(2)))//'/'//trim(str(idate(3)))//' '//trim(str(idate(4)))//':'//trim(str(idate(5)))//':'//trim(str(idate(6)))//')'

   simParam%it_vtk = simParam%it_vtk + 1 
   ! 1. Define File Headers
    Write(FileName,'(a,a,i8.8,a)') trim(simParam%OutputFile),'_',simParam%it_vtk
    Call unix2c(simParam%time, idate)
    Write(OutputTime,'(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,f13.10)') floor(idate(1)),'/',floor(idate(2)),'/',floor(idate(3)),' ',floor(idate(4)),':',floor(idate(5)),':',idate(6)
    
    E_IO = VTK_INI('ascii', &
                    trim(simParam%OutputPath)//'\'//trim(filename)//'.vtk', &
                    trim(OutputTime), &
                    'UNSTRUCTURED_GRID')
    
    NearZero = 0.001
    
    ! 2. Define Geometry
    ! VTK coordinates
    Do iElem = 1, MeshParam%nElem
        Do iLayer = 1,MeshParam%KMax
            MeshParam%xPoint(MeshParam%Quadri(1,iElem)+(iLayer-1)*MeshParam%nPoint+1) = MeshParam%xNode( MeshParam%Quadri(1,iElem)+1 )
            MeshParam%xPoint(MeshParam%Quadri(2,iElem)+(iLayer-1)*MeshParam%nPoint+1) = MeshParam%xNode( MeshParam%Quadri(2,iElem)+1 )
            MeshParam%xPoint(MeshParam%Quadri(3,iElem)+(iLayer-1)*MeshParam%nPoint+1) = MeshParam%xNode( MeshParam%Quadri(3,iElem)+1 )
            MeshParam%xPoint(MeshParam%Quadri(4,iElem)+(iLayer-1)*MeshParam%nPoint+1) = MeshParam%xNode( MeshParam%Quadri(4,iElem)+1 )
            MeshParam%xPoint(MeshParam%Quadri(1,iElem)+iLayer*MeshParam%nPoint+1) = MeshParam%xNode( MeshParam%Quadri(1,iElem)+1 )
            MeshParam%xPoint(MeshParam%Quadri(2,iElem)+iLayer*MeshParam%nPoint+1) = MeshParam%xNode( MeshParam%Quadri(2,iElem)+1 )
            MeshParam%xPoint(MeshParam%Quadri(3,iElem)+iLayer*MeshParam%nPoint+1) = MeshParam%xNode( MeshParam%Quadri(3,iElem)+1 )
            MeshParam%xPoint(MeshParam%Quadri(4,iElem)+iLayer*MeshParam%nPoint+1) = MeshParam%xNode( MeshParam%Quadri(4,iElem)+1 )
            MeshParam%yPoint(MeshParam%Quadri(1,iElem)+(iLayer-1)*MeshParam%nPoint+1) = MeshParam%yNode( MeshParam%Quadri(1,iElem)+1 )
            MeshParam%yPoint(MeshParam%Quadri(2,iElem)+1) = MeshParam%yNode( MeshParam%Quadri(2,iElem)+1 )
            MeshParam%yPoint(MeshParam%Quadri(3,iElem)+(iLayer-1)*MeshParam%nPoint+1) = MeshParam%yNode( MeshParam%Quadri(3,iElem)+1 )
            MeshParam%yPoint(MeshParam%Quadri(4,iElem)+(iLayer-1)*MeshParam%nPoint+1) = MeshParam%yNode( MeshParam%Quadri(4,iElem)+1 )
            MeshParam%yPoint(MeshParam%Quadri(1,iElem)+iLayer*MeshParam%nPoint+1) = MeshParam%yNode( MeshParam%Quadri(1,iElem)+1 )
            MeshParam%yPoint(MeshParam%Quadri(2,iElem)+iLayer*MeshParam%nPoint+1) = MeshParam%yNode( MeshParam%Quadri(2,iElem)+1 )
            MeshParam%yPoint(MeshParam%Quadri(3,iElem)+iLayer*MeshParam%nPoint+1) = MeshParam%yNode( MeshParam%Quadri(3,iElem)+1 )
            MeshParam%yPoint(MeshParam%Quadri(4,iElem)+iLayer*MeshParam%nPoint+1) = MeshParam%yNode( MeshParam%Quadri(4,iElem)+1 )
            If (iLayer == 1) Then
                MeshParam%zPoint(MeshParam%Quadri(1,iElem)+(iLayer-1)*MeshParam%nPoint+1) = HydroParam%Ze(MeshParam%Kmax-iLayer+2,iElem) !zR
                MeshParam%zPoint(MeshParam%Quadri(2,iElem)+(iLayer-1)*MeshParam%nPoint+1) = HydroParam%Ze(MeshParam%Kmax-iLayer+2,iElem) !zR
                MeshParam%zPoint(MeshParam%Quadri(3,iElem)+(iLayer-1)*MeshParam%nPoint+1) = HydroParam%Ze(MeshParam%Kmax-iLayer+2,iElem) !zR
                MeshParam%zPoint(MeshParam%Quadri(4,iElem)+(iLayer-1)*MeshParam%nPoint+1) = HydroParam%Ze(MeshParam%Kmax-iLayer+2,iElem) !zR
            Else
                MeshParam%zPoint(MeshParam%Quadri(1,iElem)+(iLayer-1)*MeshParam%nPoint+1) = HydroParam%Ze(MeshParam%Kmax-iLayer+2,iElem) !LIMCAMAUX(iLayer-1)
                MeshParam%zPoint(MeshParam%Quadri(2,iElem)+(iLayer-1)*MeshParam%nPoint+1) = HydroParam%Ze(MeshParam%Kmax-iLayer+2,iElem) !LIMCAMAUX(iLayer-1)
                MeshParam%zPoint(MeshParam%Quadri(3,iElem)+(iLayer-1)*MeshParam%nPoint+1) = HydroParam%Ze(MeshParam%Kmax-iLayer+2,iElem) !LIMCAMAUX(iLayer-1)
                MeshParam%zPoint(MeshParam%Quadri(4,iElem)+(iLayer-1)*MeshParam%nPoint+1) = HydroParam%Ze(MeshParam%Kmax-iLayer+2,iElem) !LIMCAMAUX(iLayer-1)
            EndIf
            MeshParam%zPoint(MeshParam%Quadri(1,iElem)+iLayer*MeshParam%nPoint+1) = HydroParam%Ze(MeshParam%Kmax-iLayer+1,iElem) !LIMCAMAUX(iLayer)
            MeshParam%zPoint(MeshParam%Quadri(2,iElem)+iLayer*MeshParam%nPoint+1) = HydroParam%Ze(MeshParam%Kmax-iLayer+1,iElem) !LIMCAMAUX(iLayer)
            MeshParam%zPoint(MeshParam%Quadri(3,iElem)+iLayer*MeshParam%nPoint+1) = HydroParam%Ze(MeshParam%Kmax-iLayer+1,iElem) !LIMCAMAUX(iLayer)
            MeshParam%zPoint(MeshParam%Quadri(4,iElem)+iLayer*MeshParam%nPoint+1) = HydroParam%Ze(MeshParam%Kmax-iLayer+1,iElem) !LIMCAMAUX(iLayer)
            
            If (iLayer == MeshParam%KMax .and. HydroParam%Ze(MeshParam%Kmax-iLayer+2,iElem) > HydroParam%eta(iElem)) Then
                MeshParam%zPoint(MeshParam%Quadri(1,iElem)+(iLayer-1)*MeshParam%nPoint+1) = Max(HydroParam%eta(iElem), NearZero)  !HydroParam%eta(iElem) !LIMCAMAUX(iLayer)
                MeshParam%zPoint(MeshParam%Quadri(2,iElem)+(iLayer-1)*MeshParam%nPoint+1) = Max(HydroParam%eta(iElem), NearZero) !HydroParam%eta(iElem) !LIMCAMAUX(iLayer)
                MeshParam%zPoint(MeshParam%Quadri(3,iElem)+(iLayer-1)*MeshParam%nPoint+1) = Max(HydroParam%eta(iElem), NearZero) !HydroParam%eta(iElem) !LIMCAMAUX(iLayer)
                MeshParam%zPoint(MeshParam%Quadri(4,iElem)+(iLayer-1)*MeshParam%nPoint+1) = Max(HydroParam%eta(iElem), NearZero) !HydroParam%eta(iElem) !LIMCAMAUX(iLayer)               
            EndIf
            
        EndDo
    EndDo
    
    E_IO = VTK_GEO(MeshParam%nPoint*(MeshParam%KMax+1),MeshParam%xPoint,MeshParam%yPoint,MeshParam%zPoint)
    
    ! 3. Define Connectivity

    ! VTK connectivity between cells
    
    Do iElem = 1,MeshParam%nElem
        Do iLayer = 1,MeshParam%KMax
            MeshParam%Connect((iElem-1)*9+(MeshParam%KMax-iLayer)*MeshParam%nElem*9+1) = 8
            MeshParam%Connect((iElem-1)*9+(MeshParam%KMax-iLayer)*MeshParam%nElem*9+2) = MeshParam%Quadri(1,iElem) + iLayer*MeshParam%nPoint  
            MeshParam%Connect((iElem-1)*9+(MeshParam%KMax-iLayer)*MeshParam%nElem*9+3) = MeshParam%Quadri(2,iElem) + iLayer*MeshParam%nPoint  
            MeshParam%Connect((iElem-1)*9+(MeshParam%KMax-iLayer)*MeshParam%nElem*9+4) = MeshParam%Quadri(3,iElem) + iLayer*MeshParam%nPoint  
            MeshParam%Connect((iElem-1)*9+(MeshParam%KMax-iLayer)*MeshParam%nElem*9+5) = MeshParam%Quadri(4,iElem) + iLayer*MeshParam%nPoint  
            MeshParam%Connect((iElem-1)*9+(MeshParam%KMax-iLayer)*MeshParam%nElem*9+6) = MeshParam%Quadri(1,iElem) + (iLayer-1)*MeshParam%nPoint  
            MeshParam%Connect((iElem-1)*9+(MeshParam%KMax-iLayer)*MeshParam%nElem*9+7) = MeshParam%Quadri(2,iElem) + (iLayer-1)*MeshParam%nPoint 
            MeshParam%Connect((iElem-1)*9+(MeshParam%KMax-iLayer)*MeshParam%nElem*9+8) = MeshParam%Quadri(3,iElem) + (iLayer-1)*MeshParam%nPoint 
            MeshParam%Connect((iElem-1)*9+(MeshParam%KMax-iLayer)*MeshParam%nElem*9+9) = MeshParam%Quadri(4,iElem) + (iLayer-1)*MeshParam%nPoint 
        EndDo    
    EndDo    
    MeshParam%cell_type = 12
    E_IO = VTK_CON(MeshParam%nElem*MeshParam%KMax,MeshParam%connect,MeshParam%cell_type)  !
    
    ! 4. Initializes the saving of data variables
    E_IO = VTK_DAT(MeshParam%nElem*MeshParam%KMax,'cell')
    
    ! 5. Saves the data variables related to geometric mesh
    ! 5.1 Velocity Vector
    If (simParam%OutputHydro(1)==1) Then
        !Subsurface case:
        If (MeshParam%iBedrock == 1) Then
            Do iElem = 1,MeshParam%nElem
                Do iLayer = 1,MeshParam%KMax
                    HydroParam%SVector(iElem + (iLayer-1)*MeshParam%nElem,1) = HydroParam%ubsub(iLayer,1,iElem)
                    HydroParam%SVector(iElem + (iLayer-1)*MeshParam%nElem,2) = HydroParam%ubsub(iLayer,2,iElem)
                    HydroParam%SVector(iElem + (iLayer-1)*MeshParam%nElem,3) = HydroParam%ubsub(iLayer,3,iElem)
                EndDo
            EndDo
        Else
        !Superficial case:
           Do iElem = 1,MeshParam%nElem
                Do iLayer = 1,MeshParam%KMax
                    HydroParam%SVector(iElem + (iLayer-1)*MeshParam%nElem,1) = HydroParam%ub(iLayer,1,iElem)
                    HydroParam%SVector(iElem + (iLayer-1)*MeshParam%nElem,2) = HydroParam%ub(iLayer,2,iElem)
                    HydroParam%SVector(iElem + (iLayer-1)*MeshParam%nElem,3) = HydroParam%ub(iLayer,3,iElem)
                EndDo
            EndDo          
        EndIf
        E_IO = VTK_VAR('vect',MeshParam%nElem*MeshParam%KMax,'VelocityField',HydroParam%SVector(:,1),HydroParam%SVector(:,2),HydroParam%SVector(:,3))
    EndIf
    ! 5.2 Surface Water Elevation
    If (simParam%OutputHydro(2)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = Max(HydroParam%eta(iElem), NearZero) 
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'SurfaceWaterElevation',HydroParam%SScalar)
    EndIf
    
    ! 5.3 Depth
    If (simParam%OutputHydro(3)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                !HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = V(HydroParam%eta(iElem),HydroParam%hb(iElem)) 
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = Max( V(HydroParam%eta(iElem),HydroParam%sb(iElem)), NearZero)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Depth',HydroParam%SScalar)
    EndIf

    !! x.x Saturation
    !Do iElem = 1,MeshParam%nElem
    !    If (MeshParam%Kmax > 1) Then
    !        satLayers = MeshParam%Kmax
    !        Do iLayer = 1,MeshParam%KMax
    !            HydroParam%SScalarSaturation(iElem + (iLayer-1)*MeshParam%nElem) =  MeshParam%Si(iLayer,iElem)
    !        EndDo
    !    Else
    !        satLayers = MeshParam%subfactor
    !        Do iLayer = 1,MeshParam%subfactor
    !            HydroParam%SScalarSaturation(iElem + (iLayer-1)*MeshParam%nElem) =  MeshParam%Si(iLayer,iElem)
    !        EndDo            
    !    EndIf
    !EndDo
    !E_IO = VTK_VAR(MeshParam%nElem*satLayers,'Saturation',HydroParam%SScalarSaturation)
    
    ! 5.4 NonHydrostatic Pressure
    Do iElem = 1,MeshParam%nElem
        
        Do iLayer = 1,MeshParam%KMax
            HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = HydroParam%q(iLayer,iElem) 
        EndDo
    EndDo
    E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'NonHydrostaticPressure',HydroParam%SScalar)
    
    
    ! 5.2 Salinity
    !Do iElem = 1,nElem
    !    !Calculate lake fetch
    !    !DIR_VENTO = ATAN2(Windiy(iElem),Windix(iElem))-PI-3.*PI/4.
    !    DIR_VENTO = ATAN2(1.,1.)-PI-3.*PI/4.
    !    IF (DIR_VENTO < -2.*PI) THEN
    !        DIR_VENTO = DIR_VENTO + 2.*PI
    !    ELSEIF (DIR_VENTO == -2.*PI) THEN
    !        DIR_VENTO = 2.*PI
    !    ENDIF 
    !    DIR_VENTO = NINT(ABS(DIR_VENTO/(PI/4.))) !índice da matrix fetch (de acordo com a direção do vento)
    !    DIR_VENTO = MIN(DIR_VENTO,8.)
    !    DIR_VENTO = MAX(DIR_VENTO,1.)
    !    FF = fetch_m(iElem,int(DIR_VENTO)) !fetch (m)
    !    Do iLayer = 1,KMax
    !        SScalar(iElem + (iLayer-1)*nElem) = FF 
    !    EndDo
    !EndDo
    !E_IO = VTK_VAR(nElem*KMax,'Salinity',SScalar)
    
     ! 5.5 Water temperature
    If (simParam%OutputWQ(1)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDTempW(iLayer,iElem) 
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'WaterTemp',HydroParam%SScalar)
    EndIf
    
    ! 5.6 Salinity
    !If (simParam%OutputWQ(2)==1) Then
        Do iElem = 1,MeshParam%nElem
        
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = HydroParam%sDRhoW(iLayer,iElem)!LimnoParam%sDSal(iLayer,iElem) !HydroParam%sDRhoW(iLayer,iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Salinity',HydroParam%SScalar)
    !EndIf
    
    Do iElem = 1,MeshParam%nElem
        Do iLayer = 1,MeshParam%KMax
            HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) =  HydroParam%psi_edge(iLayer, MeshParam%Edge(2,iElem))!LimnoParam%sDSal(iLayer,iElem) !HydroParam%sDRhoW(iLayer,iElem)
        EndDo
    EndDo
    E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Psi_face',HydroParam%SScalar)

    Do iElem = 1,MeshParam%nElem
        Do iLayer = 1,MeshParam%KMax
            HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) =  HydroParam%psi_cell(iLayer, iElem)!LimnoParam%sDSal(iLayer,iElem) !HydroParam%sDRhoW(iLayer,iElem)
        EndDo
    EndDo
    E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Psi_cell',HydroParam%SScalar)

    ! 5.7 DO_in_Water
    If (simParam%OutputWQ(3)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sO2W(iLayer,iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'DO_in_Water',HydroParam%SScalar)
    EndIf
  
    ! 5.8 pH_in_Water
    If (simParam%OutputWQ(4)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%spHW(iLayer,iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'pH_in_Water',HydroParam%SScalar)
    EndIf
    
    ! 5.9 pH_in_Sediment
    If (simParam%OutputWQ(5)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%spHS(iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'pH_in_Sediment',HydroParam%SScalar)
    EndIf
    
    ! 5.10 Alkalinity_in_Water
    If (simParam%OutputWQ(6)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sAlkW(iLayer,iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Alkalinity_in_Water',HydroParam%SScalar)
    EndIf      
    
    ! 5.11 Alkalinity_in_Sediment
    If (simParam%OutputWQ(7)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sAlkS(iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Alkalinity_in_Sediment',HydroParam%SScalar)
    EndIf   
    
    ! 5.12 DIC_in_Water
    If (simParam%OutputWQ(8)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDicW(iLayer,iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'DIC_in_Water',HydroParam%SScalar)
    EndIf        
    
    ! 5.13 CH4_in_Water
    If (simParam%OutputWQ(9)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sCH4W(iLayer,iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'CH4_in_Water',HydroParam%SScalar)
    EndIf      
    
    ! 5.14 CH4_in_Sediment
    If (simParam%OutputWQ(10)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sCH4S(iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'CH4_in_Sediment',HydroParam%SScalar)
    EndIf      
    
    ! 5.15 PO4_in_Water
    If (simParam%OutputWQ(11)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sPO4W(iLayer,iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'PO4_in_Water',HydroParam%SScalar)
    EndIf      

    ! 5.16 Adsoved_P_in_Water
    If (simParam%OutputWQ(12)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sPAIMW(iLayer,iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Adsoved_P_in_Water',HydroParam%SScalar)
    EndIf      
    
    ! 5.17 NO3_in_Water
    If (simParam%OutputWQ(13)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sNO3W(iLayer,iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'NO3_in_Water',HydroParam%SScalar)
    EndIf   
    
    ! 5.18 NO3_in_Sediment
    If (simParam%OutputWQ(14)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sNO3S(iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'NO3_in_Sediment',HydroParam%SScalar)
    EndIf     
    
    ! 5.19 NH4_in_Water
    If (simParam%OutputWQ(15)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sNH4W(iLayer,iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'NH4_in_Water',HydroParam%SScalar)
    EndIf   
    
    ! 5.20 NH4_in_Sediment
    If (simParam%OutputWQ(16)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sNH4S(iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'NH4_in_Sediment',HydroParam%SScalar)
    EndIf  
    
    ! 5.21 SiO2_in_Water
    If (simParam%OutputWQ(17)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sSiO2W(iLayer,iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'SiO2_in_Water',HydroParam%SScalar)
    EndIf   
    
    ! 5.22 Inorganic_Matter_in_Water
    If (simParam%OutputWQ(18)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDIMW(iLayer,iElem)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Inorganic_Matter_in_Water',HydroParam%SScalar)
    EndIf   
    
    ! 5.23 Organic_Matter_in_Water
    If (simParam%OutputWQ(19)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDPomW(iLayer,iElem,1)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Organic_Matter_in_Water',HydroParam%SScalar)
    EndIf      
    
    ! 5.24 Particulated_Organic_Matter_in_Water
    If (simParam%OutputWQ(20)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDPomW(iLayer,iElem,1)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Particulated_Organic_Matter_in_Water',HydroParam%SScalar)
    EndIf  
    
    ! 5.25 Particulated_Organic_Matter_in_Sediment
    If (simParam%OutputWQ(21)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDPomS(iElem,1)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Particulated_Organic_Matter_in_Sediment',HydroParam%SScalar)
    EndIf  
    
    ! 5.26 Dissolved_Organic_Matter_in_Water
    If (simParam%OutputWQ(22)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDDomW(iLayer,iElem,1)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Dissolved_Organic_Matter_in_Water',HydroParam%SScalar)
    EndIf  
    
    ! 5.27 Dissolved_Organic_Matter_in_Sediment
    If (simParam%OutputWQ(23)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDDomS(iElem,1)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Dissolved_Organic_Matter_in_Sediment',HydroParam%SScalar)
    EndIf     
    
    ! 5.28 Labile_Particulated_Organic_Matter_in_Water
    If (simParam%OutputWQ(24)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDPomW(iLayer,iElem,1)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Labile_Particulated_Organic_Matter_in_Water',HydroParam%SScalar)
    EndIf      
    
    ! 5.29 Refractory_Particulated_Organic_Matter_in_Water
    If (simParam%OutputWQ(25)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDPomW(iLayer,iElem,2)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Refractory_Particulated_Organic_Matter_in_Water',HydroParam%SScalar)
    EndIf      
    
    ! 5.30 Labile_Particulated_Organic_Matter_in_Sediment
    If (simParam%OutputWQ(26)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDPomS(iElem,1)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Labile_Particulated_Organic_Matter_in_Sediment',HydroParam%SScalar)
    EndIf      
    
    ! 5.31 Refractory_Particulated_Organic_Matter_in_Sediment
    If (simParam%OutputWQ(27)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDPomS(iElem,2)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Refractory_Particulated_Organic_Matter_in_Sediment',HydroParam%SScalar)
    EndIf      
    
    ! 5.32 Labile_Dissolved_Organic_Matter_in_Water
    If (simParam%OutputWQ(28)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDDomW(iLayer,iElem,1)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Labile_Dissolved_Organic_Matter_in_Water',HydroParam%SScalar)
    EndIf      
    
    ! 5.33 Refractory_Dissolved__Organic_Matter_in_Water
    If (simParam%OutputWQ(29)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDDomW(iLayer,iElem,2)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Refractory_Dissolved__Organic_Matter_in_Water',HydroParam%SScalar)
    EndIf      
    
    ! 5.34 Labile_Dissolved__Organic_Matter_in_Sediment
    If (simParam%OutputWQ(30)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDDomS(iElem,1)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Labile_Dissolved__Organic_Matter_in_Sediment',HydroParam%SScalar)
    EndIf      
    
    ! 5.35 Refractory_Dissolved__Organic_Matter_in_Sediment
    If (simParam%OutputWQ(31)==1) Then
        Do iElem = 1,MeshParam%nElem
            Do iLayer = 1,MeshParam%KMax
                HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDDomS(iElem,2)
            EndDo
        EndDo
        E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,'Refractory_Dissolved__Organic_Matter_in_Sediment',HydroParam%SScalar)
    EndIf       
    
    ! 5.36 Bacterioplankton_in_Water
    If (simParam%OutputWQ(32)==1) Then
        Do bb = 1, LimnoParam%nbac
            Write(VarName,'(a,a,i8.8)') 'Bacterioplankton_in_Water','_Group',bb
            Do iElem = 1,MeshParam%nElem
                Do iLayer = 1,MeshParam%KMax
                    HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDBacW(iLayer,iElem,bb)
                EndDo
            EndDo
            E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,VarName,HydroParam%SScalar)
        EndDo
    EndIf       
    
    ! 5.37 Bacterioplankton_in_Sediment
    If (simParam%OutputWQ(33)==1) Then
        Do bb = 1, LimnoParam%nbac
            Write(VarName,'(a,a,i8.8)') 'Bacterioplankton_in_Sediment','_Group',bb
            Do iElem = 1,MeshParam%nElem
                Do iLayer = 1,MeshParam%KMax
                    HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDBacS(iElem,bb)
                EndDo
            EndDo
            E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,VarName,HydroParam%SScalar)
        EndDo
    EndIf  
    
    ! 5.38 Phytoplankton_in_Water
    If (simParam%OutputWQ(34)==1) Then
        Do gg = 1, LimnoParam%nphy
            Write(VarName,'(a,a,i8.8)') 'Phytoplankton_in_Water','_Group',gg
            Do iElem = 1,MeshParam%nElem
                Do iLayer = 1,MeshParam%KMax
                    HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDPhytW(iLayer,iElem,gg) 
                EndDo
            EndDo
            E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,VarName,HydroParam%SScalar)
        EndDo
    EndIf

    ! 5.39 Phytoplankton_in_Sediment
    If (simParam%OutputWQ(35)==1) Then
        Do gg = 1, LimnoParam%nphy
            Write(VarName,'(a,a,i8.8)') 'Phytoplankton_in_Sediment','_Group',gg
            Do iElem = 1,MeshParam%nElem
                Do iLayer = 1,MeshParam%KMax
                    HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDPhytS(iElem,gg) 
                EndDo
            EndDo
            E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,VarName,HydroParam%SScalar)
        EndDo
    EndIf    

    ! 5.40 Macrophytes
    If (simParam%OutputWQ(36)==1) Then
        Do mm = 1, LimnoParam%nmac
            Write(VarName,'(a,a,i8.8)') 'Macrophytes','_Group',mm
            Do iElem = 1,MeshParam%nElem
                Do iLayer = 1,MeshParam%KMax
                    HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDMac(iElem,mm) 
                EndDo
            EndDo
            E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,VarName,HydroParam%SScalar)
        EndDo
    EndIf      
    
    ! 5.41 Zooplankton
    If (simParam%OutputWQ(37)==1) Then
        Do zz = 1, LimnoParam%nzoo
            Write(VarName,'(a,a,i8.8)') 'Zooplankton','_Group',zz
            Do iElem = 1,MeshParam%nElem
                Do iLayer = 1,MeshParam%KMax
                    HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDZoo(iLayer,iElem,zz) 
                EndDo
            EndDo
            E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,VarName,HydroParam%SScalar)
        EndDo
    EndIf    
    
    ! 5.42 Zoobenthos
    If (simParam%OutputWQ(38)==1) Then
        Do ben = 1, LimnoParam%nben
            Write(VarName,'(a,a,i8.8)') 'Zoobenthos','_Group',zz
            Do iElem = 1,MeshParam%nElem
                Do iLayer = 1,MeshParam%KMax
                    HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDBent(iElem,ben) 
                EndDo
            EndDo
            E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,VarName,HydroParam%SScalar)
        EndDo
    EndIf    
    
    ! 5.43 Adult_Fish
    If (simParam%OutputWQ(39)==1) Then
        Do ff = 1, LimnoParam%nfish
            Write(VarName,'(a,a,i8.8)') 'Adult_Fish','_Group',zz
            Do iElem = 1,MeshParam%nElem
                Do iLayer = 1,MeshParam%KMax
                    HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDFiAd(iLayer,iElem,ff) 
                EndDo
            EndDo
            E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,VarName,HydroParam%SScalar)
        EndDo
    EndIf    

    ! 5.43 Juvenile_Fish
    If (simParam%OutputWQ(40)==1) Then
        Do ff = 1, LimnoParam%nfish
            Write(VarName,'(a,a,i8.8)') 'Juvenile_Fish','_Group',zz
            Do iElem = 1,MeshParam%nElem
                Do iLayer = 1,MeshParam%KMax
                    HydroParam%SScalar(iElem + (iLayer-1)*MeshParam%nElem) = LimnoParam%sDFiJv(iLayer,iElem,ff) 
                EndDo
            EndDo
            E_IO = VTK_VAR(MeshParam%nElem*MeshParam%KMax,VarName,HydroParam%SScalar)
        EndDo
    EndIf    
    
    
    ! 6. finalize the file opened
    E_IO = VTK_END()
    
    Return
End Subroutine VTKOutput