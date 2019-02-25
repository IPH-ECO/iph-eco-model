!> This subroutine reads the hydrodynamic parameters. 
Subroutine ReadOutputs(sim,simParam,Kmax,nEdge,nElem)
    
    Use domain_types
    Use SimulationModel
    Use iso_c_binding
    !Use MeshVars, Only: Kmax,nEdge,nElem
    
    Implicit none
    
    Integer:: i,j,Kmax,nEdge,nElem
    type(Simulation) :: sim
    type(SimulationParam) :: simParam
    character, pointer :: ModelOutputPath(:)
    type(hydroOutputParameter), pointer :: hydroOutputParameters(:)
    type(wqOutputParameter), pointer :: wqOutputParameters(:)
    character(len=200):: text
    
    
    call c_f_pointer(sim%outputDirectory, ModelOutputPath, [sim%outputDirectoryLength])
    
    simParam%OutputPath =''
    Do j=1,sim%outputDirectoryLength
        simParam%OutputPath = simParam%OutputPath(1:j)//ModelOutputPath(j)
    EndDo
    simParam%OutputPath = trim(simParam%OutputPath)
    
    call c_f_pointer(sim%hydrooutputParameters, hydroOutputParameters, [sim%hydrooutputParametersLength])
    
    call c_f_pointer(sim%wqOutputParameters, wqOutputParameters, [sim%wqOutputParametersLength])
    
    simParam%OutputHydro = 0
    Do i = 1, sim%hydrooutputParametersLength
        
        text = trim(hydroOutputParameters(i)%name)

        !Hydrodynamic module
        If (trim(text) == 'VelocityField') Then
            simParam%OutputHydro(1) = 1
        ElseIf (trim(text) == 'SurfaceWaterElevation') Then
            simParam%OutputHydro(2) = 1
        ElseIf (trim(text) == 'depth') Then
            simParam%OutputHydro(3) = 1
        ElseIf (trim(text) == 'horizontalEddyViscosity') Then
            simParam%OutputHydro(4) = 1
        ElseIf (trim(text) == 'horizontalEddyDiffusivity') Then
            simParam%OutputHydro(5) = 1
        ElseIf (trim(text) == 'verticalEddyViscosity') Then
            simParam%OutputHydro(6) = 1
        ElseIf (trim(text) == 'verticalEddyDiffusivity') Then
            simParam%OutputHydro(7) = 1
        ElseIf (trim(text) == 'barotropicPressure') Then
            simParam%OutputHydro(8) = 1
        ElseIf (trim(text) == 'baroclinicPressure') Then
            simParam%OutputHydro(9) = 1
        ElseIf (trim(text) == 'nonHydrostaticPressureComponent') Then
            simParam%OutputHydro(10) = 1
        ElseIf (trim(text) == 'convectiveAcceleration') Then
            simParam%OutputHydro(11) = 1
        ElseIf (trim(text) == 'shearStressOnWaterSurface') Then
            simParam%OutputHydro(12) = 1
        ElseIf (trim(text) == 'shearStressOnSedimentBed') Then
            simParam%OutputHydro(13) = 1
        ElseIf (trim(text) == 'coriolisComponent') Then
            simParam%OutputHydro(14) = 1
            
            
        EndIf
    EndDo
    
    simParam%OutputWQ = 0
    
Do i = 1, sim%wqoutputParametersLength
        
        text = trim(wqOutputParameters(i)%name)

        !Hydrodynamic module
        If (trim(text) == 'sDTempW') Then
            simParam%OutputWQ(1) = 1
        ElseIf (trim(text) == 'sDSal') Then
            simParam%OutputWQ(2) = 1
        ElseIf (trim(text) == 'sO2W') Then
            simParam%OutputWQ(3) = 1
        ElseIf (trim(text) == 'spHW') Then
            simParam%OutputWQ(4) = 1
        ElseIf (trim(text) == 'spHS') Then
            simParam%OutputWQ(5) = 1
        ElseIf (trim(text) == 'sAlkW') Then
            simParam%OutputWQ(6) = 1
        ElseIf (trim(text) == 'sAlkS') Then
            simParam%OutputWQ(7) = 1
        ElseIf (trim(text) == 'sDicW') Then
            simParam%OutputWQ(8) = 1
        ElseIf (trim(text) == 'sCH4W') Then
            simParam%OutputWQ(9) = 1
        ElseIf (trim(text) == 'sCH4S') Then
            simParam%OutputWQ(10) = 1
        ElseIf (trim(text) == 'sPO4W') Then
            simParam%OutputWQ(11) = 1
        ElseIf (trim(text) == 'sPAIMW') Then
            simParam%OutputWQ(12) = 1
        ElseIf (trim(text) == 'sNO3W') Then
            simParam%OutputWQ(13) = 1
        ElseIf (trim(text) == 'sNO3S') Then
            simParam%OutputWQ(14) = 1
        ElseIf (trim(text) == 'sNH4W') Then
            simParam%OutputWQ(15) = 1
        ElseIf (trim(text) == 'sNH4S') Then
            simParam%OutputWQ(16) = 1
        ElseIf (trim(text) == 'sSiO2W') Then
            simParam%OutputWQ(17) = 1
        ElseIf (trim(text) == 'sDIMW') Then
            simParam%OutputWQ(18) = 1
        ElseIf (trim(text) == 'sOMW') Then
            simParam%OutputWQ(19) = 1
        ElseIf (trim(text) == 'sDPomW') Then
            simParam%OutputWQ(20) = 1
        ElseIf (trim(text) == 'sDPomS') Then
            simParam%OutputWQ(21) = 1
        ElseIf (trim(text) == 'sDDomW') Then
            simParam%OutputWQ(22) = 1
        ElseIf (trim(text) == 'sDDomS') Then
            simParam%OutputWQ(23) = 1
        ElseIf (trim(text) == 'sDPomLW') Then
            simParam%OutputWQ(24) = 1
        ElseIf (trim(text) == 'sDPomRW') Then
            simParam%OutputWQ(25) = 1
        ElseIf (trim(text) == 'sDPomLS') Then
            simParam%OutputWQ(26) = 1
        ElseIf (trim(text) == 'sDPomRS') Then
            simParam%OutputWQ(27) = 1
        ElseIf (trim(text) == 'sDDomLW') Then
            simParam%OutputWQ(28) = 1
        ElseIf (trim(text) == 'sDDomRW') Then
            simParam%OutputWQ(29) = 1
        ElseIf (trim(text) == 'sDDomLS') Then
            simParam%OutputWQ(30) = 1
        ElseIf (trim(text) == 'sDDomRS') Then
            simParam%OutputWQ(31) = 1
        ElseIf (trim(text) == 'sDBacW') Then
            simParam%OutputWQ(32) = 1
        ElseIf (trim(text) == 'sDBacS') Then
            simParam%OutputWQ(33) = 1
        ElseIf (trim(text) == 'sDPhytW') Then
            simParam%OutputWQ(34) = 1
        ElseIf (trim(text) == 'sDPhytS') Then
            simParam%OutputWQ(35) = 1
        ElseIf (trim(text) == 'sDMac') Then
            simParam%OutputWQ(36) = 1
        ElseIf (trim(text) == 'sDZoo') Then
            simParam%OutputWQ(37) = 1
        ElseIf (trim(text) == 'sDBent') Then
            simParam%OutputWQ(38) = 1
        ElseIf (trim(text) == 'sDFiAd') Then
            simParam%OutputWQ(39) = 1
        ElseIf (trim(text) == 'sDFiJv') Then
            simParam%OutputWQ(40) = 1
        EndIf
    EndDo
    Return
   
End Subroutine ReadOutputs
    