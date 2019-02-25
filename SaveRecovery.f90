Subroutine SaveRecovery(sim,simParam,MeshParam,HydroParam,LimnoParam)
    !DIR$ NOOPTIMIZE
    Use domain_types
    Use SimulationModel
    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars
    
    type(Simulation) :: sim
    type(SimulationParam) :: simParam
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    
    character(len=200):: text
    
    If (simParam%recoveryTimeStep /= 0) Then
        If (mod(simParam%it,simParam%recoveryTimeStep)==0.or.simParam%it==1) Then
        !If (mod(simParam%time-simParam%IniTime,int(simParam%recoveryTimeStep*simParam%dt))==0.) Then
            simParam%SaveVariables%layers = MeshParam%KMax
            simParam%SaveVariables%edges = MeshParam%nEdge
            simParam%SaveVariables%elements = MeshParam%nElem
            simParam%SaveVariables%simulationTime = simParam%time
            simParam%SaveVariables%u = c_loc(HydroParam%u)
            simParam%SaveVariables%w = c_loc(HydroParam%w)
            simParam%SaveVariables%eta = c_loc(HydroParam%eta)
                
            !simParam%SaveVariables%wqo = c_loc(LimnoParam%sDSalW)
            Do iOutput = 1,sim%wqoOutputParametersLength
                text = trim(simParam%wqoOutputParameters(iOutput)%name)
                If (trim(text)=='sDTempW') Then
                    Do iElem = 1,MeshParam%nElem
                        Do iLayer = 1,MeshParam%KMax
                            simParam%wqoutput(iLayer + (iElem-1)*MeshParam%KMax +(iOutput-1)*(MeshParam%nElem*MeshParam%KMax)) = LimnoParam%sDTempW(iLayer,iElem) 
                        EndDo
                    EndDo
                ElseIf (trim(text)=='sDSal') Then
                    Do iElem = 1,MeshParam%nElem
                        Do iLayer = 1,MeshParam%KMax
                            simParam%wqoutput(iLayer + (iElem-1)*MeshParam%KMax +(iOutput-1)*(MeshParam%nElem*MeshParam%KMax)) = LimnoParam%sDSal(iLayer,iElem) 
                        EndDo
                    EndDo
                EndIf
            EndDo
                
            simParam%SaveVariables%wqo = c_loc(simParam%wqoutput)
                
            sim%recoveryVariables = c_loc(simParam%SaveVariables)
                
            simParam%SaveVariables%changed = logical(.true., kind = c_bool)
                
            do while (simParam%SaveVariables%changed==.true.) 
                continue
            end do
        EndIf
    EndIf
    
    
    
End Subroutine SaveRecovery