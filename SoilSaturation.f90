  Subroutine SoilSaturation(iLayer,iElem,HydroParam,MeshParam)
    
    Use MeshVars 
    Use Hydrodynamic
    Implicit none
    Integer :: nLayers, iElem, iLayer
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    
    Do iElem = 1,MeshParam%nElem
        MeshParam%Si(:,iElem) = 1.d0
        MeshParam%Ki(:,iElem) = MeshParam%Ksat(:,iElem)
        If(HydroParam%ElSmallms(iElem) == HydroParam%ElCapitalMs(iElem)) Then 
            continue
        ElseIf (HydroParam%eta(iElem) <= sum(HydroParam%DZsi(:HydroParam%CapitalMs))) Then
            nLayers = Max(Floor(HydroParam%ElCapitalMs(iElem)-HydroParam%ElSmallms(iElem))/2,1)
            If (HydroParam%eta(iElem) > sum(HydroParam%DZsi(:nLayers))) Then
                nLayers = nLayers + 1
                Do iLayer = nLayers,HydroParam%ElCapitalMs(iElem)
                    Call Saturation(iLayer, iElem, HydroParam, MeshParam) 
                EndDo
            Else
                Do iLayer = HydroParam%ElSmallms(iElem),nLayers
                    Call Saturation(iLayer, iElem, HydroParam, MeshParam) 
                EndDo
            EndIf
        EndIf

    EndDo
    
    End Subroutine SoilSaturation
    
    
    Subroutine Saturation(iLayer, iElem, HydroParam, MeshParam)
                 
    Use MeshVars 
    Use Hydrodynamic
            
    Implicit none
    
    Real :: SoilSat, mParam
    Integer :: nnel, jlev, id0, FuFw_flag, Interpolate_Flag, i
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam

    If(Saturation_flag == 0) Then !Darcy's Model
        MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)*MeshParam%Si(iLayer,iElem) 
    ElseIf(Saturation_flag == 1) Then ! Brooks and Corey's Model
        If(HydroParam%eta(iElem) < HydroParam%Zb(iLayer,iElem) - 1/MeshParam%alpha(iLayer,iElem)) Then
            MeshParam%Si(iLayer,iElem) = MeshParam%alpha(iLayer,iElem)*(HydroParam%Zb(iLayer,iElem) - HydroParam%eta(iElem))**(MeshParam%nSoil(iLayer,iElem))
            MeshParam%Si(iLayer,iElem) = 1/MeshParam%Si(iLayer,iElem)
        EndIf      
        MeshParam%Ki(iLayer,iElem) = MeshParam%Ksati(iLayer,iElem)*MeshParam%Si(iLayer,iElem)**(3+2/MeshParam%nSoil(iLayer,iElem))
    Else !van Genuchten's Model
        mParam = (1-1/MeshParam%nSoil(iLayer,iElem))
        If(HydroParam%eta(iElem) < HydroParam%Zb(iLayer,iElem)) Then
            MeshParam%Si(iLayer,iElem) = 1/(1 + MeshParam%alpha(iLayer,iElem)*(HydroParam%Zb(iLayer,iElem) - HydroParam%eta(iElem)))**mParam
        EndIf
        MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)*(1-(1-MeshParam%Si(iLayer,iElem)**(1/mParam))**mParam)**2
        MeshParam%Ki(iLayer,iElem) = MeshParam%Ki(iLayer,iElem)*(MeshParam%Si(iLayer,iElem))**0.5
    EndIf

            
    return
        
    End Subroutine Saturation    