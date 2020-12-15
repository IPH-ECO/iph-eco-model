  Subroutine SoilSaturation(eta,iElem,HydroParam,MeshParam)
    
    Use MeshVars 
    Use Hydrodynamic
    Implicit none
    Real :: eta
    Integer :: nLayers, iElem, iLayer
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    
    MeshParam%Si(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem),iElem) = 0.d0
    MeshParam%Ki(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem),iElem) = 0.d0
    If(HydroParam%ElSmallms(iElem) == HydroParam%ElCapitalMs(iElem)) Then !One- two-Dimensional, layer satured
        MeshParam%Si(HydroParam%ElSmallms(iElem),iElem) = 1.d0
        MeshParam%Ki(HydroParam%ElSmallms(iElem),iElem) = MeshParam%Ksat(HydroParam%ElSmallms(iElem)),iElem)
    ElseIf (eta - HydroParam%Pcri<= HydroParam%hb(iElem)) Then ! Three-Dimensional
        !If eta is above the bathmetry, the soil can be unsatured, calculate saturation by models:
        nLayers = Max(Floor(HydroParam%ElCapitalMs(iElem)-HydroParam%ElSmallms(iElem))/2,1)
        If (eta > sum(HydroParam%DZsi(:nLayers))) Then
            nLayers = nLayers + 1
            Do iLayer = nLayers,HydroParam%ElCapitalMs(iElem)
                Call Saturation(eta, iLayer, iElem, HydroParam, MeshParam) 
            EndDo
        Else
            Do iLayer = HydroParam%ElSmallms(iElem),nLayers
                Call Saturation(eta, iLayer, iElem, HydroParam, MeshParam) 
            EndDo
        EndIf
    Else
        !Soil satured:
        MeshParam%Si(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalM(iElem),iElem) = 1.d0
        MeshParam%Ki(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalMs(iElem),iElem) = MeshParam%Ksat(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalMs(iElem),iElem)
    EndIf
    
    End Subroutine SoilSaturation    

    Subroutine Saturation(eta, iLayer, iElem, HydroParam, MeshParam)
                 
    Use MeshVars 
    Use Hydrodynamic
            
    Implicit none
    
    Real :: SoilSat, mParam, eta
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam

    If(MeshParam%iSaturation == 0) Then !Darcy's Model
        MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)*MeshParam%Si(iLayer,iElem) 
    ElseIf(MeshParam%iSaturation == 1) Then ! Brooks and Corey's Model
        If(eta < HydroParam%Zb(iLayer,iElem) - 1/MeshParam%alpha(iLayer,iElem)) Then
            MeshParam%Si(iLayer,iElem) = MeshParam%alpha(iLayer,iElem)*(HydroParam%Zb(iLayer,iElem) - eta)**(MeshParam%nSoil(iLayer,iElem))
            MeshParam%Si(iLayer,iElem) = 1/MeshParam%Si(iLayer,iElem)
        EndIf
        MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)*MeshParam%Si(iLayer,iElem)**(3+2/MeshParam%nSoil(iLayer,iElem))
    Else !van Genuchten's Model
        mParam = (1-1/MeshParam%nSoil(iLayer,iElem))
        If(eta < HydroParam%Zb(iLayer,iElem)) Then
            MeshParam%Si(iLayer,iElem) = 1/(1 + MeshParam%alpha(iLayer,iElem)*(HydroParam%Zb(iLayer,iElem) - eta))**mParam
        EndIf
        MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)*(1-(1-MeshParam%Si(iLayer,iElem)**(1/mParam))**mParam)**2
        MeshParam%Ki(iLayer,iElem) = MeshParam%Ki(iLayer,iElem)*(MeshParam%Si(iLayer,iElem))**0.5
    EndIf

            
    return
        
    End Subroutine Saturation    