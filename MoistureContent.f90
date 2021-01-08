  Subroutine MoistureContent(eta,etaplus,iElem,HydroParam,MeshParam)
    
    Use MeshVars 
    Use Hydrodynamic
    Implicit none
    Real :: eta, etaplus, V
    Integer :: iElem, iLayer
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    
    HydroParam%Vol(iElem) = MeshParam%Area(iElem)*etaplus
    MeshParam%Si(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalM(iElem),iElem) = 1.d0
    MeshParam%Ki(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem),iElem) = 0.d0    
    MeshParam%Ki(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalMs(iElem),iElem) = MeshParam%Ksat(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalMs(iElem),iElem)
    If (eta <= HydroParam%hb(iElem) .and. eta > HydroParam%sb(iElem)) Then ! 
        MeshParam%Si(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem),iElem) = 0.d0
        !If eta is above the bathmetry, the soil can be unsatured, calculate saturation by models:
        Do iLayer = HydroParam%ElSmallms(iElem),HydroParam%ElCapitalMs(iElem)
            MeshParam%Si(iLayer,iElem) = 0.d0
            Call SoilSaturation(eta, iLayer, iElem, HydroParam, MeshParam) 
        EndDo
    ElseIf(eta > HydroParam%hb(iElem)) Then
        HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*(Dot_Product(MeshParam%ei(:,iElem)*MeshParam%Si(:,iElem),HydroParam%DZsi(:,iElem)) + V(eta,HydroParam%hb(iElem)) )
    Else
        MeshParam%Si(:,iElem) = 0.d0
        MeshParam%Ki(:,iElem) = 0.d0
    EndIf
    
    End Subroutine MoistureContent    

    Subroutine SoilSaturation(eta, iLayer, iElem, HydroParam, MeshParam)
                 
    Use MeshVars 
    Use Hydrodynamic
            
    Implicit none
    
    Real :: SoilSat, mParam, eta, V, Kiaux, Z 
    Integer :: iElem, iLayer,subfactor
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam

    If(MeshParam%iSaturation == 0) Then !Darcy's Model
        If(HydroParam%Ze(iLayer,iElem) < eta < HydroParam%Ze(iLayer+1,iElem)) Then
            HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*(V(eta,HydroParam%Ze(iLayer,iElem)))*MeshParam%ei(iLayer,iElem)
            MeshParam%Si(iLayer,iElem) = Min(1.d0,(eta - HydroParam%Ze(iLayer,iElem))/HydroParam%DZsi(iLayer,iElem)) !layer thickness wet percentage
            MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)
        ElseIf(eta >= HydroParam%Ze(iLayer+1,iElem)) Then
            HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)
            MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)
        Else
            MeshParam%Ki(iLayer,iElem) = 0.d0
        EndIf        
    ElseIf(MeshParam%iSaturation == 1) Then ! Brooks and Corey's Model
        If(eta < HydroParam%Zb(iLayer,iElem) - 1/MeshParam%alpha(iLayer,iElem)) Then
            MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(HydroParam%Zb(iLayer,iElem) - eta))**(MeshParam%nSoil(iLayer,iElem))
            MeshParam%Si(iLayer,iElem) = 1/MeshParam%Si(iLayer,iElem)
        Else
            MeshParam%Si(iLayer,iElem) = 1.d0
        EndIf
        HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)
        MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)*MeshParam%Si(iLayer,iElem)**(3+2/MeshParam%nSoil(iLayer,iElem))
    Else !van Genuchten's Model
        If(HydroParam%ElSmallms(iElem) == HydroParam%ElCapitalMs(iElem)) Then
            subfactor = 100
            Z = 0.d0
            mParam = (1-1/MeshParam%nSoil(iLayer,iElem))
            Do While(Z <= HydroParam%hb(iElem))
                Z = Z + HydroParam%DZsi(HydroParam%ElSmallms(iElem),iElem)/subfactor
                If(eta < Z) Then
                    MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(Z - eta))**(MeshParam%nSoil(iLayer,iElem))
                    MeshParam%Si(iLayer,iElem) = 1/(1 + MeshParam%Si(iLayer,iElem))**mParam
                Else
                    MeshParam%Si(iLayer,iElem) = 1.d0 
                EndIf
                HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*(HydroParam%DZsi(HydroParam%ElSmallms(iElem),iElem)/subfactor)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)
                Kiaux = MeshParam%Ksat(iLayer,iElem)*(1-(1-MeshParam%Si(iLayer,iElem)**(1/mParam))**mParam)**2
                MeshParam%Ki(iLayer,iElem) = MeshParam%Ki(iLayer,iElem) + Kiaux*(MeshParam%Si(iLayer,iElem))**0.5
            EndDo
            MeshParam%Ki(iLayer,iElem) = MeshParam%Ki(iLayer,iElem)/subfactor
        Else
            mParam = (1-1/MeshParam%nSoil(iLayer,iElem))
            If(eta < HydroParam%Ze(iLayer+1,iElem)) Then
                MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta))**(MeshParam%nSoil(iLayer,iElem))
                MeshParam%Si(iLayer,iElem) = 1/(1 + MeshParam%Si(iLayer,iElem))**mParam
            Else
                MeshParam%Si(iLayer,iElem) = 1.d0 
            EndIf
            HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)
            MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)*(1-(1-MeshParam%Si(iLayer,iElem)**(1/mParam))**mParam)**2
            MeshParam%Ki(iLayer,iElem) = MeshParam%Ki(iLayer,iElem)*(MeshParam%Si(iLayer,iElem))**0.5        
        EndIf
    EndIf    
    return
        
    End Subroutine SoilSaturation   