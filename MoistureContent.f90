    Subroutine MoistureContent(eta,etaplus,iElem,HydroParam,MeshParam)
    
    Use MeshVars 
    Use Hydrodynamic
    Implicit none
    Real :: eta, etaplus, V, DSmin, H1, H2, DSmax, vaux
    Integer :: iElem, iLayer
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam

    HydroParam%Vol(iElem) = MeshParam%Area(iElem)*etaplus
    HydroParam%Vol1(iElem) = 0.d0
    MeshParam%Si(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalM(iElem),iElem) = 1.d0
    MeshParam%Ki(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem),iElem)  = 0.d0    
    MeshParam%Ki(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalMs(iElem),iElem) = MeshParam%Ksat(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalMs(iElem),iElem)
    HydroParam%Ci(iElem) = 0.d0
    HydroParam%P(iElem) = MeshParam%Area(iElem)
    If (V(eta,HydroParam%sb(iElem)) >= 0.0d0 ) Then ! Wet Cell
        If (V(eta,HydroParam%hb(iElem)) > 0.0d0) Then 
            ! Superficial Layer:
            HydroParam%Vol(iElem)  = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*(Dot_Product(MeshParam%ei(:,iElem),HydroParam%DZsi(:,iElem)) + V(eta,HydroParam%hb(iElem)) ) !V1 = %Vol and V2 is zero.
            HydroParam%Vol1(iElem) =  MeshParam%Area(iElem)*(sum(HydroParam%DZsi(:,iElem)) + V(eta,HydroParam%hb(iElem)) )
            
            !HydroParam%Vol(iElem)  = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*(Dot_Product(MeshParam%ei(:,iElem),HydroParam%DZsi(:,iElem)) + V(eta,HydroParam%hb(iElem)) )
            
            !Haux = 0.d0
            !if(HydroParam%hb(iElem) > 0.d0) Then
            !    !DSmax = max(HydroParam%hj(MeshParam%Edge(1,iElem)),HydroParam%hj(MeshParam%Edge(2,iElem)),HydroParam%hj(MeshParam%Edge(3,iElem)),HydroParam%hj(MeshParam%Edge(4,iElem)))            
            !    !DSmin = min(HydroParam%hj(MeshParam%Edge(1,iElem)),HydroParam%hj(MeshParam%Edge(2,iElem)),HydroParam%hj(MeshParam%Edge(3,iElem)),HydroParam%hj(MeshParam%Edge(4,iElem)))  
            !    !
            !    DSmax = max(HydroParam%hj(MeshParam%Edge(2,iElem)),HydroParam%hb(iElem))
            !    DSmin = min(HydroParam%hj(MeshParam%Edge(2,iElem)),HydroParam%hb(iElem))
            !    H1 = 0.5*(DSmax - DSmin)*MeshParam%Area(iElem)
            !    
            !    DSmax = max(HydroParam%hj(MeshParam%Edge(4,iElem)),HydroParam%hb(iElem))
            !    DSmin = min(HydroParam%hj(MeshParam%Edge(4,iElem)),HydroParam%hb(iElem))
            !    H2 = 0.5*(DSmax - DSmin)*MeshParam%Area(iElem)             
            !    
            !endif
            !
            !H1 = 0.5*(-HydroParam%hj(MeshParam%Edge(2,iElem)) + HydroParam%hb(iElem))*MeshParam%Area(iElem)
            !if(H1<0.d0) then
            !    H1 = H1*MeshParam%ei(1,iElem)
            !endif
            !
            !H2 = 0.5*(-HydroParam%hj(MeshParam%Edge(4,iElem)) + HydroParam%hb(iElem))*MeshParam%Area(iElem)
            !if(H2<0.d0) then
            !    H2 = H2*MeshParam%ei(1,iElem)
            !endif
            ! 
            !if(H1/= 0.d0) then
            !    continue
            !endif
            !
            !
            !HydroParam%Vol(iElem)  = H1 + H2 + HydroParam%Vol(iElem) +  MeshParam%Area(iElem)*(Dot_Product(MeshParam%ei(:,iElem),HydroParam%DZsi(:,iElem)) + V(eta,HydroParam%hb(iElem)) )
            !
            !vaux = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*(Dot_Product(MeshParam%ei(:,iElem),HydroParam%DZsi(:,iElem)) + V(eta,HydroParam%hb(iElem)) )
            !HydroParam%Vol1(iElem) = MeshParam%Area(iElem)*etaplus + MeshParam%Area(iElem)*V(eta,HydroParam%sb(iElem))
            
            
            !!! Subsurface Layer:
            !Do iLayer = HydroParam%ElSmallms(iElem),HydroParam%ElCapitalMs(iElem)
            !    If(HydroParam%DZsi(iLayer,iElem) > 0 ) Then
            !        MeshParam%Si(iLayer,iElem) = 0.d0
            !        Call SoilSaturation(eta, HydroParam%Vol(iElem), HydroParam%Vol1(iElem), iLayer, iElem, HydroParam, MeshParam)
            !    EndIf
            !EndDo
            HydroParam%P(iElem) = MeshParam%Area(iElem)
            !HydroParam%Ci(iElem) = HydroParam%P(iElem) 
            HydroParam%Ci(iElem) = MeshParam%Area(iElem)
        Else
            HydroParam%Vol(iElem)  = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*MeshParam%ei(HydroParam%ElCapitalMs(iElem),iElem)*V(eta,HydroParam%sb(iElem)) !V1 = %Vol and V2 is zero.
            HydroParam%Vol1(iElem) =  MeshParam%Area(iElem)*V(eta,HydroParam%sb(iElem))
            
            HydroParam%P(iElem)  = MeshParam%Area(iElem)
            HydroParam%Ci(iElem) = MeshParam%Area(iElem)*MeshParam%ei(HydroParam%ElCapitalMs(iElem),iElem)
                            
            !H1 = 0.d0
            !if(eta - HydroParam%hj(MeshParam%Edge(2,iElem))> 0.d0) then
            !    H1 = (eta - HydroParam%hj(MeshParam%Edge(2,iElem)))*MeshParam%Area(iElem)*0.5
            !endif
            !
            !H2 = 0.d0
            !if(eta - HydroParam%hj(MeshParam%Edge(4,iElem))> 0.d0) then
            !    H2 = (eta - HydroParam%hj(MeshParam%Edge(4,iElem)))*MeshParam%Area(iElem)*0.5
            !endif
            ! 
            !if(H1/= 0.d0) then
            !    continue
            !endif    
            !
            !HydroParam%Vol(iElem)  = (1/MeshParam%ei(1,iElem)-MeshParam%ei(1,iElem))*H1 + (1/MeshParam%ei(1,iElem)-MeshParam%ei(1,iElem))*H2 + HydroParam%Vol(iElem) + MeshParam%Area(iElem)*MeshParam%ei(HydroParam%ElCapitalMs(iElem),iElem)*V(eta,HydroParam%sb(iElem)) !V1 = %Vol and V2 is zero.
            !! Only Subsurface Flow:
            !MeshParam%Si(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem),iElem) = 0.d0
            !Do iLayer = HydroParam%ElSmallms(iElem),HydroParam%ElCapitalMs(iElem)
            !    MeshParam%Si(iLayer,iElem) = 0.d0
            !    If(HydroParam%DZsi(iLayer,iElem) > 0 ) Then
            !        MeshParam%Si(iLayer,iElem) = 0.d0
            !        Call SoilSaturation(eta, HydroParam%Vol(iElem), HydroParam%Vol1(iElem), iLayer, iElem, HydroParam, MeshParam)
            !    EndIf
            !EndDo
            !HydroParam%P(iElem)  = MeshParam%Area(iElem)*MeshParam%ei(HydroParam%ElCapitalMs(iElem),iElem)*HydroParam%P(iElem)
            !HydroParam%Ci(iElem) = MeshParam%Area(iElem)*MeshParam%ei(HydroParam%ElCapitalMs(iElem),iElem)*HydroParam%Ci(iElem)
            !HydroParam%Ci(iElem) = 0.0d0
        EndIf    
    Else ! Dry Cell
        MeshParam%Si(:,iElem) = 0.d0
        MeshParam%Ki(:,iElem) = 0.d0
        HydroParam%P(iElem)   = MeshParam%Area(iElem)
        HydroParam%Ci(iElem)  = MeshParam%Area(iElem)
    EndIf

    End Subroutine MoistureContent  

    !
    !Subroutine Vol1(eta,Vol,iElem,HydroParam,MeshParam)
    !
    !Use MeshVars 
    !Use Hydrodynamic
    !Implicit none
    !Real :: eta, etaplus, V, Vol
    !Integer :: iElem, iLayer
    !type(MeshGridParam) :: MeshParam
    !type(HydrodynamicParam) :: HydroParam
    !
    !Vol = 0.d0
    !MeshParam%Si(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalM(iElem),iElem) = 1.d0
    !MeshParam%Ki(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem),iElem) = 0.d0    
    !MeshParam%Ki(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalMs(iElem),iElem) = MeshParam%Ksat(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalMs(iElem),iElem)
    !
    !If (V(eta,HydroParam%sb(iElem)) > 0.0d0 ) Then !Cell Wet
    !    If (V(eta,HydroParam%hb(iElem) > HydroParam%Pcri)) Then ! Superficial Layer
    !        HydroParam%P(iElem) = MeshParam%Area(iElem)
    !        Vol = MeshParam%Area(iElem)*(Dot_Product(MeshParam%ei(:,iElem),HydroParam%DZsi(:,iElem)) + V(eta,HydroParam%hb(iElem)) ) !V1 = %Vol and V2 is zero.
    !    Else
    !        !Subsurface Flow
    !        MeshParam%Si(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem),iElem) = 0.d0
    !        Do iLayer = HydroParam%ElSmallms(iElem),HydroParam%ElCapitalMs(iElem)
    !            MeshParam%Si(iLayer,iElem) = 0.d0
    !            Call SoilSaturation(HydroParam%hb(iElem)-HydroParam%PsiCrit(iElem), Vol, iLayer, iElem, HydroParam, MeshParam) 
    !        EndDo
    !        HydroParam%P = MeshParam%Area(iElem)*MeshParam%ei(HydroParam%ElCapitalMs(iElem),iElem)*HydroParam%P(iElem)
    !        Vol = Vol + HydroParam%P(iElem)*HydroParam%PsiCrit(iElem)
    !    EndIf    
    !        
    !EndIf
    !        
    !End Subroutine Vol1
    !
    
    Subroutine SoilSaturation(eta,Vol,V1,iLayer,iElem,HydroParam,MeshParam)
                 
    Use MeshVars 
    Use Hydrodynamic
            
    Implicit none
    
    Real :: SoilSat, mParam, eta, V, Kiaux, Z, Vol, V1, dz, iLayeraux
    !Integer :: iElem, iLayer,subfactor
    Integer :: iElem, iLayer
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam

    If(MeshParam%iSaturation == 0) Then !Darcy's Model
        
        If(HydroParam%Ze(iLayer,iElem) < eta < HydroParam%Ze(iLayer+1,iElem)) Then
            
            Vol = Vol + MeshParam%Area(iElem)*( V(eta,HydroParam%Ze(iLayer,iElem)) - V(eta,HydroParam%hb(iElem)))*MeshParam%ei(iLayer,iElem)
            V1  = V1  + MeshParam%Area(iElem)*( V(eta,HydroParam%Ze(iLayer,iElem)) - V(eta,HydroParam%hb(iElem) ))
            MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)
            MeshParam%Si(iLayer,iElem) = 1.0d0
            HydroParam%Ci(iElem) = 1.0d0
            
        ElseIf(eta >= HydroParam%Ze(iLayer+1,iElem)) Then
            
            Vol = Vol + MeshParam%Area(iElem)*(HydroParam%DZsi(iLayer,iElem))*MeshParam%ei(iLayer,iElem)
            V1  = V1  + MeshParam%Area(iElem)*(HydroParam%DZsi(iLayer,iElem))
            MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)
            MeshParam%Si(iLayer,iElem) = 1.0d0
            HydroParam%Ci(iElem) = 1.0d0/MeshParam%ei(iLayer,iElem)
            
        Else
            MeshParam%Ki(iLayer,iElem) = 0.d0
            MeshParam%Si(iLayer,iElem) = 0.d0
        EndIf
        HydroParam%P(iElem) = MeshParam%Si(iLayer,iElem)/MeshParam%ei(iLayer,iElem)
                
    ElseIf(MeshParam%iSaturation == 1) Then !Brooks and Corey's Model
        
        If(HydroParam%ElSmallms(iElem) == HydroParam%ElCapitalMs(iElem)) Then
            !subfactor = 100
            iLayeraux = 1
            Z = HydroParam%hb(iElem)
            mParam = (1-1/MeshParam%nSoil(iLayeraux,iElem))
            !dz = (HydroParam%hb(iElem)-Z)/MeshParam%subfactor 
            dz = (Z - HydroParam%sb(iElem))/MeshParam%subfactor
            !HydroParam%P(iElem)  =  1/(MeshParam%alpha(iLayeraux,iElem)*(-HydroParam%PsiCrit(iElem)))**(MeshParam%nSoil(iLayeraux,iElem))                
            If(eta - HydroParam%hb(iElem) < 0.d0 .and. eta - HydroParam%hb(iElem) < -HydroParam%PsiCrit(iElem)) Then
                Do While(Z > HydroParam%sb(iElem))
                    Z = Z - dz
                 !Z = Z + HydroParam%DZsi(HydroParam%ElSmallms(iElem),iElem)/subfactor
                    If(eta - Z + dz/2 <= -HydroParam%PsiCrit(iElem)) Then
                        MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayeraux,iElem)*(Z + dz/2 - eta))**(MeshParam%nSoil(iLayeraux,iElem))
                        MeshParam%Si(iLayer,iElem) = 1/MeshParam%Si(iLayer,iElem)
                        ! %P = Saturation(Critical Pressure)
                        HydroParam%P(iElem)  = HydroParam%P(iElem) + dz*MeshParam%nSoil(iLayeraux,iElem)*(-HydroParam%PsiCrit(iElem)/(eta-Z+dz/2))**(MeshParam%nSoil(iLayeraux,iElem) + 1)/HydroParam%PsiCrit(iElem)
                        HydroParam%Ci(iElem) = HydroParam%Ci(iElem) + dz*MeshParam%nSoil(iLayeraux,iElem)*(-HydroParam%PsiCrit(iElem)/(eta-Z+dz/2))**(MeshParam%nSoil(iLayeraux,iElem) + 1)/HydroParam%PsiCrit(iElem)
                        !V1    = V1 + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*MeshParam%Si(iLayer,iElem)
                        !Vol   = Vol + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*MeshParam%Si(iLayer,iElem)      
                        V1    = V1 + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*HydroParam%P(iElem) 
                        Vol   = Vol + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*HydroParam%P(iElem)                
                    Else
                        MeshParam%Si(iLayer,iElem) = 1.d0
                        HydroParam%Ci(iElem) = 1.d0
                        HydroParam%P(iElem)  = HydroParam%P(iElem) + dz*MeshParam%nSoil(iLayeraux,iElem)/HydroParam%PsiCrit(iElem)
                        V1    = V1 + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*HydroParam%P(iElem)
                        Vol   = Vol + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*MeshParam%Si(iLayer,iElem)
                    EndIf
                    !V1    = V1 + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*HydroParam%P(iElem)
                    !Vol   = Vol + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*MeshParam%Si(iLayer,iElem)
                    Kiaux = MeshParam%Ksat(iLayeraux,iElem)*MeshParam%Si(iLayer,iElem)**(3+2/MeshParam%nSoil(iLayeraux,iElem))
                    MeshParam%Ki(iLayeraux,iElem) = MeshParam%Ki(iLayeraux,iElem) + Kiaux
                    iLayer = iLayer + 1
                    
                EndDo
                MeshParam%Ki(iLayeraux,iElem) = MeshParam%Ki(iLayeraux,iElem)/MeshParam%subfactor
            Else
                !Soil satured:
                MeshParam%Si(:,iElem) = 1.d0
                HydroParam%P(iElem)  = 1/MeshParam%ei(iLayeraux,iElem)
                HydroParam%Ci(iElem) = 1.0d0
                
                MeshParam%Ki(iLayeraux,iElem) = MeshParam%Ksat(iLayeraux,iElem)
                V1  = V1 + MeshParam%Area(iElem)*Min(eta - HydroParam%sb(iElem),HydroParam%DZsi(iLayeraux,iElem))
                Vol = Vol + MeshParam%Area(iElem)*Min(eta - HydroParam%sb(iElem),HydroParam%DZsi(iLayeraux,iElem))*MeshParam%ei(iLayeraux,iElem)
            Endif
            !
            !If(eta - HydroParam%hb(iElem) < 0.d0 .and. eta - HydroParam%hb(iElem) < -HydroParam%PsiCrit(iElem)) Then
            !    Do While(Z <=  HydroParam%hb(iElem))
            !        !Z = Z + HydroParam%DZsi(HydroParam%ElSmallms(iElem),iElem)/subfactor
            !        Z = Z + dz
            !        If(eta - Z - dz/2 <= -HydroParam%PsiCrit(iElem)) Then
            !            MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayeraux,iElem)*(Z-dz/2 - eta))**(MeshParam%nSoil(iLayeraux,iElem))
            !            MeshParam%Si(iLayer,iElem) = 1/MeshParam%Si(iLayer,iElem)
            !            ! %P = Saturation(Critical Pressure)
            !            HydroParam%P(iElem)  =  HydroParam%P(iElem) + dz*MeshParam%ei(iLayeraux,iElem)*MeshParam%nSoil(iLayeraux,iElem)*(-HydroParam%PsiCrit(iElem)/(eta-Z-dz/2))**(MeshParam%nSoil(iLayeraux,iElem) + 1)/HydroParam%PsiCrit(iElem)
            !            HydroParam%Ci(iElem) = HydroParam%Ci(iElem) + dz*MeshParam%ei(iLayeraux,iElem)*MeshParam%nSoil(iLayeraux,iElem)*(-HydroParam%PsiCrit(iElem)/(eta-Z-dz/2))**(MeshParam%nSoil(iLayeraux,iElem) + 1)/HydroParam%PsiCrit(iElem)
            !            V1    = V1 + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*MeshParam%Si(iLayer,iElem)
            !            Vol   = Vol + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*MeshParam%Si(iLayer,iElem)                    
            !        Else
            !            MeshParam%Si(iLayer,iElem) = 1.d0
            !            !HydroParam%Ci(iElem) = 1.d0
            !            HydroParam%P(iElem)  = HydroParam%P(iElem) + dz*MeshParam%nSoil(iLayeraux,iElem)/HydroParam%PsiCrit(iElem)
            !            V1    = V1 + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*HydroParam%P(iElem)
            !            Vol   = Vol + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*MeshParam%Si(iLayer,iElem)
            !        EndIf
            !        !V1    = V1 + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*HydroParam%P(iElem)
            !        !Vol   = Vol + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*MeshParam%Si(iLayer,iElem)
            !        Kiaux = MeshParam%Ksat(iLayeraux,iElem)*MeshParam%Si(iLayer,iElem)**(3+2/MeshParam%nSoil(iLayeraux,iElem))
            !        MeshParam%Ki(iLayeraux,iElem) = MeshParam%Ki(iLayeraux,iElem) + Kiaux
            !        iLayer = iLayer + 1
            !    EndDo
            !    MeshParam%Ki(iLayeraux,iElem) = MeshParam%Ki(iLayeraux,iElem)/MeshParam%subfactor
            !Else
            !    !Soil satured:
            !    MeshParam%Si(:,iElem) = 1.d0
            !    HydroParam%P(iElem)  = 1/MeshParam%ei(iLayeraux,iElem)
            !    HydroParam%Ci(iElem) = 1.0d0
            !    
            !    MeshParam%Ki(iLayeraux,iElem) = MeshParam%Ksat(iLayeraux,iElem)
            !    V1 = MeshParam%Area(iElem)*Min(eta - HydroParam%sb(iElem),HydroParam%DZsi(iLayeraux,iElem))
            !    Vol = MeshParam%Area(iElem)*Min(eta - HydroParam%sb(iElem),HydroParam%DZsi(iLayeraux,iElem))*MeshParam%ei(iLayeraux,iElem)
            !Endif
            !
            
        Else              
            HydroParam%P(iElem) = 1.0d0
            If(eta - HydroParam%Zb(iLayer,iElem) <= -HydroParam%PsiCrit(iElem)) Then
                MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(HydroParam%Zb(iLayer+1,iElem) - eta))**(MeshParam%nSoil(iLayer,iElem))
                MeshParam%Si(iLayer,iElem) = 1/MeshParam%Si(iLayer,iElem)
                Vol = Vol + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)  
                V1  = V1  + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*HydroParam%P(iElem)
                HydroParam%Ci(iElem) = HydroParam%Ci(iElem) + HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%nSoil(iLayeraux,iElem)*(-HydroParam%PsiCrit(iElem)/(eta - HydroParam%Zb(iLayer,iElem)))**(MeshParam%nSoil(iLayeraux,iElem) + 1)/HydroParam%PsiCrit(iElem)
            !If(eta-HydroParam%Ze(iLayer+1,iElem) <= -HydroParam%PsiCrit(iElem)) Then
            !    HydroParam%P(iElem)  =  1.0d0
            !    !HydroParam%P(iElem) = 1/(MeshParam%alpha(iLayer,iElem)*(-HydroParam%PsiCrit(iElem)))**(MeshParam%nSoil(iLayer,iElem))
            !    
            !    If(eta-HydroParam%Ze(iLayer,iElem) < -HydroParam%PsiCrit(iElem)) Then 
            !        ! eta < Ze & Pressure in Cell <= Critical Pressure
            !        MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta))**(MeshParam%nSoil(iLayer,iElem))
            !        MeshParam%Si(iLayer,iElem) = 1/MeshParam%Si(iLayer,iElem)
            !        V1  = V1  + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*HydroParam%P(iElem)
            !        Vol = Vol + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)
            !        
            !        HydroParam%Ci(iElem) = HydroParam%Ci(iElem) + HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%nSoil(iLayeraux,iElem)*(-HydroParam%PsiCrit(iElem)/(eta - HydroParam%Ze(iLayer+1,iElem)))**(MeshParam%nSoil(iLayeraux,iElem) + 1)/HydroParam%PsiCrit(iElem)
            !        !HydroParam%Ci(iElem) = HydroParam%P(iElem)                
            !    Else
            !        ! Ze < eta < Ze + 1:  Cell Partially Wet
            !        !Vol to Ze - eta:
            !        V1  = V1  + MeshParam%Area(iElem)*(eta - HydroParam%Ze(iLayer,iElem) + HydroParam%PsiCrit(iElem))*MeshParam%ei(iLayer,iElem)
            !        Vol = Vol + MeshParam%Area(iElem)*(eta - HydroParam%Ze(iLayer,iElem) + HydroParam%PsiCrit(iElem))*MeshParam%ei(iLayer,iElem)
            !        
            !        !Vol to eta - Ze + 1:
            !        MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta - HydroParam%PsiCrit(iElem)))**(MeshParam%nSoil(iLayer,iElem))
            !        MeshParam%Si(iLayer,iElem) = 1/MeshParam%Si(iLayer,iElem) 
            !        Vol = Vol + MeshParam%Area(iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta - HydroParam%PsiCrit(iElem))*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)
            !        V1  = V1  + MeshParam%Area(iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta - HydroParam%PsiCrit(iElem))*MeshParam%ei(iLayer,iElem)*HydroParam%P(iElem)
            !        
            !        HydroParam%Ci(iElem) = HydroParam%Ci(iElem) + (HydroParam%Ze(iLayer+1,iElem) - eta)*MeshParam%nSoil(iLayeraux,iElem)*(-HydroParam%PsiCrit(iElem)/(eta-HydroParam%Ze(iLayer+1,iElem)))**(MeshParam%nSoil(iLayeraux,iElem) + 1)/HydroParam%PsiCrit(iElem)
            !        
            !        !Saturation Weighted to calc the hydraulic conductivity:
            !        MeshParam%Si(iLayer,iElem) = (MeshParam%Si(iLayer,iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta) + (eta - HydroParam%Ze(iLayer,iElem)))/HydroParam%DZsi(iLayer,iElem)
            !    EndIf
            Else
                ! eta > Zb + 1
                MeshParam%Si(:,iElem) = 1.d0
                HydroParam%Ci(iElem) = HydroParam%P(iElem)
                Vol = Vol + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)   
                V1  = V1  + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)
            EndIf
            !Vol = Vol + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)
            MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)*MeshParam%Si(iLayer,iElem)**(3+2/MeshParam%nSoil(iLayer,iElem))
        EndIf
    
    Else !van Genuchten's Model
        mParam = (1-1/MeshParam%nSoil(iLayer,iElem))
        If(HydroParam%ElSmallms(iElem) == HydroParam%ElCapitalMs(iElem)) Then
            !subfactor = 100
            iLayeraux = 1
            Z = HydroParam%sb(iElem)
            mParam = (1-1/MeshParam%nSoil(iLayer,iElem))
            dz = (HydroParam%hb(iElem)-Z)/MeshParam%subfactor  
            HydroParam%P(iElem)  =  1.0d0            
            Do While(Z <=  HydroParam%hb(iElem))
                Z = Z + dz
                If(eta - Z <= 0) Then
                    MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayeraux,iElem)*(Z - eta))**(MeshParam%nSoil(iLayeraux,iElem))
                    MeshParam%Si(iLayer,iElem) = 1/(1 + MeshParam%Si(iLayer,iElem))**mParam
                    If(eta - Z <= -HydroParam%PsiCrit(iElem)) Then
                        ! (Theta)zz is stricly positive and increasing for all eta - Z < -Psicrit and decreasing for eta - Z >= -Psicrit
                        ! Which means that the Psicrit is the Inflection Point of Theta
                        !HydroParam%P(iElem) = 1/( 1 + MeshParam%alpha(iLayeraux,iElem)*(-HydroParam%PsiCrit(iElem))**(MeshParam%nSoil(iLayeraux,iElem)))**mParam
                        HydroParam%Ci(iElem) = HydroParam%Ci(iElem) + dz*MeshParam%nSoil(iLayeraux,iElem)*MeshParam%alpha(iLayeraux,iElem)*mParam*abs(MeshParam%alpha(iLayeraux,iElem)*HydroParam%PsiCrit(iElem))**(MeshParam%nSoil(iLayeraux,iElem)-1)/(1 + abs(MeshParam%alpha(iLayeraux,iElem)*(Z - eta)**MeshParam%nSoil(iLayeraux,iElem)))**(mParam + 1)
                    Else
                        HydroParam%P(iElem) = MeshParam%Si(iLayer,iElem)
                    EndIf
                Else
                    MeshParam%Si(iLayer,iElem) = 1.d0
                    HydroParam%Ci(iElem) = HydroParam%P(iElem)
                EndIf
                
                Vol = Vol + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*MeshParam%Si(iLayer,iElem)
                V1  = V1  + MeshParam%Area(iElem)*dz*MeshParam%ei(iLayeraux,iElem)*HydroParam%P(iElem)
                Kiaux = MeshParam%Ksat(iLayeraux,iElem)*(1-(1-MeshParam%Si(iLayer,iElem)**(1/mParam))**mParam)**2
                MeshParam%Ki(iLayeraux,iElem) = MeshParam%Ki(iLayeraux,iElem) + Kiaux*(MeshParam%Si(iLayer,iElem))**0.5 
                iLayer = iLayer + 1
            EndDo
            
            MeshParam%Ki(iLayeraux,iElem) = MeshParam%Ki(iLayeraux,iElem)/MeshParam%subfactor
            
        Else
            HydroParam%P(iElem) = 1
            If(eta - HydroParam%Zb(iLayer,iElem) <= 0) Then
                
                MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(HydroParam%Zb(iLayer+1,iElem) - eta))**(MeshParam%nSoil(iLayer,iElem))
                MeshParam%Si(iLayer,iElem) = 1/(1 + MeshParam%Si(iLayer,iElem))**mParam
                Vol = Vol + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)  
                V1  = V1  + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*HydroParam%P(iElem)
                
                If(eta - HydroParam%Zb(iLayer,iElem) < -HydroParam%PsiCrit(iElem)) Then
                    HydroParam%Ci(iElem) = HydroParam%Ci(iElem) + HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%nSoil(iLayeraux,iElem)*MeshParam%alpha(iLayeraux,iElem)*mParam*abs(MeshParam%alpha(iLayeraux,iElem)*HydroParam%PsiCrit(iElem))**(MeshParam%nSoil(iLayeraux,iElem)-1)/(1 + abs(MeshParam%alpha(iLayeraux,iElem)*(eta - HydroParam%Zb(iLayer,iElem))**MeshParam%nSoil(iLayeraux,iElem)))**(mParam + 1)
                    !HydroParam%P(iElem) = HydroParam%P(iElem) + HydroParam%DZsi(iLayer,iElem)*MeshParam%nSoil(iLayeraux,iElem)*MeshParam%alpha(iLayeraux,iElem)*mParam*abs(MeshParam%alpha(iLayeraux,iElem)*HydroParam%PsiCrit(iElem))**(MeshParam%nSoil(iLayeraux,iElem)-1)/(1 + abs(MeshParam%alpha(iLayeraux,iElem)*(eta - HydroParam%Zb(iLayer,iElem))**MeshParam%nSoil(iLayeraux,iElem)))**(mParam + 1)
                Else
                    !HydroParam%P(iElem) = 1/( 1 + MeshParam%alpha(iLayer,iElem)*(-1.05*HydroParam%PsiCrit(iElem))**(MeshParam%nSoil(iLayer,iElem)))**mParam
                EndIf
                
 
            !If(eta - HydroParam%Ze(iLayer+1,iElem) <= 0) Then
            !    
            !    HydroParam%P(iElem) = 1/( 1 + MeshParam%alpha(iLayer,iElem)*(-HydroParam%PsiCrit(iElem))**(MeshParam%nSoil(iLayer,iElem)))**mParam
            !    
            !    If(eta-HydroParam%Ze(iLayer,iElem) <= 0) Then
            !        ! eta < Ze                    
            !        MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta))**(MeshParam%nSoil(iLayer,iElem))
            !        MeshParam%Si(iLayer,iElem) = 1/(1 + MeshParam%Si(iLayer,iElem))**mParam
            !        Vol = Vol + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)  
            !        V1  = V1  + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*HydroParam%P(iElem)
            !    Else
            !        ! Ze < eta < Ze + 1:  Cell Partially Wet
            !        !Vol to Ze - eta:
            !        V1  = V1  + MeshParam%Area(iElem)*(eta - HydroParam%Ze(iLayer,iElem) + HydroParam%PsiCrit(iElem))*MeshParam%ei(iLayer,iElem)
            !        Vol = Vol + MeshParam%Area(iElem)*(eta - HydroParam%Ze(iLayer,iElem) + HydroParam%PsiCrit(iElem))*MeshParam%ei(iLayer,iElem)
            !        
            !        ! Ze < eta < Ze + 1: 
            !        Vol = Vol + MeshParam%Area(iElem)*(eta - HydroParam%Ze(iLayer,iElem))*MeshParam%ei(iLayer,iElem)
            !        
            !        !Vol to eta - Ze + 1:
            !        MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta - HydroParam%PsiCrit(iElem)))**(MeshParam%nSoil(iLayer,iElem))
            !        MeshParam%Si(iLayer,iElem) = 1/(1 + MeshParam%Si(iLayer,iElem))**mParam
            !        Vol = Vol + MeshParam%Area(iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta - HydroParam%PsiCrit(iElem))*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)
            !        V1  = V1  + MeshParam%Area(iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta - HydroParam%PsiCrit(iElem))*MeshParam%ei(iLayer,iElem)*HydroParam%P(iElem)
            !        
            !        !Saturation Weighted to calc the hydraulic conductivity:
            !        MeshParam%Si(iLayer,iElem) = (MeshParam%Si(iLayer,iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta) + (eta - HydroParam%Ze(iLayer,iElem)))/HydroParam%DZsi(iLayer,iElem)
            !    EndIf
            !    
            !    
            !     If(eta-HydroParam%Ze(iLayer,iElem) < -HydroParam%PsiCrit(iElem)) Then
            !         
            !        HydroParam%Ci(iElem) = HydroParam%Ci(iElem) + HydroParam%DZsi(iLayer,iElem)*MeshParam%nSoil(iLayeraux,iElem)*MeshParam%alpha(iLayeraux,iElem)*mParam*abs(MeshParam%alpha(iLayeraux,iElem)*HydroParam%PsiCrit(iElem))**(MeshParam%nSoil(iLayeraux,iElem)-1)/(1 + abs(MeshParam%alpha(iLayeraux,iElem)*(HydroParam%Ze(iLayer+1,iElem)  - eta)**MeshParam%nSoil(iLayeraux,iElem)))**(mParam + 1)
            !     
            !     ElseIf(eta-HydroParam%Ze(iLayer+1,iElem) < -HydroParam%PsiCrit(iElem)) Then
            !         
            !        HydroParam%Ci(iElem) = HydroParam%Ci(iElem) + HydroParam%DZsi(iLayer,iElem)*MeshParam%nSoil(iLayeraux,iElem)*MeshParam%alpha(iLayeraux,iElem)*mParam*abs(MeshParam%alpha(iLayeraux,iElem)*HydroParam%PsiCrit(iElem))**(MeshParam%nSoil(iLayeraux,iElem)-1)/(1 + abs(MeshParam%alpha(iLayeraux,iElem)*(HydroParam%Ze(iLayer+1,iElem)  - eta)**MeshParam%nSoil(iLayeraux,iElem)))**(mParam + 1)
            !
            !     EndIf 
            !           
            Else
                MeshParam%Si(iLayer,iElem) = 1.d0 
                Vol = Vol + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)
                V1  = V1  + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)
            EndIf
            
            !Vol = Vol + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)
            MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)*(1-(1-MeshParam%Si(iLayer,iElem)**(1/mParam))**mParam)**2
            MeshParam%Ki(iLayer,iElem) = MeshParam%Ki(iLayer,iElem)*(MeshParam%Si(iLayer,iElem))**0.5        
        EndIf
    EndIf    
    return
        
    End Subroutine SoilSaturation
    
    
    !Subroutine MoistureContent(eta,etaplus,iElem,HydroParam,MeshParam)
  !  
  !  Use MeshVars 
  !  Use Hydrodynamic
  !  Implicit none
  !  Real :: eta, etaplus, V
  !  Integer :: iElem, iLayer
  !  type(MeshGridParam) :: MeshParam
  !  type(HydrodynamicParam) :: HydroParam
  !  HydroParam%PsiCrit(iElem) = 0.d0
  !  HydroParam%Vol(iElem) = MeshParam%Area(iElem)*etaplus
  !  MeshParam%Si(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalM(iElem),iElem) = 1.d0
  !  MeshParam%Ki(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem),iElem) = 0.d0    
  !  MeshParam%Ki(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalMs(iElem),iElem) = MeshParam%Ksat(HydroParam%ElSmallms(iElem):HydroParam%ElCapitalMs(iElem),iElem)
  !  If (eta <= HydroParam%hb(iElem) .and. eta > HydroParam%sb(iElem)) Then ! 
  !      MeshParam%Si(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem),iElem) = 0.d0
  !      !If eta is above the bathmetry, the soil can be unsatured, calculate saturation by models:
  !      Do iLayer = HydroParam%ElSmallms(iElem),HydroParam%ElCapitalMs(iElem)
  !          MeshParam%Si(iLayer,iElem) = 0.d0
  !          Call SoilSaturation(eta, iLayer, iElem, HydroParam, MeshParam) 
  !      EndDo
  !  ElseIf(eta > HydroParam%hb(iElem)) Then
  !      HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*(Dot_Product(MeshParam%ei(:,iElem)*MeshParam%Si(:,iElem),HydroParam%DZsi(:,iElem)) + V(eta,HydroParam%hb(iElem)) )
  !  Else
  !      MeshParam%Si(:,iElem) = 0.d0
  !      MeshParam%Ki(:,iElem) = 0.d0
  !  EndIf
  !  
  !  End Subroutine MoistureContent    
  !
  !  Subroutine SoilSaturation(eta, iLayer, iElem, HydroParam, MeshParam)
  !               
  !  Use MeshVars 
  !  Use Hydrodynamic
  !          
  !  Implicit none
  !  
  !  Real :: SoilSat, mParam, eta, V, Kiaux, Z 
  !  Integer :: iElem, iLayer,subfactor
  !  type(MeshGridParam) :: MeshParam
  !  type(HydrodynamicParam) :: HydroParam
  !
  !  If(MeshParam%iSaturation == 0) Then !Darcy's Model
  !      HydroParam%PsiCrit(iElem) = 0.d0
  !      If(HydroParam%Ze(iLayer,iElem) < eta < HydroParam%Ze(iLayer+1,iElem)) Then
  !          HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*(V(eta,HydroParam%Ze(iLayer,iElem)))*MeshParam%ei(iLayer,iElem)
  !          MeshParam%Si(iLayer,iElem) = Min(1.d0,(eta - HydroParam%Ze(iLayer,iElem))/HydroParam%DZsi(iLayer,iElem)) !layer thickness wet percentage
  !          MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)
  !      ElseIf(eta >= HydroParam%Ze(iLayer+1,iElem)) Then
  !          HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)
  !          MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)
  !      Else
  !          MeshParam%Ki(iLayer,iElem) = 0.d0
  !      EndIf        
  !  ElseIf(MeshParam%iSaturation == 1) Then ! Brooks and Corey's Model
  !      HydroParam%PsiCrit(iElem) = 1/MeshParam%alpha(HydroParam%ElCapitalMs(iElem),iElem)
  !      If(HydroParam%ElSmallms(iElem) == HydroParam%ElCapitalMs(iElem)) Then
  !          subfactor = 100
  !          Z = 0.d0
  !          mParam = (1-1/MeshParam%nSoil(iLayer,iElem))
  !          Do While(Z <= HydroParam%hb(iElem))
  !              Z = Z + HydroParam%DZsi(HydroParam%ElSmallms(iElem),iElem)/subfactor
  !              If(eta < Z) Then
  !                  MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(Z - eta))**(MeshParam%nSoil(iLayer,iElem))
  !                  MeshParam%Si(iLayer,iElem) = 1/MeshParam%Si(iLayer,iElem)
  !              Else
  !                  MeshParam%Si(iLayer,iElem) = 1.d0 
  !              EndIf
  !              HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*(HydroParam%DZsi(HydroParam%ElSmallms(iElem),iElem)/subfactor)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)
  !              Kiaux = MeshParam%Ksat(iLayer,iElem)*MeshParam%Si(iLayer,iElem)**(3+2/MeshParam%nSoil(iLayer,iElem))
  !              MeshParam%Ki(iLayer,iElem) = MeshParam%Ki(iLayer,iElem) + Kiaux
  !          EndDo
  !          MeshParam%Ki(iLayer,iElem) = MeshParam%Ki(iLayer,iElem)/subfactor
  !      Else
  !          If(eta < HydroParam%Ze(iLayer+1,iElem) - 1/MeshParam%alpha(iLayer,iElem)) Then
  !              MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta))**(MeshParam%nSoil(iLayer,iElem))
  !              MeshParam%Si(iLayer,iElem) = 1/MeshParam%Si(iLayer,iElem)
  !          Else
  !              MeshParam%Si(iLayer,iElem) = 1.d0
  !          EndIf
  !          HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)
  !          MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)*MeshParam%Si(iLayer,iElem)**(3+2/MeshParam%nSoil(iLayer,iElem))
  !      EndIf
  !  Else !van Genuchten's Model
  !      HydroParam%PsiCrit(iElem) = (MeshParam%nSoil(HydroParam%ElCapitalMs(iElem),iElem) - 1)**(1/MeshParam%nSoil(HydroParam%ElCapitalMs(iElem),iElem))/(MeshParam%alpha(HydroParam%ElCapitalMs(iElem),iElem)*MeshParam%nSoil(HydroParam%ElCapitalMs(iElem),iElem))        
  !      If(HydroParam%ElSmallms(iElem) == HydroParam%ElCapitalMs(iElem)) Then
  !          subfactor = 100
  !          Z = 0.d0
  !          mParam = (1-1/MeshParam%nSoil(iLayer,iElem))
  !          Do While(Z <= HydroParam%hb(iElem))
  !              Z = Z + HydroParam%DZsi(HydroParam%ElSmallms(iElem),iElem)/subfactor
  !              If(eta < Z) Then
  !                  MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(Z - eta))**(MeshParam%nSoil(iLayer,iElem))
  !                  MeshParam%Si(iLayer,iElem) = 1/(1 + MeshParam%Si(iLayer,iElem))**mParam
  !              Else
  !                  MeshParam%Si(iLayer,iElem) = 1.d0 
  !              EndIf
  !              HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*(HydroParam%DZsi(HydroParam%ElSmallms(iElem),iElem)/subfactor)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)
  !              Kiaux = MeshParam%Ksat(iLayer,iElem)*(1-(1-MeshParam%Si(iLayer,iElem)**(1/mParam))**mParam)**2
  !              MeshParam%Ki(iLayer,iElem) = MeshParam%Ki(iLayer,iElem) + Kiaux*(MeshParam%Si(iLayer,iElem))**0.5
  !          EndDo
  !          MeshParam%Ki(iLayer,iElem) = MeshParam%Ki(iLayer,iElem)/subfactor
  !      Else
  !          mParam = (1-1/MeshParam%nSoil(iLayer,iElem))
  !          If(eta < HydroParam%Ze(iLayer+1,iElem)) Then
  !              MeshParam%Si(iLayer,iElem) = (MeshParam%alpha(iLayer,iElem)*(HydroParam%Ze(iLayer+1,iElem) - eta))**(MeshParam%nSoil(iLayer,iElem))
  !              MeshParam%Si(iLayer,iElem) = 1/(1 + MeshParam%Si(iLayer,iElem))**mParam
  !          Else
  !              MeshParam%Si(iLayer,iElem) = 1.d0 
  !          EndIf
  !          HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*HydroParam%DZsi(iLayer,iElem)*MeshParam%ei(iLayer,iElem)*MeshParam%Si(iLayer,iElem)
  !          MeshParam%Ki(iLayer,iElem) = MeshParam%Ksat(iLayer,iElem)*(1-(1-MeshParam%Si(iLayer,iElem)**(1/mParam))**mParam)**2
  !          MeshParam%Ki(iLayer,iElem) = MeshParam%Ki(iLayer,iElem)*(MeshParam%Si(iLayer,iElem))**0.5        
  !      EndIf
  !  EndIf    
  !  return
  !      
  !  End Subroutine SoilSaturation