!> This subroutine gets the values of water quality boundary conditions. 
Subroutine GetWQBoundaryConditions(HydroParam,LimnoParam,time)
    
    Use Hydrodynamic
    Use LimnologyVars
    
    Implicit none
    
    Integer:: i,iElem,iLayer
    Real:: time
    Real:: t_interp(1)
    Real:: p_interp(1,1)
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam

    t_interp = time ! Current simulation time
    
    !1. Getting the value of water temperature in the current simulation time 
    LimnoParam%uDLoadTemp = LimnoParam%sDTempW
    Do i =1,LimnoParam%NuDLoadTemp
        iElem = LimnoParam%NTempIndex(i,2)
        Call interp_linear( 1, LimnoParam%TempnTime(i), LimnoParam%TempTime(i,1:LimnoParam%TempnTime(i)), LimnoParam%TempValue(i,1:LimnoParam%TempnTime(i)), 1, t_interp, p_interp )
        Do iLayer = LimnoParam%TempSmallm(i),LimnoParam%TempCapitalM(i)
            LimnoParam%uDLoadTemp(iLayer,iElem) = p_interp(1,1)
        EndDo
    EndDo

    !2. Getting the value of Salinity in the current simulation time 
    LimnoParam%uDLoadSal = LimnoParam%sDSal
    Do i =1,LimnoParam%NuDLoadSal
        iElem = LimnoParam%NSalIndex(i,2)
        Call interp_linear( 1, LimnoParam%SalnTime(i), LimnoParam%SalTime(i,1:LimnoParam%SalnTime(i)), LimnoParam%SalValue(i,1:LimnoParam%SalnTime(i)), 1, t_interp, p_interp )
        Do iLayer = LimnoParam%SalSmallm(i),LimnoParam%SalCapitalM(i)
            LimnoParam%uDLoadSal(iLayer,iElem) = p_interp(1,1)
        EndDo
    EndDo
    
    !3. Getting the value of Inorganic Matter in the current simulation time  
    LimnoParam%uDLoadIM = LimnoParam%sDIMW
    Do i =1,LimnoParam%NuDLoadIM
        iElem = LimnoParam%NIMIndex(i,2)
        Call interp_linear( 1, LimnoParam%IMnTime(i), LimnoParam%IMTime(i,1:LimnoParam%IMnTime(i)), LimnoParam%IMValue(i,1:LimnoParam%IMnTime(i)), 1, t_interp, p_interp )
        Do iLayer = LimnoParam%IMSmallm(i),LimnoParam%IMCapitalM(i)
            LimnoParam%uDLoadIM(iLayer,iElem) = p_interp(1,1)
        EndDo
    EndDo
    
    !4. Getting the value of  Organic Matter in the current simulation time  
    
    If (LimnoParam%InclMatOrgSplit == 1) Then
        If (LimnoParam%InclPomDomSplit == 1) Then
            LimnoParam%uDLoadDom = LimnoParam%sDDomW(:,:,2)
            Do i =1,LimnoParam%NuDLoadDom
                iElem = LimnoParam%NDomIndex(i,2)
                Call interp_linear( 1, LimnoParam%DomnTime(i), LimnoParam%DomTime(i,1:LimnoParam%DomnTime(i)), LimnoParam%DomValue(i,1:LimnoParam%DomnTime(i)), 1, t_interp, p_interp )
                Do iLayer = LimnoParam%DomSmallm(i),LimnoParam%DomCapitalM(i)
                    LimnoParam%uDLoadDom(iLayer,iElem) = p_interp(1,1)
                EndDo
            EndDo
    
            LimnoParam%uDLoadPom = LimnoParam%sDPomW(:,:,2)
            Do i =1,LimnoParam%NuDLoadPom
                iElem = LimnoParam%NPomIndex(i,2)
                Call interp_linear( 1, LimnoParam%PomnTime(i), LimnoParam%PomTime(i,1:LimnoParam%PomnTime(i)), LimnoParam%PomValue(i,1:LimnoParam%PomnTime(i)), 1, t_interp, p_interp )
                Do iLayer = LimnoParam%PomSmallm(i),LimnoParam%PomCapitalM(i)
                    LimnoParam%uDLoadPom(iLayer,iElem) = p_interp(1,1)
                EndDo
            EndDo
        Else
            LimnoParam%uDLoadDom = LimnoParam%sDDomW(:,:,1)
            Do i =1,LimnoParam%NuDLoadDom
                iElem = LimnoParam%NDomIndex(i,2)
                Call interp_linear( 1, LimnoParam%DomnTime(i), LimnoParam%DomTime(i,1:LimnoParam%DomnTime(i)), LimnoParam%DomValue(i,1:LimnoParam%DomnTime(i)), 1, t_interp, p_interp )
                Do iLayer = LimnoParam%DomSmallm(i),LimnoParam%DomCapitalM(i)
                    LimnoParam%uDLoadDom(iLayer,iElem) = p_interp(1,1)
                EndDo
            EndDo
    
            LimnoParam%uDLoadPom = LimnoParam%sDPomW(:,:,1)
            Do i =1,LimnoParam%NuDLoadPom
                iElem = LimnoParam%NPomIndex(i,2)
                Call interp_linear( 1, LimnoParam%PomnTime(i), LimnoParam%PomTime(i,1:LimnoParam%PomnTime(i)), LimnoParam%PomValue(i,1:LimnoParam%PomnTime(i)), 1, t_interp, p_interp )
                Do iLayer = LimnoParam%PomSmallm(i),LimnoParam%PomCapitalM(i)
                    LimnoParam%uDLoadPom(iLayer,iElem) = p_interp(1,1)
                EndDo
            EndDo
        EndIf
    Else
        LimnoParam%uDLoadPom = LimnoParam%sDPomW(:,:,1)
        Do i =1,LimnoParam%NuDLoadPom
            iElem = LimnoParam%NPomIndex(i,2)
            Call interp_linear( 1, LimnoParam%PomnTime(i), LimnoParam%PomTime(i,1:LimnoParam%PomnTime(i)), LimnoParam%PomValue(i,1:LimnoParam%PomnTime(i)), 1, t_interp, p_interp )
            Do iLayer = LimnoParam%PomSmallm(i),LimnoParam%PomCapitalM(i)
                LimnoParam%uDLoadPom(iLayer,iElem) = p_interp(1,1)
            EndDo
        EndDo    
    EndIf
    
        

    
    !6. Getting the value of Orthophosphate (PO4) in the current simulation time  
    LimnoParam%uPLoadPO4 = LimnoParam%sPO4W
    Do i =1,LimnoParam%NuPLoadPO4
        iElem = LimnoParam%NPO4Index(i,2)
        Call interp_linear( 1, LimnoParam%PO4nTime(i), LimnoParam%PO4Time(i,1:LimnoParam%PO4nTime(i)), LimnoParam%PO4Value(i,1:LimnoParam%PO4nTime(i)), 1, t_interp, p_interp )
        Do iLayer = LimnoParam%PO4Smallm(i),LimnoParam%PO4CapitalM(i)
            LimnoParam%uPLoadPO4(iLayer,iElem) = p_interp(1,1)
        EndDo
    EndDo
    
    !7. Getting the value of Particulate Adsorbed Inorganic Phosphorus (PAIM) in the current simulation time  
    LimnoParam%uPLoadPAIM = LimnoParam%sPAIMW
    Do i =1,LimnoParam%NuPLoadPAIM
        iElem = LimnoParam%NPAIMIndex(i,2)
        Call interp_linear( 1, LimnoParam%PAIMnTime(i), LimnoParam%PAIMTime(i,1:LimnoParam%PAIMnTime(i)), LimnoParam%PAIMValue(i,1:LimnoParam%PAIMnTime(i)), 1, t_interp, p_interp )
        Do iLayer = LimnoParam%PAIMSmallm(i),LimnoParam%PAIMCapitalM(i)
            LimnoParam%uPLoadPAIM(iLayer,iElem) = p_interp(1,1)
        EndDo
    EndDo
    
    !8. Getting the value of Ammonium (NH4) in the current simulation time  
    LimnoParam%uNLoadNH4 = LimnoParam%sNH4W
    Do i =1,LimnoParam%NuNLoadNH4
        iElem = LimnoParam%NNH4Index(i,2)
        Call interp_linear( 1, LimnoParam%NH4nTime(i), LimnoParam%NH4Time(i,1:LimnoParam%NH4nTime(i)), LimnoParam%NH4Value(i,1:LimnoParam%NH4nTime(i)), 1, t_interp, p_interp )
        Do iLayer = LimnoParam%NH4Smallm(i),LimnoParam%NH4CapitalM(i)
            LimnoParam%uNLoadNH4(iLayer,iElem) = p_interp(1,1)
        EndDo
    EndDo
    
    !9. Getting the value of Nitrate (NO3) in the current simulation time  
    LimnoParam%uNLoadNO3 = LimnoParam%sNO3W
    Do i =1,LimnoParam%NuNLoadNO3
        iElem = LimnoParam%NNO3Index(i,2)
        Call interp_linear( 1, LimnoParam%NO3nTime(i), LimnoParam%NO3Time(i,1:LimnoParam%NO3nTime(i)), LimnoParam%NO3Value(i,1:LimnoParam%NO3nTime(i)), 1, t_interp, p_interp )
        Do iLayer = LimnoParam%NO3Smallm(i),LimnoParam%NO3CapitalM(i)
            LimnoParam%uNLoadNO3(iLayer,iElem) = p_interp(1,1)
        EndDo
    EndDo
    
    !10. Getting the value of Silicate (SiO2) in the current simulation time  
    LimnoParam%uSiLoadSiO2 = LimnoParam%sSiO2W
    Do i =1,LimnoParam%NuSiLoadSiO2
        iElem = LimnoParam%NSiO2Index(i,2)
        Call interp_linear( 1, LimnoParam%SiO2nTime(i), LimnoParam%SiO2Time(i,1:LimnoParam%SiO2nTime(i)), LimnoParam%SiO2Value(i,1:LimnoParam%SiO2nTime(i)), 1, t_interp, p_interp )
        Do iLayer = LimnoParam%SiO2Smallm(i),LimnoParam%SiO2CapitalM(i)
            LimnoParam%uSiLoadSiO2(iLayer,iElem) = p_interp(1,1)
        EndDo
    EndDo
    
    !11. Getting the value of Dissolved Oxygen (O2) in the current simulation time  
    LimnoParam%uO2LoadO2 = LimnoParam%sO2W
    Do i =1,LimnoParam%NuO2LoadO2
        iElem = LimnoParam%NO2Index(i,2)
        Call interp_linear( 1, LimnoParam%O2nTime(i), LimnoParam%O2Time(i,1:LimnoParam%O2nTime(i)), LimnoParam%O2Value(i,1:LimnoParam%O2nTime(i)), 1, t_interp, p_interp )
        Do iLayer = LimnoParam%O2Smallm(i),LimnoParam%O2CapitalM(i)
            LimnoParam%uO2LoadO2(iLayer,iElem) = p_interp(1,1)
        EndDo
    EndDo
    
    !12. Getting the value of Biological Oxygen Demand (BOD) in the current simulation time  
    LimnoParam%uDLoadDBO = LimnoParam%sDBOW
    Do i =1,LimnoParam%NuDBOLoadDBO
        iElem = LimnoParam%NDBOIndex(i,2)
        Call interp_linear( 1, LimnoParam%DBOnTime(i), LimnoParam%DBOTime(i,1:LimnoParam%DBOnTime(i)), LimnoParam%DBOValue(i,1:LimnoParam%DBOnTime(i)), 1, t_interp, p_interp )
        Do iLayer = LimnoParam%DBOSmallm(i),LimnoParam%DBOCapitalM(i)
            LimnoParam%uDLoadDBO(iLayer,iElem) = p_interp(1,1)
        EndDo
    EndDo
    
    
    !13. Getting the value of Dissolved Inorganic Carbon in the current simulation time  
    LimnoParam%uDLoadDic = LimnoParam%sDicW
    Do i =1,LimnoParam%NuDicLoadDic
        iElem = LimnoParam%NDicIndex(i,2)
        Call interp_linear( 1, LimnoParam%DicnTime(i), LimnoParam%DicTime(i,1:LimnoParam%DicnTime(i)), LimnoParam%DicValue(i,1:LimnoParam%DicnTime(i)), 1, t_interp, p_interp )
        Do iLayer = LimnoParam%DicSmallm(i),LimnoParam%DicCapitalM(i)
            LimnoParam%uDLoadDic(iLayer,iElem) = p_interp(1,1)
        EndDo
    EndDo
    
   
End Subroutine GetWQBoundaryConditions
    