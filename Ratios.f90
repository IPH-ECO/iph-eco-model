Subroutine Ratios(HydroParam,MeshParam,LimnoParam)
    
    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars

    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    Integer:: index,iElem,iLayer,gg,pom,dom,bb,zz,ff,mm
    Real:: V
    Real:: NearZero = 1e-10
    
    Do iElem = 1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            !-----------------------------------------------------------------------
            !  Current POM ratios(check the curent state)
            !-----------------------------------------------------------------------
            Do pom = 1, LimnoParam%npom
                LimnoParam%rCDPomW(iLayer,iElem,pom) = LimnoParam%sCPomW(iLayer,iElem,pom)/(LimnoParam%sDPomW(iLayer,iElem,pom)+NearZero)
                LimnoParam%rPDPomW(iLayer,iElem,pom) = LimnoParam%sPPomW(iLayer,iElem,pom)/(LimnoParam%sDPomW(iLayer,iElem,pom)+NearZero)
                LimnoParam%rNDPomW(iLayer,iElem,pom) = LimnoParam%sNPomW(iLayer,iElem,pom)/(LimnoParam%sDPomW(iLayer,iElem,pom)+NearZero)
            EndDo
            !-----------------------------------------------------------------------
            !  Current DOM ratios(check the curent state)
            !-----------------------------------------------------------------------
            Do dom = 1, LimnoParam%ndom
                LimnoParam%rCDDomW(iLayer,iElem,dom) = LimnoParam%sCDomW(iLayer,iElem,dom)/(LimnoParam%sDDomW(iLayer,iElem,dom)+NearZero)
                LimnoParam%rPDDomW(iLayer,iElem,dom) = LimnoParam%sPDomW(iLayer,iElem,dom)/(LimnoParam%sDDomW(iLayer,iElem,dom)+NearZero)
                LimnoParam%rNDDomW(iLayer,iElem,dom) = LimnoParam%sNDomW(iLayer,iElem,dom)/(LimnoParam%sDDomW(iLayer,iElem,dom)+NearZero)
            EndDo
    
            !-----------------------------------------------------------------------
            !   Current Phytoplankton in the Water ratios(check the curent state)
            !-----------------------------------------------------------------------
            !  Local status
            Do gg = 1, LimnoParam%nphy
                !   C/D_ratio_of_Phytoplankton
                LimnoParam%rCDPhytW(iLayer,iElem,gg) = LimnoParam%sCPhytW(iLayer,iElem,gg) /(LimnoParam%sDPhytW(iLayer,iElem,gg)+NearZero)
                !   P/D_ratio_of_Phytoplankton
                LimnoParam%rPDPhytW(iLayer,iElem,gg) = LimnoParam%sPPhytW(iLayer,iElem,gg) /(LimnoParam%sDPhytW(iLayer,iElem,gg)+NearZero)
                !   N/D_ratio_of_Phytoplankton
                LimnoParam%rNDPhytW(iLayer,iElem,gg) = LimnoParam%sNPhytW(iLayer,iElem,gg) /(LimnoParam%sDPhytW(iLayer,iElem,gg)+NearZero)
            
                LimnoParam%rCDPhytW(iLayer,iElem,gg) = Max(LimnoParam%cCDPhytMin(gg),Min(LimnoParam%cCDPhytMax(gg),LimnoParam%rCDPhytW(iLayer,iElem,gg)))
                LimnoParam%rPDPhytW(iLayer,iElem,gg) = Max(LimnoParam%cPDPhytMin(gg),Min(LimnoParam%cPDPhytMax(gg),LimnoParam%rPDPhytW(iLayer,iElem,gg)))
                LimnoParam%rNDPhytW(iLayer,iElem,gg) = Max(LimnoParam%cNDPhytMin(gg),Min(LimnoParam%cNDPhytMax(gg),LimnoParam%rNDPhytW(iLayer,iElem,gg)))
            EndDo
            
            !-----------------------------------------------------------------------
            !   Current Macrophytes ratios(check the curent state)
            !-----------------------------------------------------------------------
            !  Local status
            Do mm = 1, LimnoParam%nmac
                !   C/D_ratio_of_Macrophytes
                LimnoParam%rCDMac(iElem,mm) = LimnoParam%sCMac(iElem,mm) /(LimnoParam%sDMac(iElem,mm)+NearZero)
                !   N/D_ratio_of_Macrophytes
                LimnoParam%rNDMac(iElem,mm) = LimnoParam%sNMac(iElem,mm) /(LimnoParam%sDMac(iElem,mm)+NearZero)
                !   P/D_ratio_of_Macrophytes
                LimnoParam%rPDMac(iElem,mm) = LimnoParam%spMac(iElem,mm) /(LimnoParam%sDMac(iElem,mm)+NearZero)
            EndDo
            
            !-----------------------------------------------------------------------
            !   Current Bacterioplankton in the Water ratios(check the curent state)
            !-----------------------------------------------------------------------
            !  Local status
            Do bb = 1, LimnoParam%nbac
                !   C/D_ratio_of_Bacterioplankton
                LimnoParam%rCDBacW(iLayer,iElem,bb) = LimnoParam%sCBacW(iLayer,iElem,bb) /(LimnoParam%sDBacW(iLayer,iElem,bb)+NearZero)
                !   P/D_ratio_of_Bacterioplankton
                LimnoParam%rPDBacW(iLayer,iElem,bb) = LimnoParam%sPBacW(iLayer,iElem,bb) /(LimnoParam%sDBacW(iLayer,iElem,bb)+NearZero)
                !   N/D_ratio_of_Bacterioplankton
                LimnoParam%rNDBacW(iLayer,iElem,bb) = LimnoParam%sNBacW(iLayer,iElem,bb) /(LimnoParam%sDBacW(iLayer,iElem,bb)+NearZero)
            EndDo

            !-----------------------------------------------------------------------
            !   Current Zooplankton ratios(check the curent state)
            !-----------------------------------------------------------------------
            Do zz = 1, LimnoParam%nzoo
                !  C/D_ratio_herb.zooplankton
                LimnoParam%rCDZoo(iLayer,iElem,zz) = LimnoParam%sCZoo(iLayer,iElem,zz) /(LimnoParam%sDZoo(iLayer,iElem,zz)+NearZero)
                !  P/D_ratio_herb.zooplankton
                LimnoParam%rPDZoo(iLayer,iElem,zz) = LimnoParam%sPZoo(iLayer,iElem,zz) /(LimnoParam%sDZoo(iLayer,iElem,zz)+NearZero)
                !  N/C_ratio_herb.zooplankton
                LimnoParam%rNDZoo(iLayer,iElem,zz) = LimnoParam%sNZoo(iLayer,iElem,zz) /(LimnoParam%sDZoo(iLayer,iElem,zz)+NearZero)
            EndDo

            !-----------------------------------------------------------------------
            !   Current Fish ratios(check the curent state)
            !-----------------------------------------------------------------------
            Do ff = 1, LimnoParam%nfish
                !  C/D_ratio_adult_fish
                LimnoParam%rCDFiAd(iLayer,iElem,ff) = LimnoParam%sCFiAd(iLayer,iElem,ff) /(LimnoParam%sDFiAd(iLayer,iElem,ff)+NearZero)
                !  P/D_ratio_adult_fish
                LimnoParam%rPDFiAd(iLayer,iElem,ff) = LimnoParam%sPFiAd(iLayer,iElem,ff) /(LimnoParam%sDFiAd(iLayer,iElem,ff)+NearZero)
                !  N/C_ratio_adult_fish
                LimnoParam%rNDFiAd(iLayer,iElem,ff) = LimnoParam%sNFiAd(iLayer,iElem,ff) /(LimnoParam%sDFiAd(iLayer,iElem,ff)+NearZero)

                !  C/D_ratio_juvenile_fish
                LimnoParam%rCDFiJv(iLayer,iElem,ff) = LimnoParam%sCFiJv(iLayer,iElem,ff) /(LimnoParam%sDFiJv(iLayer,iElem,ff)+NearZero)
                !  P/D_ratio_juvenile_fish
                LimnoParam%rPDFiJv(iLayer,iElem,ff) = LimnoParam%sPFiJv(iLayer,iElem,ff) /(LimnoParam%sDFiJv(iLayer,iElem,ff)+NearZero)
                !  N/C_ratio_juvenile_fish
                LimnoParam%rNDFiJv(iLayer,iElem,ff) = LimnoParam%sNFiJv(iLayer,iElem,ff) /(LimnoParam%sDFiJv(iLayer,iElem,ff)+NearZero)
            EndDo
            
            
        EndDo
    EndDo

Return
End
