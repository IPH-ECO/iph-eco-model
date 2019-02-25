Subroutine RatiosSediment(MeshParam,LimnoParam)
    Use MeshVars
    Use LimnologyVars

    Implicit None
    type(MeshGridParam) :: MeshParam
    type(LimnologyParam) :: LimnoParam
    Integer:: iElem,pom,dom,gg,ben
    Real:: NearZero = 1e-10
    
    Do iElem = 1,MeshParam%nElem
        !Organic Matter in the sediment
        Do pom = 1, LimnoParam%npom
            LimnoParam%rCDPomS(iElem,pom) = LimnoParam%sCPomS(iElem,pom) /(LimnoParam%sDPomS(iElem,pom)+NearZero)
            LimnoParam%rNDPomS(iElem,pom) = LimnoParam%sNPomS(iElem,pom) /(LimnoParam%sDPomS(iElem,pom)+NearZero)
            LimnoParam%rPDPomS(iElem,pom) = LimnoParam%sPPomS(iElem,pom) /(LimnoParam%sDPomS(iElem,pom)+NearZero)
            !LimnoParam%rSiDPomS(iElem,pom) = LimnoParam%sSiPomS(iElem,pom) /(LimnoParam%sDPomS(iElem,pom)+NearZero)
        EndDo
        Do dom = 1, LimnoParam%ndom
            LimnoParam%rCDDomS(iElem,dom) = LimnoParam%sCDomS(iElem,dom) /(LimnoParam%sDDomS(iElem,dom)+NearZero)
            LimnoParam%rNDDomS(iElem,dom) = LimnoParam%sNDomS(iElem,dom) /(LimnoParam%sDDomS(iElem,dom)+NearZero)
            LimnoParam%rPDDomS(iElem,dom:) = LimnoParam%sPDomS(iElem,dom) /(LimnoParam%sDDomS(iElem,dom)+NearZero)
            !LimnoParam%rSiDDomS(iElem,dom) = LimnoParam%sSiDomS(iElem,dom) /(LimnoParam%sDDomS(iElem,dom)+NearZero)
        EndDo
        !Phytoplankton in the sediment
        Do gg = 1, LimnoParam%nphy
            !  P/D_ratio_of_Phytoplankton
            LimnoParam%rCDPhytS(iElem,gg) = Min(LimnoParam%cCDPhytMin(gg),Max(LimnoParam%cCDPhytMax(gg),LimnoParam%sCPhytS(iElem,gg) /(LimnoParam%sDPhytS(iElem,gg)+NearZero)))
            !  N/D_ratio_of_Phytoplankton
            LimnoParam%rNDPhytS(iElem,gg) = Min(LimnoParam%cNDPhytMin(gg),Max(LimnoParam%cNDPhytMax(gg),LimnoParam%sNPhytS(iElem,gg) /(LimnoParam%sDPhytS(iElem,gg)+NearZero)))
            !  P/D_ratio_of_Phytoplankton
            LimnoParam%rPDPhytS(iElem,gg) = Min(LimnoParam%cPDPhytMin(gg),Max(LimnoParam%cPDPhytMax(gg),LimnoParam%sPPhytS(iElem,gg) /(LimnoParam%sDPhytS(iElem,gg)+NearZero)))
        EndDo
        !----------------------------------------------------------------------
        !  Current local nutrients ratios in zoobenthos(check the curent state)
        !----------------------------------------------------------------------
        Do ben = 1, LimnoParam%nben
            LimnoParam%rCDBent(iElem,ben)=LimnoParam%sCBent(iElem,ben)/(LimnoParam%sDBent(iElem,ben)+NearZero)
            LimnoParam%rNDBent(iElem,ben)=LimnoParam%sNBent(iElem,ben)/(LimnoParam%sDBent(iElem,ben)+NearZero)
            LimnoParam%rPDBent(iElem,ben)=LimnoParam%sPBent(iElem,ben)/(LimnoParam%sDBent(iElem,ben)+NearZero)
        EndDo
        
    EndDo !Loop Cell
      

Return
End
