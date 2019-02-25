!> This subroutine reads the simulation parameters. 
Subroutine AllocateWQVars(LimnoParam,MeshParam)
    
    Use LimnologyVars
    Use MeshVars

    
    Implicit none
    type(MeshGridParam) :: MeshParam
    type(LimnologyParam) :: LimnoParam

    ! 1. Water Quality boundary condition variables
    !Allocate(LimnoParam%IndexWQ(MeshParam%nElem,15))
    Allocate(LimnoParam%dVarEstS(MeshParam%nElem,4))
    
    ! 2. Transport Equation Variables
    
    ! 3. Water temperature, Salinity and Density
	Allocate(LimnoParam%sDTempW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sDTempWP(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sDSal(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sDSalP(MeshParam%KMax,MeshParam%nElem))
    
    
!------------------------------------------------------------------------------------------------------------------!	
!  B I O T I C O                                                                                                   !	
!------------------------------------------------------------------------------------------------------------------!	
!------------------------------------------------------------------------------------------------------------------!	
!  B A C T E R I O P L A N C T O N                                                                                 !	
!------------------------------------------------------------------------------------------------------------------!
!1)Configuracoes Iniciais
!2)Parametros	
!3)Condicoes Iniciais 
!4)Fluxos
	Allocate(LimnoParam%wDAssBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wCAssBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wNAssBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wPAssBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wSiAssBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac)) 
	Allocate(LimnoParam%tDAssBac(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tCAssBac(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tNAssBac(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tPAssBac(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wDAssBacDom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac,LimnoParam%ndom))  
	Allocate(LimnoParam%wCAssBacDom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac,LimnoParam%ndom))  
	Allocate(LimnoParam%wNAssBacDom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac,LimnoParam%ndom))  
	Allocate(LimnoParam%wPAssBacDom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac,LimnoParam%ndom))  
	Allocate(LimnoParam%wSiAssBacDom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac,LimnoParam%ndom))  
	Allocate(LimnoParam%tDAssBacDom(MeshParam%nElem,LimnoParam%nbac,LimnoParam%ndom))  
	Allocate(LimnoParam%tCAssBacDom(MeshParam%nElem,LimnoParam%nbac,LimnoParam%ndom))  
	Allocate(LimnoParam%tNAssBacDom(MeshParam%nElem,LimnoParam%nbac,LimnoParam%ndom))  
	Allocate(LimnoParam%tPAssBacDom(MeshParam%nElem,LimnoParam%nbac,LimnoParam%ndom))  
	Allocate(LimnoParam%wCRespBac2CO2(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wCRespBac2CH4(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wDRespBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wCRespBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
    Allocate(LimnoParam%wNExcrBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))
    Allocate(LimnoParam%wPExcrBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))
	Allocate(LimnoParam%tCRespBac2CO2(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tCRespBac2CH4(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tDRespBac(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tCRespBac(MeshParam%nElem,LimnoParam%nbac))  
    Allocate(LimnoParam%tNExcrBac(MeshParam%nElem,LimnoParam%nbac))
    Allocate(LimnoParam%tPExcrBac(MeshParam%nElem,LimnoParam%nbac))
    Allocate(LimnoParam%wPUptBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wNUptBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac)) 
    Allocate(LimnoParam%tPUptBac(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tNUptBac(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wDMortBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wCMortBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wNMortBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%wPMortBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tDMortBac(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tCMortBac(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tNMortBac(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tPMortBac(MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tDSetBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tCSetBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tNSetBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
	Allocate(LimnoParam%tPSetBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  
!5)Funcoes Auxiliares
	!IF (.NOT. AllocateD(aTLimMinBac))   Allocate(LimnoParam%aTLimMinBac(LimnoParam%nbac))  
	!IF (.NOT. AllocateD(apHLimMinBac))  Allocate(LimnoParam%apHLimMinBac(LimnoParam%nbac))  
	!IF (.NOT. AllocateD(aO2LimMinBac))  Allocate(LimnoParam%aO2LimMinBac(LimnoParam%nbac))  
	!IF (.NOT. AllocateD(aNutLimMinBac)) Allocate(LimnoParam%aNutLimMinBac(LimnoParam%nbac))  
	!IF (.NOT. AllocateD(aCorMinBac))    Allocate(LimnoParam%aCorMinBac(LimnoParam%nbac))  	 
	!IF (.NOT. AllocateD(oDFoodBac))     Allocate(LimnoParam%oDFoodBac(LimnoParam%nbac))  	 
	!IF (.NOT. AllocateD(aDFoodSatBac))  Allocate(LimnoParam%aDFoodSatBac(LimnoParam%nbac))  	 
	Allocate(LimnoParam%oDFoodBacDom(LimnoParam%ndom))  	 
	Allocate(LimnoParam%fDAssBacDom(LimnoParam%ndom)) 
    Allocate(LimnoParam%oDSetBacSum(LimnoParam%ndom))
	!IF (.NOT. AllocateD(aDAssBacSub))   Allocate(LimnoParam%aDAssBacSub(MeshParam%KMax,LimnoParam%nbac,LimnoParam%ndom))	 
	!IF (.NOT. AllocateD(kCorResbBac))   Allocate(LimnoParam%kCorResbBac(LimnoParam%nbac))  	 
!	IF (.NOT. AllocateD(kCorMortBac))   Allocate(LimnoParam%kCorMortBac(LimnoParam%nbac))  	 
!7)Variaveis de Estado	
	Allocate(LimnoParam%sDBacW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%sCBacW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%sNBacW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%sPBacW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  	 
	
	Allocate(LimnoParam%sDBacWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%sCBacWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%sNBacWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%sPBacWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  	 

	Allocate(LimnoParam%sDBacS(MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%sCBacS(MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%sNBacS(MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%sPBacS(MeshParam%nElem,LimnoParam%nbac))  	   
	
	Allocate(LimnoParam%sDBacSP(MeshParam%nElem,LimnoParam%nbac))     
	Allocate(LimnoParam%sCBacSP(MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%sNBacSP(MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%sPBacSP(MeshParam%nElem,LimnoParam%nbac))  	 
	
	Allocate(LimnoParam%rCDBacW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%rNDBacW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%rPDBacW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nbac))  	                          
	
	Allocate(LimnoParam%rCDBacS(MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%rNDBacS(MeshParam%nElem,LimnoParam%nbac))  	 
	Allocate(LimnoParam%rPDBacS(MeshParam%nElem,LimnoParam%nbac))  	                           

!------------------------------------------------------------------------------------------------------------------!	
!  F I T O P L A N C T O N                                                                                         !	
!------------------------------------------------------------------------------------------------------------------!
!1)Configuracoes Iniciais
!2)Parametros(pg)	
    !Produção primária
    !Consumo de Nutrientes
    !Respiração
    !Sedimentação
    !Lise
!3)Condicoes Iniciais 
!4)Fluxos - Agua(K,gg),Sed(K,gg)
	!Consumo de nutrientes
	Allocate(LimnoParam%wCUptPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%wNUptPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%wPUptPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%wSiUptPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy)) 
	Allocate(LimnoParam%wNUptNH4Phyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%wNUptNO3Phyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))                          
	!Produção Primária
	Allocate(LimnoParam%wDAssPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!Respiração 
	!IF (.NOT. AllocateD(wDRespbPhyt))  Allocate(LimnoParam%wDRespbPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(wDResppPhyt))  Allocate(LimnoParam%wDResppPhyt(MeshParam%KMax,LimnoParam%nphy))
	Allocate(LimnoParam%wDRespPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
    Allocate(LimnoParam%wCExcrPhytW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
    Allocate(LimnoParam%wPExcrPhytW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
    Allocate(LimnoParam%wNExcrPhytW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tDRespPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tCExcrPhytS(MeshParam%nElem,LimnoParam%nphy))
    Allocate(LimnoParam%tNExcrPhytS(MeshParam%nElem,LimnoParam%nphy))
    Allocate(LimnoParam%tPExcrPhytS(MeshParam%nElem,LimnoParam%nphy))
	!Resuspensão
	Allocate(LimnoParam%tDResusPhyt(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tCResusPhyt(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tPResusPhyt(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tNResusPhyt(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tSiResusPhyt(MeshParam%nElem,LimnoParam%nphy))
	!Sedimentação
    Allocate(LimnoParam%uCorVSetPhyt(LimnoParam%nphy))
    Allocate(LimnoParam%StokesVelPhyt(MeshParam%KMax,LimnoParam%nphy))
    Allocate(LimnoParam%tDSetPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(tCSetPhyt))    Allocate(LimnoParam%tCSetPhyt(MeshParam%KMax,LimnoParam%nphy))
	!IF (.NOT. AllocateD(tNSetPhyt))    Allocate(LimnoParam%tNSetPhyt(MeshParam%KMax,LimnoParam%nphy))
	!IF (.NOT. AllocateD(tPSetPhyt))    Allocate(LimnoParam%tPSetPhyt(MeshParam%KMax,LimnoParam%nphy))
	!IF (.NOT. AllocateD(tSiSetPhyt))   Allocate(LimnoParam%tSiSetPhyt(MeshParam%KMax,LimnoParam%nphy))
	!Lise
	Allocate(LimnoParam%wDLisPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%wCLisPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%wNLisPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%wPLisPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%wSiLisPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
    
	!IF (.NOT. AllocateD(wDLisPomPhyt))  Allocate(LimnoParam%wDLisPomPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(wCLisPomPhyt))  Allocate(LimnoParam%wCLisPomPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(wNLisPomPhyt))  Allocate(LimnoParam%wNLisPomPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(wPLisPomPhyt))  Allocate(LimnoParam%wPLisPomPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(wSiLisPomPhyt)) Allocate(LimnoParam%wSiLisPomPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy)) 
	!IF (.NOT. AllocateD(wDLisDomPhyt))  Allocate(LimnoParam%wDLisDomPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(wCLisDomPhyt))  Allocate(LimnoParam%wCLisDomPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(wNLisDomPhyt))  Allocate(LimnoParam%wNLisDomPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(wPLisDomPhyt))  Allocate(LimnoParam%wPLisDomPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(wSiLisDomPhyt)) Allocate(LimnoParam%wSiLisDomPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tDLisPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tCLisPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tNLisPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tPLisPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tSiLisPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tDLisPomPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tCLisPomPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tNLisPomPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tPLisPomPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tSiLisPomPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tDLisDomPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tCLisDomPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tNLisDomPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tPLisDomPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%tSiLisDomPhytS(MeshParam%nElem,LimnoParam%nphy))
	!Exudacao
	!IF (.NOT. AllocateD(wDExdPhyt))  Allocate(LimnoParam%wDExdPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(wCExdPhyt))  Allocate(LimnoParam%wCExdPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(wNExdPhyt))  Allocate(LimnoParam%wNExdPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	!IF (.NOT. AllocateD(wPExdPhyt))  Allocate(LimnoParam%wPExdPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))	
	!IF (.NOT. AllocateD(wSiExdPhyt)) Allocate(LimnoParam%wSiExdPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
!5)Funcoes Auxiliares(gg)
	!Consumo de nutrientes
	!IF (.NOT. AllocateD(uFunTmPhyt))        Allocate(LimnoParam%uFunTmPhyt(MeshParam%KMax,LimnoParam%nphy))  
	!IF (.NOT. AllocateD(aVCUptMaxCorPhyt))  Allocate(LimnoParam%aVCUptMaxCorPhyt(LimnoParam%nphy))
	!IF (.NOT. AllocateD(aVNUptMaxCorPhyt))  Allocate(LimnoParam%aVNUptMaxCorPhyt(LimnoParam%nphy))
	!IF (.NOT. AllocateD(aVPUptMaxCorPhyt))  Allocate(LimnoParam%aVPUptMaxCorPhyt(LimnoParam%nphy))
	!IF (.NOT. AllocateD(aVSiUptMaxCorPhyt)) Allocate(LimnoParam%aVSiUptMaxCorPhyt(LimnoParam%nphy))
	!IF (.NOT. AllocateD(aVCUptPhyt))        Allocate(LimnoParam%aVCUptPhyt(LimnoParam%nphy))       
	!IF (.NOT. AllocateD(aVNUptPhyt))        Allocate(LimnoParam%aVNUptPhyt(LimnoParam%nphy))       
	!IF (.NOT. AllocateD(aVPUptPhyt))		Allocate(LimnoParam%aVPUptPhyt(LimnoParam%nphy))       
	!IF (.NOT. AllocateD(aVSiUptPhyt))		Allocate(LimnoParam%aVSiUptPhyt(LimnoParam%nphy))
	!IF (.NOT. AllocateD(ahPUptPhyt))		Allocate(LimnoParam%ahPUptPhyt(LimnoParam%nphy))       
	!IF (.NOT. AllocateD(ahNUptPhyt))		Allocate(LimnoParam%ahNUptPhyt(LimnoParam%nphy))      
	!IF (.NOT. AllocateD(afNH4UptPhyt))		Allocate(LimnoParam%afNH4UptPhyt(LimnoParam%nphy))  
	!Produção Primária
	Allocate(LimnoParam%aExtIM(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%aExtPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%aExtDom(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%aExtZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%aExtPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%aExtCoefOpen(MeshParam%KMax,MeshParam%nElem))
	!IF (.NOT. AllocateD(aPACoef))	   Allocate(LimnoParam%aPACoef(MeshParam%KMax))
	Allocate(LimnoParam%aExtCoef(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%aSecchi(MeshParam%KMax,MeshParam%nElem))
	!IF (.NOT. AllocateD(aLLimPhyt))    Allocate(LimnoParam%aLLimPhyt(MeshParam%KMax,LimnoParam%nphy))
	Allocate(LimnoParam%uLPAR0(MeshParam%KMax))
	Allocate(LimnoParam%aLPARBot(MeshParam%KMax))
	!IF (.NOT. AllocateD(uhLPhyt))      Allocate(LimnoParam%uhLPhyt(LimnoParam%nphy))
	!IF (.NOT. AllocateD(aMuTmLPhyt))   Allocate(LimnoParam%aMuTmLPhyt(LimnoParam%nphy))
	!IF (.NOT. AllocateD(aCLimPhyt))    Allocate(LimnoParam%aCLimPhyt(LimnoParam%nphy))
	!IF (.NOT. AllocateD(aPLimPhyt))    Allocate(LimnoParam%aPLimPhyt(LimnoParam%nphy))
	!IF (.NOT. AllocateD(aNLimPhyt))    Allocate(LimnoParam%aNLimPhyt(LimnoParam%nphy))
	!IF (.NOT. AllocateD(aSiLimPhyt))   Allocate(LimnoParam%aSiLimPhyt(LimnoParam%nphy))
	!IF (.NOT. AllocateD(aNutLimPhyt))  Allocate(LimnoParam%aNutLimPhyt(LimnoParam%nphy))
	!IF (.NOT. AllocateD(aMuPhyt))      Allocate(LimnoParam%aMuPhyt(LimnoParam%nphy))
	Allocate(LimnoParam%rChDPhyt(LimnoParam%nphy))
	!Respiração
	!IF (.NOT. AllocateD(ukDRespbTmPhyt)) Allocate(LimnoParam%ukDRespbTmPhyt(LimnoParam%nphy))
	!Lise
	!Allocate(LimnoParam%aCLisPhyt(LimnoParam%nphy))
	!Allocate(LimnoParam%aPLisPhyt(LimnoParam%nphy))
	!Allocate(LimnoParam%aNLisPhyt(LimnoParam%nphy))
	!Allocate(LimnoParam%aSiLisPhyt(LimnoParam%nphy))
	!Allocate(LimnoParam%aCLisPhytS(LimnoParam%nphy))
	!Allocate(LimnoParam%aPLisPhytS(LimnoParam%nphy))
	!Allocate(LimnoParam%aNLisPhytS(LimnoParam%nphy))
	!Allocate(LimnoParam%aSiLisPhytS(LimnoParam%nphy))
	!Allocate(LimnoParam%aNutLisPhyt(LimnoParam%nphy))
	!Allocate(LimnoParam%aLisPomPhyt(LimnoParam%nphy))
	!Allocate(LimnoParam%aNutLisPhytS(LimnoParam%nphy))
	!Allocate(LimnoParam%aLisPomPhytS(LimnoParam%nphy))
!7)Variaveis de Estado	
	Allocate(LimnoParam%sDPhytW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sCPhytW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sNPhytW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sPPhytW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sSiPhytW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	
	Allocate(LimnoParam%sDPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sCPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sNPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sPPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sSiPhytS(MeshParam%nElem,LimnoParam%nphy))

	Allocate(LimnoParam%sDPhytWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sCPhytWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sNPhytWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sPPhytWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sSiPhytWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	
	Allocate(LimnoParam%sDPhytSP(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sCPhytSP(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sNPhytSP(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sPPhytSP(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%sSiPhytSP(MeshParam%nElem,LimnoParam%nphy))
	
	Allocate(LimnoParam%rCDPhytW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%rNDPhytW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%rPDPhytW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%rSiDPhytW(MeshParam%KMax,MeshParam%nElem,LimnoParam%nphy))  
	
	Allocate(LimnoParam%rCDPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%rNDPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%rPDPhytS(MeshParam%nElem,LimnoParam%nphy))
	Allocate(LimnoParam%rSiDPhytS(MeshParam%nElem,LimnoParam%nphy))    

!------------------------------------------------------------------------------------------------------------------!	
!  M A C R O F I T A S                                                                                             !	
!------------------------------------------------------------------------------------------------------------------!	
!1)Configuracoes Iniciais
!2)Parametros(mg)
 	!Produção primária
    !Consumo de Nutrientes
    !Respiração e Excreção
    !Mortalidade
    !Consumo por aves
!3)Condicoes Iniciais (mm)
!4)Fluxos (mm)
	!Consumo de Nutrientes
	Allocate(LimnoParam%tDProdMac(MeshParam%nElem,LimnoParam%nmac))
    Allocate(LimnoParam%tCUptMacW(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNUptMacW(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tPUptMacW(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNUptNH4MacW(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNUptNO3MacW(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tCUptMacS(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNUptMacS(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tPUptMacS(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNUptNH4MacS(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNUptNO3MacS(MeshParam%nElem,LimnoParam%nmac))
    Allocate(LimnoParam%tO2UptNO3MacW(MeshParam%nElem,LimnoParam%nmac))
    Allocate(LimnoParam%tO2RespMacW(MeshParam%nElem,LimnoParam%nmac))
    Allocate(LimnoParam%tO2RespMacS(MeshParam%nElem,LimnoParam%nmac))
    Allocate(LimnoParam%tO2ProdMac(MeshParam%nElem,LimnoParam%nmac))
    Allocate(LimnoParam%tO2ProdMacW(MeshParam%nElem,LimnoParam%nmac))
    Allocate(LimnoParam%tO2ProdMacS(MeshParam%nElem,LimnoParam%nmac))
	!Exudacao
	Allocate(LimnoParam%tCExdMac(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNExdMac(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tPExdMac(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tCExdMacW(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNExdMacW(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tPExdMacW(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tCExdMacS(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNExdMacS(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tPExdMacS(MeshParam%nElem,LimnoParam%nmac))
	!Respiração
    Allocate(LimnoParam%tDRespMac(MeshParam%nElem,LimnoParam%nmac))
    Allocate(LimnoParam%tCRespMac(MeshParam%nElem,LimnoParam%nmac))
    Allocate(LimnoParam%tDRespMacW(MeshParam%nElem,LimnoParam%nmac))
    Allocate(LimnoParam%tDRespMacS(MeshParam%nElem,LimnoParam%nmac))
	!Mortalidade
	Allocate(LimnoParam%tDMortMac(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tCMortMac(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNMortMac(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tPMortMac(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tDMortMacW(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tCMortMacW(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNMortMacW(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tPMortMacW(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tDMortMacS(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tCMortMacS(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNMortMacS(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tPMortMacS(MeshParam%nElem,LimnoParam%nmac))
    !Consumo por aves
	Allocate(LimnoParam%tDGrazMacBird(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tCGrazMacBird(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNGrazMacBird(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tPGrazMacBird(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tDEgesBird(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tCEgesBird(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tNEgesBird(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%tPEgesBird(MeshParam%nElem,LimnoParam%nmac))
!7)Variaveis de Estado (I,J,mm)	
	Allocate(LimnoParam%sDMac(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%sCMac(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%sNMac(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%sPMac(MeshParam%nElem,LimnoParam%nmac))
    
	Allocate(LimnoParam%sDMacP(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%sCMacP(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%sNMacP(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%sPMacP(MeshParam%nElem,LimnoParam%nmac))

    Allocate(LimnoParam%rCDMac(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%rNDMac(MeshParam%nElem,LimnoParam%nmac))
	Allocate(LimnoParam%rPDMac(MeshParam%nElem,LimnoParam%nmac))
		
!------------------------------------------------------------------------------------------------------------------!	
!  Z O O P L A N C T O N                                                                                           !	
!------------------------------------------------------------------------------------------------------------------!	
!1)Configuracoes Iniciais
!2)Parametros(zg)
!3)Condicoes Iniciais(zz)
!4)Fluxos(K,zz) ou (K,zz,presa)
	Allocate(LimnoParam%wDGrzZoo(MeshParam%KMax,LimnoParam%nzoo))
	Allocate(LimnoParam%wCGrzZoo(MeshParam%KMax,LimnoParam%nzoo))
	Allocate(LimnoParam%wNGrzZoo(MeshParam%KMax,LimnoParam%nzoo))
	Allocate(LimnoParam%wPGrzZoo(MeshParam%KMax,LimnoParam%nzoo))
	Allocate(LimnoParam%wDGrzZooZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nzoo))
	Allocate(LimnoParam%wCGrzZooZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nzoo))
	Allocate(LimnoParam%wNGrzZooZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nzoo))
	Allocate(LimnoParam%wPGrzZooZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nzoo))
	Allocate(LimnoParam%wDGrzZooPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nphy))
	Allocate(LimnoParam%wCGrzZooPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nphy))
	Allocate(LimnoParam%wNGrzZooPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nphy))
	Allocate(LimnoParam%wPGrzZooPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nphy))
	Allocate(LimnoParam%wSiGrzZooPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nphy))
	Allocate(LimnoParam%wDGrzZooBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nbac))
	Allocate(LimnoParam%wCGrzZooBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nbac))
	Allocate(LimnoParam%wNGrzZooBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nbac))
	Allocate(LimnoParam%wPGrzZooBac(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%nbac))
	Allocate(LimnoParam%wDGrzZooPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%npom))
	Allocate(LimnoParam%wCGrzZooPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%npom))
	Allocate(LimnoParam%wNGrzZooPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%npom))
	Allocate(LimnoParam%wPGrzZooPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo,LimnoParam%npom))
	Allocate(LimnoParam%wDRespZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wCRespZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wDEgesZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wCEgesZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wNEgesZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wPEgesZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wDMortZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wCMortZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wNMortZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wPMortZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wNExcrZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wPExcrZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wDMessZooOM(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wCMessZooOM(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wNMessZooOM(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wPMessZooOM(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wDAssZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wCAssZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wNAssZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%wPAssZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
!5)Funcoes Auxiliares(zz)	
	Allocate(LimnoParam%aDSatFoodZoo(LimnoParam%nzoo))
	Allocate(LimnoParam%fDGrzZooZoo(LimnoParam%nzoo,LimnoParam%nzoo))
	Allocate(LimnoParam%fDGrzZooPhyt(LimnoParam%nzoo,LimnoParam%nphy))
	Allocate(LimnoParam%fDGrzZooBac(LimnoParam%nzoo,LimnoParam%nbac))
	Allocate(LimnoParam%fDGrzZooPom(LimnoParam%nzoo,LimnoParam%npom))
	Allocate(LimnoParam%akCExcrZoo(LimnoParam%nzoo))
	Allocate(LimnoParam%akNExcrZoo(LimnoParam%nzoo))
	Allocate(LimnoParam%akPExcrZoo(LimnoParam%nzoo))
!7)Variaveis de Estado(I,J,K,zz)	
	Allocate(LimnoParam%sDZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%sCZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%sNZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%sPZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))                        
	
	Allocate(LimnoParam%sDZooP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%sCZooP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%sNZooP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%sPZooP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))                       
	
	Allocate(LimnoParam%rCDZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%rNDZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))
	Allocate(LimnoParam%rPDZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nzoo))

!------------------------------------------------------------------------------------------------------------------!	
!  Z O O B E N T O S                                                                                               !	
!------------------------------------------------------------------------------------------------------------------!	
!1)Configuracoes Iniciais
!2)Parametros(beng)	
!3)Condicoes Iniciais(ben) 
!4)Fluxos(ben)
	Allocate(LimnoParam%tDGrzBent(LimnoParam%nben))
	Allocate(LimnoParam%tCGrzBent(LimnoParam%nben))
	Allocate(LimnoParam%tNGrzBent(LimnoParam%nben))
	Allocate(LimnoParam%tPGrzBent(LimnoParam%nben))
	Allocate(LimnoParam%tDGrzBentPhyt(MeshParam%nElem,LimnoParam%nben,LimnoParam%nphy))
	Allocate(LimnoParam%tCGrzBentPhyt(MeshParam%nElem,LimnoParam%nben,LimnoParam%nphy))
	Allocate(LimnoParam%tNGrzBentPhyt(MeshParam%nElem,LimnoParam%nben,LimnoParam%nphy))
	Allocate(LimnoParam%tPGrzBentPhyt(MeshParam%nElem,LimnoParam%nben,LimnoParam%nphy))
	Allocate(LimnoParam%tSiGrzBentPhyt(MeshParam%nElem,LimnoParam%nben,LimnoParam%nphy))
	Allocate(LimnoParam%tDGrzBentBac(MeshParam%nElem,LimnoParam%nben,LimnoParam%nbac))
	Allocate(LimnoParam%tCGrzBentBac(MeshParam%nElem,LimnoParam%nben,LimnoParam%nbac))
	Allocate(LimnoParam%tNGrzBentBac(MeshParam%nElem,LimnoParam%nben,LimnoParam%nbac))
	Allocate(LimnoParam%tPGrzBentBac(MeshParam%nElem,LimnoParam%nben,LimnoParam%nbac))
	Allocate(LimnoParam%tDGrzBentPom(MeshParam%nElem,LimnoParam%nben,LimnoParam%npom))
	Allocate(LimnoParam%tCGrzBentPom(MeshParam%nElem,LimnoParam%nben,LimnoParam%npom))
	Allocate(LimnoParam%tNGrzBentPom(MeshParam%nElem,LimnoParam%nben,LimnoParam%npom))
	Allocate(LimnoParam%tPGrzBentPom(MeshParam%nElem,LimnoParam%nben,LimnoParam%npom))
	Allocate(LimnoParam%tDRespBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tCRespBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tDEgesBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tCEgesBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tNEgesBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tPEgesBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tDMortBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tCMortBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tNMortBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tPMortBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tNExcBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tPExcBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tDMessBentOM(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tCMessBentOM(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tNMessBentOM(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%tPMessBentOM(MeshParam%nElem,LimnoParam%nben))
!7)Variaveis de Estado(I,J,ben)	
	Allocate(LimnoParam%sDBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%sCBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%sNBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%sPBent(MeshParam%nElem,LimnoParam%nben))            
	
	Allocate(LimnoParam%sDBentP(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%sCBentP(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%sNBentP(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%sPBentP(MeshParam%nElem,LimnoParam%nben))          
	
	Allocate(LimnoParam%rCDBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%rNDBent(MeshParam%nElem,LimnoParam%nben))
	Allocate(LimnoParam%rPDBent(MeshParam%nElem,LimnoParam%nben))                            
	
!------------------------------------------------------------------------------------------------------------------!	
!  P E I X E S                                                                                                     !	
!------------------------------------------------------------------------------------------------------------------!	
!1)Configuracoes Iniciais
!2)Parametros	
!3)Condicoes Iniciais 
!4)Fluxos
    Allocate(LimnoParam%wDMigrFiAd(LimnoParam%nfish))
	Allocate(LimnoParam%wCMigrFiAd(LimnoParam%nfish))
	Allocate(LimnoParam%wNMigrFiAd(LimnoParam%nfish))
	Allocate(LimnoParam%wPMigrFiAd(LimnoParam%nfish))					
    Allocate(LimnoParam%wDMigrFiJv(LimnoParam%nfish))
	Allocate(LimnoParam%wCMigrFiJv(LimnoParam%nfish))
	Allocate(LimnoParam%wNMigrFiJv(LimnoParam%nfish))
	Allocate(LimnoParam%wPMigrFiJv(LimnoParam%nfish))					
    Allocate(LimnoParam%wDReprFish(LimnoParam%nfish))
	Allocate(LimnoParam%wCReprFish(LimnoParam%nfish))   
	Allocate(LimnoParam%wNReprFish(LimnoParam%nfish))
	Allocate(LimnoParam%wPReprFish(LimnoParam%nfish))
	Allocate(LimnoParam%wDAgeFish(LimnoParam%nfish))      
	Allocate(LimnoParam%wCAgeFish(LimnoParam%nfish))      
	Allocate(LimnoParam%wNAgeFish(LimnoParam%nfish))
	Allocate(LimnoParam%wPAgeFish(LimnoParam%nfish)) 
    Allocate(LimnoParam%wDAssFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
    Allocate(LimnoParam%wCAssFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
    Allocate(LimnoParam%wNAssFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
    Allocate(LimnoParam%wPAssFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
    Allocate(LimnoParam%wDAssFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
    Allocate(LimnoParam%wCAssFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
    Allocate(LimnoParam%wNAssFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
    Allocate(LimnoParam%wPAssFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wDGrzFiAdFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nfish))
	Allocate(LimnoParam%wCGrzFiAdFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nfish))
	Allocate(LimnoParam%wNGrzFiAdFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nfish))
	Allocate(LimnoParam%wPGrzFiAdFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nfish))
	Allocate(LimnoParam%wDGrzFiAdFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nfish))
	Allocate(LimnoParam%wCGrzFiAdFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nfish))
	Allocate(LimnoParam%wNGrzFiAdFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nfish))
	Allocate(LimnoParam%wPGrzFiAdFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nfish))
	Allocate(LimnoParam%wDGrzFiAdZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nzoo)) 
	Allocate(LimnoParam%wCGrzFiAdZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nzoo))
	Allocate(LimnoParam%wNGrzFiAdZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nzoo))
	Allocate(LimnoParam%wPGrzFiAdZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nzoo))
	Allocate(LimnoParam%wDGrzFiJvZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nzoo)) 
	Allocate(LimnoParam%wCGrzFiJvZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nzoo))
	Allocate(LimnoParam%wNGrzFiJvZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nzoo))
	Allocate(LimnoParam%wPGrzFiJvZoo(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nzoo))
	Allocate(LimnoParam%wDGrzFiAdBent(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nben)) 
	Allocate(LimnoParam%wCGrzFiAdBent(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nben))
	Allocate(LimnoParam%wNGrzFiAdBent(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nben))
	Allocate(LimnoParam%wPGrzFiAdBent(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nben))
	Allocate(LimnoParam%wDGrzFiJvBent(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nben)) 
	Allocate(LimnoParam%wCGrzFiJvBent(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nben))
	Allocate(LimnoParam%wNGrzFiJvBent(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nben))
	Allocate(LimnoParam%wPGrzFiJvBent(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nben))
    Allocate(LimnoParam%wDGrzFiAdPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nphy))
	Allocate(LimnoParam%wCGrzFiAdPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nphy))
	Allocate(LimnoParam%wNGrzFiAdPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nphy))
	Allocate(LimnoParam%wPGrzFiAdPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nphy))
	Allocate(LimnoParam%wDGrzFiJvPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nphy))
	Allocate(LimnoParam%wCGrzFiJvPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nphy))
	Allocate(LimnoParam%wNGrzFiJvPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nphy))
	Allocate(LimnoParam%wPGrzFiJvPhyt(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%nphy))
	!IF (.NOT. AllocateD(tDGrzFishBac))  Allocate(LimnoParam%tDGrzFishBac(LimnoParam%nfish,LimnoParam%nbac)) 
	!IF (.NOT. AllocateD(tCGrzFishBac))  Allocate(LimnoParam%tCGrzFishBac(LimnoParam%nfish,LimnoParam%nbac))
	!IF (.NOT. AllocateD(tNGrzFishBac))  Allocate(LimnoParam%tNGrzFishBac(LimnoParam%nfish,LimnoParam%nbac))
	!IF (.NOT. AllocateD(tPGrzFishBac))  Allocate(LimnoParam%tPGrzFishBac(LimnoParam%nfish,LimnoParam%nbac))
	!IF (.NOT. AllocateD(tDGrzFishBent)) Allocate(LimnoParam%tDGrzFishBent(LimnoParam%nfish,LimnoParam%nben))
	!IF (.NOT. AllocateD(tCGrzFishBent)) Allocate(LimnoParam%tCGrzFishBent(LimnoParam%nfish,LimnoParam%nben))
	!IF (.NOT. AllocateD(tNGrzFishBent)) Allocate(LimnoParam%tNGrzFishBent(LimnoParam%nfish,LimnoParam%nben))
	!IF (.NOT. AllocateD(tPGrzFishBent)) Allocate(LimnoParam%tPGrzFishBent(LimnoParam%nfish,LimnoParam%nben))
	Allocate(LimnoParam%wDGrzFiAdPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%npom)) 
	Allocate(LimnoParam%wCGrzFiAdPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%npom)) 
	Allocate(LimnoParam%wNGrzFiAdPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%npom))
	Allocate(LimnoParam%wPGrzFiAdPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%npom))
	Allocate(LimnoParam%wDGrzFiJvPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%npom)) 
	Allocate(LimnoParam%wCGrzFiJvPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%npom)) 
	Allocate(LimnoParam%wNGrzFiJvPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%npom))
	Allocate(LimnoParam%wPGrzFiJvPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish,LimnoParam%npom))
	Allocate(LimnoParam%wDHarvFiAd(LimnoParam%nfish))     
	Allocate(LimnoParam%wCHarvFiAd(LimnoParam%nfish))     
	Allocate(LimnoParam%wNHarvFiAd(LimnoParam%nfish))
	Allocate(LimnoParam%wPHarvFiAd(LimnoParam%nfish))
	Allocate(LimnoParam%wDHarvFiJv(LimnoParam%nfish))     
	Allocate(LimnoParam%wCHarvFiJv(LimnoParam%nfish))     
	Allocate(LimnoParam%wNHarvFiJv(LimnoParam%nfish))
	Allocate(LimnoParam%wPHarvFiJv(LimnoParam%nfish))
	Allocate(LimnoParam%wDRespFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
    Allocate(LimnoParam%wDRespFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wCRespFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
    Allocate(LimnoParam%wCRespFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wDEgesFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))     
	Allocate(LimnoParam%wCEgesFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))     
	Allocate(LimnoParam%wNEgesFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wPEgesFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wDEgesFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))     
	Allocate(LimnoParam%wCEgesFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))     
	Allocate(LimnoParam%wNEgesFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wPEgesFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wDMortFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))     
	Allocate(LimnoParam%wCMortFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))     
	Allocate(LimnoParam%wNMortFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wPMortFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wDMortFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))     
	Allocate(LimnoParam%wCMortFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))     
	Allocate(LimnoParam%wNMortFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wPMortFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))

	Allocate(LimnoParam%wNExcrFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wPExcrFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wNExcrFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wPExcrFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))

    Allocate(LimnoParam%wDMessFiAdOM(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wCMessFiAdOM(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wNMessFiAdOM(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wPMessFiAdOM(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
    Allocate(LimnoParam%wDMessFiJvOM(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wCMessFiJvOM(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wNMessFiJvOM(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%wPMessFiJvOM(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))

    !7)Variaveis de Estado	
	Allocate(LimnoParam%sDFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%sCFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%sNFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%sPFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))        

    Allocate(LimnoParam%sDFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%sCFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%sNFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%sPFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))        

	Allocate(LimnoParam%sDFiAdP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%sCFiAdP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%sNFiAdP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%sPFiAdP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))       

	Allocate(LimnoParam%sDFiJvP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%sCFiJvP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%sNFiJvP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%sPFiJvP(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))       
    
	Allocate(LimnoParam%rCDFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%rNDFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%rPDFiAd(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))

	Allocate(LimnoParam%rCDFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%rNDFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
	Allocate(LimnoParam%rPDFiJv(MeshParam%KMax,MeshParam%nElem,LimnoParam%nfish))
    
!------------------------------------------------------------------------------------------------------------------!	
!  M O D U L O    E C O L O G I C O  =>  B I O T I C O / E C O L O G I C O                                         !	
!------------------------------------------------------------------------------------------------------------------!	
!------------------------------------------------------------------------------------------------------------------!	
!  A B I O T I C O                                                                                                 !	
!------------------------------------------------------------------------------------------------------------------!	
!------------------------------------------------------------------------------------------------------------------!	
!  C O N F I G U R A C O E S     I N I C I A I S                                                                   !	
!------------------------------------------------------------------------------------------------------------------!	
!------------------------------------------------------------------------------------------------------------------!	
!  P A R A M E T R O S                                                                                             !	
!------------------------------------------------------------------------------------------------------------------!	
!Link BIOTICO<->ABIOTICO
	!Particulate Organic Matter (DOM)
	!Dissolved Organic Matter (DOM)
	!Climaticos
    !Sedimento
    !Infiltracao
	!Erosao
    !Resuspensão
    !Sedimentacao
    !Hidrolise - POM -> DOM
    !Mineração - IF InclBac=0 - DOM->CO2,NH4,PO4,SiO2
    !DBO e Coliformes
    !Desnitrificação
    !Nitrificação
    !Adsorção do P
    !Imobilizacao do P
    !Difusao Agua-Sed
    !Difusao Agua-Ar
	!O2
    !CO2
	!Geral
!------------------------------------------------------------------------------------------------------------------!	
! F L U X O S   E   V A R S   A U X I L I A R E S                                                                  !	
!------------------------------------------------------------------------------------------------------------------!	
	!Sedimento
	Allocate(LimnoParam%aDepthS(LimnoParam%nsed))
	Allocate(LimnoParam%fDepthS(LimnoParam%nsed))
	!Infiltração
	Allocate(LimnoParam%tDInfDomW(LimnoParam%ndom))
	Allocate(LimnoParam%tCInfDomW(LimnoParam%ndom))
	Allocate(LimnoParam%tNInfDomW(LimnoParam%ndom))
	Allocate(LimnoParam%tPInfDomW(LimnoParam%ndom))
	Allocate(LimnoParam%tDInfDomS(LimnoParam%ndom))
	Allocate(LimnoParam%tCInfDomS(LimnoParam%ndom))
	Allocate(LimnoParam%tNInfDomS(LimnoParam%ndom))
	Allocate(LimnoParam%tPInfDomS(LimnoParam%ndom))
    !Erosion
    Allocate(LimnoParam%uDErosIMS(MeshParam%nElem))
	Allocate(LimnoParam%tDErosPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%tCErosPom(MeshParam%KMax,LimnoParam%npom))
	Allocate(LimnoParam%tNErosPom(MeshParam%KMax,LimnoParam%npom))
	Allocate(LimnoParam%tPErosPom(MeshParam%KMax,LimnoParam%npom))
	!Enterro
	Allocate(LimnoParam%tDPomS(LimnoParam%npom))
	Allocate(LimnoParam%tDBurPom(LimnoParam%npom))
	Allocate(LimnoParam%tCBurPom(LimnoParam%npom))
	Allocate(LimnoParam%tNBurPom(LimnoParam%npom))
	Allocate(LimnoParam%tPBurPom(LimnoParam%npom))
	Allocate(LimnoParam%tDBurDom(LimnoParam%ndom))
	Allocate(LimnoParam%tCBurDom(LimnoParam%ndom))
	Allocate(LimnoParam%tNBurDom(LimnoParam%ndom))
	Allocate(LimnoParam%tPBurDom(LimnoParam%ndom))
	!Resuspensao
    Allocate(LimnoParam%tDResusIM(MeshParam%nElem))
	Allocate(LimnoParam%tDResusTot(MeshParam%nElem))
	Allocate(LimnoParam%tDResusPom(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%tCResusPom(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%tPResusPom(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%tNResusPom(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%tDResusDom(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%tCResusDom(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%tNResusDom(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%tPResusDom(MeshParam%nElem,LimnoParam%ndom))
    Allocate(LimnoParam%tPResusPO4(MeshParam%nElem))
    Allocate(LimnoParam%tPResusAIM(MeshParam%nElem))
    Allocate(LimnoParam%tNResusNH4(MeshParam%nElem))
    Allocate(LimnoParam%tNResusNO3(MeshParam%nElem))
    Allocate(LimnoParam%tCResusDic(LimnoParam%ndom))
	!Sedimentação
	Allocate(LimnoParam%StokesVelIM(MeshParam%KMax))
    Allocate(LimnoParam%StokesVelPom(MeshParam%KMax,LimnoParam%npom))
	Allocate(LimnoParam%tDSetIM(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%tPSetAIM(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%tDSetPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%tCSetPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%tNSetPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%tPSetPom(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	!Hidrolise
	Allocate(LimnoParam%aTLimHidPom(LimnoParam%npom))
	Allocate(LimnoParam%aO2LimHidPom(LimnoParam%npom))
	Allocate(LimnoParam%aBacLimHidPom(LimnoParam%npom))
	Allocate(LimnoParam%apHLimHidPom(LimnoParam%npom))
	Allocate(LimnoParam%kCorHidPom(LimnoParam%npom))	
	Allocate(LimnoParam%wDHidPomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%wCHidPomW(MeshParam%KMax,LimnoParam%npom))
	Allocate(LimnoParam%wNHidPomW(MeshParam%KMax,LimnoParam%npom))
	Allocate(LimnoParam%wPHidPomW(MeshParam%KMax,LimnoParam%npom))
	Allocate(LimnoParam%aTLimHidPomS(LimnoParam%npom))
    Allocate(LimnoParam%tDHidPomS(LimnoParam%nsed,LimnoParam%npom))
    Allocate(LimnoParam%tCHidPomS(LimnoParam%nsed,LimnoParam%npom))
    Allocate(LimnoParam%tNHidPomS(LimnoParam%nsed,LimnoParam%npom))
    Allocate(LimnoParam%tPHidPomS(LimnoParam%nsed,LimnoParam%npom))
    
    
	!Mineração - IF InclBac=0
	Allocate(LimnoParam%wDMinAerW(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%wCMinAerW(MeshParam%KMax,LimnoParam%ndom))
	Allocate(LimnoParam%wNMinAerW(MeshParam%KMax,LimnoParam%ndom))
	Allocate(LimnoParam%wPMinAerW(MeshParam%KMax,LimnoParam%ndom))
	Allocate(LimnoParam%wSiMinAerW(MeshParam%KMax,LimnoParam%ndom))
	Allocate(LimnoParam%wCMinAer2CO2W(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%wCMinAer2CH4W(MeshParam%KMax,LimnoParam%ndom))
	Allocate(LimnoParam%wDMinAnaerW(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%wCMinAnaerW(MeshParam%KMax,LimnoParam%ndom))
	Allocate(LimnoParam%wNMinAnaerW(MeshParam%KMax,LimnoParam%ndom))
	Allocate(LimnoParam%wPMinAnaerW(MeshParam%KMax,LimnoParam%ndom))
	Allocate(LimnoParam%wSiMinAnaerW(MeshParam%KMax,LimnoParam%ndom))
	Allocate(LimnoParam%wCMinAnaer2CO2W(MeshParam%KMax,LimnoParam%ndom))
	Allocate(LimnoParam%wCMinAnaer2CH4W(MeshParam%KMax,LimnoParam%ndom))
	Allocate(LimnoParam%aTLimMinAer(LimnoParam%ndom))
	Allocate(LimnoParam%aTLimMinAnaer(LimnoParam%ndom))
	Allocate(LimnoParam%aO2LimMinAer(LimnoParam%ndom))
	Allocate(LimnoParam%aO2LimMinAnaer(LimnoParam%ndom))
	Allocate(LimnoParam%apHLimMinAer(LimnoParam%ndom))
	Allocate(LimnoParam%apHLimMinAnaer(LimnoParam%ndom))
	Allocate(LimnoParam%aBacLimMinAer(LimnoParam%ndom))
	Allocate(LimnoParam%aBacLimMinAnaer(LimnoParam%ndom))
	Allocate(LimnoParam%aTLimMinDomS(LimnoParam%ndom))
	Allocate(LimnoParam%tDMinAerS(LimnoParam%nsed,LimnoParam%ndom))
	Allocate(LimnoParam%tCMinAerS(LimnoParam%nsed,LimnoParam%ndom))
	Allocate(LimnoParam%tNMinAerS(LimnoParam%nsed,LimnoParam%ndom))
	Allocate(LimnoParam%tPMinAerS(LimnoParam%nsed,LimnoParam%ndom))      
	Allocate(LimnoParam%tCMinAer2CO2S(LimnoParam%nsed,LimnoParam%ndom))
    Allocate(LimnoParam%tCMinAer2CH4S(LimnoParam%nsed,LimnoParam%ndom))
	Allocate(LimnoParam%tDMinAnaerS(LimnoParam%nsed,LimnoParam%ndom))
	Allocate(LimnoParam%tCMinAnaerS(LimnoParam%nsed,LimnoParam%ndom))
	Allocate(LimnoParam%tNMinAnaerS(LimnoParam%nsed,LimnoParam%ndom))
	Allocate(LimnoParam%tPMinAnaerS(LimnoParam%nsed,LimnoParam%ndom))      
	Allocate(LimnoParam%tCMinAnaer2CO2S(LimnoParam%nsed,LimnoParam%ndom))
    Allocate(LimnoParam%tCMinAnaer2CH4S(LimnoParam%nsed,LimnoParam%ndom))

    
	!Desnitrificação
	Allocate(LimnoParam%wNDenitW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%tNDenitS(LimnoParam%nsed))   
	!Nitrificação
	Allocate(LimnoParam%wNNitrW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%tNNitrS(LimnoParam%nsed))
    Allocate(LimnoParam%tO2NitrS(LimnoParam%nsed))
	!Adsorção do Fósforo
	Allocate(LimnoParam%wPSorpIMW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%tPSorpIMS(LimnoParam%nsed))
	!Imobilização do Fósforo
	!Difusao Agua-Sed
	Allocate(LimnoParam%tDDifDom(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%tCDifDom(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%tNDifDom(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%tPDifDom(MeshParam%nElem,LimnoParam%ndom))
    Allocate(LimnoParam%tPDifPO4(MeshParam%nElem))
    Allocate(LimnoParam%tNDifNH4(MeshParam%nElem))
    Allocate(LimnoParam%tNDifNO3(MeshParam%nElem))
    Allocate(LimnoParam%tCDifDic(MeshParam%nElem))
    Allocate(LimnoParam%tO2Dif(MeshParam%nElem))
	!Difusao Agua-Ar
	!O2
	Allocate(LimnoParam%tO2Reaer(MeshParam%KMax))
	!CO2
	Allocate(LimnoParam%tCO2Reaer(MeshParam%KMax))
    
!------------------------------------------------------------------------------------------------------------------!	
!  V A R I A V E I S    D E    E S T A D O                                                                         !	
!------------------------------------------------------------------------------------------------------------------!
!Variáveis de Estado (No tempo ATUAL)
	!Materia Organica
	Allocate(LimnoParam%sDPomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sCPomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sNPomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sPPomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sSiPomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	
	Allocate(LimnoParam%sDDomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sCDomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sNDomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sPDomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sSiDomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	
	Allocate(LimnoParam%sDPomS(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sCPomS(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sNPomS(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sPPomS(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sSiPomS(MeshParam%nElem,LimnoParam%npom))
	
	Allocate(LimnoParam%sDDomS(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sCDomS(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sNDomS(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sPDomS(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sSiDomS(MeshParam%nElem,LimnoParam%ndom))
    
	Allocate(LimnoParam%sDBOW(MeshParam%KMax,MeshParam%nElem))
	!Bioticos
    Allocate(LimnoParam%sCFW(MeshParam%KMax,MeshParam%nElem))
	!Inorganicos
	Allocate(LimnoParam%sDIMW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sPO4W(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sPAIMW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sNH4W(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sNO3W(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sSiO2W(MeshParam%KMax,MeshParam%nElem))

 	Allocate(LimnoParam%sDIMS(MeshParam%nElem))
	Allocate(LimnoParam%sPO4S(MeshParam%nElem))
	Allocate(LimnoParam%sPAIMS(MeshParam%nElem))
	Allocate(LimnoParam%sNH4S(MeshParam%nElem))
	Allocate(LimnoParam%sNO3S(MeshParam%nElem))                  
	!Gases
	Allocate(LimnoParam%sO2W(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sDicW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sCH4W(MeshParam%KMax,MeshParam%nElem))
 	Allocate(LimnoParam%sDicS(MeshParam%nElem))
	Allocate(LimnoParam%sCH4S(MeshParam%nElem))

!Variáveis de Estado (No tempo PASSADO)
	!Materia Organica
	Allocate(LimnoParam%sDPomWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sCPomWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sNPomWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sPPomWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sSiPomWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	
	Allocate(LimnoParam%sDDomWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sCDomWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sNDomWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sPDomWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sSiDomWP(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	
	Allocate(LimnoParam%sDPomSP(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sCPomSP(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sNPomSP(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sPPomSP(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%sSiPomSP(MeshParam%nElem,LimnoParam%npom))
	
	Allocate(LimnoParam%sDDomSP(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sCDomSP(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sNDomSP(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sPDomSP(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%sSiDomSP(MeshParam%nElem,LimnoParam%ndom))
    
	Allocate(LimnoParam%sDBOWP(MeshParam%KMax,MeshParam%nElem))
    Allocate(LimnoParam%sCFWP(MeshParam%KMax,MeshParam%nElem))
    !Inorganicos
	Allocate(LimnoParam%sDIMWP(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sPO4WP(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sPAIMWP(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sNH4WP(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sNO3WP(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sSiO2WP(MeshParam%KMax,MeshParam%nElem))

 	Allocate(LimnoParam%sDIMSP(MeshParam%nElem))
	Allocate(LimnoParam%sPO4SP(MeshParam%nElem))
	Allocate(LimnoParam%sPAIMSP(MeshParam%nElem))
	Allocate(LimnoParam%sNH4SP(MeshParam%nElem))
	Allocate(LimnoParam%sNO3SP(MeshParam%nElem))          
	!Gases
	Allocate(LimnoParam%sO2WP(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sDicWP(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sCH4WP(MeshParam%KMax,MeshParam%nElem))

	Allocate(LimnoParam%sDicSP(MeshParam%nElem))
	Allocate(LimnoParam%sCH4SP(MeshParam%nElem))

!Razoes de nutrientes internos Abióticas e Bioticas
	Allocate(LimnoParam%rCDPomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%rNDPomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%rPDPomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%rSiDPomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%npom))
	
	Allocate(LimnoParam%rCDDomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%rNDDomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%rPDDomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%rSiDDomW(MeshParam%KMax,MeshParam%nElem,LimnoParam%ndom))       
	
	Allocate(LimnoParam%rCDPomS(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%rNDPomS(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%rPDPomS(MeshParam%nElem,LimnoParam%npom))
	Allocate(LimnoParam%rSiDPomS(MeshParam%nElem,LimnoParam%npom))
	
	Allocate(LimnoParam%rCDDomS(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%rNDDomS(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%rPDDomS(MeshParam%nElem,LimnoParam%ndom))
	Allocate(LimnoParam%rSiDDomS(MeshParam%nElem,LimnoParam%ndom))        
	
	Allocate(LimnoParam%rPDIMW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%rPDIMS(MeshParam%nElem))

	!Variaveis Secundarias	
	Allocate(LimnoParam%aDTotS(MeshParam%nElem))  
	Allocate(LimnoParam%aRhoTotS(MeshParam%nElem))   
	Allocate(LimnoParam%aRhoSolidS(MeshParam%nElem))
	Allocate(LimnoParam%afDTotS(MeshParam%nElem))    
	Allocate(LimnoParam%afDOrgS(MeshParam%nElem))
	Allocate(LimnoParam%afPomS(MeshParam%nElem,LimnoParam%npom))  
	Allocate(LimnoParam%afPomTotS(MeshParam%nElem)) 
	
	Allocate(LimnoParam%oDPhytW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%oCPhytW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%oNPhytW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%oPPhytW(MeshParam%KMax,MeshParam%nElem))
	
	Allocate(LimnoParam%aDPhytS(MeshParam%nElem))
	Allocate(LimnoParam%aCPhytS(MeshParam%nElem))
	Allocate(LimnoParam%aNPhytS(MeshParam%nElem))
	Allocate(LimnoParam%aPPhytS(MeshParam%nElem))

	Allocate(LimnoParam%oDSestW(MeshParam%KMax,MeshParam%nElem))                                                                   
	Allocate(LimnoParam%oCSestW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%oNSestW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%oPSestW(MeshParam%KMax,MeshParam%nElem))
	
	Allocate(LimnoParam%oDOMW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%oCOMW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%oNOMW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%oPOMW(MeshParam%KMax,MeshParam%nElem))
	
	Allocate(LimnoParam%oNTotW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%oNkjW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%oNDissW(MeshParam%KMax,MeshParam%nElem))                      
	Allocate(LimnoParam%aNDissS(MeshParam%nElem)) 
	Allocate(LimnoParam%aNkjAvailS(MeshParam%nElem)) 
	Allocate(LimnoParam%aNkjS(MeshParam%nElem))     
	Allocate(LimnoParam%aNTotAvailS(MeshParam%nElem))
	Allocate(LimnoParam%aNTotS(MeshParam%nElem)) 
	Allocate(LimnoParam%afNInorS(MeshParam%nElem))
	Allocate(LimnoParam%afNTotS(MeshParam%nElem)) 
	Allocate(LimnoParam%oNO3S(MeshParam%nElem))
	Allocate(LimnoParam%oNH4S(MeshParam%nElem))
	Allocate(LimnoParam%oNDissS(MeshParam%nElem)) 
	
	Allocate(LimnoParam%oPTotW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%oPInorgW(MeshParam%KMax,MeshParam%nElem))                                 
	Allocate(LimnoParam%aPInorgS(MeshParam%nElem))
	Allocate(LimnoParam%aPTotAvailS(MeshParam%nElem))
	Allocate(LimnoParam%aPTotS(MeshParam%nElem))    
	Allocate(LimnoParam%afPInorgS(MeshParam%nElem))  
	Allocate(LimnoParam%afPO4S(MeshParam%nElem))
	Allocate(LimnoParam%oPO4S(MeshParam%nElem))                                                    
	
	Allocate(LimnoParam%oDicS(MeshParam%nElem))
	Allocate(LimnoParam%oH2CO3S(MeshParam%nElem))
    Allocate(LimnoParam%oHCO3S(MeshParam%nElem))
	Allocate(LimnoParam%oCO3S(MeshParam%nElem))
	Allocate(LimnoParam%oDDomS(MeshParam%nElem,LimnoParam%ndom))
	
 	Allocate(LimnoParam%sAlkW(MeshParam%KMax,MeshParam%nElem))
    Allocate(LimnoParam%spHW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sH2CO3W(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sHCO3W(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sCO3W(MeshParam%KMax,MeshParam%nElem))
 	Allocate(LimnoParam%pCO2W(MeshParam%KMax,MeshParam%nElem))
 	Allocate(LimnoParam%sAlkS(MeshParam%nElem))
 	Allocate(LimnoParam%spHS(MeshParam%nElem))
	Allocate(LimnoParam%sH2CO3S(MeshParam%nElem))
	Allocate(LimnoParam%sHCO3S(MeshParam%nElem))
	Allocate(LimnoParam%sCO3S(MeshParam%nElem))
    
    Allocate(LimnoParam%sO2S(LimnoParam%nsed)) 
    Allocate(LimnoParam%sGPPW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sRespW(MeshParam%KMax,MeshParam%nElem))
	Allocate(LimnoParam%sNEPW(MeshParam%KMax,MeshParam%nElem))    
    
End Subroutine AllocateWQVars
    