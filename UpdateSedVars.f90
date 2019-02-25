Subroutine UpdateSedVars(HydroParam,MeshParam,LimnoParam)
    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars

    Implicit None
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam
    type(LimnologyParam) :: LimnoParam
    Integer:: iElem,pom

    !Updating state variables in the sediment
    LimnoParam%sDPomSP = LimnoParam%sDPomS
    LimnoParam%sCPomSP = LimnoParam%sCPomS
    LimnoParam%sNPomSP = LimnoParam%sNPomS
    LimnoParam%sPPomSP = LimnoParam%sPPomS
    LimnoParam%sDDomSP = LimnoParam%sDDomS
    LimnoParam%sCDomSP = LimnoParam%sCDomS
    LimnoParam%sNDomSP = LimnoParam%sNDomS
    LimnoParam%sPDomSP = LimnoParam%sPDomS
    LimnoParam%sPO4SP = LimnoParam%sPO4S
    LimnoParam%sPAIMSP = LimnoParam%sPAIMS
    LimnoParam%sNH4SP = LimnoParam%sNH4S
    LimnoParam%sNO3SP = LimnoParam%sNO3S
    LimnoParam%sDIMSP = LimnoParam%sDIMS
    LimnoParam%sDicSP = LimnoParam%sDicS
    LimnoParam%sDBacSP = LimnoParam%sDBacS
    LimnoParam%sCBacSP = LimnoParam%sCBacS
    LimnoParam%sNBacSP = LimnoParam%sNBacS
    LimnoParam%sPBacSP = LimnoParam%sPBacS
    LimnoParam%sDPhytSP = LimnoParam%sDPhytS
    LimnoParam%sCPhytSP = LimnoParam%sCPhytS
    LimnoParam%sNPhytSP = LimnoParam%sNPhytS
    LimnoParam%sPPhytSP = LimnoParam%sPPhytS
    LimnoParam%sDBentP = LimnoParam%sDBent
    LimnoParam%sCBentP = LimnoParam%sCBent
    LimnoParam%sNBentP = LimnoParam%sNBent
    LimnoParam%sPBentP = LimnoParam%sPBent
    
    Do iElem = 1,MeshParam%nElem
	    LimnoParam%aDPhytS(iElem)=SUM(LimnoParam%sDPhytS(iElem,:))	                                    !Fitoplanton no sedimento (gP/m�)
        LimnoParam%aCPhytS(iElem)=SUM(LimnoParam%sCPhytS(iElem,:))	                
        LimnoParam%aNPhytS(iElem)=SUM(LimnoParam%sNPhytS(iElem,:))	                                    !Fitoplanton no sedimento (gP/m�)
        LimnoParam%aPPhytS(iElem)=SUM(LimnoParam%sPPhytS(iElem,:))	                                    !Fitoplancton total no sedimento (gN/m�)
        
		LimnoParam%aDTotS(iElem)=LimnoParam%sDIMS(iElem)+SUM(LimnoParam%sDPoms(iElem,:)) 	                                        !sedimento total (excluinda a biota) (g/m�)
		LimnoParam%aRhoTotS(iElem)=LimnoParam%aDTotS(iElem)/LimnoParam%cDepthS			                                        !densidade aparente do sedimento (g solido/m� sedimento)
		LimnoParam%aRhoSolidS(iElem)=(LimnoParam%sDIMS(iElem)*LimnoParam%cRhoIM+(SUM(LimnoParam%sDPoms(iElem,:)))*LimnoParam%cRhoOM)/LimnoParam%aDTotS(iElem)	        !densidade do solido m�dia (g/m�)
		LimnoParam%afDTotS(iElem)=1.0/(1.0+LimnoParam%bPorS/(1.0-LimnoParam%bPorS)*HydroParam%sDRhoW(HydroParam%ElSmallm(iElem),iElem)/LimnoParam%aRhoSolidS(iElem))	                !fra��o de peso seco do sedimento (g solido/g de sedimento)
		LimnoParam%afDOrgS(iElem)=(SUM(LimnoParam%sDPoms(iElem,:)))/LimnoParam%aDTotS(iElem)	                                    !fra��o de mat�ria organica total do sedimento (-)
		LimnoParam%afPomTotS(iElem)=SUM(LimnoParam%sDPomS(iElem,:))/(LimnoParam%sDIMS(iElem)+SUM(LimnoParam%sDPomS(iElem,:)))	                !fra��o de detritos do sedimento (-)
		
		Do pom=1,LimnoParam%npom
		    LimnoParam%afPomS(iElem,pom)=LimnoParam%sDPomS(iElem,pom)/(LimnoParam%aDTotS(iElem)*LimnoParam%afDOrgS(iElem))		                !fra��o de detritos na materia organica (-)
		EndDo
		!
		!Carbono
		LimnoParam%oH2CO3S(iElem)=LimnoParam%sH2CO3S(iElem)/LimnoParam%cDepthS/LimnoParam%bPorS
		LimnoParam%oHCO3S(iElem)=LimnoParam%sHCO3S(iElem)/LimnoParam%cDepthS/LimnoParam%bPorS
		LimnoParam%oCO3S(iElem)=LimnoParam%sCO3S(iElem)/LimnoParam%cDepthS/LimnoParam%bPorS
		LimnoParam%oDicS(iElem)=LimnoParam%sDicS(iElem)/LimnoParam%cDepthS/LimnoParam%bPorS
  !      
        !Nitrogenio
		LimnoParam%aNDissS(iElem)=LimnoParam%sNH4S(iElem)+LimnoParam%sNO3S(iElem)+SUM(LimnoParam%sNDomS(iElem,:))			                    !N dissolvido total no poro do sedimento (gN/m�)
		LimnoParam%aNkjAvailS(iElem)=SUM(LimnoParam%sNPomS(iElem,:))+SUM(LimnoParam%sNDomS(iElem,:))+LimnoParam%aNPhytS(iElem)+LimnoParam%sNH4S(iElem)	    !N kjeldahl no sedimento (excluindo humus) (gN/m�)
		LimnoParam%aNkjS(iElem)=LimnoParam%aNkjAvailS(iElem)          				                                !N kjeldahl no sedimento (gN/m�)
		LimnoParam%aNTotAvailS(iElem)=LimnoParam%aNkjAvailS(iElem)+LimnoParam%sNO3S(iElem)			                                !N total no sedimento (excluinda humus) (gN/m�)
		LimnoParam%aNTotS(iElem)=LimnoParam%aNkjS(iElem)+LimnoParam%sNO3S(iElem)        				                            !N total no sedimento (gN/m�)
		LimnoParam%afNInorS(iElem)=LimnoParam%aNDissS(iElem)/LimnoParam%aDTotS(iElem)				                                !fra��o de N inorganico no sedimento (gN/gD)
		LimnoParam%afNTotS(iElem)=LimnoParam%aNTotS(iElem)/LimnoParam%aDTotS(iElem)				                                !fra��o de N total no sedimento (gN/gD)
		LimnoParam%oNDissS(iElem)=LimnoParam%aNDissS(iElem)/LimnoParam%cDepthS/LimnoParam%bPorS				                                !conc. de N dissolvido na �gua intersticial (gN/m�)
		!
		!Fosforo			
		LimnoParam%aPInorgS(iElem)=LimnoParam%sPO4S(iElem)+LimnoParam%sPAIMS(iElem)	                                            !P inorganico no sedimento (gP/m�)
		LimnoParam%aPTotAvailS(iElem)=SUM(LimnoParam%sPPomS(iElem,:))+SUM(LimnoParam%sPDomS(iElem,:))+LimnoParam%aPInorgS(iElem)+LimnoParam%aPPhytS(iElem)	!P total no sedimento (excluindo humus, animais e vegeta��o) (gP/m�)
		LimnoParam%aPTotS(iElem)=LimnoParam%aPTotAvailS(iElem)            	                                        !P total no sedimento (excluindo animais e vegeta��o) (gP/m�)
		LimnoParam%afPInorgS(iElem)=LimnoParam%aPInorgS(iElem)/LimnoParam%aDTotS(iElem)	                                        !fra��o de fosforo inorganico no sedimento (gP/gD)
		LimnoParam%afPO4S(iElem)=LimnoParam%sPO4S(iElem)/LimnoParam%aPTotAvailS(iElem)		                                        !fra��o de f�sforo dissolvido no sedimento (-)
        
    EndDo !Loop Cell
    
    !Updating ratios in the sediment
    Call RatiosSediment(MeshParam,LimnoParam)
      

Return
End
