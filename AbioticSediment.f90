Subroutine AbioticSediment(HydroParam,MeshParam,LimnoParam,dt,dtday)
    
    Use Hydrodynamic
    Use MeshVars
    Use LimnologyVars
    Use uTempFunction

    Implicit None
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam
    type(LimnologyParam) :: LimnoParam
    Integer:: iElem,sed,dom,pom
    Real:: dt,dtday
    Real:: d1, d2
    Real:: NearZero = 1e-10
	Integer:: nphS,iter	
	Real:: aAlkS,aDicS
	Real:: aTkelvin,apH,npH,fw
	Real:: aPHS,phaux,Hion
	Real:: k1,k2,kwater
	Real:: incr,Fs
	Real:: a0,a1,a2
    
    !Abiotic processes in the sediment
    
    Do iElem = 1,MeshParam%nElem
        !-----------------------------------------------------------------------
        !  Temperature functions for sediment abiotic process
        !-----------------------------------------------------------------------
        !  Temperature_dependence_for_nitrification
           LimnoParam%uFunTmNitr=uFunTmAbio(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem),iElem),LimnoParam%cThetaNitr)
        !  temperature_function_of_diffusion
           LimnoParam%uFunTmDif= uFunTmAbio(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem),iElem),LimnoParam%cThetaDif)
        !-----------------------------------------------------------------------
        !  dissolved nutrients concentration in sediment(converting)
        !-----------------------------------------------------------------------
        !  conc._dissolved_N-NO3_in_interstitial_water
           LimnoParam%oNO3S(iElem) = LimnoParam%sNO3S(iElem) / LimnoParam%cDepthS / LimnoParam%bPorS
        !  conc._dissolved_N-NH4_in_interstitial_water
           LimnoParam%oNH4S(iElem) = LimnoParam%sNH4S(iElem) / LimnoParam%cDepthS / LimnoParam%bPorS
        !  conc._dissolved_P_in_interstitial_water
           LimnoParam%oPO4S(iElem) = LimnoParam%sPO4S(iElem) / LimnoParam%cDepthS / LimnoParam%bPorS
           
        
           
        !-----------------------------------------------------------------------
        !  Oxygen conditions in sediment
        !----------------------------------------------------------------------
	    If (LimnoParam%sO2W(HydroParam%ElSmallm(iElem),iElem)>0.1) Then
	        LimnoParam%akO2DifCor=LimnoParam%kO2Dif*LimnoParam%uFunTmDif*LimnoParam%cturbDifO2*LimnoParam%bPorCorS
	        !  sediment_oxygen_demand (only over labile fraction of DOM in the aerobic sediment)
	        d1 = LimnoParam%molO2molC*LimnoParam%cCPerDW*sum(LimnoParam%tDMinAerS(:,1)) 
	        d2 = LimnoParam%O2PerNH4*LimnoParam%molO2molN*LimnoParam%kNitrS*LimnoParam%uFunTmNitr*LimnoParam%sNH4S(iElem)
    !	    d3 = metanotrofia
	        LimnoParam%tSOD = (d1+d2)/LimnoParam%cDepthS 
            !  oxygen_penetration_depth
            LimnoParam%aDepthOxySed=min(LimnoParam%cDepthS,(2.0*LimnoParam%sO2W(HydroParam%ElSmallm(iElem),iElem)*LimnoParam%akO2DifCor/LimnoParam%tSOD)**0.5)
            !  fraction_aerobic_sediment
            LimnoParam%afOxySed=LimnoParam%aDepthOxySed/LimnoParam%cDepthS
            !  aerobic_mineralisation
            LimnoParam%tDMinOxyDetS=LimnoParam%afOxySed*sum(LimnoParam%tDMinAerS(:,1))
            !  sediment_oxygen_demand
            LimnoParam%tO2MinDetS=LimnoParam%molO2molC*LimnoParam%cCPerDW*LimnoParam%tDMinOxyDetS
            
        Else
            LimnoParam%aDepthOxySed=0.0
            LimnoParam%afOxySed=LimnoParam%aDepthOxySed/LimnoParam%cDepthS
        EndIf
        
        Do sed=1,LimnoParam%nsed
            
            !-----------------------------------------------------------------------
            !  Processes depending on Oxygen conditions in the sediment
            !----------------------------------------------------------------------
            If (LimnoParam%sO2W(HydroParam%ElSmallm(iElem),iElem)>0.1) Then !Aerobic and Anaerobic sediment
                If (sed==1) Then !Aerobic layer of the sediment
                    LimnoParam%sO2S(sed) = LimnoParam%sO2W(HydroParam%ElSmallm(iElem),iElem)
                    LimnoParam%aDepthS(sed) = LimnoParam%aDepthOxySed
                    LimnoParam%fDepthS(sed) = LimnoParam%aDepthOxySed/LimnoParam%cDepthS 
                Else  !Anaerobic layer of the sediment
                    LimnoParam%sO2S(sed) = 0.0
                    LimnoParam%aDepthS(sed) = LimnoParam%cDepthS-LimnoParam%aDepthOxySed
                    LimnoParam%fDepthS(sed) = 1.0-LimnoParam%aDepthOxySed/LimnoParam%cDepthS
                EndIf
            Else !Anaerobic sediment
                If (sed==1) Then !Aerobic layer of the sediment
                    LimnoParam%sO2S(sed) = 0.
                    LimnoParam%aDepthS(sed) = 0.
                    LimnoParam%fDepthS(sed) = 0. 
                Else  !Anaerobic layer of the sediment
                    LimnoParam%sO2S(sed) = 0.0
                    LimnoParam%aDepthS(sed) = LimnoParam%cDepthS
                    LimnoParam%fDepthS(sed) = 1.0
                EndIf
            EndIf
            
            !-----------------------------------------------------------------------
            !  Mineralization functions
            !-----------------------------------------------------------------------
            !  aerobic_mineralisation - DOM -> CO2,PO4,NH4
            Do dom=1,LimnoParam%ndom
                !  temp._function_of_aerobic_mineralisation
                LimnoParam%aTLimMinAerS  = LimnoParam%cThetaMinAerS ** (LimnoParam%sDTempW(HydroParam%ElSmallm(iElem),iElem)-20.)
                !  oxy_function_of_aerobic_mineralisation
                LimnoParam%aO2LimMinAerS = LimnoParam%sO2S(sed)/(LimnoParam%sO2S(sed)+LimnoParam%hO2Bac)

                LimnoParam%tDMinAerS(sed,dom) = LimnoParam%kMinAerDomS(dom) * LimnoParam%aO2LimMinAerS * LimnoParam%aTLimMinAerS * LimnoParam%sDDomS(iElem,dom) * LimnoParam%fDepthS(sed)	
                LimnoParam%tCMinAerS(sed,dom) = LimnoParam%rCDDomS(iElem,dom) * LimnoParam%tDMinAerS(sed,dom)
                LimnoParam%tNMinAerS(sed,dom) = LimnoParam%rNDDomS(iElem,dom) * LimnoParam%tDMinAerS(sed,dom)
                LimnoParam%tPMinAerS(sed,dom) = LimnoParam%rPDDomS(iElem,dom) * LimnoParam%tDMinAerS(sed,dom)        

     	        LimnoParam%tCMinAer2CO2S(sed,dom) = LimnoParam%tCMinAerS(sed,dom)
    	        LimnoParam%tCMinAer2CH4S(sed,dom) = 0.0
            EndDo    
            
            !  Anaerobic_mineralisation - DOM -> CO2,PO4,NH4
            Do dom=1,LimnoParam%ndom
                !  temp._function_of_anaerobic_mineralisation
    	        LimnoParam%aTLimMinAnaerS  = LimnoParam%cThetaMinAnaerS ** (LimnoParam%sDTempW(HydroParam%ElSmallm(iElem),iElem)-20.)
                !  oxy_function_of_anaerobic_mineralisation
    	        LimnoParam%aO2LimMinAnaerS = LimnoParam%kbacAn/(LimnoParam%sO2S(sed)+LimnoParam%kbacAn)
    	        
    	        LimnoParam%tDMinAnaerS(sed,dom) = LimnoParam%kMinAnaerDomS(dom) * LimnoParam%aO2LimMinAnaerS * LimnoParam%aTLimMinAnaerS * LimnoParam%sDDomS(iElem,dom) * LimnoParam%fDepthS(sed)	
    	        LimnoParam%tCMinAnaerS(sed,dom) = LimnoParam%rCDDomS(iElem,dom) * LimnoParam%tDMinAnaerS(sed,dom)
    	        LimnoParam%tNMinAnaerS(sed,dom) = LimnoParam%rNDDomS(iElem,dom) * LimnoParam%tDMinAnaerS(sed,dom)
    	        LimnoParam%tPMinAnaerS(sed,dom) = LimnoParam%rPDDomS(iElem,dom) * LimnoParam%tDMinAnaerS(sed,dom)
    	        
    	        !tCMinAnaer2CO2S(sed,dom) =     fDom2CO2(dom)  * tCMinAnaerS(sed,dom)
    	        !tCMinAnaer2CH4S(sed,dom) = (1.-fDom2CO2(dom)) * tCMinAnaerS(sed,dom)            
            EndDo
            
            !-----------------------------------------------------------------------
            !  denitrification flux
            !-----------------------------------------------------------------------
            !  mineralisation_flux_by_denitrification
            LimnoParam%tDDenitS=LimnoParam%oNO3S(iElem)*LimnoParam%oNO3S(iElem)/(LimnoParam%hNO3Denit*LimnoParam%hNO3Denit+LimnoParam%oNO3S(iElem)*LimnoParam%oNO3S(iElem))*(1.0-LimnoParam%fDepthS(sed))*(SUM(LimnoParam%tDMinAerS(sed,:))+SUM(LimnoParam%tDMinAnaerS(sed,:)))
            !  Denitrification_flux
            LimnoParam%tNDenitS(sed)=LimnoParam%NO3PerC*LimnoParam%molNmolC*LimnoParam%cCPerDW*LimnoParam%tDDenitS
            !-----------------------------------------------------------------------
            !  nitrification flux
            !-----------------------------------------------------------------------
            !  nitrification_flux
            LimnoParam%tNNitrS(sed)=LimnoParam%fDepthS(sed)*LimnoParam%kNitrS*LimnoParam%uFunTmNitr*LimnoParam%sNH4S(iElem)*LimnoParam%fDepthS(sed) 
            !  O2_flux_due_to_nitrification
            LimnoParam%tO2NitrS(sed)=LimnoParam%O2PerNH4*LimnoParam%molO2molN*LimnoParam%tNNitrS(sed)  
            
            !-----------------------------------------------------------------------
            !  absorbed P in sediment,oxygen dependent
            !-----------------------------------------------------------------------
            !  max._P_adsorption_per_g_inorg._matter_in_sediment
            LimnoParam%aPAdsMaxS = LimnoParam%cRelPAdsD+ LimnoParam%fDepthS(sed)*LimnoParam%cRelPAdsFe*LimnoParam%fFeDIM+LimnoParam%cRelPAdsAl*LimnoParam%fAlDIM
            !  P_adsorption_affinity,_corrected_for_redox_conditions
            LimnoParam%aKPAdsS=(1.0-LimnoParam%fRedMax*(1.0-LimnoParam%fDepthS(sed)))*LimnoParam%cKPAdsOx
            !  P_adsorption_isotherm_onto_inorg._matter_in_sediment
            LimnoParam%aPIsoAdsS=LimnoParam%aPAdsMaxS*LimnoParam%aKPAdsS*LimnoParam%oPO4S(iElem)/(1.0+LimnoParam%aKPAdsS*LimnoParam%oPO4S(iElem))
            !  equilibrium_amount
            LimnoParam%aPEqIMS = LimnoParam%aPIsoAdsS * LimnoParam%sDIMS(iElem)
            !  sorption
            LimnoParam%tPSorpIMS(sed) =LimnoParam%kPSorp*(LimnoParam%aPEqIMS-LimnoParam%sPAIMS(iElem))*LimnoParam%fDepthS(sed)
            
            !-----------------------------------------------------------------------
            !  Hydrolise in the sediment
            !-----------------------------------------------------------------------
            LimnoParam%aO2LimHidPomS = LimnoParam%sO2S(sed)/(LimnoParam%sO2S(sed)+LimnoParam%hO2Bac) + LimnoParam%fbacAn*(LimnoParam%kbacAn/(LimnoParam%sO2S(sed)+LimnoParam%kbacAn))
            
            Do pom=1,LimnoParam%npom
                If (LimnoParam%InclMatOrgSplit == 1) Then
		            LimnoParam%aTLimHidPomS(pom) = LimnoParam%cThetaHidS ** (LimnoParam%sDTempW(HydroParam%ElSmallm(iElem),iElem)-20.)
			        LimnoParam%kCorHidPomS   = LimnoParam%kHidPomS(pom) * LimnoParam%aO2LimHidPomS * LimnoParam%aTLimHidPomS(pom) 
			        LimnoParam%tDHidPomS(sed,pom) = LimnoParam%kCorHidPomS * LimnoParam%sDPomS(iElem,pom)
			        LimnoParam%tCHidPomS(sed,pom) = LimnoParam%rCDPomS(iElem,pom) * LimnoParam%tDHidPomS(sed,pom)
			        LimnoParam%tNHidPomS(sed,pom) = LimnoParam%rNDPomS(iElem,pom) * LimnoParam%tDHidPomS(sed,pom)
			        LimnoParam%tPHidPomS(sed,pom) = LimnoParam%rPDPomS(iElem,pom) * LimnoParam%tDHidPomS(sed,pom)    
                Else
                    LimnoParam%tDHidPomS(sed,pom) = 0.
                    LimnoParam%tCHidPomS(sed,pom) = 0.
                    LimnoParam%tNHidPomS(sed,pom) = 0.
                    LimnoParam%tPHidPomS(sed,pom) = 0.
                EndIf
            EndDo
            
            !-----------------------------------------------------------------------
            !  Methanotroph in the sediment
            !-----------------------------------------------------------------------
            !aTLimCH42CO2S  = cThetaCH4toCO2S ** (sDTempW(ElSmallm(iElem),iElem)-20.)
            !aO2LimCH42CO2S = sO2S(sed)/(sO2S(sed)+hO2CH42CO2)
            !aCH4LimCH42CO2S = sCH4S(I,J)/(sCH4S(I,J)+hCH4)
            !tCH42CO2S(sed) = kCH42CO2S(domg) * aO2LimCH42CO2S * aTLimMinDomS(dom) * aCH4LimCH42CO2S* sCH4S(I,J) * fDepthS(sed)
            
        EndDo

        !-----------------------------------------------------------------------
        !  diffusion process
        !-----------------------------------------------------------------------
        !  average_diffusion_distance
        LimnoParam%aDepthDif=LimnoParam%fDepthDifS*LimnoParam%cDepthS
        !  diffusion_flux_of_NH4_from_sediment_to_water
	    Do dom=1,LimnoParam%ndom
	        LimnoParam%tDDifDom(iElem,dom) = LimnoParam%kPDifDom(dom)*LimnoParam%uFunTmDif*LimnoParam%cTurbDifNut*LimnoParam%bPorCorS*(LimnoParam%oDDomS(iElem,dom)-LimnoParam%sDDomW(HydroParam%ElSmallm(iElem),iElem,dom))/LimnoParam%aDepthDif
	        LimnoParam%tCDifDom(iElem,dom) = LimnoParam%rCDDomS(iElem,dom) * LimnoParam%tDDifDom(iElem,dom)
	        LimnoParam%tNDifDom(iElem,dom) = LimnoParam%rNDDomS(iElem,dom) * LimnoParam%tDDifDom(iElem,dom)
	        LimnoParam%tPDifDom(iElem,dom) = LimnoParam%rPDDomS(iElem,dom) * LimnoParam%tDDifDom(iElem,dom)
	    EndDo
            
        !  diffusion_flux_of_NH4_from_sediment_to_water
        LimnoParam%tNDifNH4(iElem)=LimnoParam%kNDifNH4*LimnoParam%uFunTmDif*LimnoParam%cTurbDifNut*LimnoParam%bPorCorS*(LimnoParam%oNH4S(iElem)-LimnoParam%sNH4W(HydroParam%ElSmallm(iElem),iElem))/LimnoParam%aDepthDif
        !  diffusion_flux_of_NO3_from_sediment_to_water
        LimnoParam%tNDifNO3(iElem)=LimnoParam%kNDifNO3*LimnoParam%uFunTmDif*LimnoParam%cTurbDifNut*LimnoParam%bPorCorS*(LimnoParam%oNO3S(iElem)-LimnoParam%sNO3W(HydroParam%ElSmallm(iElem),iElem))/LimnoParam%aDepthDif
        !  diffusion_flux_of_dissolved_P_from_sediment_to_water
        LimnoParam%tPDifPO4(iElem)=LimnoParam%kPDifPO4*LimnoParam%uFunTmDif*LimnoParam%cTurbDifNut*LimnoParam%bPorCorS*(LimnoParam%oPO4S(iElem)-LimnoParam%sPO4W(HydroParam%ElSmallm(iElem),iElem))/LimnoParam%aDepthDif
        !  corrected_O2_diffusion_coefficient
        LimnoParam%akO2DifCor=LimnoParam%kO2Dif*LimnoParam%uFunTmDif*LimnoParam%cTurbDifO2*LimnoParam%bPorCorS
        !  O2_diffusion_(water_->_sediment)
        LimnoParam%tO2Dif(iElem)= LimnoParam%akO2DifCor*LimnoParam%sO2W(HydroParam%ElSmallm(iElem),iElem)/LimnoParam%aDepthDif
        !  CO2_diffusion_(water_->_sediment)
        LimnoParam%tCDifDic(iElem)=LimnoParam%kCDifDic*LimnoParam%uFunTmDif*LimnoParam%cTurbDifNut*LimnoParam%bPorCorS*(LimnoParam%oDicS(iElem)-LimnoParam%sDicW(HydroParam%ElSmallm(iElem),iElem))/LimnoParam%aDepthDif
	    !Fluxo de CH4
	    !tCDifDic=kPDifCH4*(cThetaDif**(TEMPI(I,J,NCAM(I,J))-20.))*cTurbDifNut*bPorCorS*(oCH4S(I,J)-sCH4W(I,J,NCAM(I,J)))/aDepthDif
           
        !!  chem._loss_of_dissolved_P_from_pore_water
        LimnoParam%tPChemPO4=max(0.0,LimnoParam%kPChemPO4*(LimnoParam%oPO4S(iElem)-LimnoParam%cPO4Max))
        
        !-----------------------------------------------------------------------
        !  Burial of sediment,contains erosion process
        !-----------------------------------------------------------------------
        
		LimnoParam%tDIMS =  LimnoParam%uDErosIMS(iElem)                                                                     & !Erosion
		       + LimnoParam%tDSetIM(HydroParam%ElSmallm(iElem),iElem)                                                       & !Sedimentation
		       - LimnoParam%tDResusIM(iElem)                                                                                !Resuspension
        
		Do pom=1,LimnoParam%npom
            
            If (LimnoParam%InclMatOrgSplit == 1) Then ! Case 1: Splited Organic Matter (POM and DOM)
                    
                If (LimnoParam%InclPomDomSplit == 1) Then ! Case 2: Splited Organic Matter (POM and DOM) and Splited fractions (Labile and Refratory)
                    If (pom==LimnoParam%LABIL) Then
		                LimnoParam%tDPomS(pom)=  LimnoParam%tDSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                              & ! Sedimentation POM
		                            - LimnoParam%tDResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                            - SUM(LimnoParam%tDHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                    + LimnoParam%fDLisDomPhyt*SUM(LimnoParam%tDLisPomPhytS(iElem,:))                                    & !Lise Fito
                                    + LimnoParam%fDMortDomMac*SUM(LimnoParam%tDMortMac(iElem,:))                                        & !Lise Macrophytes
                                    + LimnoParam%fDMessDomBen*SUM(LimnoParam%tDMessBentOM(iElem,:))                                     & !Messy Feeding Benthos
                                    + LimnoParam%fDEgesDomBen*SUM(LimnoParam%tDEgesBent(iElem,:))                                       & !Pellets Benthos
                                    + LimnoParam%fDMortDomBen*SUM(LimnoParam%tDMortBent(iElem,:))                                       !Mortality Benthos
                    ElseIf(pom==LimnoParam%REFRATARIA) Then 
		                LimnoParam%tDPomS(pom)=  LimnoParam%tDSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                              & ! Sedimentation POM
		                            - LimnoParam%tDResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                            - SUM(LimnoParam%tDHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                    + LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)                                        & ! Erosion
                                    + (1.-LimnoParam%fDLisDomPhyt)*SUM(LimnoParam%tDLisPomPhytS(iElem,:))                               & !Lise Fito
                                    + (1.-LimnoParam%fDMortDomMac)*SUM(LimnoParam%tDMortMac(iElem,:))                                   & !Lise Macrophytes
                                    + (1.-LimnoParam%fDMessDomBen)*SUM(LimnoParam%tDMessBentOM(iElem,:))                                & !Messy Feeding Benthos
                                    + (1.-LimnoParam%fDEgesDomBen)*SUM(LimnoParam%tDEgesBent(iElem,:))                                  & !Pellets Benthos
                                    + (1.-LimnoParam%fDMortDomBen)*SUM(LimnoParam%tDMortBent(iElem,:))                                  !Mortality Benthos
                    EndIf
                Else ! Case 3: Splited Organic Matter (POM and DOM) and non-splited fractions
		            LimnoParam%tDPomS(pom)=  LimnoParam%tDSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                              & ! Sedimentation POM
		                        - LimnoParam%tDResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                        - SUM(LimnoParam%tDHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                + LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)                                        & ! Erosion
                                + SUM(LimnoParam%tDLisPomPhytS(iElem,:))                                                            & !Lise Fito
                                + SUM(LimnoParam%tDMortMac(iElem,:))                                                                & !Lise Macrophytes
                                + SUM(LimnoParam%tDMessBentOM(iElem,:))                                                             & !Messy Feeding Benthos
                                + SUM(LimnoParam%tDEgesBent(iElem,:))                                                               & !Pellets Benthos
                                + SUM(LimnoParam%tDMortBent(iElem,:))                                                               !Mortality Benthos
                
                EndIf
            Else ! Case 4: Non-splited Organic Matter and non-splited fractions  (only one box of organic matter -  POM=Detritus - PCLake)
		        LimnoParam%tDPomS(pom)=  LimnoParam%tDSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                              & ! Sedimentation POM
		                    - LimnoParam%tDResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                    - SUM(LimnoParam%tDHidPomS(:,pom))                                                                  & ! Hydrolise POM
                            + LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)                                        & ! Erosion
                            + SUM(LimnoParam%tDMortBac(iElem,:))                                                                & !Lise Bacterioplankton
                            + SUM(LimnoParam%tDLisPomPhytS(iElem,:))                                                            & !Lise Fito
                            + SUM(LimnoParam%tDMortMac(iElem,:))                                                                & !Lise Macrophytes
                            + SUM(LimnoParam%tDMessBentOM(iElem,:))                                                             & !Messy Feeding Benthos
                            + SUM(LimnoParam%tDEgesBent(iElem,:))                                                               & !Pellets Benthos
                            + SUM(LimnoParam%tDMortBent(iElem,:))                                                               & !Mortality Benthos
                            - SUM(LimnoParam%tDMinAerS(:,pom))                                                                  & ! Mineralization
                            - SUM(LimnoParam%tDMinAnaerS(:,pom))                                                                & ! Mineralization
                            - LimnoParam%tDDifDom(iElem,pom)                                                                      ! Difusion
            EndIf
		EndDo
            
           
        LimnoParam%vDeltaS= (LimnoParam%tDIMS/LimnoParam%cRhoIM + SUM(LimnoParam%tDPomS(:))/LimnoParam%cRhoOM)/(1.0 - LimnoParam%bPorS)
            
        !  burial_flux_of_DW_in_inorganic_matter_in_lake
        if (LimnoParam%vDeltaS >= 0.0) then
            LimnoParam%tDBurIM = (SUM(LimnoParam%tDPomS(:)) +(LimnoParam%cRhoOM / LimnoParam%cRhoIM) * LimnoParam%tDIMS) / (SUM(LimnoParam%sDPomS(iElem,:)) /LimnoParam%sDIMS(iElem)  + LimnoParam%cRhoOM / LimnoParam%cRhoIM)
        else
            LimnoParam%tDBurIM = (SUM(LimnoParam%tDPomS(:)) +(LimnoParam%cRhoOM / LimnoParam%cRhoIM) * LimnoParam%tDIMS) / (LimnoParam%fDOrgSoil /(1.0 - LimnoParam%fDOrgSoil) + LimnoParam%cRhoOM / LimnoParam%cRhoIM) 
        endif
        
        !  burial_flux_in_lake
        LimnoParam%tDBurOM=SUM(LimnoParam%sDPomS(iElem,:))/LimnoParam%sDIMS(iElem)*LimnoParam%tDBurIM !OM = RPOM + LPOM
        LimnoParam%tDBurTot=LimnoParam%tDBurIM+LimnoParam%tDBurOM
        If (LimnoParam%vDeltaS >= 0.0) then
            LimnoParam%tDBurIM = (SUM(LimnoParam%tDPomS(:)) +(LimnoParam%cRhoOM / LimnoParam%cRhoIM) * LimnoParam%tDIMS) / (SUM(LimnoParam%sDPomS(iElem,:)) /LimnoParam%sDIMS(iElem)  + LimnoParam%cRhoOM / LimnoParam%cRhoIM)
            Do pom=1,LimnoParam%npom
			    If(pom==LimnoParam%LABIL)Then
			        LimnoParam%tDBurPom(pom)=LimnoParam%sDPomS(iElem,pom)/SUM(LimnoParam%sDPomS(iElem,:))*LimnoParam%tDBurOM
			        LimnoParam%tCBurPom(pom)=LimnoParam%rCDPomS(iElem,pom)*LimnoParam%tDBurPom(pom)
			        LimnoParam%tNBurPom(pom)=LimnoParam%rNDPomS(iElem,pom)*LimnoParam%tDBurPom(pom)
			        LimnoParam%tPBurPom(pom)=LimnoParam%rPDPomS(iElem,pom)*LimnoParam%tDBurPom(pom)
			    ElseIf(pom==LimnoParam%REFRATARIA)Then
			        LimnoParam%tDBurPom(pom)=LimnoParam%tDBurOM-LimnoParam%tDBurPom(pom-1)
			        LimnoParam%tCBurPom(pom)=LimnoParam%rCDPomS(iElem,pom)*LimnoParam%tDBurPom(pom)
			        LimnoParam%tNBurPom(pom)=LimnoParam%rNDPomS(iElem,pom)*LimnoParam%tDBurPom(pom)
			        LimnoParam%tPBurPom(pom)=LimnoParam%rPDPomS(iElem,pom)*LimnoParam%tDBurPom(pom)
			    EndIf
			EndDo
            Do dom=1,LimnoParam%ndom
			    LimnoParam%tDBurDom(dom)=LimnoParam%sDDomS(iElem,dom)*(LimnoParam%vDeltaS/LimnoParam%cDepthS)
			    LimnoParam%tCBurDom(dom)=LimnoParam%rCDPomS(iElem,dom)*LimnoParam%tDBurDom(dom)
			    LimnoParam%tNBurDom(dom)=LimnoParam%rNDPomS(iElem,dom)*LimnoParam%tDBurDom(dom)
			    LimnoParam%tPBurDom(dom)=LimnoParam%rPDPomS(iElem,dom)*LimnoParam%tDBurDom(dom)
			EndDo
			LimnoParam%tPBurPO4=LimnoParam%sPO4S(iElem)*(LimnoParam%vDeltaS/LimnoParam%cDepthS)
			LimnoParam%tPBurAIM=LimnoParam%sPAIMS(iElem)*(LimnoParam%vDeltaS/LimnoParam%cDepthS)
			LimnoParam%tNBurNH4=LimnoParam%sNH4S(iElem)*(LimnoParam%vDeltaS/LimnoParam%cDepthS)
			LimnoParam%tNBurNO3=LimnoParam%sNO3S(iElem)*(LimnoParam%vDeltaS/LimnoParam%cDepthS)
			LimnoParam%tCBurDic=LimnoParam%sDicS(iElem)*(LimnoParam%vDeltaS/LimnoParam%cDepthS)
        
        Else
			LimnoParam%tDBurIM=-(SUM(LimnoParam%tDPomS(:))+(LimnoParam%cRhoOM/LimnoParam%cRhoIM)*LimnoParam%tDIMS)/(LimnoParam%fDOrgSoil/(1.-LimnoParam%fDOrgSoil)+LimnoParam%cRhoOM/LimnoParam%cRhoIM)
			LimnoParam%tDBurOM=-(LimnoParam%fDOrgSoil/(1.-LimnoParam%fDOrgSoil))*LimnoParam%tDBurIM
			LimnoParam%tDBurTot=LimnoParam%tDBurIM+LimnoParam%tDBurOM
			Do pom=1,LimnoParam%npom
			    LimnoParam%tDBurPom(pom)=LimnoParam%tDBurOM
			    LimnoParam%tCBurPom(pom)=LimnoParam%rCDPomS(iElem,pom)*LimnoParam%tDBurPom(pom)
			    LimnoParam%tNBurPom(pom)=LimnoParam%rNDPomS(iElem,pom)*LimnoParam%tDBurPom(pom)
			    LimnoParam%tPBurPom(pom)=LimnoParam%rPDPomS(iElem,pom)*LimnoParam%tDBurPom(pom)
		    EndDo
			Do dom=1,LimnoParam%ndom
		        LimnoParam%tDBurDom(dom)=LimnoParam%cDomGround(dom)*(LimnoParam%vDeltaS/LimnoParam%cDepthS)
		        LimnoParam%tCBurDom(dom)=LimnoParam%rCDDomS(iElem,dom)*LimnoParam%tDBurDom(dom)
		        LimnoParam%tNBurDom(dom)=LimnoParam%rNDDomS(iElem,dom)*LimnoParam%tDBurDom(dom)
		        LimnoParam%tPBurDom(dom)=LimnoParam%rPDDomS(iElem,dom)*LimnoParam%tDBurDom(dom)
		    EndDo
			LimnoParam%tPBurPO4=LimnoParam%cPO4Ground*(LimnoParam%vDeltaS/LimnoParam%cDepthS)
			LimnoParam%tPBurAIM=0.0
			LimnoParam%tNBurNH4=LimnoParam%cNH4Ground*(LimnoParam%vDeltaS/LimnoParam%cDepthS)
			LimnoParam%tNBurNO3=LimnoParam%cNO3Ground*(LimnoParam%vDeltaS/LimnoParam%cDepthS)
			LimnoParam%tCBurDic=LimnoParam%cDicGround*(LimnoParam%vDeltaS/LimnoParam%cDepthS)			
        Endif
        
            
        !-------------------------------------------------------------------------------------!
        !                               Source/Sink Term for POM                              !
        !-------------------------------------------------------------------------------------!
        
	    Do pom=1,LimnoParam%npom
            
            If (LimnoParam%InclMatOrgSplit == 1) Then ! Case 1: Splited Organic Matter (POM and DOM)
                    
                If (LimnoParam%InclPomDomSplit == 1) Then ! Case 2: Splited Organic Matter (POM and DOM) and Splited fractions (Labile and Refratory)
                    If (pom==LimnoParam%LABIL) Then
		               LimnoParam%dVarEstS(iElem,1)=  LimnoParam%tDSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                         & ! Sedimentation POM
		                            - LimnoParam%tDResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                            - SUM(LimnoParam%tDHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                    + LimnoParam%fDLisDomPhyt*SUM(LimnoParam%tDLisPomPhytS(iElem,:))                                    & !Lise Fito
                                    + LimnoParam%fDMortDomMac*SUM(LimnoParam%tDMortMac(iElem,:))                                        & !Lise Macrophytes
                                    + LimnoParam%fDMessDomBen*SUM(LimnoParam%tDMessBentOM(iElem,:))                                     & !Messy Feeding Benthos
                                    + LimnoParam%fDEgesDomBen*SUM(LimnoParam%tDEgesBent(iElem,:))                                       & !Pellets Benthos
                                    + LimnoParam%fDMortDomBen*SUM(LimnoParam%tDMortBent(iElem,:))                                       & !Mortality Benthos
                                    - LimnoParam%tDBurPom(pom)                                                                          ! Burial
                    
		                LimnoParam%dVarEstS(iElem,2)=  LimnoParam%tCSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                        & ! Sedimentation POM
		                            - LimnoParam%tCResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                            - SUM(LimnoParam%tCHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                    + LimnoParam%fDLisDomPhyt*SUM(LimnoParam%tCLisPomPhytS(iElem,:))                                    & !Lise Fito
                                    + LimnoParam%fDMortDomMac*SUM(LimnoParam%tCMortMac(iElem,:))                                        & !Lise Macrophytes
                                    + LimnoParam%fDMessDomBen*SUM(LimnoParam%tCMessBentOM(iElem,:))                                     & !Messy Feeding Benthos
                                    + LimnoParam%fDEgesDomBen*SUM(LimnoParam%tCEgesBent(iElem,:))                                       & !Pellets Benthos
                                    + LimnoParam%fDMortDomBen*SUM(LimnoParam%tCMortBent(iElem,:))                                       & !Mortality Benthos
                                    - LimnoParam%tCBurPom(pom)                                                                          ! Burial

		                LimnoParam%dVarEstS(iElem,3)=  LimnoParam%tNSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                        & ! Sedimentation POM
		                            - LimnoParam%tNResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                            - SUM(LimnoParam%tNHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                    + LimnoParam%fDLisDomPhyt*SUM(LimnoParam%tNLisPomPhytS(iElem,:))                                    & !Lise Fito
                                    + LimnoParam%fDMortDomMac*SUM(LimnoParam%tNMortMac(iElem,:))                                        & !Lise Macrophytes
                                    + LimnoParam%fDMessDomBen*SUM(LimnoParam%tNMessBentOM(iElem,:))                                     & !Messy Feeding Benthos
                                    + LimnoParam%fDEgesDomBen*SUM(LimnoParam%tNEgesBent(iElem,:))                                       & !Pellets Benthos
                                    + LimnoParam%fDMortDomBen*SUM(LimnoParam%tNMortBent(iElem,:))                                       & !Mortality Benthos
                                    - LimnoParam%tNBurPom(pom)                                                                          ! Burial
                        
		                LimnoParam%dVarEstS(iElem,4)=  LimnoParam%tPSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                        & ! Sedimentation POM
		                            - LimnoParam%tPResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                            - SUM(LimnoParam%tPHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                    + LimnoParam%fDLisDomPhyt*SUM(LimnoParam%tPLisPomPhytS(iElem,:))                                    & !Lise Fito
                                    + LimnoParam%fDMortDomMac*SUM(LimnoParam%tPMortMac(iElem,:))                                        & !Lise Macrophytes
                                    + LimnoParam%fDMessDomBen*SUM(LimnoParam%tPMessBentOM(iElem,:))                                     & !Messy Feeding Benthos
                                    + LimnoParam%fDEgesDomBen*SUM(LimnoParam%tPEgesBent(iElem,:))                                       & !Pellets Benthos
                                    + LimnoParam%fDMortDomBen*SUM(LimnoParam%tPMortBent(iElem,:))                                       & !Mortality Benthos
                                    - LimnoParam%tPBurPom(pom)                                                                          ! Burial
                        
                    ElseIf(pom==LimnoParam%REFRATARIA) Then 
		                LimnoParam%dVarEstS(iElem,1)=  LimnoParam%tDSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                       & ! Sedimentation POM
		                            - LimnoParam%tDResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                            - SUM(LimnoParam%tDHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                    + LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)                                        & ! Erosion
                                    + (1.-LimnoParam%fDLisDomPhyt)*SUM(LimnoParam%tDLisPomPhytS(iElem,:))                               & !Lise Fito
                                    + (1.-LimnoParam%fDMortDomMac)*SUM(LimnoParam%tDMortMac(iElem,:))                                   & !Lise Macrophytes
                                    + (1.-LimnoParam%fDMessDomBen)*SUM(LimnoParam%tDMessBentOM(iElem,:))                                & !Messy Feeding Benthos
                                    + (1.-LimnoParam%fDEgesDomBen)*SUM(LimnoParam%tDEgesBent(iElem,:))                                  & !Pellets Benthos
                                    + (1.-LimnoParam%fDMortDomBen)*SUM(LimnoParam%tDMortBent(iElem,:))                                  & !Mortality Benthos
                                    - LimnoParam%tDBurPom(pom)                                                                          ! Burial
                        
		                LimnoParam%dVarEstS(iElem,2)=  LimnoParam%tCSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                       & ! Sedimentation POM
		                            - LimnoParam%tCResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                            - SUM(LimnoParam%tCHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                    + LimnoParam%cCDPomRef(pom)*LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)              & ! Erosion
                                    + (1.-LimnoParam%fDLisDomPhyt)*SUM(LimnoParam%tCLisPomPhytS(iElem,:))                               & !Lise Fito
                                    + (1.-LimnoParam%fDMortDomMac)*SUM(LimnoParam%tCMortMac(iElem,:))                                   & !Lise Macrophytes
                                    + (1.-LimnoParam%fDMessDomBen)*SUM(LimnoParam%tCMessBentOM(iElem,:))                                & !Messy Feeding Benthos
                                    + (1.-LimnoParam%fDEgesDomBen)*SUM(LimnoParam%tCEgesBent(iElem,:))                                  & !Pellets Benthos
                                    + (1.-LimnoParam%fDMortDomBen)*SUM(LimnoParam%tCMortBent(iElem,:))                                  & !Mortality Benthos
                                    - LimnoParam%tCBurPom(pom)                                                                          ! Burial
                        
		                LimnoParam%dVarEstS(iElem,3)=  LimnoParam%tNSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                       & ! Sedimentation POM
		                            - LimnoParam%tNResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                            - SUM(LimnoParam%tNHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                    + LimnoParam%cNDPomRef(pom)*LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)              & ! Erosion
                                    + (1.-LimnoParam%fDLisDomPhyt)*SUM(LimnoParam%tNLisPomPhytS(iElem,:))                               & !Lise Fito
                                    + (1.-LimnoParam%fDMortDomMac)*SUM(LimnoParam%tNMortMac(iElem,:))                                   & !Lise Macrophytes
                                    + (1.-LimnoParam%fDMessDomBen)*SUM(LimnoParam%tNMessBentOM(iElem,:))                                & !Messy Feeding Benthos
                                    + (1.-LimnoParam%fDEgesDomBen)*SUM(LimnoParam%tNEgesBent(iElem,:))                                  & !Pellets Benthos
                                    + (1.-LimnoParam%fDMortDomBen)*SUM(LimnoParam%tNMortBent(iElem,:))                                  & !Mortality Benthos
                                    - LimnoParam%tNBurPom(pom)                                                                          ! Burial
                        
		                LimnoParam%dVarEstS(iElem,4)=  LimnoParam%tPSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                       & ! Sedimentation POM
		                            - LimnoParam%tPResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                            - SUM(LimnoParam%tPHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                    + LimnoParam%cPDPomRef(pom)*LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)              & ! Erosion
                                    + (1.-LimnoParam%fDLisDomPhyt)*SUM(LimnoParam%tPLisPomPhytS(iElem,:))                               & !Lise Fito
                                    + (1.-LimnoParam%fDMortDomMac)*SUM(LimnoParam%tPMortMac(iElem,:))                                   & !Lise Macrophytes
                                    + (1.-LimnoParam%fDMessDomBen)*SUM(LimnoParam%tPMessBentOM(iElem,:))                                & !Messy Feeding Benthos
                                    + (1.-LimnoParam%fDEgesDomBen)*SUM(LimnoParam%tPEgesBent(iElem,:))                                  & !Pellets Benthos
                                    + (1.-LimnoParam%fDMortDomBen)*SUM(LimnoParam%tPMortBent(iElem,:))                                  & !Mortality Benthos
                                    - LimnoParam%tPBurPom(pom)                                                                          ! Burial
                        
                    EndIf
                Else ! Case 3: Splited Organic Matter (POM and DOM) and non-splited fractions
		            LimnoParam%dVarEstS(iElem,1)=  LimnoParam%tDSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                        & ! Sedimentation POM
		                        - LimnoParam%tDResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                        - SUM(LimnoParam%tDHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                + LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)                                        & ! Erosion
                                + SUM(LimnoParam%tDLisPomPhytS(iElem,:))                                                            & !Lise Fito
                                + SUM(LimnoParam%tDMortMac(iElem,:))                                                                & !Lise Macrophytes
                                + SUM(LimnoParam%tDMessBentOM(iElem,:))                                                             & !Messy Feeding Benthos
                                + SUM(LimnoParam%tDEgesBent(iElem,:))                                                               & !Pellets Benthos
                                + SUM(LimnoParam%tDMortBent(iElem,:))                                                               & !Mortality Benthos
                                - LimnoParam%tDBurPom(pom)                                                                          ! Burial
                    
		            LimnoParam%dVarEstS(iElem,2)=  LimnoParam%tCSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                        & ! Sedimentation POM
		                        - LimnoParam%tCResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                        - SUM(LimnoParam%tCHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                + LimnoParam%cCDPomRef(pom)*LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)              & ! Erosion
                                + SUM(LimnoParam%tCLisPomPhytS(iElem,:))                                                            & !Lise Fito
                                + SUM(LimnoParam%tCMortMac(iElem,:))                                                                & !Lise Macrophytes
                                + SUM(LimnoParam%tCMessBentOM(iElem,:))                                                             & !Messy Feeding Benthos
                                + SUM(LimnoParam%tCEgesBent(iElem,:))                                                               & !Pellets Benthos
                                + SUM(LimnoParam%tCMortBent(iElem,:))                                                               & !Mortality Benthos
                                - LimnoParam%tCBurPom(pom)                                                                          ! Burial
                    
		            LimnoParam%dVarEstS(iElem,3)=  LimnoParam%tNSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                        & ! Sedimentation POM
		                        - LimnoParam%tNResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                        - SUM(LimnoParam%tNHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                + LimnoParam%cNDPomRef(pom)*LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)              & ! Erosion
                                + SUM(LimnoParam%tNLisPomPhytS(iElem,:))                                                            & !Lise Fito
                                + SUM(LimnoParam%tNMortMac(iElem,:))                                                                & !Lise Macrophytes
                                + SUM(LimnoParam%tNMessBentOM(iElem,:))                                                             & !Messy Feeding Benthos
                                + SUM(LimnoParam%tNEgesBent(iElem,:))                                                               & !Pellets Benthos
                                + SUM(LimnoParam%tNMortBent(iElem,:))                                                               & !Mortality Benthos
                                - LimnoParam%tNBurPom(pom)                                                                          ! Burial
                    
		            LimnoParam%dVarEstS(iElem,4)=  LimnoParam%tPSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                        & ! Sedimentation POM
		                        - LimnoParam%tPResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                        - SUM(LimnoParam%tPHidPomS(:,pom))                                                                  & ! Hydrolise POM
                                + LimnoParam%cPDPomRef(pom)*LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)              & ! Erosion
                                + SUM(LimnoParam%tPLisPomPhytS(iElem,:))                                                            & !Lise Fito
                                + SUM(LimnoParam%tPMortMac(iElem,:))                                                                & !Lise Macrophytes
                                + SUM(LimnoParam%tPMessBentOM(iElem,:))                                                             & !Messy Feeding Benthos
                                + SUM(LimnoParam%tPEgesBent(iElem,:))                                                               & !Pellets Benthos
                                + SUM(LimnoParam%tPMortBent(iElem,:))                                                               & !Mortality Benthos
                                - LimnoParam%tPBurPom(pom)                                                                          ! Burial
                
                EndIf
            Else ! Case 4: Non-splited Organic Matter and non-splited fractions  (only one box of organic matter -  POM=Detritus - PCLake)
		        LimnoParam%dVarEstS(iElem,1)=  LimnoParam%tDSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                              & ! Sedimentation POM
		                    - LimnoParam%tDResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                    - SUM(LimnoParam%tDHidPomS(:,pom))                                                                  & ! Hydrolise POM
                            + LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)                                        & ! Erosion
                            + SUM(LimnoParam%tDMortBac(iElem,:))                                                                & !Lise Bacterioplankton
                            + SUM(LimnoParam%tDLisPomPhytS(iElem,:))                                                            & !Lise Fito
                            + SUM(LimnoParam%tDMortMac(iElem,:))                                                                & !Lise Macrophytes
                            + SUM(LimnoParam%tDMessBentOM(iElem,:))                                                             & !Messy Feeding Benthos
                            + SUM(LimnoParam%tDEgesBent(iElem,:))                                                               & !Pellets Benthos
                            + SUM(LimnoParam%tDMortBent(iElem,:))                                                               & !Mortality Benthos
                            - SUM(LimnoParam%tDMinAerS(:,pom))                                                                  & ! Mineralization
                            - SUM(LimnoParam%tDMinAnaerS(:,pom))                                                                & ! Mineralization
                            - LimnoParam%tDDifDom(iElem,pom)                                                                    & ! Difusion
                            - LimnoParam%tDBurPom(pom)                                                                          ! Burial
                
		        LimnoParam%dVarEstS(iElem,2)=  LimnoParam%tCSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                        & ! Sedimentation POM
		                    - LimnoParam%tCResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                    - SUM(LimnoParam%tCHidPomS(:,pom))                                                                  & ! Hydrolise POM
                            + LimnoParam%cCDPomRef(pom)*LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)              & ! Erosion
                            + SUM(LimnoParam%tCMortBac(iElem,:))                                                                & !Lise Bacterioplankton
                            + SUM(LimnoParam%tCLisPomPhytS(iElem,:))                                                            & !Lise Fito
                            + SUM(LimnoParam%tCMortMac(iElem,:))                                                                & !Lise Macrophytes
                            + SUM(LimnoParam%tCMessBentOM(iElem,:))                                                             & !Messy Feeding Benthos
                            + SUM(LimnoParam%tCEgesBent(iElem,:))                                                               & !Pellets Benthos
                            + SUM(LimnoParam%tCMortBent(iElem,:))                                                               & !Mortality Benthos
                            - SUM(LimnoParam%tCMinAerS(:,pom))                                                                  & ! Mineralization
                            - SUM(LimnoParam%tCMinAnaerS(:,pom))                                                                & ! Mineralization
                            - LimnoParam%tCDifDom(iElem,pom)                                                                    & ! Difusion
                            - LimnoParam%tCBurPom(pom)                                                                          ! Burial
                
		        LimnoParam%dVarEstS(iElem,3)=  LimnoParam%tNSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                        & ! Sedimentation POM
		                    - LimnoParam%tNResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                    - SUM(LimnoParam%tNHidPomS(:,pom))                                                                  & ! Hydrolise POM
                            + LimnoParam%cNDPomRef(pom)*LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)              & ! Erosion
                            + SUM(LimnoParam%tNMortBac(iElem,:))                                                                & !Lise Bacterioplankton
                            + SUM(LimnoParam%tNLisPomPhytS(iElem,:))                                                            & !Lise Fito
                            + SUM(LimnoParam%tNMortMac(iElem,:))                                                                & !Lise Macrophytes
                            + SUM(LimnoParam%tNMessBentOM(iElem,:))                                                             & !Messy Feeding Benthos
                            + SUM(LimnoParam%tNEgesBent(iElem,:))                                                               & !Pellets Benthos
                            + SUM(LimnoParam%tNMortBent(iElem,:))                                                               & !Mortality Benthos
                            - SUM(LimnoParam%tNMinAerS(:,pom))                                                                  & ! Mineralization
                            - SUM(LimnoParam%tNMinAnaerS(:,pom))                                                                & ! Mineralization
                            - LimnoParam%tNDifDom(iElem,pom)                                                                    & ! Difusion
                            - LimnoParam%tNBurPom(pom)                                                                          ! Burial
                
		        LimnoParam%dVarEstS(iElem,4)=  LimnoParam%tPSetPom(HydroParam%ElSmallm(iElem),iElem,pom)                        & ! Sedimentation POM
		                    - LimnoParam%tPResusPom(iElem,pom)                                                                  & ! Resuspension POM
		                    - SUM(LimnoParam%tPHidPomS(:,pom))                                                                  & ! Hydrolise POM
                            + LimnoParam%cPDPomRef(pom)*LimnoParam%tDErosPom(HydroParam%ElSmallm(iElem),iElem,pom)              & ! Erosion
                            + SUM(LimnoParam%tPMortBac(iElem,:))                                                                & !Lise Bacterioplankton
                            + SUM(LimnoParam%tPLisPomPhytS(iElem,:))                                                            & !Lise Fito
                            + SUM(LimnoParam%tPMortMac(iElem,:))                                                                & !Lise Macrophytes
                            + SUM(LimnoParam%tPMessBentOM(iElem,:))                                                             & !Messy Feeding Benthos
                            + SUM(LimnoParam%tPEgesBent(iElem,:))                                                               & !Pellets Benthos
                            + SUM(LimnoParam%tPMortBent(iElem,:))                                                               & !Mortality Benthos
                            - SUM(LimnoParam%tPMinAerS(:,pom))                                                                  & ! Mineralization
                            - SUM(LimnoParam%tPMinAnaerS(:,pom))                                                                & ! Mineralization
                            - LimnoParam%tPDifDom(iElem,pom)                                                                    & ! Difusion
                            - LimnoParam%tPBurPom(pom)                                                                          ! Burial
                
            EndIf
    	
	        LimnoParam%sDPomS(iElem,pom)=MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,1)+LimnoParam%sDPomSP(iElem,pom)))
	        LimnoParam%sCPomS(iElem,pom)=MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,2)+LimnoParam%sCPomSP(iElem,pom)))  
	        LimnoParam%sNPomS(iElem,pom)=MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,3)+LimnoParam%sNPomSP(iElem,pom)))  
	        LimnoParam%sPPomS(iElem,pom)=MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,4)+LimnoParam%sPPomSP(iElem,pom)))
	                 
        EndDo      
        
        
	    Do dom=1,LimnoParam%ndom
            
            If (LimnoParam%InclMatOrgSplit == 1) Then ! Case 1: Splited Organic Matter (POM and DOM)
                If (LimnoParam%InclPomDomSplit == 1) Then ! Case 2: Splited Organic Matter (POM and DOM) and Splited fractions (Labile and Refratory)
                    If (dom==LimnoParam%LABIL) Then
	                    LimnoParam%dVarEstS(iElem,1) = - LimnoParam%tDResusDom(iElem,dom)                                               & ! Ressusp POM
	                                - LimnoParam%tDDifDom(iElem,dom)                                                                    & ! Difusion
	                                + SUM(LimnoParam%tDHidPomS(:,dom))                                                                  & ! Hidrolise POM -> DOM
                                    - SUM(LimnoParam%tDMinAerS(:,dom))                                                                  & ! Mineralization
                                    - SUM(LimnoParam%tDMinAnaerS(:,dom))                                                                ! Mineralization
                        
	                    LimnoParam%dVarEstS(iElem,2) = - LimnoParam%tCResusDom(iElem,dom)                                               & ! Ressusp POM
	                                - LimnoParam%tCDifDom(iElem,dom)                                                                    & ! Difusion
	                                + SUM(LimnoParam%tCHidPomS(:,dom))                                                                  & ! Hidrolise POM -> DOM
                                    - SUM(LimnoParam%tCMinAerS(:,dom))                                                                  & ! Mineralization
                                    - SUM(LimnoParam%tCMinAnaerS(:,dom))                                                                ! Mineralization

	                    LimnoParam%dVarEstS(iElem,3) = - LimnoParam%tNResusDom(iElem,dom)                                               & ! Ressusp POM
	                                - LimnoParam%tNDifDom(iElem,dom)                                                                    & ! Difusion
	                                + SUM(LimnoParam%tNHidPomS(:,dom))                                                                  & ! Hidrolise POM -> DOM
                                    - SUM(LimnoParam%tNMinAerS(:,dom))                                                                  & ! Mineralization
                                    - SUM(LimnoParam%tNMinAnaerS(:,dom))                                                                ! Mineralization

	                    LimnoParam%dVarEstS(iElem,4) = - LimnoParam%tPResusDom(iElem,dom)                                               & ! Ressusp POM
	                                - LimnoParam%tPDifDom(iElem,dom)                                                                    & ! Difusion
	                                + SUM(LimnoParam%tPHidPomS(:,dom))                                                                  & ! Hidrolise POM -> DOM
                                    - SUM(LimnoParam%tPMinAerS(:,dom))                                                                  & ! Mineralization
                                    - SUM(LimnoParam%tPMinAnaerS(:,dom))                                                                ! Mineralization
                        
                    ElseIf(dom==LimnoParam%REFRATARIA) Then       
	                    LimnoParam%dVarEstS(iElem,1) = - LimnoParam%tDResusDom(iElem,dom)                                               & ! Ressusp POM
                                    + SUM(LimnoParam%tDMortBac(iElem,:))                                                                & !Lise Bacterioplankton
	                                + SUM(LimnoParam%tDHidPomS(:,dom))                                                                  & ! Hidrolise POM -> DOM
                                    - SUM(LimnoParam%tDMinAerS(:,dom))                                                                  & ! Mineralization
                                    - SUM(LimnoParam%tDMinAnaerS(:,dom))                                                                ! Mineralization
                        
	                    LimnoParam%dVarEstS(iElem,2) = - LimnoParam%tCResusDom(iElem,dom)                                               & ! Ressusp POM
                                    + SUM(LimnoParam%tCMortBac(iElem,:))                                                                & !Lise Bacterioplankton
	                                + SUM(LimnoParam%tCHidPomS(:,dom))                                                                  & ! Hidrolise POM -> DOM
                                    - SUM(LimnoParam%tCMinAerS(:,dom))                                                                  & ! Mineralization
                                    - SUM(LimnoParam%tCMinAnaerS(:,dom))                                                                ! Mineralization

	                    LimnoParam%dVarEstS(iElem,3) = - LimnoParam%tNResusDom(iElem,dom)                                               & ! Ressusp POM
                                    + SUM(LimnoParam%tNMortBac(iElem,:))                                                                & !Lise Bacterioplankton
	                                + SUM(LimnoParam%tNHidPomS(:,dom))                                                                  & ! Hidrolise POM -> DOM
                                    - SUM(LimnoParam%tNMinAerS(:,dom))                                                                  & ! Mineralization
                                    - SUM(LimnoParam%tNMinAnaerS(:,dom))                                                                ! Mineralization

	                    LimnoParam%dVarEstS(iElem,4) = - LimnoParam%tPResusDom(iElem,dom)                                               & ! Ressusp POM
                                    + SUM(LimnoParam%tPMortBac(iElem,:))                                                                & !Lise Bacterioplankton
	                                + SUM(LimnoParam%tPHidPomS(:,dom))                                                                  & ! Hidrolise POM -> DOM
                                    - SUM(LimnoParam%tPMinAerS(:,dom))                                                                  & ! Mineralization
                                    - SUM(LimnoParam%tPMinAnaerS(:,dom))                                                                ! Mineralization
                    EndIf
                Else ! Case 3: Splited Organic Matter (POM and DOM) and non-splited fractions    
	                LimnoParam%dVarEstS(iElem,1) = - LimnoParam%tDResusDom(iElem,dom)                                               & ! Ressusp POM
                                + SUM(LimnoParam%tDMortBac(iElem,:))                                                                & !Lise Bacterioplankton
	                            - LimnoParam%tDDifDom(iElem,dom)                                                                    & ! Difusion
	                            + SUM(LimnoParam%tDHidPomS(:,dom))                                                                  & ! Hidrolise POM -> DOM
                                - SUM(LimnoParam%tDMinAerS(:,dom))                                                                  & ! Mineralization
                                - SUM(LimnoParam%tDMinAnaerS(:,dom))                                                                ! Mineralization
 
	                LimnoParam%dVarEstS(iElem,2) = - LimnoParam%tCResusDom(iElem,dom)                                               & ! Ressusp POM
                                + SUM(LimnoParam%tCMortBac(iElem,:))                                                                & !Lise Bacterioplankton
	                            - LimnoParam%tCDifDom(iElem,dom)                                                                    & ! Difusion
	                            + SUM(LimnoParam%tCHidPomS(:,dom))                                                                  & ! Hidrolise POM -> DOM
                                - SUM(LimnoParam%tCMinAerS(:,dom))                                                                  & ! Mineralization
                                - SUM(LimnoParam%tCMinAnaerS(:,dom))                                                                ! Mineralization

	                LimnoParam%dVarEstS(iElem,3) = - LimnoParam%tNResusDom(iElem,dom)                                               & ! Ressusp POM
                                + SUM(LimnoParam%tNMortBac(iElem,:))                                                                & !Lise Bacterioplankton
	                            - LimnoParam%tNDifDom(iElem,dom)                                                                    & ! Difusion
	                            + SUM(LimnoParam%tNHidPomS(:,dom))                                                                  & ! Hidrolise POM -> DOM
                                - SUM(LimnoParam%tNMinAerS(:,dom))                                                                  & ! Mineralization
                                - SUM(LimnoParam%tNMinAnaerS(:,dom))                                                                ! Mineralization

	                LimnoParam%dVarEstS(iElem,4) = - LimnoParam%tPResusDom(iElem,dom)                                               & ! Ressusp POM
                                + SUM(LimnoParam%tPMortBac(iElem,:))                                                                & !Lise Bacterioplankton
	                            - LimnoParam%tPDifDom(iElem,dom)                                                                    & ! Difusion
	                            + SUM(LimnoParam%tPHidPomS(:,dom))                                                                  & ! Hidrolise POM -> DOM
                                - SUM(LimnoParam%tPMinAerS(:,dom))                                                                  & ! Mineralization
                                - SUM(LimnoParam%tPMinAnaerS(:,dom))                                                                ! Mineralization
                    
                EndIf
            Else ! Case 4: Non-splited Organic Matter and non-splited fractions  (only one box of organic matter -  POM=Detritus - PCLake)

            EndIf            
            
            If (LimnoParam%InclMatOrgSplit == 1) Then                                              
	            LimnoParam%sDDomS(iElem,dom)=MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,1)+LimnoParam%sDDomSP(iElem,dom)))
	            LimnoParam%sCDomS(iElem,dom)=MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,2)+LimnoParam%sCDomSP(iElem,dom)))  
	            LimnoParam%sNDomS(iElem,dom)=MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,3)+LimnoParam%sNDomSP(iElem,dom)))  
	            LimnoParam%sPDomS(iElem,dom)=MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,4)+LimnoParam%sPDomSP(iElem,dom)))
            Else
	            LimnoParam%sDDomS(iElem,dom)= LimnoParam%sDPomS(iElem,dom)
	            LimnoParam%sCDomS(iElem,dom)= LimnoParam%sCPomS(iElem,dom)
	            LimnoParam%sNDomS(iElem,dom)= LimnoParam%sNPomS(iElem,dom)
	            LimnoParam%sPDomS(iElem,dom)= LimnoParam%sPPomS(iElem,dom)
            EndIf
            
             
        EndDo
        
	    !PO4Sed
	    LimnoParam%dVarEstS(iElem,1) =  - LimnoParam%tPResusPO4(iElem)                  & ! Ressuspension PO4S -> PO4W
	                        - LimnoParam%tPDifPO4(iElem)                                & ! Difusion PO4S <-> PO4W
                            - LimnoParam%tPChemPO4                                      & ! Imobilization
	                        + SUM(LimnoParam%tPMinAerS(:,:))                            & ! Mineralization DOM -> PO4S
	                        + SUM(LimnoParam%tPMinAnaerS(:,:))                          & ! Mineralization DOM -> PO4S
                            - SUM(LimnoParam%tPSorpIMS(:))                              & ! Soption PO4S -> IMSed
	                        + SUM(LimnoParam%tPExdMacS(iElem,:))                        & ! Respiration MacSed -> PO4S
	                        - SUM(LimnoParam%tPUptMacS(iElem,:))                        & ! Upatake MacroSed <- PO4S
	                        - LimnoParam%tPBurPO4                                       ! Buried PO4top->PO4bottom
	       
	    LimnoParam%sPO4S(iElem)=Max(NearZero,(dtday*LimnoParam%dVarEstS(iElem,1)+LimnoParam%sPO4SP(iElem)))               
        
        !------------------------------------------------------------------!
	    !PAIMSed
	    LimnoParam%dVarEstS(iElem,1) = LimnoParam%tPSetAIM(HydroParam%ElSmallm(iElem),iElem)        & ! Sedimentation
	                        - LimnoParam%tPResusAIM(iElem)                                          & ! Ressuspension
	                        + SUM(LimnoParam%tPSorpIMS(:))                                          & ! Soption
                            - LimnoParam%tPBurAIM                                                   ! Buried   
	
	    LimnoParam%sPAIMS(iElem)=Max(NearZero,(dtday*LimnoParam%dVarEstS(iElem,1)+LimnoParam%sPAIMSP(iElem))) 
        
        !------------------------------------------------------------------!
	    !NH4Sed		
	    LimnoParam%dVarEstS(iElem,1) =  - LimnoParam%tNResusNH4(iElem)                  & ! Ressuspension NH4S -> NH4W
	                        - LimnoParam%tNDifNH4(iElem)                                & ! Difusion NH4S <-> NH4W
	                        - SUM(LimnoParam%tNNitrS(:))                                & ! Nitrification NH4 -> NO3
                            + SUM(LimnoParam%tNMinAerS(:,:))                            & ! Mineralization DOM -> NH4S
	                        + SUM(LimnoParam%tNMinAnaerS(:,:))                          & ! Mineralization DOM -> NH4S
	                        + SUM(LimnoParam%tNExdMacS(iElem,:))                        & ! Respiration MACSed -> NH4
	                        - SUM(LimnoParam%tNUptNH4MacS(iElem,:))                     & ! Uptake MacroSed NH4S
	                        - LimnoParam%tNBurNH4                                       ! Buried
	
	    LimnoParam%sNH4S(iElem)=MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,1)+LimnoParam%sNH4SP(iElem)))

        !------------------------------------------------------------------!
	    !NO3
	    LimnoParam%dVarEstS(iElem,1) =  - LimnoParam%tNResusNO3(iElem)                  & ! Ressuspension NH4S -> NH4W
	                        - LimnoParam%tNDifNO3(iElem)                                & ! Difusion NH4S <-> NH4W
                            + SUM(LimnoParam%tNNitrS(:))                                & ! Nitrification NH4 -> NO3
	                        - SUM(LimnoParam%tNDenitS(:))                               & ! Denitrification NO3 -> N2
	                        - SUM(LimnoParam%tNUptNO3MacS(iElem,:))                     & ! Uptake MacroSed NH4S
	                        - LimnoParam%tNBurNO3                                       ! Buried

	    LimnoParam%sNO3S(iElem)=MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,1)+LimnoParam%sNO3SP(iElem))) 
        
        !------------------------------------------------------------------!
	    !IMSed
		LimnoParam%dVarEstS(iElem,1) =  LimnoParam%uDErosIMS(iElem)                     & !Erosion
		                    + LimnoParam%tDSetIM(HydroParam%ElSmallm(iElem),iElem)      & !Sedimentation
		                    - LimnoParam%tDResusIM(iElem)                               & !Resuspension
	                        - LimnoParam%tDBurIM 	                                    ! Buried
	
	    LimnoParam%sDIMS(iElem)=MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,1)+LimnoParam%sDIMSP(iElem))) 
        
        !------------------------------------------------------------------!
        !DICSed - mg/L
	    LimnoParam%dVarEstS(iElem,1) = - LimnoParam%tCDifDic(iElem)                     & ! Difusion DicS -> DicW
                            !+ SUM(tCMinDom2CO2W(:,:))     & ! Mineralizacao DOM -> DicS
                            + SUM(LimnoParam%tCExcrPhytS(iElem,:))                      & ! Respiration Phyto -> DicS
                            + SUM(LimnoParam%tCRespBent(iElem,:))                       & ! Respiration Zoops -> DicS
                            + SUM(LimnoParam%tCExdMacS(iElem,:))                        & ! Respiration Mac -> DicS
                            - SUM(LimnoParam%tCUptMacS(iElem,:))                         ! Uptake MacroSed DicS
	
	    LimnoParam%sDicS(iElem)=MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,1)+LimnoParam%sDicSP(iElem))) 
        
        !------------------------------------------------------------------!
        ! pH and carbonate species
        aAlkS = LimnoParam%sAlkS(iElem) / 5.0E+04  !mg/L CaCO3 -> eq/L 
        ! mgCaCO3-- => eq 
        ! 1gCaCO3--  = 2/100 eq
        ! 1mgCaCO3-- = 2/(1E5) eq 
        ! 100 = massa molar CaCO3
        !   2 = carga
        aDicS = LimnoParam%sDicS(iElem)/12000.  !mol/L   
        ! molar(mol/L)
        ! 1molC = 12g
        ! 1g = 1mol/12g/1000
        aPHS = LimnoParam%spHS(iElem)

        k1=10**-6.3
        k2=10**-10.3
        kwater=10**-14.

        phaux = -aPHS-2.1
        If (aPH <= 0.0) phaux = -14.0
        incr = 10.0
        Do nph=1,3
        Fw    = 1.0
        incr = incr/10.0
        iter = 0
            Do While (Fw > 0.0 .and. iter < 12)
                phaux  = phaux+INCR
                Hion   = 10.0**phaux
                a0 = Hion**2 / (Hion*Hion + k1*Hion + k1*k2)
                a1 = k1*Hion / (Hion*Hion + k1*Hion + k1*k2)
                a2 = k1*k2   / (Hion*Hion + k1*Hion + k1*k2)
                Fw = (a1 + 2*a2) * aDicS + kwater/Hion - Hion - aAlkS
                iter   = iter+1
            EndDo
            phaux = phaux-INCR
        EndDo

        Hion =  10.0 ** phaux
        LimnoParam%spHS(iElem) = - phaux

        a0 = Hion*Hion/(Hion*Hion + k1*Hion + k1*k2)
        a1 = k1*Hion / (Hion*Hion + k1*Hion + k1*k2)
        a2 = k1*k2   / (Hion*Hion + k1*Hion + k1*k2)

        LimnoParam%sH2CO3S(iElem)  = LimnoParam%sDicS(iElem)  *a0 
        LimnoParam%sHCO3S(iElem) =   LimnoParam%sDicS(iElem) * a1 
        LimnoParam%sCO3S(iElem)  =   LimnoParam%sDicS(iElem) * a2
        
        !-------------------------------------------------------------------------!
        !  Alkalinity correction  		                                          !
        !-------------------------------------------------------------------------!
        
        LimnoParam%dVarEstS(iElem,1) =  -SUM(LimnoParam%tNNitrS(:))*(50000.*0.01)*2./14.007             & !Nitrification   - Diminui:  NH4+ + 2O2 -> NO3- + H2O + 2H+
                                        -SUM(LimnoParam%tNDenitS(:))*(50000.*0.01)*1./14.007            & !Denitrificacao - Aumenta:  5CH2O + 4NO3- + 4H+ -> 5CO2 + 2N2 + 7H2O
                                        +SUM(LimnoParam%tNUptNH4MacS(iElem,:))*(50000.*0.001)/14.007    & !Fotossintese - Aumenta:  106CO2 + 16NO3- + HPO4-2 + 122H2O + 18H+ -> (C106.H263.O110.N16.P1) + 138O2
                                        +SUM(LimnoParam%tNUptNO3MacS(iElem,:))*(50000.*0.001)/14.007    & !Fotossintese - Aumenta:  106CO2 + 16NO3- + HPO4-2 + 122H2O + 18H+ -> (C106.H263.O110.N16.P1) + 138O2
                                        +SUM(LimnoParam%tPUptMacS(iElem,:))*(50000.*0.001)/31.974
            
        LimnoParam%sAlkS(iElem) = MAX(NearZero,(LimnoParam%sAlkS(iElem) + DtDay*LimnoParam%dVarEstS(iElem,1)*aAlkS))
        
        
        
    EndDo !Loop Cell
      

Return
End
