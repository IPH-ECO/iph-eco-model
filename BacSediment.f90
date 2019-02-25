Subroutine BacSediment(HydroParam,MeshParam,LimnoParam,dt,dtday)
    
    Use Hydrodynamic
    Use MeshVars
    Use LimnologyVars
    Use uTempFunction

    Implicit None
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam
    type(LimnologyParam) :: LimnoParam
    Integer:: iElem,dom,bb
    Real:: dt,dtday,aDOMW
    Real:: NearZero = 1e-10
    
    
    Do bb = 1, LimnoParam%nbac
        Do iElem = 1,MeshParam%nElem
            !-------------------------------------------------------------------------------------!
            !                           Bacterioplankton rocesses in Sediment                     !
            !-------------------------------------------------------------------------------------!
            !---------------------<BAC>-------------------------!
            !            Processes								!
            !			(1) DOC => BAC							!
            !			(2) Respiration							!
            !			(3) Mortality			                !
            !			(4) Uptake/Release Nutrients            !
            !			(5) Sedimentation                       !
            !			(6) Grazing by Zoops                    !
            !---------------------------------------------------!
            !-------------------------------------------------------------------------!
			!-------------------------------------------------------------------------!
			! (1) DOM Assimilation 				                                      !

			! 1.1 Limiting Factors
		    !Temp
		    LimnoParam%aTLimMinBac = LimnoParam%cThetaBac(bb) ** (LimnoParam%sDTempW(HydroParam%ElSmallm(iElem),iElem)-20.)
			
			!pH
			If (LimnoParam%spHS(iElem).GT.LimnoParam%pHmax) Then
				LimnoParam%apHLimMinBac = exp(LimnoParam%pHmax - LimnoParam%spHS(iElem))
			ElseIf (LimnoParam%spHS(iElem).LT.LimnoParam%pHmin) Then
				LimnoParam%apHLimMinBac = exp(LimnoParam%spHS(iElem) - LimnoParam%pHmin)
			Else
				LimnoParam%apHLimMinBac = 1.0
			EndIf
			
			!O2
			If (LimnoParam%abac(bb)==LimnoParam%AEROBICA) Then
			    LimnoParam%aO2LimMinBac = LimnoParam%sO2W(HydroParam%ElSmallm(iElem),iElem)/(LimnoParam%sO2W(HydroParam%ElSmallm(iElem),iElem)+LimnoParam%hO2Bac+NearZero)
			ElseIf (LimnoParam%abac(bb)==LimnoParam%ANAEROBICA) Then
			    LimnoParam%aO2LimMinBac = LimnoParam%kbacAn/(LimnoParam%sO2W(HydroParam%ElSmallm(iElem),iElem)+LimnoParam%kbacAn+NearZero)
			Else
			    LimnoParam%aO2LimMinBac = (1.-LimnoParam%fbacAn)*LimnoParam%sO2W(HydroParam%ElSmallm(iElem),iElem)/(LimnoParam%sO2W(HydroParam%ElSmallm(iElem),iElem)+LimnoParam%hO2Bac+NearZero) + LimnoParam%fbacAn*(LimnoParam%kbacAn/(LimnoParam%sO2W(HydroParam%ElSmallm(iElem),iElem)+LimnoParam%kbacAn+NearZero))
			EndIf
			
			!Internal nutrients
			LimnoParam%aNutLimMinBac = Min(LimnoParam%rNDBacS(iElem,bb)/LimnoParam%cNDBacRef(bb) , LimnoParam%rPDBacS(iElem,bb)/LimnoParam%cPDBacRef(bb))
			
			!Substrate availability : Saturation - Pref * FOOD -> oDFood"Predador""Presa"(pred,food)	
    	    Do dom = 1, LimnoParam%ndom
    		    LimnoParam%oDFoodBacDom(dom) = LimnoParam%cPrefBacDom(bb,dom) * LimnoParam%sDDomS(iElem,dom) 
    	    EndDo
    	    LimnoParam%oDFoodBac = sum(LimnoParam%oDFoodBacDom(:))
			LimnoParam%aDFoodSatBac = LimnoParam%oDFoodBac / (LimnoParam%oDFoodBac + LimnoParam%hAssBac(bb)+NearZero)
                
			! 1.2 Uptake
			LimnoParam%tDAssBac(iElem,bb) = LimnoParam%MuBac(bb) * LimnoParam%apHLimMinBac * LimnoParam%aO2LimMinBac * LimnoParam%aTLimMinBac * LimnoParam%aNutLimMinBac * LimnoParam%aDFoodSatBac * LimnoParam%sDBacS(iElem,bb) 

            ! Grazing Function x Assimilation = G(D) x Fgrz
            Do dom = 1, LimnoParam%ndom
		        LimnoParam%fDAssBacDom(dom) = LimnoParam%oDFoodBacDom(dom) / LimnoParam%oDFoodBac
                LimnoParam%tDAssBacDom(iElem,bb,dom) = LimnoParam%fDAssBacDom(dom) * LimnoParam%tDAssBac(iElem,bb)
		        LimnoParam%tCAssBacDom(iElem,bb,dom) = LimnoParam%tDAssBacDom(iElem,bb,dom) * LimnoParam%rCDDomS(iElem,dom)
		        LimnoParam%tNAssBacDom(iElem,bb,dom) = LimnoParam%tDAssBacDom(iElem,bb,dom) * LimnoParam%rNDDomS(iElem,dom)
	            LimnoParam%tPAssBacDom(iElem,bb,dom) = LimnoParam%tDAssBacDom(iElem,bb,dom) * LimnoParam%rPDDomS(iElem,dom)
            EndDo
            ! Total assimilated nutrients of preys
            LimnoParam%tCAssBac(iElem,bb) =  SUM(LimnoParam%tCAssBacDom(iElem,bb,:)) 
            LimnoParam%tNAssBac(iElem,bb) =  SUM(LimnoParam%tNAssBacDom(iElem,bb,:))
            LimnoParam%tPAssBac(iElem,bb) =  SUM(LimnoParam%tPAssBacDom(iElem,bb,:))

                
    		!-------------------------------------------------------------------------!
			! (2) Respiration BAC => CO2 e CH4										  !
			! 2.1 basal respiration - Respb
            LimnoParam%kCorResbBac = LimnoParam%kResbBac(bb) * LimnoParam%aTLimMinBac
            LimnoParam%tDRespbBac = LimnoParam%kCorResbBac * LimnoParam%sDBacS(iElem,bb)         
			! 2.2 Parcela de perda por respiracao devido atividade metabolica (ERSEM) - Respa
            !wDRespaBac = fBacLoss(bb) * wDAssBac(iLayer,iElem,bb) 		
			! 2.3 Total Respiration -> basal + atividade metabolica(mg/m3/s)
            LimnoParam%tDRespBac(iElem,bb) = LimnoParam%tDRespbBac !+ wDRespaBac
                
            !  C_Respiration
            LimnoParam%tCRespBac(iElem,bb) = LimnoParam%rCDBacS(iElem,bb) / LimnoParam%cCDBacRef(bb) * LimnoParam%kResbBac(bb) * LimnoParam%aTLimMinBac*LimnoParam%sCBacS(iElem,bb)
   
            LimnoParam%tCRespBac2CO2(iElem,bb) =     LimnoParam%fBac2CO2(bb)  * LimnoParam%tCRespBac(iElem,bb) 
			LimnoParam%tCRespBac2CH4(iElem,bb) = (1.-LimnoParam%fBac2CO2(bb)) * LimnoParam%tCRespBac(iElem,bb) 
                
            !  N_excretion
            LimnoParam%tNExcrBac(iElem,bb) = LimnoParam%rNDBacS(iElem,bb) / LimnoParam%cNDBacRef(bb) * LimnoParam%kResbBac(bb) * LimnoParam%aTLimMinBac*LimnoParam%sNBacS(iElem,bb)
            !  P_excretion
            LimnoParam%tPExcrBac(iElem,bb) = LimnoParam%rPDBacS(iElem,bb) / LimnoParam%cPDBacRef(bb) * LimnoParam%kResbBac(bb) * LimnoParam%aTLimMinBac*LimnoParam%sPBacS(iElem,bb)
                
                
			!-------------------------------------------------------------------------!
			! (3) Mortality      												  !
            LimnoParam%kCorMortBac = LimnoParam%aTLimMinBac * LimnoParam%kMortBac(bb)
            LimnoParam%tDMortBac(iElem,bb) = LimnoParam%kCorMortBac * LimnoParam%sDBacS(iElem,bb)
            LimnoParam%tCMortBac(iElem,bb) = LimnoParam%rCDBacS(iElem,bb) * LimnoParam%tDMortBac(iElem,bb)
            LimnoParam%tNMortBac(iElem,bb) = LimnoParam%rNDBacS(iElem,bb) * LimnoParam%tDMortBac(iElem,bb)
            LimnoParam%tPMortBac(iElem,bb) = LimnoParam%rPDBacS(iElem,bb) * LimnoParam%tDMortBac(iElem,bb)
   
			!-------------------------------------------------------------------------!
			! (4) Uptake/Release PO4,NH4				                              !
			! 4.1 PO4 Uptake/Release (CAEDYM)
			If (LimnoParam%tDAssBac(iElem,bb)*LimnoParam%cPDBacRef(bb).GT.LimnoParam%tPAssBac(iElem,bb)) Then
			    LimnoParam%tPUptBac(iElem,bb) = (LimnoParam%tDAssBac(iElem,bb)*LimnoParam%cPDBacRef(bb) - LimnoParam%tPAssBac(iElem,bb)) - LimnoParam%tDRespBac(iElem,bb) * LimnoParam%rPDBacS(iElem,bb)
			Else
			    LimnoParam%tPUptBac(iElem,bb) = - LimnoParam%tDRespBac(iElem,bb) * LimnoParam%rPDBacS(iElem,bb)
			EndIf
   
			! 4.2 NH4 Uptake/Release (CAEDYM)		    
		    If (LimnoParam%tDAssBac(iElem,bb)*LimnoParam%cNDBacRef(bb).GT.LimnoParam%tNAssBac(iElem,bb)) Then
			    LimnoParam%tNUptBac(iElem,bb) = (LimnoParam%tDAssBac(iElem,bb)*LimnoParam%cNDBacRef(bb) - LimnoParam%tNAssBac(iElem,bb)) - LimnoParam%tDRespBac(iElem,bb) * LimnoParam%rNDBacS(iElem,bb)
			Else
			    LimnoParam%tNUptBac(iElem,bb) = - LimnoParam%tDRespBac(iElem,bb) * LimnoParam%rNDBacS(iElem,bb)
            EndIf
                
            
            !-------------------------------------------------------------------------------------!
            !                               Source/Sink Term for Phytoplankton in the sediment                       !
            !-------------------------------------------------------------------------------------!
            
			LimnoParam%dVarEstS(iElem,1) =  LimnoParam%tDAssBac(iElem,bb)                               & !Assimilation
                                            - LimnoParam%tDRespBac(iElem,bb)                            & !Respiration
                                            - LimnoParam%tDMortBac(iElem,bb)                            & !Mortality
                                            - LimnoParam%tDSetBac(HydroParam%ElSmallm(iElem),iElem,bb)  & !Sedimentation
                                            -SUM(LimnoParam%tDGrzBentBac(iElem,:,bb))                   !Grazing by Benthos

			LimnoParam%dVarEstS(iElem,2) =  LimnoParam%tCAssBac(iElem,bb)                               & !Assimilation
                                            - LimnoParam%tCRespBac(iElem,bb)                            & !Respiration
                                            - LimnoParam%tCMortBac(iElem,bb)                            & !Mortality
                                            - LimnoParam%tCSetBac(HydroParam%ElSmallm(iElem),iElem,bb)  & !Sedimentation
                                            -SUM(LimnoParam%tCGrzBentBac(iElem,:,bb))                   !Grazing by Benthos

			LimnoParam%dVarEstS(iElem,3) =  LimnoParam%tNAssBac(iElem,bb)                               & !Assimilation
                                            - LimnoParam%tNExcrBac(iElem,bb)                            & !Excretion
                                            - LimnoParam%tNMortBac(iElem,bb)                            & !Mortality
                                            - LimnoParam%tNSetBac(HydroParam%ElSmallm(iElem),iElem,bb)  & !Sedimentation
                                            -SUM(LimnoParam%tNGrzBentBac(iElem,:,bb))                   !Grazing by Benthos
  
			LimnoParam%dVarEstS(iElem,4) =  LimnoParam%tPAssBac(iElem,bb)                               & !Assimilation
                                            - LimnoParam%tPExcrBac(iElem,bb)                            & !Excretion
                                            - LimnoParam%tPMortBac(iElem,bb)                            & !Mortality
                                            - LimnoParam%tPSetBac(HydroParam%ElSmallm(iElem),iElem,bb)  & !Sedimentation
                                            -SUM(LimnoParam%tPGrzBentBac(iElem,:,bb))                   !Grazing by Benthos
            
            LimnoParam%sDBacS(iElem,bb) = MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,1)+LimnoParam%sDBacSP(iElem,bb)))
            LimnoParam%sCBacS(iElem,bb) = MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,2)+LimnoParam%sCBacSP(iElem,bb)))
            LimnoParam%sNBacS(iElem,bb) = MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,3)+LimnoParam%sNBacSP(iElem,bb)))
            LimnoParam%sPBacS(iElem,bb) = MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,4)+LimnoParam%sPBacSP(iElem,bb)))
            
        EndDo !Loop Cell
    EndDo !Loop Functional Group
      

Return
End
