Subroutine Bacterioplankton(HydroParam,MeshParam,LimnoParam,dt,dtday)

    ! Bacterioplankton Modelling in the Water Column
    
    ! List of Modifications:
    !   24.07.2015: Routine Implementation      (Carlos Ruberto)
    !   25.07.2015: Routine Validation      (Carlos Ruberto)
    ! Programmer: Carlos Ruberto
    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars

    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    Integer:: index,iElem,iLayer,bb,dom
    Real:: wDSetBac_Layer,wCSetBac_Layer,wNSetBac_Layer,wPSetBac_Layer
    Real:: dt,dtday
    Real:: V
    Real:: NearZero = 1e-10
    
    Do bb = 1, LimnoParam%nbac
        Do iElem = 1,MeshParam%nElem
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                !-------------------------------------------------------------------------------------!
                !                           Bacterioplankton Processes in Water                       !
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
		        LimnoParam%aTLimMinBac = LimnoParam%cThetaBac(bb) ** (LimnoParam%sDTempW(iLayer,iElem)-20.)
			
			    !pH
			    If (LimnoParam%spHW(iLayer,iElem).GT.LimnoParam%pHmax) Then
				    LimnoParam%apHLimMinBac = exp(LimnoParam%pHmax - LimnoParam%spHW(iLayer,iElem))
			    ElseIf (LimnoParam%spHW(iLayer,iElem).LT.LimnoParam%pHmin) Then
				    LimnoParam%apHLimMinBac = exp(LimnoParam%spHW(iLayer,iElem) - LimnoParam%pHmin)
			    Else
				    LimnoParam%apHLimMinBac = 1.0
			    EndIf
			
			    !O2
			    If (LimnoParam%abac(bb)==LimnoParam%AEROBICA) Then
			        LimnoParam%aO2LimMinBac = LimnoParam%sO2W(iLayer,iElem)/(LimnoParam%sO2W(iLayer,iElem)+LimnoParam%hO2Bac+NearZero)
			    ElseIf (LimnoParam%abac(bb)==LimnoParam%ANAEROBICA) Then
			        LimnoParam%aO2LimMinBac = LimnoParam%kbacAn/(LimnoParam%sO2W(iLayer,iElem)+LimnoParam%kbacAn+NearZero)
			    Else
			        LimnoParam%aO2LimMinBac = (1.-LimnoParam%fbacAn)*LimnoParam%sO2W(iLayer,iElem)/(LimnoParam%sO2W(iLayer,iElem)+LimnoParam%hO2Bac+NearZero) + LimnoParam%fbacAn*(LimnoParam%kbacAn/(LimnoParam%sO2W(iLayer,iElem)+LimnoParam%kbacAn+NearZero))
			    EndIf
			
			    !Internal nutrients
			    LimnoParam%aNutLimMinBac = Min(LimnoParam%rNDBacW(iLayer,iElem,bb)/LimnoParam%cNDBacRef(bb) , LimnoParam%rPDBacW(iLayer,iElem,bb)/LimnoParam%cPDBacRef(bb))
			
			    !Substrate availability : Saturation - Pref * FOOD -> oDFood"Predador""Presa"(pred,food)
                If (LimnoParam%InclMatOrg==1) Then
    	            Do dom = 1, LimnoParam%ndom
    		            LimnoParam%oDFoodBacDom(dom) = LimnoParam%cPrefBacDom(bb,dom) * LimnoParam%sDDomW(iLayer,iElem,dom) 
                    EndDo
                Else
                    LimnoParam%oDFoodBacDom(dom) = 0.
                EndIf
                
                
    	        LimnoParam%oDFoodBac = sum(LimnoParam%oDFoodBacDom(:))
			    LimnoParam%aDFoodSatBac = LimnoParam%oDFoodBac / (LimnoParam%oDFoodBac + LimnoParam%hAssBac(bb)+NearZero)
                
			    ! 1.2 Uptake
			    LimnoParam%wDAssBac(iLayer,iElem,bb) = LimnoParam%MuBac(bb) * LimnoParam%apHLimMinBac * LimnoParam%aO2LimMinBac * LimnoParam%aTLimMinBac * LimnoParam%aNutLimMinBac * LimnoParam%aDFoodSatBac * LimnoParam%sDBacW(iLayer,iElem,bb) 

                ! Grazing Function x Assimilation = G(D) x Fgrz
                Do dom = 1, LimnoParam%ndom
		            LimnoParam%fDAssBacDom(dom) = LimnoParam%oDFoodBacDom(dom) / LimnoParam%oDFoodBac
                    LimnoParam%wDAssBacDom(iLayer,iElem,bb,dom) = LimnoParam%fDAssBacDom(dom) * LimnoParam%wDAssBac(iLayer,iElem,bb)
		            LimnoParam%wCAssBacDom(iLayer,iElem,bb,dom) = LimnoParam%wDAssBacDom(iLayer,iElem,bb,dom) * LimnoParam%rCDDomW(iLayer,iElem,dom)
		            LimnoParam%wNAssBacDom(iLayer,iElem,bb,dom) = LimnoParam%wDAssBacDom(iLayer,iElem,bb,dom) * LimnoParam%rNDDomW(iLayer,iElem,dom)
	                LimnoParam%wPAssBacDom(iLayer,iElem,bb,dom) = LimnoParam%wDAssBacDom(iLayer,iElem,bb,dom) * LimnoParam%rPDDomW(iLayer,iElem,dom)
                    LimnoParam%wSiAssBacDom(iLayer,iElem,bb,dom) = LimnoParam%wDAssBacDom(iLayer,iElem,bb,dom) * LimnoParam%cSiDDomRef(dom)
                EndDo
                ! Total assimilated nutrients of preys
                LimnoParam%wCAssBac(iLayer,iElem,bb) =  SUM(LimnoParam%wCAssBacDom(iLayer,iElem,bb,:)) 
                LimnoParam%wNAssBac(iLayer,iElem,bb) =  SUM(LimnoParam%wNAssBacDom(iLayer,iElem,bb,:))
                LimnoParam%wPAssBac(iLayer,iElem,bb) =  SUM(LimnoParam%wPAssBacDom(iLayer,iElem,bb,:))
                LimnoParam%wSiAssBac(iLayer,iElem,bb) = SUM(LimnoParam%wSiAssBacDom(iLayer,iElem,bb,:))

                
    		    !-------------------------------------------------------------------------!
			    ! (2) Respiration BAC => CO2 e CH4										  !
			    ! 2.1 basal respiration - Respb
                LimnoParam%kCorResbBac = LimnoParam%kResbBac(bb) * LimnoParam%aTLimMinBac
                LimnoParam%wDRespbBac = LimnoParam%kCorResbBac * LimnoParam%sDBacW(iLayer,iElem,bb)         
			    ! 2.2 Parcela de perda por respiracao devido atividade metabolica (ERSEM) - Respa
                !wDRespaBac = fBacLoss(bb) * wDAssBac(iLayer,iElem,bb) 		
			    ! 2.3 Total Respiration -> basal + atividade metabolica(mg/m3/s)
                LimnoParam%wDRespBac(iLayer,iElem,bb) = LimnoParam%wDRespbBac !+ wDRespaBac
                
                !  C_Respiration
                LimnoParam%wCRespBac(iLayer,iElem,bb) = LimnoParam%rCDBacW(iLayer,iElem,bb) / LimnoParam%cCDBacRef(bb) * LimnoParam%kResbBac(bb) * LimnoParam%aTLimMinBac*LimnoParam%sCBacW(iLayer,iElem,bb)

                LimnoParam%wCRespBac2CO2(iLayer,iElem,bb) =     LimnoParam%fBac2CO2(bb)  * LimnoParam%wCRespBac(iLayer,iElem,bb) 
				LimnoParam%wCRespBac2CH4(iLayer,iElem,bb) = (1.-LimnoParam%fBac2CO2(bb)) * LimnoParam%wCRespBac(iLayer,iElem,bb) 
                
                !  N_excretion
                LimnoParam%wNExcrBac(iLayer,iElem,bb) = LimnoParam%rNDBacW(iLayer,iElem,bb) / LimnoParam%cNDBacRef(bb) * LimnoParam%kResbBac(bb) * LimnoParam%aTLimMinBac*LimnoParam%sNBacW(iLayer,iElem,bb)
                !  P_excretion
                LimnoParam%wPExcrBac(iLayer,iElem,bb) = LimnoParam%rPDBacW(iLayer,iElem,bb) / LimnoParam%cPDBacRef(bb) * LimnoParam%kResbBac(bb) * LimnoParam%aTLimMinBac*LimnoParam%sPBacW(iLayer,iElem,bb)
                
                
			    !-------------------------------------------------------------------------!
			    ! (3) Mortality      												  !
                LimnoParam%kCorMortBac = LimnoParam%aTLimMinBac * LimnoParam%kMortBac(bb)
                LimnoParam%wDMortBac(iLayer,iElem,bb) = LimnoParam%kCorMortBac * LimnoParam%sDBacW(iLayer,iElem,bb)
                LimnoParam%wCMortBac(iLayer,iElem,bb) = LimnoParam%rCDBacW(iLayer,iElem,bb) * LimnoParam%wDMortBac(iLayer,iElem,bb)
                LimnoParam%wNMortBac(iLayer,iElem,bb) = LimnoParam%rNDBacW(iLayer,iElem,bb) * LimnoParam%wDMortBac(iLayer,iElem,bb)
                LimnoParam%wPMortBac(iLayer,iElem,bb) = LimnoParam%rPDBacW(iLayer,iElem,bb) * LimnoParam%wDMortBac(iLayer,iElem,bb)

			    !-------------------------------------------------------------------------!
			    ! (4) Uptake/Release PO4,NH4				                              !
			    ! 4.1 PO4 Uptake/Release (CAEDYM)
			    If (LimnoParam%wDAssBac(iLayer,iElem,bb)*LimnoParam%cPDBacRef(bb).GT.LimnoParam%wPAssBac(iLayer,iElem,bb)) Then
			        LimnoParam%wPUptBac(iLayer,iElem,bb) = (LimnoParam%wDAssBac(iLayer,iElem,bb)*LimnoParam%cPDBacRef(bb) - LimnoParam%wPAssBac(iLayer,iElem,bb)) - LimnoParam%wDRespBac(iLayer,iElem,bb) * LimnoParam%rPDBacW(iLayer,iElem,bb)
			    Else
			        LimnoParam%wPUptBac(iLayer,iElem,bb) = - LimnoParam%wDRespBac(iLayer,iElem,bb) * LimnoParam%rPDBacW(iLayer,iElem,bb)
			    EndIf

			    ! 4.2 NH4 Uptake/Release (CAEDYM)		    
		        If (LimnoParam%wDAssBac(iLayer,iElem,bb)*LimnoParam%cNDBacRef(bb).GT.LimnoParam%wNAssBac(iLayer,iElem,bb)) Then
			        LimnoParam%wNUptBac(iLayer,iElem,bb) = (LimnoParam%wDAssBac(iLayer,iElem,bb)*LimnoParam%cNDBacRef(bb) - LimnoParam%wNAssBac(iLayer,iElem,bb)) - LimnoParam%wDRespBac(iLayer,iElem,bb) * LimnoParam%rNDBacW(iLayer,iElem,bb)
			    Else
			        LimnoParam%wNUptBac(iLayer,iElem,bb) = - LimnoParam%wDRespBac(iLayer,iElem,bb) * LimnoParam%rNDBacW(iLayer,iElem,bb)
                EndIf
                
			    !-------------------------------------------------------------------------!
			    ! (5) Sedimentation				                              !
                If ( iLayer == HydroParam%ElCapitalM(iElem) ) Then        ! Bottom layer or intermediary layer has contribuition of the up layer
                    wDSetBac_Layer = LimnoParam%tDSetBac(iLayer,iElem,bb)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wCSetBac_Layer = wDSetBac_Layer * LimnoParam%rCDBacW(iLayer,iElem,bb)
                    wNSetBac_Layer = wDSetBac_Layer * LimnoParam%rNDBacW(iLayer,iElem,bb)
                    wPSetBac_Layer = wDSetBac_Layer * LimnoParam%rPDBacW(iLayer,iElem,bb)
                Else
                    wDSetBac_Layer = LimnoParam%tDSetBac(iLayer,iElem,bb)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem)) - LimnoParam%tDSetBac(iLayer+1,iElem,bb)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer+1,iElem))
                    wCSetBac_Layer = wDSetBac_Layer * LimnoParam%rCDBacW(iLayer,iElem,bb)
                    wNSetBac_Layer = wDSetBac_Layer * LimnoParam%rNDBacW(iLayer,iElem,bb)
                    wPSetBac_Layer = wDSetBac_Layer * LimnoParam%rPDBacW(iLayer,iElem,bb)
                Endif
                
                !-------------------------------------------------------------------------------------!
                !                   Source/Sink Term for Bacterioplankton						      !
                !-------------------------------------------------------------------------------------!
                
                HydroParam%dVarEst(iLayer,iElem,1) =  LimnoParam%wDAssBac(iLayer,iElem,bb)              & !Assimilation  
                                            - LimnoParam%wDRespBac(iLayer,iElem,bb)                     & !Respiration
                                            - LimnoParam%wDMortBac(iLayer,iElem,bb)                     & !Mortality
                                            - wDSetBac_Layer                                            & !Sedimentation
                                            - Sum(LimnoParam%wDGrzZooBac(iLayer,iElem,:,bb))            !Grazing by Zoops

                HydroParam%dVarEst(iLayer,iElem,2) =  LimnoParam%wCAssBac(iLayer,iElem,bb)              & !Assimilation  
                                            - LimnoParam%wCRespBac(iLayer,iElem,bb)                     & !Respiration
                                            - LimnoParam%wCMortBac(iLayer,iElem,bb)                     & !Mortality
                                            - wCSetBac_Layer                                            & !Sedimentation
                                            - Sum(LimnoParam%wCGrzZooBac(iLayer,iElem,:,bb))            !Grazing by Zoops
                
                HydroParam%dVarEst(iLayer,iElem,3) =  LimnoParam%wNAssBac(iLayer,iElem,bb)              & !Assimilation  
                                            + LimnoParam%wNUptBac(iLayer,iElem,bb)                      & !Uptake/Release
                                            - LimnoParam%wNExcrBac(iLayer,iElem,bb)                     & !Excretion
                                            - LimnoParam%wNMortBac(iLayer,iElem,bb)                     & !Mortality
                                            - wNSetBac_Layer                                            & !Sedimentation
                                            - Sum(LimnoParam%wNGrzZooBac(iLayer,iElem,:,bb))            !Grazing by Zoops

                HydroParam%dVarEst(iLayer,iElem,4) =  LimnoParam%wPAssBac(iLayer,iElem,bb)              & !Assimilation  
                                            + LimnoParam%wPUptBac(iLayer,iElem,bb)                      & !Uptake/Release
                                            - LimnoParam%wPExcrBac(iLayer,iElem,bb)                     & !Excretion
                                            - LimnoParam%wPMortBac(iLayer,iElem,bb)                     & !Mortality
                                            - wPSetBac_Layer                                            & !Sedimentation
                                            - Sum(LimnoParam%wPGrzZooBac(iLayer,iElem,:,bb))            !Grazing by Zoops
                
                
            EndDo !Loop Layer
        EndDo !Loop Cell
      
        !-------------------------------------------------------------------------------------!
        !                   Solver Transport Equation for Bacterioplankton                    !
        !-------------------------------------------------------------------------------------!
        index = 0
        
        If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDBacW(:,:,bb),LimnoParam%sDBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCBacW(:,:,bb),LimnoParam%sCBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNBacW(:,:,bb),LimnoParam%sNBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPBacW(:,:,bb),LimnoParam%sPBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam)
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDBacW(:,:,bb),LimnoParam%sDBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCBacW(:,:,bb),LimnoParam%sCBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNBacW(:,:,bb),LimnoParam%sNBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPBacW(:,:,bb),LimnoParam%sPBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
        ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDBacW(:,:,bb),LimnoParam%sDBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCBacW(:,:,bb),LimnoParam%sCBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNBacW(:,:,bb),LimnoParam%sNBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPBacW(:,:,bb),LimnoParam%sPBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam)
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDBacW(:,:,bb),LimnoParam%sDBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCBacW(:,:,bb),LimnoParam%sCBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNBacW(:,:,bb),LimnoParam%sNBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPBacW(:,:,bb),LimnoParam%sPBacWP(:,:,bb),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
        EndiF     
        
    EndDo !Loop Functional Group

    
Return
End
