Subroutine PhytoWater(HydroParam,MeshParam,MeteoParam,LimnoParam,dt,dtday)

    
    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars
    Use Meteorological
    Use uTempFunction

    Implicit None
    type(MeshGridParam) :: MeshParam
    type(MeteorologicalParam) :: MeteoParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
!   Declaracao das variaveis da subrotina
    Integer:: index,iElem,iLayer,gg
    Real:: wDSetPhyt_Layer,wCSetPhyt_Layer,wNSetPhyt_Layer,wPSetPhyt_Layer,tDResusPhyt_bottom,tCResusPhyt_bottom,tNResusPhyt_bottom,tPResusPhyt_bottom
    Real:: V
    Real:: dt,dtday
    Real:: NearZero = 1e-10
    
    Do gg = 1, LimnoParam%nphy
        Do iElem = 1,MeshParam%nElem
            Do iLayer = HydroParam%ElCapitalM(iElem),HydroParam%ElSmallm(iElem),-1
                ! Light
                If ( HydroParam%ElSmallm(iElem) == HydroParam%ElCapitalM(iElem) ) Then        ! Only One Vertical Layer
                    ! light_at_the_surface,in FABM it's at the top of the current layer (W/m² PAR)
                    LimnoParam%uLPAR0(iLayer) = MeshParam%CREDV(iElem)*(1.-HydroParam%ALB)*LimnoParam%fPAR*MeteoParam%SolarRad(iElem)*exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer+1,iElem)))
                    ! light_at_the_bottom,in FABM it's at the bottom of the current layer  (W/m² PAR)
                    LimnoParam%aLPARBot(iLayer) = LimnoParam%uLPAR0(iLayer)*exp(-LimnoParam%aExtCoef(iLayer,iElem)*1.)
                Else
                    If ( iLayer == HydroParam%ElSmallm(iElem) ) Then !Bottom layer
                        ! light_at_the_surface,in FABM it's at the top of the current layer (W/m² PAR)
                        LimnoParam%uLPAR0(iLayer) = LimnoParam%aLPARBot(iLayer+1)
                        ! light_at_the_bottom,in FABM it's at the bottom of the current layer  (W/m² PAR)
                        LimnoParam%aLPARBot(iLayer) = LimnoParam%uLPAR0(iLayer)*exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer,iElem)))
                    Elseif ( iLayer == HydroParam%ElCapitalM(iElem) ) Then !Surface layer
                        ! light_at_the_surface,in FABM it's at the top of the current layer (W/m² PAR)
                        LimnoParam%uLPAR0(iLayer) = MeshParam%CREDV(iElem)*(1.-HydroParam%ALB)*LimnoParam%fPAR*MeteoParam%SolarRad(iElem)*exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer+1,iElem)))
                        ! light_at_the_bottom,in FABM it's at the bottom of the current layer  (W/m² PAR)
                        LimnoParam%aLPARBot(iLayer) = LimnoParam%uLPAR0(iLayer)*exp(-LimnoParam%aExtCoef(iLayer,iElem)*Max(1.,HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer,iElem)))
                    Else !Intermediary layer
                        ! light_at_the_surface,in FABM it's at the top of the current layer (W/m² PAR)
                        LimnoParam%uLPAR0(iLayer) = LimnoParam%aLPARBot(iLayer+1)
                        ! light_at_the_bottom,in FABM it's at the bottom of the current layer  (W/m² PAR)
                        LimnoParam%aLPARBot(iLayer) = LimnoParam%uLPAR0(iLayer)*exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer,iElem)))
                    EndIf
                Endif
            EndDo
            
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                !-------------------------------------------------------------------------------------!
                !                           Phytoplankton Processes in Water					      !
                !-------------------------------------------------------------------------------------!
                !-Phyto-------------------------------------------------------------------------------!
                ! 1. Primary Production/Nutrient Uptake                                               !
                ! 2. Respiration/Excretion                                                            !
                ! 3. Lise (Mortality)                                                                 !
                ! 4. Sedimentation                                                                    !
                ! 5. Resuspension                                                                     !
                !-------------------------------------------------------------------------------------!
                !-------------------------------------------------------------------------!
			    !1) Primary Production
                ! Temperature function for pelagic phytoplankton
                LimnoParam%uFunTmPhyt = uFunTmBio(LimnoParam%sDTempW(iLayer,iElem),LimnoParam%cSigTmPhyt(gg),LimnoParam%cTmOptPhyt(gg)) 
                !-----------------------------------------------------------------------
                !  Light functions for pelagic phytoplankton
                !-----------------------------------------------------------------------
                ! light limitation function for each group, each group has two different types
                If (LimnoParam%LightMethodPhyt == 0) Then ! case 1, without photoinhibition,Chalker (1980) model
                    ! half-satuation light intensity at current temperature
                    LimnoParam%uhLPhyt = LimnoParam%hLRefPhyt(gg) * LimnoParam%uFunTmPhyt
                    ! light limitation function for green algae, no photo-inhibition
                    If (HydroParam%ElSmallm(iElem)== HydroParam%ElCapitalM(iElem)) Then ! 2D - Production in the first meter
                        LimnoParam%aLLimPhyt = 1.0 / (LimnoParam%aExtCoef(iLayer,iElem) * 1.) * log((1.0 + LimnoParam%uLPAR0(iLayer) / (LimnoParam%uhLPhyt+NearZero)) / (1.0 + LimnoParam%aLPARBot(iLayer) / (LimnoParam%uhLPhyt+NearZero))) !função Lehman para luz (integração)
                    Else
                        LimnoParam%aLLimPhyt = 1.0 / (LimnoParam%aExtCoef(iLayer,iElem) * Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))) * log((1.0 + LimnoParam%uLPAR0(iLayer) / (LimnoParam%uhLPhyt+NearZero)) / (1.0 + LimnoParam%aLPARBot(iLayer) / (LimnoParam%uhLPhyt+NearZero))) !função Lehman para luz (integração)
                    EndIf
                ElseIf (LimnoParam%LightMethodPhyt == 1) Then !  case 2, Klepper et al. (1988) / Ebenhoh et al. (1997) model.
                    LimnoParam%uhLPhyt = LimnoParam%cLOptRefPhyt(gg) * LimnoParam%uFunTmPhyt
                    If (HydroParam%ElSmallm(iElem)== HydroParam%ElCapitalM(iElem)) Then ! 2D - Production in the first meter
                        LimnoParam%aLLimPhyt = exp(1.0) /(LimnoParam%aExtCoef(iLayer,iElem) * 1.) *(exp(- LimnoParam%aLPARBot(iLayer) /(LimnoParam%uhLPhyt+NearZero)) - exp(- LimnoParam%uLPAR0(iLayer) /(LimnoParam%uhLPhyt+NearZero)))
                    Else
                        LimnoParam%aLLimPhyt = exp(1.0) /(LimnoParam%aExtCoef(iLayer,iElem) * Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))) *(exp(- LimnoParam%aLPARBot(iLayer) /(LimnoParam%uhLPhyt+NearZero)) - exp(- LimnoParam%uLPAR0(iLayer) /(LimnoParam%uhLPhyt+NearZero)))
                    EndIf
                Else
                    Stop 'Light limitation function for Diatoms does not exist'
                EndIf
                ! growth_rate_at_current_light_AND_temp.
                LimnoParam%aMuTmLPhyt = LimnoParam%ufDay *(1.0 - LimnoParam%afCovSurfVeg)* LimnoParam%aLLimPhyt * LimnoParam%uFunTmPhyt*LimnoParam%cMuMaxPhyt(gg)

                !-----------------------------------------------------------------------
                !  Nutrient limitation functions
                !-----------------------------------------------------------------------
                !  Nutrient functions for diatom
                !  Droop_function(C)_for_Phytoplankton
                LimnoParam%aCLimPhyt = max(0.0,(1.0 - LimnoParam%cCDPhytMin(gg) / (LimnoParam%rCDPhytW(iLayer,iElem,gg) + NearZero)) * LimnoParam%cCDPhytMax(gg) /(LimnoParam%cCDPhytMax(gg) - LimnoParam%cCDPhytMin(gg)))
                !  Droop_function(P)_for_Phytoplankton
                LimnoParam%aPLimPhyt = max(0.0,(1.0 - LimnoParam%cPDPhytMin(gg) / (LimnoParam%rPDPhytW(iLayer,iElem,gg) + NearZero)) * LimnoParam%cPDPhytMax(gg) /(LimnoParam%cPDPhytMax(gg) - LimnoParam%cPDPhytMin(gg)))
                !  Droop_function(N)_for_Phytoplankton
                LimnoParam%aNLimPhyt = max(0.0,(1.0 - LimnoParam%cNDPhytMin(gg) / (LimnoParam%rNDPhytW(iLayer,iElem,gg) + NearZero)) * LimnoParam%cNDPhytMax(gg) /(LimnoParam%cNDPhytMax(gg) - LimnoParam%cNDPhytMin(gg)))
                !  silica_dependence_of_growth_rate
                LimnoParam%aSiLimPhyt = LimnoParam%sSiO2W(iLayer,iElem) /(LimnoParam%hSiAssDiat(gg) + LimnoParam%sSiO2W(iLayer,iElem))
                !  nutrient_limitation_function_of_Diatom
                !   aNutLimDiat = min(aPLimDiat,aNLimDiat,aSiLimDiat)  ! v5.09
                LimnoParam%aNutLimPhyt = min(LimnoParam%aPLimPhyt,LimnoParam%aNLimPhyt)   ! pl613
                !-----------------------------------------------------------------------
                !  Algae growth_DW
                !-----------------------------------------------------------------------
                !  growth_rate_Phytoplankton
                
                LimnoParam%aMuPhyt = LimnoParam%aNutLimPhyt * LimnoParam%aMuTmLPhyt
                !  assimilation_Algae
                LimnoParam%wDAssPhyt(iLayer,iElem,gg) = LimnoParam%aMuPhyt*LimnoParam%sDPhytW(iLayer,iElem,gg)
                !-------------------------------------------------------------------------------------------------------------
                !  Algae uptake_C
                !-------------------------------------------------------------------------------------------------------------
                ! maximum_P_uptake_rate_of_Algae,corrected_for_P/D_ratio
                LimnoParam%aVCUptMaxCorPhyt = max(0.0,LimnoParam%cVCUptMaxPhyt(gg) * LimnoParam%aLLimPhyt * LimnoParam%uFunTmPhyt *(LimnoParam%cCDPhytMax(gg) - LimnoParam%rCDPhytW(iLayer,iElem,gg))&
                                     &/(LimnoParam%cCDPhytMax(gg) - LimnoParam%cCDPhytMin(gg)))
                ! C_uptake_rate_of_Algae
                LimnoParam%aVCUptPhyt = LimnoParam%sDICW(iLayer,iElem) * LimnoParam%aVCUptMaxCorPhyt /(LimnoParam%cVCUptMaxPhyt(gg) / LimnoParam%cAffCUptPhyt(gg) + LimnoParam%sDICW(iLayer,iElem))
                ! C_uptake_Algae
                LimnoParam%wCUptPhyt(iLayer,iElem,gg) = LimnoParam%aVCUptPhyt * LimnoParam%sDPhytW(iLayer,iElem,gg)
                !-------------------------------------------------------------------------------------------------------------
                !  Algae uptake_P
                !-------------------------------------------------------------------------------------------------------------
                ! maximum_P_uptake_rate_of_Algae,corrected_for_P/D_ratio
                LimnoParam%aVPUptMaxCorPhyt = max(0.0,LimnoParam%cVPUptMaxPhyt(gg) * LimnoParam%uFunTmPhyt *(LimnoParam%cPDPhytMax(gg) - LimnoParam%rPDPhytW(iLayer,iElem,gg))&
                                     &/(LimnoParam%cPDPhytMax(gg) - LimnoParam%cPDPhytMin(gg)))
                ! P_uptake_rate_of_Algae
                LimnoParam%aVPUptPhyt = LimnoParam%sPO4W(iLayer,iElem) * LimnoParam%aVPUptMaxCorPhyt /(LimnoParam%cVPUptMaxPhyt(gg) / LimnoParam%cAffPUptPhyt(gg) + LimnoParam%sPO4W(iLayer,iElem))
                ! P_uptake_Algae
                LimnoParam%wPUptPhyt(iLayer,iElem,gg) = LimnoParam%aVPUptPhyt * LimnoParam%sDPhytW(iLayer,iElem,gg)
                !-------------------------------------------------------------------------------------------------------------
                !  Algae uptake_N
                !-------------------------------------------------------------------------------------------------------------
                ! maximum_N_uptake_rate_of_Algae,corrected_for_N/D_ratio
                LimnoParam%aVNUptMaxCorPhyt = max(0.0,LimnoParam%cVNUptMaxPhyt(gg) * LimnoParam%uFunTmPhyt *(LimnoParam%cNDPhytMax(gg) - LimnoParam%rNDPhytW(iLayer,iElem,gg))&
                                     &/(LimnoParam%cNDPhytMax(gg) - LimnoParam%cNDPhytMin(gg)))
                !  half-sat._NDissW_for_uptake_by_Algae
                LimnoParam%ahNUptPhyt = LimnoParam%aVNUptMaxCorPhyt / LimnoParam%cAffNUptPhyt(gg)
                ! N_uptake_rate_of_Algae
                LimnoParam%aVNUptPhyt = (LimnoParam%sNO3W(iLayer,iElem) + LimnoParam%sNH4W(iLayer,iElem)) * LimnoParam%aVNUptMaxCorPhyt /(LimnoParam%ahNUptPhyt + LimnoParam%sNO3W(iLayer,iElem) + LimnoParam%sNH4W(iLayer,iElem))
                ! N_uptake_Algae
                LimnoParam%wNUptPhyt(iLayer,iElem,gg) = LimnoParam%aVNUptPhyt * LimnoParam%sDPhytW(iLayer,iElem,gg)
                !-----------------------------------------------------------------------
                !  NH4 exchange with abiotic module
                !-----------------------------------------------------------------------
                !  fraction_ammonium_uptake_by_Algae
                LimnoParam%afNH4UptPhyt = LimnoParam%sNH4W(iLayer,iElem) * LimnoParam%sNO3W(iLayer,iElem) /((LimnoParam%ahNUptPhyt + LimnoParam%sNH4W(iLayer,iElem)) *(LimnoParam%ahNUptPhyt + LimnoParam%sNO3W(iLayer,iElem))) + LimnoParam%sNH4W(iLayer,iElem) &
                & * LimnoParam%ahNUptPhyt /((LimnoParam%sNH4W(iLayer,iElem) + LimnoParam%sNO3W(iLayer,iElem)) *(LimnoParam%ahNUptPhyt + LimnoParam%sNO3W(iLayer,iElem)))
                !  ammonium_uptake_by_Algae
                LimnoParam%wNUptNH4Phyt(iLayer,iElem,gg) = LimnoParam%afNH4UptPhyt * LimnoParam%wNUptPhyt(iLayer,iElem,gg)
                !-----------------------------------------------------------------------
                !  NO3 exchange with abiotic module
                !-----------------------------------------------------------------------
                !   nitrate_uptake_by_Algae
                LimnoParam%wNUptNO3Phyt(iLayer,iElem,gg) = LimnoParam%wNUptPhyt(iLayer,iElem,gg) - LimnoParam%wNUptNH4Phyt(iLayer,iElem,gg)
                
                !2) Algae respiration_DW and excretion of nutrients
                !  temp._corrected_respiration_constant_of_Algae
                LimnoParam%ukDRespbTmPhyt = LimnoParam%kDRespbPhyt(gg) * LimnoParam%uFunTmPhyt
                !  basal respiration_of_Algae_in_water
                LimnoParam%wDRespbPhyt = LimnoParam%ukDRespbTmPhyt * LimnoParam%sDPhytW(iLayer,iElem,gg) ! Used in O2 exchange
                !  Photorespiration_of_Algae_in_water (ERSEM, Polimene et al.(2006)) - alfap -> fracao da producao primaria bruta (PPG)
                LimnoParam%wDResppPhyt = LimnoParam%alfap(gg) * LimnoParam%wDAssPhyt(iLayer,iElem,gg)
                !  Total respiration_of_Algae_in_water
                LimnoParam%wDRespPhyt(iLayer,iElem,gg) =  LimnoParam%wDRespbPhyt + LimnoParam%wDResppPhyt
                !DIC Respiration
                
                !  C_excretion_Algae_in_water (Exudation)
                !   wPExcrDiatW = rPDDiatW /(self%cPDDiatMin + rPDDiatW) * rPDDiatW * wDRespDiatW ! v5.09
                LimnoParam%wCExcrPhytW(iLayer,iElem,gg) = (2.0 * LimnoParam%rCDPhytW(iLayer,iElem,gg)) /(LimnoParam%cCDPhytMax(gg) + LimnoParam%rCDPhytW(iLayer,iElem,gg)) * LimnoParam%rCDPhytW(iLayer,iElem,gg) * LimnoParam%wDRespPhyt(iLayer,iElem,gg)
                !  P_Exudation_Algae_in_water
                !   wPExcrDiatW = rPDDiatW /(self%cPDDiatMin + rPDDiatW) * rPDDiatW * wDRespDiatW ! v5.09
                LimnoParam%wPExcrPhytW(iLayer,iElem,gg) = (2.0 * LimnoParam%rPDPhytW(iLayer,iElem,gg)) /(LimnoParam%cPDPhytMax(gg) + LimnoParam%rPDPhytW(iLayer,iElem,gg)) * LimnoParam%rPDPhytW(iLayer,iElem,gg) * LimnoParam%wDRespPhyt(iLayer,iElem,gg)
                !  N_Exudation_Algae_in_water
                !   wNExcrDiatW = rNDDiatW /(self%cNDDiatMin + rNDDiatW) * rNDDiatW * wDRespDiatW ! V5.09
                LimnoParam%wNExcrPhytW(iLayer,iElem,gg) = (2.0 * LimnoParam%rNDPhytW(iLayer,iElem,gg)) /(LimnoParam%cNDPhytMax(gg) + LimnoParam%rNDPhytW(iLayer,iElem,gg)) * LimnoParam%rNDPhytW(iLayer,iElem,gg) * LimnoParam%wDRespPhyt(iLayer,iElem,gg) ! pl613

                !3) Lise_Algae_in_water
                LimnoParam%wDLisPhyt(iLayer,iElem,gg) = LimnoParam%cLisPhytW(gg) * LimnoParam%sDPhytW(iLayer,iElem,gg)
   				LimnoParam%wCLisPhyt(iLayer,iElem,gg) = LimnoParam%wDLisPhyt(iLayer,iElem,gg) * LimnoParam%rCDPhytW(iLayer,iElem,gg)
				LimnoParam%wNLisPhyt(iLayer,iElem,gg) = LimnoParam%wDLisPhyt(iLayer,iElem,gg) * LimnoParam%rNDPhytW(iLayer,iElem,gg)
                LimnoParam%wPLisPhyt(iLayer,iElem,gg) = LimnoParam%wDLisPhyt(iLayer,iElem,gg) * LimnoParam%rPDPhytW(iLayer,iElem,gg)
				If (LimnoParam%aphy(gg)==LimnoParam%MDIAT.or.LimnoParam%aphy(gg)==LimnoParam%FDIAT) Then
					LimnoParam%wSiLisPhyt(iLayer,iElem,gg) = LimnoParam%wDLisPhyt(iLayer,iElem,gg) * LimnoParam%rSiDPhytW(iLayer,iElem,gg)					
				Else
					LimnoParam%wSiLisPhyt(iLayer,iElem,gg) = 0.					
                EndIf
                
                !4) Algae Sedimentation
                If ( iLayer == HydroParam%ElCapitalM(iElem) ) Then        ! Bottom layer or intermediary layer has contribuition of the up layer
                    wDSetPhyt_Layer = LimnoParam%tDSetPhyt(iLayer,iElem,gg)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wCSetPhyt_Layer = wDSetPhyt_Layer * LimnoParam%rCDPhytW(iLayer,iElem,gg)
                    wNSetPhyt_Layer = wDSetPhyt_Layer * LimnoParam%rNDPhytW(iLayer,iElem,gg)
                    wPSetPhyt_Layer = wDSetPhyt_Layer * LimnoParam%rPDPhytW(iLayer,iElem,gg)
                Else
                    wDSetPhyt_Layer = LimnoParam%tDSetPhyt(iLayer,iElem,gg)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem)) - LimnoParam%tDSetPhyt(iLayer+1,iElem,gg)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer+1,iElem))
                    wCSetPhyt_Layer = wDSetPhyt_Layer * LimnoParam%rCDPhytW(iLayer,iElem,gg)
                    wNSetPhyt_Layer = wDSetPhyt_Layer * LimnoParam%rNDPhytW(iLayer,iElem,gg)
                    wPSetPhyt_Layer = wDSetPhyt_Layer * LimnoParam%rPDPhytW(iLayer,iElem,gg)
                Endif
            
                !5) Algae Resuspension 
                If (iLayer == HydroParam%ElSmallm(iElem)) Then
                    tDResusPhyt_bottom = LimnoParam%tDResusPhyt(iElem,gg)
                    tCResusPhyt_bottom = tDResusPhyt_bottom * LimnoParam%rCDPhytW(iLayer,iElem,gg)
                    tNResusPhyt_bottom = tDResusPhyt_bottom * LimnoParam%rNDPhytW(iLayer,iElem,gg)
                    tPResusPhyt_bottom = tDResusPhyt_bottom * LimnoParam%rPDPhytW(iLayer,iElem,gg)
                Else
                    tDResusPhyt_bottom = 0.
                    tCResusPhyt_bottom = 0.
                    tNResusPhyt_bottom = 0.
                    tPResusPhyt_bottom = 0.
                EndIf
                
                
           
                !-------------------------------------------------------------------------------------!
                !                               Source/Sink Term for Algae						      !
                !-------------------------------------------------------------------------------------!
                
                HydroParam%dVarEst(iLayer,iElem,1) =  LimnoParam%wDAssPhyt(iLayer,iElem,gg)                             & !Assimilation  
                                        - LimnoParam%wDRespPhyt(iLayer,iElem,gg)                                        & !Respiration
	                                    - LimnoParam%wDLisPhyt(iLayer,iElem,gg)                                         & !Lise Celular
		                                - wDSetPhyt_Layer                                                               & !Sedimentation
                                        + tDResusPhyt_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Resuspension
                                        - SUM(LimnoParam%wDGrzZooPhyt(iLayer,iElem,:,gg))                               & !Grazing by Zooplankton
                                        - SUM(LimnoParam%wDGrzFiAdPhyt(iLayer,iElem,:,gg))                              & !Grazing by Aduld Fish
                                        - SUM(LimnoParam%wDGrzFiAdPhyt(iLayer,iElem,:,gg))                              !Grazing by Juveline Fish
                
                HydroParam%dVarEst(iLayer,iElem,2) =  LimnoParam%wCUptPhyt(iLayer,iElem,gg)                             & !Nutrient Uptake  
                                        - LimnoParam%wCExcrPhytW(iLayer,iElem,gg)                                       & !Excretion (Exudation)
	                                    - LimnoParam%wCLisPhyt(iLayer,iElem,gg)                                         & !Lise Celular
		                                - wCSetPhyt_Layer                                                               & !Sedimentation
                                        + tCResusPhyt_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Resuspension
                                        - SUM(LimnoParam%wCGrzZooPhyt(iLayer,iElem,:,gg))                               & !Grazing by Zooplankton
                                        - SUM(LimnoParam%wCGrzFiAdPhyt(iLayer,iElem,:,gg))                              & !Grazing by Aduld Fish
                                        - SUM(LimnoParam%wCGrzFiAdPhyt(iLayer,iElem,:,gg))                              !Grazing by Juveline Fish
                
                HydroParam%dVarEst(iLayer,iElem,3) =  LimnoParam%wNUptPhyt(iLayer,iElem,gg)                             & !Nutrient Uptake  
                                        - LimnoParam%wNExcrPhytW(iLayer,iElem,gg)                                       & !Excretion (Exudation)
	                                    - LimnoParam%wNLisPhyt(iLayer,iElem,gg)                                         & !Lise Celular
		                                - wNSetPhyt_Layer                                                               & !Sedimentation
                                        + tNResusPhyt_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Resuspension
                                        - SUM(LimnoParam%wNGrzZooPhyt(iLayer,iElem,:,gg))                               & !Grazing by Zooplankton
                                        - SUM(LimnoParam%wNGrzFiAdPhyt(iLayer,iElem,:,gg))                              & !Grazing by Aduld Fish
                                        - SUM(LimnoParam%wNGrzFiAdPhyt(iLayer,iElem,:,gg))                              !Grazing by Juveline Fish

                HydroParam%dVarEst(iLayer,iElem,4) =  LimnoParam%wPUptPhyt(iLayer,iElem,gg)                             & !Nutrient Uptake  
                                        - LimnoParam%wPExcrPhytW(iLayer,iElem,gg)                                       & !Excretion (Exudation)
	                                    - LimnoParam%wPLisPhyt(iLayer,iElem,gg)                                         & !Lise Celular
		                                - wPSetPhyt_Layer                                                               & !Sedimentation
                                        + tPResusPhyt_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Resuspension
                                        - SUM(LimnoParam%wPGrzZooPhyt(iLayer,iElem,:,gg))                               & !Grazing by Zooplankton
                                        - SUM(LimnoParam%wPGrzFiAdPhyt(iLayer,iElem,:,gg))                              & !Grazing by Aduld Fish
                                        - SUM(LimnoParam%wPGrzFiAdPhyt(iLayer,iElem,:,gg))                              !Grazing by Juveline Fish
                        
            EndDo !Loop Layer
        EndDo !Loop Cell

        !-------------------------------------------------------------------------------------!
        !                       Solver Transport Equation for Phytoplankton               	    	      !
        !-------------------------------------------------------------------------------------!
        index = 0
    
        If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDPhytW(:,:,gg),LimnoParam%sDPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCPhytW(:,:,gg),LimnoParam%sCPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNPhytW(:,:,gg),LimnoParam%sNPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPPhytW(:,:,gg),LimnoParam%sPPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam)
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDPhytW(:,:,gg),LimnoParam%sDPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCPhytW(:,:,gg),LimnoParam%sCPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNPhytW(:,:,gg),LimnoParam%sNPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPPhytW(:,:,gg),LimnoParam%sPPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
        ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDPhytW(:,:,gg),LimnoParam%sDPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCPhytW(:,:,gg),LimnoParam%sCPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNPhytW(:,:,gg),LimnoParam%sNPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPPhytW(:,:,gg),LimnoParam%sPPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam)
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDPhytW(:,:,gg),LimnoParam%sDPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCPhytW(:,:,gg),LimnoParam%sCPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNPhytW(:,:,gg),LimnoParam%sNPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPPhytW(:,:,gg),LimnoParam%sPPhytWP(:,:,gg),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
        EndiF     
       
    EndDo !Loop Functional Group
      

Return
End
