Subroutine Macrophytes(HydroParam,MeshParam,MeteoParam,LimnoParam,dt,dtday)

    ! Macrophytes Modelling
    
    ! List of Modifications:
    !   24.07.2015: Routine Implementation      (Carlos Ruberto)
    !   25.07.2015: Routine Validation      (Carlos Ruberto)
    ! Programmer: Carlos Ruberto
    Use MeshVars
    Use Meteorological
    Use Hydrodynamic
    Use LimnologyVars
    Use uTempFunction 

    Implicit None
    type(MeshGridParam) :: MeshParam
    type(MeteorologicalParam) :: MeteoParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    Integer:: index,iElem,iLayer,mm
    Real:: dt,dtday
    Real:: d1,d2,d3,d4
    Real:: V
    Real:: NearZero = 1e-10
    
    Do mm = 1, LimnoParam%nmac
        Do iElem = 1,MeshParam%nElem
            
            !-----------------------------------------------------------------------
            !  the germination, allocation and reallocation process 
            !-----------------------------------------------------------------------
            !    setting_root_fration
            If (sum(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))/(max(1,size(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem)))) >= LimnoParam%cTmInitVeg(mm)) Then
                LimnoParam%bfRootVeg = LimnoParam%fRootVegSum(mm)
            Else
                LimnoParam%bfRootVeg = LimnoParam%fRootVegWin(mm)
            EndIf
            !-----------------------------------------------------------------------
            !  fractions of roots and shoots
            !-----------------------------------------------------------------------
            !  shoot_fraction
            LimnoParam%bfShootVeg = 1.0 - LimnoParam%bfRootVeg
            !  root_biomass
            LimnoParam%aDRootVeg = LimnoParam%bfRootVeg * LimnoParam%sDMac(iElem,mm)
            !  shoot_biomass
            LimnoParam%aDShootVeg = LimnoParam%bfShootVeg * LimnoParam%sDMac(iElem,mm)
            !  emergent_biomass
            LimnoParam%aDEmergVeg = LimnoParam%fEmergVeg(mm) * LimnoParam%aDShootVeg
            !  floating_biomass
            LimnoParam%aDFloatVeg = LimnoParam%fFloatVeg(mm) * LimnoParam%aDShootVeg
            !  submerged_fraction_of_shoot
            LimnoParam%bfSubVeg = 1.0 - LimnoParam%fFloatVeg(mm) - LimnoParam%fEmergVeg(mm)
            !  submerged_biomass
            LimnoParam%aDSubVeg = LimnoParam%bfSubVeg * LimnoParam%aDShootVeg
            !-----------------------------------------------------------------------
            !  water coverage by vegetation
            !-----------------------------------------------------------------------
            !  fraction_of_water_SURFACE_covered_by_plant_species
            LimnoParam%afCovSurfVeg = min(1.0, max(LimnoParam%aDFloatVeg / (LimnoParam%cDLayerVeg(mm) + NearZero), LimnoParam%aDEmergVeg / (LimnoParam%fEmergVeg(mm) * LimnoParam%cDCarrMac(mm) + NearZero) ))
            !  fraction_emergent_coverage
            LimnoParam%afCovEmergVeg = min(1.0, 0.01 * LimnoParam%cCovSpVeg(mm) * LimnoParam%aDEmergVeg)
            !  percent_cover
            LimnoParam%aCovVeg = min(100.0, LimnoParam%cCovSpVeg(mm) * LimnoParam%aDShootVeg)
            
            !-------------------------------------------------------------------------------------!
            !                           Macrophytes Processes                                     !
            !-------------------------------------------------------------------------------------!
            
            Do iLayer = HydroParam%ElCapitalM(iElem),HydroParam%ElSmallm(iElem),-1
                ! Light
                If ( HydroParam%ElSmallm(iElem) == HydroParam%ElCapitalM(iElem) ) Then        ! Only One Vertical Layer
                    ! light_at_the_surface,in FABM it's at the top of the current layer (W/m� PAR)
                    LimnoParam%uLPAR0(iLayer) = MeshParam%CREDV(iElem)*(1.-HydroParam%ALB)*LimnoParam%fPAR*MeteoParam%SolarRad(iElem)*exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer+1,iElem)))
                    ! light_at_the_bottom,in FABM it's at the bottom of the current layer  (W/m� PAR)
                    LimnoParam%aLPARBot(iLayer) = LimnoParam%uLPAR0(iLayer)*exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer,iElem)))
                Else
                    If ( iLayer == HydroParam%ElSmallm(iElem) ) Then !Bottom layer
                        ! light_at_the_surface,in FABM it's at the top of the current layer (W/m� PAR)
                        LimnoParam%uLPAR0(iLayer) = LimnoParam%aLPARBot(iLayer+1)
                        ! light_at_the_bottom,in FABM it's at the bottom of the current layer  (W/m� PAR)
                        LimnoParam%aLPARBot(iLayer) = LimnoParam%uLPAR0(iLayer)*exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer,iElem)))
                    Elseif ( iLayer == HydroParam%ElCapitalM(iElem) ) Then !Surface layer
                        ! light_at_the_surface,in FABM it's at the top of the current layer (W/m� PAR)
                        LimnoParam%uLPAR0(iLayer) = MeshParam%CREDV(iElem)*(1.-HydroParam%ALB)*LimnoParam%fPAR*MeteoParam%SolarRad(iElem)*exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer+1,iElem)))
                        ! light_at_the_bottom,in FABM it's at the bottom of the current layer  (W/m� PAR)
                        LimnoParam%aLPARBot(iLayer) = LimnoParam%uLPAR0(iLayer)*exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer,iElem)))
                    Else !Intermediary layer
                        ! light_at_the_surface,in FABM it's at the top of the current layer (W/m� PAR)
                        LimnoParam%uLPAR0(iLayer) = LimnoParam%aLPARBot(iLayer+1)
                        ! light_at_the_bottom,in FABM it's at the bottom of the current layer  (W/m� PAR)
                        LimnoParam%aLPARBot(iLayer) = LimnoParam%uLPAR0(iLayer)*exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer,iElem)))
                    EndIf
                Endif
            EndDo            
            
            If (LimnoParam%amac(mm)==LimnoParam%NYMP) Then
                !-----------------------------------------------------------------------
                !  temperature functions of vegetation
                !-----------------------------------------------------------------------
                !  temperature_function_of_vegetation_production
                LimnoParam%uFunTmProdMac = uFunTmVeg(LimnoParam%sDTempW(HydroParam%ElCapitalm(iElem),iElem),LimnoParam%cQ10ProdMac(mm)) 
                !  temperature_function_of_vegetation_respiration
                LimnoParam%uFunTmRespMac = uFunTmVeg(LimnoParam%sDTempW(HydroParam%ElCapitalm(iElem),iElem),LimnoParam%cQ10RespMac(mm)) 
                
                LimnoParam%uMuMaxTmMac = LimnoParam%cMuMaxMac(mm)*LimnoParam%cQ10ProdMac(mm)**(0.1*(LimnoParam%sDTempW(HydroParam%ElCapitalm(iElem),iElem)-20.))
		        LimnoParam%aLPAR1Mac = LimnoParam%uLPAR0(HydroParam%ElCapitalm(iElem))*exp(-LimnoParam%aExtCoef(HydroParam%ElCapitalm(iElem),iElem)*0.)
                LimnoParam%aLPAR2Mac =  LimnoParam%aLPAR1Mac*exp(-LimnoParam%aExtCoef(HydroParam%ElCapitalm(iElem),iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(HydroParam%ElCapitalM(iElem),iElem)))
                LimnoParam%uhLMac = LimnoParam%hLRefMac(mm)*LimnoParam%uFunTmProdMac
                LimnoParam%aFunLSubMac = 1./(-LimnoParam%aExtCoef(HydroParam%ElCapitalm(iElem),iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(HydroParam%ElCapitalM(iElem),iElem)))*log((1.+LimnoParam%aLPAR1Mac/(LimnoParam%uhLMac+NearZero))/(1.+LimnoParam%aLPAR2Mac/(LimnoParam%uhLMac+NearZero)))
                
            Else
                !-----------------------------------------------------------------------
                !  temperature functions of vegetation
                !-----------------------------------------------------------------------
                !  temperature_function_of_vegetation_production
                LimnoParam%uFunTmProdMac = uFunTmVeg(sum(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))/(max(1,size(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem)))),LimnoParam%cQ10ProdMac(mm)) 
                !  temperature_function_of_vegetation_respiration
                LimnoParam%uFunTmRespMac = uFunTmVeg(sum(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))/(max(1,size(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem)))),LimnoParam%cQ10RespMac(mm)) 

                LimnoParam%uMuMaxTmMac=LimnoParam%cMuMaxMac(mm)*LimnoParam%cQ10ProdMac(mm)**(0.1*(sum(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))/(max(1,size(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))))-20.))
		        LimnoParam%aLPAR1Mac = LimnoParam%uLPAR0(HydroParam%ElCapitalm(iElem))*exp(-LimnoParam%aExtCoef(HydroParam%ElCapitalm(iElem),iElem)*0.)
                LimnoParam%aLPAR2Mac =  LimnoParam%aLPAR1Mac*exp(-sum(LimnoParam%aExtCoef(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))/(max(1,size(LimnoParam%aExtCoef(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))))*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(HydroParam%ElSmallm(iElem),iElem)))
                LimnoParam%uhLMac = LimnoParam%hLRefMac(mm)*LimnoParam%uFunTmProdMac
                LimnoParam%aFunLSubMac = 1./(-sum(LimnoParam%aExtCoef(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))/(max(1,size(LimnoParam%aExtCoef(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))))*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(HydroParam%ElSmallm(iElem),iElem)))*log((1.+LimnoParam%aLPAR1Mac/(LimnoParam%uhLMac+NearZero))/(1.+LimnoParam%aLPAR2Mac/(LimnoParam%uhLMac+NearZero)))
            EndIf
            
            !  max._growth_rate_at_current_temp._AND_light
            LimnoParam%aMuTmLMac= LimnoParam%uMuMaxTmMac*LimnoParam%ufDay*LimnoParam%aFunLSubMac*LimnoParam%bfShootVeg
            !-----------------------------------------------------------------------
            !  Nutrient limitation functions
            !-----------------------------------------------------------------------
            !  Droop_function_(C)_for_vegetation
            LimnoParam%aCLimMac = max(0.0, (1.0 - LimnoParam%cCDMacMin(mm) / LimnoParam%rCDMac(iElem,mm)) * LimnoParam%cCDMacMax(mm) / (LimnoParam%cCDMacMax(mm) - LimnoParam%cCDMacMin(mm)) )
            !  Droop_function_(N)_for_vegetation
            LimnoParam%aNLimMac = max(0.0, (1.0 - LimnoParam%cNDMacMin(mm) / LimnoParam%rNDMac(iElem,mm)) * LimnoParam%cNDMacMax(mm) / (LimnoParam%cNDMacMax(mm) - LimnoParam%cNDMacMin(mm)) )
            !  Droop_function_(P)_for_vegetation
            LimnoParam%aPLimMac = max(0.0, (1.0 - LimnoParam%cPDMacMin(mm) / LimnoParam%rPDMac(iElem,mm)) * LimnoParam%cPDMacMax(mm) / (LimnoParam%cPDMacMax(mm) - LimnoParam%cPDMacMin(mm)) )
            !  nutrient_limitation_function_of_vegetation
            LimnoParam%aNutLimMac = min(LimnoParam%aCLimMac,LimnoParam%aPLimMac,LimnoParam%aNLimMac)
            !  actual_growth_rate_of_vegetation
            LimnoParam%aMuMac = LimnoParam%aMuTmLMac * LimnoParam%aNutLimMac  
            
            !-----------------------------------------------------------------------
            !  vegetation growth rate adjust and correction 
            !-----------------------------------------------------------------------
            !  intrinsic_net_increase_rate_of_vegetation
            LimnoParam%akDIncrMac = LimnoParam%aMuTmLMac - LimnoParam%kDRespMac(mm) * LimnoParam%uFunTmRespMac - LimnoParam%kMortMac(mm)
            !  logistic_correction_of_vegetation
            LimnoParam%tDEnvMac = max(0.0, LimnoParam%sDMac(iElem,mm)**2.0*LimnoParam%akDIncrMac / (LimnoParam%cDCarrMac(mm)+NearZero) )
            !  logistic_correction_of_production
            LimnoParam%tDEnvProdMac = LimnoParam%aMuMac / LimnoParam%cMuMaxMac(mm) * LimnoParam%tDEnvMac
            !  vegetation_production
            LimnoParam%tDProdMac(iElem,mm) = max(0.0, LimnoParam%aMuMac * LimnoParam%sDMac(iElem,mm) - LimnoParam%tDEnvProdMac)
           
            !-----------------------------------------------------------------------
            !  vegetation nutrient uptake --- Carbon
            !-----------------------------------------------------------------------
            !  fraction_of_C_uptake_from_sediment
            d1 = sum(LimnoParam%sDicW(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))/(max(1,size(LimnoParam%sDicW(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))))
            If (LimnoParam%bfRootVeg <= NearZero) Then
                LimnoParam%afCUptVegS = 0.0
            Elseif (LimnoParam%fFloatVeg(mm) + LimnoParam%bfSubVeg <= NearZero) Then
                LimnoParam%afCUptVegS = 1.0
            Else
                LimnoParam%afCUptVegS = LimnoParam%fSedUptVegMax(mm) / (1.0 + LimnoParam%fSedUptVegCoef(mm) * ((((LimnoParam%oDicS(iElem)+NearZero) / (d1+NearZero)) )** LimnoParam%fSedUptVegExp(mm)))
            EndIf
            !  maximum_C_uptake_rate_of_vegetation,_corrected_for_C/D_ratio
            LimnoParam%aVCUptMaxCorMac = max( 0.0, LimnoParam%cVCUptMaxMaC(mm) * LimnoParam%aFunLSubMac * LimnoParam%uFunTmProdMac * (LimnoParam%cCDMacMax(mm)-LimnoParam%rCDMac(iElem,mm)) / (LimnoParam%cCDMacMax(mm)-LimnoParam%cCDMacMin(mm)))
            !  C_uptake_RATE_by_subm_AND_floating_parts
            LimnoParam%aVCUptMacW = d1 * LimnoParam%aVCUptMaxCorMac / (LimnoParam%aVCUptMaxCorMac / LimnoParam%cAffCUptMac(mm) + d1)
            LimnoParam%tCUptMacW(iElem,mm) = (1.0 - LimnoParam%afCUptVegS) * LimnoParam%aVCUptMaxCorMac * d1 / (LimnoParam%aVCUptMaxCorMac / LimnoParam%cAffCUptMac(mm) + d1) * LimnoParam%sDMac(iElem,mm) 
            !  C_uptake_rate_by_roots
            LimnoParam%aVCUptMacS = LimnoParam%oDicS(iElem) * LimnoParam%aVCUptMaxCorMac / (LimnoParam%aVCUptMaxCorMac / LimnoParam%cAffCUptMac(mm) + LimnoParam%oDicS(iElem))
            LimnoParam%tCUptMacS(iElem,mm) = LimnoParam%oDicS(iElem) * LimnoParam%aVCUptMaxCorMac / (LimnoParam%aVCUptMaxCorMac / LimnoParam%cAffCUptMac(mm) + LimnoParam%oDicS(iElem))

            !-----------------------------------------------------------------------
            !  vegetation nutrient uptake --- Nitrogen
            !-----------------------------------------------------------------------
            !  fraction_of_N_uptake_from_sediment
            d1 = sum(LimnoParam%oNDissW(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))/(max(1,size(LimnoParam%oNDissW(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))))
            d2 = sum(LimnoParam%sNH4W(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))/(max(1,size(LimnoParam%sNH4W(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))))
            d3 = sum(LimnoParam%sNO3W(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))/(max(1,size(LimnoParam%sNO3W(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))))
            If (LimnoParam%bfRootVeg <= NearZero) Then
                LimnoParam%afNUptVegS = 0.0
            Elseif (LimnoParam%fFloatVeg(mm) + LimnoParam%bfSubVeg <= NearZero) Then
                LimnoParam%afNUptVegS = 1.0
            Else
                LimnoParam%afNUptVegS = LimnoParam%fSedUptVegMax(mm) / (1.0 + LimnoParam%fSedUptVegCoef(mm) * ((((LimnoParam%oNDissS(iElem)+NearZero) / (d1+NearZero)) )** LimnoParam%fSedUptVegExp(mm)))
            EndIf
            !  maximum_N_uptake_rate_of_vegetation,_corrected_for_N/D_ratio
            LimnoParam%aVNUptMaxCorMac = max( 0.0, LimnoParam%cVNUptMaxMaC(mm) *LimnoParam%uFunTmProdMac * (LimnoParam%cNDMacMax(mm)-LimnoParam%rNDMac(iElem,mm)) / (LimnoParam%cNDMacMax(mm)-LimnoParam%cNDMacMin(mm)))
            LimnoParam%ahNUptMac = LimnoParam%aVNUptMaxCorMac / LimnoParam%cAffNUptMac(mm)
            !  N_uptake_RATE_by_subm_AND_floating_parts
            LimnoParam%aVNUptMacW = d1 * LimnoParam%aVNUptMaxCorMac / (LimnoParam%aVNUptMaxCorMac / LimnoParam%cAffNUptMac(mm) + d1)
            LimnoParam%tNUptMacW(iElem,mm) = (1.0 - LimnoParam%afNUptVegS) * LimnoParam%aVNUptMaxCorMac * d1 / (LimnoParam%aVNUptMaxCorMac / LimnoParam%cAffNUptMac(mm) + d1) * LimnoParam%sDMac(iElem,mm) 
            
            LimnoParam%afNH4UptMacW = d2*d3/((LimnoParam%ahNUptMac+d2)*(LimnoParam%ahNUptMac+d3)) + d2*LimnoParam%ahNUptMac/((d2+d3)*(LimnoParam%ahNUptMac+d3))
            LimnoParam%tNUptNH4MacW(iElem,mm) =  LimnoParam%afNH4UptMacW * LimnoParam%tNUptMacW(iElem,mm) 
	        LimnoParam%tNUptNO3MacW(iElem,mm) =  (1.-LimnoParam%afNH4UptMacW) * LimnoParam%tNUptMacW(iElem,mm)
            
            !  N_uptake_rate_by_roots
            LimnoParam%aVNUptMacS = LimnoParam%oNDissS(iElem) * LimnoParam%aVNUptMaxCorMac / (LimnoParam%aVNUptMaxCorMac / LimnoParam%cAffNUptMac(mm) + LimnoParam%oNDissS(iElem))
            LimnoParam%tNUptMacS(iElem,mm) = LimnoParam%oNDissS(iElem) * LimnoParam%aVNUptMaxCorMac / (LimnoParam%aVNUptMaxCorMac / LimnoParam%cAffNUptMac(mm) + LimnoParam%oNDissS(iElem))
		    
            d2 = LimnoParam%oNH4S(iElem)
            d3 = LimnoParam%oNH4S(iElem)
            LimnoParam%afNH4UptMacS = d2*d3/((LimnoParam%ahNUptMac+d2)*(LimnoParam%ahNUptMac+d3)) + d2*LimnoParam%ahNUptMac/((d2+d3)*(LimnoParam%ahNUptMac+d3))
            LimnoParam%tNUptNH4MacS(iElem,mm) =  LimnoParam%afNH4UptMacS * LimnoParam%tNUptMacS(iElem,mm) 
	        LimnoParam%tNUptNO3MacS(iElem,mm) =  (1.-LimnoParam%afNH4UptMacS) * LimnoParam%tNUptMacS(iElem,mm)
            
            !-----------------------------------------------------------------------
            !  vegetation nutrient uptake --- Phosphorus
            !-----------------------------------------------------------------------
            !  fraction_of_P_uptake_from_sediment
            d1 = sum(LimnoParam%sPO4W(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))/(max(1,size(LimnoParam%sPO4W(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))))
            If (LimnoParam%bfRootVeg <= NearZero) Then
                LimnoParam%afPUptVegS = 0.0
            Elseif (LimnoParam%fFloatVeg(mm) + LimnoParam%bfSubVeg <= NearZero) Then
                LimnoParam%afPUptVegS = 1.0
            Else
                LimnoParam%afPUptVegS = LimnoParam%fSedUptVegMax(mm) / (1.0 + LimnoParam%fSedUptVegCoef(mm) * ((((LimnoParam%oPO4S(iElem)+NearZero) / (d1+NearZero)) )** LimnoParam%fSedUptVegExp(mm)))
            EndIf
            !  maximum_P_uptake_rate_of_vegetation,_corrected_for_P/D_ratio
            LimnoParam%aVPUptMaxCorMac = max( 0.0, LimnoParam%cVPUptMaxMaC(mm) *LimnoParam%uFunTmProdMac * (LimnoParam%cPDMacMax(mm)-LimnoParam%rPDMac(iElem,mm)) / (LimnoParam%cPDMacMax(mm)-LimnoParam%cPDMacMin(mm)))
            !  P_uptake_RATE_by_subm_AND_floating_parts
            LimnoParam%aVPUptMacW = d1 * LimnoParam%aVPUptMaxCorMac / (LimnoParam%aVPUptMaxCorMac / LimnoParam%cAffPUptMac(mm) + d1)
            LimnoParam%tPUptMacW(iElem,mm) = (1.0 - LimnoParam%afPUptVegS) * LimnoParam%aVPUptMaxCorMac * d1 / (LimnoParam%aVPUptMaxCorMac / LimnoParam%cAffPUptMac(mm) + d1) * LimnoParam%sDMac(iElem,mm) 
            !    P_uptake_rate_by_roots
            LimnoParam%aVPUptMacS = LimnoParam%oPO4S(iElem) * LimnoParam%aVPUptMaxCorMac / (LimnoParam%aVPUptMaxCorMac / LimnoParam%cAffPUptMac(mm) + LimnoParam%oPO4S(iElem))
            LimnoParam%tPUptMacS(iElem,mm) = LimnoParam%oPO4S(iElem) * LimnoParam%aVPUptMaxCorMac / (LimnoParam%aVPUptMaxCorMac / LimnoParam%cAffPUptMac(mm) + LimnoParam%oPO4S(iElem))

            !=======================================================================
            !  the Dissimilation part!
            !=======================================================================
            !-----------------------------------------------------------------------
            !  vegetation respiration
            !-----------------------------------------------------------------------
            !  dark_respiration_of_vegetation
            LimnoParam%tDRespMac(iElem,mm) = LimnoParam%kDRespMac(mm) * LimnoParam%uFunTmRespMac * LimnoParam%sDMac(iElem,mm)
            LimnoParam%tDRespMacW(iElem,mm) = (1.-LimnoParam%bfRootVeg) * LimnoParam%tDRespMac(iElem,mm)
            LimnoParam%tDRespMacS(iElem,mm) = LimnoParam%bfRootVeg * LimnoParam%tDRespMac(iElem,mm)
            !-----------------------------------------------------------------------
            !  vegetation excretion
            !-----------------------------------------------------------------------
            !  C_excretion_by_vegetation
            LimnoParam%tCExdMac(iElem,mm) = (2.0 * LimnoParam%rCDMac(iElem,mm)) / (LimnoParam%cCDMacMax(mm) + LimnoParam%rCDMac(iElem,mm)) *LimnoParam%rCDMac(iElem,mm) * LimnoParam%tDRespMac(iElem,mm) ! pl613
            !  C_excretion_by_vegetation_to_sediment
            LimnoParam%tCExdMacS(iElem,mm) = LimnoParam%bfRootVeg * LimnoParam%tCExdMac(iElem,mm)
            !  C_excretion_by_vegetation_to_water
            LimnoParam%tCExdMacW(iElem,mm) = LimnoParam%tCExdMac(iElem,mm) - LimnoParam%tCExdMacS(iElem,mm)

            !  N_excretion_by_vegetation
            LimnoParam%tNExdMac(iElem,mm) = (2.0 * LimnoParam%rNDMac(iElem,mm)) / (LimnoParam%cNDMacMax(mm) + LimnoParam%rNDMac(iElem,mm)) *LimnoParam%rNDMac(iElem,mm) * LimnoParam%tDRespMac(iElem,mm) ! pl613
            !  N_excretion_by_vegetation_to_sediment
            LimnoParam%tNExdMacS(iElem,mm) = LimnoParam%bfRootVeg * LimnoParam%tNExdMac(iElem,mm)
            !  N_excretion_by_vegetation_to_water
            LimnoParam%tNExdMacW(iElem,mm) = LimnoParam%tNExdMac(iElem,mm) - LimnoParam%tNExdMacS(iElem,mm)
            
            !  P_excretion_by_vegetation
            LimnoParam%tPExdMac(iElem,mm) = (2.0 * LimnoParam%rPDMac(iElem,mm)) / (LimnoParam%cPDMacMax(mm) + LimnoParam%rPDMac(iElem,mm)) *LimnoParam%rPDMac(iElem,mm) * LimnoParam%tDRespMac(iElem,mm) ! pl613
            !  P_excretion_by_vegetation_to_sediment
            LimnoParam%tPExdMacS(iElem,mm) = LimnoParam%bfRootVeg * LimnoParam%tPExdMac(iElem,mm)
            !  P_excretion_by_vegetation_to_water
            LimnoParam%tPExdMacW(iElem,mm) = LimnoParam%tPExdMac(iElem,mm) - LimnoParam%tPExdMacS(iElem,mm)

            !-----------------------------------------------------------------------
            !  vegetation mortality,Dry-weight
            !-----------------------------------------------------------------------
            !  logistic_correction_of_mortality
            LimnoParam%tDEnvMortMac = LimnoParam%tDEnvMac - LimnoParam%tDEnvProdMac
            !  total_mortality_flux_DW_vegetation
            LimnoParam%tDMortMac(iElem,mm) = LimnoParam%kMortMac(mm) * LimnoParam%sDMac(iElem,mm) + LimnoParam%tDEnvMortMac
            !-----------------------------------------------------------------------
            !  vegetation Mortality,Carbon
            !-----------------------------------------------------------------------
            !  C_mortality_flux_of_vegetation
            LimnoParam%tCMortMac(iElem,mm) = LimnoParam%rCDMac(iElem,mm) * LimnoParam%tDMortMac(iElem,mm)
            !  mortality_flux_of_vegetation_becoming_dissolved_C_in_sediment
            LimnoParam%tCMortMacS(iElem,mm) = LimnoParam%bfRootVeg * LimnoParam%tCMortMac(iElem,mm)
            !  mortality_flux_of_vegetation_becoming_dissolved_N_in_water
            LimnoParam%tCMortMacW(iElem,mm) = LimnoParam%tCMortMac(iElem,mm) - LimnoParam%tCMortMacS(iElem,mm)
            
            !-----------------------------------------------------------------------
            !  vegetation Mortality,Nitrogen
            !-----------------------------------------------------------------------
            !  N_mortality_flux_of_vegetation
            LimnoParam%tNMortMac(iElem,mm) = LimnoParam%rNDMac(iElem,mm) * LimnoParam%tDMortMac(iElem,mm)
            !  mortality_flux_of_vegetation_becoming_dissolved_N_in_sediment
            LimnoParam%tNMortMacS(iElem,mm) = LimnoParam%bfRootVeg * LimnoParam%tNMortMac(iElem,mm)
            !LimnoParam%tNMortMacNH4S(iElem,mm) = LimnoParam%fDissMortMac  * LimnoParam%tNMortMacS(iElem,mm)
            !  mortality_flux_of_vegetation_becoming_dissolved_N_in_water
            LimnoParam%tNMortMacW(iElem,mm) = LimnoParam%tNMortMac(iElem,mm) - LimnoParam%tNMortMacS(iElem,mm)
            !LimnoParam%tNMortMacNH4W(iElem,mm) = LimnoParam%fDissMortMac  * LimnoParam%tNMortMacW(iElem,mm)
            
            !-----------------------------------------------------------------------
            !  vegetation Mortality,Phosphrus
            !-----------------------------------------------------------------------
            !  P_mortality_flux_of_vegetation
            LimnoParam%tPMortMac(iElem,mm) = LimnoParam%rPDMac(iElem,mm) * LimnoParam%tDMortMac(iElem,mm)
            !  mortality_flux_of_vegetation_becoming_dissolved_C_in_sediment
            LimnoParam%tPMortMacS(iElem,mm) = LimnoParam%bfRootVeg * LimnoParam%tPMortMac(iElem,mm)
            !  mortality_flux_of_vegetation_becoming_dissolved_N_in_water
            LimnoParam%tPMortMacW(iElem,mm) = LimnoParam%tPMortMac(iElem,mm) - LimnoParam%tPMortMacS(iElem,mm)
            
            !-----------------------------------------------------------------------
            !  Grazing by Birds,Dry-weight
            !-----------------------------------------------------------------------
		    LimnoParam%tDGrazMacBird(iElem,mm)=LimnoParam%cPrefMacBird(mm)*LimnoParam%sDMac(iElem,mm)/(LimnoParam%hDMacBird(mm)+LimnoParam%sDMac(iElem,mm))*LimnoParam%cBirdsPerha/10000.*LimnoParam%cDGrazPerBird(mm)
		    LimnoParam%tCGrazMacBird(iElem,mm)=LimnoParam%tDGrazMacBird(iElem,mm)*LimnoParam%rCDMac(iElem,mm)
            LimnoParam%tNGrazMacBird(iElem,mm)=LimnoParam%tDGrazMacBird(iElem,mm)*LimnoParam%rNDMac(iElem,mm)
            LimnoParam%tPGrazMacBird(iElem,mm)=LimnoParam%tDGrazMacBird(iElem,mm)*LimnoParam%rPDMac(iElem,mm)
            !-----------------------------------------------------------------------
            !  Egestion by Birds,Dry-weight
            !-----------------------------------------------------------------------
		    LimnoParam%tDEgesBird(iElem,mm) = (1.-LimnoParam%fDAssBird(mm))*LimnoParam%tDGrazMacBird(iElem,mm)
		    LimnoParam%tCEgesBird(iElem,mm) = (1.-LimnoParam%fDAssBird(mm))*LimnoParam%tCGrazMacBird(iElem,mm)*LimnoParam%rCDMac(iElem,mm)
		    LimnoParam%tNEgesBird(iElem,mm) = (1.-LimnoParam%fDAssBird(mm))*LimnoParam%tNGrazMacBird(iElem,mm)*LimnoParam%rNDMac(iElem,mm)
		    LimnoParam%tPEgesBird(iElem,mm) = (1.-LimnoParam%fDAssBird(mm))*LimnoParam%tPGrazMacBird(iElem,mm)*LimnoParam%rPDMac(iElem,mm)
            
            !-----------------------------------------------------------------------
            !  Update O2 in water
            !-----------------------------------------------------------------------
            !  O2_production_to_water_due_to_NO3_uptake_by_macrophytes
            LimnoParam%tO2UptNO3MacW(iElem,mm) = LimnoParam%O2PerNO3 * LimnoParam%molO2molN * LimnoParam%bfSubVeg * LimnoParam%tNUptNO3MacW(iElem,mm)
            !  correction_of_O2_demand_in_water_at_low_oxygen_conc.
            d1 = sum(LimnoParam%sO2W(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))/(max(1,size(LimnoParam%sO2W(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalm(iElem),iElem))))
            LimnoParam%aCorO2BOD = d1 / (LimnoParam%hO2Bac + d1)
            !  submerged_O2_respiration
            LimnoParam%tO2RespMacW(iElem,mm) = LimnoParam%molO2molC * LimnoParam%cCPerDW * LimnoParam%bfSubVeg * LimnoParam%tDRespMac(iElem,mm) * LimnoParam%aCorO2BOD
            !  root_O2_respiration
            LimnoParam%tO2RespMacS(iElem,mm) = LimnoParam%molO2molC * LimnoParam%cCPerDW * LimnoParam%bfRootVeg * LimnoParam%tDRespMac(iElem,mm) * LimnoParam%afOxySed
            !  vegetation_O2_production
            LimnoParam%tO2ProdMac(iElem,mm) = LimnoParam%molO2molC * LimnoParam%cCPerDW * LimnoParam%tDProdMac(iElem,mm)
            !  O2_transport_to_roots
            LimnoParam%tO2ProdMacS(iElem,mm) = min(LimnoParam%tO2RespMacS(iElem,mm),LimnoParam%tO2ProdMac(iElem,mm))
            !  O2_used_for_vegetation_production
            LimnoParam%tO2ProdMacW(iElem,mm) = min( LimnoParam%tO2ProdMac(iElem,mm) - LimnoParam%tO2ProdMacS(iElem,mm), LimnoParam%bfSubVeg * LimnoParam%tO2ProdMac(iElem,mm))            
            
            !-------------------------------------------------------------------------------------!
            !                               Source/Sink Term for Macrophytes                      !
            !-------------------------------------------------------------------------------------!
            
			LimnoParam%dVarEstS(iElem,1) =   LimnoParam%tDProdMac(iElem,mm)                             &!Production
                                            -LimnoParam%tDRespMac(iElem,mm)                             &!Respiration
                                            -LimnoParam%tDMortMac(iElem,mm)                             &!Mortality
                                            -LimnoParam%tDGrazMacBird(iElem,mm)                         !Grazing by Birds
			
            LimnoParam%dVarEstS(iElem,2) =   LimnoParam%tCUptMacW(iElem,mm)                             &!Uptake in water
                                            -LimnoParam%tCUptMacS(iElem,mm)                             &!Uptake in sediment
                                            -LimnoParam%tCExdMac(iElem,mm)                              &!Excretion
                                            -LimnoParam%tCMortMac(iElem,mm)                             &!Mortality
                                            -LimnoParam%tCGrazMacBird(iElem,mm)                         !Grazing by Birds
            
            LimnoParam%dVarEstS(iElem,3) =   LimnoParam%tNUptMacW(iElem,mm)                             &!Uptake in water
                                            -LimnoParam%tNUptMacS(iElem,mm)                             &!Uptake in sediment
                                            -LimnoParam%tNExdMac(iElem,mm)                              &!Excretion
                                            -LimnoParam%tNMortMac(iElem,mm)                             &!Mortality
                                            -LimnoParam%tNGrazMacBird(iElem,mm)                         !Grazing by Birds
           
            LimnoParam%dVarEstS(iElem,4) =   LimnoParam%tPUptMacW(iElem,mm)                             &!Uptake in water
                                            -LimnoParam%tPUptMacS(iElem,mm)                             &!Uptake in sediment
                                            -LimnoParam%tPExdMac(iElem,mm)                              &!Excretion
                                            -LimnoParam%tPMortMac(iElem,mm)                             &!Mortality
                                            -LimnoParam%tPGrazMacBird(iElem,mm)                         !Grazing by Birds           
            
            LimnoParam%sDMac(iElem,mm) = max(NearZero,(dtday*LimnoParam%dVarEstS(iElem,1)+LimnoParam%sDMacP(iElem,mm)))
            LimnoParam%sCMac(iElem,mm) = max(NearZero,(dtday*LimnoParam%dVarEstS(iElem,2)+LimnoParam%sDMacP(iElem,mm)))
            LimnoParam%sNMac(iElem,mm) = max(NearZero,(dtday*LimnoParam%dVarEstS(iElem,3)+LimnoParam%sDMacP(iElem,mm)))
            LimnoParam%sPMac(iElem,mm) = max(NearZero,(dtday*LimnoParam%dVarEstS(iElem,4)+LimnoParam%sDMacP(iElem,mm)))
                
        EndDo !Loop Cell
    EndDo !Loop Functional Group

    
Return
End
