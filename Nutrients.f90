Subroutine Nutrients(HydroParam,MeshParam,LimnoParam,dt,dtday)   !PO4,NH4,NO3 & SiO2

    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars

    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    Integer:: index,iElem,iLayer,iter 
    Real:: dt,dtday
    Real:: V
    Real:: tPResusPO4_bottom,tPResusPAIM_bottom,tNResusNH4_bottom,tNResusNO3_bottom,tPSetAIM_Layer
    Real:: tPDifPO4_bottom,tNDifNH4_bottom,tNDifNO3_bottom,tCDifDic_bottom,tDDifO2_bottom

    !-------------------------------------------------------------------------------------!
    !                           Boundary Condition for PO4                                !
    !-------------------------------------------------------------------------------------!
    index = 6
    HydroParam%uLoadVarEst = LimnoParam%uPLoadPO4
    
    Do iElem = 1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
    
        !-PO4-----------------------------------------------------------------------------------!
        ! 1. Phosphorus sorption                                                                 !
        ! 2. Resuspension (from SedimentFluxes.f90)                                             !
        ! 3. BAC Uptake/Release (from Bacterioplakton.f90)                                      !
        ! 3. Phyto Uptake (from Phytoplankton.f90)                                              !
        ! 4. Macro Uptake (from Macrophytes.f90)                                                !
        ! 5. Mineralization (done in DOM IF InclBac=0)                                          !
        ! 6. Sediment flux (from Sediment.f90)                                                  !
        ! 7. Infiltration                                                                       !
        !---------------------------------------------------------------------------------------!
        !-------------------------------------------------------------------------------------!
        !                              PO4 Processes in Water					              !
        !-------------------------------------------------------------------------------------!
        !1) Phosphorus sorption 
            !  correction_of_O2_demand_in_water_at_low_oxygen_conc.
            LimnoParam%aCorO2BOD = LimnoParam%sO2W(iLayer,iElem)/(LimnoParam%hO2Bac+LimnoParam%sO2W(iLayer,iElem))
	        !  max._P_adsorption_per_g_inorg._matter_in_water (pP/gD))
	        LimnoParam%aPAdsMaxW=LimnoParam%cRelPAdsD + LimnoParam%aCorO2BOD*LimnoParam%cRelPAdsFe*LimnoParam%fFeDIM + LimnoParam%cRelPAdsAl*LimnoParam%fAlDIM 
	        !  P_adsorption_affinity,_corrected_for_redox_conditions (m³/gP)
	        LimnoParam%aKPAdsW=(1.-LimnoParam%fRedMax*(1.-LimnoParam%aCorO2BOD))*LimnoParam%cKPAdsOx						     
	        !!  P_adsorption_isotherm_onto_inorg._matter_in_sediment (gP/gD)
	        LimnoParam%aPIsoAdsW=LimnoParam%aPAdsMaxW*LimnoParam%aKPAdsW*LimnoParam%sPO4W(iLayer,iElem)/(1.+LimnoParam%aKPAdsW*LimnoParam%sPO4W(iLayer,iElem))                             
	        !  equilibrium_conc. (gP/m³)
	        LimnoParam%aPEqIMW=LimnoParam%aPIsoAdsW*LimnoParam%sDIMW(iLayer,iElem)											                         
	        !  sorption_flux_in_water (gP/m³/d)
            If (LimnoParam%InclMatInorg == 1) Then
	            LimnoParam%wPSorpIMW(iLayer,iElem)=LimnoParam%kPSorp*(LimnoParam%aPEqIMW-LimnoParam%sPAIMW(iLayer,iElem))							                         
            Else
                LimnoParam%wPSorpIMW(iLayer,iElem)=0.
            EndIf
            
            
            
        !2) PO4 Resuspension and difusion
            If (iLayer == HydroParam%ElSmallm(iElem)) Then
                tPDifPO4_bottom = LimnoParam%tPDifPO4(iElem)
                tPResusPO4_bottom = LimnoParam%tPResusPO4(iElem)
            Else
                tPDifPO4_bottom = 0.
                tPResusPO4_bottom = 0.
            EndIf
            

            !-------------------------------------------------------------------------------------!
            !                               Source/Sink Term for PO4 			        	      !
            !-------------------------------------------------------------------------------------!
	        HydroParam%dVarEst(iLayer,iElem,1) =  LimnoParam%cPBackLoad                                                         & ! BackLoad                      
		                            - LimnoParam%wPSorpIMW(iLayer,iElem)                                                        & ! Adsorption/Desorption PO4 <-> IM
                                    + tPResusPO4_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))                       & ! Resuspension
                                    + tPDifPO4_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))                         & ! Difusion
		                            + LimnoParam%rPDDomW(iLayer,iElem,1) * LimnoParam%wDMinAerW(iLayer,iElem,1)                 & ! Mineralization DOM Labile -> PO4 If InclBac=0
                                    + LimnoParam%rPDDomW(iLayer,iElem,1) * LimnoParam%wDMinAnaerW(iLayer,iElem,1)               & ! Mineralization DOM Labile -> PO4 If InclBac=0
                                    - SUM(LimnoParam%wPUptPhyt(iLayer,iElem,:))                                                 & ! Uptake Phyto PO4 -> PHYTO
                                    + SUM(LimnoParam%wPExcrPhytW(iLayer,iElem,:))                                               & ! Phytoplankton Excretion
		                            - SUM(LimnoParam%wPUptBac(iLayer,iElem,:))                                                  & ! Uptake/Release BAC <-> PO4 If InclBac=1
        		                    + SUM(LimnoParam%wPExcrZoo(iLayer,iElem,:))                                                 & ! zooplankton P excretion
                                    + SUM( LimnoParam%wPExcrFiAd(iLayer,iElem,:))                                               & ! adult fish excretion
                                    + SUM( LimnoParam%wPExcrFiJv(iLayer,iElem,:))                                               & ! Juveline fish excretion
                                    - SUM(LimnoParam%tPUptMacW(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))                                            ! Uptake Macro PO4 -> MACRO
   
        EndDo
    EndDo
            
    !-------------------------------------------------------------------------------------!
    !                       Solver Transport Equation for PO4              	    	      !
    !-------------------------------------------------------------------------------------!
    If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
        Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sPO4W,LimnoParam%sPO4WP,dt,dtday,HydroParam,MeshParam) 
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
        Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sPO4W,LimnoParam%sPO4WP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
    ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
        Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sPO4W,LimnoParam%sPO4WP,dt,dtday,HydroParam,MeshParam)
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
        Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sPO4W,LimnoParam%sPO4WP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
    EndiF    
    
    index = 7     
    HydroParam%uLoadVarEst = LimnoParam%uPLoadPAIM
    Do iElem = 1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
    
            !-PAIM--------------------------------------------------------------------------------!
            ! 1. Adsorption/Desorption                                                            !
            ! 2. Sedimentation (from SedimentFluxes.f90)                                          !
            ! 3. Resuspension (from SedimentFluxes.f90)                                          !
            !-------------------------------------------------------------------------------------!
            !-------------------------------------------------------------------------------------!
            !                               PAIM Processes in Water					      !
            !-------------------------------------------------------------------------------------!
		    
            !2) PAIM Sedimentation
            If ( iLayer == HydroParam%ElCapitalM(iElem) ) Then        ! Bottom layer or intermediary layer has contribuition of the up layer
                tPSetAIM_Layer = LimnoParam%tPSetAIM(iLayer,iElem)
            Else
                tPSetAIM_Layer = LimnoParam%tPSetAIM(iLayer,iElem) - LimnoParam%tPSetAIM(iLayer+1,iElem)
            Endif
            
            
            !3) PAIM Resuspension 
            If (iLayer == HydroParam%ElSmallm(iElem)) Then
                tPResusPAIM_bottom = LimnoParam%tPResusAIM(iElem)
            Else
                tPResusPAIM_bottom = 0.
            EndIf
            
    
            !-------------------------------------------------------------------------------------!
            !                               Source/Sink Term for PAIM				      !
            !-------------------------------------------------------------------------------------!
	        HydroParam%dVarEst(iLayer,iElem,1) =  LimnoParam%wPSorpIMW(iLayer,iElem)                                            & ! Adsorption/Desorption PO4 <-> IM
                                    - tPSetAIM_Layer/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))                          & ! Sedimentation PAIM
                                    + tPResusPAIM_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))                      ! Resuspension PAIM
            
            
        EndDo

    EndDo    
    !-------------------------------------------------------------------------------------!
    !                       Solver Transport Equation for PAIM              	    	      !
    !-------------------------------------------------------------------------------------!
    
    If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
        Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sPAIMW,LimnoParam%sPAIMWP,dt,dtday,HydroParam,MeshParam) 
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
        Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sPAIMW,LimnoParam%sPAIMWP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
    ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
        Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sPAIMW,LimnoParam%sPAIMWP,dt,dtday,HydroParam,MeshParam)
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
        Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sPAIMW,LimnoParam%sPAIMWP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
    EndiF   
    
    
    
    !-------------------------------------------------------------------------------------!
    !                           Boundary Condition for NH4     	    			      !
    !-------------------------------------------------------------------------------------!
    index = 8
    HydroParam%uLoadVarEst = LimnoParam%uNLoadNH4

    Do iElem = 1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            !-NH4---------------------------------------------------------------------------------!
            ! 1. Nitrification/Denitrification                                                    !
            ! 2. Resuspension                !
            ! 2. Contribution from NH4 Mineralization (from OrgMatter If InclBac=0)               !
            ! 3. Bac Contribution (from OrgMatter If InclBac=1)                                                                !
            ! 4. Phyto Contribution                                                               !
            ! 5. Macro uptake                                                                     !
            ! 6. Sediment flux (from SedimentFluxes.for)                                                       !
            ! 7. Infiltracao                                                                      !
            !-------------------------------------------------------------------------------------!
            !-------------------------------------------------------------------------------------!
            !                               NH4 Processes in Water					      !
            !-------------------------------------------------------------------------------------!
            !1) Nitrificação
            LimnoParam%aCorO2NitrW=(LimnoParam%sO2W(iLayer,iElem)**2.)/(LimnoParam%hO2Nitr**2.+LimnoParam%sO2W(iLayer,iElem)**2.)
            If (MeshParam%BANHADO(iElem)==1) Then
            	LimnoParam%wNNitrW(iLayer,iElem)=2.*LimnoParam%kNitrW*(LimnoParam%cThetaNitr**(LimnoParam%sDTempW(iLayer,iElem)-20.))*LimnoParam%aCorO2NitrW*LimnoParam%sNH4W(iLayer,iElem) 
            Else
            	LimnoParam%wNNitrW(iLayer,iElem)=   LimnoParam%kNitrW*(LimnoParam%cThetaNitr**(LimnoParam%sDTempW(iLayer,iElem)-20.))*LimnoParam%aCorO2NitrW*LimnoParam%sNH4W(iLayer,iElem) 
            EndIf
            
            !2) NH4 Resuspension and Difusion
            If (iLayer == HydroParam%ElSmallm(iElem)) Then
                tNDifNH4_bottom = LimnoParam%tNDifNH4(iElem)
                tNResusNH4_bottom = LimnoParam%tNResusNH4(iElem)
            Else
                tNDifNH4_bottom = 0.
                tNResusNH4_bottom = 0.
            EndIf
            			
            !-------------------------------------------------------------------------------------!
            !                               Source/Sink Term for NH4				      !
            !-------------------------------------------------------------------------------------!
            HydroParam%dVarEst(iLayer,iElem,1) =   LimnoParam%cNBackLoad                                                        & ! BackLoad
                                        - LimnoParam%wNNitrW(iLayer,iElem)                                                      & ! Nitrification NH4->NO3
                                        + tNResusNH4_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))                   & ! Resuspension
                                        + tNDifNH4_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))                     & ! Difusion
		                                + LimnoParam%rNDDomW(iLayer,iElem,1) * LimnoParam%wDMinAerW(iLayer,iElem,1)             & ! Mineralization DOM Labile -> PO4 If InclBac=0
                                        + LimnoParam%rNDDomW(iLayer,iElem,1) * LimnoParam%wDMinAnaerW(iLayer,iElem,1)           & ! Mineralization DOM Labile -> PO4 If InclBac=0
                                        - SUM(LimnoParam%wNUptNH4Phyt(iLayer,iElem,:))                                          & ! Uptake Phyto NH4 -> PHYTO
                                        + SUM(LimnoParam%wNExcrPhytW(iLayer,iElem,:))                                           & ! Phytoplankton Excretion
            		                    - SUM(LimnoParam%wNUptBac(iLayer,iElem,:))                                              & ! Uptake/Release BAC <-> NH4
                                        + SUM(LimnoParam%wNExcrZoo(iLayer,iElem,:))                                             & ! zooplankton N excretion
                                        + SUM( LimnoParam%wNExcrFiAd(iLayer,iElem,:))                                           & ! adult fish excretion
                                        + SUM( LimnoParam%wNExcrFiJv(iLayer,iElem,:))                                           & ! Juveline fish excretion
                                        - SUM(LimnoParam%tNUptNH4MacW(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))                                            ! Uptake Macro NH4 -> MACRO 
            
        EndDo

    EndDo    
    !-------------------------------------------------------------------------------------!
    !                       Solver Transport Equation for NH4              	    	      !
    !-------------------------------------------------------------------------------------!
    
    If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
        Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sNH4W,LimnoParam%sNH4WP,dt,dtday,HydroParam,MeshParam) 
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
        Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sNH4W,LimnoParam%sNH4WP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
    ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
        Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sNH4W,LimnoParam%sNH4WP,dt,dtday,HydroParam,MeshParam)
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
        Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sNH4W,LimnoParam%sNH4WP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
    EndiF   
    
    !-------------------------------------------------------------------------------------!
    !                           Boundary Condition for NO3     	    			      !
    !-------------------------------------------------------------------------------------!
    index = 9
    HydroParam%uLoadVarEst = LimnoParam%uNLoadNO3
    
    Do iElem = 1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            !-NO3---------------------------------------------------------------------------------!
            ! 1. Nitrification/Denitrification                                                    !
            ! 2. Phyto Uptake                                                                     !
            ! 3. Macrofitas Uptake                                                                !
            ! 4. Sediment flux (done later)                                                       !
            ! 7. Infiltracao                                                                      !
            !-------------------------------------------------------------------------------------!
            !-------------------------------------------------------------------------------------!
            !                               NO3 Processes in Water					      !
            !-------------------------------------------------------------------------------------!
            !1) Denitrification
            !  correction_of_O2_demand_in_water_at_low_oxygen_conc.
            LimnoParam%aCorO2BOD = LimnoParam%sO2W(iLayer,iElem)/(LimnoParam%hO2Bac+LimnoParam%sO2W(iLayer,iElem))
            !  mineralisation_flux_by_denitrification
            LimnoParam%wDDenitW = LimnoParam%sNO3W(iLayer,iElem)*LimnoParam%sNO3W(iLayer,iElem)/(LimnoParam%hNO3Denit*LimnoParam%hNO3Denit+LimnoParam%sNO3W(iLayer,iElem)*LimnoParam%sNO3W(iLayer,iElem))* &
               & (1.0-LimnoParam%aCorO2BOD)*(SUM(LimnoParam%wDMinAnaerW(iLayer,iElem,:))+SUM(LimnoParam%wDMinAerW(iLayer,iElem,:)))
            !  Denitrification_flux
            LimnoParam%wNDenitW(iLayer,iElem) = LimnoParam%NO3PerC*LimnoParam%molNmolC*LimnoParam%cCPerDW*LimnoParam%wDDenitW
            
            !2) NO3 Resuspension and Difusion
            If (iLayer == HydroParam%ElSmallm(iElem)) Then
                tNDifNO3_bottom = LimnoParam%tNDifNO3(iElem)
                tNResusNO3_bottom = LimnoParam%tNResusNO3(iElem)
            Else
                tNResusNO3_bottom = 0.
            EndIf
            			
            !-------------------------------------------------------------------------------------!
            !                               Source/Sink Term for NO3				      !
            !-------------------------------------------------------------------------------------!
            HydroParam%dVarEst(iLayer,iElem,1) = LimnoParam%wNNitrW(iLayer,iElem)                               & !Nitrification NH4->NO3
                                    - LimnoParam%wNDenitW(iLayer,iElem)                                         & !Denitrification
                                    + tNResusNO3_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))       & !Resuspension
                                    + tNDifNO3_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))         & ! Difusion
                                    - SUM(LimnoParam%wNUptNO3Phyt(iLayer,iElem,:))                              & !Uptake Phyto 
				                    - SUM(LimnoParam%tNUptNO3MacW(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))                                            ! Uptake Macro NO3 -> MACRO 
            
        EndDo
    EndDo    
    !-------------------------------------------------------------------------------------!
    !                       Solver Transport Equation for NO3              	    	      !
    !-------------------------------------------------------------------------------------! 
    
    If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
        Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sNO3W,LimnoParam%sNO3WP,dt,dtday,HydroParam,MeshParam) 
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
        Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sNO3W,LimnoParam%sNO3WP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
    ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
        Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sNO3W,LimnoParam%sNO3WP,dt,dtday,HydroParam,MeshParam)
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
        Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sNO3W,LimnoParam%sNO3WP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
    EndiF   

    !-------------------------------------------------------------------------------------!
    !                           Boundary Condition for SiO2     	    			      !
    !-------------------------------------------------------------------------------------!
    index = 10
    HydroParam%uLoadVarEst = LimnoParam%uSiLoadSiO2
    
    Do iElem = 1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            !-SiO2--------------------------------------------------------------------------------!
            ! 1. Contribution from Si Mineralization (done in DOM IF InclBac=0)                   !
            ! 2. Diatoms Uptake                                                                   !
            ! 3. Bac Contribution                                                                 !
            ! 4. Sediment flux (done later)                                                       !
            !-------------------------------------------------------------------------------------!
            !-------------------------------------------------------------------------------------!
            !                               Source/Sink Term for SiO2				      !
            !-------------------------------------------------------------------------------------!
            HydroParam%dVarEst(iLayer,iElem,1) =   LimnoParam%cSiDDomRef(1) * LimnoParam%wDMinAerW(iLayer,iElem,1)               !& ! Mineralization DOM Labile -> Si
                                    !+ SUM(wSiAssBac(KAUX,:))*fSiBac  & ! Mineralizacao DOM -> Si (IF InclBac=1)
		                            !- SUM(wSiUptPhyt(KAUX,:))        & ! Uptake Phyto SiO2 -> Diatoms
		                            !+ SUM(wSiLisPhyt(KAUX,:))          ! Lise Phyto Diatoms -> SiO2

        EndDo

    EndDo    
    !-------------------------------------------------------------------------------------!
    !                       Solver Transport Equation for SiO2              	    	      !
    !-------------------------------------------------------------------------------------!
    If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
        Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sSiO2W,LimnoParam%sSiO2WP,dt,dtday,HydroParam,MeshParam) 
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
        Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sSiO2W,LimnoParam%sSiO2WP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
    ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
        Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sSiO2W,LimnoParam%sSiO2WP,dt,dtday,HydroParam,MeshParam)
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
        Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sSiO2W,LimnoParam%sSiO2WP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
    EndiF   

RETURN
END

