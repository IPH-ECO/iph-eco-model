Subroutine PhytoSediment(HydroParam,MeshParam,LimnoParam,dt,dtday)
    
    Use Hydrodynamic
    Use MeshVars
    Use LimnologyVars
    Use uTempFunction

    Implicit None
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam
    type(LimnologyParam) :: LimnoParam
    Integer:: iElem,gg
    Real:: dt,dtday
    Real:: NearZero = 1e-10
    
    
    
    !-------------------------------------------------------------------------------------!
    !                              Phytoplankton in the sediment					      !
    !-------------------------------------------------------------------------------------!
    !-Phytoplankton in the sediment-------------------------------------------------------!
    ! 1. Erosion                                                                          !
    ! 2. Sedimentation                                                                    !
    ! 3. Ressuspension                                                                    !
    !-------------------------------------------------------------------------------------!
    Do gg = 1, LimnoParam%nphy
        Do iElem = 1,MeshParam%nElem
            !-----------------------------------------------------------------------
            !  Temperature functions for pelagic phytoplankton
            !-----------------------------------------------------------------------
            !  temperature_function_of_phytoplankton
            LimnoParam%uFunTmPhyt = uFunTmBio(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem),iElem),LimnoParam%cSigTmPhyt(gg),LimnoParam%cTmOptPhyt(gg))

            ! Algae respiration_DW and excretion of nutrients
            !  temp._corrected_respiration_constant_of_Algae
            LimnoParam%ukDRespbTmPhyt = LimnoParam%kDRespbPhyt(gg) * LimnoParam%uFunTmPhyt
            !  Total respiration_of_Algae_in_sediment
            LimnoParam%tDRespPhytS(iElem,gg) = LimnoParam%ukDRespbTmPhyt * LimnoParam%sDPhytS(iElem,gg) ! Used in O2 exchange
            !  C_excretion_Algae_in_sediment (Exudation)
            LimnoParam%tCExcrPhytS(iElem,gg) = (2.0 * LimnoParam%rCDPhytS(iElem,gg)) /(LimnoParam%cCDPhytMax(gg) + LimnoParam%rCDPhytS(iElem,gg)) * LimnoParam%rCDPhytS(iElem,gg) * LimnoParam%tDRespPhytS(iElem,gg)
            !  P_Exudation_Algae_in_sediment
            LimnoParam%tPExcrPhytS(iElem,gg) = (2.0 * LimnoParam%rPDPhytS(iElem,gg)) /(LimnoParam%cPDPhytMax(gg) + LimnoParam%rPDPhytS(iElem,gg)) * LimnoParam%rPDPhytS(iElem,gg) * LimnoParam%tDRespPhytS(iElem,gg)
            !  N_Exudation_Algae_in_sediment
            LimnoParam%tNExcrPhytS(iElem,gg) = (2.0 * LimnoParam%rNDPhytS(iElem,gg)) /(LimnoParam%cNDPhytMax(gg) + LimnoParam%rNDPhytS(iElem,gg)) * LimnoParam%rNDPhytS(iElem,gg) * LimnoParam%tDRespPhytS(iElem,gg) ! pl613
            
            ! Lise_Algae_in_water
            LimnoParam%tDLisPhytS(iElem,gg) = LimnoParam%cLisPhytS(gg) * LimnoParam%sDPhytS(iElem,gg)
   			LimnoParam%tCLisPhytS(iElem,gg) = LimnoParam%tDLisPhytS(iElem,gg) * LimnoParam%rCDPhytS(iElem,gg)
			LimnoParam%tNLisPhytS(iElem,gg) = LimnoParam%tDLisPhytS(iElem,gg) * LimnoParam%rNDPhytS(iElem,gg)
            LimnoParam%tPLisPhytS(iElem,gg) = LimnoParam%tDLisPhytS(iElem,gg) * LimnoParam%rPDPhytS(iElem,gg)
			If (LimnoParam%aphy(gg)==LimnoParam%MDIAT.or.LimnoParam%aphy(gg)==LimnoParam%FDIAT) Then
				LimnoParam%tSiLisPhytS(iElem,gg) = LimnoParam%tDLisPhytS(iElem,gg) * LimnoParam%rSiDPhytS(iElem,gg)					
			Else
				LimnoParam%tSiLisPhytS(iElem,gg) = 0.					
            EndIf
            
	        ! POM fraction of Lise
	        LimnoParam%tDLisPomPhytS(iElem,gg) = (1.-LimnoParam%fDLisDomPhyt) * LimnoParam%tDLisPhytS(iElem,gg)
	        LimnoParam%tCLisPomPhytS(iElem,gg) = (1.-LimnoParam%fDLisDomPhyt) * LimnoParam%tCLisPhytS(iElem,gg)
	        LimnoParam%tNLisPomPhytS(iElem,gg) = (1.-LimnoParam%fDLisDomPhyt) * LimnoParam%tNLisPhytS(iElem,gg)			
	        LimnoParam%tPLisPomPhytS(iElem,gg) = (1.-LimnoParam%fDLisDomPhyt)* LimnoParam%tPLisPhytS(iElem,gg)
	        ! DOM fraction of Lise
	        LimnoParam%tDLisDomPhytS(iElem,gg) = LimnoParam%fDLisDomPhyt * LimnoParam%tDLisPhytS(iElem,gg)
	        LimnoParam%tCLisDomPhytS(iElem,gg) = LimnoParam%fDLisDomPhyt * LimnoParam%tCLisPhytS(iElem,gg)
	        LimnoParam%tNLisDomPhytS(iElem,gg) = LimnoParam%fDLisDomPhyt * LimnoParam%tNLisPhytS(iElem,gg)		
	        LimnoParam%tPLisDomPhytS(iElem,gg) = LimnoParam%fDLisDomPhyt * LimnoParam%tPLisPhytS(iElem,gg)	
        	
	        If (LimnoParam%aphy(gg)==LimnoParam%MDIAT.or.LimnoParam%aphy(gg)==LimnoParam%FDIAT) THEN
		        LimnoParam%tSiLisPomPhytS(iElem,gg) = (1.-LimnoParam%fDLisDomPhyt) * LimnoParam%tSiLisPhytS(iElem,gg)
		        LimnoParam%tSiLisDomPhytS(iElem,gg) = LimnoParam%fDLisDomPhyt * LimnoParam%tSiLisPhytS(iElem,gg)			
	        ELSE
		        LimnoParam%tSiLisPomPhytS(iElem,gg) = 0.
		        LimnoParam%tSiLisDomPhytS(iElem,gg) = 0.					
	        ENDIF
            
            
            
                !-------------------------------------------------------------------------------------!
                !                               Source/Sink Term for Phytoplankton in the sediment                       !
                !-------------------------------------------------------------------------------------!
            
			    LimnoParam%dVarEstS(iElem,1) =   LimnoParam%tDSetPhyt(HydroParam%ElSmallm(iElem),iElem,gg)          & !Sedimentation
                                                -LimnoParam%tDResusPhyt(iElem,gg)	                                & !Resuspension
                                                -LimnoParam%tDRespPhytS(iElem,gg)	                                & !Respiration
                                                -LimnoParam%tDLisPhytS(iElem,gg)	                                & !Mortality
                                                -SUM(LimnoParam%tDGrzBentPhyt(iElem,:,gg))                            !Grazing

			    LimnoParam%dVarEstS(iElem,2) =   LimnoParam%tDSetPhyt(HydroParam%ElSmallm(iElem),iElem,gg)*LimnoParam%rCDPhytW(HydroParam%ElSmallm(iElem),iElem,gg)          & !Sedimentation
                                                -LimnoParam%tCResusPhyt(iElem,gg)                                   & !Resuspension
                                                -LimnoParam%tCExcrPhytS(iElem,gg)	                                & !Respiration
                                                -LimnoParam%tCLisPhytS(iElem,gg)	                                & !Mortality
                                                -SUM(LimnoParam%tCGrzBentPhyt(iElem,:,gg))                            !Grazing
                                                
			    LimnoParam%dVarEstS(iElem,3) =   LimnoParam%tDSetPhyt(HydroParam%ElSmallm(iElem),iElem,gg)*LimnoParam%rNDPhytW(HydroParam%ElSmallm(iElem),iElem,gg)          & !Sedimentation
                                                -LimnoParam%tNResusPhyt(iElem,gg)	                                & !Resuspension
                                                -LimnoParam%tNExcrPhytS(iElem,gg)	                                & !Respiration
                                                -LimnoParam%tNLisPhytS(iElem,gg)	                                & !Mortality
                                                -SUM(LimnoParam%tNGrzBentPhyt(iElem,:,gg))                            !Grazing
			    
                LimnoParam%dVarEstS(iElem,4) =   LimnoParam%tDSetPhyt(HydroParam%ElSmallm(iElem),iElem,gg)*LimnoParam%rPDPhytW(HydroParam%ElSmallm(iElem),iElem,gg)          & !Sedimentation
                                                -LimnoParam%tPResusPhyt(iElem,gg)	                                & !Resuspension
                                                -LimnoParam%tPExcrPhytS(iElem,gg)	                                & !Respiration
                                                -LimnoParam%tPLisPhytS(iElem,gg)	                                & !Mortality
                                                -SUM(LimnoParam%tPGrzBentPhyt(iElem,:,gg))                            !Grazing
                                                
                LimnoParam%sDPhytS(iElem,gg) = MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,1)+LimnoParam%sDPhytSP(iElem,gg)))
                LimnoParam%sCPhytS(iElem,gg) = MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,2)+LimnoParam%sCPhytSP(iElem,gg)))
                LimnoParam%sNPhytS(iElem,gg) = MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,3)+LimnoParam%sNPhytSP(iElem,gg)))
                LimnoParam%sPPhytS(iElem,gg) = MAX(NearZero,(dtday*LimnoParam%dVarEstS(iElem,4)+LimnoParam%sPPhytSP(iElem,gg)))
                
        EndDo !Loop Cell
        
        
    EndDo !Loop Functional Group
      

Return
End
