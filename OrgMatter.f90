Subroutine OrgMatter(HydroParam,MeshParam,LimnoParam,dt,dtday)

    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars

    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    Integer:: index,iElem,iLayer,iter,om,bb
    Real:: wDSetPom_Layer,wCSetPom_Layer,wNSetPom_Layer,wPSetPom_Layer
    Real:: wDDifDom_Layer,wCDifDom_Layer,wNDifDom_Layer,wPDifDom_Layer
    Real:: wDResusPom_bottom,wCResusPom_bottom,wNResusPom_bottom,wPResusPom_bottom
    Real:: wDEgesBirdPom_Layer,wCEgesBirdPom_Layer,wNEgesBirdPom_Layer,wPEgesBirdPom_Layer
    Real:: dt,dtday
    Real:: V, sumbac
    Real:: NearZero = 1e-10
    
    Do om = 1, LimnoParam%npom ! =1 -> Labile; =2 -> Refractory
        
        !-------------------------------------------------------------------------------------!
        !                               Boundary Condition for POM				                !
        !-------------------------------------------------------------------------------------!
            
        If (LimnoParam%InclPomDomSplit == 1) Then
            If (om==LimnoParam%LABIL) Then
                index = 0
            ElseIf(om==LimnoParam%REFRATARIA) Then
                index = 4
                HydroParam%uLoadVarEst = LimnoParam%uDLoadPom
            EndIf 
        Else 
            index = 4
            HydroParam%uLoadVarEst = LimnoParam%uDLoadPom
        EndIf
        
        Do iElem = 1,MeshParam%nElem
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)

                !-------------------------------------------------------------------------------------!
                !                               POM Processes in Water					      !
                !-------------------------------------------------------------------------------------!
                !-<POM>---------------------------------------------!
                !(1) Hidrólise agua - POM => DOM					!
                !(5) Sedimentacao									!
                !(6) Contribuicao mortalidade fito      			!
                !(7) Contribuicao/Consumo do Zoops 	    			!
                !(8) Contribuicao/Consumo do Peixes 				!
                !(2) Ressuspensao(Sedimento)						!
                !(3) Erosao (Sedimento)								!
                !(4) Enterro (Sedimento)							!
                !---------------------------------------------------!

			    !-------------------------------------------------------------------------!
			    ! (1) Hidrólise - POM => DOM											  !
			    !-------------------------------------------------------------------------!
                If (LimnoParam%InclMatOrgSplit == 1) Then
			        ! (1.1) Limitacao da temperatura
                    LimnoParam%aTLimHidPom(om) = LimnoParam%cThetaHidW ** (LimnoParam%sDTempW(iLayer,iElem)-20.)
        
			        ! (1.2) Limitacao de OD
                    LimnoParam%aO2LimHidPom(om) = LimnoParam%sO2W(iLayer,iElem)/(LimnoParam%sO2W(iLayer,iElem)+LimnoParam%hO2Bac+NearZero) + LimnoParam%fbacAn*(LimnoParam%kbacAn/(LimnoParam%sO2W(iLayer,iElem)+LimnoParam%kbacAn+NearZero))
    			
			        ! (1.3) Inluencia da biomassa de bacterias
			        If (LimnoParam%InclBac == 1) Then
				        LimnoParam%aBacLimHidPom(om) = SUM(LimnoParam%sDBacW(iLayer,iElem,:)) / (SUM(LimnoParam%sDBacW(iLayer,iElem,:)) + LimnoParam%hBacHid+NearZero) 
			        Else
				        LimnoParam%aBacLimHidPom(om) = 1.
			        EndIf
    	
			        ! (1.4) Influencia do pH
			        If (LimnoParam%spHW(iLayer,iElem).GT.LimnoParam%pHmax) Then
				        LimnoParam%apHLimHidPom(om) = exp(LimnoParam%pHmax - LimnoParam%spHW(iLayer,iElem))
			        ElseIf (LimnoParam%spHW(iLayer,iElem).LT.LimnoParam%pHmin) Then
				        LimnoParam%apHLimHidPom(om) = exp(LimnoParam%spHW(iLayer,iElem) - LimnoParam%pHmin)
			        Else
				        LimnoParam%apHLimHidPom(om) = 1.0
			        EndIf
    			
			        ! (1.5) Taxa de Hidrolise corrigida
                    LimnoParam%kCorHidPom(om) = LimnoParam%kHidPomW(om) * LimnoParam%aTLimHidPom(om) * LimnoParam%aO2LimHidPom(om) * LimnoParam%aBacLimHidPom(om) * LimnoParam%apHLimHidPom(om)

			        ! (1.6) Calcula a transferencia POC => DOC 
                    LimnoParam%wDHidPomW(iLayer,iElem,om) = LimnoParam%kCorHidPom(om) * LimnoParam%sDPomW(iLayer,iElem,om)	
			        ! (1.7) Calcula os compartimentos de nutrientes C, N e P  na POM 
                    LimnoParam%wCHidPomW(iLayer,om) = LimnoParam%rCDPomW(iLayer,iElem,om) * LimnoParam%wDHidPomW(iLayer,iElem,om)
                    LimnoParam%wNHidPomW(iLayer,om) = LimnoParam%rNDPomW(iLayer,iElem,om) * LimnoParam%wDHidPomW(iLayer,iElem,om)
                    LimnoParam%wPHidPomW(iLayer,om) = LimnoParam%rPDPomW(iLayer,iElem,om) * LimnoParam%wDHidPomW(iLayer,iElem,om)
                Else
                    LimnoParam%wDHidPomW(iLayer,iElem,om) = 0.0
                    LimnoParam%wCHidPomW(iLayer,om) = 0.0
                    LimnoParam%wNHidPomW(iLayer,om) = 0.0
                    LimnoParam%wPHidPomW(iLayer,om) = 0.0
                EndIf
                
                
			    
                !------------------------------------------------------------------!
                ! 2) Mineralizacao Aerobica - DOM -> CO2,PO4,NH4
                    ! (2.1) Limitacao da temperatura
                    LimnoParam%aTLimMinAer(om)  = LimnoParam%cThetaMinAerW ** (LimnoParam%sDTempW(iLayer,iElem)-20.)
                    ! (2.2) Limitacao de OD
                    LimnoParam%aO2LimMinAer(om) = LimnoParam%sO2W(iLayer,iElem)/(LimnoParam%sO2W(iLayer,iElem)+LimnoParam%hO2Bac)
		            ! (2.3) Limitacao de pH
	                If(LimnoParam%spHW(iLayer,iElem).GT.LimnoParam%pHmax) Then
		                LimnoParam%apHLimMinAer(om) = exp(LimnoParam%pHmax - LimnoParam%spHW(iLayer,iElem))
	                ElseIf (LimnoParam%spHW(iLayer,iElem).LT.LimnoParam%pHmin) Then
		                LimnoParam%apHLimMinAer(om) = exp(LimnoParam%spHW(iLayer,iElem) - LimnoParam%pHmin)
	                Else
		                LimnoParam%apHLimMinAer(om) = 1.0
                    EndIf       
                    
                    ! (2.4) Limitacao de biomassa de bacterioplankton
                    If (LimnoParam%InclBac == 0) Then
                        LimnoParam%aBacLimMinAer(om) = 1.
                    Else
                        sumbac = 0.
                        Do bb = 1, LimnoParam%nbac
                            If (LimnoParam%abac(bb)==1) Then !Aerobic bacterioplankton
                                sumbac = sumbac + LimnoParam%sDBacW(iLayer,iElem,bb)
                            EndIf
                        EndDo
                        LimnoParam%aBacLimMinAer(om) = sumbac/(sumbac+LimnoParam%hBacMin(om))
                    EndIf
                    
                    LimnoParam%wDMinAerW(iLayer,iElem,om) = LimnoParam%kMinAerDomW(om) * LimnoParam%apHLimMinAer(om) * LimnoParam%aO2LimMinAer(om) * LimnoParam%aTLimMinAer(om) * LimnoParam%aBacLimMinAer(om) * LimnoParam%sDDomW(iLayer,iElem,om)	
                    LimnoParam%wCMinAerW(iLayer,om) = LimnoParam%rCDDomW(iLayer,iElem,om) * LimnoParam%wDMinAerW(iLayer,iElem,om)
                    LimnoParam%wNMinAerW(iLayer,om) = LimnoParam%rNDDomW(iLayer,iElem,om) * LimnoParam%wDMinAerW(iLayer,iElem,om)
                    LimnoParam%wPMinAerW(iLayer,om) = LimnoParam%rPDDomW(iLayer,iElem,om) * LimnoParam%wDMinAerW(iLayer,iElem,om)         
    			    !LimnoParam%wSiMinAerW(iLayer,pom) = LimnoParam%cSiDDomRef(pom) * LimnoParam%wDMinAerW(iLayer,iElem,pom)

                    LimnoParam%wCMinAer2CO2W(iLayer,iElem,om) = LimnoParam%wCMinAerW(iLayer,om)
                    LimnoParam%wCMinAer2CH4W(iLayer,om) = 0.0

                !------------------------------------------------------------------!
                ! 2) Mineralizacao Anaerobica - DOM -> CH4,CO2,PO4,NH4
                    ! (2.1) Limitacao da temperatura
                    LimnoParam%aTLimMinAnaer(om)  = LimnoParam%cThetaMinAnaerW ** (LimnoParam%sDTempW(iLayer,iElem)-20.)
                    ! (2.2) Limitacao de OD
                    LimnoParam%aO2LimMinAnaer(om) = LimnoParam%kbacAn/(LimnoParam%sO2W(iLayer,iElem)+LimnoParam%kbacAn)
		            ! (2.3) Limitacao de pH
	                If(LimnoParam%spHW(iLayer,iElem).GT.LimnoParam%pHmax) Then
		                LimnoParam%apHLimMinAnaer(om) = exp(LimnoParam%pHmax - LimnoParam%spHW(iLayer,iElem))
	                ElseIf (LimnoParam%spHW(iLayer,iElem).LT.LimnoParam%pHmin) Then
		                LimnoParam%apHLimMinAnaer(om) = exp(LimnoParam%spHW(iLayer,iElem) - LimnoParam%pHmin)
	                Else
		                LimnoParam%apHLimMinAnaer(om) = 1.0
                    EndIf
                    ! (2.4) Limitacao de biomassa de bacterioplankton
                    If (LimnoParam%InclBac == 0) Then
                        LimnoParam%aBacLimMinAnaer(om) = 1.
                    Else
                        sumbac = 0.
                        Do bb = 1, LimnoParam%nbac
                            If (LimnoParam%abac(bb)==2) Then !Aerobic bacterioplankton
                                sumbac = sumbac + LimnoParam%sDBacW(iLayer,iElem,bb)
                            EndIf
                        EndDo
                        LimnoParam%aBacLimMinAnaer(om) = sumbac/(sumbac+LimnoParam%hBacMin(om))
                    EndIf
	                                       
                    LimnoParam%wDMinAnaerW(iLayer,iElem,om) = LimnoParam%kMinAnaerDomW(om) * LimnoParam%apHLimMinAnaer(om) * LimnoParam%aO2LimMinAnaer(om) * LimnoParam%aTLimMinAnaer(om) * LimnoParam%aBacLimMinAnaer(om) * LimnoParam%sDDomW(iLayer,iElem,om)	
                    
                    LimnoParam%wCMinAnaerW(iLayer,om) = LimnoParam%rCDDomW(iLayer,iElem,om) * LimnoParam%wDMinAnaerW(iLayer,iElem,om)
                    LimnoParam%wNMinAnaerW(iLayer,om) = LimnoParam%rNDDomW(iLayer,iElem,om) * LimnoParam%wDMinAnaerW(iLayer,iElem,om)
                    LimnoParam%wPMinAnaerW(iLayer,om) = LimnoParam%rPDDomW(iLayer,iElem,om) * LimnoParam%wDMinAnaerW(iLayer,iElem,om)

                
                
                            
			    !-------------------------------------------------------------------------!
			    ! (2) Sedimentation											  !
			    !-------------------------------------------------------------------------!
                If ( iLayer == HydroParam%ElCapitalM(iElem) ) Then        ! Bottom layer or intermediary layer has contribuition of the up layer
                    wDSetPom_Layer = LimnoParam%tDSetPom(iLayer,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wCSetPom_Layer = LimnoParam%tCSetPom(iLayer,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wPSetPom_Layer = LimnoParam%tPSetPom(iLayer,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wNSetPom_Layer = LimnoParam%tNSetPom(iLayer,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                Else
                    wDSetPom_Layer = LimnoParam%tDSetPom(iLayer,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem)) - LimnoParam%tDSetPom(iLayer+1,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer+1,iElem))
                    wCSetPom_Layer = LimnoParam%tCSetPom(iLayer,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem)) - LimnoParam%tCSetPom(iLayer+1,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer+1,iElem))
                    wPSetPom_Layer = LimnoParam%tPSetPom(iLayer,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem)) - LimnoParam%tPSetPom(iLayer+1,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer+1,iElem))
                    wNSetPom_Layer = LimnoParam%tNSetPom(iLayer,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem)) - LimnoParam%tNSetPom(iLayer+1,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer+1,iElem))
                Endif
    		
			    !-------------------------------------------------------------------------!
			    ! (3) Erosao POM, onde: cDErosTot = mg/m2/dia                             !	
			    !-------------------------------------------------------------------------!
                If (LimnoParam%iSed==1) Then
		            If (MeshParam%Neighbor(1,iElem)==0.or.MeshParam%Neighbor(2,iElem)==0.or.MeshParam%Neighbor(3,iElem)==0.or.MeshParam%Neighbor(4,iElem)==0.and.om==LimnoParam%REFRATARIA) Then !
                        !  organic_matter_input_from_banks
			            LimnoParam%tDErosPom(iLayer,iElem,om) = LimnoParam%fDOrgSoil*LimnoParam%cDErosTot
                        !  organic_C_input_from_banks
			            LimnoParam%tCErosPom(iLayer,om) = LimnoParam%cCDPomRef(om)*LimnoParam%tDErosPom(iLayer,iElem,om)
                        !  organic_N_input_from_banks
			            LimnoParam%tNErosPom(iLayer,om) = LimnoParam%cNDPomRef(om)*LimnoParam%tDErosPom(iLayer,iElem,om)
                        !  organic_P_input_from_banks
			            LimnoParam%tPErosPom(iLayer,om) = LimnoParam%cPDPomRef(om)*LimnoParam%tDErosPom(iLayer,iElem,om)
		            Else
			            LimnoParam%tDErosPom(iLayer,iElem,om) = 0.0
			            LimnoParam%tCErosPom(iLayer,om) = 0.0
			            LimnoParam%tNErosPom(iLayer,om) = 0.0
			            LimnoParam%tPErosPom(iLayer,om) = 0.0
                    EndIf
                EndIf
                
			    !-------------------------------------------------------------------------!
			    ! (4)  Egestion by Birds                                                    !	
			    !-------------------------------------------------------------------------!
                If ( iLayer == HydroParam%ElCapitalM(iElem) ) Then        ! Bottom layer or intermediary layer has contribuition of the up layer
                    wDEgesBirdPom_Layer = SUM(LimnoParam%tDEgesBird(iElem,:))/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wCEgesBirdPom_Layer = SUM(LimnoParam%tCEgesBird(iElem,:))/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wNEgesBirdPom_Layer = SUM(LimnoParam%tNEgesBird(iElem,:))/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wPEgesBirdPom_Layer = SUM(LimnoParam%tPEgesBird(iElem,:))/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                Else
                    wDEgesBirdPom_Layer = 0.
                    wCEgesBirdPom_Layer = 0.
                    wNEgesBirdPom_Layer = 0.
                    wPEgesBirdPom_Layer = 0.
                EndIf
                
			    !-------------------------------------------------------------------------!
			    ! (5)  DOM Difusion                                                    !	
			    !-------------------------------------------------------------------------!
                If ( iLayer == HydroParam%ElSmallm(iElem) ) Then
                    wDDifDom_Layer = LimnoParam%tDDifDom(iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wCDifDom_Layer = LimnoParam%tCDifDom(iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wNDifDom_Layer = LimnoParam%tNDifDom(iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wPDifDom_Layer = LimnoParam%tPDifDom(iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                Else
                    wDDifDom_Layer = 0.
                    wCDifDom_Layer = 0.
                    wNDifDom_Layer = 0.
                    wPDifDom_Layer = 0.
                EndIf
			    !-------------------------------------------------------------------------!
			    ! (6)  POM Resuspension                                                    !	
			    !-------------------------------------------------------------------------!
                If ( iLayer == HydroParam%ElSmallm(iElem) ) Then
                    wDResusPom_bottom = LimnoParam%tDResusPom(iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wCResusPom_bottom = LimnoParam%tCResusPom(iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wNResusPom_bottom = LimnoParam%tNResusPom(iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                    wPResusPom_bottom = LimnoParam%tPResusPom(iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
                Else
                    wDResusPom_bottom = 0.
                    wCResusPom_bottom = 0.
                    wNResusPom_bottom = 0.
                    wPResusPom_bottom = 0.
                EndIf
                
                
                
                If (LimnoParam%InclMatOrgSplit == 1) Then ! Case 1: Splited Organic Matter (POM and DOM)
                    
                    If (LimnoParam%InclPomDomSplit == 1) Then ! Case 2: Splited Organic Matter (POM and DOM) and Splited fractions (Labile and Refratory)
                        If (om==LimnoParam%LABIL) Then
                            !-------------------------------------------------------------------------------------!
                            !                               Source/Sink Term for POM Labile                       !
                            !-------------------------------------------------------------------------------------!
                            HydroParam%dVarEst(iLayer,iElem,1) = LimnoParam%fDLisDomPhyt*SUM(LimnoParam%wDLisPhyt(iLayer,iElem,:))              & !Lise Fito
		                                    + LimnoParam%fDMessDomZoo*SUM(LimnoParam%wDMessZooOM(iLayer,iElem,:))                               & !Messy Feeding Zoops
		                                    + LimnoParam%fDEgesDomZoo*SUM(LimnoParam%wDEgesZoo(iLayer,iElem,:))                                 & !Pellets Zoops
                                            + LimnoParam%fDMortDomZoo*SUM(LimnoParam%wDMortZoo(iLayer,iElem,:))                                 & !Mortalidade Zoops
 	                                        - SUM(LimnoParam%wDGrzZooPom(iLayer,iElem,:,om))                                                   & !Grazing Zoops
	                                        + LimnoParam%fDMessDomFish*SUM(LimnoParam%wDMessFiAdOM(iLayer,iElem,:))                             & !Messy Feeding Adult Fish
                                            + LimnoParam%fDMessDomFish*SUM(LimnoParam%wDMessFiJvOM(iLayer,iElem,:))                             & !Messy Feeding Juveline Fish
	                                        + LimnoParam%fDEgesDomFish*SUM(LimnoParam%wDEgesFiAd(iLayer,iElem,:))                               & !Aduld Fish Pellets-> POM
                                            + LimnoParam%fDEgesDomFish*SUM(LimnoParam%wDEgesFiJv(iLayer,iElem,:))                               & !Juveline Fish Pellets-> POM
	                                        + LimnoParam%fDMortDomFish*SUM(LimnoParam%wDMortFiAd(iLayer,iElem,:))                               & !Aduld Fish Mortality-> POM
                                            + LimnoParam%fDMortDomFish*SUM(LimnoParam%wDMortFiJv(iLayer,iElem,:))                               & !Juveline Fish Mortality-> POM
 	                                        - SUM(LimnoParam%wDGrzFiAdPom(iLayer,iElem,:,om))                                                  & !Grazing Aduld Fish
                                            - SUM(LimnoParam%wDGrzFiJvPom(iLayer,iElem,:,om))                                                  & !Grazing Juveline Fish
                                            + LimnoParam%fDMortDomMac*SUM(LimnoParam%tDMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                                    + LimnoParam%fDAssDomBird*wDEgesBirdPom_Layer                                                       & !Egestao Bird -> POM                    
                                            - LimnoParam%wDHidPomW(iLayer,iElem,om)                                                             & ! Hidrólise - POM => DOM
                                            + wDResusPom_bottom                                                                                 ! Resuspension
  
                            HydroParam%dVarEst(iLayer,iElem,2) = LimnoParam%fDLisDomPhyt*SUM(LimnoParam%wCLisPhyt(iLayer,iElem,:))              & !Lise Fito
		                                    + LimnoParam%fDMessDomZoo*SUM(LimnoParam%wCMessZooOM(iLayer,iElem,:))                               & !Messy Feeding Zoops
		                                    + LimnoParam%fDEgesDomZoo*SUM(LimnoParam%wCEgesZoo(iLayer,iElem,:))                                 & !Pellets Zoops
                                            + LimnoParam%fDMortDomZoo*SUM(LimnoParam%wCMortZoo(iLayer,iElem,:))                                 & !Mortalidade Zoops
 	                                        - SUM(LimnoParam%wCGrzZooPom(iLayer,iElem,:,om))                                                   & !Grazing Zoops
	                                        + LimnoParam%fDMessDomFish*SUM(LimnoParam%wCMessFiAdOM(iLayer,iElem,:))                             & !Messy Feeding Adult Fish
                                            + LimnoParam%fDMessDomFish*SUM(LimnoParam%wCMessFiJvOM(iLayer,iElem,:))                             & !Messy Feeding Juveline Fish
	                                        + LimnoParam%fDEgesDomFish*SUM(LimnoParam%wCEgesFiAd(iLayer,iElem,:))                               & !Aduld Fish Pellets-> POM
                                            + LimnoParam%fDEgesDomFish*SUM(LimnoParam%wCEgesFiJv(iLayer,iElem,:))                               & !Juveline Fish Pellets-> POM
	                                        + LimnoParam%fDMortDomFish*SUM(LimnoParam%wCMortFiAd(iLayer,iElem,:))                               & !Aduld Fish Mortality-> POM
                                            + LimnoParam%fDMortDomFish*SUM(LimnoParam%wCMortFiJv(iLayer,iElem,:))                               & !Juveline Fish Mortality-> POM
 	                                        - SUM(LimnoParam%wCGrzFiAdPom(iLayer,iElem,:,om))                                                  & !Grazing Aduld Fish
                                            - SUM(LimnoParam%wCGrzFiJvPom(iLayer,iElem,:,om))                                                  & !Grazing Juveline Fish
                                            + LimnoParam%fDMortDomMac*SUM(LimnoParam%tCMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                                    + LimnoParam%fDAssDomBird*wCEgesBirdPom_Layer                                                       & !Egestao Bird -> POM                    
                                            - LimnoParam%wCHidPomW(iLayer,om)                                                                  & ! Hidrólise - POM => DOM
                                            + wCResusPom_bottom                                                                                 ! Resuspension
 
                            HydroParam%dVarEst(iLayer,iElem,3) = LimnoParam%fDLisDomPhyt*SUM(LimnoParam%wNLisPhyt(iLayer,iElem,:))              & !Lise Fito
		                                    + LimnoParam%fDMessDomZoo*SUM(LimnoParam%wNMessZooOM(iLayer,iElem,:))                               & !Messy Feeding Zoops
		                                    + LimnoParam%fDEgesDomZoo*SUM(LimnoParam%wNEgesZoo(iLayer,iElem,:))                                 & !Pellets Zoops
                                            + LimnoParam%fDMortDomZoo*SUM(LimnoParam%wNMortZoo(iLayer,iElem,:))                                 & !Mortalidade Zoops
 	                                        - SUM(LimnoParam%wNGrzZooPom(iLayer,iElem,:,om))                                                   & !Grazing Zoops
	                                        + LimnoParam%fDMessDomFish*SUM(LimnoParam%wNMessFiAdOM(iLayer,iElem,:))                             & !Messy Feeding Adult Fish
                                            + LimnoParam%fDMessDomFish*SUM(LimnoParam%wNMessFiJvOM(iLayer,iElem,:))                             & !Messy Feeding Juveline Fish
	                                        + LimnoParam%fDEgesDomFish*SUM(LimnoParam%wNEgesFiAd(iLayer,iElem,:))                               & !Aduld Fish Pellets-> POM
                                            + LimnoParam%fDEgesDomFish*SUM(LimnoParam%wNEgesFiJv(iLayer,iElem,:))                               & !Juveline Fish Pellets-> POM
	                                        + LimnoParam%fDMortDomFish*SUM(LimnoParam%wNMortFiAd(iLayer,iElem,:))                               & !Aduld Fish Mortality-> POM
                                            + LimnoParam%fDMortDomFish*SUM(LimnoParam%wNMortFiJv(iLayer,iElem,:))                               & !Juveline Fish Mortality-> POM
 	                                        - SUM(LimnoParam%wNGrzFiAdPom(iLayer,iElem,:,om))                                                  & !Grazing Aduld Fish
                                            - SUM(LimnoParam%wNGrzFiJvPom(iLayer,iElem,:,om))                                                  & !Grazing Juveline Fish
                                            + LimnoParam%fDMortDomMac*SUM(LimnoParam%tNMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                                    + LimnoParam%fDAssDomBird*wNEgesBirdPom_Layer                                                       & !Egestao Bird -> POM                    
                                            - LimnoParam%wNHidPomW(iLayer,om)                                                                  & ! Hidrólise - POM => DOM
                                            + wNResusPom_bottom                                                                                 ! Resuspension
 
                            HydroParam%dVarEst(iLayer,iElem,4) = LimnoParam%fDLisDomPhyt*SUM(LimnoParam%wPLisPhyt(iLayer,iElem,:))              & !Lise Fito
		                                    + LimnoParam%fDMessDomZoo*SUM(LimnoParam%wPMessZooOM(iLayer,iElem,:))                               & !Messy Feeding Zoops
		                                    + LimnoParam%fDEgesDomZoo*SUM(LimnoParam%wPEgesZoo(iLayer,iElem,:))                                 & !Pellets Zoops
                                            + LimnoParam%fDMortDomZoo*SUM(LimnoParam%wPMortZoo(iLayer,iElem,:))                                 & !Mortalidade Zoops
 	                                        - SUM(LimnoParam%wPGrzZooPom(iLayer,iElem,:,om))                                                   & !Grazing Zoops
	                                        + LimnoParam%fDMessDomFish*SUM(LimnoParam%wPMessFiAdOM(iLayer,iElem,:))                             & !Messy Feeding Adult Fish
                                            + LimnoParam%fDMessDomFish*SUM(LimnoParam%wPMessFiJvOM(iLayer,iElem,:))                             & !Messy Feeding Juveline Fish
	                                        + LimnoParam%fDEgesDomFish*SUM(LimnoParam%wPEgesFiAd(iLayer,iElem,:))                               & !Aduld Fish Pellets-> POM
                                            + LimnoParam%fDEgesDomFish*SUM(LimnoParam%wPEgesFiJv(iLayer,iElem,:))                               & !Juveline Fish Pellets-> POM
	                                        + LimnoParam%fDMortDomFish*SUM(LimnoParam%wPMortFiAd(iLayer,iElem,:))                               & !Aduld Fish Mortality-> POM
                                            + LimnoParam%fDMortDomFish*SUM(LimnoParam%wPMortFiJv(iLayer,iElem,:))                               & !Juveline Fish Mortality-> POM
 	                                        - SUM(LimnoParam%wPGrzFiAdPom(iLayer,iElem,:,om))                                                  & !Grazing Aduld Fish
                                            - SUM(LimnoParam%wPGrzFiJvPom(iLayer,iElem,:,om))                                                  & !Grazing Juveline Fish
                                            + LimnoParam%fDMortDomMac*SUM(LimnoParam%tPMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                                    + LimnoParam%fDAssDomBird*wPEgesBirdPom_Layer                                                       & !Egestao Bird -> POM                    
                                            - LimnoParam%wPHidPomW(iLayer,om)                                                                  & ! Hidrólise - POM => DOM
                                            + wPResusPom_bottom                                                                                 ! Resuspension
                            
                        ElseIf(om==LimnoParam%REFRATARIA) Then
                            !-------------------------------------------------------------------------------------!
                            !                               Source/Sink Term for POM Refratory                    !
                            !-------------------------------------------------------------------------------------!
                            HydroParam%dVarEst(iLayer,iElem,1) = - wDSetPom_Layer                                                               & !Sedimentacao
		                                    + LimnoParam%tDErosPom(iLayer,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Erosao Margem
                                            !+ SUM(LimnoParam%wDMortBac(iLayer,iElem,:))                                                         & !Lise Bacterioplankton
		                                    + (1.-LimnoParam%fDLisDomPhyt)*SUM(LimnoParam%wDLisPhyt(iLayer,iElem,:))                            & !Lise Fito
		                                    + (1.-LimnoParam%fDMessDomZoo)*SUM(LimnoParam%wDMessZooOM(iLayer,iElem,:))                          & !Messy Feeding Zoops
		                                    + (1.-LimnoParam%fDEgesDomZoo)*SUM(LimnoParam%wDEgesZoo(iLayer,iElem,:))                            & !Pellets Zoops
                                            + (1.-LimnoParam%fDMortDomZoo)*SUM(LimnoParam%wDMortZoo(iLayer,iElem,:))                            & !Mortalidade Zoops
	                                        + (1.-LimnoParam%fDMessDomFish)*SUM(LimnoParam%wDMessFiAdOM(iLayer,iElem,:))                        & !Messy Feeding Adult Fish
                                            + (1.-LimnoParam%fDMessDomFish)*SUM(LimnoParam%wDMessFiJvOM(iLayer,iElem,:))                        & !Messy Feeding Juveline Fish
	                                        + (1.-LimnoParam%fDEgesDomFish)*SUM(LimnoParam%wDEgesFiAd(iLayer,iElem,:))                          & !Aduld Fish Pellets-> POM
                                            + (1.-LimnoParam%fDEgesDomFish)*SUM(LimnoParam%wDEgesFiJv(iLayer,iElem,:))                          & !Juveline Fish Pellets-> POM
	                                        + (1.-LimnoParam%fDMortDomFish)*SUM(LimnoParam%wDMortFiAd(iLayer,iElem,:))                          & !Aduld Fish Mortality-> POM
                                            + (1.-LimnoParam%fDMortDomFish)*SUM(LimnoParam%wDMortFiJv(iLayer,iElem,:))                          & !Juveline Fish Mortality-> POM
                                            + (1.-LimnoParam%fDMortDomMac)*SUM(LimnoParam%tDMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                                    + (1.-LimnoParam%fDAssDomBird)*wDEgesBirdPom_Layer                                                  & !Egestao Bird -> POM                    
                                            - LimnoParam%wDHidPomW(iLayer,iElem,om)                                                            ! Hidrólise - POM => DOM
  
                            HydroParam%dVarEst(iLayer,iElem,2) = - wCSetPom_Layer                                                               & !Sedimentacao
		                                    + LimnoParam%tCErosPom(iLayer,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Erosao Margem
                                            !+ SUM(LimnoParam%wDMortBac(iLayer,iElem,:))                                                         & !Lise Bacterioplankton
		                                    + (1.-LimnoParam%fDLisDomPhyt)*SUM(LimnoParam%wCLisPhyt(iLayer,iElem,:))                            & !Lise Fito
		                                    + (1.-LimnoParam%fDMessDomZoo)*SUM(LimnoParam%wCMessZooOM(iLayer,iElem,:))                          & !Messy Feeding Zoops
		                                    + (1.-LimnoParam%fDEgesDomZoo)*SUM(LimnoParam%wCEgesZoo(iLayer,iElem,:))                            & !Pellets Zoops
                                            + (1.-LimnoParam%fDMortDomZoo)*SUM(LimnoParam%wCMortZoo(iLayer,iElem,:))                            & !Mortalidade Zoops
	                                        + (1.-LimnoParam%fDMessDomFish)*SUM(LimnoParam%wCMessFiAdOM(iLayer,iElem,:))                        & !Messy Feeding Adult Fish
                                            + (1.-LimnoParam%fDMessDomFish)*SUM(LimnoParam%wCMessFiJvOM(iLayer,iElem,:))                        & !Messy Feeding Juveline Fish
	                                        + (1.-LimnoParam%fDEgesDomFish)*SUM(LimnoParam%wCEgesFiAd(iLayer,iElem,:))                          & !Aduld Fish Pellets-> POM
                                            + (1.-LimnoParam%fDEgesDomFish)*SUM(LimnoParam%wCEgesFiJv(iLayer,iElem,:))                          & !Juveline Fish Pellets-> POM
	                                        + (1.-LimnoParam%fDMortDomFish)*SUM(LimnoParam%wCMortFiAd(iLayer,iElem,:))                          & !Aduld Fish Mortality-> POM
                                            + (1.-LimnoParam%fDMortDomFish)*SUM(LimnoParam%wCMortFiJv(iLayer,iElem,:))                          & !Juveline Fish Mortality-> POM
                                            + (1.-LimnoParam%fDMortDomMac)*SUM(LimnoParam%tCMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                                    + (1.-LimnoParam%fDAssDomBird)*wCEgesBirdPom_Layer                                                  & !Egestao Bird -> POM                    
                                            - LimnoParam%wCHidPomW(iLayer,om)                                                                  ! Hidrólise - POM => DOM
 
                            HydroParam%dVarEst(iLayer,iElem,3) = - wNSetPom_Layer                                                               & !Sedimentacao
		                                    + LimnoParam%tNErosPom(iLayer,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Erosao Margem
                                            !+ SUM(LimnoParam%wDMortBac(iLayer,iElem,:))                                                         & !Lise Bacterioplankton
		                                    + (1.-LimnoParam%fDLisDomPhyt)*SUM(LimnoParam%wNLisPhyt(iLayer,iElem,:))                            & !Lise Fito
		                                    + (1.-LimnoParam%fDMessDomZoo)*SUM(LimnoParam%wNMessZooOM(iLayer,iElem,:))                          & !Messy Feeding Zoops
		                                    + (1.-LimnoParam%fDEgesDomZoo)*SUM(LimnoParam%wNEgesZoo(iLayer,iElem,:))                            & !Pellets Zoops
                                            + (1.-LimnoParam%fDMortDomZoo)*SUM(LimnoParam%wNMortZoo(iLayer,iElem,:))                            & !Mortalidade Zoops
	                                        + (1.-LimnoParam%fDMessDomFish)*SUM(LimnoParam%wNMessFiAdOM(iLayer,iElem,:))                        & !Messy Feeding Adult Fish
                                            + (1.-LimnoParam%fDMessDomFish)*SUM(LimnoParam%wNMessFiJvOM(iLayer,iElem,:))                        & !Messy Feeding Juveline Fish
	                                        + (1.-LimnoParam%fDEgesDomFish)*SUM(LimnoParam%wNEgesFiAd(iLayer,iElem,:))                          & !Aduld Fish Pellets-> POM
                                            + (1.-LimnoParam%fDEgesDomFish)*SUM(LimnoParam%wNEgesFiJv(iLayer,iElem,:))                          & !Juveline Fish Pellets-> POM
	                                        + (1.-LimnoParam%fDMortDomFish)*SUM(LimnoParam%wNMortFiAd(iLayer,iElem,:))                          & !Aduld Fish Mortality-> POM
                                            + (1.-LimnoParam%fDMortDomFish)*SUM(LimnoParam%wNMortFiJv(iLayer,iElem,:))                          & !Juveline Fish Mortality-> POM
                                            + (1.-LimnoParam%fDMortDomMac)*SUM(LimnoParam%tNMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                                    + (1.-LimnoParam%fDAssDomBird)*wNEgesBirdPom_Layer                                                  & !Egestao Bird -> POM                    
                                            - LimnoParam%wNHidPomW(iLayer,om)                                                                  ! Hidrólise - POM => DOM

                            HydroParam%dVarEst(iLayer,iElem,3) = - wPSetPom_Layer                                                               & !Sedimentacao
		                                    + LimnoParam%tPErosPom(iLayer,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Erosao Margem
                                            !+ SUM(LimnoParam%wDMortBac(iLayer,iElem,:))                                                         & !Lise Bacterioplankton
		                                    + (1.-LimnoParam%fDLisDomPhyt)*SUM(LimnoParam%wPLisPhyt(iLayer,iElem,:))                            & !Lise Fito
		                                    + (1.-LimnoParam%fDMessDomZoo)*SUM(LimnoParam%wPMessZooOM(iLayer,iElem,:))                          & !Messy Feeding Zoops
		                                    + (1.-LimnoParam%fDEgesDomZoo)*SUM(LimnoParam%wPEgesZoo(iLayer,iElem,:))                            & !Pellets Zoops
                                            + (1.-LimnoParam%fDMortDomZoo)*SUM(LimnoParam%wPMortZoo(iLayer,iElem,:))                            & !Mortalidade Zoops
	                                        + (1.-LimnoParam%fDMessDomFish)*SUM(LimnoParam%wPMessFiAdOM(iLayer,iElem,:))                        & !Messy Feeding Adult Fish
                                            + (1.-LimnoParam%fDMessDomFish)*SUM(LimnoParam%wPMessFiJvOM(iLayer,iElem,:))                        & !Messy Feeding Juveline Fish
	                                        + (1.-LimnoParam%fDEgesDomFish)*SUM(LimnoParam%wPEgesFiAd(iLayer,iElem,:))                          & !Aduld Fish Pellets-> POM
                                            + (1.-LimnoParam%fDEgesDomFish)*SUM(LimnoParam%wPEgesFiJv(iLayer,iElem,:))                          & !Juveline Fish Pellets-> POM
	                                        + (1.-LimnoParam%fDMortDomFish)*SUM(LimnoParam%wPMortFiAd(iLayer,iElem,:))                          & !Aduld Fish Mortality-> POM
                                            + (1.-LimnoParam%fDMortDomFish)*SUM(LimnoParam%wPMortFiJv(iLayer,iElem,:))                          & !Juveline Fish Mortality-> POM
                                            + (1.-LimnoParam%fDMortDomMac)*SUM(LimnoParam%tPMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                                    + (1.-LimnoParam%fDAssDomBird)*wPEgesBirdPom_Layer                                                  & !Egestao Bird -> POM                    
                                            - LimnoParam%wPHidPomW(iLayer,om)       
                        EndIf
                    Else ! Case 3: Splited Organic Matter (POM and DOM) and non-splited fractions
                        !-------------------------------------------------------------------------------------!
                        !                               Source/Sink Term for POM						      !
                        !-------------------------------------------------------------------------------------!
                        HydroParam%dVarEst(iLayer,iElem,1) = - wDSetPom_Layer                                                               & !Sedimentacao
		                                + LimnoParam%tDErosPom(iLayer,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Erosao Margem
                                        !+ SUM(LimnoParam%wDMortBac(iLayer,iElem,:))                                                         & !Lise Bacterioplankton
		                                + SUM(LimnoParam%wDLisPhyt(iLayer,iElem,:))                                                         & !Lise Fito
		                                + SUM(LimnoParam%wDMessZooOM(iLayer,iElem,:))                                                       & !Messy Feeding Zoops
		                                + SUM(LimnoParam%wDEgesZoo(iLayer,iElem,:))                                                         & !Pellets Zoops
                                        + SUM(LimnoParam%wDMortZoo(iLayer,iElem,:))                                                         & !Mortalidade Zoops
 	                                    - SUM(LimnoParam%wDGrzZooPom(iLayer,iElem,:,om))                                                   & !Grazing Zoops
	                                    + SUM(LimnoParam%wDMessFiAdOM(iLayer,iElem,:))                                                      & !Messy Feeding Adult Fish
                                        + SUM(LimnoParam%wDMessFiJvOM(iLayer,iElem,:))                                                      & !Messy Feeding Juveline Fish
	                                    + SUM(LimnoParam%wDEgesFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Pellets-> POM
                                        + SUM(LimnoParam%wDEgesFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Pellets-> POM
	                                    + SUM(LimnoParam%wDMortFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Mortality-> POM
                                        + SUM(LimnoParam%wDMortFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Mortality-> POM
 	                                    - SUM(LimnoParam%wDGrzFiAdPom(iLayer,iElem,:,om))                                                  & !Grazing Aduld Fish
                                        - SUM(LimnoParam%wDGrzFiJvPom(iLayer,iElem,:,om))                                                  & !Grazing Juveline Fish
                                        + SUM(LimnoParam%tDMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                                + wDEgesBirdPom_Layer                                                                               & !Egestao Bird -> POM                    
                                        - LimnoParam%wDHidPomW(iLayer,iElem,om)                                                            & ! Hidrólise - POM => DOM
                                        + wDResusPom_bottom                                                                                 ! Resuspension
                    
                        HydroParam%dVarEst(iLayer,iElem,2) = - wCSetPom_Layer                                                               & !Sedimentacao
		                                + LimnoParam%tCErosPom(iLayer,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Erosao Margem
                                        !+ SUM(LimnoParam%wDMortBac(iLayer,iElem,:))                                                         & !Lise Bacterioplankton
		                                + SUM(LimnoParam%wCLisPhyt(iLayer,iElem,:))                                                         & !Lise Fito
		                                + SUM(LimnoParam%wCMessZooOM(iLayer,iElem,:))                                                       & !Messy Feeding Zoops
		                                + SUM(LimnoParam%wCEgesZoo(iLayer,iElem,:))                                                         & !Pellets Zoops
                                        + SUM(LimnoParam%wCMortZoo(iLayer,iElem,:))                                                         & !Mortalidade Zoops
 	                                    - SUM(LimnoParam%wCGrzZooPom(iLayer,iElem,:,om))                                                   & !Grazing Zoops
	                                    + SUM(LimnoParam%wCMessFiAdOM(iLayer,iElem,:))                                                      & !Messy Feeding Adult Fish
                                        + SUM(LimnoParam%wCMessFiJvOM(iLayer,iElem,:))                                                      & !Messy Feeding Juveline Fish
	                                    + SUM(LimnoParam%wCEgesFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Pellets-> POM
                                        + SUM(LimnoParam%wCEgesFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Pellets-> POM
	                                    + SUM(LimnoParam%wCMortFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Mortality-> POM
                                        + SUM(LimnoParam%wCMortFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Mortality-> POM
 	                                    - SUM(LimnoParam%wCGrzFiAdPom(iLayer,iElem,:,om))                                                  & !Grazing Aduld Fish
                                        - SUM(LimnoParam%wCGrzFiJvPom(iLayer,iElem,:,om))                                                  & !Grazing Juveline Fish
                                        + SUM(LimnoParam%tCMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                                + wCEgesBirdPom_Layer                                                                               & !Egestao Bird -> POM                    
                                        - LimnoParam%wCHidPomW(iLayer,om)                                                                   & ! Hidrólise - POM => DOM
                                        + wCResusPom_bottom                                                                                 ! Resuspension
                        
                        HydroParam%dVarEst(iLayer,iElem,3) = - wNSetPom_Layer                                                               & !Sedimentacao
		                                + LimnoParam%tNErosPom(iLayer,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Erosao Margem
                                        !+ SUM(LimnoParam%wDMortBac(iLayer,iElem,:))                                                         & !Lise Bacterioplankton
		                                + SUM(LimnoParam%wNLisPhyt(iLayer,iElem,:))                                                         & !Lise Fito
		                                + SUM(LimnoParam%wNMessZooOM(iLayer,iElem,:))                                                       & !Messy Feeding Zoops
		                                + SUM(LimnoParam%wNEgesZoo(iLayer,iElem,:))                                                         & !Pellets Zoops
                                        + SUM(LimnoParam%wNMortZoo(iLayer,iElem,:))                                                         & !Mortalidade Zoops
 	                                    - SUM(LimnoParam%wNGrzZooPom(iLayer,iElem,:,om))                                                   & !Grazing Zoops
	                                    + SUM(LimnoParam%wNMessFiAdOM(iLayer,iElem,:))                                                      & !Messy Feeding Adult Fish
                                        + SUM(LimnoParam%wNMessFiJvOM(iLayer,iElem,:))                                                      & !Messy Feeding Juveline Fish
	                                    + SUM(LimnoParam%wNEgesFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Pellets-> POM
                                        + SUM(LimnoParam%wNEgesFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Pellets-> POM
	                                    + SUM(LimnoParam%wNMortFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Mortality-> POM
                                        + SUM(LimnoParam%wNMortFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Mortality-> POM
 	                                    - SUM(LimnoParam%wNGrzFiAdPom(iLayer,iElem,:,om))                                                  & !Grazing Aduld Fish
                                        - SUM(LimnoParam%wNGrzFiJvPom(iLayer,iElem,:,om))                                                  & !Grazing Juveline Fish
                                        + SUM(LimnoParam%tNMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                                + wNEgesBirdPom_Layer                                                                               & !Egestao Bird -> POM                    
                                        - LimnoParam%wNHidPomW(iLayer,om)                                                                   &! Hidrólise - POM => DOM
                                        + wNResusPom_bottom                                                                                 ! Resuspension
  
                        HydroParam%dVarEst(iLayer,iElem,4) = - wPSetPom_Layer                                                               & !Sedimentacao
		                                + LimnoParam%tPErosPom(iLayer,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Erosao Margem
                                        !+ SUM(LimnoParam%wDMortBac(iLayer,iElem,:))                                                         & !Lise Bacterioplankton
		                                + SUM(LimnoParam%wPLisPhyt(iLayer,iElem,:))                                                         & !Lise Fito
		                                + SUM(LimnoParam%wPMessZooOM(iLayer,iElem,:))                                                       & !Messy Feeding Zoops
		                                + SUM(LimnoParam%wPEgesZoo(iLayer,iElem,:))                                                         & !Pellets Zoops
                                        + SUM(LimnoParam%wPMortZoo(iLayer,iElem,:))                                                         & !Mortalidade Zoops
 	                                    - SUM(LimnoParam%wPGrzZooPom(iLayer,iElem,:,om))                                                   & !Grazing Zoops
	                                    + SUM(LimnoParam%wPMessFiAdOM(iLayer,iElem,:))                                                      & !Messy Feeding Adult Fish
                                        + SUM(LimnoParam%wPMessFiJvOM(iLayer,iElem,:))                                                      & !Messy Feeding Juveline Fish
	                                    + SUM(LimnoParam%wPEgesFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Pellets-> POM
                                        + SUM(LimnoParam%wPEgesFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Pellets-> POM
	                                    + SUM(LimnoParam%wPMortFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Mortality-> POM
                                        + SUM(LimnoParam%wPMortFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Mortality-> POM
 	                                    - SUM(LimnoParam%wpGrzFiAdPom(iLayer,iElem,:,om))                                                  & !Grazing Aduld Fish
                                        - SUM(LimnoParam%wPGrzFiJvPom(iLayer,iElem,:,om))                                                  & !Grazing Juveline Fish
                                        + SUM(LimnoParam%tPMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                                + wPEgesBirdPom_Layer                                                                              & !Egestao Bird -> POM                    
                                        - LimnoParam%wPHidPomW(iLayer,om)                                                                  & ! Hidrólise - POM => DOM
                                        + wPResusPom_bottom                                                                                 ! Resuspension
                    EndIf
                Else ! Case 4: Non-splited Organic Matter and non-splited fractions  (only one box of organic matter -  Detritus - PCLake)
                    !-------------------------------------------------------------------------------------!
                    !                               Source/Sink Term for POM						      !
                    !-------------------------------------------------------------------------------------!
                    HydroParam%dVarEst(iLayer,iElem,1) = - wDSetPom_Layer                                                               & !Sedimentacao
		                            + LimnoParam%tDErosPom(iLayer,iElem,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))          & !Erosao Margem
                                    + SUM(LimnoParam%wDMortBac(iLayer,iElem,:))                                                         & !Lise Bacterioplankton
		                            + SUM(LimnoParam%wDLisPhyt(iLayer,iElem,:))                                                         & !Lise Fito
		                            + SUM(LimnoParam%wDMessZooOM(iLayer,iElem,:))                                                       & !Messy Feeding Zoops
		                            + SUM(LimnoParam%wDEgesZoo(iLayer,iElem,:))                                                         & !Pellets Zoops
                                    + SUM(LimnoParam%wDMortZoo(iLayer,iElem,:))                                                         & !Mortalidade Zoops
 	                                - SUM(LimnoParam%wDGrzZooPom(iLayer,iElem,:,om))                                                   & !Grazing Zoops
	                                + SUM(LimnoParam%wDMessFiAdOM(iLayer,iElem,:))                                                      & !Messy Feeding Adult Fish
                                    + SUM(LimnoParam%wDMessFiJvOM(iLayer,iElem,:))                                                      & !Messy Feeding Juveline Fish
	                                + SUM(LimnoParam%wDEgesFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Pellets-> POM
                                    + SUM(LimnoParam%wDEgesFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Pellets-> POM
	                                + SUM(LimnoParam%wDMortFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Mortality-> POM
                                    + SUM(LimnoParam%wDMortFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Mortality-> POM
 	                                - SUM(LimnoParam%wDGrzFiAdPom(iLayer,iElem,:,om))                                                  & !Grazing Aduld Fish
                                    - SUM(LimnoParam%wDGrzFiJvPom(iLayer,iElem,:,om))                                                  & !Grazing Juveline Fish
                                    + SUM(LimnoParam%tDMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                            + wDEgesBirdPom_Layer                                                                               & !Egestao Bird -> POM                    
                                    - LimnoParam%wDMinAerW(iLayer,iElem,om)                                                             & ! Mineralization
                                    - LimnoParam%wDMinAnaerW(iLayer,iElem,om)                                                           & ! Mineralization
                                    + wDDifDom_Layer                                                                                    & ! Difusion
                                    + wDResusPom_bottom                                                                                 ! Resuspension

                    HydroParam%dVarEst(iLayer,iElem,2) = - wCSetPom_Layer                                                               & !Sedimentacao
		                            + LimnoParam%tCErosPom(iLayer,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))                & !Erosao Margem
                                    + SUM(LimnoParam%wCMortBac(iLayer,iElem,:))                                                         & !Lise Bacterioplankton
		                            + SUM(LimnoParam%wCLisPhyt(iLayer,iElem,:))                                                         & !Lise Fito
		                            + SUM(LimnoParam%wCMessZooOM(iLayer,iElem,:))                                                       & !Messy Feeding Zoops
		                            + SUM(LimnoParam%wCEgesZoo(iLayer,iElem,:))                                                         & !Pellets Zoops
                                    + SUM(LimnoParam%wCMortZoo(iLayer,iElem,:))                                                         & !Mortalidade Zoops
 	                                - SUM(LimnoParam%wCGrzZooPom(iLayer,iElem,:,om))                                                   & !Grazing Zoops
	                                + SUM(LimnoParam%wCMessFiAdOM(iLayer,iElem,:))                                                      & !Messy Feeding Adult Fish
                                    + SUM(LimnoParam%wCMessFiJvOM(iLayer,iElem,:))                                                      & !Messy Feeding Juveline Fish
	                                + SUM(LimnoParam%wCEgesFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Pellets-> POM
                                    + SUM(LimnoParam%wCEgesFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Pellets-> POM
	                                + SUM(LimnoParam%wCMortFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Mortality-> POM
                                    + SUM(LimnoParam%wCMortFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Mortality-> POM
 	                                - SUM(LimnoParam%wCGrzFiAdPom(iLayer,iElem,:,om))                                                  & !Grazing Aduld Fish
                                    - SUM(LimnoParam%wCGrzFiJvPom(iLayer,iElem,:,om))                                                  & !Grazing Juveline Fish
                                    + SUM(LimnoParam%tCMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                            + wCEgesBirdPom_Layer                                                                               & !Egestao Bird -> POM                    
                                    - LimnoParam%wCMinAerW(iLayer,om)                                                                   & ! Mineralization
                                    - LimnoParam%wCMinAnaerW(iLayer,om)                                                                 & ! Mineralization
                                    + wCDifDom_Layer                                                                                    & ! Difusion
                                    + wCResusPom_bottom                                                                                 ! Resuspension

                    HydroParam%dVarEst(iLayer,iElem,3) = - wNSetPom_Layer                                                               & !Sedimentacao
		                            + LimnoParam%tNErosPom(iLayer,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))                 & !Erosao Margem
                                    + SUM(LimnoParam%wNMortBac(iLayer,iElem,:))                                                         & !Lise Bacterioplankton
		                            + SUM(LimnoParam%wNLisPhyt(iLayer,iElem,:))                                                         & !Lise Fito
		                            + SUM(LimnoParam%wNMessZooOM(iLayer,iElem,:))                                                       & !Messy Feeding Zoops
		                            + SUM(LimnoParam%wNEgesZoo(iLayer,iElem,:))                                                         & !Pellets Zoops
                                    + SUM(LimnoParam%wNMortZoo(iLayer,iElem,:))                                                         & !Mortalidade Zoops
 	                                - SUM(LimnoParam%wNGrzZooPom(iLayer,iElem,:,om))                                                    & !Grazing Zoops
	                                + SUM(LimnoParam%wNMessFiAdOM(iLayer,iElem,:))                                                      & !Messy Feeding Adult Fish
                                    + SUM(LimnoParam%wNMessFiJvOM(iLayer,iElem,:))                                                      & !Messy Feeding Juveline Fish
	                                + SUM(LimnoParam%wNEgesFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Pellets-> POM
                                    + SUM(LimnoParam%wNEgesFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Pellets-> POM
	                                + SUM(LimnoParam%wNMortFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Mortality-> POM
                                    + SUM(LimnoParam%wNMortFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Mortality-> POM
 	                                - SUM(LimnoParam%wNGrzFiAdPom(iLayer,iElem,:,om))                                                   & !Grazing Aduld Fish
                                    - SUM(LimnoParam%wNGrzFiJvPom(iLayer,iElem,:,om))                                                   & !Grazing Juveline Fish
                                    + SUM(LimnoParam%tNMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                            + wNEgesBirdPom_Layer                                                                               & !Egestao Bird -> POM                    
                                    - LimnoParam%wNMinAerW(iLayer,om)                                                                   & ! Mineralization
                                    - LimnoParam%wNMinAnaerW(iLayer,om)                                                                 & ! Mineralization
                                    + wNDifDom_Layer                                                                                    & ! Difusion
                                    + wNResusPom_bottom                                                                                 ! Resuspension

                    HydroParam%dVarEst(iLayer,iElem,4) = - wPSetPom_Layer                                                               & !Sedimentacao
		                            + LimnoParam%tPErosPom(iLayer,om)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))                 & !Erosao Margem
                                    + SUM(LimnoParam%wPMortBac(iLayer,iElem,:))                                                         & !Lise Bacterioplankton
		                            + LimnoParam%fDLisDomPhyt*SUM(LimnoParam%wPLisPhyt(iLayer,iElem,:))                                                         & !Lise Fito
		                            + SUM(LimnoParam%wPMessZooOM(iLayer,iElem,:))                                                       & !Messy Feeding Zoops
		                            + SUM(LimnoParam%wPEgesZoo(iLayer,iElem,:))                                                         & !Pellets Zoops
                                    + SUM(LimnoParam%wPMortZoo(iLayer,iElem,:))                                                         & !Mortalidade Zoops
 	                                - SUM(LimnoParam%wPGrzZooPom(iLayer,iElem,:,om))                                                    & !Grazing Zoops
	                                + SUM(LimnoParam%wPMessFiAdOM(iLayer,iElem,:))                                                      & !Messy Feeding Adult Fish
                                    + SUM(LimnoParam%wPMessFiJvOM(iLayer,iElem,:))                                                      & !Messy Feeding Juveline Fish
	                                + SUM(LimnoParam%wPEgesFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Pellets-> POM
                                    + SUM(LimnoParam%wPEgesFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Pellets-> POM
	                                + SUM(LimnoParam%wPMortFiAd(iLayer,iElem,:))                                                        & !Aduld Fish Mortality-> POM
                                    + SUM(LimnoParam%wPMortFiJv(iLayer,iElem,:))                                                        & !Juveline Fish Mortality-> POM
 	                                - SUM(LimnoParam%wPGrzFiAdPom(iLayer,iElem,:,om))                                                   & !Grazing Aduld Fish
                                    - SUM(LimnoParam%wPGrzFiJvPom(iLayer,iElem,:,om))                                                   & !Grazing Juveline Fish
                                    + SUM(LimnoParam%tPMortMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))  & !Mortallidade Mac Agua -> POM
		                            + wPEgesBirdPom_Layer                                                                               & !Egestao Bird -> POM                    
                                    - LimnoParam%wPMinAerW(iLayer,om)                                                                   & ! Mineralization
                                    - LimnoParam%wPMinAnaerW(iLayer,om)                                                                 & ! Mineralization
                                    + wPDifDom_Layer                                                                                    & ! Difusion
                                    + wPResusPom_bottom                                                                                 ! Resuspension
                
                EndIf
            EndDo !Loop Layer
        EndDo !Loop Cell
        

        !-------------------------------------------------------------------------------------!
        !                       Solver Transport Equation for POM   	    	      !
        !-------------------------------------------------------------------------------------!
    
        If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDPomW(:,:,om),LimnoParam%sDPomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCPomW(:,:,om),LimnoParam%sCPomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNPomW(:,:,om),LimnoParam%sNPomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPPomW(:,:,om),LimnoParam%sPPomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDPomW(:,:,om),LimnoParam%sDPomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCPomW(:,:,om),LimnoParam%sCPomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNPomW(:,:,om),LimnoParam%sNPomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPPomW(:,:,om),LimnoParam%sPPomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
        ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDPomW(:,:,om),LimnoParam%sDPomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCPomW(:,:,om),LimnoParam%sCPomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNPomW(:,:,om),LimnoParam%sNPomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPPomW(:,:,om),LimnoParam%sPPomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDPomW(:,:,om),LimnoParam%sDPomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCPomW(:,:,om),LimnoParam%sCPomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNPomW(:,:,om),LimnoParam%sNPomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPPomW(:,:,om),LimnoParam%sPPomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
        EndiF     

    EndDo !Loop DOM
    
    Do om = 1, LimnoParam%npom ! =1 -> Labile; =2 -> Refractory
        
        !-------------------------------------------------------------------------------------!
        !                               Boundary Condition for DOM				      !
        !-------------------------------------------------------------------------------------!
        If (LimnoParam%InclMatOrgSplit == 1) Then
            If (LimnoParam%InclPomDomSplit == 1) Then
                If (om==LimnoParam%LABIL) Then
                    index = 0
                ElseIf(om==LimnoParam%REFRATARIA) Then
                    index = 5
                    HydroParam%uLoadVarEst = LimnoParam%uDLoadDom
                EndIf 
            Else 
                index = 5
                HydroParam%uLoadVarEst = LimnoParam%uDLoadDom
            EndIf
        Else
            index = 0
        EndIf       
        
        Do iElem = 1,MeshParam%nElem
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)

            !-------------------------------------------------------------------------------------!
            !                            DOM Processes in Water							  !
            !-------------------------------------------------------------------------------------!
            !-<DOM>----------------------------------------------!
            !(1) Fotolise - DOCR=>DOCL				             !
            !(1) Miner.(BAC=0)-DOC->CO2/NH4/PO4/Si               !
            !(3) Miner.(BAC=1)-DOC->CO2/NH4/PO4/Si*              !
            !(2) Hidrólise agua - POC->DOC*		                 !
            !(4) Mortalidade BAC 					             !
            !(5) Excrecao FITO          			             !
            !(6) Contribuicao da biota (mort e exc)	             !
            !(7) Erosao								             !
            !(8) Fluxo do sedimento					             !
            !----------------------------------------------------!      
                
            !-------------------------------------------------------------------------------------!
            !                               Source/Sink Term for DOM						      !
            !-------------------------------------------------------------------------------------!
                If (LimnoParam%InclMatOrgSplit == 1) Then ! Case 1: Splited Organic Matter (POM and DOM)
                    
                    If (LimnoParam%InclPomDomSplit == 1) Then ! Case 2: Splited Organic Matter (POM and DOM) and Splited fractions (Labile and Refratory)
                        If (om==LimnoParam%LABIL) Then
		                    HydroParam%dVarEst(iLayer,iElem,1) = - LimnoParam%wDMinAerW(iLayer,iElem,1)                                 & ! Mineralization (DOM Labile -> Inorganic Components)
                                        - LimnoParam%wDMinAnaerW(iLayer,iElem,1)                                                        & ! Mineralization (DOM Labile -> Inorganic Components)
                                        + LimnoParam%wDMinAerW(iLayer,iElem,2)                                                          & !Mineralization (DOM Refratory -> DOM Labile)
                                        + LimnoParam%wDMinAnaerW(iLayer,iElem,2)                                                        & !Mineralization (DOM Refratory -> DOM Labile)
		                                + LimnoParam%wDHidPomW(iLayer,iElem,om)                                                         & !Hidrolise POML
                                        + wDDifDom_Layer                                                                                ! Difusion

		                    HydroParam%dVarEst(iLayer,iElem,2) = - LimnoParam%wCMinAerW(iLayer,1)                                       & ! Mineralization (DOM Labile -> Inorganic Components)
                                        - LimnoParam%wCMinAnaerW(iLayer,1)                                                              & ! Mineralization (DOM Labile -> Inorganic Components)
                                        + LimnoParam%wCMinAerW(iLayer,2)                                                                & !Mineralization (DOM Refratory -> DOM Labile)
                                        + LimnoParam%wCMinAnaerW(iLayer,2)                                                              & !Mineralization (DOM Refratory -> DOM Labile)
		                                + LimnoParam%wCHidPomW(iLayer,om)                                                               & !Hidrolise POML
                                        + wCDifDom_Layer                                                                                ! Difusion

		                    HydroParam%dVarEst(iLayer,iElem,3) = - LimnoParam%wNMinAerW(iLayer,1)                                       & ! Mineralization (DOM Labile -> Inorganic Components)
                                        - LimnoParam%wNMinAnaerW(iLayer,1)                                                              & ! Mineralization (DOM Labile -> Inorganic Components)
                                        + LimnoParam%wNMinAerW(iLayer,2)                                                                & !Mineralization (DOM Refratory -> DOM Labile)
                                        + LimnoParam%wNMinAnaerW(iLayer,2)                                                              & !Mineralization (DOM Refratory -> DOM Labile)
		                                + LimnoParam%wNHidPomW(iLayer,om)                                                               & !Hidrolise POML
                                        + wNDifDom_Layer                                                                                ! Difusion

		                    HydroParam%dVarEst(iLayer,iElem,4) = - LimnoParam%wPMinAerW(iLayer,1)                                       & ! Mineralization (DOM Labile -> Inorganic Components)
                                        - LimnoParam%wPMinAnaerW(iLayer,1)                                                              & ! Mineralization (DOM Labile -> Inorganic Components)
                                        + LimnoParam%wPMinAerW(iLayer,2)                                                                & !Mineralization (DOM Refratory -> DOM Labile)
                                        + LimnoParam%wPMinAnaerW(iLayer,2)                                                              & !Mineralization (DOM Refratory -> DOM Labile)
		                                + LimnoParam%wPHidPomW(iLayer,om)                                                               & !Hidrolise POML
                                        + wPDifDom_Layer                                                                                ! Difusion
                            
                        ElseIf(om==LimnoParam%REFRATARIA) Then       
		                    HydroParam%dVarEst(iLayer,iElem,1) = LimnoParam%wDHidPomW(iLayer,iElem,om)                                  & !Hidrolise POML
                                        - LimnoParam%wDMinAerW(iLayer,iElem,2)                                                          & !Mineralization (DOM Refratory -> DOM Labile)
                                        - LimnoParam%wDMinAnaerW(iLayer,iElem,2)                                                        & !Mineralization (DOM Refratory -> DOM Labile)
                                        + SUM(LimnoParam%wDMortBac(iLayer,iElem,:))                                                     !Lise Bacterioplankton

		                    HydroParam%dVarEst(iLayer,iElem,2) = LimnoParam%wCHidPomW(iLayer,om)                                        & !Hidrolise POML
                                        - LimnoParam%wCMinAerW(iLayer,2)                                                                & !Mineralization (DOM Refratory -> DOM Labile)
                                        - LimnoParam%wCMinAnaerW(iLayer,2)                                                              & !Mineralization (DOM Refratory -> DOM Labile)
                                        + SUM(LimnoParam%wCMortBac(iLayer,iElem,:))                                                     !Lise Bacterioplankton

		                    HydroParam%dVarEst(iLayer,iElem,3) = LimnoParam%wNHidPomW(iLayer,om)                                        & !Hidrolise POML
                                        - LimnoParam%wNMinAerW(iLayer,2)                                                                & !Mineralization (DOM Refratory -> DOM Labile)
                                        - LimnoParam%wNMinAnaerW(iLayer,2)                                                              & !Mineralization (DOM Refratory -> DOM Labile)
                                        + SUM(LimnoParam%wNMortBac(iLayer,iElem,:))                                                     !Lise Bacterioplankton

		                    HydroParam%dVarEst(iLayer,iElem,4) = LimnoParam%wPHidPomW(iLayer,om)                                        & !Hidrolise POML
                                        - LimnoParam%wPMinAerW(iLayer,2)                                                                & !Mineralization (DOM Refratory -> DOM Labile)
                                        - LimnoParam%wPMinAnaerW(iLayer,2)                                                              & !Mineralization (DOM Refratory -> DOM Labile)
                                        + SUM(LimnoParam%wPMortBac(iLayer,iElem,:))                                                     !Lise Bacterioplankton
                            
                        EndIf
                        
                    Else ! Case 3: Splited Organic Matter (POM and DOM) and non-splited fractions    
		                HydroParam%dVarEst(iLayer,iElem,1) = - LimnoParam%wDMinAerW(iLayer,iElem,om)                                    & ! Mineralization (DOM Labile -> Inorganic Components)
                                        - LimnoParam%wDMinAnaerW(iLayer,iElem,om)                                                       & ! Mineralization (DOM Labile -> Inorganic Components)
		                                + LimnoParam%wDHidPomW(iLayer,iElem,om)                                                         & ! Hidrólise - POM => DOM
                                        + SUM(LimnoParam%wDMortBac(iLayer,iElem,:))                                                     & !Lise Bacterioplankton
                                        + wDDifDom_Layer                                                                                ! Difusion

		                HydroParam%dVarEst(iLayer,iElem,2) = - LimnoParam%wCMinAerW(iLayer,om)                                          & ! Mineralization (DOM Labile -> Inorganic Components)
                                        - LimnoParam%wCMinAnaerW(iLayer,om)                                                             & ! Mineralization (DOM Labile -> Inorganic Components)
		                                + LimnoParam%wCHidPomW(iLayer,om)                                                               & ! Hidrólise - POM => DOM
                                        + SUM(LimnoParam%wCMortBac(iLayer,iElem,:))                                                     & !Lise Bacterioplankton
                                        + wCDifDom_Layer                                                                                ! Difusion

		                HydroParam%dVarEst(iLayer,iElem,3) = - LimnoParam%wNMinAerW(iLayer,om)                                          & ! Mineralization (DOM Labile -> Inorganic Components)
                                        - LimnoParam%wNMinAnaerW(iLayer,om)                                                             & ! Mineralization (DOM Labile -> Inorganic Components)
		                                + LimnoParam%wNHidPomW(iLayer,om)                                                               & ! Hidrólise - POM => DOM
                                        + SUM(LimnoParam%wNMortBac(iLayer,iElem,:))                                                     & !Lise Bacterioplankton
                                        + wNDifDom_Layer                                                                                ! Difusion

		                HydroParam%dVarEst(iLayer,iElem,4) = - LimnoParam%wPMinAerW(iLayer,om)                                          & ! Mineralization (DOM Labile -> Inorganic Components)
                                        - LimnoParam%wPMinAnaerW(iLayer,om)                                                             & ! Mineralization (DOM Labile -> Inorganic Components)
		                                + LimnoParam%wPHidPomW(iLayer,om)                                                               & ! Hidrólise - POM => DOM
                                        + SUM(LimnoParam%wPMortBac(iLayer,iElem,:))                                                     & !Lise Bacterioplankton
                                        + wPDifDom_Layer                                                                                ! Difusion
                        
                    EndIf
                    
                        
                Else ! Case 4: Non-splited Organic Matter and non-splited fractions  (only one box of organic matter -  POM=Detritus - PCLake)
                        
                EndIf
            EndDo !Loop Layer
        EndDo !Loop Cell
        
      
        !-------------------------------------------------------------------------------------!
        !                       Solver Transport Equation for DOM   	    	      !
        !-------------------------------------------------------------------------------------!
        If (LimnoParam%InclMatOrgSplit == 1) Then
            If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
                Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDDomW(:,:,om),LimnoParam%sDDomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
                Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCDomW(:,:,om),LimnoParam%sCDomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
                Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNDomW(:,:,om),LimnoParam%sNDomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
                Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPDomW(:,:,om),LimnoParam%sPDomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
            ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
                Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDDomW(:,:,om),LimnoParam%sDDomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
                Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCDomW(:,:,om),LimnoParam%sCDomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
                Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNDomW(:,:,om),LimnoParam%sNDomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
                Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPDomW(:,:,om),LimnoParam%sPDomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
                Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDDomW(:,:,om),LimnoParam%sDDomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
                Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCDomW(:,:,om),LimnoParam%sCDomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
                Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNDomW(:,:,om),LimnoParam%sNDomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
                Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPDomW(:,:,om),LimnoParam%sPDomWP(:,:,om),dt,dtday,HydroParam,MeshParam)
            ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
                Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDDomW(:,:,om),LimnoParam%sDDomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
                Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCDomW(:,:,om),LimnoParam%sCDomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
                Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNDomW(:,:,om),LimnoParam%sNDomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
                Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPDomW(:,:,om),LimnoParam%sPDomWP(:,:,om),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            EndIf
        Else
            LimnoParam%sDDomW(:,:,om) = LimnoParam%sDPomW(:,:,om)
            LimnoParam%sCDomW(:,:,om) = LimnoParam%sCPomW(:,:,om)
            LimnoParam%sNDomW(:,:,om) = LimnoParam%sNPomW(:,:,om)
            LimnoParam%sPDomW(:,:,om) = LimnoParam%sPPomW(:,:,om)
        EndIf
        
    EndDo !Loop Labile and Refratory
    
Return
End
