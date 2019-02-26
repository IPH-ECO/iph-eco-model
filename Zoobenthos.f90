Subroutine Zoobenthos(HydroParam,MeshParam,LimnoParam,dt,dtday)
    
    Use Hydrodynamic
    Use MeshVars
    Use LimnologyVars
    Use uTempFunction

    Implicit None
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam
    type(LimnologyParam) :: LimnoParam
    Integer:: iElem,ben,gg,pom,bb
    Real:: dt,dtday,aDOMW
    Real:: NearZero = 1e-10
    
    
    !Zoobenthos processes in the sediment
    
    
    
    !-------------------------------------------------------------------------------------!
    !                               Zoobenthos Processes in Water					      !
    !-------------------------------------------------------------------------------------!
    !-Zoobenthos--------------------------------------------------------------------------!
    ! 1. Grazing                                                                          !
    ! 2. Respiration and Excretion                                                        !
    ! 3. Mortality                                                                        !
    ! 4. Messy                                                                            !
    ! 5. Egestion                                                                         !
    !-------------------------------------------------------------------------------------!
    Do ben = 1, LimnoParam%nben 
        Do iElem = 1,MeshParam%nElem
            !-----------------------------------------------------------------------
            !  temperature function
            !-----------------------------------------------------------------------
            !  temp._function_of_zoobenthos
            LimnoParam%uFunTmBent = uFunTmBio(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem),iElem),LimnoParam%cSigTmBent(ben),LimnoParam%cTmOptBent(ben))
        
            !---------------------------------------------------------------------------
            !  zoobenthos assimilation,DW
            !---------------------------------------------------------------------------
            !  total_food_for_zoobenthos
            LimnoParam%oDFoodBent = 0.
            LimnoParam%oCFoodBent = 0.
            LimnoParam%oNFoodBent = 0.
            LimnoParam%oPFoodBent = 0.
            !  total_phytoplankton_for_zoobenthos
            LimnoParam%oDFoodBentPhyt = 0.
            !  total_bacterioplankton_for_zoobenthos
            LimnoParam%oDFoodBentBac = 0.
            !  total_POM_for_zoobenthos
            LimnoParam%oDFoodBentPom = 0.
            !  organic_seston
            aDOMW = 0.
            ! Phyto_for_zoobenthos
            If (LimnoParam%aben(ben)==1) Then !Micro-zoobenthos
                Do gg = 1, LimnoParam%nphy
                    If (LimnoParam%aphy(gg)==LimnoParam%DINOF.or.LimnoParam%aphy(gg)==LimnoParam%MDIAT.or.LimnoParam%aphy(gg)==LimnoParam%FDIAT) Then
                        LimnoParam%oDFoodBent = LimnoParam%oDFoodBent + LimnoParam%cPrefBentPhyt(ben,gg)*LimnoParam%sDPhytS(iElem,gg)
                        LimnoParam%oCFoodBent = LimnoParam%oCFoodBent + LimnoParam%cPrefBentPhyt(ben,gg)*LimnoParam%sCPhytS(iElem,gg)
                        LimnoParam%oNFoodBent = LimnoParam%oNFoodBent + LimnoParam%cPrefBentPhyt(ben,gg)*LimnoParam%sNPhytS(iElem,gg)
                        LimnoParam%oPFoodBent = LimnoParam%oPFoodBent + LimnoParam%cPrefBentPhyt(ben,gg)*LimnoParam%sPPhytS(iElem,gg)
                        LimnoParam%oDFoodBentPhyt = LimnoParam%oDFoodBentPhyt + LimnoParam%cPrefBentPhyt(ben,gg)*LimnoParam%sDPhytS(iElem,gg)
                        aDOMW = aDOMW + LimnoParam%sDPhytS(iElem,gg)
                    EndIf
                EndDo
            Else ! Meso and Macro-zoobenthos
                Do gg = 1, LimnoParam%nphy
                    If (LimnoParam%aphy(gg)==LimnoParam%DINOF.or.LimnoParam%aphy(gg)==LimnoParam%CYANO.or.LimnoParam%aphy(gg)==LimnoParam%NODUL.or.LimnoParam%aphy(gg)==LimnoParam%CHLOR.or.LimnoParam%aphy(gg)==LimnoParam%CRYPT.or.LimnoParam%aphy(gg)==LimnoParam%MDIAT.or.LimnoParam%aphy(gg)==LimnoParam%FDIAT) Then
                        LimnoParam%oDFoodBent = LimnoParam%oDFoodBent + LimnoParam%cPrefBentPhyt(ben,gg)*LimnoParam%sDPhytS(iElem,gg)
                        LimnoParam%oCFoodBent = LimnoParam%oCFoodBent + LimnoParam%cPrefBentPhyt(ben,gg)*LimnoParam%sCPhytS(iElem,gg)
                        LimnoParam%oNFoodBent = LimnoParam%oNFoodBent + LimnoParam%cPrefBentPhyt(ben,gg)*LimnoParam%sNPhytS(iElem,gg)
                        LimnoParam%oPFoodBent = LimnoParam%oPFoodBent + LimnoParam%cPrefBentPhyt(ben,gg)*LimnoParam%sPPhytS(iElem,gg)
                        LimnoParam%oDFoodBentPhyt = LimnoParam%oDFoodBentPhyt + LimnoParam%cPrefBentPhyt(ben,gg)*LimnoParam%sDPhytS(iElem,gg)
                        aDOMW = aDOMW + LimnoParam%sDPhytS(iElem,gg)
                    EndIf
                EndDo
            EndIf
            
            ! Bacterioplankton_for_zoobenthos
            Do bb = 1, LimnoParam%nbac
                If (LimnoParam%abac(bb)==LimnoParam%AEROBICA.or.LimnoParam%abac(bb)==LimnoParam%ANAEROBICA.or.LimnoParam%abac(bb)==LimnoParam%BACTERIA) Then
                    LimnoParam%oDFoodBent = LimnoParam%oDFoodBent + LimnoParam%cPrefBentBac(ben,bb)*LimnoParam%sDBacS(iElem,bb)
                    LimnoParam%oCFoodBent = LimnoParam%oCFoodBent + LimnoParam%cPrefBentBac(ben,bb)*LimnoParam%sCBacS(iElem,bb)
                    LimnoParam%oNFoodBent = LimnoParam%oNFoodBent + LimnoParam%cPrefBentBac(ben,bb)*LimnoParam%sNBacS(iElem,bb)
                    LimnoParam%oPFoodBent = LimnoParam%oPFoodBent + LimnoParam%cPrefBentBac(ben,bb)*LimnoParam%sPBacS(iElem,bb)
                    LimnoParam%oDFoodBentBac = LimnoParam%oDFoodBentBac + LimnoParam%cPrefBentBac(ben,bb)*LimnoParam%sDBacS(iElem,bb)
                    aDOMW = aDOMW + LimnoParam%sDBacS(iElem,bb)
                EndIf
            EndDo
            ! POM_for_zoobenthos
            Do pom = 1, LimnoParam%npom
                If (pom==LimnoParam%LABIL.or.pom==LimnoParam%REFRATARIA) Then
                    LimnoParam%oDFoodBent = LimnoParam%oDFoodBent + LimnoParam%cPrefBentPom(ben,pom)*LimnoParam%sDPomS(iElem,pom)
                    LimnoParam%oCFoodBent = LimnoParam%oCFoodBent + LimnoParam%cPrefBentPom(ben,pom)*LimnoParam%sCPomS(iElem,pom)
                    LimnoParam%oNFoodBent = LimnoParam%oNFoodBent + LimnoParam%cPrefBentPom(ben,pom)*LimnoParam%sNPomS(iElem,pom)
                    LimnoParam%oPFoodBent = LimnoParam%oPFoodBent + LimnoParam%cPrefBentPom(ben,pom)*LimnoParam%sPPomS(iElem,pom)
                    LimnoParam%oDFoodBentPom = LimnoParam%oDFoodBentPom + LimnoParam%cPrefBentPom(ben,pom)*LimnoParam%sDPomS(iElem,pom)
                    aDOMW = aDOMW + LimnoParam%sDPomS(iElem,pom)
                EndIf
            EndDo            
            
            
            !  food_limitation_function_of_zoobenthos
            LimnoParam%aDSatBent = LimnoParam%oDFoodBent /(LimnoParam%hDFoodBent(ben) + aDOMW)
            !  intrinsic_net_increase_rate_of_zoobenthos
            LimnoParam%ukDIncrBent = (LimnoParam%kDAssBent(ben) - LimnoParam%kDRespBent(ben)) * LimnoParam%uFunTmBent - LimnoParam%kMortBent(ben)
            !  environmental_correction_of_zoobenthos
            LimnoParam%tDEnvBent = max(0.0,LimnoParam%ukDIncrBent / LimnoParam%cDCarrBent(ben) * LimnoParam%sDBent(iElem,ben)*LimnoParam%sDBent(iElem,ben))
            !  assimilation_of_zoobenthos
            LimnoParam%tDAssBent = LimnoParam%aDSatBent *(LimnoParam%kDAssBent(ben) * LimnoParam%uFunTmBent * LimnoParam%sDBent(iElem,ben) - LimnoParam%tDEnvBent)
            !  consumption_of_zoobenthos
            LimnoParam%tDConsBent = LimnoParam%tDAssBent / LimnoParam%fDAssBent(ben)
            
            !  DW_phytoplankton_consumption_by_zoobenthos
            Do gg = 1, LimnoParam%nphy
                LimnoParam%tDGrzBentPhyt(iElem,ben,gg) = LimnoParam%cPrefBentPhyt(ben,gg)*LimnoParam%sDPhytS(iElem,gg)/LimnoParam%oDFoodBent * LimnoParam%tDConsBent
                LimnoParam%tCGrzBentPhyt(iElem,ben,gg) = LimnoParam%rCDPhytS(iElem,gg)*LimnoParam%tDGrzBentPhyt(iElem,ben,gg)
                LimnoParam%tNGrzBentPhyt(iElem,ben,gg) = LimnoParam%rNDPhytS(iElem,gg)*LimnoParam%tDGrzBentPhyt(iElem,ben,gg)
                LimnoParam%tPGrzBentPhyt(iElem,ben,gg) = LimnoParam%rPDPhytS(iElem,gg)*LimnoParam%tDGrzBentPhyt(iElem,ben,gg)
            EndDo  
            !  DW_bacterioplankton_consumption_by_zoobenthos
            Do bb = 1, LimnoParam%nbac
                LimnoParam%tDGrzBentBac(iElem,ben,bb) = LimnoParam%cPrefBentBac(ben,bb)*LimnoParam%sDBacS(iElem,bb)/LimnoParam%oDFoodBent * LimnoParam%tDConsBent
                LimnoParam%tCGrzBentBac(iElem,ben,bb) = LimnoParam%rCDBacS(iElem,bb)*LimnoParam%tDGrzBentBac(iElem,ben,bb)
                LimnoParam%tNGrzBentBac(iElem,ben,bb) = LimnoParam%rNDBacS(iElem,bb)*LimnoParam%tDGrzBentBac(iElem,ben,bb)
                LimnoParam%tPGrzBentBac(iElem,ben,bb) = LimnoParam%rPDBacS(iElem,bb)*LimnoParam%tDGrzBentBac(iElem,ben,bb)
            EndDo
            !  DW_detritus_consumption_by_zoobenthos
            Do pom = 1, LimnoParam%npom
                LimnoParam%tDGrzBentPom(iElem,ben,pom) = LimnoParam%cPrefBentPom(ben,pom)*LimnoParam%sDPomS(iElem,pom)/LimnoParam%oDFoodBent * LimnoParam%tDConsBent
                LimnoParam%tCGrzBentPom(iElem,ben,pom) = LimnoParam%rCDPomS(iElem,pom)*LimnoParam%tDGrzBentPom(iElem,ben,pom)
                LimnoParam%tNGrzBentPom(iElem,ben,pom) = LimnoParam%rNDPomS(iElem,pom)*LimnoParam%tDGrzBentPom(iElem,ben,pom)
                LimnoParam%tPGrzBentPom(iElem,ben,pom) = LimnoParam%rPDPomS(iElem,pom)*LimnoParam%tDGrzBentPom(iElem,ben,pom)
            EndDo
            !-----------------------------------------------------------------------
            !  zoobenthos assimilation C
            !-----------------------------------------------------------------------
            !  C/DW_ratio_of_zoobenthos_food
            LimnoParam%rCDFoodBent = LimnoParam%oCFoodBent /(LimnoParam%oDFoodBent+NearZero)
            !  total_C_consumption
            LimnoParam%tCConsBent = (sum(LimnoParam%tCGrzBentPhyt(iElem,ben,:)) + sum(LimnoParam%tCGrzBentBac(iElem,ben,:)) + sum(LimnoParam%tCGrzBentPom(iElem,ben,:)))
            !  C_assimilation_efficiency_of_zoobenthos
            LimnoParam%afCAssBent = min(1.0,LimnoParam%cCDBentRef(ben) / LimnoParam%rCDFoodBent * LimnoParam%fDAssBent(ben))
            !  assimilation_by_zoobenthos
            LimnoParam%tCAssBent = LimnoParam%afCAssBent*LimnoParam%tCConsBent
            !-----------------------------------------------------------------------
            !  zoobenthos assimilation N
            !-----------------------------------------------------------------------
            !  N/DW_ratio_of_zoobenthos_food
            LimnoParam%rNDFoodBent = LimnoParam%oNFoodBent /(LimnoParam%oDFoodBent+NearZero)
            !  total_N_consumption
            LimnoParam%tNConsBent = (sum(LimnoParam%tNGrzBentPhyt(iElem,ben,:)) + sum(LimnoParam%tNGrzBentBac(iElem,ben,:)) + sum(LimnoParam%tNGrzBentPom(iElem,ben,:)))
            !  N_assimilation_efficiency_of_zoobenthos
            LimnoParam%afNAssBent = min(1.0,LimnoParam%cNDBentRef(ben) / LimnoParam%rNDFoodBent * LimnoParam%fDAssBent(ben))
            !  assimilation_by_zoobenthos
            LimnoParam%tNAssBent = LimnoParam%afNAssBent*LimnoParam%tNConsBent
            !-----------------------------------------------------------------------
            !  zoobenthos assimilation P
            !-----------------------------------------------------------------------
            !  P/DW_ratio_of_zoobenthos_food
            LimnoParam%rPDFoodBent = LimnoParam%oPFoodBent /(LimnoParam%oDFoodBent+NearZero)
            !  total_P_consumption
            LimnoParam%tPConsBent = (sum(LimnoParam%tPGrzBentPhyt(iElem,ben,:)) + sum(LimnoParam%tPGrzBentBac(iElem,ben,:)) + sum(LimnoParam%tPGrzBentPom(iElem,ben,:)))
            !  P_assimilation_efficiency_of_zoobenthos
            LimnoParam%afPAssBent = min(1.0,LimnoParam%cPDBentRef(ben) / LimnoParam%rPDFoodBent * LimnoParam%fDAssBent(ben))
            !  assimilation_by_zoobenthos
            LimnoParam%tPAssBent = LimnoParam%afPAssBent*LimnoParam%tPConsBent
            
            !-----------------------------------------------------------------------
            !  zoobenthos respiration and excretion, DW,P and N
            !-----------------------------------------------------------------------
            !  corr._factor_of_zoobenthos_respiration_for_P_and_N_content
            LimnoParam%aCorDRespBent = max(LimnoParam%cPDBentRef(ben) / LimnoParam%rPDBent(iElem,ben),LimnoParam%cNDBentRef(ben) / LimnoParam%rNDBent(iElem,ben))
            !  respiration_of_zoobenthos
            LimnoParam%tDRespBent(iElem,ben) = LimnoParam%aCorDRespBent * LimnoParam%kDRespBent(ben) * LimnoParam%uFunTmBent * LimnoParam%sDBent(iElem,ben)
            !  C_excretion
            LimnoParam%tCRespBent(iElem,ben) = LimnoParam%rCDBent(iElem,ben) / LimnoParam%cCDBentRef(ben) * LimnoParam%kDRespBent(ben) * LimnoParam%uFunTmBent * LimnoParam%sCBent(iElem,ben)
            !  N_excretion
            LimnoParam%tNExcBent(iElem,ben) = LimnoParam%rNDBent(iElem,ben) / LimnoParam%cNDBentRef(ben) * LimnoParam%kDRespBent(ben) * LimnoParam%uFunTmBent * LimnoParam%sNBent(iElem,ben)
            !  P_excretion
            LimnoParam%tPExcBent(iElem,ben) = LimnoParam%rPDBent(iElem,ben) / LimnoParam%cPDBentRef(ben) * LimnoParam%kDRespBent(ben) * LimnoParam%uFunTmBent * LimnoParam%sPBent(iElem,ben)
            
            !-----------------------------------------------------------------------
            !  zoobenthos mortality
            !-----------------------------------------------------------------------
            !  zoobenthos_mortality_incl._environmental_correction
            LimnoParam%tDMortBent(iElem,ben) = LimnoParam%kMortBent(ben) * LimnoParam%sDBent(iElem,ben) +(1.0 - LimnoParam%aDSatBent) * LimnoParam%tDEnvBent
            !  zoopl._mortality_C
            LimnoParam%tCMortBent(iElem,ben) = LimnoParam%rCDBent(iElem,ben) * LimnoParam%tDMortBent(iElem,ben)
            !  zoopl._mortality_N
            LimnoParam%tNMortBent(iElem,ben) = LimnoParam%rNDBent(iElem,ben) * LimnoParam%tDMortBent(iElem,ben)
            !  zoopl._mortality_P
            LimnoParam%tPMortBent(iElem,ben) = LimnoParam%rPDBent(iElem,ben) * LimnoParam%tDMortBent(iElem,ben)
            
            
            !-----------------------------------------------------------------------
            !  zoobenthos Messy -> Assimilation => POM
            !-----------------------------------------------------------------------
            !  Messy_of_zoobenthos_DW
            LimnoParam%tDMessBentOM(iElem,ben) = LimnoParam%tDConsBent - LimnoParam%tDAssBent
            !  Messy_of_zoobenthos_C
            LimnoParam%tDMessBentOM(iElem,ben) = LimnoParam%tCConsBent - LimnoParam%tCAssBent
            !  Messy_of_zoobenthos_N
            LimnoParam%tDMessBentOM(iElem,ben) = LimnoParam%tNConsBent - LimnoParam%tNAssBent
            !  Messy_of_zoobenthos_P
            LimnoParam%tDMessBentOM(iElem,ben) = LimnoParam%tPConsBent - LimnoParam%tPAssBent
            
            !-----------------------------------------------------------------------
            !  zoobenthos Egestion/Pellets (-) -> POM
            !-----------------------------------------------------------------------
            !  egestion_of_zoobenthos_DW
            LimnoParam%tDEgesBent(iElem,ben) = LimnoParam%kPelBent(ben) * LimnoParam%sDBent(iElem,ben)
            !  egestion_of_zoobenthos_C
            LimnoParam%tCEgesBent(iElem,ben) = LimnoParam%rCDBent(iElem,ben) * LimnoParam%tDEgesBent(iElem,ben)
            !  egestion_of_zoobenthos_N
            LimnoParam%tNEgesBent(iElem,ben) = LimnoParam%rNDBent(iElem,ben) * LimnoParam%tDEgesBent(iElem,ben)
            !  egestion_of_zoobenthos_P
            LimnoParam%tPEgesBent(iElem,ben) = LimnoParam%rPDBent(iElem,ben) * LimnoParam%tDEgesBent(iElem,ben)
            
            !-------------------------------------------------------------------------------------!
            !                               Source/Sink Term for Zoobenthos                       !
            !-------------------------------------------------------------------------------------!
            
			LimnoParam%dVarEstS(iElem,1) =   LimnoParam%tDAssBent                                                               & !Assimilation
                                            -LimnoParam%tDRespBent(iElem,ben)                                                   & !Respiration
                                            -LimnoParam%tDEgesBent(iElem,ben)                                                   & !Egestion/Pellets
                                            -LimnoParam%tDMortBent(iElem,ben)                                                   & !Mortality
                                            - SUM(LimnoParam%wDGrzFiAdBent(HydroParam%ElSmallm(iElem),iElem,:,ben))             & !Predation by Adult Fish
                                            - SUM(LimnoParam%wDGrzFiJvBent(HydroParam%ElSmallm(iElem),iElem,:,ben))             !Predation by Juvenile Fish
            
			LimnoParam%dVarEstS(iElem,2) =   LimnoParam%tCAssBent                                                               & !Assimilation
                                            -LimnoParam%tCRespBent(iElem,ben)                                                   & !Respiration
                                            -LimnoParam%tCEgesBent(iElem,ben)                                                   & !Egestion/Pellets
                                            -LimnoParam%tCMortBent(iElem,ben)                                                   & !Mortality
                                            - SUM(LimnoParam%wCGrzFiAdBent(HydroParam%ElSmallm(iElem),iElem,:,ben))             & !Predation by Adult Fish
                                            - SUM(LimnoParam%wCGrzFiJvBent(HydroParam%ElSmallm(iElem),iElem,:,ben))             !Predation by Juvenile Fish

			LimnoParam%dVarEstS(iElem,3) =   LimnoParam%tNAssBent                                                               & !Assimilation
                                            -LimnoParam%tNExcBent(iElem,ben)                                                    & !Excretion
                                            -LimnoParam%tNEgesBent(iElem,ben)                                                   & !Egestion/Pellets
                                            -LimnoParam%tNMortBent(iElem,ben)                                                   & !Mortality
                                            - SUM(LimnoParam%wNGrzFiAdBent(HydroParam%ElSmallm(iElem),iElem,:,ben))             & !Predation by Adult Fish
                                            - SUM(LimnoParam%wNGrzFiJvBent(HydroParam%ElSmallm(iElem),iElem,:,ben))             !Predation by Juvenile Fish
			
            LimnoParam%dVarEstS(iElem,4) =   LimnoParam%tPAssBent                                                               & !Assimilation
                                            -LimnoParam%tPExcBent(iElem,ben)                                                    & !Excretion
                                            -LimnoParam%tPEgesBent(iElem,ben)                                                   & !Egestion/Pellets
                                            -LimnoParam%tPMortBent(iElem,ben)                                                   & !Mortality
                                            - SUM(LimnoParam%wNGrzFiAdBent(HydroParam%ElSmallm(iElem),iElem,:,ben))             & !Predation by Adult Fish
                                            - SUM(LimnoParam%wNGrzFiJvBent(HydroParam%ElSmallm(iElem),iElem,:,ben))             !Predation by Juvenile Fish
            
            LimnoParam%sDBent(iElem,ben) = max(NearZero,(dtday*LimnoParam%dVarEstS(iElem,1)+LimnoParam%sDBentP(iElem,ben)))
            LimnoParam%sCBent(iElem,ben) = max(NearZero,(dtday*LimnoParam%dVarEstS(iElem,2)+LimnoParam%sCBentP(iElem,ben)))
            LimnoParam%sNBent(iElem,ben) = max(NearZero,(dtday*LimnoParam%dVarEstS(iElem,3)+LimnoParam%sNBentP(iElem,ben)))
            LimnoParam%sPBent(iElem,ben) = max(NearZero,(dtday*LimnoParam%dVarEstS(iElem,4)+LimnoParam%sPBentP(iElem,ben)))
            
        EndDo !Loop Cell
    EndDo !Loop Functional Group
      

Return
End
