Subroutine Zooplankton(HydroParam,MeshParam,LimnoParam,dt,dtday)

    ! Zooplankton Modelling in the Water Column
    
    ! List of Modifications:
    !   24.07.2015: Routine Implementation      (Carlos Ruberto)
    !   25.07.2015: Routine Validation      (Carlos Ruberto)
    ! Programmer: Carlos Ruberto
    
    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars
    Use uTempFunction
    
    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    Integer:: index,iElem,iLayer,pom,zz,bb,gg,zg
    Real:: aDOMW
    Real:: V
    Real:: dt,dtday
    Real:: NearZero = 1e-10
    
    Do zz = 1, LimnoParam%nzoo
        Do iElem = 1,MeshParam%nElem
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                !-------------------------------------------------------------------------------------!
                !                               Zooplankton Processes in Water					      !
                !-------------------------------------------------------------------------------------!
                !-Zoops-------------------------------------------------------------------------------!
                !1. Grazing (+)                                                                       !
                !2. Messy Fedding (Egestion) (-) -> POM                                               !
                !3. Respiration/Excretion(-) -> CO2, N(DOM), P(DOM)                                   !
                !4. Egestion/Pellets(-) -> POM                                                        !
                !5. Mortality (-) -> POM                                                              !
                !6. Predation (-)                                                                     !
                !-------------------------------------------------------------------------------------!			
                !-------------------------------------------------------------------------------------!
                !-------------------------------------------------------------------------!
			    !1) Grazing                  
                !  temp._function_of_zooplankton
                LimnoParam%uFunTmZoo = uFunTmBio(LimnoParam%sDTempW(iLayer,iElem),LimnoParam%cSigTmZoo(zz),LimnoParam%cTmOptZoo(zz))
                !-----------------------------------------------------------------------
                !  zooplankton assimilation DW
                !-----------------------------------------------------------------------
                If (LimnoParam%azoo(zz)==1) Then !Micro-zooplankton
                    !  total_food_for_zooplankton
                    LimnoParam%oDFoodZoo = 0.
                    LimnoParam%oCFoodZoo = 0.
                    LimnoParam%oNFoodZoo = 0.
                    LimnoParam%oPFoodZoo = 0.
                    !  total_phytoplankton_for_zooplankton
                    LimnoParam%oDFoodZooPhyt = 0.
                    !  total_bacterioplankton_for_zooplankton
                    LimnoParam%oDFoodZooBac = 0.
                    !  total_POM_for_zooplankton
                    LimnoParam%oDFoodZooPom = 0.
                    !  organic_seston
                    aDOMW = 0.
                    ! Phyto_for_zooplankton
                    Do gg = 1, LimnoParam%nphy
                        If (LimnoParam%aphy(gg)==LimnoParam%DINOF.or.LimnoParam%aphy(gg)==LimnoParam%MDIAT.or.LimnoParam%aphy(gg)==LimnoParam%FDIAT) Then
                            LimnoParam%oDFoodZoo = LimnoParam%oDFoodZoo + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            LimnoParam%oCFoodZoo = LimnoParam%oCFoodZoo + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sCPhytW(iLayer,iElem,gg)
                            LimnoParam%oNFoodZoo = LimnoParam%oNFoodZoo + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sNPhytW(iLayer,iElem,gg)
                            LimnoParam%oPFoodZoo = LimnoParam%oPFoodZoo + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sPPhytW(iLayer,iElem,gg)
                            LimnoParam%oDFoodZooPhyt = LimnoParam%oDFoodZooPhyt + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            aDOMW = aDOMW + LimnoParam%sDPhytW(iLayer,iElem,gg)
                        EndIf
                    EndDo
                    ! Bacterioplankton_for_zooplankton
                    Do bb = 1, LimnoParam%nbac
                        If (LimnoParam%abac(bb)==LimnoParam%AEROBICA.or.LimnoParam%abac(bb)==LimnoParam%ANAEROBICA.or.LimnoParam%abac(bb)==LimnoParam%BACTERIA) Then
                            LimnoParam%oDFoodZoo = LimnoParam%oDFoodZoo + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sDBacW(iLayer,iElem,bb)
                            LimnoParam%oCFoodZoo = LimnoParam%oCFoodZoo + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sCBacW(iLayer,iElem,bb)
                            LimnoParam%oNFoodZoo = LimnoParam%oNFoodZoo + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sNBacW(iLayer,iElem,bb)
                            LimnoParam%oPFoodZoo = LimnoParam%oPFoodZoo + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sPBacW(iLayer,iElem,bb)
                            LimnoParam%oDFoodZooBac = LimnoParam%oDFoodZooBac + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sDBacW(iLayer,iElem,bb)
                            aDOMW = aDOMW + LimnoParam%sDBacW(iLayer,iElem,bb)
                        EndIf
                    EndDo
                    ! POM_for_zooplankton
                    Do pom = 1, LimnoParam%npom
                        If (pom==LimnoParam%LABIL.or.pom==LimnoParam%REFRATARIA) Then
                            LimnoParam%oDFoodZoo = LimnoParam%oDFoodZoo + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            LimnoParam%oCFoodZoo = LimnoParam%oCFoodZoo + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sCPomW(iLayer,iElem,pom)
                            LimnoParam%oNFoodZoo = LimnoParam%oNFoodZoo + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sNPomW(iLayer,iElem,pom)
                            LimnoParam%oPFoodZoo = LimnoParam%oPFoodZoo + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sPPomW(iLayer,iElem,pom)
                            LimnoParam%oDFoodZooPom = LimnoParam%oDFoodZooPom + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            aDOMW = aDOMW + LimnoParam%sDPomW(iLayer,iElem,pom)
                        EndIf
                    EndDo
                ElseIf (LimnoParam%azoo(zz)==2.or.LimnoParam%azoo(zz)==4) Then !Meso-zooplankton or General zooplankton
                    !  food_for_zooplankton
                    LimnoParam%oDFoodZoo = 0.
                    !  total_phytoplankton_for_zooplankton
                    LimnoParam%oDFoodZooPhyt = 0.
                    !  total_bacterioplankton_for_zooplankton
                    LimnoParam%oDFoodZooBac = 0.
                    !  total_POM_for_zooplankton
                    LimnoParam%oDFoodZooPom = 0.
                    !  organic_seston
                    LimnoParam%oDOMW = 0.
                    ! Phyto_for_zooplankton
                    Do gg = 1, LimnoParam%nphy
                        If (LimnoParam%aphy(gg)==LimnoParam%DINOF.or.LimnoParam%aphy(gg)==LimnoParam%CYANO.or.LimnoParam%aphy(gg)==LimnoParam%NODUL.or.LimnoParam%aphy(gg)==LimnoParam%CHLOR.or.LimnoParam%aphy(gg)==LimnoParam%CRYPT.or.LimnoParam%aphy(gg)==LimnoParam%MDIAT.or.LimnoParam%aphy(gg)==LimnoParam%FDIAT) Then
                            LimnoParam%oDFoodZoo = LimnoParam%oDFoodZoo + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            LimnoParam%oCFoodZoo = LimnoParam%oCFoodZoo + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sCPhytW(iLayer,iElem,gg)
                            LimnoParam%oNFoodZoo = LimnoParam%oNFoodZoo + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sNPhytW(iLayer,iElem,gg)
                            LimnoParam%oPFoodZoo = LimnoParam%oPFoodZoo + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sPPhytW(iLayer,iElem,gg)
                            LimnoParam%oDFoodZooPhyt = LimnoParam%oDFoodZooPhyt + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            aDOMW = aDOMW + LimnoParam%sDPhytW(iLayer,iElem,gg)
                        EndIf
                    EndDo
                    ! Bacterioplankton_for_zooplankton
                    Do bb = 1, LimnoParam%nbac
                        If (LimnoParam%abac(bb)==LimnoParam%AEROBICA.or.LimnoParam%abac(bb)==LimnoParam%ANAEROBICA.or.LimnoParam%abac(bb)==LimnoParam%BACTERIA) Then
                            LimnoParam%oDFoodZoo = LimnoParam%oDFoodZoo + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sDBacW(iLayer,iElem,bb)
                            LimnoParam%oCFoodZoo = LimnoParam%oCFoodZoo + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sCBacW(iLayer,iElem,bb)
                            LimnoParam%oNFoodZoo = LimnoParam%oNFoodZoo + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sNBacW(iLayer,iElem,bb)
                            LimnoParam%oPFoodZoo = LimnoParam%oPFoodZoo + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sPBacW(iLayer,iElem,bb)
                            LimnoParam%oDFoodZooBac = LimnoParam%oDFoodZooBac + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sDBacW(iLayer,iElem,bb)
                            aDOMW = aDOMW + LimnoParam%sDBacW(iLayer,iElem,bb)
                        EndIf
                    EndDo
                    ! POM_for_zooplankton
                    Do pom = 1, LimnoParam%npom
                        If (pom==LimnoParam%LABIL.or.pom==LimnoParam%REFRATARIA) Then
                            LimnoParam%oDFoodZoo = LimnoParam%oDFoodZoo + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            LimnoParam%oCFoodZoo = LimnoParam%oCFoodZoo + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sCPomW(iLayer,iElem,pom)
                            LimnoParam%oNFoodZoo = LimnoParam%oNFoodZoo + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sNPomW(iLayer,iElem,pom)
                            LimnoParam%oPFoodZoo = LimnoParam%oPFoodZoo + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sPPomW(iLayer,iElem,pom)
                            LimnoParam%oDFoodZooPom = LimnoParam%oDFoodZooPom + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            aDOMW = aDOMW + LimnoParam%sDPomW(iLayer,iElem,pom)
                        EndIf
                    EndDo
                    !! Zooplankton_for_zooplankton
                    !Do zg = 1, nzoo
                    !    If (azoo(zg)==1) Then
                    !        oDFoodZoo = oDFoodZoo + cPrefZooZoo(zz,zg)*sDZoo(iLayer,iElem,zg)
                    !         oDOMW = oDOMW + sDZoo(iLayer,iElem,zg)
                    !    EndIf
                    !EndDo

                ElseIf (LimnoParam%azoo(zz)==3) Then !Macro-zooplankton
                    !  food_for_zooplankton
                    LimnoParam%oDFoodZoo = 0.
                    !  total_phytoplankton_for_zooplankton
                    LimnoParam%oDFoodZooPhyt = 0.
                    !  total_bacterioplankton_for_zooplankton
                    LimnoParam%oDFoodZooBac = 0.
                    !  total_POM_for_zooplankton
                    LimnoParam%oDFoodZooPom = 0.
                    !  organic_seston
                    aDOMW = 0.
                    ! Phyto_for_zooplankton
                    Do gg = 1, LimnoParam%nphy
                        If (LimnoParam%aphy(gg)==LimnoParam%DINOF.or.LimnoParam%aphy(gg)==LimnoParam%CYANO.or.LimnoParam%aphy(gg)==LimnoParam%NODUL.or.LimnoParam%aphy(gg)==LimnoParam%CHLOR.or.LimnoParam%aphy(gg)==LimnoParam%CRYPT.or.LimnoParam%aphy(gg)==LimnoParam%MDIAT.or.LimnoParam%aphy(gg)==LimnoParam%FDIAT) Then
                            LimnoParam%oDFoodZoo = LimnoParam%oDFoodZoo + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            LimnoParam%oCFoodZoo = LimnoParam%oCFoodZoo + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sCPhytW(iLayer,iElem,gg)
                            LimnoParam%oNFoodZoo = LimnoParam%oNFoodZoo + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sNPhytW(iLayer,iElem,gg)
                            LimnoParam%oPFoodZoo = LimnoParam%oPFoodZoo + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sPPhytW(iLayer,iElem,gg)
                            LimnoParam%oDFoodZooPhyt = LimnoParam%oDFoodZooPhyt + LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            aDOMW = aDOMW + LimnoParam%sDPhytW(iLayer,iElem,gg)
                        EndIf
                    EndDo
                    ! Bacterioplankton_for_zooplankton
                    Do bb = 1, LimnoParam%nbac
                        If (LimnoParam%abac(bb)==LimnoParam%AEROBICA.or.LimnoParam%abac(bb)==LimnoParam%ANAEROBICA.or.LimnoParam%abac(bb)==LimnoParam%BACTERIA) Then
                            LimnoParam%oDFoodZoo = LimnoParam%oDFoodZoo + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sDBacW(iLayer,iElem,bb)
                            LimnoParam%oCFoodZoo = LimnoParam%oCFoodZoo + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sCBacW(iLayer,iElem,bb)
                            LimnoParam%oNFoodZoo = LimnoParam%oNFoodZoo + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sNBacW(iLayer,iElem,bb)
                            LimnoParam%oPFoodZoo = LimnoParam%oPFoodZoo + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sPBacW(iLayer,iElem,bb)
                            LimnoParam%oDFoodZooBac = LimnoParam%oDFoodZooBac + LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sDBacW(iLayer,iElem,bb)
                            aDOMW = aDOMW + LimnoParam%sDBacW(iLayer,iElem,bb)
                        EndIf
                    EndDo
                    ! POM_for_zooplankton
                    Do pom = 1, LimnoParam%npom
                        If (pom==LimnoParam%LABIL.or.pom==LimnoParam%REFRATARIA) Then
                            LimnoParam%oDFoodZoo = LimnoParam%oDFoodZoo + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            LimnoParam%oCFoodZoo = LimnoParam%oCFoodZoo + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sCPomW(iLayer,iElem,pom)
                            LimnoParam%oNFoodZoo = LimnoParam%oNFoodZoo + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sNPomW(iLayer,iElem,pom)
                            LimnoParam%oPFoodZoo = LimnoParam%oPFoodZoo + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sPPomW(iLayer,iElem,pom)
                            LimnoParam%oDFoodZooPom = LimnoParam%oDFoodZooPom + LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            aDOMW = aDOMW + LimnoParam%sDPomW(iLayer,iElem,pom)
                        EndIf
                    EndDo
                    ! Zooplankton_for_zooplankton
                    Do zg = 1, LimnoParam%nzoo
                        If (LimnoParam%azoo(zg)==1.or.LimnoParam%azoo(zg)==2) Then
                            LimnoParam%oDFoodZoo = LimnoParam%oDFoodZoo + LimnoParam%cPrefZooZoo(zz,zg)*LimnoParam%sDZoo(iLayer,iElem,zg)
                            LimnoParam%oCFoodZoo = LimnoParam%oCFoodZoo + LimnoParam%cPrefZooZoo(zz,zg)*LimnoParam%sCZoo(iLayer,iElem,zg)
                            LimnoParam%oNFoodZoo = LimnoParam%oNFoodZoo + LimnoParam%cPrefZooZoo(zz,zg)*LimnoParam%sNZoo(iLayer,iElem,zg)
                            LimnoParam%oPFoodZoo = LimnoParam%oPFoodZoo + LimnoParam%cPrefZooZoo(zz,zg)*LimnoParam%sPZoo(iLayer,iElem,zg)
                            aDOMW = aDOMW + LimnoParam%sDZoo(iLayer,iElem,zg)
                        EndIf
                    EndDo
                Else
                    Stop 'Zooplankton Group does not exist'
                EndIf
                
                !  max._assimilation_rate_of_zooplankton,temp._corrected
                LimnoParam%ukDAssTmZoo = LimnoParam%fDAssZoo(zz) * LimnoParam%cFiltMax(zz) * LimnoParam%uFunTmZoo * LimnoParam%hFilt(zz)
                !  food_saturation_function_of_zooplankton
                LimnoParam%aDSatZoo = LimnoParam%oDFoodZoo /(LimnoParam%hFilt(zz) + aDOMW)
                !  respiration_constant_of_zooplankton
                LimnoParam%ukDRespTmZoo = LimnoParam%kDRespZoo(zz) * LimnoParam%uFunTmZoo
                !  intrinsic_rate_of_increase_of_zooplankton
                LimnoParam%ukDIncrZoo = LimnoParam%ukDAssTmZoo - LimnoParam%ukDRespTmZoo - LimnoParam%kMortZoo(zz)
                !  environmental_correction_of_zooplankton
                LimnoParam%wDEnvZoo = max(0.0,LimnoParam%ukDIncrZoo / LimnoParam%cDCarrZoo(zz) * LimnoParam%sDZoo(iLayer,iElem,zz)*LimnoParam%sDZoo(iLayer,iElem,zz))
                !  assimilation_of_zooplankton
                LimnoParam%wDAssZoo(iLayer,iElem,zz) = LimnoParam%aDSatZoo *(LimnoParam%ukDAssTmZoo * LimnoParam%sDZoo(iLayer,iElem,zz) - LimnoParam%wDEnvZoo)
                
                LimnoParam%wDConsZoo = LimnoParam%wDAssZoo(iLayer,iElem,zz) / LimnoParam%fDAssZoo(zz) 
                !  DW_phytoplankton_consumption_by_zooplankton
                Do gg = 1, LimnoParam%nphy
                    LimnoParam%wDGrzZooPhyt(iLayer,iElem,zz,gg) = LimnoParam%cPrefZooPhyt(zz,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)/LimnoParam%oDFoodZoo * LimnoParam%wDConsZoo
                    LimnoParam%wCGrzZooPhyt(iLayer,iElem,zz,gg) = LimnoParam%rCDPhytW(iLayer,iElem,gg)*LimnoParam%wDGrzZooPhyt(iLayer,iElem,zz,gg)
                    LimnoParam%wNGrzZooPhyt(iLayer,iElem,zz,gg) = LimnoParam%rNDPhytW(iLayer,iElem,gg)*LimnoParam%wDGrzZooPhyt(iLayer,iElem,zz,gg)
                    LimnoParam%wPGrzZooPhyt(iLayer,iElem,zz,gg) = LimnoParam%rPDPhytW(iLayer,iElem,gg)*LimnoParam%wDGrzZooPhyt(iLayer,iElem,zz,gg)
                EndDo
                !  DW_bacterioplankton_consumption_by_zooplankton
                Do bb = 1, LimnoParam%nbac
                    LimnoParam%wDGrzZooBac(iLayer,iElem,zz,bb) = LimnoParam%cPrefZooBac(zz,bb)*LimnoParam%sDBacW(iLayer,iElem,bb)/LimnoParam%oDFoodZoo * LimnoParam%wDConsZoo
                    LimnoParam%wCGrzZooBac(iLayer,iElem,zz,bb) = LimnoParam%rCDBacW(iLayer,iElem,bb)*LimnoParam%wDGrzZooBac(iLayer,iElem,zz,bb)
                    LimnoParam%wNGrzZooBac(iLayer,iElem,zz,bb) = LimnoParam%rNDBacW(iLayer,iElem,bb)*LimnoParam%wDGrzZooBac(iLayer,iElem,zz,bb)
                    LimnoParam%wPGrzZooBac(iLayer,iElem,zz,bb) = LimnoParam%rPDBacW(iLayer,iElem,bb)*LimnoParam%wDGrzZooBac(iLayer,iElem,zz,bb)
                EndDo
                !  DW_detritus_consumption_by_zooplankton
                Do pom = 1, LimnoParam%npom
                    LimnoParam%wDGrzZooPom(iLayer,iElem,zz,pom) = LimnoParam%cPrefZooPom(zz,pom)*LimnoParam%sDPomW(iLayer,iElem,pom) / LimnoParam%oDFoodZoo * LimnoParam%wDConsZoo
                    LimnoParam%wCGrzZooPom(iLayer,iElem,zz,pom) = LimnoParam%rCDPomW(iLayer,iElem,pom)*LimnoParam%wDGrzZooPom(iLayer,iElem,zz,pom)
                    LimnoParam%wNGrzZooPom(iLayer,iElem,zz,pom) = LimnoParam%rNDPomW(iLayer,iElem,pom)*LimnoParam%wDGrzZooPom(iLayer,iElem,zz,pom)
                    LimnoParam%wPGrzZooPom(iLayer,iElem,zz,pom) = LimnoParam%rPDPomW(iLayer,iElem,pom)*LimnoParam%wDGrzZooPom(iLayer,iElem,zz,pom)
                EndDo
                !  DW_zooplankton(zg)_consumption_by_zooplankton(zz)
                Do zg = 1, LimnoParam%nzoo
                    LimnoParam%wDGrzZooZoo(iLayer,iElem,zz,zg) = LimnoParam%cPrefZooZoo(zz,zg)*LimnoParam%sDZoo(iLayer,iElem,zg) / LimnoParam%oDFoodZoo * LimnoParam%wDConsZoo
                    LimnoParam%wCGrzZooZoo(iLayer,iElem,zz,zg) = LimnoParam%rCDZoo(iLayer,iElem,zg)*LimnoParam%wDGrzZooZoo(iLayer,iElem,zz,zg)
                    LimnoParam%wNGrzZooZoo(iLayer,iElem,zz,zg) = LimnoParam%rNDZoo(iLayer,iElem,zg)*LimnoParam%wDGrzZooZoo(iLayer,iElem,zz,zg)
                    LimnoParam%wPGrzZooZoo(iLayer,iElem,zz,zg) = LimnoParam%rPDZoo(iLayer,iElem,zg)*LimnoParam%wDGrzZooZoo(iLayer,iElem,zz,zg)
                EndDo
                !-----------------------------------------------------------------------
                !  zooplankton assimilation C
                !-----------------------------------------------------------------------
                !  C/DW_ratio_of_zooplankton_food
                LimnoParam%rCDFoodZoo = LimnoParam%oCFoodZoo /(LimnoParam%oDFoodZoo+NearZero)
                !  total_C_consumption
                LimnoParam%wCConsZoo = (sum(LimnoParam%wCGrzZooPhyt(iLayer,iElem,zz,:)) + sum(LimnoParam%wCGrzZooBac(iLayer,iElem,zz,:)) + sum(LimnoParam%wCGrzZooPom(iLayer,iElem,zz,:)) + sum(LimnoParam%wCGrzZooZoo(iLayer,iElem,zz,:)))
                !  C_assimilation_efficiency_of_herbivores
                LimnoParam%afCAssZoo = min(1.0,LimnoParam%cCDZooRef(zz) / LimnoParam%rCDFoodZoo * LimnoParam%fDAssZoo(zz))
                !  assimilation_by_herbivores
                LimnoParam%wCAssZoo(iLayer,iElem,zz) = LimnoParam%afCAssZoo*LimnoParam%wCConsZoo
                
                !-----------------------------------------------------------------------
                !  zooplankton assimilation N
                !-----------------------------------------------------------------------
                !  N/DW_ratio_of_zooplankton_food
                LimnoParam%rNDFoodZoo = LimnoParam%oNFoodZoo /(LimnoParam%oDFoodZoo+NearZero)
                !  total_N_consumption
                LimnoParam%wNConsZoo = (sum(LimnoParam%wNGrzZooPhyt(iLayer,iElem,zz,:)) + sum(LimnoParam%wNGrzZooBac(iLayer,iElem,zz,:)) + sum(LimnoParam%wNGrzZooPom(iLayer,iElem,zz,:)) + sum(LimnoParam%wNGrzZooZoo(iLayer,iElem,zz,:)))
                !  N_assimilation_efficiency_of_herbivores
                LimnoParam%afNAssZoo = min(1.0,LimnoParam%cNDZooRef(zz) / LimnoParam%rNDFoodZoo * LimnoParam%fDAssZoo(zz))
                !  assimilation_by_herbivores
                LimnoParam%wNAssZoo(iLayer,iElem,zz) = LimnoParam%afNAssZoo*LimnoParam%wNConsZoo

                !-----------------------------------------------------------------------
                !  zooplankton assimilation P
                !-----------------------------------------------------------------------
                !  N/DW_ratio_of_zooplankton_food
                LimnoParam%rPDFoodZoo = LimnoParam%oPFoodZoo /(LimnoParam%oDFoodZoo+NearZero)
                !  total_N_consumption
                LimnoParam%wPConsZoo = (sum(LimnoParam%wPGrzZooPhyt(iLayer,iElem,zz,:)) + sum(LimnoParam%wPGrzZooBac(iLayer,iElem,zz,:)) + sum(LimnoParam%wPGrzZooPom(iLayer,iElem,zz,:)) + sum(LimnoParam%wPGrzZooZoo(iLayer,iElem,zz,:)))
                !  N_assimilation_efficiency_of_herbivores
                LimnoParam%afPAssZoo = min(1.0,LimnoParam%cPDZooRef(zz) / LimnoParam%rPDFoodZoo * LimnoParam%fDAssZoo(zz))
                !  assimilation_by_herbivores
                LimnoParam%wPAssZoo(iLayer,iElem,zz) = LimnoParam%afPAssZoo*LimnoParam%wPConsZoo
                
                
                !-------------------------------------------------------------------------------------!
			    ! 2) Messy -> Assimilation => POM
                !-------------------------------------------------------------------------------------!
                 LimnoParam%wDMessZooOM(iLayer,iElem,zz) = LimnoParam%wDConsZoo - LimnoParam%wDAssZoo(iLayer,iElem,zz)
                 LimnoParam%wCMessZooOM(iLayer,iElem,zz) = LimnoParam%wCConsZoo - LimnoParam%wCAssZoo(iLayer,iElem,zz)
                 LimnoParam%wNMessZooOM(iLayer,iElem,zz) = LimnoParam%wNConsZoo - LimnoParam%wNAssZoo(iLayer,iElem,zz)
                 LimnoParam%wPMessZooOM(iLayer,iElem,zz) = LimnoParam%wPConsZoo - LimnoParam%wPAssZoo(iLayer,iElem,zz)
                
                
                !-----------------------------------------------------------------------
                !  3) zooplankton respiration and excretion
                !-----------------------------------------------------------------------
                !  corr._factor_of_zoopl._respiration_for_P_and_N_content
                LimnoParam%aCorDRespZoo = max(LimnoParam%cPDZooRef(zz) / LimnoParam%rPDZoo(iLayer,iElem,zz),LimnoParam%cNDZooRef(zz) / LimnoParam%rNDZoo(iLayer,iElem,zz))
                !  zoopl._respiration_DW (DIC)
                LimnoParam%wDRespZoo(iLayer,iElem,zz) = LimnoParam%aCorDRespZoo * LimnoParam%kDRespZoo(zz) * LimnoParam%uFunTmZoo * LimnoParam%sDZoo(iLayer,iElem,zz)
                !  C_excretion
                LimnoParam%wCRespZoo(iLayer,iElem,zz) = LimnoParam%rCDZoo(iLayer,iElem,zz) / LimnoParam%cCDZooRef(zz) * LimnoParam%kDRespZoo(zz) * LimnoParam%uFunTmZoo*LimnoParam%sCZoo(iLayer,iElem,zz)
                !  N_excretion
                LimnoParam%wNExcrZoo(iLayer,iElem,zz) = LimnoParam%rNDZoo(iLayer,iElem,zz) / LimnoParam%cNDZooRef(zz) * LimnoParam%kDRespZoo(zz) * LimnoParam%uFunTmZoo*LimnoParam%sNZoo(iLayer,iElem,zz)
                !  P_excretion
                LimnoParam%wPExcrZoo(iLayer,iElem,zz) = LimnoParam%rPDZoo(iLayer,iElem,zz) / LimnoParam%cPDZooRef(zz) * LimnoParam%kDRespZoo(zz) * LimnoParam%uFunTmZoo*LimnoParam%sPZoo(iLayer,iElem,zz)
                !-----------------------------------------------------------------------
                !  4) zooplankton Egestion/Pellets (-) -> POM
                !-----------------------------------------------------------------------
				LimnoParam%wDEgesZoo(iLayer,iElem,zz) = LimnoParam%kPelZoo(zz) * LimnoParam%sDZoo(iLayer,iElem,zz) 
			    LimnoParam%wCEgesZoo(iLayer,iElem,zz) = LimnoParam%rCDZoo(iLayer,iElem,zz) * LimnoParam%wDEgesZoo(iLayer,iElem,zz) 
			    LimnoParam%wNEgesZoo(iLayer,iElem,zz) = LimnoParam%rNDZoo(iLayer,iElem,zz) * LimnoParam%wDEgesZoo(iLayer,iElem,zz)
			    LimnoParam%wPEgesZoo(iLayer,iElem,zz) = LimnoParam%rPDZoo(iLayer,iElem,zz) * LimnoParam%wDEgesZoo(iLayer,iElem,zz)
                
                !-----------------------------------------------------------------------
                !  5) zooplankton mortality
                !-----------------------------------------------------------------------
                !  zoopl._mortality,incl._environmental_correction_DW
                LimnoParam%wDMortZoo(iLayer,iElem,zz) = LimnoParam%kMortZoo(zz) * LimnoParam%sDZoo(iLayer,iElem,zz) +(1.0 - LimnoParam%aDSatZoo) * LimnoParam%wDEnvZoo
                !  zoopl._mortality_C
                LimnoParam%wCMortZoo(iLayer,iElem,zz) = LimnoParam%rCDZoo(iLayer,iElem,zz) * LimnoParam%wDMortZoo(iLayer,iElem,zz)
                !  zoopl._mortality_N
                LimnoParam%wNMortZoo(iLayer,iElem,zz) = LimnoParam%rNDZoo(iLayer,iElem,zz) * LimnoParam%wDMortZoo(iLayer,iElem,zz)
                !  zoopl._mortality_P
                LimnoParam%wPMortZoo(iLayer,iElem,zz) = LimnoParam%rPDZoo(iLayer,iElem,zz) * LimnoParam%wDMortZoo(iLayer,iElem,zz)
                
            
                !-------------------------------------------------------------------------------------!
                !                       Source/Sink Term for Zooplankton						      !
                !-------------------------------------------------------------------------------------!
              
                HydroParam%dVarEst(iLayer,iElem,1) =  LimnoParam%wDAssZoo(iLayer,iElem,zz)                          & !Assimilation
                                        - LimnoParam%wDRespZoo(iLayer,iElem,zz)                                     & !Respiration
                                        - LimnoParam%wDEgesZoo(iLayer,iElem,zz)                                     & !Egestion
                                        - LimnoParam%wDMortZoo(iLayer,iElem,zz)					                    & !Mortality
                                        - SUM(LimnoParam%wDGrzZooZoo(iLayer,iElem,:,zz))                            & !Predation by Zoops
                                        - SUM(LimnoParam%wDGrzFiAdZoo(iLayer,iElem,:,zz))                           & !Predation by Adult Fish
                                        - SUM(LimnoParam%wDGrzFiJvZoo(iLayer,iElem,:,zz))                           !Predation by Juvenile Fish
             
                HydroParam%dVarEst(iLayer,iElem,2) =  LimnoParam%wCAssZoo(iLayer,iElem,zz)                          & !Assimilation
                                        - LimnoParam%wCRespZoo(iLayer,iElem,zz)                                     & !Respiration (DIC)
                                        - LimnoParam%wCEgesZoo(iLayer,iElem,zz)                                     & !Egestion
                                        - LimnoParam%wCMortZoo(iLayer,iElem,zz)					                    & !Mortality
                                        - SUM(LimnoParam%wCGrzZooZoo(iLayer,iElem,:,zz))                            & !Predation by Zoops
                                        - SUM(LimnoParam%wCGrzFiAdZoo(iLayer,iElem,:,zz))                           & !Predation by Adult Fish
                                        - SUM(LimnoParam%wCGrzFiJvZoo(iLayer,iElem,:,zz))                           !Predation by Juvenile Fish
  
                HydroParam%dVarEst(iLayer,iElem,3) =  LimnoParam%wNAssZoo(iLayer,iElem,zz)                          & !Assimilation
                                        - LimnoParam%wNExcrZoo(iLayer,iElem,zz)                                     & !Excretion
                                        - LimnoParam%wNEgesZoo(iLayer,iElem,zz)                                     & !Egestion
                                        - LimnoParam%wNMortZoo(iLayer,iElem,zz)					                    & !Mortality
                                        - SUM(LimnoParam%wNGrzZooZoo(iLayer,iElem,:,zz))                            & !Predation by Zoops
                                        - SUM(LimnoParam%wNGrzFiAdZoo(iLayer,iElem,:,zz))                           & !Predation by Adult Fish
                                        - SUM(LimnoParam%wNGrzFiJvZoo(iLayer,iElem,:,zz))                           !Predation by Juvenile Fish

                HydroParam%dVarEst(iLayer,iElem,4) =  LimnoParam%wPAssZoo(iLayer,iElem,zz)                          & !Assimilation
                                        - LimnoParam%wPExcrZoo(iLayer,iElem,zz)                                     & !Excretion
                                        - LimnoParam%wPEgesZoo(iLayer,iElem,zz)                                     & !Egestion
                                        - LimnoParam%wPMortZoo(iLayer,iElem,zz)					                    & !Mortality
                                        - SUM(LimnoParam%wPGrzZooZoo(iLayer,iElem,:,zz))                            & !Predation by Zoops
                                        - SUM(LimnoParam%wPGrzFiAdZoo(iLayer,iElem,:,zz))                           & !Predation by Adult Fish
                                        - SUM(LimnoParam%wPGrzFiJvZoo(iLayer,iElem,:,zz))                           !Predation by Juvenile Fish
                
            EndDo !Loop Layer

        EndDo !Loop Cell
        
      
        !-------------------------------------------------------------------------------------!
        !                   Solver Transport Equation for Zooplankton                         !
        !-------------------------------------------------------------------------------------!
        Index = 0
        If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDZoo(:,:,zz),LimnoParam%sDZooP(:,:,zz),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCZoo(:,:,zz),LimnoParam%sCZooP(:,:,zz),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNZoo(:,:,zz),LimnoParam%sNZooP(:,:,zz),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPZoo(:,:,zz),LimnoParam%sPZooP(:,:,zz),dt,dtday,HydroParam,MeshParam)
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDZoo(:,:,zz),LimnoParam%sDZooP(:,:,zz),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCZoo(:,:,zz),LimnoParam%sCZooP(:,:,zz),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNZoo(:,:,zz),LimnoParam%sNZooP(:,:,zz),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPZoo(:,:,zz),LimnoParam%sPZooP(:,:,zz),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
        ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDZoo(:,:,zz),LimnoParam%sDZooP(:,:,zz),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCZoo(:,:,zz),LimnoParam%sCZooP(:,:,zz),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNZoo(:,:,zz),LimnoParam%sNZooP(:,:,zz),dt,dtday,HydroParam,MeshParam)
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPZoo(:,:,zz),LimnoParam%sPZooP(:,:,zz),dt,dtday,HydroParam,MeshParam)
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDZoo(:,:,zz),LimnoParam%sDZooP(:,:,zz),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sCZoo(:,:,zz),LimnoParam%sCZooP(:,:,zz),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,3),LimnoParam%sNZoo(:,:,zz),LimnoParam%sNZooP(:,:,zz),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,4),LimnoParam%sPZoo(:,:,zz),LimnoParam%sPZooP(:,:,zz),dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
        EndiF     
        

    EndDo !Loop Functional Groups
    
Return
End
