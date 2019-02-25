SUBROUTINE Fishes(HydroParam,MeshParam,LimnoParam,dt,dtday,Julday)

    ! Fishes Modelling in the Water Column
    
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
    Integer:: index,iElem,iLayer,pom,zg,ben,gg,ff,fg,Julday
    Real:: V
    Real:: dt,dtday
    Real:: aDOMWAd,aDOMWJv
    Real:: NearZero = 1e-10
    
    Do ff = 1, LimnoParam%nfish
        Do iElem = 1,MeshParam%nElem
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                !-------------------------------------------------------------------------------------!
                !       Fishes (Piscivorous, Planktivorous and Omnivorous) Processes in Water         !
                !-------------------------------------------------------------------------------------!
                !1) Interannual Migration (+/-)                                                       !
                !2) Reprodution/Aging(+/-) JUV->AD                                                    !
                !3) Harvest(-)                                                                        !
                !4) Predation (+)                                                                     !
                !5) Messy Fedding (-) -> POM                                                          !
                !6) Respiration/Excretion(-) -> CO2, N(DOM),P(DOM)                                    !
                !7) Pellets(-) -> POM                                                                 !
                !8) Mortality (-) -> POM                                                              !
                !9) Predation (-)                                                                     !
                !-------------------------------------------------------------------------------------!			
                !-------------------------------------------------------------------------------------!
                
                LimnoParam%uFunTmFish = uFunTmBio(LimnoParam%sDTempW(iLayer,iElem),LimnoParam%cSigTmFish(ff),LimnoParam%cTmOptFish(ff))
            
			    !1) Interannual Migration (gD/m2/d)  
                !-----------------------------------------------------------------------
                !  adult fish migration
                !-----------------------------------------------------------------------
                !  migration_flux of adult fish,DW
                !  PCLake_osis:sDFiAd,in g/m^2
                LimnoParam%wDMigrFiAd(ff) = LimnoParam%kMigrFiAd(ff) *(LimnoParam%cDFiAdIn(ff) - LimnoParam%sDFiAd(iLayer,iElem,ff))  
                !  net_migration_flux of adult fish, C
                !  PCLake_osis:sCFiAd,in g/m^2
                LimnoParam%wCMigrFiAd(ff) = LimnoParam%cCDFishRef(ff)*LimnoParam%wDMigrFiAd(ff) !kMigrFish(ff) *(rCDFiAd(ff) * cDFishIn(ff) - sCFiAd(iLayer,iElem,ff))   
                !  net_migration_flux of adult fish, N
                !  PCLake_osis:sPFiAd,in g/m^2
                LimnoParam%wNMigrFiAd(ff) = LimnoParam%cNDFishRef(ff)*LimnoParam%wDMigrFiAd(ff) !kMigrFish(ff) *(rNDFiAd(ff) * cDFishIn(ff) - sNFiAd(iLayer,iElem,ff))   
                !  net_migration_flux of adult fish, P
                !  PCLake_osis:sPFiAd,in g/m^2
                LimnoParam%wPMigrFiAd(ff) = LimnoParam%cPDFishRef(ff)*LimnoParam%wDMigrFiAd(ff) !kMigrFish(ff) *(rPDFiAd(ff) * cDFishIn(ff) - sPFiAd(iLayer,iElem,ff))
                !-----------------------------------------------------------------------
                !  young fish migration
                !-----------------------------------------------------------------------
                !  migration_flux of young fish, DW
                !  PCLake_osis:sDFiJv,in g/m^2
                If (LimnoParam%afish(ff) == 1) Then !Piscivorous fish (afish = 1)
                    LimnoParam%wDMigrFiJv(ff) = 0.
                Else
                    LimnoParam%wDMigrFiJv(ff) = LimnoParam%kMigrFiJv(ff) *(LimnoParam%cDFiJvIn(ff) - LimnoParam%sDFiJv(iLayer,iElem,ff))
                EndIf
                !  net_migration_flux of young fish,C
                !  PCLake_osis:sCFiJv,in g/m^2
                LimnoParam%wCMigrFiJv(ff) = LimnoParam%cCDFishRef(ff)*LimnoParam%wDMigrFiJv(ff) !kMigrFish *(self%cPDFishRef * self%cDFiJvIn - sPFiJv)
                !  net_migration_flux of young fish,N
                !  PCLake_osis:sNFiJv,in g/m^2
                LimnoParam%wNMigrFiJv(ff) = LimnoParam%cNDFishRef(ff)*LimnoParam%wDMigrFiJv(ff) !kMigrFish *(self%cNDFishRef * self%cDFiJvIn - sNFiJv)
                !  net_migration_flux of young fish,P
                !  PCLake_osis:sPFiJv,in g/m^2
                LimnoParam%wPMigrFiJv(ff) = LimnoParam%cPDFishRef(ff)*LimnoParam%wDMigrFiJv(ff) !kMigrFish *(self%cPDFishRef * self%cDFiJvIn - sPFiJv)
                
                ! 2) Reprodution/Aging(+/-) JUV<->AD
                !-----------------------------------------------------------------------
                !  fish reproduction
                !-----------------------------------------------------------------------
                !  Reproduction_flux_DW
                If (Julday >= LimnoParam%cDayReprFish(ff) .and. Julday < LimnoParam%cDayReprFish(ff) + 1.0) Then
                !  PCLake_osis:sDFiAd,in g/m^2
                   LimnoParam%wDReprFish(ff) = LimnoParam%fReprFish(ff) * LimnoParam%sDFiAd(iLayer,iElem,ff)
                Else
                   LimnoParam%wDReprFish(ff) =0.0
                EndIf
                !  Reproduction_flux_C
                LimnoParam%wCReprFish(ff) = LimnoParam%rCDFiAd(iLayer,iElem,ff) * LimnoParam%wDReprFish(ff)
                !  Reproduction_flux_N
                LimnoParam%wNReprFish(ff) = LimnoParam%rNDFiAd(iLayer,iElem,ff) * LimnoParam%wDReprFish(ff)
                !  Reproduction_flux_P
                LimnoParam%wPReprFish(ff) = LimnoParam%rPDFiAd(iLayer,iElem,ff) * LimnoParam%wDReprFish(ff)
                !-----------------------------------------------------------------------
                !  fish aging
                !-----------------------------------------------------------------------
                !  Ageing_DW
                If (Julday >=  365.0 .AND. Julday <= 365.0) Then
                !  PCLake_osis:sDFiAd,in g/m^2
                   LimnoParam%wDAgeFish(ff) = LimnoParam%fAgeFish(ff) * LimnoParam%sDFiJv(iLayer,iElem,ff)
                Else
                   LimnoParam%wDAgeFish(ff) = 0.0
                EndIf
                !  Ageing_C
                LimnoParam%wCAgeFish(ff) = LimnoParam%rCDFiJv(iLayer,iElem,ff) * LimnoParam%wDAgeFish(ff)
                !  Ageing_N
                LimnoParam%wNAgeFish(ff) = LimnoParam%rNDFiJv(iLayer,iElem,ff) * LimnoParam%wDAgeFish(ff)
                !  Ageing_P
                LimnoParam%wPAgeFish(ff) = LimnoParam%rPDFiJv(iLayer,iElem,ff) * LimnoParam%wDAgeFish(ff)
                
                !3) Harvest(-)
                !-----------------------------------------------------------------------
                !  adult fish Harvest
                !-----------------------------------------------------------------------
				LimnoParam%wDHarvFiAd(ff) = LimnoParam%kHarvFiAd(ff) * LimnoParam%sDFiAd(iLayer,iElem,ff)
				LimnoParam%wCHarvFiAd(ff) = LimnoParam%wDHarvFiAd(ff)*LimnoParam%rCDFiAd(iLayer,iElem,ff) 
				LimnoParam%wNHarvFiAd(ff) = LimnoParam%wDHarvFiAd(ff)*LimnoParam%rNDFiAd(iLayer,iElem,ff) 
				LimnoParam%wPHarvFiAd(ff) = LimnoParam%wDHarvFiAd(ff)*LimnoParam%rPDFiAd(iLayer,iElem,ff) 
                !-----------------------------------------------------------------------
                !  juvenile fish Harvest
                !-----------------------------------------------------------------------
                If (LimnoParam%afish(ff) == 1) Then !Piscivorous fish (afish = 1)
                    LimnoParam%wDHarvFiJv(ff) = 0.
                Else
                    LimnoParam%wDHarvFiJv(ff) = LimnoParam%kHarvFiJv(ff) * LimnoParam%sDFiJv(iLayer,iElem,ff)
                EndIf
				LimnoParam%wCHarvFiJv(ff) = LimnoParam%wDHarvFiJv(ff)*LimnoParam%rCDFiJv(iLayer,iElem,ff) 
				LimnoParam%wNHarvFiJv(ff) = LimnoParam%wDHarvFiJv(ff)*LimnoParam%rNDFiJv(iLayer,iElem,ff) 
				LimnoParam%wPHarvFiJv(ff) = LimnoParam%wDHarvFiJv(ff)*LimnoParam%rPDFiJv(iLayer,iElem,ff) 
                
                !4) Predation (Assimilation)  
                If (LimnoParam%afish(ff) == 1) Then !Piscivorous fish (afish = 1)
                    !!  total_food_for_fish
                    LimnoParam%oDFoodFiAd = 0.
                    LimnoParam%oCFoodFiAd = 0.
                    LimnoParam%oNFoodFiAd = 0.
                    LimnoParam%oPFoodFiAd = 0.
                    LimnoParam%oDFoodFiJv = 0.
                    LimnoParam%oCFoodFiJv = 0.
                    LimnoParam%oNFoodFiJv = 0.
                    LimnoParam%oPFoodFiJv = 0.

                    !  total_fish_for_fish
                    LimnoParam%oDFoodFiAdFish = 0.
                    LimnoParam%oDFoodFiAdZoo = 0.
                    LimnoParam%oDFoodFiJvZoo = 0.
                    LimnoParam%oDFoodFiAdBent = 0.
                    LimnoParam%oDFoodFiJvBent = 0.
                    LimnoParam%oDFoodFiAdPhyt = 0.
                    LimnoParam%oDFoodFiJvPhyt = 0.
                    LimnoParam%oDFoodFiAdPom = 0.
                    LimnoParam%oDFoodFiJvPom = 0.
                    !  organic_seston
                    aDOMWAd = 0.
                    ! Fish_for_fish
                    Do fg = 1, LimnoParam%nfish
                        If (LimnoParam%afish(fg)==2.or.LimnoParam%afish(fg)==3) Then ! It can predate omnivorous and planktivorous fish
                            LimnoParam%oDFoodFiAd = LimnoParam%oDFoodFiAd + LimnoParam%cPrefFiAdFiJv(ff,fg)*LimnoParam%sDFiJv(iLayer,iElem,fg) + LimnoParam%cPrefFiAdFiAd(ff,fg)*LimnoParam%sDFiAd(iLayer,iElem,fg)
                            LimnoParam%oCFoodFiAd = LimnoParam%oCFoodFiAd + LimnoParam%cPrefFiAdFiJv(ff,fg)*LimnoParam%sCFiJv(iLayer,iElem,fg) + LimnoParam%cPrefFiAdFiAd(ff,fg)*LimnoParam%sCFiAd(iLayer,iElem,fg)
                            LimnoParam%oNFoodFiAd = LimnoParam%oNFoodFiAd + LimnoParam%cPrefFiAdFiJv(ff,fg)*LimnoParam%sNFiJv(iLayer,iElem,fg) + LimnoParam%cPrefFiAdFiAd(ff,fg)*LimnoParam%sNFiAd(iLayer,iElem,fg)
                            LimnoParam%oPFoodFiAd = LimnoParam%oPFoodFiAd + LimnoParam%cPrefFiAdFiJv(ff,fg)*LimnoParam%sPFiJv(iLayer,iElem,fg) + LimnoParam%cPrefFiAdFiAd(ff,fg)*LimnoParam%sPFiAd(iLayer,iElem,fg)
                            LimnoParam%oDFoodFiAdFish = LimnoParam%oDFoodFiAdFish + LimnoParam%cPrefFiAdFiJv(ff,fg)*LimnoParam%sDFiJv(iLayer,iElem,fg) + LimnoParam%cPrefFiAdFiAd(ff,fg)*LimnoParam%sDFiAd(iLayer,iElem,fg)
                            aDOMWAd = aDOMWAd + LimnoParam%sDFiJv(iLayer,iElem,fg) + LimnoParam%sDFiAd(iLayer,iElem,fg)
                        EndIf
                    EndDo
                ElseIf  (LimnoParam%afish(ff) == 2) Then !Omnivorous fish (afish = 2)
                    !!  total_food_for_fish
                    LimnoParam%oDFoodFiAd = 0.
                    LimnoParam%oCFoodFiAd = 0.
                    LimnoParam%oNFoodFiAd = 0.
                    LimnoParam%oPFoodFiAd = 0.
                    !  total_fish_for_fish
                    LimnoParam%oDFoodFiAdFish = 0.
                    LimnoParam%oDFoodFiAdZoo = 0.
                    LimnoParam%oDFoodFiJvZoo = 0.
                    LimnoParam%oDFoodFiAdBent = 0.
                    LimnoParam%oDFoodFiJvBent = 0.
                    LimnoParam%oDFoodFiAdPhyt = 0.
                    LimnoParam%oDFoodFiJvPhyt = 0.
                    LimnoParam%oDFoodFiAdPom = 0.
                    LimnoParam%oDFoodFiJvPom = 0.
                    !  organic_seston
                    aDOMWAd = 0.
                    aDOMWJv = 0.
                    ! Fish_for_fish
                    Do fg = 1, LimnoParam%nfish
                        If (LimnoParam%afish(fg)==3) Then ! It can predate planktivorous fish
                            !Adult Fish
                            LimnoParam%oDFoodFiAd = LimnoParam%oDFoodFiAd + LimnoParam%cPrefFiAdFiJv(ff,fg)*LimnoParam%sDFiJv(iLayer,iElem,fg) + LimnoParam%cPrefFiAdFiAd(ff,fg)*LimnoParam%sDFiAd(iLayer,iElem,fg)
                            LimnoParam%oCFoodFiAd = LimnoParam%oCFoodFiAd + LimnoParam%cPrefFiAdFiJv(ff,fg)*LimnoParam%sCFiJv(iLayer,iElem,fg) + LimnoParam%cPrefFiAdFiAd(ff,fg)*LimnoParam%sCFiAd(iLayer,iElem,fg)
                            LimnoParam%oNFoodFiAd = LimnoParam%oNFoodFiAd + LimnoParam%cPrefFiAdFiJv(ff,fg)*LimnoParam%sNFiJv(iLayer,iElem,fg) + LimnoParam%cPrefFiAdFiAd(ff,fg)*LimnoParam%sNFiAd(iLayer,iElem,fg)
                            LimnoParam%oPFoodFiAd = LimnoParam%oPFoodFiAd + LimnoParam%cPrefFiAdFiJv(ff,fg)*LimnoParam%sPFiJv(iLayer,iElem,fg) + LimnoParam%cPrefFiAdFiAd(ff,fg)*LimnoParam%sPFiAd(iLayer,iElem,fg)
                            LimnoParam%oDFoodFiAdFish = LimnoParam%oDFoodFiAdFish + LimnoParam%cPrefFiAdFiJv(ff,fg)*LimnoParam%sDFiJv(iLayer,iElem,fg) + LimnoParam%cPrefFiAdFiAd(ff,fg)*LimnoParam%sDFiAd(iLayer,iElem,fg)
                            aDOMWAd = aDOMWAd + LimnoParam%sDFiJv(iLayer,iElem,fg) + LimnoParam%sDFiAd(iLayer,iElem,fg)
                        EndIf
                    EndDo
                    ! Zooplankton_for_fish
                    Do zg = 1, LimnoParam%nzoo
                        If (LimnoParam%azoo(zg)>=2) Then ! It can predate meso, macro and general zooplankton
                            !Adult Fish
                            LimnoParam%oDFoodFiAd = LimnoParam%oDFoodFiAd + LimnoParam%cPrefFiAdZoo(ff,zg)*LimnoParam%sDZoo(iLayer,iElem,zg)
                            LimnoParam%oCFoodFiAd = LimnoParam%oCFoodFiAd + LimnoParam%cPrefFiAdZoo(ff,zg)*LimnoParam%sCZoo(iLayer,iElem,zg)
                            LimnoParam%oNFoodFiAd = LimnoParam%oNFoodFiAd + LimnoParam%cPrefFiAdZoo(ff,zg)*LimnoParam%sNZoo(iLayer,iElem,zg)
                            LimnoParam%oPFoodFiAd = LimnoParam%oPFoodFiAd + LimnoParam%cPrefFiAdZoo(ff,zg)*LimnoParam%sPZoo(iLayer,iElem,zg)
                            LimnoParam%oDFoodFiAdZoo = LimnoParam%oDFoodFiAdZoo + LimnoParam%cPrefFiAdZoo(ff,zg)*LimnoParam%sDZoo(iLayer,iElem,zg)
                            aDOMWAd = aDOMWAd + LimnoParam%sDZoo(iLayer,iElem,zg)
                            !Juveline Fish
                            LimnoParam%oDFoodFiJv = LimnoParam%oDFoodFiJv + LimnoParam%cPrefFiJvZoo(ff,zg)*LimnoParam%sDZoo(iLayer,iElem,zg)
                            LimnoParam%oCFoodFiJv = LimnoParam%oCFoodFiJv + LimnoParam%cPrefFiJvZoo(ff,zg)*LimnoParam%sCZoo(iLayer,iElem,zg)
                            LimnoParam%oNFoodFiJv = LimnoParam%oNFoodFiJv + LimnoParam%cPrefFiJvZoo(ff,zg)*LimnoParam%sNZoo(iLayer,iElem,zg)
                            LimnoParam%oPFoodFiJv = LimnoParam%oPFoodFiJv + LimnoParam%cPrefFiJvZoo(ff,zg)*LimnoParam%sPZoo(iLayer,iElem,zg)
                            LimnoParam%oDFoodFiJvZoo = LimnoParam%oDFoodFiJvZoo + LimnoParam%cPrefFiJvZoo(ff,zg)*LimnoParam%sDZoo(iLayer,iElem,zg)
                            aDOMWJv = aDOMWJv + LimnoParam%sDZoo(iLayer,iElem,zg)
                        EndIf
                    EndDo
                    ! Zoobenthos_for_fish
                    Do ben = 1, LimnoParam%nben
                        If (iLayer == HydroParam%ElSmallm(iElem)) Then
                            If (LimnoParam%aben(ben)>=2) Then ! It can predate meso, macro and general zoobenthos
                                !Adult Fish
                                LimnoParam%oDFoodFiAd = LimnoParam%oDFoodFiAd + LimnoParam%cPrefFiAdBent(ff,ben)*LimnoParam%sDBent(iElem,ben)
                                LimnoParam%oCFoodFiAd = LimnoParam%oCFoodFiAd + LimnoParam%cPrefFiAdBent(ff,ben)*LimnoParam%sCBent(iElem,ben)
                                LimnoParam%oNFoodFiAd = LimnoParam%oNFoodFiAd + LimnoParam%cPrefFiAdBent(ff,ben)*LimnoParam%sNBent(iElem,ben)
                                LimnoParam%oPFoodFiAd = LimnoParam%oPFoodFiAd + LimnoParam%cPrefFiAdBent(ff,ben)*LimnoParam%sPBent(iElem,ben)
                                LimnoParam%oDFoodFiAdBent = LimnoParam%oDFoodFiAdBent + LimnoParam%cPrefFiAdBent(ff,ben)*LimnoParam%sDBent(iElem,ben)
                                aDOMWAd = aDOMWAd + LimnoParam%sDBent(iElem,ben)
                                !Juveline Fish
                                LimnoParam%oDFoodFiJv = LimnoParam%oDFoodFiJv + LimnoParam%cPrefFiJvBent(ff,ben)*LimnoParam%sDBent(iElem,ben)
                                LimnoParam%oCFoodFiJv = LimnoParam%oCFoodFiJv + LimnoParam%cPrefFiJvBent(ff,ben)*LimnoParam%sCBent(iElem,ben)
                                LimnoParam%oNFoodFiJv = LimnoParam%oNFoodFiJv + LimnoParam%cPrefFiJvBent(ff,ben)*LimnoParam%sNBent(iElem,ben)
                                LimnoParam%oPFoodFiJv = LimnoParam%oPFoodFiJv + LimnoParam%cPrefFiJvBent(ff,ben)*LimnoParam%sPBent(iElem,ben)
                                LimnoParam%oDFoodFiJvBent = LimnoParam%oDFoodFiJvBent + LimnoParam%cPrefFiJvBent(ff,ben)*LimnoParam%sDBent(iElem,ben)
                                aDOMWJv = aDOMWJv + LimnoParam%sDBent(iElem,zg)
                            EndIf
                        EndIf
                    EndDo
                    ! Phyto_for_fish
                    Do gg = 1, LimnoParam%nphy
                        If (LimnoParam%aphy(gg)==LimnoParam%DINOF.or.LimnoParam%aphy(gg)==LimnoParam%CYANO.or.LimnoParam%aphy(gg)==LimnoParam%NODUL.or.LimnoParam%aphy(gg)==LimnoParam%CHLOR.or.LimnoParam%aphy(gg)==LimnoParam%CRYPT.or.LimnoParam%aphy(gg)==LimnoParam%MDIAT.or.LimnoParam%aphy(gg)==LimnoParam%FDIAT) Then
                            !Adult Fish
                            LimnoParam%oDFoodFiAd = LimnoParam%oDFoodFiAd + LimnoParam%cPrefFiAdPhyt(ff,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            LimnoParam%oCFoodFiAd = LimnoParam%oCFoodFiAd + LimnoParam%cPrefFiAdPhyt(ff,gg)*LimnoParam%sCPhytW(iLayer,iElem,gg)
                            LimnoParam%oNFoodFiAd = LimnoParam%oNFoodFiAd + LimnoParam%cPrefFiAdPhyt(ff,gg)*LimnoParam%sNPhytW(iLayer,iElem,gg)
                            LimnoParam%oPFoodFiAd = LimnoParam%oPFoodFiAd + LimnoParam%cPrefFiAdPhyt(ff,gg)*LimnoParam%sPPhytW(iLayer,iElem,gg)
                            LimnoParam%oDFoodFiAdPhyt = LimnoParam%oDFoodFiAdPhyt + LimnoParam%cPrefFiAdPhyt(ff,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            aDOMWAd = aDOMWAd + LimnoParam%sDPhytW(iLayer,iElem,gg)
                            !Juveline Fish
                            LimnoParam%oDFoodFiJv = LimnoParam%oDFoodFiJv + LimnoParam%cPrefFiJvPhyt(ff,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            LimnoParam%oCFoodFiJv = LimnoParam%oCFoodFiJv + LimnoParam%cPrefFiJvPhyt(ff,gg)*LimnoParam%sCPhytW(iLayer,iElem,gg)
                            LimnoParam%oNFoodFiJv = LimnoParam%oNFoodFiJv + LimnoParam%cPrefFiJvPhyt(ff,gg)*LimnoParam%sNPhytW(iLayer,iElem,gg)
                            LimnoParam%oPFoodFiJv = LimnoParam%oPFoodFiJv + LimnoParam%cPrefFiJvPhyt(ff,gg)*LimnoParam%sPPhytW(iLayer,iElem,gg)
                            LimnoParam%oDFoodFiJvPhyt = LimnoParam%oDFoodFiJvPhyt + LimnoParam%cPrefFiJvPhyt(ff,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            aDOMWJv = aDOMWJv + LimnoParam%sDPhytW(iLayer,iElem,gg)
                        EndIf
                    EndDo
                    ! POM_for_fish
                    Do pom = 1, LimnoParam%npom
                        If (pom==LimnoParam%LABIL.or.pom==LimnoParam%REFRATARIA) Then
                            !Adult Fish
                            LimnoParam%oDFoodFiAd = LimnoParam%oDFoodFiAd + LimnoParam%cPrefFiAdPom(ff,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            LimnoParam%oCFoodFiAd = LimnoParam%oCFoodFiAd + LimnoParam%cPrefFiAdPom(ff,pom)*LimnoParam%sCPomW(iLayer,iElem,pom)
                            LimnoParam%oNFoodFiAd = LimnoParam%oNFoodFiAd + LimnoParam%cPrefFiAdPom(ff,pom)*LimnoParam%sNPomW(iLayer,iElem,pom)
                            LimnoParam%oPFoodFiAd = LimnoParam%oPFoodFiAd + LimnoParam%cPrefFiAdPom(ff,pom)*LimnoParam%sPPomW(iLayer,iElem,pom)
                            LimnoParam%oDFoodFiAdPom = LimnoParam%oDFoodFiAdPom + LimnoParam%cPrefFiAdPom(ff,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            aDOMWAd = aDOMWAd + LimnoParam%sDPomW(iLayer,iElem,pom)
                            !Juveline Fish
                            LimnoParam%oDFoodFiJv = LimnoParam%oDFoodFiJv + LimnoParam%cPrefFiJvPom(ff,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            LimnoParam%oCFoodFiJv = LimnoParam%oCFoodFiJv + LimnoParam%cPrefFiJvPom(ff,pom)*LimnoParam%sCPomW(iLayer,iElem,pom)
                            LimnoParam%oNFoodFiJv = LimnoParam%oNFoodFiJv + LimnoParam%cPrefFiJvPom(ff,pom)*LimnoParam%sNPomW(iLayer,iElem,pom)
                            LimnoParam%oPFoodFiJv = LimnoParam%oPFoodFiJv + LimnoParam%cPrefFiJvPom(ff,pom)*LimnoParam%sPPomW(iLayer,iElem,pom)
                            LimnoParam%oDFoodFiJvPom = LimnoParam%oDFoodFiJvPom + LimnoParam%cPrefFiJvPom(ff,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            aDOMWJv = aDOMWJv + LimnoParam%sDPomW(iLayer,iElem,pom)
                        EndIf
                    EndDo
                    
                ElseIf  (LimnoParam%afish(ff) == 3) Then !Planktivorous fish (afish = 3)
                    !!  total_food_for_fish
                    LimnoParam%oDFoodFiAd = 0.
                    LimnoParam%oCFoodFiAd = 0.
                    LimnoParam%oNFoodFiAd = 0.
                    LimnoParam%oPFoodFiAd = 0.
                    
                    LimnoParam%oDFoodFiAdZoo = 0.
                    LimnoParam%oDFoodFiJvZoo = 0.
                    LimnoParam%oDFoodFiAdBent = 0.
                    LimnoParam%oDFoodFiJvBent = 0.
                    LimnoParam%oDFoodFiAdPhyt = 0.
                    LimnoParam%oDFoodFiJvPhyt = 0.
                    LimnoParam%oDFoodFiAdPom = 0.
                    LimnoParam%oDFoodFiJvPom = 0.
                    !  organic_seston
                    aDOMWAd = 0.
                    aDOMWJv = 0.
                    ! Zooplankton_for_fish
                    Do zg = 1, LimnoParam%nzoo
                        If (LimnoParam%azoo(zg)>=2) Then ! It can predate meso, macro and general zooplankton
                            !Adult Fish
                            LimnoParam%oDFoodFiAd = LimnoParam%oDFoodFiAd + LimnoParam%cPrefFiAdZoo(ff,zg)*LimnoParam%sDZoo(iLayer,iElem,zg)
                            LimnoParam%oCFoodFiAd = LimnoParam%oCFoodFiAd + LimnoParam%cPrefFiAdZoo(ff,zg)*LimnoParam%sCZoo(iLayer,iElem,zg)
                            LimnoParam%oNFoodFiAd = LimnoParam%oNFoodFiAd + LimnoParam%cPrefFiAdZoo(ff,zg)*LimnoParam%sNZoo(iLayer,iElem,zg)
                            LimnoParam%oPFoodFiAd = LimnoParam%oPFoodFiAd + LimnoParam%cPrefFiAdZoo(ff,zg)*LimnoParam%sPZoo(iLayer,iElem,zg)
                            LimnoParam%oDFoodFiAdZoo = LimnoParam%oDFoodFiAdZoo + LimnoParam%cPrefFiAdZoo(ff,zg)*LimnoParam%sDZoo(iLayer,iElem,zg)
                            aDOMWAd = aDOMWAd + LimnoParam%sDZoo(iLayer,iElem,zg)
                            !Juveline Fish
                            LimnoParam%oDFoodFiJv = LimnoParam%oDFoodFiJv + LimnoParam%cPrefFiJvZoo(ff,zg)*LimnoParam%sDZoo(iLayer,iElem,zg)
                            LimnoParam%oCFoodFiJv = LimnoParam%oCFoodFiJv + LimnoParam%cPrefFiJvZoo(ff,zg)*LimnoParam%sCZoo(iLayer,iElem,zg)
                            LimnoParam%oNFoodFiJv = LimnoParam%oNFoodFiJv + LimnoParam%cPrefFiJvZoo(ff,zg)*LimnoParam%sNZoo(iLayer,iElem,zg)
                            LimnoParam%oPFoodFiJv = LimnoParam%oPFoodFiJv + LimnoParam%cPrefFiJvZoo(ff,zg)*LimnoParam%sPZoo(iLayer,iElem,zg)
                            LimnoParam%oDFoodFiJvZoo = LimnoParam%oDFoodFiJvZoo + LimnoParam%cPrefFiJvZoo(ff,zg)*LimnoParam%sDZoo(iLayer,iElem,zg)
                            aDOMWJv = aDOMWJv + LimnoParam%sDZoo(iLayer,iElem,zg)
                        EndIf
                    EndDo
                    ! Zoobenthos_for_fish
                    Do ben = 1, LimnoParam%nben
                        If (iLayer == HydroParam%ElSmallm(iElem)) Then
                            If (LimnoParam%aben(ben)>=2) Then ! It can predate meso, macro and general zoobenthos
                                !Adult Fish
                                LimnoParam%oDFoodFiAd = LimnoParam%oDFoodFiAd + LimnoParam%cPrefFiAdBent(ff,ben)*LimnoParam%sDBent(iElem,ben)
                                LimnoParam%oCFoodFiAd = LimnoParam%oCFoodFiAd + LimnoParam%cPrefFiAdBent(ff,ben)*LimnoParam%sCBent(iElem,ben)
                                LimnoParam%oNFoodFiAd = LimnoParam%oNFoodFiAd + LimnoParam%cPrefFiAdBent(ff,ben)*LimnoParam%sNBent(iElem,ben)
                                LimnoParam%oPFoodFiAd = LimnoParam%oPFoodFiAd + LimnoParam%cPrefFiAdBent(ff,ben)*LimnoParam%sPBent(iElem,ben)
                                LimnoParam%oDFoodFiAdBent = LimnoParam%oDFoodFiAdBent + LimnoParam%cPrefFiAdBent(ff,ben)*LimnoParam%sDBent(iElem,ben)
                                aDOMWAd = aDOMWAd + LimnoParam%sDBent(iElem,ben)
                                !Juveline Fish
                                LimnoParam%oDFoodFiJv = LimnoParam%oDFoodFiJv + LimnoParam%cPrefFiJvBent(ff,ben)*LimnoParam%sDBent(iElem,ben)
                                LimnoParam%oCFoodFiJv = LimnoParam%oCFoodFiJv + LimnoParam%cPrefFiJvBent(ff,ben)*LimnoParam%sCBent(iElem,ben)
                                LimnoParam%oNFoodFiJv = LimnoParam%oNFoodFiJv + LimnoParam%cPrefFiJvBent(ff,ben)*LimnoParam%sNBent(iElem,ben)
                                LimnoParam%oPFoodFiJv = LimnoParam%oPFoodFiJv + LimnoParam%cPrefFiJvBent(ff,ben)*LimnoParam%sPBent(iElem,ben)
                                LimnoParam%oDFoodFiJvBent = LimnoParam%oDFoodFiJvBent + LimnoParam%cPrefFiJvBent(ff,ben)*LimnoParam%sDBent(iElem,ben)
                                aDOMWJv = aDOMWJv + LimnoParam%sDBent(iElem,zg)
                            EndIf
                        EndIf
                    EndDo
                    ! Phyto_for_fish
                    Do gg = 1, LimnoParam%nphy
                        If (LimnoParam%aphy(gg)==LimnoParam%DINOF.or.LimnoParam%aphy(gg)==LimnoParam%CYANO.or.LimnoParam%aphy(gg)==LimnoParam%NODUL.or.LimnoParam%aphy(gg)==LimnoParam%CHLOR.or.LimnoParam%aphy(gg)==LimnoParam%CRYPT.or.LimnoParam%aphy(gg)==LimnoParam%MDIAT.or.LimnoParam%aphy(gg)==LimnoParam%FDIAT) Then
                            !Adult Fish
                            LimnoParam%oDFoodFiAd = LimnoParam%oDFoodFiAd + LimnoParam%cPrefFiAdPhyt(ff,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            LimnoParam%oCFoodFiAd = LimnoParam%oCFoodFiAd + LimnoParam%cPrefFiAdPhyt(ff,gg)*LimnoParam%sCPhytW(iLayer,iElem,gg)
                            LimnoParam%oNFoodFiAd = LimnoParam%oNFoodFiAd + LimnoParam%cPrefFiAdPhyt(ff,gg)*LimnoParam%sNPhytW(iLayer,iElem,gg)
                            LimnoParam%oPFoodFiAd = LimnoParam%oPFoodFiAd + LimnoParam%cPrefFiAdPhyt(ff,gg)*LimnoParam%sPPhytW(iLayer,iElem,gg)
                            LimnoParam%oDFoodFiAdPhyt = LimnoParam%oDFoodFiAdPhyt + LimnoParam%cPrefFiAdPhyt(ff,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            aDOMWAd = aDOMWAd + LimnoParam%sDPhytW(iLayer,iElem,gg)
                            !Juveline Fish
                            LimnoParam%oDFoodFiJv = LimnoParam%oDFoodFiJv + LimnoParam%cPrefFiJvPhyt(ff,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            LimnoParam%oCFoodFiJv = LimnoParam%oCFoodFiJv + LimnoParam%cPrefFiJvPhyt(ff,gg)*LimnoParam%sCPhytW(iLayer,iElem,gg)
                            LimnoParam%oNFoodFiJv = LimnoParam%oNFoodFiJv + LimnoParam%cPrefFiJvPhyt(ff,gg)*LimnoParam%sNPhytW(iLayer,iElem,gg)
                            LimnoParam%oPFoodFiJv = LimnoParam%oPFoodFiJv + LimnoParam%cPrefFiJvPhyt(ff,gg)*LimnoParam%sPPhytW(iLayer,iElem,gg)
                            LimnoParam%oDFoodFiJvPhyt = LimnoParam%oDFoodFiJvPhyt + LimnoParam%cPrefFiJvPhyt(ff,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg)
                            aDOMWJv = aDOMWJv + LimnoParam%sDPhytW(iLayer,iElem,gg)
                        EndIf
                    EndDo
                    ! POM_for_fish
                    Do pom = 1, LimnoParam%npom
                        If (pom==LimnoParam%LABIL.or.pom==LimnoParam%REFRATARIA) Then
                            !Adult Fish
                            LimnoParam%oDFoodFiAd = LimnoParam%oDFoodFiAd + LimnoParam%cPrefFiAdPom(ff,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            LimnoParam%oCFoodFiAd = LimnoParam%oCFoodFiAd + LimnoParam%cPrefFiAdPom(ff,pom)*LimnoParam%sCPomW(iLayer,iElem,pom)
                            LimnoParam%oNFoodFiAd = LimnoParam%oNFoodFiAd + LimnoParam%cPrefFiAdPom(ff,pom)*LimnoParam%sNPomW(iLayer,iElem,pom)
                            LimnoParam%oPFoodFiAd = LimnoParam%oPFoodFiAd + LimnoParam%cPrefFiAdPom(ff,pom)*LimnoParam%sPPomW(iLayer,iElem,pom)
                            LimnoParam%oDFoodFiAdPom = LimnoParam%oDFoodFiAdPom + LimnoParam%cPrefFiAdPom(ff,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            aDOMWAd = aDOMWAd + LimnoParam%sDPomW(iLayer,iElem,pom)
                            !Juveline Fish
                            LimnoParam%oDFoodFiJv = LimnoParam%oDFoodFiJv + LimnoParam%cPrefFiJvPom(ff,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            LimnoParam%oCFoodFiJv = LimnoParam%oCFoodFiJv + LimnoParam%cPrefFiJvPom(ff,pom)*LimnoParam%sCPomW(iLayer,iElem,pom)
                            LimnoParam%oNFoodFiJv = LimnoParam%oNFoodFiJv + LimnoParam%cPrefFiJvPom(ff,pom)*LimnoParam%sNPomW(iLayer,iElem,pom)
                            LimnoParam%oPFoodFiJv = LimnoParam%oPFoodFiJv + LimnoParam%cPrefFiJvPom(ff,pom)*LimnoParam%sPPomW(iLayer,iElem,pom)
                            LimnoParam%oDFoodFiJvPom = LimnoParam%oDFoodFiJvPom + LimnoParam%cPrefFiJvPom(ff,pom)*LimnoParam%sDPomW(iLayer,iElem,pom)
                            aDOMWJv = aDOMWJv + LimnoParam%sDPomW(iLayer,iElem,pom)
                        EndIf
                    EndDo
                EndIf
                
                !  food_limitation_function_of_adult_fish
                LimnoParam%aDSatFiAd = LimnoParam%oDFoodFiAd /(LimnoParam%hDFiAd(ff) + aDOMWAd)
                !  intrinsic_net_increase_rate_of_fish
                LimnoParam%ukDIncrFiAd = (LimnoParam%kDAssFiAd(ff) - LimnoParam%kDRespFiAd(ff)) * LimnoParam%uFunTmFish - LimnoParam%kMortFiAd(ff)
                !  environmental_correction_of_fish
                LimnoParam%tDEnvFiAd = max(0.0,LimnoParam%ukDIncrFiAd /(LimnoParam%cDCarrFish(ff)-LimnoParam%sDFiJv(iLayer,iElem,ff)) * LimnoParam%sDFiAd(iLayer,iElem,ff)*LimnoParam%sDFiAd(iLayer,iElem,ff))
                !  assimilation_of_fish
                !  PCLake_Osis: 
                LimnoParam%wDAssFiAd(iLayer,iElem,ff) = LimnoParam%aDSatFiAd *(LimnoParam%kDAssFiAd(ff) * LimnoParam%uFunTmFish * LimnoParam%sDFiAd(iLayer,iElem,ff) - LimnoParam%tDEnvFiAd)
                
                LimnoParam%wDConsFiAd = LimnoParam%wDAssFiAd(iLayer,iElem,ff) / LimnoParam%fDAssFiAd(ff)

                !  food_limitation_function_of_young_fish
                LimnoParam%aDSatFiJv = LimnoParam%oDFoodFiJv /(LimnoParam%hDFiJv(ff) + aDOMWJv)
                !  intrinsic_net_increase_rate_of_fish
                LimnoParam%ukDIncrFiJv = (LimnoParam%kDAssFiJv(ff) - LimnoParam%kDRespFiJv(ff)) * LimnoParam%uFunTmFish - LimnoParam%kMortFiJv(ff)
                !  environmental_correction_of_fish
                LimnoParam%tDEnvFiJv = max(0.0,LimnoParam%ukDIncrFiJv /(LimnoParam%cDCarrFish(ff)-LimnoParam%sDFiAd(iLayer,iElem,ff)) * LimnoParam%sDFiJv(iLayer,iElem,ff)*LimnoParam%sDFiJv(iLayer,iElem,ff))
                !  assimilation_of_fish
                !  PCLake_Osis: 
                LimnoParam%wDAssFiJv(iLayer,iElem,ff) = LimnoParam%aDSatFiJv *(LimnoParam%kDAssFiJv(ff) * LimnoParam%uFunTmFish * LimnoParam%sDFiJv(iLayer,iElem,ff) - LimnoParam%tDEnvFiJv)
                
                LimnoParam%wDConsFiJv = LimnoParam%wDAssFiJv(iLayer,iElem,ff) / LimnoParam%fDAssFiJv(ff)
                
                ! DW_fish(fg)_consumption_by_Aduld_fish(ff). Obs. Sum of wDGrz should be wDCons
                Do fg = 1, LimnoParam%nfish
                    LimnoParam%wDGrzFiAdFiAd(iLayer,iElem,ff,fg) = (LimnoParam%cPrefFiAdFiAd(ff,fg)*LimnoParam%sDFiAd(iLayer,iElem,fg)) / LimnoParam%oDFoodFiAd * LimnoParam%wDConsFiAd
                    LimnoParam%wCGrzFiAdFiAd(iLayer,iElem,ff,fg) = LimnoParam%rCDFiAd(iLayer,iElem,fg)*LimnoParam%wDGrzFiAdFiAd(iLayer,iElem,ff,fg)
                    LimnoParam%wNGrzFiAdFiAd(iLayer,iElem,ff,fg) = LimnoParam%rNDFiAd(iLayer,iElem,fg)*LimnoParam%wDGrzFiAdFiAd(iLayer,iElem,ff,fg)
                    LimnoParam%wPGrzFiAdFiAd(iLayer,iElem,ff,fg) = LimnoParam%rPDFiAd(iLayer,iElem,fg)*LimnoParam%wDGrzFiAdFiAd(iLayer,iElem,ff,fg)
                    
                    LimnoParam%wDGrzFiAdFiJv(iLayer,iElem,ff,fg) = (LimnoParam%cPrefFiAdFiJv(ff,fg)*LimnoParam%sDFiJv(iLayer,iElem,fg)) / LimnoParam%oDFoodFiAd * LimnoParam%wDConsFiAd
                    LimnoParam%wCGrzFiAdFiJv(iLayer,iElem,ff,fg) = LimnoParam%rCDFiAd(iLayer,iElem,fg)*LimnoParam%wDGrzFiAdFiJv(iLayer,iElem,ff,fg)
                    LimnoParam%wNGrzFiAdFiJv(iLayer,iElem,ff,fg) = LimnoParam%rNDFiAd(iLayer,iElem,fg)*LimnoParam%wDGrzFiAdFiJv(iLayer,iElem,ff,fg)
                    LimnoParam%wPGrzFiAdFiJv(iLayer,iElem,ff,fg) = LimnoParam%rPDFiAd(iLayer,iElem,fg)*LimnoParam%wDGrzFiAdFiJv(iLayer,iElem,ff,fg)
                EndDo                
               
                !  DW_zooplankton(zg)_consumption_by_Aduld_fish(ff)
                Do zg = 1, LimnoParam%nzoo
                    LimnoParam%wDGrzFiAdZoo(iLayer,iElem,ff,zg) = LimnoParam%cPrefFiAdZoo(ff,zg)*LimnoParam%sDZoo(iLayer,iElem,zg) / LimnoParam%oDFoodFiAd * LimnoParam%wDConsFiAd
                    LimnoParam%wCGrzFiAdZoo(iLayer,iElem,ff,zg) = LimnoParam%rCDZoo(iLayer,iElem,zg)*LimnoParam%wDGrzFiAdZoo(iLayer,iElem,ff,zg)
                    LimnoParam%wNGrzFiAdZoo(iLayer,iElem,ff,zg) = LimnoParam%rNDZoo(iLayer,iElem,zg)*LimnoParam%wDGrzFiAdZoo(iLayer,iElem,ff,zg)
                    LimnoParam%wPGrzFiAdZoo(iLayer,iElem,ff,zg) = LimnoParam%rPDZoo(iLayer,iElem,zg)*LimnoParam%wDGrzFiAdZoo(iLayer,iElem,ff,zg)
                EndDo

                !  DW_zooplankton(zg)_consumption_by_Juveline_fish(ff)
                Do zg = 1, LimnoParam%nzoo
                    LimnoParam%wDGrzFiJvZoo(iLayer,iElem,ff,zg) = LimnoParam%cPrefFiJvZoo(ff,zg)*LimnoParam%sDZoo(iLayer,iElem,zg) / LimnoParam%oDFoodFiJv * LimnoParam%wDConsFiJv
                    LimnoParam%wCGrzFiJvZoo(iLayer,iElem,ff,zg) = LimnoParam%rCDZoo(iLayer,iElem,zg)*LimnoParam%wDGrzFiJvZoo(iLayer,iElem,ff,zg)
                    LimnoParam%wNGrzFiJvZoo(iLayer,iElem,ff,zg) = LimnoParam%rNDZoo(iLayer,iElem,zg)*LimnoParam%wDGrzFiJvZoo(iLayer,iElem,ff,zg)
                    LimnoParam%wPGrzFiJvZoo(iLayer,iElem,ff,zg) = LimnoParam%rPDZoo(iLayer,iElem,zg)*LimnoParam%wDGrzFiJvZoo(iLayer,iElem,ff,zg)
                EndDo
                If (iLayer == HydroParam%ElSmallm(iElem)) Then
                    !  DW_zoobenthos(ben)_consumption_by_Aduld_fish(ff)
                    Do ben = 1, LimnoParam%nben
                        LimnoParam%wDGrzFiAdBent(iLayer,iElem,ff,ben) = LimnoParam%cPrefFiAdBent(ff,ben)*LimnoParam%sDBent(iElem,ben) / LimnoParam%oDFoodFiAd * LimnoParam%wDConsFiAd
                        LimnoParam%wCGrzFiAdBent(iLayer,iElem,ff,ben) = LimnoParam%rCDBent(iElem,ben)*LimnoParam%wDGrzFiAdBent(iLayer,iElem,ff,ben)
                        LimnoParam%wNGrzFiAdBent(iLayer,iElem,ff,ben) = LimnoParam%rNDBent(iElem,ben)*LimnoParam%wDGrzFiAdBent(iLayer,iElem,ff,ben)
                        LimnoParam%wPGrzFiAdBent(iLayer,iElem,ff,ben) = LimnoParam%rPDBent(iElem,ben)*LimnoParam%wDGrzFiAdBent(iLayer,iElem,ff,ben)
                    EndDo

                    !  DW_zoobenthos(ben)_consumption_by_Juveline_fish(ff)
                    Do ben = 1, LimnoParam%nben
                        LimnoParam%wDGrzFiJvBent(iLayer,iElem,ff,ben) = LimnoParam%cPrefFiJvBent(ff,ben)*LimnoParam%sDBent(iElem,zg) / LimnoParam%oDFoodFiJv * LimnoParam%wDConsFiJv
                        LimnoParam%wCGrzFiJvBent(iLayer,iElem,ff,ben) = LimnoParam%rCDBent(iElem,ben)*LimnoParam%wDGrzFiJvBent(iLayer,iElem,ff,zg)
                        LimnoParam%wNGrzFiJvBent(iLayer,iElem,ff,ben) = LimnoParam%rNDBent(iElem,ben)*LimnoParam%wDGrzFiJvBent(iLayer,iElem,ff,zg)
                        LimnoParam%wPGrzFiJvBent(iLayer,iElem,ff,ben) = LimnoParam%rPDBent(iElem,ben)*LimnoParam%wDGrzFiJvBent(iLayer,iElem,ff,zg)
                    EndDo
                Else
                    LimnoParam%wDGrzFiAdBent(iLayer,iElem,ff,:) = 0.
                    LimnoParam%wCGrzFiAdBent(iLayer,iElem,ff,:) = 0.
                    LimnoParam%wNGrzFiAdBent(iLayer,iElem,ff,:) = 0.
                    LimnoParam%wPGrzFiAdBent(iLayer,iElem,ff,:) = 0.
                    
                    LimnoParam%wDGrzFiJvBent(iLayer,iElem,ff,:) = 0.
                    LimnoParam%wCGrzFiJvBent(iLayer,iElem,ff,:) = 0.
                    LimnoParam%wNGrzFiJvBent(iLayer,iElem,ff,:) = 0.
                    LimnoParam%wPGrzFiJvBent(iLayer,iElem,ff,:) = 0.
                EndIf
                
                !  DW_phytoplankton(gg)_consumption_by_Aduld_fish(ff)
                Do gg = 1, LimnoParam%nphy
                    LimnoParam%wDGrzFiAdPhyt(iLayer,iElem,ff,gg) = LimnoParam%cPrefFiAdPhyt(ff,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg) / LimnoParam%oDFoodFiAd * LimnoParam%wDConsFiAd
                    LimnoParam%wCGrzFiAdPhyt(iLayer,iElem,ff,gg) = LimnoParam%rCDPhytW(iLayer,iElem,gg)*LimnoParam%wDGrzFiAdPhyt(iLayer,iElem,ff,gg)
                    LimnoParam%wNGrzFiAdPhyt(iLayer,iElem,ff,gg) = LimnoParam%rNDPhytW(iLayer,iElem,gg)*LimnoParam%wDGrzFiAdPhyt(iLayer,iElem,ff,gg)
                    LimnoParam%wPGrzFiAdPhyt(iLayer,iElem,ff,gg) = LimnoParam%rPDPhytW(iLayer,iElem,gg)*LimnoParam%wDGrzFiAdPhyt(iLayer,iElem,ff,gg)
                EndDo

                !  DW_phytoplankton(gg)_consumption_by_Juveline_fish(ff)
                Do gg = 1, LimnoParam%nphy
                    LimnoParam%wDGrzFiJvPhyt(iLayer,iElem,ff,gg) = LimnoParam%cPrefFiJvPhyt(ff,gg)*LimnoParam%sDPhytW(iLayer,iElem,gg) / LimnoParam%oDFoodFiJv * LimnoParam%wDConsFiJv
                    LimnoParam%wCGrzFiJvPhyt(iLayer,iElem,ff,gg) = LimnoParam%rCDPhytW(iLayer,iElem,gg)*LimnoParam%wDGrzFiJvPhyt(iLayer,iElem,ff,gg)
                    LimnoParam%wNGrzFiJvPhyt(iLayer,iElem,ff,gg) = LimnoParam%rNDPhytW(iLayer,iElem,gg)*LimnoParam%wDGrzFiJvPhyt(iLayer,iElem,ff,gg)
                    LimnoParam%wPGrzFiJvPhyt(iLayer,iElem,ff,gg) = LimnoParam%rPDPhytW(iLayer,iElem,gg)*LimnoParam%wDGrzFiJvPhyt(iLayer,iElem,ff,gg)
                EndDo
                
                !  DW_detritus_consumption_by_Aduld_fish(ff)
                Do pom = 1, LimnoParam%npom
                    LimnoParam%wDGrzFiAdPom(iLayer,iElem,ff,pom) = LimnoParam%cPrefFiAdPom(ff,pom)*LimnoParam%sDPomW(iLayer,iElem,pom) / LimnoParam%oDFoodFiAd * LimnoParam%wDConsFiAd
                    LimnoParam%wCGrzFiAdPom(iLayer,iElem,ff,pom) = LimnoParam%rCDPomW(iLayer,iElem,pom)*LimnoParam%wDGrzFiAdPom(iLayer,iElem,ff,pom)
                    LimnoParam%wNGrzFiAdPom(iLayer,iElem,ff,pom) = LimnoParam%rNDPomW(iLayer,iElem,pom)*LimnoParam%wDGrzFiAdPom(iLayer,iElem,ff,pom)
                    LimnoParam%wPGrzFiAdPom(iLayer,iElem,ff,pom) = LimnoParam%rPDPomW(iLayer,iElem,pom)*LimnoParam%wDGrzFiAdPom(iLayer,iElem,ff,pom)
                EndDo

                !  DW_detritus_consumption_by_Juveline_fish(ff)
                Do pom = 1, LimnoParam%npom
                    LimnoParam%wDGrzFiJvPom(iLayer,iElem,ff,pom) = LimnoParam%cPrefFiJvPom(ff,pom)*LimnoParam%sDPomW(iLayer,iElem,pom) / LimnoParam%oDFoodFiJv * LimnoParam%wDConsFiJv
                    LimnoParam%wCGrzFiJvPom(iLayer,iElem,ff,pom) = LimnoParam%rCDPomW(iLayer,iElem,pom)*LimnoParam%wDGrzFiJvPom(iLayer,iElem,ff,pom)
                    LimnoParam%wNGrzFiJvPom(iLayer,iElem,ff,pom) = LimnoParam%rNDPomW(iLayer,iElem,pom)*LimnoParam%wDGrzFiJvPom(iLayer,iElem,ff,pom)
                    LimnoParam%wPGrzFiJvPom(iLayer,iElem,ff,pom) = LimnoParam%rPDPomW(iLayer,iElem,pom)*LimnoParam%wDGrzFiJvPom(iLayer,iElem,ff,pom)
                EndDo
                
                !-----------------------------------------------------------------------
                !  Aduld Fish assimilation C
                !-----------------------------------------------------------------------
                !  C/DW_ratio_of_adult_fish_food
                LimnoParam%rCDFoodFiAd = LimnoParam%oCFoodFiAd /(LimnoParam%oDFoodFiAd+NearZero)
                !  total_C_consumption
                LimnoParam%wCConsFiAd = (sum(LimnoParam%wCGrzFiAdFiAd(iLayer,iElem,ff,:)) + sum(LimnoParam%wCGrzFiAdFiJv(iLayer,iElem,ff,:)) + sum(LimnoParam%wCGrzFiAdZoo(iLayer,iElem,ff,:)) + sum(LimnoParam%wCGrzFiAdBent(iLayer,iElem,ff,:)) + sum(LimnoParam%wCGrzFiAdPhyt(iLayer,iElem,ff,:)) + sum(LimnoParam%wCGrzFiAdPom(iLayer,iElem,ff,:)))
                !  C_assimilation_efficiency_of_fish
                LimnoParam%afCAssFiAd = min(1.0,LimnoParam%cCDFishRef(ff) / (LimnoParam%rCDFoodFiAd+NearZero) * LimnoParam%fDAssFiAd(ff))
                !  assimilation_by_fish
                LimnoParam%wCAssFiAd(iLayer,iElem,ff) = LimnoParam%afCAssFiAd*LimnoParam%wCConsFiAd

                !-----------------------------------------------------------------------
                !  Aduld Fish assimilation N
                !-----------------------------------------------------------------------
                !  N/DW_ratio_of_adult_fish_food
                LimnoParam%rNDFoodFiAd = LimnoParam%oNFoodFiAd /(LimnoParam%oDFoodFiAd+NearZero)
                !  total_C_consumption
                LimnoParam%wNConsFiAd = (sum(LimnoParam%wNGrzFiAdFiAd(iLayer,iElem,ff,:)) + sum(LimnoParam%wNGrzFiAdFiJv(iLayer,iElem,ff,:)) + sum(LimnoParam%wNGrzFiAdZoo(iLayer,iElem,ff,:)) + sum(LimnoParam%wNGrzFiAdBent(iLayer,iElem,ff,:)) + sum(LimnoParam%wNGrzFiAdPhyt(iLayer,iElem,ff,:)) + sum(LimnoParam%wNGrzFiAdPom(iLayer,iElem,ff,:)))
                !  C_assimilation_efficiency_of_fish
                LimnoParam%afNAssFiAd = min(1.0,LimnoParam%cNDFishRef(ff) / (LimnoParam%rNDFoodFiAd+NearZero) * LimnoParam%fDAssFiAd(ff))
                !  assimilation_by_fish
                LimnoParam%wNAssFiAd(iLayer,iElem,ff) = LimnoParam%afNAssFiAd*LimnoParam%wNConsFiAd

                !-----------------------------------------------------------------------
                !  Aduld Fish assimilation P
                !-----------------------------------------------------------------------
                !  P/DW_ratio_of_adult_fish_food
                LimnoParam%rPDFoodFiAd = LimnoParam%oPFoodFiAd /(LimnoParam%oDFoodFiAd+NearZero)
                !  total_C_consumption
                LimnoParam%wPConsFiAd = (sum(LimnoParam%wPGrzFiAdFiAd(iLayer,iElem,ff,:)) + sum(LimnoParam%wPGrzFiAdFiJv(iLayer,iElem,ff,:)) + sum(LimnoParam%wPGrzFiAdZoo(iLayer,iElem,ff,:)) + sum(LimnoParam%wPGrzFiAdBent(iLayer,iElem,ff,:)) + sum(LimnoParam%wPGrzFiAdPhyt(iLayer,iElem,ff,:)) + sum(LimnoParam%wPGrzFiAdPom(iLayer,iElem,ff,:)))
                !  C_assimilation_efficiency_of_fish
                LimnoParam%afPAssFiAd = min(1.0,LimnoParam%cPDFishRef(ff) / (LimnoParam%rPDFoodFiAd+NearZero) * LimnoParam%fDAssFiAd(ff))
                !  assimilation_by_fish
                LimnoParam%wPAssFiAd(iLayer,iElem,ff) = LimnoParam%afPAssFiAd*LimnoParam%wPConsFiAd
                
                !-----------------------------------------------------------------------
                !  Juvenile Fish assimilation C
                !-----------------------------------------------------------------------
                !  C/DW_ratio_of_Juveline_fish_food
                LimnoParam%rCDFoodFiJv = LimnoParam%oCFoodFiJv /(LimnoParam%oDFoodFiJv+NearZero)
                !  total_C_consumption
                LimnoParam%wCConsFiJv = (sum(LimnoParam%wCGrzFiJvZoo(iLayer,iElem,ff,:)) + sum(LimnoParam%wCGrzFiJvBent(iLayer,iElem,ff,:)) + sum(LimnoParam%wCGrzFiJvPhyt(iLayer,iElem,ff,:)) + sum(LimnoParam%wCGrzFiJvPom(iLayer,iElem,ff,:)))
                !  C_assimilation_efficiency_of_fish
                LimnoParam%afCAssFiJv = min(1.0,LimnoParam%cCDFishRef(ff) / (LimnoParam%rCDFoodFiJv+NearZero) * LimnoParam%fDAssFiJv(ff))
                !  assimilation_by_fish
                LimnoParam%wCAssFiJv(iLayer,iElem,ff) = LimnoParam%afCAssFiJv*LimnoParam%wCConsFiJv

                !-----------------------------------------------------------------------
                !  Juvenile Fish assimilation N
                !-----------------------------------------------------------------------
                !  C/DW_ratio_of_Juveline_fish_food
                LimnoParam%rNDFoodFiJv = LimnoParam%oNFoodFiJv /(LimnoParam%oDFoodFiJv+NearZero)
                !  total_C_consumption
                LimnoParam%wNConsFiJv = (sum(LimnoParam%wNGrzFiJvZoo(iLayer,iElem,ff,:)) + sum(LimnoParam%wNGrzFiJvBent(iLayer,iElem,ff,:)) + sum(LimnoParam%wNGrzFiJvPhyt(iLayer,iElem,ff,:)) + sum(LimnoParam%wNGrzFiJvPom(iLayer,iElem,ff,:)))
                !  C_assimilation_efficiency_of_fish
                LimnoParam%afNAssFiJv = min(1.0,LimnoParam%cNDFishRef(ff) / (LimnoParam%rNDFoodFiJv+NearZero) * LimnoParam%fDAssFiJv(ff))
                !  assimilation_by_fish
                LimnoParam%wNAssFiJv(iLayer,iElem,ff) = LimnoParam%afNAssFiJv*LimnoParam%wNConsFiJv

                !-----------------------------------------------------------------------
                !  Juvenile Fish assimilation P
                !-----------------------------------------------------------------------
                !  P/DW_ratio_of_Juveline_fish_food
                LimnoParam%rPDFoodFiJv = LimnoParam%oPFoodFiJv /(LimnoParam%oDFoodFiJv+NearZero)
                !  total_C_consumption
                LimnoParam%wPConsFiJv = (sum(LimnoParam%wPGrzFiJvZoo(iLayer,iElem,ff,:)) + sum(LimnoParam%wPGrzFiJvBent(iLayer,iElem,ff,:)) + sum(LimnoParam%wPGrzFiJvPhyt(iLayer,iElem,ff,:)) + sum(LimnoParam%wPGrzFiJvPom(iLayer,iElem,ff,:)))
                !  C_assimilation_efficiency_of_fish
                LimnoParam%afPAssFiJv = min(1.0,LimnoParam%cPDFishRef(ff) / (LimnoParam%rPDFoodFiJv+NearZero) * LimnoParam%fDAssFiJv(ff))
                !  assimilation_by_fish
                LimnoParam%wPAssFiJv(iLayer,iElem,ff) = LimnoParam%afPAssFiJv*LimnoParam%wPConsFiJv
                
                !5) Messy Fedding (-) -> POM 
                !-----------------------------------------------------------------------
                !  adult fish Messy Fedding
                !-----------------------------------------------------------------------
                 LimnoParam%wDMessFiAdOM(iLayer,iElem,ff) = LimnoParam%wDConsFiAd - LimnoParam%wDAssFiAd(iLayer,iElem,ff)
                 LimnoParam%wCMessFiAdOM(iLayer,iElem,ff) = LimnoParam%wDMessFiAdOM(iLayer,iElem,ff) * LimnoParam%wCConsFiAd/(LimnoParam%wDConsFiAd+NearZero)
                 LimnoParam%wNMessFiAdOM(iLayer,iElem,ff) = LimnoParam%wDMessFiAdOM(iLayer,iElem,ff) * LimnoParam%wNConsFiAd/(LimnoParam%wDConsFiAd+NearZero)
                 LimnoParam%wPMessFiAdOM(iLayer,iElem,ff) = LimnoParam%wDMessFiAdOM(iLayer,iElem,ff) * LimnoParam%wPConsFiAd/(LimnoParam%wDConsFiAd+NearZero)
                !-----------------------------------------------------------------------
                !  juvenile fish Messy Fedding
                !-----------------------------------------------------------------------
                If (LimnoParam%afish(ff) == 1) Then !Piscivorous fish (afish = 1)
                    LimnoParam%wDMessFiJvOM(iLayer,iElem,ff) = 0.
                Else
                    LimnoParam%wDMessFiJvOM(iLayer,iElem,ff) = LimnoParam%wDConsFiJv - LimnoParam%wDAssFiJv(iLayer,iElem,ff)
                EndIf
                 LimnoParam%wCMessFiJvOM(iLayer,iElem,ff) = LimnoParam%wDMessFiJvOM(iLayer,iElem,ff) * LimnoParam%wCConsFiJv/(LimnoParam%wDConsFiJv+NearZero)
                 LimnoParam%wNMessFiJvOM(iLayer,iElem,ff) = LimnoParam%wDMessFiJvOM(iLayer,iElem,ff) * LimnoParam%wNConsFiJv/(LimnoParam%wDConsFiJv+NearZero)
                 LimnoParam%wPMessFiJvOM(iLayer,iElem,ff) = LimnoParam%wDMessFiJvOM(iLayer,iElem,ff) * LimnoParam%wPConsFiJv/(LimnoParam%wDConsFiJv+NearZero)
                
                
                !6) Respiration/Excretion(-) -> CO2,N,P 
                !-----------------------------------------------------------------------
                !  adult fish respiration and excretion
                !-----------------------------------------------------------------------
                !  respiration_of_fish
                !  PCLake_osis:sDFiAd,in g/m^2
                LimnoParam%wDRespFiAd(iLayer,iElem,ff) = LimnoParam%kDRespFiAd(ff) * LimnoParam%uFunTmFish * LimnoParam%sDFiAd(iLayer,iElem,ff)
                !  C_respiration_of_FiAd
                !  PCLake_osis:sCFiAd,in g/m^2
                LimnoParam%wCRespFiAd(iLayer,iElem,ff) = (LimnoParam%rCDFiAd(iLayer,iElem,ff) / (LimnoParam%cCDFishRef(ff)+NearZero)) * LimnoParam%kDRespFiAd(ff) * LimnoParam%uFunTmFish * LimnoParam%sCFiAd(iLayer,iElem,ff)
                !  N_excretion_of_FiAd
                !  PCLake_osis:sNFiAd,in g/m^2
                LimnoParam%wNExcrFiAd(iLayer,iElem,ff) = (LimnoParam%rNDFiAd(iLayer,iElem,ff) / (LimnoParam%cNDFishRef(ff)+NearZero)) * LimnoParam%kDRespFiAd(ff) * LimnoParam%uFunTmFish * LimnoParam%sNFiAd(iLayer,iElem,ff)       
                !  P_excretion_of_FiAd
                !  PCLake_osis:sPFiAd,in g/m^2
                LimnoParam%wPExcrFiAd(iLayer,iElem,ff) = (LimnoParam%rPDFiAd(iLayer,iElem,ff) / (LimnoParam%cPDFishRef(ff)+NearZero)) * LimnoParam%kDRespFiAd(ff) * LimnoParam%uFunTmFish * LimnoParam%sPFiAd(iLayer,iElem,ff)
                !-----------------------------------------------------------------------
                !  young fish respiration and excretion
                !-----------------------------------------------------------------------
                !  respiration_of_fish_DW
                !  PCLake_osis:sDFiAd,in g/m^2
                If (LimnoParam%afish(ff) == 1) Then !Piscivorous fish (afish = 1)
                    LimnoParam%wDRespFiJv(iLayer,iElem,ff) = 0.
                    !  C_respiration_of_FiJv
                    !  PCLake_osis:sCFiJv,in g/m^2
                    LimnoParam%wCRespFiJv(iLayer,iElem,ff) = 0.
                    !  N_excretion_of_FiJv
                    !  PCLake_osis:sNFiAd,in g/m^2
                    LimnoParam%wNExcrFiJv(iLayer,iElem,ff) = 0.   
                    !  P_excretion_of_FiJv
                    !  PCLake_osis:sPFiAd,in g/m^2
                    LimnoParam%wPExcrFiJv(iLayer,iElem,ff) = 0.
                Else
                    LimnoParam%wDRespFiJv(iLayer,iElem,ff) = LimnoParam%kDRespFiJv(ff) * LimnoParam%uFunTmFish * LimnoParam%sDFiJv(iLayer,iElem,ff)
                    !  C_respiration_of_FiJv
                    !  PCLake_osis:sCFiJv,in g/m^2
                    LimnoParam%wCRespFiJv(iLayer,iElem,ff) = (LimnoParam%rCDFiJv(iLayer,iElem,ff) / (LimnoParam%cCDFishRef(ff)+NearZero)) * LimnoParam%kDRespFiJv(ff) * LimnoParam%uFunTmFish * LimnoParam%sCFiJv(iLayer,iElem,ff)
                    !  N_excretion_of_FiJv
                    !  PCLake_osis:sNFiAd,in g/m^2
                    LimnoParam%wNExcrFiJv(iLayer,iElem,ff) = (LimnoParam%rNDFiJv(iLayer,iElem,ff) / (LimnoParam%cNDFishRef(ff)+NearZero)) * LimnoParam%kDRespFiJv(ff) * LimnoParam%uFunTmFish * LimnoParam%sNFiJv(iLayer,iElem,ff)       
                    !  P_excretion_of_FiJv
                    !  PCLake_osis:sPFiAd,in g/m^2
                    LimnoParam%wPExcrFiJv(iLayer,iElem,ff) = (LimnoParam%rPDFiJv(iLayer,iElem,ff) / (LimnoParam%cPDFishRef(ff)+NearZero)) * LimnoParam%kDRespFiJv(ff) * LimnoParam%uFunTmFish * LimnoParam%sPFiJv(iLayer,iElem,ff)
                EndIf
  
                !7) Pellets(-) -> POM
                !-----------------------------------------------------------------------
                !  adult fish Pellets
                !-----------------------------------------------------------------------
				LimnoParam%wDEgesFiAd(iLayer,iElem,ff) = LimnoParam%kPelFiAd(ff) * LimnoParam%sDFiAd(iLayer,iElem,ff) 
			    LimnoParam%wCEgesFiAd(iLayer,iElem,ff) = LimnoParam%rCDFiAd(iLayer,iElem,ff) * LimnoParam%wDEgesFiAd(iLayer,iElem,ff)
			    LimnoParam%wNEgesFiAd(iLayer,iElem,ff) = LimnoParam%rNDFiAd(iLayer,iElem,ff) * LimnoParam%wDEgesFiAd(iLayer,iElem,ff)
			    LimnoParam%wPEgesFiAd(iLayer,iElem,ff) = LimnoParam%rPDFiAd(iLayer,iElem,ff) * LimnoParam%wDEgesFiAd(iLayer,iElem,ff)
                !-----------------------------------------------------------------------
                !  juvenile fish Pellets
                !-----------------------------------------------------------------------
                If (LimnoParam%afish(ff) == 1) Then !Piscivorous fish (afish = 1)
                    LimnoParam%wDEgesFiJv(iLayer,iElem,ff) = 0.
                Else
                    LimnoParam%wDEgesFiJv(iLayer,iElem,ff) = LimnoParam%kPelFiJv(ff) * LimnoParam%sDFiJv(iLayer,iElem,ff) 
                EndIf
			    LimnoParam%wCEgesFiJv(iLayer,iElem,ff) = LimnoParam%rCDFiJv(iLayer,iElem,ff) * LimnoParam%wDEgesFiJv(iLayer,iElem,ff)
			    LimnoParam%wNEgesFiJv(iLayer,iElem,ff) = LimnoParam%rNDFiJv(iLayer,iElem,ff) * LimnoParam%wDEgesFiJv(iLayer,iElem,ff)
			    LimnoParam%wPEgesFiJv(iLayer,iElem,ff) = LimnoParam%rPDFiJv(iLayer,iElem,ff) * LimnoParam%wDEgesFiJv(iLayer,iElem,ff)
                
                !8) Mortality -> POM 
                !-----------------------------------------------------------------------
                !  adult fish mortality
                !-----------------------------------------------------------------------
                !  fish_mortality_incl._environmental_correction
                !  PCLake_osis:sDFiAd,in g/m^2
                LimnoParam%wDMortFiAd(iLayer,iElem,ff) = LimnoParam%kMortFiAd(ff) * LimnoParam%sDFiAd(iLayer,iElem,ff) +(1.0 - LimnoParam%aDSatFiAd) * LimnoParam%tDEnvFiAd
                !  mortality_of_FiAd
                !  PCLake_osis:sPFiAd,in g/m^2
                LimnoParam%wCMortFiAd(iLayer,iElem,ff) = LimnoParam%rCDFiAd(iLayer,iElem,ff) * LimnoParam%wDMortFiAd(iLayer,iElem,ff)
                !  PCLake_osis:sNFiAd,in g/m^2
                LimnoParam%wNMortFiAd(iLayer,iElem,ff) = LimnoParam%rNDFiAd(iLayer,iElem,ff) * LimnoParam%wDMortFiAd(iLayer,iElem,ff)
                !  PCLake_osis:sPFiAd,in g/m^2
                LimnoParam%wPMortFiAd(iLayer,iElem,ff) = LimnoParam%rPDFiAd(iLayer,iElem,ff) * LimnoParam%wDMortFiAd(iLayer,iElem,ff)
                !-----------------------------------------------------------------------
                !  young fish mortality
                !-----------------------------------------------------------------------
                !  fish_mortality_incl._environmental_correction
                !  PCLake_osis:sDFiAd,in g/m^2
                If (LimnoParam%afish(ff) == 1) Then !Piscivorous fish (afish = 1)
                    LimnoParam%wDMortFiJv(iLayer,iElem,ff) = 0.
                Else
                    LimnoParam%wDMortFiJv(iLayer,iElem,ff) = LimnoParam%kMortFiJv(ff) * LimnoParam%sDFiJv(iLayer,iElem,ff) +(1.0 - LimnoParam%aDSatFiJv) * LimnoParam%tDEnvFiJv
                EndIf
                !  mortality_of_FiAd
                !  PCLake_osis:sPFiAd,in g/m^2
                LimnoParam%wCMortFiJv(iLayer,iElem,ff) = LimnoParam%rCDFiJv(iLayer,iElem,ff) * LimnoParam%wDMortFiJv(iLayer,iElem,ff)
                !  PCLake_osis:sNFiAd,in g/m^2
                LimnoParam%wNMortFiJv(iLayer,iElem,ff) = LimnoParam%rNDFiJv(iLayer,iElem,ff) * LimnoParam%wDMortFiJv(iLayer,iElem,ff)
                !  PCLake_osis:sPFiAd,in g/m^2
                LimnoParam%wPMortFiJv(iLayer,iElem,ff) = LimnoParam%rPDFiJv(iLayer,iElem,ff) * LimnoParam%wDMortFiJv(iLayer,iElem,ff)
                   
                   
                !-------------------------------------------------------------------------------------!
                !                               Source/Sink Term for adult fish                       !
                !-------------------------------------------------------------------------------------!
                
			    HydroParam%dVarEst(iLayer,iElem,1) =   -LimnoParam%wDReprFish(ff)	                            & !Reprodution
				                            +LimnoParam%wDMigrFiAd(ff)	                            & !Migration
				                            +LimnoParam%wDAgeFish(ff)  	                        & !Aging JV -> AD
				                            -LimnoParam%wDHarvFiAd(ff)                             & !Harvest
                                            +LimnoParam%wDAssFiAd(iLayer,iElem,ff)                 & !Assimilation
						                    -LimnoParam%wDRespFiAd(iLayer,iElem,ff)                & !Respiration/Excretion
						                    -LimnoParam%wDEgesFiAd(iLayer,iElem,ff)                & !Egestion/Pellets 
						                    -LimnoParam%wDMortFiAd(iLayer,iElem,ff)                & !Mortality
				                            -SUM(LimnoParam%wDGrzFiAdFiAd(iLayer,iElem,:,ff))      !Predation

			    HydroParam%dVarEst(iLayer,iElem,2) =   -LimnoParam%wCReprFish(ff)	                            & !Reprodution
				                            +LimnoParam%wCMigrFiAd(ff)	                            & !Migration
				                            +LimnoParam%wCAgeFish(ff)  	                        & !Aging JV -> AD
				                            -LimnoParam%wCHarvFiAd(ff)                             & !Harvest
                                            +LimnoParam%wCAssFiAd(iLayer,iElem,ff)                 & !Assimilation
						                    -LimnoParam%wCRespFiAd(iLayer,iElem,ff)                & !Respiration/Excretion
						                    -LimnoParam%wCEgesFiAd(iLayer,iElem,ff)                & !Egestion/Pellets 
						                    -LimnoParam%wCMortFiAd(iLayer,iElem,ff)                & !Mortality
				                            -SUM(LimnoParam%wCGrzFiAdFiAd(iLayer,iElem,:,ff))      !Predation
                
			    HydroParam%dVarEst(iLayer,iElem,3) =   -LimnoParam%wNReprFish(ff)	                            & !Reprodution
				                            +LimnoParam%wNMigrFiAd(ff)	                            & !Migration
				                            +LimnoParam%wNAgeFish(ff)  	                        & !Aging JV -> AD
				                            -LimnoParam%wNHarvFiAd(ff)                             & !Harvest
                                            +LimnoParam%wNAssFiAd(iLayer,iElem,ff)                 & !Assimilation
						                    -LimnoParam%wNExcrFiAd(iLayer,iElem,ff)                & !Respiration/Excretion
						                    -LimnoParam%wNEgesFiAd(iLayer,iElem,ff)                & !Egestion/Pellets 
						                    -LimnoParam%wNMortFiAd(iLayer,iElem,ff)                & !Mortality
				                            -SUM(LimnoParam%wNGrzFiAdFiAd(iLayer,iElem,:,ff))      !Predation

			    HydroParam%dVarEst(iLayer,iElem,4) =   -LimnoParam%wPReprFish(ff)	                            & !Reprodution
				                            +LimnoParam%wPMigrFiAd(ff)	                            & !Migration
				                            +LimnoParam%wPAgeFish(ff)  	                        & !Aging JV -> AD
				                            -LimnoParam%wPHarvFiAd(ff)                             & !Harvest
                                            +LimnoParam%wPAssFiAd(iLayer,iElem,ff)                 & !Assimilation
						                    -LimnoParam%wPExcrFiAd(iLayer,iElem,ff)                & !Respiration/Excretion
						                    -LimnoParam%wPEgesFiAd(iLayer,iElem,ff)                & !Egestion/Pellets 
						                    -LimnoParam%wPMortFiAd(iLayer,iElem,ff)                & !Mortality
				                            -SUM(LimnoParam%wPGrzFiAdFiAd(iLayer,iElem,:,ff))      !Predation
                
                LimnoParam%sDFiAd(iLayer,iElem,ff) = MAX(NearZero,(dtday*HydroParam%dVarEst(iLayer,iElem,1)+LimnoParam%sDFiAdP(iLayer,iElem,ff)))
                LimnoParam%sCFiAd(iLayer,iElem,ff) = MAX(NearZero,(dtday*HydroParam%dVarEst(iLayer,iElem,2)+LimnoParam%sCFiAdP(iLayer,iElem,ff)))
                LimnoParam%sNFiAd(iLayer,iElem,ff) = MAX(NearZero,(dtday*HydroParam%dVarEst(iLayer,iElem,3)+LimnoParam%sNFiAdP(iLayer,iElem,ff)))
                LimnoParam%sPFiAd(iLayer,iElem,ff) = MAX(NearZero,(dtday*HydroParam%dVarEst(iLayer,iElem,4)+LimnoParam%sPFiAdP(iLayer,iElem,ff)))

                !-------------------------------------------------------------------------------------!
                !                               Source/Sink Term for juvenile fish                       !
                !-------------------------------------------------------------------------------------!
                
			    HydroParam%dVarEst(iLayer,iElem,1) =   LimnoParam%wDReprFish(ff)	                            & !Reprodution
				                            +LimnoParam%wDMigrFiJv(ff)	                            & !Migration
				                            -LimnoParam%wDAgeFish(ff)  	                        & !Aging JV -> AD
				                            -LimnoParam%wDHarvFiJv(ff)                             & !Harvest
                                            +LimnoParam%wDAssFiJv(iLayer,iElem,ff)                 & !Assimilation
						                    -LimnoParam%wDRespFiJv(iLayer,iElem,ff)                & !Respiration/Excretion
						                    -LimnoParam%wDEgesFiJv(iLayer,iElem,ff)                & !Egestion/Pellets 
						                    -LimnoParam%wDMortFiJv(iLayer,iElem,ff)                & !Mortality
				                            -SUM(LimnoParam%wDGrzFiAdFiJv(iLayer,iElem,:,ff))      !Predation

			    HydroParam%dVarEst(iLayer,iElem,2) =   LimnoParam%wCReprFish(ff)	                            & !Reprodution
				                            +LimnoParam%wCMigrFiJv(ff)	                            & !Migration
				                            -LimnoParam%wCAgeFish(ff)  	                        & !Aging JV -> AD
				                            -LimnoParam%wCHarvFiJv(ff)                             & !Harvest
                                            +LimnoParam%wCAssFiJv(iLayer,iElem,ff)                 & !Assimilation
						                    -LimnoParam%wCRespFiJv(iLayer,iElem,ff)                & !Respiration/Excretion
						                    -LimnoParam%wCEgesFiJv(iLayer,iElem,ff)                & !Egestion/Pellets 
						                    -LimnoParam%wCMortFiJv(iLayer,iElem,ff)                & !Mortality
				                            -SUM(LimnoParam%wCGrzFiAdFiJv(iLayer,iElem,:,ff))      !Predation
                
			    HydroParam%dVarEst(iLayer,iElem,3) =   LimnoParam%wNReprFish(ff)	                            & !Reprodution
				                            +LimnoParam%wNMigrFiJv(ff)	                            & !Migration
				                            -LimnoParam%wNAgeFish(ff)  	                        & !Aging JV -> AD
				                            -LimnoParam%wNHarvFiJv(ff)                             & !Harvest
                                            +LimnoParam%wNAssFiJv(iLayer,iElem,ff)                 & !Assimilation
						                    -LimnoParam%wNExcrFiJv(iLayer,iElem,ff)                & !Respiration/Excretion
						                    -LimnoParam%wNEgesFiJv(iLayer,iElem,ff)                & !Egestion/Pellets 
						                    -LimnoParam%wNMortFiJv(iLayer,iElem,ff)                & !Mortality
				                            -SUM(LimnoParam%wNGrzFiAdFiJv(iLayer,iElem,:,ff))      !Predation

			    HydroParam%dVarEst(iLayer,iElem,4) =   LimnoParam%wPReprFish(ff)	                            & !Reprodution
				                            +LimnoParam%wPMigrFiJv(ff)	                            & !Migration
				                            -LimnoParam%wPAgeFish(ff)  	                        & !Aging JV -> AD
				                            -LimnoParam%wPHarvFiJv(ff)                             & !Harvest
                                            +LimnoParam%wPAssFiJv(iLayer,iElem,ff)                 & !Assimilation
						                    -LimnoParam%wPExcrFiJv(iLayer,iElem,ff)                & !Respiration/Excretion
						                    -LimnoParam%wPEgesFiJv(iLayer,iElem,ff)                & !Egestion/Pellets 
						                    -LimnoParam%wPMortFiJv(iLayer,iElem,ff)                & !Mortality
				                            -SUM(LimnoParam%wPGrzFiAdFiJv(iLayer,iElem,:,ff))      !Predation
                
                LimnoParam%sDFiJv(iLayer,iElem,ff) = MAX(NearZero,(dtday*HydroParam%dVarEst(iLayer,iElem,1)+LimnoParam%sDFiJvP(iLayer,iElem,ff)))
                LimnoParam%sCFiJv(iLayer,iElem,ff) = MAX(NearZero,(dtday*HydroParam%dVarEst(iLayer,iElem,2)+LimnoParam%sCFiJvP(iLayer,iElem,ff)))
                LimnoParam%sNFiJv(iLayer,iElem,ff) = MAX(NearZero,(dtday*HydroParam%dVarEst(iLayer,iElem,3)+LimnoParam%sNFiJvP(iLayer,iElem,ff)))
                LimnoParam%sPFiJv(iLayer,iElem,ff) = MAX(NearZero,(dtday*HydroParam%dVarEst(iLayer,iElem,4)+LimnoParam%sPFiJvP(iLayer,iElem,ff)))
                
            EndDo !Loop Layer

        EndDo !Loop Cell
    
    EndDo !Loop Fishes
      
    !-------------------------------------------------------------------------------------!
    !                       Solver Transport for Fishes                    	    	      !
    !-------------------------------------------------------------------------------------!
    ! Call MovFishes
    
Return
End
