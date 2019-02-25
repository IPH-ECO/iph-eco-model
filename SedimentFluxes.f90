SUBROUTINE SedimentFluxes(HydroParam,MeshParam,LimnoParam,dt,dtday)
! !DESCRIPTION:
!  SedimentFluxes is created for the purpose of computing resuspension and sedimentation,
!  sediment burial processes. 
!  No local state variable is registed here.
!  resuspension and sedimentation involve: sDIMW<==>sDIMS,sD/N/PDetW<==>sD/N/PDetS,sPAIMW<==>sPAIMS
!  aD/N/PPhytW<==>aD/N/PPhytS,sDDiatW<==>sDDiatS,sNH4S==>sNH4W,sNO3S==>sNO3W,sPO4S==>sPO4W
!  Burial process involve: sDIMS==>,sD/N/PDetS==>,sPAIMS==>
!  feh: Sep.8
!  Diatom Si sedimentation and resuspension can't ben handled here, since SiDiat is not state variable both 
!  in water column and sediment. Something could be further considered.
! !USES:
    Use Hydrodynamic
    Use MeshVars
    Use LimnologyVars
    Use uTempFunction

    Implicit None
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam
    type(LimnologyParam) :: LimnoParam
    Integer:: iElem,iLayer,sed,ff,mm,gg,dom,bb
    Real:: dt,dtday
    Real:: NearZero = 1e-10
    Real:: sDFiAdBottom,sDVeg
    Real:: DIR_VENTO,fetchm
    Real:: V

    Do iElem = 1,MeshParam%nElem
        !-----------------------------------------------------------------------
        !  Current nutrients ratios(check the curent state)
        !-----------------------------------------------------------------------
        !rPDDetS=sPDetS(iElem)/(sDDetS(iElem)+NearZero)
        !rNDDetS=sNDetS(iElem)/(sDDetS(iElem)+NearZero)
        !rSiDDetS=sSiDetS(iElem)/(sDDetS(iElem)+NearZero)
        !-----------------------------------------------------------------------
        !  Temperature functions for sediment abiotic process
        !-----------------------------------------------------------------------
        !  temperature_correction_of_sedimentation
        LimnoParam%uFunTmSet= uFunTmAbio(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem),iElem),LimnoParam%cThetaSet)
        !------------------------------------------------------------------------------------------------------------
        !  resuspension and sedimentation 
        !------------------------------------------------------------------------------------------------------------
        !------------------------------------------------------------------------------------------------------------
        !   The resuspended matter in the water column calculation
        !------------------------------------------------------------------------------------------------------------
        !  bioturbation_by_fish
        LimnoParam%tDTurbFish = 0.
        Do ff = 1, LimnoParam%nfish
            If (LimnoParam%afish(ff) == 2 .or. LimnoParam%afish(ff) == 3) Then !Omnivorous adult (afish = 2) or Planktivous adult (afish = 3) fish in the bottom layer
                sDFiAdBottom = LimnoParam%sDFiAd(HydroParam%ElSmallm(iElem),iElem,ff)
                LimnoParam%uFunTmFish = uFunTmBio(LimnoParam%sDTempW(HydroParam%ElSmallm(iElem),iElem),LimnoParam%cSigTmFish(ff),LimnoParam%cTmOptFish(ff))
                LimnoParam%tDTurbFish = LimnoParam%tDTurbFish + LimnoParam%kTurbFish*LimnoParam%uFunTmFish*sDFiAdBottom
            EndIf
        EndDo
        
        !Calculate lake fetch for each cell
        DIR_VENTO = ATAN2(HydroParam%Windiy(iElem),HydroParam%Windix(iElem))-HydroParam%PI-3.*HydroParam%PI/4.
        IF (DIR_VENTO < -2.*HydroParam%PI) THEN
            DIR_VENTO = DIR_VENTO + 2.*HydroParam%PI
        ELSEIF (DIR_VENTO == -2.*HydroParam%PI) THEN
            DIR_VENTO = 2.*HydroParam%PI
        ENDIF 
        DIR_VENTO = NINT(ABS(DIR_VENTO/(HydroParam%PI/4.))) !índice da matrix fetch (de acordo com a direção do vento)
        DIR_VENTO = MIN(DIR_VENTO,8.)
        DIR_VENTO = MAX(DIR_VENTO,1.)
        fetchm = HydroParam%fetch_m(iElem,int(DIR_VENTO)) !fetch (m)
            
        !  calculate resuspension rate, two methods
        If (LimnoParam%ResuspMethodFlag == 0) Then ! No resuspension
        
        ElseIf (LimnoParam%ResuspMethodFlag == 1) Then !PCLake
            if (LimnoParam%sDTempW(HydroParam%ElSmallm(iElem),iElem) >= 0.1) then
                LimnoParam%aFunDimSusp=LimnoParam%cSuspRef*((LimnoParam%cSuspMin+LimnoParam%cSuspMax/(1.0+exp(LimnoParam%cSuspSlope*&
                & (V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem)) -LimnoParam%hDepthSusp))))*((((fetchm)/ (LimnoParam%cFetchRef+NearZero)) )** (0.5)))
            else
                LimnoParam%aFunDimSusp=0.0
            endif
            LimnoParam%tDResusTauDead=min(LimnoParam%aFunDimSusp, ((LimnoParam%aFunDimSusp +NearZero )**(0.5))) &
            &*((LimnoParam%fLutum/ LimnoParam%fLutumRef )** (0.5))*LimnoParam%bPorS

            !  vegetation_dependence_of_resuspension (see Later)
            sDVeg = 0
            Do mm = 1, LimnoParam%nmac
                If (LimnoParam%amac(mm) == 1 .or. LimnoParam%amac(mm) == 5) Then !Rooted (amac = 1) or general (amac = 5) macrophytes
                    sDVeg = sDVeg + LimnoParam%sDMac(iElem,mm)
                EndIf
            EndDo
            LimnoParam%aFunVegResus=max(0.0,1.0-LimnoParam%kVegResus*sDVeg)

            !  resuspension_due_to_shear_stress_AND_fish,_corrected_for_vegetation_effect, (1/day)
            LimnoParam%tDResusDead=(LimnoParam%tDResusTauDead+LimnoParam%tDTurbFish)*LimnoParam%aFunVegResus !/secs_pr_day
        ElseIf (LimnoParam%ResuspMethodFlag == 2) Then !Garcia & Parker
            
        EndIf

        !------------------------------------------------------------------------------------------------------------
        !  Resuspension Fluxes
        !------------------------------------------------------------------------------------------------------------
        If (LimnoParam%ResuspMethodFlag == 0) Then ! No resuspension
            !  The inorganic matter resuspension
            LimnoParam%tDResusIM(iElem)=0.
            !  detrital_resuspension_DW
            LimnoParam%tDResusPom(iElem,:)=0.
            !  detrital_resuspension_C
            LimnoParam%tCResusPom(iElem,:)=0.
            !  detrital_resuspension_P
            LimnoParam%tPResusPom(iElem,:)=0.
            !  detrital_resuspension_N
            LimnoParam%tNResusPom(iElem,:)=0.
            
            LimnoParam%tDResusDom(iElem,:) = 0.
            LimnoParam%tCResusDom(iElem,:) = 0.
            LimnoParam%tNResusDom(iElem,:) = 0.
            LimnoParam%tPResusDom(iElem,:) = 0.
            !  detrital_resuspension_SI
    !        tSiResusPom(iElem,:)=Dot_Product(rSiDPomS(iElem,:),tDResusPom(iElem,:))
            !  resuspension_nutrient_P
            LimnoParam%tPResusPO4(iElem)=0.
            !  resuspension_absorbed_PAIM
            LimnoParam%tPResusAIM(iElem)=0.
            !  resuspension_nutrient_NO3
            LimnoParam%tNResusNO3(iElem)=0.
            !  resuspension_nutrient_NH4
            LimnoParam%tNResusNH4(iElem)=0.
        
            !!  Algae group resuspension
            !  resuspension of phytoplankton,DW
            LimnoParam%tDResusPhyt(iElem,:) = 0.
            !  resuspension of phytoplankton,C
            LimnoParam%tCResusPhyt(iElem,:) = 0.
            !  resuspension of phytoplankton,N
            LimnoParam%tNResusPhyt(iElem,:) = 0.
            !  resuspension of phytoplankton,P
            LimnoParam%tPResusPhyt(iElem,:) = 0.
        
        ElseIf (LimnoParam%ResuspMethodFlag == 1) Then !PCLake
            !  The inorganic matter resuspension
            LimnoParam%tDResusIM(iElem)=LimnoParam%fLutum*LimnoParam%sDIMS(iElem)/(LimnoParam%fLutum*LimnoParam%sDIMS(iElem)+SUM(LimnoParam%sDPomS(iElem,:))+NearZero)*LimnoParam%tDResusDead
            !  detrital_resuspension_DW
            LimnoParam%tDResusPom(iElem,:)=LimnoParam%sDPomS(iElem,:)/(LimnoParam%fLutum*LimnoParam%sDIMS(iElem)+LimnoParam%sDPomS(iElem,:)+NearZero)*LimnoParam%tDResusDead
            !  detrital_resuspension_C
            LimnoParam%tCResusPom(iElem,:)=Dot_Product(LimnoParam%rCDPomS(iElem,:),LimnoParam%tDResusPom(iElem,:))
            !  detrital_resuspension_P
            LimnoParam%tPResusPom(iElem,:)=Dot_Product(LimnoParam%rPDPomS(iElem,:),LimnoParam%tDResusPom(iElem,:))
            !  detrital_resuspension_N
            LimnoParam%tNResusPom(iElem,:)=Dot_Product(LimnoParam%rNDPomS(iElem,:),LimnoParam%tDResusPom(iElem,:))
            
            LimnoParam%tDResusDom(iElem,:) = 0. !LimnoParam%sDDomS(iElem,:)/SUM(LimnoParam%sDPomS(iElem,:))*SUM(LimnoParam%tDResusPom(iElem,:))
            LimnoParam%tCResusDom(iElem,:) = 0. !Dot_Product(LimnoParam%rCDDomS(iElem,:) , LimnoParam%tDResusDom(iElem,:))
            LimnoParam%tNResusDom(iElem,:) = 0. !Dot_Product(LimnoParam%rNDDomS(iElem,:) , LimnoParam%tDResusDom(iElem,:))
            LimnoParam%tPResusDom(iElem,:) = 0. !Dot_Product(LimnoParam%rPDDomS(iElem,:) , LimnoParam%tDResusDom(iElem,:))

            
            
            
            !  detrital_resuspension_SI
    !        tSiResusPom(iElem,:)=Dot_Product(rSiDPomS(iElem,:),tDResusPom(iElem,:))
            !  resuspension_nutrient_P
            LimnoParam%tPResusPO4(iElem)=LimnoParam%sPO4S(iElem)*Dot_Product(1./(LimnoParam%sDPomS(iElem,:)+NearZero),LimnoParam%tDResusPom(iElem,:))
            !  resuspension_absorbed_PAIM
            LimnoParam%tPResusAIM(iElem)=LimnoParam%sPAIMS(iElem)/(LimnoParam%sDIMS(iElem)+NearZero)*LimnoParam%tDResusIM(iElem)
            !  resuspension_nutrient_NO3
            LimnoParam%tNResusNO3(iElem)=LimnoParam%sNO3S(iElem)*Dot_Product(1./(LimnoParam%sDPomS(iElem,:)+NearZero),LimnoParam%tDResusPom(iElem,:))
            !  resuspension_nutrient_NH4
            LimnoParam%tNResusNH4(iElem)=LimnoParam%sNH4S(iElem)*Dot_Product(1./(LimnoParam%sDPomS(iElem,:)+NearZero),LimnoParam%tDResusPom(iElem,:))   
        
            !  phytoplankton_resuspension_rate_constant, in day
            LimnoParam%akResusPhytRef = LimnoParam%kResusPhytMax * (1.0 - exp(LimnoParam%cResusPhytExp * LimnoParam%tDResusDead))
            !!  Algae group resuspension
            !  resuspension of phytoplankton,DW
            LimnoParam%tDResusPhyt(iElem,:) = LimnoParam%akResusPhytRef*LimnoParam%sDPhytS(iElem,:)
            !  resuspension of phytoplankton,C
            LimnoParam%tCResusPhyt(iElem,:) = LimnoParam%rCDPhytS(iElem,:) * LimnoParam%tDResusPhyt(iElem,:)
            !  resuspension of phytoplankton,N
            LimnoParam%tNResusPhyt(iElem,:) = LimnoParam%rNDPhytS(iElem,:) * LimnoParam%tDResusPhyt(iElem,:)
            !  resuspension of phytoplankton,P
            LimnoParam%tPResusPhyt(iElem,:) = LimnoParam%rPDPhytS(iElem,:) * LimnoParam%tDResusPhyt(iElem,:)
            
        EndIf
        
        
        !-----------------------------------------------------------------------
        !  The sedimentation calculation
        !-----------------------------------------------------------------------
        If (LimnoParam%SetMethodFlag == 0) Then ! No sedimentation

        ElseIf (LimnoParam%SetMethodFlag==1) Then !PCLake (sedimentation calculation, based on resuspension)
            !  correction_factor_for_settling_rate_(<=_1),basic settling rate, in day
            LimnoParam%aFunTauSet=min(1.0,1.0/((LimnoParam%aFunDimSusp + NearZero)**(0.5)))
            !-----------------------------------------------------------------------
            !  Different matter sedimentation based on the basic settling rate
            !-----------------------------------------------------------------------
            !  sedimentation_velocity_of_IM, in day
            LimnoParam%uCorVSetIM=LimnoParam%aFunTauSet*((LimnoParam%fLutumRef/(LimnoParam%fLutum + NearZero))**(0.5))*LimnoParam%uFunTmSet*LimnoParam%cVSetIM
            !  sedimentation_velocity_of_detritus, in day
            LimnoParam%uCorVSetDet=LimnoParam%cVSetPom*LimnoParam%aFunTauSet*LimnoParam%uFunTmSet
            !  corrected_sedimentation_velocity_of_Algae, in day
            Do gg = 1, LimnoParam%nphy
                LimnoParam%uCorVSetPhyt(gg) = LimnoParam%cVSetPhyt(gg) * LimnoParam%aFunTauSet * LimnoParam%uFunTmSet
            EndDo
            
            
        ElseIf (LimnoParam%SetMethodFlag==2) Then !Stokes
			!-------------------------------------------------------------------------!
			! (2) Sedimentacao POM (Stokes Equation)			                      !	
			!-------------------------------------------------------------------------!
			! (2.1) Get the viscosity of the water and particle density                 

		    !------------------------------------------------------------------!
		    ! Calculates the molecular viscosity of water for a given          !  
		    ! temperature                                                      !  
		    !                                                                  !  
		    ! NOTE: to use an array of value pass temp as -1.0, actually	   !  
		    !       any realnumber will do since temp is ignored if            !  
		    !       temp_array ispresent                                       !  
		    ! I had to fit lines to tabulated values for viscosity:            !  
		    !                                                                  !  
		    !  Temp (C)     Viscosity (N s / m^2) x 10^3                       !  
		    !  --------     ---------                                          !  
		    !      0          1.781                                            !  
		    !      5          1.518                                            !  
		    !     10          1.307                                            ! 
		    !     15          1.139                                            ! 
		    !     20          1.002                                            !  
		    !     25          0.890                                            !  
		    !     30          0.798                                            !  
		    !     40          0.653                                            ! 
		    !     50          0.547                                            !  
		    !     60          0.466                                            !  
		    !     70          0.404                                            !  
		    !     80          0.354                                            !  
		    !     90          0.315                                            !  
		    !    100          0.282                                            !  
		    !                                                                  !  
		    ! From Table A.1b, _FLUID_MECHANICS_With_Engineering_Applications_ !  
		    ! by Robert L. Daugherty and Joseph B. Franzini, however, these    !  
		    ! values are common in most fluid mechanics texts.                 !  
		    ! NOTE: N s / m^2  = kg / m / s                                    !  
		    !------------------------------------------------------------------!
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
			    If (LimnoParam%sDTempW(iLayer,iElem) <= 20.) Then
				    ! 0C to 20C
				    ! y = 0.0008 * x^2 - 0.0556 * x + 1.7789
				    ! r^2 = 0.9999
				    LimnoParam%mu =  0.0008 * LimnoParam%sDTempW(iLayer,iElem) * LimnoParam%sDTempW(iLayer,iElem) - 0.0556 * LimnoParam%sDTempW(iLayer,iElem) + 1.7789
    		
			    ElseIf(LimnoParam%sDTempW(iLayer,iElem) <= 60.) Then 
				    ! 20C to 60C
				    ! y = 0.0002 * x^2 - 0.0323 * x + 1.5471
				    ! r^2 = 0.9997
				    LimnoParam%mu =  0.0002 * LimnoParam%sDTempW(iLayer,iElem) * LimnoParam%sDTempW(iLayer,iElem) - 0.0323 * LimnoParam%sDTempW(iLayer,iElem) + 1.5471
			    Else
				    ! 60C to 100C
				    ! y = 0.00006 * x^2 - 0.0141 * x + 1.1026
				    ! r^2 = 0.9995
				    LimnoParam%mu = 0.00006 * LimnoParam%sDTempW(iLayer,iElem) * LimnoParam%sDTempW(iLayer,iElem) - 0.0141 * LimnoParam%sDTempW(iLayer,iElem) + 1.1026
			    EndIf

			    ! Now divide it to the proper units (N s / m^2)
                LimnoParam%mu = LimnoParam%mu / 1000.00
    	
			    ! (2.2) Get the vertical velocities via Stokes settling formula (m/s)
                LimnoParam%StokesVelPom(iLayer,:) = MAX(0.0 , 1.0 * (9.8 * (LimnoParam%ppOM - HydroParam%sDRhoW(iLayer,iElem)) * LimnoParam%diOM * LimnoParam%diOM / (18.0 * LimnoParam%mu)))
                LimnoParam%StokesVelIM(iLayer) =       MAX(0.0 , 1.0 * (9.8 * (LimnoParam%ppIM - HydroParam%sDRhoW(iLayer,iElem)) * LimnoParam%diIM * LimnoParam%diIM / (18.0 * LimnoParam%mu)))
                Do gg = 1, LimnoParam%nphy
                    LimnoParam%StokesVelPhyt(iLayer,gg) = MAX(0.0 , 1.0 * (9.8 * (LimnoParam%ppPhyt(gg) - HydroParam%sDRhoW(iLayer,iElem)) * LimnoParam%diPhyt(gg) * LimnoParam%diPhyt(gg) / (18.0 * LimnoParam%mu)))
                EndDo
                
            EndDo
            
        EndIf
        !
        !------------------------------------------------------------------------------------------------------------
        !  Sedimentation Fluxes
        !------------------------------------------------------------------------------------------------------------
        If (LimnoParam%SetMethodFlag == 0) Then ! No sedimentation
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                !  sedimentation_IM
                LimnoParam%tDSetIM(iLayer,iElem)=0.
                !  sedimentation_PAIM
                LimnoParam%tPSetAIM(iLayer,iElem)=0.
                !  sedimentation_flux_of_detritus
                LimnoParam%tDSetPom(iLayer,iElem,:)=0.
                !  sedimentation_detrital_C
                LimnoParam%tCSetPom(iLayer,iElem,:)=0.
                !  sedimentation_detrital_P
                LimnoParam%tPSetPom(iLayer,iElem,:)=0.
                !  sedimentation_detrital_N
                LimnoParam%tNSetPom(iLayer,iElem,:)=0.
                !  sedimentation_detrital_Si
                !tSiSetPom(iLayer,iElem,:)=uCorVSetDet*sSiPomW(ElSmallm(iElem),iElem,:)
                !  sedimentation_flux_of_Algae
                LimnoParam%tDSetPhyt(iLayer,iElem,:)=0.
                !  sedimentation_flux_of_Bacteria
                LimnoParam%tDSetBac(iLayer,iElem,:)=0.
            EndDo

        ElseIf (LimnoParam%SetMethodFlag==1) Then !PCLake (sedimentation calculation, based on resuspension)
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                !  sedimentation_IM
                LimnoParam%tDSetIM(iLayer,iElem)=LimnoParam%uCorVSetIM*LimnoParam%sDIMW(iLayer,iElem)
                !  sedimentation_PAIM
                LimnoParam%tPSetAIM(iLayer,iElem)=LimnoParam%sPAIMW(iLayer,iElem)/(LimnoParam%sDIMW(iLayer,iElem) +NearZero)*LimnoParam%tDSetIM(iLayer,iElem)
                !  sedimentation_flux_of_detritus
                LimnoParam%tDSetPom(iLayer,iElem,:)=LimnoParam%uCorVSetDet*LimnoParam%sDPomW(iLayer,iElem,:)
                !  sedimentation_detrital_C
                LimnoParam%tCSetPom(iLayer,iElem,:)=LimnoParam%uCorVSetDet*LimnoParam%sCPomW(iLayer,iElem,:)
                !  sedimentation_detrital_P
                LimnoParam%tPSetPom(iLayer,iElem,:)=LimnoParam%uCorVSetDet*LimnoParam%sPPomW(iLayer,iElem,:)
                !  sedimentation_detrital_N
                LimnoParam%tNSetPom(iLayer,iElem,:)=LimnoParam%uCorVSetDet*LimnoParam%sNPomW(iLayer,iElem,:)
                !  sedimentation_detrital_Si
                !tSiSetPom(iLayer,iElem,:)=uCorVSetDet*sSiPomW(ElSmallm(iElem),iElem,:)
                !  sedimentation_flux_of_Algae
                LimnoParam%tDSetPhyt(iLayer,iElem,:)=LimnoParam%uCorVSetPhyt(:)*LimnoParam%sDPhytW(iLayer,iElem,:)
                !  sedimentation_flux_of_Bacteria
                If (LimnoParam%InclBac == 1) Then
    	            Do dom = 1, LimnoParam%ndom
    		            LimnoParam%oDFoodBacDom(dom) = LimnoParam%cPrefBacDom(bb,dom) * LimnoParam%sDDomW(iLayer,iElem,dom)
                        LimnoParam%oDSetBacSum(dom) = LimnoParam%cPrefBacDom(bb,dom) * LimnoParam%tDSetPom(iLayer,iElem,dom)
    	            EndDo
    	            LimnoParam%oDFoodBac = Sum(LimnoParam%oDFoodBacDom(:))
                    LimnoParam%oDSetBac = Sum(LimnoParam%oDSetBacSum(:))/Sum(LimnoParam%cPrefBacDom(bb,:))
                    LimnoParam%tDSetBac(iLayer,iElem,:)= (LimnoParam%oDFoodBac / (LimnoParam%oDFoodBac + LimnoParam%hSedBac(:) +NearZero))*LimnoParam%oDSetBac
                Else
                    LimnoParam%tDSetBac(iLayer,iElem,:) = 0.
                EndIf

            EndDo
        ElseIf (LimnoParam%SetMethodFlag==2) Then !Stokes
			! Fluxo de Sedimentacao (mg/m2/d)
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
			    ! Estande de macrófitas e vvel < 0 = fator que aumenta sedimentacao da particula		
			    If (MeshParam%BANHADO(iElem)==1) Then
				    LimnoParam%tDSetPom(iLayer,iElem,:) = 2.  * LimnoParam%StokesVelPom(iLayer,:)*86400.* LimnoParam%sDPomW(iLayer,iElem,:)   
                    LimnoParam%tDSetIM(iLayer,iElem)      = 2.  * LimnoParam%StokesVelIM(iLayer) * 86400. * LimnoParam%sDIMW(iLayer,iElem)
                    LimnoParam%tDSetPhyt(iLayer,iElem,:) = 2. * LimnoParam%StokesVelPhyt(iLayer,:)*86400.* LimnoParam%sDPhytW(iLayer,iElem,:)
			    ! Ponto nao esta localizado em estande de macrófitas = nao há reducao no fluxo		
			    Else
				    LimnoParam%tDSetPom(iLayer,iElem,:) =       LimnoParam%StokesVelPom(iLayer,:)*86400.* LimnoParam%sDPomW(iLayer,iElem,:)
                    LimnoParam%tDSetIM(iLayer,iElem)      =       LimnoParam%StokesVelIM(iLayer) * 86400. * LimnoParam%sDIMW(iLayer,iElem)
                    LimnoParam%tDSetPhyt(iLayer,iElem,:) =      LimnoParam%StokesVelPhyt(iLayer,:)*86400.* LimnoParam%sDPhytW(iLayer,iElem,:)
                EndIf
			    !tDSetPom(KAUX,pom) = 0.0
			    ! (2.4) Calcula os compartimentos de nutrientes C, N e P  na POM no fluxo de sedimentacao (mg/m3/s)
                !  sedimentation_PAIM
                LimnoParam%tPSetAIM(iLayer,iElem)=LimnoParam%sPAIMW(iLayer,iElem)/(LimnoParam%sDIMW(iLayer,iElem) +NearZero)*LimnoParam%tDSetIM(iLayer,iElem)
                LimnoParam%tCSetPom(iLayer,iElem,:) = LimnoParam%rCDPomW(iLayer,iElem,:) * LimnoParam%tDSetPom(iLayer,iElem,:)
                LimnoParam%tNSetPom(iLayer,iElem,:) = LimnoParam%rNDPomW(iLayer,iElem,:) * LimnoParam%tDSetPom(iLayer,iElem,:)
                LimnoParam%tPSetPom(iLayer,iElem,:) = LimnoParam%rPDPomW(iLayer,iElem,:) * LimnoParam%tDSetPom(iLayer,iElem,:)
                !  sedimentation_flux_of_Bacteria
                If (LimnoParam%InclBac == 1) Then
    	            Do dom = 1, LimnoParam%ndom
    		            LimnoParam%oDFoodBacDom(dom) = LimnoParam%cPrefBacDom(bb,dom) * LimnoParam%sDDomW(iLayer,iElem,dom)
                        LimnoParam%oDSetBacSum(dom) = LimnoParam%cPrefBacDom(bb,dom) * LimnoParam%tDSetPom(iLayer,iElem,dom)
    	            EndDo
    	            LimnoParam%oDFoodBac = Sum(LimnoParam%oDFoodBacDom(:))
                    LimnoParam%oDSetBac = Sum(LimnoParam%oDSetBacSum(:))/Sum(LimnoParam%cPrefBacDom(bb,:))
                    LimnoParam%tDSetBac(iLayer,iElem,:)= (LimnoParam%oDFoodBac / (LimnoParam%oDFoodBac + LimnoParam%hSedBac(:) +NearZero))*LimnoParam%oDSetBac
                Else
                    LimnoParam%tDSetBac(iLayer,iElem,:)=0.
                EndIf
                
                 
            EndDo
        EndIf
    EndDo
Return
End

