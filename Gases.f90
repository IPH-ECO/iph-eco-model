SUBROUTINE Gases(HydroParam,MeshParam,LimnoParam,dt,dtday) !OD,DIC & CH4

    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars
    
    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
	Integer:: nph,iter
	Real:: pCO2Sat	
	Real:: aAlkW,aDicW
	Real:: aTkelvin
	Real:: Scnumber	
	Real:: aPH,phaux,Hion
	Real:: k1,k2,kwater
	Real:: incr,Fw
	Real:: a0,a1,a2
	Real:: uDLoadDicaux
	Real:: uO2LoadO2aux
	Real:: dAlkW
    Real:: ScCO2
    Real:: kHenryCO2
    Real:: V
    Integer:: index,iElem,iEdge,iLayer
    Real:: tPDifO2_bottom,tCDifDic_bottom
    Real:: dt,dtday
    Real:: NearZero = 1e-10
    
!   Inicializa Variaveis da Rotina    
    nph=0.
    aAlkW=0.
    aDicW=0.
    aTkelvin=0.
    aPH=0.
    phaux=0.
    Hion=0.
    k1=0.
    k2=0.
    kwater=0.
    incr=0.
    Fw=0.
    a0=0.
    a1=0.
    a2=0.
	uDLoadDicaux=0.
	uO2LoadO2aux=0.
	
!-OD----------------------------------------------------------------------------------!
! 1. Exchange to/from the air-water interface                                         !
! 2. Photosynthetic oxygen production and respiratory oxygen                          !
!    consumption by phytoplankton.                                                    !
! 3. Photosynthetic oxygen production and respiratory oxygen                          !
!    consumption by macrophytes.                                                       !
! 4. Utilization of oxygen due to the action of bacteria on                           !
!    organic matter.                                                                  !
! 5. Utilization of oxygen in the process of nitrification.                           !
! 6. Exchange of oxygen at the sediment/water interface.                              !
!-------------------------------------------------------------------------------------!
  If (LimnoParam%iOxyBOD==1) Then  
      
    !-------------------------------------------------------------------------------------!
    !                               Boundary Condition for O2				      !
    !-------------------------------------------------------------------------------------!
    index = 11
    HydroParam%uLoadVarEst = LimnoParam%uO2LoadO2
            
    Do iElem = 1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
!-------------------------------------------------------------------------------------!
!                               O2 Processes in Water					      !
!-------------------------------------------------------------------------------------!
            !1. Exchange to/from the air-water interface
            If (iLayer==HydroParam%ElCapitalM(iElem)) Then !Surface Layer
                
			    !Oxygen saturation concentration (mgO2/l)
			    LimnoParam%uO2Sat=14.652-0.41022*LimnoParam%sDTempW(iLayer,iElem)+7.991e-3*LimnoParam%sDTempW(iLayer,iElem)**2.-7.7774e-5*LimnoParam%sDTempW(iLayer,iElem)**3. 
                
			    LimnoParam%kReaer=0.727*SQRT(MeshParam%CREDV(iElem)*sqrt(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.)     &
			                -0.371*MeshParam%CREDV(iElem)*sqrt(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.)     &
			                +0.0376*(MeshParam%CREDV(iElem)*sqrt(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.)))**2.
    			
			    LimnoParam%tO2Reaer(iLayer)=LimnoParam%kReaer*(LimnoParam%uO2Sat-LimnoParam%sO2W(iLayer,iElem))*(LimnoParam%cThetaReaer**(LimnoParam%sDTempW(iLayer,iElem)-20.))
            Else
                LimnoParam%tO2Reaer(iLayer)=0.0
            EndIf
            !2. Difusion
            If (iLayer == HydroParam%ElSmallm(iElem)) Then
                tPDifO2_bottom = LimnoParam%tO2Dif(iElem)
            Else
                tPDifO2_bottom = 0.
            EndIf
            
            ! 6. Oxygen consumption at the sediment/water interface.
            If (iLayer==HydroParam%ElSmallm(iElem)) Then !Bottom Layer
	            !akO2DifCor=kO2Dif*(cThetaDif**(TEMPI(I,J,NCAM(I,J))-20.))*cturbDifO2*bPorCorS	  !Coeficiente de difusão de O2 corrigido
	            !tSOD=(molO2molC*cCPerDW*(1.-fRefrDetS)*tDMinDetS+O2PerNH4*molO2molN*kNitrS*(cThetaNitr**(TEMPI(I,J,NCAM(I,J))-20.))*sNH4S(I,J))/cDepthS !Demanda de oxigenio no sedimento
	            !aDepthOxySed=(2.0*sO2W(I,J,NCAM(I,J))*akO2DifCor/tSOD)**0.5 !Profundidade de penetração de oxigenio
	            !afOxySed=aDepthOxySed/cDepthS !Propoção aerobica no sedimento (-)
	            !tO2MinDetS=molO2molC*cCPerDW*afOxySed*(1.0-fRefrDetS)*tDMinDetS !Consumo de O2 devido a mineração de detritos no sedimento
            Else
                !tO2MinDetS = 0.
            EndIf
            
            
!-------------------------------------------------------------------------------------!
!                               Source/Sink Term for O2				      !
!-------------------------------------------------------------------------------------!
            
            HydroParam%dVarEst(iLayer,iElem,1) = LimnoParam%tO2Reaer(iLayer)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))              & ! Exchange to/from the air-water interface
                    -LimnoParam%molO2molC*LimnoParam%cCPerDW*SUM(LimnoParam%wDMinAerW(iLayer,iElem,:))                                      & ! Mineralization (InclBac=0)
                    +LimnoParam%molO2molC*LimnoParam%cCPerDW*SUM(LimnoParam%wDAssPhyt(iLayer,iElem,:))                                      & ! Phytoplankton Production
                    -LimnoParam%molO2molC*LimnoParam%cCPerDW*LimnoParam%aCorO2BOD*SUM(LimnoParam%wDRespPhyt(iLayer,iElem,:))                & ! Phytoplankton Respiration
                    +LimnoParam%molO2molN*LimnoParam%O2PerNO3*SUM(LimnoParam%wNUptNO3Phyt(iLayer,iElem,:))                                  & !  O2_production_due_to_NO3_uptake_by_phytopl.
                    +LimnoParam%molO2molC*LimnoParam%cCPerDW*SUM(LimnoParam%tDProdMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))           & ! PP Mac
                    -LimnoParam%molO2molC*LimnoParam%cCPerDW*SUM(LimnoParam%tDRespMacW(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))          & ! Resp Macro
                    -LimnoParam%molO2molC*SUM(LimnoParam%wCRespZoo(iLayer,iElem,:))                                                         & ! Zooplankton Respiration
                    -LimnoParam%molO2molC*SUM(LimnoParam%wCRespFiAd(iLayer,iElem,:))                                                        & ! Adult Fish Respiration
                    -LimnoParam%molO2molC*SUM(LimnoParam%wCRespFiJv(iLayer,iElem,:))                                                        & ! Juveline Fish Respiration
                    -LimnoParam%molO2molN*LimnoParam%O2PerNH4*LimnoParam%wNNitrW(iLayer,iElem)                                              & ! Nitrification
                    +tPDifO2_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))                                                       ! Difusion
                                                                                                                ! Oxygen consumption at the sediment/water interface.
      
!-------------------------------------------------------------------------------------!
!                               Source/Sink Term for BOD				      !
!-------------------------------------------------------------------------------------!
            
            HydroParam%dVarEst(iLayer,iElem,2) =   -LimnoParam%molO2molC*LimnoParam%cCPerDW*SUM(LimnoParam%wDMinAerW(iLayer,iElem,:))                       & ! Mineralization in the water
                                        -LimnoParam%molO2molN*LimnoParam%O2PerNH4*LimnoParam%wNNitrW(iLayer,iElem)                                          ! Nitrification in the water
                                                                                                                ! Mineralization and Nitrification in the sediment

        EndDo
        
    EndDo
            

!-------------------------------------------------------------------------------------!
!                       Solver Transport Equation for O2   	    	      !
!-------------------------------------------------------------------------------------!
    If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
        Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sO2W,LimnoParam%sO2WP,dt,dtday,HydroParam,MeshParam) 
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
        Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sO2W,LimnoParam%sO2WP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
    ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
        Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sO2W,LimnoParam%sO2WP,dt,dtday,HydroParam,MeshParam)
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
        Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sO2W,LimnoParam%sO2WP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
    EndiF        

    
    !-------------------------------------------------------------------------------------!
    !                               Boundary Condition for BOD				      !
    !-------------------------------------------------------------------------------------!
    index = 12
    HydroParam%uLoadVarEst = LimnoParam%uDLoadDBO

    
!-------------------------------------------------------------------------------------!
!                       Solver Transport Equation for BOD   	    	      !
!-------------------------------------------------------------------------------------!
    If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
        Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sDBOW,LimnoParam%sDBOWP,dt,dtday,HydroParam,MeshParam) 
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
        Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sDBOW,LimnoParam%sDBOWP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
    ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
        Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sDBOW,LimnoParam%sDBOWP,dt,dtday,HydroParam,MeshParam)
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
        Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,2),LimnoParam%sDBOW,LimnoParam%sDBOWP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
    EndIf
  EndIf
  
    
!-DIC---------------------------------------------------------------------------------!
! 1. Atmospheric Exchange                                                             !
! 2. Carbonate Kinetics                                                               !
! 3. CO2  Phyto*			    			                                          !
! 4. CO2  Bacteria*			    		                                              !
! 5. CO2  Zooplankton*		    	        	                                      !
! 6. CO2  Macrophytes *                                                                !
! 7. CO2  Fishes *                                                                    !
! 8. Sediment fluxes*    		        		                                  !
!-------------------------------------------------------------------------------------!
  If (LimnoParam%iCarbon==1) Then
      
    !-------------------------------------------------------------------------------------!
    !                               Boundary Condition for DIC				      !
    !-------------------------------------------------------------------------------------!
    index = 13
    HydroParam%uLoadVarEst = LimnoParam%uDLoadDic
      
    Do iElem = 1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            !-------------------------------------------------------------------------!
            ! 1. Surface CO2 flux                                            !
            !-------------------------------------------------------------------------!
            ! METODO 1) Q U A 2 K
            If (iLayer==HydroParam%ElCapitalM(iElem)) Then !Surface Layer
                ! Coeficiente de Troca Gasosa
			    LimnoParam%kReaer=0.727*SQRT(MeshParam%CREDV(iElem)*sqrt(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.)     &
			                -0.371*MeshParam%CREDV(iElem)*sqrt(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.)     &
			                +0.0376*(MeshParam%CREDV(iElem)*sqrt(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.)))**2.
                LimnoParam%kCO2Exch = (32./44.)**0.25 * LimnoParam%kReaer * (LimnoParam%cThetaReaer**(LimnoParam%sDTempW(iLayer,iElem)-20.)) !m/dia
            
                ! Constante de Henry
                !pkHenry = -2385.73/TEMPI(I,J,KAUX) - 0.0152642*TEMPI(I,J,KAUX) + 14.0184 ![mole (L.atm)-1] 	Convert to g/m3/atm:  pkHenry = pkHenry * 12.0 * 1000.0
                !kHenry = 10**-(pkHenry)
                LimnoParam%kHenry = 10**-1.46 ![moleCO2/L.atm]
                ! Saturacao de CO2 na agua
                LimnoParam%pCO2 = 10**-3.416  ![atm]
                pCO2Sat = LimnoParam%kHenry * LimnoParam%pCO2 * 12000. ! molCO2/L.atm * atm = molCO2/L = molCO2/L*[12gC/molCO2]*1000 -> gC/m3 , pCO2 = Values in 2007 are approximately 10^-3.416 atm (= 383.7 ppm)
    			
                ! Fluxo de CO2 na interface
                LimnoParam%tCO2Reaer(iLayer)=LimnoParam%kCO2Exch*(pCO2Sat - LimnoParam%sH2CO3W(iLayer,iElem)) ! CO2 == H2CO3*, H2CO3* = CO2 + H2CO3, CO2>>H2CO3
			    
                !-------------------------------------------------------------------------!
                ! METODO 2) C A E D Y M
                ! Schmidt, Sc:
                ScCO2 = 1911.1 - 118.11*LimnoParam%sDTempW(iLayer,iElem) + 3.4527*LimnoParam%sDTempW(iLayer,iElem)**2. - 0.043132*LimnoParam%sDTempW(iLayer,iElem)**3.
                
                ! Gas transfer velocity, kCO2 (cm/hr):
                !kCO2 = 0.31 u^2 (Sc/600)^-0.5
                LimnoParam%kCO2Exch = 0.31 * (HydroParam%Windix(iElem)**2.          &
    			                +  HydroParam%Windiy(iElem)**2.)         &
    			                * (ScCO2/660.)**(-0.5)
                ! convert to m/day:
                LimnoParam%kCO2Exch = LimnoParam%kCO2Exch * 24.0 / 100.0
    			
                ! Henry - molCO2/L.atm
                aTkelvin = LimnoParam%sDTempW(iLayer,iElem) + 273.15
                kHenryCO2= exp(-58.0931 + 90.5069*(100./aTkelvin) + 22.294*log(aTkelvin/100.))
			    
                ! Dioxido de Carbono Dissolvido: 
                !atm            =         gC/m3    /   molCO2/L.atm*(12gC/molCO2)*(1000L/1m3)
                LimnoParam%pCO2w(iLayer,iElem) = LimnoParam%sH2CO3W(iLayer,iElem)/(LimnoParam%kHenry*12.*1000.)
    			
                ! Fluxo de CO2 na interface
                ! FCO2 = kCO2 * Ko * (pCO2 - PCO2a):
                ! g/m2/day      =  m/day * molCO2/L.atm * atm * (12gC/molCO2) * (1000L/1m3)
                LimnoParam%tCO2Reaer(iLayer) = LimnoParam%kCO2Exch*kHenryCO2*(LimnoParam%pCO2 - LimnoParam%pCO2w(iLayer,iElem))* 12.* 1000.
            Else
                LimnoParam%tCO2Reaer(iLayer) = 0.
            EndIf
            
            !2. Difusion
            If (iLayer == HydroParam%ElSmallm(iElem)) Then
                tCDifDic_bottom = LimnoParam%tCDifDic(iElem)
            Else
                tCDifDic_bottom = 0.
            EndIf
            
!-------------------------------------------------------------------------------------!
!                               Source/Sink Term for DIC				      !
!-------------------------------------------------------------------------------------!
            
            HydroParam%dVarEst(iLayer,iElem,1) =   +LimnoParam%tCO2Reaer(iLayer)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))  & ! Surface CO2 flux
                                        +SUM(LimnoParam%wCMinAer2CO2W(iLayer,iElem,:))                                              & ! Mineralization in the water (InclBac=0)
                                        +SUM(LimnoParam%wCRespBac2CO2(iLayer,iElem,:))                                              & ! Bacterioplankton Respiration (InclBac=1)
                                        +LimnoParam%cCPerDW*SUM(LimnoParam%wDRespPhyt(iLayer,iElem,:))                              & ! Phytoplankton Respiration
                                        -LimnoParam%cCPerDW*SUM(LimnoParam%wDAssPhyt(iLayer,iElem,:))                               & ! Phytoplankton Production
                                        +SUM(LimnoParam%tCExdMac(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))                      & !Resp Macro
                                        -SUM(LimnoParam%tCUptMacW(iElem,:))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))   & !PP Macro
                                        +SUM(LimnoParam%wCRespZoo(iLayer,iElem,:))                                                  & ! Zooplankton Respiration
                                        +SUM(LimnoParam%wCRespFiAd(iLayer,iElem,:))                                                 & !Adult Fish Respiration
                                        +SUM(LimnoParam%wCRespFiJv(iLayer,iElem,:))                                                 & !Juveline Fish Respiration
                                        +tCDifDic_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))                          ! Difusion
                                                                                                            !+ Metanotrofia
        EndDo
    EndDo
    
!-------------------------------------------------------------------------------------!
!                       Solver Transport Equation for DIC   	    	                !
!-------------------------------------------------------------------------------------!
    If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
        Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDicW,LimnoParam%sDicWP,dt,dtday,HydroParam,MeshParam) 
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
        Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDicW,LimnoParam%sDicWP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
    ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
        Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDicW,LimnoParam%sDicWP,dt,dtday,HydroParam,MeshParam)
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
        Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDicW,LimnoParam%sDicWP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
    EndiF        
    
!-------------------------------------------------------------------------!
! Carbonate equilibrium  		                                          !
!-------------------------------------------------------------------------!
    Do iElem = 1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            ! pH and carbonate species
            aAlkW = LimnoParam%sAlkW(iLayer,iElem) / 5.0E+04  !mg/L CaCO3 -> eq/L 
            ! mgCaCO3-- => eq 
            ! 1gCaCO3--  = 2/100 eq
            ! 1mgCaCO3-- = 2/(1E5) eq 
            ! 100 = massa molar CaCO3
            !   2 = carga
            aDicW = LimnoParam%sDicW(iLayer,iElem)/12000.  !mol/L   
            ! molar(mol/L)
            ! 1molC = 12g
            ! 1g = 1mol/12g/1000
            aPH = LimnoParam%spHW(iLayer,iElem)
	
            ! Temperature adjustment (sem considerar forca ionica) MUDAR!!!!!!!!!!!!!!!
            !k1 = 10**-(0.000142121212 * aTkelvin * aTkelvin - 0.012648181818 * aTkelvin + 6.577539393939)                 
            !k2 = 10**-(0.000113679654 * aTkelvin * aTkelvin - 0.014687186147 * aTkelvin + 10.625769696970)               
            !kwater = 10**-(0.000201991342 * aTkelvin * aTkelvin - 0.043419653680 * aTkelvin + 14.949709090909)            
            k1=10**-6.3
            k2=10**-10.3
            kwater=10**-14.
    
            phaux = -aPH-2.1
            If (aPH <= 0.0) phaux = -14.0
            incr = 10.0
            Do nph=1,3
            Fw    = 1.0
            incr = incr/10.0
            iter = 0
                Do While (Fw > 0.0 .and. iter < 12)
                    phaux  = phaux+INCR
                    Hion   = 10.0**phaux
                    a0 = Hion**2 / (Hion*Hion + k1*Hion + k1*k2)
                    a1 = k1*Hion / (Hion*Hion + k1*Hion + k1*k2)
                    a2 = k1*k2   / (Hion*Hion + k1*Hion + k1*k2)
                    Fw = (a1 + 2*a2) * aDicW + kwater/Hion - Hion - aAlkW
                    iter   = iter+1
                EndDo
                phaux = phaux-INCR
            EndDo

            Hion =  10.0 ** phaux
            LimnoParam%spHW(iLayer,iElem) = - phaux
	
            a0 = Hion*Hion/(Hion*Hion + k1*Hion + k1*k2)
            a1 = k1*Hion / (Hion*Hion + k1*Hion + k1*k2)
            a2 = k1*k2   / (Hion*Hion + k1*Hion + k1*k2)

            LimnoParam%sH2CO3W(iLayer,iElem)  = LimnoParam%sDicW(iLayer,iElem) * a0  !gC/m3
            LimnoParam%sHCO3W(iLayer,iElem) =   LimnoParam%sDicW(iLayer,iElem) * a1  !gC/m3
            LimnoParam%sCO3W(iLayer,iElem)  =   LimnoParam%sDicW(iLayer,iElem) * a2  !gC/m3    
            
!-------------------------------------------------------------------------!
!  Alkalinity correction  		                                          !
!-------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------!
!                               Source/Sink Term for DIC				      !
!-------------------------------------------------------------------------------------!
            
            HydroParam%dVarEst(iLayer,iElem,1) =   +SUM(LimnoParam%wNUptNO3Phyt(iLayer,iElem,:))*(50000.*0.001)/14.007                    & !Fotossintese - Aumenta:  106CO2 + 16NO3- + HPO4-2 + 122H2O + 18H+ -> (C106.H263.O110.N16.P1) + 138O2
                                        -SUM(LimnoParam%wNUptNH4Phyt(iLayer,iElem,:))*(50000.*0.001)/14.007                    & !Fotossintese
                                        +SUM(LimnoParam%wPUptPhyt(iLayer,iElem,:))*(50000.*0.001)/31.974                       & !Fotossintese
                                        +SUM(LimnoParam%tNUptNO3MacW(iElem,:))*(50000.*0.001)/14.007/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))                               & !Fotossintese - Aumenta:  106CO2 + 16NO3- + HPO4-2 + 122H2O + 18H+ -> (C106.H263.O110.N16.P1) + 138O2
                                        +SUM(LimnoParam%tNUptNH4MacW(iElem,:))*(50000.*0.001)/14.007/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))                                & !Fotossintese
                                        +SUM(LimnoParam%tPUptMacW(iElem,:))*(50000.*0.001)/31.974/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))                                   & !Fotossintese
                                        -LimnoParam%wNNitrW(iLayer,iElem)*(50000.*0.01)*2./14.007                              & !Nitrificacao   - Diminui:  NH4+ + 2O2 -> NO3- + H2O + 2H+
                                        -LimnoParam%wNDenitW(iLayer,iElem)*(50000.*0.01)*1./14.007                             !Denitrificacao - Aumenta:  5CH2O + 4NO3- + 4H+ -> 5CO2 + 2N2 + 7H2O
            
           LimnoParam%sAlkW(iLayer,iElem) = MAX(NearZero,(LimnoParam%sAlkW(iLayer,iElem) + DtDay*HydroParam%dVarEst(iLayer,iElem,1)*aAlkW))
        EndDo

    
    EndDo
EndIf
	
Return
End















!IF(TEMPODIA.GT.TEMPODIA0)THEN
!!-------------------------------------------------------------------------------------!
!!                                     PROCESSOS									      !
!!-------------------------------------------------------------------------------------!			
!!-CH4---------------------------------------------------------------------------------!
!! 1. Atmospheric Exchange                                                             !
!! 2. Metanotrofia                                                                     !
!!-------------------------------------------------------------------------------------!
!
!        IF (KAUX==1) THEN
!
!	        ! Schmidt, Sc:
!		    ScCO2 = 1911.1 - 118.11*TEMPI(I,J,KAUX) + 3.4527*TEMPI(I,J,KAUX)**2. - 0.043132*TEMPI(I,J,KAUX)**3.
!            
!            ! Gas transfer velocity, kCO2 (cm/hr):
!	        !kCO2 = 0.31 u^2 (Sc/600)^-0.5
!		    kCO2Exch = 0.31 * (VENTOX(int((TEMPODIA-TEMPODIA0)/DTDIA)+1)**2.          &
!		                    +  VENTOY(int((TEMPODIA-TEMPODIA0)/DTDIA)+1)**2.)         &
!		    	            * (ScCO2/660.)**(-0.5)
!	        ! convert to m/day:
!		    kCO2Exch = kCO2Exch * 24.0 / 100.0
!    		
!	        ! Henry - molCO2/L.atm
!            aTkelvin = TEMPI(I,J,KAUX) + 273.15
!	        kHenryCO2= exp(-58.0931 + 90.5069*(100./aTkelvin) + 22.294*log(aTkelvin/100.))
!    	    
!	        ! Dioxido de Carbono Dissolvido: 
!	        !atm            =         gC/m3    /   molCO2/L.atm*(12gC/molCO2)*(1000L/1m3)
!		    pCO2w(I,J,KAUX) = sH2CO3W(I,J,KAUX)/(kHenry*12.*1000.)
!    		
!	        ! Fluxo de CO2 na interface
!	        ! FCO2 = kCO2 * Ko * (pCO2 - PCO2a):
!	        ! g/m2/day      =  m/day * molCO2/L.atm * atm * (12gC/molCO2) * (1000L/1m3)
!		    tCO2Reaer(KAUX) = kCO2Exch*kHenryCO2*(pCO2 - pCO2w(I,J,KAUX))* 12.* 1000.
!			
!        ENDIF
!			
!		!-------------------------------------------------------------------------!
!		! 2. Metanotrofia                                                         !
!            aTLimCH42CO2W = cThetaCH42CO2W ** (TEMPI(I,J,KAUX)-20.)
!            aO2LimCH42CO2W = sO2W(I,J,KAUX)/(sO2W(I,J,KAUX)+hO2CH42CO2) 
!            aCH4LimCH42CO2W = sCH4W(I,J,KAUX)/(sCH4W(I,J,KAUX)+hCH4)
!            wCH4toCO2S(KAUX) = kCH42CO2W * aO2LimCH42CO2W * aTLimCH42CO2W * aCH4LimCH42CO2W * sCH4W(I,J,KAUX)
!		
!		!-------------------------------------------------------------------------!
!        ! 3. Difusao Bolhas                                                       !
!        
!ENDIF
!-------------------------------------------------------------------------------------!
!                               BALANCOS ZONA PELAGICA							      !
!-------------------------------------------------------------------------------------!

!            dCH4W = - Fluxos interface ar-agua
!                    - Metonotrofia
!                    + Bubbling losses - tCH4BubS * fracao /PROF
!                    + Respiracao anaerobica

!-------------------------------------------------------------------------------------!
!                                       SOLVER									      !
!-------------------------------------------------------------------------------------!
!IF (K==NCAMAUX(I,J).AND.K.NE.1)THEN
!    CALL solver2 (KAUX,dCH4W,sCH4W(I,J,KAUX),sCH4WP(I,J,KAUX),sCH4WP(I-1,J,KAUX),sCH4WP(I+1,J,KAUX),sCH4WP(I,J-1,KAUX),sCH4WP(I,J+1,KAUX),sCH4WP(I,J,KAUX-1),0,4,uDLoadCH4aux)
!ELSEIF(K==1.AND.K.NE.NCAMAUX(I,J))THEN
!    CALL solver2 (KAUX,dCH4W,sCH4W(I,J,KAUX),sCH4WP(I,J,KAUX),sCH4WP(I-1,J,KAUX),sCH4WP(I+1,J,KAUX),sCH4WP(I,J-1,KAUX),sCH4WP(I,J+1,KAUX),0,sCH4WP(I,J,KAUX+1),4,uDLoadCH4aux)
!ELSEIF(K==NCAMAUX(I,J).AND.K==1)THEN
!    CALL solver2 (KAUX,dCH4W,sCH4W(I,J,KAUX),sCH4WP(I,J,KAUX),sCH4WP(I-1,J,KAUX),sCH4WP(I+1,J,KAUX),sCH4WP(I,J-1,KAUX),sCH4WP(I,J+1,KAUX),0,0,4,uDLoadCH4aux)
!ELSE
!    CALL solver2 (KAUX,dCH4W,sCH4W(I,J,KAUX),sCH4WP(I,J,KAUX),sCH4WP(I-1,J,KAUX),sCH4WP(I+1,J,KAUX),sCH4WP(I,J-1,KAUX),sCH4WP(I,J+1,KAUX),sCH4WP(I,J,KAUX-1),sCH4WP(I,J,KAUX+1),4,uDLoadCH4aux)
!ENDIF
!-------------------------------------------------------------------------------------!	