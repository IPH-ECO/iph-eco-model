Subroutine UpdateLimnVars(HydroParam,MeshParam,LimnoParam,simParam)

    ! This routine updates the Past-Time value of Water Quality variables with the Real-Time value
    ! Called in routines: 
    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars
    Use SimulationModel

    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    type(SimulationParam) :: simParam
    Integer:: index,iElem,iLayer,om,gg
    Integer:: ano2,diajulcont1,dayoftheyear
    Real:: delta,cos_hss,sunSetHourAngle,dayPhotoFration
    Real:: V,idate(6)
    
    ! Abiotic Variables in Water
    LimnoParam%sO2WP       = LimnoParam%sO2W               ! Dissolve Oxygen
    LimnoParam%sDicWP       = LimnoParam%sDicW
    LimnoParam%sCH4WP      = LimnoParam%sCH4W
    LimnoParam%sDPomWP     = LimnoParam%sDPomW
    LimnoParam%sCPomWP     = LimnoParam%sCPomW
    LimnoParam%sNPomWP     = LimnoParam%sNPomW
    LimnoParam%sPPomWP     = LimnoParam%sPPomW
    LimnoParam%sDDomWP     = LimnoParam%sDDomW
    LimnoParam%sCDomWP     = LimnoParam%sCDomW
    LimnoParam%sNDomWP     = LimnoParam%sNDomW
    LimnoParam%sPDomWP     = LimnoParam%sPDomW
    LimnoParam%sSiDomWP    = LimnoParam%sSiDomW
	LimnoParam%sDIMWP      = LimnoParam%sDIMW              ! Inorganic Matter in Water
	LimnoParam%sPO4WP     = LimnoParam%sPO4W              ! PO4 in Water
	LimnoParam%sNO3WP     = LimnoParam%sNO3W              ! NO3 in Water
	LimnoParam%sNH4WP     = LimnoParam%sNH4W              ! NH4 in Water
	LimnoParam%sSiO2WP    = LimnoParam%sSiO2W             ! Silica in Water
	LimnoParam%sPAIMWP    = LimnoParam%sPAIMW             ! Adsorved P in Inorganic Matter in Water
    LimnoParam%sDBacWP    = LimnoParam%sDBacW
    LimnoParam%sCBacWP    = LimnoParam%sCBacW
    LimnoParam%sNBacWP    = LimnoParam%sNBacW
    LimnoParam%sPBacWP    = LimnoParam%sPBacW
    LimnoParam%sDPhytWP   = LimnoParam%sDPhytW
    LimnoParam%sCPhytWP   = LimnoParam%sCPhytW
    LimnoParam%sNPhytWP   = LimnoParam%sNPhytW
    LimnoParam%sPPhytWP   = LimnoParam%sPPhytW
	LimnoParam%sDMacP     = LimnoParam%sDMac              ! Total Macrophyte Biomass
    LimnoParam%sCMacP     = LimnoParam%sCMac              ! Total Macrophyte (C Fraction) Biomass
	LimnoParam%sPMacP     = LimnoParam%sPMac              ! Total Macrophyte (P Fraction) Biomass
	LimnoParam%sNMacP     = LimnoParam%sNMac              ! Total Macrophyte (N Fraction) Biomass
	LimnoParam%sDZooP     = LimnoParam%sDZoo              ! Zooplankton Biomass
    LimnoParam%sCZooP     = LimnoParam%sCZoo              ! Zooplankton (P Fraction) Biomass
	LimnoParam%sPZooP     = LimnoParam%sPZoo              ! Zooplankton (P Fraction) Biomass
	LimnoParam%sNZooP     = LimnoParam%sNZoo              ! Zooplankton (N Fraction) Biomass
    LimnoParam%sDFiAdP     = LimnoParam%sDFiAd
    LimnoParam%sCFiAdP     = LimnoParam%sCFiAd
    LimnoParam%sNFiAdP     = LimnoParam%sNFiAd
    LimnoParam%sPFiAdP     = LimnoParam%sPFiAd
    LimnoParam%sDFiJvP     = LimnoParam%sDFiJv
    LimnoParam%sCFiJvP     = LimnoParam%sCFiJv
    LimnoParam%sNFiJvP     = LimnoParam%sNFiJv
    LimnoParam%sPFiJvP     = LimnoParam%sPFiJv

    
    Do iElem = 1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            
            !Organic Phytoplankton  (mgD/l)
			LimnoParam%oDPhytW(iLayer,iElem)=sum(LimnoParam%sDPhytW(iLayer,iElem,:))	                   
			LimnoParam%oCPhytW(iLayer,iElem)=sum(LimnoParam%sCPhytW(iLayer,iElem,:))	                    
			LimnoParam%oNPhytW(iLayer,iElem)=sum(LimnoParam%sNPhytW(iLayer,iElem,:))	                    
			LimnoParam%oPPhytW(iLayer,iElem)=sum(LimnoParam%sPPhytW(iLayer,iElem,:))	                
            
            !Organic seston  (mgD/l)
            LimnoParam%oDOMW(iLayer,iElem) = sum(LimnoParam%sDPomW(iLayer,iElem,:))+LimnoParam%oDPhytW(iLayer,iElem)
            LimnoParam%oCOMW(iLayer,iElem) = sum(LimnoParam%sCPomW(iLayer,iElem,:))+LimnoParam%oCPhytW(iLayer,iElem)
            LimnoParam%oNOMW(iLayer,iElem) = sum(LimnoParam%sNPomW(iLayer,iElem,:))+LimnoParam%oNPhytW(iLayer,iElem)
            LimnoParam%oPOMW(iLayer,iElem) = sum(LimnoParam%sPPomW(iLayer,iElem,:))+LimnoParam%oPPhytW(iLayer,iElem)
            
			LimnoParam%oDSestW(iLayer,iElem)=LimnoParam%oDOMW(iLayer,iElem)+LimnoParam%sDIMW(iLayer,iElem)			        !seston total (mgD/l)
			LimnoParam%oNSestW(iLayer,iElem)=LimnoParam%oNOMW(iLayer,iElem)					                !seston total (mgN/l)
			LimnoParam%oPSestW(iLayer,iElem)=LimnoParam%oPOMW(iLayer,iElem)+LimnoParam%sPAIMW(iLayer,iElem)			        !seston total (mgP/l)
            
            !Nitrogenio
			LimnoParam%oNDissW(iLayer,iElem)=LimnoParam%sNO3W(iLayer,iElem)+LimnoParam%sNH4W(iLayer,iElem)			        !liberação lenta de nitrogenio na água (mgN/l)
			LimnoParam%oNkjW(iLayer,iElem)=LimnoParam%oNSestW(iLayer,iElem)+LimnoParam%sNH4W(iLayer,iElem)		        !N Kjeldahl na agua (mgN/l)
			LimnoParam%oNTotW(iLayer,iElem)=LimnoParam%oNkjW(iLayer,iElem)+LimnoParam%sNO3W(iLayer,iElem)			        !N total na agua-excluindo animais e vegetação (mgP/l)
			
			!Fosforo
			LimnoParam%oPInorgW(iLayer,iElem)=LimnoParam%sPO4W(iLayer,iElem)+LimnoParam%sPAIMW(iLayer,iElem)			        !P inorganico na agua (mgP/l)
			LimnoParam%oPTotW(iLayer,iElem)=LimnoParam%oPSestW(iLayer,iElem)+LimnoParam%sPO4W(iLayer,iElem)			        !P total na agua-excluindo animais e vegetação (mgP/l)            
        EndDo
    EndDo
            
   	! Total Light extintion atenuation coefficient
    Do iElem = 1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            ! Light extintion atenuation coefficient due to IM
            LimnoParam%aExtIM(iLayer,iElem) = LimnoParam%cExtSpIM * LimnoParam%sDIMW(iLayer,iElem)
            ! Light extintion atenuation coefficient due to DOM
            Do om = 1, LimnoParam%npom ! =1 -> Labile; =2 -> Refractory
                LimnoParam%aExtPom(iLayer,iElem,om) = LimnoParam%cExtSpPom(om)* LimnoParam%sDPomW(iLayer,iElem,om)  
                LimnoParam%aExtDom(iLayer,iElem,om) = LimnoParam%cExtSpDom(om)* LimnoParam%sDDomW(iLayer,iElem,om)    
            EndDo
            ! Light extintion atenuation coefficient due to Phytoplankton
            Do gg = 1, LimnoParam%nphy
                LimnoParam%aExtPhyt(iLayer,iElem,gg) = LimnoParam%cExtSpPhyt(gg) * LimnoParam%sDPhytW(iLayer,iElem,gg)
            EndDo
        
            
            ! Total Light extintion atenuation coefficient (1/m)
        	LimnoParam%aExtCoef(iLayer,iElem) = LimnoParam%cExtWat + LimnoParam%aExtIM(iLayer,iElem) + SUM(LimnoParam%aExtDom(iLayer,iElem,:)) + SUM(LimnoParam%aExtPom(iLayer,iElem,:)) + SUM(LimnoParam%aExtPhyt(iLayer,iElem,:)) ! + SUM(aExtZoo(KAUX,:)) + aExtVeg 
            !Light extiction coefficient without vegetation (1/m)
            LimnoParam%aExtCoefOpen(iLayer,iElem) = LimnoParam%cExtWat + LimnoParam%aExtIM(iLayer,iElem) + SUM(LimnoParam%aExtDom(iLayer,iElem,:)) + SUM(LimnoParam%aExtPom(iLayer,iElem,:)) + SUM(LimnoParam%aExtPhyt(iLayer,iElem,:))           !Coeficiente de extinção da luz-sem contribuição da vegetação (1/m)
            !Poole-Atkins coefficient
            LimnoParam%aPACoef = LimnoParam%cPACoefMin + (LimnoParam%cPACoefMax - LimnoParam%cPACoefMin) * LimnoParam%hPACoef / (LimnoParam%hPACoef + LimnoParam%oDOMW(iLayer,iElem))	    
            !Secchi Depth (m)
            LimnoParam%aSecchi(iLayer,iElem) = MIN(V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem)) , LimnoParam%aPACoef / LimnoParam%aExtCoefOpen(iLayer,iElem))							!Profundidade de Secchi (m)
        EndDo
    EndDo

				
	!  Day_length
    Call unix2c(simParam%time, idate)
    !print *,idate
    !Day0*19                                        ! Initial time of simulation (AAAA/MM/DD HH:MM:SS)
	CALL JULIANDAY(floor(idate(1)),1,1,diajulcont1) !January 1st if the current year
    
    CALL JULIANDAY(floor(idate(1)),floor(idate(2)),floor(idate(3)),dayoftheyear)
    simParam%Simday = dayoftheyear
	dayoftheyear = int(dayoftheyear-diajulcont1)+1
                
	IF (dayoftheyear<1.OR.dayoftheyear>366) THEN
        write(*,*) ' ERROR: Day of year is out of bounds.'
        write(*,*) ' dayoftheyear =',dayoftheyear
        STOP 
    ENDIF
    simParam%Julday = dayoftheyear
    
				
	delta = 23.45*(HydroParam%PI/180)*COS(2*HydroParam%PI/365*(172-dayoftheyear))
	cos_hss = -tan(HydroParam%LAT*HydroParam%PI/180.)*tan(delta)
                
	IF (abs(cos_hss) > 1.0) THEN
        write(*,*) ' ERROR: cos_hss (cosine of the sunset hour angle) is out of range.'
        write(*,*) ' cos_hss =',cos_hss
        STOP 
    ENDIF
					    
	sunSetHourAngle = acos(cos_hss)
	dayPhotoFration = sunSetHourAngle/HydroParam%PI
	LimnoParam%ufDay = dayPhotoFration  !cfDayAve-cfDayVar*(-cos(2*PI*(DiaJul+TempLag)/365.)) !Fotoperíodo (comprimento do dia) (h/24h)
    
    Call Ratios(HydroParam,MeshParam,LimnoParam)
				
	
	Return
    
End Subroutine UpdateLimnVars