Subroutine VerticalWB(HydroParam,MeshParam,MeteoParam,dt,SimTime)

    ! This routine calculates the vertical water balance on the free-surface water by Evaporation and Precipitation
    ! Called in routine 0-MAIN
    !$ use omp_lib
    Use MeshVars 
    Use Hydrodynamic
    Use Meteorological
    
	Implicit none
    
	Real:: Patm0,Temp0,cp,Lv,Cw,Patm,Lambda,GAMA,e0,DELTA,Q,Es,fu_evap,es_air,e_evap,Ea
    Real:: e_water,Sigma,Fsw,DayOfTheYear,Cos_hss,SunSetHourAngle,DayPhotoFration,fc
    Real:: W_evap,Epslon_s,Epslon_c,Epslon,Fir,H_evap, Er, e_sat, B_EVAP
    Integer:: iElem
    Real:: NearZero = 1e-10
    Real:: dt,SimTime,etaplus0
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam

    etaplus0 = 0.d0
    If(SimTime <= 5400 ) Then !Bench 02 Superficial
        etaplus0 = 10.8d0/1000/3600*dt
    EndIf
    
    !etaplus0 = 0.d0
    !If(SimTime <= 3024000 ) Then !Bench 02 Subsurface
    !    etaplus0 = 10.8d0/1000/3600*dt
    !EndIf
    
    !!
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam,MeteoParam) private(iElem,Er,e_sat,B_Evap,e_evap,Ea,DELTA,NearZero,dt)
    Do iElem=1,MeshParam%nElem
        HydroParam%etaplus(iElem) = etaplus0
        
        !If(MeshParam%xb(iElem) == 810.0d0) Then  !Bench02
        !    HydroParam%etaplus(iElem) = etaplus0
        !EndIf 
      !  If (MeteoParam%iReadEvap==0) Then
      !      HydroParam%etaplus(iElem) = -MeteoParam%Evap(iElem)/86400./1000.
      !  ElseIf (MeteoParam%iReadEvap==1) Then
	     !   HydroParam%etaplus(iElem) = -MeteoParam%Evap(iElem)/86400./1000.
      !  ElseIf (MeteoParam%iReadEvap==2) Then     !Requires AirTemp, SolarRad, Wind, Humidity   
      !       ! Energy Budget
      !       Er = MeteoParam%SolarRad(iElem)*86.4*1.E6/( ( 2.501*1.E6-2370*MeteoParam%AirTemp(iElem) ) * 977 )
      !       ! Aerodinamic Budget
      !       e_sat = 611.*EXP( ( 17.27*MeteoParam%AirTemp(iElem) )/( 273.3 + MeteoParam%AirTemp(iElem) ) )
      !       B_Evap = 0.102*( sqrt(MeteoParam%WindX(iElem)**2.+MeteoParam%WindY(iElem)**2.) )/( log(2./0.04) )**2.
      !       e_evap = MeteoParam%RelHum(iElem)/100 * e_sat
      !       Ea = B_Evap*( e_sat - e_evap )
      !       ! DELTA 
      !       DELTA = 4098*e_sat/( 273.3 + MeteoParam%AirTemp(iElem) )**2.
      !       HydroParam%etaplus(iElem) = -MAX(0.0, ( ( DELTA/( DELTA+66.8 ) )*Er + ( 66.8/( DELTA + 66.8 ) )*Ea )/86400./1000.)         
      !  EndIf
	     !
      !  If (MeteoParam%iReadRail==0) Then       ! Constant Value
		    !HydroParam%etaplus(iElem) = HydroParam%etaplus(iElem) + MeteoParam%Precip(iElem)/86400./1000.
      !  ElseIf (MeteoParam%iReadRail==1) Then
      !      !etaplus = etaplus + Precipmm(int((TEMPODIA-TEMPODIA0)/DTDIA)+1)
      !  EndIf
      !
      !  HydroParam%etaplus(iElem) = HydroParam%etaplus(iElem)*dt
      !  
      !If (HydroParam%eta(iElem) - HydroParam%hb(iElem) <= HydroParam%Pcri+NearZero) HydroParam%etaplus(iElem) = 0.
    
    EndDo
    !!$OMP end parallel do

    Return
End Subroutine VerticalWB
	
