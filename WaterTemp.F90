Subroutine WaterTemp(HydroParam,MeshParam,MeteoParam,LimnoParam,dt,dtday)
    
    ! This routine calculates the Water Temperature based on heat-budget
    ! Called in routine 0-MAIN
    ! Call the following routine: UpWindCWC or LTS

    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars
    Use Meteorological

    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam
    type(LimnologyParam) :: LimnoParam
    Integer:: index,iElem,iLayer,count,iRange
    Real:: eair,dt,dtday
    Real, Dimension(MeshParam%KMax,MeshParam%nELem):: Source
    
    
    ! Water Temperature Option -> ReadTemp = 1 --> Water Temperature Model
    If (LimnoParam%iTempW==0) Then
        LimnoParam%sDTempW = HydroParam%WtempRef
    ElseIf (LimnoParam%iTempW==1) Then
        index = 1              ! Boundary Condition Index for Water Temperature (If there is not Boundary condition index = 0)
        LimnoParam%sDTempWP = LimnoParam%sDTempW      ! Save the previous Water Temperature Field
        ! 1. Water Temperature Boundary Condition
        HydroParam%uLoadVarEst = LimnoParam%uDLoadTemp
        
        Call WTemp(HydroParam,MeshParam,MeteoParam,LimnoParam,Source)
            
        HydroParam%dVarEst(:,:,1) = (dt/dtday)*Source
        LimnoParam%iTranspFlag=0
        ! 3. Solver Transport equation    
        If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
            Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDTempW,LimnoParam%sDTempWP,dt,dtday,HydroParam,MeshParam) 
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
            Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDTempW,LimnoParam%sDTempWP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
        ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
            Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDTempW,LimnoParam%sDTempWP,dt,dtday,HydroParam,MeshParam)
        ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
            Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDTempW,LimnoParam%sDTempWP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
        EndiF        
    EndIf
      
    !If (simParam%TIME == simParam%IniTime) Then
   ! Open(90,FILE=trim('Temp-type')//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
    !EndIf
            
    !< Saving position on text file
    !Do i=1,PartParam%nPart
            
    !Write(90,'(I10,100F30.20)') 1,HydroParam%dVarEst(HydroParam%ElCapitalM(12884),12884,1),LimnoParam%sDTempW(HydroParam%ElCapitalM(12884),12884) !HydroParam%ElSmallM(icell):HydroParam%ElSmallM(icell):
    
    	Return
End Subroutine WaterTemp
    
Subroutine WTemp(HydroParam,MeshParam,MeteoParam,LimnoParam,Source)
    
    ! Compute Water Temperature Source Term
    
    ! Based on: 
    ! [1] Hodges, Ben. Heat Budget and Thermodynamics at a free surface: Some theory and numerical implementation.
    !   Centre for Water Research, The University of Western Australia. ED 1300 BH. 1999
    ! [2] Chapra, Steve. Surface Water Quality Modeling.
    !   McGraw-Hill, 1997. ISBN 0-07-011364-5.
    ! [3] Ji, Zhen-Gang. Hydrodynamics and Water Quality: Modeling Rivers, Lakes, and Estuaries.
    !   John Wiley & Sons, 2008. ISBN: 978-0-470-13543-3.
    ! [4] SEML.m, a MatLab code from Sally MacIntyre. Available at Marice Science Institute, UCSB.
    
    ! Input:
    ! SolarRad -> Solar Radiation (W/m²)
    ! AirTemp  -> Air Temperature (°C)
    ! Humidity -> Relative Humidity (%)
    ! Wind     -> Wind Intensity (m/s)
    ! sDTempWP -> Water Temperature (°C)
    ! sDRhoWP  -> Water Density (Kg/m³)
    
    ! Output:
    ! Source -> Water Temperature Source Term (°C/s)
    
    ! List of Modifications:
    !   19.12.2016: Routine Implementation       (J. Rafael Cavalcanti)
    !   20.05.2017: Longwave Radiation from File (J. Rafael Cavalcanti)
    !   14.08.2017: SEML Heat Budget             (J. Rafael Cavalcanti)
    !   13.08.2017: Dynamic Bulk Transfer        (J. Rafael Cavalcanti)
    
    ! Programmer: J. Rafael Cavalcanti
    
    Use Hydrodynamic
    Use MeshVars
    Use LimnologyVars
    Use Meteorological
    !Use Environmental, Only: SolarRad, LWin, AirTemp, Humidity, Wind, ReadSolarRad, AirPressure, CloudCover
    !Use Processes, Only: ShortWave, AtmLongWave, WtrLongWave, Conduction, Evaporation
    !Use Processes, Only: CHE, LMO, WTemp_A, WTemp_e, WTemp_Cp, WExtCoef
    !Use Transport, Only: sDRhoWP
    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    type(MeteorologicalParam) :: MeteoParam
    Integer:: iElem, iLayer
    Integer, Parameter:: iHeat = 2
    Real, Dimension(MeshParam%KMax,MeshParam%nELem):: Source
    Real, Dimension(MeshParam%KMax+1,MeshParam%nELem):: QSw, QS
    Real, Dimension(MeshParam%nELem):: ShortWave,AtmLongWave, WtrLongWave, Conduction, Evaporation,CHE,LMO
    Real:: StefanBoltz       ! Stefan-Boltzmann Constant (W/m²/K)
    Real, Parameter:: Hum100 = 100. 
    Real, Parameter:: VonKarman = 0.4 
    Real, Parameter:: Cp_air = 1004.0                 ! Specific Heat of Air at 25°C (J/kg/°C)
    Real, Parameter:: Ch = 1.4*1.E-3                  ! Bulk Transfer Coefficient of Sensible Heat Transfer
    Real, Parameter:: Cw = 1.4*1.E-3                  ! Bulk Transfer Coefficient of Latent Heat Transfer
    Real, Parameter:: CDN10m  = 1.30E-3               ! Drag Coefficient at 10m
    Real, Parameter:: CHEN10m = 1.35E-3               ! Bulk Transfer Coefficient at 10m
    Real:: eair, esat, FWind, V, A_Wind, B_Wind
    Real:: RhoAir, SpecHum, LatHeat, AirEmiss
    Real:: SphAir, SphWtr, AirDens
    Real:: CDN, CHEN, CD, Tv, WS
    !Real:: CHE
    
    ! ----- Chapra's Heat Budget: -----
    !A_Wind = 6.9; B_Wind = 0.345 ! [3]
    A_Wind = 19.0; B_Wind = 0.95 ! [2]
    StefanBoltz = 0.484490740740741*LimnoParam%tau_param !5.67*1.E-8

    ! -----
    Source = 0.; QSw = 0.; QS = 0.; CHE = 0.; LMO = 0.
    ! 1. Computing the Heat Fluxes
    Do iElem = 1,MeshParam%nElem
        ! 1.1 Absorbed LongWave Radiation - Atmospheric
        !If (ReadSolarRad == 2) Then
        !    AtmLongWave(iElem) = LWin(iElem) !*(1 - Albedo)
        !Else
        AtmLongWave(iElem) = StefanBoltz*(1 - HydroParam%ALB)*( LimnoParam%a_param + 0.031*Sqrt( eair(MeteoParam%AirTemp(iElem),MeteoParam%RelHum(iElem)) ) )*( MeteoParam%AirTemp(iElem) + 273.15 )**4.
        !AtmLongWave(iElem) = AirEmiss(MeteoParam%AirTemp(iElem))*StefanBoltz*(1. - HydroParam%ALB)*(1 + 0.17*CloudCover(iElem)**2. )*( MeteoParam%AirTemp(iElem) + 273.15 )**4.
        !AtmLongWave(iElem) = AirEmiss(MeteoParam%AirTemp(iElem))*StefanBoltz*(1. - HydroParam%ALB)*( LimnoParam%sDTempW(HydroParam%ElCapitalM(iElem),iElem) + 273.15 )**4.
        !End If
        ! 1.2 Emitted LongWave Radiation - Water
        WtrLongWave(iElem) = LimnoParam%e_param*StefanBoltz*( MeteoParam%AirTemp(iElem) + 273.15 )**4.
        Select Case (iHeat)
        Case(1)         ! Chapra Based [2,3]
            ! ----- Compute the Wind Velocity Function -----
            FWind = A_Wind + B_Wind*( HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2. )
            ! 1.3 Latent Heat (Evaporation)
            Evaporation(iElem) = ( esat(MeteoParam%AirTemp(iElem)) - eair(MeteoParam%AirTemp(iElem),MeteoParam%RelHum(iElem)) )*FWind
            ! 1.4 Sensible Heat (Conduction)
            Conduction(iElem) = LimnoParam%c1_param*FWind*( LimnoParam%sDTempW(HydroParam%ElCapitalM(iElem),iElem) - MeteoParam%AirTemp(iElem) )
        Case(2)         ! Sally Based [1,4]
            ! ----- Update Mass Transfer Coefficients [1] -----
            CDN  = CDN10m /( 1 - CDN10m **0.5/VonKarman*log(10./2.) )**2.
            CHEN = CHEN10m/( 1 - CHEN10m**0.5/VonKarman*log(10./2.) )**2.
            SphAir  = SpecHum(MeteoParam%AirTemp(iElem),101325*MeteoParam%AtmPressure(iElem),MeteoParam%RelHum(iElem))
            SphWtr  = SpecHum(LimnoParam%sDTempW(HydroParam%ElCapitalM(iElem),iElem),101325*MeteoParam%AtmPressure(iElem),Hum100)
            AirDens = RhoAir(101325*MeteoParam%AtmPressure(iElem),MeteoParam%AirTemp(iElem),SphAir)
            LatHeat = 2.5008*1.E+6 - 2.3*1.E+3*( MeteoParam%AirTemp(iElem) )
            Tv      = ( MeteoParam%AirTemp(iElem) + 273.15 )*( 1 + 0.6078*SphAir )
            WS      = HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.
            Call Hicks(CDN,CHEN,LimnoParam%sDTempW(HydroParam%ElCapitalM(iElem),iElem),HydroParam%sDRhoW(HydroParam%ElCapitalM(iElem),iElem),SphAir,SphWtr,AirDens,LatHeat,Tv,WS,MeteoParam%AirTemp(iElem),HydroParam%g,HydroParam%pi,CD,CHE(iElem),LMO(iElem))
            ! 1.3 Latent Heat (Evaporation)
            Evaporation(iElem) = AirDens*LatHeat*CHE(iElem)*DSqrt( HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2. )*( SphWtr - SphAir )
            ! 1.4 Sensible Heat (Conduction) 
            Conduction(iElem) = AirDens*Cp_air*CHE(iElem)*DSqrt( HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2. )*( LimnoParam%sDTempW(HydroParam%ElCapitalM(iElem),iElem) - MeteoParam%AirTemp(iElem) )
        Case Default
            Print*, '**************************************************************'
            Print*, 'Something Wrong...'
            Print*, 'UnIPHECO could not Recognize the Heat Balance Option...'
            Print*, 'It is not your fault. Talk to the development team to solve it.'
            Pause
            Stop
        End Select
        ! 1.5 ShortWave Radiation
        ShortWave(iElem) = MeteoParam%SolarRad(iElem)*(1. - HydroParam%ALB)
        ! 1.6 Heat Balance and Source Term
        ! Obs.:  Exponential of zero is equal to 1
        QSw(HydroParam%ElCapitalM(iElem)+1,iElem) = ShortWave(iElem)
        QS(HydroParam%ElCapitalM(iElem)+1,iElem)  = AtmLongWave(iElem) - WtrLongWave(iElem) - Evaporation(iElem) - Conduction(iElem)
        Do iLayer = HydroParam%ElCapitalM(iElem), HydroParam%ElSmallm(iElem), - 1
            QSw(iLayer, iElem) = QSw(iLayer+1,iElem)*EXP( -LimnoParam%aExtCoef(iLayer,iElem)*HydroParam%Dzit(iLayer,iElem) )
            !QS(iLayer, iElem)  = QS (iLayer+1,iElem)*EXP( -WExtCoef(iLayer,iElem)*HydroParam%Dzin(iLayer,iElem) )
            If ( V(HydroParam%etan(iElem),HydroParam%hb(iElem)) <= HydroParam%PCri ) Then      ! The Cell is Dry
                Source(iLayer,iElem) = 0.
            Else
                If (HydroParam%ElSmallm(iElem)==HydroParam%ElCapitalM(iElem)) Then
                    Source(iLayer,iElem) = (QSw(iLayer+1,iElem) + QS(iLayer+1,iElem))/( HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param*Min(1.,Max(1.,HydroParam%Dzit(iLayer,iElem))) )
                Else
                    If (iLayer == HydroParam%ElCapitalM(iElem)) Then
                        Source(iLayer,iElem) = ( ( QSw(iLayer+1,iElem) - QSw(iLayer,iElem) ) + ( QS(iLayer+1,iElem) - QS(iLayer,iElem) ) )/( HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param*Min(1.,Max(1.,HydroParam%Dzit(iLayer,iElem))) )
                    Else
                        Source(iLayer,iElem) = ( ( QSw(iLayer+1,iElem) - QSw(iLayer,iElem) ) + ( QS(iLayer+1,iElem) - QS(iLayer,iElem) ) )/( HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param*HydroParam%Dzit(iLayer,iElem) )
                    EndIf
                EndIf
            End If
        EndDo       ! Layer Loop
    End Do       ! Element Loop
        
    Return
End Subroutine WTemp
    
! ----------------------------------------------------------------------------------------    
! ---------------------------------- SUBROUTINES -----------------------------------------
! ----------------------------------------------------------------------------------------
Subroutine Hicks(CDN,CHEN,SfrT,SfrRho,SphAir,SphWtr,AirDens,LatHeat,Tv,WS,AirT,g,pi,CD,CHE,LMO)
    
    ! Update the Bulk Transfer Coefficients based on an iterative Procedure
    
    ! Based on: 
    ! [1] Hicks,B.B. A procedure for the formulation of bulk transfer coefficients over water.
    !   Boundary-Layer Meteorology, vol. 8, pp. 151-524, 1975.
    ! [2] Paulson,C.A. The mathematical representation of Wind Speed and Temperature profiles in the Unstable atmospheric surface layer.
    !   Journal of Applied Meteorology, vol. 9, pp. 857-861, 1970.
    ! [3] SEML.m, a MatLab code from Sally MacIntyre. Available at Marice Science Institute, UCSB.
    ! [4] MacIntyre,S.; Romero,J.R.; Kling,G.W. Spatial-temporal variability in surface layer deepening and lateral advection in an embayment of lake Victoria, East Africa.
    !   Limnology and Oceanography, vol. 47(3), pp. 656-671, 2002.

    ! Input:
    ! CDN  -> Momentum Drag in Neutral Atmospheric Conditions
    ! CHEN -> Bulk Transfer Coefficient in Neutral Atmospheric Conditions
    ! SfrT -> Surface Temperature

    ! Output:
    ! CD  -> Momentum Drag in Unstable Atmospheric Conditions
    ! CHE -> Bulk Transfer Coefficient in Unstable Atmospheric Conditions
    ! LMO -> Monin-Obukhov Length Scale

    ! List of Modifications:
    !   13.08.2017: Routine Implementation       (J. Rafael Cavalcanti)
    
    ! Programmer: J. Rafael Cavalcanti

    Implicit None
    Real, Intent(In):: CDN, CHEN, SfrRho, SfrT, SphAir, SphWtr, AirDens, LatHeat, Tv, WS, AirT
    Real, Intent(Out):: CD, CHE, LMO
    Real, Parameter:: Cp_air = 1004.0               ! Specific Heat of Air at 25°C (J/kg/°C)
    Real, Parameter:: Lcnst  = 0.61                 ! Monin-Obukhov Constant
    Real, Parameter:: Hconst  = 2.0                 !
    Logical:: Init, Cont
    Integer:: Count
    Real:: Tau, ustrAir, ustrWtr, SE, LE, OldL
    Real, Parameter:: VonKarman = 0.4
    Real::g,pi
        
    ! 1. Initialize Variables
    Init = .TRUE.; Cont = .TRUE.
    CHE = CHEN; CD = CDN
    Count = 0
    ! 2. Begin Hicks Iteration 
    Do While (Cont == .TRUE.)
        Tau     = AirDens*CD*WS
        ustrAir = Sqrt( Tau/AirDens )
        ustrWtr = Sqrt( Tau/SfrRho )
        ! ----- According to [1] -----
        SE      = AirDens*Cp_air *CHE*WS*( AirT - SfrT )
        LE      = AirDens*LatHeat*CHE*WS*( SphAir - SphWtr )        
        LMO     = AirDens*(ustrAir**3.)*Tv/VonKarman/g
        If ( SE == 0. .AND. LE == 0. ) Then
            LMO = LMO
        Else
            LMO = LMO/( ( SE/Cp_air ) + ( Lcnst*( AirT + 273.15 )*LE/LatHeat ) )      ! Monin-Obukhov Length
        End If      ! If ( SE == 0. .AND. LE == 0. ) Then
        Call BulkTransfer(Hconst,pi,LMO,CDN,CHEN,CD,CHE)
        ! ----- Check Convergence -----
        If ( Init ) Then
            OldL = LMO
            Init = .FALSE.
        End If      ! Check Init
        Count = Count + 1
        If ( ( OldL - LMO )/OldL > 1.E-5 ) Then
            OldL = LMO
        Else
            Cont = .FALSE.
        End If      ! Check Convergence
    End Do      ! End Hicks Interation
    
    !Print*, 'The Hicks iteration took: ', Count, 'iterations'
    
    Return
End Subroutine Hicks

! ----------------------------------------------------------------------------------------
    
Subroutine BulkTransfer(H,pi,L,CDN,CEN,CD,CE)
    
    ! Compute Bulk Transfer Coefficients
    
    ! Based on: 
    ! [1] Paulson,C.A. The mathematical representation of Wind Speed and Temperature profiles in the Unstable atmospheric surface layer.
    !   Journal of Applied Meteorology, vol. 9, pp. 857-861, 1970.
    ! [2] SEML.m, a MatLab code from Sally MacIntyre. Available at Marice Science Institute, UCSB.

    ! Input:
    ! H   -> Sensor Height
    ! LMO -> Monin-Obukhov Length Scale
    ! CDN -> Momentum Drag in Neutral Atmospheric Conditions
    ! CEN -> Bulk Transfer Coefficient in Neutral Atmospheric Conditions

    ! Output:
    ! CD -> Momentum Drag in Unstable Atmospheric Conditions
    ! CE -> Bulk Transfer Coefficient in Unstable Atmospheric Conditions

    ! List of Modifications:
    !   13.08.2017: Routine Implementation       (J. Rafael Cavalcanti)
    
    ! Programmer: J. Rafael Cavalcanti
    
    Implicit None
    Real, Intent(In):: H, CDN, CEN
    Real, Intent(Out):: CD, CE
    Real:: sm, se, x, L
    Real, Parameter:: VonKarman = 0.4
    Real:: pi
    
    ! 1. Compute Integrated Similarity Functions
    If ( (H/L) > 15. ) Then
        L = H/15.
    Else
        Continue
    End If      ! If ( (H/L) > 15. ) Then
    If ( (H/L) < -15. ) Then
        L = -H/15.
    Else
        Continue
    End If      ! If ( (H/L) < -15. ) Then
    If ( (H/L) < 0. ) Then
        x = ( 1 - 16*(H/L) )**0.25
        sm  = 2.*log( ( 1. + x )/2. ) + log( (1. + x**2.)/2. ) - 2.*Atan(x)
        sm  = sm + Pi/2.
        se  = 2.*log( ( 1 + x**2. )/2. )
    ElseIf ( (H/L) >= 0. .AND. (H/L) < 0.5 ) Then
        sm = -5.*(H/L)
        se = sm
    ElseIf ( (H/L) >= 0.5 .AND. (H/L) < 10.0 ) Then
        sm = 0.5/(H/L)**2. - 4.25/(H/L) - 7.*log( H/L ) - 0.852
        se = sm
    Else
        sm = log( H/L ) - 0.76*( H/L ) - 12.093
        se = sm;
    End If      ! If ( (H/L) < 0. ) Then
    
    CD = CDN/( 1 + CDN/VonKarman**2.*( sm*sm - VonKarman*sm/CDN**0.5 - VonKarman*sm*CDN**0.5/CDN ) )
    CE = CEN/( 1 + CEN/VonKarman**2.*( sm*se - VonKarman*se/CDN**0.5 - VonKarman*sm*CDN**0.5/CEN ) ) 
    If ( isNaN(CD) ) CD = CDN
    If ( isNaN(CE) ) CE = CEN
    
    Return
End Subroutine BulkTransfer    
