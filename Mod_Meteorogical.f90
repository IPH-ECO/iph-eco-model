!> This module declares the simulation parameters. 
Module Meteorological
    
    use domain_types
    
    Implicit None
    
     type MeteorologicalParam

        ! Meterological parameters
	    Real TempMax					!<Temperatura máxima do lago (oC)
	    Real TempMin					!<Temperatura minima do lago (oC)
	    Real TempLag					!<Dia Juliano que ocorre a temperatura minima (Dia Juliano)
	    Real TempDev					!<fator de desvio da temperatura (-)
	    Real MaxI						!<Irradiação solar máxima do lago (W/m²/d)
	    Real MinI						!<Irradiação solar minima do lago (W/m²/d)
	    Real ILag						!<Dia Juliano que ocorre a Irradiação solar minima (Dia Juliano)
	
	    Real ConsRad					!<Solar radiation constant  read from simulation preference 
        Integer:: Nstation                          !<Number of meteorological stations 
    
        Real, Allocatable:: rhoair(:) !<Air density (kg/m³)
    
        Real, Allocatable:: CoordMetStations(:,:)
        Real, Allocatable:: ValueMetStations(:)
	    !Atmospheric Pressure
        Integer, Allocatable:: AtmPressurenTime(:)
        Real, Allocatable:: AtmPressureTime(:,:)
        Real, Allocatable:: AtmPressureValue(:,:)
        Real, Allocatable:: AtmPressure(:)

        !Air temperature
        Integer, Allocatable:: AirTempnTime(:)
        Real, Allocatable:: AirTempTime(:,:)
        Real, Allocatable:: AirTempValue(:,:)
        Real, Allocatable:: AirTemp(:)
    
        !Solar radiation
        Integer, Allocatable:: SolarRadnTime(:)
        Real, Allocatable:: SolarRadTime(:,:)
        Real, Allocatable:: SolarRadValue(:,:)
        Real, Allocatable:: SolarRad(:)
    
        !Relative Humidity
        Integer, Allocatable:: RelHumnTime(:)
        Real, Allocatable:: RelHumTime(:,:)
        Real, Allocatable:: RelHumValue(:,:)
        Real, Allocatable:: RelHum(:)

        !Wind
        Integer, Allocatable:: WindnTime(:)
        Real, Allocatable:: WindTime(:,:)
        Real, Allocatable:: WindValuex(:,:)
        Real, Allocatable:: WindValuey(:,:)
        Real, Allocatable:: WindX(:),WindY(:)
        !Precipitation
        Integer, Allocatable:: PrecipnTime(:)
        Real, Allocatable:: PrecipTime(:,:)
        Real, Allocatable:: PrecipValue(:,:)
        Real, Allocatable:: Precip(:)
    
        !Evaporation
        Integer, Allocatable:: EvapnTime(:)
        Real, Allocatable:: EvapTime(:,:)
        Real, Allocatable:: EvapValue(:,:)
        Real, Allocatable:: Evap(:)
    
    
        !Flags on simulation 
	    INTEGER iReadAtmP       !< Atmospheric pressure option: 0 = Constant Temp; 1 = time-series
        INTEGER iReadAirTemp    !< Air Temperature option: 0 = Constant Temp; 1 = time-series
        INTEGER iReadHum        !< Relative Humidity option: 0 = Constant Temp; 1 = time-series
	    INTEGER iReadRad        !< Solar radiation option: 0 = Constant Temp; 1 = time-series
	    INTEGER iReadWind       !< Wind option: 0 = Constant Temp; 1 = time-series (X and Y components); 2 = time-series (intensity and direction)
	    INTEGER iReadRail       !< Precipitation option: 0 = Constant Temp; 1 = time-series
	    INTEGER iReadEvap       !< Evaporation option: 0 = Constant Temp; 1 = time-series; 2 = Modelling
    
    contains
        procedure :: initializeMeteoParam
    end type
    
    contains
    
    subroutine initializeMeteoParam(this,MeteorConfiguration,IniTime,FinalTime,AirtempRef,pi)
    
        Integer:: NObj,nTime, i,j,k
        !Integer:: i,j,iElem,iEdge,iLayer,r,l,Face
        !character, pointer :: hydroParametersName(:)
        class(MeteorologicalParam) :: this
        type(MeteorologicalConfiguration) :: MeteorConfiguration
        type(MeteorologicalStation), pointer :: MeteorStation(:)
        type(MeteorologicalParameter), pointer :: MeteorParameter(:)
        character, pointer :: MeteorParametersName(:)
        type(TimeSeries), pointer :: MeteortimeSeries(:)
        character(len=200):: text
        Real:: idate(6)
        Integer:: IniTime !< Initial time of simulation (unix time)
        Integer:: FinalTime !< Final time of simulation (unix time)
        Real::AirtempRef,pi
    
    
        call c_f_pointer(MeteorConfiguration%stations, MeteorStation,[MeteorConfiguration%stationsLength])
    
    
        this%Nstation = MeteorConfiguration%stationsLength 
    
        Allocate (this%CoordMetStations(2,this%Nstation))
        Allocate (this%ValueMetStations(this%Nstation))
    
        If (this%Nstation==0) Then
            Return
        Else
            Do i = 1, this%Nstation
                this%CoordMetStations(1,i) = MeteorStation(i)%utmX
                this%CoordMetStations(2,i) = MeteorStation(i)%utmY
                call c_f_pointer(MeteorStation(i)%parameters, MeteorParameter, [MeteorStation(i)%parametersLength])
                Do j = 1, MeteorStation(i)%parametersLength
                    text =''
                    call c_f_pointer(MeteorParameter(j)%name, MeteorParametersName, [MeteorParameter(j)%nameLength])
            
                    Do k=1,MeteorParameter(j)%nameLength
                        text = trim(text)//MeteorParametersName(k)
                    EndDo
            
                    If (trim(text) == 'AtmosphericPressure') Then 
                        If (MeteorParameter(j)%functionMet == 1) Then ! Constant
                            Allocate (this%AtmPressureValue(i,2))
                            Allocate (this%AtmPressureTime(i,2))
                            Allocate (this%AtmPressurenTime(i))
                            this%AtmPressurenTime(i) = 2
                            this%AtmPressureTime(i,1) = IniTime
                            this%AtmPressureTime(i,2) = FinalTime
		                    this%AtmPressureValue(i,:) = MeteorParameter(j)%constantValue
                        ElseIf (MeteorParameter(j)%functionMet == 2) Then ! Time-series
                            nTime = MeteorParameter(j)%timeSizeListLength
                            If (nTime<2) Then
                                print *,'Number of time intervals is lesser then two. Please check ' //trim(text) // 'time series.'// 'Station '// char(i)//'.' ! 
                                pause
                                Stop 
                            EndIf
                            Allocate (this%AtmPressureValue(i,nTime))
                            Allocate (this%AtmPressureTime(i,nTime))
                            Allocate (this%AtmPressurenTime(i))
                            this%AtmPressurenTime(i) = nTime
                            call c_f_pointer(MeteorParameter(j)%timeSeriesList, MeteortimeSeries, [MeteorParameter(j)%timeSizeListLength])
                            Do k = 1, this%AtmPressurenTime(i)
                                ! Storing stamp time and values of time series 
                                this%AtmPressureTime(i,k) = MeteortimeSeries(k)%timeStamp
                                Call unix2c(this%AirTempTime(i,k), idate)
                                print *,idate
                                this%AtmPressureValue(i,k) = MeteortimeSeries(k)%value1
                            EndDo
                        EndIf

                    ElseIf (trim(text) == 'AirTemperature') Then
                        If (MeteorParameter(j)%functionMet == 1) Then ! Constant
                            Allocate (this%AirTempValue(i,2))
                            Allocate (this%AirTempTime(i,2))
                            Allocate (this%AirTempnTime(i))
                            this%AirTempnTime(i) = 2
                            this%AirTempTime(i,1) = IniTime
                            this%AirTempTime(i,2) = FinalTime
		                    this%AirTempValue(i,:) = MeteorParameter(j)%constantValue
                        ElseIf (MeteorParameter(j)%functionMet == 2) Then ! Time-series
                            nTime = MeteorParameter(j)%timeSizeListLength
                            If (nTime<2) Then
                                print *,'Number of time intervals is lesser then two. Please check ' //trim(text) // 'time series.'// 'Station '// char(i)//'.' ! 
                                pause
                                Stop 
                            EndIf
                            Allocate (this%AirTempValue(i,nTime))
                            Allocate (this%AirTempTime(i,nTime))
                            Allocate (this%AirTempnTime(i))
                            this%AirTempnTime(i) = nTime
                            call c_f_pointer(MeteorParameter(j)%timeSeriesList, MeteortimeSeries, [MeteorParameter(j)%timeSizeListLength])
                            Do k = 1, this%AirTempnTime(i)
                                ! Storing stamp time and values of time series 
                                this%AirTempTime(i,k) = MeteortimeSeries(k)%timeStamp
                                this%AirTempValue(i,k) = MeteortimeSeries(k)%value1
                            EndDo
                        EndIf
                    
                    ElseIf (trim(text) == 'RelativeHumidity') Then
                        If (MeteorParameter(j)%functionMet == 1) Then ! Constant
                            Allocate (this%RelHumValue(i,2))
                            Allocate (this%RelHumTime(i,2))
                            Allocate (this%RelHumnTime(i))
                            this%RelHumnTime(i) = 2
                            this%RelHumTime(i,1) = IniTime
                            this%RelHumTime(i,2) = FinalTime
		                    this%RelHumValue(i,:) = MeteorParameter(j)%constantValue
                        ElseIf (MeteorParameter(j)%functionMet == 2) Then ! Time-series
                            nTime = MeteorParameter(j)%timeSizeListLength
                            If (nTime<2) Then
                                print *,'Number of time intervals is lesser then two. Please check ' //trim(text) // 'time series.'// 'Station '// char(i)//'.' ! 
                                pause
                                Stop 
                            EndIf
                            Allocate (this%RelHumValue(i,nTime))
                            Allocate (this%RelHumTime(i,nTime))
                            Allocate (this%RelHumnTime(i))
                            this%RelHumnTime(i) = nTime
                            call c_f_pointer(MeteorParameter(j)%timeSeriesList, MeteortimeSeries, [MeteorParameter(j)%timeSizeListLength])
                            Do k = 1, this%RelHumnTime(i)
                                ! Storing stamp time and values of time series 
                                this%RelHumTime(i,k) = MeteortimeSeries(k)%timeStamp
                                this%RelHumValue(i,k) = MeteortimeSeries(k)%value1
                            EndDo
                        EndIf
                
                    ElseIf (trim(text) == 'SolarRadiation') Then
                        If (MeteorParameter(j)%functionMet == 1) Then ! Constant
                            Allocate (this%SolarRadValue(i,2))
                            Allocate (this%SolarRadTime(i,2))
                            Allocate (this%SolarRadnTime(i))
                            this%SolarRadnTime(i) = 2
                            this%SolarRadTime(i,1) = IniTime
                            this%SolarRadTime(i,2) = FinalTime
		                    this%SolarRadValue(i,:) = MeteorParameter(j)%constantValue
                        ElseIf (MeteorParameter(j)%functionMet == 2) Then ! Time-series
                            nTime = MeteorParameter(j)%timeSizeListLength
                            If (nTime<2) Then
                                print *,'Number of time intervals is lesser then two. Please check ' //trim(text) // 'time series.'// 'Station '// char(i)//'.' ! 
                                pause
                                Stop 
                            EndIf
                            Allocate (this%SolarRadValue(i,nTime))
                            Allocate (this%SolarRadTime(i,nTime))
                            Allocate (this%SolarRadnTime(i))
                            this%SolarRadnTime(i) = nTime
                            call c_f_pointer(MeteorParameter(j)%timeSeriesList, MeteortimeSeries, [MeteorParameter(j)%timeSizeListLength])
                            Do k = 1, this%SolarRadnTime(i)
                                ! Storing stamp time and values of time series 
                                this%SolarRadTime(i,k) = MeteortimeSeries(k)%timeStamp
                                this%SolarRadValue(i,k) = MeteortimeSeries(k)%value1
                            EndDo
                        EndIf
                
                    ElseIf (trim(text) == 'Wind') Then
                        If (MeteorParameter(j)%functionMet == 1) Then ! Constant
                            Allocate (this%WindValuex(i,2))
                            Allocate (this%WindValuey(i,2))
                            Allocate (this%WindTime(i,2))
                            Allocate (this%WindnTime(i))
                            this%WindnTime(i) = 2
                            this%WindTime(i,1) = IniTime
                            this%WindTime(i,2) = FinalTime
                            !Precisa converter para componentes x e y se vier em int/dir
                            If (MeteorParameter(j)%useXYComponent) Then
		                        this%WindValuex(i,:) = MeteorParameter(j)%xComponent
                                this%WindValuey(i,:) = MeteorParameter(j)%yComponent
                            Else
		                        this%WindValuex(i,:) = -MeteorParameter(j)%intensity*sin(pi*MeteorParameter(j)%direction/180.)
                                this%WindValuey(i,:) = -MeteorParameter(j)%intensity*cos(pi*MeteorParameter(j)%direction/180.)
                            EndIf
                        
                        ElseIf (MeteorParameter(j)%functionMet == 2) Then ! Time-series
                            nTime = MeteorParameter(j)%timeSizeListLength
                            If (nTime<2) Then
                                print *,'Number of time intervals is lesser then two. Please check ' //trim(text) // 'time series.'// 'Station '// char(i)//'.' ! 
                                pause
                                Stop 
                            EndIf
                            Allocate (this%WindValuex(i,nTime))
                            Allocate (this%WindValuey(i,nTime))
                            Allocate (this%WindTime(i,nTime))
                            Allocate (this%WindnTime(i))
                            this%WindnTime(i) = nTime
                            call c_f_pointer(MeteorParameter(j)%timeSeriesList, MeteortimeSeries, [MeteorParameter(j)%timeSizeListLength])
                            Do k = 1, this%WindnTime(i)
                                ! Storing stamp time and values of time series 
                                this%WindTime(i,k) = MeteortimeSeries(k)%timeStamp
                                If (MeteorParameter(j)%useXYComponent) Then
                                    this%WindValuex(i,k) = MeteortimeSeries(k)%value1
                                    this%WindValuey(i,k) = MeteortimeSeries(k)%value2
                                Else
		                            this%WindValuex(i,k) = -MeteortimeSeries(k)%value1*sin(pi*MeteortimeSeries(k)%value2/180.)
                                    this%WindValuey(i,k) = -MeteortimeSeries(k)%value1*cos(pi*MeteortimeSeries(k)%value2/180.)
                                EndIf
                            EndDo
                        EndIf
                    
                    ElseIf (trim(text) == 'Precipitation') Then
                        If (MeteorParameter(j)%functionMet == 1) Then ! Constant
                            Allocate (this%PrecipValue(i,2))
                            Allocate (this%PrecipTime(i,2))
                            Allocate (this%PrecipnTime(i))
                            this%PrecipnTime(i) = 2
                            this%PrecipTime(i,1) = IniTime
                            this%PrecipTime(i,2) = FinalTime
		                    this%PrecipValue(i,:) = MeteorParameter(j)%constantValue
                        ElseIf (MeteorParameter(j)%functionMet == 2) Then ! Time-series
                            nTime = MeteorParameter(j)%timeSizeListLength
                            If (nTime<2) Then
                                print *,'Number of time intervals is lesser then two. Please check ' //trim(text) // 'time series.'// 'Station '// char(i)//'.' ! 
                                pause
                                Stop 
                            EndIf
                            Allocate (this%PrecipValue(i,nTime))
                            Allocate (this%PrecipTime(i,nTime))
                            Allocate (this%PrecipnTime(i))
                            this%PrecipnTime(i) = nTime
                            call c_f_pointer(MeteorParameter(j)%timeSeriesList, MeteortimeSeries, [MeteorParameter(j)%timeSizeListLength])
                            Do k = 1, this%PrecipnTime(i)
                                ! Storing stamp time and values of time series 
                                this%PrecipTime(i,k) = MeteortimeSeries(k)%timeStamp
                                this%PrecipValue(i,k) = MeteortimeSeries(k)%value1
                            EndDo
                        EndIf
                    
                    ElseIf (trim(text) == 'Evaporation') Then
                        If (MeteorParameter(j)%functionMet == 1) Then ! Constant
                            Allocate (this%EvapValue(i,2))
                            Allocate (this%EvapTime(i,2))
                            Allocate (this%EvapnTime(i))
                            this%EvapnTime(i) = 2
                            this%EvapTime(i,1) = IniTime
                            this%EvapTime(i,2) = FinalTime
		                    this%EvapValue(i,:) = MeteorParameter(j)%constantValue
                        ElseIf (MeteorParameter(j)%functionMet == 2) Then ! Time-series
                            nTime = MeteorParameter(j)%timeSizeListLength
                            If (nTime<2) Then
                                print *,'Number of time intervals is lesser then two. Please check ' //trim(text) // 'time series.'// 'Station '// char(i)//'.' ! 
                                pause
                                Stop 
                            EndIf
                            Allocate (this%EvapValue(i,nTime))
                            Allocate (this%EvapTime(i,nTime))
                            Allocate (this%EvapnTime(i))
                            this%EvapnTime(i) = nTime
                            call c_f_pointer(MeteorParameter(j)%timeSeriesList, MeteortimeSeries, [MeteorParameter(j)%timeSizeListLength])
                            Do k = 1, this%EvapnTime(i)
                                ! Storing stamp time and values of time series 
                                this%EvapTime(i,k) = MeteortimeSeries(k)%timeStamp
                                this%EvapValue(i,k) = MeteortimeSeries(k)%value1
                            EndDo
                        EndIf
                    
                    EndIf
            
                
            
                EndDo
            EndDo
        EndIf
    
    End subroutine
        
        
End Module Meteorological
    
    