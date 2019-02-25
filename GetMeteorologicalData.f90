!> This subroutine reads the hydrodynamic parameters. 
Subroutine GetMeteorologicalData(MeshParam,MeteoParam,time,AirtempRef)
    
    Use Meteorological
    Use MeshVars
    
    Implicit none
    
    Integer:: i,iElem,iEdge
    Real:: t_interp(1)
    Real:: p_interp(1,1)
    Real:: weit
    Real:: t_interp_double
    Real:: NearZero = 1e-10
    Real:: AirtempRef
    Real:: time
    type(MeshGridParam) :: MeshParam
    type(MeteorologicalParam) :: MeteoParam
    Integer:: leftt(MeteoParam%NStation),rightt(MeteoParam%NStation)
    Real :: staint(MeteoParam%NStation)

    t_interp = time
    t_interp_double = time
    
    
    !1. Atmosferic pressure data
    If (MeteoParam%NStation == 0) Then
        MeteoParam%AtmPressure = 1.01325e5 
    ElseIf (MeteoParam%Nstation==1) Then
        
        call r8vec_bracket ( MeteoParam%AtmPressurenTime(1), MeteoParam%AtmPressureTime(1,1:MeteoParam%AtmPressurenTime(1)), t_interp_double, leftt(1),rightt(1) )
        
        Do iElem=1,MeshParam%nElem
            !Atmosferic pressure for each cell
            MeteoParam%AtmPressure(iElem) = &
                ( ( MeteoParam%AtmPressureTime(1,rightt(1)) - time                )  * MeteoParam%AtmPressureValue(1,leftt(1))    &
                + (                 time - MeteoParam%AtmPressureTime(1,leftt(1)) ) * MeteoParam%AtmPressureValue(1,rightt(1)) ) &
                / ( MeteoParam%AtmPressureTime(1,rightt(1))     - MeteoParam%AtmPressureTime(1,leftt(1)) )
            !Call interp_linear( 1, AtmPressurenTime(1), AtmPressureTime(1,1:AtmPressurenTime(1)), AtmPressureValue(1,1:AtmPressurenTime(1)), 1, t_interp, p_interp )
            !AtmPressure(iElem) = p_interp(1,1)
        EndDo
    Else
        Do i=1,MeteoParam%Nstation
            call r8vec_bracket ( MeteoParam%AtmPressurenTime(i), MeteoParam%AtmPressureTime(i,1:MeteoParam%AtmPressurenTime(i)), t_interp_double, leftt(i),rightt(i) )
        EndDo
        Do iElem=1,MeshParam%nElem
            !IDW interpolation for each cell
            weit = 0.
            Do i=1,MeteoParam%Nstation
                staint(i) = 1./(dsqrt((MeteoParam%CoordMetStations(1,i)-MeshParam%xb(iElem))**2.+(MeteoParam%CoordMetStations(2,i)-MeshParam%yb(iElem))**2.)+NearZero)
                weit = weit + staint(i)
            EndDo
            MeteoParam%AtmPressure(iElem) = 0.
            Do i=1,MeteoParam%Nstation
                MeteoParam%AtmPressure(iElem) = &
                    ( ( MeteoParam%AtmPressureTime(i,rightt(i)) - time                )  * MeteoParam%AtmPressureValue(i,leftt(i))    &
                    + (                 time - MeteoParam%AtmPressureTime(i,leftt(i)) ) * MeteoParam%AtmPressureValue(i,rightt(i)) ) &
                    / ( MeteoParam%AtmPressureTime(i,rightt(i))     - MeteoParam%AtmPressureTime(i,leftt(i)) )
                !Call interp_linear( 1, AtmPressurenTime(i), AtmPressureTime(i,1:AtmPressurenTime(i)), AtmPressureValue(i,1:AtmPressurenTime(i)), 1, t_interp, p_interp )
                !AtmPressure(iElem) = AtmPressure(iElem) + p_interp(1,1)*staint(i)/weit
            EndDo
        EndDo
    EndIf
    

    !2. Reading Air Temperature 
    If (MeteoParam%NStation == 0) Then
        MeteoParam%AirTemp = AirtempRef
    ElseIf (MeteoParam%Nstation==1) Then
        call r8vec_bracket ( MeteoParam%AirTempnTime(1), MeteoParam%AirTempTime(1,1:MeteoParam%AirTempnTime(1)), t_interp_double, leftt(1),rightt(1) )
        Do iElem=1,MeshParam%nElem
            !Value for each cell
            MeteoParam%AirTemp(iElem) = &
                ( ( MeteoParam%AirTempTime(1,rightt(1)) - time                )  * MeteoParam%AirTempValue(1,leftt(1))    &
                + (                 time - MeteoParam%AirTempTime(1,leftt(1)) ) * MeteoParam%AirTempValue(1,rightt(1)) ) &
                / ( MeteoParam%AirTempTime(1,rightt(1))     - MeteoParam%AirTempTime(1,leftt(1)) )
            !Call interp_linear( 1, AirTempnTime(1), AirTempTime(1,1:AirTempnTime(1)), AirTempValue(1,1:AirTempnTime(1)), 1, t_interp, p_interp )
            !AirTemp(iElem) = p_interp(1,1)
        EndDo
    Else
        Do i=1,MeteoParam%Nstation
            call r8vec_bracket ( MeteoParam%AirTempnTime(i), MeteoParam%AirTempTime(i,1:MeteoParam%AirTempnTime(i)), t_interp_double, leftt(i),rightt(i) )
        EndDo
        Do iElem=1,MeshParam%nElem
            !IDW interpolation for each cell
            weit = 0.
            Do i=1,MeteoParam%Nstation
                staint(i) = 1./(dsqrt((MeteoParam%CoordMetStations(1,i)-MeshParam%xb(iElem))**2.+(MeteoParam%CoordMetStations(2,i)-MeshParam%yb(iElem))**2.)+NearZero)
                weit = weit + staint(i)
            EndDo
            MeteoParam%AirTemp(iElem) = 0.
            Do i=1,MeteoParam%Nstation
                MeteoParam%AirTemp(iElem) = &
                    ( ( MeteoParam%AirTempTime(i,rightt(i)) - time                )  * MeteoParam%AirTempValue(i,leftt(i))    &
                    + (                 time - MeteoParam%AirTempTime(i,leftt(i)) ) * MeteoParam%AirTempValue(i,rightt(i)) ) &
                    / ( MeteoParam%AirTempTime(i,rightt(i))     - MeteoParam%AirTempTime(i,leftt(i)) )
                !Call interp_linear( 1, AirTempnTime(i), AirTempTime(i,1:AirTempnTime(i)), AirTempValue(i,1:AirTempnTime(i)), 1, t_interp, p_interp )
                !AirTemp(iElem) = AirTemp(iElem) + p_interp(1,1)*staint(i)/weit
            EndDo
        EndDo
    EndIf
    
    !3. Reading Relative Humidity
    If (MeteoParam%NStation == 0) Then
        MeteoParam%RelHum = 100.
    ElseIf (MeteoParam%Nstation==1) Then
        call r8vec_bracket ( MeteoParam%RelHumnTime(1), MeteoParam%RelHumTime(1,1:MeteoParam%RelHumnTime(1)), t_interp_double, leftt(1),rightt(1) )
        Do iElem=1,MeshParam%nElem
            !Value for each cell
            MeteoParam%RelHum(iElem) = &
                ( ( MeteoParam%RelHumTime(1,rightt(1)) - time                )  * MeteoParam%RelHumValue(1,leftt(1))    &
                + (                 time - MeteoParam%RelHumTime(1,leftt(1)) ) * MeteoParam%RelHumValue(1,rightt(1)) ) &
                / ( MeteoParam%RelHumTime(1,rightt(1))     - MeteoParam%RelHumTime(1,leftt(1)) )
            !Call interp_linear( 1, RelHumnTime(1), RelHumTime(1,1:RelHumnTime(1)), RelHumValue(1,1:RelHumnTime(1)), 1, t_interp, p_interp )
            !RelHum(iElem) = p_interp(1,1)
        EndDo
    Else
        Do i=1,MeteoParam%Nstation
            call r8vec_bracket ( MeteoParam%RelHumnTime(i), MeteoParam%RelHumTime(i,1:MeteoParam%RelHumnTime(i)), t_interp_double, leftt(i),rightt(i) )
        EndDo
        Do iElem=1,MeshParam%nElem
            !IDW interpolation for each cell
            weit = 0.
            Do i=1,MeteoParam%Nstation
                staint(i) = 1./(dsqrt((MeteoParam%CoordMetStations(1,i)-MeshParam%xb(iElem))**2.+(MeteoParam%CoordMetStations(2,i)-MeshParam%yb(iElem))**2.)+NearZero)
                weit = weit + staint(i)
            EndDo
            MeteoParam%RelHum(iElem) = 0.
            Do i=1,MeteoParam%Nstation
                MeteoParam%RelHum(iElem) = &
                    ( ( MeteoParam%RelHumTime(i,rightt(i)) - time                )  * MeteoParam%RelHumValue(i,leftt(i))    &
                    + (                 time - MeteoParam%RelHumTime(i,leftt(i)) ) * MeteoParam%RelHumValue(i,rightt(i)) ) &
                    / ( MeteoParam%RelHumTime(i,rightt(i))     - MeteoParam%RelHumTime(i,leftt(i)) )
                !Call interp_linear( 1, RelHumnTime(i), RelHumTime(i,1:RelHumnTime(i)), RelHumValue(i,1:RelHumnTime(i)), 1, t_interp, p_interp )
                !RelHum(iElem) = RelHum(iElem) + p_interp(1,1)*staint(i)/weit
            EndDo
        EndDo
    EndIf
    

    !4. Reading solar radiation data
    If (MeteoParam%NStation == 0) Then
        MeteoParam%SolarRad = 0.
    ElseIf (MeteoParam%Nstation==1) Then
        call r8vec_bracket ( MeteoParam%SolarRadnTime(1), MeteoParam%SolarRadTime(1,1:MeteoParam%SolarRadnTime(1)), t_interp_double, leftt(1),rightt(1) )
        Do iElem=1,MeshParam%nElem
            !Atmosferic pressure for each cell
            MeteoParam%SolarRad(iElem) = &
                ( ( MeteoParam%SolarRadTime(1,rightt(1)) - time                )  * MeteoParam%SolarRadValue(1,leftt(1))    &
                + (                 time - MeteoParam%SolarRadTime(1,leftt(1)) ) * MeteoParam%SolarRadValue(1,rightt(1)) ) &
                / ( MeteoParam%SolarRadTime(1,rightt(1))     - MeteoParam%SolarRadTime(1,leftt(1)) )
            !Call interp_linear( 1, SolarRadnTime(1), SolarRadTime(1,1:SolarRadnTime(1)), SolarRadValue(1,1:SolarRadnTime(1)), 1, t_interp, p_interp )
            !SolarRad(iElem) = p_interp(1,1)
        EndDo
    Else
        Do i=1,MeteoParam%Nstation
            call r8vec_bracket ( MeteoParam%SolarRadnTime(i), MeteoParam%SolarRadTime(i,1:MeteoParam%SolarRadnTime(i)), t_interp_double, leftt(i),rightt(i) )
        EndDo
        Do iElem=1,MeshParam%nElem
            !IDW interpolation for each cell
            weit = 0.
            Do i=1,MeteoParam%Nstation
                staint(i) = 1./(dsqrt((MeteoParam%CoordMetStations(1,i)-MeshParam%xb(iElem))**2.+(MeteoParam%CoordMetStations(2,i)-MeshParam%yb(iElem))**2.)+NearZero)
                weit = weit + staint(i)
            EndDo
            MeteoParam%SolarRad(iElem) = 0.
            Do i=1,MeteoParam%Nstation
                MeteoParam%SolarRad(iElem) = &
                    ( ( MeteoParam%SolarRadTime(i,rightt(i)) - time                )  * MeteoParam%SolarRadValue(i,leftt(i))    &
                    + (                 time - MeteoParam%SolarRadTime(i,leftt(i)) ) * MeteoParam%SolarRadValue(i,rightt(i)) ) &
                    / ( MeteoParam%SolarRadTime(i,rightt(i))     - MeteoParam%SolarRadTime(i,leftt(i)) )
                !Call interp_linear( 1, SolarRadnTime(i), SolarRadTime(i,1:SolarRadnTime(i)), SolarRadValue(i,1:SolarRadnTime(i)), 1, t_interp, p_interp )
                !SolarRad(iElem) = SolarRad(iElem) + p_interp(1,1)*staint(i)/weit
            EndDo
        EndDo
    EndIf
    
    !5. Reading wind data
    If (MeteoParam%NStation == 0) Then
        MeteoParam%WindX = 0.
        MeteoParam%WindY = 0.
    ElseIf (MeteoParam%Nstation==1) Then
        
        call r8vec_bracket ( MeteoParam%WindnTime(1), MeteoParam%WindTime(1,1:MeteoParam%WindnTime(1)), t_interp_double, leftt(1),rightt(1) )
        
        Do iEdge=1,MeshParam%nEdge
            !Value for each cell
            
            MeteoParam%WindX(iEdge) = &
                ( ( MeteoParam%WindTime(1,rightt(1)) - time                )  * MeteoParam%WindValuex(1,leftt(1))    &
                + (                 time - MeteoParam%WindTime(1,leftt(1)) ) * MeteoParam%WindValuex(1,rightt(1)) ) &
                / ( MeteoParam%WindTime(1,rightt(1))     - MeteoParam%WindTime(1,leftt(1)) )
            
            !Call interp_linear( 1, WindnTime(1), WindTime(1,1:WindnTime(1)), WindValuex(1,1:WindnTime(1)), 1, t_interp, p_interp )
            !WindX(iEdge) = p_interp(1,1)
            !Call interp_linear( 1, WindnTime(1), WindTime(1,1:WindnTime(1)), WindValuey(1,1:WindnTime(1)), 1, t_interp, p_interp )
            !WindY(iEdge) = p_interp(1,1)
           MeteoParam%WindY(iEdge) = &
                ( ( MeteoParam%WindTime(1,rightt(1)) - time                )  * MeteoParam%WindValuey(1,leftt(1))    &
                + (                 time - MeteoParam%WindTime(1,leftt(1)) ) * MeteoParam%WindValuey(1,rightt(1)) ) &
                / ( MeteoParam%WindTime(1,rightt(1))     - MeteoParam%WindTime(1,leftt(1)) )
            
        EndDo
    Else
        Do i=1,MeteoParam%Nstation
            call r8vec_bracket ( MeteoParam%WindnTime(i), MeteoParam%WindTime(i,1:MeteoParam%WindnTime(i)), t_interp_double, leftt(i),rightt(i) )
        EndDo
        Do iEdge=1,MeshParam%nEdge
            !IDW interpolation for each cell
            weit = 0.
            Do i=1,MeteoParam%Nstation
                staint(i) = 1./(dsqrt((MeteoParam%CoordMetStations(1,i)-MeshParam%EdgeBary(1,iEdge))**2.+(MeteoParam%CoordMetStations(2,i)-MeshParam%EdgeBary(2,iEdge))**2.)+NearZero)
                weit = weit + staint(i)
            EndDo
            MeteoParam%WindX(iEdge) = 0.
            MeteoParam%WindY(iEdge) = 0.
            Do i=1,MeteoParam%Nstation
               MeteoParam%WindX(iEdge) = &
                    ( ( MeteoParam%WindTime(i,rightt(i)) - time                )  * MeteoParam%WindValuex(i,leftt(i))    &
                    + (                 time - MeteoParam%WindTime(i,leftt(i)) ) * MeteoParam%WindValuex(i,rightt(i)) ) &
                    / ( MeteoParam%WindTime(i,rightt(i))     - MeteoParam%WindTime(i,leftt(i)) )
                !Call interp_linear( 1, WindnTime(i), WindTime(i,1:WindnTime(i)), WindValuex(i,1:WindnTime(i)), 1, t_interp, p_interp )
                !WindX(iEdge) = WindX(iEdge) + p_interp(1,1)*staint(i)/weit
                MeteoParam%WindY(iEdge) = &
                    ( ( MeteoParam%WindTime(i,rightt(i)) - time                )  * MeteoParam%WindValuey(i,leftt(i))    &
                    + (                 time - MeteoParam%WindTime(i,leftt(i)) ) * MeteoParam%WindValuey(i,rightt(i)) ) &
                    / ( MeteoParam%WindTime(i,rightt(i))     - MeteoParam%WindTime(i,leftt(i)) )
                !Call interp_linear( 1, WindnTime(i), WindTime(i,1:WindnTime(i)), WindValuex(i,1:WindnTime(i)), 1, t_interp, p_interp )
                !WindY(iEdge) = WindY(iEdge) + p_interp(1,1)*staint(i)/weit
            EndDo
        EndDo
    EndIf
   
    !6. Reading precipitation data
    If (MeteoParam%NStation == 0) Then
        MeteoParam%Precip = 0.
    ElseIf (MeteoParam%Nstation==1) Then
        call r8vec_bracket ( MeteoParam%PrecipnTime(1), MeteoParam%PrecipTime(1,1:MeteoParam%PrecipnTime(1)), t_interp_double, leftt(1),rightt(1) )
        Do iElem=1,MeshParam%nElem
            !Precipitation for each cell
            MeteoParam%Precip(iElem) = &
                ( ( MeteoParam%PrecipTime(1,rightt(1)) - time                )  * MeteoParam%PrecipValue(1,leftt(1))    &
                + (                 time - MeteoParam%PrecipTime(1,leftt(1)) ) * MeteoParam%PrecipValue(1,rightt(1)) ) &
                / ( MeteoParam%PrecipTime(1,rightt(1))     - MeteoParam%PrecipTime(1,leftt(1)) )
            !Call interp_linear( 1, PrecipnTime(1), PrecipTime(1,1:PrecipnTime(1)), PrecipValue(1,1:PrecipnTime(1)), 1, t_interp, p_interp )
            !Precip(iElem) = p_interp(1,1)
        EndDo
    Else
        Do i=1,MeteoParam%Nstation
            call r8vec_bracket ( MeteoParam%PrecipnTime(i), MeteoParam%PrecipTime(i,1:MeteoParam%PrecipnTime(i)), t_interp_double, leftt(i),rightt(i) )
        EndDo
        Do iElem=1,MeshParam%nElem
            !IDW interpolation for each cell
            weit = 0.
            Do i=1,MeteoParam%Nstation
                staint(i) = 1./(dsqrt((MeteoParam%CoordMetStations(1,i)-MeshParam%xb(iElem))**2.+(MeteoParam%CoordMetStations(2,i)-MeshParam%yb(iElem))**2.)+NearZero)
                weit = weit + staint(i)
            EndDo
            MeteoParam%Precip(iElem) = 0.
            Do i=1,MeteoParam%Nstation
                MeteoParam%Precip(iElem) = &
                    ( ( MeteoParam%PrecipTime(i,rightt(i)) - time                )  * MeteoParam%PrecipValue(i,leftt(i))    &
                    + (                 time - MeteoParam%PrecipTime(i,leftt(i)) ) * MeteoParam%PrecipValue(i,rightt(i)) ) &
                    / ( MeteoParam%PrecipTime(i,rightt(i))     - MeteoParam%PrecipTime(i,leftt(i)) )

                !Call interp_linear( 1, PrecipnTime(i), PrecipTime(i,1:PrecipnTime(i)), PrecipValue(i,1:PrecipnTime(i)), 1, t_interp, p_interp )
                !Precip(iElem) = Precip(iElem) + p_interp(1,1)*staint(i)/weit
            EndDo
        EndDo
    EndIf
    
    !7. Reading evaporation data
    If (MeteoParam%NStation == 0) Then
        MeteoParam%Evap = 0.
    ElseIf (MeteoParam%Nstation==1) Then
        call r8vec_bracket ( MeteoParam%EvapnTime(1), MeteoParam%EvapTime(1,1:MeteoParam%EvapnTime(1)), t_interp_double, leftt(1),rightt(1) )
        Do iElem=1,MeshParam%nElem
            !Value for each cell
            MeteoParam%Evap(iElem) = &
                ( ( MeteoParam%EvapTime(1,rightt(1)) - time                )  * MeteoParam%EvapValue(1,leftt(1))    &
                + (                 time - MeteoParam%EvapTime(1,leftt(1)) ) * MeteoParam%EvapValue(1,rightt(1)) ) &
                / ( MeteoParam%EvapTime(1,rightt(1))     - MeteoParam%EvapTime(1,leftt(1)) )
            !Call interp_linear( 1, EvapnTime(1), EvapTime(1,1:EvapnTime(1)), EvapValue(1,1:EvapnTime(1)), 1, t_interp, p_interp )
            !Evap(iElem) = p_interp(1,1)
        EndDo
    Else
        Do i=1,MeteoParam%Nstation
            call r8vec_bracket ( MeteoParam%EvapnTime(i), MeteoParam%EvapTime(i,1:MeteoParam%EvapnTime(i)), t_interp_double, leftt(i),rightt(i) )
        EndDo
        Do iElem=1,MeshParam%nElem
            !IDW interpolation for each cell
            weit = 0.
            Do i=1,MeteoParam%Nstation
                staint(i) = 1./(dsqrt((MeteoParam%CoordMetStations(1,i)-MeshParam%xb(iElem))**2.+(MeteoParam%CoordMetStations(2,i)-MeshParam%yb(iElem))**2.)+NearZero)
                weit = weit + staint(i)
            EndDo
            MeteoParam%Evap(iElem) = 0.
            Do i=1,MeteoParam%Nstation
            MeteoParam%Evap(iElem) = &
                ( ( MeteoParam%EvapTime(i,rightt(i)) - time                )  * MeteoParam%EvapValue(i,leftt(i))    &
                + (                 time - MeteoParam%EvapTime(i,leftt(i)) ) * MeteoParam%EvapValue(i,rightt(i)) ) &
                / ( MeteoParam%EvapTime(i,rightt(i))     - MeteoParam%EvapTime(i,leftt(i)) )
                !Call interp_linear( 1, EvapnTime(i), EvapTime(i,1:EvapnTime(i)), EvapValue(i,1:EvapnTime(i)), 1, t_interp, p_interp )
                !Evap(iElem) = Evap(iElem) + p_interp(1,1)*staint(i)/weit
            EndDo
        EndDo
    EndIf

    Return
    
End Subroutine GetMeteorologicalData
    