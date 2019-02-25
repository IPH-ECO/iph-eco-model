Function AirEmiss(AirTemp) 
    ! Compute the air emissivity based on Swinbank's Equation
    ! Based on:
    ! [1] Swinbank, W.C. Longwave radiation from clear skies.
    !   Quartely Journal of the Royal Meteorological Society, Vol. 89, pp. 339-448. 1963.
    ! [2] Hodges, Ben. Heat Budget and Thermodynamics at a Free-Surface.
    !   CWR, University of Western Australia, 1999.
    ! Input:
    ! Air Temp -> Air Temperature (2m above the surface)
    ! Output:
    ! AirEmiss -> Emissivity of Air
    ! List of Modifications: 
    !   -> 14.08.2017: Routine Implementation (J. Rafael Cavalcanti)
    ! Programmer: J. Rafael Cavalcanti
    
    Implicit None
    Real:: AirTemp, AirEmiss
    Real, Parameter:: C_e = 0.937*1.E-5         ! Swinbank Coefficient (deg^-2)
    
    AirEmiss = C_e*( 273.15 + AirTemp )**2.
    
End Function AirEmiss    
! ------------------------------------------------------------------------------------    
Function Det(X1,X2,X3,Y1,Y2,Y3) 
    ! Compute the area of the triangle defined by three points with coordinates (X1,Y1), (X2,Y2) and (X3,Y3) using determinant formula 
    ! Based on:
    ! [1] Sloan, Scott W. A point-in-polygon program. 
    !   Advances in Engineering Software, Vol. 7 (1), pp. 45-47. 1985
    ! Input:
    ! X1, Y1 -> Coordinates of point 1
    ! X2, Y2 -> Coordinates of point 2
    ! X3, Y3 -> Coordinates of point 3
    ! Output:
    ! Det -> Area of the triangle defined by the three points
    ! List of Modifications: 
    !   -> 17.09.2015: Routine Implementation (J. Rafael Cavalcanti)
    !   -> 18.09.2015: Compute the Area       (J. Rafael Cavalcanti)
    ! Notes:
    ! Det is positive if points 1, 2, and 3 define triangle in anticlockwise order
    ! Det is negative if points 1, 2, and 3 define triangle in clockwise order
    ! Det is zero if at least two of the points are coincident of if all three points are collinear
    ! Programmer: J. Rafael Cavalcanti
    
    Implicit None
    Real:: X1, X2, X3, Y1, Y2, Y3, Det
    
    Det = 0.5*((X1 - X3)*(Y2 - Y3) - (X2 - X3)*(Y1 - Y3))
    
End Function Det    
! ------------------------------------------------------------------------------------    
Function DV(eta,h) 
    ! Flag for Wet or Dry Cell
    ! Casulli, V. A high-resolution wetting and drying algorithm for free-surfacehydrodynamics. 
    !   International Journal for Numerical Methods in Fluids, v. 60 (4), p. 391-408, 2009.
    ! Input:
    ! eta -> Free-Surface Elevation
    ! h   -> Bottom Elevation
    ! Output:
    ! DV  -> Flag for Wet and Dry Cell
    ! List of Modifications: 
    !   -> 10.03.2014: Routine Implementation (J. Rafael Cavalcanti)
    ! Programmer: J. Rafael Cavalcanti
    
    Implicit None
    Real:: DV, eta, h
    
    If (eta - h > 0) Then      ! Wet Cell
        DV = 1
    Else                     ! Dry Cell
        DV = 0
    End If
    
End Function DV
! ------------------------------------------------------------------------------------        
Function eair(Tair,RelHum) 
    ! Compute the Atmospheric Vapour Pressure
    ! Based on:
    ! [1] Chapra, Steve. Surface Water Quality Modeling.
    !   McGraw-Hill, 1997. ISBN 0-07-011364-5.
    ! Input:
    ! Tair   -> Air Temperature
    ! RelHum -> Air Relative Humidity 
    ! Output:
    ! eair -> Atmospheric Vapour Pressure (mmHg)
    ! List of Modifications: 
    !   -> 19.12.2016: Routine Implementation (J. Rafael Cavalcanti)
    ! Programmer: J. Rafael Cavalcanti

    Implicit None
    Real:: Aux, Tair, RelHum, eair
    ! 1. Compute the Saturation Vapour Pressure
    Aux = 4.596*EXP( (17.27*Tair)/(237.3 + Tair) )
    ! 2. Compute the Atmospheric Vapour Pressure
    eair = RelHum*Aux/100.
        
End Function eair
! ------------------------------------------------------------------------------------
Function esat(Tair) 
    ! Compute the Saturation Atmospheric Vapour Pressure
    ! Based on:
    ! [1] Chapra, Steve. Surface Water Quality Modeling.
    !   McGraw-Hill, 1997. ISBN 0-07-011364-5.
    ! [2] Hodges, Ben. Heat Budget and Thermodynamics at a Free-Surface.
    !   CWR, University of Western Australia, 1999.
    ! Input:
    ! Tair   -> Air Temperature (degC)
    ! Output:
    ! esat -> Atmospheric Vapour Pressure (mmHg)
    ! List of Modifications: 
    !   -> 19.12.2016: Routine Implementation           (J. Rafael Cavalcanti)
    !   -> 10.08.2017: Changes in Pressure Calculation  (J. Rafael Cavalcanti)
    ! Programmer: J. Rafael Cavalcanti
    ! Note:
    !   -> 760 mmHg = 101325     Pa
    !   ->   1 mmHg = 133.322365 Pa
    !   ->   1 Pa   = 0.007501   mmHg

    Implicit None
    Real:: esat, Tair, Aux
    
    ! 1. Compute the Saturation Vapour Pressure
    Aux = 4.596*EXP( (17.27*Tair)/(237.3 + Tair) )
    esat = Aux                                                         ! mmHg - [1]
    
End Function esat
! ------------------------------------------------------------------------------------
Function FWall(l,Bottom,Surface) 
    ! Compute the Wall Proximity Function for Mellor and Yamada Turbulence Scheme
    ! Based on:
    ! [1] Warner, John C.; Sherwood, Christopher R.; Arango, Hernan G.; Signell, Richard P. 2005. Performance of four turbulence closure models implemented using a generic length scale method.
    !   Ocean Modelling, 8, 81-113.
    ! [2] Mellor, George L.; Yamada, Tetsuji. 1982. Development of a Turbulence Closure Model for Geophysical Fluid Problems.
    !   Reviews of Geophysics and Space Physics, 20 (4), 851-875.
    ! Input:
    ! Bottom  -> Distance from Bottom
    ! Surface -> Distance from Free-Surface
    ! Output:
    ! FWall   -> Proximity Function
    ! List of Modifications: 
    !   -> 24.01.2017: Routine Implementation (J. Rafael Cavalcanti)
    ! Programmer: J. Rafael Cavalcanti
    
    Implicit None
    Real:: Bottom, Surface, l, FWall
    Real, Parameter:: E1 = 1.33 ! [1]
    Real, Parameter:: vonKarman = 0.4
    
    FWall = 1 + E1*( (l/vonKarman)*( (Bottom+Surface)/(Bottom*Surface) ) )**2.
    
End Function FWall    
! ------------------------------------------------------------------------------------
Function RhoAir(Pair,Tair,Sph) 
    ! Compute the Air Density based on equation of state for air: RhoAir = P/(R*Tv)
    ! Where R if the universal gas constant: 287.04 J/kg/K and Tv if the Virtual Temperature: Tv = T(K)*(1 + 0.6078*SpH)
    ! Based on:
    ! [1] SEML.m, a MatLab code from Sally MacIntyre. Available at Marice Science Institute, UCSB.
    ! Input:
    ! Pair   -> Air Pressure
    ! Tair   -> Air Temperature
    ! Sph    -> Specific Humidity
    ! Output:
    ! RhoAir -> Air Density (kg/m3)
    ! List of Modifications: 
    !   -> 03.08.2017: Routine Implementation (J. Rafael Cavalcanti)
    !   -> 13.09.2017: Bug Corrections        (J. Rafael Cavalcanti)
    ! Notes:
    !   -> Must Convert Temperature to units of K
    ! Programmer: J. Rafael Cavalcanti
    
    Implicit None
    Real:: Pair, Tair, RelHum, RhoAir
    Real:: P_mb, RH, svp, SpHSat, SpH, Tv
    Real, Parameter:: A1 = 0.7859, A2 = 0.03477, A3 = 0.00412
    Real, Parameter:: B1 = 1.E-6, B2 = 4.5, B3 = 0.006
    Real, Parameter:: R_gas = 287.04
    
    ! 1. Compute the Virtual Temperature
    Tv = (Tair + 273.15)*( 1 + 0.6078*Sph )
    ! 2. Compute the Air Density
    RhoAir = Pair/(R_gas*Tv)
        
End Function RhoAir    
! ------------------------------------------------------------------------------------
Function Sh(Gh) 
    ! Compute the Stability Function for Turbulent Eddy Diffusivity
    ! Based on:
    ! [1] Gross, Edward S.; Koseff, Jeffrey R.; Monismith, Stephen G. Three-Dimensional Salinity Simulations of South San Francisco Bay.
    !   Journal of Hydraulic Engineering, v. 125 (11), pp. 1999-1209, 1999.
    ! [2] Galperin, B.; Kantha, L.H.; Hassid, S.; Rosati, A. A quasi-equilibrium turbulent energy model for geophysical flows.
    !   Journal of Atmospheric Sciences, v. 45, pp. 55-62, 1998.
    ! Input:
    ! Gh -> Parameter
    ! Output:
    ! Sh   -> Stability Function
    ! List of Modifications: 
    !   -> 23.01.2017: Routine Implementation (J. Rafael Cavalcanti)
    ! Programmer: J. Rafael Cavalcanti
    
    Implicit None
    Real:: Sh, Gh
    Real, Parameter:: A1=0.92, A2=0.74, B1=16.6, B2=10.1, C3=0.2
    
    Sh = ( A2*(1-6*A1/B1) )/( 1 - 3*A2*Gh*(6*A1+B2*(1-C3)) )
    
End Function Sh    
! ------------------------------------------------------------------------------------
Function Sig(i,Right,Left)
    ! Signal Function
    ! Casulli, V.; Walters, R. An Unstructured Grid, Three-Dimensional Model based on the Shallow Water Equations. 
    !   International Journal for Numerical Methods in Fluids, v. 60 (4), p. 391-408, 2009.
    ! Input:
    ! i     -> Element
    ! Right -> Element on the Right of the Edge
    ! Left  -> Element on the Left  of the Edge
    ! Output:
    ! Sig   -> Flux Signal (< 0 - OutFlow; > 0 - InFlow) 
    ! List of Modifications: 
    !   -> 10.03.2014: Routine Implementation (J. Rafael Cavalcanti)
    ! Programmer: J. Rafael Cavalcanti
    
    Implicit None
    Integer:: i, Right, Left, Sig
    
    Sig = ( Right - 2*i + Left )/( Right - Left )
    
End Function Sig    
! ------------------------------------------------------------------------------------
Function Sm(Gh) 
    ! Compute the Stability Function for Turbulent Eddy Viscosity
    ! Based on:
    ! [1] Gross, Edward S.; Koseff, Jeffrey R.; Monismith, Stephen G. Three-Dimensional Salinity Simulations of South San Francisco Bay.
    !   Journal of Hydraulic Engineering, v. 125 (11), pp. 1999-1209, 1999.
    ! [2] Galperin, B.; Kantha, L.H.; Hassid, S.; Rosati, A. A quasi-equilibrium turbulent energy model for geophysical flows.
    !   Journal of Atmospheric Sciences, v. 45, 55-62, 1988.
    ! Input:
    ! Gh -> Parameter
    ! Output:
    ! Sm   -> Stability Function
    ! List of Modifications: 
    !   -> 23.01.2017: Routine Implementation (J. Rafael Cavalcanti)
    ! Programmer: J. Rafael Cavalcanti
    
    Implicit None
    Real:: Sm, Sh,Gh
    Real, Parameter:: A1=0.92, A2=0.74, B1=16.6, B2=10.1, C2=0.08, C3=0.2
    
    Sh = ( A2*(1-6*A1/B1) )/( 1 - 3*A2*Gh*(6*A1+B2*(1-C3)) )
    Sm = ( B1**(-1/3.) + (18*A1*A1 + 9*A1*A2*(1*C2))*Sh*Gh )/( 1 - 9*A1*A2*Gh )
    
End Function Sm    
! ------------------------------------------------------------------------------------
Function SpecHum(Temp,Pair,RelHum) 
    ! Compute the Specific Humidity    
    ! Based on:
    ! [1] SEML.m, a MatLab code from Sally MacIntyre. Available at Marice Science Institute, UCSB.
    ! Input:
    ! Temp   -> Air or Water Temperature (deg C)
    ! Pair   -> Air Pressure (Pa)
    ! RelHum -> Air Relative Humidity (%)
    ! Output:
    ! RhoAir -> Air Density (kg/m3)
    ! List of Modifications: 
    !   -> 03.08.2017: Routine Implementation (J. Rafael Cavalcanti)
    !   -> 13.09.2017: Bug Corrections        (J. Rafael Cavalcanti)
    ! Notes:
    !   -> Must Convert Pressure to units of mb (= 100 Pa)
    ! Programmer: J. Rafael Cavalcanti
    
    Implicit None
    Real:: Pair, Temp, RelHum, SpecHum
    Real:: P_mb, RH, svp, SpHSat, SpH
    Real, Parameter:: A1 = 0.7859, A2 = 0.03477, A3 = 0.00412
    Real, Parameter:: B1 = 1.E-6, B2 = 4.5, B3 = 0.006
    
    P_mb = Pair/100.; RH = RelHum/100.
    ! 1. Compute the Saturation Vapor Pressure
    svp = ( 10**( ( A1 + A2*Temp )/( 1 + A3*Temp ) ) )*( 1 + B1*P_mb*( B2 + B3*Temp**2. ) )
    ! 2. Compute the Saturation Specific Humidity
    SpHSat = ( 0.62197*svp )/( P_mb - ( 1 - 0.62197 )*svp )
    ! 3. Compute the Specific Humidity
    SpecHum = RH*SpHSat/( 1 - SpHSat + RH*SpHSat )
    
End Function SpecHum    
! ------------------------------------------------------------------------------------    
Function V(eta,h) 
    ! Compute the Water Elevation
    ! Based on:
    ! [1] Casulli, V. A high-resolution wetting and drying algorithm for free-surfacehydrodynamics. 
    !   International Journal for Numerical Methods in Fluids, v. 60 (4), p. 391-408, 2009.
    ! Input:
    ! eta -> Free-Surface Elevation
    ! h   -> Bottom Elevation
    ! Output:
    ! V   -> Water Elevation (Bottom to Free-Surface)
    ! List of Modifications: 
    !   -> 10.03.2014: Routine Implementation (J. Rafael Cavalcanti)
    ! Programmer: J. Rafael Cavalcanti
        
    Implicit None
    Real:: V, eta, h
    
    V = Max( 0., eta - h )
    
End Function V