Function V(eta,h) 

    ! Compute the Water Elevation
    ! Casulli, V. A high-resolution wetting and drying algorithm for free-surfacehydrodynamics. 
    ! INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN FLUIDS, v. 60, n. 4 (2009), p. 391-408
    
    ! Input:
    ! eta -> Free-Surface Elevation
    ! h   -> Bottom Elevation
    ! Output:
    ! V   -> Water Elevation (Bottom to Free-Surface)
    
    ! List of Modifications: 
    !   -> 10.03.2014: Routine Implementation (Rafael Cavalcanti)
    ! Programmer: Rafael Cavalcanti
    
    
    Implicit None
    Double Precision:: V, eta, h
    
    V = Max( 0., eta - h )
    
End Function V