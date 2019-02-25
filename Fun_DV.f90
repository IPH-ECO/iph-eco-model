Function DV(eta,h) 

    ! Flag for Wet or Dry Cell
    ! Casulli, V. A high-resolution wetting and drying algorithm for free-surfacehydrodynamics. 
    ! INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN FLUIDS, v. 60, n. 4 (2009), p. 391-408
    
    ! Input:
    ! eta -> Free-Surface Elevation
    ! h   -> Bottom Elevation
    ! Output:
    ! DV  -> Flag for Wet and Dry Cell
    
    ! List of Modifications: 
    !   -> 10.03.2014: Routine Implementation (Rafael Cavalcanti)
    ! Programmer: Rafael Cavalcanti
    
    Implicit None
    Double Precision:: DV, eta, h
    
    If (eta - h > 0) Then      ! Wet Cell
        DV = 1
    Else                     ! Dry Cell
        DV = 0
    EndIf
    
End Function DV