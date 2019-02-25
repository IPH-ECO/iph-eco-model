Function Sig(i,Right,Left)
    
    ! Signal Function
    ! Casulli, V.; Walters, R. An Unstructured Grid, Three-Dimensional Model based on the Shallow Water Equations. 
    ! INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN FLUIDS, v. 32, n. 3 (2000), p. 331-348
    
    ! Input:
    ! i     -> Element
    ! Right -> Element on the Right of the Edge
    ! Left  -> Element on the Left  of the Edge
    ! Output:
    ! Sig   -> Flux Signal (< 0 - OutFlow; > 0 - InFlow) 
    
    ! List of Modifications: 
    !   -> 10.03.2014: Routine Implementation (Rafael Cavalcanti)
    ! Programmer: Rafael Cavalcanti
    
    Implicit None
    Integer:: i, Right, Left, Sig
    
    Sig = ( Right - 2*i + Left )/( Right - Left )
    
    Return
End Function Sig