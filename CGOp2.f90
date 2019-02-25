Subroutine CGOp(a,b)
    
    ! Matrix Free Conjugate Gradient Method
    ! Input:
    ! a -> Input Matrix
    ! Output:
    ! b -> Solution
    
    ! List of Modifications: 
    !   -> 10.03.2014: Routine Implementation (Rafael Cavalcanti)
    !   -> 10.03.2014: Fortran Sintax         (Rafael Cavalcanti)
    ! Programmer: Michael Dumbser
    
    Use MeshVars, Only: nElem
    Implicit None
    Double Precision:: a(nElem), b(nElem), r(nElem), Ab(nElem), pCG(nElem), v(nElem)
    Double Precision:: alpha, tol, lambda, alphanew
    Integer:: k, N
    
    N = nElem
    b = a                           ! Initial guess
    Call MatOp(b,Ab)
    r = b - Ab
    pCG = r                         ! Steepest Descent Direction
    alpha = Dot_Product(r,r)        ! Square of the norm of r
    tol = 1e-14                     ! Tolerance
    
    Do k = 1,N
       If ( sqrt(alpha) < tol ) Then
           ! System has been solved
           !Print*, 'The system has been solved with: ',k, 'CGOp iterations'
           Return
       EndIf
       Call MatOp(pCG,v)
       lambda = alpha/Dot_Product(pCG,v)
       b = b + lambda*pCG
       r = r - lambda*v
       alphanew = Dot_Product(r,r)
       pCG = r + alphanew/alpha*pCG
       alpha = alphanew
    EndDo
        
    Return
End Subroutine CGOp