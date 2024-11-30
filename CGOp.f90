Subroutine CGOp(a,b,dt,HydroParam,MeshParam)
    
    ! Matrix Free Conjugate Gradient Method
    ! Input:
    ! a -> Input Matrix
    ! Output:
    ! b -> Solution
    
    ! List of Modifications: 
    !   -> 10.03.2014: Routine Implementation (Rafael Cavalcanti)
    !   -> 10.03.2014: Fortran Sintax         (Rafael Cavalcanti)
    ! Programmer: Michael Dumbser
    
    Use MeshVars !, Only:  nElem, Edge, Neighbor, EdgeLength, Cirdistance
    Use Hydrodynamic
    !Use SimulationModel, Only: dt,NearZero
    Implicit None
    Real:: alpha, tol, lambda, alphanew 
    Integer:: k, N
    Integer:: iElem, iEdge, Pij, Face
    Real:: Soma, Coef
    Real:: NearZero = 1e-10
    Real:: dt
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Real, intent(in) :: a(MeshParam%nElem)
    Real, intent(out) :: b(MeshParam%nElem)
    Real:: r(MeshParam%nElem), Ab(MeshParam%nElem), pCG(MeshParam%nElem), v(MeshParam%nElem)
    
    N = MeshParam%nElem
    b = a                           ! Initial guess
    Call MatOp(b,Ab,dt,HydroParam,MeshParam)
    ! b - Ab = %Deta - (P+T).eta(k,0)
    r = b - Ab
    pCG = r                         ! Steepest Descent Direction
    alpha = Dot_Product(r,r)        ! Square of the norm of r
    tol = 1e-14                     ! Tolerance
    
    Do k = 1,1000*N
       If ( sqrt(alpha) < tol ) Then
           ! System has been solved
           !Print*, 'The system has been solved with: ',k, 'CGOp iterations'
           Return
       EndIf

       Call MatOp(pCG,v,dt,HydroParam,MeshParam)
       
       lambda = alpha/Dot_Product(pCG,v)
       b = b + lambda*pCG
       r = r - lambda*v
       alphanew = Dot_Product(r,r)
       pCG = r + alphanew/alpha*pCG
       alpha = alphanew
       
    EndDo
        
    Return
End Subroutine CGOp