Subroutine CGOpNH(a,b,dt,HydroParam,MeshParam)
    
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
    Double Precision:: alpha, tol, lambda, alphanew 
    Integer:: k, N
    Integer:: iElem, iEdge,iLayer, Pij, Face
    Double Precision:: Soma, Coef
    Double Precision:: NearZero = 1e-10
    Double Precision:: dt
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    
    Double Precision, intent(in) :: a(MeshParam%kMax+1,MeshParam%nElem)
    Double Precision, intent(out) :: b(MeshParam%kMax+1,MeshParam%nElem)
    Double Precision:: r(MeshParam%kMax+1,MeshParam%nElem), Ab(MeshParam%kMax+1,MeshParam%nElem), pCG(MeshParam%kMax+1,MeshParam%nElem), v(MeshParam%kMax+1,MeshParam%nElem)
    
    N = MeshParam%nElem
    b = a                           ! Initial guess
    Call MatOpNH(b,Ab,dt,HydroParam,MeshParam)
    Do iElem = 1, MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)+1        
            r(iLayer,iElem) = b(iLayer,iElem) - Ab(iLayer,iElem)
        EndDo
    EndDo
    pCG = r                         ! Steepest Descent Direction
    alpha = 0.
    Do iElem = 1, MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)+1        
            alpha = alpha + r(iLayer,iElem)*r(iLayer,iElem)
        EndDo
    EndDo
    tol = 1e-20                      ! Tolerance
    
    Do k = 1,1000*N
       If ( sqrt(alpha) < tol ) Then
           ! System has been solved
           !Print*, 'The system has been solved with: ',k, 'CGOp iterations'
           Return
       EndIf
       Call MatOpNH(pCG,v,dt,HydroParam,MeshParam)
       Soma = 0.
       Do iElem = 1, MeshParam%nElem
           Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)+1        
               Soma = Soma + pCG(iLayer,iElem)*v(iLayer,iElem)
           EndDo
       EndDo
       lambda = alpha/Soma
       Do iElem = 1, MeshParam%nElem
           Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)+1        
               b(iLayer,iElem) = b(iLayer,iElem) + lambda*pCG(iLayer,iElem)
           EndDo
       EndDo
       Do iElem = 1, MeshParam%nElem
           Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)+1        
               r(iLayer,iElem) = r(iLayer,iElem) - lambda*v(iLayer,iElem)
           EndDo
       EndDo
       alphanew = 0.
       Do iElem = 1, MeshParam%nElem
           Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)+1        
               alphanew = alphanew + r(iLayer,iElem)*r(iLayer,iElem)
           EndDo
       EndDo
       Do iElem = 1, MeshParam%nElem
           Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)+1        
               pCG(iLayer,iElem) = r(iLayer,iElem) + alphanew/alpha*pCG(iLayer,iElem)
           EndDo
       EndDo
       alpha = alphanew
    EndDo
        
    Return
End Subroutine CGOpNH