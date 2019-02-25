Subroutine MatOp_OPENMP(a,b,dt,g,Theta,PCRI,H,P,IndexWaterLevel,DZiADZ,nElem,Edge,Neighbor,EdgeLength,CirDistance,nEdge)
    
    ! Compute the Semi-Implicit Method Coefficient Matrix
    ! Casulli, V. A high-resolution wetting and drying algorithm for free-surfacehydrodynamics. 
    ! INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN FLUIDS, v. 60, n. 4 (2009), p. 391-408
    
    ! Input:
    ! a -> Free-Surface Elevation Vector
    ! Output:
    ! b -> Solution
    
    ! List of Modifications: 
    !   -> 10.03.2014: Routine Implementation (Rafael Cavalcanti)
    ! Programmer: Rafael Cavalcanti
    
    !$ use omp_lib
    Use MeshVars
    Use Hydrodynamic
    
    Implicit None
    Integer:: iElem, iEdge, Pij, Face
    Double Precision:: Sum, Coef
    Double Precision:: NearZero = 1e-10
    Double Precision:: dt
    !type(MeshGridParam) :: MeshParam
    !type(HydrodynamicParam) :: HydroParam
    Double Precision, intent(in) :: a(nElem)
    Double Precision, intent(out) :: b(nElem)
    Integer:: nElem,nEdge
    Double Precision:: g,Theta,PCRI
    Double Precision:: H(nEdge)
    Double Precision:: P(nElem)
    Double Precision:: IndexWaterLevel(nElem)
    Double Precision:: DZiADZ(nEdge)
    Integer:: Edge(4,nElem)
    Integer:: Neighbor(4,nElem)
    Double Precision:: EdgeLength(nEdge)
    Double Precision:: CirDistance(nEdge)
    
    
    Coef = g*(Theta*dt)**2
    !call omp_set_num_threads(2)        
    ! 1. Compute T Matrix (Casulli, 2009)
    !$OMP parallel do default(none) shared(a,b,P,Edge,Neighbor,EdgeLength,CirDistance,IndexWaterLevel,DZiADZ,nElem,H,Pcri,NearZero,Coef) private(iElem,iEdge,Face,Pij, Sum)
    Do iElem = 1, nElem
        b(iElem) = P(iElem)*a(iElem) ! Initializing b
        Sum = 0.
        Do iEdge = 1, 4
            Face = Edge(iEdge,iElem)
            Pij = Neighbor(iEdge,iElem)
            If (Pij == 0.or.H(Face) <= PCRI+NearZero) Then
                If (IndexWaterLevel(iElem)>0) Then
                    Sum = Sum + Coef*( EdgeLength(Face)/CirDistance(Face) )*( (  - a(iElem) ) )*DZiADZ(Face) !( EdgeLength(Face)/CirDistance(Face) )*
                EndIf
            Else
                Sum = Sum + Coef*( EdgeLength(Face)/CirDistance(Face) )*( ( a(Pij) - a(iElem) ) )*DZiADZ(Face)
            EndIf
        EndDo
        b(iElem) = b(iElem) - Sum
        
    EndDo
    !$OMP end parallel do
    
    Return    
End Subroutine MatOp_OPENMP