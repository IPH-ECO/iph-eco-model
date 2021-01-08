Subroutine MatOp(a,b,dt,HydroParam,MeshParam)
    
    ! Compute the Semi-Implicit Method Coefficient Matrix
    !
    ! [1] Casulli, V. A high-resolution wetting and drying algorithm for free-surfacehydrodynamics. 
    !     International Journal for Numerical Methods in Fluids, v. 60, n. 4 (2009), p. 391-408
    ! [2] Casulli, Vincenzo. A conservative semi‐implicit method for coupled surface–subsurface flows in regional scale. 
    !     International Journal for Numerical Methods in Fluids, v. 79, n. 4 (2015), p. 199-214.
    
    ! Input:
    ! a -> Free-Surface Elevation Vector
    ! Output:
    ! b -> Solution
    
    ! List of Modifications: 
    !   -> 10.03.2014: Routine Implementation (Rafael Cavalcanti)
    !   -> 07.05.2019: Update - Subsurface Flow Term (Cayo Lopes)
    ! Programmer: Rafael Cavalcanti
    
    !$ use omp_lib
    Use MeshVars
    Use Hydrodynamic
     
    Implicit None
    Integer:: iElem, iEdge, Pij, Face
    Real:: Sum1, Coef
    Real:: NearZero = 1e-10
    Real:: dt
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Real, intent(in) :: a(MeshParam%nElem)
    Real, intent(out) :: b(MeshParam%nElem)
    Real:: aGhost
    
    !Coef = HydroParam%g*(HydroParam%Theta*dt)**2
    !call omp_set_num_threads(1)        
    ! 1. Compute T Matrix (Casulli, 2009)
    
    !!$OMP parallel do default(none) shared(a,b,MeshParam,HydroParam,NearZero,Coef) private(iElem,iEdge,Face,Pij,Sum)
    Do iElem = 1, MeshParam%nElem       
        
        b(iElem) = HydroParam%P(iElem)*a(iElem) !Initializing b vector
        Sum1 = 0.d0
        Do iEdge = 1,4
            Face = MeshParam%Edge(iEdge,iElem)
            Pij = MeshParam%Neighbor(iEdge,iElem)
            If (HydroParam%IndexWaterLevelEdge(Face)>0) Then
                !Face with pressure head boundary condition:
                !!Casulli,2015:
                Sum1 = Sum1 + ( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( - a(iElem) )*(HydroParam%Theta*dt)*(HydroParam%Theta*HydroParam%g*dt*HydroParam%DZiADZ(Face) + HydroParam%DZK(Face))
                !Sum1 = Sum1 + ( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( - a(iElem) )*(HydroParam%Theta*dt)*(HydroParam%Theta*HydroParam%g*dt*HydroParam%DZiADZ(Face) + HydroParam%Theta*HydroParam%DZK(Face))
            Else
                If (Pij == 0) Then
                    Sum1 = Sum1
                    continue
                Else
                    !Casulli,2015:
                    Sum1 = Sum1 + ( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( ( a(Pij) - a(iElem) ) )*(HydroParam%Theta*dt)*(HydroParam%Theta*HydroParam%g*dt*HydroParam%DZiADZ(Face) + HydroParam%DZK(Face))
                    !Sum1 = Sum1 + ( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( ( a(Pij) - a(iElem) ) )*(HydroParam%Theta*dt)*(HydroParam%Theta*HydroParam%g*dt*HydroParam%DZiADZ(Face) + HydroParam%Theta*HydroParam%DZK(Face))
                EndIf
            EndIf            
            !If (Pij == 0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero) Then
            !    If (HydroParam%IndexWaterLevel(iElem)>0) Then
            !        Sum1 = Sum1 + Coef*( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( ( - a(iElem) ) )*HydroParam%DZiADZ(Face) !( EdgeLength(Face)/CirDistance(Face) )*
            !    EndIf
            !Else
            !    Sum1 = Sum1 + Coef*( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( ( a(Pij) - a(iElem) ) )*HydroParam%DZiADZ(Face)
            !EndIf
        EndDo
        !b(iElem) = (P+T).eta
        b(iElem) = b(iElem) - Sum1 
        
    EndDo 
    !!$OMP end parallel do
   
    Return    
End Subroutine MatOp