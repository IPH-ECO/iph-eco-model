Subroutine MatOpNH(a,b,dt,Hsublayer,HydroParam,MeshParam)
    
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
    Integer:: iElem,iLayer, iEdge, Pij, Face
    Real:: Sum1, Coef, test
    Real:: NearZero = 1e-10
    Real:: dt, dzp,dzm,aGhost, Hsublayer
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Real, intent(in) :: a(MeshParam%kMax+1,MeshParam%nElem)
    Real, intent(out) :: b(MeshParam%kMax+1,MeshParam%nElem)
    
    Coef = (HydroParam%Theta*dt)
    !call omp_set_num_threads(1)        
    ! 1. Compute T Matrix (Casulli, 2009)
    b = 0.
    Do iElem = 1, MeshParam%nElem
        If (HydroParam%ElSmallm(iElem)==HydroParam%ElCapitalM(iElem)) Then
            b(HydroParam%ElSmallm(iElem),iElem) = 0.
        Else
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)+1      
                Sum1 = 0.
                If (iLayer==HydroParam%ElSmallm(iElem)) Then
                    Do iEdge = 1, 4
                        Face = MeshParam%Edge(iEdge,iElem)
                        Pij = MeshParam%Neighbor(iEdge,iElem)
                        If (Pij == 0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero) Then
                            !Sum1 = Sum1 + Coef*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer,Face)*(ai(Layer,iElem) - aGhost)/MeshParam%CirDistance(Face)                
                        Else
                            If (iLayer>=HydroParam%ElSmallm(MeshParam%Neighbor(iEdge,iElem))) Then
                                Sum1 = Sum1 + (Coef/MeshParam%Area(iElem))*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer,Face)*(a(iLayer,iElem) - a(iLayer,Pij))/MeshParam%CirDistance(Face)
                            EndIf
                        EndIf
                    EndDo 
                    dzp = 0.5d0*( HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer+1,iElem) )       ! Dz at Upper Interface
                    b(iLayer,iElem) = Sum1 + Coef*((a(iLayer,iElem)-a(iLayer+1,iElem))/dzp)
                    !dzp = 0.5*( HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer+1,iElem) )       ! Dz at Upper Interface
                    !dzm = 0.5*( HydroParam%DZit(iLayer,iElem) )       ! Dz at Lower Interface
                    !b(iLayer,iElem) = Sum1 + Coef*((a(iLayer,iElem)-a(iLayer+1,iElem))/dzp - (-a(iLayer,iElem))/dzm)
                 ElseIf (iLayer==HydroParam%ElCapitalM(iElem)+1) Then
                    
                    Do iEdge = 1, 4 
                        Face = MeshParam%Edge(iEdge,iElem)
                        Pij = MeshParam%Neighbor(iEdge,iElem)
                        If (Pij == 0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero) Then
                            !Sum1 = Sum1 + Coef*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer,Face)*(a(iLayer,iElem) - aGhost)/MeshParam%CirDistance(Face)                
                        Else
                            Sum1 = Sum1 + (Coef/MeshParam%Area(iElem))*MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face)*Hsublayer*(a(iLayer,iElem) - a(iLayer,Pij))     
                            !Sum1 = 0.
                        EndIf
                    
                    EndDo
                    dzm = 0.5*(HydroParam%DZit(iLayer-1,iElem))       ! Dz at Lower Interface
                    b(iLayer,iElem) = Sum1 - Coef*((a(iLayer-1,iElem)-a(iLayer,iElem))/dzm) + (a(iLayer,iElem)/(HydroParam%g*HydroParam%Theta*dt))
                    !b(iLayer,iElem) = Sum1 - Coef*(a(iLayer-1,iElem)/dzm) !+ (a(iLayer,iElem)/(HydroParam%g*HydroParam%Theta*dt))
                    
                ElseIf (iLayer==HydroParam%ElCapitalM(iElem)) Then
                    
                    !Do iEdge = 1, 4 
                    !    Face = MeshParam%Edge(iEdge,iElem)
                    !    Pij = MeshParam%Neighbor(iEdge,iElem)
                    !    If (Pij == 0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero) Then
                    !        !Sum1 = Sum1 + Coef*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer,Face)*(a(iLayer,iElem) - aGhost)/MeshParam%CirDistance(Face)                
                    !    Else
                    !        If (iLayer>=HydroParam%ElSmallm(MeshParam%Neighbor(iEdge,iElem))) Then
                    !            Sum1 = Sum1 + (Coef/MeshParam%Area(iElem))*MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face)*HydroParam%Dzjt(iLayer,Face)*(0.5d0*(a(iLayer,iElem)+a(iLayer+1,iElem)) - 0.5d0*(a(iLayer,Pij)+a(iLayer+1,Pij)))     
                    !        EndIf
                    !    EndIf
                    !
                    !EndDo
                    !dzm = 0.5*( HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer-1,iElem) )       ! Dz at Lower Interface
                    !b(iLayer,iElem) = Sum1 - Coef*((0.5d0*(a(iLayer-1,iElem)+a(iLayer,iElem))-0.5d0*(a(iLayer,iElem)+a(iLayer+1,iElem)))/dzm) + (0.5d0*(a(iLayer,iElem)+a(iLayer+1,iElem))/(HydroParam%g*HydroParam%Theta*dt))
                    !b(iLayer,iElem) = Sum1 - Coef*(a(iLayer-1,iElem)/dzm) !+ (a(iLayer,iElem)/(HydroParam%g*HydroParam%Theta*dt))
                    Do iEdge = 1, 4
                        Face = MeshParam%Edge(iEdge,iElem)
                        Pij = MeshParam%Neighbor(iEdge,iElem)
                        If (Pij == 0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero) Then
                            !Sum1 = Sum1 + Coef*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer,Face)*(a(iLayer,iElem) - aGhost)/MeshParam%CirDistance(Face)                
                        Else
                            !If (iLayer>=HydroParam%ElSmallm(MeshParam%Neighbor(iEdge,iElem))) Then
                                Sum1 = Sum1 + (Coef/MeshParam%Area(iElem))*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer,Face)*(a(iLayer,iElem) - a(iLayer,Pij))/MeshParam%CirDistance(Face) 
                            !EndIf
                        EndIf
                    EndDo
                    dzp = 0.5*( HydroParam%DZit(iLayer,iElem))       ! Dz at Upper Interface
                    dzm = 0.5*( HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer-1,iElem) )       ! Dz at Lower Interface
                    b(iLayer,iElem) = Sum1 + Coef*((a(iLayer,iElem)-a(iLayer+1,iElem))/dzp - (a(iLayer-1,iElem)-a(iLayer,iElem))/dzm)

                Else
                    Do iEdge = 1, 4
                        Face = MeshParam%Edge(iEdge,iElem)
                        Pij = MeshParam%Neighbor(iEdge,iElem)
                        If (Pij == 0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero) Then
                            !Sum1 = Sum1 + Coef*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer,Face)*(a(iLayer,iElem) - aGhost)/MeshParam%CirDistance(Face)                
                        Else
                            If (iLayer>=HydroParam%ElSmallm(MeshParam%Neighbor(iEdge,iElem))) Then
                                Sum1 = Sum1 + (Coef/MeshParam%Area(iElem))*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer,Face)*(a(iLayer,iElem) - a(iLayer,Pij))/MeshParam%CirDistance(Face) 
                            EndIf
                        EndIf
                    EndDo
                    dzp = 0.5*( HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer+1,iElem) )       ! Dz at Upper Interface
                    dzm = 0.5*( HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer-1,iElem) )       ! Dz at Lower Interface
                    b(iLayer,iElem) = Sum1 + Coef*((a(iLayer,iElem)-a(iLayer+1,iElem))/dzp - (a(iLayer-1,iElem)-a(iLayer,iElem))/dzm)
                EndIf
            EndDo
        EndIf
    EndDo
    
    Return
    
    End Subroutine MatOpNH
    
    
