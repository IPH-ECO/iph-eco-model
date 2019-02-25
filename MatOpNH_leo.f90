Subroutine MatOpNH(a,b,dt,HydroParam,MeshParam)
    
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
    Double Precision:: Sum1, Coef, test, aGhost
    Double Precision:: NearZero = 1e-10
    Double Precision:: dt, dzp,dzm
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Double Precision, intent(in) :: a(MeshParam%kMax+1,MeshParam%nElem)
    Double Precision, intent(out) :: b(MeshParam%kMax+1,MeshParam%nElem)
    !HydroParam%Theta*
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

                        Else
                            Sum1 = Sum1 + dt*HydroParam%Dzjt(iLayer,Face)*(a(iLayer,Pij) - a(iLayer,iElem))/(MeshParam%CirDistance(Face)**2.)
                        EndIf
                    EndDo
                    dzp = 0.5*( HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer+1,iElem) )       ! Dz at Upper Interface
                    b(iLayer,iElem) = Sum1 + dt*((a(iLayer+1,iElem)-a(iLayer,iElem))/dzp)
                ElseIf (iLayer==HydroParam%ElCapitalM(iElem)+1) Then
                    Do iEdge = 1, 4 
                        Face = MeshParam%Edge(iEdge,iElem)
                        Pij = MeshParam%Neighbor(iEdge,iElem)
                        If (Pij == 0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero) Then

                        Else
                            Sum1 = Sum1 + dt*HydroParam%PCRI*(a(iLayer,Pij) - a(iLayer,iElem))/(MeshParam%CirDistance(Face)**2.)
                        EndIf
                    
                    EndDo
                    dzm = HydroParam%DZit(iLayer-1,iElem)/2. !0.5*( HydroParam%PCRI + HydroParam%DZit(iLayer-1,iElem) )       ! Dz at Lower Interface
                    b(iLayer,iElem) = Sum1 - dt*((a(iLayer,iElem)-a(iLayer-1,iElem))/dzm) - a(iLayer,iElem)/(HydroParam%g*HydroParam%Theta*dt) !dt*a(iLayer,iElem)/ HydroParam%DZit(iLayer,iElem)
                Else
                    If (iLayer==HydroParam%ElCapitalM(iElem)) Then
                        Do iEdge = 1, 4
                            Face = MeshParam%Edge(iEdge,iElem)
                            Pij = MeshParam%Neighbor(iEdge,iElem)
                            If (Pij == 0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero) Then

                            Else
                                Sum1 = Sum1 + dt*(HydroParam%Dzjt(iLayer,Face)-HydroParam%PCRI)*(a(iLayer,Pij) - a(iLayer,iElem))/(MeshParam%CirDistance(Face)**2.)
                            EndIf
                        EndDo
                        dzp = HydroParam%DZit(iLayer,iElem)/2. !0.5*( HydroParam%DZit(iLayer,iElem) + HydroParam%PCRI )       ! Dz at Upper Interface
                        dzm = 0.5*( HydroParam%DZit(iLayer,iElem)-HydroParam%PCRI + HydroParam%DZit(iLayer-1,iElem) )       ! Dz at Lower Interface
                        b(iLayer,iElem) = Sum1 + dt*((a(iLayer+1,iElem)-a(iLayer,iElem))/dzp - (a(iLayer,iElem)-a(iLayer-1,iElem))/dzm)
                    Else
                        Do iEdge = 1, 4
                            Face = MeshParam%Edge(iEdge,iElem)
                            Pij = MeshParam%Neighbor(iEdge,iElem)
                            If (Pij == 0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero) Then

                            Else
                                Sum1 = Sum1 + dt*HydroParam%Dzjt(iLayer,Face)*(a(iLayer,Pij) - a(iLayer,iElem))/(MeshParam%CirDistance(Face)**2.)
                            EndIf
                        EndDo
                        dzp = 0.5*( HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer+1,iElem) )       ! Dz at Upper Interface
                        dzm = 0.5*( HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer-1,iElem) )       ! Dz at Lower Interface
                        b(iLayer,iElem) = Sum1 + dt*((a(iLayer+1,iElem)-a(iLayer,iElem))/dzp - (a(iLayer,iElem)-a(iLayer-1,iElem))/dzm)
                    EndIf
                    
                    
                EndIf
            EndDo
        EndIf
    EndDo
    
    Return    
End Subroutine MatOpNH