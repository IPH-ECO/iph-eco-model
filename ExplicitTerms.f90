Subroutine ExplicitTerms(HydroParam,MeshParam,dt)
    
    !$ use omp_lib
    Use MeshVars !, Only: nElem,nEdge,Edge,Left,Right,NormalVectorEdge,TangentVector,CirDistance,EdgeNodes,EdgeLength
    Use Hydrodynamic
    
    Implicit None
    Integer:: iElem,iEdge,iLayer,Face,Sig
    Integer:: l,r,iNode1,iNode2
    Real:: fc,DuuDnn,DuuDtt,Fu1,Fu2,DwwDxx,DwwDyy,Fw1,Fw2
    Real:: NearZero = 1e-10
    Real:: dt
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    
   
    ! Turbulence and Coriolis effect
    If (HydroParam%iConv == 0) Then
        !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iLayer,Face,l,r,fc,Fu1,Fu2,DuuDnn,DuuDtt,iNode1,iNode2,dt)
        Do Face = 1,MeshParam%nEdge
            l = MeshParam%Left(Face)
            r = MeshParam%Right(Face)
            
            Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face)
                HydroParam%Wu(iLayer,Face) = 0.
                If (r == 0) Then
                    HydroParam%Wu(iLayer,Face) = 0.
                Else
                    If (iLayer>=HydroParam%ElSmallm(l).and.iLayer>=HydroParam%ElSmallm(r)) Then
                        If (HydroParam%iCoriolis == 1) Then
                            fc = 2.*HydroParam%Omega*sin(HydroParam%Pi*HydroParam%Lat/180.)
                            HydroParam%Wu(iLayer,Face) = fc*HydroParam%Fv(iLayer,Face)
                        EndIf
                        Fu2 = -Sig(l,r,l)*HydroParam%Fub(iLayer,1,l)*MeshParam%NormalVectorEdge(1,Face)  - Sig(l,r,l)*HydroParam%Fub(iLayer,2,l)*MeshParam%NormalVectorEdge(2,Face)
                        Fu1 = -Sig(l,r,l)*HydroParam%Fub(iLayer,1,r)*MeshParam%NormalVectorEdge(1,Face)  - Sig(l,r,l)*HydroParam%Fub(iLayer,2,r)*MeshParam%NormalVectorEdge(2,Face)
                        DuuDnn = (Fu1 - 2.*HydroParam%Fu(iLayer,Face) + Fu2)/((MeshParam%CirDistance(Face)/2.)**2.)       ! Second derivative of U related to X in (I-1/2,J)
        
                        iNode1 = MeshParam%EdgeNodes(1,Face)
                        iNode2 = MeshParam%EdgeNodes(2,Face)
                        Fu1 = -Sig(l,r,l)*HydroParam%FuxyNode(iLayer,1,iNode1)*MeshParam%NormalVectorEdge(1,Face)  - Sig(l,r,l)*HydroParam%FuxyNode(iLayer,2,iNode1)*MeshParam%NormalVectorEdge(2,Face)
                        Fu2 = -Sig(l,r,l)*HydroParam%FuxyNode(iLayer,1,iNode2)*MeshParam%NormalVectorEdge(1,Face)  - Sig(l,r,l)*HydroParam%FuxyNode(iLayer,2,iNode2)*MeshParam%NormalVectorEdge(2,Face)
                        DuuDtt = (Fu1 - 2.*HydroParam%Fu(iLayer,Face) + Fu2)/(MeshParam%EdgeLength(Face)/2.)**2. 
                        HydroParam%Wu(iLayer,Face) = HydroParam%Wu(iLayer,Face) + 0.5*(HydroParam%HorViscosity(1,iLayer,l)+HydroParam%HorViscosity(1,iLayer,r))*( DuuDnn + DuuDtt )

                        HydroParam%Fu(iLayer,Face) = HydroParam%Fu(iLayer,Face) + dt*HydroParam%Wu(iLayer,Face)
                    EndIf
                EndIf
            EndDo
        EndDo
        !!$OMP end parallel do
        
        !If (HydroParam%iNonHydro==0.and.MeshParam%Kmax>1) Then ! Non-hydrostatic pressure
            Do iElem = 1,MeshParam%nElem 
                Do iLayer = HydroParam%ElSmallm(iElem)+1, HydroParam%ElCapitalM(iElem)
                    If (MeshParam%Neighbor(4,iElem)==0) Then
                        Fw1 = HydroParam%w(iLayer,iElem)
                    Else
                        Fw1 = HydroParam%w(iLayer,MeshParam%Neighbor(4,iElem))
                    EndIf
                    If (MeshParam%Neighbor(2,iElem)==0) Then
                        Fw2 = HydroParam%w(iLayer,iElem)
                    Else
                        Fw2 = HydroParam%w(iLayer,MeshParam%Neighbor(2,iElem))
                    EndIf
                    DwwDxx = (Fw1 - 2.*HydroParam%w(iLayer,iElem) + Fw2)/(MeshParam%DX**2)
                
                    If (MeshParam%Neighbor(1,iElem)==0) Then
                        Fw1 = HydroParam%w(iLayer,iElem)
                    Else
                        Fw1 = HydroParam%w(iLayer,MeshParam%Neighbor(1,iElem))
                    EndIf
                    If (MeshParam%Neighbor(3,iElem)==0) Then
                        Fw2 = HydroParam%w(iLayer,iElem)
                    Else
                        Fw2 = HydroParam%w(iLayer,MeshParam%Neighbor(3,iElem))
                    EndIf
                    DwwDyy = (Fw1 - 2.*HydroParam%w(iLayer,iElem) + Fw2)/(MeshParam%DY**2)
                
                    HydroParam%Fw(iLayer,iElem) = HydroParam%Fw(iLayer,iElem) + 0.5*(HydroParam%HorViscosity(1,iLayer,iElem)+HydroParam%HorViscosity(1,iLayer-1,iElem))*(DwwDxx+DwwDxx)
                EndDo
            EndDo
        !EndIf
    
    Else
        !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iLayer,iElem,Face,l,r,fc,Fu1,Fu2,DuuDnn,DuuDtt,iNode1,iNode2,dt)
        Do iElem = 1,MeshParam%nElem 
            Do iEdge = 1,4
                Face = MeshParam%Edge(iEdge,iElem)
                l = MeshParam%Left(Face)
                r = MeshParam%Right(Face)
                Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face)
                    HydroParam%Wu(iLayer,Face) = 0.
                    If (r == 0) Then
                        HydroParam%Wu(iLayer,Face) = 0.
                    Else
                        If (HydroParam%iCoriolis == 1) Then
                            fc = 2.*HydroParam%Omega*sin(HydroParam%Pi*HydroParam%Lat/180.)
                            if (iEdge==1) then ! North face
                                HydroParam%Wu(iLayer,Face) = -fc*Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%uxy(iLayer,1,MeshParam%Edge(iEdge,iElem))
                            elseif (iEdge== 2) then !West face
                                HydroParam%Wu(iLayer,Face) = fc*Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%uxy(iLayer,2,MeshParam%Edge(iEdge,iElem))
                            elseif (iEdge== 3) then !South face
                                HydroParam%Wu(iLayer,Face) = -fc*Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%uxy(iLayer,1,MeshParam%Edge(iEdge,iElem))
                            elseif (iEdge== 4) then !East face
                                HydroParam%Wu(iLayer,Face) = fc*Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%uxy(iLayer,2,MeshParam%Edge(iEdge,iElem))
                            endif
                        EndIf
                        
            
                        Fu2 = -Sig(l,r,l)*HydroParam%ub(iLayer,1,l)*MeshParam%NormalVectorEdge(1,Face)  - Sig(l,r,l)*HydroParam%ub(iLayer,2,l)*MeshParam%NormalVectorEdge(2,Face)
                        Fu1 = -Sig(l,r,l)*HydroParam%ub(iLayer,1,r)*MeshParam%NormalVectorEdge(1,Face)  - Sig(l,r,l)*HydroParam%ub(iLayer,2,r)*MeshParam%NormalVectorEdge(2,Face)
                        DuuDnn = (Fu1 - 2.*HydroParam%u(iLayer,Face) + Fu2)/((MeshParam%CirDistance(Face)/2.)**2.)       ! Second derivative of U related to X in (I-1/2,J)
        
                        iNode1 = MeshParam%EdgeNodes(1,Face)
                        iNode2 = MeshParam%EdgeNodes(2,Face)
                        Fu1 = -Sig(l,r,l)*HydroParam%uNode(iLayer,1,iNode1)*MeshParam%NormalVectorEdge(1,Face)  - Sig(l,r,l)*HydroParam%uNode(iLayer,2,iNode1)*MeshParam%NormalVectorEdge(2,Face)
                        Fu2 = -Sig(l,r,l)*HydroParam%uNode(iLayer,1,iNode2)*MeshParam%NormalVectorEdge(1,Face)  - Sig(l,r,l)*HydroParam%uNode(iLayer,2,iNode2)*MeshParam%NormalVectorEdge(2,Face)
                        DuuDtt = (Fu1 - 2.*HydroParam%u(iLayer,Face) + Fu2)/(MeshParam%EdgeLength(Face)/2.)**2. 
                        HydroParam%Wu(iLayer,Face) = HydroParam%Wu(iLayer,Face) + 0.5*(HydroParam%HorViscosity(1,iLayer,l)+HydroParam%HorViscosity(1,iLayer,r))*( DuuDnn + DuuDtt )
            
                        HydroParam%Fu(iLayer,Face) = HydroParam%Fu(iLayer,Face) + dt*HydroParam%Wu(iLayer,Face)
                    EndIf
                EndDo
            EndDo
        EndDo
        
        !If (HydroParam%iNonHydro==0.and.MeshParam%Kmax>1) Then ! Non-hydrostatic pressure
            Do iElem = 1,MeshParam%nElem 
                Do iLayer = HydroParam%ElSmallm(iElem)+1, HydroParam%ElCapitalM(iElem)
                    If (MeshParam%Neighbor(4,iElem)==0) Then
                        Fw1 = HydroParam%Fw(iLayer,iElem)
                    Else
                        Fw1 = HydroParam%Fw(iLayer,MeshParam%Neighbor(4,iElem))
                    EndIf
                    If (MeshParam%Neighbor(2,iElem)==0) Then
                        Fw2 = HydroParam%Fw(iLayer,iElem)
                    Else
                        Fw2 = HydroParam%Fw(iLayer,MeshParam%Neighbor(2,iElem))
                    EndIf
                    DwwDxx = (Fw1 - 2.*HydroParam%Fw(iLayer,iElem) + Fw2)/(MeshParam%DX**2)
                
                    If (MeshParam%Neighbor(1,iElem)==0) Then
                        Fw1 = HydroParam%Fw(iLayer,iElem)
                    Else
                        Fw1 = HydroParam%Fw(iLayer,MeshParam%Neighbor(1,iElem))
                    EndIf
                    If (MeshParam%Neighbor(3,iElem)==0) Then
                        Fw2 = HydroParam%Fw(iLayer,iElem)
                    Else
                        Fw2 = HydroParam%Fw(iLayer,MeshParam%Neighbor(3,iElem))
                    EndIf
                    DwwDyy = (Fw1 - 2.*HydroParam%Fw(iLayer,iElem) + Fw2)/(MeshParam%DY**2)
                
                    HydroParam%Fw(iLayer,iElem) = HydroParam%Fw(iLayer,iElem) + 0.5*(HydroParam%HorViscosity(1,iLayer,iElem)+HydroParam%HorViscosity(1,iLayer-1,iElem))*(DwwDxx+DwwDxx)
                EndDo
            EndDo
        !EndIf        
        !!$OMP end parallel do
    EndIf

    !write(*,*) HorViscosity(:,1,1),VerEddyVisc(2,3)
    Return
End Subroutine ExplicitTerms