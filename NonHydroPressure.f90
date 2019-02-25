Subroutine NonHydroPressure(HydroParam,MeshParam,dt)
    
    !$ use omp_lib
    Use MeshVars !, Only: nEdge,Left,Right,CirDistance
    Use Hydrodynamic
    Use ELM
    
    Implicit None
    Integer:: iElem, iEdge, iNode, iLayer, Face, iLayer_bar
    Integer:: l, r, Sig, iNewton, j, ie
    Real:: NearZero = 1e-10
    Real:: dt, SumRHS,SumRHSn, SumLayer, V, Coef,Coefn, res, SumW, SumWb, Hsublayer,sum1,sum0
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    !real :: start, finish
   
    
    Hsublayer = 0. !(MeshParam%zR - MeshParam%LIMCAMAUX(1))/1000.
    !! 1.0 Vertical velocities
    !HydroParam%wt = HydroParam%w
    
    Call wvelocity(HydroParam,MeshParam,dt)
    !Compute fictitious vertical velocity in top layer (for transport equation)
    Do iElem = 1, MeshParam%nElem
        SumW = 0.
        SumWb =0.
        Do iEdge = 1,4
            Face = MeshParam%Edge(iEdge,iElem)
            SumW = SumW + (Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*(HydroParam%DZj(HydroParam%ElCapitalM(iElem),Face)-Hsublayer)*HydroParam%u(HydroParam%ElCapitalM(iElem),Face))
            !SumWb = SumWb + (Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*(HydroParam%DZj(HydroParam%ElSmallm(iElem),Face)-Hsublayer)*HydroParam%u(HydroParam%ElSmallm(iElem),Face))
        EndDo
        HydroParam%w(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%w(HydroParam%ElCapitalM(iElem),iElem) - SumW/MeshParam%Area(iElem)         ! w(k+1/2,iElem) [1]
    EndDo
    
    ! 1.1 Assemble the Right Hand Side (RHS)
    !Coef = HydroParam%g*HydroParam%Theta*dt
    !Coefn = HydroParam%g*(1.-HydroParam%Theta)*dt
    Coefn = (1.-HydroParam%Theta)/(HydroParam%Theta)
    Do iElem = 1, MeshParam%nElem
        If (HydroParam%ElSmallm(iElem)==HydroParam%ElCapitalM(iElem)) Then
            HydroParam%rhsnonHydro(HydroParam%ElSmallm(iElem),iElem) = 0.
        Else
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)+1
                SumRHS = 0.
                SumRHSn = 0.
                If (iLayer==HydroParam%ElCapitalM(iElem)+1) Then
                    Do iEdge = 1,4
                        Face = MeshParam%Edge(iEdge,iElem)
                        SumLayer = 0.
                        Do iLayer_bar = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                            SumLayer = SumLayer + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer_bar,Face)*HydroParam%ut(iLayer_bar,Face)
                        EndDo
                        SumRHS  = SumRHS  + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*Hsublayer*HydroParam%u(iLayer-1,Face)
                        !SumRHS  = 0.
                        SumRHSn = SumRHSn + SumLayer
                    EndDo
                    HydroParam%rhsnonHydro(iLayer,iElem) = (HydroParam%w(iLayer,iElem)) - (SumRHS/MeshParam%Area(iElem)) + (HydroParam%etan(iElem) - HydroParam%eta(iElem))/(HydroParam%Theta*dt) - (Coefn*SumRHSn/MeshParam%Area(iElem))
                    !HydroParam%q(iLayer,iElem) = 0.
                ElseIf (iLayer==HydroParam%ElCapitalM(iElem)) Then
                    !Do iEdge = 1,4
                    !    Face = MeshParam%Edge(iEdge,iElem)
                    !    SumLayer = 0.
                    !    Do iLayer_bar = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                    !        SumLayer = SumLayer + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer_bar,Face)*HydroParam%ut(iLayer_bar,Face)
                    !    EndDo
                    !    SumRHS  = SumRHS  + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer,Face)*HydroParam%u(iLayer,Face)
                    !    !SumRHS  = 0.
                    !    SumRHSn = SumRHSn + SumLayer
                    !EndDo
                    !HydroParam%rhsnonHydro(iLayer,iElem) = (HydroParam%w(iLayer,iElem)) - (SumRHS/MeshParam%Area(iElem)) + (HydroParam%etan(iElem) - HydroParam%eta(iElem))/(HydroParam%Theta*dt) - (Coefn*SumRHSn/MeshParam%Area(iElem))
                    Do iEdge = 1,4
                        Face = MeshParam%Edge(iEdge,iElem)
                        SumRHS = SumRHS + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer,Face)*HydroParam%u(iLayer,Face)
                    EndDo
                    HydroParam%rhsnonHydro(iLayer,iElem) = (HydroParam%w(iLayer,iElem)-HydroParam%w(iLayer+1,iElem)) - (SumRHS/MeshParam%Area(iElem)) 

                    !HydroParam%q(iLayer+1,iElem) = 0.
                Else
                    Do iEdge = 1,4
                        Face = MeshParam%Edge(iEdge,iElem)
                        SumRHS = SumRHS + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*HydroParam%Dzjt(iLayer,Face)*HydroParam%u(iLayer,Face)
                    EndDo
                    HydroParam%rhsnonHydro(iLayer,iElem) = (HydroParam%w(iLayer,iElem)-HydroParam%w(iLayer+1,iElem)) - (SumRHS/MeshParam%Area(iElem)) 
                EndIf
            EndDo
        EndIf
    EndDo  
    ! 1.2 Non-hydrostatic correction
    Call CGOpNH(HydroParam%rhsnonHydro,HydroParam%q,dt, Hsublayer,HydroParam,MeshParam)
    !Do iElem = 1, MeshParam%nElem
    !    HydroParam%q(HydroParam%ElCapitalM(iElem)+1,:) = 0.
    !EndDo
    
    !! 7.3 Compute nodal dynamic pressure for tang. vel.
    !Do iNode=1,MeshParam%nNode
    !    Do iLayer = 1, MeshParam%Kmax
    !        sum1=0 !sum of elemental elevations
    !        sum0=0 !sum of areas
    !        do j=1,MeshParam%nVertexElem(iNode)
    !            ie=MeshParam%VertexElem(j,iNode)
    !            sum1=sum1+MeshParam%Area(ie)*HydroParam%q(iLayer,ie)
    !            sum0=sum0+MeshParam%Area(ie)
    !        Enddo !j=1,nne(i)
    !        HydroParam%pq(iLayer,iNode)=sum1/sum0
    !    EndDo
    !EndDo !i=1,np   
    
    ! 1.3 Correcting the provisional values of the horizontal velocities
    Do iEdge = 1,MeshParam%nEdge
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
        Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
            If (r==0) Then
                HydroParam%u(iLayer,iEdge) = HydroParam%u(iLayer,iEdge)
                !HydroParam%utang(iLayer,iEdge) = HydroParam%utang(iLayer,iEdge) - HydroParam%Theta*dt*(HydroParam%pq(iLayer,MeshParam%EdgeNodes(2,iEdge))-HydroParam%pq(iLayer,MeshParam%EdgeNodes(1,iEdge)))/MeshParam%EdgeLength(iEdge)
            Else
                HydroParam%u(iLayer,iEdge) = HydroParam%u(iLayer,iEdge) - HydroParam%Theta*dt*(HydroParam%q(iLayer,r)-HydroParam%q(iLayer,l))/MeshParam%CirDistance(iEdge)
                !HydroParam%utang(iLayer,iEdge) = HydroParam%utang(iLayer,iEdge) - HydroParam%Theta*dt*(HydroParam%pq(iLayer,MeshParam%EdgeNodes(2,iEdge))-HydroParam%pq(iLayer,MeshParam%EdgeNodes(1,iEdge)))/MeshParam%EdgeLength(iEdge)
                !HydroParam%u(iLayer,iEdge) = HydroParam%u(iLayer,iEdge) - dt*(HydroParam%q(iLayer,r)-HydroParam%q(iLayer,l))/MeshParam%CirDistance(iEdge)
            EndIf
        EndDo
    EndDo
    ! 1.4 Correcting the provisional values of the Free-Surface Elevation and vertical velocities
    Do iElem = 1, MeshParam%nElem

        HydroParam%eta(iElem) = HydroParam%eta(iElem) + HydroParam%q(HydroParam%ElCapitalM(iElem)+1,iElem)/HydroParam%g
        
            
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)-1 
            SumW = 0.
            Do iEdge = 1,4
                Face = MeshParam%Edge(iEdge,iElem)
                SumW = SumW + (Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*HydroParam%DZjt(iLayer,Face)*HydroParam%u(iLayer,Face))
            EndDo
            HydroParam%w(iLayer+1,iElem) = HydroParam%w(iLayer,iElem) - SumW/MeshParam%Area(iElem)         ! w(k+1/2,iElem) [1]
        EndDo
        
        !Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)-1
        !    HydroParam%w(iLayer+1,iElem) = HydroParam%w(iLayer+1,iElem) - HydroParam%Theta*dt*(HydroParam%q(iLayer+1,iElem)-HydroParam%q(iLayer,iElem))/(0.5*( HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer+1,iElem) ))
        !EndDo
        
        !Compute fictitious vertical velocity in top layer (for transport equation)
        SumW = 0.
        Do iEdge = 1,4
            Face = MeshParam%Edge(iEdge,iElem)
            SumW = SumW + (Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*HydroParam%DZjt(iLayer,Face)*HydroParam%u(iLayer,Face))
        EndDo
        HydroParam%w(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%w(HydroParam%ElCapitalM(iElem),iElem) - SumW/MeshParam%Area(iElem)         ! w(k+1/2,iElem) [1]
    EndDo
    
    ! 1.5 Compute nodal elevations for tang. vel.
    Do iNode=1,MeshParam%nNode
        sum1=0 !sum of elemental elevations
        sum0=0 !sum of areas
        do j=1,MeshParam%nVertexElem(iNode)
            ie=MeshParam%VertexElem(j,iNode)
            sum1=sum1+MeshParam%Area(ie)*HydroParam%eta(ie)
            sum0=sum0+MeshParam%Area(ie)
        Enddo !j=1,nne(i)
        HydroParam%peta(iNode)=sum1/sum0
    EndDo !i=1,np   
    
    ! 1.5 Calculating the non-hydrostatic pressure component
    Do iElem = 1, MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)+1
            HydroParam%q(iLayer,iElem) = HydroParam%q(iLayer,iElem) - HydroParam%q(HydroParam%ElCapitalM(iElem)+1,iElem)
        EndDo
    EndDo 
    
    !! 7.3 Compute nodal dynamic pressure for tang. vel.
    !Do iNode=1,MeshParam%nNode
    !    Do iLayer = 1, MeshParam%Kmax
    !        sum1=0 !sum of elemental elevations
    !        sum0=0 !sum of areas
    !        do j=1,MeshParam%nVertexElem(iNode)
    !            ie=MeshParam%VertexElem(j,iNode)
    !            sum1=sum1+MeshParam%Area(ie)*HydroParam%q(iLayer,ie)
    !            sum0=sum0+MeshParam%Area(ie)
    !        Enddo !j=1,nne(i)
    !        HydroParam%pq(iLayer,iNode)=sum1/sum0
    !    EndDo
    !EndDo !i=1,np   
    !call cpu_time(start)
    HydroParam%SumVer = 0.
	Do iElem = 1, MeshParam%nElem
		!Balanço de volume baseado na equação da continuidade
		Do iLayer =  HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) !ElSmallm(iElem), ElSmallm(iElem) !ElCapitalM(iElem), ElCapitalM(iElem) !
			If (iLayer==HydroParam%ElCapitalM(iElem)) Then
				SumW = 0.
				SumWb = 0.
				Do iEdge = 1,4
					Face = MeshParam%Edge(iEdge,iElem)
					Do iLayer_bar =  HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
						SumW = SumW + ((1.-HydroParam%theta)*dt*Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*HydroParam%DZjt(iLayer_bar,Face)*HydroParam%ut(iLayer_bar,Face))
					EndDo
					SumWb = SumWb + (HydroParam%theta*dt*Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*HydroParam%DZjt(iLayer,Face)*HydroParam%u(iLayer,Face))
				EndDo
				HydroParam%SumVer = HydroParam%SumVer + SumWb + SumW + (HydroParam%eta(iElem)-HydroParam%etan(iElem))*MeshParam%Area(iElem)
			Else
				SumW = 0.
				Do iEdge = 1,4
					Face = MeshParam%Edge(iEdge,iElem)
					SumW = SumW + (Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*HydroParam%DZjt(iLayer,Face)*HydroParam%u(iLayer,Face))
				EndDo
				HydroParam%SumVer = HydroParam%SumVer + SumW + (HydroParam%w(iLayer+1,iElem)-HydroParam%w(iLayer,iElem))*MeshParam%Area(iElem)
			EndIf
		EndDo
		!!Balanço de volume baseado na equação da continuidade integrada
		!SumVer = -Sig(iElem,Right(Edge(1,iElem)),MeshParam%Left(Edge(1,iElem)))*(Theta*u(1,Edge(1,iElem))*DZjt(1,Edge(1,iElem))+(1.-Theta)*ut(1,Edge(1,iElem))*DZjt(1,Edge(1,iElem)))*DX*DT-Sig(iElem,Right(Edge(2,iElem)),MeshParam%Left(Edge(2,iElem)))*(Theta*u(1,Edge(2,iElem))*DZjt(1,Edge(2,iElem))+ (1.-Theta)*ut(1,Edge(2,iElem))*DZjt(1,Edge(2,iElem)))*DY*simParam%DT-Sig(iElem,Right(Edge(3,iElem)),MeshParam%Left(Edge(3,iElem)))*(Theta*u(1,Edge(3,iElem))*DZjt(1,Edge(3,iElem))+ (1.-Theta)*ut(1,Edge(3,iElem))*DZjt(1,Edge(3,iElem)))*DY*simParam%DT-Sig(iElem,Right(Edge(4,iElem)),MeshParam%Left(Edge(4,iElem)))*(Theta*u(1,Edge(4,iElem))*DZjt(1,Edge(4,iElem))+(1.-Theta)*ut(1,Edge(4,iElem))*DZjt(1,Edge(4,iElem)))*DY*simParam%DT - (eta(iElem)-etaplus(ielem)-etan(iElem))*DX*DY 
		!SumVerAcum = SumVerAcum + SumVer
		!Write(*,*) '1-',SumVer,eta(iElem),sDSal(1,iElem)!SumVer, SumVerAcum!, VT2(I,J+1,1), VT1(I,J+1,1)
		!Print*,HydroParam%SumVer,HydroParam%SumVerAcum
		!pause
		!EndIf
	EndDo
	HydroParam%SumVerAcum = HydroParam%SumVerAcum + HydroParam%SumVer  
    !call cpu_time(finish)
    !print '("Time = ",f20.15," seconds.")',finish-start    
    
    Do iEdge = 1,MeshParam%nEdge
        
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
        ! 9.1 Compute Index Smallm(j) and CapitalM(j)
        ! Smallm and CapitalM are related to mj and Mj in [1]
        If (r == 0) Then
            HydroParam%H(iEdge) = Max( HydroParam%PCRI, -HydroParam%hj(iEdge) + HydroParam%eta(l) )
        Else
            HydroParam%H(iEdge) = Max( HydroParam%PCRI,0.5*((-HydroParam%hj(iEdge) + HydroParam%eta(l)) + (-HydroParam%hj(iEdge) + HydroParam%eta(r))) ) !Max( HydroParam%PCRI,-HydroParam%hj(iEdge) + HydroParam%eta(l), -HydroParam%hj(iEdge) + HydroParam%eta(r) ) !Max( PCRI,-hj(iEdge) + 0.5*(eta(l) + eta(r)) ) !
        EndIf

        !9.2 Compute Elevation in the Edges
        If (r == 0) Then
            HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%eta(l) !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
        Else
            If ( HydroParam%H(iEdge) <= HydroParam%PCRI+NearZero ) Then
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%eta(l),HydroParam%eta(r))
            Else
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = 0.5*(HydroParam%eta(l)+HydroParam%eta(r)) !Max(HydroParam%eta(l),HydroParam%eta(r)) !0.5*(eta(l)+eta(r)) !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
            EndIf
        EndIf
        
        ! 9.3 Compute the Vertical Mesh Spacing
        HydroParam%DZj(HydroParam%CapitalM(iEdge),iEdge) = HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) - HydroParam%Z(HydroParam%CapitalM(iEdge),iEdge)
        HydroParam%DZj(HydroParam%CapitalM(iEdge),iEdge) = Max(HydroParam%Pcri,HydroParam%DZj(HydroParam%CapitalM(iEdge),iEdge))
    EndDo    
    Do iElem = 1, MeshParam%nElem
       
        ! 11.1 Define the range between Bottom and Top Layers. (Tricky part - Not explained in [1])
        ! In a real world application, the bathymetry might change between the faces of the prism.
        ! We need to find "Max(Smallm(:))" and "Min(CapitalM(:))" in each face. 
        ! Then we can calculate the Vertical Velocity only in those Layers
        ! *Obs: The velocity in faces with Smallm < Max(Smallm(:)) == 0. The same goes to faces with CapitalM > Min(CapitalM(:))
        ! 4.1.1 Find "mi" and "Mi", the Bottom and Top Layers range for the Element, not for the faces
        ! Set mi
        
        ! Upper Index
        Do iLayer = 1,MeshParam%KMax
            If (iLayer == MeshParam%KMax) Then
                If (HydroParam%eta(iElem)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
                    HydroParam%ElCapitalM(iElem) = iLayer
                    exit
                EndIf
            Else
                If (HydroParam%eta(iElem)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%eta(iElem)<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
                    HydroParam%ElCapitalM(iElem) = iLayer
                    exit
                EndIf
            EndIf
        EndDo        
        
        
        !If ( Smallm(Edge(1,iElem))==Smallm(Edge(2,iElem)).AND.Smallm(Edge(2,iElem))==Smallm(Edge(3,iElem)).AND.Smallm(Edge(3,iElem))==Smallm(Edge(4,iElem)) ) Then
        !    ElSmallm(iElem) = Smallm(Edge(1,iElem))         ! All "mj" are equal. We can choose whatever we want!
        !Else
        !    ! If we have different Bottom Layers, we merge all in one cell!
        !    !ElSmallm(iElem) = Max( Smallm(Edge(1,iElem)),Smallm(Edge(2,iElem)),Smallm(Edge(3,iElem)),Smallm(Edge(4,iElem)) ) - 1
        !    !ElSmallm(iElem) = Max( Smallm(Edge(1,iElem)),Smallm(Edge(2,iElem)),Smallm(Edge(3,iElem)),Smallm(Edge(4,iElem)) )
        !    ElSmallm(iElem) = Min( Smallm(Edge(1,iElem)),Smallm(Edge(2,iElem)),Smallm(Edge(3,iElem)),Smallm(Edge(4,iElem)) ) 
        !EndIf
        !! Set Mi
        !If ( CapitalM(Edge(1,iElem))==CapitalM(Edge(2,iElem)).AND.CapitalM(Edge(2,iElem))==CapitalM(Edge(3,iElem)).AND.CapitalM(Edge(3,iElem))==CapitalM(Edge(4,iElem)) ) Then
        !    ElCapitalM(iElem) = CapitalM(Edge(1,iElem))         ! All "Mj" are equal. We can choose whatever we want!
        !Else
        !    ! If we have different Top Layers, we merge all in one cell!
        !    !ElCapitalM(iElem) = Min( CapitalM(Edge(1,iElem)),CapitalM(Edge(2,iElem)),CapitalM(Edge(3,iElem)),CapitalM(Edge(4,iElem)) ) + 1
        !    !ElCapitalM(iElem) = Min( CapitalM(Edge(1,iElem)),CapitalM(Edge(2,iElem)),CapitalM(Edge(3,iElem)),CapitalM(Edge(4,iElem)) )
        !    ElCapitalM(iElem) = Min( CapitalM(Edge(1,iElem)),CapitalM(Edge(2,iElem)),CapitalM(Edge(3,iElem)),CapitalM(Edge(4,iElem)) )
        !EndIf
        ! 11.2 Update the Element Vertical Spacing
        HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%eta(iElem) !- hb(iElem)       ! Free-Surface
        Do iLayer = HydroParam%ElCapitalM(iElem)+2, MeshParam%KMax+1
            HydroParam%Ze(iLayer,iElem) = HydroParam%eta(iElem)
        EndDo
        
        ! 11.3 Compute the Vertical Mesh Spacing
        HydroParam%DZi(HydroParam%ElCapitalM(iElem),iElem) = HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(HydroParam%ElCapitalM(iElem),iElem)
        HydroParam%DZi(HydroParam%ElCapitalM(iElem),iElem) = Max(HydroParam%Pcri,HydroParam%DZi(HydroParam%ElCapitalM(iElem),iElem))
        HydroParam%Zb(HydroParam%ElCapitalM(iElem),iElem) = 0.5*( HydroParam%Ze(HydroParam%ElCapitalM(iElem),iElem) + HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) )



    EndDo
   
    
    Return
End Subroutine NonHydroPressure