!> This subroutine reads the hydrodynamic parameters. 
Subroutine ReadHydroIniCond(HydroParam,hydroConfiguration,simParam,MeshParam)
    
    Use domain_types
    Use SimulationModel
    Use MeshVars
    Use Hydrodynamic
    Use Sorting
    Use iso_c_binding
    
    Implicit none
    
    Integer:: i,j,iElem,iEdge,iNode,iLayer,r,l,Face,ie,Sig
    !character, pointer :: hydroParametersName(:)
    type(HydrodynamicConfiguration) :: hydroConfiguration  
    type(HydrodynamicParameter), pointer :: hydroParameters(:)
    type(SimulationParam) :: simParam
    character(len=200):: text
    Real:: sum1,sum0,AuxVel,V
    Real:: NearZero = 1e-10
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Real, Dimension(4,MeshParam%nElem):: HNeighbor
    Integer, Dimension(MeshParam%nElem):: nNeighbor
    Integer:: npit
    
    call c_f_pointer(hydroConfiguration%parameters, hydroParameters, [hydroConfiguration%numberOfParameters])
    
    Do i = 1, hydroConfiguration%numberOfParameters
        text = trim(hydroParameters(i)%name)

        !Flags
        If (trim(text) == 'envFreeSurfaceElevation') Then
            HydroParam%Zini = hydroParameters(i)%value
        ElseIf (trim(text) == 'envWaterVelocityX') Then
            HydroParam%Uini  = hydroParameters(i)%value
        ElseIf (trim(text) == 'envWaterVelocityY') Then
            HydroParam%Vini = hydroParameters(i)%value
        ElseIf (trim(text) == 'envWaterVelocityW') Then
            HydroParam%Wini = hydroParameters(i)%value
        EndIf
        !print *, hydroParametersName, hydroParameters(i)%value
    EndDo
    HydroParam%uIniVec(:,1) = HydroParam%Uini
    HydroParam%uIniVec(:,2) = HydroParam%Vini
    HydroParam%uIniVec(:,3) = HydroParam%Wini
    HydroParam%uxyback = 0.d0
    HydroParam%uxy = 0.d0
    HydroParam%Fu = 0.d0
    HydroParam%Fv = 0.d0
    HydroParam%Fub = 0.d0
    HydroParam%FuxyNode = 0.d0
    HydroParam%u = 0.d0
    HydroParam%utang = 0.d0
    HydroParam%ub = 0.d0
    HydroParam%uNode = 0.d0
    HydroParam%HorViscosity = 0.d0
    HydroParam%HorDiffusivity = 0.d0
    If ( HydroParam%iVTurb == 0 ) Then 
        HydroParam%VerEddyVisc = HydroParam%VerEddyVisc_Cte
        HydroParam%VerEddyDiff = HydroParam%VerEddyDiff_Cte
        HydroParam%VerEddyViscCell = HydroParam%VerEddyVisc_Cte
        HydroParam%VerEddyDiffCell = HydroParam%VerEddyDiff_Cte
    Elseif ( HydroParam%iVTurb == 1 ) Then     
        HydroParam%VerEddyVisc = HydroParam%vmin
        HydroParam%VerEddyDiff = HydroParam%tdmin_pp
        HydroParam%VerEddyViscCell = HydroParam%vmin
        HydroParam%VerEddyDiffCell = HydroParam%tdmin_pp
    EndIf
    HydroParam%iAG = 0.d0
    HydroParam%iADZ = 1.
    HydroParam%DZj = HydroParam%Pcri
    HydroParam%DZjt = HydroParam%Pcri
    HydroParam%DZi = HydroParam%Pcri
    HydroParam%DZit = HydroParam%Pcri
    HydroParam%Gu = 0.d0
    
    !CAYO
    HydroParam%DZK = 0.d0
    MeshParam%Kj = 0.d0
    HydroParam%DZsi = HydroParam%Pcri
    HydroParam%DZsit = HydroParam%Pcri
    HydroParam%DZhi = HydroParam%Pcri
    HydroParam%DZhit = HydroParam%Pcri
    HydroParam%us = 0.d0
    HydroParam%um = 0.d0    
    
    HydroParam%rhsnonHydro = 0.d0
    HydroParam%q = 0.d0
    HydroParam%pq = 0.d0
     
    ! Set Initial Conditions of Free-Surface Elevation
    If (MeshParam%ieta0 == 1) Then
        Do iElem = 1, MeshParam%nElem
            HydroParam%eta(iElem) = MeshParam%eta0(iElem)
        EndDo
    Else
        HydroParam%eta = HydroParam%Zini
    EndIf
    
    !HydroParam%eta(1) = 0.01    
    
    
    If (simParam%it > 0) Then
        HydroParam%eta = simParam%etasave
    EndIf
    
    HydroParam%etan = HydroParam%eta
    
    !!Filling sinks
    !npit = 1 
    !Do While (npit /= 0)   
    !    
    !    Do iElem = 1, MeshParam%nElem
    !        nNeighbor(iElem) = 0
    !        Do iEdge = 1, 4
    !            Face = MeshParam%Edge(iEdge,iElem)
    !            l = MeshParam%Left(Face) 
    !            r = MeshParam%Right(Face)
    !            If (r == 0) Then 
    !                nNeighbor(iElem) = nNeighbor(iElem)
    !            Else
    !                nNeighbor(iElem) = nNeighbor(iElem) +1
    !                HNeighbor(nNeighbor(iElem),iElem) = HydroParam%hb(MeshParam%Neighbor(iEdge,iElem))
    !            EndIf
    !        EndDo
    !        
    !    EndDo
    !    npit = 0
    !    Do iElem = 1, MeshParam%nElem
    !        Call Sort(HNeighbor(:,iElem),nNeighbor(iElem))
    !        If (nNeighbor(iElem) == 0) Then
    !            HydroParam%hb(iElem) = HydroParam%hb(iElem)
    !        Else
    !            If (HydroParam%hb(iElem)<MinVal(HNeighbor(1:nNeighbor(iElem),iElem))) Then
    !                HydroParam%hb(iElem) = MinVal(HNeighbor(1:nNeighbor(iElem),iElem))
    !                npit = npit + 1
    !            Elseif (HydroParam%hb(iElem)==MinVal(HNeighbor(1:nNeighbor(iElem),iElem))) Then
    !                If (HydroParam%hb(iElem)<MinVal(HNeighbor(2:nNeighbor(iElem),iElem)).and.nNeighbor(iElem)>=2) Then
    !                    HydroParam%hb(iElem) = MinVal(HNeighbor(2:nNeighbor(iElem),iElem))
    !                    npit = npit + 1
    !                Else
    !                    HydroParam%hb(iElem) = HydroParam%hb(iElem)
    !                EndIf
    !            Else
    !                HydroParam%hb(iElem) = HydroParam%hb(iElem)
    !            EndIf
    !        EndIf
    !        
    !    EndDo
    !EndDo
    !
    
    Do iElem = 1, MeshParam%nElem
        
        Do iLayer = 1,MeshParam%KMax
            If (abs(HydroParam%hb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<HydroParam%PCRI) then
                HydroParam%hb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
            EndIf
       EndDo
        
        If (MeshParam%KMax==1) Then ! Two-dimensional
            Do iLayer = 1,MeshParam%KMax
                If (abs(HydroParam%hb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<HydroParam%PCRI) then
                    HydroParam%hb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
                EndIf
            EndDo
        ElseIf (MeshParam%KMax>1) Then ! Three-dimensional
            Do iLayer = 1,MeshParam%KMax-1
                
                If (HydroParam%hb(iElem)+0.001>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%hb(iElem)+0.001<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) Then
                    If (abs(HydroParam%hb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<abs(HydroParam%hb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer))) Then ! 
                        HydroParam%hb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
                    Else
                        HydroParam%hb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)
                    EndIf
                EndIf
            EndDo
        EndIf
        
        If (HydroParam%eta(iElem) - HydroParam%hb(iElem) <= HydroParam%Pcri.and.HydroParam%eta(iElem) - HydroParam%hb(iElem) >= 0.) Then
            HydroParam%hb(iElem) = HydroParam%eta(iElem) - HydroParam%Pcri
        EndIf
        
        If (HydroParam%eta(iElem) - HydroParam%hb(iElem) < HydroParam%Pcri) Then
            HydroParam%eta(iElem) = HydroParam%hb(iElem)+HydroParam%Pcri !*** Verificar sinal posteriormente, se já entrar como cota não precisa do sinal 
        EndIf
        
    EndDo
    
    
    !Compute nodal elevations for tang. vel.
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
    HydroParam%petan=HydroParam%peta !store for transport eqs.
    
    ! Set Bathymetry in th edges
    Do iElem = 1, MeshParam%nElem
        
        Do iEdge = 1, 4
            Face = MeshParam%Edge(iEdge,iElem)
            l = MeshParam%Left(Face) 
            r = MeshParam%Right(Face)
            If (r == 0) Then
                HydroParam%hj(Face) = HydroParam%hb(l) !*** Verificar sinal posteriormente, se já entrar como cota não precisa do sinal 
                HydroParam%sj(Face) = HydroParam%sb(l)
                Do iLayer = 1,MeshParam%KMax
                    If (abs(HydroParam%hj(Face)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<HydroParam%PCRI) then
                        HydroParam%hj(Face)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
                    EndIf
                    If (abs(HydroParam%sj(Face)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<HydroParam%PCRI) then
                        HydroParam%sj(Face)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
                    EndIf

                EndDo
            Else
                HydroParam%hj(Face) = Max(HydroParam%hb(l),HydroParam%hb(r))  !  (0.5d0)*( Hbat(I,J) + HBat(IViz,JViz) ) !Max(Hbat(I,J),HBat(IViz,JViz))!*** !(0.5d0)*( Hbat(I,J) + HBat(IViz,JViz) ) Verificar sinal posteriormente, se já entrar como cota não precisa do sinal 
                HydroParam%sj(Face) = Max(HydroParam%sb(l),HydroParam%sb(r))
                Do iLayer = 1,MeshParam%KMax
                    If (abs(HydroParam%hj(Face)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<HydroParam%PCRI) then
                        HydroParam%hj(Face)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
                    EndIf
                    If (abs(HydroParam%sj(Face)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<HydroParam%PCRI) then
                        HydroParam%sj(Face)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
                    EndIf
                EndDo

            EndIf
        EndDo
        
    EndDo
    
    
    !Compute Index Smallm(j) and CapitalM(j)
    ! Smallm and CapitalM are related to mj and Mj in [1]
    Do iEdge = 1,MeshParam%nEdge

        ! 1.1 Compute Water Depth - Superficial Flow Layer
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
        If (r == 0) Then
            HydroParam%H(iEdge) = Max( HydroParam%PCRI, -HydroParam%hj(iEdge) + HydroParam%eta(l) )
        Else
            HydroParam%H(iEdge) = Max( HydroParam%PCRI,-HydroParam%hj(iEdge) + HydroParam%eta(l), -HydroParam%hj(iEdge) + HydroParam%eta(r) ) !Max( PCRI,-hj(iEdge) + (0.5d0)*(eta(l) + eta(r)) ) !
        EndIf
 
        
        ! Lower Index - Superficial Flow Layer
        Do iLayer = 1,MeshParam%KMax
            If (iLayer == MeshParam%KMax) Then
                If (HydroParam%hj(iEdge)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
                    HydroParam%Smallm(iEdge) = iLayer
                    exit
                EndIf
            Else
                If (HydroParam%hj(iEdge)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%hj(iEdge)<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
                    HydroParam%Smallm(iEdge) = iLayer
                    exit
                EndIf
            EndIf
        EndDo

        ! Upper Index - Superficial Flow Layer
        Do iLayer = 1,MeshParam%KMax
            If (iLayer == MeshParam%KMax) Then
                If (HydroParam%H(iEdge) + HydroParam%hj(iEdge)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
                    HydroParam%CapitalM(iEdge) = iLayer
                    exit
                EndIf
            Else
                If (HydroParam%H(iEdge) + HydroParam%hj(iEdge)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%H(iEdge) + HydroParam%hj(iEdge)<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
                    HydroParam%CapitalM(iEdge) = iLayer
                    exit
                EndIf
            EndIf
        EndDo

        !Do iLayer = HydroParam%Smallm(iEdge) + 1, HydroParam%CapitalM(iEdge)
        !    HydroParam%Z(iLayer,iEdge) = MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1) !zL + (iLayer - 1)*dz                  ! Equidistant Core Grid
        !EndDo
        !HydroParam%Z(HydroParam%Smallm(iEdge),iEdge)     = HydroParam%hj(iEdge)                 ! Bottom
        !If (r == 0) Then
        !    If ( HydroParam%eta(l) - HydroParam%hj(iEdge) <= HydroParam%PCRI+NearZero) Then
        !        HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%hj(iEdge) + HydroParam%PCRI
        !    Else
        !        HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%eta(l)
        !    EndIf
        !     !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
        !Else
        !    If ( HydroParam%H(iEdge) <= HydroParam%PCRI+NearZero ) Then
        !        HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%eta(l),HydroParam%eta(r))
        !    Else
        !        HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = (0.5d0)*(HydroParam%eta(l)+HydroParam%eta(r)) !Max(HydroParam%eta(l),HydroParam%eta(r)) !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
        !    EndIf
        !EndIf

        !! 1.4 Compute the Vertical Mesh Spacing
        !Do iLayer = HydroParam%Smallms(iEdge), HydroParam%CapitalM(iEdge) !CAYO
        !    HydroParam%DZj(iLayer,iEdge) = HydroParam%Z(iLayer+1,iEdge) - HydroParam%Z(iLayer,iEdge)
        !    HydroParam%DZj(iLayer,iEdge) = Max(HydroParam%Pcri,HydroParam%DZj(iLayer,iEdge))
        !EndDo

        ! Compute Water Depth - Subsuperficial Flow Layer  !CAYO  
        If (r == 0) Then
            HydroParam%Hs(iEdge) = Max( HydroParam%PCRI, -HydroParam%sj(iEdge) + HydroParam%eta(l) ) !CAYO
        Else
            HydroParam%Hs(iEdge) = Max( HydroParam%PCRI,-HydroParam%sj(iEdge) + HydroParam%eta(l), -HydroParam%sj(iEdge) + HydroParam%eta(r) )!CAYO
        EndIf
 
        
        ! Lower Index - Subsuperficial Flow Layer  !CAYO
        Do iLayer = 1,MeshParam%KMax
            If (iLayer == MeshParam%KMax) Then
                If (HydroParam%sj(iEdge)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
                    HydroParam%Smallms(iEdge) = iLayer
                    exit
                EndIf
            Else
                If (HydroParam%sj(iEdge)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%sj(iEdge)<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
                    HydroParam%Smallms(iEdge) = iLayer
                    exit
                EndIf
            EndIf
        EndDo

        ! Upper Index - Subsuperficial Flow Layer  !CAYO
        If ( HydroParam%Smallm(iEdge) > 1) Then   ! Cayo 
            HydroParam%CapitalMs(iEdge) = HydroParam%Smallm(iEdge) - 1
            if (HydroParam%CapitalMs(iEdge)<HydroParam%Smallms(iEdge)) then
                HydroParam%CapitalMs(iEdge) = HydroParam%Smallms(iEdge)
            endif
        Else
            HydroParam%CapitalMs(iEdge) = HydroParam%Smallms(iEdge)
        EndIf
        
        !9.2 Compute Elevation in the Edges !CAYO
        Do iLayer = HydroParam%Smallms(iEdge) + 1, HydroParam%CapitalM(iEdge)
            HydroParam%Z(iLayer,iEdge) = MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1) !zL + (iLayer - 1)*dz                  ! Equidistant Core Grid
        EndDo
        HydroParam%Z(HydroParam%Smallms(iEdge),iEdge)     = HydroParam%sj(iEdge)                 ! Bottom
        If (r == 0) Then
            If ( HydroParam%eta(l) - HydroParam%sj(iEdge) <= HydroParam%PCRI+NearZero) Then
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%sj(iEdge) + HydroParam%PCRI
            Else
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%eta(l)
            EndIf
             !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
        Else
            If ( HydroParam%Hs(iEdge) <= HydroParam%PCRI+NearZero ) Then
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%eta(l),HydroParam%eta(r))
            Else
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = (0.5d0)*(HydroParam%eta(l)+HydroParam%eta(r)) !Max(HydroParam%eta(l),HydroParam%eta(r)) !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
            EndIf
        EndIf        
        
        
        ! 1.4 Compute the Vertical Mesh Spacing
        Do iLayer = HydroParam%Smallms(iEdge), HydroParam%CapitalM(iEdge) !CAYO
            HydroParam%DZj(iLayer,iEdge) = HydroParam%Z(iLayer+1,iEdge) - HydroParam%Z(iLayer,iEdge)
            HydroParam%DZj(iLayer,iEdge) = Max(HydroParam%Pcri,HydroParam%DZj(iLayer,iEdge))
        EndDo        
        
        !1.5 Compute Kj and DZhj/DZsj - Cayo (Loop adicionado)        
        Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge) ! 
            If (iLayer >= HydroParam%Smallms(iEdge) .or. HydroParam%Smallms(iEdge) == HydroParam%Smallm(iEdge) ) Then
                MeshParam%Kj(iLayer,iEdge) = 0.d0
            Else
                MeshParam%Kj(iLayer,iEdge) = 0.01
            EndIf
        EndDo  

    EndDo
    

    
    ! 4. Compute the New Vertical Velocity Field
    ! w(Smallm-1,:) = 0.d0 -> No flux through the Bottom
    ! w(CapitalM,:) = 0.d0 -> No flux through the Surface
    ! The velocity is defined in the barycenter of the Top Face of each Cell
    ! *Obs: The outward direction is from Bottom to Top.
    Do iElem = 1, MeshParam%nElem
        ! 4.1 Define the range between Bottom and Top Layers. (Tricky part - Not explained in [1])
        ! In a real world application, the bathymetry might change between the faces of the prism.
        ! We need to find "Max(Smallm(:))" and "Min(CapitalM(:))" in each face. 
        ! Then we can calculate the Vertical Velocity only in those Layers
        ! *Obs: The velocity in faces with Smallm < Max(Smallm(:)) == 0. The same goes to faces with CapitalM > Min(CapitalM(:))
        ! 4.1.1 Find "mi" and "Mi", the Bottom and Top Layers range for the Element, not for the faces
        ! Set mi
        !If ( Smallm(Edge(1,iElem))==Smallm(Edge(2,iElem)).AND.Smallm(Edge(2,iElem))==Smallm(Edge(3,iElem)).AND.Smallm(Edge(3,iElem))==Smallm(Edge(4,iElem)) ) Then
        !    ElSmallm(iElem) = Smallm(Edge(1,iElem))         ! All "mj" are equal. We can choose whatever we want!
        !Else
        !    ! If we have different Bottom Layers, we merge all in one cell!
        !    !ElSmallm(iElem) = Max( Smallm(Edge(1,iElem)),Smallm(Edge(2,iElem)),Smallm(Edge(3,iElem)),Smallm(Edge(4,iElem)) ) !- 1
        !    ElSmallm(iElem) = Min( Smallm(Edge(1,iElem)),Smallm(Edge(2,iElem)),Smallm(Edge(3,iElem)),Smallm(Edge(4,iElem)) )
        !EndIf
        !! Set Mi
        !If ( CapitalM(Edge(1,iElem))==CapitalM(Edge(2,iElem)).AND.CapitalM(Edge(2,iElem))==CapitalM(Edge(3,iElem)).AND.CapitalM(Edge(3,iElem))==CapitalM(Edge(4,iElem)) ) Then
        !    ElCapitalM(iElem) = CapitalM(Edge(1,iElem))         ! All "Mj" are equal. We can choose whatever we want!
        !Else
        !    ! If we have different Top Layers, we merge all in one cell!
        !    !ElCapitalM(iElem) = Min( CapitalM(Edge(1,iElem)),CapitalM(Edge(2,iElem)),CapitalM(Edge(3,iElem)),CapitalM(Edge(4,iElem)) ) !+ 1
        !    ElCapitalM(iElem) = Max( CapitalM(Edge(1,iElem)),CapitalM(Edge(2,iElem)),CapitalM(Edge(3,iElem)),CapitalM(Edge(4,iElem)) )
        !EndIf
        !If (HydroParam%hb(iElem)<-0.4) Then
        !    HydroParam%hb(iElem) = -0.4
        !EndIf
        
        ! Lower Index - Superficial Flow Layer       
        Do iLayer = 1,MeshParam%KMax
            If (iLayer == MeshParam%KMax) Then
                If (HydroParam%hb(iElem)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
                    HydroParam%ElSmallm(iElem) = iLayer
                    exit
                EndIf
            Else
                If (HydroParam%hb(iElem)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%hb(iElem)<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
                    HydroParam%ElSmallm(iElem) = iLayer
                    exit
                EndIf
            EndIf
        EndDo
        
        ! Upper Index - Superficial Flow Layer
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
        
              
        ! 4.1.2.1 Compute the Vertical Mesh Spacing
        !Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem) !CAYO
        !    HydroParam%DZi(iLayer,iElem) = HydroParam%Ze(iLayer+1,iElem) - HydroParam%Ze(iLayer,iElem)
        !    HydroParam%DZi(iLayer,iElem) = Max(HydroParam%Pcri,HydroParam%DZi(iLayer,iElem))
        !EndDo
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            If (iLayer == HydroParam%ElSmallm(iElem)) Then
                HydroParam%Zb(iLayer,iElem) = (0.5d0)*( HydroParam%Ze(iLayer,iElem) + HydroParam%Ze(iLayer+1,iElem) ) !(0.5d0)*( 0. + Ze(iLayer+1,iElem) )   (verificar com Rafael)             
            Else
                HydroParam%Zb(iLayer,iElem) = (0.5d0)*( HydroParam%Ze(iLayer,iElem) + HydroParam%Ze(iLayer+1,iElem) )
            EndIf
        EndDo       
        
        ! Lower Index - Subsuperficial Flow Layer !CAYO        
        Do iLayer = 1,MeshParam%KMax
            If (iLayer == MeshParam%KMax) Then
                If (HydroParam%sb(iElem)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
                    HydroParam%ElSmallms(iElem) = iLayer
                    exit
                EndIf
            Else
                If (HydroParam%sb(iElem)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%sb(iElem)<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
                    HydroParam%ElSmallms(iElem) = iLayer
                    exit
                EndIf
            EndIf
        EndDo

        ! Upper Index - Superficial Flow Layer        
        If ( HydroParam%ElSmallm(iElem) > 1) Then   ! Cayo 
            HydroParam%ElCapitalMs(iElem) = HydroParam%ElSmallm(iElem) - 1
            if (HydroParam%ElCapitalMs(iElem)<HydroParam%ElSmallms(iElem)) then
                HydroParam%ElCapitalMs(iElem) = HydroParam%ElSmallms(iElem)
            endif
        Else
            HydroParam%ElCapitalMs(iElem) = HydroParam%ElSmallms(iElem)
        EndIf        
        
        ! 4.1.2 Update the Element Vertical Spacing
        
        Do iLayer = HydroParam%ElSmallms(iElem)+1, HydroParam%ElCapitalM(iElem) ! Cayo 
            HydroParam%Ze(iLayer,iElem) = MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)      ! Equidistant Core Grid
        EndDo
        HydroParam%Ze(HydroParam%ElSmallm(iElem),iElem)     = HydroParam%sb(iElem)                    ! Bottom
        If ( HydroParam%eta(iElem) - HydroParam%sb(iElem) <= HydroParam%PCRI+NearZero ) Then
            HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%hb(iElem) + HydroParam%PCRI !- hb(iElem)       ! Free-Surface (verificar com Rafael)
        Else
            HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%eta(iElem) !- hb(iElem)       ! Free-Surface (verificar com Rafael)    
        EndIf
        Do iLayer = 1, HydroParam%ElSmallm(iElem) - 1
            HydroParam%Ze(iLayer,iElem) = HydroParam%sb(iElem)
        EndDo
        
        ! 4.1.2.1 Compute the Vertical Mesh Spacing
        Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem) !CAYO
            HydroParam%DZi(iLayer,iElem) = HydroParam%Ze(iLayer+1,iElem) - HydroParam%Ze(iLayer,iElem)
            HydroParam%DZi(iLayer,iElem) = Max(HydroParam%Pcri,HydroParam%DZi(iLayer,iElem))
        EndDo
        
        
        !Porosity
        Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem) ! Cayo (Loop adicionado)
            If (iLayer >= HydroParam%ElSmallm(iElem) .or.  HydroParam%ElSmallms(iElem) ==  HydroParam%ElSmallm(iElem)) Then
                MeshParam%ei(iLayer,iElem) = 1
            Else
                MeshParam%ei(iLayer,iElem) = 0.3
            EndIf
        EndDo        
        
    EndDo     
    
    ! Set Initial Conditions of velocity components 
	Do iElem = 1, MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            Do iEdge = 1,4
                r = MeshParam%Right(MeshParam%Edge(iEdge,iElem))
                !uIniVec(iLayer,1) = -6*(EdgeBary(1,Edge(iEdge,iElem))/100.)**2. + 6*EdgeBary(1,Edge(iEdge,iElem))/100.
                !uIniVec(iLayer,2) = 0.d0
                If (r==0.or.HydroParam%H(MeshParam%Edge(iEdge,iElem))<=HydroParam%Pcri) Then ! No slip condition
                    HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem)) = 0.d0
                Else
                    If (simParam%it > 0) Then
                        HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem)) = simParam%usave(iLayer,MeshParam%Edge(iEdge,iElem))
                    Else
                        HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem))  = Dot_Product(HydroParam%uIniVec(iLayer,1:2),MeshParam%NormalVector(:,iEdge,iElem))
                    EndIf
                    
                EndIf
            EndDo
        EndDo
    EndDo
    If (simParam%it > 0) Then
        HydroParam%w = simParam%wsave
        Call Velocities(HydroParam,MeshParam)
    Else
        HydroParam%w = HydroParam%Wini
    EndIf
    
    If (HydroParam%iVTurb ==3 ) Then
        If(.NOT.Allocated(HydroParam%TKE))             Allocate(HydroParam%TKE(MeshParam%Kmax+1,MeshParam%nElem))
        If(.NOT.Allocated(HydroParam%LengthScale))     Allocate(HydroParam%LengthScale(MeshParam%Kmax+1,MeshParam%nElem))
        If(.NOT.Allocated(HydroParam%TKEP))            Allocate(HydroParam%TKEP(MeshParam%Kmax+1,MeshParam%nElem))
        If(.NOT.Allocated(HydroParam%LengthScaleP))    Allocate(HydroParam%LengthScaleP(MeshParam%Kmax+1,MeshParam%nElem))
        If(.NOT.Allocated(HydroParam%DissipRate))      Allocate(HydroParam%DissipRate(MeshParam%Kmax+1,MeshParam%nElem))
        If(.NOT.Allocated(HydroParam%ShearProd))       Allocate(HydroParam%ShearProd(MeshParam%Kmax+1,MeshParam%nElem))
        If(.NOT.Allocated(HydroParam%BuoyancyProd))    Allocate(HydroParam%BuoyancyProd(MeshParam%Kmax+1,MeshParam%nElem))
        
    EndIf
    
    
    !Open(1001, File ='VelocityEdge.txt', Action='Write', Status='Replace')
    !Do iEdge = 1, nEdge
    !    AuxVel = -6*(EdgeBary(1,iEdge)/100.)**2. + 6*EdgeBary(1,iEdge)/100.
    !    Write(1001,'(i8.10,4F20.10)') iEdge, uxy(1,1,iEdge), AuxVel, uxy(1,2,iEdge), 0.
    !EndDo
    !Close(1001)
    
    
    
End Subroutine ReadHydroIniCond
    