!> This subroutine reads the hydrodynamic parameters. 
Subroutine ReadHydroIniCond(HydroParam,hydroConfiguration,simParam,MeshParam)
       
    Use domain_types
    Use SimulationModel
    Use MeshVars
    Use Hydrodynamic
    Use Sorting
    Use iso_c_binding
    Use Meteorological
    Implicit none
    
    Integer:: i,j,iElem,iEdge,iNode,iLayer,r,l,Face,ie,Sig
    !character, pointer :: hydroParametersName(:)
    type(HydrodynamicConfiguration) :: hydroConfiguration  
    type(HydrodynamicParameter), pointer :: hydroParameters(:)
    type(SimulationParam) :: simParam
    character(len=200):: text
    Real:: sum1,sum0,AuxVel,V, SimTime, e0, k0
    Real:: NearZero = 1e-10
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam    
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
    HydroParam%CFL = 0.d0
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
    HydroParam%ug = 0.d0    
    HydroParam%vg = 0.d0
    HydroParam%wg = 0.d0
    HydroParam%ubV = 0.d0
    HydroParam%uxyL = 0.d0
    HydroParam%wfc = 0.d0
    MeshParam%Si = 0.d0
    MeshParam%ei = 1.d0
    HydroParam%DZK = 0.d0
    MeshParam%Kj = 0.d0
    HydroParam%DZsi = 0.d0
    HydroParam%DZhi = 0.d0
    HydroParam%DZhit = 0.d0
    HydroParam%DZhj = 0.d0
    HydroParam%DZsj = 0.d0
    HydroParam%us = 0.d0
    HydroParam%ustang = 0.d0
    HydroParam%um = 0.d0  
    HydroParam%wm = 0.d0 
    HydroParam%uxysub = 0.d0
    HydroParam%ubsub = 0.d0
    HydroParam%etaInfn = 0.d0
    HydroParam%etaInf = 0.d0
    MeshParam%ei = 1
    
    HydroParam%rhsnonHydro = 0.d0
    HydroParam%q = 0.d0
    HydroParam%pq = 0.d0
    SimTime = 0.d0
    
    MeshParam%iSaturation = 0
    
    !!!Bench 01:
    e0 = 0.3 !e0 0.1 b01 0.3
    MeshParam%Ksat = 0.01 !k0 0.01 b01 
    MeshParam%alpha = 29.0!Casulli == 1 Panday == 2.25
    MeshParam%nSoil = 4.0 !Casulli == 1.4 Panday == 1.89
    !
    !!!!Bench 02:
    !e0 = 0.9
    !MeshParam%Ksat  = 0.005 !k0 0.01 b01 0.00005    
    !MeshParam%alpha = 2.25 !Casulli == 1 Panday == 2.25
    !MeshParam%nSoil = 1.89 !Casulli == 1.4 Panday == 1.89    
    !
    !!!!Bench 03:
    !e0 = 0.5!0.2 !e0 0.1 b01 0.3
    !MeshParam%Ksat = 0.005 !0.0005 !k0 0.01 b01
    
    !k0 = 0.01 !k0 0.01 b01
    HydroParam%PsiCrit = 0.0d0
     
    ! Set Initial Conditions of Free-Surface Elevation
    If (MeshParam%ieta0 == 1) Then
        Do iElem = 1, MeshParam%nElem
            HydroParam%eta(iElem) = MeshParam%eta0(iElem)
        EndDo
    Else
        HydroParam%eta = HydroParam%Zini
    EndIf
    
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
            If (abs(HydroParam%sb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<HydroParam%PCRI) then
                HydroParam%sb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
            EndIf
        EndDo
        
      ! LIMCAMAUX stores bottom level from each layer defined by user. The following "for" set a correction in the hb/sb level in order to your values
      ! match with the bottom level from this layers in cases where the difference between hb/sb and LIMCAUX is lower PCRI.
        If (MeshParam%KMax==1) Then ! Two-dimensional
            Do iLayer = 1,MeshParam%KMax
                If (abs(HydroParam%hb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<HydroParam%PCRI) then
                    HydroParam%hb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
                    !HydroParam%hb(iElem) = 0.0d0
                    !HydroParam%hb(iElem)= Max(-0.33,0.5*abs(MeshParam%xb(iElem)-1.40 - 2.64) - 0.35)
                    !HydroParam%hb(iElem)= -dmin1(0.d0,HydroParam%hb(iElem))
                    if(HydroParam%hb(iElem)<=HydroParam%Pcri) then
                        HydroParam%hb(iElem) = 0.0d0
                    endif
                EndIf
                If (abs(HydroParam%sb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<HydroParam%PCRI) then
                    HydroParam%sb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
                EndIf   
            EndDo
            
        ElseIf (MeshParam%KMax>1) Then ! Three-dimensional
            Do iLayer = 1,MeshParam%KMax-1
                
                ! Check if hb is between two layers
                If (HydroParam%hb(iElem)+0.001>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%hb(iElem)+0.001<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) Then
                    ! If hb is closest the lower LIMCAMAUX (Kmax - 1 - ilayer) than upper LIMCAUX (Kmax - ilayer): hb = LIMCAUX lower
                    If (abs(HydroParam%hb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<abs(HydroParam%hb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer))) Then ! 
                        HydroParam%hb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
                    Else ! Else, hb is closest the upper layer
                        HydroParam%hb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)
                    EndIf
                EndIf
                ! Same check as hb:
                If (HydroParam%sb(iElem)+0.001>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%sb(iElem)+0.001<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) Then
                    If (abs(HydroParam%sb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<abs(HydroParam%sb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer))) Then ! 
                        HydroParam%sb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
                    Else
                        HydroParam%sb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)
                    EndIf
                EndIf
            EndDo
        EndIf
        
        !If (HydroParam%eta(iElem) - HydroParam%hb(iElem) <= HydroParam%Pcri+NearZero.and.HydroParam%eta(iElem) - HydroParam%hb(iElem) >= 0.) Then
        !    !HydroParam%hb(iElem) = HydroParam%eta(iElem) - HydroParam%Pcri/2.d0
        !    HydroParam%hb(iElem) = HydroParam%eta(iElem)
        !EndIf
        
        !If (HydroParam%eta(iElem) - HydroParam%sb(iElem) <= HydroParam%Pcri+NearZero) Then
        If (HydroParam%eta(iElem) - HydroParam%sb(iElem) <=NearZero) Then
            !HydroParam%eta(iElem) = HydroParam%sb(iElem) + HydroParam%Pcri/2.d0 !*** Verificar sinal posteriormente, se já entrar como cota não precisa do sinal 
            HydroParam%eta(iElem) = HydroParam%sb(iElem)!*** Verificar sinal posteriormente, se já entrar como cota não precisa do sinal 
        EndIf
    EndDo


    If (MeshParam%iBedrock == 0) Then
        HydroParam%sb =  HydroParam%hb
    EndIf
    
    !Compute nodal elevations for tang. vel.
    Do iNode=1,MeshParam%nNode
        sum1 = 0.d0 !sum of elemental elevations
        sum0 = 0.d0 !sum of areas
        do j=1,MeshParam%nVertexElem(iNode)
            ie=MeshParam%VertexElem(j,iNode)
            sum1 = sum1 + MeshParam%Area(ie)*HydroParam%eta(ie)
            sum0 = sum0 + MeshParam%Area(ie)
        Enddo !j=1,nne(i)
        HydroParam%peta(iNode) = sum1/sum0
    EndDo !i=1,np   
    HydroParam%petan = HydroParam%peta !store for transport eqs.
    
    ! Set Bathymetry in th edges
    Do iElem = 1, MeshParam%nElem 
        Do iEdge = 1, 4
            Face = MeshParam%Edge(iEdge,iElem)
            l = MeshParam%Left(Face) 
            r = MeshParam%Right(Face)
            If (r == 0) Then
                HydroParam%hj(Face) = HydroParam%hb(l) !*** Verificar sinal posteriormente, se j� entrar como cota n�o precisa do sinal 
                HydroParam%sj(Face) = HydroParam%sb(l)
                Do iLayer = 1,MeshParam%KMax
                    If (abs(HydroParam%hj(Face)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<=HydroParam%PCRI+NearZero) then
                        HydroParam%hj(Face)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
                    EndIf
                    If (abs(HydroParam%sj(Face)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<=HydroParam%PCRI+NearZero) then
                        HydroParam%sj(Face)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
                    EndIf
                EndDo        
            Else
                HydroParam%hj(Face) = Max(HydroParam%hb(l),HydroParam%hb(r))  !  (0.5d0)*( Hbat(I,J) + HBat(IViz,JViz) ) !Max(Hbat(I,J),HBat(IViz,JViz))!*** !(0.5d0)*( Hbat(I,J) + HBat(IViz,JViz) ) Verificar sinal posteriormente, se j� entrar como cota n�o precisa do sinal 
                HydroParam%sj(Face) = Max(HydroParam%sb(l),HydroParam%sb(r))              
                !HydroParam%hj(Face) = 0.5*(HydroParam%hb(l)+HydroParam%hb(r))  !  (0.5d0)*( Hbat(I,J) + HBat(IViz,JViz) ) !Max(Hbat(I,J),HBat(IViz,JViz))!*** !(0.5d0)*( Hbat(I,J) + HBat(IViz,JViz) ) Verificar sinal posteriormente, se j� entrar como cota n�o precisa do sinal 
                !HydroParam%sj(Face) = 0.5*(HydroParam%sb(l)+HydroParam%sb(r))
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
        HydroParam%ElCapitalM(iElem) = HydroParam%ElSmallm(iElem)
        If(HydroParam%eta(iElem) - HydroParam%hb(iElem) > 0.d0) Then
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
        EndIf                
        
        ! 4.1.2.1 Compute the Vertical Mesh Spacing
        !Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem) !CAYO
        !    HydroParam%DZi(iLayer,iElem) = HydroParam%Ze(iLayer+1,iElem) - HydroParam%Ze(iLayer,iElem)
        !    HydroParam%DZi(iLayer,iElem) = Max(HydroParam%Pcri,HydroParam%DZi(iLayer,iElem))
        !EndDo

        HydroParam%ElSmallms(iElem)     = HydroParam%ElSmallm(iElem)
        HydroParam%ElCapitalMs(iElem)   = HydroParam%ElSmallm(iElem)
		If (MeshParam%iBedrock == 1) Then
        ! Lower Index - Subsuperficial Flow Layer       
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
            If ( HydroParam%ElSmallm(iElem) > 1) Then
                HydroParam%ElCapitalMs(iElem) = HydroParam%ElSmallm(iElem) - 1
                if (HydroParam%ElCapitalMs(iElem)<HydroParam%ElSmallms(iElem)) then
                    HydroParam%ElCapitalMs(iElem) = HydroParam%ElSmallms(iElem)
                endif
            Else
                HydroParam%ElCapitalMs(iElem) = HydroParam%ElSmallms(iElem)
            EndIf
		EndIf
                      
       
        !CAYO
        ! 4.1.2 Update the Element Vertical Spacing
        !If ( HydroParam%eta(iElem) - HydroParam%hb(iElem) <= HydroParam%PCRI+NearZero ) Then
        If ( HydroParam%eta(iElem) - HydroParam%hb(iElem) <= NearZero ) Then
            !HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%hb(iElem) + HydroParam%PCRI !- hb(iElem)       ! Free-Surface (verificar com Rafael)
            HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%hb(iElem) !- hb(iElem)       ! Free-Surface (verificar com Rafael)
        Else
            HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%eta(iElem) !- hb(iElem)       ! Free-Surface (verificar com Rafael)    
            !HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = Max(HydroParam%eta(iElem), HydroParam%Ze(HydroParam%ElCapitalM(iElem),iElem)     
        EndIf
        HydroParam%Ze(:,iElem) = HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem)
        
        Do iLayer = HydroParam%ElSmallms(iElem)+1, HydroParam%ElCapitalM(iElem)
            HydroParam%Ze(iLayer,iElem) = MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)      ! Equidistant Core Grid
        EndDo
        HydroParam%Ze(HydroParam%ElSmallms(iElem),iElem)     = HydroParam%sb(iElem)                    ! Bottom
        
        Do iLayer = 1, HydroParam%ElSmallms(iElem) - 1
            HydroParam%Ze(iLayer,iElem) = HydroParam%sb(iElem)
        EndDo
        
        Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem)
            HydroParam%Zb(iLayer,iElem) = (0.5d0)*( HydroParam%Ze(iLayer,iElem) + HydroParam%Ze(iLayer+1,iElem) )
        EndDo   

        ! 4.1.2.1 Compute the Vertical Mesh Spacing
        Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem)
            HydroParam%DZi(iLayer,iElem) = HydroParam%Ze(iLayer+1,iElem) - HydroParam%Ze(iLayer,iElem)
            HydroParam%DZi(iLayer,iElem) = Max(0.0d0,HydroParam%DZi(iLayer,iElem))
            !HydroParam%DZi(iLayer,iElem) = Max(HydroParam%Pcri,HydroParam%DZi(iLayer,iElem))
        EndDo
        
        !xx. Porosity and and DZhi/DZsi
        HydroParam%DZhi(:,iElem) = HydroParam%DZi(:,iElem)
		If (MeshParam%iBedrock == 1) Then
            
            !If(MeshParam%iSaturation == 0) Then !Darcy's Model
            !    HydroParam%PsiCrit(iElem) = 0.0
            !ElseIf(MeshParam%iSaturation == 1) Then ! Brooks and Corey's Model
            !    HydroParam%PsiCrit(iElem) = 1/MeshParam%alpha(HydroParam%ElCapitalMs(iElem),iElem)
            !Else !van Genuchten's Model
            !    HydroParam%PsiCrit(iElem) = ((MeshParam%nSoil(HydroParam%ElCapitalMs(iElem),iElem) - 1)/MeshParam%nSoil(HydroParam%ElCapitalMs(iElem),iElem))**(1/MeshParam%nSoil(HydroParam%ElCapitalMs(iElem),iElem))/(MeshParam%alpha(HydroParam%ElCapitalMs(iElem),iElem)*MeshParam%nSoil(HydroParam%ElCapitalMs(iElem),iElem))        
            !    !HydroParam%PsiCrit(iElem) = etatol
            !EndIf            
            
            Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem)
            
                !If (iLayer >= HydroParam%ElSmallm(iElem) .or.  HydroParam%ElSmallms(iElem) ==  HydroParam%ElSmallm(iElem)) Then
                !    continue
                !Else
                !    MeshParam%ei(iLayer,iElem) = e0
                !EndIf

                If (iLayer < HydroParam%ElSmallm(iElem)) Then
                    HydroParam%DZhi(iLayer,iElem) = 0.d0
                    HydroParam%DZsi(iLayer,iElem) = HydroParam%DZi(iLayer,iElem)               
                ElseIf (iLayer > HydroParam%ElSmallm(iElem)) Then
                    HydroParam%DZhi(iLayer,iElem) = HydroParam%DZi(iLayer,iElem)
                    HydroParam%DZsi(iLayer,iElem) = 0.d0
                Else
                    If (HydroParam%Ze(iLayer+1,iElem) > HydroParam%hb(iElem) ) Then
                        HydroParam%DZhi(iLayer,iElem) = HydroParam%Ze(iLayer+1,iElem) - HydroParam%hb(iElem)
                        HydroParam%DZsi(iLayer,iElem) = HydroParam%hb(iElem) - HydroParam%Ze(iLayer,iElem)
                    Else
                        HydroParam%DZhi(iLayer,iElem) = 0.d0
                        HydroParam%DZsi(iLayer,iElem) = HydroParam%Ze(iLayer+1,iElem) - HydroParam%Ze(iLayer,iElem)
                        !If(HydroParam%ElSmallms(iElem) == HydroParam%ElCapitalM(iElem)) Then
                        !    If(HydroParam%eta(iElem) < HydroParam%hb(iElem) ) Then
                        !        HydroParam%DZsi(iLayer,iElem) = HydroParam%eta(iElem)
                        !    ElseIf (HydroParam%eta(iElem) <= HydroParam%sb(iElem)) Then
                        !        HydroParam%DZsi(iLayer,iElem) = 0.d0
                        !    EndIf
                        !EndIf                        
                    EndIf
                EndIf                                    
                
                If(HydroParam%DZsi(iLayer,iElem)>0) Then
                    MeshParam%ei(iLayer,iElem) = e0
                EndIf                  
            EndDo
        EndIf
        
        Call MoistureContent(HydroParam%eta(iElem),0.d0,iElem,HydroParam,MeshParam)  
        
    EndDo  
    
    
    !Compute Index Smallm(j) and CapitalM(j)
    ! Smallm and CapitalM are related to mj and Mj in [1]
    Do iEdge = 1,MeshParam%nEdge

        ! 1.1 Compute Water Depth
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
        
        If (r==0) Then
            r = l
        EndIf
        !
        !If (r == 0) Then
        !    HydroParam%H(iEdge) = Max( HydroParam%PCRI, -HydroParam%sj(iEdge) + HydroParam%eta(l) )
        !Else
        !    HydroParam%H(iEdge) = Max( HydroParam%PCRI,-HydroParam%sj(iEdge) + HydroParam%eta(l), -HydroParam%sj(iEdge) + HydroParam%eta(r) )!CAYO
        !EndIf
        !HydroParam%H(iEdge) = Max( HydroParam%PCRI,-HydroParam%sj(iEdge) + HydroParam%eta(l), -HydroParam%sj(iEdge) + HydroParam%eta(r) )!CAYO
        HydroParam%H(iEdge) = Max(0.d0,-HydroParam%sj(iEdge) + HydroParam%eta(l), -HydroParam%sj(iEdge) + HydroParam%eta(r) )!CAYO
        
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
        HydroParam%CapitalM(iEdge) = HydroParam%Smallm(iEdge)
        Do iLayer = 1,MeshParam%KMax
            If (iLayer == MeshParam%KMax) Then
                !If (Max(HydroParam%PCRI, HydroParam%eta(l), HydroParam%eta(r))>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
                If (Max(0.d0, HydroParam%eta(l), HydroParam%eta(r))>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
                    HydroParam%CapitalM(iEdge) = Max(iLayer,HydroParam%CapitalM(iEdge))
                    exit
                EndIf
            Else
                !If (Max(HydroParam%PCRI, HydroParam%eta(l), HydroParam%eta(r))>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.Max(HydroParam%PCRI, HydroParam%eta(l), HydroParam%eta(r))<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
                If (Max(0.d0, HydroParam%eta(l), HydroParam%eta(r))>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.Max(HydroParam%PCRI, HydroParam%eta(l), HydroParam%eta(r))<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
                    HydroParam%CapitalM(iEdge) = Max(iLayer,HydroParam%CapitalM(iEdge))
                    exit
                EndIf
            EndIf
        EndDo
        !
        !Do iLayer = 1,MeshParam%KMax
        !    If (iLayer == MeshParam%KMax) Then
        !        If (HydroParam%H(iEdge) + HydroParam%hj(iEdge)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
        !            HydroParam%CapitalM(iEdge) = iLayer
        !            exit
        !        EndIf
        !    Else
        !        If (HydroParam%H(iEdge) + HydroParam%hj(iEdge)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%H(iEdge) -HydroParam%sj(iEdge)+ HydroParam%hj(iEdge)<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
        !            HydroParam%CapitalM(iEdge) = iLayer
        !            exit
        !        EndIf
        !    EndIf
        !EndDo
        
        HydroParam%Smallms(iEdge) = HydroParam%Smallm(iEdge)
        HydroParam%CapitalMs(iEdge) = HydroParam%Smallm(iEdge)   
        If (MeshParam%iBedrock == 1) Then
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

            ! Upper Index - Subsuperficial Flow Layer
            If ( HydroParam%Smallm(iEdge) > 1) Then 
                HydroParam%CapitalMs(iEdge) = HydroParam%Smallm(iEdge) - 1
                if (HydroParam%CapitalMs(iEdge)<HydroParam%Smallms(iEdge)) then
                    HydroParam%CapitalMs(iEdge) = HydroParam%Smallms(iEdge)
                endif
            Else
                HydroParam%CapitalMs(iEdge) = HydroParam%Smallms(iEdge)
            EndIf
        EndIf  
        !CAYO
        !9.2 Compute Elevation in the Edges
        Do iLayer = HydroParam%Smallms(iEdge) + 1, HydroParam%CapitalM(iEdge)
            HydroParam%Z(iLayer,iEdge) = MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1) !zL + (iLayer - 1)*dz                  ! Equidistant Core Grid
        EndDo
        HydroParam%Z(HydroParam%Smallms(iEdge),iEdge)     = HydroParam%sj(iEdge)                 ! Bottom
        If (r == 0) Then
            !If ( HydroParam%eta(l) - HydroParam%hj(iEdge) <= HydroParam%PCRI+NearZero) Then
            If ( HydroParam%eta(l) - HydroParam%hj(iEdge) <= NearZero) Then
                !HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%hj(iEdge) + HydroParam%PCRI
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%hj(iEdge)
            Else
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%eta(l)
            EndIf
             !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
        Else
            !If ( HydroParam%H(iEdge) - HydroParam%hj(iEdge) <= HydroParam%PCRI+NearZero ) Then
            !If ( Max(HydroParam%PCRI, HydroParam%eta(l), HydroParam%eta(r)) - HydroParam%hj(iEdge) <= HydroParam%PCRI+NearZero ) Then !( HydroParam%H(iEdge) <= HydroParam%PCRI/2 + NearZero ) Then
            If ( Max(0.0d0, HydroParam%eta(l), HydroParam%eta(r)) - HydroParam%hj(iEdge) <= NearZero ) Then !( HydroParam%H(iEdge) <= HydroParam%PCRI/2 + NearZero ) Then
                !HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%eta(l),HydroParam%eta(r))
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%hb(l),HydroParam%hb(r))                
            Else
                !HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = (0.5d0)*(HydroParam%eta(l)+HydroParam%eta(r)) !Max(HydroParam%eta(l),HydroParam%eta(r)) !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%eta(l),HydroParam%eta(r))             
            EndIf
        EndIf        
                
        ! 1.4 Compute the Vertical Mesh Spacing
        Do iLayer = HydroParam%Smallms(iEdge), HydroParam%CapitalM(iEdge)
            HydroParam%DZj(iLayer,iEdge) = HydroParam%Z(iLayer+1,iEdge) - HydroParam%Z(iLayer,iEdge)
            HydroParam%DZj(iLayer,iEdge) = Max(HydroParam%Pcri,HydroParam%DZj(iLayer,iEdge))
        EndDo        
        
        !1.5 Compute Kj and DZhj/DZsj
		HydroParam%DZhj(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge) = HydroParam%DZj(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge)
        HydroParam%DZsj(:,iEdge) = 0.d0
        MeshParam%Kj(:,iEdge) = 0.d0
        HydroParam%DZK(iEdge) = 0.d0
        If (MeshParam%iBedrock == 1) Then
            Do iLayer = HydroParam%Smallms(iEdge), HydroParam%CapitalM(iEdge) ! 
                    
                If (HydroParam%H(iEdge) > HydroParam%Pcri) Then
                    If (iLayer < HydroParam%Smallm(iEdge)) Then
                        HydroParam%DZhj(iLayer,iEdge) = 0.
                        HydroParam%DZsj(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)
                    ElseIf (iLayer > HydroParam%Smallm(iEdge)) Then
                        HydroParam%DZhj(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)
                        HydroParam%DZsj(iLayer,iEdge) = 0.
                        continue
                    Else
                        If ( HydroParam%Z(iLayer+1,iEdge) > HydroParam%hj(iEdge) ) Then
                            HydroParam%DZhj(iLayer,iEdge) = HydroParam%Z(iLayer+1,iEdge) - HydroParam%hj(iEdge)
                            HydroParam%DZsj(iLayer,iEdge) = HydroParam%hj(iEdge) - HydroParam%Z(iLayer,iEdge)
                        Else
                            HydroParam%DZhj(iLayer,iEdge) = 0.
                            HydroParam%DZsj(iLayer,iEdge) = Max(HydroParam%Pcri,HydroParam%Z(iLayer+1,iEdge) - HydroParam%Z(iLayer,iEdge))
                            !If(HydroParam%Smallms(iEdge) == HydroParam%CapitalMs(iEdge) .and. MeshParam%iSaturation == 0) Then
                            !    If(HydroParam%H(iEdge) < HydroParam%hj(iEdge) ) Then
                            !        HydroParam%DZsj(iLayer,iEdge) = HydroParam%H(iEdge)
                            !    ElseIf(HydroParam%H(iEdge) <= HydroParam%sj(iEdge))Then
                            !        HydroParam%DZsj(iLayer,iEdge) = 0.d0
                            !    EndIf   
                            !EndIf                       
                        EndIf
                    EndIf    
                
                    If(HydroParam%Smallms(iEdge) == HydroParam%CapitalMs(iEdge) .and. MeshParam%iSaturation == 0) Then
                        If(HydroParam%H(iEdge) < HydroParam%hj(iEdge) ) Then
                            HydroParam%DZsj(iLayer,iEdge) = HydroParam%H(iEdge)
                        ElseIf(HydroParam%H(iEdge) <= HydroParam%sj(iEdge))Then
                            HydroParam%DZsj(iLayer,iEdge) = 0.d0
                        EndIf   
                    EndIf                   
                
                    If (HydroParam%DZsj(iLayer,iEdge) > 0) Then 
                        !If Neighbour lower layer is above Elem layer, the lower Element layer is 
                        If (iLayer < HydroParam%ElSmallms(r) .and. MeshParam%Ki(iLayer,l) > 0) Then
                            MeshParam%Kj(iLayer,iEdge) = MeshParam%Ki(iLayer,l)
                        ElseIf (iLayer < HydroParam%ElSmallms(l) .and. MeshParam%Ki(iLayer,r) > 0) Then     
                            MeshParam%Kj(iLayer,iEdge) = MeshParam%Ki(iLayer,r)                
                        Else
                            MeshParam%Kj(iLayer,iEdge) = Max(MeshParam%Ki(iLayer,l),MeshParam%Ki(iLayer,r))
                        EndIf
                    EndIf
                
                    HydroParam%DZK(iEdge) = HydroParam%DZK(iEdge) + HydroParam%DZsj(iLayer,iEdge)*MeshParam%Kj(iLayer,iEdge) !Sediment Layer
                EndIf
            EndDo
		EndIf
        
    EndDo
   
    
    ! Set Initial Conditions of velocity components 
	Do iElem = 1, MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            Do iEdge = 1,4
                r = MeshParam%Right(MeshParam%Edge(iEdge,iElem))
                !uIniVec(iLayer,1) = -6*(EdgeBary(1,Edge(iEdge,iElem))/100.)**2. + 6*EdgeBary(1,Edge(iEdge,iElem))/100.
                !uIniVec(iLayer,2) = 0.d0
                !If (r==0.or.HydroParam%H(MeshParam%Edge(iEdge,iElem))<=HydroParam%Pcri+NearZero) Then ! No slip condition
                If (r==0.or.HydroParam%H(MeshParam%Edge(iEdge,iElem))<=NearZero) Then ! No slip condition
                    HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem)) = 0.d0
                Else
                    If (simParam%it > 0) Then
                        HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem)) = simParam%usave(iLayer,MeshParam%Edge(iEdge,iElem))
                        !HydroParam%utang(iLayer,MeshParam%Edge(iEdge,iElem)) = simParam%utangsave(iLayer,MeshParam%Edge(iEdge,iElem))
                    Else
                        HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem))  = Dot_Product(HydroParam%uIniVec(iLayer,1:2),MeshParam%NormalVector(:,iEdge,iElem))
                        HydroParam%utang(iLayer,MeshParam%Edge(iEdge,iElem))  = Dot_Product(HydroParam%uIniVec(iLayer,1:2),MeshParam%TangentVector(:,iEdge,iElem))
                    EndIf
                    
                EndIf
            EndDo
        EndDo
    EndDo
    If (simParam%it > 0) Then
        HydroParam%w = simParam%wsave
        Call VelocitiesSub(HydroParam,MeshParam)
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
    
!    
!    
!    
!!> This subroutine reads the hydrodynamic parameters. 
!Subroutine ReadHydroIniCond(HydroParam,hydroConfiguration,simParam,MeshParam)
!    
!    Use domain_types
!    Use SimulationModel
!    Use MeshVars
!    Use Hydrodynamic
!    Use Sorting
!    Use iso_c_binding
!    
!    Implicit none
!    
!    Integer:: i,j,iElem,iEdge,iNode,iLayer,r,l,Face,ie,Sig
!    !character, pointer :: hydroParametersName(:)
!    type(HydrodynamicConfiguration) :: hydroConfiguration  
!    type(HydrodynamicParameter), pointer :: hydroParameters(:)
!    type(SimulationParam) :: simParam
!    character(len=200):: text
!    Real:: sum1,sum0,AuxVel,V, SimTime, e0, k0
!    Real:: NearZero = 1e-10
!    type(MeshGridParam) :: MeshParam
!    type(HydrodynamicParam) :: HydroParam
!    Real, Dimension(4,MeshParam%nElem):: HNeighbor
!    Integer, Dimension(MeshParam%nElem):: nNeighbor
!    Integer:: npit
!    
!    call c_f_pointer(hydroConfiguration%parameters, hydroParameters, [hydroConfiguration%numberOfParameters])
!    
!    Do i = 1, hydroConfiguration%numberOfParameters
!        text = trim(hydroParameters(i)%name)
!
!        !Flags
!        If (trim(text) == 'envFreeSurfaceElevation') Then
!            HydroParam%Zini = hydroParameters(i)%value
!        ElseIf (trim(text) == 'envWaterVelocityX') Then
!            HydroParam%Uini  = hydroParameters(i)%value
!        ElseIf (trim(text) == 'envWaterVelocityY') Then
!            HydroParam%Vini = hydroParameters(i)%value
!        ElseIf (trim(text) == 'envWaterVelocityW') Then
!            HydroParam%Wini = hydroParameters(i)%value
!        EndIf
!        !print *, hydroParametersName, hydroParameters(i)%value
!    EndDo
!    HydroParam%uIniVec(:,1) = HydroParam%Uini
!    HydroParam%uIniVec(:,2) = HydroParam%Vini
!    HydroParam%uIniVec(:,3) = HydroParam%Wini
!    HydroParam%uxyback = 0.d0
!    HydroParam%uxy = 0.d0
!    HydroParam%Fu = 0.d0
!    HydroParam%Fv = 0.d0
!    HydroParam%Fub = 0.d0
!    HydroParam%FuxyNode = 0.d0
!    HydroParam%u = 0.d0
!    HydroParam%utang = 0.d0
!    HydroParam%ub = 0.d0
!    HydroParam%uNode = 0.d0
!    HydroParam%HorViscosity = 0.d0
!    HydroParam%HorDiffusivity = 0.d0
!    If ( HydroParam%iVTurb == 0 ) Then 
!        HydroParam%VerEddyVisc = HydroParam%VerEddyVisc_Cte
!        HydroParam%VerEddyDiff = HydroParam%VerEddyDiff_Cte
!        HydroParam%VerEddyViscCell = HydroParam%VerEddyVisc_Cte
!        HydroParam%VerEddyDiffCell = HydroParam%VerEddyDiff_Cte
!    Elseif ( HydroParam%iVTurb == 1 ) Then     
!        HydroParam%VerEddyVisc = HydroParam%vmin
!        HydroParam%VerEddyDiff = HydroParam%tdmin_pp
!        HydroParam%VerEddyViscCell = HydroParam%vmin
!        HydroParam%VerEddyDiffCell = HydroParam%tdmin_pp
!    EndIf
!    HydroParam%iAG = 0.d0
!    HydroParam%iADZ = 1.
!    HydroParam%DZj = HydroParam%Pcri
!    HydroParam%DZjt = HydroParam%Pcri
!    HydroParam%DZi = HydroParam%Pcri
!    HydroParam%DZit = HydroParam%Pcri
!    HydroParam%Gu = 0.d0
!    
!    !CAYO
!    HydroParam%ug = 0.d0    
!    HydroParam%vg = 0.d0
!    HydroParam%wg = 0.d0
!    HydroParam%ubV = 0.d0
!    HydroParam%uxyL = 0.d0
!    HydroParam%wfc = 0.d0
!    
!    MeshParam%ei = 1.d0
!    HydroParam%DZK = 0.d0
!    MeshParam%Kj = 0.d0
!    HydroParam%DZsi = HydroParam%Pcri
!    HydroParam%DZhi = HydroParam%Pcri
!    HydroParam%DZhit = HydroParam%Pcri
!    HydroParam%DZhj = HydroParam%Pcri
!    HydroParam%DZsj = HydroParam%Pcri
!    HydroParam%us = 0.d0
!    HydroParam%ustang = 0.d0
!    HydroParam%um = 0.d0  
!    HydroParam%wm = 0.d0 
!    HydroParam%uxysub = 0.d0
!    HydroParam%ubsub = 0.d0
!    HydroParam%etaInfn = 0.d0
!    HydroParam%etaInf = 0.d0
!    
!    HydroParam%rhsnonHydro = 0.d0
!    HydroParam%q = 0.d0
!    HydroParam%pq = 0.d0
!    SimTime = 0.d0
!    
!    e0 = 0.1 !e0 0.3 b01
!    k0 = 0.00005 !k0 0.01 b01
!    !e0 = 0.3 !e0 0.3 b01
!    !k0 = 0.01 !k0 0.01 b01
!     
!    ! Set Initial Conditions of Free-Surface Elevation
!    If (MeshParam%ieta0 == 1) Then
!        Do iElem = 1, MeshParam%nElem
!            HydroParam%eta(iElem) = MeshParam%eta0(iElem)
!        EndDo
!    Else
!        HydroParam%eta = HydroParam%Zini
!    EndIf
!    
!    If (simParam%it > 0) Then
!        HydroParam%eta = simParam%etasave
!    EndIf
!    HydroParam%etan = HydroParam%eta
!    
!    !!Filling sinks
!    !npit = 1 
!    !Do While (npit /= 0)   
!    !    
!    !    Do iElem = 1, MeshParam%nElem
!    !        nNeighbor(iElem) = 0
!    !        Do iEdge = 1, 4
!    !            Face = MeshParam%Edge(iEdge,iElem)
!    !            l = MeshParam%Left(Face) 
!    !            r = MeshParam%Right(Face)
!    !            If (r == 0) Then 
!    !                nNeighbor(iElem) = nNeighbor(iElem)
!    !            Else
!    !                nNeighbor(iElem) = nNeighbor(iElem) +1
!    !                HNeighbor(nNeighbor(iElem),iElem) = HydroParam%hb(MeshParam%Neighbor(iEdge,iElem))
!    !            EndIf
!    !        EndDo
!    !        
!    !    EndDo
!    !    npit = 0
!    !    Do iElem = 1, MeshParam%nElem
!    !        Call Sort(HNeighbor(:,iElem),nNeighbor(iElem))
!    !        If (nNeighbor(iElem) == 0) Then
!    !            HydroParam%hb(iElem) = HydroParam%hb(iElem)
!    !        Else
!    !            If (HydroParam%hb(iElem)<MinVal(HNeighbor(1:nNeighbor(iElem),iElem))) Then
!    !                HydroParam%hb(iElem) = MinVal(HNeighbor(1:nNeighbor(iElem),iElem))
!    !                npit = npit + 1
!    !            Elseif (HydroParam%hb(iElem)==MinVal(HNeighbor(1:nNeighbor(iElem),iElem))) Then
!    !                If (HydroParam%hb(iElem)<MinVal(HNeighbor(2:nNeighbor(iElem),iElem)).and.nNeighbor(iElem)>=2) Then
!    !                    HydroParam%hb(iElem) = MinVal(HNeighbor(2:nNeighbor(iElem),iElem))
!    !                    npit = npit + 1
!    !                Else
!    !                    HydroParam%hb(iElem) = HydroParam%hb(iElem)
!    !                EndIf
!    !            Else
!    !                HydroParam%hb(iElem) = HydroParam%hb(iElem)
!    !            EndIf
!    !        EndIf
!    !        
!    !    EndDo
!    !EndDo
!    !
!    
!    Do iElem = 1, MeshParam%nElem
!        
!        Do iLayer = 1,MeshParam%KMax
!            If (abs(HydroParam%hb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<HydroParam%PCRI) then
!                HydroParam%hb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
!            EndIf
!        EndDo
!        
!        Do iLayer = 1,MeshParam%KMaxsub
!            If (abs(HydroParam%sb(iElem)-MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer+1))<HydroParam%PCRI) then
!                HydroParam%sb(iElem)=MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer+1)
!            EndIf
!        EndDo        
!        
!      ! LIMCAMAUX stores bottom level from each layer defined by user. The following statements set a correction in the hb/sb level in order to your values
!      ! match with the bottom level from this layers in cases where the difference between hb/sb and LIMCAUX is lower PCRI.
!        !Surface layers:
!        If (MeshParam%KMax==1) Then ! Two-dimensional
!            Do iLayer = 1,MeshParam%KMax
!                If (abs(HydroParam%hb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<HydroParam%PCRI) then
!                    HydroParam%hb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
!                EndIf      
!        ElseIf (MeshParam%KMax>1) Then ! Three-dimensional
!            Do iLayer = 1,MeshParam%KMax-1 
!                ! Check if hb is between two layers
!                If (HydroParam%hb(iElem)+0.001>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%hb(iElem)+0.001<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) Then
!                    ! If hb is closest the lower LIMCAMAUX (Kmax - 1 - ilayer) than upper LIMCAUX (Kmax - ilayer): hb = LIMCAUX lower
!                    If (abs(HydroParam%hb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<abs(HydroParam%hb(iElem)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer))) Then ! 
!                        HydroParam%hb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
!                    Else ! Else, hb is closest the upper layer
!                        HydroParam%hb(iElem)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)
!                    EndIf
!                EndIf
!            EndDo
!        EndIf
!        
!        !Subsurface layers:
!        If (MeshParam%KMaxsub==1) Then ! Two-dimensional or non-existent subsurface flow
!            Do iLayer = 1,MeshParam%KMaxsub
!                If (abs(HydroParam%sb(iElem)-MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer+1))<HydroParam%PCRI) then
!                    HydroParam%sb(iElem)=MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer+1)
!                EndIf   
!            EndDo 
!        ElseIf (MeshParam%KMaxsub>1) Then ! Three-dimensional
!            Do iLayer = 1,MeshParam%KMaxsub-1
!                ! Same check for sb:
!                If (HydroParam%sb(iElem)+0.001>=MeshParam%LIMCAMAUXsub(MeshParam%KMax-iLayer+1).and.HydroParam%sb(iElem)+0.001<MeshParam%LIMCAMAUXsub(MeshParam%KMax-iLayer)) Then
!                    If (abs(HydroParam%sb(iElem)-MeshParam%LIMCAMAUXsub(MeshParam%KMax-iLayer+1))<abs(HydroParam%sb(iElem)-MeshParam%LIMCAMAUXsub(MeshParam%KMax-iLayer))) Then ! 
!                        HydroParam%sb(iElem)=MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer+1)
!                    Else
!                        HydroParam%sb(iElem)=MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer)
!                    EndIf
!                EndIf
!            EndDo
!        EndIf
!        
!        If (HydroParam%eta(iElem) - HydroParam%hb(iElem) <= HydroParam%Pcri+NearZero.and.HydroParam%eta(iElem) - HydroParam%hb(iElem) >= 0.) Then
!            HydroParam%hb(iElem) = HydroParam%eta(iElem) - HydroParam%Pcri/2.d0
!        EndIf
!        
!        If (HydroParam%eta(iElem) - HydroParam%sb(iElem) <= HydroParam%Pcri+NearZero) Then
!            HydroParam%eta(iElem) = HydroParam%sb(iElem) + HydroParam%Pcri/2.d0 !*** Verificar sinal posteriormente, se j� entrar como cota n�o precisa do sinal 
!        EndIf
!
!    EndDo
!
!
!    If (MeshParam%iBedrock == 0) Then
!            HydroParam%sb =  HydroParam%hb
!    EndIf
!    
!    !Compute nodal elevations for tang. vel.
!    Do iNode=1,MeshParam%nNode
!        sum1 = 0.d0 !sum of elemental elevations
!        sum0 = 0.d0 !sum of areas
!        do j=1,MeshParam%nVertexElem(iNode)
!            ie=MeshParam%VertexElem(j,iNode)
!            sum1 = sum1 + MeshParam%Area(ie)*HydroParam%eta(ie)
!            sum0 = sum0 + MeshParam%Area(ie)
!        Enddo !j=1,nne(i)
!        HydroParam%peta(iNode) = sum1/sum0
!    EndDo !i=1,np   
!    HydroParam%petan = HydroParam%peta !store for transport eqs.
!    
!    ! Set Bathymetry in th edges
!    Do iElem = 1, MeshParam%nElem
!        Do iEdge = 1, 4
!            Face = MeshParam%Edge(iEdge,iElem)
!            l = MeshParam%Left(Face) 
!            r = MeshParam%Right(Face)
!            If (r == 0) Then
!                HydroParam%hj(Face) = HydroParam%hb(l) !*** Verificar sinal posteriormente, se j� entrar como cota n�o precisa do sinal 
!                HydroParam%sj(Face) = HydroParam%sb(l)
!                Do iLayer = 1,MeshParam%KMax
!                    If (abs(HydroParam%hj(Face)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<=HydroParam%PCRI+NearZero) then
!                        HydroParam%hj(Face)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
!                    EndIf
!                EndDo   
!                Do iLayer = 1,MeshParam%KMaxsub
!                    If (abs(HydroParam%sj(Face)-MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer+1))<=HydroParam%PCRI+NearZero) then
!                        HydroParam%sj(Face)=MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer+1)
!                    EndIf
!                EndDo                  
!            Else
!                HydroParam%hj(Face) = Max(HydroParam%hb(l),HydroParam%hb(r))  !  (0.5d0)*( Hbat(I,J) + HBat(IViz,JViz) ) !Max(Hbat(I,J),HBat(IViz,JViz))!*** !(0.5d0)*( Hbat(I,J) + HBat(IViz,JViz) ) Verificar sinal posteriormente, se j� entrar como cota n�o precisa do sinal 
!                HydroParam%sj(Face) = Max(HydroParam%sb(l),HydroParam%sb(r))
!                Do iLayer = 1,MeshParam%KMax
!                    If (abs(HydroParam%hj(Face)-MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1))<HydroParam%PCRI) then
!                        HydroParam%hj(Face)=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)
!                    EndIf
!                EndDo
!                Do iLayer = 1,MeshParam%KMaxsub
!                    If (abs(HydroParam%sj(Face)-MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer+1))<HydroParam%PCRI) then
!                        HydroParam%sj(Face)=MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer+1)
!                    EndIf
!                EndDo
!            EndIf
!        EndDo   
!    EndDo
!    
!    
!    !Compute Index Smallm(j) and CapitalM(j)
!    ! Smallm and CapitalM are related to mj and Mj in [1]
!    Do iEdge = 1,MeshParam%nEdge
!
!        ! 1.1 Compute Water Depth
!        l = MeshParam%Left(iEdge)
!        r = MeshParam%Right(iEdge)
!
!        If (r == 0) Then
!            HydroParam%H(iEdge) = Max( HydroParam%PCRI, -HydroParam%sj(iEdge) + HydroParam%eta(l) )
!        Else
!            HydroParam%H(iEdge) = Max( HydroParam%PCRI,-HydroParam%sj(iEdge) + HydroParam%eta(l), -HydroParam%sj(iEdge) + HydroParam%eta(r) )!CAYO
!        EndIf
! 
!        ! Lower Index - Superficial Flow Layer
!        Do iLayer = 1,MeshParam%KMax
!            If (iLayer == MeshParam%KMax) Then
!                If (HydroParam%hj(iEdge)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
!                    HydroParam%Smallm(iEdge) = iLayer
!                    exit
!                EndIf
!            Else
!                If (HydroParam%hj(iEdge)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%hj(iEdge)<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
!                    HydroParam%Smallm(iEdge) = iLayer
!                    exit
!                EndIf
!            EndIf
!        EndDo
!
!        ! Upper Index - Superficial Flow Layer
!        Do iLayer = 1,MeshParam%KMax
!            If (iLayer == MeshParam%KMax) Then
!                If (HydroParam%H(iEdge) + HydroParam%hj(iEdge)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
!                    HydroParam%CapitalM(iEdge) = iLayer
!                    exit
!                EndIf
!            Else
!                If (HydroParam%H(iEdge) + HydroParam%hj(iEdge)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%H(iEdge) + HydroParam%hj(iEdge)<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
!                    HydroParam%CapitalM(iEdge) = iLayer
!                    exit
!                EndIf
!            EndIf
!        EndDo
!
!        HydroParam%Smallms(iEdge) = HydroParam%Smallm(iEdge)
!        HydroParam%CapitalMs(iEdge) = HydroParam%Smallm(iEdge)   
!        If (MeshParam%iBedrock == 1) Then
!            ! Lower Index - Subsuperficial Flow Layer  !CAYO
!            Do iLayer = 1,MeshParam%KMaxsub
!                If (iLayer == MeshParam%KMaxsub) Then
!                    If (HydroParam%sj(iEdge)>=MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer+1)) then
!                        HydroParam%Smallms(iEdge) = iLayer
!                        exit
!                    EndIf
!                Else
!                    If (HydroParam%sj(iEdge)>=MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer+1).and.HydroParam%sj(iEdge)<MeshParam%LIMCAMAUXsub(MeshParam%KMaxsub-iLayer)) then
!                        HydroParam%Smallms(iEdge) = iLayer
!                        exit
!                    EndIf
!                EndIf
!            EndDo
!
!            ! Upper Index - Subsuperficial Flow Layer
!            If ( HydroParam%Smallm(iEdge) > 1) Then 
!                HydroParam%CapitalMs(iEdge) = HydroParam%Smallm(iEdge) - 1
!                if (HydroParam%CapitalMs(iEdge)<HydroParam%Smallms(iEdge)) then
!                    HydroParam%CapitalMs(iEdge) = HydroParam%Smallms(iEdge)
!                endif
!            Else
!                HydroParam%CapitalMs(iEdge) = HydroParam%Smallms(iEdge)
!            EndIf
!        EndIf  
!        
!        !9.2 Compute Elevation in the Edges
!        Do iLayer = HydroParam%Smallms(iEdge) + 1, HydroParam%CapitalM(iEdge)
!            HydroParam%Z(iLayer,iEdge) = MeshParam%LIMCAMAUX(MeshParam%KMax+MeshParam%KMaxsub-iLayer) !MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1) !zL + (iLayer - 1)*dz                  ! Equidistant Core Grid
!        EndDo
!        
!        HydroParam%Z(HydroParam%Smallms(iEdge),iEdge)     = HydroParam%sj(iEdge)                 ! Bottom
!        If (r == 0) Then
!            If ( HydroParam%eta(l) - HydroParam%sj(iEdge) <= HydroParam%PCRI+NearZero) Then
!                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%sj(iEdge) + HydroParam%PCRI
!            Else
!                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%eta(l)
!            EndIf
!             !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
!        Else
!            If ( HydroParam%H(iEdge) <= HydroParam%PCRI+NearZero ) Then
!                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%eta(l),HydroParam%eta(r))
!            Else
!                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = (0.5d0)*(HydroParam%eta(l)+HydroParam%eta(r)) !Max(HydroParam%eta(l),HydroParam%eta(r)) !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
!            EndIf
!        EndIf        
!        
!        
!        ! 1.4 Compute the Vertical Mesh Spacing
!        Do iLayer = HydroParam%Smallms(iEdge), HydroParam%CapitalM(iEdge)
!            HydroParam%DZj(iLayer,iEdge) = HydroParam%Z(iLayer+1,iEdge) - HydroParam%Z(iLayer,iEdge)
!            HydroParam%DZj(iLayer,iEdge) = Max(HydroParam%Pcri,HydroParam%DZj(iLayer,iEdge))
!        EndDo        
!        
!        !1.5 Compute Kj and DZhj/DZsj
!		HydroParam%DZhj(:,iEdge) = HydroParam%DZj(:,iEdge)
!        If (MeshParam%iBedrock == 1) Then
!            Do iLayer = HydroParam%Smallms(iEdge), HydroParam%CapitalM(iEdge) ! 
!        
!                !If (iLayer >= HydroParam%Smallm(iEdge) .or. HydroParam%Smallms(iEdge) == HydroParam%Smallm(iEdge) ) Then
!                !    continue
!                !Else
!                !    MeshParam%Kj(iLayer,iEdge) = 0.01
!                !EndIf
!                
!                If (iLayer < HydroParam%Smallm(iEdge)) Then
!                    HydroParam%DZhj(iLayer,iEdge) = 0.
!                    HydroParam%DZsj(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)
!                ElseIf (iLayer > HydroParam%Smallm(iEdge)) Then
!                    HydroParam%DZhj(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)
!                    HydroParam%DZsj(iLayer,iEdge) = 0.
!                Else
!                    If ( HydroParam%Z(iLayer+1,iEdge) >= HydroParam%hj(iEdge) ) Then
!                        HydroParam%DZhj(iLayer,iEdge) = HydroParam%Z(iLayer+1,iEdge) - HydroParam%hj(iEdge)
!                        HydroParam%DZsj(iLayer,iEdge) = HydroParam%hj(iEdge) - HydroParam%Z(iLayer,iEdge)
!                    Else
!                        HydroParam%DZhj(iLayer,iEdge) = 0.
!                        HydroParam%DZsj(iLayer,iEdge) = HydroParam%Z(iLayer+1,iEdge) - HydroParam%Z(iLayer,iEdge)
!                    EndIf
!                EndIf    
!                
!                If (HydroParam%DZsj(iLayer,iEdge) > 0) Then
!                    MeshParam%Kj(iLayer,iEdge) = k0
!                EndIf
!                
!            EndDo
!		EndIf
!        
!    EndDo
!     
!    ! 4. Compute the New Vertical Velocity Field
!    ! w(Smallm-1,:) = 0.d0 -> No flux through the Bottom
!    ! w(CapitalM,:) = 0.d0 -> No flux through the Surface
!    ! The velocity is defined in the barycenter of the Top Face of each Cell
!    ! *Obs: The outward direction is from Bottom to Top.
!    Do iElem = 1, MeshParam%nElem
!        ! 4.1 Define the range between Bottom and Top Layers. (Tricky part - Not explained in [1])
!        ! In a real world application, the bathymetry might change between the faces of the prism.
!        ! We need to find "Max(Smallm(:))" and "Min(CapitalM(:))" in each face. 
!        ! Then we can calculate the Vertical Velocity only in those Layers
!        ! *Obs: The velocity in faces with Smallm < Max(Smallm(:)) == 0. The same goes to faces with CapitalM > Min(CapitalM(:))
!        ! 4.1.1 Find "mi" and "Mi", the Bottom and Top Layers range for the Element, not for the faces
!        ! Set mi
!        !If ( Smallm(Edge(1,iElem))==Smallm(Edge(2,iElem)).AND.Smallm(Edge(2,iElem))==Smallm(Edge(3,iElem)).AND.Smallm(Edge(3,iElem))==Smallm(Edge(4,iElem)) ) Then
!        !    ElSmallm(iElem) = Smallm(Edge(1,iElem))         ! All "mj" are equal. We can choose whatever we want!
!        !Else
!        !    ! If we have different Bottom Layers, we merge all in one cell!
!        !    !ElSmallm(iElem) = Max( Smallm(Edge(1,iElem)),Smallm(Edge(2,iElem)),Smallm(Edge(3,iElem)),Smallm(Edge(4,iElem)) ) !- 1
!        !    ElSmallm(iElem) = Min( Smallm(Edge(1,iElem)),Smallm(Edge(2,iElem)),Smallm(Edge(3,iElem)),Smallm(Edge(4,iElem)) )
!        !EndIf
!        !! Set Mi
!        !If ( CapitalM(Edge(1,iElem))==CapitalM(Edge(2,iElem)).AND.CapitalM(Edge(2,iElem))==CapitalM(Edge(3,iElem)).AND.CapitalM(Edge(3,iElem))==CapitalM(Edge(4,iElem)) ) Then
!        !    ElCapitalM(iElem) = CapitalM(Edge(1,iElem))         ! All "Mj" are equal. We can choose whatever we want!
!        !Else
!        !    ! If we have different Top Layers, we merge all in one cell!
!        !    !ElCapitalM(iElem) = Min( CapitalM(Edge(1,iElem)),CapitalM(Edge(2,iElem)),CapitalM(Edge(3,iElem)),CapitalM(Edge(4,iElem)) ) !+ 1
!        !    ElCapitalM(iElem) = Max( CapitalM(Edge(1,iElem)),CapitalM(Edge(2,iElem)),CapitalM(Edge(3,iElem)),CapitalM(Edge(4,iElem)) )
!        !EndIf
!        !If (HydroParam%hb(iElem)<-0.4) Then
!        !    HydroParam%hb(iElem) = -0.4
!        !EndIf
!        
!        ! Lower Index - Superficial Flow Layer       
!        Do iLayer = 1,MeshParam%KMax
!            If (iLayer == MeshParam%KMax) Then
!                If (HydroParam%hb(iElem)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
!                    HydroParam%ElSmallm(iElem) = iLayer
!                    exit
!                EndIf
!            Else
!                If (HydroParam%hb(iElem)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%hb(iElem)<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
!                    HydroParam%ElSmallm(iElem) = iLayer
!                    exit
!                EndIf
!            EndIf
!        EndDo
!        
!        ! Upper Index - Superficial Flow Layer
!        Do iLayer = 1,MeshParam%KMax
!            If (iLayer == MeshParam%KMax) Then
!                If (HydroParam%eta(iElem)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
!                    HydroParam%ElCapitalM(iElem) = iLayer
!                    exit
!                EndIf
!            Else
!                If (HydroParam%eta(iElem)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%eta(iElem)<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
!                    HydroParam%ElCapitalM(iElem) = iLayer
!                    exit
!                EndIf
!            EndIf
!        EndDo                
!        
!        ! 4.1.2.1 Compute the Vertical Mesh Spacing
!        !Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem) !CAYO
!        !    HydroParam%DZi(iLayer,iElem) = HydroParam%Ze(iLayer+1,iElem) - HydroParam%Ze(iLayer,iElem)
!        !    HydroParam%DZi(iLayer,iElem) = Max(HydroParam%Pcri,HydroParam%DZi(iLayer,iElem))
!        !EndDo
!        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
!            If (iLayer == HydroParam%ElSmallm(iElem)) Then
!                HydroParam%Zb(iLayer,iElem) = (0.5d0)*( HydroParam%Ze(iLayer,iElem) + HydroParam%Ze(iLayer+1,iElem) ) !(0.5d0)*( 0. + Ze(iLayer+1,iElem) )   (verificar com Rafael)             
!            Else
!                HydroParam%Zb(iLayer,iElem) = (0.5d0)*( HydroParam%Ze(iLayer,iElem) + HydroParam%Ze(iLayer+1,iElem) )
!            EndIf
!        EndDo       
!        
!        HydroParam%ElSmallms(iElem)     = HydroParam%ElSmallm(iElem)
!        HydroParam%ElCapitalMs(iElem)   = HydroParam%ElSmallm(iElem)
!		If (MeshParam%iBedrock == 1) Then
!        ! Lower Index - Subsuperficial Flow Layer       
!            Do iLayer = 1,MeshParam%KMax
!                If (iLayer == MeshParam%KMax) Then
!                    If (HydroParam%sb(iElem)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
!                        HydroParam%ElSmallms(iElem) = iLayer
!                        exit
!                    EndIf
!                Else
!                    If (HydroParam%sb(iElem)>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.HydroParam%sb(iElem)<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
!                        HydroParam%ElSmallms(iElem) = iLayer
!                        exit
!                    EndIf
!                EndIf
!            EndDo
!
!            ! Upper Index - Superficial Flow Layer        
!            If ( HydroParam%ElSmallm(iElem) > 1) Then
!                HydroParam%ElCapitalMs(iElem) = HydroParam%ElSmallm(iElem) - 1
!                if (HydroParam%ElCapitalMs(iElem)<HydroParam%ElSmallms(iElem)) then
!                    HydroParam%ElCapitalMs(iElem) = HydroParam%ElSmallms(iElem)
!                endif
!            Else
!                HydroParam%ElCapitalMs(iElem) = HydroParam%ElSmallms(iElem)
!            EndIf
!		EndIf
!        
!        ! 4.1.2 Update the Element Vertical Spacing
!        
!        Do iLayer = HydroParam%ElSmallms(iElem)+1, HydroParam%ElCapitalM(iElem) ! Cayo 
!            HydroParam%Ze(iLayer,iElem) = MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)      ! Equidistant Core Grid
!        EndDo
!        HydroParam%Ze(HydroParam%ElSmallms(iElem),iElem)     = HydroParam%sb(iElem)                    ! Bottom
!        If ( HydroParam%eta(iElem) - HydroParam%sb(iElem) <= HydroParam%PCRI+NearZero ) Then
!            HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%sb(iElem) + HydroParam%PCRI !- hb(iElem)       ! Free-Surface (verificar com Rafael)
!        Else
!            HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%eta(iElem) !- hb(iElem)       ! Free-Surface (verificar com Rafael)    
!        EndIf
!        Do iLayer = 1, HydroParam%ElSmallms(iElem) - 1
!            HydroParam%Ze(iLayer,iElem) = HydroParam%sb(iElem)
!        EndDo
!        
!        ! 4.1.2.1 Compute the Vertical Mesh Spacing
!        Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem)
!            HydroParam%DZi(iLayer,iElem) = HydroParam%Ze(iLayer+1,iElem) - HydroParam%Ze(iLayer,iElem)
!            HydroParam%DZi(iLayer,iElem) = Max(HydroParam%Pcri,HydroParam%DZi(iLayer,iElem))
!        EndDo
!        
!        !xx. Porosity and and DZhi/DZsi
!        HydroParam%DZhi(:,iElem) = HydroParam%DZi(:,iElem)
!		If (MeshParam%iBedrock == 1) Then
!            Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem)
!            
!                If (iLayer >= HydroParam%ElSmallm(iElem) .or.  HydroParam%ElSmallms(iElem) ==  HydroParam%ElSmallm(iElem)) Then
!                    continue
!                Else
!                    MeshParam%ei(iLayer,iElem) = e0
!                EndIf
!
!                If (iLayer < HydroParam%ElSmallm(iElem)) Then
!                    HydroParam%DZhi(iLayer,iElem) = 0.d0
!                    HydroParam%DZsi(iLayer,iElem) = HydroParam%DZi(iLayer,iElem)               
!                ElseIf (iLayer > HydroParam%ElSmallm(iElem)) Then
!                    HydroParam%DZhi(iLayer,iElem) = HydroParam%DZi(iLayer,iElem)
!                    HydroParam%DZsi(iLayer,iElem) = 0.d0
!                Else
!                    If (HydroParam%Ze(iLayer+1,iElem) >= HydroParam%hb(iElem) ) Then
!                        HydroParam%DZhi(iLayer,iElem) = HydroParam%Ze(iLayer+1,iElem) - HydroParam%hb(iElem)
!                        HydroParam%DZsi(iLayer,iElem) = HydroParam%hb(iElem) - HydroParam%Ze(iLayer,iElem)
!                    Else
!                        HydroParam%DZhi(iLayer,iElem) = 0.d0
!                        HydroParam%DZsi(iLayer,iElem) = HydroParam%Ze(iLayer+1,iElem) - HydroParam%Ze(iLayer,iElem)
!                    EndIf
!                EndIf                                    
!                
!                If(HydroParam%DZsi(iLayer,iElem)>0) Then
!                    MeshParam%ei(iLayer,iElem) = e0
!                EndIf
!                                
!            EndDo
!		EndIf
!        
!    EndDo     
!    
!    ! Set Initial Conditions of velocity components 
!	Do iElem = 1, MeshParam%nElem
!        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
!            Do iEdge = 1,4
!                r = MeshParam%Right(MeshParam%Edge(iEdge,iElem))
!                !uIniVec(iLayer,1) = -6*(EdgeBary(1,Edge(iEdge,iElem))/100.)**2. + 6*EdgeBary(1,Edge(iEdge,iElem))/100.
!                !uIniVec(iLayer,2) = 0.d0
!                If (r==0.or.HydroParam%H(MeshParam%Edge(iEdge,iElem))<=HydroParam%Pcri+NearZero) Then ! No slip condition
!                    HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem)) = 0.d0
!                Else
!                    If (simParam%it > 0) Then
!                        HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem)) = simParam%usave(iLayer,MeshParam%Edge(iEdge,iElem))
!                        !HydroParam%utang(iLayer,MeshParam%Edge(iEdge,iElem)) = simParam%utangsave(iLayer,MeshParam%Edge(iEdge,iElem))
!                    Else
!                        HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem))  = Dot_Product(HydroParam%uIniVec(iLayer,1:2),MeshParam%NormalVector(:,iEdge,iElem))
!                        HydroParam%utang(iLayer,MeshParam%Edge(iEdge,iElem))  = Dot_Product(HydroParam%uIniVec(iLayer,1:2),MeshParam%TangentVector(:,iEdge,iElem))
!                    EndIf
!                    
!                EndIf
!            EndDo
!        EndDo
!    EndDo
!    If (simParam%it > 0) Then
!        HydroParam%w = simParam%wsave
!        Call VelocitiesSub(HydroParam,MeshParam)
!    Else
!        HydroParam%w = HydroParam%Wini
!    EndIf
!    
!    If (HydroParam%iVTurb ==3 ) Then
!        If(.NOT.Allocated(HydroParam%TKE))             Allocate(HydroParam%TKE(MeshParam%Kmax+1,MeshParam%nElem))
!        If(.NOT.Allocated(HydroParam%LengthScale))     Allocate(HydroParam%LengthScale(MeshParam%Kmax+1,MeshParam%nElem))
!        If(.NOT.Allocated(HydroParam%TKEP))            Allocate(HydroParam%TKEP(MeshParam%Kmax+1,MeshParam%nElem))
!        If(.NOT.Allocated(HydroParam%LengthScaleP))    Allocate(HydroParam%LengthScaleP(MeshParam%Kmax+1,MeshParam%nElem))
!        If(.NOT.Allocated(HydroParam%DissipRate))      Allocate(HydroParam%DissipRate(MeshParam%Kmax+1,MeshParam%nElem))
!        If(.NOT.Allocated(HydroParam%ShearProd))       Allocate(HydroParam%ShearProd(MeshParam%Kmax+1,MeshParam%nElem))
!        If(.NOT.Allocated(HydroParam%BuoyancyProd))    Allocate(HydroParam%BuoyancyProd(MeshParam%Kmax+1,MeshParam%nElem))
!        
!    EndIf
!    
!    
!    !Open(1001, File ='VelocityEdge.txt', Action='Write', Status='Replace')
!    !Do iEdge = 1, nEdge
!    !    AuxVel = -6*(EdgeBary(1,iEdge)/100.)**2. + 6*EdgeBary(1,iEdge)/100.
!    !    Write(1001,'(i8.10,4F20.10)') iEdge, uxy(1,1,iEdge), AuxVel, uxy(1,2,iEdge), 0.
!    !EndDo
!    !Close(1001)
!    
!    
!    
!End Subroutine ReadHydroIniCond
!