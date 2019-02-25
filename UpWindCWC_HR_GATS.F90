Subroutine UpWindCWC_HR_GATS(index,uLoadVar,dVar,VarEst,VarEstP,dt,dtday,HydroParam,MeshParam,psi_flag)

    ! This routine solves the Transport Equation usind a centered finite difference scheme for variables with Boundary Condition allowed.
    ! Called in routines: ABIOTICO, SALINITY

    ! Example:
    ! sDSal(I,J,K)       == VarEst
    ! dVarSal            == dVar
    ! sDSalP(I,J,K)      == VarEstP
    ! sDSalP(I-1,J,K)    == VarEstPIM
    ! sDSalP(I+1,J,K)    == VarEstPIP	
    ! sDSalP(I,J-1,K)    == VarEstPJM
    ! sDSalP(I,J+1,K)    == VarEstPJP
    ! sDSalP(I,J,K-1)    == VarEstPKM
    ! sDSalP(I,J,K+1)    == VarEstPKP
    ! INDEXCONT(ICEL,14) == index => 14
    ! uLoadVar        == uDLoadSal
    ! sDSal0             == VarEstINI
    
    ! Obs: A dry cell receive the value of -99. The code recognizes this value and the concetration is not showed in the display (Intel - Array Visualizer)
    Use MeshVars
    Use Hydrodynamic
    
    Implicit None
    Real, Dimension(4):: ConcFace
	Real:: CONCU,CONCD,DDCFace,DCDFace,Dispersion,Advection
    Real:: DCDFace_psi,Dispersion_psi,fi_small,r_face,fi_big,psi_face
	Real:: DDCZ,DC2DZ2,DCDZU,DCDZD,DCDZU_psi,DCDZD_psi,DC2DZ2_psi
    Real:: V,SumQj,Sumdj,dtau
	Integer:: index,iElem,iEdge,Face,iLayer,Sig,l,r,nSubCycle,iSubCycle,psi_flag,iEdgein, rr
    Real:: NearZero = 1e-10
    Real:: dt,dtday,Courant,qj
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Real, Dimension(MeshParam%KMax,MeshParam%nElem):: dVar,VarEst,VarEstP
    Real, Dimension(MeshParam%KMax,MeshParam%nElem):: uLoadVar
    Real, Dimension(MeshParam%KMax):: r_num_plus,r_den_plus,r_num_minus,r_den_minus,r_num,r_den
    !Boundary tratment variables
    Real :: a,b,c,ad,bd,cd,num,den,numconc,denconc,numdist,dendist, Neighbor
    
    HydroParam%psi_edge = 0.d0
    HydroParam%psi_cell = 0.d0
    Do iElem = 1, MeshParam%nElem
        
        
        ! 1. Define Sub-Cycle Time Step    
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            SumQj = 0.
            Sumdj = 0.
            Do iEdge = 1,4
                Face = MeshParam%Edge(iEdge,iElem)
                l = MeshParam%Left(Face)
                r = MeshParam%Right(Face)
                
                If ( MeshParam%Neighbor(iEdge,iElem)==0.or.HydroParam%H(MeshParam%Edge(iEdge,iElem)) <= HydroParam%PCRI+NearZero) Then 
		            Sumdj = Sumdj 
                    SumQj = SumQj
                Else
                    If (iLayer<HydroParam%ElSmallm(MeshParam%Neighbor(iEdge,iElem))) Then
                        DDCFace = 0.
		                Sumdj = Sumdj 
                        SumQj = SumQj
                    Else
                        If (Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)*HydroParam%DZjt(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)*HydroParam%DZjt(iLayer,Face))>=0) Then ! Water leaving
                            SumQj = SumQj + abs(MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))
                        Else
                            SumQj = SumQj
                        EndIf
                        DDCFace = (0.5*(HydroParam%HorDiffusivity(2,iLayer,MeshParam%Neighbor(iEdge,iElem))+HydroParam%HorDiffusivity(2,iLayer,iElem)))
		                Sumdj = Sumdj + max(0.,abs(MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face)*HydroParam%DZjt(iLayer,Face)*DDCFace - 0.5*abs(MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))))
                                
                    EndIf
                EndIf
            EndDo
            HydroParam%Locdt(iLayer,iElem) = MeshParam%Area(iElem)*HydroParam%DZit(iLayer,iElem)/(2*SumQj+2*abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem)))+Sumdj+NearZero)
        EndDo       ! Layer Loop

        
        !! 2. Define Psi Values
        !r_num_plus = 0.
        !r_den_plus = 0.
        !
        !Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
        !    
        !    Do iEdge = 1,4
        !        Face = MeshParam%Edge(iEdge,iElem)
        !        l = MeshParam%Left(Face)
        !        r = MeshParam%Right(Face)
        !        If ( r==0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero ) Then 
        !            !r_num_plus(iLayer) = r_num_plus(iLayer)
        !            !r_den_plus(iLayer) = r_den_plus(iLayer)
        !        Else
        !            If (Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))<=0) Then !inflow faces
        !                If (iLayer<HydroParam%ElSmallm(MeshParam%Neighbor(iEdge,iElem))) Then 
        !                    r_num_plus(iLayer) = r_num_plus(iLayer)
        !                    r_den_plus(iLayer) = r_den_plus(iLayer)
        !                Else
        !                    r_num_plus(iLayer) = r_num_plus(iLayer) + abs((HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*HydroParam%DZjt(iLayer,Face)*MeshParam%EdgeLength(Face))*abs(VarEstP(iLayer,iElem)-VarEstP(iLayer,MeshParam%Neighbor(iEdge,iElem)))
        !                    r_den_plus(iLayer) = r_den_plus(iLayer) + abs((HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*HydroParam%DZjt(iLayer,Face)*MeshParam%EdgeLength(Face))
        !                EndIf
        !            EndIf
        !        EndIf
        !    EndDo
        !        
        !    If (HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem)>0) Then !Bottom layer inflow
        !        If (iLayer>HydroParam%ElSmallm(iElem)) Then
        !            r_num_plus(iLayer) = r_num_plus(iLayer) + abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem)))*abs(VarEstP(iLayer,iElem)-VarEstP(iLayer-1,iElem))
        !            r_den_plus(iLayer) = r_den_plus(iLayer) + abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem)))
        !        EndIf
        !    EndIf
        !        
        !    If (HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem)<0) Then !top layer inflow
        !        If (iLayer< HydroParam%ElCapitalM(iElem)) Then
        !            r_num_plus(iLayer) = r_num_plus(iLayer) + abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem)))*abs(VarEstP(iLayer,iElem)-VarEstP(iLayer+1,iElem))
        !            r_den_plus(iLayer) = r_den_plus(iLayer) + abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem)))
        !        EndIf
        !    EndIf
        !        
        !    Do iEdge = 1,4
        !        Face = MeshParam%Edge(iEdge,iElem)
        !        l = MeshParam%Left(Face)
        !        r = MeshParam%Right(Face)
        !        If ( r==0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero ) Then 
        !            HydroParam%psi_edge(iLayer,Face) = 0.
        !        Else
        !            If (Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))>0) Then !outflow faces
        !                fi_small = Min(1., 2.*0.5*(HydroParam%HorDiffusivity(2,iLayer,MeshParam%Neighbor(iEdge,iElem))+HydroParam%HorDiffusivity(2,iLayer,iElem))*HydroParam%DZjt(iLayer,Face)/(abs((HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face))+NearZero))
        !                r_face = (r_num_plus(iLayer)/(abs(VarEstP(iLayer,MeshParam%Neighbor(iEdge,iElem))-VarEstP(iLayer,iElem))*(r_den_plus(iLayer)) + NearZero))
        !                Courant = 0.5*(HydroParam%Theta*HydroParam%u(iLayer,Face)*+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*dt/MeshParam%CirDistance(Face)
        !                Call Psi_value(psi_flag,r_face,Courant,fi_small,psi_face)
        !                HydroParam%psi_edge(iLayer,Face) = psi_face
        !            EndIf
        !        EndIf
        !    EndDo
        !        
        !EndDo
        !
        !Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)    
        !
        !    If (HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem)<0) Then !Bottom layer outflow
        !            If (iLayer>HydroParam%ElSmallm(iElem)) Then
        !                DDCZ = HydroParam%VerEddyDiffCell(iLayer,iElem)
        !                fi_small = Min(1., 2.*MeshParam%Area(iElem)*DDCZ/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))/2.)/abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))+NearZero))
        !                r_face = (r_num_plus(iLayer)/(abs(VarEstP(iLayer-1,iElem)-VarEstP(iLayer,iElem))*(r_den_plus(iLayer)) + NearZero))
        !                Courant = 0.5*(HydroParam%Theta*HydroParam%u(iLayer,Face)*+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*dt/MeshParam%CirDistance(Face)
        !                Call Psi_value(psi_flag,r_face,Courant,fi_small,psi_face)
        !                HydroParam%psi_cell(iLayer,iElem) = psi_face
        !            EndIf
        !        EndIf
        !        
        !        If (HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem)>0) Then !top layer outflow
        !            If (iLayer< HydroParam%ElCapitalM(iElem)) Then
        !                DDCZ = HydroParam%VerEddyDiffCell(iLayer+1,iElem)
        !                fi_small = Min(1., 2.*MeshParam%Area(iElem)*DDCZ/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))/2.)/abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))+NearZero))
        !                r_face = (r_num_plus(iLayer)/(abs(VarEstP(iLayer+1,iElem)-VarEstP(iLayer,iElem))*(r_den_plus(iLayer)) + NearZero))
        !                Courant = 0.5*(HydroParam%Theta*HydroParam%u(iLayer,Face)*+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*dt/MeshParam%CirDistance(Face)
        !                Call Psi_value(psi_flag,r_face,Courant,fi_small,psi_face)
        !                HydroParam%psi_cell(iLayer+1,iElem) = psi_face
        !            EndIf
        !        EndIf
        !
        !EndDo ! Layer Loop
        
    EndDo       ! Element Loop
    dtau = Min( dt, MinVal( HydroParam%Locdt ) )
    nSubCycle     = Ceiling(dt/dtau)
    If (nSubCycle>1) Then
        print*,nSubCycle
    EndIf
    
    HydroParam%DZitau = HydroParam%DZit
   ! 2. Start Global SubCycling Time Stepping
    
    Do iSubCycle = 1, nSubCycle
        
        HydroParam%DZistau = HydroParam%DZit+(HydroParam%DZi-HydroParam%DZit)*iSubCycle/nSubCycle
        
        Do iElem = 1,MeshParam%nElem
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                !if (HydroParam%ElSmallm(iElem) == HydroParam%ElCapitalM(iElem)) then
                !    !HydroParam%psi_edge(iLayer,Face) = 0.d0
                !    !HydroParam%psi_cell(iLayer,iElem) = 0.d0
                !    cycle
                !endif

                !Define Psi Values
                Do iEdge = 1,4
                    Face = MeshParam%Edge(iEdge,iElem)
                    l = MeshParam%Left(Face)
                    r = MeshParam%Right(Face)
                    If (iEdge==1) Then
                        iEdgein = 3
                    ElseIf (iEdge==2) Then
                        iEdgein = 4
                    ElseIf (iEdge==3) Then
                        iEdgein = 1
                    ElseIf (iEdge==4) Then
                        iEdgein = 2
                    EndIf
                    rr = MeshParam%Neighbor(iEdgein,iElem)
                    
                    !set boundary conditions for psi_edge
                    b= VarEstP(iLayer,iElem)
                    if (rr==0) then
                        a = b! VarEstP(iLayer,iElem)
                    else
                        if (iLayer<HydroParam%ElSmallm(rr)) then ! r==0 to the bottom
                            a = b !VarEstP(iLayer,iElem)
                        else
                            a= VarEstP(iLayer,MeshParam%Neighbor(iEdgein,iElem))
                        endif
                    endif
                    if (MeshParam%Neighbor(iEdge,iElem)==0) then
                        c = b !VarEstP(iLayer,iElem)
                        Neighbor = iElem
                    else
                        if (iLayer<HydroParam%ElSmallm(MeshParam%Neighbor(iEdge,iElem))) then ! r==0 to the bottom
                            c = b !VarEstP(iLayer,iElem)
                            Neighbor = iElem
                        else
                            c= VarEstP(iLayer,MeshParam%Neighbor(iEdge,iElem))
                            Neighbor = MeshParam%Neighbor(iEdge,iElem)
                        endif
                    endif
                    !r_face
                    num = (b-a)
                    den= (c-b)+NearZero
                    
                    if (abs(c-b)<=NearZero.or.abs(b-a)<=NearZero) then
                        HydroParam%psi_edge(iLayer,Face) = 0
                        cycle
                    endif
                    
                    !if (HydroParam%IndexInflowEdge(Face) > 0) then
                    !    HydroParam%psi_edge(iLayer,Face) = 0
                    !    cycle
                    !endif
                    
                    If ((HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))>0) Then !outflow faces 
                        fi_small = Min(1., 2.*0.5*(HydroParam%HorDiffusivity(2,iLayer,Neighbor)+HydroParam%HorDiffusivity(2,iLayer,iElem))*HydroParam%DZjt(iLayer,Face)/(abs((HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face))+NearZero))
                        r_face = num/den 
                        Courant = 0.5*(HydroParam%Theta*HydroParam%u(iLayer,Face)*+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*dt/MeshParam%CirDistance(Face)
                        Call Psi_value(psi_flag,r_face,Courant,fi_small,psi_face)
                        HydroParam%psi_edge(iLayer,Face) = psi_face
                        
                    EndIf
                    
                EndDo
                
                If (HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem)<0) Then !Bottom layer outflow 
                    !set boundary condition for psi_cel bottom layer outflow
                    b = VarEstP(iLayer,iElem)
                    bd = HydroParam%DZit(iLayer,iElem)
                    if (iLayer == HydroParam%ElCapitalM(iElem)) then
                        a = b
                        ad = bd
                    else
                        a = VarEstP(iLayer+1,iElem)
                        ad = HydroParam%DZit(iLayer+1,iElem)
                    endif
                    if (iLayer == HydroParam%ElSmallm(iElem)) then
                        c = b
                        cd = bd
                    else
                        c = VarEstP(iLayer-1,iElem)
                        cd = HydroParam%DZit(iLayer-1,iElem)
                    endif
                    numconc = b - a !(VarEstP(iLayer,iElem)-VarEstP(iLayer+1,iElem)
                    denconc = (c-b) !(VarEstP(iLayer-1,iElem)-VarEstP(iLayer,iElem)+NearZero))
                    numdist = 0.5d0*(bd+cd)
                    dendist = (0.5d0*(bd+ad))
                    
                    if (abs(b-a)<=NearZero.or.abs(c-b)<=NearZero) then
                        HydroParam%psi_cell(iLayer,iElem) = 0
                        cycle
                    endif
                    if (abs(numdist)<=NearZero.or.abs(dendist)<=NearZero) then
                        HydroParam%psi_cell(iLayer,iElem) = 0
                        cycle
                    endif
                    
                    DDCZ = HydroParam%VerEddyDiffCell(iLayer,iElem)
                    fi_small = Min(1., 2.*MeshParam%Area(iElem)*DDCZ/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))/2.)/abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))+NearZero))
                    r_face =  (numconc/denconc)*(numdist/dendist) !((VarEstP(iLayer,iElem)-VarEstP(iLayer+1,iElem))/(VarEstP(iLayer-1,iElem)-VarEstP(iLayer,iElem)+NearZero))*(0.5d0*(HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))/(0.5d0*(HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer+1,iElem))+NearZero))
                    Courant = 0.5*(HydroParam%Theta*HydroParam%u(iLayer,Face)*+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*dt/MeshParam%CirDistance(Face)
                    Call Psi_value(psi_flag,r_face,Courant,fi_small,psi_face)
                    HydroParam%psi_cell(iLayer,iElem) = psi_face
                    
                EndIf
                
                If (HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem)>0) Then !top layer outflow 
                    !if (iLayer < HydroParam%ElCapitalM(iElem)) then
                        !set boundary condition for psi_cel top layer outflow
                        b = VarEstP(iLayer,iElem)
                        bd = HydroParam%DZit(iLayer,iElem)
                        if (iLayer == HydroParam%ElCapitalM(iElem)) then
                            a = b
                            ad = HydroParam%DZit(iLayer-1,iElem) !equal to cd
                        else
                            a = VarEstP(iLayer+1,iElem)
                            ad = HydroParam%DZit(iLayer+1,iElem)
                        endif
                        if (iLayer == HydroParam%ElSmallm(iElem)) then
                            c = b
                            cd = bd
                        else
                            c = VarEstP(iLayer-1,iElem)
                            cd = HydroParam%DZit(iLayer-1,iElem)
                        endif
                        numconc = b - c !((VarEstP(iLayer,iElem)-VarEstP(iLayer-1,iElem))
                        denconc = (a-b) !(VarEstP(iLayer+1,iElem)-VarEstP(iLayer,iElem)+NearZero))
                        numdist = 0.5d0*(bd+ad)  !(0.5d0*(HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer+1,iElem))
                        dendist = (0.5d0*(bd+cd))  !(0.5d0*(HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))+NearZero))  
                    
                        if (abs(a-b)<=NearZero.or.abs(b-c)<=NearZero) then
                            HydroParam%psi_cell(iLayer,iElem) = 0
                            cycle
                        endif
                        if (abs(numdist)<=NearZero.or.abs(dendist)<=NearZero) then
                            HydroParam%psi_cell(iLayer,iElem) = 0
                            cycle
                        endif
                    
                        DDCZ = HydroParam%VerEddyDiffCell(iLayer+1,iElem)
                        fi_small = Min(1., 2.*MeshParam%Area(iElem)*DDCZ/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))/2.)/abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))+NearZero))
                        r_face =  ((VarEstP(iLayer,iElem)-VarEstP(iLayer-1,iElem))/(VarEstP(iLayer+1,iElem)-VarEstP(iLayer,iElem)+NearZero))*(0.5d0*(HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer+1,iElem))/(0.5d0*(HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))+NearZero)) !r_face = (r_num_plus(iLayer)/(abs(VarEstP(iLayer+1,iElem)-VarEstP(iLayer,iElem))*(r_den_plus(iLayer)) + NearZero))
                        Courant = 0.5*(HydroParam%Theta*HydroParam%u(iLayer,Face)*+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*dt/MeshParam%CirDistance(Face)
                        Call Psi_value(psi_flag,r_face,Courant,fi_small,psi_face)
                        HydroParam%psi_cell(iLayer,iElem) = psi_face
                    
                        !endif
                    EndIf  
                
                !Do iEdge = 1,4
                !    Face = MeshParam%Edge(iEdge,iElem)
                !    if (HydroParam%psi_edge(iLayer,Face)==0) then
                !        
                !        If (HydroParam%w(iLayer,iElem)<0) Then !Bottom layer outflow
                !            if (HydroParam%IndexInflowEdge(Face) > 0) then
                !                HydroParam%psi_cell(iLayer,iElem) = 0.d0
                !            endif
                !        EndIf
                !        
                !        If (HydroParam%w(iLayer+1,iElem)>0) Then !top layer outflow
                !            if (HydroParam%IndexInflowEdge(Face) > 0) then
                !                HydroParam%psi_cell(iLayer+1,iElem) = 0.d0
                !            endif
                !        EndIf
                !        
                !    endif
                !enddo                    
                 
            EndDo
        EndDo
        
        
    
        Do iElem = 1,MeshParam%nElem
        
            If ( HydroParam%eta(iElem)-HydroParam%hb(iElem) <= HydroParam%PCRI+NearZero.or.HydroParam%etan(iElem)-HydroParam%hb(iElem)<= HydroParam%PCRI+NearZero) Then    
                VarEst(HydroParam%ElCapitalM(iElem),iElem) = VarEstP(HydroParam%ElCapitalM(iElem),iElem)  
                Cycle
            EndIf

                Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                    
                    ! Horizontal Eddy-Diffusivity (Dispersion coefficient)
                    ! Concentration in each face
                    DCDFace = 0.
                    Do iEdge = 1,4
                        Face = MeshParam%Edge(iEdge,iElem)
                        l = MeshParam%Left(Face)
                        r = MeshParam%Right(Face)
                        If ( MeshParam%Neighbor(iEdge,iElem)==0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero ) Then 
		                    ConcFace(iEdge) = VarEstP(iLayer,iElem)
                            DDCFace = 0.
		                    DCDFace = DCDFace 
                        Else
                            If (iLayer<HydroParam%ElSmallm(MeshParam%Neighbor(iEdge,iElem))) Then
                                ConcFace(iEdge) = VarEstP(iLayer,iElem)
                                DDCFace = 0.
		                        DCDFace = DCDFace
                            else
                                If (Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)*HydroParam%DZjt(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)*HydroParam%DZjt(iLayer,Face))>=0) Then !outflow faces
                                    ConcFace(iEdge) = VarEstP(iLayer,iElem)
                                Else !inflow faces
                                    ConcFace(iEdge) = VarEstP(iLayer,MeshParam%Neighbor(iEdge,iElem))
                                endif
                            EndIf
		                    DCDFace = DCDFace + (MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face))*HydroParam%DZjt(iLayer,Face)*DDCFace*(ConcFace(iEdge)-VarEstP(iLayer,iElem)) 
                        EndIf
                    EndDo
                    
                    Dispersion = DCDFace/MeshParam%Area(iElem)
                    
                    
                    !flux limiter to reduce Horizontal numerical diffusion
                    DCDFace_psi = 0.
                    Do iEdge = 1,4
                        Face = MeshParam%Edge(iEdge,iElem)
                        l = MeshParam%Left(Face)
                        r = MeshParam%Right(Face)
                        If ( MeshParam%Neighbor(iEdge,iElem)==0.or.HydroParam%H(Face) <= HydroParam%PCRI+NearZero ) Then 
		                    DCDFace_psi = DCDFace_psi
                        Else
                            If (iLayer<HydroParam%ElSmallm(MeshParam%Neighbor(iEdge,iElem))) Then
                                DCDFace_psi = DCDFace_psi
                            else
                                DCDFace_psi = DCDFace_psi + HydroParam%psi_edge(iLayer,Face)*abs((HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*HydroParam%DZjt(iLayer,Face)*MeshParam%EdgeLength(Face))*(ConcFace(iEdge)-VarEstP(iLayer,iElem))                        
                            endif
                        EndIf
                    EndDo
                    Dispersion_psi = DCDFace_psi/MeshParam%Area(iElem)
                    
                
                    If ( HydroParam%ElSmallm(iElem) == HydroParam%ElCapitalM(iElem) ) Then        ! Only One Vertical Layer
                        
                        If (abs(-Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*(HydroParam%Theta*HydroParam%u(iLayer,MeshParam%Edge(1,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(1,iElem))+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,MeshParam%Edge(1,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(1,iElem)))*MeshParam%dx*DT-Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*(HydroParam%Theta*HydroParam%u(iLayer,MeshParam%Edge(2,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(2,iElem))+ (1.-HydroParam%Theta)*HydroParam%ut(iLayer,MeshParam%Edge(2,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(2,iElem)))*MeshParam%dy*dt-Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*(HydroParam%Theta*HydroParam%u(iLayer,MeshParam%Edge(3,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(3,iElem))+ (1.-HydroParam%Theta)*HydroParam%ut(iLayer,MeshParam%Edge(3,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(3,iElem)))*MeshParam%dy*dt-Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*(HydroParam%Theta*HydroParam%u(iLayer,MeshParam%Edge(4,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(4,iElem))+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,MeshParam%Edge(4,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(4,iElem)))*MeshParam%dy*dt- (HydroParam%eta(iElem)-HydroParam%etaplus(ielem)/nSubCycle-HydroParam%etan(iElem))*MeshParam%dx*MeshParam%dy)>0.1) Then
                            VarEst(HydroParam%ElCapitalM(iElem),iElem) = VarEstP(HydroParam%ElCapitalM(iElem),iElem)  
                            Cycle
                        EndIf
                        
                        ! Top
                        CONCU = VarEstP(iLayer,iElem)
                        ! Bottom
		                CONCD = VarEstP(iLayer,iElem)
                    
                        DC2DZ2 = 0.
        
                        ! Turbulence Model (pegar referência)  
                        DDCZ = 0.
        
                        ! Transport Equation Solution for 2D
                        If (index== 0) Then !No boundary condition
                            Advection = 0.
                            Do iEdge = 1,4
                                Face = MeshParam%Edge(iEdge,iElem)
                                l = MeshParam%Left(Face)
                                r = MeshParam%Right(Face)
                                Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                            EndDo
                                VarEst(iLayer,iElem) = MAX(NearZero,(dt/nSubCycle*(Dispersion-Advection)-dt/nSubCycle/2.*Dispersion_psi+dtday/nSubCycle*dVar(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem))/(HydroParam%DZistau(iLayer,iElem)-HydroParam%etaplus(iElem)/nSubCycle)) 
                        Else
                            Advection = 0.
                            num=0
                            Do iEdge = 1,4
                                Face = MeshParam%Edge(iEdge,iElem)
                                l = MeshParam%Left(Face)
                                r = MeshParam%Right(Face)
                                If (HydroParam%IndexInflowEdge(Face) > 0.or.HydroParam%IndexWaterLevelEdge(Face) > 0) Then !Inflow/Outflow boundary condition
                                    If (HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)>=0) Then ! Outflow Loading
                                        Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)    
                                    Else !InFlow Loading
                                        Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*uLoadVar(iLayer,iElem)/MeshParam%Area(iElem) 
                                    EndIf
                                    num=1
                                Else
                                    Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                                EndIf
                            EndDo
                            if (num==1) then
                                VarEst(iLayer,iElem) = MAX(NearZero,(dt/nSubCycle*(Dispersion-Advection)+dtday/nSubCycle*dVar(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem))/(HydroParam%DZistau(iLayer,iElem)-HydroParam%etaplus(iElem)/nSubCycle)) 
                            else
                                VarEst(iLayer,iElem) = MAX(NearZero,(dt/nSubCycle*(Dispersion-Advection)-dt/nSubCycle/2.*Dispersion_psi+dtday/nSubCycle*dVar(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem))/(HydroParam%DZistau(iLayer,iElem)-HydroParam%etaplus(iElem)/nSubCycle)) 
                            endif
                        EndIf
                
                    Else
                        If ( iLayer == HydroParam%ElSmallm(iElem) ) Then !Bottom layer

                            ! Concentration in the Top layer
                            If ((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))>=0) Then
                                CONCU = VarEstP(iLayer,iElem)                         ! Top face concentration is equal to the concentration in the layer    
                            Else
                                CONCU = VarEstP(iLayer+1,iElem)                         ! Top face concentration is equal to the concentration in the layer
                            EndIF
                            ! Concentration in the Top Bottom layer
		                    CONCD = VarEstP(iLayer,iElem) 

		                    ! Vertical Diffusivity Flux in the top layer
                            DDCZ = HydroParam%VerEddyDiffCell(iLayer+1,iElem)

                            !DCDZU = DDCZ*( VarEstP(iLayer+1,iElem) - VarEstP(iLayer,iElem) )/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer+1,iElem))/2.)
                            DCDZU = (VarEstP(iLayer+1,iElem) - VarEstP(iLayer,iElem))*max(0., (MeshParam%Area(iElem)*DDCZ/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer+1,iElem))/2.))-((1./2.)*abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem)))))
                        
                            !DCDZU = Max(0., DCDZU - 0.5*abs(MeshParam%Area(iElem)*HydroParam%w(iLayer+1,iElem)))
                        
                            !flux limiter to reduce vertical numerical diffusion in the top layer
                            DCDZU_psi = HydroParam%psi_cell(iLayer+1,iElem)*abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem)))*(VarEstP(iLayer+1,iElem) - VarEstP(iLayer,iElem))
                        
                            ! Vertical Diffusivity Flux in the bottom layer
                            DCDZD = 0.
                            DCDZD_psi = 0.
                        
                            !DCDZD = Max(0., DCDZD - 0.5*abs(Area(iElem)*w(iLayer,iElem)))
                            ! Vertical gradient of Diffusivity
                            DC2DZ2 = (DCDZU - DCDZD)/MeshParam%Area(iElem) 
                        
                            DC2DZ2_psi = (DCDZU_psi - DCDZD_psi)/MeshParam%Area(iElem)
        
        
                            ! Transport Equation Solution
                            If (index== 0) Then !No boundary condition
                                Advection = 0.
                                Do iEdge = 1,4
                                    Face = MeshParam%Edge(iEdge,iElem)
                                    l = MeshParam%Left(Face)
                                    r = MeshParam%Right(Face)
                                    Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                                EndDo
                                    VarEst(iLayer,iElem) = MAX(NearZero,(dt/nSubCycle*(Dispersion-Advection+DC2DZ2)-dt/nSubCycle/2.*Dispersion_psi-dt/nSubCycle/2.*DC2DZ2_psi+dtday/nSubCycle*dVar(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)-((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*dt/nSubCycle)/HydroParam%DZitau(iLayer,iElem))
                            Else
                                Advection = 0.
                                num=0
                                Do iEdge = 1,4
                                    Face = MeshParam%Edge(iEdge,iElem)
                                    l = MeshParam%Left(Face)
                                    r = MeshParam%Right(Face)
                                    If (HydroParam%IndexInflowEdge(Face) > 0.or.HydroParam%IndexWaterLevelEdge(Face) > 0) Then !Inflow/Outflow boundary condition
                                        If (HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)>=0) Then ! Outflow Loading
                                            Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)    
                                        Else !InFlow Loading
                                            Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*uLoadVar(iLayer,iElem)/MeshParam%Area(iElem) 
                                        EndIf
                                        num=1
                                    Else
                                        Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                                    EndIf
                                EndDo
                                if (num==1) then
                                    VarEst(iLayer,iElem) = MAX(NearZero,(dt/nSubCycle*(Dispersion-Advection+DC2DZ2)+dtday/nSubCycle*dVar(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)-((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*dt/nSubCycle)/HydroParam%DZitau(iLayer,iElem))
                                else
                                    VarEst(iLayer,iElem) = MAX(NearZero,(dt/nSubCycle*(Dispersion-Advection+DC2DZ2)-dt/nSubCycle/2.*Dispersion_psi-dt/nSubCycle/2.*DC2DZ2_psi+dtday/nSubCycle*dVar(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)-((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*dt/nSubCycle)/HydroParam%DZitau(iLayer,iElem))
                                endif
                            EndIf
                    
                        Elseif ( iLayer == HydroParam%ElCapitalM(iElem) ) Then !Surface layer  
                        
                            ! Concentration in the Top layer
                            CONCU = VarEstP(iLayer,iElem)                      ! Up
                            ! Concentration in the bottom layer
                            If ((HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))>=0) Then
		                        CONCD = VarEstP(iLayer-1,iElem) 
                            Else
                                CONCD = VarEstP(iLayer,iElem)
                            EndIf

		                    ! Vertical Diffusivity Flux in the top layer
                            DCDZU = 0.
                        
                            DCDZU_psi = 0.
                        
                            ! Vertical Diffusivity Flux in the bottom layer
                            DDCZ = HydroParam%VerEddyDiffCell(iLayer,iElem)
                            !DCDZD = DDCZ*( VarEstP(iLayer,iElem) - VarEstP(iLayer-1,iElem) )/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))/2.)
                            DCDZD = ( VarEstP(iLayer,iElem) - VarEstP(iLayer-1,iElem) )* max(0., (MeshParam%Area(iElem)*DDCZ/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))/2.))-((1./2.)*abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem)))))
                                                    
                            !DCDZD = Max(0., DCDZD - 0.5*abs(MeshParam%Area(iElem)*HydroParam%w(iLayer,iElem)))
                        
                            !flux limiter to reduce vertical numerical diffusion in the bottom layer
                            DCDZD_psi = HydroParam%psi_cell(iLayer,iElem)*abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem)))*(VarEstP(iLayer,iElem) - VarEstP(iLayer-1,iElem))
                        
                            !DCDZU = Max(0., DCDZU - 0.5*abs(Area(iElem)*w(iLayer+1,iElem)))
                            !DCDZD = Max(0., DCDZD - 0.5*abs(MeshParam%Area(iElem)*HydroParam%w(iLayer,iElem)))
                            ! Vertical gradient of Diffusivity
                            DC2DZ2 = (DCDZU - DCDZD)/MeshParam%Area(iElem) 
                        
                            DC2DZ2_psi = (DCDZU_psi - DCDZD_psi)/MeshParam%Area(iElem)
                        
                            ! Transport Equation Solution
                            If (index== 0) Then !No boundary condition
                                Advection = 0.
                                Do iEdge = 1,4
                                    Face = MeshParam%Edge(iEdge,iElem)
                                    l = MeshParam%Left(Face)
                                    r = MeshParam%Right(Face)
                                    Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                                EndDo
                                    VarEst(iLayer,iElem) = MAX(NearZero,(dt/nSubCycle*(Dispersion-Advection+DC2DZ2)-dt/nSubCycle/2.*Dispersion_psi-dt/nSubCycle/2.*DC2DZ2_psi+dtday/nSubCycle*dVar(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem) - ((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*dt/nSubCycle)/HydroParam%DZitau(iLayer,iElem))
                            Else
                                Advection = 0.
                                num=0
                                Do iEdge = 1,4
                                    Face = MeshParam%Edge(iEdge,iElem)
                                    l = MeshParam%Left(Face)
                                    r = MeshParam%Right(Face)
                                    If (HydroParam%IndexInflowEdge(Face) > 0.or.HydroParam%IndexWaterLevelEdge(Face) > 0) Then !Inflow/Outflow boundary condition
                                        If (HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)>=0) Then ! Outflow Loading
                                            Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)    
                                        Else !InFlow Loading
                                            Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*uLoadVar(iLayer,iElem)/MeshParam%Area(iElem) 
                                        EndIf
                                        num=1
                                    Else
                                        Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                                    EndIf
                                EndDo
                                if (num==1) then
                                    VarEst(iLayer,iElem) = MAX(NearZero,(dt/nSubCycle*(Dispersion-Advection+DC2DZ2)+dtday/nSubCycle*dVar(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem) - ((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*dt/nSubCycle)/HydroParam%DZitau(iLayer,iElem))
                                else
                                    VarEst(iLayer,iElem) = MAX(NearZero,(dt/nSubCycle*(Dispersion-Advection+DC2DZ2)-dt/nSubCycle/2.*Dispersion_psi-dt/nSubCycle/2.*DC2DZ2_psi+dtday/nSubCycle*dVar(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem) - ((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*dt/nSubCycle)/HydroParam%DZitau(iLayer,iElem))
                                endif
                            EndIf

                        Else !Intermediary layer

                            ! Calcula a Concentração na camada
                            ! Top
                            If ((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))>=0) Then
                                CONCU = VarEstP(iLayer,iElem)                         ! Top face concentration is equal to the concentration in the layer    
                            Else
                                CONCU = VarEstP(iLayer+1,iElem)                         ! Top face concentration is equal to the concentration in the layer
                            EndIF
                        
                            ! Bottom
                            If ((HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))>=0) Then
		                        CONCD = VarEstP(iLayer-1,iElem) 
                            Else
                                CONCD = VarEstP(iLayer,iElem)
                            EndIf
        
		                    ! Vertical Diffusivity Flux in the top layer
                            DDCZ = HydroParam%VerEddyDiffCell(iLayer+1,iElem)
                            !DCDZU = DDCZ*( VarEstP(iLayer+1,iElem) - VarEstP(iLayer,iElem) )/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer+1,iElem))/2.)
                            !DCDZU = Max(0., DCDZU - 0.5*abs(MeshParam%Area(iElem)*HydroParam%w(iLayer+1,iElem)))
                            DCDZU = (VarEstP(iLayer+1,iElem) - VarEstP(iLayer,iElem))*max(0., (MeshParam%Area(iElem)*DDCZ/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer+1,iElem))/2.))-((1./2.)*abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem)))))
                        
                            !flux limiter to reduce vertical numerical diffusion in the top layer
                            DCDZU_psi = HydroParam%psi_cell(iLayer+1,iElem)*abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem)))*(VarEstP(iLayer+1,iElem) - VarEstP(iLayer,iElem))
                        
                        
                            ! Vertical Diffusivity Flux in the bottom layer
                            DDCZ = HydroParam%VerEddyDiffCell(iLayer,iElem)
                            !DCDZD = DDCZ*( VarEstP(iLayer,iElem) - VarEstP(iLayer-1,iElem) )/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))/2.)
                            !DCDZD = Max(0., DCDZD - 0.5*abs(MeshParam%Area(iElem)*HydroParam%w(iLayer,iElem)))
                            DCDZD = ( VarEstP(iLayer,iElem) - VarEstP(iLayer-1,iElem) )* max(0., (MeshParam%Area(iElem)*DDCZ/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))/2.))-((1./2.)*abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem)))))
                        
                            !flux limiter to reduce vertical numerical diffusion in the bottom layer
                            DCDZD_psi = HydroParam%psi_cell(iLayer,iElem)*abs(MeshParam%Area(iElem)*(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.-HydroParam%Theta)*HydroParam%wt(iLayer,iElem)))*(VarEstP(iLayer,iElem) - VarEstP(iLayer-1,iElem))
                        
                            ! Vertical gradient of Diffusivity
                            DC2DZ2 = (DCDZU - DCDZD)/MeshParam%Area(iElem) 
                        
                            DC2DZ2_psi = (DCDZU_psi - DCDZD_psi)/MeshParam%Area(iElem)
                        
                        
                            ! Transport Equation Solution
                            If (index== 0) Then !No boundary condition
                                Advection = 0.
                                Do iEdge = 1,4
                                    Face = MeshParam%Edge(iEdge,iElem)
                                    l = MeshParam%Left(Face)
                                    r = MeshParam%Right(Face)
                                    Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                                EndDo
                                    VarEst(iLayer,iElem) = MAX(NearZero,(dt/nSubCycle*(Dispersion-Advection+DC2DZ2)-dt/nSubCycle/2.*Dispersion_psi-dt/nSubCycle/2.*DC2DZ2_psi+dtday/nSubCycle*dVar(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem) - ((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*dt/nSubCycle)/HydroParam%DZitau(iLayer,iElem))
                            Else
                                Advection = 0.
                                num=0
                                Do iEdge = 1,4
                                    Face = MeshParam%Edge(iEdge,iElem)
                                    l = MeshParam%Left(Face)
                                    r = MeshParam%Right(Face)
                                    If (HydroParam%IndexInflowEdge(Face) > 0.or.HydroParam%IndexWaterLevelEdge(Face) > 0) Then !Inflow/Outflow boundary condition
                                        If (HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)>=0) Then ! Outflow Loading
                                            Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)    
                                        Else !InFlow Loading
                                            Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*uLoadVar(iLayer,iElem)/MeshParam%Area(iElem) 
                                        EndIf
                                        num=1
                                    Else
                                        Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                                    EndIf
                                EndDo
                                if (num==1) then
                                    VarEst(iLayer,iElem) = MAX(NearZero,(dt/nSubCycle*(Dispersion-Advection+DC2DZ2)+dtday/nSubCycle*dVar(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem) - ((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*dt/nSubCycle)/HydroParam%DZitau(iLayer,iElem))
                                else
                                    VarEst(iLayer,iElem) = MAX(NearZero,(dt/nSubCycle*(Dispersion-Advection+DC2DZ2)-dt/nSubCycle/2.*Dispersion_psi-dt/nSubCycle/2.*DC2DZ2_psi+dtday/nSubCycle*dVar(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZitau(iLayer,iElem) - ((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*dt/nSubCycle)/HydroParam%DZitau(iLayer,iElem))
                                endif
                            EndIf

                        EndIf
                    EndIf
                EndDo
            
                ! Copy concentration Above the Free-Surface
                Do iLayer = HydroParam%ElCapitalM(iElem) + 1, MeshParam%KMAX
                    VarEst(iLayer,iElem) = VarEst(HydroParam%ElCapitalM(iElem),iElem)
                EndDo
                ! Copy concentration Below the Bottom
                Do iLayer = 1, HydroParam%ElSmallm(iElem) - 1
                    VarEst(iLayer,iElem) = VarEst(HydroParam%ElSmallm(iElem),iElem)
                EndDo
            
        EndDo !Loop Elem
        ! preencher valores fora dos limites da camada
        VarEstP = VarEst
        HydroParam%DZitau = HydroParam%DZistau
    EndDo ! Subcycling time stepping

    Return
End Subroutine UpWindCWC_HR_GATS