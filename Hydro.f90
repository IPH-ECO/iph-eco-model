Subroutine Hydro(HydroParam,MeshParam,MeteoParam,dt,time,Simtime)
    
    ! Semi-Implicit Solution for Shallow Water Equations in Structured Grid
    ! Based on: 
    ! [1] Casulli, V.; Walters, R.A. An unstructured grid, three-dimensional model based on the shallow water equations.
    !     International Journal for Numerical Methods in Fluids, 32 (2000), p. 331 – 348
    ! [2] Casulli, V. A high-resolution wetting and drying algorithm for free-surface hydrodynamics.
    !     International Journal for Numerical Methods in Fluids, 60 (2009), p. 391-408
    ! [3] Casulli, V. A conservative semi-implicit method for coupled surface–subsurface flows in regional scale
    !   International Journal for Numerical Methods in Fluids, v. 79, p. 199-214, 2015.
    
    ! Input:
    ! Initial Conditions and Bathymetric Information
    ! Output:
    ! eta -> Free-Surface Elevation
    ! u   -> Water Velocity in Normal Direction in Each Edge
    
    ! List of Modifications:
    !   16.12.2014: Routine Implementation      (Rafael Cavalcanti)
    !   16.12.2014: Routine Implementation      (Carlos Ruberto)
    !   07.05.2020: Subsurface flows adaptation (Cayo Lopes)
    ! Programmer: Rafael Cavalcanti
    
    !$ use omp_lib
    Use MeshVars
    Use Hydrodynamic
    Use Meteorological
    Use ELM
    
    Implicit None
    Integer:: iElem, iEdge, lEdge, iNode, iNewton, iLayer, INFO, iLayer_bar
    Integer:: r, l, Sig, Face, Pij, DIM, ie, j, k, iNode1,iNode2, cont
    Real:: V, DV, SumRHS, SumH, res, SumLhS,gamma,teste
    Real:: dzp, dzm, SumW, rAux, Aux, VerEddyViscUp,VerEddyViscDown, vel
    Real:: Chezy,sum1,sum0,rhoairCell
    Real:: NearZero = 1e-10
    Real:: dt, man,raioh,slope,SimTime,time
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam
    type(MeteorologicalParam) :: MeteoParam
    Real:: aTh(MeshParam%KMax), bTh(MeshParam%KMax), cTh(MeshParam%KMax), iBGhost(3,3,MeshParam%KMax,MeshParam%nElem)
    
    Do iEdge = 1,MeshParam%nEdge
        HydroParam%epson(:,iEdge) = 0.
        !!Sponge layer
        !If (iEdge>2641) Then !2641 for dx=dy=0.025 !5281 for dx=dy=0.0125
        !    Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)  
        !        HydroParam%epson(iLayer,iEdge) = 0.5*(((MeshParam%EdgeBary(1,iEdge)-MeshParam%EdgeBary(1,2641))/(MeshParam%EdgeBary(1,MeshParam%nEdge)-MeshParam%EdgeBary(1,2641)))**2.)*((HydroParam%Z(HydroParam%Smallm(iEdge),iEdge)-HydroParam%Z(iLayer,iEdge))/(HydroParam%Z(HydroParam%Smallm(iEdge),iEdge)-HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge)))
        !    EndDo
        !Else
        !     HydroParam%epson(:,iEdge) = 0.
        !EndIf
    EndDo
    !
    ! 0. Compute turbulence
    Call Turbulence(HydroParam,MeshParam,MeteoParam,dt)
    !HydroParam%iConv = 3
    
    
    HydroParam%DZK(:)   = 0.d0
        
    !! 1. Convective Term
    If (HydroParam%iConv == 0) Then
        HydroParam%Fw = HydroParam%w
        HydroParam%Fu = HydroParam%u
        call FuFv(HydroParam,MeshParam,dt)
        Do iEdge = 1,MeshParam%nEdge
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge) 
                HydroParam%Fu(iLayer,iEdge) = (1.-HydroParam%epson(iLayer,iEdge))*HydroParam%Fu(iLayer,iEdge) 
            EndDo
            
            if (HydroParam%IndexWaterLevelEdge(iEdge)>0) then
                l = MeshParam%Left(iEdge)
                Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                    cont=0
                    vel=0
                    Do lEdge=1,4 
                        If (HydroParam%Fu(iLayer,MeshParam%Edge(lEdge,l))>0.and.HydroParam%IndexWaterLevelEdge(MeshParam%Edge(lEdge,l))==0) Then
                            cont=cont+1
                            vel = vel + HydroParam%Fu(iLayer,MeshParam%Edge(lEdge,l))
                        EndIf
                    EndDo
                    if (cont>0) then    
                        HydroParam%Fu(iLayer,iEdge) = vel/cont
                    endif
                EndDo    
            endif 
        EndDo
        
    ElseIf (HydroParam%iConv == 1) Then
        
    ElseIf (HydroParam%iConv == 2) Then
        !Call Convective
    ElseIf (HydroParam%iConv == 3) Then ! Neglect Nonlinear Convection
        !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iEdge)
        Do iEdge = 1,MeshParam%nEdge
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge) 
                HydroParam%Fu(iLayer,iEdge) = (1.-HydroParam%epson(iLayer,iEdge))*HydroParam%u(iLayer,iEdge)
            EndDo
        EndDo
        HydroParam%Fw = HydroParam%w
        !!$OMP end parallel do
    EndIf
    
    !Getting hydrodynamic boundary condition values at current step time
    Call GetHydroBoundaryConditions(HydroParam,MeshParam,time)

    
    ! 2. Baroclinic pressure effect
    Call Pressure(HydroParam,MeshParam,dt)
    
    ! 3. Viscosity and Coriolis effect
    Call ExplicitTerms(HydroParam,MeshParam,dt)
    
    ! 4. Get wind velocity components
    Call WindVelocity(HydroParam,MeshParam,MeteoParam)
    
    ! 5. Vertical Water Balance (Precipitation and Evaporation)
    Call VerticalWB(HydroParam,MeshParam,MeteoParam,dt)
    
    
    ! 6. Assemble Matrix
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam,MeteoParam) private(iEdge,iLayer,l,r,Chezy,rhoairCell,aTh,bTh,cTh,dzp,dzm,DIM,NearZero,dt)
    Do iEdge = 1,MeshParam%nEdge
        
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
       
        ! 6.1 Get roughness 
        If (r == 0) Then
            If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
                Chezy = HydroParam%Rug(l)
            ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
                Chezy = Max(HydroParam%Pcri,HydroParam%H(iEdge))**(1./6.)/(HydroParam%Rug(l)+NearZero)
            ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
                Chezy = 18.*log10(12.*Max(HydroParam%Pcri,HydroParam%H(iEdge))/(HydroParam%Rug(l)/30.+NearZero))
            EndIf
            rhoairCell = MeteoParam%rhoair(l)
        Else
            If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
                Chezy = 0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))
            ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
                Chezy = Max(HydroParam%Pcri,HydroParam%H(iEdge))**(1./6.)/(0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))+NearZero)
            ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
                Chezy = 18.*log10(12.*Max(HydroParam%Pcri,HydroParam%H(iEdge))/(0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))/30.+NearZero))
            EndIf
            rhoairCell = 0.5*(MeteoParam%rhoair(l) + MeteoParam%rhoair(r))
        EndIf

        
        ! 6.2 Get Inflow/Outflow and Normal Depth Boundary Condition 
        ! If a face is dry, set the normal velocity to zero
        If (HydroParam%H(iEdge) <= HydroParam%PCRI+NearZero) Then
            HydroParam%u(:,iEdge)  = 0.
            HydroParam%Fu(:,iEdge) = 0.
        EndIf
        
        ! 6.3 Get Shear stresses in the edges
        Call Tension(HydroParam%iWindStress,HydroParam%BottomTensionFlag,iEdge,HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),HydroParam%g,HydroParam%uxy(HydroParam%Smallm(iEdge),1,iEdge),HydroParam%uxy(HydroParam%Smallm(iEdge),2,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge),Chezy,rhoairCell,HydroParam%rho0,HydroParam%windDragConstant,HydroParam%WindVel(:,iEdge),HydroParam%WindXY(:,iEdge),HydroParam%GammaT,HydroParam%GammaB)

        ! 6.4 Assemble Matrix G 
        ! Different from the article [1]
        ! Here the multiplication by DZj is done when we calculate iAG
        ! Only at the First Layer the Wind force takes place
        ! When there is no neighbour, there is no flux ( (1-Theta)*g*(dt/CirDistance(iEdge))*(eta(r) - eta(l)) == 0 )
        Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)      
            If ( iLayer == HydroParam%CapitalM(iEdge) ) Then
                If ( r == 0 .or. HydroParam%H(iEdge) <= HydroParam%PCRI+NearZero) Then
                    HydroParam%Gu(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)*HydroParam%Fu(iLayer,iEdge) !+ dt*GammaT*WindVel
                Else
                    HydroParam%Gu(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)*(HydroParam%Fu(iLayer,iEdge) - (1-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*( HydroParam%g*(HydroParam%eta(r) - HydroParam%eta(l)) + (HydroParam%q(iLayer,r) - HydroParam%q(iLayer,l)))) + dt*HydroParam%GammaT*HydroParam%WindVel(1,iEdge)
                EndIf
            Else
                If ( r == 0 .or. HydroParam%H(iEdge) <= HydroParam%PCRI+NearZero) Then
                    HydroParam%Gu(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)*HydroParam%Fu(iLayer,iEdge)
                Else
                    HydroParam%Gu(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)*(HydroParam%Fu(iLayer,iEdge) - (1-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*( HydroParam%g*(HydroParam%eta(r) - HydroParam%eta(l)) + (HydroParam%q(iLayer,r) - HydroParam%q(iLayer,l))))
                EndIf
            EndIf
        EndDo
        
        ! 6.5 Assemble the TriDiagonal System in Vertical Direction
        If ( HydroParam%Smallm(iEdge) == HydroParam%CapitalM(iEdge) ) Then        ! Only One Vertical Layer
            If ( r == 0 .or. HydroParam%H(iEdge) <= HydroParam%PCRI+NearZero) Then
                HydroParam%GammaB = 0. 
            EndIf
            
            HydroParam%iADZ(HydroParam%Smallm(iEdge),iEdge) = HydroParam%H(iEdge)/( HydroParam%H(iEdge) + dt*HydroParam%GammaB + NearZero)
            HydroParam%iAG(HydroParam%Smallm(iEdge),iEdge)  = HydroParam%Gu(HydroParam%Smallm(iEdge),iEdge)/( HydroParam%H(iEdge) + dt*HydroParam%GammaB + NearZero )
            HydroParam%DZiADZ(iEdge)             = HydroParam%H(iEdge)*HydroParam%H(iEdge)/( HydroParam%H(iEdge) + dt*HydroParam%GammaB  + NearZero )
            HydroParam%DZiAG (iEdge)             = HydroParam%H(iEdge)*HydroParam%Gu(HydroParam%Smallm(iEdge),iEdge)/( HydroParam%H(iEdge) + dt*HydroParam%GammaB + NearZero )
            
        Else
            ! 6.5.1 Assemble Matrix A
            If ( r == 0 .or. HydroParam%H(iEdge) <= HydroParam%PCRI+NearZero) Then
                HydroParam%GammaB = 0. 
            EndIf
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                If ( iLayer == HydroParam%Smallm(iEdge) ) Then
                    dzp         = 0.5*( HydroParam%DZj(iLayer,iEdge) + HydroParam%DZj(iLayer+1,iEdge) )       ! Dz at Upper Interface
                    aTh(iLayer) = 0. !-dt*VerEddyVisc(iLayer+1,iEdge)/dzp !-dt*nuz/dzp   !                                        ! Lower Diagonal 
                    bTh(iLayer) = HydroParam%DZj(iLayer,iEdge) + dt*HydroParam%VerEddyVisc(iLayer+1,iEdge)*( 1/dzp ) + HydroParam%GammaB*dt + NearZero!DZj(iLayer,iEdge) + dt*nuz*( 1/dzp ) + GammaB*dt      ! Diagonal
                    cTh(iLayer) = -dt*HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzp !0.                                                    ! Upper Diagonal 

                Elseif ( iLayer == HydroParam%CapitalM(iEdge) ) Then
                    dzm         = 0.5*( HydroParam%DZj(iLayer,iEdge) + HydroParam%DZj(iLayer-1,iEdge) )       ! Dz at Lower Interface
                    aTh(iLayer) = -dt*HydroParam%VerEddyVisc(iLayer,iEdge)/dzm !0.                                                    ! Lower Diagonal 
                    bTh(iLayer) = Max(HydroParam%Pcri,HydroParam%DZj(iLayer,iEdge) + dt*HydroParam%VerEddyVisc(iLayer,iEdge)*( 1/dzm ) + HydroParam%GammaT*dt + NearZero) !DZj(iLayer,iEdge) + dt*nuz*( 1/dzm ) + GammaT*dt      ! Diagonal
                    cTh(iLayer) = 0. !-dt*VerEddyVisc(iLayer,iEdge)/dzm !-dt*nuz/dzm                                           ! Upper Diagonal

                Else
                    dzp         = 0.5*( HydroParam%DZj(iLayer,iEdge) + HydroParam%DZj(iLayer+1,iEdge) )       ! Dz at Upper Interface 
                    dzm         = 0.5*( HydroParam%DZj(iLayer,iEdge) + HydroParam%DZj(iLayer-1,iEdge) )       ! Dz at Lower Interface 
                    aTh(iLayer) = -dt*HydroParam%VerEddyVisc(iLayer,iEdge)/dzm !-dt*VerEddyVisc(iLayer+1,iEdge)/dzp !-dt*nuz/dzm                                           ! Lower Diagonal 
                    bTh(iLayer) = HydroParam%DZj(iLayer,iEdge) + dt*HydroParam%VerEddyVisc(iLayer,iEdge)*( 1/dzm ) +  dt*HydroParam%VerEddyVisc(iLayer+1,iEdge)*( 1/dzp ) + NearZero         ! Diagonal
                    cTh(iLayer) = -dt*HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzp !-dt*VerEddyVisc(iLayer,iEdge)/dzm !-dt*nuz/dzp                                           ! Upper Diagonal 

                EndIf
            EndDo

            ! 6.5.2 Assemble Matrix: iAG, iADZ, DZiADZ, DZiAG
            DIM = size(HydroParam%DZj(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),1)
            Call solve_tridiag(cTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),bTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),aTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),HydroParam%DZj(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),HydroParam%iADZ(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),DIM)
            Call solve_tridiag(cTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),bTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),aTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),HydroParam%Gu(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),HydroParam%iAG(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),DIM)

            HydroParam%DZiADZ(iEdge) = Dot_Product(HydroParam%DZj(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),HydroParam%iADZ(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iedge),iEdge) )
            HydroParam%DZiAG (iEdge) = Dot_Product(HydroParam%DZj(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),HydroParam%iAG (HydroParam%Smallm(iEdge):HydroParam%CapitalM(iedge),iEdge) )
        EndIf       ! If ( Smallm(iEdge) == CapitalM(iEdge) ) Then
    EndDo
    !!$OMP end parallel do
    
    
    !6.6 Assemble Matrix DZK    !CAYO Flag

    If (MeshParam%iBedrock == 1) Then
        Do iEdge = 1,MeshParam%nEdge       
            If (HydroParam%Smallms(iEdge) == HydroParam%CapitalMs(iEdge).and. HydroParam%Smallms(iEdge) == HydroParam%CapitalM(iEdge)) Then !Rever 
                HydroParam%DZK(iEdge)   = HydroParam%DZj(HydroParam%Smallms(iEdge),iEdge)*MeshParam%Kj(HydroParam%Smallms(iEdge),iEdge) !Sediment Layer
            Else
                HydroParam%DZK(iEdge)   = Dot_Product(HydroParam%DZj(HydroParam%Smallms(iEdge):HydroParam%CapitalMs(iEdge),iEdge),MeshParam%Kj(HydroParam%Smallms(iEdge):HydroParam%CapitalMs(iEdge),iEdge) )
            EndIf
        EndDo
    EndIf
        
    Call Volume(HydroParam,MeshParam) !CAYO
      
    ! 7. Compute the New Free-Surface Elevation
    ! 7.1 Assemble the Right Hand Side (RHS)
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam,NearZero,dt) private(iElem,iEdge,SumRHS,SumLhS,gamma,Face,Pij)
    !HydroParam%etan = HydroParam%eta
    Do iElem = 1, MeshParam%nElem
        
        SumRHS = 0d0
        SumLhS = 0.d0
        gamma = 0.d0
        HydroParam%etaInf(iElem) = 0.d0
        
        Do iEdge = 1,4
 
            Face = MeshParam%Edge(iEdge,iElem)
            SumRHS = SumRHS + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*((1.d0-HydroParam%Theta)*(Dot_product( HydroParam%DZj(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face),HydroParam%u(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face)) + Dot_product(HydroParam%DZj(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face),HydroParam%us(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face)))+ HydroParam%Theta*HydroParam%DZiAG(Face))      

            ! 6.1.1 If there is a Pressure Boundary Condition
            Pij = MeshParam%Neighbor(iEdge,iElem) 
            
            If (Pij == 0.and.HydroParam%IndexWaterLevelEdge(Face)>0) Then
				If ((HydroParam%WaterLevel(HydroParam%IndexWaterLevelEdge(Face))-HydroParam%hj(Face))<HydroParam%PCRI+NearZero) THEN
                    HydroParam%etaInf(iElem) = HydroParam%hj(Face) + HydroParam%PCRI 
				Else
				    HydroParam%etaInf(iElem) = HydroParam%WaterLevel(HydroParam%IndexWaterLevelEdge(Face))
                    !HydroParam%etaInf(iElem) = 3.00*sin(((0.025*(1-1))*(0.67d0/0.4d0))+(2.0d0*HydroParam%pi*(Simtime)/43200.0d0-HydroParam%pi/2))
                EndIf
                SumLHS = SumLHS + ( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( HydroParam%etaInf(iElem) )*(HydroParam%g*HydroParam%Theta*dt*HydroParam%DZiADZ(Face) + HydroParam%DZK(Face))
                !SumLHS = SumLHS + ( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( ( HydroParam%etaInf(iElem) )*HydroParam%H(Face)**2 )/( HydroParam%H(Face) + dt*gamma )  ! 
				
            EndIf
            
        EndDo
        !HydroParam%rhs(iElem) =  MeshParam%Area(iElem)*V(HydroParam%eta(iElem)+HydroParam%etaplus(iElem),HydroParam%hb(iElem)) - dt*SumRHS + ( HydroParam%g*(HydroParam%Theta*dt)**2. )*SumLHS
        HydroParam%rhs(iElem) = HydroParam%Vol(iElem) - dt*SumRHS + (HydroParam%Theta*dt)*SumLHS
    EndDo
    
    !Do iElem = 1, MeshParam%nElem
    !    If (HydroParam%IndexWaterLevel(iElem)>0) Then
    !            HydroParam%rhs(iElem) = MeshParam%Area(iElem)*V(0.01*cos(2.0*HydroParam%pi*(Simtime-2.02d0/4.0d0)/2.02d0) ,HydroParam%hb(iElem))
    !    EndIf
    !EndDo
    ! 
    !!$OMP end parallel do
    
    ! 7.2 Newton Loop for Non-Linear Wet- and Dry-ing Algorithm [2]
    HydroParam%etan = HydroParam%eta

    Do iNewton = 1,100
        ! 7.2.1 Assemble Matrices F (System Matrix - Newton Method) and P (Derivative Matrix)
        HydroParam%P = 0.
        HydroParam%F = 0.
        Call MatOp(HydroParam%eta,HydroParam%Aeta,dt,HydroParam,MeshParam)
        Call Volume(HydroParam,MeshParam) !CAYO

        Do iElem = 1,MeshParam%nElem
            !HydroParam%F(iElem) = MeshParam%Area(iElem)*V(HydroParam%eta(iElem),HydroParam%hb(iElem)) + HydroParam%Aeta(iElem) - HydroParam%rhs(iElem)
            HydroParam%F(iElem) = HydroParam%Vol(iElem) + HydroParam%Aeta(iElem) - HydroParam%rhs(iElem) !CAYO
        EndDo
        
        !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iElem,sumH)
        Do iElem = 1, MeshParam%nElem
            HydroParam%P(iElem) = MeshParam%Area(iElem)*dV(HydroParam%eta(iElem),HydroParam%sb(iElem)) !CAYO
            sumH = Sum( HydroParam%H(MeshParam%Edge(:,iElem)) )
            If (V(HydroParam%eta(iElem),HydroParam%hb(iElem)) < HydroParam%PCRI+NearZero) Then
                ! The Cell is Dry. The Solution is V(eta^(n+1)) = V(eta^(n))
                HydroParam%P(iElem) = MeshParam%Area(iElem)       ! We can't allow the row to have zeros in all elements - If P = 1, the water level remains the same
            EndIf
        EndDo
        !!$OMP end parallel do
            
        res = sqrt(sum(HydroParam%f**2))      ! Residual of the Method
        !Print*, 'iNewton = ',iNewton , 'res = ',res
        If ( res < 1e-8 ) Then
            continue
            exit
        EndIf
        ! 7.2.2 Compute the New Free-Surface Elevation
        Call CGOp(HydroParam%F,HydroParam%Deta,dt,HydroParam,MeshParam)
        !!$OMP parallel do default(none) shared(HydroParam,MeshParam) private(iElem)
        Do iElem = 1, MeshParam%nElem
            HydroParam%eta(iElem) = HydroParam%eta(iElem) - HydroParam%Deta(iElem)
        EndDo
        !!$OMP end parallel do
    EndDo
    
    ! 7.3 Compute nodal elevations for tang. vel.
    HydroParam%petan=HydroParam%peta !store for transport eqs. 
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iNode,sum1,sum0,ie,j)
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
    !!$OMP end parallel do
    
    ! 8. Water Level Corrective Step
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iElem)
    Do iElem = 1,MeshParam%nElem
        If ( HydroParam%eta(iElem) - HydroParam%sb(iElem) < HydroParam%PCRI/2. ) Then
            HydroParam%eta(iElem) = HydroParam%sb(iElem) + HydroParam%PCRI/2.
        EndIf
    EndDo
    !!$OMP end parallel do
    ! 10. Updating free-surface water level and  Vertical Mesh Spacing
    HydroParam%DZjt = HydroParam%DZj 
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iLayer,iEdge,l,r,NearZero)
    Do iEdge = 1,MeshParam%nEdge
        
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
        ! 9.1 Compute Index Smallm(j) and CapitalM(j)
        ! Smallm and CapitalM are related to mj and Mj in [1]
        If (r == 0) Then
            HydroParam%H(iEdge) = Max( HydroParam%PCRI, -HydroParam%sj(iEdge) + HydroParam%eta(l) ) !CAYO
        Else
            HydroParam%H(iEdge) = Max( HydroParam%PCRI,-HydroParam%sj(iEdge) + HydroParam%eta(l), -HydroParam%sj(iEdge) + HydroParam%eta(r) )!CAYO
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
        !9.2 Compute Elevation in the Edges
        !Do iLayer = HydroParam%Smallm(iEdge) + 1, HydroParam%CapitalM(iEdge)
        !    HydroParam%Z(iLayer,iEdge) = MeshParam%LIMCAMAUX( MeshParam%KMax-iLayer+1) !zL + (iLayer - 1)*dz                  ! Equidistant Core Grid
        !EndDo
        !HydroParam%Z(HydroParam%Smallm(iEdge),iEdge)     = HydroParam%hj(iEdge)                 ! Bottom
        !If (r == 0) Then
        !    HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%eta(l) !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
        !Else
        !    If ( HydroParam%H(iEdge) <= HydroParam%PCRI+NearZero ) Then
        !        HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%eta(l),HydroParam%eta(r))
        !    Else
        !        HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = 0.5*(HydroParam%eta(l)+HydroParam%eta(r)) !Max(HydroParam%eta(l),HydroParam%eta(r)) !0.5*(eta(l)+eta(r)) !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
        !    EndIf
        !EndIf
        
        !! 9.3 Compute the Vertical Mesh Spacing
        !Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
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
                MeshParam%Kj(iLayer,iEdge) = 0.001
            EndIf
        EndDo        
    EndDo
    !!$OMP end parallel do
    
    ! 11. Compute the New Vertical Velocity Field
    ! w(Smallm-1,:) = 0. -> No flux through the Bottom
    ! w(CapitalM,:) = 0. -> No flux through the Surface
    ! The velocity is defined in the barycenter of the Top Face of each Cell
    ! *Obs: The outward direction is from Bottom to Top.
    HydroParam%DZit = HydroParam%DZi
    HydroParam%wt = HydroParam%w
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iLayer,iElem,iEdge,Face,SumW,NearZero)
    Do iElem = 1, MeshParam%nElem
       
        ! 11.1 Define the range between Bottom and Top Layers. (Tricky part - Not explained in [1])
        ! In a real world application, the bathymetry might change between the faces of the prism.
        ! We need to find "Max(Smallm(:))" and "Min(CapitalM(:))" in each face. 
        ! Then we can calculate the Vertical Velocity only in those Layers
        ! *Obs: The velocity in faces with Smallm < Max(Smallm(:)) == 0. The same goes to faces with CapitalM > Min(CapitalM(:))
        ! 4.1.1 Find "mi" and "Mi", the Bottom and Top Layers range for the Element, not for the faces
        ! Set mi
        
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
        !! 11.2 Update the Element Vertical Spacing
        !
        !Do iLayer = HydroParam%ElSmallm(iElem)+1, HydroParam%ElCapitalM(iElem)
        !    HydroParam%Ze(iLayer,iElem) = MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)      ! Equidistant Core Grid
        !EndDo
        !HydroParam%Ze(HydroParam%ElSmallm(iElem),iElem)     = HydroParam%hb(iElem)                    ! Bottom
        !HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%eta(iElem) !- hb(iElem)       ! Free-Surface
        !Do iLayer = 1, HydroParam%ElSmallm(iElem) - 1
        !    HydroParam%Ze(iLayer,iElem) = HydroParam%hb(iElem)
        !EndDo
        !Do iLayer = HydroParam%ElCapitalM(iElem)+2, MeshParam%KMax+1
        !    HydroParam%Ze(iLayer,iElem) = HydroParam%eta(iElem)
        !EndDo
        !
        !! 11.3 Compute the Vertical Mesh Spacing
        !Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
        !    HydroParam%DZi(iLayer,iElem) = HydroParam%Ze(iLayer+1,iElem) - HydroParam%Ze(iLayer,iElem)
        !    HydroParam%DZi(iLayer,iElem) = Max(HydroParam%Pcri,HydroParam%DZi(iLayer,iElem))
        !EndDo        
        !
        !
        !Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem)
        !    If (iLayer == HydroParam%ElSmallms(iElem)) Then
        !        HydroParam%Zb(iLayer,iElem) = 0.5*( HydroParam%Ze(iLayer,iElem) + HydroParam%Ze(iLayer+1,iElem) ) !0.5*( 0. + Ze(iLayer+1,iElem) )   (verificar com Rafael)               
        !    Else
        !        HydroParam%Zb(iLayer,iElem) = 0.5*( HydroParam%Ze(iLayer,iElem) + HydroParam%Ze(iLayer+1,iElem) )
        !    EndIf
        !EndDo    
        
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
        If ( HydroParam%ElSmallm(iElem) > 1 ) Then   ! Cayo 
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
        HydroParam%Ze(HydroParam%ElSmallms(iElem),iElem)     = HydroParam%sb(iElem)                    ! Bottom
        If ( HydroParam%eta(iElem) - HydroParam%sb(iElem) <= HydroParam%PCRI+NearZero ) Then
            HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%sb(iElem) + HydroParam%PCRI !- hb(iElem)       ! Free-Surface (verificar com Rafael)
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
        Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem) !CAYO
            If (iLayer == HydroParam%ElSmallms(iElem)) Then
                HydroParam%Zb(iLayer,iElem) = 0.5*( HydroParam%Ze(iLayer,iElem) + HydroParam%Ze(iLayer+1,iElem) ) !0.5*( 0. + Ze(iLayer+1,iElem) )   (verificar com Rafael)               
            Else
                HydroParam%Zb(iLayer,iElem) = 0.5*( HydroParam%Ze(iLayer,iElem) + HydroParam%Ze(iLayer+1,iElem) )
            EndIf 
        EndDo
        !Porosity
        Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem) ! Cayo (Loop adicionado)
            If (iLayer >= HydroParam%ElSmallm(iElem) .or.  HydroParam%ElSmallms(iElem) ==  HydroParam%ElSmallm(iElem)) Then
                MeshParam%ei(iLayer,iElem) = 1
            Else
                MeshParam%ei(iLayer,iElem) = 0.5
            EndIf
        EndDo            
    EndDo
    
    
    ! 9. Compute the New Horizontal and Vertical Velocity Field 
    HydroParam%ut = HydroParam%u
    Call uvelocity(HydroParam,MeshParam,dt)
    
    !Call utangvelocity(HydroParam,MeshParam,dt)
    
    HydroParam%wt = HydroParam%w
    !Call wvelocity(HydroParam,MeshParam,dt)
    
    Do iElem = 1, MeshParam%nElem
         !11.4 Evaluate the new Vertical Velocity Field
        HydroParam%w(HydroParam%ElSmallm(iElem),iElem)       = 0.       ! No Flux through the Bottom
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
            SumW = 0.
            Do iEdge = 1,4
                Face = MeshParam%Edge(iEdge,iElem)
                SumW = SumW + (Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*HydroParam%DZjt(iLayer,Face)*HydroParam%u(iLayer,Face))
            EndDo
            HydroParam%w(iLayer+1,iElem) = HydroParam%w(iLayer,iElem) - SumW/MeshParam%Area(iElem)         ! w(k+1/2,iElem) [1]
        EndDo
    EndDo
    
	!Do iElem = 1, MeshParam%nElem
 !       If (iElem==5) Then
 !           !Balanço de volume baseado na equação da continuidade
 !           HydroParam%SumVer = 0.
 !           Do iLayer =  HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) !ElSmallm(iElem), ElSmallm(iElem) !ElCapitalM(iElem), ElCapitalM(iElem) !
 !               If (iLayer==HydroParam%ElCapitalM(iElem)) Then
 !                   SumW = 0.
 !                   Do iEdge = 1,4
 !                       Face = MeshParam%Edge(iEdge,iElem)
 !                       Do iLayer_bar =  HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
 !                           SumW = SumW + (HydroParam%theta*dt*Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*HydroParam%DZjt(iLayer_bar,Face)*HydroParam%u(iLayer_bar,Face)) + ((1.-HydroParam%theta)*dt*Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*HydroParam%DZjt(iLayer_bar,Face)*HydroParam%ut(iLayer_bar,Face))
 !                       EndDo
 !                   EndDo
 !                   HydroParam%SumVer = HydroParam%SumVer + SumW + (HydroParam%eta(iElem)-HydroParam%etan(iElem))*MeshParam%Area(iElem)
 !               Else
 !                   SumW = 0.
 !                   Do iEdge = 1,4
 !                       Face = MeshParam%Edge(iEdge,iElem)
 !                       SumW = SumW + (Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*HydroParam%DZjt(iLayer,Face)*HydroParam%u(iLayer,Face))
 !                   EndDo
 !                   HydroParam%SumVer = HydroParam%SumVer + SumW + (HydroParam%w(iLayer+1,iElem)-HydroParam%w(iLayer,iElem))*MeshParam%Area(iElem)
 !               EndIf
 !           EndDo
 !           HydroParam%SumVerAcum = HydroParam%SumVerAcum + HydroParam%SumVer
 !           !!Balanço de volume baseado na equação da continuidade integrada
 !           !SumVer = -Sig(iElem,Right(Edge(1,iElem)),MeshParam%Left(Edge(1,iElem)))*(Theta*u(1,Edge(1,iElem))*DZjt(1,Edge(1,iElem))+(1.-Theta)*ut(1,Edge(1,iElem))*DZjt(1,Edge(1,iElem)))*DX*DT-Sig(iElem,Right(Edge(2,iElem)),MeshParam%Left(Edge(2,iElem)))*(Theta*u(1,Edge(2,iElem))*DZjt(1,Edge(2,iElem))+ (1.-Theta)*ut(1,Edge(2,iElem))*DZjt(1,Edge(2,iElem)))*DY*simParam%DT-Sig(iElem,Right(Edge(3,iElem)),MeshParam%Left(Edge(3,iElem)))*(Theta*u(1,Edge(3,iElem))*DZjt(1,Edge(3,iElem))+ (1.-Theta)*ut(1,Edge(3,iElem))*DZjt(1,Edge(3,iElem)))*DY*simParam%DT-Sig(iElem,Right(Edge(4,iElem)),MeshParam%Left(Edge(4,iElem)))*(Theta*u(1,Edge(4,iElem))*DZjt(1,Edge(4,iElem))+(1.-Theta)*ut(1,Edge(4,iElem))*DZjt(1,Edge(4,iElem)))*DY*simParam%DT - (eta(iElem)-etaplus(ielem)-etan(iElem))*DX*DY 
 !           !SumVerAcum = SumVerAcum + SumVer
 !           !Write(*,*) '1-',SumVer,eta(iElem),sDSal(1,iElem)!SumVer, SumVerAcum!, VT2(I,J+1,1), VT1(I,J+1,1)
 !           Print*,HydroParam%SumVer,HydroParam%SumVerAcum
 !           pause
 !       EndIf
 !   EndDo      
    
    
    
    !print*, HydroParam%w
    !!$OMP end parallel do
    !10. Non-hydrostatic correction
    
    If (HydroParam%iNonHydro==0.and.MeshParam%Kmax>1) Then 
        Call NonHydroPressure(HydroParam,MeshParam,dt)
    EndIf
    !HydroParam%um = HydroParam%u !CAYO
    !HydroParam%u = HydroParam%u + HydroParam%us !CAYO
    
    ! 12. Compute average and tangential velocities
    !Call Velocities(HydroParam,MeshParam)
    Call VelocitiesSUB(HydroParam,MeshParam) !CAYO
    
    !HydroParam%u = HydroParam%um !CAYO
    !HydroParam%um =  HydroParam%u + HydroParam%us !CAYO
    
   ! Call ExchangeTime(HydroParam,MeshParam,MeteoParam,dt)
    
    Return
End Subroutine Hydro