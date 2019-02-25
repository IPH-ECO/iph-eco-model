Subroutine Hydro(HydroParam,MeshParam,MeteoParam,dt,time)
    
    ! Semi-Implicit Solution for Shallow Water Equations in Structured Grid
    ! Based on: 
    ! [1] Casulli, V.; Walters, R.A. An unstructured grid, three-dimensional model based on the shallow water equations.
    !     International Journal for Numerical Methods in Fluids, 32 (2000), p. 331 – 348
    ! [2] Casulli, V. A high-resolution wetting and drying algorithm for free-surface hydrodynamics.
    !     International Journal for Numerical Methods in Fluids, 60 (2009), p. 391-408
    
    ! Input:
    ! Initial Conditions and Bathymetric Information
    ! Output:
    ! eta -> Free-Surface Elevation
    ! u   -> Water Velocity in Normal Direction in Each Edge
    
    ! List of Modifications:
    !   16.12.2014: Routine Implementation      (Rafael Cavalcanti)
    !   16.12.2014: Routine Implementation      (Carlos Ruberto)
    ! Programmer: Rafael Cavalcanti
    
    !$ use omp_lib
    Use MeshVars
    Use Hydrodynamic
    Use Meteorological
    Use ELM
    
    Implicit None
    Integer:: iElem, iEdge, iNode, iNewton, iLayer, INFO
    Integer:: r, l, Sig, Face, Pij, DIM, ie, j, k, iNode1,iNode2,time
    Double Precision:: V, DV, SumRHS, SumH, res, SumLhS,gamma,teste
    Double Precision:: dzp, dzm, SumW, rAux, Aux, VerEddyViscUp,VerEddyViscDown
    Double Precision:: Chezy,sum1,sum0,rhoairCell
    Double Precision:: NearZero = 1e-10
    Double Precision:: dt, man,raioh,slope
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam
    Double Precision:: aTh(MeshParam%KMax), bTh(MeshParam%KMax), cTh(MeshParam%KMax)

    
    ! 0. Compute turbulence
    Call Turbulence(HydroParam,MeshParam,MeteoParam,dt)
    !HydroParam%iConv = 3
    ! 1. Convective Term
    If (HydroParam%iConv == 0) Then
        Call FuFv(HydroParam,MeshParam,dt)
    ElseIf (HydroParam%iConv == 1) Then
        
    ElseIf (HydroParam%iConv == 2) Then
        !Call Convective
    ElseIf (HydroParam%iConv == 3) Then ! Neglect Nonlinear Convection
        !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iEdge)
        Do iEdge = 1,MeshParam%nEdge
            HydroParam%Fu(:,iEdge) = HydroParam%u(:,iEdge)              ! If Fu == u -> Neglect Nonlinear Convection
            
        EndDo
        !!$OMP end parallel do
    EndIf
    HydroParam%Fu = HydroParam%u
    
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

        If (l==1326) Then
            Continue
        EndIf
        
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
                    HydroParam%Gu(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)*(HydroParam%Fu(iLayer,iEdge) - (1-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*( HydroParam%g*(HydroParam%eta(r) - HydroParam%eta(l)))) + dt*HydroParam%GammaT*HydroParam%WindVel(1,iEdge)
                EndIf
            Else
                If ( r == 0 .or. HydroParam%H(iEdge) <= HydroParam%PCRI+NearZero) Then
                    HydroParam%Gu(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)*HydroParam%Fu(iLayer,iEdge)
                Else
                    HydroParam%Gu(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)*(HydroParam%Fu(iLayer,iEdge) - (1-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*( HydroParam%g*(HydroParam%eta(r) - HydroParam%eta(l)) ))
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
               
    ! 7. Compute the New Free-Surface Elevation
    ! 7.1 Assemble the Right Hand Side (RHS)
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam,NearZero,dt) private(iElem,iEdge,SumRHS,SumLhS,gamma,Face,Pij)
    Do iElem = 1, MeshParam%nElem
        
        SumRHS = 0
        SumLhS = 0.
        gamma = 0.
        
        Do iEdge = 1,4

            Face = MeshParam%Edge(iEdge,iElem)
            SumRHS = SumRHS + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*( (1-HydroParam%Theta)*Dot_product(HydroParam%Dzj(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face),HydroParam%u(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face)) + HydroParam%Theta*HydroParam%DZiAG(Face) )      

            ! 6.1.1 If there is a Pressure Boundary Condition
            Pij = MeshParam%Neighbor(iEdge,iElem)
            
            If (Pij == 0.and.HydroParam%IndexWaterLevel(iElem)>0) Then
				If ((HydroParam%WaterLevel(HydroParam%IndexWaterLevel(iElem))-HydroParam%hj(Face))<HydroParam%PCRI+NearZero) THEN
                    HydroParam%etaInf(iElem) = HydroParam%hj(Face) + HydroParam%PCRI 
				Else
				    HydroParam%etaInf(iElem) = HydroParam%WaterLevel(HydroParam%IndexWaterLevel(iElem))
				EndIf
                SumLHS = SumLHS + ( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( ( HydroParam%etaInf(iElem) )*HydroParam%H(Face)**2 )/( HydroParam%H(Face) + dt*gamma )  ! 
            EndIf
            
        EndDo
        HydroParam%rhs(iElem) = MeshParam%Area(iElem)*V(HydroParam%eta(iElem)+HydroParam%etaplus(iElem),HydroParam%hb(iElem)) - dt*SumRHS + ( HydroParam%g*(HydroParam%Theta*dt)**2. )*SumLHS

    EndDo
    !!$OMP end parallel do
    
    ! 7.2 Newton Loop for Non-Linear Wet- and Dry-ing Algorithm [2]
    HydroParam%etan = HydroParam%eta
    Do iNewton = 1,100
        ! 7.2.1 Assemble Matrices F (System Matrix - Newton Method) and P (Derivative Matrix)
        HydroParam%P = 0.
        HydroParam%F = 0.
        Call MatOp(HydroParam%eta,HydroParam%Aeta,dt,HydroParam,MeshParam)

        Do iElem = 1,MeshParam%nElem
            HydroParam%F(iElem) = MeshParam%Area(iElem)*V(HydroParam%eta(iElem),HydroParam%hb(iElem)) + HydroParam%Aeta(iElem) - HydroParam%rhs(iElem)
        EndDo
        
        !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iElem,sumH)
        Do iElem = 1, MeshParam%nElem
            HydroParam%P(iElem) = MeshParam%Area(iElem)*dV(HydroParam%eta(iElem),HydroParam%hb(iElem))
            sumH = Sum( HydroParam%H(MeshParam%Edge(:,iElem)) )
            If (sumH < 1e-3) Then
                ! The Cell is Dry. The Solution is V(eta^(n+1)) = V(eta^(n))
                HydroParam%P(iElem) = MeshParam%Area(iElem)       ! We can't allow the row to have zeros in all elements - If P = 1, the water level remains the same
            EndIf
        EndDo
        !!$OMP end parallel do
            
        res = sqrt(sum(HydroParam%f**2))      ! Residual of the Method
        !Print*, 'iNewton = ',iNewton , 'res = ',res
        If ( res < 1e-14 ) Then
            Exit
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
        If ( HydroParam%eta(iElem) - HydroParam%hb(iElem) < HydroParam%PCRI/2. ) Then
            HydroParam%eta(iElem) = HydroParam%hb(iElem) + HydroParam%PCRI/2.
        EndIf
    EndDo
    !!$OMP end parallel do

    ! 9. Compute the New Horizontal Velocity Field 
    
    !!!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iElem,iLayer,iEdge,Face,l,r,NearZero,dt)
    !Do iElem = 1, MeshParam%nElem
    !    Do iEdge = 1,4
    !        Face = MeshParam%Edge(iEdge,iElem)
    !        l = MeshParam%Left(Face) 
    !        r = MeshParam%Right(Face)
    !        Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face)
    !            If (r == 0) Then
    !                If (HydroParam%IndexInflowEdge(Face) > 0) Then        ! Boundary Condition
    !                    HydroParam%u(iLayer,Face)  = HydroParam%Fu(iLayer,Face)
    !                Else
    !                    HydroParam%u(iLayer,Face) = 0.
    !                    If (HydroParam%IndexWaterLevel(iElem)>0.and.-HydroParam%hj(Face) + HydroParam%eta(l)>HydroParam%PCRI+NearZero) Then
    !                        HydroParam%u(iLayer,Face)  = HydroParam%iAG(iLayer,Face) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(Face))*(HydroParam%etaInf(iElem) - HydroParam%eta(l))*HydroParam%iADZ(iLayer,Face) 
    !                    EndIf
    !                   ! !Vertedor -----------------------------------------
    !                   ! If ( Face == 8048.AND.l == 3777 ) Then
    !                   !     HydroParam%u(iLayer,Face)  = HydroParam%Fu(iLayer,Face)
    !                   ! Else
    !                   !     Continue
    !                   ! End If  
    !                   !!-----------------------------------------
    !                    
    !                EndIf
    !            Else
    !                ! If a face is dry, set the velocity to zero
    !                If ( Max( HydroParam%PCRI,-HydroParam%hj(Face) + HydroParam%etan(l), -HydroParam%hj(Face) + HydroParam%etan(r) ) <= HydroParam%PCRI+NearZero.or.Max( HydroParam%PCRI,-HydroParam%hj(Face) + HydroParam%eta(l), -HydroParam%hj(Face) + HydroParam%eta(r) ) <= HydroParam%PCRI+NearZero) Then
    !                    !If ( Max( HydroParam%PCRI,-HydroParam%hj(Face) + HydroParam%eta(l), -HydroParam%hj(Face) + HydroParam%eta(r) ) <= HydroParam%PCRI+NearZero) Then
    !                        HydroParam%u(iLayer,Face)  = 0.
    !                    !Else
    !                        !HydroParam%u(iLayer,Face)  = HydroParam%iAG(iLayer,Face) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(Face))*(HydroParam%eta(r) - HydroParam%eta(l))*HydroParam%iADZ(iLayer,Face)
    !                        !HydroParam%ut(iLayer,Face) = HydroParam%u(iLayer,Face)
    !                    !Endif
    !                Else
    !                    HydroParam%u(iLayer,Face)  = HydroParam%iAG(iLayer,Face) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(Face))*(HydroParam%eta(r) -HydroParam%eta(l))*HydroParam%iADZ(iLayer,Face)
    !                EndIf
    !            EndIf
    !        EndDo
    !        ! 10.1 Copy Velocities Above the Free-Surface (du/dz=0) 
    !        Do iLayer = HydroParam%CapitalM(Face) + 1, MeshParam%KMAX
    !            HydroParam%u(iLayer,Face) = 0. !u(CapitalM(Face),Face) 
    !        EndDo
    !        ! 10.2 Nullify Velocities Below the Bottom (u=0) 
    !        Do iLayer = 1, HydroParam%Smallm(Face) - 1
    !            HydroParam%u(iLayer,Face) = 0.
    !        EndDo
    !        HydroParam%Hu(Face) = Sum( HydroParam%u(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face) )/(HydroParam%CapitalM(Face)-HydroParam%Smallm(Face)+1)
    !    EndDo
    !EndDo
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
            HydroParam%H(iEdge) = Max( HydroParam%PCRI, -HydroParam%hj(iEdge) + HydroParam%eta(l) )
        Else
            HydroParam%H(iEdge) = Max( HydroParam%PCRI,-HydroParam%hj(iEdge) + HydroParam%eta(l), -HydroParam%hj(iEdge) + HydroParam%eta(r) ) !Max( PCRI,-hj(iEdge) + 0.5*(eta(l) + eta(r)) ) !
        EndIf
        ! Lower Index
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

        ! Upper Index
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
        Do iLayer = HydroParam%Smallm(iEdge) + 1, HydroParam%CapitalM(iEdge)
            HydroParam%Z(iLayer,iEdge) = MeshParam%LIMCAMAUX( MeshParam%KMax-iLayer+1) !zL + (iLayer - 1)*dz                  ! Equidistant Core Grid
        EndDo
        HydroParam%Z(HydroParam%Smallm(iEdge),iEdge)     = HydroParam%hj(iEdge)                 ! Bottom
        If (r == 0) Then
            HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%eta(l) !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
        Else
            If ( HydroParam%H(iEdge) <= HydroParam%PCRI+NearZero ) Then
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%eta(l),HydroParam%eta(r))
            Else
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%eta(l),HydroParam%eta(r)) !0.5*(HydroParam%eta(l)+HydroParam%eta(r)) !Max(HydroParam%eta(l),HydroParam%eta(r)) !0.5*(eta(l)+eta(r)) !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
            EndIf
        EndIf
        
        ! 9.3 Compute the Vertical Mesh Spacing
        Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
            HydroParam%DZj(iLayer,iEdge) = HydroParam%Z(iLayer+1,iEdge) - HydroParam%Z(iLayer,iEdge)
            HydroParam%DZj(iLayer,iEdge) = Max(HydroParam%Pcri,HydroParam%DZj(iLayer,iEdge))
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
        
        ! Lower Index
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
        Do iLayer = HydroParam%ElSmallm(iElem)+1, HydroParam%ElCapitalM(iElem)
            HydroParam%Ze(iLayer,iElem) = MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)      ! Equidistant Core Grid
        EndDo
        HydroParam%Ze(HydroParam%ElSmallm(iElem),iElem)     = HydroParam%hb(iElem)                    ! Bottom
        HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%eta(iElem) !- hb(iElem)       ! Free-Surface
        Do iLayer = 1, HydroParam%ElSmallm(iElem) - 1
            HydroParam%Ze(iLayer,iElem) = HydroParam%hb(iElem)
        EndDo
        Do iLayer = HydroParam%ElCapitalM(iElem)+2, MeshParam%KMax+1
            HydroParam%Ze(iLayer,iElem) = HydroParam%eta(iElem)
        EndDo
        
        ! 11.3 Compute the Vertical Mesh Spacing
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            HydroParam%DZi(iLayer,iElem) = HydroParam%Ze(iLayer+1,iElem) - HydroParam%Ze(iLayer,iElem)
            HydroParam%DZi(iLayer,iElem) = Max(HydroParam%Pcri,HydroParam%DZi(iLayer,iElem))
        EndDo
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            If (iLayer == HydroParam%ElSmallm(iElem)) Then
                HydroParam%Zb(iLayer,iElem) = 0.5*( HydroParam%Ze(iLayer,iElem) + HydroParam%Ze(iLayer+1,iElem) ) !0.5*( 0. + Ze(iLayer+1,iElem) )   (verificar com Rafael)               
            Else
                HydroParam%Zb(iLayer,iElem) = 0.5*( HydroParam%Ze(iLayer,iElem) + HydroParam%Ze(iLayer+1,iElem) )
            EndIf
        EndDo

        ! 11.4 Evaluate the new Vertical Velocity Field
        !HydroParam%w(HydroParam%ElSmallm(iElem),iElem)       = 0.       ! No Flux through the Bottom
        !Do iLayer = HydroParam%ElSmallm(iElem) + 1, HydroParam%ElCapitalM(iElem)
        !    SumW = 0.
        !    Do iEdge = 1,4
        !        Face = MeshParam%Edge(iEdge,iElem)
        !        SumW = SumW + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*HydroParam%DZjt(iLayer-1,Face)*(HydroParam%Theta*HydroParam%u(iLayer-1,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer-1,Face))
        !    EndDo
        !    HydroParam%w(iLayer,iElem) = HydroParam%w(iLayer-1,iElem) - (1./MeshParam%Area(iElem))*SumW         ! w(k+1/2,iElem) [1]
        !    
        !EndDo

    EndDo
    
    ! 9. Compute the New Horizontal and Vertical Velocity Field 
    HydroParam%ut = HydroParam%u
    Call uvelocity(HydroParam,MeshParam,dt)
    
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

    !print*, HydroParam%w
    !!$OMP end parallel do
    !10. Non-hydrostatic correction
    If (HydroParam%iNonHydro==0.and.MeshParam%Kmax>1) Then 
        Call NonHydroPressure(HydroParam,MeshParam,dt)
    EndIf
    
    ! 12. Compute average and tangential velocicities
    Call Velocities(HydroParam,MeshParam)
    
   ! Call ExchangeTime(HydroParam,MeshParam,MeteoParam,dt)
    
    Return
End Subroutine Hydro