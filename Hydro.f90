Subroutine Hydro(HydroParam,MeshParam,MeteoParam,dt,time,Simtime,iNewton,innerNewton)
    
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
    Use ELMConservative
    !Use AdvConservativeScheme
    
    Implicit None
    Integer:: iElem, iEdge, lEdge, iNode, iNewton, iLayer, INFO, iLayer_bar, Small, FlagLayer, innerNewton
    Integer:: llElem, rrElem, psi_flag,  indEdge
    Real:: Courant, fi_small, DZjAcum, Futn, Fvtn, dummy
    Integer:: r, l, Sig, Face, Pij, DIM, ie, j, k, iNode1,iNode2, cont
    Real:: V, DV, SumRHS, SumH, res, SumLhS,gamma, teste
    Real:: dzp, dzm, SumW, rAux, Aux, VerEddyViscUp,VerEddyViscDown, vel, e0, k0, H
    Real:: Chezy,sum1,sum2,sum3,sum4,sum0,rhoairCell
    Real:: NearZero = 1e-10
    Real:: dt, man,raioh,slope,SimTime,time
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam
    type(MeteorologicalParam) :: MeteoParam
    Real:: aTh(MeshParam%KMax), bTh(MeshParam%KMax), cTh(MeshParam%KMax), iBGhost(3,3,MeshParam%KMax,MeshParam%nElem)

    !Do iEdge = 1,MeshParam%nEdge
    !    HydroParam%epson(:,iEdge) = 0.
    !    !Sponge layer
    !    If (iEdge>2641) Then !2641 for dx=dy=0.025 !5281 for dx=dy=0.0125
    !        Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)  
    !            HydroParam%epson(iLayer,iEdge) = 0.5*(((MeshParam%EdgeBary(1,iEdge)-MeshParam%EdgeBary(1,2641))/(MeshParam%EdgeBary(1,MeshParam%nEdge)-MeshParam%EdgeBary(1,2641)))**2.)*((HydroParam%Z(HydroParam%Smallm(iEdge),iEdge)-HydroParam%Z(iLayer,iEdge))/(HydroParam%Z(HydroParam%Smallm(iEdge),iEdge)-HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge)))
    !        EndDo
    !    Else
    !         HydroParam%epson(:,iEdge) = 0.
    !    EndIf
    !EndDo
    
    !MeshParam%ei = 0.1 !e0 0.3 b01
    !MeshParam%Ksat = 0.00005!k0 0.01 b01 0.00005
    !Bench 01:    
    e0 = 0.3 !e0 0.3 b01
    !Bench 02:
    !e0 = 0.1 !e0 0.1 b02
    !!Bench 03:
    !e0 = 0.5 !0.2 !e0 0.1 b01 0.3

    HydroParam%iConv = 0

    !HydroParam%iNonHydro=0
    
    ! 0. Compute turbulence
    Call Turbulence(HydroParam,MeshParam,MeteoParam,dt)
        
    ! 1. Convective Term
    If (HydroParam%iConv == 0.or.HydroParam%iConv == 4) Then
        
        HydroParam%Fw = HydroParam%w
        HydroParam%Fu = HydroParam%u
        HydroParam%Fv = HydroParam%utang
        Call FuFv(HydroParam,MeshParam,dt)

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
    ElseIf (HydroParam%iConv == 5) Then
        
        !HydroParam%Fw = HydroParam%w
        HydroParam%Fu = HydroParam%u
        HydroParam%Fv = HydroParam%utang
        Call FuFvConservative(HydroParam,MeshParam,dt)
        
    EndIf
    
    ! 2. Getting hydrodynamic boundary condition values at current step time
    Call GetHydroBoundaryConditions(HydroParam,MeshParam,dt,time,SimTime)
    
    !if(simtime == dt*202) then
    !    continue
    !endif

    ! 3. Baroclinic pressure effect
    Call Pressure(HydroParam,MeshParam,dt)
    
    ! 4. Viscosity and Coriolis effect
    Call ExplicitTerms(HydroParam,MeshParam,dt)
    
    ! 5. Get wind velocity components
    Call WindVelocity(HydroParam,MeshParam,MeteoParam)
    
    !! 6. Vertical Water Balance (Precipitation and Evaporation)
    Call VerticalWB(HydroParam,MeshParam,MeteoParam,dt,SimTime)
    
    ! Calculate bed friction coefficient, in this point the velocity field is time tn. This coefficient is using to calculate new velocity field (in tn + 1):
    !Call BedFriction(HydroParam,MeshParam,dt)
    
    ! 7. Assemble Matrix
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam,MeteoParam) private(iEdge,iLayer,l,r,Chezy,rhoairCell,aTh,bTh,cTh,dzp,dzm,DIM,NearZero,dt)
    !If (MeshParam%iBedrock == 1) Then
    Do iEdge = 1,MeshParam%nEdge
        
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
        
        If (r == 0) Then
            rhoairCell = MeteoParam%rhoair(l)
        Else
            rhoairCell = 0.5*(MeteoParam%rhoair(l) + MeteoParam%rhoair(r))
        EndIf
        !
        ! 7.1 Get roughness 
        H = HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge) !Surface Water Height
        If (r == 0) Then
            If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
                Chezy = HydroParam%Rug(l)
            ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
                Chezy = Max(HydroParam%Pcri,H)**(1./6.)/(HydroParam%Rug(l)+NearZero)
            ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
                Chezy = 18.*log10(12.*Max(HydroParam%Pcri,H)/(HydroParam%Rug(l)/30.+NearZero))
            EndIf
            rhoairCell = MeteoParam%rhoair(l)
        Else
            If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
                Chezy = 0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))
            ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
                Chezy = Max(HydroParam%Pcri,H)**(1./6.)/(0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))+NearZero)
            ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
                Chezy = 18.*log10(12.*Max(HydroParam%Pcri,H)/(0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))/30.+NearZero))
            EndIf
            rhoairCell = 0.5*(MeteoParam%rhoair(l) + MeteoParam%rhoair(r))
        EndIf
    
        ! 7.2 Get Inflow/Outflow and Normal Depth Boundary Condition 
        ! If a face is dry, set the normal velocity to zero
        If (HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge) <= HydroParam%PCRI+NearZero) Then
            HydroParam%u(:,iEdge)  = 0.
            HydroParam%Fu(:,iEdge) = 0.
        EndIf
        
        !If(r==0)Then
        !    r = l
        !EndIf
        !If (HydroParam%IndexWaterLevelEdge(iEdge)>0.and. H >HydroParam%PCRI+NearZero) Then
        !    Futn = HydroParam%Fu(HydroParam%Smallm(iEdge),iEdge) - (dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etaInf(l) - HydroParam%eta(l)))
        !Else
        !    Futn = HydroParam%Fu(HydroParam%Smallm(iEdge),iEdge) - (dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%eta(r) - HydroParam%eta(l)))
        !EndIf
        !Fvtn = HydroParam%Fv(HydroParam%Smallm(iEdge),iEdge) - dt/MeshParam%CirDistance(iEdge)*HydroParam%g*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge)) )
        !!Fvtn = HydroParam%uxy(HydroParam%Smallm(iEdge),2,iEdge) - dt/MeshParam%CirDistance(iEdge)*HydroParam%g*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge)) )
        !

        !Call Tension(HydroParam%iWindStress,HydroParam%BottomTensionFlag,iEdge,HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),HydroParam%g,HydroParam%uxy(HydroParam%Smallm(iEdge),1,iEdge),HydroParam%uxy(HydroParam%Smallm(iEdge),2,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge),Chezy,rhoairCell,HydroParam%rho0,HydroParam%windDragConstant,HydroParam%WindVel(:,iEdge),HydroParam%WindXY(:,iEdge),HydroParam%GammaT,HydroParam%GammaB)
        !Call Tension(HydroParam%iWindStress,HydroParam%BottomTensionFlag,iEdge,HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),HydroParam%g,HydroParam%uxy(HydroParam%Smallm(iEdge),1,iEdge),HydroParam%uxy(HydroParam%Smallm(iEdge),2,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge),Chezy,rhoairCell,HydroParam%rho0,HydroParam%windDragConstant,HydroParam%WindVel(:,iEdge),HydroParam%WindXY(:,iEdge),HydroParam%GammaT,dummy)
        Call Tension(HydroParam%iWindStress,HydroParam%BottomTensionFlag,iEdge,HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),HydroParam%g,HydroParam%uxy(HydroParam%Smallm(iEdge),1,iEdge),HydroParam%uxy(HydroParam%Smallm(iEdge),2,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge),Chezy,rhoairCell,HydroParam%rho0,HydroParam%windDragConstant,HydroParam%WindVel(:,iEdge),HydroParam%WindXY(:,iEdge),HydroParam%GammaT,HydroParam%GammaB(iEdge))

!        Call Tension(HydroParam%iWindStress,HydroParam%BottomTensionFlag,iEdge,HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),HydroParam%g,Futn,Futn,HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge),Chezy,rhoairCell,HydroParam%rho0,HydroParam%windDragConstant,HydroParam%WindVel(:,iEdge),HydroParam%WindXY(:,iEdge),HydroParam%GammaT,HydroParam%GammaB)
        !r = MeshParam%Right(iEdge)
        
        !! 7.3 Get Shear stresses in the edges
        !Call Tension(HydroParam%iWindStress,HydroParam%BottomTensionFlag,iEdge,HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),HydroParam%g,HydroParam%uxy(HydroParam%Smallm(iEdge),1,iEdge),HydroParam%uxy(HydroParam%Smallm(iEdge),2,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge),Chezy,rhoairCell,HydroParam%rho0,HydroParam%windDragConstant,HydroParam%WindVel(:,iEdge),HydroParam%WindXY(:,iEdge),HydroParam%GammaT,HydroParam%GammaB)

        ! 7.4 Assemble Matrix G 
        ! Different from the article [1]
        ! Here the multiplication by DZj is done when we calculate iAG
        ! Only at the First Layer the Wind force takes place
        ! When there is no neighbour, there is no flux ( (1-Theta)*g*(dt/CirDistance(iEdge))*(eta(r) - eta(l)) == 0 )
        Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)      
            If ( iLayer == HydroParam%CapitalM(iEdge) ) Then
                If ( r == 0 .or. HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge)<= HydroParam%PCRI+NearZero) Then
                    HydroParam%Gu(iLayer,iEdge) = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)*HydroParam%Fu(iLayer,iEdge) !+ dt*GammaT*WindVel
                Else
                    HydroParam%Gu(iLayer,iEdge) = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)*(HydroParam%Fu(iLayer,iEdge) - (1-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*( HydroParam%g*(HydroParam%eta(r) - HydroParam%eta(l)) + (HydroParam%q(iLayer,r) - HydroParam%q(iLayer,l)))) + dt*HydroParam%GammaT*HydroParam%WindVel(1,iEdge)
                EndIf
            Else
                If ( r == 0 .or. HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge)<= HydroParam%PCRI+NearZero) Then
                    HydroParam%Gu(iLayer,iEdge) = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)*HydroParam%Fu(iLayer,iEdge)
                Else
                    HydroParam%Gu(iLayer,iEdge) = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)*(HydroParam%Fu(iLayer,iEdge) - (1-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*( HydroParam%g*(HydroParam%eta(r) - HydroParam%eta(l)) + (HydroParam%q(iLayer,r) - HydroParam%q(iLayer,l))))
                EndIf
            EndIf
        EndDo
        
        ! 7.5 Assemble the TriDiagonal System in Vertical Direction
        If ( HydroParam%Smallm(iEdge) == HydroParam%CapitalM(iEdge) ) Then        ! Only One Vertical Layer
            !If no neighbour or the edge's surface layer is dry:
            !If ( r == 0 .or. HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge) <= HydroParam%PCRI+NearZero) Then
            !If (HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge) <= HydroParam%PCRI+NearZero) Then
            !    HydroParam%GammaB = 0. 
            !EndIf
            !HydroParam%GammaT = 0.d0
            !DZhj is surface water thickness, in following formulation DZhj is equivalent H in [2]:
            !HydroParam%iADZ(HydroParam%Smallm(iEdge),iEdge) = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)/( HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge) + dt*(HydroParam%GammaB + HydroParam%GammaT) + NearZero)
            !HydroParam%iAG(HydroParam%Smallm(iEdge),iEdge)  = HydroParam%Gu(HydroParam%Smallm(iEdge),iEdge)/(HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge) + dt*(HydroParam%GammaB + HydroParam%GammaT) + NearZero )
            !HydroParam%DZiADZ(iEdge)             = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)*HydroParam%iADZ(HydroParam%Smallm(iEdge),iEdge) 
            !!In 2D Model iADZ = DZiA, which implies that DZiAG = iADZG:
            !HydroParam%DZiAG(iEdge)             = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)*HydroParam%iAG(HydroParam%Smallm(iEdge),iEdge)

            !DZhj is surface water thickness, in following formulation DZhj is equivalent H in [2]:
            HydroParam%iADZ(HydroParam%Smallm(iEdge),iEdge) = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)/( HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge) + dt*(HydroParam%GammaB(iEdge) + HydroParam%GammaT) + NearZero)
            HydroParam%iAG(HydroParam%Smallm(iEdge),iEdge)  = HydroParam%Gu(HydroParam%Smallm(iEdge),iEdge)/(HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge) + dt*(HydroParam%GammaB(iEdge) + HydroParam%GammaT) + NearZero )
            HydroParam%DZiADZ(iEdge)             = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)*HydroParam%iADZ(HydroParam%Smallm(iEdge),iEdge) 
            !In 2D Model iADZ = DZiA, which implies that DZiAG = iADZG:
            HydroParam%DZiAG(iEdge)             = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)*HydroParam%iAG(HydroParam%Smallm(iEdge),iEdge)
        Else
            ! 7.5.1 Assemble Matrix A
            If ( r == 0 .or. HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge)<= HydroParam%PCRI+NearZero) Then
                !HydroParam%GammaB = 0. 
                HydroParam%GammaB(iEdge) = 0.                 
            EndIf
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                If ( iLayer == HydroParam%Smallm(iEdge) ) Then
                    dzp         = 0.5*( HydroParam%DZhj(iLayer,iEdge) + HydroParam%DZhj(iLayer+1,iEdge) )       ! Dz at Upper Interface
                    aTh(iLayer) = 0. !-dt*VerEddyVisc(iLayer+1,iEdge)/dzp !-dt*nuz/dzp   !                                        ! Lower Diagonal 
                    !bTh(iLayer) = MeshParam%EdgeLength(iEdge)*(HydroParam%DZhj(iLayer,iEdge) + dt*HydroParam%VerEddyVisc(iLayer+1,iEdge)*( 1/dzp ) + HydroParam%GammaB*dt) + NearZero!DZj(iLayer,iEdge) + dt*nuz*( 1/dzp ) + GammaB*dt      ! Diagonal
                    bTh(iLayer) = MeshParam%EdgeLength(iEdge)*(HydroParam%DZhj(iLayer,iEdge) + dt*HydroParam%VerEddyVisc(iLayer+1,iEdge)*( 1/dzp ) + HydroParam%GammaB(iEdge)*dt) + NearZero!DZj(iLayer,iEdge) + dt*nuz*( 1/dzp ) + GammaB*dt      ! Diagonal
                    cTh(iLayer) = -MeshParam%EdgeLength(iEdge)*dt*HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzp !0.                                                    ! Upper Diagonal 

                Elseif ( iLayer == HydroParam%CapitalM(iEdge) ) Then
                    dzm         = 0.5*( HydroParam%DZhj(iLayer,iEdge) + HydroParam%DZhj(iLayer-1,iEdge) )       ! Dz at Lower Interface
                    aTh(iLayer) = -MeshParam%EdgeLength(iEdge)*dt*HydroParam%VerEddyVisc(iLayer,iEdge)/dzm !0.                                                    ! Lower Diagonal 
                    !bTh(iLayer) = Max(HydroParam%Pcri,MeshParam%EdgeLength(iEdge)*(HydroParam%DZhj(iLayer,iEdge) + dt*HydroParam%VerEddyVisc(iLayer,iEdge)*( 1/dzm )) + HydroParam%GammaT*dt + NearZero) !DZj(iLayer,iEdge) + dt*nuz*( 1/dzm ) + GammaT*dt      ! Diagonal
                    bTh(iLayer) = Max(NearZero,MeshParam%EdgeLength(iEdge)*(HydroParam%DZhj(iLayer,iEdge) + dt*HydroParam%VerEddyVisc(iLayer,iEdge)*( 1/dzm )) + HydroParam%GammaT*dt + NearZero) !DZj(iLayer,iEdge) + dt*nuz*( 1/dzm ) + GammaT*dt      ! Diagonal
                    cTh(iLayer) = 0. !-dt*VerEddyVisc(iLayer,iEdge)/dzm !-dt*nuz/dzm                                           ! Upper Diagonal

                Else
                    dzp         = 0.5*( HydroParam%DZhj(iLayer,iEdge) + HydroParam%DZhj(iLayer+1,iEdge) )       ! Dz at Upper Interface 
                    dzm         = 0.5*( HydroParam%DZhj(iLayer,iEdge) + HydroParam%DZhj(iLayer-1,iEdge) )       ! Dz at Lower Interface 
                    aTh(iLayer) = -MeshParam%EdgeLength(iEdge)*dt*HydroParam%VerEddyVisc(iLayer,iEdge)/dzm !-dt*VerEddyVisc(iLayer+1,iEdge)/dzp !-dt*nuz/dzm                                           ! Lower Diagonal 
                    bTh(iLayer) = MeshParam%EdgeLength(iEdge)*(HydroParam%DZhj(iLayer,iEdge) + dt*HydroParam%VerEddyVisc(iLayer,iEdge)*( 1/dzm ) +  dt*HydroParam%VerEddyVisc(iLayer+1,iEdge)*( 1/dzp )) + NearZero         ! Diagonal
                    cTh(iLayer) = -MeshParam%EdgeLength(iEdge)*dt*HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzp !-dt*VerEddyVisc(iLayer,iEdge)/dzm !-dt*nuz/dzp                                           ! Upper Diagonal 

                EndIf
            EndDo

            ! 7.5.2 Assemble Matrix: iAG, iADZ, DZiADZ, DZiAG
            DIM = size(HydroParam%DZhj(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),1)
            Call solve_tridiag(cTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),bTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),aTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),HydroParam%DZj(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),HydroParam%iADZ(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),DIM)
            Call solve_tridiag(cTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),bTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),aTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),HydroParam%Gu(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),HydroParam%iAG(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),DIM)

            HydroParam%DZiADZ(iEdge) = Dot_Product(HydroParam%DZhj(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),HydroParam%iADZ(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iedge),iEdge) )
            HydroParam%DZiAG (iEdge) = Dot_Product(HydroParam%DZhj(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),HydroParam%iAG(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iedge),iEdge) )
        EndIf
    EndDo
    !!$OMP end parallel do
    
    ! 8. Compute the New Free-Surface Elevation
    ! 8.1 Assemble the Right Hand Side (RHS)
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam,NearZero,dt,simtime) private(iElem,iEdge,SumRHS,SumLhS,gamma,Face,Pij)
    Do iElem = 1, MeshParam%nElem
        SumRHS = 0d0
        SumLhS = 0.d0
        gamma = 0.d0
        !If(MeshParam%xb(iElem) < 0.05) Then
        !    HydroParam%eta(iElem) =   5.0
        !EndIf        
        !
        Do iEdge = 1,4
            Face = MeshParam%Edge(iEdge,iElem)
            SumRHS = SumRHS + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*((1.d0-HydroParam%Theta)*(Dot_product( HydroParam%DZhj(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face),HydroParam%u(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face)) + Dot_product(HydroParam%DZsj(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face),HydroParam%us(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face))) + HydroParam%Theta*HydroParam%DZiAG(Face))
            
            !SumRHS = SumRHS + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*((1.d0-HydroParam%Theta)*(Dot_product( HydroParam%DZhj(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face),HydroParam%u(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face))) + HydroParam%Theta*HydroParam%DZiAG(Face))
            !SumRHS = SumRHS + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*((1.d0-HydroParam%Theta)*(Dot_product( HydroParam%DZhj(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face),HydroParam%u(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face)) + Dot_product(HydroParam%DZsj(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face),HydroParam%us(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face))) + HydroParam%Theta*(HydroParam%DZiAG(Face) + Dot_product(HydroParam%DZsj(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face),HydroParam%Gusub(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face))))
            ! 8.1.1 If there is a Pressure Boundary Condition

            Pij = MeshParam%Neighbor(iEdge,iElem) 
            If (Pij == 0.and.HydroParam%IndexWaterLevelEdge(Face)>0) Then   
                HydroParam%etaInfn(iElem) =  HydroParam%etaInf(iElem)
				If ((HydroParam%WaterLevel(HydroParam%IndexWaterLevelEdge(Face))-HydroParam%hj(Face))<HydroParam%PCRI  +NearZero) Then
                    !HydroParam%etaInf(iElem) = HydroParam%hj(Face) + HydroParam%PCRI/2.d0 !Rever
                    HydroParam%etaInf(iElem) = HydroParam%hj(Face)
                Else
                    HydroParam%etaInf(iElem) = HydroParam%WaterLevel(HydroParam%IndexWaterLevelEdge(Face))
        !            If(iElem > 3) Then
        !                HydroParam%etaInf(iElem) = 0.50d0
        !            Else
				    !    HydroParam%etaInf(iElem) = HydroParam%WaterLevel(HydroParam%IndexWaterLevelEdge(Face))
				    !EndIf
                EndIf
                SumLHS = SumLHS + ( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( HydroParam%etaInf(iElem) )*(HydroParam%g*HydroParam%Theta*dt*HydroParam%DZiADZ(Face) + HydroParam%DZK(Face))
                !SumLHS = SumLHS + ( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( HydroParam%etaInf(iElem) )*(HydroParam%Theta*dt*HydroParam%g*HydroParam%Theta*dt*HydroParam%DZiADZ(Face) + dt*HydroParam%DZK(Face))
                !SumLHS = SumLHS + ( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( HydroParam%etaInf(iElem) )*(HydroParam%g*HydroParam%Theta*dt*HydroParam%DZiADZ(Face) + HydroParam%Theta*HydroParam%DZK(Face))
            EndIf
        EndDo
        HydroParam%rhs(iElem) = HydroParam%Vol(iElem) - dt*SumRHS + (HydroParam%Theta*dt)*SumLHS
        !HydroParam%rhs(iElem) = HydroParam%Vol(iElem) - dt*SumRHS + SumLHS
        
    EndDo
    !!$OMP end parallel do  
    
!!!################################################################################################# INICIO ALGORITMO NEWTON NOVO #################    
    !! 8.2 Newton-Casulli-Zanolli Non-Linear System Solver Algorithm:
    !! The Volume Sistem is rewritten such that: V(eta) + T.eta = b -> V1(eta) - V2(eta) + T.eta = b
    !! The Volume function V(eta) can be defined by the difference between two nonnegative, bounded functions of eta.   
    HydroParam%etan = HydroParam%eta   
    !HydroParam%eta = Max(HydroParam%eta - HydroParam%Pcri,0.0d0)
    !HydroParam%eta = HydroParam%eta
    
    Do iNewton = 1,200
        !x.x.x Outer iteration - V2(eta) is linearized such that:
        !V1(k) - V2(k-1) - Q(k-1).[eta(k) - eta(k-1)] + T.eta(k) = rhs -> V1(k) + [T-Q(k-1)].eta(k) = rhs + V2(k-1) + Q(k-1).eta(k-1)
        
        !x.x.x Set eta(k-1) = eta(k,m):
        HydroParam%etak = HydroParam%eta
        
        !x.x.x Compute Q(k-1) and rhs(k) Matrices:
        HydroParam%P  = 0.0d0
        HydroParam%Qk = 0.d0 !always in k-1 (%etak)
        HydroParam%Ci = 0.d0
        HydroParam%Vol1 = 0.d0
        HydroParam%Vol2 = 0.d0
        
        !x.x.x Compute T.eta Matrix:
        Call MatOp(HydroParam%etak,HydroParam%Aeta,dt,HydroParam,MeshParam)       
        
        Do iElem = 1, MeshParam%nElem
    
            !x.x.x. Compute Q(k-1) and V2(k-1):
            Call MoistureContent(HydroParam%etak(iElem),0.d0,iElem,HydroParam,MeshParam)
            
            !V2 = V1 - V
            HydroParam%Vol2(iElem) = HydroParam%Vol1(iElem) - HydroParam%Vol(iElem)
            ! Q = P - A.e.S, Q >=0
            !HydroParam%Qk(iElem) = Max(0.d0, HydroParam%P(iElem) - HydroParam%Vol(iElem)/HydroParam%eta(iElem))
            !!HydroParam%Qk(iElem) = Max(0.d0, HydroParam%P(iElem) - MeshParam%Area(iElem)*MeshParam%ei(HydroParam%ElCapitalMs(iElem),iElem))
            !If(HydroParam%etak(iElem) - HydroParam%hb(iElem) < 0.d0) Then
            !    HydroParam%Qk(iElem) = HydroParam%Vol2(iElem)/HydroParam%etak(iElem)
            !Else
            !    HydroParam%Qk(iElem) = Max(0.d0, HydroParam%P(iElem) - HydroParam%Ci(iElem))
            !EndIf
            If (MeshParam%iBedrock == 0) Then
                HydroParam%Qk(iElem) = 0.d0
            Else
                HydroParam%Qk(iElem) = Max(0.d0, HydroParam%P(iElem) - HydroParam%Ci(iElem))
            EndIf
            !HydroParam%Qk(iElem) = Max(0.d0, HydroParam%P(iElem) - MeshParam%Area(iElem)*MeshParam%ei(HydroParam%ElCapitalMs(iElem),iElem))
            !x.x.x. Compute d(k) == rhs(k):
            !d(k) == rhs(k) = rhs + V2(k-1) - Q(k-1).eta(k-1)
            HydroParam%d(iElem) = HydroParam%rhs(iElem) + HydroParam%Vol2(iElem) - HydroParam%Qk(iElem)*V(HydroParam%etak(iElem),HydroParam%sb(iElem))
            
            !!x.x.x. outer residual
            !HydroParam%F(iElem) = HydroParam%Vol1(iElem) - HydroParam%Vol2(iElem) + HydroParam%Aeta(iElem) - HydroParam%rhs(iElem)
        EndDo    
        
        !res = sqrt(sum(HydroParam%F**2))      ! Residual of the Method
        !!Print*, 'iNewton = ',iNewton , 'res = ',res
        !If ( res < 1e-8 ) Then
        !    !x.x.x Set eta(k) = eta(k,m)
        !    continue
        !    exit
        !EndIf
        
        HydroParam%eta = HydroParam%etak
        !HydroParam%etam = HydroParam%etak
        Do innerNewton = 1,200
        !x.x.x Inner iteration - V1(eta) is linearized such that:
        !V1(m-1) + P(k,m-1).[eta(k,m) - eta(k,m-1)] + [T-Q(k-1)].eta(k,m) = rhs(k) -> [T + P(m,k-1) - Q(k-1)].eta(k,m) = rhs(k) - V1(m-1) + P(k,m-1).eta(k,m-1)
        !From here, we get the system: [T + P(k,m-1) - Q(k-1)].eta(k,m) = rhs(k,m-1) 
             
            !x.x.x Compute (T + P(k,m-1) - Q(k)).eta(k,m) Matrix:
            !(T + P(k,m-1) - Q(k)).eta(k,m) == A.eta -> A == Jacobian Matrix
            Call MatOp(HydroParam%eta,HydroParam%Aeta,dt,HydroParam,MeshParam)
            
            !x.x.x Set A.eta - rhs(k,m) = 0 to conjugate gradient method:
            !HydroParam%Aeta =  HydroParam%Aeta - HydroParam%d + HydroParam%Vol1 - HydroParam%P*V(HydroParam%eta,HydroParam%sb)
            Do iElem = 1,MeshParam%nElem
                HydroParam%Aeta(iElem) =  HydroParam%Aeta(iElem) - HydroParam%d(iElem) + HydroParam%Vol1(iElem) - HydroParam%P(iElem)*V(HydroParam%eta(iElem),HydroParam%sb(iElem))    
            EndDo
            !x.x.x Compute the New Free-Surface Elevation eta(k,m):
            Call CGOp(HydroParam%Aeta,HydroParam%Deta,dt,HydroParam,MeshParam)
                   
            !$OMP parallel do default(none) shared(HydroParam,MeshParam) private(iElem)
            Do iElem = 1, MeshParam%nElem
                HydroParam%F(iElem) = HydroParam%P(iElem)*HydroParam%Deta(iElem) - HydroParam%Vol1(iElem)
                !x.x.x New Free Surface Elevation:
                !HydroParam%eta(iElem) = max(0.d0,HydroParam%eta(iElem) - HydroParam%Deta(iElem))
                HydroParam%eta(iElem) = HydroParam%eta(iElem) - HydroParam%Deta(iElem)
                Call MoistureContent(HydroParam%eta(iElem),0.d0,iElem,HydroParam,MeshParam)
                HydroParam%F(iElem) = HydroParam%F(iElem) + HydroParam%Vol1(iElem)
            EndDo
            !$OMP end parallel do  
            
            res = sqrt(sum(HydroParam%F**2))      ! Residual of the Method
            !Print*, 'iNewton = ',iNewton , 'res = ',res
            If ( res < 1e-5 ) Then
                !x.x.x Set eta(k) = eta(k,m)
                continue
                exit
            EndIf
            
        EndDo
        
        HydroParam%F = HydroParam%Qk*(HydroParam%etak-HydroParam%eta) - HydroParam%Vol2 + (HydroParam%Vol1 - HydroParam%Vol)
        res = sqrt(sum(HydroParam%F**2))      ! Residual of the Method
        !Print*, 'iNewton = ',iNewton , 'res = ',res
        If ( res < 1e-5 ) Then
            !x.x.x Set eta(k) = eta(k,m)
            continue
            exit
        EndIf
                
    EndDo
    !
    !Do iNewton = 1,200
    !    !x.x.x Outer iteration - V2(eta) is linearized such that:
    !    !V1(k) - V2(k-1) - Q(k-1).[eta(k) - eta(k-1)] + T.eta(k) = rhs -> V1(k) + [T-Q(k-1)].eta(k) = rhs + V2(k-1) + Q(k-1).eta(k-1)
    !    
    !    !x.x.x Set eta(k-1) = eta(k,m):
    !    HydroParam%etak = HydroParam%eta
    !    
    !    !x.x.x Compute Q(k-1) and rhs(k) Matrices:
    !    HydroParam%P  = 0.0d0
    !    HydroParam%Qk = 0.d0 !always in k-1 (%etak)
    !    HydroParam%Ci = 0.d0
    !    HydroParam%Vol1 = 0.d0
    !    HydroParam%Vol2 = 0.d0
    !    
    !    !x.x.x Compute T.eta Matrix:
    !    Call MatOp(HydroParam%etak,HydroParam%Aeta,dt,HydroParam,MeshParam)       
    !    
    !    Do iElem = 1, MeshParam%nElem
    !
    !        !x.x.x. Compute Q(k-1) and V2(k-1):
    !        Call MoistureContent(HydroParam%etak(iElem),0.d0,iElem,HydroParam,MeshParam)
    !        
    !        !V2 = V1 - V
    !        HydroParam%Vol2(iElem) = HydroParam%Vol1(iElem) - HydroParam%Vol(iElem)
    !        
    !        HydroParam%Qk(iElem) = Max(0.d0, HydroParam%P(iElem) - HydroParam%Ci(iElem))
    !
    !        !x.x.x. Compute d(k) == rhs(k):
    !        !d(k) == rhs(k) = rhs + V2(k-1) - Q(k-1).eta(k-1)
    !        HydroParam%d(iElem) = HydroParam%rhs(iElem) + HydroParam%Vol2(iElem) - HydroParam%Qk(iElem)*V(HydroParam%etak(iElem),HydroParam%sb(iElem))
    !        
    !        !!x.x.x. outer residual
    !        !HydroParam%F(iElem) = HydroParam%Vol1(iElem) - HydroParam%Vol2(iElem) + HydroParam%Aeta(iElem) - HydroParam%rhs(iElem)
    !    EndDo    
    !    
    !    HydroParam%eta = HydroParam%etak
    !    
    !    !HydroParam%etam = HydroParam%etak
    !    Do innerNewton = 1,200
    !    !x.x.x Inner iteration - V1(eta) is linearized such that:
    !    !V1(m-1) + P(k,m-1).[eta(k,m) - eta(k,m-1)] + [T-Q(k-1)].eta(k,m) = rhs(k) -> [T + P(m,k-1) - Q(k-1)].eta(k,m) = rhs(k) - V1(m-1) + P(k,m-1).eta(k,m-1)
    !    !From here, we get the system: [T + P(k,m-1) - Q(k-1)].eta(k,m) = rhs(k,m-1) 
    !         
    !        !x.x.x Compute (T + P(k,m-1) - Q(k)).eta(k,m) Matrix:
    !        !(T + P(k,m-1) - Q(k)).eta(k,m) == A.eta -> A == Jacobian Matrix
    !        HydroParam%P  = 0.0d0
    !        HydroParam%Aeta = 0.0d0
    !        Call MatOp(HydroParam%eta,HydroParam%Aeta,dt,HydroParam,MeshParam)
    !        
    !        !x.x.x Set A.eta - rhs(k,m) = 0 to conjugate gradient method:
    !        !HydroParam%Aeta =  HydroParam%Aeta - HydroParam%d + HydroParam%Vol1 - HydroParam%P*V(HydroParam%eta,HydroParam%sb)
    !        !Do iElem = 1,MeshParam%nElem
    !        !    HydroParam%Aeta(iElem) =  HydroParam%Aeta(iElem) - HydroParam%d(iElem) + HydroParam%Vol1(iElem) - HydroParam%P(iElem)*V(HydroParam%eta(iElem),HydroParam%sb(iElem))    
    !        !EndDo
    !
    !        !!x.x.x Compute the New Free-Surface Elevation eta(k,m):
    !        !Call CGOp(HydroParam%Aeta,HydroParam%Deta,dt,HydroParam,MeshParam)
    !               
    !        !!$OMP parallel do default(none) shared(HydroParam,MeshParam) private(iElem)
    !        !Do iElem = 1, MeshParam%nElem
    !        !    HydroParam%F(iElem) = HydroParam%P(iElem)*HydroParam%Deta(iElem) - HydroParam%Vol1(iElem)
    !        !    !x.x.x New Free Surface Elevation:
    !        !    HydroParam%eta(iElem) = HydroParam%eta(iElem) - HydroParam%Deta(iElem)    
    !        !    Call MoistureContent(HydroParam%eta(iElem),0.d0,iElem,HydroParam,MeshParam)
    !        !    HydroParam%F(iElem) = HydroParam%F(iElem) + HydroParam%Vol1(iElem)
    !        !EndDo
    !        !!$OMP end parallel do  
    !        !
    !        !res = sqrt(sum(HydroParam%F**2))      ! Residual of the Method
    !        !!Print*, 'iNewton = ',iNewton , 'res = ',res
    !        !If ( res < 1e-5 ) Then
    !        !    !x.x.x Set eta(k) = eta(k,m)
    !        !    continue
    !        !    exit
    !        !EndIf
    !        
    !        Do iElem = 1,MeshParam%nElem
    !            Call MoistureContent(HydroParam%eta(iElem),0.d0,iElem,HydroParam,MeshParam)
    !            HydroParam%Aeta(iElem) =  HydroParam%Vol1(iElem) + HydroParam%Aeta(iElem) - HydroParam%d(iElem)
    !        EndDo
    !        
    !        res = sqrt(sum(HydroParam%Aeta**2))      ! Residual of the Method
    !        !Print*, 'iNewton = ',iNewton , 'res = ',res
    !        If ( res < 1e-5 ) Then
    !            !x.x.x Set eta(k) = eta(k,m)
    !            continue
    !            exit
    !        EndIf
    !    
    !        !x.x.x Compute the New Free-Surface Elevation eta(k,m):
    !        Call CGOp(HydroParam%Aeta,HydroParam%Deta,dt,HydroParam,MeshParam)            
    !        
    !        !$OMP parallel do default(none) shared(HydroParam,MeshParam) private(iElem)
    !        Do iElem = 1, MeshParam%nElem
    !            !x.x.x New Free Surface Elevation:
    !            HydroParam%eta(iElem) = HydroParam%eta(iElem) - HydroParam%Deta(iElem)     
    !        EndDo
    !        !$OMP end parallel do  
    !        
    !    EndDo
    !    
    !    HydroParam%F = HydroParam%Qk*(HydroParam%etak-HydroParam%eta) - HydroParam%Vol2 + (HydroParam%Vol1 - HydroParam%Vol)
    !    res = sqrt(sum(HydroParam%F**2))      ! Residual of the Method
    !    !Print*, 'iNewton = ',iNewton , 'res = ',res
    !    If ( res < 1e-5 ) Then
    !        !x.x.x Set eta(k) = eta(k,m)
    !        continue
    !        exit
    !    EndIf
    !            
    !EndDo    !
    !Do iNewton = 1,200
    !    !x.x.x Outer iteration - V2(eta) is linearized such that:
    !    !V1(k) - V2(k-1) - Q(k-1).[eta(k) - eta(k-1)] + T.eta(k) = rhs -> V1(k) + [T-Q(k-1)].eta(k) = rhs + V2(k-1) + Q(k-1).eta(k-1)
    !    
    !    !x.x.x Set eta(k-1) = eta(k,m):
    !    HydroParam%etak = HydroParam%eta
    !    
    !    !x.x.x Compute Q(k-1) and rhs(k) Matrices:
    !    HydroParam%P  = 0.0d0
    !    HydroParam%Qk = 0.d0 !always in k-1 (%etak)
    !    HydroParam%Ci = 0.d0
    !    HydroParam%Vol1 = 0.d0
    !    HydroParam%Vol2 = 0.d0
    !    
    !    !x.x.x Compute T.eta Matrix:
    !    Call MatOp(HydroParam%etak,HydroParam%Aeta,dt,HydroParam,MeshParam)       
    !    
    !    Do iElem = 1, MeshParam%nElem
    !
    !        !x.x.x. Compute Q(k-1) and V2(k-1):
    !        Call MoistureContent(HydroParam%etak(iElem),0.d0,iElem,HydroParam,MeshParam)
    !        
    !        !V2 = V1 - V
    !        HydroParam%Vol2(iElem) = HydroParam%Vol1(iElem) - HydroParam%Vol(iElem)
    !        
    !        HydroParam%Qk(iElem) = Max(0.d0, HydroParam%P(iElem) - HydroParam%Ci(iElem))
    !
    !        !x.x.x. Compute d(k) == rhs(k):
    !        !d(k) == rhs(k) = rhs + V2(k-1) - Q(k-1).eta(k-1)
    !        HydroParam%d(iElem) = HydroParam%rhs(iElem) + HydroParam%Vol2(iElem) - HydroParam%Qk(iElem)*V(HydroParam%etak(iElem),HydroParam%sb(iElem))
    !        
    !        !!x.x.x. outer residual
    !        !HydroParam%F(iElem) = HydroParam%Vol1(iElem) - HydroParam%Vol2(iElem) + HydroParam%Aeta(iElem) - HydroParam%rhs(iElem)
    !    EndDo    
    !    
    !    HydroParam%eta = HydroParam%etak
    !    
    !    !HydroParam%etam = HydroParam%etak
    !    Do innerNewton = 1,200
    !    !x.x.x Inner iteration - V1(eta) is linearized such that:
    !    !V1(m-1) + P(k,m-1).[eta(k,m) - eta(k,m-1)] + [T-Q(k-1)].eta(k,m) = rhs(k) -> [T + P(m,k-1) - Q(k-1)].eta(k,m) = rhs(k) - V1(m-1) + P(k,m-1).eta(k,m-1)
    !    !From here, we get the system: [T + P(k,m-1) - Q(k-1)].eta(k,m) = rhs(k,m-1) 
    !         
    !        !x.x.x Compute (T + P(k,m-1) - Q(k)).eta(k,m) Matrix:
    !        !(T + P(k,m-1) - Q(k)).eta(k,m) == A.eta -> A == Jacobian Matrix
    !        HydroParam%P  = 0.0d0
    !        HydroParam%Aeta = 0.0d0
    !        Call MatOp(HydroParam%eta,HydroParam%Aeta,dt,HydroParam,MeshParam)
    !        
    !        !x.x.x Set A.eta - rhs(k,m) = 0 to conjugate gradient method:
    !        !HydroParam%Aeta =  HydroParam%Aeta - HydroParam%d + HydroParam%Vol1 - HydroParam%P*V(HydroParam%eta,HydroParam%sb)
    !        !Do iElem = 1,MeshParam%nElem
    !        !    HydroParam%Aeta(iElem) =  HydroParam%Aeta(iElem) - HydroParam%d(iElem) + HydroParam%Vol1(iElem) - HydroParam%P(iElem)*V(HydroParam%eta(iElem),HydroParam%sb(iElem))    
    !        !EndDo
    !
    !        !!x.x.x Compute the New Free-Surface Elevation eta(k,m):
    !        !Call CGOp(HydroParam%Aeta,HydroParam%Deta,dt,HydroParam,MeshParam)
    !               
    !        !!$OMP parallel do default(none) shared(HydroParam,MeshParam) private(iElem)
    !        !Do iElem = 1, MeshParam%nElem
    !        !    HydroParam%F(iElem) = HydroParam%P(iElem)*HydroParam%Deta(iElem) - HydroParam%Vol1(iElem)
    !        !    !x.x.x New Free Surface Elevation:
    !        !    HydroParam%eta(iElem) = HydroParam%eta(iElem) - HydroParam%Deta(iElem)    
    !        !    Call MoistureContent(HydroParam%eta(iElem),0.d0,iElem,HydroParam,MeshParam)
    !        !    HydroParam%F(iElem) = HydroParam%F(iElem) + HydroParam%Vol1(iElem)
    !        !EndDo
    !        !!$OMP end parallel do  
    !        !
    !        !res = sqrt(sum(HydroParam%F**2))      ! Residual of the Method
    !        !!Print*, 'iNewton = ',iNewton , 'res = ',res
    !        !If ( res < 1e-5 ) Then
    !        !    !x.x.x Set eta(k) = eta(k,m)
    !        !    continue
    !        !    exit
    !        !EndIf
    !        
    !        Do iElem = 1,MeshParam%nElem
    !            Call MoistureContent(HydroParam%eta(iElem),0.d0,iElem,HydroParam,MeshParam)
    !            HydroParam%Aeta(iElem) =  HydroParam%Vol1(iElem) + HydroParam%Aeta(iElem) - HydroParam%d(iElem)
    !        EndDo
    !        
    !        res = sqrt(sum(HydroParam%Aeta**2))      ! Residual of the Method
    !        !Print*, 'iNewton = ',iNewton , 'res = ',res
    !        If ( res < 1e-5 ) Then
    !            !x.x.x Set eta(k) = eta(k,m)
    !            continue
    !            exit
    !        EndIf
    !    
    !        !x.x.x Compute the New Free-Surface Elevation eta(k,m):
    !        Call CGOp(HydroParam%Aeta,HydroParam%Deta,dt,HydroParam,MeshParam)            
    !        
    !        !$OMP parallel do default(none) shared(HydroParam,MeshParam) private(iElem)
    !        Do iElem = 1, MeshParam%nElem
    !            !x.x.x New Free Surface Elevation:
    !            HydroParam%eta(iElem) = HydroParam%eta(iElem) - HydroParam%Deta(iElem)     
    !        EndDo
    !        !$OMP end parallel do  
    !        
    !    EndDo
    !    
    !    HydroParam%F = HydroParam%Qk*(HydroParam%etak-HydroParam%eta) - HydroParam%Vol2 + (HydroParam%Vol1 - HydroParam%Vol)
    !    res = sqrt(sum(HydroParam%F**2))      ! Residual of the Method
    !    !Print*, 'iNewton = ',iNewton , 'res = ',res
    !    If ( res < 1e-5 ) Then
    !        !x.x.x Set eta(k) = eta(k,m)
    !        continue
    !        exit
    !    EndIf
    !            
    !EndDo


!!!################################################################################################# FIM ALGORITMO NOVO NEWTON #####################
!!
!!    ! 8!.2 Newton Loop for Non-Linear Wet- and Dry-ing Algorithm [2]
!    HydroParam%etan = HydroParam%eta   
!    HydroParam%Qk = 0.d0
!    !HydroParam%etak = HydroParam%eta + HydroParam%etaplus   
!    Do iNewton = 1,200
!        ! 8.2.1 Assemble Matrices F (System Matrix - Newton Method) and P (Derivative Matrix)
!        HydroParam%P = 0.
!        HydroParam%F = 0.
!        Call MatOp(HydroParam%eta,HydroParam%Aeta,dt,HydroParam,MeshParam)
!        !Call Volume(HydroParam,MeshParam)
!        !$OMP parallel do default(none) shared(HydroParam,MeshParam) private(iElem)      
!        Do iElem = 1, MeshParam%nElem
!            !8.2.2 Volume with Eta at tn+1:
!            Call MoistureContent(HydroParam%eta(iElem),0.d0,iElem,HydroParam,MeshParam)            
!    
!            !8.2.3 Compute Newton Method's Residue for each iElem:
!            !In this point, MatOp Output = (T-Matrix)Eta
!            !We want F = 0 condition satisfied: V(nt+1) + T*n(t+1) = rhs :: F  = V(n(t+1)) + T*n(t+1) - rhs -> 0
!            !HydroParam%Vol(iElem) = MeshParam%Area(iElem)*V(HydroParam%eta(iElem),0.d0)
!            
!            HydroParam%F(iElem) = HydroParam%Vol(iElem) + HydroParam%Aeta(iElem) - HydroParam%rhs(iElem)
!            
!            !8.2.4 Fill P values in dry cell for CGOp computation:
!            !HydroP = A*dn/dn, if cell is dry dn/dn = 0, else dn/dn = 1
!            !HydroParam%P(iElem) = MeshParam%Area(iElem)
!            !if (HydroParam%Vol(iElem) > 0.d0) Then
!            !    HydroParam%P(iElem) = HydroParam%Vol(iElem)/HydroParam%eta(iElem)
!            !endif
!    
!        EndDo   
!        
!        !$OMP end parallel do
!        res = sqrt(sum(HydroParam%F**2))      ! Residual of the Method
!        !Print*, 'iNewton = ',iNewton , 'res = ',res
!        If ( res < 1e-5 ) Then
!            continue
!            exit
!        EndIf
!    
!        !! 8.2.5 Compute the New Free-Surface Elevation
!        !!CGOp is used to minimize F value :: F = V(n(t+1)) + T*n(t+1) - rhs == 0
!        Call CGOp(HydroParam%F,HydroParam%Deta,dt,HydroParam,MeshParam)
!    
!        !$OMP parallel do default(none) shared(HydroParam,MeshParam) private(iElem)
!        Do iElem = 1, MeshParam%nElem
!            HydroParam%eta(iElem) = HydroParam%eta(iElem) - HydroParam%Deta(iElem)
!        EndDo
!        !$OMP end parallel do
!    EndDo
!    
!    
    ! 9. Water Level Corrective Step
    
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam,NearZero) private(iElem)
    !Do iElem = 1,MeshParam%nElem
    !    If ( HydroParam%eta(iElem) - HydroParam%sb(iElem) <= HydroParam%PCRI + NearZero ) Then
    !        HydroParam%eta(iElem) = HydroParam%sb(iElem) + HydroParam%PCRI/2.0d0
    !    EndIf
    !EndDo
    
    HydroParam%eta = Max(HydroParam%eta + HydroParam%sb, HydroParam%sb)
    !!$OMP end parallel do
    
    ! 9.1 Compute nodal elevations for tang. vel.
    HydroParam%petan=HydroParam%peta !store for transport eqs. 
    !$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iNode,sum1,sum0,ie,j)
    !Compute nodal elevations for tang. vel. MAX VALUE
    !Do iNode=1,MeshParam%nNode
    !    sum1 = 0.d0 !sum of elemental elevations
    !    sum0 = 0.d0 !sum of areas
    !    do j=1,MeshParam%nVertexElem(iNode)
    !        ie=MeshParam%VertexElem(j,iNode)
    !        sum1 = sum1 + MeshParam%Area(ie)*HydroParam%eta(ie)
    !        sum0 = sum0 + MeshParam%Area(ie)
    !    Enddo !j=1,nne(i)
    !    HydroParam%peta(iNode) = sum1/sum0
    !EndDo !i=1,np
    
    Do iNode=1,MeshParam%nNode  
        HydroParam%peta(iNode) = HydroParam%eta(MeshParam%VertexElem(1,iNode))
        do j=1,MeshParam%nVertexElem(iNode)-1
            ie=MeshParam%VertexElem(j+1,iNode)
            ! Nodal water elevation is the maximum surround value:
            HydroParam%peta(iNode) = Max(HydroParam%peta(iNode),HydroParam%eta(ie))
        Enddo !j=1,nne(i)     
    EndDo !i=1,np   
    !$OMP end parallel do
    
    ! 10. Updating Vertical Mesh Spacing
    HydroParam%DZjt = HydroParam%DZj
    HydroParam%DZhjt = HydroParam%DZhj
    HydroParam%DZsjt = HydroParam%DZsj 

    ! 11. Compute the New Horizontal Normal/Tangetial Velocity Field
    HydroParam%ut = HydroParam%u
    HydroParam%umt = HydroParam%um
    Call uvelocity(HydroParam,MeshParam,MeteoParam,dt)
    Call utangvelocity(HydroParam,MeshParam,MeteoParam,dt) 
 
    ! 12. Updating free-surface water level and  Vertical Mesh Spacing
        
    ! xx. Compute the New Vertical Velocity Field
    ! w(Smallm-1,:) = 0. -> No flux through the Bottom
    ! w(CapitalM,:) = 0. -> No flux through the Surface
    ! The velocity is defined in the barycenter of the Top Face of each Cell
    ! *Obs: The outward direction is from Bottom to Top.
    
    ! xx. Vertical Water Balance (Precipitation and Evaporation)
    Call VerticalWB(HydroParam,MeshParam,MeteoParam,dt,SimTime+dt)
    
    HydroParam%DZit = HydroParam%DZi
    HydroParam%DZhit = HydroParam%DZhi
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam,e0) private(iLayer,iElem,iEdge,Face,SumW,NearZero)
    Do iElem = 1, MeshParam%nElem
       
        ! 12.1 Define the range between Bottom and Top Layers. (Tricky part - Not explained in [1])
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
            
        
        !CAYO
        ! 4.1.2 Update the Element Vertical Spacing
        !If ( HydroParam%eta(iElem) - HydroParam%hb(iElem) <= HydroParam%PCRI+NearZero ) Then
        If ( HydroParam%eta(iElem) - HydroParam%hb(iElem) <= HydroParam%PCRI+NearZero ) Then
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
            !HydroParam%DZi(iLayer,iElem) = Max(HydroParam%Pcri,HydroParam%DZi(iLayer,iElem))
            HydroParam%DZi(iLayer,iElem) = Max(0.0d0,HydroParam%DZi(iLayer,iElem))
        EndDo
        
        !xx. Porosity and and DZhi/DZsi
        HydroParam%DZhi(:,iElem) = HydroParam%DZi(:,iElem)
		If (MeshParam%iBedrock == 1) Then
            Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem)
                If (HydroParam%H(iEdge) > HydroParam%Pcri) Then

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
                EndIf
            EndDo
        EndIf
        Call MoistureContent(HydroParam%eta(iElem),HydroParam%etaplus(iElem),iElem,HydroParam,MeshParam)
    EndDo
    !!$OMP end parallel do    
   
    
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam,k0) private(iLayer,Small,iEdge,l,r,NearZero)
    Do iEdge = 1,MeshParam%nEdge

        ! 1.1 Compute Water Depth
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
        
        If (r==0) Then
            r = l
        EndIf

        !HydroParam%H(iEdge) = Max( HydroParam%PCRI,-HydroParam%sj(iEdge) + HydroParam%eta(l), -HydroParam%sj(iEdge) + HydroParam%eta(r) )!CAYO
        If (Max( 0.d0,-HydroParam%sj(iEdge) + HydroParam%eta(l), -HydroParam%sj(iEdge) + HydroParam%eta(r) ) <= HydroParam%Pcri) Then
            HydroParam%H(iEdge) = 0.d0
        Else
            HydroParam%H(iEdge) = Max( 0.d0,-HydroParam%sj(iEdge) + HydroParam%eta(l), -HydroParam%sj(iEdge) + HydroParam%eta(r) )!CAY
            !HydroParam%H(iEdge) = Max( 0.d0,-HydroParam%sj(iEdge) + 0.5*(HydroParam%eta(l) + HydroParam%eta(r)))!CAY
        EndIf
        
        If (HydroParam%IndexWaterLevelEdge(iEdge) > 0 .and. HydroParam%H(iEdge) > HydroParam%Pcri) Then 
            HydroParam%H(iEdge) = Max( HydroParam%H(iEdge), HydroParam%WaterLevel(HydroParam%IndexWaterLevelEdge(iEdge)))
        endif
        !
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
                If (Max(HydroParam%PCRI, HydroParam%eta(l), HydroParam%eta(r))>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1)) then
                    HydroParam%CapitalM(iEdge) = Max(iLayer,HydroParam%CapitalM(iEdge))
                    exit
                EndIf
            Else
                If (Max(HydroParam%PCRI, HydroParam%eta(l), HydroParam%eta(r))>=MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer+1).and.Max(HydroParam%PCRI, HydroParam%eta(l), HydroParam%eta(r))<MeshParam%LIMCAMAUX(MeshParam%KMax-iLayer)) then
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
            If ( HydroParam%eta(l) - HydroParam%hj(iEdge) <= HydroParam%PCRI+NearZero) Then
                !HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%hj(iEdge) + HydroParam%PCRI
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%hj(iEdge)
            Else
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%eta(l)
            EndIf
        Else
            If ( Max(HydroParam%PCRI, HydroParam%eta(l), HydroParam%eta(r)) - HydroParam%hj(iEdge) <= HydroParam%PCRI+NearZero ) Then !( HydroParam%H(iEdge) <= HydroParam%PCRI/2 + NearZero ) Then
                !HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%eta(l),HydroParam%eta(r))
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%hb(l),HydroParam%hb(r))                
            Else
                !HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = (0.5d0)*(HydroParam%eta(l)+HydroParam%eta(r)) !Max(HydroParam%eta(l),HydroParam%eta(r)) !H(iEdge) + hj(iEdge)       ! Free-Surface (verificar com Rafael)
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = Max(HydroParam%eta(l),HydroParam%eta(r))             
            EndIf
        EndIf        
                
        ! 9.3 Compute the Vertical Mesh Spacing
        Do iLayer = HydroParam%Smallms(iEdge), HydroParam%CapitalM(iEdge)
            HydroParam%DZj(iLayer,iEdge) = HydroParam%Z(iLayer+1,iEdge) - HydroParam%Z(iLayer,iEdge)
            !HydroParam%DZj(iLayer,iEdge) = Max(HydroParam%Pcri,HydroParam%DZj(iLayer,iEdge))
            HydroParam%DZj(iLayer,iEdge) = Max(0.0d0,HydroParam%DZj(iLayer,iEdge))
        EndDo        
        
        !9.4 Compute Subsurface Parameters:
		HydroParam%DZhj(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge) = HydroParam%DZj(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge)
        HydroParam%DZsj(:,iEdge) = 0.d0
        MeshParam%Kj(:,iEdge) = 0.d0
        HydroParam%DZK(iEdge) = 0.d0
        If (MeshParam%iBedrock == 1) Then
            Do iLayer = HydroParam%Smallms(iEdge), HydroParam%CapitalM(iEdge) ! 
                !xx. Layers thickness:
                If (iLayer < HydroParam%Smallm(iEdge)) Then
                    HydroParam%DZhj(iLayer,iEdge) = 0.
                    HydroParam%DZsj(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)
                ElseIf (iLayer > HydroParam%Smallm(iEdge)) Then
                    HydroParam%DZhj(iLayer,iEdge) = HydroParam%DZj(iLayer,iEdge)
                    HydroParam%DZsj(iLayer,iEdge) = 0.
                    continue
                Else
                    If ( HydroParam%Z(iLayer+1,iEdge) > HydroParam%hj(iEdge) .and.  Max(HydroParam%PCRI, HydroParam%eta(l), HydroParam%eta(r)) - HydroParam%hj(iEdge) > 0.d0) Then
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
                        !     EndIf   
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
                
                If(HydroParam%DZhj(iLayer,iEdge) + HydroParam%DZsj(iLayer,iEdge) /= HydroParam%H(iEdge)) Then
                    continue
                EndIF
                
                If(HydroParam%DZhj(iLayer,iEdge) <= HydroParam%Pcri) Then
                    HydroParam%DZhj(iLayer,iEdge) = 0.d0  
                EndIf
                If(HydroParam%DZsj(iLayer,iEdge) <= HydroParam%Pcri) Then
                    !HydroParam%DZhj(iLayer,iEdge) = HydroParam%DZhj(iLayer,iEdge) + HydroParam%DZsj(iLayer,iEdge)
                    HydroParam%DZsj(iLayer,iEdge) = 0.d0  
                EndIf                
                
                !xx. Hydrauli   c Conductivity:
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
                
                !xx. Subsurface Wet Area:
                HydroParam%DZK(iEdge) = HydroParam%DZK(iEdge) + HydroParam%DZsj(iLayer,iEdge)*MeshParam%Kj(iLayer,iEdge) !Sediment Layer
            EndDo
		EndIf
    EndDo
    !!$OMP end parallel do 
    
    ! 13. Vertical Velocity Field:
    HydroParam%wt = HydroParam%w
    !Call wvelocity(HydroParam,MeshParam,dt)
    
    ! 13.1 Vertical Velocity Field to Surface Flow:
    Do iElem = 1, MeshParam%nElem
        ! 13.1.1 Evaluate the new Vertical Velocity Field
        HydroParam%w(HydroParam%ElSmallm(iElem),iElem)       = 0.       ! No Flux through the Bottom
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) 
            SumW = 0.
            Do iEdge = 1,4
                Face = MeshParam%Edge(iEdge,iElem)
                SumW = SumW + (Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*HydroParam%DZhjt(iLayer,Face)*HydroParam%u(iLayer,Face))
            EndDo
            HydroParam%w(iLayer+1,iElem) = HydroParam%w(iLayer,iElem) - SumW/MeshParam%Area(iElem)         ! w(k+1/2,iElem) [1]        
        EndDo
    EndDo

    HydroParam%wmt = HydroParam%wm
    ! 13.2 Averaged Vertical Velocity Field:
    If (MeshParam%iBedrock == 1) Then
        Do iElem = 1, MeshParam%nElem
            !13.2.1 Evaluate the Averaged Vertical Velocity Field
            HydroParam%wm(HydroParam%ElSmallms(iElem),iElem)       = 0.       ! No Flux through the Bottom
            Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem) 
                SumW = 0.
                Do iEdge = 1,4
                    Face = MeshParam%Edge(iEdge,iElem)
                    HydroParam%um(iLayer,Face) = (HydroParam%u(iLayer,Face)*HydroParam%DZhjt(iLayer,Face) + HydroParam%us(iLayer,Face)*HydroParam%DZsjt(iLayer,Face))/(HydroParam%DZsjt(iLayer,Face) + HydroParam%DZhjt(iLayer,Face))                    
                    SumW = SumW + (Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%Edgelength(Face)*(HydroParam%DZsjt(iLayer,Face) + HydroParam%DZhjt(iLayer,Face))*HydroParam%um(iLayer,Face))      
                EndDo
                HydroParam%wm(iLayer+1,iElem) = HydroParam%wm(iLayer,iElem) - SumW/MeshParam%Area(iElem)         ! w(k+1/2,iElem) [1]
            EndDo
        EndDo
    Else
        HydroParam%wm = HydroParam%w
    EndIf
    

    ! 14. Non-hydrostatic correction:
    If (HydroParam%iNonHydro==0.and.MeshParam%Kmax>1) Then 
        Call NonHydroPressure(HydroParam,MeshParam,dt)
    EndIf

    ! 15. Compute Nodal Velocities:
    !Call Velocities(HydroParam,MeshParam)
    Call VelocitiesSUB(HydroParam,MeshParam,dt)
    !HydroParam%CFL = HydroParam%CFL*dt
    
   ! Call ExchangeTime(HydroParam,MeshParam,MeteoParam,dt)
    
    Return
End Subroutine Hydro