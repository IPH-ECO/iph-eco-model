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
    Use ELMConservative
    !Use AdvConservativeScheme
    
    Implicit None
    Integer:: iElem, iEdge, lEdge, iNode, iNewton, iLayer, INFO, iLayer_bar, Small, FlagLayer, maxF,minF,maxElem,minElem
    Integer:: llElem, rrElem, psi_flag,  indEdge
    Real:: Courant,fi_small, DZjAcum
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
    MeshParam%Ksat = 0.00005!k0 0.01 b01 0.00005
    e0 = 0.1 !e0 0.3 b01
    !k0 = 0.01 !k0 0.01 b01
    
    HydroParam%iConv = 0
    !HydroParam%iConv = 5
    !HydroParam%iNonHydro=0
    
    ! 0. Compute turbulence
    Call Turbulence(HydroParam,MeshParam,MeteoParam,dt)
        
    ! 1. Convective Term
    If (HydroParam%iConv == 0.or.HydroParam%iConv == 4) Then
        HydroParam%Fw = HydroParam%w
        HydroParam%Fu = HydroParam%u
        Call FuFv(HydroParam,MeshParam,dt)
        !Do iEdge = 1,MeshParam%nEdge
            
            !Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge) 
            !    HydroParam%Fu(iLayer,iEdge) = (1.-HydroParam%epson(iLayer,iEdge))*HydroParam%Fu(iLayer,iEdge) 
            !EndDo
            
            !if (HydroParam%IndexWaterLevelEdge(iEdge)>0) then
            !    l = MeshParam%Left(iEdge)
            !    Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
            !        cont=0
            !        vel=0
            !        Do lEdge=1,4 
            !            If (HydroParam%Fu(iLayer,MeshParam%Edge(lEdge,l))>0.and.HydroParam%IndexWaterLevelEdge(MeshParam%Edge(lEdge,l))==0) Then
            !                cont=cont+1
            !                vel = vel + HydroParam%Fu(iLayer,MeshParam%Edge(lEdge,l))
            !            EndIf
            !        EndDo
            !        if (cont>0) then    
            !            HydroParam%Fu(iLayer,iEdge) = vel/cont
            !        endif
            !    EndDo    
            !Endif 
        !EndDo
        
    ElseIf (HydroParam%iConv == 1) Then
        
    ElseIf (HydroParam%iConv == 2) Then
        !Call Convective
    ElseIf (HydroParam%iConv == 3) Then ! Neglect Nonlinear Convection
        !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iEdge)
        Do iEdge = 1,MeshParam%nEdge
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge) 
                HydroParam%Fu(iLayer,iEdge) = (1.-HydroParam%epson(iLayer,iEdge))*HydroParam%u(iLayer,iEdge)
                HydroParam%Fu(iLayer,iEdge) = (1.-HydroParam%epson(iLayer,iEdge))*HydroParam%um(iLayer,iEdge)
            EndDo
        EndDo
        HydroParam%Fw = HydroParam%w
        !!$OMP end parallel do     
    ElseIf (HydroParam%iConv == 5) Then
        
        HydroParam%Fw = HydroParam%w
        HydroParam%Fu = HydroParam%u
        Call FuFvConservative(HydroParam,MeshParam,dt)
        !Do iEdge = 1,MeshParam%nEdge
        !    
        !    !Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge) 
        !    !    HydroParam%Fu(iLayer,iEdge) = (1.-HydroParam%epson(iLayer,iEdge))*HydroParam%Fu(iLayer,iEdge) 
        !    !EndDo
        !    
        !    if (HydroParam%IndexWaterLevelEdge(iEdge)>0) then
        !        l = MeshParam%Left(iEdge)
        !        Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
        !            cont=0
        !            vel=0
        !            Do lEdge=1,4 
        !                If (HydroParam%Fu(iLayer,MeshParam%Edge(lEdge,l))>0.and.HydroParam%IndexWaterLevelEdge(MeshParam%Edge(lEdge,l))==0) Then
        !                    cont=cont+1
        !                    vel = vel + HydroParam%Fu(iLayer,MeshParam%Edge(lEdge,l))
        !                EndIf
        !            EndDo
        !            If (cont>0) then    
        !                HydroParam%Fu(iLayer,iEdge) = vel/cont
        !            Endif
        !        EndDo    
        !    endif 
        !EndDo        
        
    EndIf
    
    ! 2. Getting hydrodynamic boundary condition values at current step time
    Call GetHydroBoundaryConditions(HydroParam,MeshParam,dt,time,SimTime)
    
    if(simtime == dt*72) then
        continue
    endif

    ! 3. Baroclinic pressure effect
    Call Pressure(HydroParam,MeshParam,dt)
    
    ! 4. Viscosity and Coriolis effect
    Call ExplicitTerms(HydroParam,MeshParam,dt)
    
    ! 5. Get wind velocity components
    Call WindVelocity(HydroParam,MeshParam,MeteoParam)
    
    !! 6. Vertical Water Balance (Precipitation and Evaporation)
    !Call VerticalWB(HydroParam,MeshParam,MeteoParam,dt,SimTime)
      
    ! 7. Assemble Matrix
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam,MeteoParam) private(iEdge,iLayer,l,r,Chezy,rhoairCell,aTh,bTh,cTh,dzp,dzm,DIM,NearZero,dt)
    !If (MeshParam%iBedrock == 1) Then
    Do iEdge = 1,MeshParam%nEdge
        
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
  
        !!bench 02:
        !if ( 100-NearZero <= MeshParam%EdgeBary(1,iEdge)<= 100 + NearZero) then
        !    if (MeshParam%EdgeBary(2,iEdge) >= 110 - NearZero) then
        !        HydroParam%u(:,iEdge)  = 0.
        !        HydroParam%Fu(:,iEdge) = 0.
        !        r = 0
        !    elseif(MeshParam%EdgeBary(2,iEdge) <= 90 + NearZero) then
        !        HydroParam%u(:,iEdge)  = 0.
        !        HydroParam%Fu(:,iEdge) = 0.   
        !        r = 0
        !    endif
        !endif        
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
        If (HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge)<= HydroParam%PCRI+NearZero) Then
            HydroParam%u(:,iEdge)  = 0.
            HydroParam%Fu(:,iEdge) = 0.
        EndIf
        
        ! 7.3 Get Shear stresses in the edges
        Call Tension(HydroParam%iWindStress,HydroParam%BottomTensionFlag,iEdge,HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),HydroParam%g,HydroParam%uxy(HydroParam%Smallm(iEdge),1,iEdge),HydroParam%uxy(HydroParam%Smallm(iEdge),2,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge),Chezy,rhoairCell,HydroParam%rho0,HydroParam%windDragConstant,HydroParam%WindVel(:,iEdge),HydroParam%WindXY(:,iEdge),HydroParam%GammaT,HydroParam%GammaB)

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
            If ( r == 0 .or. HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge) <= HydroParam%PCRI+NearZero) Then
                HydroParam%GammaB = 0. 
            EndIf
            
            !DZhj is surface water thickness, in following formulation DZhj is equivalent H in [2]:
            HydroParam%iADZ(HydroParam%Smallm(iEdge),iEdge) = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)/( HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge) + dt*(HydroParam%GammaB - HydroParam%GammaT) + NearZero)
            HydroParam%iAG(HydroParam%Smallm(iEdge),iEdge)  = HydroParam%Gu(HydroParam%Smallm(iEdge),iEdge)/(HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge) + dt*(HydroParam%GammaB - HydroParam%GammaT) + NearZero )
            HydroParam%DZiADZ(iEdge)             = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)*HydroParam%iADZ(HydroParam%Smallm(iEdge),iEdge) 
            !In 2D Model iADZ = DZiA, which implies DZiAG = iADZG:
            HydroParam%DZiAG(iEdge)             = HydroParam%DZhj(HydroParam%Smallm(iEdge),iEdge)*HydroParam%iAG(HydroParam%Smallm(iEdge),iEdge)

        Else
            ! 7.5.1 Assemble Matrix A
            If ( r == 0 .or. HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge)<= HydroParam%PCRI+NearZero) Then
                HydroParam%GammaB = 0. 
            EndIf
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                If ( iLayer == HydroParam%Smallm(iEdge) ) Then
                    dzp         = 0.5*( HydroParam%DZhj(iLayer,iEdge) + HydroParam%DZhj(iLayer+1,iEdge) )       ! Dz at Upper Interface
                    aTh(iLayer) = 0. !-dt*VerEddyVisc(iLayer+1,iEdge)/dzp !-dt*nuz/dzp   !                                        ! Lower Diagonal 
                    bTh(iLayer) = MeshParam%EdgeLength(iEdge)*(HydroParam%DZhj(iLayer,iEdge) + dt*HydroParam%VerEddyVisc(iLayer+1,iEdge)*( 1/dzp ) + HydroParam%GammaB*dt) + NearZero!DZj(iLayer,iEdge) + dt*nuz*( 1/dzp ) + GammaB*dt      ! Diagonal
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
    
 !   !7.6 Assemble Matrix DZK
	!HydroParam%DZK(:)   = 0.d0
 !   If (MeshParam%iBedrock == 1) Then
 !       Do iEdge = 1,MeshParam%nEdge       
 !           If (HydroParam%Smallms(iEdge) == HydroParam%CapitalMs(iEdge).and.HydroParam%Smallms(iEdge) == HydroParam%CapitalM(iEdge)) Then
 !               HydroParam%DZK(iEdge)   = HydroParam%DZsj(HydroParam%Smallms(iEdge),iEdge)*MeshParam%Kj(HydroParam%Smallms(iEdge),iEdge) !Sediment Layer
 !           Else
 !               HydroParam%DZK(iEdge)   = Dot_Product(HydroParam%DZsj(HydroParam%Smallms(iEdge):HydroParam%CapitalMs(iEdge),iEdge),MeshParam%Kj(HydroParam%Smallms(iEdge):HydroParam%CapitalMs(iEdge),iEdge) )
 !           EndIf
 !       EndDo
 !   EndIf   
    
    !Call Volume(HydroParam,MeshParam)
    !Do iElem = 1,MeshParam%nElem
    !    Call MoistureContent(HydroParam%eta(iElem),HydroParam%etaplus(iElem),iElem,HydroParam,MeshParam)
    !EndDo
    ! 8. Compute the New Free-Surface Elevation
    ! 8.1 Assemble the Right Hand Side (RHS)
    !$OMP parallel do default(none) shared(MeshParam,HydroParam,NearZero,dt) private(iElem,iEdge,SumRHS,SumLhS,gamma,Face,Pij)
    Do iElem = 1, MeshParam%nElem
        SumRHS = 0d0
        SumLhS = 0.d0
        gamma = 0.d0
        Do iEdge = 1,4
            Face = MeshParam%Edge(iEdge,iElem)
            SumRHS = SumRHS + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*((1.d0-HydroParam%Theta)*(Dot_product( HydroParam%DZhj(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face),HydroParam%u(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face)) + Dot_product(HydroParam%DZsj(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face),HydroParam%us(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face))) + HydroParam%Theta*HydroParam%DZiAG(Face))
            !SumRHS = SumRHS + Sig(iElem,MeshParam%Right(Face),MeshParam%Left(Face))*MeshParam%EdgeLength(Face)*((1.d0-HydroParam%Theta)*(Dot_product( HydroParam%DZhj(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face),HydroParam%u(HydroParam%Smallm(Face):HydroParam%CapitalM(Face),Face)) + Dot_product(HydroParam%DZsj(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face),HydroParam%us(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face))) + HydroParam%Theta*(HydroParam%DZiAG(Face) + Dot_product(HydroParam%DZsj(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face),HydroParam%Gusub(HydroParam%Smallms(Face):HydroParam%CapitalMs(Face),Face))))
            ! 8.1.1 If there is a Pressure Boundary Condition
            Pij = MeshParam%Neighbor(iEdge,iElem) 
            If (Pij == 0.and.HydroParam%IndexWaterLevelEdge(Face)>0) Then   
                HydroParam%etaInfn(iElem) =  HydroParam%etaInf(iElem)
				If ((HydroParam%WaterLevel(HydroParam%IndexWaterLevelEdge(Face))-HydroParam%hj(Face))<HydroParam%PCRI/2.d0+NearZero) Then
                    HydroParam%etaInf(iElem) = HydroParam%hj(Face) + HydroParam%PCRI/2.d0 !Rever
                Else
				    HydroParam%etaInf(iElem) = HydroParam%WaterLevel(HydroParam%IndexWaterLevelEdge(Face))
                EndIf
                SumLHS = SumLHS + ( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( HydroParam%etaInf(iElem) )*(HydroParam%g*HydroParam%Theta*dt*HydroParam%DZiADZ(Face) + HydroParam%DZK(Face))
                !SumLHS = SumLHS + ( MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face) )*( HydroParam%etaInf(iElem) )*(HydroParam%g*HydroParam%Theta*dt*HydroParam%DZiADZ(Face) + HydroParam%Theta*HydroParam%DZK(Face))
            EndIf
        EndDo
        HydroParam%rhs(iElem) = HydroParam%Vol(iElem) - dt*SumRHS + (HydroParam%Theta*dt)*SumLHS
    EndDo
    !$OMP end parallel do

   
    ! 8.2 Newton Loop for Non-Linear Wet- and Dry-ing Algorithm [2]
    HydroParam%etan = HydroParam%eta   
    HydroParam%eta = HydroParam%eta + HydroParam%etaplus
    !HydroParam%etak = HydroParam%eta + HydroParam%etaplus   
    Do iNewton = 1,200
        ! 8.2.1 Assemble Matrices F (System Matrix - Newton Method) and P (Derivative Matrix)
        HydroParam%P = 0.
        HydroParam%F = 0.
        Call MatOp(HydroParam%eta,HydroParam%Aeta,dt,HydroParam,MeshParam)
        !Call Volume(HydroParam,MeshParam)
        !$OMP parallel do default(none) shared(HydroParam,MeshParam) private(iElem)      
        Do iElem = 1, MeshParam%nElem
            !8.2.2 Volume with Eta at tn+1:
            Call MoistureContent(HydroParam%eta(iElem),0.d0,iElem,HydroParam,MeshParam)            

            !8.2.3 Compute Newton Method's Residue for each iElem:
            !In this point, MatOp Output = (T-Matrix)Eta
            !We want F = 0 condition satisfied: V(nt+1) + T*n(t+1) = rhs :: F  = V(n(t+1)) + T*n(t+1) - rhs -> 0
            HydroParam%F(iElem) = HydroParam%Vol(iElem) + HydroParam%Aeta(iElem) - HydroParam%rhs(iElem)
            
            !8.2.4 Fill P values in dry cell for CGOp computation:
            !HydroP = A*dn/dn, if cell is dry dn/dn = 0, else dn/dn = 1
            If(HydroParam%eta(iElem) > HydroParam%hb(iElem)) Then
                HydroParam%P(iElem) = MeshParam%Area(iElem) ! Vol = eta.Area -> Vol = (eta - hb)*Area + (hb-sb)*Area*e = Area*(eta-hb +hb*e -sb*e)
                !HydroParam%P(iElem) = (HydroParam%Vol(iElem)/(HydroParam%eta(iElem) + HydroParam%hb(iElem)*(MeshParam%epson(HydroParam%ElCapitalMs(iElem)) - 1) - HydroParam%sb(iElem)*MeshParam%epson(HydroParam%ElCapitalMs(iElem))
            Else
                HydroParam%P(iElem) = (HydroParam%Vol(iElem)/HydroParam%eta(iElem)) ! Vol = eta*e*S*Area -> Vol/eta = e*S*Area
            EndIf
        EndDo   
        !$OMP end parallel do
        res = sqrt(sum(HydroParam%F**2))      ! Residual of the Method
        !Print*, 'iNewton = ',iNewton , 'res = ',res
        If ( res < 1e-6 ) Then
            continue
            exit
        EndIf
        
        ! V(n(t+1)) + T*n(t+1) = b
        
        ! (P+T).eta = rhs
        
        !! 8.2.5 Compute the New Free-Surface Elevation
        !!CGOp is used to minimize F value :: F = V(n(t+1)) + T*n(t+1) - rhs == 0
        Call CGOp(HydroParam%F,HydroParam%Deta,dt,HydroParam,MeshParam)

        !$OMP parallel do default(none) shared(HydroParam,MeshParam) private(iElem)
        Do iElem = 1, MeshParam%nElem
            HydroParam%eta(iElem) = HydroParam%eta(iElem) - HydroParam%Deta(iElem)
        EndDo
        !$OMP end parallel do
    EndDo
    
    ! 9. Water Level Corrective Step
    !$OMP parallel do default(none) shared(MeshParam,HydroParam,NearZero) private(iElem)
    Do iElem = 1,MeshParam%nElem
        If ( HydroParam%eta(iElem) - HydroParam%sb(iElem) <= HydroParam%PCRI/2.0d0 + NearZero ) Then
            HydroParam%eta(iElem) = HydroParam%sb(iElem) + HydroParam%PCRI/2.0d0
        EndIf
    EndDo
    !$OMP end parallel do
    
    ! 9.1 Compute nodal elevations for tang. vel.
    HydroParam%petan=HydroParam%peta !store for transport eqs. 
    !$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iNode,sum1,sum0,ie,j)
    Do iNode=1,MeshParam%nNode
        sum1=0 !sum of elemental elevations
        sum0=0 !sum of Areas
        do j=1,MeshParam%nVertexElem(iNode)
            ie=MeshParam%VertexElem(j,iNode)
            sum1=sum1+MeshParam%Area(ie)*HydroParam%eta(ie)
            sum0=sum0+MeshParam%Area(ie)
        Enddo !j=1,nne(i)
        HydroParam%peta(iNode)=sum1/sum0
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
            HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) = HydroParam%hb(iElem) + HydroParam%PCRI !- hb(iElem)       ! Free-Surface (verificar com Rafael)
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
            HydroParam%DZi(iLayer,iElem) = Max(HydroParam%Pcri,HydroParam%DZi(iLayer,iElem))
        EndDo
        
        !xx. Porosity and and DZhi/DZsi
        HydroParam%DZhi(:,iElem) = HydroParam%DZi(:,iElem)
		If (MeshParam%iBedrock == 1) Then
            Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem)
            
                If (iLayer >= HydroParam%ElSmallm(iElem).or.HydroParam%ElSmallms(iElem) ==  HydroParam%ElSmallm(iElem)) Then
                    continue
                Else
                    MeshParam%ei(iLayer,iElem) = e0
                EndIf

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
        
        HydroParam%H(iEdge) = Max( HydroParam%PCRI,-HydroParam%sj(iEdge) + HydroParam%eta(l), -HydroParam%sj(iEdge) + HydroParam%eta(r) )!CAYO
        
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
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%hj(iEdge) + HydroParam%PCRI
            Else
                HydroParam%Z(HydroParam%CapitalM(iEdge)+1,iEdge) = HydroParam%eta(l)
            EndIf
        Else
            If ( HydroParam%H(iEdge) - HydroParam%hj(iEdge) <= HydroParam%PCRI+NearZero ) Then !( HydroParam%H(iEdge) <= HydroParam%PCRI/2 + NearZero ) Then
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
            HydroParam%DZj(iLayer,iEdge) = Max(HydroParam%Pcri,HydroParam%DZj(iLayer,iEdge))
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
                    If ( HydroParam%Z(iLayer+1,iEdge) > HydroParam%hj(iEdge) ) Then
                        HydroParam%DZhj(iLayer,iEdge) = HydroParam%Z(iLayer+1,iEdge) - HydroParam%hj(iEdge)
                        HydroParam%DZsj(iLayer,iEdge) = HydroParam%hj(iEdge) - HydroParam%Z(iLayer,iEdge)
                    Else
                        HydroParam%DZhj(iLayer,iEdge) = 0.
                        HydroParam%DZsj(iLayer,iEdge) = Max(HydroParam%Pcri,HydroParam%Z(iLayer+1,iEdge) - HydroParam%Z(iLayer,iEdge))
                        If(HydroParam%Smallms(iEdge) == HydroParam%CapitalMs(iEdge) .and. MeshParam%iSaturation == 0) Then
                            If(HydroParam%H(iEdge) < HydroParam%hj(iEdge) ) Then
                                HydroParam%DZsj(iLayer,iEdge) = HydroParam%H(iEdge)
                            ElseIf(HydroParam%H(iEdge) <= HydroParam%sj(iEdge))Then
                                HydroParam%DZsj(iLayer,iEdge) = 0.d0
                             EndIf   
                        EndIf                      
                    EndIf
                EndIf    
                
                !xx. Hydraulic Conductivity:
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
    Call VelocitiesSUB(HydroParam,MeshParam)
    
   ! Call ExchangeTime(HydroParam,MeshParam,MeteoParam,dt)
    
    Return
End Subroutine Hydro