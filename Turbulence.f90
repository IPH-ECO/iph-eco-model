Subroutine Turbulence(HydroParam,MeshParam,MeteoParam,dt)
    
    ! Calculate the eddy viscosity using different Turbulence Models
    ! Based on:
    ! [1] Ji, Z-G. Hydrodynamics and Water Quality: modeling rivers, lakes, and estuaries. John Wiley & Sons, 2008.
    ! [2] Smagorinsky, J. 1963. General circulation experiments with the primitive equations, part I: the basic experiment. 
    !       Monthly Weather Review, (91), 99-152.
    ! [3] Pacanowski, R.C.; Philander, S.G.H. 1981. Parameterization of vertical mixing in numerical models of tropical oceans.
    !       Journal of Physical Oceanography, (11), 1443-1451.
    ! [4] Umlauf, L.; Burchard, H. 2003. A generic length-scale equation for geophysical turbulence models.
    !       Journal of Marine Research, (61), 235-265.
    ! [5] Zhang Y.; Baptista, A.M.; Myers, E.P. 2004. A cross-scale model for 3D baroclinic circulation in estuary-plume-shelf systems: I. Formulation and skill assessment.
    !       Continental Shelf Research, (24), 2187-2214.
    ! [6] Mellor, G.L.; Yamada, T. 1982. Development of a turbulence closure model for geophysical fluid problems.
    !       Reviews of Geophysics and Space Physics, (20), 851-875.
    ! [7] Galperin, B.; Kantha, L.H.; Hassid, S.; Rosati, A. 1988. A quasi-equilibrium turbulent energy model for geophysical flows.
    !       Journal of Atmospheric Sciences, (45), 55-62.
    ! [8] Blumberg, A.F.; Galperin, B.; O'Connor, D.J. 1992. Modeling vertical structure of open-channel flow. 
    !       Jouornal of Hydraulic Engineering, (118), 1119-1134.
    ! [9] Hodges, B.R.; Imberger, J.; Saggio, A.; Winters, K.B. 2000. Modeling basin-scale internal waves in a stratified lake.
    !       Limnology and Oceanography, 45 (7), 1603-1620.
    ! Input:
    ! Hydrodynamic Features
    ! Output:
    ! Eddy Viscosity Coefficients
    
    ! List of Modifications:
    !   02.01.2015: Routine Implementation       (Rafael Cavalcanti)
    ! Programmer: Rafael Cavalcanti
    !$ use omp_lib
    Use MeshVars !, Only: nElem,Edge,Area,NEdge,Left,Right, dx, dy
    Use Hydrodynamic
    Use Meteorological
    
    Implicit None
    Integer:: iElem, iEdge, iLayer
    Integer:: l, r, Pij
    !Real:: lC(0:5), ln(2,4), llambda(4), lw(2)
    !Real:: lEdgeBary(2,4), lNeigBary(2,4), lxNode(4), lyNode(4)
    Real:: rich,vdiff,tdiff,drhodz,bvf,dundz,dutdz,shear2
    Real:: dudx,dvdy,dvdx,dudy
    Real:: dt
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam
    
    
    ! Note: HydroParam%DZit == HydroParam%DZhit, in subsurface coupled case DZhit is equivalent to Surface Water thickness in cell
    ! in only surface case the DZhit = DZit.
    ! Same case for HydroParm%DZj == HydroParam%DZhj, with DZhj for Edge.
    
    ! 1. Define the Cell centered Horizontal Eddy-Viscosity:
    If (HydroParam%iHTurb == 0) Then         ! User-defined Horizontal Diffusion
        !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iLayer,iElem)
        Do iElem = 1, MeshParam%nElem
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                HydroParam%HorViscosity(:,iLayer,iElem) = (/ HydroParam%HorEddyViscX_Cte, HydroParam%HorEddyViscY_Cte /)
                HydroParam%HorDiffusivity(:,iLayer,iElem) = (/ HydroParam%HorEddyDiffX_Cte, HydroParam%HorEddyDiffY_Cte /)
            EndDo       ! Layer Loop
        EndDo       ! Element Loop
        
        
        !!$OMP end parallel do
    Elseif (HydroParam%iHTurb == 1) Then        ! Smagorinsky Model
        ! Horizontal Eddy-Viscosity Calculation
        !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iLayer,iElem,dudx,dvdy,dvdx,dudy)
        Do iElem = 1, MeshParam%nElem
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                dudx = (HydroParam%uxy(iLayer,1,MeshParam%Edge(4,iElem))-HydroParam%uxy(iLayer,1,MeshParam%Edge(2,iElem)))/MeshParam%DX
                dvdy = (HydroParam%uxy(iLayer,2,MeshParam%Edge(1,iElem))-HydroParam%uxy(iLayer,2,MeshParam%Edge(3,iElem)))/MeshParam%DY
                dvdx = (HydroParam%uxy(iLayer,2,MeshParam%Edge(4,iElem))-HydroParam%uxy(iLayer,2,MeshParam%Edge(2,iElem)))/MeshParam%DX
                dudy = (HydroParam%uxy(iLayer,1,MeshParam%Edge(1,iElem))-HydroParam%uxy(iLayer,1,MeshParam%Edge(3,iElem)))/MeshParam%DY
                HydroParam%HorViscosity(1,iLayer,iElem) = HydroParam%CSmag*MeshParam%Area(iElem)*sqrt(dudx**2+dvdy**2+(dvdx+dudy)**2/2)
                HydroParam%HorViscosity(2,iLayer,iElem) = HydroParam%HorViscosity(1,iLayer,iElem)
                HydroParam%HorDiffusivity(:,iLayer,iElem) = (/ HydroParam%HorViscosity(1,iLayer,iElem)+HydroParam%HorEddyDiffY_Back, HydroParam%HorViscosity(1,iLayer,iElem)+HydroParam%HorEddyDiffY_Back /) ! According to Ji(2008), hvisc = hdiff 
            EndDo       ! Layer Loop
        EndDo       ! Element Loop
        !!$OMP end parallel do
    Else
        Print*, 'Wrong Option of Horizontal Turbulence Model.'
        Print*, 'Program Ended.'
        Stop
    EndIf
    
    ! 2. Define the Vertical Eddy-Viscosity and Eddy-Diffusivity
    ! Obs.: Eddy-Viscosity: In the Momentum Equations
    !       Eddy-Diffusivity: In the Transport Equations
    If ( HydroParam%iVTurb == 0 ) Then         ! User-defined Vertical Viscosity and Diffusivity
        ! Obs. Each Layer define the Turbulent Diffusivity at the Bottom of the cell
        Do iEdge = 1, MeshParam%nEdge
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                HydroParam%VerEddyVisc(iLayer,iEdge) = HydroParam%VerEddyVisc_Cte
                HydroParam%VerEddyDiff(iLayer,iEdge) = HydroParam%VerEddyDiff_Cte
            EndDo       ! Layer Loop
        EndDo       ! Element Loop
        
        Do iElem = 1, MeshParam%nElem
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                HydroParam%VerEddyViscCell(iLayer,iElem) = HydroParam%VerEddyVisc_Cte
                HydroParam%VerEddyDiffCell(iLayer,iElem) = HydroParam%VerEddyDiff_Cte
            EndDo       ! Layer Loop
        EndDo       ! Element Loop
        
    Elseif ( HydroParam%iVTurb == 1 ) Then        ! Zero-Equation Model [3]
        !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iLayer,iEdge,l,r,drhodz,bvf,dundz,dutdz,shear2,rich)
        Do iEdge = 1, MeshParam%nEdge
            l = MeshParam%Left(iEdge)
            r = MeshParam%Right(iEdge)
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                If ( HydroParam%Smallm(iEdge) == HydroParam%CapitalM(iEdge) ) Then        ! Only One Vertical Layer
                    HydroParam%VerEddyVisc(iLayer,iEdge) = 0.
                    HydroParam%VerEddyDiff(iLayer,iEdge) = 0.
                Else
                    If ( iLayer == HydroParam%Smallm(iEdge) ) Then !Bottom layer (m-1/2)
                        HydroParam%VerEddyVisc(iLayer,iEdge) = 0. 
                        HydroParam%VerEddyDiff(iLayer,iEdge) = 0. 
                    Else !Intermediary (k-1/2) or Surface layer (M-1/2)
                        If (r == 0) Then
                            drhodz = (HydroParam%sDRhoW(iLayer,l)-HydroParam%sDRhoW(iLayer-1,l))/((HydroParam%DZhj(iLayer,iEdge)+HydroParam%DZhj(iLayer-1,iEdge))/2.)
                        Else
                            drhodz = (0.5*(HydroParam%sDRhoW(iLayer,l)+HydroParam%sDRhoW(iLayer,r))-0.5*(HydroParam%sDRhoW(iLayer-1,l)+HydroParam%sDRhoW(iLayer-1,r)))/((HydroParam%DZhj(iLayer,iEdge)+HydroParam%DZhj(iLayer-1,iEdge))/2.)
                        EndIf
		                bvf    = -HydroParam%G*(drhodz/HydroParam%rho0+HydroParam%G/1.5e3**2.)
                        dundz  = (HydroParam%uxy(iLayer,1,iEdge)-HydroParam%uxy(iLayer-1,1,iEdge))/((HydroParam%DZhj(iLayer,iEdge)+HydroParam%DZhj(iLayer-1,iEdge))/2.)
		                dutdz  = (HydroParam%uxy(iLayer,2,iEdge)-HydroParam%uxy(iLayer-1,2,iEdge))/((HydroParam%DZhj(iLayer,iEdge)+HydroParam%DZhj(iLayer-1,iEdge))/2.)
                        shear2 = dundz**2. + dutdz**2.
		                shear2 = max(shear2,1.0d-10)
		                rich   = max(bvf/shear2,0.0d0)
                        HydroParam%VerEddyVisc(iLayer,iEdge) = HydroParam%vref/((1.+HydroParam%alfa_turbmodel1*rich)**HydroParam%n_turbmodel1)+HydroParam%vmin
                        HydroParam%VerEddyDiff(iLayer,iEdge) = HydroParam%VerEddyVisc(iLayer,iEdge)/((1.+HydroParam%alfa_turbmodel1*rich)**HydroParam%n_turbmodel1)+HydroParam%tdmin_pp
                    EndIf
                EndIf
            EndDo       ! Layer Loop
        EndDo       ! Element Loop
        !!$OMP end parallel do
         
        
        Do iElem = 1, MeshParam%nElem
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                If ( HydroParam%ElSmallm(iElem) == HydroParam%ElCapitalM(iElem) ) Then        ! Only One Vertical Layer
                    HydroParam%VerEddyViscCell(iLayer,iElem) = HydroParam%vmin
                    HydroParam%VerEddyDiffCell(iLayer,iElem) = HydroParam%tdmin_pp
                Else
                    If ( iLayer == HydroParam%ElSmallm(iElem) ) Then !Bottom layer (m-1/2)
                        HydroParam%VerEddyViscCell(iLayer,iElem) = HydroParam%vmin 
                        HydroParam%VerEddyDiffCell(iLayer,iElem) = HydroParam%tdmin_pp 
                    Else !Intermediary (k-1/2) or Surface layer (M-1/2)
                        drhodz = (HydroParam%sDRhoW(iLayer,iElem)-HydroParam%sDRhoW(iLayer-1,iElem))/((HydroParam%DZi(iLayer,iElem)+HydroParam%DZi(iLayer-1,iElem))/2.)
		                bvf    = -HydroParam%G*(drhodz/HydroParam%rho0+HydroParam%G/1.5e3**2.)
                        dundz  = (HydroParam%ub(iLayer,1,iElem)-HydroParam%ub(iLayer-1,1,iElem))/((HydroParam%DZi(iLayer,iElem)+HydroParam%DZi(iLayer-1,iElem))/2.)
		                dutdz  = (HydroParam%ub(iLayer,2,iElem)-HydroParam%ub(iLayer-1,2,iElem))/((HydroParam%DZi(iLayer,iElem)+HydroParam%DZi(iLayer-1,iElem))/2.)
                        shear2 = dundz**2. + dutdz**2.
		                shear2 = max(shear2,1.0d-10)
		                rich   = max(bvf/shear2,0.0d0)
                        HydroParam%VerEddyViscCell(iLayer,iElem) = HydroParam%vref/((1.+HydroParam%alfa_turbmodel1*rich)**HydroParam%n_turbmodel1)+HydroParam%vmin
                        HydroParam%VerEddyDiffCell(iLayer,iElem) = HydroParam%VerEddyViscCell(iLayer,iElem)/((1.+HydroParam%alfa_turbmodel1*rich)**HydroParam%n_turbmodel1)+HydroParam%tdmin_pp
                    EndIf
                EndIf
            EndDo       ! Layer Loop
        EndDo       ! Element Loop
        !!$OMP end parallel do        
    Elseif ( HydroParam%iVTurb == 2 ) Then        ! Generic Length Scale Model [4,5]
        Continue
    Elseif ( HydroParam%iVTurb == 3 ) Then        ! Mellor and Yamada Turbulence Model [6,7,8]
        Call TwoEq_KC(HydroParam,MeshParam,MeteoParam,dt)
    Elseif ( HydroParam%iVTurb == 4 ) Then        ! ELCOM Mixing Model [9]
        Continue
    Else
        Print*, 'Wrong Option of Vertical Turbulence Model.'
        Print*, 'Program Ended.'
        Stop
    EndIf
    !write(*,*) HorViscosity(:,1,1),VerEddyVisc(2,3)
    Return
End Subroutine Turbulence
 
Subroutine TwoEq_KC(HydroParam,MeshParam,MeteoParam,dt)

    ! Compute the Turbulent Vertical Eddy Viscosity and Diffusivity Using a Two Equation Closure Scheme
    
    ! Based on: 
    ! [1] Kantha,L.H.; Clayson,C.A. An improved mixed layer model for geophysical applications.
    !   Journal of Geophysical Research, v. 99 (C12), pp. 25235-25266, 1994.
    ! [2] Rueda,F.J.; Schladow,S.G. Dynamics of Large Polymictic Lake. II: Numerical Simulations.
    !   Journal of Hydraulic Engineering, v. 129(2), pp. 92-101, 2003.
    ! [3] Galperin,B.; Kantha,L.H.; Hassid,S.; Rosati,A. A Quasi-equilibrium Turbulent Energy Model for Geophysical Flows.
    !   Journal of the Atmospheric Sciences, v. 45(1), pp. 55-62, 1988.

    ! Input:
    !   Water Temperature
    !   Water Density
    !   Horizontal Velocity
    
    ! Output:
    !   Vertical Eddy Viscosity
    !   Vertical Eddy Diffusivity
    !   Richardson Number
    !   Friction Velocity
  
    ! List of Modifications:
    !   23.08.2017: Routine Implementation       (J. Rafael Cavalcanti)
    
    ! Programmer: J. Rafael Cavalcanti
    Use Hydrodynamic !, Only: nElem, ElCapitalMOld, ElSmallmOld, Dzin, ub, KMax
    Use MeshVars
    Use Meteorological
    !Use Hydrodynamic, Only: Ri, LengthScale, VerEddyVisc, VerEddyDiff
    !Use Hydrodynamic, Only: TKE, LengthScale, TKEP, LengthScaleP, DissipRate
    !Use Transport, Only: sDRhoW
    !Use Param, Only: g, BkGndVisc, BkGndDiff
    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam
    Integer:: iElem, iLayer
    Real:: dzp, RhoAv, drhodz, ShearSq, N2, rich
    Real, Dimension(MeshParam%KMax+1,MeshParam%nElem):: kpp, kp, k
    Real, Dimension(MeshParam%KMax+1,MeshParam%nElem):: klpp, klp, kl
    Real:: Gh, Gh_Min_KC, Gh_Max_KC, Sh, Sm
    Real:: t_1, t_2, t_3, t_4, t_5, Rt
    Real, Parameter:: A_1 = 0.92, A_2 = 0.74, B_1 = 16.6, B_2 = 10.1, C_3 = 0.20, C_2 = 0.70
    Real:: dt
    HydroParam%VerEddyViscCell = 0.
    HydroParam%VerEddyDiffCell = 0.
    
        
    ! ----- Set Coefficients -----
    Gh_Max_KC = 2.9E-2; Gh_Min_KC = -2.8E-1
    t_1 = A_2 * (1.-6.*A_1/B_1)
    t_2 = 3.*A_2*(6.*A_1+B_2*(1.- C_3))
    t_3 = B_1**(-1./3.)
    t_4 = 9.*A_1*(2.*A_1+A_2*(1.-C_2)) 
    t_5 = 9.*A_1*A_2    
    
    ! 1. Solve the System for TKE and LengthScale
    !DissipRate = 0.
    kp = 2.*HydroParam%TKE; kpp = 2.*HydroParam%TKEP;
    klp = HydroParam%LengthScale*HydroParam%TKE; klpp = HydroParam%LengthScaleP*HydroParam%TKEP
    Call SolveTKE_LengthScale(HydroParam,MeshParam,MeteoParam,dt,kpp,klpp,kp,klp,k,kl)
    
    ! 2. Update Vertical Eddy Viscosity and Diffusivity
    Do iElem = 1, MeshParam%nElem
        If ( HydroParam%ElSmallm(iElem) == HydroParam%ElCapitalM(iElem) ) Then      ! Bidimensional Model or Bidimensional Cell
            Cycle
        Else
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) - 1
                dzp     = 0.5*( HydroParam%DZhit(iLayer+1,iElem) + HydroParam%DZhit(iLayer,iElem) )
                RhoAv   = 0.5*( HydroParam%sDRhoW(iLayer+1,iElem) + HydroParam%sDRhoW(iLayer,iElem) ) !+ 1000.
                drhodz  = ( HydroParam%sDRhoW(iLayer+1,iElem) - HydroParam%sDRhoW(iLayer,iElem) )/dzp
                ShearSq = ( (HydroParam%ub(iLayer+1,1,iElem) - HydroParam%ub(iLayer,1,iElem))/dzp )**2. + ( (HydroParam%ub(iLayer+1,2,iElem) - HydroParam%ub(iLayer,2,iElem))/dzp )**2.
                N2      = -HydroParam%g*drhodz/rhoAv
                rich = N2/Max( ShearSq, 1.E-14 )
                ! 2.1 Update TKE, Lenght Scale, and Dissipation
                HydroParam%TKE(iLayer+1,iElem)         = k(iLayer+1,iElem)/2.
                HydroParam%LengthScale(iLayer+1,iElem) = kl(iLayer+1,iElem)/k(iLayer+1,iElem)
                ! ----- Length Scale Limitation from [3] -----
                If (N2 < 1.E-6) N2 = 1.E-12
                If (N2 > 0.0) HydroParam%LengthScale(iLayer+1,iElem) = Min(HydroParam%LengthScale(iLayer+1,iElem),0.53*Sqrt(k(iLayer+1,iElem)/N2))
                !DissipRate(iLayer+1,iElem) = k(iLayer+1,iElem)*DSqrt( kp(iLayer+1,iElem) )/( B_1*LengthScale(iLayer+1,iElem) )
                ! 2.2 Evaluate Transfer Coefficients at time n
                Gh = -N2*HydroParam%LengthScale(iLayer+1,iElem)**2./k(iLayer+1,iElem)
                Gh = Max( Gh_min_KC, Min( Gh_max_KC, Gh ) )
                Sh = t_1/( 1. - Gh*t_2 )
                Sm = ( t_3 + Sh*Gh*t_4 )/( 1. - t_5*Gh )
                Rt = DSqrt(k(iLayer+1,iElem))*HydroParam%LengthScale(iLayer+1,iElem)
                If ( k(iLayer+1,iElem) >= 1.E-6 ) Then      ! Mixed Layer
                    HydroParam%VerEddyViscCell(iLayer+1,iElem) = Rt*Sm + HydroParam%vmin
                    HydroParam%VerEddyDiffCell(iLayer+1,iElem) = Rt*Sh + HydroParam%tdmin_pp
                Else        ! Bellow the Mixed Layer
                    ! ----- Kantha and Clayson Model -----
                    If ( rich >= 0.7 ) Then
                        HydroParam%VerEddyViscCell(iLayer+1,iElem) = 1.E-4
                        HydroParam%VerEddyDiffCell(iLayer+1,iElem) = 5.E-5
                    ElseIf ( rich >= 0.0 .AND. rich < 0.7 ) Then
                        HydroParam%VerEddyViscCell(iLayer+1,iElem) = 1.E-4 + 5.E-3*( 1 - (rich/0.7)**2. )**3.
                        HydroParam%VerEddyDiffCell(iLayer+1,iElem) = 5.E-5 + 5.E-3*( 1 - (rich/0.7)**2. )**3.
                    Else
                        HydroParam%VerEddyViscCell(iLayer+1,iElem) = 1.E-4 + 5.E-3
                        HydroParam%VerEddyDiffCell(iLayer+1,iElem) = 5.E-5 + 5.E-3
                    End If
                End If
            End Do      ! Layer Loop            
        End If      ! Bidimensional Cell
        ! 2.3 Correction for Shallow Surface Elements
        ! Obs.: This is how SI3D does it. I don't agree with this. Need to double check!
        If ( HydroParam%DZhit(HydroParam%ElCapitalM(iElem),iElem) < 1.E-2 ) Then
            HydroParam%VerEddyViscCell(HydroParam%ElCapitalM(iElem),iElem) = Max(HydroParam%VerEddyViscCell(HydroParam%ElCapitalM(iElem),iElem),1.E-4)
            HydroParam%VerEddyDiffCell(HydroParam%ElCapitalM(iElem),iElem) = Max(HydroParam%VerEddyDiffCell(HydroParam%ElCapitalM(iElem),iElem),1.E-4)
        Else
            Continue
        End If
        ! 2.4 Assign Transfer Coefficients at the Top and Bottom Surfaces
        HydroParam%VerEddyViscCell(1:HydroParam%ElSmallm(iElem),iElem) = 0.0; HydroParam%VerEddyViscCell(HydroParam%ElCapitalM(iElem)+1:MeshParam%KMax,iElem) = 0.0
        HydroParam%VerEddyDiffCell(1:HydroParam%ElSmallm(iElem),iElem) = 0.0; HydroParam%VerEddyDiffCell(HydroParam%ElCapitalM(iElem)+1:MeshParam%KMax,iElem) = 0.0
    End Do      ! Element Loop
        
    HydroParam%TKEP = HydroParam%TKE
    HydroParam%LengthScaleP = HydroParam%LengthScale
    
    Return
End Subroutine TwoEq_KC
    
Subroutine SolveTKE_LengthScale(HydroParam,MeshParam,MeteoParam,dt,q2pp,q2lpp,q2p,q2lp,q2,q2l)
        
    ! Solve a 2 Equation Model For Turbulence Quantities (2*TKE and 2*TKE*LengthScale)
    
    ! Based on: 
    ! [1] Kantha,L.H.; Clayson,C.A. An improved mixed layer model for geophysical applications.
    !   Journal of Geophysical Research, v. 99 (C12), pp. 25235-25266, 1994.
    ! [2] Rueda,F.J.; Schladow,S.G. Dynamics of Large Polymictic Lake. II: Numerical Simulations.
    !   Journal of Hydraulic Engineering, v. 129(2), pp. 92-101, 2003.
    ! [3] Galperin,B.; Kantha,L.H.; Hassid,S.; Rosati,A. A Quasi-equilibrium Turbulent Energy Model for Geophysical Flows.
    !   Journal of the Atmospheric Sciences, v. 45(1), pp. 55-62, 1988.
    ! [4] Fletcher, C.A.J. Computational Techniques for Fluid Dynamics.
    !   Springer, New York, 1991.

    ! Input:
    
    ! Output:
  
    ! List of Modifications:
    !   05.09.2017: Routine Implementation       (J. Rafael Cavalcanti)
    
    ! Programmer: J. Rafael Cavalcanti

    Use Hydrodynamic !, Only: nElem, KMax, ElSmallmOld, ElCapitalMOld, Dzin, etan, hb, ub, Rough, Area
    !Use Hydrodynamic, Only: VerEddyViscP, VerEddyDiffP, TensionT, TensionB, uStarB, uStarT, ShearProd, BuoyancyProd
    Use MeshVars
    Use Meteorological
    !Use Environmental, Only: Wind
    !Use Transport, Only: sDRhoW
    !Use Param, Only: g, dt, Rho0, VonKarman
    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam
    Integer:: iElem, iLayer, BasinScaleFlag, DIM, iTurb
    Real, Parameter:: Sq = 0.20, q2_min = 2.e-6, q2l_min = 6.80236723501E-009
    Real, Dimension(MeshParam%KMax+1,MeshParam%nElem), Intent(In):: q2p, q2pp, q2lp, q2lpp
    Real, Dimension(MeshParam%KMax+1,MeshParam%nElem), Intent(Out):: q2, q2l
    Real, Dimension(MeshParam%KMax+1):: xl, zFromB, zFromT
    Real, Dimension(MeshParam%KMax+1):: dCoeff, dsTh, SolTh
    Real, Dimension(3,MeshParam%KMax+1):: aaTh
    !Real:: ShearProduction, BuoyancyProduction, Dissipation, TKinE
    Real:: V, zFromB0, zFromT0, VCell(2,2),RhoAir
    Real:: dzp, RhoAv, drhodz, ShearSq, N2, rhoxh
    Real, Parameter:: B_1 = 16.6, E_1 = 1.80, E_2 = 1.33, E_3 = 0.25
    Real:: GammaB,GammaT,Chezy,AirDens,SphAir
    Real:: dt
    Real:: NearZero = 1e-10
    Real, Parameter:: VonKarman = 0.4
    Real:: SpecHum,uStarB,uStarT
    
    ! ----- Initialize Basin Scale Variables -----
    uStarB = 0.; uStarT = 0.
    q2 = 0.; q2l = 0.; HydroParam%ShearProd = 0.; HydroParam%BuoyancyProd = 0.
    dCoeff = 0.
    aaTh = 0.; dsTh = 0.; SolTh = 0.
    BasinScaleFlag = 0; iTurb = 1
    !ShearProduction = 0.; BuoyancyProduction = 0.; Dissipation = 0.; TKinE = 0.
    !If ( BasinScaleFlag == 1 ) Then
    !    ShearProduction = 0.; BuoyancyProduction = 0.; Dissipation = 0.; TKinE = 0.
    !Else
    !    Continue
    !End If      ! If ( BasinScaleFlag == 1 ) Then
    
    ! 1. Compute the New TKE and LengthScale Field
    Do iElem = 1, MeshParam%nElem
        If ( HydroParam%ElSmallm(iElem) == HydroParam%ElCapitalM(iElem) ) Then      ! Bidimensional Model or Bidimensional Cell
            Cycle
        Else
            ! 2. Compute the Array of Vertical Distances from the Bottom of the Water Column to the Top of each layer
            zfromb0 = 0.
            Do iLayer = HydroParam%ElSmallm(iElem) + 1, HydroParam%ElCapitalM(iElem)
                zfromb0 = zfromb0 + HydroParam%DZhit(iLayer,iElem)
                zfromb(iLayer) = zfromb0
            End Do      ! Layer Loop
            zfromb(HydroParam%ElSmallm(iElem))     = 0.
            zfromb(HydroParam%ElCapitalM(iElem)+1) = V(HydroParam%etan(iElem),HydroParam%hb(iElem))
            
            ! 3. Compute the Array of Vertical Distances from the Top of the Water Column to the Top of each layer
            zfromt0 = 0.
            Do iLayer = HydroParam%ElCapitalM(iElem), HydroParam%ElSmallm(iElem) + 1, -1
                zfromt0 = zfromt0 + HydroParam%DZhit(iLayer-1,iElem)
                zfromt(iLayer) = zfromt0
            End Do      ! Layer Loop
            zfromt(HydroParam%ElSmallm(iElem))     = V(HydroParam%etan(iElem),HydroParam%hb(iElem))
            zfromt(HydroParam%ElCapitalM(iElem)+1) = 0.
            
            ! 4. Compute the Lenght Scale at time n
            xl(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)) = q2lp(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem)/q2p(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem)
            xl(HydroParam%ElSmallm(iElem)) = 0.; xl(HydroParam%ElCapitalM(iElem)+1) = 0.
            
            ! ----- Turbulent Kinetic Energy (q2) Equations -----
            ! 5. Compute Shear and Buoyancy source Terms            
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem) - 1
                ! ----- Density Gradients -----
                dzp     = 0.5*( HydroParam%DZhit(iLayer+1,iElem) + HydroParam%DZhit(iLayer,iElem) )
                RhoAv   = 0.5*( HydroParam%sDRhoW(iLayer+1,iElem) + HydroParam%sDRhoW(iLayer,iElem) ) !+ 1000.
                drhodz  = ( HydroParam%sDRhoW(iLayer+1,iElem) - HydroParam%sDRhoW(iLayer,iElem) )/dzp
                ShearSq = ( (HydroParam%ub(iLayer+1,1,iElem) - HydroParam%ub(iLayer,1,iElem))/dzp )**2. + ( (HydroParam%ub(iLayer+1,2,iElem) - HydroParam%ub(iLayer,2,iElem))/dzp )**2.
                N2      = -HydroParam%g*drhodz/rhoAv
                ! ----- Shear -----
                ShearSq = ( (HydroParam%ub(iLayer+1,1,iElem) - HydroParam%ub(iLayer,1,iElem))/dzp )**2. + ( (HydroParam%ub(iLayer+1,2,iElem) - HydroParam%ub(iLayer,2,iElem))/dzp )**2.
                ! 5.1. Corrective step: Length Scale Limitation from [3]
                If (N2 > 0.0) xl(iLayer+1) = Min(xl(iLayer+1),0.53*SQRT(q2p(iLayer+1,iElem)/N2))
                ! 5.2. Calculate TKE Production at the Interfaces
                HydroParam%ShearProd(iLayer+1,iElem) = HydroParam%VerEddyViscCell(iLayer+1,iElem)*ShearSq
                HydroParam%BuoyancyProd(iLayer+1,iElem) = -HydroParam%VerEddyDiffCell(iLayer+1,iElem)*N2
            End Do
            
            ! 6. Diffusion Coefficients for q2 and q2l
            dCoeff(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)+1) = Sq*0.5*( xl(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)+1) + xl(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)))*sqrt( 0.5*(q2p(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)+1,iElem) + q2p(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem),iElem)) )
            
            ! 7. Form Tridiagonal Matrix
            ! ----- Define Upper Diagonal Terms -----
            aaTh(3,HydroParam%ElSmallm(iElem)) = 0.
            aaTh(3,HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)) = -dCoeff(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem))/( HydroParam%Dzi(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem)*0.5*( HydroParam%DZhit(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem) + HydroParam%DZhit(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1,iElem) ) )
            aaTh(3,HydroParam%ElCapitalM(iElem)+1) = 0.
            ! ----- Define Lower Diagonal Terms -----
            aaTh(1,HydroParam%ElSmallm(iElem)) = 0.
            aaTh(1,HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)) = -dCoeff(HydroParam%ElSmallm(iElem)+2:HydroParam%ElCapitalM(iElem)+1)/( HydroParam%DZhit(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1,iElem)*0.5*( HydroParam%DZhit(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem) + HydroParam%DZhit(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1,iElem) ) )
            aaTh(1,HydroParam%ElCapitalM(iElem)+1) = 0.
            ! ----- Define Main Diagonal Terms -----
            aaTh(2,HydroParam%ElSmallm(iElem)) = 3./(2.*dt)
            aaTh(2,HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)) = 3./(2.*dt) - aaTh(1,HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)) - aaTh(3,HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)) + 2.*( Sqrt(q2p(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem))/( B_1*xl(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)) ) )
            aaTh(2,HydroParam%ElCapitalM(iElem)+1) = 3./(2.*dt)
            ! ----- Define RHS Matrix -----
            dsTh(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)) = 4.*q2p(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem)/(2.*dt) - q2pp(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem)/(2.*dt)
            Do iLayer = HydroParam%ElSmallm(iElem)+1, HydroParam%ElCapitalM(iElem)
                If ( HydroParam%BuoyancyProd(iLayer,iElem) >= 0.0 ) Then         ! Stable Stratification
                    dsTh(iLayer) = dsTh(iLayer) + 2*( HydroParam%ShearProd(iLayer,iElem) + HydroParam%BuoyancyProd(iLayer,iElem) )
                Else        ! Unstable Stratification
                    dsTh(iLayer) = dsTh(iLayer) + 2*( HydroParam%ShearProd(iLayer,iElem) )
                    aaTh(2,iLayer) = aaTh(2,iLayer) - HydroParam%BuoyancyProd(iLayer,iElem)*2./q2p(iLayer,iElem)
                End If      ! If ( BuoyancyProd(iLayer) >= 0.0 ) Then
            End Do      ! Layer Loop
            
            ! 8. Boundary Conditions
            ! ----- Compute Surface and Bottom Tensions -----
            VCell(1:2,1) = (/ HydroParam%ub(HydroParam%ElCapitalM(iElem),1,iElem), HydroParam%ub(HydroParam%ElCapitalM(iElem),1,iElem) /)
            VCell(1:2,2) = (/ HydroParam%ub(HydroParam%ElCapitalM(iElem),2,iElem), HydroParam%ub(HydroParam%ElCapitalM(iElem),2,iElem) /)
            
            ! 1. Shear stress in water surface
            If (HydroParam%iWindStress==1) Then ! Formulation 1 - Based on Wind Drag coefficient
                GammaT = HydroParam%windDragConstant*(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.)**0.5
            ElseIf (HydroParam%iWindStress==0) Then ! Formulation 2 - Based on Air density
                SphAir  = SpecHum(MeteoParam%AirTemp(iElem),101325*MeteoParam%AtmPressure(iElem),MeteoParam%RelHum(iElem))
                AirDens = RhoAir(101325*MeteoParam%AtmPressure(iElem),MeteoParam%AirTemp(iElem),SphAir)
                GammaT = AirDens/HydroParam%rho0*(0.63+0.066*sqrt(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.))*0.001*sqrt((HydroParam%Windix(iElem)-HydroParam%ub(HydroParam%ElCapitalM(iElem),1,iElem))**2. + (HydroParam%Windiy(iElem)-HydroParam%ub(HydroParam%ElCapitalM(iElem),2,iElem))**2.)
            EndIf
            
            ! 2. Shear stress in the bottom
            If (HydroParam%BottomTensionFlag==0) Then
                If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
                    Chezy = HydroParam%Rug(iElem)
                ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
                    Chezy = Max(HydroParam%Pcri,V(HydroParam%etan(iElem),HydroParam%hb(iElem)))**(1./6.)/(HydroParam%Rug(iElem)+NearZero)
                ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
                    Chezy = 18.*log10(12.*Max(HydroParam%Pcri,V(HydroParam%etan(iElem),HydroParam%hb(iElem)))/(HydroParam%Rug(iElem)/30.+NearZero))
                EndIf
                GammaB = (HydroParam%g*(HydroParam%ub(HydroParam%ElSmallm(iElem),1,iElem)**2.+HydroParam%ub(HydroParam%ElSmallm(iElem),2,iElem)**2.)**0.5)/(Chezy**2.)
            EndIf
            
            uStarT = Sqrt( GammaT/HydroParam%rho0 )
            uStarB = Sqrt( GammaB/HydroParam%rho0 )
            ! ----- Bottom Boundary Condition -----
            dsTh(HydroParam%ElSmallm(iElem)) = 3.*B_1**(2./3.)*uStarB/(2.*dt)
            ! ----- Surface Boundary Condition -----
            dsTh(HydroParam%ElCapitalM(iElem)+1) = 3.*B_1**(2./3.)*uStarT/(2.*dt)
            
            ! 9. Solve Tridiagonal System for Vertical Distribution of TKE
            DIM = size(dsTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)+1),1)
            Call TridT(DIM,aaTh(1:3,HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)+1),dsTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)+1),SolTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)+1))
            
            ! 10. Define TKE at Next Time Step
            q2(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)+1,iElem) = Max( q2_min, SolTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)+1) )
            q2(HydroParam%ElCapitalM(iElem)+2:MeshParam%KMax,iElem) = q2(HydroParam%ElCapitalM(iElem)+1,iElem)
            !! ----- Define Variables for Basin Scale Energy Balances
            !If ( BasinScaleFlag == 1 ) Then
            !    Do iLayer = ElSmallmOld(iElem) + 1, ElCapitalMOld(iElem)
            !        rhoxh = sDRhoW(iLayer-1,iElem)*DZin(iLayer-1,iElem)
            !        ShearProduction    = ShearProduction    + ShearProd(iLayer,iElem)*rhoxh*Area(iElem)
            !        BuoyancyProduction = BuoyancyProduction + BuoyancyProd(iLayer,iElem)*rhoxh*Area(iElem)
            !        Dissipation        = Dissipation        + q2(iLayer,iElem)*DSqrt( q2p(iLayer,iElem) )/( B_1*xl(iLayer) )*rhoxh*Area(iElem)
            !        TKinE              = TKinE              + q2(iLayer,iElem)*rhoxh*Area(iElem)
            !    End Do      ! Layer Loop
            !Else
            !    Continue
            !End If      ! If ( BasinScaleFlag == 1 ) Then
            
            ! ----- Turbulent Length Scale (q2l) Equations -----
            ! 11. Form Tridiagonal Matrix
            ! ----- Define Upper Diagonal Terms -----
            ! -> Same as for q2
            ! ----- Define Lower Diagonal Terms -----
            ! -> Same as for q2  
            ! ----- Define Main Diagonal Terms -----
            aaTh(2,HydroParam%ElSmallm(iElem)) = 3./(2.*dt)
            Do iLayer = HydroParam%ElSmallm(iElem)+1, HydroParam%ElCapitalM(iElem)
                aaTh(2,HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)) = 3./(2.*dt) - aaTh(1,iLayer) - aaTh(3,iLayer) + Sqrt(q2p(iLayer,iElem))/( B_1*xl(iLayer) )*( 1 + E_2*( xl(iLayer)/( VonKarman*zfromb(iLayer) ) )**2. + E_3*( xl(iLayer)/( VonKarman*zfromT(iLayer) ) )**2. )
            End Do      ! Layer Loop
            aaTh(2,HydroParam%ElCapitalM(iElem)+1) = 3./(2.*dt)
            ! ----- Define RHS Matrix -----
            Select Case (iTurb)
            Case(1)         ! Original MY2.5 Formulation
                Do iLayer = HydroParam%ElSmallm(iElem)+1, HydroParam%ElCapitalM(iElem)
                    If ( HydroParam%BuoyancyProd(iLayer,iElem) >= 0. ) Then         ! Stable Stratification
                        dsTh(iLayer) = xl(iLayer)*( E_1*HydroParam%ShearProd(iLayer,iElem) + E_1*HydroParam%BuoyancyProd(iLayer,iElem) )
                    Else        ! Unstable Stratification
                        dsTh(iLayer) = xl(iLayer)*( E_1*HydroParam%ShearProd(iLayer,iElem) )
                        aaTh(2,iLayer) = aaTh(2,iLayer) - E_1*HydroParam%BuoyancyProd(iLayer,iElem)/q2p(iLayer,iElem)
                    End If      ! If ( BuoyancyProd(iLayer) >= 0. ) Then
                End Do      ! Layer Loop
            Case(2:)         ! k-kl Implementation of [4]
                Do iLayer = HydroParam%ElSmallm(iElem)+1, HydroParam%ElCapitalM(iElem)
                    If ( HydroParam%BuoyancyProd(iLayer,iElem) >= 0. ) Then         ! Stable Stratification
                        dsTh(iLayer) = xl(iLayer)*( E_1*HydroParam%ShearProd(iLayer,iElem) + 2.00*HydroParam%BuoyancyProd(iLayer,iElem) )
                    Else        ! Unstable Stratification
                        dsTh(iLayer) = xl(iLayer)*( E_1*HydroParam%ShearProd(iLayer,iElem) )
                        aaTh(2,iLayer) = aaTh(2,iLayer) - 5.06*HydroParam%BuoyancyProd(iLayer,iElem)/q2p(iLayer,iElem)
                    End If      ! If ( BuoyancyProd(iLayer) >= 0. ) Then
                End Do      ! Layer Loop
            End Select
            dsTh(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)) = dsTh(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem)) + 4.*q2lp(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem)/(2.*dt) - q2lpp(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem)/(2.*dt)
            
            ! 12. Boundary Conditions
            ! ----- Bottom Boundary Condition -----
            dsTh(HydroParam%ElSmallm(iElem)) = 0.0
            ! ----- Surface Boundary Condition -----
            dsTh(HydroParam%ElCapitalM(iElem)+1) = 0.0
            
            ! 13. Solve Tridiagonal System for Vertical Distribution of LengthScale
            DIM = size(dsTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)+1),1)
            Call TridT(DIM,aaTh(1:3,HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)+1),dsTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)+1),SolTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)+1))
            
            ! 14. Define LengthScale at Next Time Step
            q2l(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)+1,iElem) = Max( q2l_min, SolTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)+1) )
            q2l(HydroParam%ElCapitalM(iElem)+2:MeshParam%KMax,iElem) = q2l(HydroParam%ElCapitalM(iElem)+1,iElem)
            
        End If      ! If ( ElSmallmOld(iElem) == ElCapitalMOld(iElem) ) Then
    End Do      ! Element Loop
    
    !If ( BasinScaleFlag == 1 ) Then
    !    ShearProduction    = ShearProduction*dt
    !    BuoyancyProduction = BuoyancyProduction*dt
    !    Dissipation        = Dissipation*dt
    !    TKinE              = TKinE*0.5
    !Else
    !    Continue
    !End If      ! If ( BasinScaleFlag == 1 ) Then
            
    Return
End Subroutine SolveTKE_LengthScale  