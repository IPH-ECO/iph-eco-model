
!>@brief Variable declation for hydrodynamic computation
!>@author Rafael Cavalcanti
  !>@attention List of modification
  !>@parblock
  !! -> 10.03.2014: Routine Implementation (Rafael Cavalcanti)
  !!
  !! -> 15.09.2015: Routine Update         (Carlos Ruberto Fragoso)
  !>@endparblock
 
  Module Hydrodynamic
    
   use domain_types
    
    Implicit None
    
    type HydrodynamicParam
        Real:: Pi = 3.141592653589793238462d0 !< Pi number
        Real:: SmallNumber = 1e-5 !< Small Number
    
        ! 2. Vertical Discretization Variables
        !Real:: zL !< lowest reference elevation (needs to be lower than the max elevation found in geometry mesh)
        !Real:: zR !< greatest reference elevation (needs to be greater than the max elevation found in geometry mesh)
        Real, Allocatable:: Z(:,:) !< Elevation in the bottom of the layer in each edge
        Real, Allocatable:: Ze(:,:) !< Elevation in the bottom of the layer in each element  
        Real, Allocatable:: Zb(:,:)  !< Elevation in the center of the layer in each element  
        Real, Allocatable:: DZj(:,:) !< Height of the layer in each edge in the time n+1
        Real, Allocatable:: DZjt(:,:) !< Height of the layer in each edge in the time n
        Real, Allocatable:: DZi(:,:) !< Height of the layer in each element in the time n+1
        Real, Allocatable:: DZit(:,:) !< Height of the layer in each element in the time n
        Real, Allocatable:: iADZ(:,:)
        Real, Allocatable:: iAG(:,:)
        Real, Allocatable:: DZiADZ(:)
        Real, Allocatable:: DZiAG(:)
        Integer, Allocatable::  Smallm(:)   !<Lower layer index at the edge of the element
        Integer, Allocatable::  Smallms(:)   !<Lower layer index at the edge of the element
        Integer, Allocatable :: CapitalM(:) !<greatest layer index at the edge of the element
        Integer, Allocatable :: CapitalMs(:) !<greatest layer index at the edge of the element !CAYO
        Integer, Allocatable :: ElSmallm(:) !<Lower layer index at the center of the element
        Integer, Allocatable :: ElSmallms(:) !<Lower layer index at the center of the element
        Integer, Allocatable :: ElCapitalM(:) !<greatest layer index at the center of the element
        Integer, Allocatable :: ElCapitalMs(:) !<greatest layer index at the center of the element !CAYO
        
       
        Real, Allocatable:: DZK(:) !Sediment Layer
        
        Real, Allocatable:: DZsj(:,:)!CAYO
        Real, Allocatable:: DZsjt(:,:)!CAYO
        Real, Allocatable:: DZhj(:,:)!CAYO
        Real, Allocatable:: DZhjt(:,:) !CAYO
        Real, Allocatable:: DZsi(:,:)!CAYO
        Real, Allocatable:: DZsit(:,:)!CAYO
        Real, Allocatable:: DZhi(:,:)!CAYO
        Real, Allocatable:: DZhit(:,:) !CAYO   
    
        ! 3. Hydrodynamic        
        Real:: SumVerAcum =0
        Real:: SumVer
        ! 3.1. Velocities
        Real, Allocatable:: u(:,:) !< Normal velocity at the edges in each layer, dimension: Kmax,nEdge
        Real, Allocatable:: ut(:,:) !< Normal velocity at previous time step, dimension: Kmax,nEdge
        Real, Allocatable:: utang(:,:) !< Normal velocity at the edges in each layer, dimension: Kmax,nEdge
        Real, Allocatable:: uxyback(:,:,:) !< horizontal backtracking velocity components at the edges in each layer, dimension: Kmax,2,nEdge
        Real, Allocatable:: uNode(:,:,:) !< Nodal velocity , dimension: Kmax,2,nNode
        Real, Allocatable:: uxy(:,:,:) !< horizontal velocity components at current time step ,dimension: (Kmax,2,nEdge)
        Real, Allocatable:: uxyL(:,:,:) !< horizontal velocity components in the center of each k+1/2 Layer for eack Element ,dimension: (Kmax+1,2,nElem)
        Real, Allocatable:: Wu(:,:) !< horizontal velocity components at current time step ,dimension: (Kmax,nEdge)
        Real, Allocatable:: ug(:,:),vg(:,:),wg(:,:) !< Velocities in each face in k+1/2 Layer dimension: (nEdge,Kmax+1)
        Real, Allocatable:: ub(:,:,:)  !< cell-centered three components of velocity, dimension: ub(Kmax,3,nElem)
        Real, Allocatable:: ubV(:,:,:)  !< Vertice center three components of velocity, dimension: ubV(Kmax,3,nNode)
        Real, Allocatable:: ubBack(:,:,:) !< backtracking velocity components at the center of each Element in each layer, dimension: Kmax,2,nElem
        Real, Allocatable:: w(:,:) !< vertical velocity dimension: (Kmax+1,nElem)
        Real, Allocatable:: wt(:,:) !< vertical velocity dimension at previous time step: (Kmax+1,nElem)
        Real, Allocatable:: wfc(:,:) !< Face-centered vertical velocity dimension: (Kmax,nEdge)
        Real, Allocatable:: FuxyNode(:,:,:)
        Real, Allocatable:: Fub(:,:,:)
        Real, Allocatable:: Fw(:,:) !< vertical velocity dimension at previous time step: (Kmax+1,nElem)
        Real, Allocatable:: epson(:,:) !< Normal velocity at the edges in each layer, dimension: Kmax,nEdge
        Real, Allocatable:: psi_edge(:,:)
        Real, Allocatable:: psi_cell(:,:)
        
        Real, Allocatable:: us(:,:) !< Normal superficial flow velocity at the edges in each layer, dimension: Kmax,nEdge !CAYO
        Real, Allocatable:: ust(:,:) !< Normal superficial flow  velocity at previous time step, dimension: Kmax,nEdge !CAYO
        Real, Allocatable:: um(:,:) !<  Kmax,nEdge !CAYO
        Real, Allocatable:: umt(:,:) !< Kmax,nEdge    !CAYO
        ! 3.1. Others Variables
        Real, Allocatable:: etaInf(:) !< Tidal boundary condition
        Real, Allocatable:: etaplus(:) !< Vertical water balance at current time step  dimension: nElem
        Real, Allocatable:: peta(:) !< Nodal free-Surface Elevation at current time step  dimension: nNode
        Real, Allocatable:: petan(:)     ! Nodal Free-Surface Elevation from previous timestep n (Time Step N+1; Time Step N) dimension: nElem
        Real, Allocatable:: eta(:) !< Cell-centered Free-Surface Elevation at current time step  dimension: nElem
        Real, Allocatable:: etan(:)     ! Cell-centered Free-Surface Elevation from previous timestep n (Time Step N+1; Time Step N) dimension: nElem
        Real, Allocatable:: hb(:) !< elevation at the center of each element dimension: nElem
        Real, Allocatable:: sb(:)
        Real, Allocatable:: H(:) !< depth of the Edge dimension: nEdge
        Real, Allocatable:: hj(:)  !< elevation of the edge dimension: nEdge
        Real, Allocatable:: sj(:)  !< elevation of the edge dimension: nEdge
        Real, Allocatable:: Hu(:) !< vertically integrated velocity, dimension: nEdge
        Real, Allocatable:: P(:)  !< area (term P in casulli 2000), dimension: nElem
        Real, Allocatable:: Aeta(:) !< Volume in the element
        Real, Allocatable:: f(:) !< vector P in casulli 2000
        Real, Allocatable:: Deta(:), rhs(:), Gu(:,:), Fu(:,:), Fv(:,:), Fvu(:,:), Fuv(:,:)
        Real, Allocatable:: rhsnonHydro(:,:)
        Real, Allocatable:: q(:,:),pq(:,:)
        Real, Allocatable:: Rug(:) !< Roughness coefficient
        Real, Allocatable:: uIniVec(:,:) !< Initial condition of velocity components
        Real, Allocatable:: sDRhoW(:,:) !<Water density
        Real, Allocatable:: sDRhoWt(:,:) !<Water density
    
        Real, Allocatable:: Hs(:) !CAYO
        !Real, Allocatable:: Kj(:,:), ei(:,:) !Porosity and Hydraulic Conductivity   !CAYO !MOD_Mesh
        Real, Allocatable:: Vol(:) !<Water Volume in element with porous region             !CAYO 
        
        ! 4. Turbulence Model    
        Real,Allocatable:: HorViscosity(:,:,:)      !< Horizontal Eddy Viscosity
        Real,Allocatable:: HorDiffusivity(:,:,:)    !< Horizontal Eddy Diffusivity
        Real,Allocatable:: VerEddyVisc(:,:)         !< Vertical Eddy Viscosity Edge
        Real,Allocatable:: VerEddyDiff(:,:)         !< Vertical Eddy Diffusivity Edge
        Real,Allocatable:: VerEddyViscCell(:,:)         !< Vertical Eddy Viscosity Cell-centered
        Real,Allocatable:: VerEddyDiffCell(:,:)         !< Vertical Eddy Diffusivity Cell-centered
        Real, Allocatable:: TKE(:,:), TKEP(:,:)                  ! Total Kinetic Energy
        Real, Allocatable:: LengthScale(:,:), LengthScaleP(:,:)  ! Turbulent Eddy Length Scale
        Real, Allocatable:: DissipRate(:,:)                      ! Dissipation Rate of Total Kinetic Energy
        Real, Allocatable:: ShearProd(:,:), BuoyancyProd(:,:)    ! Dissipation Rate of Total Kinetic Energy
        
    
        ! 5. Pressure
        Real,Allocatable:: PBarc(:,:)         !< baroclinic pressure contribution
    
        ! 6. Boundary Conditions
        Real, Allocatable:: WindVel(:,:)
        Real, Allocatable:: WindXY(:,:)
        Real, Allocatable:: Windix(:),Windiy(:)
        Real, Allocatable:: InFlowValue(:,:),WaterLevelValue(:,:)
        Real, Allocatable:: InFlowTime(:,:),WaterLevelTime(:,:)
        Integer, Allocatable:: InFlownTime(:),WaterLevelnTime(:)
        Integer, Allocatable:: InFlowSmallm(:),InFlowCapitalM(:)
        Integer, Allocatable:: IndexWaterLevel(:,:)
        Integer, Allocatable:: IndexWaterLevelEdge(:)
        Integer, Allocatable:: IndexInflow(:,:)
        Integer, Allocatable:: IndexInflowEdge(:)
        Real, Allocatable:: WaterLevel(:)
        Integer:: NInflow, NWaterLevel
        Integer, Allocatable:: NRange(:)
    
        ! 7. Hydrodynamic Parameters 
        Integer:: NFUT !<Numbers of Sub-time steps for ELM (used in FUFV)
        Integer:: NTRASP !<Maximum Numbers of Sub-time steps for transport solution (used in FUFV)
        Real:: Lat !<average latitude
        Real:: Altit !<average altitude 
        Real:: ALB  !<Albedo
        Real:: OMEGA !<Earth angular velocity (in rad/sec)
        Real:: g !<Acceleration due to gravity (m/s²)
        Real:: CFL  !<Courant–Friedrich–Lewy Number
        Real:: nux, nuy, nuz
        Real:: GammaB !<Tension in the bottom layer
        Real:: GammaT !<Tension in the surface layer
        Real:: Theta !<Implicitness coefficient (from 0.5 to 1.00)
        Real:: CSmag !<Smagorinsky coefficient (from 0.10 to 0.20)
        Real:: HorEddyDiffY_Back !<Background horizontal Diffusivity
        Real:: HorEddyViscX_Cte !<Horizontal Viscosity (x-direction)
        Real:: HorEddyViscY_Cte !<Horizontal Viscosity (y-direction)
        Real:: VerEddyVisc_Cte !<Vertical Viscosity (constant)
        Real:: HorEddyDiffX_Cte !<Horizontal Diffusivity (x-direction)
        Real:: HorEddyDiffY_Cte !<Horizontal Diffusivity (y-direction)
        Real:: VerEddyDiff_Cte !< Vertical Diffusivity (constant)
        Real:: alfa_turbmodel1, n_turbmodel1 !<Vertical mixing turbulence model (VerTurbFlag == 1)
	    Real:: vref !<Vertical eddy Viscosity of Reference (VerTurbFlag == 1)
	    Real:: vmin !<Background Vertical eddy Viscosity (VerTurbFlag == 1)
        Real:: tdmin_pp !<Background Vertical eddy diffuvsivity (VerTurbFlag == 1)
	    Real:: rho0 !<Water density of reference (kg/m³)
        Real:: WtempRef !<Water temperature of reference (oC)
        Real:: AirtempRef !<Air temperature of reference (oC)
        Real:: Pcri !<threshold Depth for dry/wet algorithm
        Real:: smoothingTime !<Smoothing period to eliminate effects of the initial conditions 
        Real:: windDragConstant  !<wind Drag coefficient (constant)
        Real:: windDragCoefficient(3) !<wind Drag coefficient (linear function)
        Real:: windDragWindSpeed(3) !<Wind Speed threshold for changing wind Drag coefficient (linear function)
        Real:: RugChezyConst !<Roughness coeficient of Chezy (Constant)
        Real:: RugManConst !<Roughness coeficient of Manning (Constant)
        Real:: RugWCConst !<Roughness coeficient of White-Colebrook (Constant)
        !Real, Allocatable:: rhoair(:) !<Air density (kg/m³)
        !Fetch
        Real, Allocatable:: fetch_m(:,:) !<Fetch distance (m)
    
    
        !8. Initial Condition
        Real:: Uini !<Initial condition of x velocity component (m/s)
        Real:: Vini !<Initial condition of y velocity component (m/s)
        Real:: Wini !<Initial condition of z velocity component (m/s)
        Real:: Zini !<Initial condition of elevation surface water (m)
    

        ! 9. Flags
        Integer:: OutPutFlag !< Output flag: OutPutFlag = 0 (to be defined); OutPutFlag = 1 (VKT output)
        !Integer:: ConvectiveFlag (iConv)
        !Integer:: HorTurbFlag (iHTurb)
        !Integer:: VerTurbFlag (iVTurb)
        !Integer:: SurfTensionFlag
        Integer:: BottomTensionFlag    
        !Integer:: ResuspMethodFlag (vai para modulo de limnologia)
        !Integer:: SetMethodFlag (vai para modulo de limnologia)
        Integer:: iWindStress ! (SurfTensionFlag)< Wind Stress on Water Surface formulation: iWindStress = 0 (Based on Wind Drag coefficient); iWindStress = 1 (Based on Air density)
        Integer:: iWindDrag !< Wind Drag coefficient formulation: iWindDrag = 0 (Constant); iWindDrag = 1 (Linear function)
        Integer:: iRoughForm !< Roughness Formulation: iRoughForm = 0 (roughnessChezyConstant); iRoughForm = 1 (roughnessManningConstant); iRoughForm = 2 (roughnessWhiteColebrookConstant); iRoughForm = 3 (roughnessChezyUseGridData); iRoughForm = 4 (roughnessManningUseGridData); iRoughForm = 5 (roughnessWhiteColebrookUseGridData);
        Integer:: iHTurb  !< Horizontal turbulence Formulation: iHTurb = 0 (constant); iHTurb = 1 (Smagorinsky model)
        Integer:: iVTurb !< Vertical turbulence Formulation: iVTurb = 0 (constant); iVTurb = 1 (zero model)
        Integer:: iBTurb !< Bottom turbulence Formulation (not implemented yet)
        Integer:: iConv !< Formulation for convective terms: iConv = 0 (Neglect Nonlinear Convection); iConv = 1 (Eulerian-Lagragean Method)
        Integer:: iBarot !< Barotropic effect: iBarot = 0 (no); iBarot = 1 (yes)
        Integer:: iNonHydro !< Non hydrostatic pressure effect: iNonHydro = 0 (no); iNonHydro = 1 (yes)
        Integer:: iCoriolis !< Coriolis effect: iCoriolis = 0 (no); iCoriolis = 1 (yes)
        
        Real, Allocatable:: SScalar(:)               ! Salinity Concentration - Current Time Step
        Real, Allocatable:: SScalar2D(:)               ! Salinity Concentration - Current Time Step
        Real, Allocatable:: SVector(:,:)               ! Salinity Concentration - Current Time Step
       
        
        ! Transport Equation Variables
        Real, Allocatable:: uLoadVarEst(:,:)
        Real, Allocatable:: dVarEst(:,:,:)
        Real, Allocatable:: Locdt(:,:)
        Real, Allocatable:: DZitau(:,:)
        Real, Allocatable:: DZistau(:,:)
        
    contains
        procedure :: initializeHydroParam
    end type
    
    contains
    
    subroutine initializeHydroParam(this, hydroConfiguration)
    
        Integer:: i,j
        Real:: LATR
        !character, pointer :: hydroParametersName(:)
        type(HydrodynamicConfiguration) :: hydroConfiguration  
        type(HydrodynamicParameter), pointer :: hydroParameters(:)
        character(len=200):: text
        class(HydrodynamicParam) :: this
    
        call c_f_pointer(hydroConfiguration%parameters, hydroParameters, [hydroConfiguration%numberOfParameters])
    
        Do i = 1, hydroConfiguration%numberOfParameters
            text = trim(hydroParameters(i)%name)

            !Flags
            If (trim(text) == 'iWindStress') Then
                this%iWindStress = hydroParameters(i)%value
            ElseIf (trim(text) == 'iWindDrag') Then
                this%iWindDrag  = hydroParameters(i)%value
            ElseIf (trim(text) == 'bottomRoughness') Then
                this%iRoughForm = hydroParameters(i)%value
            ElseIf (trim(text) == 'horizontalEddyVD') Then
                this%iHTurb = hydroParameters(i)%value
            ElseIf (trim(text) == 'verticalEddyVD') Then
                this%iVTurb = hydroParameters(i)%value
            ElseIf (trim(text) == 'iBTurb') Then
                this%iBTurb = hydroParameters(i)%value
            ElseIf (trim(text) == 'iConv') Then
                this%iConv = hydroParameters(i)%value
            ElseIf (trim(text) == 'pressure') Then
                this%iBarot = hydroParameters(i)%value
            ElseIf (trim(text) == 'hydrostaticComponent') Then
                this%iNonHydro = hydroParameters(i)%value
            ElseIf (trim(text) == 'iCoriolis') Then
                this%iCoriolis = hydroParameters(i)%value
            EndIf
            !Parameters
            If (trim(text) == 'waterDensity') Then
                this%rho0 = hydroParameters(i)%value
            ElseIf (trim(text) == 'waterTemperature') Then
                this%WtempRef  = hydroParameters(i)%value
            ElseIf (trim(text) == 'airTemperature') Then
                this%AirtempRef  = hydroParameters(i)%value
            ElseIf (trim(text) == 'thetaCoefficient') Then
                this%theta  = hydroParameters(i)%value
            ElseIf (trim(text) == 'thresholdDepth') Then
                this%Pcri  = hydroParameters(i)%value
            ElseIf (trim(text) == 'smoothingTime') Then
                this%smoothingTime  = hydroParameters(i)%value
            ElseIf (trim(text) == 'windDragConstant') Then
                this%windDragConstant  = hydroParameters(i)%value
            ElseIf (trim(text) == 'windDragCoefficientC1') Then
                this%windDragCoefficient(1)  = hydroParameters(i)%value
            ElseIf (trim(text) == 'windDragCoefficientC2') Then
                this%windDragCoefficient(2)  = hydroParameters(i)%value
            ElseIf (trim(text) == 'windDragCoefficientC3') Then
                this%windDragCoefficient(3)  = hydroParameters(i)%value
            ElseIf (trim(text) == 'windDragWindSpeedW1') Then
                this%windDragWindSpeed(1)  = hydroParameters(i)%value
            ElseIf (trim(text) == 'windDragWindSpeedW2') Then
                this%windDragWindSpeed(2)  = hydroParameters(i)%value
            ElseIf (trim(text) == 'windDragWindSpeedW3') Then
                this%windDragWindSpeed(3)  = hydroParameters(i)%value
            ElseIf (trim(text) == 'roughnessChezyConstant') Then
                this%RugChezyConst  = hydroParameters(i)%value
            ElseIf (trim(text) == 'roughnessManningConstant') Then
                this%RugManConst  = hydroParameters(i)%value
            ElseIf (trim(text) == 'roughnessWhiteColebrookConstant') Then
                this%RugWCConst  = hydroParameters(i)%value
            ElseIf (trim(text) == 'horizontalEddyViscosity') Then
                this%HorEddyViscX_Cte  = hydroParameters(i)%value
                this%HorEddyViscY_Cte = this%HorEddyViscX_Cte
            ElseIf (trim(text) == 'horizontalEddyDiffusivity') Then
                this%HorEddyDiffX_Cte  = hydroParameters(i)%value
                this%HorEddyDiffY_Cte  = this%HorEddyDiffX_Cte
            ElseIf (trim(text) == 'smagorinskyCoefficient') Then
                this%CSmag  = hydroParameters(i)%value
            ElseIf (trim(text) == 'backgroundHorizontalEddyDiffusivity') Then
                this%HorEddyDiffY_Back  = hydroParameters(i)%value
            ElseIf (trim(text) == 'verticalEddyViscosity') Then
                this%VerEddyVisc_Cte  = hydroParameters(i)%value
            ElseIf (trim(text) == 'verticalEddyDiffusivity') Then
                this%VerEddyDiff_Cte  = hydroParameters(i)%value
            ElseIf (trim(text) == 'refVerticalEddyViscosity') Then
                this%vref  = hydroParameters(i)%value
            ElseIf (trim(text) == 'backgroundVerticalEddyViscosity') Then
                this%vmin  = hydroParameters(i)%value
            ElseIf (trim(text) == 'backgroundVerticalEddyDiffusivity') Then
                this%tdmin_pp  = hydroParameters(i)%value
            ElseIf (trim(text) == 'numberOfSubTimeSteps') Then
                this%NFUT  = hydroParameters(i)%value
            ElseIf (trim(text) == 'earthsRotationSpeed') Then
                this%OMEGA  = hydroParameters(i)%value
            ElseIf (trim(text) == 'latitude') Then
                this%Lat  = hydroParameters(i)%value
            ElseIf (trim(text) == 'atitude') Then
                this%Altit  = hydroParameters(i)%value
            EndIf
            !print *, hydroParametersName, hydroParameters(i)%value
        EndDo
    
        this%alfa_turbmodel1 = 5.d0
        this%n_turbmodel1 = 1.d0
        this%NTRASP = 10.d0
    
        !Gravity of the Earth (m/s²)
        this%g = 9.810665d0 !9.780327*(1. + 0.0053024*(sin(this%Lat*this%Pi/180))**2. - 0.0000058*sin(2*this%Lat*this%Pi/180)**2. ) - 3.086e-6*this%Altit
    
        !Reference water density (kg/m³)
        this%rho0 = 1000.d0
    
        !Albedo
	    LATR = 180.d0*(ASIN(SIN(ABS(this%LAT)/1.33d0*this%PI/180.d0)))/this%PI
	    this%ALB  = (((COS(ABS(this%LAT)*this%PI/180.d0)-COS(LATR*this%PI/180.d0)*1.33d0)/(COS(ABS(this%LAT)*this%PI/180.d0)+COS(LATR*this%PI/180.d0)*1.33d0))**2.+((COS(LATR*this%PI/180.d0)-COS(ABS(this%LAT)*this%PI/180.d0)*1.33d0)/(COS(LATR*this%PI/180.d0)+COS(ABS(this%LAT)*this%PI/180.d0)*1.33d0))**2.d0)/2.d0    
        
    end subroutine
        
End Module Hydrodynamic