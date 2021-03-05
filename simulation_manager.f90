!>\mainpage 
  !>\section introsecIPH Introduction
  !! IPH-ECO is a three-dimensional complex dynamic model of aquatic ecosystems such as lakes,
  !! reservoirs and estuaries. The model consists of several partial differential equations being
  !! divided in three modules: (a) a detailed hydrodynamic module, describing quantitative flows
  !! and water level; (b) a nutrient module, which deals with nutrient transport mechanisms;
  !! and (c) a biological module, which describes whole aquatic food-web interactions.
  !! The differential equations are solved numerically by applying an efficient semi-implicit finite differences
  !! method in a three-dimensional regular grid.
  !! For an extended model like IPH-ECO it is important think about how to control the complexity with the numerous parameters.
  !! An effective way to control complexity is to implement flexibility in the model design in such way that it is possible to
  !! switch between different modes of complexity. Therefore IPH-ECO includes a graphical user-friendly interface (GUI)
  !! for MS Windows environment with a flexible design to vary the complexity of the model. The GUI was developed in Visual Basic
  !! which is an event-driven programming language with an integrated development environment (IDE) for its Component Object Models (COM),
  !! making it easier the building of interfaces for applications in any OS. As the simulation model was implemented in Visual FORTRAN,
  !! we connected the two modules through a FORTRAN Dynamic Link Library (DLL).
  !>\section sec_IPH_ECO_MOD Model principles
  !! text 
  !> \section sec_IPH_ECO_MOD_REF Principal references
  !> \subsection sec_IPH_ECO_MOD_ref_hydro Hydrodynamics
  !! text
  !> \subsection sec_IPH_ECO_MOD_transp Transport
  !! text
  !> \subsection sec_IPH_ECO_MOD_wq Water Quality
  !! text
  !>\image html IPH_ECo_principles.jpg "Figure 1: Distinction between lake and wetlands" width=10cm
  
  !>@file simulation_manager.f90
  !> @brief IPH-ECO: A dynamic complex tridimensional water quality model
  !> @details This model is based on 
  !> @author Development team: Ecotecnology and Ecological Modelling group (IPH/UFRGS)
  !> @author main contributors:  Carlos Ruberto Fragoso Júnior and  J. Rafael Cavalcanti
  !> @version (Current version: 1.2)  
  !> \note Main differences from Interface Official code:
  !> @parblock
  !> -#   Re-organization of the code
    !> -#  English language in comments
    !> -#  Wet-and Dry-ing algorithm for water quality model
    !! -#  Simplification of main 
  !>\attention List of modification
  !> -# Implementation Ruberto C Fragoso
  !> -# ...
    !> @endparblock
  !>\section sec_usage Usage
  !! If launched from command line, the iph_ECO program requires 7 arguments which are the paths and file names that will be read before starting the time loop
  !! The list of arguments is:
  !> -# name of the simulation reference file
  !> -# name of the boundary condition matrix file
  !> -# name of the geometry file
  !> -# path to the boundary condition file
  !> -# path to the environmental variable file
  !> -# name of the palette file
  !> -# path were output files will be saved 
  !>\section Flowchart
  !>\image html main_flowchart.jpg "Main flowchart version 1.3 " width = 6cm
  !> \section sec_dep Dependencies:
  !> -# MODULE Simulation
  !> -# MODULE Hydrodynamics
  !> -# MODULE Mesh

    
    subroutine startSimulation(sim) bind(C, name="startSimulation")
	!DEC$ ATTRIBUTES DLLEXPORT :: startSimulation
	Use domain_types
    Use SimulationModel
    Use Hydrodynamic
    Use MeshVars
    Use Meteorological
    Use LimnologyVars
    Use ParticleTracking
    !$ use omp_lib
    
    type(Simulation) :: sim
    type(HydrodynamicConfiguration), pointer :: hydroConfiguration
    type(WaterQualityConfiguration), pointer :: wqConfiguration
    type(StructuredMesh), pointer :: StructuredMeshFeatures
    type(GridDataConfiguration), pointer :: GridDataConfig
    type(GridData), pointer :: GridDataValues
	type(BoundaryCondition), pointer :: boundaryConditions
    type(MeteorologicalConfiguration), pointer :: MeteorConfiguration
    type(SimulationParam) :: simParam
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam
    type(LimnologyParam) :: LimnoParam
    
    
    !type(ParticleParam) :: PartParam
    

    character, pointer :: OutputNameFile(:)
    Integer i, icell, iOutput, iNewton,innerNewton
    Character(200):: FileName,Basename
    Real::  stime,ftime,t1, Simtime
    real :: start_time, finish_time
    
    !call cpu_time(start)
    
	call c_f_pointer(sim%hydrodynamicConfiguration, hydroConfiguration) 
    call c_f_pointer(hydroConfiguration%gridDataConfiguration, GridDataConfig)
    call c_f_pointer(hydroConfiguration%boundaryConditions, boundaryConditions)
    call c_f_pointer(GridDataConfig%structuredMesh, StructuredMeshFeatures)
    call c_f_pointer(sim%meteorologicalConfiguration, MeteorConfiguration)
    call c_f_pointer(sim%label, OutputNameFile, [sim%labelLength])
    call c_f_pointer(sim%waterQualityConfiguration, wqConfiguration)
    !call c_f_pointer(sim%recoveryVariables, SaveVariables)
    
    simParam%OutputFile = ""
    Do i=1,sim%labelLength
        simParam%OutputFile = trim(simParam%OutputFile)//OutputNameFile(i)
    EndDo
    simParam%outputTimeStep = sim%outputTimeInterval
    simParam%recoveryTimeStep = sim%autosaveTimeInterval
    simParam%nOutput = sim%wqoOutputParametersLength
        
    isStructuredMesh = GridDataConfig%isStructured
    If (isStructuredMesh==.false.) Stop 'Non-structured grid. Please check the grid.'

    !1. Reading simulation parameters
    Call simParam%initializeSimulationParam(sim)
    
    !2. Reading hydrodynamic parameters 
    !Call ReadHydroParameters(hydroConfiguration)
    Call HydroParam%initializeHydroParam(hydroConfiguration) 
    
    !3. Reading Mesh Features
    !Call ReadMesh(StructuredMeshFeatures)
    Call MeshParam%initializeMeshParam(sim, StructuredMeshFeatures)
    
    !4. Allocating hydrodynamic, Water Quality and Meteorologic variables
    Call AllocateHydroVars(HydroParam,MeshParam)
    
    If (simParam%iSimWQ == 1) Call LimnoParam%initializeLimnoParam(wqConfiguration)
    
    If (simParam%iSimWQ == 1) Call AllocateWQVars(LimnoParam,MeshParam)
    
    Call AllocateMeteoVars(MeshParam,MeteoParam) 
    
    !5. Reading Grid Data
    Call ReadGridData(GridDataConfig,MeshParam,HydroParam)
    
    !6. Reading hydrodynamic and water quality initial conditions
    Call ReadHydroIniCond(HydroParam,hydroConfiguration,simParam,MeshParam)
    HydroParam%g = 9.80665
    
    If (simParam%iSimWQ == 1) Then
        Call ReadWQIniCond(HydroParam,MeshParam,LimnoParam,simParam)
    Else
        LimnoParam%sDTempW = HydroParam%WtempRef
        LimnoParam%sDSal  = 0.
    EndIf  
    
    !7. Reading hydrodynamic boundary conditions
    Call ReadHydroBoundaryCondition(HydroParam,hydroConfiguration,simParam%IniTime,simParam%FinalTime,MeshParam)
    
    !7.1 Get Pressure Boundary Condition intial values to initialize etaInf in time n-1:
    Do i = 1,HydroParam%NWaterLevel
        HydroParam%WaterLevel(i) = HydroParam%WaterLevelValue(i,1)
    EndDo  
    Do iElem = 1,MeshParam%nElem  
        Do iEdge = 1,4
            Face = MeshParam%Edge(iEdge,iElem)
            l = MeshParam%Left(Face) 
            r = MeshParam%Right(Face)
            If (r == 0) Then
                !If the Face has a pressure boundary condition, the etaInf in time n-1 is set:
                If (HydroParam%IndexWaterLevelEdge(Face)>0) Then
                    If ((HydroParam%WaterLevel(HydroParam%IndexWaterLevelEdge(Face))-HydroParam%hj(Face))<HydroParam%PCRI/2.d0+NearZero) Then 
                        HydroParam%etaInfn(iElem) = HydroParam%eta(iElem)
                        HydroParam%etaInf(iElem) = HydroParam%eta(iElem)
		            Else
			            HydroParam%etaInfn(iElem) = HydroParam%WaterLevel(HydroParam%IndexWaterLevelEdge(Face))
			            HydroParam%etaInf(iElem) = HydroParam%WaterLevel(HydroParam%IndexWaterLevelEdge(Face))
                    EndIf  
                EndIf                
            EndIf
        EndDo
    EndDo
    
    
    If (simParam%iSimWQ == 1)  Call ReadWQBoundaryCondition(MeshParam,HydroParam,LimnoParam,hydroConfiguration,wqConfiguration,simParam%IniTime,simParam%FinalTime)
    
    !8. Reading meteorological data    
    !Call ReadMeteorologicalData(MeteorConfiguration,simParam%IniTime,simParam%FinalTime,HydroParam%AirtempRef,HydroParam%pi)
    Call MeteoParam%initializeMeteoParam(MeteorConfiguration,simParam%IniTime,simParam%FinalTime,HydroParam%AirtempRef,HydroParam%pi)

    !9. Reading model outcomes    
    Call ReadOutputs(sim,simParam,MeshParam%Kmax,MeshParam%nEdge,MeshParam%nElem) 
    
    !OPEN(1,FILE= trim(simParam%OutputPath)//'\'//'CondIniSal.txt')
    !READ (1, *)  (LimnoParam%sDSal(:,i), i = 1, MeshParam%nElem)
    !LimnoParam%sDSalP = LimnoParam%sDSal
    
    !10. Calculates the lake fetch for each computational cell for structured grids (eight directions)
    Call Fetch(MeshParam,HydroParam)
    
    !Call PartParam%InitializeParticle(HydroParam,MeshParam) 
            !icell = 1 !PAC
            !If (simParam%TIME == simParam%IniTime) Then
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(90,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
    !11. Dynamic block
    !Open(1001, File='level.txt', Action='Write', Status='Replace')
    !stime = omp_get_wtime()
    
    
    Allocate(simParam%wqoutput(MeshParam%nElem*MeshParam%KMax*simParam%nOutput))
    
    Do While (simParam%time < simParam%NumOfDays*86400+dble(simParam%IniTime))
    !Do While (simParam%it <= simParam%FinalTime/simParam%dt)
        simParam%it = simParam%it + 1 
        call cpu_time(start_time)
        !If (simParam%it>=413) Then
        !    Continue
        !EndIf        
       
        !Getting meteorogical values at current step time
        Call GetMeteorologicalData(MeshParam,MeteoParam,simParam%time,HydroParam%AirtempRef)
        Simtime = simParam%dt*simParam%it
        ! Use Casulli's Semi-Implicit Solution (TRIM) to solve the Depth-Averaged Navier-Stokes Equation
        If (simParam%iSimHydro == 1) Call Hydro(HydroParam,MeshParam,MeteoParam,simParam%dt,simParam%time,Simtime,iNewton,innerNewton)
        
        If (simParam%iSimWQ == 1) Call GetWQBoundaryConditions(HydroParam,LimnoParam,simParam%time)        
        
        If (simParam%iSimWQ == 1) Call Limnology(HydroParam,MeshParam,MeteoParam,LimnoParam,simParam)
        
        ! Water Density as a Function of Water Temperature and Salinity
        If (simParam%iSimWQ == 1) Then
        HydroParam%sDRhoWt = HydroParam%sDRhoW 
        Call WaterDensity(HydroParam,MeshParam,LimnoParam)
        
        Else
            HydroParam%sDRhoW    = HydroParam%rho0
        EndIf
        
        !Call ParticlePosition(PartParam,HydroParam,MeshParam,simParam%dt)
        
        !Call ResidenceTime(PartParam,HydroParam,MeshParam,simParam%dt,simParam%time,simParam%IniTime)
           
        ! Spatial Model Outcomes (VTK)
        call cpu_time(finish_time) 
        simParam%start=start_time
        simParam%finish=finish_time
        
        !Call WriteOutputs(simParam,HydroParam,MeshParam,LimnoParam,MeteoParam,iNewton,innerNewton)
        Call WriteOutputs(simParam,HydroParam,MeshParam,LimnoParam,MeteoParam,iNewton,innerNewton)

        Call SaveRecovery(sim,simParam,MeshParam,HydroParam,LimnoParam)
        
        simParam%time = simParam%time + simParam%dt 
        !Write(1001,'(i20.10,4F30.20)') time, eta(1), SumVer, SumVerAcum
        
        !Update progress of the simulation
        sim%progress = int(((simParam%time-simParam%IniTime)/(simParam%NumOfDays*86400))*100)
        
        !Finish or abort a simulation
        If (sim%statusCode == 4.or.sim%statusCode == 5) exit
        
        !Pause and resume a simulation according to simulation status
        do while (sim%statusCode == 3)
            call sleep(1)
            continue
        end do  
        
        
    
    End Do
    
    
    !ftime = omp_get_wtime()
    !t1 = ftime-stime
            !< Saving position on text file
            !Do i=1,PartParam%nPart
    !Write(90,'(100F30.20)') t1

    !call DestroyHydro(HydroParam)
        
    !print *, GridDataConfig%isStructured
    
end subroutine
