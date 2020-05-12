!> This module declares the simulation parameters. 
Module SimulationModel
    use domain_types   
    
    implicit none
    
    type SimulationParam
        Integer:: iSim !< Type of simulation: iCod = 0 (Default); iCod = 1 (Sensitivity analysis) 
        Integer:: iSimHydro !< Include hydrodynamic model: iSimHydro = 0 (no); iSimHydro = 1 (yes)
        Integer:: iSimWQ !< Include water quality model: iSimWQ = 0 (no); iSimWQ = 1 (yes)
        Integer:: iSimSed !< Include sediment model: iSimSed = 0 (no); iSimSed = 1 (yes)
        Integer:: iRestart !< Restart a simulation: iRestart = 0 (no); iCod = 1 (yes) 
        Logical:: isStructuredMesh !< IIs it structured Mesh: isStructuredMesh = False (no); isStructuredMesh = True (yes)
        Character:: Day0*19	!< Initial time of simulation (AAAA/MM/DD HH:MM:SS)
        Integer:: IniTime !< Initial time of simulation (unix time)
        Integer:: FinalTime !< Final time of simulation (unix time)
        Real:: NumOfDays !< Number of simulated (days)
        Real:: dt !< Step time of the simulation (seg)
        Real:: dtday !< Step time of the simulation (day)
        Real:: time !< time of simulation (unix time)
        Integer:: Julday !< Julian day of simulation
        Integer:: Simday !< day of simulation
        Integer:: it !< time step unit
        Integer:: it_vtk !< time step unit for vtk files
        Integer:: i34=4 !< Element shape: i34 = 3 (triangle); i34 = 4 (square)
        Real:: NearZero = 1e-10 !< Small Number
        Integer:: OutputHydro(14) !< Hydrodynamic output flags
        Integer:: OutputWQ(40) !< WQ output flags
        Character:: OutputPath*500	!< Directory of model outputs
        Character:: OutputFile*200	!< Filename output
        Integer:: outputTimeStep !< Time step to save model outputs
        Integer:: recoveryTimeStep !< Time step to save recovery variables
        Integer:: RestartTime !< Time to restart a simulation
        real(c_double), dimension(:, :), pointer :: usave
        real(c_double), dimension(:, :), pointer :: wsave
        real(c_double), dimension(:), pointer :: etasave
        real(c_double), dimension(:), pointer :: start
        real(c_double), dimension(:), pointer :: finish
        real(c_double), dimension(:, :, :), pointer :: wqosave
        type(RecoveryVariables), pointer :: SaveVariables
        Real, Allocatable:: wqoutput(:)
        Integer:: nOutput !< Type of simulation: iCod = 0 (Default); iCod = 1 (Sensitivity analysis) 
        type(wqoOutputParameter), pointer :: wqoOutputParameters(:)
    contains
        procedure :: initializeSimulationParam
    end type
    
    contains
    
        subroutine initializeSimulationParam(this,sim)
            
            Integer i
            Integer:: idate(6)
            !character, pointer :: simulationLabel(:)
            !double precision, pointer :: layers(:)
            class(SimulationParam) :: this
            type(Simulation), intent(in) :: sim
            character(len=200):: text
            
            
            !call c_f_pointer(sim%label, simulationLabel, [sim%labelLength])
            !call c_f_pointer(sim%layers, layers, [sim%layersLength])
            call c_f_pointer(sim%recoveryVariables, this%SaveVariables)
    
            this%iSim = sim%simulationType                       ! Type of simulation: iCod = 1 (Default); iCod = 2 (Parameter Calibration); iCod = 3 (Sensitivity analysis) 
            this%iSimHydro = sim%hydrodynamic                    ! Include hydrodynamic model
            this%iSimWQ = sim%waterQuality                       ! Include water quality model
            this%iSimSed = 0                          ! Include sediment model
            this%iRestart = 0
            this%IniTime = sim%initialTime
            !Call unix2c(sim%initialTime, idate)
            !print *,idate
            !Day0*19                                        ! Initial time of simulation (AAAA/MM/DD HH:MM:SS)
            this%NumOfDays = sim%period !0.02435*800  !(1.1574074074074074074074074074074e-5)*1  ! Number of simulated (days)
            this%FinalTime = this%IniTime + int(this%NumOfDays*86400) !this%IniTime + 5!int(this%NumOfDays*86400)      ! Final time of simulation (AAAA/MM/DD HH:MM:SS)
            this%dt = sim%stepTime  !0.001                            ! Step time of the simulation (seg)
            this%dtday = this%dt/86400.                               ! Step time of the simulation (day)
            call c_f_pointer(sim%wqoOutputParameters, this%wqoOutputParameters, [sim%wqoOutputParametersLength])
             
            allocate(this%start(int(this%FinalTime/this%dt)+1))
            allocate(this%finish(int(this%FinalTime/this%dt)+1))
                
            If (.not.associated(this%SaveVariables)) then
                Allocate(this%SaveVariables)
                this%RestartTime = 0
            else
                allocate(this%usave(this%SaveVariables%layers, this%SaveVariables%edges))
                allocate(this%wsave(this%SaveVariables%layers + 1, this%SaveVariables%elements))
                allocate(this%etasave(this%SaveVariables%elements))
                call c_f_pointer(this%SaveVariables%u, this%usave, [this%SaveVariables%layers, this%SaveVariables%edges])
                call c_f_pointer(this%SaveVariables%w, this%wsave, [this%SaveVariables%layers+1, this%SaveVariables%elements])
                call c_f_pointer(this%SaveVariables%eta, this%etasave, [this%SaveVariables%elements])
                
                call c_f_pointer(this%SaveVariables%wqo, this%wqosave, [this%SaveVariables%layers,this%SaveVariables%elements, sim%wqoOutputParametersLength])
                
                this%RestartTime = this%SaveVariables%simulationTime
            Endif      
             
            If (this%RestartTime<this%IniTime) Then
                this%time = this%IniTime
                this%it = 0
                this%it_vtk = 0
            Else
                this%time = this%RestartTime + this%dt
                this%it = int((this%RestartTime-this%IniTime)/this%dt) + 1
                this%it_vtk = int((this%RestartTime-this%IniTime)/(this%outputTimeStep*this%dt)) + 1
            EndIf 
            
            !transferir para Module Mesh
        end subroutine 
        
        
    
end Module SimulationModel
    
    