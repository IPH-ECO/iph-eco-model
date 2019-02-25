!> This subroutine reads the simulation parameters. 
Subroutine AllocateMeteoVars(MeshParam,MeteoParam)
    Use MeshVars !, Only: nElem,nEdge
    Use Meteorological
    
    Implicit none
    type(MeshGridParam) :: MeshParam
    type(MeteorologicalParam) :: MeteoParam

    ! 1. Meteorological variables
    Allocate(MeteoParam%AtmPressure(MeshParam%nElem))
    Allocate(MeteoParam%AirTemp(MeshParam%nElem))
    Allocate(MeteoParam%SolarRad(MeshParam%nElem))
    Allocate(MeteoParam%RelHum(MeshParam%nElem))
    Allocate(MeteoParam%WindX(MeshParam%nEdge))
    Allocate(MeteoParam%WindY(MeshParam%nEdge))
    Allocate(MeteoParam%Precip(MeshParam%nElem))
    Allocate(MeteoParam%Evap(MeshParam%nElem))
    Allocate(MeteoParam%rhoair(MeshParam%nElem))
    
    
End Subroutine AllocateMeteoVars
  