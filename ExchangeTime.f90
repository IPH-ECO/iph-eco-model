Subroutine ExchangeTime(HydroParam,MeshParam,MeteoParam,dt)

    ! This routine calculates the vertical water balance on the free-surface water by Evaporation and Precipitation
    ! Called in routine 0-MAIN
    !$ use omp_lib
    Use MeshVars 
    Use Hydrodynamic
    Use Meteorological
    
	Implicit none
    
    Integer:: iElem,iEdge,iLayer,Face,l,r
    Real:: NearZero = 1e-10
    Real:: dt,SumOutflow,VolCell,ExTime
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam

    !!$OMP parallel do default(none) shared(MeshParam,HydroParam,MeteoParam) private(iElem,Er,e_sat,B_Evap,e_evap,Ea,DELTA,NearZero,dt)
    Do iElem=1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            SumOutflow = 0.
            Do iEdge = 1,4
                Face = MeshParam%Edge(iEdge,iElem)
                l = MeshParam%Left(Face)
                r = MeshParam%Right(Face)
                If (HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)>=0) Then ! Outflow
                    SumOutflow = SumOutflow + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))
                EndIf
            EndDo
            VolCell = MeshParam%Area(iElem)*HydroParam%DZit(iLayer,iElem)
            ExTime = VolCell/(SumOutflow+NearZero)
        EndDo
    
    EndDo
    !!$OMP end parallel do

    Return
End Subroutine ExchangeTime
	
