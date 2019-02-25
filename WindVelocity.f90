Subroutine WindVelocity(HydroParam,MeshParam,MeteoParam)

    ! Get Wind Velocity in the Edges
    
    ! List of Modifications:
    !   15.06.2015: Routine Implementation       (Carlos Ruberto)
    ! Programmer: Carlos Ruberto    
    
    
    !$ use omp_lib
    Use Hydrodynamic
    Use MeshVars
    Use Meteorological
    !Use LimnologyVars
    Implicit None
    Integer:: iElem,iEdge,Face,l,r,iter
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam

    ! 1. Wind velocity
    !!$OMP parallel do default(none) shared(MeshParam,HydroParam,MeteoParam) private(iElem,iEdge,Face,l,r,iter)
    Do iElem = 1, MeshParam%nElem 
        HydroParam%Windix(iElem) = 0.
        HydroParam%Windiy(iElem) = 0.
        iter = 0.
        Do iEdge = 1,4
            Face = MeshParam%Edge(iEdge,iElem)
            l = MeshParam%Left(Face) 
            r = MeshParam%Right(Face)
            !If (r == 0) Then        ! Boundary Condition
            !    HydroParam%WindVel(:,Face) = 0.
            !    HydroParam%WindXY(:,Face) = 0.
            !Else
            If (MeteoParam%iReadWind == 0) Then !Constant Wind
                HydroParam%WindXY(:,Face) = (/ MeteoParam%WindX(Face), MeteoParam%WindY(Face) /)
            ElseIf (MeteoParam%iReadWind == 1) Then !Wind Time series (quando tiver mais de uma estação meteorologica fazer a interpolação em WindXY)
                HydroParam%WindXY(:,Face) = (/ MeteoParam%WindX(Face), MeteoParam%WindY(Face) /)
            EndIf
            HydroParam%WindVel(1,Face)= (-1.)*Dot_Product(HydroParam%WindXY(:,Face),MeshParam%NormalVector(:,iEdge,iElem))
            HydroParam%WindVel(2,Face)= (-1.)*Dot_Product(HydroParam%WindXY(:,Face),MeshParam%TangentVector(:,iEdge,iElem))
            iter = iter +1
            HydroParam%Windix(iElem) = HydroParam%Windix(iElem) + HydroParam%WindXY(1,Face)
            HydroParam%Windiy(iElem) = HydroParam%Windiy(iElem) + HydroParam%WindXY(2,Face)
            !EndIf
        EndDo
        ! Wind components in cell-centered
        If (iter == 0) Then ! Closed cell
            HydroParam%Windix(iElem) = 0.
            HydroParam%Windiy(iElem) = 0.
        Else
            HydroParam%Windix(iElem) = HydroParam%Windix(iElem)/iter
            HydroParam%Windiy(iElem) = HydroParam%Windiy(iElem)/iter
        EndIf
        ! 2. Air density
        MeteoParam%rhoair(iElem) = ((MeteoParam%AirTemp(iElem)+273)**(5.26)*4/(2.718281828459046*10000000000000))**(1./1.23419)
    EndDo
    !!$OMP end parallel do
    
    
    Return
End Subroutine WindVelocity