Subroutine Pressure(HydroParam,MeshParam,dt)
    
    !$ use omp_lib
    Use MeshVars !, Only: nEdge,Left,Right,CirDistance
    Use Hydrodynamic
    
    Implicit None
    Integer:: iElem, iEdge, iLayer
    Integer:: l, r, iLayer_bar
    Real:: NearZero = 1e-10
    Real:: dt
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam

    !HydroParam%iBarot = 0
    HydroParam%PBarc = 0.d0
    ! 1. Define the Cell centered Horizontal Eddy-Viscosity:
    If (HydroParam%iBarot == 1) Then         ! User-defined Horizontal Diffusion
        !!$OMP parallel do default(none) shared(MeshParam,HydroParam) private(iLayer,iEdge,l,r,iLayer_bar,NearZero,dt)
        
        
        Do iEdge = 1, MeshParam%nEdge
            l = MeshParam%Left(iEdge)
            r = MeshParam%Right(iEdge)
            
            !if (iEdge==118.or.iEdge==124) Then
            !    Continue
            !Endif
            !
            
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                If ( r == 0 .or. HydroParam%H(iEdge)-HydroParam%hj(iEdge) <= HydroParam%PCRI+NearZero) Then
                    HydroParam%PBarc(iLayer,iEdge) = 0.d0
                Else
                    HydroParam%PBarc(iLayer,iEdge) = 0.d0
                    If (iLayer>=HydroParam%ElSmallm(l).and.iLayer>=HydroParam%ElSmallm(r)) Then 
                        Do iLayer_bar=iLayer,HydroParam%CapitalM(iEdge)
                            If (iLayer_bar==iLayer) Then 
                                HydroParam%PBarc(iLayer,iEdge) = HydroParam%PBarc(iLayer,iEdge) + HydroParam%g*dt/(MeshParam%CirDistance(iEdge)*HydroParam%rho0)*0.5d0*((HydroParam%Theta*HydroParam%sDRhoW(iLayer_bar,r)+(1.-HydroParam%Theta)*HydroParam%sDRhoWt(iLayer_bar,r))-(HydroParam%Theta*HydroParam%sDRhoW(iLayer_bar,l)+(1.-HydroParam%Theta)*HydroParam%sDRhoWt(iLayer_bar,l)))!(HydroParam%sDRhoW(iLayer_bar,r)-HydroParam%sDRhoW(iLayer_bar,l))*HydroParam%DZj(iLayer_bar,iEdge) !HydroParam%g*dt/(MeshParam%CirDistance(iEdge)*HydroParam%rho0)*0.5*(HydroParam%sDRhoW(iLayer_bar,r)-HydroParam%sDRhoW(iLayer_bar,l))*HydroParam%DZj(iLayer_bar,iEdge)
                            Else
                                HydroParam%PBarc(iLayer,iEdge) = HydroParam%PBarc(iLayer,iEdge) + HydroParam%g*dt/(MeshParam%CirDistance(iEdge)*HydroParam%rho0)*1.d0*((HydroParam%Theta*HydroParam%sDRhoW(iLayer_bar,r)+(1.-HydroParam%Theta)*HydroParam%sDRhoWt(iLayer_bar,r))-(HydroParam%Theta*HydroParam%sDRhoW(iLayer_bar,l)+(1.-HydroParam%Theta)*HydroParam%sDRhoWt(iLayer_bar,l)))!(HydroParam%sDRhoW(iLayer_bar,r)-HydroParam%sDRhoW(iLayer_bar,l))*HydroParam%DZj(iLayer_bar,iEdge) !HydroParam%g*dt/(MeshParam%CirDistance(iEdge)*HydroParam%rho0)*1.0*(HydroParam%sDRhoW(iLayer_bar,r)-HydroParam%sDRhoW(iLayer_bar,l))*HydroParam%DZj(iLayer_bar,iEdge)
                            EndIf
                        EndDo
                    EndIf
                EndIf
                
                !if (abs( HydroParam%PBarc(iLayer,iEdge))<1.e-8) then
                !     HydroParam%PBarc(iLayer,iEdge)=0.d0
                !endif
                
                HydroParam%Fu(iLayer,iEdge) = HydroParam%Fu(iLayer,iEdge) - HydroParam%PBarc(iLayer,iEdge)
            EndDo       ! Layer Loop
        EndDo       ! Edge Loop
        !!$OMP end parallel do
        
    EndIf
    Return
End Subroutine Pressure