Subroutine uvelocity(HydroParam,MeshParam,dt)
    
    !$ use omp_lib
    Use MeshVars
    Use Hydrodynamic
    
    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Integer:: iEdge,iLayer, DIM, l, r
    Double Precision:: aTh(MeshParam%KMax), bTh(MeshParam%KMax), cTh(MeshParam%KMax), fTh(MeshParam%KMax)
    Double Precision:: dt,dzk,dzm
    Double Precision:: NearZero = 1e-10
    
    Do iEdge = 1, MeshParam%nEdge
        l = MeshParam%Left(iEdge) 
        r = MeshParam%Right(iEdge)
        
        If (r == 0) Then
            If (HydroParam%IndexInflowEdge(iEdge) > 0) Then        ! Boundary Condition
                HydroParam%u(:,iEdge)  = HydroParam%Fu(:,iEdge)
            Else
                HydroParam%u(:,iEdge) = 0.
                If (HydroParam%IndexWaterLevel(l)>0.and.-HydroParam%hj(iEdge) + HydroParam%eta(l)>HydroParam%PCRI+NearZero) Then
                    Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                        HydroParam%u(iLayer,iEdge)  = HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%etaInf(l) - HydroParam%eta(l)) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%etaInf(l) - HydroParam%eta(l))
                    EndDo
                EndIf
            EndIf
        Else
            ! If a face is dry, set the velocity to zero
            If ( Max( HydroParam%PCRI,-HydroParam%hj(iEdge) + HydroParam%etan(l), -HydroParam%hj(iEdge) + HydroParam%etan(r) ) <= HydroParam%PCRI+NearZero.or.Max( HydroParam%PCRI,-HydroParam%hj(iEdge) + HydroParam%eta(l), -HydroParam%hj(iEdge) + HydroParam%eta(r) ) <= HydroParam%PCRI+NearZero) Then
                HydroParam%u(:,iEdge)  = 0.
            Else
                Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)      
                    If (iLayer == HydroParam%Smallm(iEdge)) Then
                        dzk         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer+1,iEdge)) 
                        aTh(iLayer) = 0.
                        bTh(iLayer) = 1. + dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                        cTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                        fTh(iLayer) = HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etan(r) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%eta(r) - HydroParam%eta(l))
                    Elseif (iLayer == HydroParam%CapitalM(iEdge)) Then
                        dzk         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer-1,iEdge)) 
                        aTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk)
                        bTh(iLayer) = 1. + dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk)
                        cTh(iLayer) = 0.
                        fTh(iLayer) = HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etan(r) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%eta(r) - HydroParam%eta(l))
                    Else
                        dzk         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer+1,iEdge)) 
                        dzm         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer-1,iEdge)) 
                        aTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm)
                        bTh(iLayer) = 1. + dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm+HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                        cTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                        fTh(iLayer) = HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etan(r) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%eta(r) - HydroParam%eta(l))
                    EndIf
                EndDo
                DIM = size(HydroParam%u(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),1)
                Call solve_tridiag(cTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),bTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),aTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),fTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),HydroParam%u(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),DIM)
        
            EndIf
        EndIf  
        
        ! 10.1 Copy Velocities Above the Free-Surface (du/dz=0) 
        Do iLayer = HydroParam%CapitalM(iEdge) + 1, MeshParam%KMAX
            HydroParam%u(iLayer,iEdge) = 0. !u(CapitalM(Face),Face) 
        EndDo
        ! 10.2 Nullify Velocities Below the Bottom (u=0) 
        Do iLayer = 1, HydroParam%Smallm(iEdge) - 1
            HydroParam%u(iLayer,iEdge) = 0.
        EndDo
        HydroParam%Hu(iEdge) = Sum( HydroParam%u(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge) )/(HydroParam%CapitalM(iEdge)-HydroParam%Smallm(iEdge)+1)
        
    EndDo
    
    Return    
End Subroutine uvelocity