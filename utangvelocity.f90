Subroutine utangvelocity(HydroParam,MeshParam,dt)
    
    !$ use omp_lib
    Use MeshVars
    Use Hydrodynamic
    
    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Integer:: iEdge,iLayer, DIM, l, r
    Real:: aTh(MeshParam%KMax), bTh(MeshParam%KMax), cTh(MeshParam%KMax), fTh(MeshParam%KMax)
    Real:: dt,dzk,dzm
    Real:: NearZero = 1e-10
    
    Do iEdge = 1, MeshParam%nEdge
        l = MeshParam%Left(iEdge) 
        r = MeshParam%Right(iEdge)
        
        ! If a face is dry, set the velocity to zero
        
        If (HydroParam%H(iEdge) <= HydroParam%PCRI + NearZero) Then
            HydroParam%utang(:,iEdge)  = 0.
        Else !MeshParam%EdgeNodes(1,j)
            If (HydroParam%Smallm(iEdge)== HydroParam%CapitalM(iEdge)) Then
                HydroParam%utang(HydroParam%Smallm(iEdge),iEdge)  = HydroParam%Fv(HydroParam%Smallm(iEdge),iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%petan(MeshParam%EdgeNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge)))
            Else
                Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                    If (iLayer == HydroParam%Smallm(iEdge)) Then
                        dzk         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer+1,iEdge)) 
                        aTh(iLayer) = 0.
                        bTh(iLayer) = 1. + dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                        cTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                        fTh(iLayer) = HydroParam%Fv(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%petan(MeshParam%EdgeNodes(1,iEdge))) + (HydroParam%pq(iLayer,MeshParam%EdgeNodes(2,iEdge))-HydroParam%pq(iLayer,MeshParam%EdgeNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge)))
                    Elseif (iLayer == HydroParam%CapitalM(iEdge)) Then
                        dzk         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer-1,iEdge)) 
                        aTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk)
                        bTh(iLayer) = 1. + dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk)
                        cTh(iLayer) = 0.
                        fTh(iLayer) = HydroParam%Fv(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%petan(MeshParam%EdgeNodes(1,iEdge))) + (HydroParam%pq(iLayer,MeshParam%EdgeNodes(2,iEdge))-HydroParam%pq(iLayer,MeshParam%EdgeNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge)))
                    Else
                        dzk         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer+1,iEdge)) 
                        dzm         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer-1,iEdge)) 
                        aTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm)
                        bTh(iLayer) = 1. + dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm+HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                        cTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                        fTh(iLayer) = HydroParam%Fv(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%petan(MeshParam%EdgeNodes(1,iEdge))) + (HydroParam%pq(iLayer,MeshParam%EdgeNodes(2,iEdge))-HydroParam%pq(iLayer,MeshParam%EdgeNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge)))
                    EndIf
                EndDo
                DIM = size(HydroParam%utang(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),1)
                Call solve_tridiag(cTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),bTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),aTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),fTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),HydroParam%utang(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),DIM)
            EndIf
        EndIf
        
        ! 10.1 Copy Velocities Above the Free-Surface (du/dz=0) 
        Do iLayer = HydroParam%CapitalM(iEdge) + 1, MeshParam%KMAX
            HydroParam%utang(iLayer,iEdge) = 0. !u(CapitalM(Face),Face) 
        EndDo
        ! 10.2 Nullify Velocities Below the Bottom (u=0) 
        Do iLayer = 1, HydroParam%Smallm(iEdge) - 1
            HydroParam%utang(iLayer,iEdge) = 0.
        EndDo
        
    EndDo
    
    Return    
End Subroutine utangvelocity