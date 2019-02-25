Subroutine uvelocity(HydroParam,MeshParam,dt)
    
    !$ use omp_lib
    Use MeshVars
    Use Hydrodynamic
    
    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Integer:: iEdge,iLayer, DIM, l, r, i
    Real:: aTh(MeshParam%KMax), bTh(MeshParam%KMax), cTh(MeshParam%KMax), fTh(MeshParam%KMax)
    Real:: dt,dzk,dzm
    Real :: bp(MeshParam%KMax),vp(MeshParam%KMax),a_aux(MeshParam%KMax),b_aux(MeshParam%KMax),c_aux(MeshParam%KMax),v_aux(MeshParam%KMax),x_aux(MeshParam%KMax),d_aux(MeshParam%KMax),x(MeshParam%KMax)
    Real :: w
    Real:: NearZero = 1e-10
    
    Do iEdge = 1, MeshParam%nEdge
        l = MeshParam%Left(iEdge) 
        r = MeshParam%Right(iEdge)
        
        If (r == 0) Then
            If (HydroParam%IndexInflowEdge(iEdge) > 0) Then        ! Boundary Condition
                HydroParam%u(:,iEdge) = HydroParam%Fu(:,iEdge)
            Else
                HydroParam%u(:,iEdge) = 0.d0
                If (HydroParam%IndexWaterLevelEdge(iEdge)>0.and.-HydroParam%hj(iEdge) + HydroParam%eta(l)>HydroParam%PCRI+NearZero) Then
                    If (HydroParam%Smallm(iEdge)== HydroParam%CapitalM(iEdge)) Then
                        HydroParam%u(HydroParam%Smallm(iEdge),iEdge)  = HydroParam%Fu(HydroParam%Smallm(iEdge),iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etaInf(l) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%etaInf(l) - HydroParam%eta(l))
                    Else
                        Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                            If (iLayer == HydroParam%Smallm(iEdge)) Then
                                dzk         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer+1,iEdge)) 
                                aTh(iLayer) = 0. 
                                bTh(iLayer) = 1. + dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                                cTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                                fTh(iLayer) = HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etaInf(l) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%etaInf(l) - HydroParam%eta(l))
                            Elseif (iLayer == HydroParam%CapitalM(iEdge)) Then
                                dzk         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer-1,iEdge)) 
                                aTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk)
                                bTh(iLayer) = 1. + dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk)
                                cTh(iLayer) = 0.
                                fTh(iLayer) = HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etaInf(l) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%etaInf(l) - HydroParam%eta(l))
                            Else
                                dzk         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer+1,iEdge)) 
                                dzm         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer-1,iEdge)) 
                                aTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm)
                                bTh(iLayer) = 1. + dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm+HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                                cTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                                fTh(iLayer) = HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etaInf(l) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%etaInf(l) - HydroParam%eta(l))
                            EndIf
                        EndDo
                        !n = size(HydroParam%u(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),1)
                        !Call solve_tridiag(cTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),bTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),aTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),fTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),HydroParam%u(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),DIM)
                    
                        !Tridiagonal Matrix solver
                        !Flipping matrix
                        do i = HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),-1
                         a_aux(HydroParam%CapitalM(iEdge)-i+HydroParam%Smallm(iEdge)) = cTh(i)
                         b_aux(HydroParam%CapitalM(iEdge)-i+HydroParam%Smallm(iEdge)) = bTh(i)
                         c_aux(HydroParam%CapitalM(iEdge)-i+HydroParam%Smallm(iEdge)) = aTh(i)
                         v_aux(HydroParam%CapitalM(iEdge)-i+HydroParam%Smallm(iEdge)) = fTh(i)
                        enddo
 
                        d_aux = 0   
                        x = d_aux
                        w = b_aux(HydroParam%Smallm(iEdge))
                        x(HydroParam%Smallm(iEdge)) = v_aux(HydroParam%Smallm(iEdge))/w
                        do i=HydroParam%Smallm(iEdge)+1,HydroParam%CapitalM(iEdge)
                            d_aux(i-1) = c_aux(i-1)/w;
                            w = b_aux(i) - a_aux(i)*d_aux(i-1);
                            x(i) = ( v_aux(i) - a_aux(i)*x(i-1) )/w;
                        enddo
                        do i=HydroParam%CapitalM(iEdge)-1, HydroParam%Smallm(iEdge), -1
                           x(i) = x(i) - d_aux(i)*x(i+1);
                        enddo

                        !Flipping matrix
                        do i = HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),-1
                         x_aux(HydroParam%CapitalM(iEdge)-i+HydroParam%Smallm(iEdge)) = x(i)
                        enddo
                        HydroParam%u(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge) = x_aux(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge))
                    EndIf
                EndIf
            EndIf
        Else
            ! If a face is dry, set the velocity to zero
            If ( Max( HydroParam%PCRI,-HydroParam%hj(iEdge) + HydroParam%etan(l), -HydroParam%hj(iEdge) + HydroParam%etan(r) ) <= HydroParam%PCRI+NearZero.or.Max( HydroParam%PCRI,-HydroParam%hj(iEdge) + HydroParam%eta(l), -HydroParam%hj(iEdge) + HydroParam%eta(r) ) <= HydroParam%PCRI+NearZero) Then
                HydroParam%u(:,iEdge)  = 0.d0
            Else
                If (HydroParam%Smallm(iEdge)== HydroParam%CapitalM(iEdge)) Then
                    HydroParam%u(HydroParam%Smallm(iEdge),iEdge)  = HydroParam%Fu(HydroParam%Smallm(iEdge),iEdge) - (1.d0-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etan(r) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%eta(r) - HydroParam%eta(l))
                Else
                    Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                        If (iLayer == HydroParam%Smallm(iEdge)) Then
                            dzk         = 0.5d0*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer+1,iEdge)) 
                            aTh(iLayer) = 0.d0
                            bTh(iLayer) = 1.d0 + dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                            cTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                            fTh(iLayer) = HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etan(r) - HydroParam%etan(l)) + ( HydroParam%q(iLayer,r)-HydroParam%q(iLayer,l) )) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%eta(r) - HydroParam%eta(l))
                        Elseif (iLayer == HydroParam%CapitalM(iEdge)) Then
                            dzk         = 0.5d0*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer-1,iEdge)) 
                            aTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk)
                            bTh(iLayer) = 1.d0 + dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk)
                            cTh(iLayer) = 0.d0
                            fTh(iLayer) = HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etan(r) - HydroParam%etan(l)) + ( HydroParam%q(iLayer,r)-HydroParam%q(iLayer,l) )) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%eta(r) - HydroParam%eta(l))
                        Else
                            dzk         = 0.5d0*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer+1,iEdge)) 
                            dzm         = 0.5d0*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer-1,iEdge)) 
                            aTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm)
                            bTh(iLayer) = 1.d0 + dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm+HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                            cTh(iLayer) = -dt/HydroParam%DZjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                            fTh(iLayer) = HydroParam%Fu(iLayer,iEdge) - (1.d0-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etan(r) - HydroParam%etan(l)) + ( HydroParam%q(iLayer,r)-HydroParam%q(iLayer,l) )) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%eta(r) - HydroParam%eta(l))
                        EndIf
                    EndDo
                    !n = size(HydroParam%u(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),1)
                    !Call solve_tridiag(cTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),bTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),aTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),fTh(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge)),HydroParam%u(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge),DIM)
                    
                    !Tridiagonal Matrix solver
                    !Flipping matrix
                    do i = HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),-1
                     a_aux(HydroParam%CapitalM(iEdge)-i+HydroParam%Smallm(iEdge)) = cTh(i)
                     b_aux(HydroParam%CapitalM(iEdge)-i+HydroParam%Smallm(iEdge)) = bTh(i)
                     c_aux(HydroParam%CapitalM(iEdge)-i+HydroParam%Smallm(iEdge)) = aTh(i)
                     v_aux(HydroParam%CapitalM(iEdge)-i+HydroParam%Smallm(iEdge)) = fTh(i)
                    enddo
 
                    d_aux = 0.d0   
                    x = d_aux
                    w = b_aux(HydroParam%Smallm(iEdge))
                    x(HydroParam%Smallm(iEdge)) = v_aux(HydroParam%Smallm(iEdge))/w
                    do i=HydroParam%Smallm(iEdge)+1.d0,HydroParam%CapitalM(iEdge)
                        d_aux(i-1.d0) = c_aux(i-1.d0)/w;
                        w = b_aux(i) - a_aux(i)*d_aux(i-1.d0);
                        x(i) = ( v_aux(i) - a_aux(i)*x(i-1.d0) )/w;
                    enddo
                    do i=HydroParam%CapitalM(iEdge)-1.d0, HydroParam%Smallm(iEdge), -1.d0
                       x(i) = x(i) - d_aux(i)*x(i+1.d0);
                    enddo

                    !Flipping matrix
                    do i = HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),-1.d0
                     x_aux(HydroParam%CapitalM(iEdge)-i+HydroParam%Smallm(iEdge)) = x(i)
                    enddo
                    HydroParam%u(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge) = x_aux(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge))
                    
                EndIf
            EndIf
        EndIf  
        
        ! 10.1 Copy Velocities Above the Free-Surface (du/dz=0) 
        Do iLayer = HydroParam%CapitalM(iEdge) + 1.d0, MeshParam%KMAX
            HydroParam%u(iLayer,iEdge) = 0.d0 !u(CapitalM(Face),Face) 
        EndDo
        ! 10.2 Nullify Velocities Below the Bottom (u=0) 
        Do iLayer = 1, HydroParam%Smallm(iEdge) - 1.d0
            HydroParam%u(iLayer,iEdge) = 0.d0
        EndDo
        HydroParam%Hu(iEdge) = Sum( HydroParam%u(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge) )/(HydroParam%CapitalM(iEdge)-HydroParam%Smallm(iEdge)+1.d0)
        
    EndDo
    
    Return    
End Subroutine uvelocity