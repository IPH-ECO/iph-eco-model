Subroutine uvelocity(HydroParam,MeshParam,MeteoParam,dt)
    ! List of Modifications: 
    !   -> 20.08.2019: Routine Update               (Cayo Lopes)
    
    
    !$ use omp_lib
    Use MeshVars
    Use Hydrodynamic
    Use Meteorological
    
    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam
    Integer:: iEdge,iLayer, DIM, l, r, i
    Real:: aTh(MeshParam%KMax), bTh(MeshParam%KMax), cTh(MeshParam%KMax), fTh(MeshParam%KMax)
    Real:: dt,dzk,dzm,DZsjAcum,GammaB,GammaT,rhoair,chezy,rhoaircell
    Real :: bp(MeshParam%KMax),vp(MeshParam%KMax),a_aux(MeshParam%KMax),b_aux(MeshParam%KMax),c_aux(MeshParam%KMax),v_aux(MeshParam%KMax),x_aux(MeshParam%KMax),d_aux(MeshParam%KMax),x(MeshParam%KMax)
    Real :: w
    Real:: NearZero = 1e-10
    
    ! 12.1 Surface Velocity    
    Do iEdge = 1, MeshParam%nEdge
        l = MeshParam%Left(iEdge) 
        r = MeshParam%Right(iEdge)
        
        
        !!bench 02:
        !if ( 100-NearZero <= MeshParam%EdgeBary(1,iEdge)<= 100 + NearZero) then
        !    if (MeshParam%EdgeBary(2,iEdge) >= 110 - NearZero) then
        !        HydroParam%u(:,iEdge)  = 0.
        !        HydroParam%Fu(:,iEdge) = 0.
        !        r = 0
        !    elseif(MeshParam%EdgeBary(2,iEdge) <= 90 + NearZero) then
        !        HydroParam%u(:,iEdge)  = 0.
        !        HydroParam%Fu(:,iEdge) = 0.   
        !        r = 0
        !    endif
        !endif        

        ! 6.1 Get roughness 
        If (r == 0) Then
            If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
                Chezy = HydroParam%Rug(l)
            ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
                Chezy = Max(HydroParam%Pcri,HydroParam%H(iEdge))**(1./6.)/(HydroParam%Rug(l)+NearZero)
            ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
                Chezy = 18.*log10(12.*Max(HydroParam%Pcri,HydroParam%H(iEdge))/(HydroParam%Rug(l)/30.+NearZero))
            EndIf
            rhoairCell = MeteoParam%rhoair(l)
        Else
            If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
                Chezy = 0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))
            ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
                Chezy = Max(HydroParam%Pcri,HydroParam%H(iEdge))**(1./6.)/(0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))+NearZero)
            ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
                Chezy = 18.*log10(12.*Max(HydroParam%Pcri,HydroParam%H(iEdge))/(0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))/30.+NearZero))
            EndIf
            rhoairCell = 0.5*(MeteoParam%rhoair(l) + MeteoParam%rhoair(r))
        EndIf
        
        Call Tension(HydroParam%iWindStress,HydroParam%BottomTensionFlag,iEdge,HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),HydroParam%g,HydroParam%uxy(HydroParam%Smallm(iEdge),1,iEdge),HydroParam%uxy(HydroParam%Smallm(iEdge),2,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge),Chezy,rhoairCell,HydroParam%rho0,HydroParam%windDragConstant,HydroParam%WindVel(:,iEdge),HydroParam%WindXY(:,iEdge),HydroParam%GammaT,HydroParam%GammaB)

        ! If not has neighbour cell:
        If (r == 0) Then
            ! When has flow boundary condition:
            If (HydroParam%IndexInflowEdge(iEdge) > 0) Then        ! Boundary Condition
                HydroParam%u(:,iEdge) = HydroParam%Fu(:,iEdge)
            Else
                HydroParam%u(:,iEdge) = 0.d0
                ! In which case has Water Level Bound condition:
                If (HydroParam%IndexWaterLevelEdge(iEdge)>0.and.-HydroParam%hj(iEdge) + HydroParam%eta(l)>HydroParam%PCRI/2+NearZero) Then
                    ! 2D Case:
                    If (HydroParam%Smallm(iEdge) == HydroParam%CapitalM(iEdge)) Then
                        HydroParam%u(HydroParam%Smallm(iEdge),iEdge)  = (HydroParam%DZhjt(iLayer,iEdge)/(HydroParam%DZhjt(iLayer,iEdge) + HydroParam%GammaB*dt + NearZero))*(HydroParam%Fu(HydroParam%Smallm(iEdge),iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etaInf(l) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%etaInf(l) - HydroParam%eta(l)))
                    Else
                        !3D Case:
                        Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                            If (iLayer == HydroParam%Smallm(iEdge)) Then
                                
                                dzk         = 0.5*(HydroParam%DZhjt(iLayer,iEdge) + HydroParam%DZhjt(iLayer+1,iEdge)) 
                                aTh(iLayer) = 0. 
                                bTh(iLayer) = MeshParam%EdgeLength(iEdge)*(HydroParam%DZhjt(iLayer,iEdge) + dt*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk) + HydroParam%GammaB*dt)
                                cTh(iLayer) = -MeshParam%EdgeLength(iEdge)*dt*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                                fTh(iLayer) = HydroParam%DZhjt(iLayer,iEdge)*MeshParam%EdgeLength(iEdge)*(HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etaInf(l) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%etaInf(l) - HydroParam%eta(l)))
                            Elseif (iLayer == HydroParam%CapitalM(iEdge)) Then
                                dzk         = 0.5*(HydroParam%DZhjt(iLayer,iEdge) + HydroParam%DZhjt(iLayer-1,iEdge)) 
                                aTh(iLayer) = -MeshParam%EdgeLength(iEdge)*dt/HydroParam%DZhjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk)
                                bTh(iLayer) = MeshParam%EdgeLength(iEdge)*(HydroParam%DZhjt(iLayer,iEdge) + dt*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk))
                                cTh(iLayer) = 0.
                                fTh(iLayer) = HydroParam%DZhjt(iLayer,iEdge)*MeshParam%EdgeLength(iEdge)*(HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etaInf(l) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%etaInf(l) - HydroParam%eta(l)))
                            Else
                                dzk         = 0.5*(HydroParam%DZhjt(iLayer,iEdge) + HydroParam%DZhjt(iLayer+1,iEdge)) 
                                dzm         = 0.5*(HydroParam%DZhjt(iLayer,iEdge) + HydroParam%DZhjt(iLayer-1,iEdge)) 
                                aTh(iLayer) = -MeshParam%EdgeLength(iEdge)*dt*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm)
                                bTh(iLayer) = MeshParam%EdgeLength(iEdge)*(HydroParam%DZhjt(iLayer,iEdge) + dt*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm+HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk))
                                cTh(iLayer) = -MeshParam%EdgeLength(iEdge)*dt*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                                fTh(iLayer) = HydroParam%DZhjt(iLayer,iEdge)*MeshParam%EdgeLength(iEdge)*(HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etaInf(l) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%etaInf(l) - HydroParam%eta(l)))
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
            !Cell without bound conditions:
            ! If a face is dry, set the velocity to zero !CAYO must be check for subsurface
            If ( Max( HydroParam%PCRI,-HydroParam%hj(iEdge) + HydroParam%etan(l), -HydroParam%hj(iEdge) + HydroParam%etan(r) ) <= HydroParam%PCRI/2+NearZero.or.Max( HydroParam%PCRI,-HydroParam%hj(iEdge) + HydroParam%eta(l), -HydroParam%hj(iEdge) + HydroParam%eta(r) ) <= HydroParam%PCRI/2+NearZero) Then
                HydroParam%u(:,iEdge)  = 0.d0
            Else
                !2D Case:
                If (HydroParam%Smallm(iEdge) == HydroParam%CapitalM(iEdge)) Then
                    HydroParam%u(HydroParam%Smallm(iEdge),iEdge)  = (HydroParam%DZhjt(iLayer,iEdge)/(HydroParam%DZhjt(iLayer,iEdge) + HydroParam%GammaB*dt + NearZero))*(HydroParam%Fu(HydroParam%Smallm(iEdge),iEdge) - (1.d0-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etan(r) - HydroParam%etan(l))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%eta(r) - HydroParam%eta(l)))
                !3D Case:
                Else
                    
                    !DZjt must be DZjht for subsurface case
                    Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                        If (iLayer == HydroParam%Smallm(iEdge)) Then
                            dzk         = 0.5d0*(HydroParam%DZhjt(iLayer,iEdge) + HydroParam%DZhjt(iLayer+1,iEdge)) 
                            aTh(iLayer) = 0.d0
                            bTh(iLayer) = MeshParam%EdgeLength(iEdge)*(HydroParam%DZhjt(iLayer,iEdge) + dt*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk) + HydroParam%GammaB*dt)
                            cTh(iLayer) = -MeshParam%EdgeLength(iEdge)*dt*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                            fTh(iLayer) = HydroParam%DZhjt(iLayer,iEdge)*MeshParam%EdgeLength(iEdge)*(HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etan(r) - HydroParam%etan(l)) + ( HydroParam%q(iLayer,r)-HydroParam%q(iLayer,l) )) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%eta(r) - HydroParam%eta(l)))
                        Elseif (iLayer == HydroParam%CapitalM(iEdge)) Then
                            dzk         = 0.5d0*(HydroParam%DZhjt(iLayer,iEdge) + HydroParam%DZhjt(iLayer-1,iEdge)) 
                            aTh(iLayer) = -MeshParam%EdgeLength(iEdge)*dt*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk)
                            bTh(iLayer) = MeshParam%EdgeLength(iEdge)*(HydroParam%DZhjt(iLayer,iEdge) + dt*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk))
                            cTh(iLayer) = 0.d0
                            fTh(iLayer) = HydroParam%DZhjt(iLayer,iEdge)*MeshParam%EdgeLength(iEdge)*(HydroParam%Fu(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etan(r) - HydroParam%etan(l)) + ( HydroParam%q(iLayer,r)-HydroParam%q(iLayer,l) )) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%eta(r) - HydroParam%eta(l)))
                        Else
                            dzk         = 0.5d0*(HydroParam%DZhjt(iLayer,iEdge) + HydroParam%DZhjt(iLayer+1,iEdge)) 
                            dzm         = 0.5d0*(HydroParam%DZhjt(iLayer,iEdge) + HydroParam%DZhjt(iLayer-1,iEdge)) 
                            aTh(iLayer) = - MeshParam%EdgeLength(iEdge)*dt*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm)
                            bTh(iLayer) =  MeshParam%EdgeLength(iEdge)*(HydroParam%DZhjt(iLayer,iEdge) + dt*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm+HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk))
                            cTh(iLayer) = - MeshParam%EdgeLength(iEdge)*dt*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                            fTh(iLayer) = HydroParam%DZhjt(iLayer,iEdge)*MeshParam%EdgeLength(iEdge)*(HydroParam%Fu(iLayer,iEdge) - (1.d0-HydroParam%Theta)*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etan(r) - HydroParam%etan(l)) + ( HydroParam%q(iLayer,r)-HydroParam%q(iLayer,l) )) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%CirDistance(iEdge))*(HydroParam%eta(r) - HydroParam%eta(l)))
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
                
        ! 12.1 Copy Velocities Above the Free-Surface (du/dz=0) 
        Do iLayer = HydroParam%CapitalM(iEdge) + 1.d0, MeshParam%KMAX
            HydroParam%u(iLayer,iEdge) = 0.d0 !u(CapitalM(Face),Face) 
        EndDo
        ! 12.2 Nullify Velocities Below the Bottom (u=0) 
        Do iLayer = 1, HydroParam%Smallms(iEdge) - 1.d0
            HydroParam%u(iLayer,iEdge) = 0.d0
        EndDo
        HydroParam%Hu(iEdge) = Sum(HydroParam%u(HydroParam%Smallm(iEdge):HydroParam%CapitalM(iEdge),iEdge))/(HydroParam%CapitalM(iEdge)-HydroParam%Smallm(iEdge)+1.d0)
        HydroParam%um(:,iEdge) =  HydroParam%u(:,iEdge)!CAYO
    EndDo
    
            
    !12.3 Subsuperficial Velocities
    HydroParam%us(:,:) = 0.d0
    If (MeshParam%iBedrock == 1) Then
        Do iEdge = 1, MeshParam%nEdge
            DZsjAcum = 0.d0
            l = MeshParam%Left(iEdge) 
            r = MeshParam%Right(iEdge)
            !HydroParam%us(:,iEdge) = 0.d0
			If(HydroParam%DZsj(HydroParam%Smallms(iEdge),iEdge) > 0) Then ! if edge lower layer have sediment thickness
                If (r == 0) Then
                    If (HydroParam%IndexWaterLevelEdge(iEdge)>0.and.-HydroParam%sj(iEdge) + HydroParam%eta(l)>HydroParam%PCRI/2.d0+NearZero) Then
                        If ((HydroParam%Smallms(iEdge) == HydroParam%CapitalMs(iEdge).and.HydroParam%Smallms(iEdge) == HydroParam%CapitalM(iEdge) )) Then
                                HydroParam%us(HydroParam%Smallms(iEdge),iEdge)  =   -MeshParam%Kj(HydroParam%Smallms(iEdge),iEdge)*(HydroParam%etaInf(l) - HydroParam%etan(l))/MeshParam%CirDistance(iEdge)
                                !HydroParam%us(HydroParam%Smallms(iEdge),iEdge)  =   (1-HydroParam%Theta)*HydroParam%us(HydroParam%Smallms(iEdge),iEdge)-HydroParam%Theta*MeshParam%Kj(HydroParam%Smallms(iEdge),iEdge)*(HydroParam%etaInf(l) - HydroParam%eta(l))/MeshParam%CirDistance(iEdge)                                           
                                !HydroParam%u(HydroParam%Smallm(iEdge),iEdge)  =   -MeshParam%Kj(HydroParam%Smallm(iEdge),iEdge)*(HydroParam%etaInf(l) - HydroParam%eta(l))/MeshParam%CirDistance(iEdge)
                        Else
                            Do iLayer = HydroParam%Smallms(iEdge),HydroParam%CapitalMs(iEdge)
                                If(HydroParam%DZsj(iLayer,iEdge) > 0) Then  
                                    HydroParam%us(iLayer,iEdge)    =   -MeshParam%Kj(iLayer,iEdge)*(HydroParam%etaInf(l) - HydroParam%etan(l))/MeshParam%CirDistance(iEdge)                                   
                                    !HydroParam%us(iLayer,iEdge)    =  (1-HydroParam%Theta)*HydroParam%us(iLayer,iEdge) - MeshParam%Kj(iLayer,iEdge)*HydroParam%Theta*(HydroParam%etaInf(l) - HydroParam%eta(l))/MeshParam%CirDistance(iEdge)
                                    !HydroParam%u(iLayer,iEdge)    =   -MeshParam%Kj(iLayer,iEdge)*(HydroParam%etaInf(l) - HydroParam%eta(l))/MeshParam%CirDistance(iEdge)
                                EndIf
                            EndDo
                        EndIf
                    EndIf
                Else
                    Do iLayer = HydroParam%Smallms(iEdge),HydroParam%CapitalM(iEdge)
                        
                        If(HydroParam%DZsj(iLayer,iEdge) > 0) Then
                            DZsjAcum = DZsjAcum + HydroParam%DZsj(iLayer,iEdge)
                            If(HydroParam%Smallms(iEdge) == HydroParam%CapitalMs(iEdge)) Then
                                HydroParam%us(iLayer,iEdge)    =   -MeshParam%Kj(iLayer,iEdge)*(HydroParam%etan(r) - HydroParam%etan(l))/MeshParam%CirDistance(iEdge)            
                            ElseIf(HydroParam%eta(r) - DZsjAcum > HydroParam%PCRI+NearZero .or. HydroParam%eta(l) - DZsjAcum >  HydroParam%PCRI/2.d0+NearZero) Then
                                HydroParam%us(iLayer,iEdge)    =   -MeshParam%Kj(iLayer,iEdge)*(HydroParam%etan(r) - HydroParam%etan(l))/MeshParam%CirDistance(iEdge)            
                                !HydroParam%us(iLayer,iEdge)    =  (1-HydroParam%Theta)*HydroParam%us(iLayer,iEdge) - MeshParam%Kj(iLayer,iEdge)*HydroParam%Theta*(HydroParam%etan(r) - HydroParam%etan(l))/MeshParam%CirDistance(iEdge)                         
                                !HydroParam%u(iLayer,iEdge)    =   -MeshParam%Kj(iLayer,iEdge)*(HydroParam%etan(r) - HydroParam%etan(l))/MeshParam%CirDistance(iEdge)
                            EndIf
                        EndIf

				    EndDo   
                EndIf
			EndIf
	        ! 10.4 Nullify Velocities Below the Bottom (u=0) 
            Do iLayer = 1, HydroParam%Smallms(iEdge) - 1.d0
                HydroParam%us(iLayer,iEdge) = 0.d0
            EndDo            
        EndDo   
        !Superficial + Subsuperficial Velocities !CAYO
        Do iEdge = 1, MeshParam%nEdge
            HydroParam%um(:,iEdge) =  HydroParam%u(:,iEdge)
            HydroParam%um(:,iEdge)   = (HydroParam%u(:,iEdge)*HydroParam%DZhjt(:,iEdge) + HydroParam%us(:,iEdge)*HydroParam%DZsjt(:,iEdge))/HydroParam%DZjt(:,iEdge)
        EndDo

    EndIf
    

    
    Return    
End Subroutine uvelocity