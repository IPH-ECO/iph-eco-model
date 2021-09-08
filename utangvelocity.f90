Subroutine utangvelocity(HydroParam,MeshParam,MeteoParam,dt)
    
    !$ use omp_lib
    Use MeshVars
    Use Hydrodynamic
    Use Meteorological
    
    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(MeteorologicalParam) :: MeteoParam
    Integer:: iEdge,iLayer, DIM, l, r
    Real:: aTh(MeshParam%KMax), bTh(MeshParam%KMax), cTh(MeshParam%KMax), fTh(MeshParam%KMax)
    Real:: dt,dzk,dzm,DZsjAcum,GammaB,GammaT,rhoair,chezy,rhoaircell,dummy
    Real :: w, H, Futn, Fvtn
    Real:: NearZero = 1e-10
    
    Do iEdge = 1, MeshParam%nEdge
        l = MeshParam%Left(iEdge) 
        r = MeshParam%Right(iEdge)
        
        ! 11.1.1 Get roughness 
        !H = HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge) !Surface Water Height
        !If (r == 0) Then
        !    If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
        !        Chezy = HydroParam%Rug(l)
        !    ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
        !        Chezy = Max(HydroParam%Pcri,HydroParam%H(iEdge))**(1./6.)/(HydroParam%Rug(l)+NearZero)
        !    ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
        !        Chezy = 18.*log10(12.*Max(HydroParam%Pcri,HydroParam%H(iEdge))/(HydroParam%Rug(l)/30.+NearZero))
        !    EndIf
        !    rhoairCell = MeteoParam%rhoair(l)
        !Else
        !    If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
        !        Chezy = 0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))
        !    ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
        !        Chezy = Max(HydroParam%Pcri,HydroParam%H(iEdge))**(1./6.)/(0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))+NearZero)
        !    ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
        !        Chezy = 18.*log10(12.*Max(HydroParam%Pcri,HydroParam%H(iEdge))/(0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))/30.+NearZero))
        !    EndIf
        !    rhoairCell = 0.5*(MeteoParam%rhoair(l) + MeteoParam%rhoair(r))
        !EndIf
        
        !If(r==0)Then
        !    r = l
        !EndIf
        !If (HydroParam%IndexWaterLevelEdge(iEdge)>0.and. H >HydroParam%PCRI+NearZero) Then
        !    Futn = HydroParam%Fu(HydroParam%Smallm(iEdge),iEdge) - (dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etaInfn(l) - HydroParam%etan(l)))
        !Else
        !    Futn = HydroParam%Fu(HydroParam%Smallm(iEdge),iEdge) - (dt/MeshParam%CirDistance(iEdge))*(HydroParam%g*(HydroParam%etan(r) - HydroParam%etan(l)))
        !EndIf
        !Fvtn = HydroParam%Fv(HydroParam%Smallm(iEdge),iEdge) - dt/MeshParam%CirDistance(iEdge)*HydroParam%g*(HydroParam%petan(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%petan(MeshParam%EdgeNodes(1,iEdge)))
        !
        !Call Tension(HydroParam%iWindStress,HydroParam%BottomTensionFlag,iEdge,HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),HydroParam%g,HydroParam%uxy(HydroParam%Smallm(iEdge),1,iEdge),HydroParam%uxy(HydroParam%Smallm(iEdge),2,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge),Chezy,rhoairCell,HydroParam%rho0,HydroParam%windDragConstant,HydroParam%WindVel(:,iEdge),HydroParam%WindXY(:,iEdge),HydroParam%GammaT,HydroParam%GammaB)
        Call Tension(HydroParam%iWindStress,HydroParam%BottomTensionFlag,iEdge,HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),HydroParam%g,HydroParam%uxy(HydroParam%Smallm(iEdge),1,iEdge),HydroParam%uxy(HydroParam%Smallm(iEdge),2,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge),Chezy,rhoairCell,HydroParam%rho0,HydroParam%windDragConstant,HydroParam%WindVel(:,iEdge),HydroParam%WindXY(:,iEdge),HydroParam%GammaT,dummy)
        !Call Tension(HydroParam%iWindStress,HydroParam%BottomTensionFlag,iEdge,HydroParam%CapitalM(iEdge),HydroParam%Smallm(iEdge),HydroParam%g,Futn,Fvtn,HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge),HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge),Chezy,rhoairCell,HydroParam%rho0,HydroParam%windDragConstant,HydroParam%WindVel(:,iEdge),HydroParam%WindXY(:,iEdge),HydroParam%GammaT,HydroParam%GammaB)

        ! If a face is dry, set the velocity to zero
        If (HydroParam%H(iEdge) - HydroParam%hj(iEdge) <= HydroParam%PCRI + NearZero) Then
            HydroParam%utang(:,iEdge)  = 0.d0
        Else !MeshParam%EdgeNodes(1,j)
            If (HydroParam%Smallm(iEdge) == HydroParam%CapitalM(iEdge)) Then
                !HydroParam%utang(HydroParam%Smallm(iEdge),iEdge)  = (HydroParam%DZhjt(HydroParam%Smallm(iEdge),iEdge)/(HydroParam%DZhjt(HydroParam%Smallm(iEdge),iEdge) + HydroParam%GammaB*dt + NearZero))*(HydroParam%Fv(HydroParam%Smallm(iEdge),iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%petan(MeshParam%EdgeNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge))))
                HydroParam%utang(HydroParam%Smallm(iEdge),iEdge)  = (HydroParam%DZhjt(HydroParam%Smallm(iEdge),iEdge)/(HydroParam%DZhjt(HydroParam%Smallm(iEdge),iEdge) + HydroParam%GammaB(iEdge)*dt + NearZero))*(HydroParam%Fv(HydroParam%Smallm(iEdge),iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%petan(MeshParam%EdgeNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge))))
                !HydroParam%utang(HydroParam%Smallm(iEdge),iEdge)  = (HydroParam%DZhjt(HydroParam%Smallm(iEdge),iEdge)/(HydroParam%DZhjt(HydroParam%Smallm(iEdge),iEdge) + HydroParam%GammaB*dt + NearZero))*(HydroParam%Fv(HydroParam%Smallm(iEdge),iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(HydroParam%utangNodes(2,iEdge)) - HydroParam%petan(HydroParam%utangNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(HydroParam%utangNodes(2,iEdge)) - HydroParam%peta(HydroParam%utangNodes(1,iEdge))))
            Else
                Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                    If (iLayer == HydroParam%Smallm(iEdge)) Then
                        dzk         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer+1,iEdge)) 
                        aTh(iLayer) = 0.
                        !bTh(iLayer) = MeshParam%CirDistance(iEdge)*(HydroParam%DZhjt(iLayer,iEdge) + dt*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk) + HydroParam%GammaB*dt)
                        bTh(iLayer) = MeshParam%CirDistance(iEdge)*(HydroParam%DZhjt(iLayer,iEdge) + dt*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk) + HydroParam%GammaB(iEdge)*dt)
                        cTh(iLayer) = -MeshParam%CirDistance(iEdge)*dt*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                        fTh(iLayer) = HydroParam%DZhjt(iLayer,iEdge)*MeshParam%CirDistance(iEdge)*(HydroParam%Fv(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%petan(MeshParam%EdgeNodes(1,iEdge))) + (HydroParam%pq(iLayer,MeshParam%EdgeNodes(2,iEdge))-HydroParam%pq(iLayer,MeshParam%EdgeNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge))))
                        !fTh(iLayer) = HydroParam%DZhjt(iLayer,iEdge)*MeshParam%CirDistance(iEdge)*(HydroParam%Fv(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(HydroParam%utangNodes(2,iEdge)) - HydroParam%petan(HydroParam%utangNodes(1,iEdge))) + (HydroParam%pq(iLayer,HydroParam%utangNodes(2,iEdge))-HydroParam%pq(iLayer,HydroParam%utangNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(HydroParam%utangNodes(2,iEdge)) - HydroParam%peta(HydroParam%utangNodes(1,iEdge))))
                    Elseif (iLayer == HydroParam%CapitalM(iEdge)) Then
                        dzk         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer-1,iEdge)) 
                        aTh(iLayer) = -MeshParam%CirDistance(iEdge)*dt/HydroParam%DZhjt(iLayer,iEdge)*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk)
                        bTh(iLayer) = MeshParam%CirDistance(iEdge)*(HydroParam%DZhjt(iLayer,iEdge) + dt*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzk))
                        cTh(iLayer) = 0.
                        fTh(iLayer) = HydroParam%DZhjt(iLayer,iEdge)*MeshParam%CirDistance(iEdge)*(HydroParam%Fv(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%petan(MeshParam%EdgeNodes(1,iEdge))) + (HydroParam%pq(iLayer,MeshParam%EdgeNodes(2,iEdge))-HydroParam%pq(iLayer,MeshParam%EdgeNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge))))
                        !fTh(iLayer) = HydroParam%DZhjt(iLayer,iEdge)*MeshParam%CirDistance(iEdge)*(HydroParam%Fv(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(HydroParam%utangNodes(2,iEdge)) - HydroParam%petan(HydroParam%utangNodes(1,iEdge))) + (HydroParam%pq(iLayer,HydroParam%utangNodes(2,iEdge))-HydroParam%pq(iLayer,HydroParam%utangNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(HydroParam%utangNodes(2,iEdge)) - HydroParam%peta(HydroParam%utangNodes(1,iEdge))))
                    Else
                        dzk         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer+1,iEdge)) 
                        dzm         = 0.5*(HydroParam%DZjt(iLayer,iEdge) + HydroParam%DZjt(iLayer-1,iEdge)) 
                        aTh(iLayer) = -MeshParam%CirDistance(iEdge)*dt*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm)
                        bTh(iLayer) = MeshParam%CirDistance(iEdge)*MeshParam%EdgeLength(iEdge)*(HydroParam%DZhjt(iLayer,iEdge) + dt*(HydroParam%VerEddyVisc(iLayer,iEdge)/dzm+HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk))
                        cTh(iLayer) = -MeshParam%CirDistance(iEdge)*dt*(HydroParam%VerEddyVisc(iLayer+1,iEdge)/dzk)
                        fTh(iLayer) = MeshParam%CirDistance(iEdge)*MeshParam%EdgeLength(iEdge)*(HydroParam%Fv(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%petan(MeshParam%EdgeNodes(1,iEdge))) + (HydroParam%pq(iLayer,MeshParam%EdgeNodes(2,iEdge))-HydroParam%pq(iLayer,MeshParam%EdgeNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge))))
                        !fTh(iLayer) = MeshParam%CirDistance(iEdge)*MeshParam%EdgeLength(iEdge)*(HydroParam%Fv(iLayer,iEdge) - (1.-HydroParam%Theta)*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%g*(HydroParam%petan(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%petan(MeshParam%EdgeNodes(1,iEdge))) + (HydroParam%pq(iLayer,HydroParam%utangNodes(2,iEdge))-HydroParam%pq(iLayer,HydroParam%utangNodes(1,iEdge)))) - HydroParam%Theta*HydroParam%g*(dt/MeshParam%EdgeLength(iEdge))*(HydroParam%peta(HydroParam%utangNodes(2,iEdge)) - HydroParam%peta(HydroParam%utangNodes(1,iEdge))))
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
    
    !xx.x Subsuperficial Tangencial Velocities:
    HydroParam%ustang(:,:) = 0.d0
    If (MeshParam%iBedrock == 1) Then
        Do iEdge = 1, MeshParam%nEdge
            DZsjAcum = 0.d0
            l = MeshParam%Left(iEdge) 
            r = MeshParam%Right(iEdge)
            If(r==0) Then
                r=l
            EndIf
			If(HydroParam%DZsj(HydroParam%Smallms(iEdge),iEdge) > 0) Then  ! if edge lower layer have sediment thickness
                Do iLayer = HydroParam%Smallms(iEdge),HydroParam%CapitalM(iEdge)
                    If(HydroParam%DZsj(iLayer,iEdge) > 0) Then 
                        DZsjAcum = DZsjAcum + HydroParam%DZsj(iLayer,iEdge) 
                        If(HydroParam%Smallms(iEdge) == HydroParam%CapitalMs(iEdge)) Then ! It has only one sediment layer
                            HydroParam%ustang(iLayer,iEdge)    =   -MeshParam%Kj(iLayer,iEdge)*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge))/MeshParam%CirDistance(iEdge))   
                            !HydroParam%ustang(iLayer,iEdge)    =   -MeshParam%Kj(iLayer,iEdge)*(HydroParam%peta(HydroParam%utangNodes(2,iEdge)) - HydroParam%peta(HydroParam%utangNodes(1,iEdge))/MeshParam%CirDistance(iEdge))   
                        !Check if eta is above the iLayer:
                        ElseIf(HydroParam%eta(r) - DZsjAcum > HydroParam%PCRI+NearZero .or. HydroParam%eta(l) - DZsjAcum >  HydroParam%PCRI/2.d0+NearZero) Then
                            HydroParam%ustang(iLayer,iEdge)    =   -MeshParam%Kj(iLayer,iEdge)*(HydroParam%peta(MeshParam%EdgeNodes(2,iEdge)) - HydroParam%peta(MeshParam%EdgeNodes(1,iEdge))/MeshParam%CirDistance(iEdge))
                            !HydroParam%ustang(iLayer,iEdge)    =   -MeshParam%Kj(iLayer,iEdge)*(HydroParam%peta(HydroParam%utangNodes(2,iEdge)) - HydroParam%peta(HydroParam%utangNodes(1,iEdge))/MeshParam%CirDistance(iEdge))
                        EndIf
                    EndIf
				EndDo   
            EndIf
	        ! xx.x Nullify Velocities Below the Bottom (u=0) 
            Do iLayer = 1, HydroParam%Smallms(iEdge) - 1.d0
                HydroParam%ustang(iLayer,iEdge) = 0.d0
            EndDo            
        EndDo   
        !xx.x Mean Velocity in each iEdge Layer (Weighted Superficial + Subsuperficial Velocities)
        Do iEdge = 1, MeshParam%nEdge
            HydroParam%umtang(:,iEdge)  =  HydroParam%utang(:,iEdge)
            HydroParam%umtang(:,iEdge)  = (HydroParam%utang(:,iEdge)*HydroParam%DZhjt(:,iEdge) + HydroParam%ustang(:,iEdge)*HydroParam%DZsjt(:,iEdge))/HydroParam%DZjt(:,iEdge)
        EndDo
    EndIf    
    
    
    Return    
End Subroutine utangvelocity