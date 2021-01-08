!> This subroutine reads the hydrodynamic boundaru conditions. 
Subroutine GetHydroBoundaryConditions(HydroParam,MeshParam,dt,time,SimTime)
    
    Use Hydrodynamic
    Use MeshVars

    
    Implicit none
    
    Integer:: i,iElem,iEdge,iLayer
    Real:: t_interp(1)
    Real:: p_interp(1,1)
    Real:: NearZero = 1e-10
    Real:: Q_aux,SimTime,dt
    Integer:: Sig,l,Small
    Real:: time, H
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam

    t_interp = time
    !1. Reading water level 
    Do i =1,HydroParam%NWaterLevel
        Call interp_linear( 1, HydroParam%WaterLevelnTime(i), HydroParam%WaterLevelTime(i,1:HydroParam%WaterLevelnTime(i)), HydroParam%WaterLevelValue(i,1:HydroParam%WaterLevelnTime(i)), 1, t_interp, p_interp )
        HydroParam%WaterLevel(i) = p_interp(1,1)
    EndDo

    !HydroParam%WaterLevel(1) = 0.214 + 0.06*cos(2*HydroParam%pi*(Simtime-dt)/(355.0d0)) !CAYO
    !HydroParam%WaterLevel(1) = 2.14 + 0.6*cos(2*HydroParam%pi*(Simtime-dt)/(355.0d0)) !CAYO
    !2. Reading inflow/outflow
    Do i =1,HydroParam%NInflow
        iElem = HydroParam%IndexInflow(i,2)
        iEdge = HydroParam%IndexInflow(i,1)
        Small = Max(HydroParam%InFlowSmallm(i),HydroParam%SmallM(iEdge))
        Do iLayer = Small, HydroParam%InFlowCapitalM(i)
            If (HydroParam%InflownTime(i)==-999) Then
                H = HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge) !Surface Water Height
                If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
                    HydroParam%Fu(iLayer,iEdge) = (Sig(iElem,MeshParam%Right(iEdge),MeshParam%Left(iEdge)))*(Max(HydroParam%Pcri,H)**(2./3.))*(HydroParam%InflowValue(i,1)**(1./2.))/(Max(HydroParam%Pcri,H)**(1./6.)/(HydroParam%Rug(iElem)+NearZero))                            
                    HydroParam%Fu(iLayer,iEdge) = (Sig(iElem,MeshParam%Right(iEdge),MeshParam%Left(iEdge)))*sqrt(HydroParam%g*Max(HydroParam%Pcri,H)) !Critical Depth Cayo                                                            
                ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
                    HydroParam%Fu(iLayer,iEdge) = (Sig(iElem,MeshParam%Right(iEdge),MeshParam%Left(iEdge)))*(Max(HydroParam%Pcri,H)**(2./3.))*(HydroParam%InflowValue(i,1)**(1./2.))/(HydroParam%Rug(iElem)+NearZero) !Gradient/Zero-Depth Gradient Condition
                    HydroParam%Fu(iLayer,iEdge) = (Sig(iElem,MeshParam%Right(iEdge),MeshParam%Left(iEdge)))*sqrt(HydroParam%g*Max(HydroParam%Pcri,H)) !Critical Depth Cayo                                                                               
                ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
                    HydroParam%Fu(iLayer,iEdge) = (Sig(iElem,MeshParam%Right(iEdge),MeshParam%Left(iEdge)))*(Max(HydroParam%Pcri,H)**(2./3.))*(HydroParam%InflowValue(i,1)**(1./2.))/(Max(HydroParam%Pcri,H)**(1./6.)/((18.*log10(12.*Max(HydroParam%Pcri,H)/(HydroParam%Rug(iElem)/30.+NearZero)))+NearZero)) 
                EndIf
            Else
                Call interp_linear( 1, HydroParam%InflownTime(i), HydroParam%InflowTime(i,1:HydroParam%InflownTime(i)), HydroParam%InflowValue(i,1:HydroParam%InflownTime(i)), 1, t_interp, p_interp )
                HydroParam%Fu(iLayer,iEdge) = -(Sig(iElem,MeshParam%Right(iEdge),MeshParam%Left(iEdge)))*(p_interp(1,1))/(sum(HydroParam%DZj(HydroParam%InFlowSmallm(i):HydroParam%InFlowCapitalM(i),iEdge))*MeshParam%EdgeLength(iEdge))
            EndIf
        EndDo
    EndDo
    
    !!------------------------------------------------------------------------- 
    !! Vertedor
    !If ( iEdge == 8048.AND.l == 3777 ) Then
    !    If ( HydroParam%eta(l) > (698.00)  ) Then     ! If ( HydroParam%eta(l) > (700.00 + HydroParam%PCri)  ) Then
    !        Q_aux = 1.6 * 25.0 * ( HydroParam%eta(l) - 698.00 )**(3/2.)
    !        HydroParam%Fu(HydroParam%CapitalM(iEdge),iEdge) = Q_aux/( 25.0 * ( HydroParam%eta(l) - 698.00 ) )
    !    Else
    !        HydroParam%Fu(HydroParam%CapitalM(iEdge),iEdge) = 0.
    !    End If
    !Else
    !    Continue
    !End If        
    !!-------------------------------------------------------------------------
        

   Return 
End Subroutine GetHydroBoundaryConditions
    