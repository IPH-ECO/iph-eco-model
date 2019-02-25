!> This subroutine reads the hydrodynamic boundaru conditions. 
Subroutine GetHydroBoundaryConditions(HydroParam,MeshParam,time)
    
    Use Hydrodynamic
    Use MeshVars
    
    Implicit none
    
    Integer:: i,iElem,iEdge,iLayer
    Real:: t_interp(1)
    Real:: p_interp(1,1)
    Real:: NearZero = 1e-10
    Real:: Q_aux
    Integer:: Sig,l
    Real:: time
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam

    t_interp = time
    !1. Reading water level 
    Do i =1,HydroParam%NWaterLevel
        Call interp_linear( 1, HydroParam%WaterLevelnTime(i), HydroParam%WaterLevelTime(i,1:HydroParam%WaterLevelnTime(i)), HydroParam%WaterLevelValue(i,1:HydroParam%WaterLevelnTime(i)), 1, t_interp, p_interp )
        HydroParam%WaterLevel(i) = p_interp(1,1)
    EndDo
    
    !2. Reading inflow/outflow
    Do i =1,HydroParam%NInflow
        iElem = HydroParam%IndexInflow(i,2)
        iEdge = HydroParam%IndexInflow(i,1)
        Do iLayer = HydroParam%InFlowSmallm(i), HydroParam%InFlowCapitalM(i)
            If (HydroParam%InflownTime(i)==-999) Then
                If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
                    HydroParam%Fu(iLayer,iEdge) = (Sig(iElem,MeshParam%Right(iEdge),MeshParam%Left(iEdge)))*(Max(HydroParam%Pcri,HydroParam%H(iEdge))**(2./3.))*(HydroParam%InflowValue(i,1)**(1./2.))/(Max(HydroParam%Pcri,HydroParam%H(iEdge))**(1./6.)/(HydroParam%Rug(iElem)+NearZero))
                ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
                    HydroParam%Fu(iLayer,iEdge) = (Sig(iElem,MeshParam%Right(iEdge),MeshParam%Left(iEdge)))*(Max(HydroParam%Pcri,HydroParam%H(iEdge))**(2./3.))*(HydroParam%InflowValue(i,1)**(1./2.))/(HydroParam%Rug(iElem)+NearZero)
                ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
                    HydroParam%Fu(iLayer,iEdge) = (Sig(iElem,MeshParam%Right(iEdge),MeshParam%Left(iEdge)))*(Max(HydroParam%Pcri,HydroParam%H(iEdge))**(2./3.))*(HydroParam%InflowValue(i,1)**(1./2.))/(Max(HydroParam%Pcri,HydroParam%H(iEdge))**(1./6.)/((18.*log10(12.*Max(HydroParam%Pcri,HydroParam%H(iEdge))/(HydroParam%Rug(iElem)/30.+NearZero)))+NearZero)) 
                EndIf
                
            Else
                Call interp_linear( 1, HydroParam%InflownTime(i), HydroParam%InflowTime(i,1:HydroParam%InflownTime(i)), HydroParam%InflowValue(i,1:HydroParam%InflownTime(i)), 1, t_interp, p_interp )
                HydroParam%Fu(iLayer,iEdge) = -(Sig(iElem,MeshParam%Right(iEdge),MeshParam%Left(iEdge)))*(p_interp(1,1))/(sum(HydroParam%DZj(HydroParam%InFlowSmallm(i):HydroParam%InFlowCapitalM(i),iEdge))*MeshParam%EdgeLength(iEdge))
            EndIf
        EndDo
    EndDo
    
    !!------------------------------------------------------------------------- 
    !    ! Vertedor
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
    