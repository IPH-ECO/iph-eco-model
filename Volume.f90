Subroutine Volume(HydroParam,MeshParam)

    ! Calculate the Element Volume
    ! Based on:
    ! [1] Casulli, V. A high-resolution wetting and drying algorithm for free-surfacehydrodynamics. 
    !   International Journal for Numerical Methods in Fluids, v. 60 (4), p. 391-408, 2009.
    ! [2] Casulli, V. A conservative semi-implicit method for coupled surface–subsurface flows in regional scale
    !   International Journal for Numerical Methods in Fluids, v. 79, p. 199-214, 2015.
    
    ! Input:
    ! Hydrodynamic Features
    ! Output:
    ! Element Volume
    
    ! List of Modifications:
    !   15.06.2015: Routine Implementation       (Carlos Ruberto)
    !   21.05.2020: Routine Update               (Cayo Lopes)
    ! Programmer: Cayo Lopes
  
    Use MeshVars !, Only: nElem,Edge,Area
    Use Hydrodynamic
    Implicit None
    Integer:: iElem, iLayer
    Real:: V
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam
    
    Do iElem = 1, MeshParam%nElem
        HydroParam%Vol(iElem) = 0.d0
        If (V(HydroParam%eta(iElem)+HydroParam%etaplus(iElem),HydroParam%hb(iElem)) > 0) Then
            !HydroParam%Vol(iElem) = MeshParam%Area(iElem)*(HydroParam%eta(iElem) + HydroParam%etaplus(iElem))
            HydroParam%Vol(iElem) = MeshParam%Area(iElem)*(HydroParam%eta(iElem) + HydroParam%etaplus(iElem) - HydroParam%hb(iElem))
            If (HydroParam%DZsi(HydroParam%Smallms(iElem),iElem) > 0) Then
                HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*Dot_Product(MeshParam%ei(:,iElem),HydroParam%DZsi(:,iElem))
            EndIf
        ElseIf (V(HydroParam%eta(iElem) + HydroParam%etaplus(iElem),HydroParam%sb(iElem)) > 0) Then
            HydroParam%Vol(iElem) = MeshParam%Area(iElem)*(V(HydroParam%eta(iElem)+HydroParam%etaplus(iElem),HydroParam%sb(iElem)))*MeshParam%ei(HydroParam%Smallms(iElem),iElem)
        EndIf
    EndDo
    
    !Do iElem = 1, MeshParam%nElem
    !    HydroParam%Vol(iElem) = 0
    !    If (V(HydroParam%eta(iElem)+HydroParam%etaplus(iElem),HydroParam%hb(iElem)) > 0) Then
    !        HydroParam%Vol(iElem) = MeshParam%Area(iElem)*(HydroParam%eta(iElem) + HydroParam%etaplus(iElem) - HydroParam%hb(iElem))
    !        If (MeshParam%iBedrock == 1) Then
    !            Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalMs(iElem)
    !                If (HydroParam%DZsi(iLayer,iElem) > 0) Then
    !                    HydroParam%Vol(iElem) = HydroParam%Vol(iElem) + MeshParam%Area(iElem)*MeshParam%ei(iLayer,iElem)*HydroParam%DZsi(iLayer,iElem)
    !                EndIf
    !            EndDo
    !        EndIf
    !    ElseIf (V(HydroParam%eta(iElem)+HydroParam%etaplus(iElem),HydroParam%sb(iElem)) > 0) Then
    !        HydroParam%Vol(iElem) = MeshParam%Area(iElem)*(V(HydroParam%eta(iElem)+HydroParam%etaplus(iElem),HydroParam%sb(iElem)))*MeshParam%ei(HydroParam%Smallms(iElem),iElem)
    !    EndIf
    !EndDo
    
    Return
End Subroutine Volume