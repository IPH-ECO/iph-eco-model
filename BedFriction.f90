Subroutine BedFriction(HydroParam,MeshParam,dt)
    
    !$ use omp_lib
    Use MeshVars !, Only: nEdge,Left,Right,CirDistance
    Use Hydrodynamic
    
    Implicit None
    Integer:: iElem, iEdge, iLayer
    Integer:: l, r, n1, n2 , n3, n4
    Real:: NearZero = 1e-10
    Real:: dt, ub, vb, GammaB
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam

  
    Do iElem = 1, MeshParam%nElem
           
        n1 = MeshParam%Edge(1,iElem)
        n2 = MeshParam%Edge(2,iElem)
        n3 = MeshParam%Edge(3,iElem)
        n4 = MeshParam%Edge(4,iElem)
        
        ub = 0.0d0
        r = MeshParam%right(n4)
        If (r/=0) Then
            ub = HydroParam%Fu(HydroParam%Smallm(n4),n4) - (dt/MeshParam%CirDistance(n4))*(HydroParam%g*(HydroParam%eta(r) - HydroParam%eta(iElem)))
        ElseIf (HydroParam%IndexWaterLevelEdge(n4) > 0) Then
            ub = HydroParam%Fu(HydroParam%Smallm(n4),n4) - (dt/MeshParam%CirDistance(n4))*(HydroParam%g*(HydroParam%etaInf(n4) - HydroParam%eta(iElem)))
        EndIf

        vb = 0.0d0
        r = MeshParam%right(n1)
        If (r/=0) Then
            vb = HydroParam%Fu(HydroParam%Smallm(n1),n1) - (dt/MeshParam%CirDistance(n1))*(HydroParam%g*(HydroParam%eta(r) - HydroParam%eta(iElem)))
        ElseIf (HydroParam%IndexWaterLevelEdge(n1)>0) Then
            vb = HydroParam%Fu(HydroParam%Smallm(n1),n1) - (dt/MeshParam%CirDistance(n1))*(HydroParam%g*(HydroParam%etaInf(n1) - HydroParam%eta(iElem)))
        EndIf      
        
        call EdgeGammab(n1, ub, vb, HydroParam, MeshParam)
        call EdgeGammab(n2, ub, vb, HydroParam, MeshParam)

        r = MeshParam%right(n2)
        If (r==0) Then
            HydroParam%GammaB(n2) = 0.5*(HydroParam%GammaB(n1) + HydroParam%GammaB(n4))
        EndIf
        
        r = MeshParam%right(n3)
        If (r==0) Then
            HydroParam%GammaB(n3) = 0.5*(HydroParam%GammaB(n1) + HydroParam%GammaB(n4))
        EndIf
        
    EndDo
    
    End Subroutine BedFriction
    
    
    
    Subroutine EdgeGammab(iEdge, ub, vb, HydroParam, MeshParam)
       
    Use MeshVars !, Only: 
    Use Hydrodynamic ! Only:    
    Implicit none
    
    Integer:: iElem, iEdge, iLayer
    Integer:: l, r
    Real :: H, Chezy, rhoair, ub, vb
    Real:: NearZero = 1e-10
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam

    l = MeshParam%Left(iEdge)
    r = MeshParam%Right(iEdge)    
    HydroParam%GammaB(iEdge) =  0.d0
    
    ! 7.1 Get roughness 
    H = HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge) !Surface Water Height
    If (r == 0) Then
        If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
            Chezy = HydroParam%Rug(l)
        ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
            Chezy = Max(HydroParam%Pcri,H)**(1./6.)/(HydroParam%Rug(l)+NearZero)
        ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
            Chezy = 18.*log10(12.*Max(HydroParam%Pcri,H)/(HydroParam%Rug(l)/30.+NearZero))
        EndIf
    Else
        If (HydroParam%iRoughForm == 0.or.HydroParam%iRoughForm == 3) Then ! roughnessChezyConstant
            Chezy = 0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))
        ElseIf (HydroParam%iRoughForm == 1.or.HydroParam%iRoughForm == 4) Then ! roughnessManningConstant
            Chezy = Max(HydroParam%Pcri,H)**(1./6.)/(0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))+NearZero)
        ElseIf (HydroParam%iRoughForm == 2.or.HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookConstant
            Chezy = 18.*log10(12.*Max(HydroParam%Pcri,H)/(0.5*(HydroParam%Rug(l) + HydroParam%Rug(r))/30.+NearZero))
        EndIf
    EndIf
    
    If (Chezy > 0 .and. HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge) > HydroParam%PCRI+NearZero) Then 
        HydroParam%GammaB(iEdge) = (HydroParam%g*(ub**2.+vb**2.)**0.5)/(Chezy**2.)
    EndIf        
    
    End Subroutine EdgeGammab