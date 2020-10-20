!> This subroutine reads the simulation parameters. 
Subroutine AllocateHydroVars(HydroParam,MeshParam)
    
    Use Hydrodynamic
    Use MeshVars
    !Use Transport
    
    Implicit none
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam

    ! 1. GridData variables
    Allocate(HydroParam%hb(MeshParam%nElem))
    Allocate(HydroParam%sb(MeshParam%nElem))
    Allocate(HydroParam%Rug(MeshParam%nElem))
    Allocate(MeshParam%CREDV(MeshParam%nElem))
    Allocate(MeshParam%BANHADO(MeshParam%nElem))
    Allocate(MeshParam%d50(MeshParam%nElem))
    Allocate(MeshParam%OMfraction(MeshParam%nElem))
    Allocate(MeshParam%eta0(MeshParam%nElem))
    
    ! 2. Hydrodynamic variables
    ! 2.1. Velocities
    Allocate (HydroParam%ug(MeshParam%nEdge,MeshParam%kMax+1))
    Allocate (HydroParam%vg(MeshParam%nEdge,MeshParam%kMax+1))
    Allocate (HydroParam%wg(MeshParam%nEdge,MeshParam%kMax+1))
    Allocate(HydroParam%Fu(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%FuxyNode(MeshParam%Kmax,3,MeshParam%nNode))
    Allocate(HydroParam%Fub(MeshParam%Kmax,3,MeshParam%nElem))
    Allocate(HydroParam%Fv(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%Fvu(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%Fuv(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%u(MeshParam%Kmax,MeshParam%nEdge))
    
    Allocate(HydroParam%psij(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%rj(MeshParam%Kmax,MeshParam%nEdge))
    
    Allocate(HydroParam%utang(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%Wu(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%uxyback(MeshParam%Kmax,2,MeshParam%nEdge))
    Allocate(HydroParam%uArrow(MeshParam%Kmax,3,MeshParam%nEdge))    
    Allocate(HydroParam%ubBack(MeshParam%Kmax,3,MeshParam%nElem))
    Allocate(HydroParam%uxy(MeshParam%Kmax,2,MeshParam%nEdge))
    Allocate(HydroParam%uxyL(MeshParam%Kmax+1,2,MeshParam%nElem))
    Allocate(HydroParam%uNode(MeshParam%Kmax+1,3,MeshParam%nNode))
    Allocate(HydroParam%ut(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%w(MeshParam%Kmax+1,MeshParam%nElem))
    Allocate(HydroParam%wt(MeshParam%Kmax+1,MeshParam%nElem))
    Allocate(HydroParam%wfc(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%Fw(MeshParam%Kmax+1,MeshParam%nElem))
    Allocate(HydroParam%ub(MeshParam%Kmax,3,MeshParam%nElem))
    Allocate(HydroParam%ubV(MeshParam%Kmax,3,MeshParam%nNode))
    Allocate(HydroParam%epson(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%psi_edge(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%psi_cell(MeshParam%Kmax,MeshParam%nElem))

    Allocate (HydroParam%ugt(MeshParam%nEdge,MeshParam%kMax+1))
    Allocate (HydroParam%vgt(MeshParam%nEdge,MeshParam%kMax+1))
    Allocate (HydroParam%wgt(MeshParam%nEdge,MeshParam%kMax+1)) 
    Allocate(HydroParam%ubt(MeshParam%Kmax,3,MeshParam%nElem))
    Allocate(HydroParam%ubVt(MeshParam%Kmax,3,MeshParam%nNode))
    Allocate(HydroParam%uxyt(MeshParam%Kmax,2,MeshParam%nEdge))
    Allocate(HydroParam%uxyLt(MeshParam%Kmax+1,2,MeshParam%nElem))
    Allocate(HydroParam%uNodet(MeshParam%Kmax+1,3,MeshParam%nNode))
    Allocate(HydroParam%wfct(MeshParam%Kmax,MeshParam%nEdge))
    
    
    ! 2.2. Others Variables
    Allocate(HydroParam%etaInf(MeshParam%nElem))
    Allocate(HydroParam%etaplus(MeshParam%nElem))
    Allocate(HydroParam%peta(MeshParam%nNode))
    Allocate(HydroParam%petan(MeshParam%nNode))
    Allocate(HydroParam%eta(MeshParam%nElem))
    Allocate(HydroParam%etan(MeshParam%nElem))
    Allocate(HydroParam%P(MeshParam%nElem))
    Allocate(HydroParam%Aeta(MeshParam%nElem))
    Allocate(HydroParam%f(MeshParam%nElem))
    Allocate(HydroParam%Deta(MeshParam%nElem))
    Allocate(HydroParam%rhs(MeshParam%nElem))
    Allocate(HydroParam%Gu(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%uIniVec(MeshParam%Kmax,3))
    Allocate(HydroParam%DZj(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%DZjt(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%DZi(MeshParam%Kmax,MeshParam%nElem))
    Allocate(HydroParam%DZit(MeshParam%Kmax,MeshParam%nElem))
    Allocate(HydroParam%iADZ(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%iAG(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%Z(MeshParam%KMax+1,MeshParam%nEdge))
    Allocate(HydroParam%Ze(MeshParam%KMax+1,MeshParam%nElem))    
    Allocate(HydroParam%Zb(MeshParam%KMax,MeshParam%nElem))  
    Allocate(HydroParam%hj(MeshParam%nEdge))
    Allocate(HydroParam%sj(MeshParam%nEdge))
    Allocate(HydroParam%H(MeshParam%nEdge))
    Allocate(HydroParam%Hu(MeshParam%nEdge))
    Allocate(HydroParam%Smallm(MeshParam%nEdge)) !lower vertical index in the edge 
    Allocate(HydroParam%Smallms(MeshParam%nEdge)) !lower vertical index in the edge 
    Allocate(HydroParam%CapitalM(MeshParam%nEdge)) !upper vertical index in the edge
    Allocate(HydroParam%CapitalMs(MeshParam%nEdge)) !upper vertical index in the edge !CAYO
    Allocate(HydroParam%ElSmallm(MeshParam%nElem)) !lower vertical index in the cell
    Allocate(HydroParam%ElSmallms(MeshParam%nElem)) !lower vertical index in the cell
    Allocate(HydroParam%ElCapitalM(MeshParam%nElem)) !upper vertical index in the cell
    Allocate(HydroParam%ElCapitalMs(MeshParam%nElem)) !upper vertical index in the cell !CAYO
    Allocate(HydroParam%DZiADZ(MeshParam%nEdge))
    Allocate(HydroParam%DZiAG(MeshParam%nEdge))
    Allocate(HydroParam%HorViscosity(2,MeshParam%KMax,MeshParam%nElem))
    Allocate(HydroParam%HorDiffusivity(2,MeshParam%KMax,MeshParam%nElem))
    Allocate(HydroParam%VerEddyDiff(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%VerEddyVisc(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%VerEddyDiffCell(MeshParam%Kmax,MeshParam%nElem))
    Allocate(HydroParam%VerEddyViscCell(MeshParam%Kmax,MeshParam%nElem))
    Allocate(HydroParam%WindVel(2,MeshParam%nEdge))
    Allocate(HydroParam%WindXY(2,MeshParam%nEdge))
    Allocate(HydroParam%Windix(MeshParam%nElem))
    Allocate(HydroParam%Windiy(MeshParam%nElem))
    Allocate(HydroParam%fetch_m(MeshParam%nElem,8))
    Allocate(HydroParam%PBarc(MeshParam%Kmax,MeshParam%nEdge))
    Allocate(HydroParam%DZitau(MeshParam%Kmax,MeshParam%nElem))
    Allocate(HydroParam%DZistau(MeshParam%Kmax,MeshParam%nElem))
    Allocate(HydroParam%Locdt(MeshParam%Kmax,MeshParam%nElem))
    Allocate(HydroParam%uLoadVarEst(MeshParam%KMax,MeshParam%nElem))
    Allocate(HydroParam%dVarEst(MeshParam%KMax,MeshParam%nElem,4))
    Allocate(HydroParam%sDRhoW(MeshParam%Kmax,MeshParam%nElem))
    Allocate(HydroParam%sDRhoWt(MeshParam%Kmax,MeshParam%nElem))
    
    
    Allocate(HydroParam%rhsnonHydro(MeshParam%Kmax+1,MeshParam%nElem))
    Allocate(HydroParam%q(MeshParam%Kmax+1,MeshParam%nElem))
    Allocate(HydroParam%pq(MeshParam%Kmax,MeshParam%nNode))
    !Allocate(HydroParam%q0(MeshParam%Kmax-1,MeshParam%nElem))
    !Allocate(HydroParam%Dq(MeshParam%Kmax,MeshParam%nElem))
    !Allocate(HydroParam%Fq(MeshParam%Kmax,MeshParam%nElem))
    
    Allocate(HydroParam%Vol(MeshParam%nElem)) !CAYO
    Allocate(MeshParam%ei(MeshParam%Kmax,MeshParam%nElem))!CAYO
    Allocate(MeshParam%Kj(MeshParam%Kmax,MeshParam%nEdge))!CAYO
    Allocate(HydroParam%us(MeshParam%Kmax,MeshParam%nEdge)) !CAYO
    Allocate(HydroParam%ust(MeshParam%Kmax,MeshParam%nEdge)) !CAYO
    Allocate(HydroParam%um(MeshParam%Kmax,MeshParam%nEdge))!CAYO
    Allocate(HydroParam%umt(MeshParam%Kmax,MeshParam%nEdge))!CAYO
    Allocate(HydroParam%wm(MeshParam%Kmax+1,MeshParam%nElem))
    Allocate(HydroParam%wmt(MeshParam%Kmax+1,MeshParam%nElem))
    
    Allocate(HydroParam%uxysub(MeshParam%Kmax,2,MeshParam%nEdge))
    Allocate(HydroParam%ubsub(MeshParam%Kmax,3,MeshParam%nElem))
    
    Allocate(HydroParam%DZsj(MeshParam%Kmax,MeshParam%nEdge)) !CAYO   
    Allocate(HydroParam%DZsjt(MeshParam%Kmax,MeshParam%nEdge)) !CAYO
    Allocate(HydroParam%DZhj(MeshParam%Kmax,MeshParam%nEdge)) !CAYO   
    Allocate(HydroParam%DZhjt(MeshParam%Kmax,MeshParam%nEdge)) !CAYO
    
    Allocate(HydroParam%DZsi(MeshParam%Kmax,MeshParam%nElem)) !CAYO   
    Allocate(HydroParam%DZsit(MeshParam%Kmax,MeshParam%nElem)) !CAYO
    Allocate(HydroParam%DZhi(MeshParam%Kmax,MeshParam%nElem)) !CAYO   
    Allocate(HydroParam%DZhit(MeshParam%Kmax,MeshParam%nElem)) !CAYO 
    
    Allocate(HydroParam%DZK(MeshParam%nEdge)) !Sediment Layer !CAYO 
    
    ! 3. Hydrodynamic output variables (VTK)
    Allocate(MeshParam%xPoint(MeshParam%nPoint*(MeshParam%Kmax+1)))
    Allocate(MeshParam%yPoint(MeshParam%nPoint*(MeshParam%Kmax+1)))
    Allocate(MeshParam%zPoint(MeshParam%nPoint*(MeshParam%Kmax+1)))
    Allocate(MeshParam%Connect(9*MeshParam%nElem*MeshParam%Kmax))
    Allocate(MeshParam%cell_type(MeshParam%nElem*MeshParam%Kmax))
    Allocate(HydroParam%SScalar(MeshParam%nElem*MeshParam%Kmax))
    Allocate(HydroParam%SScalar2D(MeshParam%nElem))
    Allocate(HydroParam%SVector(MeshParam%nElem*MeshParam%Kmax,3))
    
    
    ! 4. Hydrodynamic boundary condition variables
    Allocate(HydroParam%IndexInflowEdge(MeshParam%nEdge))
    Allocate(HydroParam%IndexWaterLevelEdge(MeshParam%nEdge))
    
   
End Subroutine AllocateHydroVars
    