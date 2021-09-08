!>\brief Routines related to the Eulerian-Lagragean Method (FU, FV and FW) of TRIM model - Casulli
!>\details 
! References:
! [1]
! [2] Cunha, A.H.F.; Fragoso, C.R.; Tavares, M.H.; Cavalcanti, J.R.; Bonnet, M.-P.; Motta-Marques, D. Combined Use of High-Resolution Numerical 
!     Schemes to Reduce Numerical Diffusion in Coupled Nonhydrostatic Hydrodynamic and Solute Transport Model. Water 2019, 11, 2288.
!>\author Carlos Ruberto Fragoso and Rafael Cavalcanti
  Module ELM
  
  private
  public:: FuFv
  public:: ComputeFW
  public:: signa
  !-----------------------------------------------------------------
  Contains
  
  !> Main subroutine of Eulerian-Lagragean Method
  !>\return Fu, Fv and Fw the resultant vectors for each edge
  !>\author  Carlos Ruberto Fragoso
  !>\attention  List of Modifications: \n
  !> - 15.09.2015: Routine Implementation (Carlos Ruberto Fragoso)
    
    Subroutine FuFv(HydroParam,MeshParam,dt)
    
    Use MeshVars !, Only: Left,Right,nElem,nEdge,EdgeBary,NormalVector,TangentVector,Edge,EdgeLength,Quadri,xNode,yNode,Area,xb,yb,nNode,Kmax,nEgdesatNode,EgdesatNode
    Use Hydrodynamic !, Only: Smallm,CapitalM,Z,u,utang,Fu,Fv,H,Pcri,w,Nfut,uxy,uxyback,ElSmallm,ElCapitalM,lat,OMEGA,Pi,Wu,HorViscosity,Fub,FuxyNode
    !Use SimulationModel 
    
    Implicit None
    Integer:: iElem,iEdge,lEdge,iLayer,jlev,l,r,nnel,ndels,Sig,Face,n1,n2,n3,n4,iNode,icount,j,fac,kin,TrajectoryFlag,FuFw_flag,jump, ntals, iEdge0
    Real::x0,y0,z0,xt,yt,zt,uuint,vvint,wwint,hhint,wdown,wup,vmag,dtb,aa(4),weit,sinal,dtbCFL
    Real:: NearZero = 1e-10
    Real::dt, tal
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam

    ! 1.1 Find out uxy velocities components, that are equivalents to normal and tangencial velocities components in particle position after
    ! the backtrack process:
    Do iEdge = 1,MeshParam%nEdge
        !Get right/left iElements that share iEdge
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
        
        ! If is a dry cell, the backtrack velocities is equal to zero.
        ! Note that in subsurface coupled case %H(iEdge) can takes into account the freatic water level in Cell.
        ! DZsj represents the Edge freatic component thickness, thus H - DZsj is equivalent only surface water level
        ! In Cell where no subsufarce layer exist or in only surface flow simulation DZsj is set equal 0.0d0 (ReadHydroIniCond Module)
        If (HydroParam%H(iEdge) + HydroParam%sj(iEdge)-HydroParam%hj(iEdge)<=HydroParam%Pcri+NearZero) Then !CAYO
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                HydroParam%uxyback(iLayer,1:2,iEdge) = (/ 0., 0. /)
            EndDo
            Cycle
        EndIf
        
        ! In the case that Cell is Wet:
        Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
            ! The particle is set in iEdge barycenter (x0,y0,z0) on Left side Element (nnel):
            nnel = l 
            jlev = iLayer
            x0 = MeshParam%EdgeBary(1,iEdge)
            y0 = MeshParam%EdgeBary(2,iEdge)
            z0 = 0.5d0*(HydroParam%Z(iLayer+1,iEdge) + HydroParam%Z(iLayer,iEdge) + sum(HydroParam%DZsj(:,iEdge)))
            
            ! Velocities in the initial position:            
            uuint = HydroParam%uxy(iLayer,1,iEdge)
            vvint = HydroParam%uxy(iLayer,2,iEdge)
            wwint = HydroParam%wfc(iLayer,iEdge)
            hhint = HydroParam%H(iEdge)-HydroParam%hj(iEdge)
            ! Velocity magnitude:
            vmag = dsqrt(uuint**2+vvint**2+wwint**2)           

            ! In cases that the velocity magnitude in layer is very low or the iEdge no has neighbour Element (r==0), the backtracking 
            ! velocity in layer is set equal the iEdge:
            If(vmag.le.1.e-6) Then ! No activity or Water Level Condition
                HydroParam%uxyback(iLayer,1:2,iEdge) = (/ uuint, vvint /)    
            Elseif (HydroParam%eta(nnel) < HydroParam%hb(nnel) .and. HydroParam%eta(r) > HydroParam%hb(r)) Then
                HydroParam%uxyback(iLayer,1:2,iEdge) = (/ uuint, vvint /)   
            Else
            ! Else, the backtracking process is initialized:
                If (r==0) Then !.and.HydroParam%u(iLayer,iEdge)>0
                    HydroParam%uxyback(iLayer,1:2,iEdge) = (/ uuint, vvint /)
                Else
                    If (iLayer < Hydroparam%ElSmallm(r)) Then
                        HydroParam%uxyback(iLayer,1,iEdge) =  uuint
                        HydroParam%uxyback(iLayer,2,iEdge) =  vvint
                    Else 
                        ! If TrajectoryFlag == 0 System Defined Intervals to Integrate Particle Trajectory,  
                        ! else User Defined Intervals to Integrate Particle Trajectory. 
                        !TrajectoryFlag = 0
                        ! Find the Local Time Step (dtb) to Integrate Trajectory:
                        If ( TrajectoryFlag == 0 ) Then         ! System Calculation - [7,8]
                            dtb   = MinVal(MeshParam%InCircle)/vmag ! ( InCircle(lElem)/Sqrt( Veloc(1)**2. + Veloc(2)**2. ), InCircle(r)/Sqrt( Veloc(1)**2. + Veloc(2)**2. ) )    
                            dtbCFL = Max((1/HydroParam%CFL(l)),(1/HydroParam%CFL(r)))
                            ndels = Max(Floor(dt/dtb),HydroParam%NFUT,Floor(dt/dtbCFL))
                            !ndels = Max(Floor(dt/dtb),HydroParam%NFUT)
                            dtb = dt/ndels
                        ElseIf ( TrajectoryFlag == 1 ) Then
                            dtbCFL = Max((1/HydroParam%CFL(l)),(1/HydroParam%CFL(r)))
                            ndels = Max(HydroParam%NFUT,Floor(dt/dtbCFL))
                            !ndels = HydroParam%NFUT
                            dtb = dt/ndels
                        EndIf

                        FuFw_flag = 0
                        iEdge0 = iEdge
                        Call btrack(ndels,dtb,dt,uuint,vvint,wwint,hhint,&
                        &x0,y0,z0,xt,yt,zt,nnel,jlev,iEdge0,HydroParam,MeshParam, FuFw_flag, iLayer, iEdge)

                        HydroParam%uxyback(iLayer,1:2,iEdge) = (/ uuint, vvint /)
                        
                    EndIf
                EndIf
            Endif
            
        EndDo
        
    EndDo
    
    If (HydroParam%iNonHydro==0.and.MeshParam%Kmax>1) Then
        Call ComputeFW(MeshParam,HydroParam,dt)
    EndIf
    
    !Explicit terms
   
    Do iEdge=1,MeshParam%nEdge
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
        Do lEdge=1,4
            If (MeshParam%Edge(lEdge,l)==iEdge) Then
                Exit
            EndIf
        EndDo

        !If (lEdge == 1) Then
        !    Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
        !        HydroParam%Fv(iLayer,iEdge) = Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*HydroParam%uxyback(iLayer,1,iEdge)*MeshParam%TangentVector(1,lEdge,l) + Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*HydroParam%uxyback(iLayer,2,iEdge)*MeshParam%TangentVector(2,lEdge,l)
        !    EndDo        
        !ElseIf(lEdge == 2) Then
        !    Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
        !        HydroParam%Fv(iLayer,iEdge) = Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*HydroParam%uxyback(iLayer,1,iEdge)*MeshParam%TangentVector(1,lEdge,l) - Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*HydroParam%uxyback(iLayer,2,iEdge)*MeshParam%TangentVector(2,lEdge,l)
        !    EndDo   
        !Else ! For iEdge == 3 and 4
        !    Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
        !        HydroParam%Fv(iLayer,iEdge) = Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*(HydroParam%uxyback(iLayer,1,iEdge)*MeshParam%TangentVector(1,lEdge,l) + HydroParam%uxyback(iLayer,2,iEdge)*MeshParam%TangentVector(2,lEdge,l))
        !    EndDo                         
        !EndIf
        
        If (lEdge == 1) Then
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                HydroParam%Fv(iLayer,iEdge) = HydroParam%uxyback(iLayer,1,iEdge)*MeshParam%TangentVector(1,lEdge,l) + HydroParam%uxyback(iLayer,2,iEdge)*MeshParam%TangentVector(2,lEdge,l)
            EndDo        
        ElseIf(lEdge == 2) Then
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                HydroParam%Fv(iLayer,iEdge) = HydroParam%uxyback(iLayer,1,iEdge)*MeshParam%TangentVector(1,lEdge,l) - HydroParam%uxyback(iLayer,2,iEdge)*MeshParam%TangentVector(2,lEdge,l)
            EndDo   
        Else ! For iEdge == 3 and 4
            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
                HydroParam%Fv(iLayer,iEdge) = HydroParam%uxyback(iLayer,1,iEdge)*MeshParam%TangentVector(1,lEdge,l) + HydroParam%uxyback(iLayer,2,iEdge)*MeshParam%TangentVector(2,lEdge,l)
            EndDo                         
        EndIf        
        !
        ! If the face has a boundary condition:
        If (HydroParam%IndexInflowEdge(iEdge)>0.or.HydroParam%IndexWaterLevelEdge(iEdge)>0) then
            HydroParam%Fu(:,iEdge) = HydroParam%u(:,iEdge)
            HydroParam%Fv(:,iEdge) = HydroParam%utang(:,iEdge)
            cycle
        EndIf

        !MeshParam%Neighbor(:,l)
        Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
            If (r==0) Then !Velocity in lateral boundaries is equal to cell center velocity
                HydroParam%Fu(iLayer,iEdge) = HydroParam%u(iLayer,iEdge) !Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*(HydroParam%Fub(iLayer,1,l))*MeshParam%NormalVector(1,lEdge,l)  + Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*(HydroParam%Fub(iLayer,2,l))*MeshParam%NormalVector(2,lEdge,l)
            Else
                !sinal = Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))
                ! Do in the Edges - > sig func always return 1, because right and left elements in Edge doesn't change (l == left always) - > check if is necessary
                if (iLayer < HydroParam%ElSmallm(r)) Then
                    HydroParam%Fu(iLayer,iEdge) = HydroParam%u(iLayer,iEdge)
                else
                    HydroParam%Fu(iLayer,iEdge) = Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*HydroParam%uxyback(iLayer,1,iEdge)*MeshParam%NormalVector(1,lEdge,l)  + Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*HydroParam%uxyback(iLayer,2,iEdge)*MeshParam%NormalVector(2,lEdge,l)
                    !HydroParam%Fu(iLayer,iEdge) = Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*HydroParam%uxyback(iLayer,1,iEdge)*MeshParam%NormalVector(1,lEdge,l)  + (Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*HydroParam%uxyback(iLayer,2,iEdge)*MeshParam%NormalVector(2,lEdge,l)*dt)*((HydroParam%uxyback(iLayer,2,iEdge-3)-2*HydroParam%uxyback(iLayer,2,iEdge)+HydroParam%uxyback(iLayer,2,iEdge+3))/MeshParam%EdgeLength(iEdge))
                endif
            EndIf
            HydroParam%Fv(iLayer,iEdge) = HydroParam%uxyback(iLayer,1,iEdge)*MeshParam%TangentVector(1,lEdge,l) + HydroParam%uxyback(iLayer,2,iEdge)*MeshParam%TangentVector(2,lEdge,l)
        EndDo
    EndDo

    !! When we try to calculating Fu and Fw with the particle starting from the center of the cell
    !Do iElem = 1,MeshParam%nElem
    !    n1=MeshParam%Quadri(1,iElem) + 1
    !    n2=MeshParam%Quadri(2,iElem) + 1
    !    n3=MeshParam%Quadri(3,iElem) + 1
    !    n4=MeshParam%Quadri(4,iElem) + 1
    !    aa(1)=signa(MeshParam%xNode(n1),MeshParam%xNode(n2),MeshParam%xb(iElem),MeshParam%yNode(n1),MeshParam%yNode(n2),MeshParam%yb(iElem))
    !    aa(2)=signa(MeshParam%xNode(n2),MeshParam%xNode(n3),MeshParam%xb(iElem),MeshParam%yNode(n2),MeshParam%yNode(n3),MeshParam%yb(iElem))
    !    aa(3)=signa(MeshParam%xNode(n3),MeshParam%xNode(n4),MeshParam%xb(iElem),MeshParam%yNode(n3),MeshParam%yNode(n4),MeshParam%yb(iElem))
    !    aa(4)=signa(MeshParam%xNode(n4),MeshParam%xNode(n1),MeshParam%xb(iElem),MeshParam%yNode(n4),MeshParam%yNode(n1),MeshParam%yb(iElem))
    !    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
    !        HydroParam%Fub(iLayer,1:2,iElem) = 0.
    !        Do iEdge = 1,4
    !            HydroParam%Fu(iLayer,MeshParam%Edge(iEdge,iElem)) = Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%uxyback(iLayer,1,MeshParam%Edge(iEdge,iElem))*MeshParam%NormalVector(1,iEdge,iElem)  + Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%uxyback(iLayer,2,MeshParam%Edge(iEdge,iElem))*MeshParam%NormalVector(2,iEdge,iElem)
    !            HydroParam%Fv(iLayer,MeshParam%Edge(iEdge,iElem)) = Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%uxyback(iLayer,1,MeshParam%Edge(iEdge,iElem))*MeshParam%TangentVector(1,iEdge,iElem) + Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%uxyback(iLayer,2,MeshParam%Edge(iEdge,iElem))*MeshParam%TangentVector(2,iEdge,iElem)
    !            HydroParam%Fub(iLayer,1,iElem) = HydroParam%ubBack(iLayer,1,iElem) !HydroParam%Fub(iLayer,1,iElem) + HydroParam%uxyback(iLayer,1,MeshParam%Edge(iEdge,iElem))*aa(iEdge)/MeshParam%Area(iElem)
    !            HydroParam%Fub(iLayer,2,iElem) =  HydroParam%ubBack(iLayer,2,iElem) !HydroParam%Fub(iLayer,2,iElem) + HydroParam%uxyback(iLayer,2,MeshParam%Edge(iEdge,iElem))*aa(iEdge)/MeshParam%Area(iElem)
    !        EndDo
    !    EndDo
    !    
    !EndDo
     
    HydroParam%FuxyNode = 0.
    Do iNode=1,MeshParam%nNode
        Do iLayer = 1,MeshParam%Kmax
            weit=0
            icount=0
            
            Do j=1,MeshParam%nEgdesatNode(iNode)
                Face = MeshParam%EgdesatNode(j,iNode)
                l = MeshParam%Left(Face)
                r = MeshParam%Right(Face)
                If (r==0) Then
                    fac=2
                Else
                    fac=1
                EndIf
                If(iLayer>=HydroParam%Smallm(Face)) Then
                    kin=min(iLayer,HydroParam%CapitalM(Face))
                    HydroParam%FuxyNode(kin,1,iNode)=HydroParam%FuxyNode(kin,1,iNode)+HydroParam%uxyback(kin,1,Face)/MeshParam%EdgeLength(Face)*fac
                    HydroParam%FuxyNode(kin,2,iNode)=HydroParam%FuxyNode(kin,2,iNode)+HydroParam%uxyback(kin,2,Face)/MeshParam%EdgeLength(Face)*fac
                    icount=icount+1 
                EndIf
                weit=weit+fac/MeshParam%EdgeLength(Face)
            EndDo
            If(icount.ne.0) Then
                HydroParam%FuxyNode(iLayer,1,iNode) = HydroParam%FuxyNode(iLayer,1,iNode)/weit
                HydroParam%FuxyNode(iLayer,2,iNode) = HydroParam%FuxyNode(iLayer,2,iNode)/weit
            EndIf
            
        EndDo       
    EndDo
    
    Return
    End Subroutine FuFv
 

    Subroutine btrack(ndelt,dtb,dt,uuint,vvint,wwint,hhint,x0,y0,z0,&
    &xt, yt, zt, nnel, jlev, id0, HydroParam, MeshParam,FuFw_flag, iLayer, iEdge)  
    !> Routine for backtracking.
    
    !> Input: point P0(x0,y0,z0) and nnel that encloses it, and initial vel.,
    !> initial level (for quicksearch), and a flag indicating 1st or 2nd tracking.\n
    !> Output: destination point Pt(xt,yt,zt), element nnel and level jlev,
    !> and vel. there (uuint,vvint,wwint).
    !>\param imode mode of backtracking (imode = 1: hydrodynamic)
    !>\param ndelt number of sub step time (input)
    !>\param dtb sub step time in seg (input)
    !>\param uuint x velocity component in the center of the edge (input)
    !>\param vvint y velocity component in the center of the edge (input)
    !>\param wwint z velocity component in the center of the edge (input)
    !>\param nnel number of the element that enclosed the point P0 (input)
    !>\param jlev number of the layer that enclosed the point P0 (input)
    !>\param x0 x initial position of the point P0 (input)
    !>\param y0 y initial position of the point P0(input)
    !>\param z0 z initial position of the point P0 (input)
    !>\return xt: x final position of the point Pt (output)
    !>\return yt: y final position of the point Pt (output)
    !>\return zt: z final position of the point Pt (output)
    !>\return uuint: x velocity component in the Edge (output)
    !>\return vvint: y velocity component in the Edge (output)
    !>\return wwint: z velocity component in the Edge (output)
    !>\return nnel: number of the element that enclosed the point Pt (output)
    !>\return jlev: number of the layer that enclosed the point Pt (output)
    
    !>\author  Carlos Ruberto Fragoso
    !>\attention  List of Modifications: \n
    !> - 15.09.2015: Routine Implementation (Carlos Ruberto Fragoso)

    Use MeshVars !, Only: Quadri,EdgeDef,xNode,yNode,Area,Edge,EdgeNodes,EdgeBary,Left,Right
    Use Hydrodynamic !, Only: Smallm,ElSmallm,ElCapitalM,CapitalM,Z,DZj,Ze,DZi,uNode,uxy,Ze,DZi,w
    
    Implicit none
    
    Real, parameter :: small1=1e-6
    integer, intent(in) :: ndelt, FuFw_flag, iLayer
    Real, intent(inout) :: dtb, dt
    integer, intent(inout) :: nnel,jlev, id0
    Real, intent(inout) :: uuint,vvint,wwint,x0,y0,z0, hhint
    Real, intent(out) :: xt,yt,zt
    Real :: uuNode(3,9), vvNode(3,9), wwNode(3,9), uuNodet(3,9), vvNodet(3,9), wwNodet(3,9), xxNode(3,9), yyNode(3,9), zzNode(3,9), hhNode(3,9)
    Real ::  uuBtrack, vvBtrack, wwBtrack, hhBTrack, uuBtrack2, vvBtrack2, wwBtrack2, theta, Umax, Umin
    Real :: staint(8),t_xi(4),s_xi(4),sig(4),subrat(4)
    Integer:: IntersectFlag,nel,nnelIni,idt,iflqs1,i,j,l,nd,nn,lev, jjlev,n1,n2,n3,n4, n5, n6, n7, n8, n9,iNode1,iNode2,iEdge,iNode, Nodes(9), addLayer, N, S, Face, r, ll, nnel0,psi_flag,nel_j
    Real:: trat,zup,zrat,aa1,aa2,aa3,aa4,aa,csi,etta,wdown,wup,weit, xtaux, ytaux, ztaux, dtaux
    Real:: NearZero = 1e-10 !< Small Number
    Real:: talx, taly, talz, tal, timeAcum, dtin,dtbCFL
    Integer:: i34 = 4
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Integer:: ELM_flag = 1
    Integer :: Interpolate_Flag = 1
    Real:: xaux,yaux,zaux
    
    jjlev = jlev
    nnelIni = nnel
    nnel0 = nnel
    
    dtaux = dtb
    timeAcum = 0.d0
    IntersectFlag = 0
    
    Do While (timeAcum < dt)
        !Posição da partícula no primeiro subpasso de tempo
        xt=x0-dtb*uuint
        yt=y0-dtb*vvint
        zt=z0-dtb*wwint 
        
        !Call RK4order(uuint, vvint, wwint, dtb, nnel, id0, jlev, xt, yt, zt, x0, y0, z0, dt, timeAcum + dtb, Interpolate_Flag, HydroParam, MeshParam)
        dtin = dtb
        Call quicksearch(1,nnel,jlev,dtb,dtin,x0,y0,z0,xt,yt,zt,iflqs1,idt,id0,i34,uuint,vvint,wwint,IntersectFlag,nel_j,HydroParam,MeshParam)   
        dtb = dtin
        
        if(isnan(dtb)) Then
            continue
        endif
        
        uuBtrack = 0.d0
        vvBtrack = 0.d0
        wwBtrack = 0.d0
        If (IntersectFlag == 0) Then
            If (ELM_flag==0) Then !iBilinear Interpolation
            
                timeAcum = timeAcum + dtb
                If (FuFw_flag==0) Then !grid for U velocitys
                    !Call FuVelocities2(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nnel,jlev, xt,yt,zt,x0,y0,z0, id0, HydroParam,MeshParam)
                    Call FuVelocities3(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nnel, jlev, xt,yt,zt,x0,y0,z0, id0, HydroParam,MeshParam)
                Else !grid for W velocitys
                    Call FwVelocities2(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nnel, jlev, jjlev, xt,yt,zt,x0,y0,z0, HydroParam,MeshParam)
                EndIf
                Call iBilinear2(uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt, zt, x0, y0, z0, id0, nnel, FuFw_flag, MeshParam)
        
            ElseIf (ELM_flag==1) Then !iQuadratic Interpolation
            
                timeAcum = timeAcum + dtb
                ! Get nodals' velocities and positions:
                Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), uuNodet(:,:), vvNodet(:,:), wwNodet(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nnel, jlev, HydroParam, MeshParam)     
                If (HydroParam%eta(nnel) - HydroParam%hb(nnel) < HydroParam%PCRI + NearZero) Then 
                    timeAcum = dt                
                Else
                    Call iQuadratic(uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt, zt)                        
                EndIf
            
            EndIf
        Else
            timeAcum = dt   
        EndIf

        uuint = uuBtrack
        vvint = vvBtrack
        wwint = wwBtrack           

        If(iflqs1.eq.1) exit

        x0=xt
        y0=yt
        z0=zt
                   
        !Adaptative sub-time step:
        If (nnel /= nnel0) Then
            
            dtb = dtaux
            tal = min(MeshParam%dx/abs(uuint), MeshParam%dy/abs(vvint), (zzNode(3,5) - zzNode(1,5))/abs(wwint))
            dtbCFL = 1/HydroParam%CFL(nnel)
            if (tal > 0 ) then
                !dtb = min(tal,dt,dtb)
                dtb = min(tal,dt,dtb,dtbCFL) 
            endif
            nnel0 = nnel
            If (timeAcum < dt .and. dtb + timeAcum > dt) Then
                dtb = dt - timeAcum
            EndIf
            
        Elseif( uuint**2 + vvint**2 + wwint**2 == 0.0d0 ) Then
            timeAcum = dt
        EndIf  
        
        if (IntersectFlag == 1) Then
            timeAcum = dt
            uuint = HydroParam%uxy(jlev,1,MeshParam%Edge(nel_j,nel))
            vvint = HydroParam%uxy(jlev,2,MeshParam%Edge(nel_j,nel))
        EndIf
        
    EndDo 
      
    Return
    End Subroutine btrack    

    Subroutine quicksearch(iloc,nnel,jlev,dtb,dtin,x0,y0,z0,&
    &xt,yt,zt,nfl,idt,id0,i34,uuint,vvint,wwint,IntersectFlag,nel_j,HydroParam,MeshParam)

      !> Straightline search algorithm.
      !>\note Initially nnel is an element that encompasses the point P0(x0,y0).\n
      !> Input: iloc,nnel,x0,y0,z0,xt,yt,zt,jlev, time, and vt2, ww2 for abnormal cases;\n
      !> Output: the updated end pt (xt,yt,zt) (if so), nnel, jlev, a flag nfl.\n
      !> Exit btrack if a bnd or dry element is hit and vel. there is small, or death trap is reached.
      !>\param iloc nudge initial pt (iloc=0: do not nudge initial pt; iloc=1: nudge)
      !>\param dtb sub step time in seg (input)
      !>\param nnel number of the element that enclosed the point P0 (input)
      !>\param jlev number of the layer that enclosed the point P0 (input)
      !>\param x0 x initial position of the point P0 (input)
      !>\param y0 y initial position of the point P0(input)
      !>\param z0 z initial position of the point P0 (input)
      !>\param xt x final position of the point Pt (input)
      !>\param yt y final position of the point Pt (input)
      !>\param zt z final position of the point Pt (input)
      !>\return nfl: The element that enclosed the point Pt was founded (output)
      !>\return nnel: number of the element that enclosed the point Pt (output)
      !>\return jlev: number of the layer that enclosed the point Pt (output)
      !>\author  Carlos Ruberto Fragoso
      !>\attention  List of Modifications: \n
      !> - 15.09.2015: Routine Implementation (Carlos Ruberto Fragoso)
    
    Use MeshVars !, Only:Quadri,Edge,EdgeDef,xNode,yNode,Area,xb,yb,Left,Right,Neighbor,EdgeBary
    Use Hydrodynamic !, Only: ElSmallm,ElCapitalM,Ze,H,Pcri,uxy,uNode
    Implicit none
    
    Real, parameter :: small1=1e-5
    Integer, intent(in) :: iloc,idt,i34
    Real, intent(in) :: dtb,x0,y0,z0
    Integer, intent(out) :: nfl
    Integer, intent(inout) :: nnel,jlev,id0
    Real:: xpoly(4),ypoly(4)
    Real, intent(inout) :: xt,yt,zt
    Real:: trm,aa,aa1,ae,xcg,ycg,pathl,xin,yin,zin,tt1,tt2,dist,xvel,yvel,zvel,hvel
    Real:: uuint,vvint,wwint,dtin
    Integer:: IntersectFlag
    Integer:: nel,i,j,k,n1,n2,nel_j,iflag,it,md1,md2,lit,k1,k2,jd1,jd2,r,isd,nel0,INOUT
    Real:: NearZero = 1e-10 !< Small Number
    Integer :: NWater
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam

    nel0 = nnel
    nfl=0
    trm=dtb !time remaining
    nel_j = id0
    nel = nnel  
    aa=0
    aa1=0
    !The ideia is that any polygon area can be represents by sum of multiples triangules.
    !Next are calculated the multiples triangules area formed by P0 and iEdges nodes points (acumulated in aa variable) 
    !and Pt and iEdges points (acumulated in aa1 variable).
    Do i=1,i34
        !Get the nodes for each edge
        n1=MeshParam%Quadri(MeshParam%EdgeDef(1,i),nel) + 1
        n2=MeshParam%Quadri(MeshParam%EdgeDef(2,i),nel) + 1
        
        aa=aa+dabs(signa(MeshParam%xNode(n1),MeshParam%xNode(n2),x0,MeshParam%yNode(n1),MeshParam%yNode(n2),y0))
        aa1=aa1+dabs(signa(MeshParam%xNode(n1),MeshParam%xNode(n2),xt,MeshParam%yNode(n1),MeshParam%yNode(n2),yt))
    EndDo !i
    
    !If P0 is inside on the Element(nnel) iEdge barycenter, then the acumalated area calculate in previous step is equal to %Area(iElement).
    !Check if P0(x0,y0,z0) is enclosed in nel:
    ae=dabs(aa-MeshParam%Area(nel))/MeshParam%Area(nel)
    If(ae.gt.small1) Then
        print*,'(x0,y0) not in nnel initially',ae,nnel,id0
        pause
        stop
    EndIf
    
    !If PT is inside the Element(nnel), then the acumalated area calculate in previous step is lower than %Area(iElement).
    !In this case, nnel was finded and need to find iLayer inside nnel (go to 400)
    ae=dabs(aa1-MeshParam%Area(nel))/MeshParam%Area(nel)
    if(ae.lt.small1) Then
        nnel=nel
        go to 400
    endif
     
    !(xt,yt) not in nel, and thus (x0,y0) and (xt,yt) are distinctive
    !An interior pt close to (x0,y0) to prevent underflow for iloc >=1.
    If(iloc.eq.0) Then
        !xcg= MeshParam%EdgeBary(1,id0) !x0
        !ycg= MeshParam%EdgeBary(2,id0) !y0
        xcg = x0 !MeshParam%EdgeBary(1,id0) !x0
        ycg = y0 !MeshParam%EdgeBary(2,id0) !y0
    Elseif(iloc.eq.1) Then
        !xcg=(1-1.0d-4)*MeshParam%EdgeBary(1,id0)+1.0d-4*MeshParam%xb(nel)
        !ycg=(1-1.0d-4)*MeshParam%EdgeBary(2,id0)+1.0d-4*MeshParam%yb(nel)
        xcg = (1-1.0d-4)*x0+1.0d-4*MeshParam%xb(nel)
        ycg = (1-1.0d-4)*y0+1.0d-4*MeshParam%yb(nel)
    !	|x0-xcg|/|x0-xctr|=1.0d-3
    endif
      
    pathl=dsqrt((xt-xcg)**2+(yt-ycg)**2)
    If(xcg.eq.xt.and.ycg.eq.yt.or.pathl.eq.0) Then
        print*,'Zero path',x0,y0,xt,yt,xcg,ycg
        pause
        stop
    endif
    
    ! Check if the particle position crosses some Edge from Element(nel):
    ! Starting edge nel_j
    Do i=1,i34
        
        !Get the nodes for each edge
        n1=MeshParam%Quadri(MeshParam%EdgeDef(1,i),nel) + 1
        n2=MeshParam%Quadri(MeshParam%EdgeDef(2,i),nel) + 1
        !Check intersection between lines segments, one formed by nodes points (Edge) and another formed by Pt(xt,yt) and Pcg(xcg,ycg):
        call intersect2(xcg,xt,MeshParam%xNode(n1),MeshParam%xNode(n2),ycg,yt,MeshParam%yNode(n1),MeshParam%yNode(n2),iflag,xin,yin,tt1,tt2)
        !If has a intersection point PI(xin, yin), the Edge i is crosses in particle trajectory:
        If(iflag.eq.1) Then
            nel_j=i
            !Find the moment which the particle crosses the edge:
            if(sqrt((x0-xin)**0.5+(y0-yin)**0.5)>small1 .and. abs(sqrt(uuint**0.5+vvint**0.5))>0 )then
                 dtin = sqrt((x0-xin)**0.5+(y0-yin)**0.5)/sqrt(uuint**0.5+vvint**0.5)
            endif
            Go to 399
        Endif
    Enddo !i=1,3
    
    !No one Edge in element (nel) was crossed by particle:
    If(iflag.eq.0) Then
        xt = x0   
        yt = y0
        zt = z0
        return
    EndIf
    
    print*,'Found no intersecting edges'
    print*,'xcg',xcg
    print*,'xt',xt
    print*,'ycg',ycg
    print*,'yt',yt
    print*,'Elem',nel
    pause
    stop

399 continue
    
    zin=z0 !initialize
    it=0
    
    loop4: Do
!----------------------------------------------------------------------------------------
        it=it+1
        If(it.gt.1000) Then
            !write(12,*)'Death trap reached',idt,id0
            nfl=1
            xt=xin
            yt=yin
            zt=zin
            nnel=nel
            Exit loop4
        Endif
        md1 = MeshParam%Quadri(MeshParam%EdgeDef(1,nel_j),nel) + 1 
        md2 = MeshParam%Quadri(MeshParam%EdgeDef(2,nel_j),nel) + 1
        
        !Compute z position 
        dist=dsqrt((xin-xt)**2+(yin-yt)**2)
        if(dist/pathl.gt.1+1.0d-4) Then
            print*,'Path overshot'
            pause
            stop
        endif
        zin=zt-dist/pathl*(zt-zin)
        trm=trm*dist/pathl !time remaining
        
        pathl=dsqrt((xin-xt)**2+(yin-yt)**2)
        If(pathl.eq.0.or.trm.eq.0) Then
            print*,'Target reached'
            pause
            stop
        endif
        
        lit=0 !flag
        !For horizontal exit and dry elements, compute tangential vel.,
        !update target (xt,yt,zt) and continue.
        isd = MeshParam%Edge(nel_j,nel)
        id0 = isd
        r = MeshParam%Right(isd)
            
        if (nel==179) then
            continue
        endif

        !If particle cross iEdge from iElement(nel), so the new iElement is the neighbour Element which share iEdge(nel_j).
        !However, it can a abnormal case which either iEdge(isd) no has neighbour (horizontal exit or wall) or is dry:
        If (r == 0 .or. HydroParam%H(isd)-HydroParam%hj(isd)<=HydroParam%Pcri+NearZero) Then
            lit=1
            !Nudge intersect (xin,yin), and update starting pt
            xin=(1-1.0d-4)*xin+1.0d-4*MeshParam%xb(nel)
            yin=(1-1.0d-4)*yin+1.0d-4*MeshParam%yb(nel)
            xcg=xin
            ycg=yin
            
            !Set tang. velocities:
            xvel= 0.!uxy(jlev,1,isd)
            yvel= 0.!uxy(jlev,2,isd)
            zvel=0.5*((HydroParam%uNode(jlev,3,md1)+HydroParam%uNode(jlev,3,md2))/2. + (HydroParam%uNode(jlev+1,3,md1)+HydroParam%uNode(jlev+1,3,md2))/2.)
            
            !Update Pt:
            xt=xin-xvel*trm
            yt=yin-yvel*trm
            zt=zin-zvel*trm
            
            !Horizontal velocity magnitude:
            hvel=dsqrt(xvel**2+yvel**2)
            If(hvel.lt.1.e-4) Then !Checar essa condi��o, todos os casos entram aqui CAYO
                nfl=1
                xt=xin
                yt=yin
                zt=zin
                nnel=nel
                IntersectFlag = 1
                exit loop4
            EndIf
            pathl=hvel*trm    
            
        ElseIf( MeshParam%Neighbor(nel_j,nnel) .ne. 0 ) Then
                If (HydroParam%eta(MeshParam%Neighbor(nel_j,nnel)) - HydroParam%hb(MeshParam%Neighbor(nel_j,nnel)) <= HydroParam%PCRI + NearZero) Then !CAYO
                    lit=1
                    !Nudge intersect (xin,yin), and update starting pt
                    xin=(1-1.0d-4)*xin+1.0d-4*MeshParam%xb(nel)
                    yin=(1-1.0d-4)*yin+1.0d-4*MeshParam%yb(nel)
                    xcg=xin
                    ycg=yin
            
                    !Set tang. velocities:
                    xvel = 0.!uxy(jlev,1,isd)
                    yvel = 0.!uxy(jlev,2,isd)
                    zvel = 0.5*((HydroParam%uNode(jlev,3,md1)+HydroParam%uNode(jlev,3,md2))/2. + (HydroParam%uNode(jlev+1,3,md1)+HydroParam%uNode(jlev+1,3,md2))/2.)
            
                    !Update Pt:
                    xt=xin-xvel*trm
                    yt=yin-yvel*trm
                    zt=zin-zvel*trm
            
                    !Horizontal velocity magnitude:
                    hvel=dsqrt(xvel**2+yvel**2)
                    If(hvel.lt.1.e-4) Then !Checar essa condi��o, todos os casos entram aqui CAYO
                        nfl=1
                        xt=xin
                        yt=yin
                        zt=zin
                        nnel=nel
                        IntersectFlag = 1
                        exit loop4
                    EndIf
                pathl=hvel*trm
            ElseIf(dmax1(zt,HydroParam%hb(nel0)) < HydroParam%hb(MeshParam%Neighbor(nel_j,nnel)) ) Then
                nnel = nel0 
                xt=(1-1.0d-4)*MeshParam%EdgeBary(1,MeshParam%Edge(nel_j,nel0)) + 1.0d-4*MeshParam%xb(nel0)
                yt=(1-1.0d-4)*MeshParam%EdgeBary(2,MeshParam%Edge(nel_j,nel0)) + 1.0d-4*MeshParam%yb(nel0)
                nfl=1
                IntersectFlag = 1
                exit loop4
            EndIf
        EndIf

        !Else in normal cases, we need get the neighbour element which shares the iEdge(nel_j):
        If(lit.eq.0) Then
            !next front element
            nel = MeshParam%Neighbor(nel_j,nel)
        EndIf
        
        !With the element updated, we check again if Pt(xt,yt) is inside the element:
        aa=0
        Do i=1,i34
            k1=MeshParam%Quadri(MeshParam%EdgeDef(1,i),nel) + 1
            k2=MeshParam%Quadri(MeshParam%EdgeDef(2,i),nel) + 1
            aa=aa+dabs(signa(MeshParam%xNode(k1),MeshParam%xNode(k2),xt,MeshParam%yNode(k1),MeshParam%yNode(k2),yt))
        EndDo !i
        
        !If is inside, the acumulated area (aa) is lower than %Area(iElement):
        ae = dabs(aa-MeshParam%Area(nel))/MeshParam%Area(nel)
        If(ae.lt.small1) Then
            nnel=nel
            ! Element find, xt and yt position defined. Go to find zt positon (Go to 400).
            Exit loop4
        EndIf
        
        !The particle is out the element, find next intersecting edge:
        Do j=1,i34
            
            !Get the nodes for each edge:
            jd1=MeshParam%Quadri(MeshParam%EdgeDef(1,j),nel) + 1 !nm(nel,nx(i34(nel),j,1))
            jd2=MeshParam%Quadri(MeshParam%EdgeDef(2,j),nel) + 1 !nm(nel,nx(i34(nel),j,2))
           
            if(jd1.eq.md1.and.jd2.eq.md2.or.jd2.eq.md1.and.jd1.eq.md2) cycle !iEdge shared with previous element that was checked, skip  to next iEdge.
            
            !Check path intersection:
            call intersect2(xcg,xt,MeshParam%xNode(jd1),MeshParam%xNode(jd2),ycg,yt,MeshParam%yNode(jd1),MeshParam%yNode(jd2),iflag,xin,yin,tt1,tt2)
            
            !If has a intersection point PI(xin, yin), the Edge i is crosses in particle trajectory:
            if(iflag.eq.1) then
                nel_j=j !next front edge     
                
                if(sqrt((x0-xin)**0.5+(y0-yin)**0.5)>small1 .and. abs(sqrt(uuint**0.5+vvint**0.5))>0 )then
                     dtin = sqrt((x0-xin)**0.5+(y0-yin)**0.5)/sqrt(uuint**0.5+vvint**0.5)
                endif                
                
                cycle loop4
            endif
            
        EndDo !j
        !print*,'Failed to find next edge',lit,xin,yin,xt,yt,nel,&
        !&md1,md2,idt
        xt = x0
        yt = y0
        Exit loop4
        pause
        stop
    
    EndDo loop4

400 Continue
    !The element (nnel) was found, now we are looking for the iLayer (jlev) that contains Pt(xt,yt,zt).
    !First zt is set as the maximum value between Ze (smallm)lower layer and zt, this ensure that the particle no set in position
    !under the bounds of vertical grid discretization. In cases which the zt is set as Ze(smallm) layers, this implies that the
    !particle reach to bottom.
    !Next zt is set as the minimum between it and  Ze(CapitalM+1) upper layer, this ensure that particle not set in a position
    !above the bounds of vertical grid discretization. 
    !zt is set as the minimum value between max[zt calculed and Lower Layer Level in the Element (nnel)] and Upper Layer in Element(nnel)
    !We want find the iLayer that includes this value, not specfic point.

    !If (nnel /= nel0) Then
    !    !zt = dmax1(zt,HydroParam%Ze(HydroParam%ElSmallm(nel0),nel0)+sum(HydroParam%DZsi(:,nel0)))
    !    zt = dmax1(zt,HydroParam%hb(nel0))
    !    !If (zt < HydroParam%Ze(HydroParam%ElSmallm(nnel),nnel) - HydroParam%hb(nnel) ) then
    !    If (zt < HydroParam%hb(nnel) ) then
    !        nnel = nel0
    !        xt = MeshParam%EdgeBary(1,MeshParam%Edge(nel_j,nel0))
    !        yt = MeshParam%EdgeBary(2,MeshParam%Edge(nel_j,nel0))
    !    EndIf
    !EndIf   
    !zt = dmin1(dmax1(zt,HydroParam%Ze(HydroParam%ElSmallm(nnel),nnel)+sum(HydroParam%DZsi(:,nnel))),HydroParam%Ze(HydroParam%ElCapitalM(nnel)+1,nnel))
    zt = dmin1(dmax1(zt,HydroParam%Ze(HydroParam%ElSmallm(nnel),nnel)+sum(HydroParam%DZsi(:,nnel))),HydroParam%eta(nnel))
    Do k = HydroParam%ElSmallm(nnel), HydroParam%ElCapitalM(nnel)
        If (zt.gt.HydroParam%Ze(k,nnel).and.zt.le.HydroParam%Ze(k+1,nnel)) Then
           jlev = k
        EndIf
    EndDo
    
    Return
    End Subroutine quicksearch
    
    
    Subroutine ComputeFW(MeshParam,HydroParam,dt)

        ! Compute the Explicit Finite Difference Operator Fw

        ! Based on: 
        ! [1] Wang,B.; Zhao,G.; Fringer,O.B. Reconstruction of vector fields for semi-Lagrangian advection on unstructured, staggered grids.
        !   Ocean Modelling, 40, p. 52-71, 2011.
        ! [2] Zhang, Y.; Baptista, A.M. SELFE: A semi-implicit Eulerian�Lagrangian finite-element model for cross-scale ocean circulation.
        !   Ocean Modelling, 21, p. 71-96, 2008.
        ! [3] Fringer, O.B.; Gerritsen, M.; Street, R.L. An unstructured-grid, finite-volume, nonhydrostatic, parallel coastal ocean simulator.
        !   Ocean Modelling, 14, p. 139-173, 2006.
        ! [4] Cheng, Ralph T.; Casulli, Vincenzo; Gartner, Jeffrey W. Tidal, Residual, Intertidal Mudflat (TRIM) Model and its Applications to San Francisco Bay, California.
        !   Estuarine, Coastal and Shelf Sciences, 36, p. 235-280. 1993.
    
        ! Input:
        ! Edge-Based Velocity Vector
        ! Vertical Velocity Field

        ! Output:
        ! Advected Vertical Velocity Field    
    
        ! List of Modifications:
        !   28.06.2017: Routine Implementation              (J. Rafael Cavalcanti)
    
        ! Note:
        ! -> I'm neglecting the viscosity effects in the Vertical Momentum transport

        ! Programmer: J. Rafael Cavalcanti
    
        Use Hydrodynamic !, Only: nElem, KMax, ElSmallmOld, ElCapitalMOld, hb, nNode, Ze, nEdge
        Use MeshVars
        !Use Hydrodynamic, Only: InCircle, Edge, xb, yb, tri, x, y, Edgelength, Area, Dzin
        !Use Hydrodynamic, Only: eta, etan, w, u, un, uEdge
        !Use Param, Only: Pcri, dt, Theta
        Implicit None
        Integer:: iElem,iLayer,jlev, jjlev,l,r,nnel,iNode,j,f,FuFw_flag, TrajectoryFlag, ndels, iEdge, Face
        Real:: x0,y0,z0,xt,yt,zt,uuint,vvint,wwint,hhint,wdown,wup,vmag,dtb
        Real:: NearZero = 1e-10
        Real:: dt
        type(MeshGridParam) :: MeshParam
        type(HydrodynamicParam) :: HydroParam
        Integer:: id0=1
         
        HydroParam%Fw = HydroParam%w
    
        Do iElem = 1,MeshParam%nElem
        
            If ( HydroParam%eta(iElem) - HydroParam%hb(iElem) < HydroParam%PCRI/2. ) Then
                 HydroParam%Fw(:,iElem) = 0.
                Cycle
            EndIf
            
            do iEdge = 1,4
            Face = MeshParam%Edge(iEdge,iElem)
            if (HydroParam%IndexInflowEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) then
                HydroParam%Fw(:,iElem) = HydroParam%w(:,iElem)
                cycle
            endif
        enddo

            Do iLayer = HydroParam%ElSmallm(iElem)+1, HydroParam%ElCapitalM(iElem)
                 
                !Do iEdge = 1,4
                !    if (MeshParam%Neighbor(iEdge,iElem)/=0) then
                !        If (iLayer<=HydroParam%ElSmallm(MeshParam%Neighbor(iEdge,iElem))) Then
                !            HydroParam%Fw(iLayer,iElem) = HydroParam%w(iLayer,iElem)
                !            continue
                !            return
                !        EndIf
                !    else
                !        !HydroParam%Fw(iLayer,iElem) = HydroParam%w(iLayer,iElem)
                !        !continue
                !        !return
                !    Endif
                !EndDo
                
                nnel = iElem 
                jlev = iLayer
                x0 = MeshParam%xb(iElem)
                y0 = MeshParam%yb(iElem)
                z0 = HydroParam%Ze(iLayer,iElem)
                uuint = HydroParam%uxyL(iLayer,1,iElem)
                vvint = HydroParam%uxyL(iLayer,2,iElem)
                wwint = HydroParam%w(iLayer,iElem)
                vmag = dsqrt(uuint**2+vvint**2+wwint**2)
                If(vmag.le.1.e-6) Then ! No activity
                    HydroParam%Fw(iLayer,iElem) = HydroParam%w(iLayer,iElem)
                Else !do backtrack
                    TrajectoryFlag = 0 ! = 0 -> System Defined Intervals to Integrate Particle Trajectory; = 1 -> User Defined Intervals to Integrate Particle Trajectory
                    ! 3.5 Find the Local Time Step to Integrate Trajectory
                    If ( TrajectoryFlag == 0 ) Then         ! System Calculation - [7,8]
                        dtb   = MinVal(MeshParam%InCircle)/vmag ! ( InCircle(lElem)/Sqrt( Veloc(1)**2. + Veloc(2)**2. ), InCircle(r)/Sqrt( Veloc(1)**2. + Veloc(2)**2. ) )
                        ndels = Max(Floor(dt/dtb),HydroParam%NFUT)
                        dtb = dt/ndels !sub-step in backtracking
                    Else If ( TrajectoryFlag == 1 ) Then       ! User Defined sub-steps numbers
                        ndels = HydroParam%NFUT 
                        dtb = dt/ndels !sub-step in backtracking
                    EndIf
                    FuFw_flag = 1
                    Call btrack(ndels,dtb,dt,uuint,vvint,wwint,hhint,&
                    &x0,y0,z0,xt,yt,zt,nnel,jlev,id0,HydroParam,MeshParam,FuFw_flag, iLayer, iEdge)
                
                    HydroParam%Fw(iLayer,iElem) = wwint
                    
                Endif
            
            EndDo
        EndDo    
        
        Return
    End Subroutine ComputeFW        
    
  
    
    Subroutine intersect2(x1,x2,x3,x4,y1,y2,y3,y4,iflag,xin,yin,tt1,tt2)
    
  !> Subroutine to detect if two segments (1,2) and (3,4) have a common point
  !>\note The 4 pts are distinctive.\n
  !> The eqs. for the 2 lines are: X=X1+(X2-X1)*tt1 and X=X3+(X4-X3)*tt2.\n
  !> Output: iflag: 0: no intersection or colinear; 1: exactly 1 intersection..\n
  !> If iflag=1, (xin,yin) is the intersection.
  !>\param iflag A common point was founded (iflag=0: no; iflag=1: yes)
  !>\param x1 x position of the point P1 (input)
  !>\param x2 x position of the point P2 (input)
  !>\param x3 x position of the point P3 (input)
  !>\param x4 x position of the point P4 (input)
  !>\param y1 y position of the point P1 (input)
  !>\param y2 y position of the point P2 (input)
  !>\param y3 y position of the point P3 (input)
  !>\param y4 y position of the point P4 (input)
  !>\return xin: x position in the intersection (output)
  !>\return yin: y position in the intersection (output)
  !>\return tt1: angular coefficient of the segment formmed by P1 and P2
  !>\return tt2: angular coefficient of the segment formmed by P3 and P4
  !>\author  Carlos Ruberto Fragoso
  !>\attention  List of Modifications: \n
  !> - 15.09.2015: Routine Implementation (Carlos Ruberto Fragoso)
    
    Implicit none
    Real, parameter :: small2=0.0 !small positive number or 0
    Real, intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4
    Integer, intent(out) :: iflag
    Real, intent(out) :: xin,yin,tt1,tt2
    Real:: delta,delta1,delta2

    tt1=-1000
    tt2=-1000
    iflag=0
    delta=(x2-x1)*(y3-y4)-(y2-y1)*(x3-x4)
    delta1=(x3-x1)*(y3-y4)-(y3-y1)*(x3-x4)
    delta2=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

    If(delta.ne.0.0d0) Then
        tt1=delta1/delta
        tt2=delta2/delta
        If(tt1.ge.-small2.and.tt1.le.1+small2.and.tt2.ge.-small2.and.tt2.le.1+small2) Then
            iflag=1
            xin=x1+(x2-x1)*tt1
            yin=y1+(y2-y1)*tt1
        Endif
    Endif

    Return
    End Subroutine intersect2    
    
    Subroutine ibilinear(elem,x1,x2,x3,x4,y1,y2,y3,y4,x,y,xi,eta,shapef)

      !> Inverse bilinear mapping for quadrangles
      !>\note Convexity of the quad must have been checked, and (x,y) must not be outside the quad.
      !>\param elem Number of element (input)
      !>\param x1 x position of the point P1 (input)
      !>\param x2 x position of the point P2 (input)
      !>\param x3 x position of the point P3 (input)
      !>\param x4 x position of the point P4 (input)
      !>\param y1 y position of the point P1 (input)
      !>\param y2 y position of the point P2 (input)
      !>\param y3 y position of the point P3 (input)
      !>\param y4 y position of the point P4 (input)
      !>\param x x position in the interested point (input)
      !>\param y y position in the interested point (input)
      !>\return shapef: weighted coefficients in each vertices (output)
      !>\author  Carlos Ruberto Fragoso
      !>\attention  List of Modifications: \n
      !> - 15.09.2015: Routine Implementation (Carlos Ruberto Fragoso)
    
    Implicit None
    
    Real, parameter:: small1=1e-6
    Real, parameter:: small3=1.e-5

    Real, intent(in) :: elem,x1,x2,x3,x4,y1,y2,y3,y4,x,y !elem for debugging only
    Real, intent(out) :: xi,eta,shapef(4)
    Real:: x0,y0,axi,aet,bxy,root_xi,root_et,dxi,deta,dd,beta,gamma,delta
    Integer:: icaseno,icount,i,j

    dimension axi(2),aet(2),bxy(2),root_xi(2),root_et(2)

    !Consts.
    x0=(x1+x2+x3+x4)/4
    y0=(y1+y2+y3+y4)/4
    axi(1)=x2-x1+x3-x4      
    axi(2)=y2-y1+y3-y4      
    aet(1)=x3+x4-x1-x2
    aet(2)=y3+y4-y1-y2
    bxy(1)=x1-x2+x3-x4
    bxy(2)=y1-y2+y3-y4

    dxi=2*((x3-x4)*(y1-y2)-(y3-y4)*(x1-x2))
    deta=2*((x4-x1)*(y3-y2)-(y4-y1)*(x3-x2))

    !Inverse mapping
    If(dabs(bxy(1))<small3.and.dabs(bxy(2))<small3.or.dabs(dxi)<small3.and.dabs(deta)<small3) Then
        icaseno=1      
        !print*, 'Entering case 1'
        dd=axi(1)*aet(2)-axi(2)*aet(1)
        if(dd==0) then
	        print*,'Case 1 error:',dd
            pause
	        stop
        EndIf
        xi=4*(aet(2)*(x-x0)-aet(1)*(y-y0))/dd
        eta=4*(axi(1)*(y-y0)-axi(2)*(x-x0))/dd

    Elseif(dabs(dxi)<small3.and.dabs(deta)>=small3) Then   
        icaseno=2      
        !print*, 'Entering case 2'
        eta=4*(bxy(2)*(x-x0)-bxy(1)*(y-y0))/deta
        dd=(axi(1)+eta*bxy(1))**2+(axi(2)+eta*bxy(2))**2
        If(dd==0) then
	        print*,'Case 2 error:',dd
            pause
            stop
        EndIf
        xi=((4*(x-x0)-eta*aet(1))*(axi(1)+eta*bxy(1))+(4*(y-y0)-eta*aet(2))*(axi(2)+eta*bxy(2)))/dd

    Elseif(dabs(dxi)>=small3.and.dabs(deta)<small3) Then   
        icaseno=3      
    !	print*, 'Entering case 3'
        xi=4*(bxy(2)*(x-x0)-bxy(1)*(y-y0))/dxi
        dd=(aet(1)+xi*bxy(1))**2+(aet(2)+xi*bxy(2))**2
        If(dd==0) Then
            print*,'Case 3 error:',dd
            pause
            stop
        EndIf
        eta=((4*(x-x0)-xi*axi(1))*(aet(1)+xi*bxy(1))+(4*(y-y0)-xi*axi(2))*(aet(2)+xi*bxy(2)))/dd

    Else !General case
        icaseno=4      
        !print*, 'Entering case 4'
        beta=aet(2)*axi(1)-aet(1)*axi(2)-4*(bxy(2)*(x-x0)-bxy(1)*(y-y0))
        gamma=4*(aet(1)*(y-y0)-aet(2)*(x-x0))
        delta=beta*beta-4*gamma*dxi
        If(delta==0) Then
	        xi=-beta/2/dxi
	        eta=(4*(bxy(2)*(x-x0)-bxy(1)*(y-y0))-xi*dxi)/deta
        Elseif(delta>0) Then
            !print*, 'Entering case 4.2'
	        root_xi(1)=(-beta+dsqrt(delta))/2/dxi
	        root_xi(2)=(-beta-dsqrt(delta))/2/dxi
	        icount=0
	        Do i=1,2
	            root_et(i)=(4*(bxy(2)*(x-x0)-bxy(1)*(y-y0))-root_xi(i)*dxi)/deta
	            If(dabs(root_xi(i))<=1.1.and.dabs(root_et(i))<=1.1) Then
	                xi=root_xi(i)
	                eta=root_et(i)
	                icount=icount+1
	            EndIf
	        EndDo !i
            If(icount==2.and.dabs(root_xi(1)-root_xi(2)).lt.small1) Then
        !	     Do nothing
        !	     xi=root_xi(1)
        !            eta=root_et(1)
	        Elseif(icount/=1) then
	            print*,'Abnormal instances',(root_xi(j),root_et(j),j=1,2),icount,elem
	            print*,x,y,x1,x2,x3,x4,y1,y2,y3,y4
	            print*,dxi,deta,bxy(1),bxy(2)
                pause
	            stop
	        endif

        Else
	        print*,'No roots',delta,elem
            pause
	        stop
        EndIf
    EndIf

    !If(dabs(xi)>1.1.or.dabs(eta)>1.1) Then
    !    !print*,'Out of bound in ibilinear:',xi,eta,elem,icaseno
    !    !print*,x,y
    !    !pause
    !    !stop
    !endif

    xi=dmin1(1.d0,dmax1(xi,-1.d0))
    eta=dmin1(1.d0,dmax1(eta,-1.d0))
    shapef(1)=(1-xi)*(1-eta)/4
    shapef(2)=(1+xi)*(1-eta)/4
    shapef(3)=(1+xi)*(1+eta)/4
    shapef(4)=(1-xi)*(1+eta)/4

    Return
    End Subroutine ibilinear    
    
    Subroutine iQuadratic(uuBtrack, vvBtrack, wwBtrack, uuN, vvN, wwN, xxN, yyN, zzN, xp, yp, zp )
    !( uuNode(3,9), vvNode(3,9), wwNode(3,9), xxNode(3,9), yyNode(3,9), zzNode(3,9) )

    Implicit None

    Real, intent(in) :: uuN(3,9), vvN(3,9), wwN(3,9), xxN(3,9), yyN(3,9), zzN(3,9), xp, yp, zp
    Real, intent(out) :: uuBtrack, vvBtrack, wwBtrack
    Real:: LZ(3,9), LY(3,3), LX(3), Zuu(9), Zvv(9), Zww(9), Yuu(3), Yvv(3), Yww(3), Xuu, Xvv, Xww, P1, P2, resto, Uresto, Vresto, Wresto, soma
    Integer:: m, iNode, iTimes, cont
    Real:: NearZero = 1e-5
    !1. Interpolatting in Z direction 9 times, one for each node. (e.g. Hodges, 2000)
    !1.1. - find the  Lagrange coefficient formula 
    !for vertical interpolations
 
    Do iNode=1,9
        ! Check if Nodes in this vertical Line is dry or wet:
        LZ(:,iNode) = 0.0d0
        If(zzN(1,iNode) < zzN(3,iNode)) Then
            If(zzN(1,iNode)==zzN(2,iNode)) Then
                ! First Order Lagrange polynomial approach:
                LZ(1,iNode) = 0.0d0
                LZ(2,iNode) = ((zp-zzN(3,iNode))/(zzN(2,iNode)-zzN(3,iNode)))
                LZ(3,iNode) = ((zp-zzN(2,iNode))/(zzN(3,iNode)-zzN(2,iNode)))
                soma = LZ(1,iNode)+LZ(2,iNode)+LZ(3,iNode)
                If(soma /= 1.0d0) Then
                    LZ(2,iNode) = LZ(2,iNode)/soma
                    LZ(3,iNode) = LZ(3,iNode)/soma
                EndIf 
            Else
                ! Second Order Lagrange polynomial approach:
                Do m=1,3
                    if (m==1) then
                        P1 = ((zp-zzN(2,iNode))/(zzN(1,iNode)-zzN(2,iNode)))
                        P2 = ((zp-zzN(3,iNode))/(zzN(1,iNode)-zzN(3,iNode)))
                    elseif (m==2) then
                        P1 = ((zp-zzN(1,iNode))/(zzN(2,iNode)-zzN(1,iNode)))
                        P2 = ((zp-zzN(3,iNode))/(zzN(2,iNode)-zzN(3,iNode)))
                    else
                        P1 = ((zp-zzN(1,iNode))/(zzN(m,iNode)-zzN(1,iNode)))
                        P2 = ((zp-zzN(2,iNode))/(zzN(m,iNode)-zzN(2,iNode)))
                    Endif
                    LZ(m,iNode)=P1*P2
                EndDo
                soma = LZ(1,iNode)+LZ(2,iNode)+LZ(3,iNode)
                if (isnan(P1).or.isnan(P2)) then
                    continue
                endif
            EndIf
        Endif
    EndDo
    !1.2. - Interpolating velocties (u, v, w) to the btrack particle cota 
    Do iNode=1,9        
        Zuu(iNode) = LZ(1,iNode)*uuN(1,iNode)+LZ(2,iNode)*uuN(2,iNode)+LZ(3,iNode)*uuN(3,iNode)
        Zvv(iNode) = LZ(1,iNode)*vvN(1,iNode)+LZ(2,iNode)*vvN(2,iNode)+LZ(3,iNode)*vvN(3,iNode)
        Zww(iNode) = LZ(1,iNode)*wwN(1,iNode)+LZ(2,iNode)*wwN(2,iNode)+LZ(3,iNode)*wwN(3,iNode)
    EndDo
    
    !2. Interpolate in Y direction 3 times
    !2.1. - find the  Lagrange coefficient formula
    cont=0
    Do iNode=1,7,3
        cont=cont+1
        Do m=1,3
            if (m==1) then
                P1 = ((yp-yyN(1,iNode+1))/(yyN(1,iNode)-yyN(1,iNode+1)))
                P2 = ((yp-yyN(1,iNode+2))/(yyN(1,iNode)-yyN(1,iNode+2)))
            elseif (m==2) then
                P1 = ((yp-yyN(1,iNode))/(yyN(1,iNode+1)-yyN(1,iNode)))
                P2 = ((yp-yyN(1,iNode+2))/(yyN(1,iNode+1)-yyN(1,iNode+2)))
            else
                P1 = ((yp-yyN(1,iNode))/(yyN(1,iNode+2)-yyN(1,iNode)))
                P2 = ((yp-yyN(1,iNode+1))/(yyN(1,iNode+2)-yyN(1,iNode+1)))
            Endif
            LY(m,cont)=P1*P2
        EndDo
        soma = LY(1,cont)+LY(2,cont)+LY(3,cont)
        if (isnan(P1).or.isnan(P2)) then
            continue
        endif
    EndDo
    !2.2. - Interpolating velocties (u, v, w) to the btrack particle yt
    cont=0
    
    Do iNode=1,7,3
        cont=cont+1
        Yuu(cont) = LY(1,cont)*Zuu(iNode)+LY(2,cont)*Zuu(iNode+1)+LY(3,cont)*Zuu(iNode+2)
        Yvv(cont) = LY(1,cont)*Zvv(iNode)+LY(2,cont)*Zvv(iNode+1)+LY(3,cont)*Zvv(iNode+2)
        Yww(cont) = LY(1,cont)*Zww(iNode)+LY(2,cont)*Zww(iNode+1)+LY(3,cont)*Zww(iNode+2)
    EndDo
    
    !3. Interpolate in X direction one times
    !3.1. - find the  Lagrange coefficient formula
    Do m=1,3
        if (m==1) then
            P1 = ((xp-xxN(1,4))/(xxN(1,1)-xxN(1,4)))
            P2 = ((xp-xxN(1,7))/(xxN(1,1)-xxN(1,7)))
        elseif (m==2) then
            P1 = ((xp-xxN(1,1))/(xxN(1,4)-xxN(1,1)))
            P2 = ((xp-xxN(1,7))/(xxN(1,4)-xxN(1,7)))
        else
            P1 = ((xp-xxN(1,1))/(xxN(1,7)-xxN(1,1)))
            P2 = ((xp-xxN(1,4))/(xxN(1,7)-xxN(1,4)))
        Endif
        LX(m) = P1*P2
    EndDo 
    soma = LX(1)+LX(2)+LX(3)
    if (isnan(P1).or.isnan(P2)) then
            continue
    endif
    !3.2. - Interpolating velocties (u, v, w) to the btrack particle xt
    Xuu = Lx(1)*Yuu(1)+Lx(2)*Yuu(2)+Lx(3)*Yuu(3)
    Xvv = Lx(1)*Yvv(1)+Lx(2)*Yvv(2)+Lx(3)*Yvv(3)
    Xww = Lx(1)*Yww(1)+Lx(2)*Yww(2)+Lx(3)*Yww(3)
    
    !4. Setting the Btrack velocities
    uuBtrack = Xuu
    vvBtrack = Xvv
    wwBtrack = Xww

    Return
    End Subroutine iQuadratic   
    
    Subroutine iBilinear2(uuBtrack, vvBtrack, wwBtrack, uuN, vvN, wwN, xxN, yyN, zzN, xp, yp, zp, x0, y0, z0, id0, nnel, FuFw_flag, MeshParam)
    !( uuNode(3,9), vvNode(3,9), wwNode(3,9), xxNode(3,9), yyNode(3,9), zzNode(3,9) )
    
    Use MeshVars
    
    Implicit None

    Real, intent(in) :: uuN(3,9), vvN(3,9), wwN(3,9), xxN(3,9), yyN(3,9), zzN(3,9), xp, yp, zp, x0, y0, z0
    Real, intent(out) :: uuBtrack, vvBtrack, wwBtrack
    Real::  a, b, d, dx, dy, dz, soma,xr,yr,zr, Nref
    Integer:: m, iNode, iTimes, cont, nnel, id0, FuFw_flag, F1,F2,F3,F4
    type(MeshGridParam) :: MeshParam

    !1. Interpolatting in Z direction 9 times, one for each node. (e.g. Hodges, 2000)
    !1.1. - find the  Lagrange coefficient formula 
    !for vertical interpolations

	dx = abs(xp-x0)
	dy = abs(yp-y0)
	dz = abs(zp-z0)
    
    if (FuFw_flag==0) then
        F1=MeshParam%Edge(1,nnel)
	    F2=MeshParam%Edge(2,nnel)
	    F3=MeshParam%Edge(3,nnel)
	    F4=MeshParam%Edge(4,nnel)
        if (id0==F2) then
			if (yp>=y0) then !North side
				a =  dx/abs(xxN(1,1)-xxN(1,2))
				b =  dy/abs(yyN(1,1)-yyN(1,4))
                d =  dz/abs(zzN(1,1)-zzN(2,1))
				if(abs(zzN(1,1)-zzN(2,1)) == 0.d0) then
                    d = 0.d0
                endif
				uuBtrack = (1-d)*((1-a)*((1-b)*uuN(1,1)+b*uuN(1,4))+a*((1-b)*uuN(1,2)+b*uuN(1,3)))+d*((1-a)*((1-b)*uuN(2,1)+b*uuN(2,4))+a*((1-b)*uuN(2,2)+b*uuN(2,3)))
				vvBtrack = (1-d)*((1-a)*((1-b)*vvN(1,1)+b*vvN(1,4))+a*((1-b)*vvN(1,2)+b*vvN(1,3)))+d*((1-a)*((1-b)*vvN(2,1)+b*vvN(2,4))+a*((1-b)*vvN(2,2)+b*vvN(2,3)))
				wwBtrack = (1-d)*((1-a)*((1-b)*wwN(1,1)+b*wwN(1,4))+a*((1-b)*wwN(1,2)+b*wwN(1,3)))+d*((1-a)*((1-b)*wwN(2,1)+b*wwN(2,4))+a*((1-b)*wwN(2,2)+b*wwN(2,3)))
			else !south side
				a =  dx/abs(xxN(1,1)-xxN(1,4))
				b =  dy/abs(yyN(1,1)-yyN(1,2))
				d =  dz/abs(zzN(1,1)-zzN(2,1))
				if(abs(zzN(1,1)-zzN(2,1)) == 0.d0) then
                    d = 0.d0
                endif
				uuBtrack = (1-d)*((1-a)*((1-b)*uuN(1,1)+b*uuN(1,2))+a*((1-b)*uuN(1,4)+b*uuN(1,3)))+d*((1-a)*((1-b)*uuN(2,1)+b*uuN(2,2))+a*((1-b)*uuN(2,4)+b*uuN(2,3)))
				vvBtrack = (1-d)*((1-a)*((1-b)*vvN(1,1)+b*vvN(1,2))+a*((1-b)*vvN(1,4)+b*vvN(1,3)))+d*((1-a)*((1-b)*vvN(2,1)+b*vvN(2,2))+a*((1-b)*vvN(2,4)+b*vvN(2,3)))
				wwBtrack = (1-d)*((1-a)*((1-b)*wwN(1,1)+b*wwN(1,2))+a*((1-b)*wwN(1,4)+b*wwN(1,3)))+d*((1-a)*((1-b)*wwN(2,1)+b*wwN(2,2))+a*((1-b)*wwN(2,4)+b*wwN(2,3)))
			endif
		elseif (id0==F4) then
			if (yp>=y0) then !North side
				a =  dx/abs(xxN(1,1)-xxN(1,4))
				b =  dy/abs(yyN(1,1)-yyN(1,2))
				d =  dz/abs(zzN(1,1)-zzN(2,1))
				if(abs(zzN(1,1)-zzN(2,1)) == 0.d0) then
                    d = 0.d0
                endif                
				uuBtrack = (1-d)*((1-a)*((1-b)*uuN(1,1)+b*uuN(1,2))+a*((1-b)*uuN(1,4)+b*uuN(1,3)))+d*((1-a)*((1-b)*uuN(2,1)+b*uuN(2,2))+a*((1-b)*uuN(2,4)+b*uuN(2,3)))
				vvBtrack = (1-d)*((1-a)*((1-b)*vvN(1,1)+b*vvN(1,2))+a*((1-b)*vvN(1,4)+b*vvN(1,3)))+d*((1-a)*((1-b)*vvN(2,1)+b*vvN(2,2))+a*((1-b)*vvN(2,4)+b*vvN(2,3)))
				wwBtrack = (1-d)*((1-a)*((1-b)*wwN(1,1)+b*wwN(1,2))+a*((1-b)*wwN(1,4)+b*wwN(1,3)))+d*((1-a)*((1-b)*wwN(2,1)+b*wwN(2,2))+a*((1-b)*wwN(2,4)+b*wwN(2,3)))
			else !south side
				a =  dx/abs(xxN(1,1)-xxN(1,2))
				b =  dy/abs(yyN(1,1)-yyN(1,4))
				d =  dz/abs(zzN(1,1)-zzN(2,1))
				if(abs(zzN(1,1)-zzN(2,1)) == 0.d0) then
                    d = 0.d0
                endif                
				uuBtrack = (1-d)*((1-a)*((1-b)*uuN(1,1)+b*uuN(1,4))+a*((1-b)*uuN(1,2)+b*uuN(1,3)))+d*((1-a)*((1-b)*uuN(2,1)+b*uuN(2,4))+a*((1-b)*uuN(2,2)+b*uuN(2,3)))
				vvBtrack = (1-d)*((1-a)*((1-b)*vvN(1,1)+b*vvN(1,4))+a*((1-b)*vvN(1,2)+b*vvN(1,3)))+d*((1-a)*((1-b)*vvN(2,1)+b*vvN(2,4))+a*((1-b)*vvN(2,2)+b*vvN(2,3)))
				wwBtrack = (1-d)*((1-a)*((1-b)*wwN(1,1)+b*wwN(1,4))+a*((1-b)*wwN(1,2)+b*wwN(1,3)))+d*((1-a)*((1-b)*wwN(2,1)+b*wwN(2,4))+a*((1-b)*wwN(2,2)+b*wwN(2,3)))
			endif
		elseif (id0==F3) then
			if (xp>=x0) then !Right side
				a =  dx/abs(xxN(1,1)-xxN(1,2))
				b =  dy/abs(yyN(1,1)-yyN(1,4))
				d =  dz/abs(zzN(1,1)-zzN(2,1))
				if(abs(zzN(1,1)-zzN(2,1)) == 0.d0) then
                    d = 0.d0
                endif                                
				uuBtrack = (1-d)*((1-a)*((1-b)*uuN(1,1)+b*uuN(1,4))+a*((1-b)*uuN(1,2)+b*uuN(1,3)))+d*((1-a)*((1-b)*uuN(2,1)+b*uuN(2,4))+a*((1-b)*uuN(2,2)+b*uuN(2,3)))
				vvBtrack = (1-d)*((1-a)*((1-b)*vvN(1,1)+b*vvN(1,4))+a*((1-b)*vvN(1,2)+b*vvN(1,3)))+d*((1-a)*((1-b)*vvN(2,1)+b*vvN(2,4))+a*((1-b)*vvN(2,2)+b*vvN(2,3)))
				wwBtrack = (1-d)*((1-a)*((1-b)*wwN(1,1)+b*wwN(1,4))+a*((1-b)*wwN(1,2)+b*wwN(1,3)))+d*((1-a)*((1-b)*wwN(2,1)+b*wwN(2,4))+a*((1-b)*wwN(2,2)+b*wwN(2,3)))
			else !Left side
				a =  dx/abs(xxN(1,1)-xxN(1,4))
				b =  dy/abs(yyN(1,1)-yyN(1,2))
				d =  dz/abs(zzN(1,1)-zzN(2,1))
				if(abs(zzN(1,1)-zzN(2,1)) == 0.d0) then
                    d = 0.d0
                endif                
				uuBtrack = (1-d)*((1-a)*((1-b)*uuN(1,1)+b*uuN(1,2))+a*((1-b)*uuN(1,4)+b*uuN(1,3)))+d*((1-a)*((1-b)*uuN(2,1)+b*uuN(2,2))+a*((1-b)*uuN(2,4)+b*uuN(2,3)))
				vvBtrack = (1-d)*((1-a)*((1-b)*vvN(1,1)+b*vvN(1,2))+a*((1-b)*vvN(1,4)+b*vvN(1,3)))+d*((1-a)*((1-b)*vvN(2,1)+b*vvN(2,2))+a*((1-b)*vvN(2,4)+b*vvN(2,3)))
				wwBtrack = (1-d)*((1-a)*((1-b)*wwN(1,1)+b*wwN(1,2))+a*((1-b)*wwN(1,4)+b*wwN(1,3)))+d*((1-a)*((1-b)*wwN(2,1)+b*wwN(2,2))+a*((1-b)*wwN(2,4)+b*wwN(2,3)))
			endif
		elseif (id0==F1) then
			if (xp>=x0) then !Right side
				a =  dx/abs(xxN(1,1)-xxN(1,4))
				b =  dy/abs(yyN(1,1)-yyN(1,2))
				d =  dz/abs(zzN(1,1)-zzN(2,1))
				if(abs(zzN(1,1)-zzN(2,1)) == 0.d0) then
                    d = 0.d0
                endif                
				uuBtrack = (1-d)*((1-a)*((1-b)*uuN(1,1)+b*uuN(1,2))+a*((1-b)*uuN(1,4)+b*uuN(1,3)))+d*((1-a)*((1-b)*uuN(2,1)+b*uuN(2,2))+a*((1-b)*uuN(2,4)+b*uuN(2,3)))
				vvBtrack = (1-d)*((1-a)*((1-b)*vvN(1,1)+b*vvN(1,2))+a*((1-b)*vvN(1,4)+b*vvN(1,3)))+d*((1-a)*((1-b)*vvN(2,1)+b*vvN(2,2))+a*((1-b)*vvN(2,4)+b*vvN(2,3)))
				wwBtrack = (1-d)*((1-a)*((1-b)*wwN(1,1)+b*wwN(1,2))+a*((1-b)*wwN(1,4)+b*wwN(1,3)))+d*((1-a)*((1-b)*wwN(2,1)+b*wwN(2,2))+a*((1-b)*wwN(2,4)+b*wwN(2,3)))
			else !Left side
				a =  dx/abs(xxN(1,1)-xxN(1,2))
				b =  dy/abs(yyN(1,1)-yyN(1,4))
				d =  dz/abs(zzN(1,1)-zzN(2,1))
				if(abs(zzN(1,1)-zzN(2,1)) == 0.d0) then
                    d = 0.d0
                endif                
				uuBtrack = (1-d)*((1-a)*((1-b)*uuN(1,1)+b*uuN(1,4))+a*((1-b)*uuN(1,2)+b*uuN(1,3)))+d*((1-a)*((1-b)*uuN(2,1)+b*uuN(2,4))+a*((1-b)*uuN(2,2)+b*uuN(2,3)))
				vvBtrack = (1-d)*((1-a)*((1-b)*vvN(1,1)+b*vvN(1,4))+a*((1-b)*vvN(1,2)+b*vvN(1,3)))+d*((1-a)*((1-b)*vvN(2,1)+b*vvN(2,4))+a*((1-b)*vvN(2,2)+b*vvN(2,3)))
				wwBtrack = (1-d)*((1-a)*((1-b)*wwN(1,1)+b*wwN(1,4))+a*((1-b)*wwN(1,2)+b*wwN(1,3)))+d*((1-a)*((1-b)*wwN(2,1)+b*wwN(2,4))+a*((1-b)*wwN(2,2)+b*wwN(2,3)))
			endif
		endif
    else
        if (xp<=x0.and.yp<y0) then !left south side
            a =  dx/abs(xxN(1,1)-xxN(1,2))
	        b =  dy/abs(yyN(1,1)-yyN(1,4))
	        d =  dz/abs(zzN(1,1)-zzN(2,1))
		    if(abs(zzN(1,1)-zzN(2,1)) == 0.d0) then
                d = 0.d0
            endif        
		    uuBtrack = (1-d)*((1-a)*((1-b)*uuN(1,1)+b*uuN(1,4))+a*((1-b)*uuN(1,2)+b*uuN(1,3)))+d*((1-a)*((1-b)*uuN(2,1)+b*uuN(2,4))+a*((1-b)*uuN(2,2)+b*uuN(2,3)))
		    vvBtrack = (1-d)*((1-a)*((1-b)*vvN(1,1)+b*vvN(1,4))+a*((1-b)*vvN(1,2)+b*vvN(1,3)))+d*((1-a)*((1-b)*vvN(2,1)+b*vvN(2,4))+a*((1-b)*vvN(2,2)+b*vvN(2,3)))
		    wwBtrack = (1-d)*((1-a)*((1-b)*wwN(1,1)+b*wwN(1,4))+a*((1-b)*wwN(1,2)+b*wwN(1,3)))+d*((1-a)*((1-b)*wwN(2,1)+b*wwN(2,4))+a*((1-b)*wwN(2,2)+b*wwN(2,3)))
        elseif (xp<=x0.and.yp>=y0) then !left North side
            a =  dx/abs(xxN(1,1)-xxN(1,4))
	        b =  dy/abs(yyN(1,1)-yyN(1,2))
	        d =  dz/abs(zzN(1,1)-zzN(2,1))
	        if(abs(zzN(1,1)-zzN(2,1)) == 0.d0) then
                d = 0.d0
            endif        
		    uuBtrack = (1-d)*((1-a)*((1-b)*uuN(1,1)+b*uuN(1,2))+a*((1-b)*uuN(1,4)+b*uuN(1,3)))+d*((1-a)*((1-b)*uuN(2,1)+b*uuN(2,2))+a*((1-b)*uuN(2,4)+b*uuN(2,3)))
		    vvBtrack = (1-d)*((1-a)*((1-b)*vvN(1,1)+b*vvN(1,2))+a*((1-b)*vvN(1,4)+b*vvN(1,3)))+d*((1-a)*((1-b)*vvN(2,1)+b*vvN(2,2))+a*((1-b)*vvN(2,4)+b*vvN(2,3)))
		    wwBtrack = (1-d)*((1-a)*((1-b)*wwN(1,1)+b*wwN(1,2))+a*((1-b)*wwN(1,4)+b*wwN(1,3)))+d*((1-a)*((1-b)*wwN(2,1)+b*wwN(2,2))+a*((1-b)*wwN(2,4)+b*wwN(2,3)))
        elseif (xp>x0.and.yp<y0) then !Right south side
            a =  dx/abs(xxN(1,1)-xxN(1,4))
	        b =  dy/abs(yyN(1,1)-yyN(1,2))
	        d =  dz/abs(zzN(1,1)-zzN(2,1))
		    if(abs(zzN(1,1)-zzN(2,1)) == 0.d0) then
                d = 0.d0
            endif        
		    uuBtrack = (1-d)*((1-a)*((1-b)*uuN(1,1)+b*uuN(1,2))+a*((1-b)*uuN(1,4)+b*uuN(1,3)))+d*((1-a)*((1-b)*uuN(2,1)+b*uuN(2,2))+a*((1-b)*uuN(2,4)+b*uuN(2,3)))
		    vvBtrack = (1-d)*((1-a)*((1-b)*vvN(1,1)+b*vvN(1,2))+a*((1-b)*vvN(1,4)+b*vvN(1,3)))+d*((1-a)*((1-b)*vvN(2,1)+b*vvN(2,2))+a*((1-b)*vvN(2,4)+b*vvN(2,3)))
		    wwBtrack = (1-d)*((1-a)*((1-b)*wwN(1,1)+b*wwN(1,2))+a*((1-b)*wwN(1,4)+b*wwN(1,3)))+d*((1-a)*((1-b)*wwN(2,1)+b*wwN(2,2))+a*((1-b)*wwN(2,4)+b*wwN(2,3)))
        elseif (xp>x0.and.yp>=y0) then !Right North side
            a =  dx/abs(xxN(1,1)-xxN(1,2))
	        b =  dy/abs(yyN(1,1)-yyN(1,4))
	        d =  dz/abs(zzN(1,1)-zzN(2,1))
		    if(abs(zzN(1,1)-zzN(2,1)) == 0.d0) then
                d = 0.d0
            endif        
		    uuBtrack = (1-d)*((1-a)*((1-b)*uuN(1,1)+b*uuN(1,4))+a*((1-b)*uuN(1,2)+b*uuN(1,3)))+d*((1-a)*((1-b)*uuN(2,1)+b*uuN(2,4))+a*((1-b)*uuN(2,2)+b*uuN(2,3)))
		    vvBtrack = (1-d)*((1-a)*((1-b)*vvN(1,1)+b*vvN(1,4))+a*((1-b)*vvN(1,2)+b*vvN(1,3)))+d*((1-a)*((1-b)*vvN(2,1)+b*vvN(2,4))+a*((1-b)*vvN(2,2)+b*vvN(2,3)))
		    wwBtrack = (1-d)*((1-a)*((1-b)*wwN(1,1)+b*wwN(1,4))+a*((1-b)*wwN(1,2)+b*wwN(1,3)))+d*((1-a)*((1-b)*wwN(2,1)+b*wwN(2,4))+a*((1-b)*wwN(2,2)+b*wwN(2,3)))
        endif
    endif
    
    Return
    End Subroutine iBilinear2 
    
  !> Compute signed area formed by vertices 1,2,3    
    Function signa(x1,x2,x3,y1,y2,y3)

    Implicit none
    Real:: x1,x2,x3,y1,y2,y3,signa

    signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
      
    
    End Function signa     
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
 Subroutine FuVelocities3(uuNode, vvNode, wwNode, xxNode, yyNode, zzNode, nnel,jlev, xt,yt,zt,x0,y0,z0, id0, HydroParam,MeshParam)
    
    Use MeshVars !, Only: Quadri,EdgeDef,xNode,yNode,Area,Edge,EdgeNodes,EdgeBary,Left,Right
    Use Hydrodynamic !, Only: Smallm,ElSmallm,ElCapitalM,CapitalM,Z,DZj,Ze,DZi,uNode,uxy,Ze,DZi,w
    
    Implicit none
    
    integer, intent(in) :: nnel,jlev
    Real, intent(in) ::  xt,yt,zt,x0,y0,z0
    Real, intent(out) ::  uuNode(3,9), vvNode(3,9), wwNode(3,9), xxNode(3,9), yyNode(3,9), zzNode(3,9)
    Integer:: idt,iflqs1,i,j,l,nd,nn,lev,n1,n2,n3,n4, iEdge,iNode, Nodes(4), F1, F2, F3, F4, id0, NodeFlag
    Real:: NearZero = 1e-10 !< Small Number
    Integer:: i34 = 4
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    !Integer:: ELM_flag = 1
    
    F1=MeshParam%Edge(1,nnel)
	F2=MeshParam%Edge(2,nnel)
	F3=MeshParam%Edge(3,nnel)
	F4=MeshParam%Edge(4,nnel)
    !if (jlev==HydroParam%ElSmallm(nnel)) then and.jlev/=HydroParam%ElSmallm(nnel).and.jlev/=HydroParam%ElCapitalM(nnel)
    !MeshParam%Neighbor(1,nnel)
	if (id0==F1) then
        if (xt>=x0) then !Right side
            n1=F1
            n2=F3
            n3=MeshParam%Quadri(4,nnel) + 1
            n4=MeshParam%Quadri(1,nnel) + 1
			NodeFlag = 1
        else !left side
            n1=F1
            n2=MeshParam%Quadri(2,nnel) + 1
            n3=MeshParam%Quadri(3,nnel) + 1
            n4=F3
			NodeFlag = 2
        endif
        Nodes(1:4)= (/ n1, n2, n3, n4 /)
    elseif (id0==F2) then
        if (yt>=y0) then !North side
            n1=F2
            n2=F4
            n3=MeshParam%Quadri(1,nnel) + 1
            n4=MeshParam%Quadri(2,nnel) + 1
			NodeFlag = 1
        else !South side
            n1=F2
            n2=MeshParam%Quadri(3,nnel) + 1
            n3=MeshParam%Quadri(4,nnel) + 1
            n4=F4
			NodeFlag = 2
        endif
        Nodes(1:4)= (/ n1, n2, n3, n4 /)
    elseif (id0==F3) then
        if (xt>=x0) then !Right side
            n1=F3
            n2=MeshParam%Quadri(4,nnel) + 1
            n3=MeshParam%Quadri(1,nnel) + 1
            n4=F1
			NodeFlag = 2
        else !left side
            n1=F3
            n2=F1
            n3=MeshParam%Quadri(2,nnel) + 1
            n4=MeshParam%Quadri(3,nnel) + 1
			NodeFlag = 1
        endif
        Nodes(1:4)= (/ n1, n2, n3, n4 /)
	else !(id0==F4)
        if (yt>=y0) then !North side
            n1=F4
            n2=MeshParam%Quadri(1,nnel) + 1
            n3=MeshParam%Quadri(2,nnel) + 1
            n4=F2
			NodeFlag = 2
        else !South side
            n1=F4
            n2=F2
            n3=MeshParam%Quadri(3,nnel) + 1
            n4=MeshParam%Quadri(4,nnel) + 1	
			NodeFlag = 1
        endif
        Nodes(1:4)= (/ n1, n2, n3, n4 /)
    endif
	 if (NodeFlag == 1) then !n1 and n2 = Face; n3 and n4 = Nodes;
        If (HydroParam%ElSmallm(nnel)==HydroParam%ElCapitalM(nnel)) Then !2D/ 1 Layer
            lev=jlev
                
			uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			uuNode(1,2)=HydroParam%uxy(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%uxy(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%wfc(jlev,Nodes(2))
			uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			uuNode(1,4)=HydroParam%ubv(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%ubv(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%ubv(jlev,3,Nodes(4))
			if (zt<z0) then !Down
                lev=jlev
				uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				uuNode(2,2)=HydroParam%ug(Nodes(2),lev);        vvNode(2,2)=HydroParam%vg(Nodes(2),lev);        wwNode(2,2)=HydroParam%wg(Nodes(2),lev)
				uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				uuNode(2,4)=HydroParam%uNode(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uNode(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%uNode(lev,3,Nodes(4))
            else !up
                lev=jlev+1
				uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				uuNode(2,2)=HydroParam%ug(Nodes(2),lev);        vvNode(2,2)=HydroParam%vg(Nodes(2),lev);        wwNode(2,2)=HydroParam%wg(Nodes(2),lev)
				uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				uuNode(2,4)=HydroParam%uNode(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uNode(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%uNode(lev,3,Nodes(4))
                lev=jlev
			endif
			xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(1,1)=HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			xxNode(1,3) = MeshParam%xNode(Nodes(3)); yyNode(1,3) = MeshParam%yNode(Nodes(3)); zzNode(1,3)=HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			xxNode(1,4) = MeshParam%xNode(Nodes(4)); yyNode(1,4) = MeshParam%yNode(Nodes(4)); zzNode(1,4)=HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			if (zt<z0) then !Down
                lev=jlev
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel)
				xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel)
				xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel)
				xxNode(2,4) = MeshParam%xNode(Nodes(4)); yyNode(2,4) = MeshParam%yNode(Nodes(4)); zzNode(2,4)= HydroParam%Ze(lev,nnel)
            else !Up
                lev=jlev+1
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=HydroParam%Z(lev,Nodes(1))
				xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=HydroParam%Z(lev,Nodes(2))
				xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%peta(Nodes(3))
				xxNode(2,4) = MeshParam%xNode(Nodes(4)); yyNode(2,4) = MeshParam%yNode(Nodes(4)); zzNode(2,4)= HydroParam%peta(Nodes(4))
                lev=jlev
			endif
        else
		    if (jlev==HydroParam%ElCapitalM(nnel)) then

			    if (zt>=z0) then !up
				    lev=jlev+1
			    else !(zt<z0) down
				    lev=jlev-1
			    endif
			    uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			    uuNode(1,2)=HydroParam%uxy(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%uxy(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%wfc(jlev,Nodes(2))
			    uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			    uuNode(1,4)=HydroParam%ubv(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%ubv(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%ubv(jlev,3,Nodes(4))
			    if (zt>=z0) then !up
				    uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				    uuNode(2,2)=HydroParam%ug(Nodes(2),lev);        vvNode(2,2)=HydroParam%vg(Nodes(2),lev);        wwNode(2,2)=HydroParam%wg(Nodes(2),lev)
				    uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				    uuNode(2,4)=HydroParam%uNode(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uNode(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%uNode(lev,3,Nodes(4))
			    else !down
				    uuNode(2,1)=HydroParam%uxy(lev,1,Nodes(1));   vvNode(2,1)=HydroParam%uxy(lev,2,Nodes(1));   wwNode(2,1)=HydroParam%wfc(lev,Nodes(1))
				    uuNode(2,2)=HydroParam%uxy(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%uxy(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%wfc(lev,Nodes(2))
				    uuNode(2,3)=HydroParam%ubv(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%ubv(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%ubv(lev,3,Nodes(3))
				    uuNode(2,4)=HydroParam%ubv(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%ubv(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%ubv(lev,3,Nodes(4))
			    endif
			    xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(1,1)=0.5d0*(HydroParam%Z(jlev,Nodes(1)) + HydroParam%hj(Nodes(1)) + HydroParam%Z(jlev+1,Nodes(1)))
			    xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2));   yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2));   zzNode(1,2)=0.5d0*(HydroParam%Z(jlev,Nodes(2)) + HydroParam%hj(Nodes(2)) + HydroParam%Z(jlev+1,Nodes(2)))
			    xxNode(1,3) = MeshParam%xNode(Nodes(3));        yyNode(1,3) = MeshParam%yNode(Nodes(3));        zzNode(1,3)=0.5d0*(HydroParam%peta(Nodes(3)) + HydroParam%hb(nnel) + HydroParam%Ze(jlev,nnel))
			    xxNode(1,4) = MeshParam%xNode(Nodes(4));        yyNode(1,4) = MeshParam%yNode(Nodes(4));        zzNode(1,4)=0.5d0*(HydroParam%peta(Nodes(4)) + HydroParam%hb(nnel) + HydroParam%Ze(jlev,nnel))
			    if (zt>=z0) then !up
				    xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = HydroParam%Z(lev,Nodes(1))
				    xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2));   yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2));   zzNode(2,2) = HydroParam%Z(lev,Nodes(2))
				    xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = HydroParam%peta(Nodes(3))
				    xxNode(2,4) = MeshParam%xNode(Nodes(4));        yyNode(2,4) = MeshParam%yNode(Nodes(4));        zzNode(2,4) = HydroParam%peta(Nodes(4))
                else !down
                    xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1)= HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
                    xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2));   yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2));   zzNode(2,2)= HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
                    xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3)= HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
                    xxNode(2,4) = MeshParam%xNode(Nodes(4));        yyNode(2,4) = MeshParam%yNode(Nodes(4));        zzNode(2,4)= HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
                endif
           
		    ElseIf (jlev==HydroParam%ElSmallm(nnel)) Then
			    If (zt>=z0) Then !up
				    lev=jlev+1
			    Else !(zt<z0) down
				    lev=jlev
			    EndIf
			    uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			    uuNode(1,2)=HydroParam%uxy(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%uxy(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%wfc(jlev,Nodes(2))
			    uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			    uuNode(1,4)=HydroParam%ubv(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%ubv(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%ubv(jlev,3,Nodes(4))
			    if (zt<z0) then !Down
                
				    uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				    uuNode(2,2)=HydroParam%ug(Nodes(2),lev);        vvNode(2,2)=HydroParam%vg(Nodes(2),lev);        wwNode(2,2)=HydroParam%wg(Nodes(2),lev)
				    uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				    uuNode(2,4)=HydroParam%uNode(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uNode(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%uNode(lev,3,Nodes(4))
                Else !up
				    uuNode(2,1)=HydroParam%uxy(lev,1,Nodes(1));   vvNode(2,1)=HydroParam%uxy(lev,2,Nodes(1));   wwNode(2,1)=HydroParam%wfc(lev,Nodes(1))
				    uuNode(2,2)=HydroParam%uxy(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%uxy(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%wfc(lev,Nodes(2))
				    uuNode(2,3)=HydroParam%ubv(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%ubv(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%ubv(lev,3,Nodes(3))
				    uuNode(2,4)=HydroParam%ubv(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%ubv(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%ubv(lev,3,Nodes(4))
                
			    Endif
			    xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(1,1) = HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			    xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2));   yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2));   zzNode(1,2) = HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			    xxNode(1,3) = MeshParam%xNode(Nodes(3));        yyNode(1,3) = MeshParam%yNode(Nodes(3));        zzNode(1,3) = HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			    xxNode(1,4) = MeshParam%xNode(Nodes(4));        yyNode(1,4) = MeshParam%yNode(Nodes(4));        zzNode(1,4) = HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			    If (zt<z0) Then !Down
				    xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel)
				    xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2));   yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2));   zzNode(2,2) = HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel)
				    xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel)
				    xxNode(2,4) = MeshParam%xNode(Nodes(4));        yyNode(2,4) = MeshParam%yNode(Nodes(4));        zzNode(2,4) = HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel)
                Else !Up
                    If (lev == HydroParam%ElCapitalM(nnel)) Then
                        xxNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = 0.5d0*(HydroParam%Z(lev,Nodes(1)) + HydroParam%hj(Nodes(1)) + HydroParam%Z(lev+1,Nodes(1)))
			            xxNode(2,2) = MeshParam%EdgeBary(2,Nodes(2));   yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2));   zzNode(2,2) = 0.5d0*(HydroParam%Z(lev,Nodes(2)) + HydroParam%hj(Nodes(2)) + HydroParam%Z(lev+1,Nodes(2)))
			            xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = 0.5d0*(HydroParam%peta(Nodes(3)) + HydroParam%hb(nnel) + HydroParam%Ze(lev,nnel))
			            xxNode(2,4) = MeshParam%xNode(Nodes(4));        yyNode(2,4) = MeshParam%yNode(Nodes(4));        zzNode(2,4) = 0.5d0*(HydroParam%peta(Nodes(4)) + HydroParam%hb(nnel) + HydroParam%Ze(lev,nnel))
                    Else
				        xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
				        xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2));   yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2));   zzNode(2,2) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
				        xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
				        xxNode(2,4) = MeshParam%xNode(Nodes(4));        yyNode(2,4) = MeshParam%yNode(Nodes(4));        zzNode(2,4) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
			        EndIf
                EndIf
		    Else 
			    If (zt>=z0) Then !up
				    lev=jlev+1
			    Else !(zt<z0) down
				    lev=jlev-1
			    EndIf
			    uuNode(1,1) = HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1) = HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1) = HydroParam%wfc(jlev,Nodes(1))
			    uuNode(1,2) = HydroParam%uxy(jlev,1,Nodes(2));   vvNode(1,2) = HydroParam%uxy(jlev,2,Nodes(2));   wwNode(1,2) = HydroParam%wfc(jlev,Nodes(2))
			    uuNode(1,3) = HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3) = HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3) = HydroParam%ubv(jlev,3,Nodes(3))
			    uuNode(1,4) = HydroParam%ubv(jlev,1,Nodes(4));   vvNode(1,4) = HydroParam%ubv(jlev,2,Nodes(4));   wwNode(1,4) = HydroParam%ubv(jlev,3,Nodes(4))

			    uuNode(2,1) = HydroParam%uxy(lev,1,Nodes(1));   vvNode(2,1) = HydroParam%uxy(lev,2,Nodes(1));   wwNode(2,1) = HydroParam%wfc(lev,Nodes(1))
			    uuNode(2,2) = HydroParam%uxy(lev,1,Nodes(2));   vvNode(2,2) = HydroParam%uxy(lev,2,Nodes(2));   wwNode(2,2) = HydroParam%wfc(lev,Nodes(2))
			    uuNode(2,3) = HydroParam%ubv(lev,1,Nodes(3));   vvNode(2,3) = HydroParam%ubv(lev,2,Nodes(3));   wwNode(2,3) = HydroParam%ubv(lev,3,Nodes(3))
			    uuNode(2,4) = HydroParam%ubv(lev,1,Nodes(4));   vvNode(2,4) = HydroParam%ubv(lev,2,Nodes(4));   wwNode(2,4) = HydroParam%ubv(lev,3,Nodes(4))

			    xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(1,1) = HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			    xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2));   yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2));   zzNode(1,2) = HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			    xxNode(1,3) = MeshParam%xNode(Nodes(3));        yyNode(1,3) = MeshParam%yNode(Nodes(3));        zzNode(1,3) = HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			    xxNode(1,4) = MeshParam%xNode(Nodes(4));        yyNode(1,4) = MeshParam%yNode(Nodes(4));        zzNode(1,4) = HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			    If (lev== HydroParam%ElCapitalM(nnel)) Then
                    xxNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = 0.5d0*(HydroParam%Z(lev,Nodes(1)) + HydroParam%hj(Nodes(1)) + HydroParam%Z(lev+1,Nodes(1)))
			        xxNode(2,2) = MeshParam%EdgeBary(2,Nodes(2));   yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2));   zzNode(2,2) = 0.5d0*(HydroParam%Z(lev,Nodes(2)) + HydroParam%hj(Nodes(2)) + HydroParam%Z(lev+1,Nodes(2)))
			        xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = 0.5d0*(HydroParam%peta(Nodes(3)) + HydroParam%hb(nnel) + HydroParam%Ze(lev,nnel))
			        xxNode(2,4) = MeshParam%xNode(Nodes(4));        yyNode(2,4) = MeshParam%yNode(Nodes(4));        zzNode(2,4) = 0.5d0*(HydroParam%peta(Nodes(4)) + HydroParam%hb(nnel) + HydroParam%Ze(lev,nnel))
                Else
		            xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
		            xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2));   yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2));   zzNode(2,2) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
		            xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
		            xxNode(2,4) = MeshParam%xNode(Nodes(4));        yyNode(2,4) = MeshParam%yNode(Nodes(4));        zzNode(2,4) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
               EndIf
            Endif
        EndIf
     Else !NodeFlag == 2 -> n1 and n4 = Face; n2 and n3 = Nodes;
		 If (HydroParam%ElSmallm(nnel)==HydroParam%ElCapitalM(nnel)) Then !2D Case/ 1 Layer
            lev=jlev
            
			uuNode(1,1) = HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1) = HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1) = HydroParam%wfc(jlev,Nodes(1))
			uuNode(1,4) = HydroParam%uxy(jlev,1,Nodes(4));   vvNode(1,4) = HydroParam%uxy(jlev,2,Nodes(4));   wwNode(1,4) = HydroParam%wfc(jlev,Nodes(4))
			uuNode(1,3) = HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3) = HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3) = HydroParam%ubv(jlev,3,Nodes(3))
			uuNode(1,2) = HydroParam%ubv(jlev,1,Nodes(2));   vvNode(1,2) = HydroParam%ubv(jlev,2,Nodes(2));   wwNode(1,2) = HydroParam%ubv(jlev,3,Nodes(2))
			If (zt>=z0) Then !up
                lev=jlev+1
				uuNode(2,1) = HydroParam%ug(Nodes(1),lev);        vvNode(2,1) = HydroParam%vg(Nodes(1),lev);        wwNode(2,1) = HydroParam%wg(Nodes(1),lev)
				uuNode(2,4) = HydroParam%ug(Nodes(4),lev);        vvNode(2,4) = HydroParam%vg(Nodes(4),lev);        wwNode(2,4) = HydroParam%wg(Nodes(4),lev)
				uuNode(2,3) = HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3) = HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3) = HydroParam%uNode(lev,3,Nodes(3))
				uuNode(2,2) = HydroParam%uNode(lev,1,Nodes(2));   vvNode(2,2) = HydroParam%uNode(lev,2,Nodes(2));   wwNode(2,2) = HydroParam%uNode(lev,3,Nodes(2))
                lev=jlev
            Else !down
                lev=jlev
                uuNode(2,1) = HydroParam%ug(Nodes(1),lev);        vvNode(2,1) = HydroParam%vg(Nodes(1),lev);        wwNode(2,1) = HydroParam%wg(Nodes(1),lev)
				uuNode(2,4) = HydroParam%ug(Nodes(4),lev);        vvNode(2,4) = HydroParam%vg(Nodes(4),lev);        wwNode(2,4) = HydroParam%wg(Nodes(4),lev)
				uuNode(2,3) = HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3) = HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3) = HydroParam%uNode(lev,3,Nodes(3))
				uuNode(2,2) = HydroParam%uNode(lev,1,Nodes(2));   vvNode(2,2) = HydroParam%uNode(lev,2,Nodes(2));   wwNode(2,2) = HydroParam%uNode(lev,3,Nodes(2))
            EndIf
			xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(1,1)=0.5d0*(HydroParam%Z(jlev,Nodes(1)) + HydroParam%hj(Nodes(1)) + HydroParam%Z(jlev+1,Nodes(1)))
			xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4));   yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(1,4)=0.5d0*(HydroParam%Z(jlev,Nodes(4)) + HydroParam%hj(Nodes(4)) + HydroParam%Z(jlev+1,Nodes(4)))
			xxNode(1,3) = MeshParam%xNode(Nodes(3));        yyNode(1,3) = MeshParam%yNode(Nodes(3));        zzNode(1,3)=0.5d0*(HydroParam%peta(Nodes(3)) + HydroParam%hb(nnel) + HydroParam%Ze(jlev,nnel))
			xxNode(1,2) = MeshParam%xNode(Nodes(2));        yyNode(1,2) = MeshParam%yNode(Nodes(2));        zzNode(1,2)=0.5d0*(HydroParam%peta(Nodes(2)) + HydroParam%hb(nnel) + HydroParam%Ze(jlev,nnel))
			If (zt>=z0) then !up
                lev=jlev+1
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = HydroParam%Z(lev,Nodes(1))
				xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4));   yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(2,4) = HydroParam%Z(lev,Nodes(4))
				xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = HydroParam%peta(Nodes(3))
				xxNode(2,2) = MeshParam%xNode(Nodes(2));        yyNode(2,2) = MeshParam%yNode(Nodes(2));        zzNode(2,2) = HydroParam%peta(Nodes(2))
                lev=jlev
            Else !down
                lev=jlev
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = HydroParam%Ze(lev,nnel)+ HydroParam%hb(nnel)
				xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4));   yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(2,4) = HydroParam%Ze(lev,nnel)+ HydroParam%hb(nnel)
				xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = HydroParam%Ze(lev,nnel)+ HydroParam%hb(nnel)
				xxNode(2,2) = MeshParam%xNode(Nodes(2));        yyNode(2,2) = MeshParam%yNode(Nodes(2));        zzNode(2,2) = HydroParam%Ze(lev,nnel)+ HydroParam%hb(nnel)
			EndIf    
        Else
		    if (jlev==HydroParam%ElCapitalM(nnel)) then
			    if (zt>=z0) then !up
				    lev=jlev+1
			    else !(zt<z0) down
				    lev=jlev-1
			    endif
			    uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			    uuNode(1,4)=HydroParam%uxy(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%uxy(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%wfc(jlev,Nodes(4))
			    uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			    uuNode(1,2)=HydroParam%ubv(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%ubv(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%ubv(jlev,3,Nodes(2))
			    if (zt>=z0) then !up
				    uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				    uuNode(2,4)=HydroParam%ug(Nodes(4),lev);        vvNode(2,4)=HydroParam%vg(Nodes(4),lev);        wwNode(2,4)=HydroParam%wg(Nodes(4),lev)
				    uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				    uuNode(2,2)=HydroParam%uNode(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%uNode(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%uNode(lev,3,Nodes(2))
			    else !down
				    uuNode(2,1)=HydroParam%uxy(lev,1,Nodes(1));   vvNode(2,1)=HydroParam%uxy(lev,2,Nodes(1));   wwNode(2,1)=HydroParam%wfc(lev,Nodes(1))
				    uuNode(2,4)=HydroParam%uxy(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uxy(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%wfc(lev,Nodes(4))
				    uuNode(2,3)=HydroParam%ubv(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%ubv(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%ubv(lev,3,Nodes(3))
				    uuNode(2,2)=HydroParam%ubv(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%ubv(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%ubv(lev,3,Nodes(2))
			    endif
			    xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(1,1) = 0.5d0*(HydroParam%Z(jlev,Nodes(1)) + HydroParam%hj(Nodes(1)) + HydroParam%Z(jlev+1,Nodes(1)))
			    xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4));   yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(1,4) = 0.5d0*(HydroParam%Z(jlev,Nodes(4)) + HydroParam%hj(Nodes(4)) + HydroParam%Z(jlev+1,Nodes(4)))
			    xxNode(1,3) = MeshParam%xNode(Nodes(3));        yyNode(1,3) = MeshParam%yNode(Nodes(3));        zzNode(1,3) = 0.5d0*(HydroParam%peta(Nodes(3)) + HydroParam%hb(nnel) + HydroParam%Ze(jlev,nnel))
			    xxNode(1,2) = MeshParam%xNode(Nodes(2));        yyNode(1,2) = MeshParam%yNode(Nodes(2));        zzNode(1,2) = 0.5d0*(HydroParam%peta(Nodes(2)) + HydroParam%hb(nnel) + HydroParam%Ze(jlev,nnel))
			    if (zt>=z0) then !up
				    xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = HydroParam%Z(lev,Nodes(1))
				    xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4));   yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(2,4) = HydroParam%Z(lev,Nodes(4))
				    xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = HydroParam%peta(Nodes(3))
				    xxNode(2,2) = MeshParam%xNode(Nodes(2));        yyNode(2,2) = MeshParam%yNode(Nodes(2));        zzNode(2,2) = HydroParam%peta(Nodes(2))
			    else !down
				    xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
				    xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4));   yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(2,4) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
				    xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
				    xxNode(2,2) = MeshParam%xNode(Nodes(2));        yyNode(2,2) = MeshParam%yNode(Nodes(2));        zzNode(2,2) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
			    endif
		    elseif (jlev==HydroParam%ElSmallm(nnel)) then
			    if (zt>=z0) then !up
				    lev=jlev+1
			    else !(zt<z0) down
				    lev=jlev
			    endif
			    uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			    uuNode(1,4)=HydroParam%uxy(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%uxy(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%wfc(jlev,Nodes(4))
			    uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			    uuNode(1,2)=HydroParam%ubv(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%ubv(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%ubv(jlev,3,Nodes(2))
			    if (zt<z0) then !Down
				    uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				    uuNode(2,4)=HydroParam%ug(Nodes(4),lev);        vvNode(2,4)=HydroParam%vg(Nodes(4),lev);        wwNode(2,4)=HydroParam%wg(Nodes(4),lev)
				    uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				    uuNode(2,2)=HydroParam%uNode(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%uNode(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%uNode(lev,3,Nodes(2))
			    else !up
				    uuNode(2,1)=HydroParam%uxy(lev,1,Nodes(1));   vvNode(2,1)=HydroParam%uxy(lev,2,Nodes(1));   wwNode(2,1) = HydroParam%wfc(lev,Nodes(1))
				    uuNode(2,4)=HydroParam%uxy(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uxy(lev,2,Nodes(4));   wwNode(2,4) = HydroParam%wfc(lev,Nodes(4))
				    uuNode(2,3)=HydroParam%ubv(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%ubv(lev,2,Nodes(3));   wwNode(2,3) = HydroParam%ubv(lev,3,Nodes(3))
				    uuNode(2,2)=HydroParam%ubv(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%ubv(lev,2,Nodes(2));   wwNode(2,2) = HydroParam%ubv(lev,3,Nodes(2))
			    endif
			    xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(1,1) = HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			    xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4));   yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(1,4) = HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			    xxNode(1,3) = MeshParam%xNode(Nodes(3));        yyNode(1,3) = MeshParam%yNode(Nodes(3));        zzNode(1,3) = HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			    xxNode(1,2) = MeshParam%xNode(Nodes(2));        yyNode(1,2) = MeshParam%yNode(Nodes(2));        zzNode(1,2) = HydroParam%Zb(jlev,nnel) + HydroParam%hb(nnel)/2
			    if (zt<z0) then !Down
				    xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel)
				    xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4));   yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(2,4) = HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel)
				    xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel)
				    xxNode(2,2) = MeshParam%xNode(Nodes(2));        yyNode(2,2) = MeshParam%yNode(Nodes(2));        zzNode(2,2) = HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel)
                else !Up
                    if (lev== HydroParam%ElCapitalM(nnel))then
                        xxNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = 0.5d0*(HydroParam%Z(lev,Nodes(1)) + HydroParam%hj(Nodes(1)) + HydroParam%Z(lev+1,Nodes(1)))
			            xxNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));   yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(2,4) = 0.5d0*(HydroParam%Z(lev,Nodes(4)) + HydroParam%hj(Nodes(4)) + HydroParam%Z(lev+1,Nodes(4)))
			            xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = 0.5d0*(HydroParam%peta(Nodes(3)) + HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel))
			            xxNode(2,2) = MeshParam%xNode(Nodes(2));        yyNode(2,2) = MeshParam%yNode(Nodes(2));        zzNode(2,2) = 0.5d0*(HydroParam%peta(Nodes(2)) + HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel))
                    else
				        xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
				        xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4));   yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(2,4) = HydroParam%Zb(lev,nnel)+ HydroParam%hb(nnel)/2
				        xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
				        xxNode(2,2) = MeshParam%xNode(Nodes(2));        yyNode(2,2) = MeshParam%yNode(Nodes(2));        zzNode(2,2) = HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
			        endif
                endif
		    else
			    if (zt>=z0) then !up
				    lev=jlev+1
			    else !(zt<z0) down
				    lev=jlev-1
			    endif
			    uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			    uuNode(1,4)=HydroParam%uxy(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%uxy(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%wfc(jlev,Nodes(4))
			    uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			    uuNode(1,2)=HydroParam%ubv(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%ubv(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%ubv(jlev,3,Nodes(2))

			    uuNode(2,1)=HydroParam%uxy(lev,1,Nodes(1));   vvNode(2,1)=HydroParam%uxy(lev,2,Nodes(1));   wwNode(2,1)=HydroParam%wfc(lev,Nodes(1))
			    uuNode(2,4)=HydroParam%uxy(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uxy(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%wfc(lev,Nodes(4))
			    uuNode(2,3)=HydroParam%ubv(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%ubv(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%ubv(lev,3,Nodes(3))
			    uuNode(2,2)=HydroParam%ubv(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%ubv(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%ubv(lev,3,Nodes(2))

			    xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(1,1)=HydroParam%Zb(jlev,nnel)
			    xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4));   yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(1,4)=HydroParam%Zb(jlev,nnel)
			    xxNode(1,3) = MeshParam%xNode(Nodes(3));        yyNode(1,3) = MeshParam%yNode(Nodes(3));        zzNode(1,3)=HydroParam%Zb(jlev,nnel)
			    xxNode(1,2) = MeshParam%xNode(Nodes(2));        yyNode(1,2) = MeshParam%yNode(Nodes(2));        zzNode(1,2)=HydroParam%Zb(jlev,nnel)
			
                if (lev == HydroParam%ElCapitalM(nnel))then 
                    xxNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1)=0.5d0*(HydroParam%Z(lev,Nodes(1)) + HydroParam%hj(Nodes(1)) + HydroParam%Z(lev+1,Nodes(1)))
			        xxNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));   yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(2,4)=0.5d0*(HydroParam%Z(lev,Nodes(4)) + HydroParam%hj(Nodes(4)) + HydroParam%Z(lev+1,Nodes(4)))
			        xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3)=0.5d0*(HydroParam%peta(Nodes(3)) + HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel))
			        xxNode(2,2) = MeshParam%xNode(Nodes(2));        yyNode(2,2) = MeshParam%yNode(Nodes(2));        zzNode(2,2)=0.5d0*(HydroParam%peta(Nodes(2)) + HydroParam%Ze(lev,nnel) + HydroParam%hb(nnel))
                else
			        xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1));   yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1));   zzNode(2,1)= HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
			        xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4));   yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));   zzNode(2,4)= HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
			        xxNode(2,3) = MeshParam%xNode(Nodes(3));        yyNode(2,3) = MeshParam%yNode(Nodes(3));        zzNode(2,3)= HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
			        xxNode(2,2) = MeshParam%xNode(Nodes(2));        yyNode(2,2) = MeshParam%yNode(Nodes(2));        zzNode(2,2)= HydroParam%Zb(lev,nnel) + HydroParam%hb(nnel)/2
                endif
            
            endif
        endif
	endif
    
    return
 End Subroutine FuVelocities3
 
 
    Subroutine FuVelocities2(uuNode, vvNode, wwNode, xxNode, yyNode, zzNode, nnel,jlev, xt,yt,zt,x0,y0,z0, id0, HydroParam,MeshParam)
    
    Use MeshVars !, Only: Quadri,EdgeDef,xNode,yNode,Area,Edge,EdgeNodes,EdgeBary,Left,Right
    Use Hydrodynamic !, Only: Smallm,ElSmallm,ElCapitalM,CapitalM,Z,DZj,Ze,DZi,uNode,uxy,Ze,DZi,w
    
    Implicit none
    
    integer, intent(in) :: nnel,jlev
    Real, intent(in) ::  xt,yt,zt,x0,y0,z0
    Real, intent(out) ::  uuNode(3,9), vvNode(3,9), wwNode(3,9), xxNode(3,9), yyNode(3,9), zzNode(3,9)
    Integer:: idt,iflqs1,i,j,l,nd,nn,lev,n1,n2,n3,n4, iEdge,iNode, Nodes(4), F1, F2, F3, F4, id0, NodeFlag
    Real:: NearZero = 1e-10 !< Small Number
    Integer:: i34 = 4
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    !Integer:: ELM_flag = 1
    
    F1=MeshParam%Edge(1,nnel)
	F2=MeshParam%Edge(2,nnel)
	F3=MeshParam%Edge(3,nnel)
	F4=MeshParam%Edge(4,nnel)
    !if (jlev==HydroParam%ElSmallm(nnel)) then and.jlev/=HydroParam%ElSmallm(nnel).and.jlev/=HydroParam%ElCapitalM(nnel)
    !MeshParam%Neighbor(1,nnel)
	if (id0==F1) then
        if (xt>=x0) then !Right side
            n1=F1
            n2=F3
            n3=MeshParam%Quadri(4,nnel) + 1
            n4=MeshParam%Quadri(1,nnel) + 1
			NodeFlag = 1
        else !left side
            n1=F1
            n2=MeshParam%Quadri(2,nnel) + 1
            n3=MeshParam%Quadri(3,nnel) + 1
            n4=F3
			NodeFlag = 2
        endif
        Nodes(1:4)= (/ n1, n2, n3, n4 /)
    elseif (id0==F2) then
        if (yt>=y0) then !North side
            n1=F2
            n2=F4
            n3=MeshParam%Quadri(1,nnel) + 1
            n4=MeshParam%Quadri(2,nnel) + 1
			NodeFlag = 1
        else !South side
            n1=F2
            n2=MeshParam%Quadri(3,nnel) + 1
            n3=MeshParam%Quadri(4,nnel) + 1
            n4=F4
			NodeFlag = 2
        endif
        Nodes(1:4)= (/ n1, n2, n3, n4 /)
    elseif (id0==F3) then
        if (xt>=x0) then !Right side
            n1=F3
            n2=MeshParam%Quadri(4,nnel) + 1
            n3=MeshParam%Quadri(1,nnel) + 1
            n4=F1
			NodeFlag = 2
        else !left side
            n1=F3
            n2=F1
            n3=MeshParam%Quadri(2,nnel) + 1
            n4=MeshParam%Quadri(3,nnel) + 1
			NodeFlag = 1
        endif
        Nodes(1:4)= (/ n1, n2, n3, n4 /)
	else !(id0==F4)
        if (yt>=y0) then !North side
            n1=F4
            n2=MeshParam%Quadri(1,nnel) + 1
            n3=MeshParam%Quadri(2,nnel) + 1
            n4=F2
			NodeFlag = 2
        else !South side
            n1=F4
            n2=F2
            n3=MeshParam%Quadri(3,nnel) + 1
            n4=MeshParam%Quadri(4,nnel) + 1	
			NodeFlag = 1
        endif
        Nodes(1:4)= (/ n1, n2, n3, n4 /)
    endif
	 if (NodeFlag == 1) then !n1 and n2 = Face; n3 and n4 = Nodes;
        If (HydroParam%ElSmallm(nnel)==HydroParam%ElCapitalM(nnel)) Then
            lev=jlev
                
			uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			uuNode(1,2)=HydroParam%uxy(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%uxy(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%wfc(jlev,Nodes(2))
			uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			uuNode(1,4)=HydroParam%ubv(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%ubv(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%ubv(jlev,3,Nodes(4))
			if (zt<z0) then !Down
                lev=jlev
				uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				uuNode(2,2)=HydroParam%ug(Nodes(2),lev);        vvNode(2,2)=HydroParam%vg(Nodes(2),lev);        wwNode(2,2)=HydroParam%wg(Nodes(2),lev)
				uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				uuNode(2,4)=HydroParam%uNode(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uNode(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%uNode(lev,3,Nodes(4))
            else !up
                lev=jlev+1
				uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				uuNode(2,2)=HydroParam%ug(Nodes(2),lev);        vvNode(2,2)=HydroParam%vg(Nodes(2),lev);        wwNode(2,2)=HydroParam%wg(Nodes(2),lev)
				uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				uuNode(2,4)=HydroParam%uNode(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uNode(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%uNode(lev,3,Nodes(4))
                lev=jlev
			endif
			xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(1,1)=HydroParam%Zb(jlev,nnel)
			xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=HydroParam%Zb(jlev,nnel)
			xxNode(1,3) = MeshParam%xNode(Nodes(3)); yyNode(1,3) = MeshParam%yNode(Nodes(3)); zzNode(1,3)=HydroParam%Zb(jlev,nnel)
			xxNode(1,4) = MeshParam%xNode(Nodes(4)); yyNode(1,4) = MeshParam%yNode(Nodes(4)); zzNode(1,4)=HydroParam%Zb(jlev,nnel)
			if (zt<z0) then !Down
                lev=jlev
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=HydroParam%Ze(lev,nnel)
				xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=HydroParam%Ze(lev,nnel)
				xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%Ze(lev,nnel)
				xxNode(2,4) = MeshParam%xNode(Nodes(4)); yyNode(2,4) = MeshParam%yNode(Nodes(4)); zzNode(2,4)= HydroParam%Ze(lev,nnel)
            else !Up
                lev=jlev+1
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=HydroParam%Z(lev,Nodes(1))
				xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=HydroParam%Z(lev,Nodes(2))
				xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%peta(Nodes(3))
				xxNode(2,4) = MeshParam%xNode(Nodes(4)); yyNode(2,4) = MeshParam%yNode(Nodes(4)); zzNode(2,4)= HydroParam%peta(Nodes(4))
                lev=jlev
			endif
        else
		if (jlev==HydroParam%ElCapitalM(nnel)) then

			if (zt>=z0) then !up
				lev=jlev+1
			else !(zt<z0) down
				lev=jlev-1
			endif
			uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			uuNode(1,2)=HydroParam%uxy(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%uxy(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%wfc(jlev,Nodes(2))
			uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			uuNode(1,4)=HydroParam%ubv(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%ubv(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%ubv(jlev,3,Nodes(4))
			if (zt>=z0) then !up
				uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				uuNode(2,2)=HydroParam%ug(Nodes(2),lev);        vvNode(2,2)=HydroParam%vg(Nodes(2),lev);        wwNode(2,2)=HydroParam%wg(Nodes(2),lev)
				uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				uuNode(2,4)=HydroParam%uNode(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uNode(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%uNode(lev,3,Nodes(4))
			else !down
				uuNode(2,1)=HydroParam%uxy(lev,1,Nodes(1));   vvNode(2,1)=HydroParam%uxy(lev,2,Nodes(1));   wwNode(2,1)=HydroParam%wfc(lev,Nodes(1))
				uuNode(2,2)=HydroParam%uxy(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%uxy(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%wfc(lev,Nodes(2))
				uuNode(2,3)=HydroParam%ubv(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%ubv(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%ubv(lev,3,Nodes(3))
				uuNode(2,4)=HydroParam%ubv(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%ubv(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%ubv(lev,3,Nodes(4))
			endif
			xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(1,1)=0.5d0*(HydroParam%Z(jlev,Nodes(1))+HydroParam%Z(jlev+1,Nodes(1)))
			xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=0.5d0*(HydroParam%Z(jlev,Nodes(2))+HydroParam%Z(jlev+1,Nodes(2)))
			xxNode(1,3) = MeshParam%xNode(Nodes(3)); yyNode(1,3) = MeshParam%yNode(Nodes(3)); zzNode(1,3)=0.5d0*(HydroParam%peta(Nodes(3))+HydroParam%Ze(jlev,nnel))
			xxNode(1,4) = MeshParam%xNode(Nodes(4)); yyNode(1,4) = MeshParam%yNode(Nodes(4)); zzNode(1,4)=0.5d0*(HydroParam%peta(Nodes(4))+HydroParam%Ze(jlev,nnel))
			if (zt>=z0) then !up
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=HydroParam%Z(lev,Nodes(1))
				xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=HydroParam%Z(lev,Nodes(2))
				xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%peta(Nodes(3))
				xxNode(2,4) = MeshParam%xNode(Nodes(4)); yyNode(2,4) = MeshParam%yNode(Nodes(4)); zzNode(2,4)= HydroParam%peta(Nodes(4))
            else !down
				    xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)= HydroParam%Zb(lev,nnel)
				    xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)= HydroParam%Zb(lev,nnel)
				    xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%Zb(lev,nnel)
				    xxNode(2,4) = MeshParam%xNode(Nodes(4)); yyNode(2,4) = MeshParam%yNode(Nodes(4)); zzNode(2,4)= HydroParam%Zb(lev,nnel)
            endif
           
		elseif (jlev==HydroParam%ElSmallm(nnel)) then
			if (zt>=z0) then !up
				lev=jlev+1
			else !(zt<z0) down
				lev=jlev
			endif
			uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			uuNode(1,2)=HydroParam%uxy(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%uxy(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%wfc(jlev,Nodes(2))
			uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			uuNode(1,4)=HydroParam%ubv(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%ubv(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%ubv(jlev,3,Nodes(4))
			if (zt<z0) then !Down
                
				uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				uuNode(2,2)=HydroParam%ug(Nodes(2),lev);        vvNode(2,2)=HydroParam%vg(Nodes(2),lev);        wwNode(2,2)=HydroParam%wg(Nodes(2),lev)
				uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				uuNode(2,4)=HydroParam%uNode(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uNode(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%uNode(lev,3,Nodes(4))
            else !up
				uuNode(2,1)=HydroParam%uxy(lev,1,Nodes(1));   vvNode(2,1)=HydroParam%uxy(lev,2,Nodes(1));   wwNode(2,1)=HydroParam%wfc(lev,Nodes(1))
				uuNode(2,2)=HydroParam%uxy(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%uxy(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%wfc(lev,Nodes(2))
				uuNode(2,3)=HydroParam%ubv(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%ubv(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%ubv(lev,3,Nodes(3))
				uuNode(2,4)=HydroParam%ubv(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%ubv(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%ubv(lev,3,Nodes(4))
                
			endif
			xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(1,1)=HydroParam%Zb(jlev,nnel)
			xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=HydroParam%Zb(jlev,nnel)
			xxNode(1,3) = MeshParam%xNode(Nodes(3)); yyNode(1,3) = MeshParam%yNode(Nodes(3)); zzNode(1,3)=HydroParam%Zb(jlev,nnel)
			xxNode(1,4) = MeshParam%xNode(Nodes(4)); yyNode(1,4) = MeshParam%yNode(Nodes(4)); zzNode(1,4)=HydroParam%Zb(jlev,nnel)
			if (zt<z0) then !Down
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=HydroParam%Ze(lev,nnel)
				xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=HydroParam%Ze(lev,nnel)
				xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%Ze(lev,nnel)
				xxNode(2,4) = MeshParam%xNode(Nodes(4)); yyNode(2,4) = MeshParam%yNode(Nodes(4)); zzNode(2,4)= HydroParam%Ze(lev,nnel)
            else !Up
                if (lev== HydroParam%ElCapitalM(nnel)) then
                    xxNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=0.5d0*(HydroParam%Z(lev,Nodes(1))+HydroParam%Z(lev+1,Nodes(1)))
			        xxNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=0.5d0*(HydroParam%Z(lev,Nodes(2))+HydroParam%Z(lev+1,Nodes(2)))
			        xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)=0.5d0*(HydroParam%peta(Nodes(3))+HydroParam%Ze(lev,nnel))
			        xxNode(2,4) = MeshParam%xNode(Nodes(4)); yyNode(2,4) = MeshParam%yNode(Nodes(4)); zzNode(2,4)=0.5d0*(HydroParam%peta(Nodes(4))+HydroParam%Ze(lev,nnel))
                else
				    xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)= HydroParam%Zb(lev,nnel)
				    xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)= HydroParam%Zb(lev,nnel)
				    xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%Zb(lev,nnel)
				    xxNode(2,4) = MeshParam%xNode(Nodes(4)); yyNode(2,4) = MeshParam%yNode(Nodes(4)); zzNode(2,4)= HydroParam%Zb(lev,nnel)
			    endif
            endif
		else 
			if (zt>=z0) then !up
				lev=jlev+1
			else !(zt<z0) down
				lev=jlev-1
			endif
			uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			uuNode(1,2)=HydroParam%uxy(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%uxy(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%wfc(jlev,Nodes(2))
			uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			uuNode(1,4)=HydroParam%ubv(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%ubv(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%ubv(jlev,3,Nodes(4))

			uuNode(2,1)=HydroParam%uxy(lev,1,Nodes(1));   vvNode(2,1)=HydroParam%uxy(lev,2,Nodes(1));   wwNode(2,1)=HydroParam%wfc(lev,Nodes(1))
			uuNode(2,2)=HydroParam%uxy(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%uxy(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%wfc(lev,Nodes(2))
			uuNode(2,3)=HydroParam%ubv(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%ubv(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%ubv(lev,3,Nodes(3))
			uuNode(2,4)=HydroParam%ubv(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%ubv(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%ubv(lev,3,Nodes(4))

			xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(1,1)=HydroParam%Zb(jlev,nnel)
			xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=HydroParam%Zb(jlev,nnel)
			xxNode(1,3) = MeshParam%xNode(Nodes(3)); yyNode(1,3) = MeshParam%yNode(Nodes(3)); zzNode(1,3)=HydroParam%Zb(jlev,nnel)
			xxNode(1,4) = MeshParam%xNode(Nodes(4)); yyNode(1,4) = MeshParam%yNode(Nodes(4)); zzNode(1,4)=HydroParam%Zb(jlev,nnel)
			if (lev== HydroParam%ElCapitalM(nnel))then
                xxNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=0.5d0*(HydroParam%Z(lev,Nodes(1))+HydroParam%Z(lev+1,Nodes(1)))
			    xxNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=0.5d0*(HydroParam%Z(lev,Nodes(2))+HydroParam%Z(lev+1,Nodes(2)))
			    xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)=0.5d0*(HydroParam%peta(Nodes(3))+HydroParam%Ze(lev,nnel))
			    xxNode(2,4) = MeshParam%xNode(Nodes(4)); yyNode(2,4) = MeshParam%yNode(Nodes(4)); zzNode(2,4)=0.5d0*(HydroParam%peta(Nodes(4))+HydroParam%Ze(lev,nnel))
            else
		        xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)= HydroParam%Zb(lev,nnel)
		        xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)= HydroParam%Zb(lev,nnel)
		        xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%Zb(lev,nnel)
		        xxNode(2,4) = MeshParam%xNode(Nodes(4)); yyNode(2,4) = MeshParam%yNode(Nodes(4)); zzNode(2,4)= HydroParam%Zb(lev,nnel)
           endif 
        endif
        endif
     else !NodeFlag == 2 -> n1 and n4 = Face; n2 and n3 = Nodes;
		 If (HydroParam%ElSmallm(nnel)==HydroParam%ElCapitalM(nnel)) Then
            lev=jlev
            
			uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			uuNode(1,4)=HydroParam%uxy(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%uxy(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%wfc(jlev,Nodes(4))
			uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			uuNode(1,2)=HydroParam%ubv(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%ubv(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%ubv(jlev,3,Nodes(2))
			if (zt>=z0) then !up
                lev=jlev+1
				uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				uuNode(2,4)=HydroParam%ug(Nodes(4),lev);        vvNode(2,4)=HydroParam%vg(Nodes(4),lev);        wwNode(2,4)=HydroParam%wg(Nodes(4),lev)
				uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				uuNode(2,2)=HydroParam%uNode(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%uNode(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%uNode(lev,3,Nodes(2))
                lev=jlev
            else !down
                lev=jlev
                uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				uuNode(2,4)=HydroParam%ug(Nodes(4),lev);        vvNode(2,4)=HydroParam%vg(Nodes(4),lev);        wwNode(2,4)=HydroParam%wg(Nodes(4),lev)
				uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				uuNode(2,2)=HydroParam%uNode(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%uNode(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%uNode(lev,3,Nodes(2))
            endif
			xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(1,1)=0.5d0*(HydroParam%Z(jlev,Nodes(1))+HydroParam%Z(jlev+1,Nodes(1)))
			xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4)=0.5d0*(HydroParam%Z(jlev,Nodes(4))+HydroParam%Z(jlev+1,Nodes(4)))
			xxNode(1,3) = MeshParam%xNode(Nodes(3)); yyNode(1,3) = MeshParam%yNode(Nodes(3)); zzNode(1,3)=0.5d0*(HydroParam%peta(Nodes(3))+HydroParam%Ze(jlev,nnel))
			xxNode(1,2) = MeshParam%xNode(Nodes(2)); yyNode(1,2) = MeshParam%yNode(Nodes(2)); zzNode(1,2)=0.5d0*(HydroParam%peta(Nodes(2))+HydroParam%Ze(jlev,nnel))
			if (zt>=z0) then !up
                lev=jlev+1
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=HydroParam%Z(lev,Nodes(1))
				xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(2,4)=HydroParam%Z(lev,Nodes(4))
				xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%peta(Nodes(3))
				xxNode(2,2) = MeshParam%xNode(Nodes(2)); yyNode(2,2) = MeshParam%yNode(Nodes(2)); zzNode(2,2)= HydroParam%peta(Nodes(2))
                lev=jlev
            else !down
                lev=jlev
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=HydroParam%Ze(lev,nnel)
				xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(2,4)=HydroParam%Ze(lev,nnel)
				xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%Ze(lev,nnel)
				xxNode(2,2) = MeshParam%xNode(Nodes(2)); yyNode(2,2) = MeshParam%yNode(Nodes(2)); zzNode(2,2)= HydroParam%Ze(lev,nnel)
			endif    
        else
         
		if (jlev==HydroParam%ElCapitalM(nnel)) then
			if (zt>=z0) then !up
				lev=jlev+1
			else !(zt<z0) down
				lev=jlev-1
			endif
			uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			uuNode(1,4)=HydroParam%uxy(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%uxy(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%wfc(jlev,Nodes(4))
			uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			uuNode(1,2)=HydroParam%ubv(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%ubv(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%ubv(jlev,3,Nodes(2))
			if (zt>=z0) then !up
				uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				uuNode(2,4)=HydroParam%ug(Nodes(4),lev);        vvNode(2,4)=HydroParam%vg(Nodes(4),lev);        wwNode(2,4)=HydroParam%wg(Nodes(4),lev)
				uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				uuNode(2,2)=HydroParam%uNode(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%uNode(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%uNode(lev,3,Nodes(2))
			else !down
				uuNode(2,1)=HydroParam%uxy(lev,1,Nodes(1));   vvNode(2,1)=HydroParam%uxy(lev,2,Nodes(1));   wwNode(2,1)=HydroParam%wfc(lev,Nodes(1))
				uuNode(2,4)=HydroParam%uxy(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uxy(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%wfc(lev,Nodes(4))
				uuNode(2,3)=HydroParam%ubv(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%ubv(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%ubv(lev,3,Nodes(3))
				uuNode(2,2)=HydroParam%ubv(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%ubv(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%ubv(lev,3,Nodes(2))
			endif
			xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(1,1)=0.5d0*(HydroParam%Z(jlev,Nodes(1))+HydroParam%Z(jlev+1,Nodes(1)))
			xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4)=0.5d0*(HydroParam%Z(jlev,Nodes(4))+HydroParam%Z(jlev+1,Nodes(4)))
			xxNode(1,3) = MeshParam%xNode(Nodes(3)); yyNode(1,3) = MeshParam%yNode(Nodes(3)); zzNode(1,3)=0.5d0*(HydroParam%peta(Nodes(3))+HydroParam%Ze(jlev,nnel))
			xxNode(1,2) = MeshParam%xNode(Nodes(2)); yyNode(1,2) = MeshParam%yNode(Nodes(2)); zzNode(1,2)=0.5d0*(HydroParam%peta(Nodes(2))+HydroParam%Ze(jlev,nnel))
			if (zt>=z0) then !up
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=HydroParam%Z(lev,Nodes(1))
				xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(2,4)=HydroParam%Z(lev,Nodes(4))
				xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%peta(Nodes(3))
				xxNode(2,2) = MeshParam%xNode(Nodes(2)); yyNode(2,2) = MeshParam%yNode(Nodes(2)); zzNode(2,2)= HydroParam%peta(Nodes(2))
			else !down
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)= HydroParam%Zb(lev,nnel)
				xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(2,4)= HydroParam%Zb(lev,nnel)
				xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%Zb(lev,nnel)
				xxNode(2,2) = MeshParam%xNode(Nodes(2)); yyNode(2,2) = MeshParam%yNode(Nodes(2)); zzNode(2,2)= HydroParam%Zb(lev,nnel)
			endif
		elseif (jlev==HydroParam%ElSmallm(nnel)) then
			if (zt>=z0) then !up
				lev=jlev+1
			else !(zt<z0) down
				lev=jlev
			endif
			uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			uuNode(1,4)=HydroParam%uxy(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%uxy(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%wfc(jlev,Nodes(4))
			uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			uuNode(1,2)=HydroParam%ubv(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%ubv(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%ubv(jlev,3,Nodes(2))
			if (zt<z0) then !Down
				uuNode(2,1)=HydroParam%ug(Nodes(1),lev);        vvNode(2,1)=HydroParam%vg(Nodes(1),lev);        wwNode(2,1)=HydroParam%wg(Nodes(1),lev)
				uuNode(2,4)=HydroParam%ug(Nodes(4),lev);        vvNode(2,4)=HydroParam%vg(Nodes(4),lev);        wwNode(2,4)=HydroParam%wg(Nodes(4),lev)
				uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))
				uuNode(2,2)=HydroParam%uNode(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%uNode(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%uNode(lev,3,Nodes(2))
			else !up
				uuNode(2,1)=HydroParam%uxy(lev,1,Nodes(1));   vvNode(2,1)=HydroParam%uxy(lev,2,Nodes(1));   wwNode(2,1)=HydroParam%wfc(lev,Nodes(1))
				uuNode(2,4)=HydroParam%uxy(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uxy(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%wfc(lev,Nodes(4))
				uuNode(2,3)=HydroParam%ubv(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%ubv(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%ubv(lev,3,Nodes(3))
				uuNode(2,2)=HydroParam%ubv(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%ubv(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%ubv(lev,3,Nodes(2))
			endif
			xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(1,1)=HydroParam%Zb(jlev,nnel)
			xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4)=HydroParam%Zb(jlev,nnel)
			xxNode(1,3) = MeshParam%xNode(Nodes(3)); yyNode(1,3) = MeshParam%yNode(Nodes(3)); zzNode(1,3)=HydroParam%Zb(jlev,nnel)
			xxNode(1,2) = MeshParam%xNode(Nodes(2)); yyNode(1,2) = MeshParam%yNode(Nodes(2)); zzNode(1,2)=HydroParam%Zb(jlev,nnel)
			if (zt<z0) then !Down
				xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=HydroParam%Ze(lev,nnel)
				xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(2,4)=HydroParam%Ze(lev,nnel)
				xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%Ze(lev,nnel)
				xxNode(2,2) = MeshParam%xNode(Nodes(2)); yyNode(2,2) = MeshParam%yNode(Nodes(2)); zzNode(2,2)= HydroParam%Ze(lev,nnel)
            else !Up
                if (lev== HydroParam%ElCapitalM(nnel))then
                    xxNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=0.5d0*(HydroParam%Z(lev,Nodes(1))+HydroParam%Z(lev+1,Nodes(1)))
			        xxNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(2,4)=0.5d0*(HydroParam%Z(lev,Nodes(4))+HydroParam%Z(lev+1,Nodes(4)))
			        xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)=0.5d0*(HydroParam%peta(Nodes(3))+HydroParam%Ze(lev,nnel))
			        xxNode(2,2) = MeshParam%xNode(Nodes(2)); yyNode(2,2) = MeshParam%yNode(Nodes(2)); zzNode(2,2)=0.5d0*(HydroParam%peta(Nodes(2))+HydroParam%Ze(lev,nnel))
                else
				    xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)= HydroParam%Zb(lev,nnel)
				    xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(2,4)= HydroParam%Zb(lev,nnel)
				    xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%Zb(lev,nnel)
				    xxNode(2,2) = MeshParam%xNode(Nodes(2)); yyNode(2,2) = MeshParam%yNode(Nodes(2)); zzNode(2,2)= HydroParam%Zb(lev,nnel)
			    endif
            endif
		else
			if (zt>=z0) then !up
				lev=jlev+1
			else !(zt<z0) down
				lev=jlev-1
			endif
			uuNode(1,1)=HydroParam%uxy(jlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxy(jlev,2,Nodes(1));   wwNode(1,1)=HydroParam%wfc(jlev,Nodes(1))
			uuNode(1,4)=HydroParam%uxy(jlev,1,Nodes(4));   vvNode(1,4)=HydroParam%uxy(jlev,2,Nodes(4));   wwNode(1,4)=HydroParam%wfc(jlev,Nodes(4))
			uuNode(1,3)=HydroParam%ubv(jlev,1,Nodes(3));   vvNode(1,3)=HydroParam%ubv(jlev,2,Nodes(3));   wwNode(1,3)=HydroParam%ubv(jlev,3,Nodes(3))
			uuNode(1,2)=HydroParam%ubv(jlev,1,Nodes(2));   vvNode(1,2)=HydroParam%ubv(jlev,2,Nodes(2));   wwNode(1,2)=HydroParam%ubv(jlev,3,Nodes(2))

			uuNode(2,1)=HydroParam%uxy(lev,1,Nodes(1));   vvNode(2,1)=HydroParam%uxy(lev,2,Nodes(1));   wwNode(2,1)=HydroParam%wfc(lev,Nodes(1))
			uuNode(2,4)=HydroParam%uxy(lev,1,Nodes(4));   vvNode(2,4)=HydroParam%uxy(lev,2,Nodes(4));   wwNode(2,4)=HydroParam%wfc(lev,Nodes(4))
			uuNode(2,3)=HydroParam%ubv(lev,1,Nodes(3));   vvNode(2,3)=HydroParam%ubv(lev,2,Nodes(3));   wwNode(2,3)=HydroParam%ubv(lev,3,Nodes(3))
			uuNode(2,2)=HydroParam%ubv(lev,1,Nodes(2));   vvNode(2,2)=HydroParam%ubv(lev,2,Nodes(2));   wwNode(2,2)=HydroParam%ubv(lev,3,Nodes(2))

			xxNode(1,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(1,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(1,1)=HydroParam%Zb(jlev,nnel)
			xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4)=HydroParam%Zb(jlev,nnel)
			xxNode(1,3) = MeshParam%xNode(Nodes(3)); yyNode(1,3) = MeshParam%yNode(Nodes(3)); zzNode(1,3)=HydroParam%Zb(jlev,nnel)
			xxNode(1,2) = MeshParam%xNode(Nodes(2)); yyNode(1,2) = MeshParam%yNode(Nodes(2)); zzNode(1,2)=HydroParam%Zb(jlev,nnel)
			
            if (lev == HydroParam%ElCapitalM(nnel))then 
                xxNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)=0.5d0*(HydroParam%Z(lev,Nodes(1))+HydroParam%Z(lev+1,Nodes(1)))
			    xxNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(2,4)=0.5d0*(HydroParam%Z(lev,Nodes(4))+HydroParam%Z(lev+1,Nodes(4)))
			    xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)=0.5d0*(HydroParam%peta(Nodes(3))+HydroParam%Ze(lev,nnel))
			    xxNode(2,2) = MeshParam%xNode(Nodes(2)); yyNode(2,2) = MeshParam%yNode(Nodes(2)); zzNode(2,2)=0.5d0*(HydroParam%peta(Nodes(2))+HydroParam%Ze(lev,nnel))
            else
			    xxNode(2,1) = MeshParam%EdgeBary(1,Nodes(1)); yyNode(2,1) = MeshParam%EdgeBary(2,Nodes(1)); zzNode(2,1)= HydroParam%Zb(lev,nnel)
			    xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(2,4)= HydroParam%Zb(lev,nnel)
			    xxNode(2,3) = MeshParam%xNode(Nodes(3)); yyNode(2,3) = MeshParam%yNode(Nodes(3)); zzNode(2,3)= HydroParam%Zb(lev,nnel)
			    xxNode(2,2) = MeshParam%xNode(Nodes(2)); yyNode(2,2) = MeshParam%yNode(Nodes(2)); zzNode(2,2)= HydroParam%Zb(lev,nnel)
            endif
            
        endif
        endif
	endif
    
    return
    End Subroutine FuVelocities2

    Subroutine FwVelocities2(uuNode, vvNode, wwNode, xxNode, yyNode, zzNode, nnel,jlev, jjlev, xt,yt,zt,x0,y0,z0, HydroParam,MeshParam)
    
    Use MeshVars !, Only: Quadri,EdgeDef,xNode,yNode,Area,Edge,EdgeNodes,EdgeBary,Left,Right
    Use Hydrodynamic !, Only: Smallm,ElSmallm,ElCapitalM,CapitalM,Z,DZj,Ze,DZi,uNode,uxy,Ze,DZi,w
    
    Implicit none
    
    integer, intent(in) :: nnel,jlev, jjlev
    Real, intent(in) ::  xt,yt,zt,x0,y0,z0 
    Real, intent(out) ::  uuNode(3,9), vvNode(3,9), wwNode(3,9), xxNode(3,9), yyNode(3,9), zzNode(3,9)
    Integer:: idt,iflqs1,i,j,l,nd,nn,lev,n1,n2,n3,n4, iEdge,iNode, Nodes(4), N, S
    Real:: NearZero = 1e-10 !< Small Number
    Integer:: i34 = 4
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    !Integer:: ELM_flag = 1
    
    if (xt<=x0) then !west side
	    if (yt>=y0) then!North Side in Y
		    n1=nnel
		    n2=MeshParam%Edge(1,nnel)
		    n3=MeshParam%Quadri(2,nnel) + 1
		    n4=MeshParam%Edge(2,nnel)
		    Nodes(1:4)= (/ n1, n2, n3, n4 /)
	    else !South Side in Y
		    n1=nnel
		    n2=MeshParam%Edge(2,nnel)
		    n3=MeshParam%Quadri(3,nnel) + 1
		    n4=MeshParam%Edge(3,nnel)
		    Nodes(1:4)= (/ n1, n2, n3, n4 /)
	    endif !end Y analyses for (O==0.and.L/=0)
    else !east side
	    if (yt>=y0) then!North Side in Y
		    n1=nnel
		    n2=MeshParam%Edge(4,nnel)
		    n3=MeshParam%Quadri(1,nnel) + 1
		    n4=MeshParam%Edge(1,nnel)
		    Nodes(1:4)= (/ n1, n2, n3, n4 /)
	    else !South Side in Y
		    n1=nnel
		    n2=MeshParam%Edge(3,nnel)
		    n3=MeshParam%Quadri(4,nnel) + 1
		    n4=MeshParam%Edge(4,nnel)
		    Nodes(1:4)= (/ n1, n2, n3, n4 /)
	    endif !end Y analyses for (O==0.and.L/=0)
    endif !end X analyses for (O==0.and.L/=0
    
    If (HydroParam%ElSmallm(nnel)==HydroParam%ElCapitalM(nnel)) Then
        continue
        if (zt<z0) then !up
            continue
        endif
    endif
           
    if (zt>=z0) then !up
		lev=jjlev+1
    else !(zt<z0) down
		lev=jjlev-1
	endif

    uuNode(1,1)=HydroParam%uxyL(jjlev,1,Nodes(1));   vvNode(1,1)=HydroParam%uxyL(jjlev,2,Nodes(1));   wwNode(1,1)=HydroParam%w(jjlev,Nodes(1))
    uuNode(1,2)=HydroParam%ug(Nodes(2),jjlev);       vvNode(1,2)=HydroParam%vg(Nodes(2),jjlev);       wwNode(1,2)=HydroParam%wg(Nodes(2),jjlev)
    uuNode(1,4)=HydroParam%ug(Nodes(4),jjlev);       vvNode(1,4)=HydroParam%vg(Nodes(4),jjlev);       wwNode(1,4)=HydroParam%wg(Nodes(4),jjlev)
    uuNode(1,3)=HydroParam%uNode(jjlev,1,Nodes(3));  vvNode(1,3)=HydroParam%uNode(jjlev,2,Nodes(3));  wwNode(1,3)=HydroParam%uNode(jjlev,3,Nodes(3))

    uuNode(2,1)=HydroParam%uxyL(lev,1,Nodes(1));   vvNode(2,1)=HydroParam%uxyL(lev,2,Nodes(1));   wwNode(2,1)=HydroParam%w(lev,Nodes(1))
    uuNode(2,2)=HydroParam%ug(Nodes(2),lev);       vvNode(2,2)=HydroParam%vg(Nodes(2),lev);       wwNode(2,2)=HydroParam%wg(Nodes(2),lev)
    uuNode(2,4)=HydroParam%ug(Nodes(4),lev);       vvNode(2,4)=HydroParam%vg(Nodes(4),lev);       wwNode(2,4)=HydroParam%wg(Nodes(4),lev)
    uuNode(2,3)=HydroParam%uNode(lev,1,Nodes(3));  vvNode(2,3)=HydroParam%uNode(lev,2,Nodes(3));  wwNode(2,3)=HydroParam%uNode(lev,3,Nodes(3))

    if (jlev==HydroParam%ElCapitalM(nnel)) then !Top Layer
	    xxNode(1,1) = MeshParam%xb(Nodes(1)); 		  yyNode(1,1) = MeshParam%yb(Nodes(1)); 		zzNode(1,1)=HydroParam%Ze(jjlev,nnel)
	    xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=HydroParam%Ze(jjlev,nnel)
	    xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4)=HydroParam%Ze(jjlev,nnel)
	    xxNode(1,3) = MeshParam%xNode(Nodes(3));      yyNode(1,3) = MeshParam%yNode(Nodes(3));      zzNode(1,3)=HydroParam%Ze(jjlev,nnel)

	    xxNode(2,1) = MeshParam%xb(Nodes(1)); 		  yyNode(2,1) = MeshParam%yb(Nodes(1)); 		zzNode(2,1)=HydroParam%eta(nnel)
	    xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=HydroParam%Z(jlev+1,Nodes(2))
	    xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(2,4)=HydroParam%Z(jlev+1,Nodes(4))
	    xxNode(2,3) = MeshParam%xNode(Nodes(3));      yyNode(2,3) = MeshParam%yNode(Nodes(3));      zzNode(2,3)=HydroParam%peta(Nodes(3))
    else
	    xxNode(1,1) = MeshParam%xb(Nodes(1)); 		  yyNode(1,1) = MeshParam%yb(Nodes(1)); 		zzNode(1,1)=HydroParam%Ze(jjlev,nnel)
	    xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=HydroParam%Ze(jjlev,nnel)
	    xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4)=HydroParam%Ze(jjlev,nnel)
	    xxNode(1,3) = MeshParam%xNode(Nodes(3));      yyNode(1,3) = MeshParam%yNode(Nodes(3));      zzNode(1,3)=HydroParam%Ze(jjlev,nnel)

	    xxNode(2,1) = MeshParam%xb(Nodes(1)); 		  yyNode(2,1) = MeshParam%yb(Nodes(1)); 		zzNode(2,1)=HydroParam%Ze(lev,nnel)
	    xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=HydroParam%Ze(lev,nnel)
	    xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(2,4)=HydroParam%Ze(lev,nnel)
	    xxNode(2,3) = MeshParam%xNode(Nodes(3));      yyNode(2,3) = MeshParam%yNode(Nodes(3));      zzNode(2,3)=HydroParam%Ze(lev,nnel)
    endif
    return
    End Subroutine FwVelocities2

    Subroutine iQuadraticNodes(uuNode, vvNode, wwNode, uuNodet, vvNodet, wwNodet, xxNode, yyNode, zzNode, bbElem, bbLayer, HydroParam, MeshParam)
    
    Use MeshVars !, Only: 
    Use Hydrodynamic ! Only:
    
    Implicit none
    
    Integer, intent(in) :: bbElem, bbLayer
    Real, intent(inout) ::  uuNode(3,9), vvNode(3,9), wwNode(3,9), uuNodet(3,9), vvNodet(3,9), wwNodet(3,9), xxNode(3,9), yyNode(3,9), zzNode(3,9)
    Integer :: Nodes(9)
    Real:: NearZero = 1e-5
    Integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
            
    !iQuadratic Interpolation Nodes in bbLayer from bbElem
    ! For nodes positions in the vector, the Standard is this:
    !        n6
    !  n3 .--.--. n9
    !     |     |
    !  n2 .  .  . n8
    !     | n5  |      
    !  n1 .--.--. n7
    !        n4            
    n1 = MeshParam%Quadri(3,bbElem) + 1
    n2 = MeshParam%Edge(2,bbElem)
    n3 = MeshParam%Quadri(2,bbElem) + 1
    n4 = MeshParam%Edge(3,bbElem)
    n5 = bbElem
    n6 = MeshParam%Edge(1,bbElem)
    n7 = MeshParam%Quadri(4,bbElem) + 1
    n8 = MeshParam%Edge(4,bbElem)
    n9 = MeshParam%Quadri(1,bbElem) + 1
    Nodes(1:9)= (/n1, n2, n3, n4, n5, n6, n7, n8, n9 /)
    ! Velocities by vertical section and by layers, in  bottom to top (k-1/2 to k + 1/2). See Figure 3 in [2]):
    ! Vertical 1 = West Edge
    ! uuNode(Layer,Node Position), Layer = [k-1/2, k , k +1/2] == [1,2,3]
    uuNode(1,1)=HydroParam%uNode(bbLayer,1,Nodes(1));   vvNode(1,1)=HydroParam%uNode(bbLayer,2,Nodes(1));   wwNode(1,1)=HydroParam%uNode(bbLayer,3,Nodes(1))
    uuNode(1,2)=HydroParam%ug(Nodes(2),bbLayer);        vvNode(1,2)=HydroParam%vg(Nodes(2),bbLayer);        wwNode(1,2)=HydroParam%wg(Nodes(2),bbLayer)
    uuNode(1,3)=HydroParam%uNode(bbLayer,1,Nodes(3));   vvNode(1,3)=HydroParam%uNode(bbLayer,2,Nodes(3));   wwNode(1,3)=HydroParam%uNode(bbLayer,3,Nodes(3))
    uuNode(2,1)=HydroParam%ubV(bbLayer,1,Nodes(1));     vvNode(2,1)=HydroParam%ubV(bbLayer,2,Nodes(1));     wwNode(2,1)=HydroParam%ubV(bbLayer,3,Nodes(1))
    uuNode(2,2)=HydroParam%uxy(bbLayer,1,Nodes(2));     vvNode(2,2)=HydroParam%uxy(bbLayer,2,Nodes(2));     wwNode(2,2)=HydroParam%wfc(bbLayer,Nodes(2))
    uuNode(2,3)=HydroParam%ubV(bbLayer,1,Nodes(3));     vvNode(2,3)=HydroParam%ubV(bbLayer,2,Nodes(3));     wwNode(2,3)=HydroParam%ubV(bbLayer,3,Nodes(3))
    uuNode(3,1)=HydroParam%uNode(bbLayer+1,1,Nodes(1)); vvNode(3,1)=HydroParam%uNode(bbLayer+1,2,Nodes(1)); wwNode(3,1)=HydroParam%uNode(bbLayer+1,3,Nodes(1))
    uuNode(3,2)=HydroParam%ug(Nodes(2),bbLayer+1);      vvNode(3,2)=HydroParam%vg(Nodes(2),bbLayer+1);      wwNode(3,2)=HydroParam%wg(Nodes(2),bbLayer+1)
    uuNode(3,3)=HydroParam%uNode(bbLayer+1,1,Nodes(3)); vvNode(3,3)=HydroParam%uNode(bbLayer+1,2,Nodes(3)); wwNode(3,3)=HydroParam%uNode(bbLayer+1,3,Nodes(3))
    ! Vertical 2 = Cell Centered Section
    uuNode(1,4)=HydroParam%ug(Nodes(4),bbLayer);      vvNode(1,4)=HydroParam%vg(Nodes(4),bbLayer);      	  wwNode(1,4)=HydroParam%wg(Nodes(4),bbLayer)
    uuNode(1,5)=HydroParam%uxyL(bbLayer,1,Nodes(5));  vvNode(1,5)=HydroParam%uxyL(bbLayer,2,Nodes(5));  	  wwNode(1,5)=HydroParam%w(bbLayer,Nodes(5))
    uuNode(1,6)=HydroParam%ug(Nodes(6),bbLayer);      vvNode(1,6)=HydroParam%vg(Nodes(6),bbLayer);      	  wwNode(1,6)=HydroParam%wg(Nodes(6),bbLayer)
    uuNode(2,4)=HydroParam%uxy(bbLayer,1,Nodes(4));     vvNode(2,4)=HydroParam%uxy(bbLayer,2,Nodes(4));     wwNode(2,4)=HydroParam%wfc(bbLayer,Nodes(4))
    uuNode(2,5)=HydroParam%ub(bbLayer,1,Nodes(5));     vvNode(2,5)=HydroParam%ub(bbLayer,2,Nodes(5));       wwNode(2,5)=HydroParam%ub(bbLayer,3,Nodes(5))
    uuNode(2,6)=HydroParam%uxy(bbLayer,1,Nodes(6));     vvNode(2,6)=HydroParam%uxy(bbLayer,2,Nodes(6));     wwNode(2,6)=HydroParam%wfc(bbLayer,Nodes(6))
    uuNode(3,4)=HydroParam%ug(Nodes(4),bbLayer+1);      vvNode(3,4)=HydroParam%vg(Nodes(4),bbLayer+1);      wwNode(3,4)=HydroParam%wg(Nodes(4),bbLayer+1)
    uuNode(3,5)=HydroParam%uxyL(bbLayer+1,1,Nodes(5));  vvNode(3,5)=HydroParam%uxyL(bbLayer+1,2,Nodes(5));  wwNode(3,5)=HydroParam%w(bbLayer+1,Nodes(5))
    uuNode(3,6)=HydroParam%ug(Nodes(6),bbLayer+1);      vvNode(3,6)=HydroParam%vg(Nodes(6),bbLayer+1);      wwNode(3,6)=HydroParam%wg(Nodes(6),bbLayer+1)
    ! Vertical 3 = East Edge
    uuNode(1,7)=HydroParam%uNode(bbLayer,1,Nodes(7));   vvNode(1,7)=HydroParam%uNode(bbLayer,2,Nodes(7));   wwNode(1,7)=HydroParam%uNode(bbLayer,3,Nodes(7))
    uuNode(1,8)=HydroParam%ug(Nodes(8),bbLayer);        vvNode(1,8)=HydroParam%vg(Nodes(8),bbLayer);        wwNode(1,8)=HydroParam%wg(Nodes(8),bbLayer)
    uuNode(1,9)=HydroParam%uNode(bbLayer,1,Nodes(9));   vvNode(1,9)=HydroParam%uNode(bbLayer,2,Nodes(9));   wwNode(1,9)=HydroParam%uNode(bbLayer,3,Nodes(9))
    uuNode(2,7)=HydroParam%ubV(bbLayer,1,Nodes(7));     vvNode(2,7)=HydroParam%ubV(bbLayer,2,Nodes(7));     wwNode(2,7)=HydroParam%ubV(bbLayer,3,Nodes(7))
    uuNode(2,8)=HydroParam%uxy(bbLayer,1,Nodes(8));     vvNode(2,8)=HydroParam%uxy(bbLayer,2,Nodes(8));     wwNode(2,8)=HydroParam%wfc(bbLayer,Nodes(8))
    uuNode(2,9)=HydroParam%ubV(bbLayer,1,Nodes(9));     vvNode(2,9)=HydroParam%ubV(bbLayer,2,Nodes(9));     wwNode(2,9)=HydroParam%ubV(bbLayer,3,Nodes(9))
    uuNode(3,7)=HydroParam%uNode(bbLayer+1,1,Nodes(7)); vvNode(3,7)=HydroParam%uNode(bbLayer+1,2,Nodes(7)); wwNode(3,7)=HydroParam%uNode(bbLayer+1,3,Nodes(7))
    uuNode(3,8)=HydroParam%ug(Nodes(8),bbLayer+1);      vvNode(3,8)=HydroParam%vg(Nodes(8),bbLayer+1);      wwNode(3,8)=HydroParam%wg(Nodes(8),bbLayer+1)
    uuNode(3,9)=HydroParam%uNode(bbLayer+1,1,Nodes(9)); vvNode(3,9)=HydroParam%uNode(bbLayer+1,2,Nodes(9)); wwNode(3,9)=HydroParam%uNode(bbLayer+1,3,Nodes(9))
    
    ! uuNode(Layer,Node Position), Layer = [k-1/2, k , k +1/2] == [1,2,3]
    uuNodet(1,1)=HydroParam%uNodet(bbLayer,1,Nodes(1));   vvNodet(1,1)=HydroParam%uNodet(bbLayer,2,Nodes(1));   wwNodet(1,1)=HydroParam%uNodet(bbLayer,3,Nodes(1))
    uuNodet(1,2)=HydroParam%ugt(Nodes(2),bbLayer);        vvNodet(1,2)=HydroParam%vgt(Nodes(2),bbLayer);        wwNodet(1,2)=HydroParam%wgt(Nodes(2),bbLayer)
    uuNodet(1,3)=HydroParam%uNodet(bbLayer,1,Nodes(3));   vvNodet(1,3)=HydroParam%uNodet(bbLayer,2,Nodes(3));   wwNodet(1,3)=HydroParam%uNodet(bbLayer,3,Nodes(3))
    uuNodet(2,1)=HydroParam%ubVt(bbLayer,1,Nodes(1));     vvNodet(2,1)=HydroParam%ubVt(bbLayer,2,Nodes(1));     wwNodet(2,1)=HydroParam%ubVt(bbLayer,3,Nodes(1))
    uuNodet(2,2)=HydroParam%uxyt(bbLayer,1,Nodes(2));     vvNodet(2,2)=HydroParam%uxyt(bbLayer,2,Nodes(2));     wwNodet(2,2)=HydroParam%wfct(bbLayer,Nodes(2))
    uuNodet(2,3)=HydroParam%ubVt(bbLayer,1,Nodes(3));     vvNodet(2,3)=HydroParam%ubVt(bbLayer,2,Nodes(3));     wwNodet(2,3)=HydroParam%ubVt(bbLayer,3,Nodes(3))
    uuNodet(3,1)=HydroParam%uNodet(bbLayer+1,1,Nodes(1)); vvNodet(3,1)=HydroParam%uNodet(bbLayer+1,2,Nodes(1)); wwNodet(3,1)=HydroParam%uNodet(bbLayer+1,3,Nodes(1))
    uuNodet(3,2)=HydroParam%ugt(Nodes(2),bbLayer+1);      vvNodet(3,2)=HydroParam%vgt(Nodes(2),bbLayer+1);      wwNodet(3,2)=HydroParam%wgt(Nodes(2),bbLayer+1)
    uuNodet(3,3)=HydroParam%uNodet(bbLayer+1,1,Nodes(3)); vvNodet(3,3)=HydroParam%uNodet(bbLayer+1,2,Nodes(3)); wwNodet(3,3)=HydroParam%uNodet(bbLayer+1,3,Nodes(3))
    ! Vertical 2 = Cell Centered Section
    uuNodet(1,4)=HydroParam%ugt(Nodes(4),bbLayer);      vvNodet(1,4)=HydroParam%vgt(Nodes(4),bbLayer);      	  wwNodet(1,4)=HydroParam%wgt(Nodes(4),bbLayer)
    uuNodet(1,5)=HydroParam%uxyLt(bbLayer,1,Nodes(5));  vvNodet(1,5)=HydroParam%uxyLt(bbLayer,2,Nodes(5));  	  wwNodet(1,5)=HydroParam%wt(bbLayer,Nodes(5))
    uuNodet(1,6)=HydroParam%ugt(Nodes(6),bbLayer);      vvNodet(1,6)=HydroParam%vgt(Nodes(6),bbLayer);      	  wwNodet(1,6)=HydroParam%wgt(Nodes(6),bbLayer)
    uuNodet(2,4)=HydroParam%uxyt(bbLayer,1,Nodes(4));     vvNodet(2,4)=HydroParam%uxyt(bbLayer,2,Nodes(4));     wwNodet(2,4)=HydroParam%wfct(bbLayer,Nodes(4))
    uuNodet(2,5)=HydroParam%ubt(bbLayer,1,Nodes(5));     vvNodet(2,5)=HydroParam%ubt(bbLayer,2,Nodes(5));       wwNodet(2,5)=HydroParam%ubt(bbLayer,3,Nodes(5))
    uuNodet(2,6)=HydroParam%uxyt(bbLayer,1,Nodes(6));     vvNodet(2,6)=HydroParam%uxyt(bbLayer,2,Nodes(6));     wwNodet(2,6)=HydroParam%wfct(bbLayer,Nodes(6))
    uuNodet(3,4)=HydroParam%ugt(Nodes(4),bbLayer+1);      vvNodet(3,4)=HydroParam%vgt(Nodes(4),bbLayer+1);      wwNodet(3,4)=HydroParam%wgt(Nodes(4),bbLayer+1)
    uuNodet(3,5)=HydroParam%uxyLt(bbLayer+1,1,Nodes(5));  vvNodet(3,5)=HydroParam%uxyLt(bbLayer+1,2,Nodes(5));  wwNodet(3,5)=HydroParam%wt(bbLayer+1,Nodes(5))
    uuNodet(3,6)=HydroParam%ugt(Nodes(6),bbLayer+1);      vvNodet(3,6)=HydroParam%vgt(Nodes(6),bbLayer+1);      wwNodet(3,6)=HydroParam%wgt(Nodes(6),bbLayer+1)
    ! Vertical 3 = East Edge
    uuNodet(1,7)=HydroParam%uNode(bbLayer,1,Nodes(7));   vvNodet(1,7)=HydroParam%uNodet(bbLayer,2,Nodes(7));   wwNodet(1,7)=HydroParam%uNodet(bbLayer,3,Nodes(7))
    uuNodet(1,8)=HydroParam%ug(Nodes(8),bbLayer);        vvNodet(1,8)=HydroParam%vgt(Nodes(8),bbLayer);        wwNodet(1,8)=HydroParam%wgt(Nodes(8),bbLayer)
    uuNodet(1,9)=HydroParam%uNode(bbLayer,1,Nodes(9));   vvNodet(1,9)=HydroParam%uNodet(bbLayer,2,Nodes(9));   wwNodet(1,9)=HydroParam%uNodet(bbLayer,3,Nodes(9))
    uuNodet(2,7)=HydroParam%ubV(bbLayer,1,Nodes(7));     vvNodet(2,7)=HydroParam%ubVt(bbLayer,2,Nodes(7));     wwNodet(2,7)=HydroParam%ubVt(bbLayer,3,Nodes(7))
    uuNodet(2,8)=HydroParam%uxy(bbLayer,1,Nodes(8));     vvNodet(2,8)=HydroParam%uxyt(bbLayer,2,Nodes(8));     wwNodet(2,8)=HydroParam%wfct(bbLayer,Nodes(8))
    uuNodet(2,9)=HydroParam%ubV(bbLayer,1,Nodes(9));     vvNodet(2,9)=HydroParam%ubVt(bbLayer,2,Nodes(9));     wwNodet(2,9)=HydroParam%ubVt(bbLayer,3,Nodes(9))
    uuNodet(3,7)=HydroParam%uNode(bbLayer+1,1,Nodes(7)); vvNodet(3,7)=HydroParam%uNodet(bbLayer+1,2,Nodes(7)); wwNodet(3,7)=HydroParam%uNodet(bbLayer+1,3,Nodes(7))
    uuNodet(3,8)=HydroParam%ug(Nodes(8),bbLayer+1);      vvNodet(3,8)=HydroParam%vgt(Nodes(8),bbLayer+1);      wwNodet(3,8)=HydroParam%wgt(Nodes(8),bbLayer+1)
    uuNodet(3,9)=HydroParam%uNode(bbLayer+1,1,Nodes(9)); vvNodet(3,9)=HydroParam%uNodet(bbLayer+1,2,Nodes(9)); wwNodet(3,9)=HydroParam%uNodet(bbLayer+1,3,Nodes(9))
    !  
    !If (bbLayer==HydroParam%ElCapitalM(bbElem)) Then
    !    !xxNode(Layer,Node Position), Layer = [k-1/2, k , k +1/2] == [1,2,3]
    !    xxNode(1,1) = MeshParam%xNode(Nodes(1));      yyNode(1,1) = MeshParam%yNode(Nodes(1));      zzNode(1,1)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,3) = MeshParam%xNode(Nodes(3));      yyNode(1,3) = MeshParam%yNode(Nodes(3));      zzNode(1,3)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(2,1) = MeshParam%xNode(Nodes(1));      yyNode(2,1) = MeshParam%yNode(Nodes(1));      zzNode(2,1)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%peta(Nodes(1)))*0.5d0
    !    xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=Max((HydroParam%Ze(bbLayer,bbElem) + HydroParam%Z(bbLayer+1,Nodes(2)))*0.5d0,HydroParam%Ze(bbLayer,bbElem) + NearZero/2)
    !    xxNode(2,3) = MeshParam%xNode(Nodes(3));      yyNode(2,3) = MeshParam%yNode(Nodes(3));      zzNode(2,3)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%peta(Nodes(3)))*0.5d0
    !    xxNode(3,1) = MeshParam%xNode(Nodes(1));      yyNode(3,1) = MeshParam%yNode(Nodes(1));      zzNode(3,1)=HydroParam%peta(Nodes(1))
    !    !xxNode(3,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(3,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(3,2)=HydroParam%Z(bbLayer+1,Nodes(2))
    !    xxNode(3,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(3,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(3,2)= Max(HydroParam%Z(bbLayer+1,Nodes(2)),zzNode(2,2) + NearZero)
    !    xxNode(3,3) = MeshParam%xNode(Nodes(3));      yyNode(3,3) = MeshParam%yNode(Nodes(3));      zzNode(3,3)=HydroParam%peta(Nodes(3))
    !    
    !    xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,5) = MeshParam%xb(Nodes(5));         yyNode(1,5) = MeshParam%yb(Nodes(5));         zzNode(1,5)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(1,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(1,6)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));  zzNode(2,4)=Max((HydroParam%Ze(bbLayer,bbElem) + HydroParam%Z(bbLayer+1,Nodes(4)))*0.5d0, HydroParam%Ze(bbLayer,bbElem) + NearZero/2)
    !    xxNode(2,5) = MeshParam%xb(Nodes(5));         yyNode(2,5) = MeshParam%yb(Nodes(5));          zzNode(2,5)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%eta(Nodes(5)))*0.5d0
    !    xxNode(2,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(2,6) = MeshParam%EdgeBary(2,Nodes(6));  zzNode(2,6)=Max((HydroParam%Ze(bbLayer,bbElem) + HydroParam%Z(bbLayer+1,Nodes(6)))*0.5d0, HydroParam%Ze(bbLayer,bbElem) + NearZero/2)
    !    !xxNode(3,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(3,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(3,4)=HydroParam%Z(bbLayer+1,Nodes(4))
    !    xxNode(3,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(3,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(3,4)=Max(HydroParam%Z(bbLayer+1,Nodes(4)),zzNode(2,4) + NearZero)
    !    xxNode(3,5) = MeshParam%xb(Nodes(5));         yyNode(3,5) = MeshParam%yb(Nodes(5));         zzNode(3,5)=HydroParam%eta(Nodes(5))
    !    !xxNode(3,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(3,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(3,6)=HydroParam%Z(bbLayer+1,Nodes(6))
    !    xxNode(3,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(3,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(3,6)=Max(HydroParam%Z(bbLayer+1,Nodes(6)),zzNode(2,6) + NearZero)
    !    
    !    xxNode(1,7) = MeshParam%xNode(Nodes(7));      yyNode(1,7) = MeshParam%yNode(Nodes(7));      zzNode(1,7)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(1,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(1,8)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,9) = MeshParam%xNode(Nodes(9));      yyNode(1,9) = MeshParam%yNode(Nodes(9));      zzNode(1,9)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(2,7) = MeshParam%xNode(Nodes(7));      yyNode(2,7) = MeshParam%yNode(Nodes(7));      zzNode(2,7)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%peta(Nodes(7)))*0.5d0 
    !    xxNode(2,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(2,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(2,8)= Max( (HydroParam%Ze(bbLayer,bbElem) + HydroParam%Z(bbLayer+1,Nodes(8)) )*0.5d0, HydroParam%Ze(bbLayer,bbElem) + NearZero/2)
    !    xxNode(2,9) = MeshParam%xNode(Nodes(9));      yyNode(2,9) = MeshParam%yNode(Nodes(9));      zzNode(2,9)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%peta(Nodes(9)))*0.5d0
    !    xxNode(3,7) = MeshParam%xNode(Nodes(7));      yyNode(3,7) = MeshParam%yNode(Nodes(7));      zzNode(3,7)=HydroParam%peta(Nodes(7))
    !    !xxNode(3,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(3,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(3,8)=HydroParam%Z(bbLayer+1,Nodes(8))
    !    xxNode(3,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(3,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(3,8)= Max(HydroParam%Z(bbLayer+1,Nodes(8)),zzNode(2,8) + NearZero)
    !    xxNode(3,9) = MeshParam%xNode(Nodes(9));      yyNode(3,9) = MeshParam%yNode(Nodes(9));      zzNode(3,9)=HydroParam%peta(Nodes(9))      
    !Else
    !    xxNode(1,1) = MeshParam%xNode(Nodes(1));      yyNode(1,1) = MeshParam%yNode(Nodes(1));      zzNode(1,1)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=HydroParam%Ze(bbLayer,bbElem)   
    !    xxNode(1,3) = MeshParam%xNode(Nodes(3));      yyNode(1,3) = MeshParam%yNode(Nodes(3));      zzNode(1,3)=HydroParam%Ze(bbLayer,bbElem)    
    !    xxNode(2,1) = MeshParam%xNode(Nodes(1));      yyNode(2,1) = MeshParam%yNode(Nodes(1));      zzNode(2,1)=HydroParam%Zb(bbLayer,bbElem)    
    !    xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=HydroParam%Zb(bbLayer,bbElem)
    !    xxNode(2,3) = MeshParam%xNode(Nodes(3));      yyNode(2,3) = MeshParam%yNode(Nodes(3));      zzNode(2,3)=HydroParam%Zb(bbLayer,bbElem)  
    !    xxNode(3,1) = MeshParam%xNode(Nodes(1));      yyNode(3,1) = MeshParam%yNode(Nodes(1));      zzNode(3,1)=HydroParam%Ze(bbLayer+1,bbElem)
    !    xxNode(3,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(3,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(3,2)=HydroParam%Ze(bbLayer+1,bbElem)    
    !    xxNode(3,3) = MeshParam%xNode(Nodes(3));      yyNode(3,3) = MeshParam%yNode(Nodes(3));      zzNode(3,3)=HydroParam%Ze(bbLayer+1,bbElem)     
    !    
    !    xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4)=HydroParam%Ze(bbLayer,bbElem)   
    !    xxNode(1,5) = MeshParam%xb(Nodes(5));         yyNode(1,5) = MeshParam%yb(Nodes(5));         zzNode(1,5)=HydroParam%Ze(bbLayer,bbElem)     
    !    xxNode(1,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(1,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(1,6)=HydroParam%Ze(bbLayer,bbElem)   
    !    xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));  zzNode(2,4)=HydroParam%Zb(bbLayer,bbElem)
    !    xxNode(2,5) = MeshParam%xb(Nodes(5));         yyNode(2,5) = MeshParam%yb(Nodes(5));          zzNode(2,5)=HydroParam%Zb(bbLayer,bbElem)
    !    xxNode(2,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(2,6) = MeshParam%EdgeBary(2,Nodes(6));  zzNode(2,6)=HydroParam%Zb(bbLayer,bbElem)
    !    xxNode(3,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(3,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(3,4)=HydroParam%Ze(bbLayer+1,bbElem)   
    !    xxNode(3,5) = MeshParam%xb(Nodes(5));         yyNode(3,5) = MeshParam%yb(Nodes(5));         zzNode(3,5)=HydroParam%Ze(bbLayer+1,bbElem)    
    !    xxNode(3,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(3,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(3,6)=HydroParam%Ze(bbLayer+1,bbElem)    
    !    
    !    xxNode(1,7) = MeshParam%xNode(Nodes(7));      yyNode(1,7) = MeshParam%yNode(Nodes(7));      zzNode(1,7)=HydroParam%Ze(bbLayer,bbElem)   
    !    xxNode(1,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(1,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(1,8)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,9) = MeshParam%xNode(Nodes(9));      yyNode(1,9) = MeshParam%yNode(Nodes(9));      zzNode(1,9)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(2,7) = MeshParam%xNode(Nodes(7));      yyNode(2,7) = MeshParam%yNode(Nodes(7));      zzNode(2,7)=HydroParam%Zb(bbLayer,bbElem)    
    !    xxNode(2,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(2,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(2,8)=HydroParam%Zb(bbLayer,bbElem)  
    !    xxNode(2,9) = MeshParam%xNode(Nodes(9));      yyNode(2,9) = MeshParam%yNode(Nodes(9));      zzNode(2,9)=HydroParam%Zb(bbLayer,bbElem)
    !    xxNode(3,7) = MeshParam%xNode(Nodes(7));      yyNode(3,7) = MeshParam%yNode(Nodes(7));      zzNode(3,7)=HydroParam%Ze(bbLayer+1,bbElem)  
    !    xxNode(3,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(3,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(3,8)=HydroParam%Ze(bbLayer+1,bbElem)    
    !    xxNode(3,9) = MeshParam%xNode(Nodes(9));      yyNode(3,9) = MeshParam%yNode(Nodes(9));      zzNode(3,9)=HydroParam%Ze(bbLayer+1,bbElem)  
    !EndIf
    !        
    If (bbLayer==HydroParam%ElCapitalM(bbElem)) Then
        !xxNode(Layer,Node Position), Layer = [k-1/2, k , k +1/2] == [1,2,3]
        xxNode(1,1) = MeshParam%xNode(Nodes(1));      yyNode(1,1) = MeshParam%yNode(Nodes(1));      zzNode(1,1) = HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
        xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2) = HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
        xxNode(1,3) = MeshParam%xNode(Nodes(3));      yyNode(1,3) = MeshParam%yNode(Nodes(3));      zzNode(1,3) = HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
        xxNode(2,1) = MeshParam%xNode(Nodes(1));      yyNode(2,1) = MeshParam%yNode(Nodes(1));      zzNode(2,1) = Max((HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%peta(Nodes(1)))*0.5d0,HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)))
        xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2) = Max((HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%Z(bbLayer+1,Nodes(2)))*0.5d0,HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)))
        xxNode(2,3) = MeshParam%xNode(Nodes(3));      yyNode(2,3) = MeshParam%yNode(Nodes(3));      zzNode(2,3) = Max((HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%peta(Nodes(3)))*0.5d0,HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)))
        xxNode(3,1) = MeshParam%xNode(Nodes(1));      yyNode(3,1) = MeshParam%yNode(Nodes(1));      zzNode(3,1) = Max(HydroParam%peta(Nodes(1)), zzNode(2,1))
        xxNode(3,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(3,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(3,2) = Max(HydroParam%Z(bbLayer+1,Nodes(2)),zzNode(2,2))
        xxNode(3,3) = MeshParam%xNode(Nodes(3));      yyNode(3,3) = MeshParam%yNode(Nodes(3));      zzNode(3,3) = Max(HydroParam%peta(Nodes(3)), zzNode(2,3))
        
        xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4) = HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
        xxNode(1,5) = MeshParam%xb(Nodes(5));         yyNode(1,5) = MeshParam%yb(Nodes(5));         zzNode(1,5) = HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
        xxNode(1,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(1,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(1,6) = HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
        xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(2,4) = Max((HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%Z(bbLayer+1,Nodes(4)))*0.5d0, HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)))
        xxNode(2,5) = MeshParam%xb(Nodes(5));         yyNode(2,5) = MeshParam%yb(Nodes(5));         zzNode(2,5) = Max((HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%eta(Nodes(5)))*0.5d0, HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)))
        xxNode(2,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(2,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(2,6) = Max((HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%Z(bbLayer+1,Nodes(6)))*0.5d0, HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)))
        xxNode(3,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(3,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(3,4) = Max(HydroParam%Z(bbLayer+1,Nodes(4)),zzNode(2,4))
        xxNode(3,5) = MeshParam%xb(Nodes(5));         yyNode(3,5) = MeshParam%yb(Nodes(5));         zzNode(3,5) = Max(HydroParam%eta(Nodes(5)),zzNode(2,5))
        xxNode(3,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(3,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(3,6) = Max(HydroParam%Z(bbLayer+1,Nodes(6)),zzNode(2,6))
        
        xxNode(1,7) = MeshParam%xNode(Nodes(7));      yyNode(1,7) = MeshParam%yNode(Nodes(7));      zzNode(1,7) = HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
        xxNode(1,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(1,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(1,8) = HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
        xxNode(1,9) = MeshParam%xNode(Nodes(9));      yyNode(1,9) = MeshParam%yNode(Nodes(9));      zzNode(1,9) = HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
        xxNode(2,7) = MeshParam%xNode(Nodes(7));      yyNode(2,7) = MeshParam%yNode(Nodes(7));      zzNode(2,7) = Max((HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%peta(Nodes(7)))*0.5d0, zzNode(1,7)) 
        xxNode(2,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(2,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(2,8) = Max((HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%Z(bbLayer+1,Nodes(8)) )*0.5d0, HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)))
        xxNode(2,9) = MeshParam%xNode(Nodes(9));      yyNode(2,9) = MeshParam%yNode(Nodes(9));      zzNode(2,9) = Max((HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%peta(Nodes(9)))*0.5d0, zzNode(1,9))
        xxNode(3,7) = MeshParam%xNode(Nodes(7));      yyNode(3,7) = MeshParam%yNode(Nodes(7));      zzNode(3,7) = Max(HydroParam%peta(Nodes(7)), zzNode(2,7)) 
        xxNode(3,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(3,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(3,8) = Max(HydroParam%Z(bbLayer+1,Nodes(8)), zzNode(2,8))
        xxNode(3,9) = MeshParam%xNode(Nodes(9));      yyNode(3,9) = MeshParam%yNode(Nodes(9));      zzNode(3,9) = Max(HydroParam%peta(Nodes(9)), zzNode(2,9)) 
    Else
        xxNode(1,1) = MeshParam%xNode(Nodes(1));      yyNode(1,1) = MeshParam%yNode(Nodes(1));      zzNode(1,1)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
        xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
        xxNode(1,3) = MeshParam%xNode(Nodes(3));      yyNode(1,3) = MeshParam%yNode(Nodes(3));      zzNode(1,3)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
        xxNode(2,1) = MeshParam%xNode(Nodes(1));      yyNode(2,1) = MeshParam%yNode(Nodes(1));      zzNode(2,1)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2     
        xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2   
        xxNode(2,3) = MeshParam%xNode(Nodes(3));      yyNode(2,3) = MeshParam%yNode(Nodes(3));      zzNode(2,3)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2   
        xxNode(3,1) = MeshParam%xNode(Nodes(1));      yyNode(3,1) = MeshParam%yNode(Nodes(1));      zzNode(3,1)=HydroParam%Ze(bbLayer+1,bbElem)
        xxNode(3,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(3,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(3,2)=HydroParam%Ze(bbLayer+1,bbElem)    
        xxNode(3,3) = MeshParam%xNode(Nodes(3));      yyNode(3,3) = MeshParam%yNode(Nodes(3));      zzNode(3,3)=HydroParam%Ze(bbLayer+1,bbElem)     
        
        xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
        xxNode(1,5) = MeshParam%xb(Nodes(5));         yyNode(1,5) = MeshParam%yb(Nodes(5));         zzNode(1,5)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
        xxNode(1,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(1,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(1,6)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
        xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));  zzNode(2,4)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2    
        xxNode(2,5) = MeshParam%xb(Nodes(5));         yyNode(2,5) = MeshParam%yb(Nodes(5));          zzNode(2,5)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2  
        xxNode(2,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(2,6) = MeshParam%EdgeBary(2,Nodes(6));  zzNode(2,6)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2   
        xxNode(3,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(3,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(3,4)=HydroParam%Ze(bbLayer+1,bbElem)   
        xxNode(3,5) = MeshParam%xb(Nodes(5));         yyNode(3,5) = MeshParam%yb(Nodes(5));         zzNode(3,5)=HydroParam%Ze(bbLayer+1,bbElem)    
        xxNode(3,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(3,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(3,6)=HydroParam%Ze(bbLayer+1,bbElem)    
        
        xxNode(1,7) = MeshParam%xNode(Nodes(7));      yyNode(1,7) = MeshParam%yNode(Nodes(7));      zzNode(1,7)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
        xxNode(1,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(1,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(1,8)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
        xxNode(1,9) = MeshParam%xNode(Nodes(9));      yyNode(1,9) = MeshParam%yNode(Nodes(9));      zzNode(1,9)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
        xxNode(2,7) = MeshParam%xNode(Nodes(7));      yyNode(2,7) = MeshParam%yNode(Nodes(7));      zzNode(2,7)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2     
        xxNode(2,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(2,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(2,8)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2   
        xxNode(2,9) = MeshParam%xNode(Nodes(9));      yyNode(2,9) = MeshParam%yNode(Nodes(9));      zzNode(2,9)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2   
        xxNode(3,7) = MeshParam%xNode(Nodes(7));      yyNode(3,7) = MeshParam%yNode(Nodes(7));      zzNode(3,7)=HydroParam%Ze(bbLayer+1,bbElem)  
        xxNode(3,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(3,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(3,8)=HydroParam%Ze(bbLayer+1,bbElem)    
        xxNode(3,9) = MeshParam%xNode(Nodes(9));      yyNode(3,9) = MeshParam%yNode(Nodes(9));      zzNode(3,9)=HydroParam%Ze(bbLayer+1,bbElem)  
    EndIf
            
    return
    End Subroutine iQuadraticNodes  
    
    EndModule ELM