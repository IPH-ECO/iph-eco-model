    !  
    !Subroutine ELMConservative4(uuBtrack, vvBtrack, wwBtrack, bbElem, bbLayer, iLayer, iEdge, xt, yt, zt, dt, psi_flag, HydroParam, MeshParam)
    !   
    !Use MeshVars !, Only: 
    !Use Hydrodynamic ! Only:
    !
    !Implicit none
    !
    !Real :: xt, yt, zt, dt
    !Integer:: bbElem, bbLayer, psi_flag, iEdge, iLayer
    !Real,intent(inout) ::  uuBtrack, vvBtrack, wwBtrack
    !Real :: uuNode(3,9), vvNode(3,9), wwNode(3,9), xxNode(3,9), yyNode(3,9), zzNode(3,9)
    !Real :: uuNodeBT(9), vvNodeBT(9), wwNodeBT(9), xxNodeBT(9), yyNodeBT(9), zzNodeBT(9), hNodeBT(9)
    !Real :: Yuu(3), Yvv(3), Yww(3), Xuu(3), Xvv(3), Xww(3)
    !Real :: uuint, vvint, wwint
    !Integer :: rElem, uElem, dElem, lElem, nElem, ElFlag
    !type(MeshGridParam) :: MeshParam
    !type(HydrodynamicParam) :: HydroParam
    !Real :: fi_small = 0.0d0
    !Real :: epsGrad = 10  !CAYO
    !Real :: rj(2,4), psi(2,4), ru, soma
    !Real :: Courant
    !
    !psi(:,:) = 0.0d0
    !rj(:,:) = 0.0d0
    !HydroParam%uArrow(iLayer,:,iEdge) = 0.d0
    !Courant = uuNode(3)*dt/MeshParam%CirDistance(iEdge) !7 - Super-C/ 8 - Ultimate-Quickest/ 9 - Hyper-C
    !    
    !If(IEdge==568) Then
    !    continue
    !EndIf
    !
    !! For nodes positions in the vector, the Standard is this:
    !!
    !!               .----.----.
    !!               |      n9 | 
    !!               |       x |
    !!               |       | |
    !!               |     n8| |
    !!     .----.----.----.--x-.----.----. 
    !!     |         |       | |         |
    !!     |       n2|    n3 | |n4       |
    !!     | n1 x----x-------x-x----x n5 | 
    !!     |         |       | |         |
    !!     .----.----.----.--x-.----.----.
    !!               |     n7| |
    !!               |       | |         
    !!               |       x | 
    !!               |      n6 |
    !!               .----.----.
    !
    !! 1 - Neighbour Elements:    
    !uElem = MeshParam%Right(MeshParam%Edge(1,bbElem))
    !lElem = MeshParam%Right(MeshParam%Edge(2,bbElem))
    !dElem = MeshParam%Right(MeshParam%Edge(3,bbElem))
    !rElem = MeshParam%Right(MeshParam%Edge(4,bbElem))
    !
    !! 2 - Nodes Velocities:
    !! 2.1 - Velocities for node n3:
    !Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), bbElem, bbLayer, HydroParam, MeshParam)
    !Call iQuadraticCons(uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), Xuu(:), Xvv(:), Xww(:), Yuu(:), Yvv(:), Yww(:), xt, yt, zt)
    !uuNodeBT(3) = uuBtrack; vvNodeBT(3) = vvBtrack; wwNodeBT(3) = wwBtrack; xxNodeBT(3) = xt; yyNodeBT(3) = yt; zzNodeBT(3) = zt; hNodeBT(3) = Max(HydroParam%eta(bbElem)-sum(HydroParam%DZsi(:,bbElem)),0.d0)
    !
    !!2.2 - Velocities for node n2:
    !nElem = bbElem 
    !Call pointInElem(nElem, lElem, xt - MeshParam%dx/2, yt, MeshParam)
    !! If nElem =/ lElem, this implies that lElem =/ 0 (condition checked in function pointInElem)
    !If (nElem == lElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(lElem)) Then
    !        !u(2) == u(iEdge==2,bbElem)
    !        nElem = bbElem
    !        uuNodeBT(2) = HydroParam%uxy(bbLayer,1,MeshParam%Edge(2,bbElem)); vvNodeBT(2) = HydroParam%uxy(bbLayer,2,MeshParam%Edge(2,bbElem)); wwNodeBT(2) = HydroParam%wfc(bbLayer,MeshParam%Edge(2,bbElem)); xxNodeBT(2) = xxNodeBT(3) - MeshParam%dx/2; yyNodeBT(2) = yyNodeBT(3); zzNodeBT(2) = zt; hNodeBT(2) = Max(HydroParam%H(MeshParam%Edge(2,bbElem))-sum(HydroParam%DZsj(:,MeshParam%Edge(2,bbElem))),0.d0)          
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt - MeshParam%dx/2, yt, zt )
    !        uuNodeBT(2) = uuBtrack; vvNodeBT(2) = vvBtrack; wwNodeBT(2) = wwBtrack; xxNodeBT(2) =  xt - MeshParam%dx/2; yyNodeBT(2) = yt; zzNodeBT(2) = zt; hNodeBT(2) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)    
    !    EndIf
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt - MeshParam%dx/2, yt, zt )
    !    uuNodeBT(2) = uuBtrack; vvNodeBT(2) = vvBtrack; wwNodeBT(2) = wwBtrack; xxNodeBT(2) =  xt - MeshParam%dx/2; yyNodeBT(2) = yt; zzNodeBT(2) = zt; hNodeBT(2) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf
    !
    !!2.3 - Velocities for node n1:
    !If (nElem == lElem) Then
    !    lElem = MeshParam%Right(MeshParam%Edge(2,lElem))
    !    Call pointInElem(nElem, lElem, xt - 3*MeshParam%dx/2, yt, MeshParam)   
    !Else
    !    Call pointInElem(nElem, lElem, xt - 3*MeshParam%dx/2, yt, MeshParam)
    !EndIf    
    !
    !If (nElem == lElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(lElem)) Then
    !        ! n1 == n2:
    !        uuNodeBT(1) = uuNodeBT(2); vvNodeBT(1) = uuNodeBT(2); wwNodeBT(1) = uuNodeBT(2); xxNodeBT(1) = xxNodeBT(2) - MeshParam%dx/2; yyNodeBT(1) = yyNodeBT(2); zzNodeBT(1) = zt; hNodeBT(1) = hNodeBT(2)
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt - 3*MeshParam%dx/2, yt, zt )
    !        uuNodeBT(1) = uuBtrack; vvNodeBT(1) = vvBtrack; wwNodeBT(1) = wwBtrack; xxNodeBT(1) =  xt - 3*MeshParam%dx/2; yyNodeBT(1) = yt; zzNodeBT(1) = zt; hNodeBT(1) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)    
    !    EndIf
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt - 3*MeshParam%dx/2, yt, zt )
    !    uuNodeBT(1) = uuBtrack; vvNodeBT(1) = vvBtrack; wwNodeBT(1) = wwBtrack; xxNodeBT(1) =  xt - 3*MeshParam%dx/2; yyNodeBT(1) = yt; zzNodeBT(1) = zt; hNodeBT(1) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf
    !
    !!2.4 - Velocities for node n4:
    !nElem = bbElem 
    !Call pointInElem(nElem, rElem, xt + MeshParam%dx/2, yt, MeshParam)
    !! If nElem == rElem, this implies that lElem =/ 0 (condition checked in function pointInElem)
    !If (nElem == rElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(rElem)) Then
    !        !u(4) == u(iEdge==4,bbElem)
    !        nElem = bbElem
    !        uuNodeBT(4) = HydroParam%uxy(bbLayer,1,MeshParam%Edge(4,bbElem)); vvNodeBT(4) = HydroParam%uxy(bbLayer,2,MeshParam%Edge(4,bbElem)); wwNodeBT(4) = HydroParam%wfc(bbLayer,MeshParam%Edge(4,bbElem)); xxNodeBT(4) = xxNodeBT(3) + MeshParam%dx/2; yyNodeBT(4) = yyNodeBT(3); zzNodeBT(4) = zt; hNodeBT(4) = Max(HydroParam%H(MeshParam%Edge(4,bbElem))-sum(HydroParam%DZsj(:,MeshParam%Edge(4,bbElem))),0.d0)
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt + MeshParam%dx/2, yt, zt )
    !        uuNodeBT(4) = uuBtrack; vvNodeBT(4) = vvBtrack; wwNodeBT(4) = wwBtrack; xxNodeBT(4) =  xt + MeshParam%dx/2; yyNodeBT(4) = yt; zzNodeBT(4) = zt; hNodeBT(4) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !    EndIf        
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt + MeshParam%dx/2, yt, zt )
    !    uuNodeBT(4) = uuBtrack; vvNodeBT(4) = vvBtrack; wwNodeBT(4) = wwBtrack; xxNodeBT(4) =  xt + MeshParam%dx/2; yyNodeBT(4) = yt; zzNodeBT(4) = zt; hNodeBT(4) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf
    ! 
    !!Energy conservation x - direction:
    !If ((uuNode(3) - uuNode(1))/MeshParam%CirDistance(iEdge) > epsGrad > 0) Then    
    !    If (uuNode(3) + uuNode(1) >= 0) Then  
    !        HydroParam%uArrow(iLayer,1,iEdge) = 0.5*(uuNode(3) + uuNode(1))
    !
    !        rj(1,1) = ru(uuNode(3),uuNode(2),uuNode(1),xxNode(3),xxNode(2),xxNode(1),yyNode(3),yyNode(2),yyNode(1))
    !        Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
    !            
    !        rj(1,2) = ru(uuNode(4),uuNode(3),uuNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
    !                
    !        rj(2,1) = ru(vvNode(3),vvNode(2),vvNode(1),xxNode(3),xxNode(2),xxNode(1),yyNode(3),yyNode(2),yyNode(1))
    !        Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
    !            
    !        rj(2,2) = ru(vvNode(4),vvNode(3),vvNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2)) 
    !    Else
    !        
    !        !2.5 - Velocities for node n5:
    !        If (nElem == rElem) Then
    !            rElem = MeshParam%Right(MeshParam%Edge(2,rElem))
    !            Call pointInElem(nElem, lElem, xt + 3*MeshParam%dx/2, yt, MeshParam)    
    !        Else
    !            Call pointInElem(nElem, lElem, xt + 3*MeshParam%dx/2, yt, MeshParam)  
    !        EndIf    
    !
    !        ! If nElem == rElem, this implies that lElem =/ 0 (condition checked in function pointInElem)
    !        If (nElem == rElem) Then
    !            If (bbLayer > HydroParam%ElSmallms(rElem)) Then
    !                ! n5 == n4:
    !                uuNodeBT(5) = uuNodeBT(4); vvNodeBT(5) = uuNodeBT(4); wwNodeBT(5) = uuNodeBT(4); xxNodeBT(5) = xxNodeBT(4) + MeshParam%dx/2; yyNodeBT(5) = yyNodeBT(4); zzNodeBT(5) = zt; hNodeBT(5) = hNodeBT(4)
    !            Else
    !                Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !                Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt + 3*MeshParam%dx/2, yt, zt )
    !                uuNodeBT(5) = uuBtrack; vvNodeBT(5) = vvBtrack; wwNodeBT(5) = wwBtrack; xxNodeBT(5) =  xt + 3*MeshParam%dx/2; yyNodeBT(5) = yt; zzNodeBT(5) = zt; hNodeBT(5) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)          
    !            EndIf
    !        Else
    !            Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !            Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt + 3*MeshParam%dx/2, yt, zt )
    !            uuNodeBT(5) = uuBtrack; vvNodeBT(5) = vvBtrack; wwNodeBT(5) = wwBtrack; xxNodeBT(5) =  xt + 3*MeshParam%dx/2; yyNodeBT(5) = yt; zzNodeBT(5) = zt; hNodeBT(5) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !        EndIf 
    !
    !        HydroParam%uArrow(iLayer,1,iEdge) = 0.5*(uuNode(3) + uuNode(5))
    !         
    !        rj(1,1) = ru(uuNode(4),uuNode(3),uuNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
    !            
    !        rj(1,2) = ru(uuNode(5),uuNode(4),uuNode(3),xxNode(5),xxNode(4),xxNode(3),yyNode(5),yyNode(4),yyNode(3))
    !        Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
    !                
    !        rj(2,1) = ru(vvNode(4),vvNode(3),vvNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
    !            
    !        rj(2,2) = ru(vvNode(5),vvNode(4),vvNode(3),xxNode(5),xxNode(4),xxNode(3),yyNode(5),yyNode(4),yyNode(3))
    !        Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))               
    !    EndIf
    !Else
    !    ! Momentum Conservation x - direction:
    !    If (uuNode(3)*dzNode(3) + uuNode(1)*dzNode(1) >= 0 ) Then
    !        !HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(3)*dzNode(3) + uuNode(2)*dzNode(2))/(dzNode(4) + dzNode(2))                                
    !        HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(3)*dzNode(3) + uuNode(1)*dzNode(1))/(dzNode(4) + dzNode(2)) 
    !        
    !        rj(1,1) = ru(uuNode(3),uuNode(2),uuNode(1),xxNode(3),xxNode(2),xxNode(1),yyNode(3),yyNode(2),yyNode(1))
    !        Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
    !            
    !        rj(1,2) = ru(uuNode(4),uuNode(3),uuNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
    !                
    !        rj(2,1) = ru(vvNode(3),vvNode(2),vvNode(1),xxNode(3),xxNode(2),xxNode(1),yyNode(3),yyNode(2),yyNode(1))
    !        Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
    !            
    !        rj(2,2) = ru(vvNode(4),vvNode(3),vvNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))                                                      
    !    Else          
    !        !2.5 - Velocities for node n5:
    !        If (nElem == rElem) Then
    !            rElem = MeshParam%Right(MeshParam%Edge(2,rElem))
    !            Call pointInElem(nElem, lElem, xt + 3*MeshParam%dx/2, yt, MeshParam)    
    !        Else
    !            Call pointInElem(nElem, lElem, xt + 3*MeshParam%dx/2, yt, MeshParam)  
    !        EndIf    
    !
    !        ! If nElem == rElem, this implies that lElem =/ 0 (condition checked in function pointInElem)
    !        If (nElem == rElem) Then
    !            If (bbLayer > HydroParam%ElSmallms(rElem)) Then
    !                ! n5 == n4:
    !                uuNodeBT(5) = uuNodeBT(4); vvNodeBT(5) = uuNodeBT(4); wwNodeBT(5) = uuNodeBT(4); xxNodeBT(5) = xxNodeBT(4) + MeshParam%dx/2; yyNodeBT(5) = yyNodeBT(4); zzNodeBT(5) = zt; hNodeBT(5) = hNodeBT(4)
    !            Else
    !                Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !                Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt + 3*MeshParam%dx/2, yt, zt )
    !                uuNodeBT(5) = uuBtrack; vvNodeBT(5) = vvBtrack; wwNodeBT(5) = wwBtrack; xxNodeBT(5) =  xt + 3*MeshParam%dx/2; yyNodeBT(5) = yt; zzNodeBT(5) = zt; hNodeBT(5) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)          
    !            EndIf
    !        Else
    !            Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !            Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt + 3*MeshParam%dx/2, yt, zt )
    !            uuNodeBT(5) = uuBtrack; vvNodeBT(5) = vvBtrack; wwNodeBT(5) = wwBtrack; xxNodeBT(5) =  xt + 3*MeshParam%dx/2; yyNodeBT(5) = yt; zzNodeBT(5) = zt; hNodeBT(5) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !        EndIf 
    ! 
    !        !HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(4)*dzNode(4) + uuNode(3)*dzNode(3))/(dzNode(4) + dzNode(2))                                
    !        HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(3)*dzNode(3) + uuNode(5)*dzNode(5))/(dzNode(4) + dzNode(2))
    !        
    !        rj(1,1) = ru(uuNode(4),uuNode(3),uuNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
    !            
    !        rj(1,2) = ru(uuNode(5),uuNode(4),uuNode(3),xxNode(5),xxNode(4),xxNode(3),yyNode(5),yyNode(4),yyNode(3))
    !        Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
    !                
    !        rj(2,1) = ru(vvNode(4),vvNode(3),vvNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
    !            
    !        rj(2,2) = ru(vvNode(5),vvNode(4),vvNode(3),xxNode(5),xxNode(4),xxNode(3),yyNode(5),yyNode(4),yyNode(3))
    !        Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))   
    !    EndIf        
    !    
    !EndIf
    !
    !
    !! 2.6 - Velocities for node n7: 
    !nElem = bbElem 
    !Call pointInElem(nElem, dElem, xt, yt - MeshParam%dy/2, MeshParam)
    !
    !If (nElem == dElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(dElem)) Then
    !        !u(7) == u(iEdge == 3, bbElem)
    !        nElem = bbElem
    !        uuNodeBT(7) = HydroParam%uxy(bbLayer,1,MeshParam%Edge(3,bbElem)); vvNodeBT(7) = HydroParam%uxy(bbLayer,2,MeshParam%Edge(3,bbElem)); wwNodeBT(7) =  HydroParam%wfc(bbLayer,MeshParam%Edge(3,bbElem)); xxNodeBT(7) = xxNodeBT(3); yyNodeBT(7) = yyNodeBT(3)  - MeshParam%dy/2; zzNodeBT(7) = zt; hNodeBT(7) = Max(HydroParam%H(MeshParam%Edge(3,bbElem))-sum(HydroParam%DZsj(:,MeshParam%Edge(3,bbElem))),0.d0)
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt - MeshParam%dy/2, zt )
    !        uuNodeBT(7) = uuBtrack; vvNodeBT(7) = vvBtrack; wwNodeBT(7) = wwBtrack; xxNodeBT(7) =  xt; yyNodeBT(7) = yt - MeshParam%dy/2; zzNodeBT(7) = zt; hNodeBT(7) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !    EndIf
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt - MeshParam%dy/2, zt )
    !    uuNodeBT(7) = uuBtrack; vvNodeBT(7) = vvBtrack; wwNodeBT(7) = wwBtrack; xxNodeBT(7) =  xt; yyNodeBT(7) = yt - MeshParam%dy/2; zzNodeBT(7) = zt; hNodeBT(7) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf
    !
    !!2.7 - Velocities for node n6:
    !If (nElem == dElem) Then
    !    dElem = MeshParam%Right(MeshParam%Edge(3,dElem))
    !    Call pointInElem(nElem, dElem, xt, yt - 3*MeshParam%dy/2, MeshParam)
    !Else
    !    Call pointInElem(nElem, dElem, xt, yt - 3*MeshParam%dy/2, MeshParam)
    !EndIf    
    !
    !! If nElem == dElem, this implies that dElem =/ 0 (condition checked in function pointInElem)
    !If (nElem == dElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(dElem)) Then
    !        ! n6 == n7:
    !        uuNodeBT(6) = uuNodeBT(7); vvNodeBT(6) = vvNodeBT(7); wwNodeBT(6) = wwNodeBT(7); xxNodeBT(6) = xxNodeBT(7); yyNodeBT(6) = yyNodeBT(7) - MeshParam%dy/2; zzNodeBT(6) = zt; hNodeBT(6) = hNodeBT(7)
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt  - 3*MeshParam%dy/2, zt )
    !        uuNodeBT(6) = uuBtrack; vvNodeBT(6) = vvBtrack; wwNodeBT(6) = wwBtrack; xxNodeBT(6) =  xt; yyNodeBT(6) = yt  - 3*MeshParam%dy/2; zzNodeBT(6) = zt; hNodeBT(6) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)  
    !    EndIf        
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt  - 3*MeshParam%dy/2, zt )
    !    uuNodeBT(6) = uuBtrack; vvNodeBT(6) = vvBtrack; wwNodeBT(6) = wwBtrack; xxNodeBT(6) =  xt; yyNodeBT(6) = yt  - 3*MeshParam%dy/2; zzNodeBT(6) = zt; hNodeBT(6) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf
    !    
    !! 2.8 - Velocities for node n8: 
    !nElem = bbElem 
    !Call pointInElem(nElem, uElem, xt, yt + MeshParam%dy/2, MeshParam)
    !! If nElem == uElem, this implies that uElem =/ 0 (condition checked in function pointInElem)
    !If (nElem == uElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(uElem)) Then
    !        !u(8) == u(iEdge == 1, bbElem)
    !        nElem = bbElem        
    !        uuNodeBT(8) = HydroParam%uxy(bbLayer,1,MeshParam%Edge(1,bbElem)); vvNodeBT(8) = HydroParam%uxy(bbLayer,2,MeshParam%Edge(1,bbElem)); wwNodeBT(8) =  HydroParam%wfc(bbLayer,MeshParam%Edge(1,bbElem)); xxNodeBT(8) = xxNodeBT(3); yyNodeBT(8) = yyNodeBT(3) + MeshParam%dy/2; zzNodeBT(8) = zt; hNodeBT(8) = Max(HydroParam%H(MeshParam%Edge(1,bbElem))-sum(HydroParam%DZsj(:,MeshParam%Edge(1,bbElem))),0.d0)
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt + MeshParam%dy/2, zt )
    !        uuNodeBT(8) = uuBtrack; vvNodeBT(8) = vvBtrack; wwNodeBT(8) = wwBtrack; xxNodeBT(8) =  xt; yyNodeBT(8) = yt + MeshParam%dy/2; zzNodeBT(8) = zt; hNodeBT(8) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !    EndIf
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt + MeshParam%dy/2, zt )
    !    uuNodeBT(8) = uuBtrack; vvNodeBT(8) = vvBtrack; wwNodeBT(8) = wwBtrack; xxNodeBT(8) =  xt; yyNodeBT(8) = yt + MeshParam%dy/2; zzNodeBT(8) = zt; hNodeBT(8) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf
    !
    !!Energy conservation y-direction:
    !If ((vvNode(3) - vvNode(6))/MeshParam%CirDistance(iEdge) > epsGrad > 0) Then
    !        
    !    If (vvNode(3) + vvNode(6) >= 0) Then
    !        HydroParam%uArrow(iLayer,2,iEdge) =  0.5*(vvNode(3) + vvNode(6))
    !        
    !        rj(1,3) = ru(uuNode(3),uuNode(7),uuNode(6),xxNode(3),xxNode(7),xxNode(6),yyNode(3),yyNode(7),yyNode(6))
    !        Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
    !            
    !        rj(1,4) = ru(uuNode(8),uuNode(3),uuNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                     
    !                             
    !        rj(2,3) = ru(vvNode(3),vvNode(7),vvNode(6),xxNode(3),xxNode(7),xxNode(6),yyNode(3),yyNode(7),yyNode(6))
    !        Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
    !            
    !        rj(2,4) = ru(vvNode(8),vvNode(3),vvNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))               
    !    Else
    !
    !        !2.9 - Velocities for node n9:
    !        If (nElem == uElem) Then
    !            uElem = MeshParam%Right(MeshParam%Edge(1,uElem))
    !            Call pointInElem(nElem, uElem, xt, yt + 3*MeshParam%dy/2, MeshParam)
    !        Else
    !            Call pointInElem(nElem, uElem, xt, yt + 3*MeshParam%dy/2, MeshParam)
    !        EndIf    
    !
    !        ! If nElem == uElem, this implies that uElem =/ 0 (condition checked in function pointInElem)
    !        If (nElem == uElem) Then
    !            If (bbLayer > HydroParam%ElSmallms(uElem)) Then
    !                ! n9 == n8:
    !                uuNodeBT(9) = uuNodeBT(8); vvNodeBT(9) = vvNodeBT(8); wwNodeBT(9) = wwNodeBT(8); xxNodeBT(9) = xxNodeBT(8); yyNodeBT(9) = yyNodeBT(8) + MeshParam%dy/2; zzNodeBT(9) = zt; hNodeBT(9) = hNodeBT(8)
    !            Else
    !                Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !                Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt + 3*MeshParam%dy/2, zt )
    !                uuNodeBT(9) = uuBtrack; vvNodeBT(9) = vvBtrack; wwNodeBT(9) = wwBtrack; xxNodeBT(9) =  xt; yyNodeBT(9) = yt + 3*MeshParam%dy/2; zzNodeBT(9) = zt; hNodeBT(9) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !            EndIf
    !        Else
    !            Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !            Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt + 3*MeshParam%dy/2, zt )
    !            uuNodeBT(9) = uuBtrack; vvNodeBT(9) = vvBtrack; wwNodeBT(9) = wwBtrack; xxNodeBT(9) =  xt; yyNodeBT(9) = yt + 3*MeshParam%dy/2; zzNodeBT(9) = zt; hNodeBT(9) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !        EndIf
    !        
    !        HydroParam%uArrow(iLayer,2,iEdge) =  0.5*(vvNode(3) + vvNode(9))
    !       
    !        rj(1,3) = ru(uuNode(8),uuNode(3),uuNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
    !            
    !        rj(1,4) = ru(uuNode(9),uuNode(8),uuNode(3),xxNode(9),xxNode(8),xxNode(3),yyNode(9),yyNode(8),yyNode(3))
    !        Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                    
    !                
    !        rj(2,3) = ru(vvNode(8),vvNode(3),vvNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
    !
    !        rj(2,4) = ru(vvNode(9),vvNode(8),vvNode(3),xxNode(9),xxNode(8),xxNode(3),yyNode(9),yyNode(8),yyNode(3))
    !        Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))                                                    
    !    EndIf
    !Else
    !    ! Momentum Conservation y - direction:
    !    If (vvNode(3)*dzNode(3) + vvNode(6)*dzNode(6)>= 0) Then
    !        !HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(3)*dzNode(3) + vvNode(7)*dzNode(7))/(dzNode(8) + dzNode(7))
    !        HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(3)*dzNode(3) + vvNode(6)*dzNode(6))/(dzNode(8) + dzNode(7))
    !        
    !        rj(1,3) = ru(uuNode(3),uuNode(7),uuNode(6),xxNode(3),xxNode(7),xxNode(6),yyNode(3),yyNode(7),yyNode(6))
    !        Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
    !            
    !        rj(1,4) = ru(uuNode(8),uuNode(3),uuNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                     
    !                             
    !        rj(2,3) = ru(vvNode(3),vvNode(7),vvNode(6),xxNode(3),xxNode(7),xxNode(6),yyNode(3),yyNode(7),yyNode(6))
    !        Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
    !            
    !        rj(2,4) = ru(vvNode(8),vvNode(3),vvNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))         
    !    Else
    !        !2.9 - Velocities for node n9:
    !        If (nElem == uElem) Then
    !            uElem = MeshParam%Right(MeshParam%Edge(1,uElem))
    !            Call pointInElem(nElem, uElem, xt, yt + 3*MeshParam%dy/2, MeshParam)
    !        Else
    !            Call pointInElem(nElem, uElem, xt, yt + 3*MeshParam%dy/2, MeshParam)
    !        EndIf    
    !
    !        ! If nElem == uElem, this implies that uElem =/ 0 (condition checked in function pointInElem)
    !        If (nElem == uElem) Then
    !            If (bbLayer > HydroParam%ElSmallms(uElem)) Then
    !                ! n9 == n8:
    !                uuNodeBT(9) = uuNodeBT(8); vvNodeBT(9) = vvNodeBT(8); wwNodeBT(9) = wwNodeBT(8); xxNodeBT(9) = xxNodeBT(8); yyNodeBT(9) = yyNodeBT(8) + MeshParam%dy/2; zzNodeBT(9) = zt; hNodeBT(9) = hNodeBT(8)
    !            Else
    !                Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !                Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt + 3*MeshParam%dy/2, zt )
    !                uuNodeBT(9) = uuBtrack; vvNodeBT(9) = vvBtrack; wwNodeBT(9) = wwBtrack; xxNodeBT(9) =  xt; yyNodeBT(9) = yt + 3*MeshParam%dy/2; zzNodeBT(9) = zt; hNodeBT(9) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !            EndIf
    !        Else
    !            Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !            Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt + 3*MeshParam%dy/2, zt )
    !            uuNodeBT(9) = uuBtrack; vvNodeBT(9) = vvBtrack; wwNodeBT(9) = wwBtrack; xxNodeBT(9) =  xt; yyNodeBT(9) = yt + 3*MeshParam%dy/2; zzNodeBT(9) = zt; hNodeBT(9) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !        EndIf              
    !        
    !        !HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(7)*dzNode(7) + vvNode(3)*dzNode(3))/(dzNode(8) + dzNode(7))
    !        HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(3)*dzNode(3) + vvNode(9)*dzNode(9))/(dzNode(8) + dzNode(7))
    !        
    !        rj(1,3) = ru(uuNode(8),uuNode(3),uuNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
    !            
    !        rj(1,4) = ru(uuNode(9),uuNode(8),uuNode(3),xxNode(9),xxNode(8),xxNode(3),yyNode(9),yyNode(8),yyNode(3))
    !        Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                    
    !                
    !        rj(2,3) = ru(vvNode(8),vvNode(3),vvNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
    !
    !        rj(2,4) = ru(vvNode(9),vvNode(8),vvNode(3),xxNode(9),xxNode(8),xxNode(3),yyNode(9),yyNode(8),yyNode(3))
    !        Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))                                                      
    !    EndIf
    !EndIf
    !
    !uuint = HydroParam%uArrow(iLayer,1,iEdge)
    !vvint = HydroParam%uArrow(iLayer,2,iEdge)
    !wwint = wwNodeBT(3)
    !
    !return
    !End Subroutine ELMConservative4 
    !
    !Subroutine ELMConservative3(uuBtrack, vvBtrack, wwBtrack, bbElem, bbLayer, iLayer, iEdge, xt, yt, zt, dt, psi_flag, HydroParam, MeshParam)
    !   
    !Use MeshVars !, Only: 
    !Use Hydrodynamic ! Only:
    !
    !Implicit none
    !
    !Real :: xt, yt, zt, dt
    !Integer:: bbElem, bbLayer, psi_flag, iEdge, iLayer
    !Real,intent(inout) ::  uuBtrack, vvBtrack, wwBtrack
    !Real :: uuNode(3,9), vvNode(3,9), wwNode(3,9), xxNode(3,9), yyNode(3,9), zzNode(3,9)
    !Real :: uuNodeBT(9), vvNodeBT(9), wwNodeBT(9), xxNodeBT(9), yyNodeBT(9), zzNodeBT(9), hNodeBT(9)
    !Real :: Yuu(3), Yvv(3), Yww(3), Xuu(3), Xvv(3), Xww(3)
    !Real :: uuint, vvint, wwint
    !Integer :: rElem, uElem, dElem, lElem, nElem, ElFlag
    !type(MeshGridParam) :: MeshParam
    !type(HydrodynamicParam) :: HydroParam
    !
    !If(IEdge==568) Then
    !    continue
    !EndIf
    !
    !! For nodes positions in the vector, the Standard is this:
    !!
    !!               .----.----.
    !!               |      n9 | 
    !!               |       x |
    !!               |       | |
    !!               |     n8| |
    !!     .----.----.----.--x-.----.----. 
    !!     |         |       | |         |
    !!     |       n2|    n3 | |n4       |
    !!     | n1 x----x-------x-x----x n5 | 
    !!     |         |       | |         |
    !!     .----.----.----.--x-.----.----.
    !!               |     n7| |
    !!               |       | |         
    !!               |       x | 
    !!               |      n6 |
    !!               .----.----.
    !
    !! 1 - Neighbour Elements:    
    !uElem = MeshParam%Right(MeshParam%Edge(1,bbElem))
    !lElem = MeshParam%Right(MeshParam%Edge(2,bbElem))
    !dElem = MeshParam%Right(MeshParam%Edge(3,bbElem))
    !rElem = MeshParam%Right(MeshParam%Edge(4,bbElem))
    !
    !! 2 - Nodes Velocities:
    !! 2.1 - Velocities for node n3:
    !Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), bbElem, bbLayer, HydroParam, MeshParam)
    !Call iQuadraticCons(uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), Xuu(:), Xvv(:), Xww(:), Yuu(:), Yvv(:), Yww(:), xt, yt, zt)
    !uuNodeBT(3) = uuBtrack; vvNodeBT(3) = vvBtrack; wwNodeBT(3) = wwBtrack; xxNodeBT(3) = xt; yyNodeBT(3) = yt; zzNodeBT(3) = zt; hNodeBT(3) = Max(HydroParam%eta(bbElem)-sum(HydroParam%DZsi(:,bbElem)),0.d0)
    !
    !!2.2 - Velocities for node n2:
    !nElem = bbElem 
    !Call pointInElem(nElem, lElem, xt - MeshParam%dx/2, yt, MeshParam)
    !! If nElem =/ lElem, this implies that lElem =/ 0 (condition checked in function pointInElem)
    !If (nElem == lElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(lElem)) Then
    !        nElem = bbElem
    !        uuNodeBT(2) = HydroParam%uxyback(bbLayer,1,MeshParam%Edge(2,bbElem)); vvNodeBT(2) = HydroParam%uxyback(bbLayer,2,MeshParam%Edge(2,bbElem)); wwNodeBT(2) = HydroParam%wfc(bbLayer,MeshParam%Edge(2,bbElem)); xxNodeBT(2) = xxNodeBT(3); yyNodeBT(2) = yyNodeBT(3); zzNodeBT(2) = zt; hNodeBT(2) = hNodeBT(3)            
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt - MeshParam%dx/2, yt, zt )
    !        uuNodeBT(2) = uuBtrack; vvNodeBT(2) = vvBtrack; wwNodeBT(2) = wwBtrack; xxNodeBT(2) =  xt - MeshParam%dx/2; yyNodeBT(2) = yt; zzNodeBT(2) = zt; hNodeBT(2) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)    
    !    EndIf
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt - MeshParam%dx/2, yt, zt )
    !    uuNodeBT(2) = uuBtrack; vvNodeBT(2) = vvBtrack; wwNodeBT(2) = wwBtrack; xxNodeBT(2) =  xt - MeshParam%dx/2; yyNodeBT(2) = yt; zzNodeBT(2) = zt; hNodeBT(2) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf
    !
    !!2.3 - Velocities for node n1:
    !If (nElem == lElem) Then
    !    lElem = MeshParam%Right(MeshParam%Edge(2,lElem))
    !    Call pointInElem(nElem, lElem, xt - 3*MeshParam%dx/2, yt, MeshParam)   
    !Else
    !    Call pointInElem(nElem, lElem, xt - 3*MeshParam%dx/2, yt, MeshParam)
    !EndIf    
    !
    !If (nElem == lElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(lElem)) Then
    !        ! n1 == n2:
    !        uuNodeBT(1) = uuNodeBT(2); vvNodeBT(1) = uuNodeBT(2); wwNodeBT(1) = uuNodeBT(2); xxNodeBT(1) = xxNodeBT(2); yyNodeBT(1) = yyNodeBT(2); zzNodeBT(1) = zt; hNodeBT(1) = hNodeBT(2)
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt - 3*MeshParam%dx/2, yt, zt )
    !        uuNodeBT(1) = uuBtrack; vvNodeBT(1) = vvBtrack; wwNodeBT(1) = wwBtrack; xxNodeBT(1) =  xt - 3*MeshParam%dx/2; yyNodeBT(1) = yt; zzNodeBT(1) = zt; hNodeBT(1) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)    
    !    EndIf
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt - 3*MeshParam%dx/2, yt, zt )
    !    uuNodeBT(1) = uuBtrack; vvNodeBT(1) = vvBtrack; wwNodeBT(1) = wwBtrack; xxNodeBT(1) =  xt - 3*MeshParam%dx/2; yyNodeBT(1) = yt; zzNodeBT(1) = zt; hNodeBT(1) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf
    !
    !!2.4 - Velocities for node n4:
    !nElem = bbElem 
    !Call pointInElem(nElem, rElem, xt + MeshParam%dx/2, yt, MeshParam)
    !! If nElem == rElem, this implies that lElem =/ 0 (condition checked in function pointInElem)
    !If (nElem == rElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(rElem)) Then
    !        !n4 == n3:
    !        !u(3) == u(iEdge==4,bbElem)
    !        nElem = bbElem
    !        uuNodeBT(4) = HydroParam%uxyback(bbLayer,1,MeshParam%Edge(4,bbElem)); vvNodeBT(4) = HydroParam%uxyback(bbLayer,2,MeshParam%Edge(4,bbElem)); wwNodeBT(4) = HydroParam%wfc(bbLayer,MeshParam%Edge(4,bbElem)); xxNodeBT(4) = xxNodeBT(3); yyNodeBT(4) = yyNodeBT(3); zzNodeBT(4) = zt; hNodeBT(4) = hNodeBT(3)
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt + MeshParam%dx/2, yt, zt )
    !        uuNodeBT(4) = uuBtrack; vvNodeBT(4) = vvBtrack; wwNodeBT(4) = wwBtrack; xxNodeBT(4) =  xt + MeshParam%dx/2; yyNodeBT(4) = yt; zzNodeBT(4) = zt; hNodeBT(4) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !    EndIf        
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt + MeshParam%dx/2, yt, zt )
    !    uuNodeBT(4) = uuBtrack; vvNodeBT(4) = vvBtrack; wwNodeBT(4) = wwBtrack; xxNodeBT(4) =  xt + MeshParam%dx/2; yyNodeBT(4) = yt; zzNodeBT(4) = zt; hNodeBT(4) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf
    !
    !!2.5 - Velocities for node n5:
    !If (nElem == rElem) Then
    !    rElem = MeshParam%Right(MeshParam%Edge(2,rElem))
    !    Call pointInElem(nElem, lElem, xt + 3*MeshParam%dx/2, yt, MeshParam)    
    !Else
    !    Call pointInElem(nElem, lElem, xt + 3*MeshParam%dx/2, yt, MeshParam)  
    !EndIf    
    !
    !! If nElem == rElem, this implies that lElem =/ 0 (condition checked in function pointInElem)
    !If (nElem == rElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(rElem)) Then
    !        ! n5 == n4:
    !        uuNodeBT(5) = uuNodeBT(4); vvNodeBT(5) = uuNodeBT(4); wwNodeBT(5) = uuNodeBT(4); xxNodeBT(5) = xxNodeBT(4); yyNodeBT(5) = yyNodeBT(4); zzNodeBT(5) = zt; hNodeBT(5) = hNodeBT(4)
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt + 3*MeshParam%dx/2, yt, zt )
    !        uuNodeBT(5) = uuBtrack; vvNodeBT(5) = vvBtrack; wwNodeBT(5) = wwBtrack; xxNodeBT(5) =  xt + 3*MeshParam%dx/2; yyNodeBT(5) = yt; zzNodeBT(5) = zt; hNodeBT(5) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)          
    !    EndIf
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt + 3*MeshParam%dx/2, yt, zt )
    !    uuNodeBT(5) = uuBtrack; vvNodeBT(5) = vvBtrack; wwNodeBT(5) = wwBtrack; xxNodeBT(5) =  xt + 3*MeshParam%dx/2; yyNodeBT(5) = yt; zzNodeBT(5) = zt; hNodeBT(5) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf 
    !
    !! 2.6 - Velocities for node n7: 
    !nElem = bbElem 
    !Call pointInElem(nElem, dElem, xt, yt - MeshParam%dy/2, MeshParam)
    !
    !If (nElem == dElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(dElem)) Then
    !        !n7 == n3:
    !        !u(7) == u(iEdge == 3, bbElem)
    !        nElem = bbElem        
    !        uuNodeBT(7) = HydroParam%uxyback(bbLayer,1,MeshParam%Edge(3,bbElem)); vvNodeBT(7) = HydroParam%uxyback(bbLayer,2,MeshParam%Edge(3,bbElem)); wwNodeBT(7) =  HydroParam%wfc(bbLayer,MeshParam%Edge(3,bbElem)); xxNodeBT(7) = xxNodeBT(3); yyNodeBT(7) = yyNodeBT(3); zzNodeBT(7) = zt; hNodeBT(7) = hNodeBT(3)
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt - MeshParam%dy/2, zt )
    !        uuNodeBT(7) = uuBtrack; vvNodeBT(7) = vvBtrack; wwNodeBT(7) = wwBtrack; xxNodeBT(7) =  xt; yyNodeBT(7) = yt - MeshParam%dy/2; zzNodeBT(7) = zt; hNodeBT(7) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !    EndIf
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt - MeshParam%dy/2, zt )
    !    uuNodeBT(7) = uuBtrack; vvNodeBT(7) = vvBtrack; wwNodeBT(7) = wwBtrack; xxNodeBT(7) =  xt; yyNodeBT(7) = yt - MeshParam%dy/2; zzNodeBT(7) = zt; hNodeBT(7) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf
    !
    !!2.7 - Velocities for node n6:
    !If (nElem == dElem) Then
    !    dElem = MeshParam%Right(MeshParam%Edge(3,dElem))
    !    Call pointInElem(nElem, dElem, xt, yt - 3*MeshParam%dy/2, MeshParam)
    !Else
    !    Call pointInElem(nElem, dElem, xt, yt - 3*MeshParam%dy/2, MeshParam)
    !EndIf    
    !
    !! If nElem == dElem, this implies that dElem =/ 0 (condition checked in function pointInElem)
    !If (nElem == dElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(dElem)) Then
    !        ! n6 == n7:
    !        uuNodeBT(6) = uuNodeBT(7); vvNodeBT(6) = vvNodeBT(7); wwNodeBT(6) = wwNodeBT(7); xxNodeBT(6) = xxNodeBT(7); yyNodeBT(6) = yyNodeBT(7); zzNodeBT(6) = zt; hNodeBT(6) = hNodeBT(7)
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt  - 3*MeshParam%dy/2, zt )
    !        uuNodeBT(6) = uuBtrack; vvNodeBT(6) = vvBtrack; wwNodeBT(6) = wwBtrack; xxNodeBT(6) =  xt; yyNodeBT(6) = yt  - 3*MeshParam%dy/2; zzNodeBT(6) = zt; hNodeBT(6) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)  
    !    EndIf        
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt  - 3*MeshParam%dy/2, zt )
    !    uuNodeBT(6) = uuBtrack; vvNodeBT(6) = vvBtrack; wwNodeBT(6) = wwBtrack; xxNodeBT(6) =  xt; yyNodeBT(6) = yt  - 3*MeshParam%dy/2; zzNodeBT(6) = zt; hNodeBT(6) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf  
    !
    !! 2.8 - Velocities for node n8: 
    !nElem = bbElem 
    !Call pointInElem(nElem, uElem, xt, yt + MeshParam%dy/2, MeshParam)
    !! If nElem == uElem, this implies that uElem =/ 0 (condition checked in function pointInElem)
    !If (nElem == uElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(uElem)) Then
    !        !n7 == n3:
    !        !u(8) == u(iEdge == 1, bbElem)
    !        nElem = bbElem        
    !        uuNodeBT(8) = HydroParam%uxyback(bbLayer,1,MeshParam%Edge(1,bbElem)); vvNodeBT(8) = HydroParam%uxyback(bbLayer,2,MeshParam%Edge(1,bbElem)); wwNodeBT(8) =  HydroParam%wfc(bbLayer,MeshParam%Edge(1,bbElem)); xxNodeBT(8) = xxNodeBT(3); yyNodeBT(8) = yyNodeBT(3); zzNodeBT(8) = zt; hNodeBT(8) = hNodeBT(3)
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt + MeshParam%dy/2, zt )
    !        uuNodeBT(8) = uuBtrack; vvNodeBT(8) = vvBtrack; wwNodeBT(8) = wwBtrack; xxNodeBT(8) =  xt; yyNodeBT(8) = yt + MeshParam%dy/2; zzNodeBT(8) = zt; hNodeBT(8) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !    EndIf
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt + MeshParam%dy/2, zt )
    !    uuNodeBT(8) = uuBtrack; vvNodeBT(8) = vvBtrack; wwNodeBT(8) = wwBtrack; xxNodeBT(8) =  xt; yyNodeBT(8) = yt + MeshParam%dy/2; zzNodeBT(8) = zt; hNodeBT(8) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf
    !
    !!2.9 - Velocities for node n9:
    !If (nElem == uElem) Then
    !    uElem = MeshParam%Right(MeshParam%Edge(1,uElem))
    !    Call pointInElem(nElem, uElem, xt, yt + 3*MeshParam%dy/2, MeshParam)
    !Else
    !    Call pointInElem(nElem, uElem, xt, yt + 3*MeshParam%dy/2, MeshParam)
    !EndIf    
    !
    !! If nElem == uElem, this implies that uElem =/ 0 (condition checked in function pointInElem)
    !If (nElem == uElem) Then
    !    If (bbLayer > HydroParam%ElSmallms(uElem)) Then
    !        ! n9 == n8:
    !        uuNodeBT(9) = uuNodeBT(8); vvNodeBT(9) = vvNodeBT(8); wwNodeBT(9) = wwNodeBT(8); xxNodeBT(9) = xxNodeBT(8); yyNodeBT(9) = yyNodeBT(8); zzNodeBT(9) = zt; hNodeBT(9) = hNodeBT(8)
    !    Else
    !        Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !        Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt + 3*MeshParam%dy/2, zt )
    !        uuNodeBT(9) = uuBtrack; vvNodeBT(9) = vvBtrack; wwNodeBT(9) = wwBtrack; xxNodeBT(9) =  xt; yyNodeBT(9) = yt + 3*MeshParam%dy/2; zzNodeBT(9) = zt; hNodeBT(9) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !    EndIf
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), nElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadratic (uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), xt, yt + 3*MeshParam%dy/2, zt )
    !    uuNodeBT(9) = uuBtrack; vvNodeBT(9) = vvBtrack; wwNodeBT(9) = wwBtrack; xxNodeBT(9) =  xt; yyNodeBT(9) = yt + 3*MeshParam%dy/2; zzNodeBT(9) = zt; hNodeBT(9) = Max(HydroParam%eta(nElem)-sum(HydroParam%DZsi(:,nElem)),0.d0)
    !EndIf
    !
    !! 3 - Characteristic Line End's Velocities by Conservative Scheme: 
    !Call ConservativeScheme(uuBtrack, vvBtrack, uuNodeBT(:), vvNodeBT(:), xxNodeBT(:), yyNodeBT(:), hNodeBT(:), dt, psi_flag, iLayer, iEdge, HydroParam, MeshParam)
    !
    !uuint = uuBtrack
    !vvint = vvBtrack
    !wwint = wwNodeBT(3)
    !
    !return
    !End Subroutine ELMConservative3
    !    
    ! 
    !Subroutine ELMConservative2(uuBtrack, vvBtrack, wwBtrack, bbElem, bbLayer, iLayer, iEdge, xt, yt, zt, dt, psi_flag, HydroParam, MeshParam)
    !   
    !Use MeshVars !, Only: 
    !Use Hydrodynamic ! Only:
    !
    !Implicit none
    !
    !Real :: xt, yt, zt, dt
    !Integer:: bbElem, bbLayer, psi_flag, iEdge, iLayer
    !Real,intent(inout) ::  uuBtrack, vvBtrack, wwBtrack
    !Real :: uuNode(3,9), vvNode(3,9), wwNode(3,9), xxNode(3,9), yyNode(3,9), zzNode(3,9)
    !Real :: uuNodeBT(9), vvNodeBT(9), wwNodeBT(9), xxNodeBT(9), yyNodeBT(9), zzNodeBT(9), hNodeBT(9)
    !Real :: Yuu(3), Yvv(3), Yww(3), Xuu(3), Xvv(3), Xww(3)
    !Real :: uuint, vvint, wwint
    !Integer :: rElem, uElem, dElem, lElem
    !type(MeshGridParam) :: MeshParam
    !type(HydrodynamicParam) :: HydroParam
    !
    !If(IEdge==568) Then
    !    continue
    !EndIf
    !
    !!uuNodeBT(:) = 0.d0; vvNodeBT(:) = 0.d0; wwNodeBT(:) = 0.d0; xxNodeBT(:) = 0.d0; yyNodeBT(:) = 0.d0; zzNodeBT(:) = 0.d0; hNodeBT(:) = 0.d0
    !
    !! For nodes positions in the vector, the Standard is this:
    !!
    !!               .----.----.
    !!               |      n9 | 
    !!               |       x |
    !!               |       | |
    !!               |     n8| |
    !!     .----.----.----.--x-.----.----. 
    !!     |         |       | |         |
    !!     |       n2|    n3 | |n4       |
    !!     | n1 x----x-------x-x----x n5 | 
    !!     |         |       | |         |
    !!     .----.----.----.--x-.----.----.
    !!               |     n7| |
    !!               |       | |         
    !!               |       x | 
    !!               |      n6 |
    !!               .----.----.
    !
    !! 1 - Neighbour Elements:    
    !uElem = MeshParam%Right(MeshParam%Edge(1,bbElem))
    !lElem = MeshParam%Right(MeshParam%Edge(2,bbElem))
    !dElem = MeshParam%Right(MeshParam%Edge(3,bbElem))
    !rElem = MeshParam%Right(MeshParam%Edge(4,bbElem))
    !
    !! 2 - Nodes Velocities:
    !! 2.1 - Velocities for nodes n2, n3, n4, n7 and n8:
    !Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), bbElem, bbLayer, HydroParam, MeshParam)
    !Call iQuadraticCons(uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), Xuu(:), Xvv(:), Xww(:), Yuu(:), Yvv(:), Yww(:), xt, yt, zt)
    !
    !uuNodeBT(2) = Yuu(1); vvNodeBT(2) = Yvv(1); wwNodeBT(2) = Yww(1); xxNodeBT(2) = xxNode(1,1); yyNodeBT(2) = yt; zzNodeBT(2) = zt; hNodeBT(2) = Max(HydroParam%H(MeshParam%Edge(2,bbElem))-HydroParam%hj(MeshParam%Edge(2,bbElem)),0.d0)
    !uuNodeBT(3) = uuBtrack; vvNodeBT(3) = vvBtrack; wwNodeBT(3) = wwBtrack; xxNodeBT(3) = xt; yyNodeBT(3) = yt; zzNodeBT(3) = zt; hNodeBT(3) = Max(HydroParam%eta(bbElem)-sum(HydroParam%DZsi(:,bbElem)),0.d0)
    !uuNodeBT(4) = Yuu(3); vvNodeBT(4) = Yvv(3); wwNodeBT(4) = Yww(3); xxNodeBT(4) = xxNode(1,7); yyNodeBT(4) = yt; zzNodeBT(4) = zt; hNodeBT(4) = Max(HydroParam%H(MeshParam%Edge(4,bbElem))-HydroParam%hj(MeshParam%Edge(4,bbElem)),0.d0)
    !uuNodeBT(7) = Xuu(1); vvNodeBT(7) = Xvv(1); wwNodeBT(7) = Xww(1); xxNodeBT(7) = xt; yyNodeBT(7) = yyNode(1,4); zzNodeBT(7) = zt; hNodeBT(7) = Max(HydroParam%H(MeshParam%Edge(3,bbElem))-HydroParam%hj(MeshParam%Edge(3,bbElem)),0.d0)
    !uuNodeBT(8) = Xuu(3); vvNodeBT(8) = Xvv(3); wwNodeBT(8) = Xww(3); xxNodeBT(8) = xt; yyNodeBT(8) = yyNode(1,6); zzNodeBT(8) = zt; hNodeBT(8) = Max(HydroParam%H(MeshParam%Edge(1,bbElem))-HydroParam%hj(MeshParam%Edge(1,bbElem)),0.d0)
    !
    !! 2.2 - Velocities for node n1:
    !If(lElem == 0) Then        
    !    uuNodeBT(1) = uuNodeBT(2); vvNodeBT(1) = uuNodeBT(2); wwNodeBT(1) = uuNodeBT(2); xxNodeBT(1) = xxNodeBT(2); yyNodeBT(1) = yyNodeBT(2); zzNodeBT(1) = zt; hNodeBT(1) = hNodeBT(2)
    !ElseIf (bbLayer < HydroParam%ElSmallms(lElem)) Then
    !    uuNodeBT(1) = uuNodeBT(2); vvNodeBT(1) = uuNodeBT(2); wwNodeBT(1) = uuNodeBT(2); xxNodeBT(1) = xxNodeBT(2); yyNodeBT(1) = yyNodeBT(2); zzNodeBT(1) = zt; hNodeBT(1) = hNodeBT(2)
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), lElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadraticCons(uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), Yuu(:), Yvv(:), Yww(:), Xuu(:), Xvv(:), Xww(:), xxNode(1,5), yt, zt)
    !    uuNodeBT(1) = Yuu(2); vvNodeBT(1) = Yvv(2); wwNodeBT(1) = Yww(2); xxNodeBT(1) = xxNode(1,5); yyNodeBT(1) = yt; zzNodeBT(1) = zt;  hNodeBT(1) = Max(HydroParam%eta(lElem)-sum(HydroParam%DZsi(:,lElem)),0.d0)    
    !EndIf
    !
    !! 2.3 - Velocities for node n5:
    !If(rElem == 0) Then
    !    uuNodeBT(5) = uuNodeBT(4); vvNodeBT(5) = uuNodeBT(4); wwNodeBT(5) = uuNodeBT(4); xxNodeBT(5) = xxNodeBT(4); yyNodeBT(5) = yyNodeBT(4); zzNodeBT(5) = zt; hNodeBT(5) = hNodeBT(4)
    !ElseIf(bbLayer < HydroParam%ElSmallms(rElem)) Then
    !    uuNodeBT(5) = uuNodeBT(4); vvNodeBT(5) = uuNodeBT(4); wwNodeBT(5) = uuNodeBT(4); xxNodeBT(5) = xxNodeBT(4); yyNodeBT(5) = yyNodeBT(4); zzNodeBT(5) = zt; hNodeBT(5) = hNodeBT(4)        
    !Else        
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), rElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadraticCons(uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), Xuu(:), Xvv(:), Xww(:), Yuu(:), Yvv(:), Yww(:), xxNode(1,5), yt, zt)
    !    uuNodeBT(5) = Yuu(2); vvNodeBT(5) = Yvv(2); wwNodeBT(5) = Yww(2); xxNodeBT(5) = xxNode(1,5); yyNodeBT(5) = yt; zzNodeBT(5) = zt;  hNodeBT(5) = max(HydroParam%eta(rElem)-sum(HydroParam%DZsi(:,rElem)),0.d0)
    !EndIf
    !
    !! 2.4 - Velocities for node n6:
    !If (dElem == 0) Then
    !    uuNodeBT(6) = uuNodeBT(7); vvNodeBT(6) = vvNodeBT(7); wwNodeBT(6) = wwNodeBT(7); xxNodeBT(6) = xxNodeBT(7); yyNodeBT(6) = yyNodeBT(7); zzNodeBT(6) = zt; hNodeBT(6) = hNodeBT(7)
    !ElseIf (bbLayer < HydroParam%ElSmallms(dElem)) Then
    !    uuNodeBT(6) = uuNodeBT(7); vvNodeBT(6) = vvNodeBT(7); wwNodeBT(6) = wwNodeBT(7); xxNodeBT(6) = xxNodeBT(7); yyNodeBT(6) = yyNodeBT(7); zzNodeBT(6) = zt; hNodeBT(6) = hNodeBT(7)        
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), dElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadraticCons(uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), Xuu(:), Xvv(:), Xww(:), Yuu(:), Yvv(:), Yww(:), xt, yyNode(1,5), zt)
    !    uuNodeBT(6) = Xuu(2); vvNodeBT(6) = Xvv(2); wwNodeBT(6) = Xww(2); xxNodeBT(6) = xt; yyNodeBT(6) = yyNode(1,5); zzNodeBT(6) = zt;  hNodeBT(6) = max(HydroParam%eta(dElem)-sum(HydroParam%DZsi(:,dElem)),0.d0)    
    !EndIf
    !
    !! 2.5 - Velocities for node n9:
    !If (uElem == 0) Then
    !    uuNodeBT(9) = uuNodeBT(8); vvNodeBT(9) = vvNodeBT(8); wwNodeBT(9) = wwNodeBT(8); xxNodeBT(9) = xxNodeBT(8); yyNodeBT(9) = yyNodeBT(8); zzNodeBT(9) = zt; hNodeBT(9) = hNodeBT(8)
    !ElseIf (bbLayer < HydroParam%ElSmallms(uElem)) Then
    !    uuNodeBT(9) = uuNodeBT(8); vvNodeBT(9) = vvNodeBT(8); wwNodeBT(9) = wwNodeBT(8); xxNodeBT(9) = xxNodeBT(8); yyNodeBT(9) = yyNodeBT(8); zzNodeBT(9) = zt; hNodeBT(9) = hNodeBT(8)        
    !Else
    !    Call iQuadraticNodes(uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), uElem, bbLayer, HydroParam, MeshParam)
    !    Call iQuadraticCons(uuBtrack, vvBtrack, wwBtrack, uuNode(:,:), vvNode(:,:), wwNode(:,:), xxNode(:,:), yyNode(:,:), zzNode(:,:), Xuu(:), Xvv(:), Xww(:), Yuu(:), Yvv(:), Yww(:), xt, yyNode(1,5), zt)
    !    uuNodeBT(9) = Xuu(2); vvNodeBT(9) = Xvv(2); wwNodeBT(9) = Xww(2); xxNodeBT(9) = xt; yyNodeBT(9) = yyNode(1,5); zzNodeBT(9) = zt;  hNodeBT(9) = max(HydroParam%eta(uElem)-sum(HydroParam%DZsi(:,uElem)),0.d0)    
    !EndIf    
    !
    !! 3 - Characteristic Line End's Velocities by Conservative Scheme: 
    !Call ConservativeScheme(uuBtrack, vvBtrack, uuNodeBT(:), vvNodeBT(:), xxNodeBT(:), yyNodeBT(:), hNodeBT(:), dt, psi_flag, iLayer, iEdge, HydroParam, MeshParam)
    !
    !
    !uuint = uuBtrack
    !vvint = vvBtrack
    !wwint = wwNodeBT(3)
    !
    !
    !return
    !End Subroutine ELMConservative2
    !
    !
    !   
    !Subroutine ConservativeScheme(uuBtrack, vvBtrack, uuNode, vvNode, xxNode, yyNode, dzNode, dt, psi_flag, iLayer, iEdge, HydroParam, MeshParam)
    !
    !Use MeshVars
    !Use Hydrodynamic
    !
    !Implicit None
    !
    !Real, intent(in) :: uuNode(9), vvNode(9), xxNode(9), yyNode(9), dzNode(9), dt
    !Real, intent(inout) :: uuBtrack, vvBtrack
    !Real :: fi_small = 0.0d0
    !Real :: epsGrad = 10  !CAYO
    !Real :: rj(2,4), psi(2,4), ru, soma
    !Real :: Courant
    !Integer :: psi_flag, iEdge, iLayer
    !type(MeshGridParam) :: MeshParam
    !type(HydrodynamicParam) :: HydroParam
    !
    !psi(:,:) = 0.0d0
    !rj(:,:) = 0.0d0
    !HydroParam%uArrow(iLayer,:,iEdge) = 0.d0
    !Courant = uuNode(3)*dt/MeshParam%CirDistance(iEdge) !7 - Super-C/ 8 - Ultimate-Quickest/ 9 - Hyper-C
    !
    !if (iEdge==568) then
    !    continue
    !EndIf
    !
    !!x - direction:
    !If ((uuNode(3) - uuNode(1))/MeshParam%CirDistance(iEdge) > epsGrad > 0) Then
    !                          
    !    If (uuNode(3) + uuNode(1) >= 0) Then  
    !        HydroParam%uArrow(iLayer,1,iEdge) = 0.5*(uuNode(3) + uuNode(1))
    !          
    !        rj(1,1) = ru(uuNode(3),uuNode(2),uuNode(1),xxNode(3),xxNode(2),xxNode(1),yyNode(3),yyNode(2),yyNode(1))
    !        Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
    !            
    !        rj(1,2) = ru(uuNode(4),uuNode(3),uuNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
    !                
    !        rj(2,1) = ru(vvNode(3),vvNode(2),vvNode(1),xxNode(3),xxNode(2),xxNode(1),yyNode(3),yyNode(2),yyNode(1))
    !        Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
    !            
    !        rj(2,2) = ru(vvNode(4),vvNode(3),vvNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))                                                                                                       
    !    Else                                                
    !        HydroParam%uArrow(iLayer,1,iEdge) = 0.5*(uuNode(3) + uuNode(5))
    !         
    !        rj(1,1) = ru(uuNode(4),uuNode(3),uuNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
    !            
    !        rj(1,2) = ru(uuNode(5),uuNode(4),uuNode(3),xxNode(5),xxNode(4),xxNode(3),yyNode(5),yyNode(4),yyNode(3))
    !        Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
    !                
    !        rj(2,1) = ru(vvNode(4),vvNode(3),vvNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
    !            
    !        rj(2,2) = ru(vvNode(5),vvNode(4),vvNode(3),xxNode(5),xxNode(4),xxNode(3),yyNode(5),yyNode(4),yyNode(3))
    !        Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))                            
    !    EndIf                                                       
    !    ! y-direction:
    !    If (vvNode(3) + vvNode(6) >= 0) Then
    !        HydroParam%uArrow(iLayer,2,iEdge) =  0.5*(vvNode(3) + vvNode(6))
    !        
    !        rj(1,3) = ru(uuNode(3),uuNode(7),uuNode(6),xxNode(3),xxNode(7),xxNode(6),yyNode(3),yyNode(7),yyNode(6))
    !        Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
    !            
    !        rj(1,4) = ru(uuNode(8),uuNode(3),uuNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                     
    !                             
    !        rj(2,3) = ru(vvNode(3),vvNode(7),vvNode(6),xxNode(3),xxNode(7),xxNode(6),yyNode(3),yyNode(7),yyNode(6))
    !        Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
    !            
    !        rj(2,4) = ru(vvNode(8),vvNode(3),vvNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))               
    !    Else
    !        HydroParam%uArrow(iLayer,2,iEdge) =  0.5*(vvNode(3) + vvNode(9))
    !       
    !        rj(1,3) = ru(uuNode(8),uuNode(3),uuNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
    !            
    !        rj(1,4) = ru(uuNode(9),uuNode(8),uuNode(3),xxNode(9),xxNode(8),xxNode(3),yyNode(9),yyNode(8),yyNode(3))
    !        Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                    
    !                
    !        rj(2,3) = ru(vvNode(8),vvNode(3),vvNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
    !
    !        rj(2,4) = ru(vvNode(9),vvNode(8),vvNode(3),xxNode(9),xxNode(8),xxNode(3),yyNode(9),yyNode(8),yyNode(3))
    !        Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))                                                    
    !    EndIf                            
    !Else
    !    ! x-direction:
    !    If (uuNode(3)*dzNode(3) + uuNode(1)*dzNode(1) >= 0 ) Then
    !        !HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(3)*dzNode(3) + uuNode(2)*dzNode(2))/(dzNode(4) + dzNode(2))                                
    !        HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(3)*dzNode(3) + uuNode(1)*dzNode(1))/(dzNode(4) + dzNode(2)) 
    !        
    !        rj(1,1) = ru(uuNode(3),uuNode(2),uuNode(1),xxNode(3),xxNode(2),xxNode(1),yyNode(3),yyNode(2),yyNode(1))
    !        Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
    !            
    !        rj(1,2) = ru(uuNode(4),uuNode(3),uuNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
    !                
    !        rj(2,1) = ru(vvNode(3),vvNode(2),vvNode(1),xxNode(3),xxNode(2),xxNode(1),yyNode(3),yyNode(2),yyNode(1))
    !        Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
    !            
    !        rj(2,2) = ru(vvNode(4),vvNode(3),vvNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))                                                      
    !    Else
    !        !HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(4)*dzNode(4) + uuNode(3)*dzNode(3))/(dzNode(4) + dzNode(2))                                
    !        HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(3)*dzNode(3) + uuNode(5)*dzNode(5))/(dzNode(4) + dzNode(2))
    !        
    !        rj(1,1) = ru(uuNode(4),uuNode(3),uuNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
    !            
    !        rj(1,2) = ru(uuNode(5),uuNode(4),uuNode(3),xxNode(5),xxNode(4),xxNode(3),yyNode(5),yyNode(4),yyNode(3))
    !        Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
    !                
    !        rj(2,1) = ru(vvNode(4),vvNode(3),vvNode(2),xxNode(4),xxNode(3),xxNode(2),yyNode(4),yyNode(3),yyNode(2))
    !        Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
    !            
    !        rj(2,2) = ru(vvNode(5),vvNode(4),vvNode(3),xxNode(5),xxNode(4),xxNode(3),yyNode(5),yyNode(4),yyNode(3))
    !        Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))   
    !    EndIf
    !                                         
    !! y-direction:                                                                             
    !    If (vvNode(3)*dzNode(3) + vvNode(6)*dzNode(6)>= 0) Then
    !        !HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(3)*dzNode(3) + vvNode(7)*dzNode(7))/(dzNode(8) + dzNode(7))
    !        HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(3)*dzNode(3) + vvNode(6)*dzNode(6))/(dzNode(8) + dzNode(7))
    !        
    !        rj(1,3) = ru(uuNode(3),uuNode(7),uuNode(6),xxNode(3),xxNode(7),xxNode(6),yyNode(3),yyNode(7),yyNode(6))
    !        Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
    !            
    !        rj(1,4) = ru(uuNode(8),uuNode(3),uuNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                     
    !                             
    !        rj(2,3) = ru(vvNode(3),vvNode(7),vvNode(6),xxNode(3),xxNode(7),xxNode(6),yyNode(3),yyNode(7),yyNode(6))
    !        Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
    !            
    !        rj(2,4) = ru(vvNode(8),vvNode(3),vvNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))         
    !    Else
    !        !HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(7)*dzNode(7) + vvNode(3)*dzNode(3))/(dzNode(8) + dzNode(7))
    !        HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(3)*dzNode(3) + vvNode(9)*dzNode(9))/(dzNode(8) + dzNode(7))
    !        
    !        rj(1,3) = ru(uuNode(8),uuNode(3),uuNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
    !            
    !        rj(1,4) = ru(uuNode(9),uuNode(8),uuNode(3),xxNode(9),xxNode(8),xxNode(3),yyNode(9),yyNode(8),yyNode(3))
    !        Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                    
    !                
    !        rj(2,3) = ru(vvNode(8),vvNode(3),vvNode(7),xxNode(8),xxNode(3),xxNode(7),yyNode(8),yyNode(3),yyNode(7))
    !        Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
    !
    !        rj(2,4) = ru(vvNode(9),vvNode(8),vvNode(3),xxNode(9),xxNode(8),xxNode(3),yyNode(9),yyNode(8),yyNode(3))
    !        Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))                                                      
    !    EndIf                                                                                                            
    !EndIf
    !
    !
    !If (HydroParam%uArrow(iLayer,1,iEdge)*dt/(xxNode(4)-xxNode(2)) <0 .or. HydroParam%uArrow(iLayer,2,iEdge)*dt/(yyNode(8)-yyNode(7)) <0) Then
    !    soma = HydroParam%uArrow(iLayer,1,iEdge) + HydroParam%uArrow(iLayer,2,iEdge)
    !    continue
    !EndIf
    !!HydroParam%uxyback(iLayer,1,iEdge) = uuNode(3) - dt*(0.5*Max(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(1,2))*uuNode(3) - (2 + psi(1,1) + psi(1,2))*uuNode(2) + psi(1,1)*uuNode(1))/MeshParam%dx + 0.5*Min(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(1,2))*uuNode(4) - (2 + psi(1,1) + psi(1,2))*uuNode(3) + psi(1,1)*uuNode(2))/MeshParam%dx + 0.5*Max(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(1,4))*uuNode(3) - (2 + psi(1,3) + psi(1,4))*uuNode(7) + psi(1,3)*uuNode(6))/MeshParam%dy + 0.5*Min(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(1,4))*uuNode(8) - (2 + psi(1,3) + psi(1,4))*uuNode(3) + psi(1,3)*uuNode(7))/MeshParam%dy)
    !!HydroParam%uxyback(iLayer,2,iEdge) = uuNode(3) - dt*(0.5*Max(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(2,2))*vvNode(3) - (2 + psi(2,1) + psi(2,2))*vvNode(2) + psi(2,1)*vvNode(1))/MeshParam%dx + 0.5*Min(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(2,2))*vvNode(4) - (2 + psi(2,1) + psi(2,2))*vvNode(3) + psi(2,1)*vvNode(2))/MeshParam%dx  + 0.5*Max(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(2,4))*vvNode(3) - (2 + psi(2,3) + psi(2,4))*vvNode(7) + psi(2,3)*vvNode(6))/MeshParam%dy + 0.5*Min(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(2,4))*vvNode(8) - (2 + psi(2,3) + psi(2,4))*vvNode(3) + psi(2,3)*vvNode(7))/MeshParam%dy)                
    !!
    !!HydroParam%uxyback(iLayer,1,iEdge) = uuNode(3) - dt*(0.5*Max(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(1,2))*uuNode(3) - (2 + psi(1,1) + psi(1,2))*uuNode(2) + psi(1,1)*uuNode(1))/(xxNode(3)-xxNode(2)) + 0.5*Min(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(1,2))*uuNode(4) - (2 + psi(1,1) + psi(1,2))*uuNode(3) + psi(1,1)*uuNode(2))/(xxNode(4)-xxNode(3)) + 0.5*Max(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(1,4))*uuNode(3) - (2 + psi(1,3) + psi(1,4))*uuNode(7) + psi(1,3)*uuNode(6))/(yyNode(3)-yyNode(7)) + 0.5*Min(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(1,4))*uuNode(8) - (2 + psi(1,3) + psi(1,4))*uuNode(3) + psi(1,3)*uuNode(7))/(yyNode(8)-yyNode(3)))
    !!HydroParam%uxyback(iLayer,2,iEdge) = uuNode(3) - dt*(0.5*Max(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(2,2))*vvNode(3) - (2 + psi(2,1) + psi(2,2))*vvNode(2) + psi(2,1)*vvNode(1))/(xxNode(3)-xxNode(2)) + 0.5*Min(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(2,2))*vvNode(4) - (2 + psi(2,1) + psi(2,2))*vvNode(3) + psi(2,1)*vvNode(2))/(xxNode(4)-xxNode(3))  + 0.5*Max(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(2,4))*vvNode(3) - (2 + psi(2,3) + psi(2,4))*vvNode(7) + psi(2,3)*vvNode(6))/(yyNode(3)-yyNode(7)) + 0.5*Min(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(2,4))*vvNode(8) - (2 + psi(2,3) + psi(2,4))*vvNode(3) + psi(2,3)*vvNode(7))/(yyNode(8)-yyNode(3)))                
    !!
    !
    !uuBtrack = uuNode(3) - dt*(0.5*Max(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(1,2))*uuNode(3) - (2 + psi(1,1) + psi(1,2))*uuNode(2) + psi(1,1)*uuNode(1))/MeshParam%dx + 0.5*Min(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(1,2))*uuNode(4) - (2 + psi(1,1) + psi(1,2))*uuNode(3) + psi(1,1)*uuNode(2))/MeshParam%dx + 0.5*Max(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(1,4))*uuNode(3) - (2 + psi(1,3) + psi(1,4))*uuNode(7) + psi(1,3)*uuNode(6))/MeshParam%dy + 0.5*Min(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(1,4))*uuNode(8) - (2 + psi(1,3) + psi(1,4))*uuNode(3) + psi(1,3)*uuNode(7))/MeshParam%dy)
    !vvBtrack = vvNode(3) - dt*(0.5*Max(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(2,2))*vvNode(3) - (2 + psi(2,1) + psi(2,2))*vvNode(2) + psi(2,1)*vvNode(1))/MeshParam%dx + 0.5*Min(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(2,2))*vvNode(4) - (2 + psi(2,1) + psi(2,2))*vvNode(3) + psi(2,1)*vvNode(2))/MeshParam%dx + 0.5*Max(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(2,4))*vvNode(3) - (2 + psi(2,3) + psi(2,4))*vvNode(7) + psi(2,3)*vvNode(6))/MeshParam%dy + 0.5*Min(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(2,4))*vvNode(8) - (2 + psi(2,3) + psi(2,4))*vvNode(3) + psi(2,3)*vvNode(7))/MeshParam%dy)
    !
    !!uuBtrack = uuNode(3) - dt*(0.5*Max(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(1,2))*uuNode(3) - (2 + psi(1,1) + psi(1,2))*uuNode(2) + psi(1,1)*uuNode(1))/(xxNode(4)-xxNode(2)) + 0.5*Min(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(1,2))*uuNode(4) - (2 + psi(1,1) + psi(1,2))*uuNode(3) + psi(1,1)*uuNode(2))/(xxNode(4)-xxNode(2)) + 0.5*Max(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(1,4))*uuNode(3) - (2 + psi(1,3) + psi(1,4))*uuNode(7) + psi(1,3)*uuNode(6))/(yyNode(8)-yyNode(7)) + 0.5*Min(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(1,4))*uuNode(8) - (2 + psi(1,3) + psi(1,4))*uuNode(3) + psi(1,3)*uuNode(7))/(yyNode(8)-yyNode(7)))
    !!vvBtrack = vvNode(3) - dt*(0.5*Max(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(2,2))*vvNode(3) - (2 + psi(2,1) + psi(2,2))*vvNode(2) + psi(2,1)*vvNode(1))/(xxNode(4)-xxNode(2)) + 0.5*Min(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2 + psi(2,2))*vvNode(4) - (2 + psi(2,1) + psi(2,2))*vvNode(3) + psi(2,1)*vvNode(2))/(xxNode(4)-xxNode(2)) + 0.5*Max(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(2,4))*vvNode(3) - (2 + psi(2,3) + psi(2,4))*vvNode(7) + psi(2,3)*vvNode(6))/(yyNode(8)-yyNode(7)) + 0.5*Min(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2 + psi(2,4))*vvNode(8) - (2 + psi(2,3) + psi(2,4))*vvNode(3) + psi(2,3)*vvNode(7))/(yyNode(8)-yyNode(7)))
    !!
    !End Subroutine ConservativeScheme
    !
    
    !
    !Subroutine iQuadraticNodes(uuNode, vvNode, wwNode, xxNode, yyNode, zzNode, bbElem, bbLayer, HydroParam, MeshParam)
    !
    !Use MeshVars !, Only: 
    !Use Hydrodynamic ! Only:
    !
    !Implicit none
    !
    !Integer, intent(in) :: bbElem, bbLayer
    !Real, intent(inout) ::  uuNode(3,9), vvNode(3,9), wwNode(3,9), xxNode(3,9), yyNode(3,9), zzNode(3,9)
    !Real :: Nodes(9)
    !Integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9
    !type(MeshGridParam) :: MeshParam
    !type(HydrodynamicParam) :: HydroParam
    !        
    !!iQuadratic Interpolation Nodes in bbLayer from bbElem
    !! For nodes positions in the vector, the Standard is this:
    !!        n6
    !!  n3 .--.--. n9
    !!     |     |
    !!  n2 .  .  . n8
    !!     | n5  |      
    !!  n1 .--.--. n7
    !!        n4            
    !n1 = MeshParam%Quadri(3,bbElem) + 1
    !n2 = MeshParam%Edge(2,bbElem)
    !n3 = MeshParam%Quadri(2,bbElem) + 1
    !n4 = MeshParam%Edge(3,bbElem)
    !n5 = bbElem
    !n6 = MeshParam%Edge(1,bbElem)
    !n7 = MeshParam%Quadri(4,bbElem) + 1
    !n8 = MeshParam%Edge(4,bbElem)
    !n9 = MeshParam%Quadri(1,bbElem) + 1
    !Nodes(1:9)= (/n1, n2, n3, n4, n5, n6, n7, n8, n9 /)
    !! Velocities by vertical section and by layers, in  bottom to top (k-1/2 to k + 1/2). See Figure 3 in [2]):
    !! Vertical 1 = West Edge
    !! uuNode(Layer,Node Position), Layer = [k-1/2, k , k +1/2] == [1,2,3]
    !uuNode(1,1)=HydroParam%uNode(bbLayer,1,Nodes(1));   vvNode(1,1)=HydroParam%uNode(bbLayer,2,Nodes(1));   wwNode(1,1)=HydroParam%uNode(bbLayer,3,Nodes(1))
    !uuNode(1,2)=HydroParam%ug(Nodes(2),bbLayer);        vvNode(1,2)=HydroParam%vg(Nodes(2),bbLayer);        wwNode(1,2)=HydroParam%wg(Nodes(2),bbLayer)
    !uuNode(1,3)=HydroParam%uNode(bbLayer,1,Nodes(3));   vvNode(1,3)=HydroParam%uNode(bbLayer,2,Nodes(3));   wwNode(1,3)=HydroParam%uNode(bbLayer,3,Nodes(3))
    !uuNode(2,1)=HydroParam%ubV(bbLayer,1,Nodes(1));     vvNode(2,1)=HydroParam%ubV(bbLayer,2,Nodes(1));     wwNode(2,1)=HydroParam%ubV(bbLayer,3,Nodes(1))
    !uuNode(2,2)=HydroParam%uxy(bbLayer,1,Nodes(2));     vvNode(2,2)=HydroParam%uxy(bbLayer,2,Nodes(2));     wwNode(2,2)=HydroParam%wfc(bbLayer,Nodes(2))
    !uuNode(2,3)=HydroParam%ubV(bbLayer,1,Nodes(3));     vvNode(2,3)=HydroParam%ubV(bbLayer,2,Nodes(3));     wwNode(2,3)=HydroParam%ubV(bbLayer,3,Nodes(3))
    !uuNode(3,1)=HydroParam%uNode(bbLayer+1,1,Nodes(1)); vvNode(3,1)=HydroParam%uNode(bbLayer+1,2,Nodes(1)); wwNode(3,1)=HydroParam%uNode(bbLayer+1,3,Nodes(1))
    !uuNode(3,2)=HydroParam%ug(Nodes(2),bbLayer+1);      vvNode(3,2)=HydroParam%vg(Nodes(2),bbLayer+1);      wwNode(3,2)=HydroParam%wg(Nodes(2),bbLayer+1)
    !uuNode(3,3)=HydroParam%uNode(bbLayer+1,1,Nodes(3)); vvNode(3,3)=HydroParam%uNode(bbLayer+1,2,Nodes(3)); wwNode(3,3)=HydroParam%uNode(bbLayer+1,3,Nodes(3))
    !! Vertical 2 = Cell Centered Section
    !uuNode(1,4)=HydroParam%ug(Nodes(4),bbLayer);      vvNode(1,4)=HydroParam%vg(Nodes(4),bbLayer);      	  wwNode(1,4)=HydroParam%wg(Nodes(4),bbLayer)
    !uuNode(1,5)=HydroParam%uxyL(bbLayer,1,Nodes(5));  vvNode(1,5)=HydroParam%uxyL(bbLayer,2,Nodes(5));  	  wwNode(1,5)=HydroParam%w(bbLayer,Nodes(5))
    !uuNode(1,6)=HydroParam%ug(Nodes(6),bbLayer);      vvNode(1,6)=HydroParam%vg(Nodes(6),bbLayer);      	  wwNode(1,6)=HydroParam%wg(Nodes(6),bbLayer)
    !uuNode(2,4)=HydroParam%uxy(bbLayer,1,Nodes(4));     vvNode(2,4)=HydroParam%uxy(bbLayer,2,Nodes(4));     wwNode(2,4)=HydroParam%wfc(bbLayer,Nodes(4))
    !uuNode(2,5)=HydroParam%ub(bbLayer,1,Nodes(5));     vvNode(2,5)=HydroParam%ub(bbLayer,2,Nodes(5));       wwNode(2,5)=HydroParam%ub(bbLayer,3,Nodes(5))
    !uuNode(2,6)=HydroParam%uxy(bbLayer,1,Nodes(6));     vvNode(2,6)=HydroParam%uxy(bbLayer,2,Nodes(6));     wwNode(2,6)=HydroParam%wfc(bbLayer,Nodes(6))
    !uuNode(3,4)=HydroParam%ug(Nodes(4),bbLayer+1);      vvNode(3,4)=HydroParam%vg(Nodes(4),bbLayer+1);      wwNode(3,4)=HydroParam%wg(Nodes(4),bbLayer+1)
    !uuNode(3,5)=HydroParam%uxyL(bbLayer+1,1,Nodes(5));  vvNode(3,5)=HydroParam%uxyL(bbLayer+1,2,Nodes(5));  wwNode(3,5)=HydroParam%w(bbLayer+1,Nodes(5))
    !uuNode(3,6)=HydroParam%ug(Nodes(6),bbLayer+1);      vvNode(3,6)=HydroParam%vg(Nodes(6),bbLayer+1);      wwNode(3,6)=HydroParam%wg(Nodes(6),bbLayer+1)
    !! Vertical 3 = East Edge
    !uuNode(1,7)=HydroParam%uNode(bbLayer,1,Nodes(7));   vvNode(1,7)=HydroParam%uNode(bbLayer,2,Nodes(7));   wwNode(1,7)=HydroParam%uNode(bbLayer,3,Nodes(7))
    !uuNode(1,8)=HydroParam%ug(Nodes(8),bbLayer);        vvNode(1,8)=HydroParam%vg(Nodes(8),bbLayer);        wwNode(1,8)=HydroParam%wg(Nodes(8),bbLayer)
    !uuNode(1,9)=HydroParam%uNode(bbLayer,1,Nodes(9));   vvNode(1,9)=HydroParam%uNode(bbLayer,2,Nodes(9));   wwNode(1,9)=HydroParam%uNode(bbLayer,3,Nodes(9))
    !uuNode(2,7)=HydroParam%ubV(bbLayer,1,Nodes(7));     vvNode(2,7)=HydroParam%ubV(bbLayer,2,Nodes(7));     wwNode(2,7)=HydroParam%ubV(bbLayer,3,Nodes(7))
    !uuNode(2,8)=HydroParam%uxy(bbLayer,1,Nodes(8));     vvNode(2,8)=HydroParam%uxy(bbLayer,2,Nodes(8));     wwNode(2,8)=HydroParam%wfc(bbLayer,Nodes(8))
    !uuNode(2,9)=HydroParam%ubV(bbLayer,1,Nodes(9));     vvNode(2,9)=HydroParam%ubV(bbLayer,2,Nodes(9));     wwNode(2,9)=HydroParam%ubV(bbLayer,3,Nodes(9))
    !uuNode(3,7)=HydroParam%uNode(bbLayer+1,1,Nodes(7)); vvNode(3,7)=HydroParam%uNode(bbLayer+1,2,Nodes(7)); wwNode(3,7)=HydroParam%uNode(bbLayer+1,3,Nodes(7))
    !uuNode(3,8)=HydroParam%ug(Nodes(8),bbLayer+1);      vvNode(3,8)=HydroParam%vg(Nodes(8),bbLayer+1);      wwNode(3,8)=HydroParam%wg(Nodes(8),bbLayer+1)
    !uuNode(3,9)=HydroParam%uNode(bbLayer+1,1,Nodes(9)); vvNode(3,9)=HydroParam%uNode(bbLayer+1,2,Nodes(9)); wwNode(3,9)=HydroParam%uNode(bbLayer+1,3,Nodes(9))
    !!    
    !!If (bbLayer==HydroParam%ElCapitalM(bbElem)) Then
    !!    !xxNode(Layer,Node Position), Layer = [k-1/2, k , k +1/2] == [1,2,3]
    !!    xxNode(1,1) = MeshParam%xNode(Nodes(1));      yyNode(1,1) = MeshParam%yNode(Nodes(1));      zzNode(1,1)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
    !!    xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
    !!    xxNode(1,3) = MeshParam%xNode(Nodes(3));      yyNode(1,3) = MeshParam%yNode(Nodes(3));      zzNode(1,3)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
    !!    xxNode(2,1) = MeshParam%xNode(Nodes(1));      yyNode(2,1) = MeshParam%yNode(Nodes(1));      zzNode(2,1)=(HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%peta(Nodes(1)))*0.5d0
    !!    xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=(HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%Z(bbLayer+1,Nodes(2)))*0.5d0
    !!    xxNode(2,3) = MeshParam%xNode(Nodes(3));      yyNode(2,3) = MeshParam%yNode(Nodes(3));      zzNode(2,3)=(HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%peta(Nodes(3)))*0.5d0
    !!    xxNode(3,1) = MeshParam%xNode(Nodes(1));      yyNode(3,1) = MeshParam%yNode(Nodes(1));      zzNode(3,1)=HydroParam%peta(Nodes(1))
    !!    xxNode(3,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(3,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(3,2)=HydroParam%Z(bbLayer+1,Nodes(2))
    !!    xxNode(3,3) = MeshParam%xNode(Nodes(3));      yyNode(3,3) = MeshParam%yNode(Nodes(3));      zzNode(3,3)=HydroParam%peta(Nodes(3))
    !!    
    !!    xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
    !!    xxNode(1,5) = MeshParam%xb(Nodes(5));         yyNode(1,5) = MeshParam%yb(Nodes(5));         zzNode(1,5)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
    !!    xxNode(1,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(1,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(1,6)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
    !!    xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));  zzNode(2,4)=(HydroParam%Ze(bbLayer,bbElem)  + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%Z(bbLayer+1,Nodes(4)))*0.5d0
    !!    xxNode(2,5) = MeshParam%xb(Nodes(5));         yyNode(2,5) = MeshParam%yb(Nodes(5));          zzNode(2,5)=(HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%eta(Nodes(5)))*0.5d0
    !!    xxNode(2,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(2,6) = MeshParam%EdgeBary(2,Nodes(6));  zzNode(2,6)=(HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%Z(bbLayer+1,Nodes(6)))*0.5d0
    !!    xxNode(3,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(3,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(3,4)=HydroParam%Z(bbLayer+1,Nodes(4))
    !!    xxNode(3,5) = MeshParam%xb(Nodes(5));         yyNode(3,5) = MeshParam%yb(Nodes(5));         zzNode(3,5)=HydroParam%eta(Nodes(5))
    !!    xxNode(3,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(3,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(3,6)=HydroParam%Z(bbLayer+1,Nodes(6))
    !!    
    !!    xxNode(1,7) = MeshParam%xNode(Nodes(7));      yyNode(1,7) = MeshParam%yNode(Nodes(7));      zzNode(1,7)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
    !!    xxNode(1,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(1,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(1,8)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
    !!    xxNode(1,9) = MeshParam%xNode(Nodes(9));      yyNode(1,9) = MeshParam%yNode(Nodes(9));      zzNode(1,9)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
    !!    xxNode(2,7) = MeshParam%xNode(Nodes(7));      yyNode(2,7) = MeshParam%yNode(Nodes(7));      zzNode(2,7)=(HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%peta(Nodes(7)))*0.5d0 
    !!    xxNode(2,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(2,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(2,8)=(HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%Z(bbLayer+1,Nodes(8)) )*0.5d0
    !!    xxNode(2,9) = MeshParam%xNode(Nodes(9));      yyNode(2,9) = MeshParam%yNode(Nodes(9));      zzNode(2,9)=(HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem)) + HydroParam%peta(Nodes(9)))*0.5d0
    !!    xxNode(3,7) = MeshParam%xNode(Nodes(7));      yyNode(3,7) = MeshParam%yNode(Nodes(7));      zzNode(3,7)=HydroParam%peta(Nodes(7))
    !!    xxNode(3,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(3,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(3,8)=HydroParam%Z(bbLayer+1,Nodes(8))
    !!    xxNode(3,9) = MeshParam%xNode(Nodes(9));      yyNode(3,9) = MeshParam%yNode(Nodes(9));      zzNode(3,9)=HydroParam%peta(Nodes(9))      
    !!Else
    !!    xxNode(1,1) = MeshParam%xNode(Nodes(1));      yyNode(1,1) = MeshParam%yNode(Nodes(1));      zzNode(1,1)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))
    !!    xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
    !!    xxNode(1,3) = MeshParam%xNode(Nodes(3));      yyNode(1,3) = MeshParam%yNode(Nodes(3));      zzNode(1,3)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
    !!    xxNode(2,1) = MeshParam%xNode(Nodes(1));      yyNode(2,1) = MeshParam%yNode(Nodes(1));      zzNode(2,1)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2     
    !!    xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2   
    !!    xxNode(2,3) = MeshParam%xNode(Nodes(3));      yyNode(2,3) = MeshParam%yNode(Nodes(3));      zzNode(2,3)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2   
    !!    xxNode(3,1) = MeshParam%xNode(Nodes(1));      yyNode(3,1) = MeshParam%yNode(Nodes(1));      zzNode(3,1)=HydroParam%Ze(bbLayer+1,bbElem)
    !!    xxNode(3,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(3,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(3,2)=HydroParam%Ze(bbLayer+1,bbElem)    
    !!    xxNode(3,3) = MeshParam%xNode(Nodes(3));      yyNode(3,3) = MeshParam%yNode(Nodes(3));      zzNode(3,3)=HydroParam%Ze(bbLayer+1,bbElem)     
    !!    
    !!    xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
    !!    xxNode(1,5) = MeshParam%xb(Nodes(5));         yyNode(1,5) = MeshParam%yb(Nodes(5));         zzNode(1,5)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
    !!    xxNode(1,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(1,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(1,6)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
    !!    xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));  zzNode(2,4)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2    
    !!    xxNode(2,5) = MeshParam%xb(Nodes(5));         yyNode(2,5) = MeshParam%yb(Nodes(5));          zzNode(2,5)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2  
    !!    xxNode(2,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(2,6) = MeshParam%EdgeBary(2,Nodes(6));  zzNode(2,6)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2   
    !!    xxNode(3,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(3,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(3,4)=HydroParam%Ze(bbLayer+1,bbElem)   
    !!    xxNode(3,5) = MeshParam%xb(Nodes(5));         yyNode(3,5) = MeshParam%yb(Nodes(5));         zzNode(3,5)=HydroParam%Ze(bbLayer+1,bbElem)    
    !!    xxNode(3,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(3,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(3,6)=HydroParam%Ze(bbLayer+1,bbElem)    
    !!    
    !!    xxNode(1,7) = MeshParam%xNode(Nodes(7));      yyNode(1,7) = MeshParam%yNode(Nodes(7));      zzNode(1,7)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
    !!    xxNode(1,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(1,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(1,8)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
    !!    xxNode(1,9) = MeshParam%xNode(Nodes(9));      yyNode(1,9) = MeshParam%yNode(Nodes(9));      zzNode(1,9)=HydroParam%Ze(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))     
    !!    xxNode(2,7) = MeshParam%xNode(Nodes(7));      yyNode(2,7) = MeshParam%yNode(Nodes(7));      zzNode(2,7)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2     
    !!    xxNode(2,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(2,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(2,8)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2   
    !!    xxNode(2,9) = MeshParam%xNode(Nodes(9));      yyNode(2,9) = MeshParam%yNode(Nodes(9));      zzNode(2,9)=HydroParam%Zb(bbLayer,bbElem) + sum(HydroParam%DZsi(:,bbElem))/2   
    !!    xxNode(3,7) = MeshParam%xNode(Nodes(7));      yyNode(3,7) = MeshParam%yNode(Nodes(7));      zzNode(3,7)=HydroParam%Ze(bbLayer+1,bbElem)  
    !!    xxNode(3,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(3,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(3,8)=HydroParam%Ze(bbLayer+1,bbElem)    
    !!    xxNode(3,9) = MeshParam%xNode(Nodes(9));      yyNode(3,9) = MeshParam%yNode(Nodes(9));      zzNode(3,9)=HydroParam%Ze(bbLayer+1,bbElem)  
    !!EndIf
    !        
    !If (bbLayer==HydroParam%ElCapitalM(bbElem)) Then
    !    !xxNode(Layer,Node Position), Layer = [k-1/2, k , k +1/2] == [1,2,3]
    !    xxNode(1,1) = MeshParam%xNode(Nodes(1));      yyNode(1,1) = MeshParam%yNode(Nodes(1));      zzNode(1,1)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(1,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(1,2)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,3) = MeshParam%xNode(Nodes(3));      yyNode(1,3) = MeshParam%yNode(Nodes(3));      zzNode(1,3)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(2,1) = MeshParam%xNode(Nodes(1));      yyNode(2,1) = MeshParam%yNode(Nodes(1));      zzNode(2,1)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%peta(Nodes(1)))*0.5d0
    !    xxNode(2,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(2,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(2,2)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%Z(bbLayer+1,Nodes(2)))*0.5d0
    !    xxNode(2,3) = MeshParam%xNode(Nodes(3));      yyNode(2,3) = MeshParam%yNode(Nodes(3));      zzNode(2,3)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%peta(Nodes(3)))*0.5d0
    !    xxNode(3,1) = MeshParam%xNode(Nodes(1));      yyNode(3,1) = MeshParam%yNode(Nodes(1));      zzNode(3,1)=HydroParam%peta(Nodes(1))
    !    xxNode(3,2) = MeshParam%EdgeBary(1,Nodes(2)); yyNode(3,2) = MeshParam%EdgeBary(2,Nodes(2)); zzNode(3,2)=HydroParam%Z(bbLayer+1,Nodes(2))
    !    xxNode(3,3) = MeshParam%xNode(Nodes(3));      yyNode(3,3) = MeshParam%yNode(Nodes(3));      zzNode(3,3)=HydroParam%peta(Nodes(3))
    !    
    !    xxNode(1,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(1,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(1,4)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,5) = MeshParam%xb(Nodes(5));         yyNode(1,5) = MeshParam%yb(Nodes(5));         zzNode(1,5)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(1,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(1,6)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(2,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(2,4) = MeshParam%EdgeBary(2,Nodes(4));  zzNode(2,4)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%Z(bbLayer+1,Nodes(4)))*0.5d0
    !    xxNode(2,5) = MeshParam%xb(Nodes(5));         yyNode(2,5) = MeshParam%yb(Nodes(5));          zzNode(2,5)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%eta(Nodes(5)))*0.5d0
    !    xxNode(2,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(2,6) = MeshParam%EdgeBary(2,Nodes(6));  zzNode(2,6)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%Z(bbLayer+1,Nodes(6)))*0.5d0
    !    xxNode(3,4) = MeshParam%EdgeBary(1,Nodes(4)); yyNode(3,4) = MeshParam%EdgeBary(2,Nodes(4)); zzNode(3,4)=HydroParam%Z(bbLayer+1,Nodes(4))
    !    xxNode(3,5) = MeshParam%xb(Nodes(5));         yyNode(3,5) = MeshParam%yb(Nodes(5));         zzNode(3,5)=HydroParam%eta(Nodes(5))
    !    xxNode(3,6) = MeshParam%EdgeBary(1,Nodes(6)); yyNode(3,6) = MeshParam%EdgeBary(2,Nodes(6)); zzNode(3,6)=HydroParam%Z(bbLayer+1,Nodes(6))
    !    
    !    xxNode(1,7) = MeshParam%xNode(Nodes(7));      yyNode(1,7) = MeshParam%yNode(Nodes(7));      zzNode(1,7)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(1,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(1,8)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(1,9) = MeshParam%xNode(Nodes(9));      yyNode(1,9) = MeshParam%yNode(Nodes(9));      zzNode(1,9)=HydroParam%Ze(bbLayer,bbElem)
    !    xxNode(2,7) = MeshParam%xNode(Nodes(7));      yyNode(2,7) = MeshParam%yNode(Nodes(7));      zzNode(2,7)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%peta(Nodes(7)))*0.5d0 
    !    xxNode(2,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(2,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(2,8)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%Z(bbLayer+1,Nodes(8)) )*0.5d0
    !    xxNode(2,9) = MeshParam%xNode(Nodes(9));      yyNode(2,9) = MeshParam%yNode(Nodes(9));      zzNode(2,9)=(HydroParam%Ze(bbLayer,bbElem) + HydroParam%peta(Nodes(9)))*0.5d0
    !    xxNode(3,7) = MeshParam%xNode(Nodes(7));      yyNode(3,7) = MeshParam%yNode(Nodes(7));      zzNode(3,7)=HydroParam%peta(Nodes(7))
    !    xxNode(3,8) = MeshParam%EdgeBary(1,Nodes(8)); yyNode(3,8) = MeshParam%EdgeBary(2,Nodes(8)); zzNode(3,8)=HydroParam%Z(bbLayer+1,Nodes(8))
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
    !return
    !End Subroutine iQuadraticNodes  