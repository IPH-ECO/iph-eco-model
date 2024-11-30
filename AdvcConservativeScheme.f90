!Module AdvcConservativeScheme
!    private
!    public:: ConservativeFuFv
!    !public:: ComputeFW
!    public:: signa
!    !-----------------------------------------------------------------
!    Contains   
!    
!        Subroutine ConservativeFuFv(HydroParam,MeshParam,dt)
!    
!        Use MeshVars !, Only: Left,Right,nElem,nEdge,EdgeBary,NormalVector,TangentVector,Edge,EdgeLength,Quadri,xNode,yNode,Area,xb,yb,nNode,Kmax,nEgdesatNode,EgdesatNode
!        Use Hydrodynamic !, Only: Smallm,CapitalM,Z,u,utang,Fu,Fv,H,Pcri,w,Nfut,uxy,uxyback,ElSmallm,ElCapitalM,lat,OMEGA,Pi,Wu,HorViscosity,Fub,FuxyNode
!    
!        Implicit None
!        Integer :: iElem,iEdge,iEdgein,lEdge,iLayer,l,r,j,Sig,Face,n1,n2,n3,n4,n5,n6,n7,n8,n9,iNode,iCount,fac,kin,nu,nu1,nl,rElem,lElem,uElem,psi_flag,nr
!        Real :: NearZero = 1e-10
!        Real :: fi_small = 0.0d0
!        Real :: epsGrad = 100000 !CAYO
!        Real :: dt, weit
!        Real :: uuNode(9), vvNode(9), wwNode(9), xxNode(9), yyNode(9)
!        Real :: Nodes(9), dzNode(9)
!        Real :: rj(2,4), psi(2,4), ru
!        Real :: Courant
!        type(MeshGridParam) :: MeshParam
!        type(HydrodynamicParam) :: HydroParam
!
!
!        Do iElem = 1,MeshParam%nElem
!            
!            Do iEdge = 1,4
!                Face = MeshParam%Edge(iEdge,iElem)
!                l = MeshParam%Left(Face)
!                r = MeshParam%Right(Face)
!            
!                ! In the case that Cell is dry, the backtrack tang velocities is null.
!                ! Note that in subsurface coupled case %H(iEdge) can take the freatic water level in Cell.
!                ! DZsj represents the Edge freatic component thickness, thus H - DZsj is equivalent only surface water level
!                ! In Cell where no subsufarce layer exist or in only surface flow simulation DZsj is set equal 0.0d0 (ReadHydroIniCond Module)
!                If (HydroParam%H(iEdge) + HydroParam%sj(iEdge)-HydroParam%hj(iEdge)<=HydroParam%Pcri+NearZero) Then !CAYO
!                    Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
!                        HydroParam%uxyback(iLayer,1:2,iEdge) = (/ 0., 0. /)
!                    EndDo
!                    Cycle
!                EndIf
!            
!                If (iEdge==1) Then
!                    iEdgein = 3
!                    n7 = 2
!                    rElem = MeshParam%Right(MeshParam%Edge(4,l))
!                    lElem = MeshParam%Right(MeshParam%Edge(2,l))
!                    
!                    ! In the case that Cell is Wet:
!                    Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face)
!
!                        psi(:,:) = 0.0d0
!                        rj(:,:) = 0.0d0
!                        NodeValues(iEdgein, iEdge, Face, rElem,lElem,n7,r,l,uuNode(:),vvNode(:),xxNode(:), yyNode(:), dzNode(:))
!    
!                        Courant = 0.5*(HydroParam%Theta*HydroParam%u(iLayer,iEdge)*+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,iEdge))*dt/MeshParam%CirDistance(iEdge)
!
!                        !x - direction:
!                        If ((uuNode(3) - uuNode(7))/MeshParam%CirDistance(iEdge) > epsGrad > 0) Then
!                              
!                            If (uuNode(3) + uuNode(7) >= 0) Then  
!                                HydroParam%uArrow(iLayer,1,iEdge) = 0.5*(uuNode(7) + uuNode(3))
!                                
!                                rj(1,1) = ru(uuNode(3),uuNode(7),uuNode(6),xxNode(3),xxNode(7),xxNode(6))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(8),uuNode(3),uuNode(7),xxNode(8),xxNode(3),xxNode(7))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(3),vvNode(7),vvNode(6),xxNode(3),xxNode(7),xxNode(6))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(8),vvNode(3),vvNode(7),xxNode(8),xxNode(3),xxNode(7))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))  
!                                                        
!                            Else                                                
!                                HydroParam%uArrow(iLayer,1,iEdge) = 0.5*(uuNode(8) + uuNode(3))
!                                                              
!                                rj(1,1) = ru(uuNode(8),uuNode(3),uuNode(7),xxNode(8),xxNode(3),xxNode(7))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(9),uuNode(8),uuNode(3),xxNode(9),xxNode(8),xxNode(3))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(8),vvNode(3),vvNode(7),xxNode(8),xxNode(3),xxNode(7))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(9),vvNode(8),vvNode(3),xxNode(9),xxNode(8),xxNode(3))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))                                 
!                            EndIf                                                       
!                            ! y-direction:
!                            If (vvNode(3) + vvNode(2)>= 0) Then
!                                HydroParam%uArrow(iLayer,2,iEdge) =  0.5*(vvNode(3) + vvNode(2))
!                    
!                                rj(1,3) = ru(uuNode(3),uuNode(2),uuNode(1),yyNode(3),yyNode(2),yyNode(1))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(4),uuNode(3),uuNode(2),yyNode(4),yyNode(3),yyNode(2))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                     
!                                 
!                                rj(2,3) = ru(vvNode(3),vvNode(2),vvNode(1),yyNode(3),yyNode(2),yyNode(1))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!                
!                                rj(2,4) = ru(vvNode(4),vvNode(3),vvNode(2),yyNode(4),yyNode(3),yyNode(2))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))  
!                
!                            Else
!                                HydroParam%uArrow(iLayer,2,iEdge) =  0.5*(vvNode(4) + vvNode(3))
!                
!                                rj(1,3) = ru(uuNode(4),uuNode(3),uuNode(2),yyNode(4),yyNode(3),yyNode(2))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(5),uuNode(4),uuNode(3),yyNode(5),yyNode(4),yyNode(3))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                    
!                    
!                                rj(2,3) = ru(vvNode(4),vvNode(3),vvNode(2),yyNode(4),yyNode(3),yyNode(2))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!
!                                rj(2,4) = ru(vvNode(5),vvNode(4),vvNode(3),yyNode(5),yyNode(4),yyNode(3))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))                                                   
!                            EndIf                            
!                        Else
!                            ! x-direction:
!                            If (uuNode(7)*dzNode(7) + uuNode(3)*dzNode(3) >= 0 ) Then
!                                HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(7)*dzNode(7) + uuNode(3)*dzNode(3))/(dzNode(7) + dzNode(3))                                
!
!                                rj(1,1) = ru(uuNode(3),uuNode(7),uuNode(6),xxNode(3),xxNode(7),xxNode(6))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(8),uuNode(3),uuNode(7),xxNode(8),xxNode(3),xxNode(7))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(3),vvNode(7),vvNode(6),xxNode(3),xxNode(7),xxNode(6))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(8),vvNode(3),vvNode(7),xxNode(8),xxNode(3),xxNode(7))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))                                                     
!                            Else
!                                HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(8)*dzNode(8) + uuNode(3)*dzNode(3))/(dzNode(8) + dzNode(3))                                
!                                                         
!                                rj(1,1) = ru(uuNode(8),uuNode(3),uuNode(7),xxNode(8),xxNode(3),xxNode(7))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(9),uuNode(8),uuNode(3),xxNode(9),xxNode(8),xxNode(3))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(8),vvNode(3),vvNode(7),xxNode(8),xxNode(3),xxNode(7))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(9),vvNode(8),vvNode(3),xxNode(9),xxNode(8),xxNode(3))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))  
!                            EndIf
!                                             
!                        ! y-direction:                                                                             
!                            If (vvNode(3)*dzNode(3) + vvNode(2)*dzNode(2)>= 0) Then
!                                HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(3)*dzNode(3) + vvNode(2)*dzNode(2))/(dzNode(3) + dzNode(2))
!                    
!                                rj(1,3) = ru(uuNode(3),uuNode(2),uuNode(1),yyNode(3),yyNode(2),yyNode(1))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(4),uuNode(3),uuNode(2),yyNode(4),yyNode(3),yyNode(2))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                     
!                                 
!                                rj(2,3) = ru(vvNode(3),vvNode(2),vvNode(1),yyNode(3),yyNode(2),yyNode(1))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!                
!                                rj(2,4) = ru(vvNode(4),vvNode(3),vvNode(2),yyNode(4),yyNode(3),yyNode(2))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))            
!                            Else
!                                HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(4)*dzNode(4) + vvNode(3)*dzNode(3))/(dzNode(4) + dzNode(3))
!                
!                                rj(1,3) = ru(uuNode(4),uuNode(3),uuNode(2),yyNode(4),yyNode(3),yyNode(2))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(5),uuNode(4),uuNode(3),yyNode(5),yyNode(4),yyNode(3))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                    
!                    
!                                rj(2,3) = ru(vvNode(4),vvNode(3),vvNode(2),yyNode(4),yyNode(3),yyNode(2))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!
!                                rj(2,4) = ru(vvNode(5),vvNode(4),vvNode(3),yyNode(5),yyNode(4),yyNode(3))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))                                                   
!                            EndIf                                                                                                            
!                        EndIf                 
!                    EndDo
!                    
!                ElseIf (iEdge==2) Then
!                    iEdgein = 4
!                    n9 = 3
!                    rElem = MeshParam%Right(MeshParam%Edge(1,l)) 
!                    lElem = MeshParam%Right(MeshParam%Edge(3,l))   
!                    
!                    Courant = 0.5*(HydroParam%Theta*HydroParam%u(iLayer,iEdge)*+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,iEdge))*dt/MeshParam%CirDistance(iEdge)
!             
!                    ! In the case that Cell is Wet:
!                    Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face)
!
!                        psi(:,:) = 0.0d0
!                        rj(:,:) = 0.0d0
!                        NodeValues(iEdgein,iEdge,Face,rElem,lElem,n7,r,l,uuNode(:),vvNode(:),xxNode(:), yyNode(:), dzNode(:))
!
!                        !x - direction:
!                        If ((uuNode(3) - uuNode(4))/MeshParam%CirDistance(iEdge) > epsGrad > 0) Then
!                              
!                            If (uuNode(3) + uuNode(4)>= 0) Then  
!                                HydroParam%uArrow(iLayer,1,iEdge) = 0.5*(uuNode(3) + uuNode(4))
!                                
!                                rj(1,1) = ru(uuNode(3),uuNode(4),uuNode(5),xxNode(3),xxNode(4),xxNode(5))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(2),uuNode(3),uuNode(4),xxNode(2),xxNode(3),xxNode(4))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(3),vvNode(4),vvNode(5),xxNode(3),xxNode(4),xxNode(5))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(2),vvNode(3),vvNode(4),xxNode(2),xxNode(3),xxNode(4))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))  
!                                                        
!                            Else                                                
!                                HydroParam%uArrow(iLayer,1,iEdge) = 0.5*(uuNode(2) + uuNode(3))
!                                                              
!                                rj(1,1) = ru(uuNode(2),uuNode(3),uuNode(4),xxNode(2),xxNode(3),xxNode(4))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(1),uuNode(2),uuNode(3),xxNode(1),xxNode(2),xxNode(3))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(2),vvNode(3),vvNode(4),xxNode(2),xxNode(3),xxNode(4))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(1),vvNode(2),vvNode(3),xxNode(1),xxNode(2),xxNode(3))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))                                 
!                            EndIf                                                       
!                            ! y-direction:
!                            If (vvNode(3) + vvNode(7) >= 0) Then
!                                HydroParam%uArrow(iLayer,2,iEdge) =  0.5*(vvNode(3) + vvNode(7))
!                    
!                                rj(1,3) = ru(uuNode(3),uuNode(7),uuNode(6),yyNode(3),yyNode(7),yyNode(6))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(8),uuNode(3),uuNode(7),yyNode(8),yyNode(3),yyNode(7))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                     
!                                 
!                                rj(2,3) = ru(vvNode(3),vvNode(7),vvNode(6),yyNode(3),yyNode(7),yyNode(6))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!                
!                                rj(2,4) = ru(vvNode(8),vvNode(3),vvNode(7),yyNode(8),yyNode(3),yyNode(7))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))  
!                
!                            Else
!                                HydroParam%uArrow(iLayer,2,iEdge) =  0.5*(vvNode(8) + vvNode(3))
!                
!                                rj(1,3) = ru(uuNode(8),uuNode(3),uuNode(7),yyNode(8),yyNode(3),yyNode(7))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(9),uuNode(8),uuNode(3),yyNode(9),yyNode(8),yyNode(3))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                    
!                    
!                                rj(2,3) = ru(vvNode(8),vvNode(3),vvNode(7),yyNode(8),yyNode(3),yyNode(7))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!
!                                rj(2,4) = ru(vvNode(9),vvNode(8),vvNode(3),yyNode(9),yyNode(8),yyNode(3))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))                                                   
!                            EndIf                            
!                        Else
!                            ! x-direction:
!                            If (uuNode(4)*dzNode(4) + uuNode(3)*dzNode(3) >= 0 ) Then
!                                HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(4)*dzNode(4) + uuNode(3)*dzNode(3))/(dzNode(4) + dzNode(3))                                
!
!                                rj(1,1) = ru(uuNode(3),uuNode(4),uuNode(5),xxNode(3),xxNode(4),xxNode(5))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(2),uuNode(3),uuNode(4),xxNode(2),xxNode(3),xxNode(4))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(3),vvNode(4),vvNode(5),xxNode(3),xxNode(4),xxNode(5))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(2),vvNode(3),vvNode(4),xxNode(2),xxNode(3),xxNode(4))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))                                                     
!                            Else
!                                HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(2)*dzNode(2) + uuNode(3)*dzNode(3))/(dzNode(2) + dzNode(3))                                
!                                                         
!                                rj(1,1) = ru(uuNode(2),uuNode(3),uuNode(4),xxNode(2),xxNode(3),xxNode(4))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(1),uuNode(2),uuNode(3),xxNode(1),xxNode(2),xxNode(3))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(2),vvNode(3),vvNode(4),xxNode(2),xxNode(3),xxNode(4))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(1),vvNode(2),vvNode(3),xxNode(1),xxNode(2),xxNode(3))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))  
!                            EndIf
!                                             
!                        ! y-direction:                                                                             
!                            If (vvNode(3)*dzNode(3) + vvNode(7)*dzNode(7)>= 0) Then
!                                HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(3)*dzNode(3) + vvNode(7)*dzNode(7))/(dzNode(3) + dzNode(7))
!                    
!                                rj(1,3) = ru(uuNode(3),uuNode(7),uuNode(6),yyNode(3),yyNode(7),yyNode(6))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(8),uuNode(3),uuNode(7),yyNode(8),yyNode(3),yyNode(7))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                     
!                                 
!                                rj(2,3) = ru(vvNode(3),vvNode(7),vvNode(6),yyNode(3),yyNode(7),yyNode(6))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!                
!                                rj(2,4) = ru(vvNode(8),vvNode(3),vvNode(7),yyNode(8),yyNode(3),yyNode(7))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))            
!                            Else
!                                HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(8)*dzNode(8) + vvNode(3)*dzNode(3))/(dzNode(8) + dzNode(3))
!                
!                                rj(1,3) = ru(uuNode(8),uuNode(3),uuNode(7),yyNode(8),yyNode(3),yyNode(7))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(9),uuNode(8),uuNode(3),yyNode(9),yyNode(8),yyNode(3))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                    
!                    
!                                rj(2,3) = ru(vvNode(8),vvNode(3),vvNode(7),yyNode(8),yyNode(3),yyNode(7))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!
!                                rj(2,4) = ru(vvNode(9),vvNode(8),vvNode(3),yyNode(9),yyNode(8),yyNode(3))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))                                                   
!                            EndIf                                                                                                            
!                        EndIf
!                    EndDo
!                                                       
!                ElseIf (iEdge==3) Then
!                    iEdgein = 1
!                    n9 = 4
!                    rElem = MeshParam%Right(MeshParam%Edge(2,l)) 
!                    lElem = MeshParam%Right(MeshParam%Edge(4,l))
!                    
!                    Courant = 0.5*(HydroParam%Theta*HydroParam%u(iLayer,iEdge)*+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,iEdge))*dt/MeshParam%CirDistance(iEdge)
!                    
!                   ! In the case that Cell is Wet:
!                    Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face)
!
!                        psi(:,:) = 0.0d0
!                        rj(:,:) = 0.0d0
!                        NodeValues(iEdgein,iEdge,Face,rElem,lElem,n7,r,l,uuNode(:),vvNode(:),xxNode(:), yyNode(:), dzNode(:))                    
!  
!                        !x - direction:
!                        If ((uuNode(3) - uuNode(8))/MeshParam%CirDistance(iEdge) > epsGrad > 0) Then
!                              
!                            If (uuNode(3) + uuNode(8)>= 0) Then  
!                                HydroParam%uArrow(iLayer,1,iEdge) = 0.5*(uuNode(3) + uuNode(8))
!                                
!                                rj(1,1) = ru(uuNode(3),uuNode(8),uuNode(9),xxNode(3),xxNode(8),xxNode(9))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(7),uuNode(3),uuNode(8),xxNode(7),xxNode(3),xxNode(8))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(3),vvNode(8),vvNode(9),xxNode(3),xxNode(8),xxNode(9))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(7),vvNode(3),vvNode(8),xxNode(7),xxNode(3),xxNode(8))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))  
!                                                        
!                            Else                                                
!                                HydroParam%uArrow(iLayer,1,iEdge) = 0.5*(uuNode(7) + uuNode(3))
!                                                              
!                                rj(1,1) = ru(uuNode(7),uuNode(3),uuNode(8),xxNode(7),xxNode(3),xxNode(8))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(6),uuNode(7),uuNode(3),xxNode(6),xxNode(7),xxNode(3))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(7),vvNode(3),vvNode(8),xxNode(7),xxNode(3),xxNode(8))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(6),vvNode(7),vvNode(3),xxNode(6),xxNode(7),xxNode(3))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))                                 
!                            EndIf                                                       
!                            ! y-direction:
!                            If (vvNode(3) + vvNode(4) >= 0) Then
!                                HydroParam%uArrow(iLayer,2,iEdge) =  0.5*(vvNode(3) + vvNode(4))
!                    
!                                rj(1,3) = ru(uuNode(3),uuNode(4),uuNode(5),yyNode(3),yyNode(4),yyNode(5))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(2),uuNode(3),uuNode(4),yyNode(2),yyNode(3),yyNode(4))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                     
!                                 
!                                rj(2,3) = ru(vvNode(3),vvNode(4),vvNode(5),yyNode(3),yyNode(4),yyNode(5))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!                
!                                rj(2,4) = ru(vvNode(2),vvNode(3),vvNode(4),yyNode(2),yyNode(3),yyNode(4))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))  
!                
!                            Else
!                                HydroParam%uArrow(iLayer,2,iEdge) =  0.5*(vvNode(2) + vvNode(3))
!                
!                                rj(1,3) = ru(uuNode(2),uuNode(3),uuNode(4),yyNode(2),yyNode(3),yyNode(4))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(1),uuNode(2),uuNode(3),yyNode(1),yyNode(2),yyNode(3))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                    
!                    
!                                rj(2,3) = ru(vvNode(2),vvNode(3),vvNode(4),yyNode(2),yyNode(3),yyNode(4))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!
!                                rj(2,4) = ru(vvNode(1),vvNode(2),vvNode(3),yyNode(1),yyNode(2),yyNode(3))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))                                                   
!                            EndIf                            
!                        Else
!                            ! x-direction:
!                            If (uuNode(3)*dzNode(3) + uuNode(8)*dzNode(8) >= 0 ) Then
!                                HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(3)*dzNode(3) + uuNode(8)*dzNode(8))/(dzNode(3) + dzNode(8))                                
!
!                                rj(1,1) = ru(uuNode(3),uuNode(8),uuNode(9),xxNode(3),xxNode(8),xxNode(9))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(7),uuNode(3),uuNode(8),xxNode(7),xxNode(3),xxNode(8))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(3),vvNode(8),vvNode(9),xxNode(3),xxNode(8),xxNode(9))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(7),vvNode(3),vvNode(8),xxNode(7),xxNode(3),xxNode(8))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))                                                      
!                            Else
!                                HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(7)*dzNode(7) + uuNode(3)*dzNode(3))/(dzNode(7) + dzNode(3))                                
!                                                                                                                 
!                                rj(1,1) = ru(uuNode(7),uuNode(3),uuNode(8),xxNode(7),xxNode(3),xxNode(8))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(6),uuNode(7),uuNode(3),xxNode(6),xxNode(7),xxNode(3))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(7),vvNode(3),vvNode(8),xxNode(7),xxNode(3),xxNode(8))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(6),vvNode(7),vvNode(3),xxNode(6),xxNode(7),xxNode(3))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))   
!                            EndIf
!                                             
!                        ! y-direction:                                                                             
!                            If (vvNode(3)*dzNode(3) + vvNode(4)*dzNode(4)>= 0) Then
!                                HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(3)*dzNode(3) + vvNode(4)*dzNode(4))/(dzNode(3) + dzNode(4))
!                    
!                                rj(1,3) = ru(uuNode(3),uuNode(4),uuNode(5),yyNode(3),yyNode(4),yyNode(5))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(2),uuNode(3),uuNode(4),yyNode(2),yyNode(3),yyNode(4))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                     
!                                 
!                                rj(2,3) = ru(vvNode(3),vvNode(4),vvNode(5),yyNode(3),yyNode(4),yyNode(5))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!                
!                                rj(2,4) = ru(vvNode(2),vvNode(3),vvNode(4),yyNode(2),yyNode(3),yyNode(4))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))         
!                            Else
!                                HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(2)*dzNode(2) + vvNode(3)*dzNode(3))/(dzNode(2) + dzNode(3))
!                                
!                                rj(1,3) = ru(uuNode(2),uuNode(3),uuNode(4),yyNode(2),yyNode(3),yyNode(4))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(1),uuNode(2),uuNode(3),yyNode(1),yyNode(2),yyNode(3))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                    
!                    
!                                rj(2,3) = ru(vvNode(2),vvNode(3),vvNode(4),yyNode(2),yyNode(3),yyNode(4))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!
!                                rj(2,4) = ru(vvNode(1),vvNode(2),vvNode(3),yyNode(1),yyNode(2),yyNode(3))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))                                                      
!                            EndIf                                                                                                            
!                        EndIf
!                    EndDo
!                                                   
!                ElseIf (iEdge==4) Then
!                    iEdgein = 2
!                    n9 = 1
!                    rElem = MeshParam%Right(MeshParam%Edge(3,l)) 
!                    lElem = MeshParam%Right(MeshParam%Edge(1,l))
!                    
!                    Courant = 0.5*(HydroParam%Theta*HydroParam%u(iLayer,iEdge)*+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,iEdge))*dt/MeshParam%CirDistance(iEdge)
!                    
!                   ! In the case that Cell is Wet:
!                    Do iLayer = HydroParam%Smallm(Face), HydroParam%CapitalM(Face)
!
!                        psi(:,:) = 0.0d0
!                        rj(:,:) = 0.0d0
!                        NodeValues(iEdgein,iEdge,Face,rElem,lElem,n7,r,l,uuNode(:),vvNode(:),xxNode(:), yyNode(:), dzNode(:))
!                        
!                        !x - direction:
!                        If ((uuNode(3) - uuNode(2))/MeshParam%CirDistance(iEdge) > epsGrad > 0) Then
!                              
!                            If (uuNode(3) + uuNode(2)>= 0) Then  
!                                HydroParam%uArrow(iLayer,1,iEdge) = 0.5*(uuNode(3) + uuNode(2))
!                                
!                                rj(1,1) = ru(uuNode(3),uuNode(2),uuNode(1),xxNode(3),xxNode(2),xxNode(1))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(4),uuNode(3),uuNode(2),xxNode(4),xxNode(3),xxNode(2))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(3),vvNode(2),vvNode(1),xxNode(3),xxNode(2),xxNode(1))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(4),vvNode(3),vvNode(2),xxNode(4),xxNode(3),xxNode(2))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))  
!                                                                                                          
!                            Else                                                
!                                HydroParam%uArrow(iLayer,1,iEdge) = 0.5*(uuNode(4) + uuNode(3))
!                                                              
!                                rj(1,1) = ru(uuNode(4),uuNode(3),uuNode(2),xxNode(4),xxNode(3),xxNode(2))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(5),uuNode(4),uuNode(3),xxNode(5),xxNode(4),xxNode(3))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(4),vvNode(3),vvNode(2),xxNode(4),xxNode(3),xxNode(2))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(5),vvNode(4),vvNode(3),xxNode(5),xxNode(4),xxNode(3))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))                                 
!                            EndIf                                                       
!                            ! y-direction:
!                            If (vvNode(3) + vvNode(8) >= 0) Then
!                                HydroParam%uArrow(iLayer,2,iEdge) =  0.5*(vvNode(3) + vvNode(8))
!                    
!                                rj(1,3) = ru(uuNode(3),uuNode(8),uuNode(9),yyNode(3),yyNode(8),yyNode(9))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(7),uuNode(3),uuNode(8),yyNode(7),yyNode(3),yyNode(8))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                     
!                                 
!                                rj(2,3) = ru(vvNode(3),vvNode(8),vvNode(9),yyNode(3),yyNode(8),yyNode(9))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!                
!                                rj(2,4) = ru(vvNode(7),vvNode(3),vvNode(8),yyNode(7),yyNode(3),yyNode(8))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))  
!                
!                            Else
!                                HydroParam%uArrow(iLayer,2,iEdge) =  0.5*(vvNode(7) + vvNode(3))
!                
!                                rj(1,3) = ru(uuNode(7),uuNode(3),uuNode(8),yyNode(7),yyNode(3),yyNode(8))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(6),uuNode(7),uuNode(3),yyNode(6),yyNode(7),yyNode(3))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                    
!                    
!                                rj(2,3) = ru(vvNode(7),vvNode(3),vvNode(8),yyNode(7),yyNode(3),yyNode(8))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!
!                                rj(2,4) = ru(vvNode(6),vvNode(7),vvNode(3),yyNode(6),yyNode(7),yyNode(3))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))                                                   
!                            EndIf                            
!                        Else
!                            ! x-direction:
!                            If (uuNode(3)*dzNode(3) + uuNode(2)*dzNode(2) >= 0 ) Then
!                                HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(3)*dzNode(3) + uuNode(2)*dzNode(2))/(dzNode(3) + dzNode(2))                                
!
!                                rj(1,1) = ru(uuNode(3),uuNode(2),uuNode(1),xxNode(3),xxNode(2),xxNode(1))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(4),uuNode(3),uuNode(2),xxNode(4),xxNode(3),xxNode(2))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(3),vvNode(2),vvNode(1),xxNode(3),xxNode(2),xxNode(1))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(4),vvNode(3),vvNode(2),xxNode(4),xxNode(3),xxNode(2))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))                                                      
!                            Else
!                                HydroParam%uArrow(iLayer,1,iEdge) = (uuNode(4)*dzNode(4) + uuNode(3)*dzNode(3))/(dzNode(4) + dzNode(3))                                
!                                                                                                                 
!                                rj(1,1) = ru(uuNode(4),uuNode(3),uuNode(2),xxNode(4),xxNode(3),xxNode(2))
!                                Call Psi_value(psi_flag,rj(1,1),Courant,fi_small,psi(1,1))
!                
!                                rj(1,2) = ru(uuNode(5),uuNode(4),uuNode(3),xxNode(5),xxNode(4),xxNode(3))
!                                Call Psi_value(psi_flag,rj(1,2),Courant,fi_small,psi(1,2))
!                    
!                                rj(2,1) = ru(vvNode(4),vvNode(3),vvNode(2),xxNode(4),xxNode(3),xxNode(2))
!                                Call Psi_value(psi_flag,rj(2,1),Courant,fi_small,psi(2,1))
!                
!                                rj(2,2) = ru(vvNode(5),vvNode(4),vvNode(3),xxNode(5),xxNode(4),xxNode(3))
!                                Call Psi_value(psi_flag,rj(2,2),Courant,fi_small,psi(2,2))   
!                            EndIf
!                                             
!                        ! y-direction:                                                                             
!                            If (vvNode(3)*dzNode(3) + vvNode(8)*dzNode(8)>= 0) Then
!                                HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(3)*dzNode(3) + vvNode(8)*dzNode(8))/(dzNode(3) + dzNode(8))
!                    
!                                rj(1,3) = ru(uuNode(3),uuNode(8),uuNode(9),yyNode(3),yyNode(8),yyNode(9))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(7),uuNode(3),uuNode(8),yyNode(7),yyNode(3),yyNode(8))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                     
!                                 
!                                rj(2,3) = ru(vvNode(3),vvNode(8),vvNode(9),yyNode(3),yyNode(8),yyNode(9))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!                
!                                rj(2,4) = ru(vvNode(7),vvNode(3),vvNode(8),yyNode(7),yyNode(3),yyNode(8))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))         
!                            Else
!                                HydroParam%uArrow(iLayer,2,iEdge) =  (vvNode(7)*dzNode(7) + vvNode(3)*dzNode(3))/(dzNode(7) + dzNode(3))
!                                
!                                rj(1,3) = ru(uuNode(7),uuNode(3),uuNode(8),yyNode(7),yyNode(3),yyNode(8))
!                                Call Psi_value(psi_flag,rj(1,3),Courant,fi_small,psi(1,3))
!                
!                                rj(1,4) = ru(uuNode(6),uuNode(7),uuNode(3),yyNode(6),yyNode(7),yyNode(3))
!                                Call Psi_value(psi_flag,rj(1,4),Courant,fi_small,psi(1,4))                    
!                    
!                                rj(2,3) = ru(vvNode(7),vvNode(3),vvNode(8),yyNode(7),yyNode(3),yyNode(8))
!                                Call Psi_value(psi_flag,rj(2,3),Courant,fi_small,psi(2,3))
!
!                                rj(2,4) = ru(vvNode(6),vvNode(7),vvNode(3),yyNode(6),yyNode(7),yyNode(3))
!                                Call Psi_value(psi_flag,rj(2,4),Courant,fi_small,psi(2,4))                                                      
!                            EndIf                                                                                                            
!                        EndIf                              
!                    EndDo            
!                EndIf
!                
!                !HydroParam%uxyback(iLayer,1,iEdge) = HydroParam%uxyback(iLayer,1,iEdge) - dt*(0.5*Max(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2+psi(1,2))*uuNode(1,8) - (2 + psi(1,1) + psi(1,2))*uuNode(1,5) + psi(1,1)*uuNode(1,2))/MeshParam%dx + 0.5*Min(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2+psi(1,2))*uuNode(1,3) - (2 + psi(1,1) + psi(1,2))*uuNode(1,8) + psi(1,1)*uuNode(1,5))/MeshParam%dx + 0.5*Max(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2+psi(1,4))*uuNode(1,8) - (2 + psi(1,3) + psi(1,4))*uuNode(1,7) + psi(1,3)*uuNode(1,6))/MeshParam%dy + 0.5*Min(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2+psi(1,4))*uuNode(1,9) - (2 + psi(1,3) + psi(1,4))*uuNode(1,8) + psi(1,3)*uuNode(1,7))/MeshParam%dy)                           
!                !HydroParam%uxyback(iLayer,2,iEdge) = HydroParam%uxyback(iLayer,1,iEdge) - dt*(0.5*Max(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2+psi(2,2))*vvNode(1,8) - (2 + psi(2,1) + psi(2,2))*vvNode(1,5) + psi(2,1)*vvNode(1,2))/MeshParam%dx + 0.5*Min(HydroParam%uArrow(iLayer,1,iEdge),0.d0)*((2+psi(2,2))*vvNode(1,3) - (2 + psi(2,1) + psi(2,2))*vvNode(1,8) + psi(2,1)*vvNode(1,5))/MeshParam%dx  + 0.5*Max(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2+psi(2,4))*vvNode(1,8) - (2 + psi(2,3) + psi(2,4))*vvNode(1,7) + psi(2,3)*vvNode(1,1))/MeshParam%dy + 0.5*Min(HydroParam%uArrow(iLayer,2,iEdge),0.d0)*((2+psi(2,4))*vvNode(1,9) - (2 + psi(2,3) + psi(2,4))*vvNode(1,8) + psi(2,3)*vvNode(1,7))/MeshParam%dy)                
!                !
!            EndDo
!        EndDo
!        
!        !If (HydroParam%iNonHydro==0.and.MeshParam%Kmax>1) Then
!        !    Call ComputeFW(MeshParam,HydroParam,dt)
!        !EndIf
!   
!        Do iEdge=1,MeshParam%nEdge
!            l = MeshParam%Left(iEdge)
!            r = MeshParam%Right(iEdge)
!            Do lEdge=1,4
!                If (MeshParam%Edge(lEdge,l)==iEdge) Then
!                    Exit
!                EndIf
!            EndDo
!            !jump=0
!            !do Face = 1,4
!            ! If the face has a boundary condition:
!            If (HydroParam%IndexInflowEdge(iEdge)>0.or.HydroParam%IndexWaterLevelEdge(iEdge)>0) then
!                HydroParam%Fu(:,iEdge)=HydroParam%u(:,iEdge)
!                cycle
!            EndIf
!            !enddo
!            !if (jump==1) then
!            !    cycle
!            !endif
!        
!            !MeshParam%Neighbor(:,l)
!            Do iLayer = HydroParam%Smallm(iEdge), HydroParam%CapitalM(iEdge)
!                If (r==0) Then !Velocity in lateral boundaries is equal to cell center velocity
!                    HydroParam%Fu(iLayer,iEdge) = HydroParam%u(iLayer,iEdge) !Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*(HydroParam%Fub(iLayer,1,l))*MeshParam%NormalVector(1,lEdge,l)  + Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*(HydroParam%Fub(iLayer,2,l))*MeshParam%NormalVector(2,lEdge,l)
!                Else
!                    if (iLayer < HydroParam%ElSmallm(r)) Then
!                        HydroParam%Fu(iLayer,iEdge) = HydroParam%u(iLayer,iEdge)
!                    else
!                        HydroParam%Fu(iLayer,iEdge) = Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*HydroParam%uxyback(iLayer,1,iEdge)*MeshParam%NormalVector(1,lEdge,l)  + Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*HydroParam%uxyback(iLayer,2,iEdge)*MeshParam%NormalVector(2,lEdge,l)
!                        !HydroParam%Fu(iLayer,iEdge) = Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*HydroParam%uxyback(iLayer,1,iEdge)*MeshParam%NormalVector(1,lEdge,l)  + (Sig(l,MeshParam%Right(iEdge),MeshParam%Left(iEdge))*HydroParam%uxyback(iLayer,2,iEdge)*MeshParam%NormalVector(2,lEdge,l)*dt)*((HydroParam%uxyback(iLayer,2,iEdge-3)-2*HydroParam%uxyback(iLayer,2,iEdge)+HydroParam%uxyback(iLayer,2,iEdge+3))/MeshParam%EdgeLength(iEdge))
!                    endif
!                EndIf
!            EndDo
!        EndDo
!   
!        HydroParam%FuxyNode = 0.
!        Do iNode = 1,MeshParam%nNode
!            Do iLayer = 1,MeshParam%Kmax
!                weit=0
!                icount=0
!            
!                Do j = 1,MeshParam%nEgdesatNode(iNode)
!                    Face = MeshParam%EgdesatNode(j,iNode)
!                    l = MeshParam%Left(Face)
!                    r = MeshParam%Right(Face)
!                    If (r==0) Then
!                        fac=2
!                    Else
!                        fac=1
!                    EndIf
!                    If(iLayer>=HydroParam%Smallm(Face)) Then
!                        kin=min(iLayer,HydroParam%CapitalM(Face))
!                        HydroParam%FuxyNode(kin,1,iNode)=HydroParam%FuxyNode(kin,1,iNode)+HydroParam%uxyback(kin,1,Face)/MeshParam%EdgeLength(Face)*fac
!                        HydroParam%FuxyNode(kin,2,iNode)=HydroParam%FuxyNode(kin,2,iNode)+HydroParam%uxyback(kin,2,Face)/MeshParam%EdgeLength(Face)*fac
!                        icount=icount+1 
!                    EndIf
!                    weit=weit+fac/MeshParam%EdgeLength(Face)
!                EndDo
!                If(icount.ne.0) Then
!                    HydroParam%FuxyNode(iLayer,1,iNode) = HydroParam%FuxyNode(iLayer,1,iNode)/weit
!                    HydroParam%FuxyNode(iLayer,2,iNode) = HydroParam%FuxyNode(iLayer,2,iNode)/weit
!                EndIf
!            
!            EndDo       
!        EndDo
!    
!    Return
!    End Subroutine ConservativeFuFv
!
!    
!    Function ru(u,ux1,ux2,x,x1,x2) 
!    ! Compute the Water Elevation
!    ! Based on:
!    ! [1] 
!    ! Input:
!    ! u -> u velocity, eta or concentration
!    ! ux1  -> u in position -1/2
!    ! ux2 -> u in position -3/2
!    ! x -> u position
!    ! x1 -> u-1/2 position
!    ! x2 -> u-3/2 position
!    ! Output:
!    ! r   -> Ratio of consecutive gradients
!    ! List of Modifications: 
!    !   -> 04.08.2020: Routine Implementation
!    ! Programmer: 
!        
!    Implicit None
!    Real:: ru, u, ux1, ux2, x, x1, x2
!    
!    If(ux1 - ux2 == 0.or.x - x1 == 0) Then
!        ru = 0.d0
!    Else
!        ru = (u - ux1)/(ux1-ux2)*abs(x1-x2)/abs(x-x1)
!    EndIf
!    
!    End Function ru    
!    
!    Subroutine NodeValues(iEdgein,iEdge,Face,rElem,lElem,n7,r,l,uuNode(:),vvNode(:),xxNode(:), yyNode(:), dzNode(:))
!    
!    Use MeshVars
!    Use Hydrodynamic
!    
!    Implicit None
!
!    Real, intent(in) :: uuNode(9), vvNode(9),  xxNode(9), yyNode(9), dzNode(9)
!    Real, intent(out) :: uuNode(9), vvNode(9), xxNode(9), yyNode(9), dzNode(9), Nodes(9)
!    Integer:: iEdgein,iEdge,Face,rElem,lElem,n7,r,l
!    Integer:: n1,n2,n3,nu,nu1,nl,n7,n8,nr
!    type(MeshGridParam) :: MeshParam
!    type(HydrodynamicParam) :: HydroParam
!    
!    
!    n1 = MeshParam%Edge(iEdgein,l)
!    n2 = l
!    n3 = Face
!    n7 = MeshParam%Quadri(n7,l) + 1                        
!    n8 = MeshParam%Quadri(iEdge,l) + 1
!
!    Nodes(1:9)= (/n1, n2, n3, 0, 0, 0, n7, n8, 0 /)               
!            
!    uuNode(1) = HydroParam%uxy(iLayer,1,Nodes(1));     vvNode(1) = HydroParam%uxy(iLayer,2,Nodes(1));
!    uuNode(2) = HydroParam%ub(iLayer,1,Nodes(2));      vvNode(2) = HydroParam%ub(iLayer,2,Nodes(2)); 
!    uuNode(8) = HydroParam%ubV(iLayer,1,Nodes(8));     vvNode(8) = HydroParam%ubV(iLayer,2,Nodes(8));
!    uuNode(3) = HydroParam%uxy(iLayer,1,Nodes(3));     vvNode(3) = HydroParam%uxy(iLayer,2,Nodes(3));
!    uuNode(7) = HydroParam%ubV(iLayer,1,Nodes(7));     vvNode(7) = HydroParam%ubV(iLayer,2,Nodes(7));
!                       
!    xxNode(1) = MeshParam%EdgeBary(Nodes(1)); yyNode(1) = MeshParam%EdgeBary(2,Nodes(1));   dzNode(1) = HydroParam%Z(iLayer+1,Nodes(1)) - HydroParam%Ze(iLayer,l) + sum(HydroParam%DZsi(:,l))
!    xxNode(2) = MeshParam%xb(Nodes(2));         yyNode(2) = MeshParam%yb(Nodes(2));           dzNode(2) = HydroParam%Ze(iLayer+1,l)  - HydroParam%Ze(iLayer,l) + sum(HydroParam%DZsi(:,l))
!    xxNode(8) = MeshParam%xNode(Nodes(8));      yyNode(8) = MeshParam%yNode(Nodes(8));        dzNode(8) = HydroParam%Ze(iLayer+1,l) - HydroParam%Ze(iLayer,l) + sum(HydroParam%DZsi(:,l))
!    xxNode(3) = MeshParam%EdgeBary(Nodes(3)); yyNode(3) = MeshParam%EdgeBary(2,Nodes(3));   dzNode(3) = HydroParam%Z(iLayer+1,Nodes(3)) - HydroParam%Ze(iLayer,l) + sum(HydroParam%DZsi(:,l))   
!    xxNode(7) = MeshParam%xNode(Nodes(7));      yyNode(7) = MeshParam%yNode(Nodes(7));        dzNode(7) = HydroParam%Ze(iLayer+1,l) - HydroParam%Ze(iLayer,l) + sum(HydroParam%DZsi(:,l))
!                
!    ! Right neighbour element:    
!    !nr = n9, nu = n4, nu1 = n5, nl = n6
!    If (rElem==0) Then
!        nr = n8
!        uuNode(9) = uuNode(8);           vvNode(9) = vvNode(8);
!        xxNode(9) =  xxNode(8);                 yyNode(9) = yyNode(8);              dzNode(9) = dzNode(8)                  
!    ElseIf (ilayer<HydroParam%ElSmallms(rElem)) Then
!        nr = n8
!        uuNode(9) = uuNode(8);      vvNode(9) = vvNode(8);
!        xxNode(9) =  xxNode(8);     yyNode(9) = yyNode(8);              dzNode(9) = dzNode(8)                              
!    Else
!        nr   = MeshParam%Edge(iEdge,rElem)
!        uuNode(9) = HydroParam%uxy(iLayer,1,nr);   vvNode(9) = HydroParam%uxy(iLayer,2,nr);
!        xxNode(9) =  MeshParam%EdgeBary(1,nr);     yyNode(9) = MeshParam%EdgeBary(2,nr);    dzNode(9) = HydroParam%Z(iLayer+1,nr) - HydroParam%Ze(iLayer,rElem) + sum(HydroParam%DZsi(:,rElem))              
!    EndIf
!            
!    ! Left neighbour element:                
!    If (lElem==0) Then
!        nl = n7
!        uuNode(6)= uuNode(7);               vvNode(6) = vvNode(7);
!        xxNode(6) =  xxNode(7);             yyNode(6) = yyNode(7);                      dzNode(6) = dzNode(7)                       
!    ElseIf (ilayer<HydroParam%ElSmallms(lElem)) Then
!        nl = n7
!        uuNode(6)= uuNode(7);               vvNode(6) = vvNode(7);
!        xxNode(6) =  xxNode(7);             yyNode(6) = yyNode(7);                dzNode(6) = dzNode(7)                       
!    Else
!        nl    = MeshParam%Edge(iEdge,lElem)
!        uuNode(6) = HydroParam%uxy(iLayer,1,nl);   vvNode(6) = HydroParam%uxy(iLayer,2,nl);
!        xxNode(6) = MeshParam%EdgeBary(1,nl);      yyNode(6) = MeshParam%EdgeBary(2,nl);        dzNode(6) = HydroParam%Z(iLayer+1,nl) - HydroParam%Ze(iLayer,lElem) + sum(HydroParam%DZsi(:,lElem))                  
!    EndIf
!                        
!    ! Upper neighbour:
!    !nr = n9, nu = n4, nu1 = n5, nl = n6
!    nu = r              
!    If (nu == 0) Then
!        nu = l                
!        nu1 = l              
!        uuNode(4) = uuNode(3);             vvNode(4) = vvNode(3);
!        xxNode(4) =  xxNode(3);             yyNode(4) = yyNode(3);                      dzNode(4) = dzNode(3)              
!        uuNode(5) = uuNode(3);              vvNode(5) = vvNode(3);
!        xxNode(5) =  xxNode(3);            yyNode(5) = yyNode(3);                      dzNode(5) = dzNode(3)     
!    ElseIf (ilayer<HydroParam%ElSmallms(nu)) Then
!        nu = l                
!        nu1 = l              
!        uuNode(4) = uuNode(3);             vvNode(4) = vvNode(3);
!        xxNode(4) =  xxNode(3);             yyNode(4) = yyNode(3);                      dzNode(4) = dzNode(3)              
!        uuNode(5) = uuNode(3);              vvNode(5) = vvNode(3);
!        xxNode(5) =  xxNode(3);            yyNode(5) = yyNode(3);                      dzNode(5) = dzNode(3)                 
!    Else
!        nu1 = MeshParam%Edge(iEdge,r)    
!        uuNode(4)  = HydroParam%ub(iLayer,1,nu);             vvNode(4)  = HydroParam%ub(iLayer,2,nu);
!        xxNode(4)  =  MeshParam%xb(nu);                      yyNode(4)  = MeshParam%yb(nu);           dzNode(4) = HydroParam%eta(nu) - HydroParam%Ze(iLayer,nu) + sum(HydroParam%DZsi(:,nu))                         
!        uuNode(5) = HydroParam%uxy(iLayer,1,nu1);           vvNode(5) = HydroParam%uxy(iLayer,2,nu1);  
!        xxNode(5) = MeshParam%EdgeBary(1,nu1);              yyNode(5) = MeshParam%EdgeBary(2,nu1); dzNode(5) = HydroParam%Z(iLayer+1,nu1) - HydroParam%Ze(iLayer,nu) + sum(HydroParam%DZsi(:,nu))            
!    EndIf
!                                     
!    Return
!    End Subroutine NodeValues
!    
!    End Module AdvcConservativeScheme