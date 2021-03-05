	Subroutine VelocitiesSUB(HydroParam,MeshParam)

	! This routine calculates the mean water velocity for each cell
    ! Called in routine 0-MAIN
    
    ! UCEL/VCEL are used to calculate the narrows for flow direction in TELA_AV
    ! UT2M/VT2M are used to highlight areas with higher velocities in TELA_AV
    
	Use MeshVars !, Only: nElem,nEdge,nNode,kMax,NormalVector,TangentVector,nEgdesatNode,EgdesatNode,EdgeLength,Left,Right,Edge,Quadri,xb,yb,Area,xNode,yNode,Edge,nVertexElem,VertexElem,Neighbor
    Use Hydrodynamic
    !Use SimulationModel !, Only:i34,dt,dx,dynEgdesatNode
     
	Implicit none 
    Integer:: r,l,iElem, Elem,iLayer,iEdge,iNode,Sig,Face,icountU, icountV, icountW,kin,j,cond, ncond, lockElem, small
    Integer:: y1Elem,x1Elem,y2Elem,x2Elem,Facexin,Faceyin,neighElem
    Real ::weit, weitW,signa, weitU, weitV, med
    Real:: aa(5), aaU, aaV, lPoint(2), lArea,luNode(2,4), luEdge(2,4), z1, z2, z3, zp, UiQuad
    Integer:: lEdge(4), lTri(4),n1,n2,n3,n4,iel
    Real:: USUL,UNORTE,VOESTE,VLESTE
    Real:: NearZero = -1!1.e-20
    Real:: NearZero2 = -1!1.e-20
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    !type(SimulationParam), intent(in) :: simParam
    Real:: error(MeshParam%Kmax,2,MeshParam%nEdge),uxyFv(MeshParam%Kmax,2,MeshParam%nEdge)
    ! 1. Finding horizontal velocity components in the edges
	!Do iElem = 1, nElem
 !       Do iLayer = ElSmallm(iElem), ElCapitalM(iElem)
 !           Do iEdge = 1,4
 !               uxy(iLayer,1:2,Edge(iEdge,iElem)) = Sig(iElem,Right(Edge(iEdge,iElem)),MeshParam%Left(Edge(iEdge,iElem)))*u(iLayer,Edge(iEdge,iElem))*NormalVector(1:2,iEdge,iElem) + Sig(iElem,Right(Edge(iEdge,iElem)),MeshParam%Left(Edge(iEdge,iElem)))*uTang(iLayer,Edge(iEdge,iElem))*TangentVector(1:2,iEdge,iElem)	
 !           EndDo
 !       EndDo
 !   EndDo
    
    
    HydroParam%uxyt = HydroParam%uxy
    HydroParam%ubt = HydroParam%ub
    HydroParam%uNodet = HydroParam%uNode 
    HydroParam%ugt = HydroParam%ug   
    HydroParam%vgt =  HydroParam%vg 
    HydroParam%wgt = HydroParam%wg 
    HydroParam%ubVt  = HydroParam%ubV
    HydroParam%uxyLt = HydroParam%uxyL
    HydroParam%wfct  = HydroParam%wfc
     
    HydroParam%ubsub=0.d0
    HydroParam%uxysub=0.d0
    
    HydroParam%ubV=0.d0
    HydroParam%uNode=0.d0
    HydroParam%uxyL=0.d0
    HydroParam%ub=0.d0
    HydroParam%wfc=0.d0
    HydroParam%ug=0.d0
    HydroParam%vg=0.d0
    HydroParam%wg=0.d0
    HydroParam%uxy=0.d0
    
	Do iElem = 1, MeshParam%nElem
		lockElem = 0

        Do iEdge = 1,4
            Face = MeshParam%Edge(iEdge,iElem)
            l = MeshParam%Left(Face) 
            r = MeshParam%Right(Face)

            If (r==0) Then
                !If r==0 the face no have neighbour, this implies that the lower layer with velocity is the Edge Smallms:
                Small = HydroParam%Smallm(Face)
            Else
                !In this case, the lower layer with velocity is the min between Edge's elements
                Small = min(HydroParam%ElSmallm(l),HydroParam%ElSmallm(r))
            EndIf        
            
            Do iLayer = Small, HydroParam%ElCapitalM(iElem) !Small
                med=0.5d0
                
                If (r/=0) Then
                    !If Neighbour lower layer is above Elem layer, the lower Element layer is 
                    If (HydroParam%ElSmallm(r)<HydroParam%ElSmallm(l)) Then
                        If (ilayer<HydroParam%ElSmallm(l)) Then
                            l = r
                        Else
                            l = MeshParam%Left(Face)
                            r = MeshParam%Right(Face)
                        EndIf
                    EndIf
                EndIf
                
     !!            !Set particle's Z position
     !           If(iEdge==1) Then
     !               HydroParam%uxy(iLayer,2,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%u(iLayer,Face)                    
     !               If (HydroParam%IndexInflowEdge(Face)>0) Then
					!	HydroParam%uxy(iLayer,1,Face) = 0.
					!Else
     !                   HydroParam%uxy(iLayer,1,Face) =  Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%utang(iLayer,Face)
     !               EndIf
     !           ElseIf(iEdge==2) Then
     !               HydroParam%uxy(iLayer,1,Face) = -Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%u(iLayer,Face)
     !               If (HydroParam%IndexInflowEdge(Face)>0) Then
     !                   HydroParam%uxy(iLayer,2,Face) = 0.d0
				 !   Else
     !                   HydroParam%uxy(iLayer,2,Face) = -Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%utang(iLayer,Face)
     !               EndIf                      
     !           ElseIf(iEdge==3) Then
     !               HydroParam%uxy(iLayer,2,Face) = -Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%u(iLayer,Face)
     !               If (HydroParam%IndexInflowEdge(Face)>0) Then
					!    HydroParam%uxy(iLayer,1,Face) = 0.d0
				 !   Else
     !                   HydroParam%uxy(iLayer,1,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%utang(iLayer,Face)
     !               EndIf                        
     !           Else 
     !               HydroParam%uxy(iLayer,1,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%u(iLayer,Face)
     !               If (HydroParam%IndexInflowEdge(Face)>0) Then
     !                   HydroParam%uxy(iLayer,2,Face) = 0.d0
				 !   Else
     !                   HydroParam%uxy(iLayer,2,Face) = -Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%utang(iLayer,Face)
     !               EndIf   
     !           EndIf
                
     !           If(iEdge==1) Then
     !               uxyFv(iLayer,2,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%u(iLayer,Face)                    
     !               If (HydroParam%IndexInflowEdge(Face)>0) Then
					!	uxyFv(iLayer,1,Face) = 0.
					!Else
     !                   uxyFv(iLayer,1,Face) =  Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%utang(iLayer,Face)
     !               EndIf
     !           ElseIf(iEdge==2) Then
     !               uxyFv(iLayer,1,Face) = -Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%u(iLayer,Face)
     !               If (HydroParam%IndexInflowEdge(Face)>0) Then
     !                   uxyFv(iLayer,2,Face) = 0.d0
				 !   Else
     !                   uxyFv(iLayer,2,Face) = -Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%utang(iLayer,Face)
     !               EndIf                      
     !           ElseIf(iEdge==3) Then
     !               uxyFv(iLayer,2,Face) = -Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%u(iLayer,Face)
     !               If (HydroParam%IndexInflowEdge(Face)>0) Then
					!    uxyFv(iLayer,1,Face) = 0.d0
				 !   Else
     !                   uxyFv(iLayer,1,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%utang(iLayer,Face)
     !               EndIf                        
     !           Else 
     !               uxyFv(iLayer,1,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%u(iLayer,Face)
     !               If (HydroParam%IndexInflowEdge(Face)>0) Then
     !                   uxyFv(iLayer,2,Face) = 0.d0
				 !   Else
     !                   uxyFv(iLayer,2,Face) = -Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%utang(iLayer,Face)
     !               EndIf   
     !           EndIf
                If (HydroParam%H(iEdge)+HydroParam%sj(iEdge)-HydroParam%hj(iEdge) > HydroParam%PCRI+NearZero) Then            
                    If (iEdge==1) Then !North Edge
                        If (r == 0.or.lockElem==1) Then
						    If (HydroParam%IndexInflowEdge(Face)>0.or.HydroParam%IndexWaterLevelEdge(Face)>0) Then
							    HydroParam%uxy(iLayer,1,Face) = 0. !
						    Else
							    If (abs (HydroParam%u(iLayer,MeshParam%Edge(4,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(2,iElem)))<NearZero2) then
								    HydroParam%uxy(iLayer,1,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(4,iElem))
							    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(2,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(4,iElem)))<NearZero2) then
								    HydroParam%uxy(iLayer,1,Face) = - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(2,iElem))
							    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(4,iElem)))>NearZero2.and.abs (HydroParam%u(iLayer,MeshParam%Edge(2,iElem)))>NearZero2) then
								    HydroParam%uxy(iLayer,1,Face) = 0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(2,iElem)) + Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(4,iElem)) )
							    Else
								    HydroParam%uxy(iLayer,1,Face) = 0.0d0 !
							    EndIf
                            EndIf
                        Else
						    Elem=MeshParam%Neighbor(1,iElem)
						    If (abs (HydroParam%u(iLayer,MeshParam%Edge(4,Elem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(2,Elem)))<NearZero2) then
							    UNORTE = Sig(Elem,MeshParam%Right(MeshParam%Edge(4,Elem)),MeshParam%Left(MeshParam%Edge(4,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(4,Elem))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(2,Elem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(4,Elem)))<NearZero2) then
							    UNORTE = - Sig(Elem,MeshParam%Right(MeshParam%Edge(2,Elem)),MeshParam%Left(MeshParam%Edge(2,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(2,Elem))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(4,Elem)))>NearZero2.and.abs (HydroParam%u(iLayer,MeshParam%Edge(2,Elem)))>NearZero2) then
							    UNORTE = 0.5d0*( - Sig(Elem,MeshParam%Right(MeshParam%Edge(2,Elem)),MeshParam%Left(MeshParam%Edge(2,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(2,Elem)) + Sig(Elem,MeshParam%Right(MeshParam%Edge(4,Elem)),MeshParam%Left(MeshParam%Edge(4,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(4,Elem)) )
						    Else
							    UNORTE = 0. !
                                med=1
						    EndIf
						
						    If (abs (HydroParam%u(iLayer,MeshParam%Edge(4,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(2,iElem)))<NearZero2) then
							    HydroParam%uxy(iLayer,1,Face) = med*(UNORTE + (Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(4,iElem))))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(2,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(4,iElem)))<NearZero2) then
							    HydroParam%uxy(iLayer,1,Face) = med*(UNORTE + (- Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(2,iElem))))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(4,iElem)))>NearZero2.and.abs (HydroParam%u(iLayer,MeshParam%Edge(2,iElem)))>NearZero2) then
							    HydroParam%uxy(iLayer,1,Face) = med*(UNORTE + (0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(2,iElem)) + Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(4,iElem)) )))
						    Else
							    HydroParam%uxy(iLayer,1,Face) = UNORTE !
						    EndIf
      
                        EndIf
                        ! Set particle's tangential velocity
                        !uxy(iLayer,1,Face) = 0.25*( - Sig(iElem,Right(Edge(2,iElem)),Left(Edge(2,iElem)))*u(iLayer,Edge(2,iElem)) + 2.*UNORTE + Sig(iElem,Right(Edge(4,iElem)),Left(Edge(4,iElem)))*u(iLayer,Edge(4,iElem)) )
                        HydroParam%uxy(iLayer,2,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(1,iElem)) 
                    ElseIf (iEdge==2) then !West Edge
                        If (r == 0.or.lockElem==1) Then
						    If (HydroParam%IndexInflowEdge(Face)>0.or.HydroParam%IndexWaterLevelEdge(Face)>0) then
							    HydroParam%uxy(iLayer,2,Face) = 0. !
						    Else
							    If (abs (HydroParam%u(iLayer,MeshParam%Edge(1,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(3,iElem)))<NearZero2) then
								    HydroParam%uxy(iLayer,2,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(1,iElem))
							    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(3,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(1,iElem)))<NearZero2) then
								    HydroParam%uxy(iLayer,2,Face) = - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(3,iElem))
							    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(1,iElem)))>NearZero2.and.abs (HydroParam%u(iLayer,MeshParam%Edge(3,iElem)))>NearZero2) then
								    HydroParam%uxy(iLayer,2,Face) = 0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(3,iElem))+ Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(1,iElem)) )
							    Else
								    HydroParam%uxy(iLayer,2,Face) = 0. !
							    EndIf
						    EndIf
                    
                        Else
						    Elem=MeshParam%Neighbor(2,iElem)
						    If (abs (HydroParam%u(iLayer,MeshParam%Edge(1,Elem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(3,Elem)))<NearZero2) then
							    VOESTE = Sig(Elem,MeshParam%Right(MeshParam%Edge(1,Elem)),MeshParam%Left(MeshParam%Edge(1,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(1,Elem))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(3,Elem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(1,Elem)))<NearZero2) then
							    VOESTE = - Sig(Elem,MeshParam%Right(MeshParam%Edge(3,Elem)),MeshParam%Left(MeshParam%Edge(3,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(3,Elem))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(1,Elem)))>NearZero2.and.abs (HydroParam%u(iLayer,MeshParam%Edge(3,Elem)))>NearZero2) then
							    VOESTE = 0.5d0*( - Sig(Elem,MeshParam%Right(MeshParam%Edge(3,Elem)),MeshParam%Left(MeshParam%Edge(3,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(3,Elem))+ Sig(Elem,MeshParam%Right(MeshParam%Edge(1,Elem)),MeshParam%Left(MeshParam%Edge(1,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(1,Elem)) )
						    Else
							    VOESTE = 0. !
                                med=1
						    EndIf
					 
						    If (abs (HydroParam%u(iLayer,MeshParam%Edge(1,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(3,iElem)))<NearZero2) then
							    HydroParam%uxy(iLayer,2,Face) = med* (VOESTE +(Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(1,iElem))))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(3,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(1,iElem)))<NearZero2) then
							    HydroParam%uxy(iLayer,2,Face) = med* (VOESTE +(- Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(3,iElem))))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(1,iElem)))>NearZero2.and.abs (HydroParam%u(iLayer,MeshParam%Edge(3,iElem)))>NearZero2) then
							    HydroParam%uxy(iLayer,2,Face) = med* (VOESTE +(0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(3,iElem))+ Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(1,iElem)) )))
						    Else
							    HydroParam%uxy(iLayer,2,Face) = VOESTE !
						    EndIf
        
                        EndIf
                        ! Set particle's tangential velocity
                        HydroParam%uxy(iLayer,1,Face) = - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(2,iElem))
                        !uxy(iLayer,2,Face) = 0.25*( - Sig(iElem,MeshParam%Right(Edge(3,iElem)),MeshParam%Left(Edge(3,iElem)))*u(iLayer,Edge(3,iElem)) + Sig(iElem,MeshParam%Right(Edge(1,iElem)),MeshParam%Left(Edge(1,iElem)))*u(iLayer,Edge(1,iElem)) + 2.*VOESTE )
                    ElseIf (iEdge==3) then !South Edge
                        If (r == 0.or.lockElem==1) Then
						    If (HydroParam%IndexInflowEdge(Face)>0.or.HydroParam%IndexWaterLevelEdge(Face)>0) then
							    HydroParam%uxy(iLayer,1,Face) = 0. !
						    Else
							    If (abs (HydroParam%u(iLayer,MeshParam%Edge(4,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(2,iElem)))<NearZero2) then
								    HydroParam%uxy(iLayer,1,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(4,iElem))
							    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(2,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(4,iElem)))<NearZero2) then
								    HydroParam%uxy(iLayer,1,Face) = - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(2,iElem))
							    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(4,iElem)))>NearZero2.and.abs (HydroParam%u(iLayer,MeshParam%Edge(2,iElem)))>NearZero2) then
								    HydroParam%uxy(iLayer,1,Face) = 0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(2,iElem)) + Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(4,iElem)) )
							    Else
								    HydroParam%uxy(iLayer,1,Face) = 0. !
							    EndIf
						    EndIf
					        !USUL = 0.5*( - Sig(iElem,MeshParam%Right(Edge(2,iElem)),MeshParam%Left(Edge(2,iElem)))*u(iLayer,Edge(2,iElem)) + Sig(iElem,MeshParam%Right(Edge(4,iElem)),MeshParam%Left(Edge(4,iElem)))*u(iLayer,Edge(4,iElem)) )
                        Else 
						    Elem=MeshParam%Neighbor(3,iElem)
                    
						    If (abs (HydroParam%u(iLayer,MeshParam%Edge(4,Elem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(2,Elem)))<NearZero2) then
							    USUL = Sig(Elem,MeshParam%Right(MeshParam%Edge(4,Elem)),MeshParam%Left(MeshParam%Edge(4,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(4,Elem))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(2,Elem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(4,Elem)))<NearZero2) then
							    USUL = - Sig(Elem,MeshParam%Right(MeshParam%Edge(2,Elem)),MeshParam%Left(MeshParam%Edge(2,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(2,Elem))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(4,Elem)))>NearZero2.and.abs (HydroParam%u(iLayer,MeshParam%Edge(2,Elem)))>NearZero2) then
							    USUL = 0.5d0*( - Sig(Elem,MeshParam%Right(MeshParam%Edge(2,Elem)),MeshParam%Left(MeshParam%Edge(2,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(2,Elem)) + Sig(Elem,MeshParam%Right(MeshParam%Edge(4,Elem)),MeshParam%Left(MeshParam%Edge(4,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(4,Elem)) )
						    Else
							    USUL = 0. !
                                med=1
						    EndIf
					 
						    If (abs (HydroParam%u(iLayer,MeshParam%Edge(4,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(2,iElem)))<NearZero2) then
							    HydroParam%uxy(iLayer,1,Face) = med*(USUL + (Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(4,iElem))))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(2,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(4,iElem)))<NearZero2) then
							    HydroParam%uxy(iLayer,1,Face) = med*(USUL + (- Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(2,iElem))))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(4,iElem)))>NearZero2.and.abs (HydroParam%u(iLayer,MeshParam%Edge(2,iElem)))>NearZero2) then
							    HydroParam%uxy(iLayer,1,Face) = med*(USUL + (0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(2,iElem)) + Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(4,iElem)) )))
						    Else
							    HydroParam%uxy(iLayer,1,Face) = USUL !
						    EndIf
					 
                        EndIf
                        ! Set particle's tangential velocity
                        !uxy(iLayer,1,Face) =0.25*( - Sig(iElem,MeshParam%Right(Edge(2,iElem)),MeshParam%Left(Edge(2,iElem)))*u(iLayer,Edge(2,iElem)) + 2.*USUL      + Sig(iElem,MeshParam%Right(Edge(4,iElem)),MeshParam%Left(Edge(4,iElem)))*u(iLayer,Edge(4,iElem)) )
                        HydroParam%uxy(iLayer,2,Face) = - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(3,iElem))
                    Else  !East Edge
                        If (r == 0.or.lockElem==1) Then
						    If (HydroParam%IndexInflowEdge(Face)>0.or.HydroParam%IndexWaterLevelEdge(Face)>0) Then
							    HydroParam%uxy(iLayer,2,Face) = 0. !
						    Else
							    If (abs(HydroParam%u(iLayer,MeshParam%Edge(1,iElem)))>NearZero2 .and. abs(HydroParam%u(iLayer,MeshParam%Edge(3,iElem)))<NearZero2) Then
								    HydroParam%uxy(iLayer,2,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(1,iElem))
							    ElseIf (abs(HydroParam%u(iLayer,MeshParam%Edge(3,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(1,iElem)))<NearZero2) Then
								    HydroParam%uxy(iLayer,2,Face) = - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(3,iElem))
							    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(1,iElem)))>NearZero2.and.abs (HydroParam%u(iLayer,MeshParam%Edge(3,iElem)))>NearZero2) Then
								    HydroParam%uxy(iLayer,2,Face) = 0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(3,iElem))+ Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(1,iElem)) )
							    Else
								    HydroParam%uxy(iLayer,2,Face) = 0. !
							    EndIf
						    EndIf
							    !VLESTE = 0.5*( - Sig(iElem,MeshParam%Right(Edge(3,iElem)),MeshParam%Left(Edge(3,iElem)))*u(iLayer,Edge(3,iElem))+ Sig(iElem,MeshParam%Right(Edge(1,iElem)),MeshParam%Left(Edge(1,iElem)))*u(iLayer,Edge(1,iElem)) )
                        Else
						    Elem=MeshParam%Neighbor(4,iElem)					
						    If (abs (HydroParam%u(iLayer,MeshParam%Edge(1,Elem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(3,Elem)))<NearZero2) Then
							    VLESTE = Sig(Elem,MeshParam%Right(MeshParam%Edge(1,Elem)),MeshParam%Left(MeshParam%Edge(1,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(1,Elem))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(3,Elem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(1,Elem)))<NearZero2) Then
							    VLESTE = - Sig(Elem,MeshParam%Right(MeshParam%Edge(3,Elem)),MeshParam%Left(MeshParam%Edge(3,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(3,Elem))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(1,Elem)))>NearZero2.and.abs (HydroParam%u(iLayer,MeshParam%Edge(3,Elem)))>NearZero2) Then
							    VLESTE = 0.5d0*( - Sig(Elem,MeshParam%Right(MeshParam%Edge(3,Elem)),MeshParam%Left(MeshParam%Edge(3,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(3,Elem))+ Sig(Elem,MeshParam%Right(MeshParam%Edge(1,Elem)),MeshParam%Left(MeshParam%Edge(1,Elem)))*HydroParam%u(iLayer,MeshParam%Edge(1,Elem)) )
						    Else
							    VLESTE = 0. !
                                med=1
						    EndIf
		    
						    If (abs (HydroParam%u(iLayer,MeshParam%Edge(1,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(3,iElem)))<NearZero2) Then
							    HydroParam%uxy(iLayer,2,Face) = med* (VLESTE +(Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(1,iElem))))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(3,iElem)))>NearZero2 .and. abs (HydroParam%u(iLayer,MeshParam%Edge(1,iElem)))<NearZero2) Then
							    HydroParam%uxy(iLayer,2,Face) = med* (VLESTE +(- Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(3,iElem))))
						    ElseIf (abs (HydroParam%u(iLayer,MeshParam%Edge(1,iElem)))>NearZero2.and.abs (HydroParam%u(iLayer,MeshParam%Edge(3,iElem)))>NearZero2) Then
							    HydroParam%uxy(iLayer,2,Face) = med* (VLESTE +(0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(3,iElem))+ Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(1,iElem)) )))
						    Else
							    HydroParam%uxy(iLayer,2,Face) = VLESTE !
						    EndIf
                        
                        EndIf
                        ! Set particle's tangential velocity
                        HydroParam%uxy(iLayer,1,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(4,iElem))
                        !uxy(iLayer,2,Face) =  0.25*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(Edge(3,iElem)))*u(iLayer,Edge(3,iElem)) + Sig(iElem,MeshParam%Right(Edge(1,iElem)),MeshParam%Left(Edge(1,iElem)))*u(iLayer,Edge(1,iElem)) + 2.*VLESTE )
                    EndIf
                EndIf
            EndDo
      !      
            ! 3.1 Copy Velocities Above the Free-Surface (du/dz=0) 
            Do iLayer = HydroParam%ElCapitalM(iElem) + 1, MeshParam%KMAX
                HydroParam%uxy(iLayer,1,Face) = 0. !u(CapitalM(Face),Face) 
                HydroParam%uxy(iLayer,2,Face) = 0.
                !Fu(iLayer,Face) = 0.
            EndDo
            ! 3.2 Nullify Velocities Below the Bottom (u=0) 
            Do iLayer = 1, HydroParam%ElSmallms(iElem) - 1
                HydroParam%uxy(iLayer,1,Face) = 0. !u(CapitalM(Face),Face) 
                HydroParam%uxy(iLayer,2,Face) = 0.
                !Fu(iLayer,Face) = 0.
            EndDo
            
            !! 3.1 Copy Velocities Above the Free-Surface (du/dz=0) 
            !Do iLayer = HydroParam%ElCapitalM(iElem) + 1, MeshParam%KMAX
            !    uxyFv(iLayer,1,Face) = 0. !u(CapitalM(Face),Face) 
            !    uxyFv(iLayer,2,Face) = 0.
            !    !Fu(iLayer,Face) = 0.
            !EndDo
            !! 3.2 Nullify Velocities Below the Bottom (u=0) 
            !Do iLayer = 1, HydroParam%ElSmallms(iElem) - 1
            !    uxyFv(iLayer,1,Face) = 0. !u(CapitalM(Face),Face) 
            !    uxyFv(iLayer,2,Face) = 0.
            !    !Fu(iLayer,Face) = 0.
            !EndDo
            
           !Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem) 
           !     error(iLayer,1,Face) = HydroParam%uxy(iLayer,1,Face) - uxyFv(iLayer,1,Face)
           !     error(iLayer,2,Face) = HydroParam%uxy(iLayer,2,Face) - uxyFv(iLayer,2,Face)
           ! EndDo
            ! If subsurface flow occurs, uxysub are calculated:
            If (MeshParam%iBedrock == 1) Then
                
                If (r==0) Then
                    !If r==0 the face no have neighbour, this implies that the lower layer with velocity is the Edge Smallms:
                    Small = HydroParam%Smallms(Face)
                Else
                    !In this case, the lower layer with velocity is the min between Edge's elements
                    Small = min(HydroParam%ElSmallms(l),HydroParam%ElSmallms(r))
                EndIf        
            
                Do iLayer = Small, HydroParam%ElCapitalM(iElem) !Small
                    med=0.5d0
                
                    If (r/=0) Then
                        !If Neighbour lower layer is above Elem layer, the lower Element layer is 
                        If (HydroParam%ElSmallms(r)<HydroParam%ElSmallms(l)) Then
                            If (ilayer<HydroParam%ElSmallms(l)) Then
                                l = r
                            Else
                                l = MeshParam%Left(Face)
                                r = MeshParam%Right(Face)
                            EndIf
                        EndIf
                    EndIf                

                    ! Set particle's Z position
                    If (iEdge==1) Then !North Edge
                        If (r == 0.or.lockElem==1) Then
					        If (HydroParam%IndexInflowEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) Then 
						        HydroParam%uxysub(iLayer,1,Face) = 0. !
					        Else
						        If (abs (HydroParam%um(iLayer,MeshParam%Edge(4,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(2,iElem)))<NearZero2) then
                                    HydroParam%uxysub(iLayer,1,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(4,iElem))
                                ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(2,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(4,iElem)))<NearZero2) then
                                    HydroParam%uxysub(iLayer,1,Face) = - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(2,iElem))
                                ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(4,iElem)))>NearZero2.and.abs (HydroParam%um(iLayer,MeshParam%Edge(2,iElem)))>NearZero2) then
                                    HydroParam%uxysub(iLayer,1,Face) = 0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(2,iElem)) + Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(4,iElem)) )
                                Else
							        HydroParam%uxysub(iLayer,1,Face) = 0. !
						        EndIf
					        EndIf
                        Else
					        Elem=MeshParam%Neighbor(1,iElem)
					        If (abs (HydroParam%um(iLayer,MeshParam%Edge(4,Elem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(2,Elem)))<NearZero2) then
						        UNORTE = Sig(Elem,MeshParam%Right(MeshParam%Edge(4,Elem)),MeshParam%Left(MeshParam%Edge(4,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(4,Elem))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(2,Elem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(4,Elem)))<NearZero2) then
						        UNORTE = - Sig(Elem,MeshParam%Right(MeshParam%Edge(2,Elem)),MeshParam%Left(MeshParam%Edge(2,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(2,Elem))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(4,Elem)))>NearZero2.and.abs (HydroParam%um(iLayer,MeshParam%Edge(2,Elem)))>NearZero2) then
						        UNORTE = 0.5d0*( - Sig(Elem,MeshParam%Right(MeshParam%Edge(2,Elem)),MeshParam%Left(MeshParam%Edge(2,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(2,Elem)) + Sig(Elem,MeshParam%Right(MeshParam%Edge(4,Elem)),MeshParam%Left(MeshParam%Edge(4,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(4,Elem)) )
					        Else
						        UNORTE = 0. !
                                med=1
					        EndIf
						
					        If (abs (HydroParam%um(iLayer,MeshParam%Edge(4,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(2,iElem)))<NearZero2) then
						        HydroParam%uxysub(iLayer,1,Face) = med*(UNORTE + (Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(4,iElem))))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(2,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(4,iElem)))<NearZero2) then
						        HydroParam%uxysub(iLayer,1,Face) = med*(UNORTE + (- Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(2,iElem))))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(4,iElem)))>NearZero2.and.abs (HydroParam%um(iLayer,MeshParam%Edge(2,iElem)))>NearZero2) then
						        HydroParam%uxysub(iLayer,1,Face) = med*(UNORTE + (0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(2,iElem)) + Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(4,iElem)) )))
					        Else
						        HydroParam%uxysub(iLayer,1,Face) = UNORTE !
					        EndIf

                        EndIf
                        ! Set particle's tangential velocity
                        HydroParam%uxysub(iLayer,2,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(1,iElem)) 
                    ElseIf (iEdge==2) then !West Edge
                        If (r == 0.or.lockElem==1) Then
					        If (HydroParam%IndexInflowEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) then
						        HydroParam%uxysub(iLayer,2,Face) = 0. !
					        Else
						        If (abs (HydroParam%um(iLayer,MeshParam%Edge(1,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(3,iElem)))<NearZero2) then
							        HydroParam%uxysub(iLayer,2,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(1,iElem))
						        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(3,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(1,iElem)))<NearZero2) then
							        HydroParam%uxysub(iLayer,2,Face) = - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(3,iElem))
						        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(1,iElem)))>NearZero2.and.abs (HydroParam%um(iLayer,MeshParam%Edge(3,iElem)))>NearZero2) then
							        HydroParam%uxysub(iLayer,2,Face) = 0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(3,iElem))+ Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(1,iElem)) )
						        Else
							        HydroParam%uxysub(iLayer,2,Face) = 0. !
						        EndIf
					        EndIf            
                        Else
					        Elem=MeshParam%Neighbor(2,iElem)
					        If (abs (HydroParam%um(iLayer,MeshParam%Edge(1,Elem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(3,Elem)))<NearZero2) then
						        VOESTE = Sig(Elem,MeshParam%Right(MeshParam%Edge(1,Elem)),MeshParam%Left(MeshParam%Edge(1,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(1,Elem))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(3,Elem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(1,Elem)))<NearZero2) then
						        VOESTE = - Sig(Elem,MeshParam%Right(MeshParam%Edge(3,Elem)),MeshParam%Left(MeshParam%Edge(3,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(3,Elem))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(1,Elem)))>NearZero2.and.abs (HydroParam%um(iLayer,MeshParam%Edge(3,Elem)))>NearZero2) then
						        VOESTE = 0.5d0*( - Sig(Elem,MeshParam%Right(MeshParam%Edge(3,Elem)),MeshParam%Left(MeshParam%Edge(3,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(3,Elem))+ Sig(Elem,MeshParam%Right(MeshParam%Edge(1,Elem)),MeshParam%Left(MeshParam%Edge(1,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(1,Elem)) )
					        Else
						        VOESTE = 0. !
                                med=1
					        EndIf
					
					        If (abs (HydroParam%um(iLayer,MeshParam%Edge(1,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(3,iElem)))<NearZero2) then
						        HydroParam%uxysub(iLayer,2,Face) = med* (VOESTE +(Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(1,iElem))))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(3,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(1,iElem)))<NearZero2) then
						        HydroParam%uxysub(iLayer,2,Face) = med* (VOESTE +(- Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(3,iElem))))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(1,iElem)))>NearZero2.and.abs (HydroParam%um(iLayer,MeshParam%Edge(3,iElem)))>NearZero2) then
						        HydroParam%uxysub(iLayer,2,Face) = med* (VOESTE +(0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(3,iElem))+ Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(1,iElem)) )))
					        Else
						        HydroParam%uxysub(iLayer,2,Face) = VOESTE !
					        EndIf
        
                        EndIf
                        ! Set particle's tangential velocity
                        HydroParam%uxysub(iLayer,1,Face) = - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(2,iElem))
                    ElseIf (iEdge==3) then !South Edge
                        If (r == 0.or.lockElem==1) Then
					        If (HydroParam%IndexInflowEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) then
						        HydroParam%uxysub(iLayer,1,Face) = 0. !
					        Else
						        If (abs (HydroParam%um(iLayer,MeshParam%Edge(4,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(2,iElem)))<NearZero2) then
							        HydroParam%uxysub(iLayer,1,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(4,iElem))
						        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(2,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(4,iElem)))<NearZero2) then
							        HydroParam%uxysub(iLayer,1,Face) = - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(2,iElem))
						        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(4,iElem)))>NearZero2.and.abs (HydroParam%um(iLayer,MeshParam%Edge(2,iElem)))>NearZero2) then
							        HydroParam%uxysub(iLayer,1,Face) = 0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(2,iElem)) + Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(4,iElem)) )
						        Else
							        HydroParam%uxysub(iLayer,1,Face) = 0. !
						        EndIf
					        EndIf
                        Else 
					        Elem=MeshParam%Neighbor(3,iElem)
					        If (abs (HydroParam%um(iLayer,MeshParam%Edge(4,Elem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(2,Elem)))<NearZero2) then
						        USUL = Sig(Elem,MeshParam%Right(MeshParam%Edge(4,Elem)),MeshParam%Left(MeshParam%Edge(4,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(4,Elem))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(2,Elem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(4,Elem)))<NearZero2) then
						        USUL = - Sig(Elem,MeshParam%Right(MeshParam%Edge(2,Elem)),MeshParam%Left(MeshParam%Edge(2,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(2,Elem))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(4,Elem)))>NearZero2.and.abs (HydroParam%um(iLayer,MeshParam%Edge(2,Elem)))>NearZero2) then
						        USUL = 0.5d0*( - Sig(Elem,MeshParam%Right(MeshParam%Edge(2,Elem)),MeshParam%Left(MeshParam%Edge(2,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(2,Elem)) + Sig(Elem,MeshParam%Right(MeshParam%Edge(4,Elem)),MeshParam%Left(MeshParam%Edge(4,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(4,Elem)) )
					        Else
						        USUL = 0. !
                                med=1
					        EndIf
					
					        If (abs (HydroParam%um(iLayer,MeshParam%Edge(4,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(2,iElem)))<NearZero2) then
						        HydroParam%uxysub(iLayer,1,Face) = med*(USUL + (Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(4,iElem))))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(2,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(4,iElem)))<NearZero2) then
						        HydroParam%uxysub(iLayer,1,Face) = med*(USUL + (- Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(2,iElem))))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(4,iElem)))>NearZero2.and.abs (HydroParam%um(iLayer,MeshParam%Edge(2,iElem)))>NearZero2) then
						        HydroParam%uxysub(iLayer,1,Face) = med*(USUL + (0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(2,iElem)) + Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(4,iElem)) )))
					        Else
						        HydroParam%uxysub(iLayer,1,Face) = USUL !
					        EndIf
					
                        EndIf
                        ! Set particle's tangential velocity
                        HydroParam%uxysub(iLayer,2,Face) = - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(3,iElem))
                    Else  !East Edge
                        If (r == 0.or.lockElem==1) Then
					        If (HydroParam%IndexInflowEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) Then
						        HydroParam%uxysub(iLayer,2,Face) = 0. !
					        Else
						        If (abs(HydroParam%um(iLayer,MeshParam%Edge(1,iElem)))>NearZero2 .and. abs(HydroParam%um(iLayer,MeshParam%Edge(3,iElem)))<NearZero2) Then
							        HydroParam%uxysub(iLayer,2,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(1,iElem))
						        ElseIf (abs(HydroParam%um(iLayer,MeshParam%Edge(3,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(1,iElem)))<NearZero2) Then
							        HydroParam%uxysub(iLayer,2,Face) = - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(3,iElem))
						        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(1,iElem)))>NearZero2.and.abs (HydroParam%um(iLayer,MeshParam%Edge(3,iElem)))>NearZero2) Then
							        HydroParam%uxysub(iLayer,2,Face) = 0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(3,iElem))+ Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(1,iElem)) )
						        Else
							        HydroParam%uxysub(iLayer,2,Face) = 0. !
						        EndIf
					        EndIf
                        Else
					        Elem=MeshParam%Neighbor(4,iElem)					
					        If (abs (HydroParam%um(iLayer,MeshParam%Edge(1,Elem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(3,Elem)))<NearZero2) Then
						        VLESTE = Sig(Elem,MeshParam%Right(MeshParam%Edge(1,Elem)),MeshParam%Left(MeshParam%Edge(1,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(1,Elem))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(3,Elem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(1,Elem)))<NearZero2) Then
						        VLESTE = - Sig(Elem,MeshParam%Right(MeshParam%Edge(3,Elem)),MeshParam%Left(MeshParam%Edge(3,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(3,Elem))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(1,Elem)))>NearZero2.and.abs (HydroParam%um(iLayer,MeshParam%Edge(3,Elem)))>NearZero2) Then
						        VLESTE = 0.5d0*( - Sig(Elem,MeshParam%Right(MeshParam%Edge(3,Elem)),MeshParam%Left(MeshParam%Edge(3,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(3,Elem))+ Sig(Elem,MeshParam%Right(MeshParam%Edge(1,Elem)),MeshParam%Left(MeshParam%Edge(1,Elem)))*HydroParam%um(iLayer,MeshParam%Edge(1,Elem)) )
					        Else
						        VLESTE = 0. !
                                med=1
					        EndIf
		
					        If (abs (HydroParam%um(iLayer,MeshParam%Edge(1,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(3,iElem)))<NearZero2) Then
						        HydroParam%uxysub(iLayer,2,Face) = med* (VLESTE +(Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(1,iElem))))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(3,iElem)))>NearZero2 .and. abs (HydroParam%um(iLayer,MeshParam%Edge(1,iElem)))<NearZero2) Then
						        HydroParam%uxysub(iLayer,2,Face) = med* (VLESTE +(- Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(3,iElem))))
					        ElseIf (abs (HydroParam%um(iLayer,MeshParam%Edge(1,iElem)))>NearZero2.and.abs (HydroParam%um(iLayer,MeshParam%Edge(3,iElem)))>NearZero2) Then
						        HydroParam%uxysub(iLayer,2,Face) = med* (VLESTE +(0.5d0*( - Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(3,iElem))+ Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(1,iElem)) )))
					        Else
						        HydroParam%uxysub(iLayer,2,Face) = VLESTE !
					        EndIf
                        
                        EndIf
                        ! Set particle's tangential velocity
                        HydroParam%uxysub(iLayer,1,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*HydroParam%um(iLayer,MeshParam%Edge(4,iElem))
                    EndIf
                EndDo           
                ! 3.2 Nullify Velocities Below the Bottom (u=0) 
                Do iLayer = 1, HydroParam%ElSmallms(iElem) - 1
                    HydroParam%uxysub(iLayer,1,Face) = 0. !u(CapitalM(Face),Face) 
                    HydroParam%uxysub(iLayer,2,Face) = 0.
                EndDo
            EndIf       
            
        EndDo
        !If (iElem==5) Then
        !    write(*,*) HydroParam%uxy(1,1,MeshParam%Edge(2,iElem)),HydroParam%uxy(1,1,MeshParam%Edge(4,iElem))
        !    Continue
        !EndIf        
    EndDo
    !
    !Do iElem = 1, MeshParam%nElem
    !    
    !    Do iEdge = 1,4
    !        y1Elem = neighElem(MeshParam, MeshParam%Edge(1,iElem))
    !        x1Elem = neighElem(MeshParam, MeshParam%Edge(2,iElem))            
    !        y2Elem = neighElem(MeshParam, MeshParam%Edge(3,iElem))
    !        x2Elem = neighElem(MeshParam, MeshParam%Edge(4,iElem))
    !    EndDo        
    !    
    !    Do iEdge = 1,4
    !
    !        Face = MeshParam%Edge(iEdge,iElem)
    !        l = MeshParam%Left(Face) 
    !        r = MeshParam%Right(Face) 
    !
    !        If (r==0) Then
    !            r = l
    !        EndIf
    !        
    !        Do iLayer = HydroParam%Smallms(Face), HydroParam%Capitalms(Face)
    !            
    !            If (iEdge==1) Then
    !                
    !                Faceyin = MeshParam%Edge(4,iElem)                       
    !                Facexin = MeshParam%Edge(iEdge,x1Elem)
    !                
    !                If (iLayer < HydroParam%Smallm(x1Elem)) Then
    !                    Facexin = Face
    !                Else
    !                    Facexin = MeshParam%Edge(iEdge,x1Elem)
    !                EndIf                    
    !   
    !                If (HydroParam%uxy(iLayer,1,Face) + HydroParam%uxy(iLayer,1,Facexin)>=0) Then                      
    !                    HydroParam%uxy(iLayer,1,Face) = (HydroParam%uxy(iLayer,1,Face)*HydroParam%H(MeshParam%Edge(4,iElem)) + HydroParam%uxy(iLayer,1,Facexin)*HydroParam%H(MeshParam%Edge(2,iElem)))/(HydroParam%H(MeshParam%Edge(4,iElem)) + HydroParam%H(MeshParam%Edge(2,iElem)))
    !                Else
    !                    If (iLayer < HydroParam%Smallm(x2Elem)) Then
    !                        Facexin = Face
    !                    Else
    !                        Facexin = MeshParam%Edge(iEdge,x2Elem)
    !                    EndIf
    !                    HydroParam%uxy(iLayer,1,Face) = (HydroParam%uxy(iLayer,1,Face)*HydroParam%H(MeshParam%Edge(2,iElem)) + HydroParam%uxy(iLayer,1,Facexin)*HydroParam%H(MeshParam%Edge(4,iElem)))/(HydroParam%H(MeshParam%Edge(2,iElem)) + HydroParam%H(MeshParam%Edge(4,iElem)))
    !                EndIf
    !    
    !    
    !                If (HydroParam%uxy(iLayer,2,Face) + HydroParam%uxy(iLayer,2,Faceyin)>=0) Then
    !                    HydroParam%uxy(iLayer,2,Face) = (HydroParam%uxy(iLayer,2,Face)*HydroParam%H(Face) + HydroParam%uxy(iLayer,2,Faceyin)*HydroParam%H(Faceyin))/(HydroParam%H(Face) + HydroParam%H(Faceyin))            
    !                Else
    !                    If (iLayer < HydroParam%Smallm(y1Elem)) Then
    !                        Faceyin = Face
    !                    Else                            
    !                        Faceyin = MeshParam%Edge(iEdge,y1Elem)
    !                    EndIf
    !                    HydroParam%uxy(iLayer,2,Face) = (HydroParam%uxy(iLayer,2,Face)*HydroParam%H(Face) + HydroParam%uxy(iLayer,2,Faceyin)*HydroParam%H(Faceyin))/(HydroParam%H(Face) + HydroParam%H(Faceyin))            
    !                                
    !                EndIf       
    !            ElseIf (iEdge==2) Then
    !    
    !                Faceyin = MeshParam%Edge(iEdge,y2Elem)                       
    !                Facexin = MeshParam%Edge(iEdge,x1Elem)                    
    !                
    !                If (iLayer < HydroParam%Smallm(x1Elem)) Then
    !                    Facexin = Face
    !                Else
    !                    Facexin = MeshParam%Edge(iEdge,x1Elem)
    !                EndIf                    
    !                ! x-direction component:
    !                If (HydroParam%uxy(iLayer,1,Face) + HydroParam%uxy(iLayer,1,Facexin)>=0) Then                         
    !                    HydroParam%uxy(iLayer,1,Face) = (HydroParam%uxy(iLayer,1,Face)*HydroParam%H(Face) + HydroParam%uxy(iLayer,1,Facexin)*HydroParam%H(Facexin))/(HydroParam%H(Face) + HydroParam%H(Facexin))
    !                Else
    !                    HydroParam%uxy(iLayer,1,Face) = (HydroParam%uxy(iLayer,1,Face)*HydroParam%H(Face) + HydroParam%uxy(iLayer,1,4)*HydroParam%H(MeshParam%Edge(4,iElem)))/(HydroParam%H(Face) + HydroParam%H(MeshParam%Edge(4,iElem)))
    !                EndIf
    !                ! y-direction component:        
    !                If (iLayer < HydroParam%Smallm(y2Elem)) Then
    !                    Faceyin = Face
    !                Else
    !                    Faceyin = MeshParam%Edge(iEdge,y2Elem)
    !                EndIf                       
    !                
    !                If (HydroParam%uxy(iLayer,2,Face) + HydroParam%uxy(iLayer,2,Faceyin)>=0) Then
    !                    HydroParam%uxy(iLayer,2,Face) = (HydroParam%uxy(iLayer,2,Face)*HydroParam%H(MeshParam%Edge(1,iElem)) + HydroParam%uxy(iLayer,2,Faceyin)*HydroParam%H(MeshParam%Edge(3,iElem)))/(HydroParam%H(MeshParam%Edge(1,iElem)) + HydroParam%H(MeshParam%Edge(3,iElem)))            
    !                Else
    !                    If (iLayer < HydroParam%Smallm(y1Elem)) Then
    !                        Faceyin = Face
    !                    Else                            
    !                        Faceyin = MeshParam%Edge(iEdge,y1Elem)
    !                    EndIf                       
    !                    HydroParam%uxy(iLayer,2,Face) = (HydroParam%uxy(iLayer,2,Face)*HydroParam%H(MeshParam%Edge(3,iElem)) + HydroParam%uxy(iLayer,2,Faceyin)*HydroParam%H(MeshParam%Edge(1,iElem)))/(HydroParam%H(MeshParam%Edge(3,iElem)) + HydroParam%H(MeshParam%Edge(1,iElem)))                                  
    !                EndIf
    !            ElseIf (iEdge == 3) Then
    !                
    !                Faceyin = MeshParam%Edge(3,y2Elem)                       
    !                Facexin = MeshParam%Edge(iEdge,x1Elem)
    !                ! x-direction component:
    !                If (iLayer < HydroParam%Smallm(x1Elem)) Then
    !                    Facexin = Face
    !                Else
    !                    Facexin = MeshParam%Edge(iEdge,x1Elem)
    !                EndIf                    
    !   
    !                If (HydroParam%uxy(iLayer,1,Face) + HydroParam%uxy(iLayer,1,Facexin)>=0) Then                      
    !                    HydroParam%uxy(iLayer,1,Face) = (HydroParam%uxy(iLayer,1,Face)*HydroParam%H(MeshParam%Edge(4,iElem)) + HydroParam%uxy(iLayer,1,Facexin)*HydroParam%H(MeshParam%Edge(2,iElem)))/(HydroParam%H(MeshParam%Edge(4,iElem)) + HydroParam%H(MeshParam%Edge(2,iElem)))
    !                Else
    !                    If (iLayer < HydroParam%Smallm(x2Elem)) Then
    !                        Facexin = Face
    !                    Else
    !                        Facexin = MeshParam%Edge(iEdge,x2Elem)
    !                    EndIf
    !                    HydroParam%uxy(iLayer,1,Face) = (HydroParam%uxy(iLayer,1,Face)*HydroParam%H(MeshParam%Edge(2,iElem)) + HydroParam%uxy(iLayer,1,Facexin)*HydroParam%H(MeshParam%Edge(4,iElem)))/(HydroParam%H(MeshParam%Edge(4,iElem)) + HydroParam%H(MeshParam%Edge(2,iElem)))
    !                EndIf
    !                ! y-direction component:
    !                If (iLayer < HydroParam%Smallm(y2Elem)) Then
    !                    Faceyin = Face
    !                Else                            
    !                    Faceyin = MeshParam%Edge(iEdge,y2Elem)
    !                EndIf
    !
    !                If (HydroParam%uxy(iLayer,2,Face) + HydroParam%uxy(iLayer,2,Faceyin)>=0) Then
    !                    HydroParam%uxy(iLayer,2,Face) = (HydroParam%uxy(iLayer,2,Face)*HydroParam%H(Face) + HydroParam%uxy(iLayer,2,Faceyin)*HydroParam%H(Faceyin))/(HydroParam%H(Face) + HydroParam%H(Faceyin))            
    !                Else
    !                    Faceyin = MeshParam%Edge(1,iElem)
    !                    HydroParam%uxy(iLayer,2,Face) = (HydroParam%uxy(iLayer,2,Face)*HydroParam%H(Face) + HydroParam%uxy(iLayer,2,Faceyin)*HydroParam%H(Faceyin))/(HydroParam%H(Face) + HydroParam%H(Faceyin))                                             
    !                EndIf                
    !            Else
    !                !iEdge = 4                    
    !                Faceyin = MeshParam%Edge(iEdge,y2Elem)                       
    !                Facexin = MeshParam%Edge(2,iElem)                    
    !           
    !                ! x-direction component:
    !                If (HydroParam%uxy(iLayer,1,Face) + HydroParam%uxy(iLayer,1,Facexin)>=0) Then                         
    !                    HydroParam%uxy(iLayer,1,Face) = (HydroParam%uxy(iLayer,1,Face)*HydroParam%H(Face) + HydroParam%uxy(iLayer,1,Facexin)*HydroParam%H(Facexin))/(HydroParam%H(Face) + HydroParam%H(Facexin))
    !                Else
    !                    If (iLayer < HydroParam%Smallm(x2Elem)) Then
    !                        Facexin = Face
    !                    Else
    !                        Facexin = MeshParam%Edge(iEdge,x2Elem)
    !                    EndIf                                                 
    !                    HydroParam%uxy(iLayer,1,Face) = (HydroParam%uxy(iLayer,1,Face)*HydroParam%H(Face) + HydroParam%uxy(iLayer,1,Facexin)*HydroParam%H(Facexin))/(HydroParam%H(Face) + HydroParam%H(Facexin))
    !                EndIf
    !
    !                ! y-direction component:        
    !                If (iLayer < HydroParam%Smallm(y2Elem)) Then
    !                    Faceyin = Face
    !                Else
    !                    Faceyin = MeshParam%Edge(iEdge,y2Elem)
    !                EndIf                       
    !                
    !                If (HydroParam%uxy(iLayer,2,Face) + HydroParam%uxy(iLayer,2,Faceyin)>=0) Then
    !                    HydroParam%uxy(iLayer,2,Face) = (HydroParam%uxy(iLayer,2,Face)*HydroParam%H(MeshParam%Edge(1,iElem)) + HydroParam%uxy(iLayer,2,Faceyin)*HydroParam%H(MeshParam%Edge(3,iElem)))/(HydroParam%H(MeshParam%Edge(1,iElem)) + HydroParam%H(MeshParam%Edge(3,iElem)))            
    !                Else
    !                    If (iLayer < HydroParam%Smallm(y1Elem)) Then
    !                        Faceyin = Face
    !                    Else                            
    !                        Faceyin = MeshParam%Edge(iEdge,y1Elem)
    !                    EndIf                       
    !                    HydroParam%uxy(iLayer,2,Face) = (HydroParam%uxy(iLayer,2,Face)*HydroParam%H(MeshParam%Edge(3,iElem)) + HydroParam%uxy(iLayer,2,Faceyin)*HydroParam%H(MeshParam%Edge(1,iElem)))/(HydroParam%H(MeshParam%Edge(3,iElem)) + HydroParam%H(MeshParam%Edge(1,iElem)))                                  
    !                EndIf       
    !            EndIf           
    !        EndDo          
    !    EndDo
    !EndDo
    !
	!Do iElem = 1, MeshParam%nElem
 !       Do iEdge = 1, 4
 !           Face = MeshParam%Edge(iEdge,iElem)
 !           l = MeshParam%Left(Face) 
 !           r = MeshParam%Right(Face)    
 !           Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
 !               If (iEdge==1) then !North Edge
 !                   HydroParam%uxy(iLayer,1,Face) = -Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%utang(iLayer,MeshParam%Edge(iEdge,iElem)) 
 !                   HydroParam%uxy(iLayer,2,Face) =  Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem)) 
 !               ElseIf (iEdge==2) then !West Edge
 !                   HydroParam%uxy(iLayer,1,Face) = -Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem)) 
 !                   HydroParam%uxy(iLayer,2,Face) = -Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%utang(iLayer,MeshParam%Edge(iEdge,iElem)) 
 !               ElseIf (iEdge==3) then !South Edge
 !                   HydroParam%uxy(iLayer,1,Face) =  Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%utang(iLayer,MeshParam%Edge(iEdge,iElem)) 
 !                   HydroParam%uxy(iLayer,2,Face) = -Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem)) 
 !               Else !East Edge
 !                   HydroParam%uxy(iLayer,1,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%u(iLayer,MeshParam%Edge(iEdge,iElem)) 
 !                   HydroParam%uxy(iLayer,2,Face) = Sig(iElem,MeshParam%Right(MeshParam%Edge(iEdge,iElem)),MeshParam%Left(MeshParam%Edge(iEdge,iElem)))*HydroParam%utang(iLayer,MeshParam%Edge(iEdge,iElem)) 
 !               EndIf
 !           EndDo
 !       EndDo
 !   EndDo
            

    !Do iElem = 1, MeshParam%nElem 
    !    icountV=MeshParam%Area(iElem)
    !    n1=MeshParam%Quadri(1,iElem) + 1
    !    n2=MeshParam%Quadri(2,iElem) + 1
    !    n3=MeshParam%Quadri(3,iElem) + 1
    !    n4=MeshParam%Quadri(4,iElem) + 1
    !    aa(1)=signa(MeshParam%xNode(n1),MeshParam%xNode(n2),MeshParam%xb(iElem),MeshParam%yNode(n1),MeshParam%yNode(n2),MeshParam%yb(iElem))
    !    aa(2)=signa(MeshParam%xNode(n2),MeshParam%xNode(n3),MeshParam%xb(iElem),MeshParam%yNode(n2),MeshParam%yNode(n3),MeshParam%yb(iElem))
    !    aa(3)=signa(MeshParam%xNode(n3),MeshParam%xNode(n4),MeshParam%xb(iElem),MeshParam%yNode(n3),MeshParam%yNode(n4),MeshParam%yb(iElem))
    !    aa(4)=signa(MeshParam%xNode(n4),MeshParam%xNode(n1),MeshParam%xb(iElem),MeshParam%yNode(n4),MeshParam%yNode(n1),MeshParam%yb(iElem))
    !    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
    !        If (iElem==293.and.iLayer==1) Then
    !            Continue
    !        EndIf
    !        aa(5)=0.d0
    !        HydroParam%ub(iLayer,1:2,iElem) = 0.d0
    !        icountU=0.d0
    !        icountV=0.d0
    !        Do iEdge=1,4
    !            Face = MeshParam%Edge(iEdge,iElem) 
    !            l = MeshParam%Left(MeshParam%Edge(iEdge,iElem))
    !            r = MeshParam%Right(MeshParam%Edge(iEdge,iElem))
    !            
    !            If (HydroParam%IndexWaterLevelEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) Then !HydroParam%IndexWaterLevelEdge(Face)>0.or.
    !                aa(5)=aa(5)+aa(iEdge) !Total area contribution.
    !            Else
    !                if (r==0) then
    !                     cycle
    !                else
    !                    If (iLayer>=HydroParam%ElSmallm(r)) Then
    !                        aa(5)=aa(5)+aa(iEdge) !Total area contribution.
    !                    else
    !                        cycle
    !                    endif
    !                
    !                endif
    !            Endif
    !        EndDo
    !        if (aa(5)==0) then
    !            HydroParam%ub(iLayer,1,iElem) = 0
    !            HydroParam%ub(iLayer,2,iElem) = 0
    !            cycle
    !        endif
    !        
    !            Do iEdge=1,4
    !                Face = MeshParam%Edge(iEdge,iElem) 
    !                l = MeshParam%Left(MeshParam%Edge(iEdge,iElem))
    !                r = MeshParam%Right(MeshParam%Edge(iEdge,iElem))
    !                
    !                If (HydroParam%IndexWaterLevelEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) Then !HydroParam%IndexWaterLevelEdge(Face)>0.or.
    !                    HydroParam%ub(iLayer,1,iElem) = HydroParam%ub(iLayer,1,iElem) + HydroParam%uxy(iLayer,1,MeshParam%Edge(iEdge,iElem))*aa(iEdge)/aa(5)
    !                    HydroParam%ub(iLayer,2,iElem) = HydroParam%ub(iLayer,2,iElem) + HydroParam%uxy(iLayer,2,MeshParam%Edge(iEdge,iElem))*aa(iEdge)/aa(5)
    !                Else
    !                   if (r==0) then
    !                        cycle
    !                   Else
    !                        If (iLayer>=HydroParam%ElSmallm(r)) then
    !                            HydroParam%ub(iLayer,1,iElem) = HydroParam%ub(iLayer,1,iElem) + HydroParam%uxy(iLayer,1,MeshParam%Edge(iEdge,iElem))*aa(iEdge)/aa(5)
    !                            HydroParam%ub(iLayer,2,iElem) = HydroParam%ub(iLayer,2,iElem) + HydroParam%uxy(iLayer,2,MeshParam%Edge(iEdge,iElem))*aa(iEdge)/aa(5)
    !                        else
    !                            cycle
    !                        endif
    !                    EndIf
    !                    
    !                EndIf
    !                
    !            EndDo
    !            HydroParam%ub(iLayer,3,iElem) = 0.5*( HydroParam%w(iLayer,iElem) + HydroParam%w(iLayer+1,iElem) )
    !    EndDo
    !EndDo 
    
    ! 1. Finding cell-centered velocity components 

    Do iElem = 1, MeshParam%nElem 
        icountV=MeshParam%Area(iElem)
        n1=MeshParam%Quadri(1,iElem) + 1
        n2=MeshParam%Quadri(2,iElem) + 1
        n3=MeshParam%Quadri(3,iElem) + 1
        n4=MeshParam%Quadri(4,iElem) + 1
        aa(1)=signa(MeshParam%xNode(n1),MeshParam%xNode(n2),MeshParam%xb(iElem),MeshParam%yNode(n1),MeshParam%yNode(n2),MeshParam%yb(iElem))
        aa(2)=signa(MeshParam%xNode(n2),MeshParam%xNode(n3),MeshParam%xb(iElem),MeshParam%yNode(n2),MeshParam%yNode(n3),MeshParam%yb(iElem))
        aa(3)=signa(MeshParam%xNode(n3),MeshParam%xNode(n4),MeshParam%xb(iElem),MeshParam%yNode(n3),MeshParam%yNode(n4),MeshParam%yb(iElem))
        aa(4)=signa(MeshParam%xNode(n4),MeshParam%xNode(n1),MeshParam%xb(iElem),MeshParam%yNode(n4),MeshParam%yNode(n1),MeshParam%yb(iElem))

        Do iLayer = HydroParam%ElSmallms(iElem), HydroParam%ElCapitalM(iElem) 

            aa(5)=0.d0
            aaU=0.d0
            aaV=0.d0
            HydroParam%ub(iLayer,1:2,iElem) = 0.d0
            icountU=0.d0
            icountV=0.d0
            
            Do iEdge=1,4
                Face = MeshParam%Edge(iEdge,iElem) 
                l = MeshParam%Left(Face)
                r = MeshParam%Right(Face)
               
                If (abs(HydroParam%uxy(iLayer,1,Face))>NearZero) then !if iEdge has flux in X direction
                    aaU=aaU+aa(iEdge) !Total area contribution.
                EndIf
                If (abs(HydroParam%uxy(iLayer,2,Face))>NearZero) then !if iEdge has flux in Y direction
                    aaV=aaV+aa(iEdge) !Total area contribution.
                EndIf
            EndDo
            
            If (aaU==0) Then
                HydroParam%ub(iLayer,1,iElem) = 0
            EndIf
            If (aaV==0) Then
                HydroParam%ub(iLayer,2,iElem) = 0
            EndIf
            
            Do iEdge=1,4
                Face = MeshParam%Edge(iEdge,iElem) 
                l = MeshParam%Left(MeshParam%Edge(iEdge,iElem))
                r = MeshParam%Right(MeshParam%Edge(iEdge,iElem))

                If (abs(HydroParam%uxy(iLayer,1,Face))>NearZero.and.aaU>NearZero) then !if iEdge has flux in Y direction
                    HydroParam%ub(iLayer,1,iElem) = HydroParam%ub(iLayer,1,iElem) + HydroParam%uxy(iLayer,1,Face)*aa(iEdge)/aaU
                EndIf
                If (abs(HydroParam%uxy(iLayer,2,Face))>NearZero.and.aaV>NearZero) then !if iEdge has flux in X direction
                    HydroParam%ub(iLayer,2,iElem) = HydroParam%ub(iLayer,2,iElem) + HydroParam%uxy(iLayer,2,Face)*aa(iEdge)/aaV
                EndIf
                
            EndDo
            HydroParam%ub(iLayer,3,iElem) = 0.5*( HydroParam%w(iLayer,iElem) + HydroParam%w(iLayer+1,iElem) )

            ! If subsurface flow occurs, ubsub calculate as output:
            If (MeshParam%iBedrock == 1) Then
                If (aaU==0) Then
                    HydroParam%ubsub(iLayer,1,iElem) = 0
                EndIf
                If (aaV==0) Then
                    HydroParam%ubsub(iLayer,2,iElem) = 0
                EndIf            
            
                Do iEdge=1,4                  

                    If (abs(HydroParam%uxysub(iLayer,1,Face))>NearZero.and.aaU>NearZero) then !if iEdge has flux in Y direction
                        HydroParam%ubsub(iLayer,1,iElem) = HydroParam%ubsub(iLayer,1,iElem) + HydroParam%uxysub(iLayer,1,Face)*aa(iEdge)/aaU
                    EndIf
                    If (abs(HydroParam%uxysub(iLayer,2,Face))>NearZero.and.aaV>NearZero) then !if iEdge has flux in X direction
                        HydroParam%ubsub(iLayer,2,iElem) = HydroParam%ubsub(iLayer,2,iElem) + HydroParam%uxysub(iLayer,2,Face)*aa(iEdge)/aaV     
                    EndIf
                EndDo
                HydroParam%ubsub(iLayer,3,iElem) = 0.5*( HydroParam%wm(iLayer,iElem) + HydroParam%wm(iLayer+1,iElem) )
            EndIf  
                   
        EndDo
    EndDo   
    
    ! 2. Finding Face-centered Vertical velocity components
    Do iEdge=1,MeshParam%nEdge
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
        If (r==0) Then
            HydroParam%wfc(:,iEdge)= HydroParam%ub(:,3,l)
        Else        
            Do iLayer = 1, HydroParam%CapitalM(iEdge) !HydroParam%Smallm(iEdge)
                If (abs(HydroParam%ub(iLayer,3,l))>NearZero.and.abs(HydroParam%ub(iLayer,3,r))<NearZero) Then !Velocity in lateral boundries is equal to cell center velocity
                    HydroParam%wfc(iLayer,iEdge)= HydroParam%ub(iLayer,3,l)
                Elseif (abs(HydroParam%ub(iLayer,3,l))<NearZero.and.abs(HydroParam%ub(iLayer,3,r))>NearZero) Then
                    HydroParam%wfc(iLayer,iEdge)= HydroParam%ub(iLayer,3,r)
                ElseIf (abs(HydroParam%ub(iLayer,3,l))>NearZero.and.abs(HydroParam%ub(iLayer,3,r))>NearZero) Then
                    HydroParam%wfc(iLayer,iEdge) = (0.5d0*(HydroParam%ub(iLayer,3,l)+HydroParam%ub(iLayer,3,r)))
                Else
                    HydroParam%wfc(iLayer,iEdge)= 0.d0
                EndIf
                If (isnan(HydroParam%wfc(iLayer,iEdge))) then
                    continue
                EndIf
            EndDo
        EndIf
        !HydroParam%wfc(iLayer,iEdge) = 0.5d0 * (0.5d0*(HydroParam%ub(iLayer,3,l)+HydroParam%ub(iLayer,3,r)) + (0.5d0*(HydroParam%wg(iEdge,iLayer)+HydroParam%wg(iEdge,iLayer+1))))
    EndDo
    
    ! 3. Finding k+1/2 Layer-centred horizontal Velocity components
    Do iElem = 1, MeshParam%nElem
        HydroParam%uxyL(HydroParam%ElSmallms(iElem),1,iElem) = HydroParam%ub(HydroParam%ElSmallms(iElem),1,iElem)
        HydroParam%uxyL(HydroParam%ElSmallms(iElem),2,iElem) = HydroParam%ub(HydroParam%ElSmallms(iElem),2,iElem)
        !HydroParam%uxyL(HydroParam%ElCapitalM(iElem)+1,1,iElem) = HydroParam%ub(HydroParam%ElCapitalM(iElem),1,iElem)
        !HydroParam%uxyL(HydroParam%ElCapitalM(iElem)+1,2,iElem) = HydroParam%ub(HydroParam%ElCapitalM(iElem),2,iElem)

        Do iLayer = HydroParam%ElSmallms(iElem)+1, HydroParam%ElCapitalM(iElem)           
            HydroParam%uxyL(iLayer,1,iElem) = HydroParam%ub(iLayer-1,1,iElem)+((HydroParam%Zb(iLayer,iElem)-HydroParam%Ze(iLayer,iElem))*((HydroParam%ub(iLayer,1,iElem)-HydroParam%ub(iLayer-1,1,iElem))/(HydroParam%Zb(iLayer,iElem)-HydroParam%Zb(iLayer-1,iElem))))
            HydroParam%uxyL(iLayer,2,iElem) = HydroParam%ub(iLayer-1,2,iElem)+((HydroParam%Zb(iLayer,iElem)-HydroParam%Ze(iLayer,iElem))*((HydroParam%ub(iLayer,2,iElem)-HydroParam%ub(iLayer-1,2,iElem))/(HydroParam%Zb(iLayer,iElem)-HydroParam%Zb(iLayer-1,iElem))))
        EndDo
        If (HydroParam%ElCapitalM(iElem)==HydroParam%ElSmallms(iElem)) then ! Only 1 Layer
            HydroParam%uxyL(HydroParam%ElCapitalM(iElem)+1,1,iElem) = HydroParam%ub(HydroParam%ElCapitalM(iElem),1,iElem)!+abs((HydroParam%eta(iElem))-(HydroParam%Zb(HydroParam%ElCapitalM(iElem),iElem)))*((HydroParam%ub(HydroParam%ElCapitalM(iElem),1,iElem)-HydroParam%ub(HydroParam%ElCapitalM(iElem)-1,1,iElem))/(abs(abs(HydroParam%Zb(HydroParam%ElCapitalM(iElem)-1,iElem))-abs(HydroParam%Zb(HydroParam%ElCapitalM(iElem),iElem)))))
            HydroParam%uxyL(HydroParam%ElCapitalM(iElem)+1,2,iElem) = HydroParam%ub(HydroParam%ElCapitalM(iElem),2,iElem)!+abs((HydroParam%eta(iElem))-(HydroParam%Zb(HydroParam%ElCapitalM(iElem),iElem)))*((HydroParam%ub(HydroParam%ElCapitalM(iElem),2,iElem)-HydroParam%ub(HydroParam%ElCapitalM(iElem)-1,2,iElem))/(abs(abs(HydroParam%Zb(HydroParam%ElCapitalM(iElem)-1,iElem))-abs(HydroParam%Zb(HydroParam%ElCapitalM(iElem),iElem)))))
        Else !If the domain has more than one layer, use lagrange polynomial interpolation
            HydroParam%uxyL(HydroParam%ElCapitalM(iElem)+1,1,iElem) = HydroParam%ub(HydroParam%ElCapitalM(iElem),1,iElem)
            HydroParam%uxyL(HydroParam%ElCapitalM(iElem)+1,2,iElem) = HydroParam%ub(HydroParam%ElCapitalM(iElem),2,iElem)
            !HydroParam%uxyL(HydroParam%ElCapitalM(iElem)+1,1,iElem) = HydroParam%ub(HydroParam%ElCapitalM(iElem),1,iElem)+abs((HydroParam%eta(iElem))-(HydroParam%Zb(HydroParam%ElCapitalM(iElem),iElem)))*((HydroParam%ub(HydroParam%ElCapitalM(iElem),1,iElem)-HydroParam%ub(HydroParam%ElCapitalM(iElem)-1,1,iElem))/(abs(abs(HydroParam%Zb(HydroParam%ElCapitalM(iElem)-1,iElem))-abs(HydroParam%Zb(HydroParam%ElCapitalM(iElem),iElem)))))
            !HydroParam%uxyL(HydroParam%ElCapitalM(iElem)+1,2,iElem) = HydroParam%ub(HydroParam%ElCapitalM(iElem),2,iElem)+abs((HydroParam%eta(iElem))-(HydroParam%Zb(HydroParam%ElCapitalM(iElem),iElem)))*((HydroParam%ub(HydroParam%ElCapitalM(iElem),2,iElem)-HydroParam%ub(HydroParam%ElCapitalM(iElem)-1,2,iElem))/(abs(abs(HydroParam%Zb(HydroParam%ElCapitalM(iElem)-1,iElem))-abs(HydroParam%Zb(HydroParam%ElCapitalM(iElem),iElem)))))
        EndIf
    EndDo    
    
    ! 4. Finding Velocities in each face in k+1/2 Layer
    Do iEdge=1,MeshParam%nEdge
        l = MeshParam%Left(iEdge)
        r = MeshParam%Right(iEdge)
        !Bottom and top horizontal velocitys are iqual to cell center velocities - layer "k" is defined in the center of each cell
        !Bottom vertical Velocity is = 0, top layer vertical velocity is = w top layer velocity
        !if (r==0) then
        !    small = HydroParam%Smallm(iEdge)
        !else
        !    small = min(HydroParam%ElSmallm(l),HydroParam%ElSmallm(r))
        !endif
        Do iLayer = HydroParam%Smallms(iEdge)+1, HydroParam%CapitalM(iEdge) !HydroParam%Smallm(iEdge)+1
        
        !if (r/=0) then
        !    if (HydroParam%ElSmallm(r)<HydroParam%ElSmallm(l)) then
        !        if (ilayer<HydroParam%ElSmallm(l)) then
        !            l=r
        !        else
        !            l = MeshParam%Left(iEdge)
        !            r = MeshParam%Right(iEdge)
        !        endif
        !    endif
        !endif
    
            !Find Velocities to K+1/2 for each edge
            If (r==0) Then !Velocity in lateral boundaries is equal to cell center velocity
                If (HydroParam%uxy(iLayer,1,iEdge)==0.and.HydroParam%uxy(iLayer,2,iEdge)/=0) then !no flux in U direction
                     HydroParam%ug(iEdge,iLayer)=0
                     HydroParam%vg(iEdge,iLayer)= HydroParam%uxy(iLayer-1,2,iEdge)+(((HydroParam%Zb(iLayer,l)-HydroParam%Ze(iLayer,l))/(HydroParam%Zb(iLayer,l)-HydroParam%Zb(iLayer-1,l)))*(HydroParam%uxy(iLayer,2,iEdge)-HydroParam%uxy(iLayer-1,2,iEdge)))
                     HydroParam%wg(iEdge,iLayer)= HydroParam%w(iLayer,l) 
                Elseif (HydroParam%uxy(iLayer,2,iEdge)==0.and.HydroParam%uxy(iLayer,1,iEdge)/=0) then !no flux in V direction
                     HydroParam%ug(iEdge,iLayer)= HydroParam%uxy(iLayer-1,1,iEdge)+(((HydroParam%Zb(iLayer,l)-HydroParam%Ze(iLayer,l))/(HydroParam%Zb(iLayer,l)-HydroParam%Zb(iLayer-1,l)))*(HydroParam%uxy(iLayer,1,iEdge)-HydroParam%uxy(iLayer-1,1,iEdge)))
                     HydroParam%vg(iEdge,iLayer)=0
                     HydroParam%wg(iEdge,iLayer)= HydroParam%w(iLayer,l)
                Elseif (HydroParam%uxy(iLayer,2,iEdge)==0.and.HydroParam%uxy(iLayer,1,iEdge)==0) then !no flux in U or V direction
                     HydroParam%ug(iEdge,iLayer)=0
                     HydroParam%vg(iEdge,iLayer)=0
                     HydroParam%wg(iEdge,iLayer)= HydroParam%w(iLayer,l) !(HydroParam%wfc(iLayer-1,iEdge)+(HydroParam%DZj(iLayer-1,iEdge)/(HydroParam%DZj(iLayer-1,iEdge)+HydroParam%DZj(iLayer,iEdge)))*(HydroParam%wfc(iLayer,iEdge)-HydroParam%wfc(iLayer-1,iEdge)))
                Else !Flux in U and V direction
                    HydroParam%ug(iEdge,iLayer)= HydroParam%uxy(iLayer-1,1,iEdge)+(((HydroParam%Zb(iLayer,l)-HydroParam%Ze(iLayer,l))/(HydroParam%Zb(iLayer,l)-HydroParam%Zb(iLayer-1,l)))*(HydroParam%uxy(iLayer,1,iEdge)-HydroParam%uxy(iLayer-1,1,iEdge)))
                    HydroParam%vg(iEdge,iLayer)= HydroParam%uxy(iLayer-1,2,iEdge)+(((HydroParam%Zb(iLayer,l)-HydroParam%Ze(iLayer,l))/(HydroParam%Zb(iLayer,l)-HydroParam%Zb(iLayer-1,l)))*(HydroParam%uxy(iLayer,2,iEdge)-HydroParam%uxy(iLayer-1,2,iEdge)))
                !HydroParam%ug(iEdge,iLayer)=HydroParam%uxyL(iLayer,1,l)
                !HydroParam%vg(iEdge,iLayer)= HydroParam%uxyL(iLayer,2,l)
                    HydroParam%wg(iEdge,iLayer)= HydroParam%w(iLayer,l) !(HydroParam%wfc(iLayer-1,iEdge)+(HydroParam%DZj(iLayer-1,iEdge)/(HydroParam%DZj(iLayer-1,iEdge)+HydroParam%DZj(iLayer,iEdge)))*(HydroParam%wfc(iLayer,iEdge)-HydroParam%wfc(iLayer-1,iEdge)))
                Endif
            Else
                If (iLayer>=HydroParam%ElSmallms(r)) Then 
                    !HydroParam%ug(iEdge,iLayer)= 0.5d0*((0.5d0*(HydroParam%uxyL(iLayer,1,r)+HydroParam%uxyL(iLayer,1,l))+HydroParam%uxy(iLayer-1,1,iEdge)+(HydroParam%DZj(iLayer-1,iEdge)/(HydroParam%DZj(iLayer-1,iEdge)+HydroParam%DZj(iLayer,iEdge)))*(HydroParam%uxy(iLayer,1,iEdge)-HydroParam%uxy(iLayer-1,1,iEdge))))
                    !HydroParam%vg(iEdge,iLayer)= 0.5d0*((0.5d0*(HydroParam%uxyL(iLayer,2,r)+HydroParam%uxyL(iLayer,2,l))+HydroParam%uxy(iLayer-1,2,iEdge)+(HydroParam%DZj(iLayer-1,iEdge)/(HydroParam%DZj(iLayer-1,iEdge)+HydroParam%DZj(iLayer,iEdge)))*(HydroParam%uxy(iLayer,2,iEdge)-HydroParam%uxy(iLayer-1,2,iEdge))))
                    !HydroParam%wg(iEdge,iLayer) = 0.5d0*((0.5d0*(HydroParam%w(iLayer,r)+HydroParam%w(iLayer,l)))+(HydroParam%wfc(iLayer-1,iEdge)+(HydroParam%DZj(iLayer-1,iEdge)/(HydroParam%DZj(iLayer-1,iEdge)+HydroParam%DZj(iLayer,iEdge)))*(HydroParam%wfc(iLayer,iEdge)-HydroParam%wfc(iLayer-1,iEdge))))
                    HydroParam%ug(iEdge,iLayer)= HydroParam%uxy(iLayer-1,1,iEdge)+(((HydroParam%Zb(iLayer,l)-HydroParam%Ze(iLayer,l))/(HydroParam%Zb(iLayer,l)-HydroParam%Zb(iLayer-1,l)))*(HydroParam%uxy(iLayer,1,iEdge)-HydroParam%uxy(iLayer-1,1,iEdge)))
                    HydroParam%vg(iEdge,iLayer)= HydroParam%uxy(iLayer-1,2,iEdge)+(((HydroParam%Zb(iLayer,l)-HydroParam%Ze(iLayer,l))/(HydroParam%Zb(iLayer,l)-HydroParam%Zb(iLayer-1,l)))*(HydroParam%uxy(iLayer,2,iEdge)-HydroParam%uxy(iLayer-1,2,iEdge)))
                    HydroParam%wg(iEdge,iLayer) = 0.5d0*(HydroParam%w(iLayer,r)+HydroParam%w(iLayer,l))!+(HydroParam%wfc(iLayer-1,iEdge)+(HydroParam%DZj(iLayer-1,iEdge)/(HydroParam%DZj(iLayer-1,iEdge)+HydroParam%DZj(iLayer,iEdge)))*(HydroParam%wfc(iLayer,iEdge)-HydroParam%wfc(iLayer-1,iEdge))))
                Else
                    If (HydroParam%uxy(iLayer,1,iEdge)==0.and.HydroParam%uxy(iLayer,2,iEdge)/=0) then
                         HydroParam%ug(iEdge,iLayer)=0
                         HydroParam%vg(iEdge,iLayer)= HydroParam%uxy(iLayer-1,2,iEdge)+(((HydroParam%Zb(iLayer,l)-HydroParam%Ze(iLayer,l))/(HydroParam%Zb(iLayer,l)-HydroParam%Zb(iLayer-1,l)))*(HydroParam%uxy(iLayer,2,iEdge)-HydroParam%uxy(iLayer-1,2,iEdge)))
                         HydroParam%wg(iEdge,iLayer)= HydroParam%w(iLayer,l) 
                    Elseif (HydroParam%uxy(iLayer,2,iEdge)==0.and.HydroParam%uxy(iLayer,1,iEdge)/=0) then
                         HydroParam%ug(iEdge,iLayer)= HydroParam%uxy(iLayer-1,1,iEdge)+(((HydroParam%Zb(iLayer,l)-HydroParam%Ze(iLayer,l))/(HydroParam%Zb(iLayer,l)-HydroParam%Zb(iLayer-1,l)))*(HydroParam%uxy(iLayer,1,iEdge)-HydroParam%uxy(iLayer-1,1,iEdge)))
                         HydroParam%vg(iEdge,iLayer)=0
                         HydroParam%wg(iEdge,iLayer)= HydroParam%w(iLayer,l)
                    Elseif (HydroParam%uxy(iLayer,2,iEdge)==0.and.HydroParam%uxy(iLayer,1,iEdge)==0) then
                         HydroParam%ug(iEdge,iLayer)=0
                         HydroParam%vg(iEdge,iLayer)=0
                         HydroParam%wg(iEdge,iLayer)= HydroParam%w(iLayer,l) !(HydroParam%wfc(iLayer-1,iEdge)+(HydroParam%DZj(iLayer-1,iEdge)/(HydroParam%DZj(iLayer-1,iEdge)+HydroParam%DZj(iLayer,iEdge)))*(HydroParam%wfc(iLayer,iEdge)-HydroParam%wfc(iLayer-1,iEdge)))
                    Else
                        HydroParam%ug(iEdge,iLayer)= HydroParam%uxy(iLayer-1,1,iEdge)+(((HydroParam%Zb(iLayer,l)-HydroParam%Ze(iLayer,l))/(HydroParam%Zb(iLayer,l)-HydroParam%Zb(iLayer-1,l)))*(HydroParam%uxy(iLayer,1,iEdge)-HydroParam%uxy(iLayer-1,1,iEdge)))
                        HydroParam%vg(iEdge,iLayer)= HydroParam%uxy(iLayer-1,2,iEdge)+(((HydroParam%Zb(iLayer,l)-HydroParam%Ze(iLayer,l))/(HydroParam%Zb(iLayer,l)-HydroParam%Zb(iLayer-1,l)))*(HydroParam%uxy(iLayer,2,iEdge)-HydroParam%uxy(iLayer-1,2,iEdge)))
                        !HydroParam%ug(iEdge,iLayer)=HydroParam%uxyL(iLayer,1,l)
                        !HydroParam%vg(iEdge,iLayer)=HydroParam%uxyL(iLayer,2,l)
                        HydroParam%wg(iEdge,iLayer)= HydroParam%w(iLayer,l) !(HydroParam%wfc(iLayer-1,iEdge)+(HydroParam%DZj(iLayer-1,iEdge)/(HydroParam%DZj(iLayer-1,iEdge)+HydroParam%DZj(iLayer,iEdge)))*(HydroParam%wfc(iLayer,iEdge)-HydroParam%wfc(iLayer-1,iEdge)))
                    EndIf
                EndIf
            EndIf
        EndDo
    ! 4.1. Finding ug,vg and wg to bottom and top layers, applying the boundary conditions
        If (HydroParam%ElCapitalM(l)==HydroParam%ElSmallms(l)) then
            HydroParam%ug(iEdge,HydroParam%CapitalM(iEdge)+1)= HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge)!+weitW*((HydroParam%uxy(HydroParam%ElCapitalM(l),1,iEdge)-HydroParam%uxy(HydroParam%ElCapitalM(l)-1,1,iEdge))/(weit))
            HydroParam%vg(iEdge,HydroParam%CapitalM(iEdge)+1)= HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)!+weitW*((HydroParam%uxy(HydroParam%ElCapitalM(l),2,iEdge)-HydroParam%uxy(HydroParam%ElCapitalM(l)-1,2,iEdge))/(weit))
        Else !If the domain has more than one layer, use lagrange polynomial interpolation
            HydroParam%ug(iEdge,HydroParam%CapitalM(iEdge)+1)= HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge)
            HydroParam%vg(iEdge,HydroParam%CapitalM(iEdge)+1)= HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)
        EndIf
 
        If (r==0) Then !Velocity in lateral boundaries is equal to cell center velocity
            HydroParam%wg(iEdge,HydroParam%Smallms(iEdge))=0
            HydroParam%wg(iEdge,HydroParam%CapitalM(iEdge)+1)=HydroParam%w(HydroParam%CapitalM(iEdge)+1,l)
                
            If (HydroParam%uxy(HydroParam%Smallms(iEdge),1,iEdge)==0.and.HydroParam%uxy(HydroParam%Smallms(iEdge),2,iEdge)/=0) then
                HydroParam%ug(iEdge,HydroParam%Smallms(iEdge))= 0.d0
                HydroParam%vg(iEdge,HydroParam%Smallms(iEdge))=HydroParam%uxy(HydroParam%Smallms(iElem),2,iEdge)
            Elseif (HydroParam%uxy(HydroParam%Smallms(iEdge),2,iEdge)==0.and.HydroParam%uxy(HydroParam%Smallms(iEdge),1,iEdge)/=0) then
                HydroParam%ug(iEdge,HydroParam%Smallms(iEdge))= HydroParam%uxy(HydroParam%Smallms(iEdge),1,iEdge) 
                HydroParam%vg(iEdge,HydroParam%Smallms(iEdge))= 0
            Elseif (HydroParam%uxy(HydroParam%Smallms(iEdge),2,iEdge)==0.and.HydroParam%uxy(HydroParam%Smallms(iEdge),1,iEdge)==0) then
                HydroParam%ug(iEdge,HydroParam%Smallms(iEdge))= 0
                HydroParam%vg(iEdge,HydroParam%Smallms(iEdge))= 0
            Else
                HydroParam%ug(iEdge,HydroParam%Smallms(iEdge))= HydroParam%uxy(HydroParam%Smallms(iEdge),1,iEdge) 
                HydroParam%vg(iEdge,HydroParam%Smallms(iEdge))= HydroParam%uxy(HydroParam%Smallms(iEdge),2,iEdge)
            EndIf
            
            If (HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge)==0.and.HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)/=0) then
                HydroParam%ug(iEdge,HydroParam%CapitalM(iEdge)+1)= 0.d0
                !HydroParam%vg(iEdge,HydroParam%CapitalM(iEdge)+1)=HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)
            Elseif (HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)==0.and.HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge)/=0) then
                !HydroParam%ug(iEdge,HydroParam%CapitalM(iEdge)+1)= HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge)
                HydroParam%vg(iEdge,HydroParam%CapitalM(iEdge)+1)= 0
            Elseif (HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)==0.and.HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge)==0) then
                HydroParam%ug(iEdge,HydroParam%CapitalM(iEdge)+1)= 0
                HydroParam%vg(iEdge,HydroParam%CapitalM(iEdge)+1)= 0
            Else
                !HydroParam%ug(iEdge,HydroParam%CapitalM(iEdge)+1)= HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge) 
                !HydroParam%vg(iEdge,HydroParam%CapitalM(iEdge)+1)= HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)
            EndIf
                        
        Else
            If (iLayer>=HydroParam%ElSmallm(r)) Then
                HydroParam%ug(iEdge,HydroParam%Smallms(iEdge))= HydroParam%uxy(HydroParam%Smallms(iEdge),1,iEdge)
                HydroParam%vg(iEdge,HydroParam%Smallms(iEdge))= HydroParam%uxy(HydroParam%Smallms(iEdge),2,iEdge)
                
                HydroParam%wg(iEdge,HydroParam%Smallms(iEdge))= 0.5d0*(HydroParam%w(HydroParam%Smallms(iEdge),l)+HydroParam%w(HydroParam%Smallms(iEdge),r))
                HydroParam%wg(iEdge,HydroParam%CapitalM(iEdge)+1) = 0.5d0*(HydroParam%w(HydroParam%CapitalM(iEdge)+1,l)+HydroParam%w(HydroParam%CapitalM(iEdge)+1,r))
            Else
                HydroParam%wg(iEdge,HydroParam%Smallms(iEdge))=0
                HydroParam%wg(iEdge,HydroParam%CapitalM(iEdge)+1)=HydroParam%w(HydroParam%CapitalM(iEdge)+1,l)
                
                If (HydroParam%uxy(HydroParam%Smallms(iEdge),1,iEdge)==0.and.HydroParam%uxy(HydroParam%Smallms(iEdge),2,iEdge)/=0) then
                    HydroParam%ug(iEdge,HydroParam%Smallms(iEdge))= 0.d0
                    HydroParam%vg(iEdge,HydroParam%Smallms(iEdge))=HydroParam%uxy(HydroParam%Smallms(iElem),2,iEdge)
                Elseif (HydroParam%uxy(HydroParam%Smallms(iEdge),2,iEdge)==0.and.HydroParam%uxy(HydroParam%Smallms(iEdge),1,iEdge)/=0) then
                    HydroParam%ug(iEdge,HydroParam%Smallms(iEdge))= HydroParam%uxy(HydroParam%Smallms(iEdge),1,iEdge) 
                    HydroParam%vg(iEdge,HydroParam%Smallms(iEdge))= 0
                Elseif (HydroParam%uxy(HydroParam%Smallms(iEdge),2,iEdge)==0.and.HydroParam%uxy(HydroParam%Smallms(iEdge),1,iEdge)==0) then
                    HydroParam%ug(iEdge,HydroParam%Smallms(iEdge))= 0
                    HydroParam%vg(iEdge,HydroParam%Smallms(iEdge))= 0
                Else
                    HydroParam%ug(iEdge,HydroParam%Smallms(iEdge))= HydroParam%uxy(HydroParam%Smallms(iEdge),1,iEdge) 
                    HydroParam%vg(iEdge,HydroParam%Smallms(iEdge))= HydroParam%uxy(HydroParam%Smallms(iEdge),2,iEdge)
                EndIf
            
                If (HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge)==0.and.HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)/=0) then
                    HydroParam%ug(iEdge,HydroParam%CapitalM(iEdge)+1)= 0.d0
                    !HydroParam%vg(iEdge,HydroParam%CapitalM(iEdge)+1)=HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)
                Elseif (HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)==0.and.HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge)/=0) then
                    !HydroParam%ug(iEdge,HydroParam%CapitalM(iEdge)+1)= HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge)
                    HydroParam%vg(iEdge,HydroParam%CapitalM(iEdge)+1)= 0
                Elseif (HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)==0.and.HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge)==0) then
                    HydroParam%ug(iEdge,HydroParam%CapitalM(iEdge)+1)= 0
                    HydroParam%vg(iEdge,HydroParam%CapitalM(iEdge)+1)= 0
                Else
                    !HydroParam%ug(iEdge,HydroParam%CapitalM(iEdge)+1)= HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge) 
                    !HydroParam%vg(iEdge,HydroParam%CapitalM(iEdge)+1)= HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)
                EndIf
            EndIf
            
        EndIf
    
    EndDo !i=1,ns
    
  !  Do iEdge=1,MeshParam%nEdge
  !      l = MeshParam%Left(iEdge)
  !      r = MeshParam%Right(iEdge)
  !      !Bottom and top horizontal velocitys are iqual to cell center velocities - layer "k" is difined in the center of each cell
  !      !Bottom vertical Velocity is = 0, top layer vertical velocity is = w top layer velocity
  !      if (r==0) then
  !          small = HydroParam%Smallms(iEdge)
  !      else
  !          small = min(HydroParam%ElSmallm(l),HydroParam%ElSmallm(r))
  !      endif
  !      
  !      !4.1. Finding ug,vg and wg to bottom and top layers, apllying the boundry conditions
  !      if (HydroParam%ElCapitalM(l)==1.and.HydroParam%ElSmallm(l)==1) then
		!	HydroParam%ug(iEdge,HydroParam%ElCapitalM(l))= HydroParam%uxy(HydroParam%ElCapitalM(l),1,iEdge)
		!	HydroParam%ug(iEdge,HydroParam%ElCapitalM(l)+1)= HydroParam%uxy(HydroParam%ElCapitalM(l),1,iEdge)
		!	
		!	HydroParam%vg(iEdge,HydroParam%ElCapitalM(l))= HydroParam%uxy(HydroParam%ElCapitalM(l),2,iEdge)
		!	HydroParam%vg(iEdge,HydroParam%ElCapitalM(l)+1)= HydroParam%uxy(HydroParam%ElCapitalM(l),2,iEdge)
		!	
		!	If (r==0.or.r==l) Then 
		!		HydroParam%wg(iEdge,HydroParam%ElCapitalM(l))= HydroParam%w(HydroParam%ElCapitalM(l),l)
		!		HydroParam%wg(iEdge,HydroParam%ElCapitalM(l)+1)= HydroParam%w(HydroParam%ElCapitalM(l)+1,l)
		!	else
		!		HydroParam%wg(iEdge,HydroParam%ElCapitalM(l)) = 0.5d0*(HydroParam%w(HydroParam%ElCapitalM(l),r)+HydroParam%w(HydroParam%ElCapitalM(l),l))
		!		HydroParam%wg(iEdge,HydroParam%ElCapitalM(l)+1) = 0.5d0*(HydroParam%w(HydroParam%ElCapitalM(l)+1,r)+HydroParam%w(HydroParam%ElCapitalM(l)+1,l))
		!	EndIf
		!cycle
  !      endif
		!
		!HydroParam%ug(iEdge,HydroParam%CapitalM(iEdge)+1)= HydroParam%uxy(HydroParam%CapitalM(iEdge),1,iEdge)!+weitW*((HydroParam%uxy(HydroParam%ElCapitalM(l),1,iEdge)-HydroParam%uxy(HydroParam%ElCapitalM(l)-1,1,iEdge))/(weit))
		!HydroParam%vg(iEdge,HydroParam%CapitalM(iEdge)+1)= HydroParam%uxy(HydroParam%CapitalM(iEdge),2,iEdge)!+weitW*((HydroParam%uxy(HydroParam%ElCapitalM(l),2,iEdge)-HydroParam%uxy(HydroParam%ElCapitalM(l)-1,2,iEdge))/(weit))
		!
		!Do iLayer = HydroParam%Smallms(iEdge), HydroParam%CapitalM(iEdge) !HydroParam%Smallms(iEdge)+1
	 !
		!	if (r/=0) then
		!		if (HydroParam%ElSmallm(r)<HydroParam%ElSmallm(l)) then
		!			if (ilayer<HydroParam%ElSmallm(l)) then
		!				l=r
		!			else
		!				l = MeshParam%Left(iEdge)
		!				r = MeshParam%Right(iEdge)
		!			endif
		!		endif
		!	endif
		!	
		!	if (HydroParam%ElSmallm(l)==1) then
		!		HydroParam%ug(iEdge,iLayer)= HydroParam%uxy(iLayer,1,iEdge)
		!		HydroParam%vg(iEdge,iLayer)= HydroParam%uxy(iLayer,2,iEdge)
		!		
		!		If (r==0.or.r==l) Then 
		!			HydroParam%wg(iEdge,iLayer)= HydroParam%w(iLayer,l)
		!		else
		!			If (iLayer>=HydroParam%ElSmallm(r)) Then
		!				HydroParam%wg(iEdge,iLayer) = 0.5d0*(HydroParam%w(iLayer,r)+HydroParam%w(iLayer,l))
		!			else
		!				HydroParam%wg(iEdge,iLayer)= HydroParam%w(iLayer,l)
		!			endif
		!		EndIf
		!		
		!	Else
		!		!Find Velocities to K+1/2 for each edge
		!		if (abs(HydroParam%uxy(iLayer,1,iEdge))>0.and.abs(HydroParam%uxy(iLayer-1,1,iEdge))>0) then !flux in U direction
		!			HydroParam%ug(iEdge,iLayer)= HydroParam%uxy(iLayer-1,1,iEdge)+(((HydroParam%Zb(iLayer,l)-HydroParam%Ze(iLayer,l))/(HydroParam%Zb(iLayer,l)-HydroParam%Zb(iLayer-1,l)))*(HydroParam%uxy(iLayer,1,iEdge)-HydroParam%uxy(iLayer-1,1,iEdge)))
		!		else
		!			HydroParam%ug(iEdge,iLayer)= HydroParam%uxy(iLayer,1,iEdge)
		!		endif	
		!		if (abs(HydroParam%uxy(iLayer,2,iEdge))>0.and.abs(HydroParam%uxy(iLayer-1,2,iEdge))>0) then !flux in V direction
		!			HydroParam%vg(iEdge,iLayer)= HydroParam%uxy(iLayer-1,2,iEdge)+(((HydroParam%Zb(iLayer,l)-HydroParam%Ze(iLayer,l))/(HydroParam%Zb(iLayer,l)-HydroParam%Zb(iLayer-1,l)))*(HydroParam%uxy(iLayer,2,iEdge)-HydroParam%uxy(iLayer-1,2,iEdge)))
		!		else
		!			HydroParam%vg(iEdge,iLayer)= HydroParam%uxy(iLayer,2,iEdge)
		!		endif
		!		
		!		If (r==0.or.r==l) Then 
		!			HydroParam%wg(iEdge,iLayer)= HydroParam%w(iLayer,l)
		!		else
		!			If (iLayer>=HydroParam%ElSmallm(r)) Then
		!				HydroParam%wg(iEdge,iLayer) = 0.5d0*(HydroParam%w(iLayer,r)+HydroParam%w(iLayer,l))
		!			else
		!				HydroParam%wg(iEdge,iLayer)= HydroParam%w(iLayer,l)
		!			endif
		!		EndIf
		!	endif
		!EndDo  
		!
		!If (r==0.or.r==l) Then 
		!	HydroParam%wg(iEdge,HydroParam%ElCapitalM(l)+1)= HydroParam%w(HydroParam%ElCapitalM(l)+1,l)
		!else
		!	HydroParam%wg(iEdge,HydroParam%ElCapitalM(l)+1) = 0.5d0*(HydroParam%w(HydroParam%ElCapitalM(l)+1,r)+HydroParam%w(HydroParam%ElCapitalM(l)+1,l))
		!EndIf
  !              
  !
  !  EndDo !i=1,ns
  !  
  !  Do iEdge=1,MeshParam%nEdge
  !      l = MeshParam%Left(iEdge)
  !      r = MeshParam%Right(iEdge)
  !      Do iLayer = 1, HydroParam%CapitalM(iEdge)
  !          if (isnan(HydroParam%ug(iEdge,iLayer)).or.isnan(HydroParam%vg(iEdge,iLayer)).or.isnan(HydroParam%wg(iEdge,iLayer))) then
  !              continue
  !          endif
  !      enddo
  !  enddo
    
    
    ! 6. Finding Vertice-centered Velocity Components
    HydroParam%ubV = 0.
  
    Do iNode=1,MeshParam%nNode
        Do iLayer = 1,MeshParam%Kmax
            weit=0
            weitU=0
            weitV=0
            weitW=0
            icountU=0
            icountV=0
            icountW=0
            ! 5.1. Find U and V Node velocities
            Do j=1,MeshParam%nEgdesatNode(iNode)
                Face = MeshParam%EgdesatNode(j,iNode)
                r = MeshParam%Right(Face)
                If(abs(HydroParam%uxy(iLayer,1,Face))>NearZero) Then
                    HydroParam%ubV(iLayer,1,iNode)=HydroParam%ubV(iLayer,1,iNode)+(HydroParam%uxy(iLayer,1,Face))/MeshParam%EdgeLength(Face)
                    icountU=icountU+1
                    weitU=weitU+1/MeshParam%EdgeLength(Face) 
                EndIf
                If(abs(HydroParam%uxy(iLayer,2,Face))>NearZero) Then
                    HydroParam%ubV(iLayer,2,iNode)=HydroParam%ubV(iLayer,2,iNode)+(HydroParam%uxy(iLayer,2,Face))/MeshParam%EdgeLength(Face)
                    icountV=icountV+1
                    weitV=weitV+1/MeshParam%EdgeLength(Face) 
                EndIf
            EndDo
                
            ! 5.3. Find W node velocity
            Do j=1,MeshParam%nEgdesatNode(iNode)
                Face = MeshParam%EgdesatNode(j,iNode)
                r = MeshParam%Right(Face)
                if(abs(HydroParam%wfc(iLayer,Face))>NearZero) Then
                    HydroParam%ubV(iLayer,3,iNode)=HydroParam%ubV(iLayer,3,iNode)+HydroParam%wfc(iLayer,Face)/MeshParam%EdgeLength(Face)
                    icountW=icountW+1 
                    weitW=weitW+1/MeshParam%EdgeLength(Face)
                EndIf
            EndDo
  
                If(icountU.ne.0) Then
                    HydroParam%ubV(iLayer,1,iNode) = HydroParam%ubV(iLayer,1,iNode)/weitU
                    !!Considering inflo/ouflow condition
                    !Do j=1,MeshParam%nEgdesatNode(iNode)
                    !    Face = MeshParam%EgdesatNode(j,iNode)
                    !    if (HydroParam%IndexWaterLevelEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) then 
                    !        HydroParam%ubV(iLayer,1,iNode)=(HydroParam%uxy(iLayer,1,Face))
                    !        exit
                    !    endif
                    !enddo
                else
                    HydroParam%ubV(iLayer,1,iNode) = 0.d0
                endif
                if(icountV.ne.0) Then
                    HydroParam%ubV(iLayer,2,iNode) = HydroParam%ubV(iLayer,2,iNode)/weitV
                    !Considering inflo/ouflow condition
                    !Do j=1,MeshParam%nEgdesatNode(iNode)
                    !    Face = MeshParam%EgdesatNode(j,iNode)
                    !    if (HydroParam%IndexWaterLevelEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) then 
                    !        HydroParam%ubV(iLayer,2,iNode)=(HydroParam%uxy(iLayer,2,Face))
                    !        exit
                    !    endif
                    !enddo
                else
                    HydroParam%ubV(iLayer,2,iNode) = 0.d0
                endif
                if(icountW.ne.0) Then
                    HydroParam%ubV(iLayer,3,iNode) = HydroParam%ubV(iLayer,3,iNode)/weitW
                    !Considering inflo/ouflow condition
                    !Do j=1,MeshParam%nEgdesatNode(iNode)
                    !    Face = MeshParam%EgdesatNode(j,iNode)
                    !    if (HydroParam%IndexWaterLevelEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) then 
                    !        HydroParam%ubV(iLayer,3,iNode)= HydroParam%wfc(iLayer,Face)
                    !        exit
                    !    endif
                    !enddo
                else
                    HydroParam%ubV(iLayer,3,iNode) = 0.d0
                EndIf
            
        EndDo
    EndDo   
    
    Do iNode=1,MeshParam%nNode
        Do iLayer = 1,MeshParam%Kmax
            if (isnan(HydroParam%ubV(iLayer,1,iNode)).or.isnan(HydroParam%ubV(iLayer,2,iNode)).or.isnan(HydroParam%ubV(iLayer,3,iNode))) then
                continue
            endif
        enddo
    enddo
    
    ! 5. Finding Nodal Velocities
    HydroParam%uNode = 0.d0
    Do iNode=1,MeshParam%nNode
        !if (iNode==487.or.iNode==511) then
        !    continue 
        !endif
        Do iLayer = 1,MeshParam%Kmax+1
            weit=0
            weitU=0
            weitV=0
            weitW=0
            icountU=0
            icountV=0
            icountW=0
            if (iLayer < MeshParam%Kmax+1) then
                ! 5.1. Find U and V Node velocities
                Do j=1,MeshParam%nEgdesatNode(iNode)
                    Face = MeshParam%EgdesatNode(j,iNode)
                    r = MeshParam%Right(Face)
                    If(abs(HydroParam%ug(Face,iLayer))>NearZero) Then
                        HydroParam%uNode(iLayer,1,iNode)=HydroParam%uNode(iLayer,1,iNode)+HydroParam%ug(Face,iLayer)/MeshParam%EdgeLength(Face)
                        icountU=icountU+1
                        weitU=weitU+1/MeshParam%EdgeLength(Face) 
                    EndIf 
                    If(abs(HydroParam%vg(Face,iLayer))>NearZero) Then
                        HydroParam%uNode(iLayer,2,iNode)=HydroParam%uNode(iLayer,2,iNode)+HydroParam%vg(Face,iLayer)/MeshParam%EdgeLength(Face)
                        icountV=icountV+1
                        weitV=weitV+1/MeshParam%EdgeLength(Face) 
                    EndIf  
                EndDo
                
                ! 5.3. Find W node velocity
                Do j=1,MeshParam%nEgdesatNode(iNode)
                    Face = MeshParam%EgdesatNode(j,iNode)
                    if(abs(HydroParam%wg(Face,iLayer))>NearZero) Then
                        HydroParam%uNode(iLayer,3,iNode)=HydroParam%uNode(iLayer,3,iNode)+HydroParam%wg(Face,iLayer)/MeshParam%EdgeLength(Face)
                        icountW=icountW+1 
                        weitW=weitW+1/MeshParam%EdgeLength(Face)
                    EndIf
                EndDo
                
                If(icountU.ne.0) Then
                    HydroParam%uNode(iLayer,1,iNode) = HydroParam%uNode(iLayer,1,iNode)/weitU
                    !Do j=1,MeshParam%nEgdesatNode(iNode)
                    !    Face = MeshParam%EgdesatNode(j,iNode)
                    !    if (HydroParam%IndexWaterLevelEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) then 
                    !        HydroParam%uNode(iLayer,1,iNode)=HydroParam%ug(Face,iLayer)
                    !        exit
                    !    endif
                    !enddo
                else
                    HydroParam%uNode(iLayer,1,iNode) = 0.d0
                endif
                if(icountV.ne.0) Then
                    HydroParam%uNode(iLayer,2,iNode) = HydroParam%uNode(iLayer,2,iNode)/weitV
                    !considering inflow/outflow condition
                    !Do j=1,MeshParam%nEgdesatNode(iNode)
                    !    Face = MeshParam%EgdesatNode(j,iNode)
                    !    if (HydroParam%IndexWaterLevelEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) then 
                    !        HydroParam%uNode(iLayer,2,iNode)=HydroParam%vg(Face,iLayer)
                    !        exit
                    !    endif
                    !enddo
                else
                    HydroParam%uNode(iLayer,2,iNode) = 0.d0
                endif
                
                if(icountW.ne.0) Then
                    HydroParam%uNode(iLayer,3,iNode) = HydroParam%uNode(iLayer,3,iNode)/weitW
                    !Do j=1,MeshParam%nEgdesatNode(iNode)
                    !    Face = MeshParam%EgdesatNode(j,iNode)
                    !    if (HydroParam%IndexWaterLevelEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) then 
                    !        HydroParam%uNode(iLayer,3,iNode)=HydroParam%wg(Face,iLayer)
                    !        exit
                    !    endif
                    !enddo
                else
                    HydroParam%uNode(iLayer,3,iNode) = 0.d0
                EndIf
            else !for kMax+1 Layer
                Do j=1,MeshParam%nEgdesatNode(iNode)
                    Face = MeshParam%EgdesatNode(j,iNode)
                    r = MeshParam%Right(Face)
                    l = MeshParam%Left(Face)
                   if (HydroParam%ElCapitalM(l)==HydroParam%ElSmallms(l)) then
                        HydroParam%uNode(iLayer,1,iNode)=HydroParam%ubV(iLayer-1,1,iNode)!+(weitW)*((HydroParam%ubV(iLayer-1,1,iNode)-HydroParam%ubV(iLayer-2,1,iNode))/(weit))
                        HydroParam%uNode(iLayer,2,iNode)=HydroParam%ubV(iLayer-1,2,iNode)!+(weitW)*((HydroParam%ubV(iLayer-1,2,iNode)-HydroParam%ubV(iLayer-2,2,iNode))/(weit))       
                   else !if the domain has more then one layer, use lagrange polinomial interpolation
                        HydroParam%uNode(iLayer,1,iNode)=HydroParam%ubV(iLayer-1,1,iNode)
                        HydroParam%uNode(iLayer,2,iNode)=HydroParam%ubV(iLayer-1,2,iNode)                        
                    endif
                EndDo
                
                weitW=0
                Do j=1,MeshParam%nEgdesatNode(iNode)
                    Face = MeshParam%EgdesatNode(j,iNode)
                    r = MeshParam%Right(Face)
                    if(abs(HydroParam%wg(Face,iLayer))>NearZero) Then
                        HydroParam%uNode(iLayer,3,iNode)=HydroParam%uNode(iLayer,3,iNode)+HydroParam%wg(Face,iLayer)/MeshParam%EdgeLength(Face)
                        icountW=icountW+1 
                        weitW=weitW+1/MeshParam%EdgeLength(Face)
                    EndIf
                    
                EndDo
                if(icountW.ne.0) Then
                    HydroParam%uNode(iLayer,3,iNode) = HydroParam%uNode(iLayer,3,iNode)/weitW
                    !Do j=1,MeshParam%nEgdesatNode(iNode)
                    !    Face = MeshParam%EgdesatNode(j,iNode)
                    !    if (HydroParam%IndexWaterLevelEdge(Face)>0.or.HydroParam%IndexInflowEdge(Face)>0) then 
                    !        HydroParam%uNode(iLayer,3,iNode)=HydroParam%wg(Face,iLayer)
                    !        exit
                    !    endif
                    !enddo
                else
                    HydroParam%uNode(iLayer,3,iNode)=0.0d0
                EndIf 
            endif
        EndDo
    EndDo   
    
    Do iNode=1,MeshParam%nNode
        Do iLayer = 1,MeshParam%Kmax+1
            if (isnan(HydroParam%uNode(iLayer,1,iNode)).or.isnan(HydroParam%uNode(iLayer,2,iNode)).or.isnan(HydroParam%uNode(iLayer,3,iNode))) then
                continue
            endif
        enddo
    enddo
    
	Return
    End Subroutine VelocitiesSUB  
    
       
    Function neighElem(MeshParam,Face)

    use MeshVars
    Implicit none
    Integer:: Face, neighElem
    type(MeshGridParam) :: MeshParam
    neighElem = MeshParam%Right(Face)
    If (neighElem == 0) Then
        neighElem = MeshParam%Left(Face)
    EndIf
    
    End Function neighElem  