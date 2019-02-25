!>  This module declares and Computes Mesh Characteristics for Structured Finite Volume Framework
!> @attention List of Modifications: 
!> @parblock
!!  -> 10.03.2014: Routine Implementation               (Michael Dumbser)
!!
!!  -> 10.03.2014: Fortran Sintax                       (Rafael Cavalcanti)
!!
!!  -> 10.03.2014: Casulli's Semi-Implicit Features     (Rafael Cavalcanti)
!! 
!!  -> 25.03.2014: Edge Barycenter Correction           (Michael Dumbser)
!> @endparblock
!> @author Rafael Cavalcanti and  Ruberto C. Fragoso
!>\remark all arguments are passed through module dependencies:\n
!> <b>variables from SimulationModel:</b>\n
!! dx, dy \n
!> <b>variables from MeshVars :</b>\n
!! Quadri,QuadriSort,indPoint, nElem,nPoint,nNode,Neighbor,Edge,Area,
!! normalVector,inCircle,xNode,yNode,nVertexElem
!! KMAX,VertexElem,EdgeBary,Left, Right, nEdge,EdgeLength,CirDistance \n
!!
!>\note Vertex are counted from bottom downleft corner and following an anti-clockwise direction:\n
!>\ For nodes numeration, the Standard is this:
        !     4     3
        !     .-----.
        !     |     |
        !     |     |
        !     .-----.
        !     1     2
!>\ For nodes positions in the vector, the Standard is this:
        !     2     1
        !     .-----.
        !     |     |
        !     |     |
        !     .-----.
        !     3     4
    
Module MeshVars
    
    use domain_types

    implicit none
    type MeshGridParam
       ! 1. Mesh Features
        Real, Allocatable:: Incircle(:) !<  Diameter of the included circle of each element \f$\frac{4*Area(iElem)}{\sum{EdgeLength(Edge(:,iElem)}}\f$
        Integer, Allocatable:: Neighbor(:,:)  !< the 4 neighbor ID of each element allocated in quadriMesh.f90
        Real, Allocatable:: EdgeLength(:) !< Edge length of each edge allocated in quadriMesh.f90
        Real, Allocatable:: Area(:)       !< Area of each element allocated in quadriMesh.f90
        Real, Allocatable:: CirDistance(:) !<distance between the centers of the included circles
        Real, Allocatable:: NormalVector(:,:,:) !< component of the normal vector of each of the 4 edges for each element allocated in quadriMesh.f90
        Real, Allocatable:: NormalVectorEdge(:,:) !< component of the normal vector of each edges allocated in quadriMesh.f90
        Real, Allocatable:: TangentVector(:,:,:) !< component of the tangential vector of each of the 4 edges for each element allocated in quadriMesh.f90
        Real, Allocatable:: EdgeBary(:,:)  !< position of the barycentre of each adge allocated in quadriMesh.f90
        Integer, Allocatable:: Edge(:,:)  !< ID of the 4 edge of each element
        Integer:: EdgeDef(2,4)
        Integer, Allocatable:: Left(:)    !< ID of the left neighbour element for each edge
        Integer, Allocatable:: Right(:)   !<ID of the right neighbour element for each edge
        Integer, Allocatable:: nVertexElem(:)
        Integer, Allocatable:: VertexElem(:,:)
        Integer:: nEdge, nNode, nElem, nPoint
        Real, Allocatable:: xb(:), yb(:)
        Integer, Allocatable:: Quadri(:,:)
        Real, Allocatable :: xNode(:), yNode(:)
        Real, Allocatable :: xPoint(:), yPoint(:), zPoint(:)
        Integer, Allocatable:: indPoint(:,:),Connect(:),cell_type(:)
        Real, Allocatable:: LIMCAM(:) !< Layer levels (m)
        Real, Allocatable:: LIMCAMAUX(:) 
        Integer:: NCAMMAX !< Number of Vertical Layers
        Real:: zL !< y-cell resolution (m)
        Real:: zR !< y-cell resolution (m)
        Integer:: KMax !< Maximum Number of Vertical Layers
        Integer, Allocatable:: EdgeNodes(:,:)
        Integer, Allocatable:: EgdesatNode(:,:)
        Integer, Allocatable:: nEgdesatNode(:)
        Real:: dx !< x-cell resolution (m)
        Real:: dy !< y-cell resolution (m)
        Integer, Allocatable :: LeftAux(:), RightAux(:) 
        Real, Allocatable:: AbsDelta(:),  EdgeBaryAux(:,:) 
        Integer:: nMaxVertexElem,nMaxEgdesatNode
    
        !2. Grid Data
        Integer:: iWindRed
        Integer:: iWetland
        Integer:: id50
        Integer:: iOMfraction
        Integer:: iBedrock
        Integer:: ieta0
        Real, Allocatable:: CREDV(:)
        Integer, Allocatable:: BANHADO(:)
        Real, Allocatable:: d50(:)
        Real, Allocatable:: OMfraction(:)
        Real, Allocatable:: eta0(:)
        
        Integer, Allocatable:: trv(:),NEIBN(:),NEIBL(:),NEIBS(:),NEIBO(:)   
        

        
        
    contains
        procedure :: initializeMeshParam
    end type
    
    contains
    
    subroutine initializeMeshParam(this, sim, StructuredMeshFeatures)
    
        Integer:: iElem,iEdge,iLayer,iNode,kNode,iNode1,iNode2
        type(Simulation), intent(in) :: sim
        class(MeshGridParam) :: this
        type(StructuredMesh) :: StructuredMeshFeatures    
        character, pointer :: simulationLabel(:)
        Real, pointer :: layers(:)
        Real, pointer :: xCoordinates(:)
        Real, pointer :: yCoordinates(:)
        integer(c_long_long), pointer :: verticeIds(:)
        Integer:: i,j,jEdge,jNode1,jNode2
        Integer:: nMaxVertexElem,nMaxEgdesatNode
        Real:: n3d(3), tempvec(3), Check(2), v(4,3), zn(3), Norm_n3d, CheckMesh, Det
        Integer:: l, r
    
        call c_f_pointer(StructuredMeshFeatures%xCoordinates, xCoordinates, [StructuredMeshFeatures%numberOfElements])
        call c_f_pointer(StructuredMeshFeatures%yCoordinates, yCoordinates, [StructuredMeshFeatures%numberOfElements])
        call c_f_pointer(StructuredMeshFeatures%verticeIds, verticeIds, [StructuredMeshFeatures%verticeIdsLength])
        call c_f_pointer(sim%label, simulationLabel, [sim%labelLength])
        call c_f_pointer(sim%layers, layers, [sim%layersLength])

        this%nElem = StructuredMeshFeatures%numberOfElements 
        
        this%dx = StructuredMeshFeatures%resolution
        !Lock horizontal resolution test
        !this%dx = 0.025d0
        !!
        this%dy = this%dx
        
    
        Allocate(this%xb(this%nElem))
        Allocate(this%yb(this%nElem))
        Allocate(this%Quadri(4,this%nElem))
        
        Do iElem = 1, this%nElem
            this%xb(iElem) = xCoordinates(iElem)
            this%yb(iElem) = yCoordinates(iElem)
            this%Quadri(3,iElem) = verticeIds(4*(iElem-1)+1)
            this%Quadri(4,iElem) = verticeIds(4*(iElem-1)+2)
            this%Quadri(1,iElem) = verticeIds(4*(iElem-1)+3)
            this%Quadri(2,iElem) = verticeIds(4*(iElem-1)+4)
        EndDo
    
    
        !transferir para Module Mesh
        this%NCAMMAX = sim%layersLength                      ! Number of Vertical Layers
        this%zL = sim%minimumVerticalLimit-0.001
        this%zR = sim%maximumVerticalLimit
    
        If (this%NCAMMAX > 0) Then
            ALLOCATE (this%LIMCAM(this%NCAMMAX))
            ALLOCATE (this%LIMCAMAUX(this%NCAMMAX+1))
            Do i = 1,this%NCAMMAX
                this%LIMCAM(i) = layers(i)                   ! Layer levels (m)
            EndDo
        Else
            ALLOCATE (this%LIMCAM(1))
            ALLOCATE (this%LIMCAMAUX(1))
            this%LIMCAM = this%zL
            this%LIMCAMAUX = this%zL
        EndIf            
        
        
        ! Number of Vertical Layers
        If (this%NCAMMAX ==0) then !Two-dimensional
            this%KMax = 1
            this%LIMCAMAUX(this%KMax)=this%zL
        ElseIf (this%NCAMMAX >= 1) then !Three-dimensional
            this%KMax = 1
            Do iLayer = 1,this%NCAMMAX
                If (iLayer==1) then
                    If (this%LIMCAM(iLayer)>this%zL.and.this%LIMCAM(iLayer)<this%zR) then
                        this%KMax = this%KMax + 1
                        this%LIMCAMAUX(this%KMax-1) = this%LIMCAM(iLayer)
                    EndIf
                Elseif (iLayer==this%NCAMMAX) then
                    If (this%LIMCAM(iLayer)>this%zL.and.this%LIMCAM(iLayer)<this%LIMCAM(iLayer-1)) then
                        this%KMax = this%KMax + 1
                        this%LIMCAMAUX(this%KMax-1) = this%LIMCAM(iLayer)
                        !LIMCAMAUX(KMax-1) = zL
                    Else
                        !KMax = KMax + 1
                        !LIMCAMAUX(KMax-1) = LIMCAM(iLayer)
                    EndIf
                Else
                    If (this%LIMCAM(iLayer)>this%zL.and.this%LIMCAM(iLayer)<this%LIMCAM(iLayer-1)) then
                        this%KMax = this%KMax + 1
                        this%LIMCAMAUX(this%KMax-1) = this%LIMCAM(iLayer)
                    EndIf
                EndIf
            EndDo
            this%LIMCAMAUX(this%KMax) = this%zL
        Else
            Stop 'Incorrect number of Layers'
        Endif    
        
        this%nNode = MAXVAL(this%Quadri) + 1
        this%nPoint = this%nNode
        Allocate (this%Neighbor(4,this%nElem))
        Allocate (this%Edge(4,this%nElem))
        Allocate (this%EdgeBaryAux(2,4*this%nElem))
        Allocate (this%LeftAux(4*this%nElem))
        Allocate (this%RightAux(4*this%nElem))
        Allocate (this%Area(this%nElem))
        Allocate (this%NormalVector(2,4,this%nElem))
        Allocate (this%TangentVector(2,4,this%nElem))
        Allocate (this%Incircle(this%nElem))
        Allocate (this%xNode(this%nNode))
        Allocate (this%yNode(this%nNode))
        Allocate (this%nVertexElem(this%nNode))

        
        ! 1. Defining Mesh Vertex
        ! Using NC and NL we can define a Mathematical Equation to save the vertex for each node without repetition
        ! For nodes numeration, the Standard is this:
        !     4     3
        !     .-----.
        !     |     |
        !     |     |
        !     .-----.
        !     1     2
        ! For nodes positions in the vector, the Standard is this:
        !     2     1
        !     .-----.
        !     |     |
        !     |     |
        !     .-----.
        !     3     4
    
        Do iElem = 1, this%nElem
            this%xNode( this%Quadri(1,iElem)+1 ) = this%xb(iElem) + this%DX/2
            this%yNode( this%Quadri(1,iElem)+1 ) = this%yb(iElem) + this%DY/2
            this%xNode( this%Quadri(2,iElem)+1 ) = this%xb(iElem) - this%DX/2
            this%yNode( this%Quadri(2,iElem)+1 ) = this%yb(iElem) + this%DY/2
            this%xNode( this%Quadri(3,iElem)+1 ) = this%xb(iElem) - this%DX/2
            this%yNode( this%Quadri(3,iElem)+1 ) = this%yb(iElem) - this%DY/2
            this%xNode( this%Quadri(4,iElem)+1 ) = this%xb(iElem) + this%DX/2
            this%yNode( this%Quadri(4,iElem)+1 ) = this%yb(iElem) - this%DY/2
        EndDo
        
        ! 2. Defining The Neighbors
        
        !Counting the Elements Attached to Each Node - Voronoi Neighbors
        this%nVertexElem = 0
        Do iElem = 1, this%nElem
            Do kNode = 1, 4
                iNode = this%Quadri(kNode,iElem) + 1
                this%nVertexElem(iNode) = this%nVertexElem(iNode) + 1 
            EndDo
        EndDo
        this%nMaxVertexElem = MaxVal(this%nVertexElem)
        Allocate (this%VertexElem(this%nMaxVertexElem,this%nNode))
        !Define Elements Attached to Each Node
        this%nVertexElem = 0
        this%VertexElem  = 0
        Do iElem = 1, this%nElem
            Do kNode = 1,4
                iNode = this%Quadri(kNode,iElem) + 1
                this%nVertexElem(iNode) = this%nVertexElem(iNode) + 1 
                this%VertexElem(this%nVertexElem(iNode),iNode) = iElem 
            EndDo
        EndDo        
        ! 2.3 Compute the Neighbors of Each Quadrilateral
        this%Neighbor = 0 
        this%EdgeBary = 0 
        this%Edge     = 0 
        this%Left     = 0 
        this%Right    = 0
        ! 2.3.1 Edge Definition
        ! Edge I   of Quadrilateral i is Composed of Local Nodes 1 and 2
        this%EdgeDef(1,1) = 1
        this%EdgeDef(2,1) = 2
        ! Edge II  of Quadrilateral i is Composed of Local Nodes 2 and 3
        this%EdgeDef(1,2) = 2
        this%EdgeDef(2,2) = 3
        ! Edge III of Quadrilateral i is Composed of Local Nodes 3 and 4
        this%EdgeDef(1,3) = 3
        this%EdgeDef(2,3) = 4
        ! Edge IV  of Quadrilateral i is Composed of Local Nodes 4 and 1
        this%EdgeDef(1,4) = 4
        this%EdgeDef(2,4) = 1
        ! 2.3.2 Start the Neighbor Search
        this%nEdge = 0 
        Do iElem = 1, this%nElem
            Do iEdge = 1,4
                iNode1 = this%Quadri(this%EdgeDef(1,iEdge),iElem) + 1         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                iNode2 = this%Quadri(this%EdgeDef(2,iEdge),iElem) + 1         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                ! Look for a Potential Neighbor Element
                Do kNode = 1,this%nVertexElem(iNode1)            ! It could be iNode2
                    j = this%VertexElem(kNode,iNode1) 
                    If(iElem == j) Then
                        Continue
                    EndIf
                    Do jEdge = 1,4
                        jNode1 = this%Quadri(this%EdgeDef(1,jEdge),j) + 1 
                        jNode2 = this%Quadri(this%EdgeDef(2,jEdge),j) + 1 
                        If (iNode1 == jNode2 .AND. iNode2 == jNode1) Then
                            this%Neighbor(iEdge,iElem) = j 
                            If ( this%Edge(iEdge,iElem)==0 ) Then
                                this%nEdge = this%nEdge + 1 
                                this%Edge(iEdge,iElem) = this%nEdge 
                                this%Edge(jEdge,j)     = this%nEdge 
                                this%LeftAux(this%nEdge)  = iElem 
                                this%RightAux(this%nEdge) = j 
                                this%EdgeBaryAux(1,this%nEdge) = (this%xNode(iNode1)+this%xNode(iNode2))/2
                                this%EdgeBaryAux(2,this%nEdge) = (this%yNode(iNode1)+this%yNode(iNode2))/2
                            EndIf
                            Cycle           ! Check it With Michael
                        EndIf
                    EndDo
                    If (this%Neighbor(iEdge,iElem) /= 0) Then
                        Exit         ! If Neighbor is Found, Exit! ==== Check it With Michael
                    EndIf
                EndDo
                If (this%Neighbor(iEdge,iElem) == 0) Then
                    ! We are at the Boundary. So Edge should be zero. Add the edge.
                    If (this%Edge(iEdge,iElem) /= 0) Then
                        Print*, 'Impossible bug!!'
                        Print*, 'TriMeshCasulli - Ln 103'
                        Pause
                        Stop
                    EndIf
                    this%nEdge = this%nEdge + 1 
                    this%Edge(iEdge,iElem) = this%nEdge 
                    this%LeftAux(this%nEdge)  = iElem 
                    this%RightAux(this%nEdge) = 0 
                    this%EdgeBaryAux(1,this%nEdge) = (this%xNode(iNode1)+this%xNode(iNode2))/2                
                    this%EdgeBaryAux(2,this%nEdge) = (this%yNode(iNode1)+this%yNode(iNode2))/2
                EndIf
            EndDo       ! Do iEdge = 1,3
        EndDo       ! Do iElem = 1, NUMCEL
    
        Allocate(this%Left(this%nEdge))
        Allocate(this%Right(this%nEdge))
        Allocate(this%EdgeLength(this%nEdge))
        Allocate(this%CirDistance(this%nEdge))
        Allocate(this%AbsDelta(this%nEdge))
        Allocate(this%EdgeBary(2,this%nEdge))
        Allocate(this%EdgeNodes(2,this%nEdge))
        Allocate(this%NormalVectorEdge(2,this%nEdge))
        
        this%Left     = this%LeftAux(1:this%nEdge) 
        this%Right    = this%RightAux(1:this%nEdge) 
        this%EdgeBary = this%EdgeBaryAux(1:2,1:this%nEdge)   
        
        ! 3. Compute Edge Lengths, Incircle Diameter, and Barycenter Positions
        zn(1) = 0
        zn(2) = 0
        zn(3) = 1
    
        Do iElem = 1, this%nElem
            check = 0.
            Do iEdge = 1, 4
                j = this%Edge(iEdge,iElem) 
                iNode1 = this%Quadri(this%EdgeDef(1,iEdge),iElem) + 1         ! Global Node Number of Node 1 of Edge "iEdge" of Element i
                iNode2 = this%Quadri(this%EdgeDef(2,iEdge),iElem) + 1         ! Global Node Number of Node 2 of Edge "iEdge" of Element i
                this%EdgeNodes(1,j)=iNode1
                this%EdgeNodes(2,j)=iNode2
                v(iEdge,1) = this%xNode(iNode2) - this%xNode(iNode1)
                v(iEdge,2) = this%yNode(iNode2) - this%yNode(iNode1)
                v(iEdge,3) = 0
                Call Norm(v(iEdge,:),3,this%EdgeLength(j))
                this%TangentVector(1:2,iEdge,iElem) = v(iEdge,1:2)/this%EdgeLength(j) 
                ! Checking the Mesh Quality
                Call Cross(v(iEdge,:),zn,n3d)
                Call Norm(n3d,3,Norm_n3d)
                this%NormalVector(1:2,iEdge,iElem) = n3d(1:2)/Norm_n3d                  ! Outward-Pointing Unit Normal Vector
                this%NormalVectorEdge(1:2,j) = this%NormalVector(1:2,iEdge,iElem)
                Check = Check + this%NormalVector(:,iEdge,iElem)*this%EdgeLength(j)          ! This Quantity has to be Zero!!
            EndDo
            Call Norm(Check,2,CheckMesh)
            If (CheckMesh > 1e-12) Then
                Print*, 'Mesh is not consistent! Elem = ', iElem
                Print*, 'TriMeshCasulli - Ln 145'
                Pause
                Stop
            EndIf
            ! Area of the Quadrilateral
            !Call Cross(v(1,:),-v(3,:),tempvec)
            !tempvec = tempvec/2.
            
            this%Area(iElem) = Det(this%xNode(this%Quadri(1,iElem)+1),this%xNode(this%Quadri(2,iElem)+1),this%xNode(this%Quadri(3,iElem)+1),this%yNode(this%Quadri(1,iElem)+1),this%yNode(this%Quadri(2,iElem)+1),this%yNode(this%Quadri(3,iElem)+1))
            this%Area(iElem) = this%Area(iElem) + Det(this%xNode(this%Quadri(3,iElem)+1),this%xNode(this%Quadri(4,iElem)+1),this%xNode(this%Quadri(1,iElem)+1),this%yNode(this%Quadri(3,iElem)+1),this%yNode(this%Quadri(4,iElem)+1),this%yNode(this%Quadri(1,iElem)+1))
            !this%Area(iElem) = this%DX*this%DY
            If (this%Area(iElem) < 0) Then
                Print*, 'The nodes are not anti-clockwise! Elem = ', iElem
                Print*,  'TriMeshCasulli - Ln 152'
                Pause
                Stop
            EndIf
            ! Circumcenter of the Triangle
            this%Incircle(iElem) = 4.*this%Area(iElem)/sum(this%EdgeLength(this%Edge(:,iElem))) 
        EndDo       ! iElem = 1: NUMCEL        
        
        this%nMaxEgdesatNode = this%nMaxVertexElem + 1
        
        Allocate(this%EgdesatNode(this%nMaxEgdesatNode,this%nNode)) 
        Allocate(this%nEgdesatNode(this%nNode))
        this%nEgdesatNode = 0
        this%EgdesatNode = 0
        Do iEdge = 1,this%nEdge
            Do j = 1,2
              this%nEgdesatNode(this%EdgeNodes(j,iEdge)) = this%nEgdesatNode(this%EdgeNodes(j,iEdge)) + 1
              this%EgdesatNode(this%nEgdesatNode(this%EdgeNodes(j,iEdge)),this%EdgeNodes(j,iEdge))  = iEdge
            EndDo
        EndDo        
        
        
        ! 4. Compute Circumcenter Distance
        Do iEdge = 1,this%nEdge
            l = this%Left(iEdge)
            r = this%Right(iEdge)
            If (r==0) Then
                this%CirDistance(iEdge)  = 2*sqrt((this%EdgeBary(1,iEdge) - this%xb(l))**2 + (this%EdgeBary(2,iEdge) - this%yb(l))**2)
            Else
                this%CirDistance(iEdge)  = sqrt((this%xb(r) - this%xb(l))**2 + (this%yb(r) - this%yb(l))**2)
            EndIf
        EndDo
    
        ! Check if the Grid is Admissible
        this%AbsDelta = abs(this%CirDistance)
        If ( MinVal(this%AbsDelta) < 1e-14 ) Then
            Print*, 'Grid is not Admissible. You may have some Right Triangles... '
            Print*, 'TriMeshCasulli - Ln 184'
            Pause
            Stop
        Else
            Print*, 'Grid Seems to be Admissible'
        EndIf        
     
        
        Allocate(this%trv(this%nElem))
        Allocate(this%NEIBN(this%nElem))
        Allocate(this%NEIBL(this%nElem))
        Allocate(this%NEIBS(this%nElem))
        Allocate(this%NEIBO(this%nElem))
        
    end subroutine

End Module MeshVars
    
    