    Module ParticleTracking

    Implicit None

    type ParticleParam
        Integer:: npart !< Number of particle
        Integer, Allocatable:: npartElem(:) !< Number of cell which the particle is located
        Integer, Allocatable:: npartElem0(:) !< Number of particle initial cell
        Integer, Allocatable:: npartLayer(:) !< Number of layer which the particle is located
        Real, Allocatable:: xpart(:),ypart(:),zpart(:) !< Particle's position

        Integer:: ntarget !< Number of target cells
        Integer, Allocatable:: ntargetElem(:) !< Number of cell which the target is located
        Real, Allocatable:: RTime(:) !< Residence Time
        Integer, Allocatable:: FlagTarget(:) !< Flag for Residence Time

    contains
    procedure :: InitializeParticle
    end type

    !-----------------------------------------------------------------
    Contains

    Subroutine InitializeParticle(this,HydroParam,MeshParam)

    Use Hydrodynamic
    Use MeshVars

    class(ParticleParam) :: this
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam
    Integer:: i

    !< Defining number of particles
    this%npart = 27483

    Allocate(this%xpart(this%npart)) !< Particle actual latitude
    Allocate(this%ypart(this%npart)) !< Particle actual longitude
    Allocate(this%zpart(this%npart)) !< Particle actual depth
    Allocate(this%npartElem(this%npart)) !< Particle actual cell
    Allocate(this%npartElem0(this%npart)) !< Particle initial cell
    Allocate(this%npartLayer(this%npart)) !< Particle layer
    Allocate(this%RTime(this%npart)) !< Calculated Residence Time

    !< Read domain file

    Open(1001, File='Domain_Cells.txt', Action='Read', Status='Old') !< Domain where the particles will be released. File Format (3 columns): Cell(#) Long(UTM) LAT(UTM)

    Do i = 1,this%npart
        Read(1001,*) this%npartElem(i),this%xpart(i),this%ypart(i) !< Reading cell number and position
        this%npartElem0(i)=this%npartElem(i) !< Record initial cell
    EndDo
    close(1001)

    !< Defining target area
    this%ntarget = 46 !< Number of target cell

    Allocate(this%FlagTarget(this%npart)) !< Target cell reached: Flag = 1
    Allocate(this%ntargetElem(this%ntarget)) !< # of target cells

    this%FlagTarget = 0
    this%RTime = 0

    !< Read Target File
    Open(1001, File='Target_cells.txt', Action='Read', Status='Old') !< RT is calculated when the particles reaches the target cells. File format (1 column): Cell(#)
    Do i = 1,this%ntarget
        Read(1001,*) this%ntargetElem(i)
    EndDo
    close(1001)


    Return
    End Subroutine InitializeParticle

    Subroutine ParticlePosition(PartParam,HydroParam,MeshParam,dt)

    !
    ! Called in routine 0-MAIN
    !$ use omp_lib
    Use MeshVars
    Use Hydrodynamic
    Use SimulationModel
    Use LimnologyVars
    !Use LIB_VTK_IO

    Implicit none

    Integer:: iPart,nnnel,nel,nnel,knel,jlev,idt,ndelt,id0,j,jnel,nd,nn,lev,INOUT, tgtcell
    Real:: dtb,dt,csi,etta,zrat,zup
    Real:: upart,vpart,wpart
    Real:: x0,y0,z0
    Real:: xt,yt,zt
    Real:: xpoly(4),ypoly(4)
    Real :: vxn(8),vyn(8),vzn(8)
    Real :: staint(8)
    Integer:: iflqs1
    Real:: NearZero = 1e-10
    Integer:: Finaltime,IniTime, iTarget


    type(ParticleParam) :: PartParam
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam
    !type(SimulationParam) :: simParam

    ndelt = 10 !<substep time



    Do iPart=1,PartParam%npart
        If (PartParam%FlagTarget(iPart)==1) Then
            Cycle
        EndIf

        nnel = PartParam%npartElem(iPart) !< Actual cell
        jlev = PartParam%npartLayer(iPart) !< Actual layer

        !< Initial Position
        x0 = PartParam%xpart(iPart)
        y0 = PartParam%ypart(iPart)
        z0 = PartParam%zpart(iPart)

        call ibilinear(dble(nnel),MeshParam%EdgeBary(1,MeshParam%Edge(1,nnel)),MeshParam%EdgeBary(1,MeshParam%Edge(2,nnel)),MeshParam%EdgeBary(1,MeshParam%Edge(3,nnel)),MeshParam%EdgeBary(1,MeshParam%Edge(4,nnel)),MeshParam%EdgeBary(2,MeshParam%Edge(1,nnel)),MeshParam%EdgeBary(2,MeshParam%Edge(2,nnel)),MeshParam%EdgeBary(2,MeshParam%Edge(3,nnel)),MeshParam%EdgeBary(2,MeshParam%Edge(4,nnel)),x0,y0,csi,etta,staint)

        call ParticleVelocity(nnel,jlev,staint,z0,upart,vpart,wpart,HydroParam,MeshParam)


        Do idt=1,ndelt
            dtb = dt/ndelt
            xt=x0+dtb*upart
            yt=y0+dtb*vpart
            zt=z0+dtb*wpart

            xpoly(1)= MeshParam%xNode( MeshParam%Quadri(1,nnel)+1 )
            xpoly(2)= MeshParam%xNode( MeshParam%Quadri(2,nnel)+1 )
            xpoly(3)= MeshParam%xNode( MeshParam%Quadri(3,nnel)+1 )
            xpoly(4)= MeshParam%xNode( MeshParam%Quadri(4,nnel)+1 )
            ypoly(1)= MeshParam%yNode( MeshParam%Quadri(1,nnel)+1 )
            ypoly(2)= MeshParam%yNode( MeshParam%Quadri(2,nnel)+1 )
            ypoly(3)= MeshParam%yNode( MeshParam%Quadri(3,nnel)+1 )
            ypoly(4)= MeshParam%yNode( MeshParam%Quadri(4,nnel)+1 )

            Call PNPOLY(xt,yt,xpoly,ypoly,4,INOUT)

            If (INOUT==1.OR.INOUT==0) Then !Inside

                call ibilinear(dble(nnel),MeshParam%EdgeBary(1,MeshParam%Edge(1,nnel)),MeshParam%EdgeBary(1,MeshParam%Edge(2,nnel)),MeshParam%EdgeBary(1,MeshParam%Edge(3,nnel)),MeshParam%EdgeBary(1,MeshParam%Edge(4,nnel)),MeshParam%EdgeBary(2,MeshParam%Edge(1,nnel)),MeshParam%EdgeBary(2,MeshParam%Edge(2,nnel)),MeshParam%EdgeBary(2,MeshParam%Edge(3,nnel)),MeshParam%EdgeBary(2,MeshParam%Edge(4,nnel)),xt,yt,csi,etta,staint)

                call ParticleVelocity(nnel,jlev,staint,zt,upart,vpart,wpart,HydroParam,MeshParam)


            Else !Outside
                Do jnel=1,4
                    nel = MeshParam%Neighbor(jnel,nnel)
                    If (nel == 0) Then
                        continue
                        cycle
                    EndIf
                    xpoly(1)= MeshParam%xNode( MeshParam%Quadri(1,nel)+1 )
                    xpoly(2)= MeshParam%xNode( MeshParam%Quadri(2,nel)+1 )
                    xpoly(3)= MeshParam%xNode( MeshParam%Quadri(3,nel)+1 )
                    xpoly(4)= MeshParam%xNode( MeshParam%Quadri(4,nel)+1 )
                    ypoly(1)= MeshParam%yNode( MeshParam%Quadri(1,nel)+1 )
                    ypoly(2)= MeshParam%yNode( MeshParam%Quadri(2,nel)+1 )
                    ypoly(3)= MeshParam%yNode( MeshParam%Quadri(3,nel)+1 )
                    ypoly(4)= MeshParam%yNode( MeshParam%Quadri(4,nel)+1 )

                    Call PNPOLY(xt,yt,xpoly,ypoly,4,INOUT)

                    If (INOUT==1.OR.INOUT==0) Then !Inside
                        call ibilinear(dble(nel),MeshParam%EdgeBary(1,MeshParam%Edge(1,nel)),MeshParam%EdgeBary(1,MeshParam%Edge(2,nel)),MeshParam%EdgeBary(1,MeshParam%Edge(3,nel)),MeshParam%EdgeBary(1,MeshParam%Edge(4,nel)),MeshParam%EdgeBary(2,MeshParam%Edge(1,nel)),MeshParam%EdgeBary(2,MeshParam%Edge(2,nel)),MeshParam%EdgeBary(2,MeshParam%Edge(3,nel)),MeshParam%EdgeBary(2,MeshParam%Edge(4,nel)),xt,yt,csi,etta,staint)
                        call ParticleVelocity(nel,jlev,staint,zt,upart,vpart,wpart,HydroParam,MeshParam)
                        nnel = nel
                        Exit
                    Else ! Outside
                        Do knel=1,4
                            nnnel = MeshParam%Neighbor(knel,nel)
                            If (nnnel==nnel.or.nnnel==0) Then
                                Cycle
                            EndIf

                            xpoly(1)= MeshParam%xNode( MeshParam%Quadri(1,nnnel)+1 )
                            xpoly(2)= MeshParam%xNode( MeshParam%Quadri(2,nnnel)+1 )
                            xpoly(3)= MeshParam%xNode( MeshParam%Quadri(3,nnnel)+1 )
                            xpoly(4)= MeshParam%xNode( MeshParam%Quadri(4,nnnel)+1 )
                            ypoly(1)= MeshParam%yNode( MeshParam%Quadri(1,nnnel)+1 )
                            ypoly(2)= MeshParam%yNode( MeshParam%Quadri(2,nnnel)+1 )
                            ypoly(3)= MeshParam%yNode( MeshParam%Quadri(3,nnnel)+1 )
                            ypoly(4)= MeshParam%yNode( MeshParam%Quadri(4,nnnel)+1 )

                            Call PNPOLY(xt,yt,xpoly,ypoly,4,INOUT)

                            If (INOUT==1.OR.INOUT==0) Then !Inside
                                call ibilinear(dble(nnnel),MeshParam%EdgeBary(1,MeshParam%Edge(1,nnnel)),MeshParam%EdgeBary(1,MeshParam%Edge(2,nnnel)),MeshParam%EdgeBary(1,MeshParam%Edge(3,nnnel)),MeshParam%EdgeBary(1,MeshParam%Edge(4,nnnel)),MeshParam%EdgeBary(2,MeshParam%Edge(1,nnnel)),MeshParam%EdgeBary(2,MeshParam%Edge(2,nnnel)),MeshParam%EdgeBary(2,MeshParam%Edge(3,nnnel)),MeshParam%EdgeBary(2,MeshParam%Edge(4,nnnel)),xt,yt,csi,etta,staint)
                                call ParticleVelocity(nnnel,jlev,staint,zt,upart,vpart,wpart,HydroParam,MeshParam)
                                nnel = nnnel
                                Exit
                            EndIf
                        EndDo
                    EndIf
                EndDo
            Endif

            x0 = xt
            y0 = yt
            z0 = zt

        EndDo

        !Saving position of particle
        
        PartParam%npartElem(iPart) = nnel

        xpoly(1)= MeshParam%xNode( MeshParam%Quadri(1,nnel)+1 )
        xpoly(2)= MeshParam%xNode( MeshParam%Quadri(2,nnel)+1 )
        xpoly(3)= MeshParam%xNode( MeshParam%Quadri(3,nnel)+1 )
        xpoly(4)= MeshParam%xNode( MeshParam%Quadri(4,nnel)+1 )
        ypoly(1)= MeshParam%yNode( MeshParam%Quadri(1,nnel)+1 )
        ypoly(2)= MeshParam%yNode( MeshParam%Quadri(2,nnel)+1 )
        ypoly(3)= MeshParam%yNode( MeshParam%Quadri(3,nnel)+1 )
        ypoly(4)= MeshParam%yNode( MeshParam%Quadri(4,nnel)+1 )

        Call PNPOLY(xt,yt,xpoly,ypoly,4,INOUT)

        !< If the particle get out the domain, back to margin
        If (INOUT==1.OR.INOUT==0) Then !Inside
            !Save actual position
            PartParam%xpart(iPart) = xt
            PartParam%ypart(iPart) = yt
            PartParam%zpart(iPart) = zt
        Else
            !back to margin
            PartParam%xpart(iPart)=MeshParam%xb(nnel)
            PartParam%ypart(iPart)=MeshParam%yb(nnel)
            
        EndIf

    EndDo

    Return
    End Subroutine ParticlePosition

    Subroutine ResidenceTime(PartParam,HydroParam,MeshParam,dt,Finaltime,IniTime)

    Use MeshVars
    Use Hydrodynamic

    Implicit none

    Integer:: iPart,iTarget,nnel,jlev,INOUT
    Real:: dt
    Real:: xpoly(4),ypoly(4)
    Integer:: Finaltime,IniTime

    type(ParticleParam) :: PartParam
    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam

    Do iPart=1,PartParam%npart

        If (PartParam%FlagTarget(iPart)==1) Then
            Cycle
        EndIf

        Do iTarget = 1,PartParam%ntarget

            nnel = PartParam%ntargetElem(iTarget)

            xpoly(1)= MeshParam%xNode( MeshParam%Quadri(1,nnel)+1 )
            xpoly(2)= MeshParam%xNode( MeshParam%Quadri(2,nnel)+1 )
            xpoly(3)= MeshParam%xNode( MeshParam%Quadri(3,nnel)+1 )
            xpoly(4)= MeshParam%xNode( MeshParam%Quadri(4,nnel)+1 )
            ypoly(1)= MeshParam%yNode( MeshParam%Quadri(1,nnel)+1 )
            ypoly(2)= MeshParam%yNode( MeshParam%Quadri(2,nnel)+1 )
            ypoly(3)= MeshParam%yNode( MeshParam%Quadri(3,nnel)+1 )
            ypoly(4)= MeshParam%yNode( MeshParam%Quadri(4,nnel)+1 )

            Call PNPOLY(PartParam%xpart(iPart),PartParam%ypart(iPart),xpoly,ypoly,4,INOUT)

            If (INOUT==1.OR.INOUT==0) Then !Inside
                PartParam%RTime(iPart) = (Finaltime-IniTime)/86400.
                PartParam%FlagTarget(iPart) = 1
            EndIf

        EndDo
    EndDo




    Return
    End Subroutine ResidenceTime


    Subroutine PNPOLY(PX,PY,XX,YY,N,INOUT)

    ! .................................................................. C
    ! SUBROUTINE PNPOLY made by W. Randolph Franklin (WRF)
    !
    ! PURPOSE
    ! TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON
    !
    ! USAGE
    ! CALL PNPOLY (PX, PY, XX, YY, N, INOUT )
    !
    ! DESCRIPTION OF THE PARAMETERS
    ! PX - X-COORDINATE OF POINT IN QUESTION.
    ! PY - Y-COORDINATE OF POINT IN QUESTION.
    ! XX - N LONG VECTOR CONTAINING X-COORDINATES OF
    ! VERTICES OF POLYGON.
    ! YY - N LONG VECTOR CONTAING Y-COORDINATES OF
    ! VERTICES OF POLYGON.
    ! N - NUMBER OF VERTICES IN THE POLYGON.
    ! INOUT - THE SIGNAL RETURNED:
    ! -1 IF THE POINT IS OUTSIDE OF THE POLYGON,
    ! 0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,
    ! 1 IF THE POINT IS INSIDE OF THE POLYGON.
    !
    ! REMARKS
    ! THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.
    ! THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY
    ! OPTIONALLY BE INCREASED BY 1.
    ! THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING
    ! OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX
    ! OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING
    ! N, THESE FIRST VERTICES MUST BE COUNTED TWICE.
    ! INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.
    ! THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM
    ! WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.
    !
    ! SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
    ! NONE
    !
    ! METHOD
    ! A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT
    ! CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE
    ! POINT IS INSIDE OF THE POLYGON.
    !
    ! ..................................................................
    Real X(200),Y(200),XX(N),YY(N),PX,PY
    LOGICAL MX,MY,NX,NY
    INTEGER O,I,J,N,MAXDIM,INOUT

    ! OUTPUT UNIT FOR PRINTED MESSAGES
    DATA O/6/
    MAXDIM=200
    IF (N.LE.MAXDIM) GO TO 6
    WRITE(O,7)

7   FORMAT('0WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY. 1RESULTS INVALID')
    Return


6   DO 1 I=1,N

        X(I)=XX(I)-PX
1   Y(I)=YY(I)-PY
    INOUT=-1

    DO 2 I=1,N
        J=1+MOD(I,N)
        MX=X(I)>0.0
        NX=X(J)>0.0
        MY=Y(I)>0.0
        NY=Y(J)>0.0
        IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2
        IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3
        INOUT=-INOUT
        GO TO 2

3       IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5
4       INOUT=0


        Return
5       INOUT=-INOUT
2   Continue
    Return
    End Subroutine PNPOLY


    Subroutine quicksearch

    Use MeshVars !, Only:Quadri,Edge,EdgeDef,xNode,yNode,Area,xb,yb,Left,Right,Neighbor,EdgeBary
    Use Hydrodynamic !, Only: ElSmallm,ElCapitalM,Ze,H,Pcri,uxy,uNode
    Implicit none



    Return
    End Subroutine quicksearch

    Subroutine ParticleVelocity(nnel,jlev,staint,zt,upart,vpart,wpart,HydroParam,MeshParam)

    Use MeshVars !, Only:Quadri,Edge,EdgeDef,xNode,yNode,Area,xb,yb,Left,Right,Neighbor,EdgeBary
    Use Hydrodynamic !, Only: ElSmallm,ElCapitalM,Ze,H,Pcri,uxy,uNode
    Implicit none

    Integer:: nnel,jlev,nd,nn,lev,j
    Real:: upart,vpart,wpart
    Real:: vxn(8),vyn(8),vzn(8)
    Real:: staint(8)
    Real:: zup,zrat,zt

    type(HydrodynamicParam) :: HydroParam
    type(MeshGridParam) :: MeshParam

    jlev=1
    upart = 0.
    vpart = 0.
    wpart = 0.
    Do j=1,4
        nd = MeshParam%Edge(j,nnel)
        nn = MeshParam%Quadri(j,nnel)+1
        lev=jlev
        If (lev<HydroParam%ElSmallm(nnel)) Then
            vxn(j) = 0.
            vyn(j) = 0.
        Else
            !x, y and z velocities components at Edges and Nodes of nnel
            If (HydroParam%ElSmallm(nnel)==HydroParam%ElCapitalM(nnel)) Then
                vxn(j) = HydroParam%uxy(lev,1,nd)
                vyn(j) = HydroParam%uxy(lev,2,nd)
            Else
                If (lev==HydroParam%ElSmallm(nnel)) Then
                    If (zt <= (HydroParam%Ze(lev,nnel)+HydroParam%Ze(lev+1,nnel))/2.) Then
                        vxn(j) = HydroParam%uxy(lev,1,nd)
                        vyn(j) = HydroParam%uxy(lev,2,nd)
                    Else
                        zup= (HydroParam%Ze(jlev+2,nnel)+HydroParam%Ze(jlev+1,nnel))/2.
                        zrat=(zup-zt)/(HydroParam%DZi(lev,nnel)/2.+HydroParam%DZi(lev+1,nnel)/2.)
                        vxn(j) = HydroParam%uxy(lev+1,1,nd)*(1-zrat)+HydroParam%uxy(lev,1,nd)*zrat
                        vyn(j) = HydroParam%uxy(lev+1,2,nd)*(1-zrat)+HydroParam%uxy(lev,2,nd)*zrat
                    EndIf
                ElseIf (lev==HydroParam%ElCapitalM(nnel)) Then
                    If (zt <= (HydroParam%Ze(lev,nnel)+HydroParam%Ze(lev+1,nnel))/2.) Then
                        zup= (HydroParam%Ze(lev-1,nnel)+HydroParam%Ze(lev,nnel))/2.
                        zrat=(zt-zup)/(HydroParam%DZi(lev,nnel)/2.+HydroParam%DZi(lev-1,nnel)/2.)
                        vxn(j) = HydroParam%uxy(lev-1,1,nd)*(1-zrat)+HydroParam%uxy(lev,1,nd)*zrat
                        vyn(j) = HydroParam%uxy(lev-1,2,nd)*(1-zrat)+HydroParam%uxy(lev,2,nd)*zrat
                    Else
                        vxn(j) = HydroParam%uxy(lev,1,nd)
                        vyn(j) = HydroParam%uxy(lev,2,nd)
                    EndIf
                Else
                    If (zt <= (HydroParam%Ze(lev,nnel)+HydroParam%Ze(lev+1,nnel))/2.) Then
                        zup= (HydroParam%Ze(lev-1,nnel)+HydroParam%Ze(lev,nnel))/2.
                        zrat=(zt-zup)/(HydroParam%DZi(lev,nnel)/2.+HydroParam%DZi(lev-1,nnel)/2.)
                        vxn(j) = HydroParam%uxy(lev-1,1,nd)*(1-zrat)+HydroParam%uxy(lev,1,nd)*zrat
                        vyn(j) = HydroParam%uxy(lev-1,2,nd)*(1-zrat)+HydroParam%uxy(lev,2,nd)*zrat
                    Else
                        zup= (HydroParam%Ze(lev+2,nnel)+HydroParam%Ze(lev+1,nnel))/2.
                        zrat=(zup-zt)/(HydroParam%DZi(lev,nnel)/2.+HydroParam%DZi(lev+1,nnel)/2.)
                        vxn(j) = HydroParam%uxy(lev+1,1,nd)*(1-zrat)+HydroParam%uxy(lev,1,nd)*zrat
                        vyn(j) = HydroParam%uxy(lev+1,2,nd)*(1-zrat)+HydroParam%uxy(lev,2,nd)*zrat
                    EndIf
                EndIf

            EndIf
        EndIf

        If (MeshParam%Right(nd)==0) Then
            vzn(j) = (HydroParam%w(jlev,MeshParam%Left(nd)))+(zt-HydroParam%Ze(jlev,nnel))*((HydroParam%w(jlev+1,MeshParam%Left(nd)))-(HydroParam%w(jlev,MeshParam%Left(nd))))/HydroParam%DZi(jlev,nnel)
        Else
            vzn(j) = 0.5*(HydroParam%w(jlev,MeshParam%Right(nd))+HydroParam%w(jlev,MeshParam%Left(nd)))+(zt-HydroParam%Ze(jlev,nnel))*(0.5*(HydroParam%w(jlev+1,MeshParam%Right(nd))+HydroParam%w(jlev+1,MeshParam%Left(nd)))-0.5*(HydroParam%w(jlev,MeshParam%Right(nd))+HydroParam%w(jlev,MeshParam%Left(nd))))/HydroParam%DZi(jlev,nnel)
        EndIf

        !Update x, y and z velocities components of the particle
        upart=upart+vxn(j)*staint(j)
        vpart=vpart+vyn(j)*staint(j)
        wpart=wpart+vzn(j)*staint(j)
    EndDo !j

    Return
    End Subroutine ParticleVelocity


    Subroutine ibilinear(elem,x1,x2,x3,x4,y1,y2,y3,y4,x,y,xi,eta,shapef)

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



    End Module ParticleTracking