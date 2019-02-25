Subroutine UpWindCWC(index,uLoadVar,dVar,VarEst,VarEstP,dt,dtday,HydroParam,MeshParam)

    ! This routine solves the Transport Equation usind a centered finite difference scheme for variables with Boundary Condition allowed.
    ! Called in routines: ABIOTICO, SALINITY

    ! Example:
    ! sDSal(I,J,K)       == VarEst
    ! dVarSal            == dVar
    ! sDSalP(I,J,K)      == VarEstP
    ! sDSalP(I-1,J,K)    == VarEstPIM
    ! sDSalP(I+1,J,K)    == VarEstPIP	
    ! sDSalP(I,J-1,K)    == VarEstPJM
    ! sDSalP(I,J+1,K)    == VarEstPJP
    ! sDSalP(I,J,K-1)    == VarEstPKM
    ! sDSalP(I,J,K+1)    == VarEstPKP
    ! INDEXCONT(ICEL,14) == index => 14
    ! uLoadVar        == uDLoadSal
    ! sDSal0             == VarEstINI
    
    Use MeshVars
    Use Hydrodynamic

    Implicit None
    Real, Dimension(4):: ConcFace
	Real:: CONCU,CONCD,DDCFace,DCDFace,Dispersion,Advection
	Real:: DDCZ,DC2DZ2,DCDZU,DCDZD
    Real:: V
	Integer:: index,iElem,iEdge,Face,iLayer,Sig,l,r
    Real:: NearZero = 1e-10
    Real:: dt,dtday
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Real, Dimension(MeshParam%KMax,MeshParam%nElem):: dVar,VarEst,VarEstP
    Real, Dimension(MeshParam%KMax,MeshParam%nElem):: uLoadVar
   
    Do iElem = 1,MeshParam%nElem
        !If (iElem==11802) Then
        !    Continue
        !EndIf
        !

        If (HydroParam%eta(iElem)-HydroParam%hb(iElem) <= HydroParam%PCRI+NearZero.or.HydroParam%etan(iElem)-HydroParam%hb(iElem)<= HydroParam%PCRI+NearZero) Then
            VarEst(HydroParam%ElCapitalM(iElem),iElem) = VarEstP(HydroParam%ElCapitalM(iElem),iElem)  
            Cycle
        EndIf

            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)

                    ! Horizontal Eddy-Diffusivity (Dispersion coefficient)
                    ! Concentration in each face
                    DCDFace = 0.
                    Do iEdge = 1,4
                        Face = MeshParam%Edge(iEdge,iElem)
                        l = MeshParam%Left(Face)
                        r = MeshParam%Right(Face)
                        If ( r==0 .or.HydroParam%H(MeshParam%Edge(iEdge,iElem)) <= HydroParam%PCRI+NearZero) Then 
		                    ConcFace(iEdge) = VarEstP(iLayer,iElem)
                            DDCFace = 0.
		                    DCDFace = DCDFace 
                        Else
                            If (Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)*HydroParam%DZjt(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)*HydroParam%DZjt(iLayer,Face))>=0) Then
                                ConcFace(iEdge) = VarEstP(iLayer,iElem)
                            Else
                                If (iLayer<HydroParam%ElSmallm(MeshParam%Neighbor(iEdge,iElem))) Then 
                                    ConcFace(iEdge) = VarEstP(iLayer,iElem)
                                Else
                                    ConcFace(iEdge) = VarEstP(iLayer,MeshParam%Neighbor(iEdge,iElem))
                                EndIf
                            EndIf
                            DDCFace = (0.5*(HydroParam%HorDiffusivity(2,iLayer,MeshParam%Neighbor(iEdge,iElem))+HydroParam%HorDiffusivity(2,iLayer,iElem)))
		                    DCDFace = DCDFace + MeshParam%EdgeLength(Face)/MeshParam%CirDistance(Face)*HydroParam%DZjt(iLayer,Face)*DDCFace*(ConcFace(iEdge)-VarEstP(iLayer,iElem)) 
                        EndIf
                    EndDo
                    Dispersion = DCDFace/MeshParam%Area(iElem)
                
                If ( HydroParam%ElSmallm(iElem) == HydroParam%ElCapitalM(iElem) ) Then        ! Only One Vertical Layer
                    
                    If (abs(-Sig(iElem,MeshParam%Right(MeshParam%Edge(1,iElem)),MeshParam%Left(MeshParam%Edge(1,iElem)))*(HydroParam%Theta*HydroParam%u(iLayer,MeshParam%Edge(1,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(1,iElem))+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,MeshParam%Edge(1,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(1,iElem)))*MeshParam%dx*DT-Sig(iElem,MeshParam%Right(MeshParam%Edge(2,iElem)),MeshParam%Left(MeshParam%Edge(2,iElem)))*(HydroParam%Theta*HydroParam%u(iLayer,MeshParam%Edge(2,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(2,iElem))+ (1.-HydroParam%Theta)*HydroParam%ut(iLayer,MeshParam%Edge(2,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(2,iElem)))*MeshParam%dy*dt-Sig(iElem,MeshParam%Right(MeshParam%Edge(3,iElem)),MeshParam%Left(MeshParam%Edge(3,iElem)))*(HydroParam%Theta*HydroParam%u(iLayer,MeshParam%Edge(3,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(3,iElem))+ (1.-HydroParam%Theta)*HydroParam%ut(iLayer,MeshParam%Edge(3,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(3,iElem)))*MeshParam%dy*dt-Sig(iElem,MeshParam%Right(MeshParam%Edge(4,iElem)),MeshParam%Left(MeshParam%Edge(4,iElem)))*(HydroParam%Theta*HydroParam%u(iLayer,MeshParam%Edge(4,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(4,iElem))+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,MeshParam%Edge(4,iElem))*HydroParam%DZjt(iLayer,MeshParam%Edge(4,iElem)))*MeshParam%dy*dt- (HydroParam%eta(iElem)-HydroParam%etaplus(ielem)-HydroParam%etan(iElem))*MeshParam%dx*MeshParam%dy)>0.1) Then
                        VarEst(HydroParam%ElCapitalM(iElem),iElem) = VarEstP(HydroParam%ElCapitalM(iElem),iElem)  
                        Cycle
                    EndIf
                        
                    ! Top
                    CONCU = VarEstP(iLayer,iElem)
                    ! Bottom
		            CONCD = VarEstP(iLayer,iElem)
                    
                    DC2DZ2 = 0.
        
                    ! Turbulence Model (pegar referência)  
                    DDCZ = 0.
        
                    ! Transport Equation Solution for 2D
                    If (index== 0) Then !No boundary condition
                        Advection = 0.
                        Do iEdge = 1,4
                            Face = MeshParam%Edge(iEdge,iElem)
                            l = MeshParam%Left(Face)
                            r = MeshParam%Right(Face)
                            Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                        EndDo
                        VarEst(iLayer,iElem) = MAX(NearZero,(DT*(Dispersion-Advection)+dtday*dVar(iLayer,iElem)*V(HydroParam%etan(iElem),HydroParam%hb(iElem))+VarEstP(iLayer,iElem)*V(HydroParam%etan(iElem),HydroParam%hb(iElem)))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))) 
                    Else
                        Advection = 0.
                        Do iEdge = 1,4
                            Face = MeshParam%Edge(iEdge,iElem)
                            l = MeshParam%Left(Face)
                            r = MeshParam%Right(Face)
                            If (HydroParam%IndexInflowEdge(Face) > 0.or.HydroParam%IndexWaterLevelEdge(Face) > 0) Then !Inflow/Outflow boundary condition
                                If (HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)>=0) Then ! Outflow Loading
                                    Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)    
                                Else !InFlow Loading
                                    Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*uLoadVar(iLayer,iElem)/MeshParam%Area(iElem) 
                                EndIf
                            Else
                                Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                            EndIf
                        EndDo
                        VarEst(iLayer,iElem) = MAX(NearZero,(DT*(Dispersion-Advection)+dtday*dVar(iLayer,iElem)*V(HydroParam%etan(iElem),HydroParam%hb(iElem))+VarEstP(iLayer,iElem)*V(HydroParam%etan(iElem),HydroParam%hb(iElem)))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))) 
                    EndIf
                    
                
                Else
                    If ( iLayer == HydroParam%ElSmallm(iElem) ) Then !Bottom layer

                        ! Concentration in the Top layer
                        If (HydroParam%w(iLayer+1,iElem)>=0) Then
                            CONCU = VarEstP(iLayer,iElem)                         ! Top face concentration is equal to the concentration in the layer    
                        Else
                            CONCU = VarEstP(iLayer+1,iElem)                         ! Top face concentration is equal to the concentration in the layer
                        EndIF
                        ! Concentration in the Top Bottom layer
		                CONCD = VarEstP(iLayer,iElem) 

		                ! Vertical Diffusivity Flux in the top layer
                        DDCZ = HydroParam%VerEddyDiffCell(iLayer+1,iElem)

                        DCDZU = DDCZ*( VarEstP(iLayer+1,iElem) - VarEstP(iLayer,iElem) )/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer+1,iElem))/2.)
                        
                        !DCDZU = Max(0., DCDZU - 0.5*abs(MeshParam%Area(iElem)*HydroParam%w(iLayer+1,iElem)))

                        ! Vertical Diffusivity Flux in the bottom layer
                        DCDZD = 0.
                        ! Vertical gradient of Diffusivity
                        DC2DZ2 = DCDZU - DCDZD 
        
        
                        ! Transport Equation Solution
                        If (index== 0) Then !No boundary condition
                            Advection = 0.
                            Do iEdge = 1,4
                                Face = MeshParam%Edge(iEdge,iElem)
                                l = MeshParam%Left(Face)
                                r = MeshParam%Right(Face)
                                Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                            EndDo
                            VarEst(iLayer,iElem) = MAX(NearZero,(DT*(Dispersion-Advection+DC2DZ2)+dtday*dVar(iLayer,iElem)*HydroParam%DZit(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZit(iLayer,iElem)-((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*DT)/HydroParam%DZit(iLayer,iElem))
                        Else
                            Advection = 0.
                            Do iEdge = 1,4
                                Face = MeshParam%Edge(iEdge,iElem)
                                l = MeshParam%Left(Face)
                                r = MeshParam%Right(Face)
                                If (HydroParam%IndexInflowEdge(Face) > 0.or.HydroParam%IndexWaterLevelEdge(Face) > 0) Then !Inflow/Outflow boundary condition
                                    If (HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)>=0) Then ! Outflow Loading
                                        Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)    
                                    Else !InFlow Loading
                                        Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*uLoadVar(iLayer,iElem)/MeshParam%Area(iElem) 
                                    EndIf
                                Else
                                    Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                                EndIf
                            EndDo
                            VarEst(iLayer,iElem) = MAX(NearZero,(DT*(Dispersion-Advection+DC2DZ2)+dtday*dVar(iLayer,iElem)*HydroParam%DZit(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZit(iLayer,iElem)-((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*DT)/HydroParam%DZit(iLayer,iElem))

                        EndIf
                    
                    Elseif ( iLayer == HydroParam%ElCapitalM(iElem) ) Then !Surface layer  
                        
                        ! Concentration in the Top layer
                        CONCU = VarEstP(iLayer,iElem)                      ! Up
                        ! Concentration in the bottom layer
                        If (HydroParam%w(iLayer,iElem)>=0) Then
		                    CONCD = VarEstP(iLayer-1,iElem) 
                        Else
                            CONCD = VarEstP(iLayer,iElem)
                        EndIf

		                ! Vertical Diffusivity Flux in the top layer
                        DCDZU = 0.
                        ! Vertical Diffusivity Flux in the bottom layer
                        DDCZ = HydroParam%VerEddyDiffCell(iLayer,iElem)
                        DCDZD = DDCZ*( VarEstP(iLayer,iElem) - VarEstP(iLayer-1,iElem) )/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))/2.)
                        
                        !DCDZD = Max(0., DCDZD - 0.5*abs(MeshParam%Area(iElem)*HydroParam%w(iLayer,iElem)))
                        ! Vertical gradient of Diffusivity
                        DC2DZ2 = DCDZU - DCDZD 
                        
                        ! Transport Equation Solution
                        If (index== 0) Then !No boundary condition
                            Advection = 0.
                            Do iEdge = 1,4
                                Face = MeshParam%Edge(iEdge,iElem)
                                l = MeshParam%Left(Face)
                                r = MeshParam%Right(Face)
                                Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                            EndDo
                            VarEst(iLayer,iElem) = MAX(NearZero,(DT*(Dispersion-Advection+DC2DZ2)+dtday*dVar(iLayer,iElem)*HydroParam%DZit(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZit(iLayer,iElem) - ((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*DT)/HydroParam%DZit(iLayer,iElem))
                        Else
                            Advection = 0.
                            Do iEdge = 1,4
                                Face = MeshParam%Edge(iEdge,iElem)
                                l = MeshParam%Left(Face)
                                r = MeshParam%Right(Face)
                                If (HydroParam%IndexInflowEdge(Face) > 0.or.HydroParam%IndexWaterLevelEdge(Face) > 0) Then !Inflow/Outflow boundary condition
                                    If (HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)>=0) Then ! Outflow Loading
                                        Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)    
                                    Else !InFlow Loading
                                        Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*uLoadVar(iLayer,iElem)/MeshParam%Area(iElem) 
                                    EndIf
                                Else
                                    Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                                EndIf
                            EndDo
                            VarEst(iLayer,iElem) = MAX(NearZero,(DT*(Dispersion-Advection+DC2DZ2)+dtday*dVar(iLayer,iElem)*HydroParam%DZit(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZit(iLayer,iElem) - ((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*DT)/HydroParam%DZit(iLayer,iElem))
                        EndIf

                    Else !Intermediary layer

                        ! Calcula a Concentração na camada
                        ! Top
                        If (HydroParam%w(iLayer+1,iElem)>=0) Then
                            CONCU = VarEstP(iLayer,iElem)                         ! Top face concentration is equal to the concentration in the layer    
                        Else
                            CONCU = VarEstP(iLayer+1,iElem)                         ! Top face concentration is equal to the concentration in the layer
                        EndIF
                        
                        ! Bottom
                        If (HydroParam%w(iLayer,iElem)>=0) Then
		                    CONCD = VarEstP(iLayer-1,iElem) 
                        Else
                            CONCD = VarEstP(iLayer,iElem)
                        EndIf
        
		                ! Vertical Diffusivity Flux in the top layer
                        DDCZ = HydroParam%VerEddyDiffCell(iLayer+1,iElem)
                        DCDZU = DDCZ*( VarEstP(iLayer+1,iElem) - VarEstP(iLayer,iElem) )/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer+1,iElem))/2.)
                        !DCDZU = Max(0., DCDZU - 0.5*abs(MeshParam%Area(iElem)*HydroParam%w(iLayer+1,iElem)))
                        ! Vertical Diffusivity Flux in the bottom layer
                        DDCZ = HydroParam%VerEddyDiffCell(iLayer,iElem)
                        DCDZD = DDCZ*( VarEstP(iLayer,iElem) - VarEstP(iLayer-1,iElem) )/((HydroParam%DZit(iLayer,iElem)+HydroParam%DZit(iLayer-1,iElem))/2.)
                        !DCDZD = Max(0., DCDZD - 0.5*abs(MeshParam%Area(iElem)*HydroParam%w(iLayer,iElem)))
                        ! Vertical gradient of Diffusivity
                        DC2DZ2 = DCDZU - DCDZD 
                        
                        
                        ! Transport Equation Solution
                        If (index== 0) Then !No boundary condition
                            Advection = 0.
                            Do iEdge = 1,4
                                Face = MeshParam%Edge(iEdge,iElem)
                                l = MeshParam%Left(Face)
                                r = MeshParam%Right(Face)
                                Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                            EndDo
                            VarEst(iLayer,iElem) = MAX(NearZero,(DT*(Dispersion-Advection+DC2DZ2)+dtday*dVar(iLayer,iElem)*HydroParam%DZit(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZit(iLayer,iElem) - ((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*DT)/HydroParam%DZit(iLayer,iElem))
                        Else
                            Advection = 0.
                            Do iEdge = 1,4
                                Face = MeshParam%Edge(iEdge,iElem)
                                l = MeshParam%Left(Face)
                                r = MeshParam%Right(Face)
                                If (HydroParam%IndexInflowEdge(Face) > 0.or.HydroParam%IndexWaterLevelEdge(Face) > 0) Then !Inflow/Outflow boundary condition
                                    If (HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)>=0) Then ! Outflow Loading
                                        Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)    
                                    Else !InFlow Loading
                                        Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*uLoadVar(iLayer,iElem)/MeshParam%Area(iElem) 
                                    EndIf
                                Else
                                    Advection = Advection + (MeshParam%EdgeLength(Face)*HydroParam%DZjt(iLayer,Face)*Sig(iElem,r,l)*(HydroParam%Theta*HydroParam%u(iLayer,Face)+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face)))*ConcFace(iEdge)/MeshParam%Area(iElem)
                                EndIf
                            EndDo
                            VarEst(iLayer,iElem) = MAX(NearZero,(DT*(Dispersion-Advection+DC2DZ2)+dtday*dVar(iLayer,iElem)*HydroParam%DZit(iLayer,iElem)+VarEstP(iLayer,iElem)*HydroParam%DZit(iLayer,iElem) - ((HydroParam%Theta*HydroParam%w(iLayer+1,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer+1,iElem))*CONCU-(HydroParam%Theta*HydroParam%w(iLayer,iElem)+(1.0d0-HydroParam%Theta)*HydroParam%wt(iLayer,iElem))*CONCD)*DT)/HydroParam%DZit(iLayer,iElem))
                        EndIf

                    EndIf
                EndIf
            EndDo
            
            ! 3.1 Copy concentration Above the Free-Surface
            Do iLayer = HydroParam%ElCapitalM(iElem) + 1, MeshParam%KMAX
                VarEst(iLayer,iElem) = VarEst(HydroParam%ElCapitalM(iElem),iElem)
            EndDo
            ! 3.2 Nullify Velocities Below the Bottom (u=0) 
            Do iLayer = 1, HydroParam%ElSmallm(iElem) - 1
                VarEst(iLayer,iElem) = VarEst(HydroParam%ElSmallm(iElem),iElem)
            EndDo
            
    EndDo
    ! preencher valores fora dos limites da camada

    Return
End Subroutine UpWindCWC