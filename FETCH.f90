SUBROUTINE Fetch(MeshParam,HydroParam)

    ! This routine calculates the "Lake Fetch" for every computational cell in eigth directions:
    ! North (N), Northeast (NE), East (E), Southeast (SE), South (S), Southwest (SW), West (W), Northwest (NW)
    ! Called in routine 0-MAIN
    
	Use MeshVars
    Use Hydrodynamic
    
	IMPLICIT NONE
	
	Real:: dDIAG
    INTEGER:: contador, N
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    
	
    ! Setting initial conditions for auxiliary variables
    MeshParam%trv     = 0
    MeshParam%NEIBN   = 0
    MeshParam%NEIBL   = 0
    MeshParam%NEIBS   = 0
    MeshParam%NEIBO   = 0
    HydroParam%fetch_m = 0.
	dDIAG   = SQRT(MeshParam%DX*MeshParam%DX+MeshParam%DY*MeshParam%DY)

    ! Winds in the Northward Direction 
    MeshParam%trv = 1        ! Save the kind of cell
    WHERE(MeshParam%trv < 0)
        MeshParam%trv = 0         ! This cell is out of the interest region
    END WHERE
    MeshParam%NEIBN    = MeshParam%Neighbor(1,:)        ! North Neighbour
    MeshParam%NEIBL    = MeshParam%Neighbor(4,:)        ! East Neighbour
    MeshParam%NEIBS    = MeshParam%Neighbor(3,:)        ! South Neighbour
    MeshParam%NEIBO    = MeshParam%Neighbor(2,:)        ! West Neighbour
    contador = 1
    DO WHILE(contador > 0)
        WHERE(MeshParam%trv > 0)
            HydroParam%fetch_m(:,1) = HydroParam%fetch_m(:,1) + MeshParam%DY
        END WHERE
        DO N = 1,MeshParam%nElem
            IF ( MeshParam%NEIBN(N) == 0.or.MeshParam%NEIBN(N) == -1 ) THEN
                MeshParam%trv(N)     = 0
            ENDIF
        ENDDO
        contador = 0
        DO N=1,MeshParam%nElem
            IF ( MeshParam%trv(N) == 0 ) THEN
                If (MeshParam%NEIBS(N)/=0) MeshParam%NEIBN(MeshParam%NEIBS(N))=0
                    
                
            ENDIF
            IF (MeshParam%trv(N) > 0) THEN
                contador = contador + 1
            ENDIF
        ENDDO
    END DO

! --------------------------------------------------------------
    
    ! Winds in the Northeastward Direction 
    MeshParam%trv      = 1
    MeshParam%NEIBN    = MeshParam%Neighbor(1,:)        ! North Neighbour
    MeshParam%NEIBL    = MeshParam%Neighbor(4,:)        ! East Neighbour
    MeshParam%NEIBS    = MeshParam%Neighbor(3,:)        ! South Neighbour
    MeshParam%NEIBO    = MeshParam%Neighbor(2,:)        ! West Neighbour
    contador = 1
    DO WHILE(contador > 0)
        WHERE(MeshParam%trv > 0)
            HydroParam%fetch_m(:,2) = HydroParam%fetch_m(:,2) + dDIAG
        END WHERE
        DO N=1,MeshParam%nElem
		    IF ( MeshParam%NEIBN(N)==0.OR.MeshParam%NEIBL(N)==0.or.MeshParam%NEIBN(N)==-1.OR.MeshParam%NEIBL(N)==-1 ) THEN
		        MeshParam%trv(N)     = 0
		    ENDIF
        ENDDO
        contador = 0
        DO N=1,MeshParam%nElem
		    IF (MeshParam%trv(N) == 0) THEN
                If (MeshParam%NEIBS(N)/=0) MeshParam%NEIBN(MeshParam%NEIBS(N))=0
                If (MeshParam%NEIBO(N)/=0) MeshParam%NEIBL(MeshParam%NEIBO(N))=0
		    ENDIF
		    IF (MeshParam%trv(N) > 0) THEN
		        contador = contador + 1
		    ENDIF
        ENDDO
    END DO

! --------------------------------------------------------------
    
    ! Winds in the Eastward Direction 
    MeshParam%trv      = 1
    MeshParam%NEIBN    = MeshParam%Neighbor(1,:)        ! North Neighbour
    MeshParam%NEIBL    = MeshParam%Neighbor(4,:)        ! East Neighbour
    MeshParam%NEIBS    = MeshParam%Neighbor(3,:)        ! South Neighbour
    MeshParam%NEIBO    = MeshParam%Neighbor(2,:)        ! West Neighbour
    contador = 1
    DO WHILE(contador > 0)
        WHERE(MeshParam%trv > 0)
            HydroParam%fetch_m(:,3) = HydroParam%fetch_m(:,3) + MeshParam%DX
        END WHERE
        DO N=1,MeshParam%nElem
		    IF ( MeshParam%NEIBL(N) ==0.or.MeshParam%NEIBL(N) ==-1 ) THEN
		        MeshParam%trv(N)     = 0
		    ENDIF
        ENDDO
        contador = 0
        DO N=1,MeshParam%nElem
		    IF (MeshParam%trv(N) == 0) THEN
                If (MeshParam%NEIBO(N)/=0) MeshParam%NEIBL(MeshParam%NEIBO(N))=0
		    ENDIF
		    IF (MeshParam%trv(N) > 0) THEN
		        contador = contador + 1
		    ENDIF
        ENDDO
    END DO
! --------------------------------------------------------------
    
    ! Winds in the Southeastward Direction 
    MeshParam%trv      = 1
    MeshParam%NEIBN    = MeshParam%Neighbor(1,:)        ! North Neighbour
    MeshParam%NEIBL    = MeshParam%Neighbor(4,:)        ! East Neighbour
    MeshParam%NEIBS    = MeshParam%Neighbor(3,:)        ! South Neighbour
    MeshParam%NEIBO    = MeshParam%Neighbor(2,:)        ! West Neighbour
    contador = 1
    DO WHILE(contador > 0)
        WHERE(MeshParam%trv > 0)
            HydroParam%fetch_m(:,4) = HydroParam%fetch_m(:,4) + dDIAG
        END WHERE
        DO N=1,MeshParam%nElem
		    IF ( MeshParam%NEIBS(N)==0.OR.MeshParam%NEIBL(N)==0.or.MeshParam%NEIBS(N)==-1.OR.MeshParam%NEIBL(N)==-1 ) THEN
		        MeshParam%trv(N)     = 0
		    ENDIF
        ENDDO
        contador = 0
        DO N=1,MeshParam%nElem
		    IF ( MeshParam%trv(N)==0 ) THEN
                If (MeshParam%NEIBN(N)/=0) MeshParam%NEIBS(MeshParam%NEIBN(N)) = 0
                If (MeshParam%NEIBO(N)/=0) MeshParam%NEIBL(MeshParam%NEIBO(N)) = 0
		    ENDIF
		    IF (MeshParam%trv(N) > 0) THEN
		        contador = contador + 1
		    ENDIF
        ENDDO
    ENDDO
! --------------------------------------------------------------
    
    ! Winds in the Southward Direction 
    MeshParam%trv      = 1
    MeshParam%NEIBN    = MeshParam%Neighbor(1,:)        ! North Neighbour
    MeshParam%NEIBL    = MeshParam%Neighbor(4,:)        ! East Neighbour
    MeshParam%NEIBS    = MeshParam%Neighbor(3,:)        ! South Neighbour
    MeshParam%NEIBO    = MeshParam%Neighbor(2,:)        ! West Neighbour
    contador = 1
    DO WHILE(contador > 0)
        WHERE(MeshParam%trv > 0)
            HydroParam%fetch_m(:,5) = HydroParam%fetch_m(:,5) + MeshParam%DY
        END WHERE
        DO N=1,MeshParam%nElem
		    IF ( MeshParam%NEIBS(N) == 0.or.MeshParam%NEIBS(N) == -1 ) THEN
		        MeshParam%trv(N)     = 0
		    ENDIF
        ENDDO
        contador = 0
        DO N=1,MeshParam%nElem
		    IF ( MeshParam%trv(N) == 0 ) THEN
                If (MeshParam%NEIBN(N)/=0) MeshParam%NEIBS(MeshParam%NEIBN(N))=0
		    ENDIF
		    IF ( MeshParam%trv(N) > 0 ) THEN
		        contador = contador + 1
		    ENDIF
        ENDDO
    ENDDO

! --------------------------------------------------------------
    
    ! Winds in the Southwestward Direction
    MeshParam%trv      = 1
    MeshParam%NEIBN    = MeshParam%Neighbor(1,:)        ! North Neighbour
    MeshParam%NEIBL    = MeshParam%Neighbor(4,:)        ! East Neighbour
    MeshParam%NEIBS    = MeshParam%Neighbor(3,:)        ! South Neighbour
    MeshParam%NEIBO    = MeshParam%Neighbor(2,:)        ! West Neighbour
    contador = 1
    DO WHILE(contador > 0)
        WHERE(MeshParam%trv > 0)
            HydroParam%fetch_m(:,6) = HydroParam%fetch_m(:,6) + dDIAG
        END WHERE
        DO N=1,MeshParam%nElem
		    IF ( MeshParam%NEIBS(N)==0.OR.MeshParam%NEIBO(N)==0.or.MeshParam%NEIBS(N)==-1.OR.MeshParam%NEIBO(N)==-1 ) THEN
		        MeshParam%trv(N)     = 0
		    ENDIF
        ENDDO
        contador = 0
        DO N=1,MeshParam%nElem
		    IF ( MeshParam%trv(N) == 0 ) THEN
                If (MeshParam%NEIBN(N)/=0) MeshParam%NEIBS(MeshParam%NEIBN(N))=0
                If (MeshParam%NEIBL(N)/=0) MeshParam%NEIBO(MeshParam%NEIBL(N))=0
		    ENDIF
		    IF ( MeshParam%trv(N) > 0 ) THEN
		        contador = contador + 1
		    ENDIF
        ENDDO
    ENDDO

! --------------------------------------------------------------
    
    ! Winds in the Westward Direction
    MeshParam%trv      = 1
    MeshParam%NEIBN    = MeshParam%Neighbor(1,:)        ! North Neighbour
    MeshParam%NEIBL    = MeshParam%Neighbor(4,:)        ! East Neighbour
    MeshParam%NEIBS    = MeshParam%Neighbor(3,:)        ! South Neighbour
    MeshParam%NEIBO    = MeshParam%Neighbor(2,:)        ! West Neighbour
    contador = 1
    DO WHILE(contador > 0)
        WHERE(MeshParam%trv > 0)
            HydroParam%fetch_m(:,7) = HydroParam%fetch_m(:,7) + MeshParam%DX
        END WHERE
        DO N=1,MeshParam%nElem
		    IF ( MeshParam%NEIBO(N) == 0.or.MeshParam%NEIBO(N) == -1 ) THEN
		        MeshParam%trv(N)     = 0
		    ENDIF
        ENDDO
        contador = 0
        DO N=1,MeshParam%nElem
		    IF ( MeshParam%trv(N)==0 ) THEN
                If (MeshParam%NEIBL(N)/=0) MeshParam%NEIBO(MeshParam%NEIBL(N))=0
		    ENDIF
		    IF ( MeshParam%trv(N) > 0 ) THEN
		        contador = contador + 1
		    ENDIF
        ENDDO
    ENDDO
! --------------------------------------------------------------
    
    ! Winds in the Northwestward Direction
    MeshParam%trv      = 1
    MeshParam%NEIBN    = MeshParam%Neighbor(1,:)        ! North Neighbour
    MeshParam%NEIBL    = MeshParam%Neighbor(4,:)        ! East Neighbour
    MeshParam%NEIBS    = MeshParam%Neighbor(3,:)        ! South Neighbour
    MeshParam%NEIBO    = MeshParam%Neighbor(2,:)        ! West Neighbour
    contador = 1
    DO WHILE(contador > 0)
        WHERE(MeshParam%trv > 0)
            HydroParam%fetch_m(:,8) = HydroParam%fetch_m(:,8) + dDIAG
        END WHERE
        DO N=1,MeshParam%nElem
		    IF ( MeshParam%NEIBN(N)==0.OR.MeshParam%NEIBO(N)==0.or.MeshParam%NEIBN(N)==-1.OR.MeshParam%NEIBO(N)==-1 ) THEN
		        MeshParam%trv(N)     = 0
		    ENDIF
        ENDDO
        contador = 0
        DO N=1,MeshParam%nElem
		    IF ( MeshParam%trv(N) == 0) THEN
                If (MeshParam%NEIBS(N)/=0) MeshParam%NEIBN(MeshParam%NEIBS(N)) = 0
                If (MeshParam%NEIBL(N)/=0) MeshParam%NEIBO(MeshParam%NEIBL(N)) = 0
		    ENDIF
		    IF ( MeshParam%trv(N) > 0 ) THEN
		        contador = contador + 1
		    ENDIF
        ENDDO
    ENDDO
! --------------------------------------------------------------
    
	RETURN
END SUBROUTINE Fetch