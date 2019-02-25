Subroutine wvelocity(HydroParam,MeshParam,dt)
    
    !$ use omp_lib
    Use MeshVars
    Use Hydrodynamic
    
    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Integer:: iElem,iLayer, i
    Real:: aTh(MeshParam%KMax-1), bTh(MeshParam%KMax-1), cTh(MeshParam%KMax-1), fTh(MeshParam%KMax-1)
    Real :: bp(MeshParam%KMax-1),vp(MeshParam%KMax-1),a_aux(MeshParam%KMax-1),b_aux(MeshParam%KMax),c_aux(MeshParam%KMax-1),v_aux(MeshParam%KMax-1),x_aux(MeshParam%KMax-1),d_aux(MeshParam%KMax-1),x(MeshParam%KMax-1)
    Real :: w
    Real:: dt,dzk
    
    Do iElem = 1, MeshParam%nElem
        If (HydroParam%ElSmallm(iElem)==HydroParam%ElCapitalM(iElem)) Then
            HydroParam%w(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem),iElem) = 0.
        Else
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)-1       
                If (iLayer == HydroParam%ElSmallm(iElem)) Then
                    dzk         = 0.5*(HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer+1,iElem)) 
                    aTh(iLayer) = 0.
                    bTh(iLayer) = 1. + dt/dzk*(HydroParam%VerEddyViscCell(iLayer,iElem)/HydroParam%DZit(iLayer,iElem) + HydroParam%VerEddyViscCell(iLayer+1,iElem)/HydroParam%DZit(iLayer+1,iElem))
                    cTh(iLayer) = -dt/dzk*(HydroParam%VerEddyViscCell(iLayer+1,iElem)/HydroParam%DZit(iLayer+1,iElem))
                    fTh(iLayer) = HydroParam%Fw(iLayer+1,iElem)- (dt/dzk*(1.-HydroParam%theta)*( HydroParam%q(iLayer+1,iElem)-HydroParam%q(iLayer,iElem) ))
                Elseif (iLayer == HydroParam%ElCapitalM(iElem)-1) Then
                    dzk         = 0.5*(HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer+1,iElem)) 
                    aTh(iLayer) = -dt/dzk*(HydroParam%VerEddyViscCell(iLayer,iElem)/HydroParam%DZit(iLayer,iElem))
                    bTh(iLayer) = 1. + dt/dzk*(HydroParam%VerEddyViscCell(iLayer,iElem)/HydroParam%DZit(iLayer,iElem) + HydroParam%VerEddyViscCell(iLayer+1,iElem)/HydroParam%DZit(iLayer+1,iElem))
                    cTh(iLayer) = 0.
                    fTh(iLayer) = HydroParam%Fw(iLayer+1,iElem)- (dt/dzk*(1.-HydroParam%theta)*( HydroParam%q(iLayer+1,iElem)-HydroParam%q(iLayer,iElem) ))
                Else
                    dzk         = 0.5*(HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer+1,iElem)) 
                    aTh(iLayer) = -dt/dzk*(HydroParam%VerEddyViscCell(iLayer,iElem)/HydroParam%DZit(iLayer,iElem))
                    bTh(iLayer) = 1. + dt/dzk*(HydroParam%VerEddyViscCell(iLayer,iElem)/HydroParam%DZit(iLayer,iElem) + HydroParam%VerEddyViscCell(iLayer+1,iElem)/HydroParam%DZit(iLayer+1,iElem))
                    cTh(iLayer) = -dt/dzk*(HydroParam%VerEddyViscCell(iLayer+1,iElem)/HydroParam%DZit(iLayer+1,iElem))
                    fTh(iLayer) = HydroParam%Fw(iLayer+1,iElem)- (dt/dzk*(1.-HydroParam%theta)*( HydroParam%q(iLayer+1,iElem)-HydroParam%q(iLayer,iElem) ))
                EndIf
            EndDo
            !DIM = size(HydroParam%w(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1,iElem),1)
            !Call solve_tridiag(cTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1),bTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1),aTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1),fTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1),HydroParam%w(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem),DIM)
            
            !Tridiagonal Matrix solver
            !Flipping matrix
            do i = HydroParam%ElCapitalM(iElem)-1,HydroParam%ElSmallm(iElem),-1
                a_aux(HydroParam%ElCapitalM(iElem)-1-i+HydroParam%ElSmallm(iElem)) = cTh(i)
                b_aux(HydroParam%ElCapitalM(iElem)-1-i+HydroParam%ElSmallm(iElem)) = bTh(i)
                c_aux(HydroParam%ElCapitalM(iElem)-1-i+HydroParam%ElSmallm(iElem)) = aTh(i)
                v_aux(HydroParam%ElCapitalM(iElem)-1-i+HydroParam%ElSmallm(iElem)) = fTh(i)
            enddo
 
            d_aux = 0   
            x = d_aux
            w = b_aux(HydroParam%ElSmallm(iElem))
            x(HydroParam%ElSmallm(iElem)) = v_aux(HydroParam%ElSmallm(iElem))/w
            do i=HydroParam%ElSmallm(iElem)+1,HydroParam%ElCapitalM(iElem)-1
                d_aux(i-1) = c_aux(i-1)/w;
                w = b_aux(i) - a_aux(i)*d_aux(i-1);
                x(i) = ( v_aux(i) - a_aux(i)*x(i-1) )/w;
            enddo
            If (HydroParam%ElCapitalM(iElem)>2) Then
                do i=HydroParam%ElCapitalM(iElem)-1-1, HydroParam%ElSmallm(iElem), -1
                    x(i) = x(i) - d_aux(i)*x(i+1);
                enddo
            EndIf

            !Flipping matrix
            do i = HydroParam%ElCapitalM(iElem)-1,HydroParam%ElSmallm(iElem),-1
                x_aux(HydroParam%ElCapitalM(iElem)-1-i+HydroParam%ElSmallm(iElem)) = x(i)
            enddo
            HydroParam%w(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem) = x_aux(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1)
            
        EndIf
    EndDo
  
    Return    
End Subroutine wvelocity