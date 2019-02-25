Subroutine wvelocity(HydroParam,MeshParam,dt)
    
    !$ use omp_lib
    Use MeshVars
    Use Hydrodynamic
    
    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    Integer:: iElem,iLayer, DIM
    Double Precision:: aTh(MeshParam%KMax-1), bTh(MeshParam%KMax-1), cTh(MeshParam%KMax-1), fTh(MeshParam%KMax-1)
    Double Precision:: dt,dzk
    
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
                    fTh(iLayer) = HydroParam%Fw(iLayer+1,iElem)
                Elseif (iLayer == HydroParam%ElCapitalM(iElem)-1) Then
                    dzk         = 0.5*(HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer+1,iElem)) 
                    aTh(iLayer) = -dt/dzk*(HydroParam%VerEddyViscCell(iLayer,iElem)/HydroParam%DZit(iLayer,iElem))
                    bTh(iLayer) = 1. + dt/dzk*(HydroParam%VerEddyViscCell(iLayer,iElem)/HydroParam%DZit(iLayer,iElem) + HydroParam%VerEddyViscCell(iLayer+1,iElem)/HydroParam%DZit(iLayer+1,iElem))
                    cTh(iLayer) = 0.
                    fTh(iLayer) = HydroParam%Fw(iLayer+1,iElem)
                Else
                    dzk         = 0.5*(HydroParam%DZit(iLayer,iElem) + HydroParam%DZit(iLayer+1,iElem)) 
                    aTh(iLayer) = -dt/dzk*(HydroParam%VerEddyViscCell(iLayer,iElem)/HydroParam%DZit(iLayer,iElem))
                    bTh(iLayer) = 1. + dt/dzk*(HydroParam%VerEddyViscCell(iLayer,iElem)/HydroParam%DZit(iLayer,iElem) + HydroParam%VerEddyViscCell(iLayer+1,iElem)/HydroParam%DZit(iLayer+1,iElem))
                    cTh(iLayer) = -dt/dzk*(HydroParam%VerEddyViscCell(iLayer+1,iElem)/HydroParam%DZit(iLayer+1,iElem))
                    fTh(iLayer) = HydroParam%Fw(iLayer+1,iElem)
                EndIf
            EndDo
            DIM = size(HydroParam%w(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1,iElem),1)
            Call solve_tridiag(cTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1),bTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1),aTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1),fTh(HydroParam%ElSmallm(iElem):HydroParam%ElCapitalM(iElem)-1),HydroParam%w(HydroParam%ElSmallm(iElem)+1:HydroParam%ElCapitalM(iElem),iElem),DIM)
        EndIf
    EndDo
    
    Return    
End Subroutine wvelocity