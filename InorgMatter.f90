SUBROUTINE InorgMatter(HydroParam,MeshParam,LimnoParam,dt,dtday)

    ! Inorganic Matter Modelling in the Water Column
    
    ! List of Modifications:
    !   24.07.2015: Routine Implementation      (Carlos Ruberto)
    !   25.07.2015: Routine Validation      (Carlos Ruberto)
    ! Programmer: Carlos Ruberto
    
    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars

    Implicit None
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    Integer:: index,iElem,iLayer,iter
    Real:: tPResusIM_bottom,wDSetIM_Layer
    Real:: dt,dtday
    Real:: V
    
    
    index = 3
    !-------------------------------------------------------------------------------------!
    !                               Boundary Condition for IM       				      !
    !-------------------------------------------------------------------------------------!
    HydroParam%uLoadVarEst = LimnoParam%uDLoadIM
    
    Do iElem = 1,MeshParam%nElem
        Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
            !-------------------------------------------------------------------------------------!
            !                               IM Processes in Water					      !
            !-------------------------------------------------------------------------------------!
            !-IM----------------------------------------------------------------------------------!
            ! 1. Erosion                                                                          !
            ! 2. Sedimentation                                                                    !
            ! 3. Ressuspension                                                                    !
            !-------------------------------------------------------------------------------------!
            !-------------------------------------------------------------------------!
			! 1) Erosion (gD/m2/d)                   
		    If (MeshParam%Neighbor(1,iElem)==0.or.MeshParam%Neighbor(2,iElem)==0.or.MeshParam%Neighbor(3,iElem)==0.or.MeshParam%Neighbor(4,iElem)==0) Then
                !  IM_input_from_banks (gD/m2/d)
			    LimnoParam%uDErosIM=(1.0-LimnoParam%fDOrgSoil)*LimnoParam%cDErosTot
                !  IM_input_to_sediment_from_banks (gD/m2/d)
			    LimnoParam%uDErosIMS(iElem)=LimnoParam%fSedErosIM*LimnoParam%uDErosIM
                !  IM_input_to_water_column_from_banks (gD/m2/d)
			    LimnoParam%uDErosIMW=LimnoParam%uDErosIM-LimnoParam%uDErosIMS(iElem)
		    Else
			    LimnoParam%uDErosIM=0.0                        
			    LimnoParam%uDErosIMS(iElem)=0.0		                
			    LimnoParam%uDErosIMW=0.0		          
            EndIf 

            !2) IM Sedimentation
            If ( iLayer == HydroParam%ElCapitalM(iElem) ) Then        ! Bottom layer or intermediary layer has contribuition of the up layer
                wDSetIM_Layer = LimnoParam%tDSetIM(iLayer,iElem)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))
            Else ! Bottom layer or intermediary layer has contribuition of the up layer
                wDSetIM_Layer = LimnoParam%tDSetIM(iLayer,iElem)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem)) - LimnoParam%tDSetIM(iLayer+1,iElem)/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer+1,iElem))
            Endif
            
            !3) IM Resuspension 
            If (iLayer == HydroParam%ElSmallm(iElem)) Then
                tPResusIM_bottom = LimnoParam%tDResusIM(iElem)
            Else
                tPResusIM_bottom = 0.
            EndIf
            
            !-------------------------------------------------------------------------------------!
            !                               Source/Sink Term for IM                               !
            !-------------------------------------------------------------------------------------!
            HydroParam%dVarEst(iLayer,iElem,1) = - wDSetIM_Layer                                                        & !Sedimentation
                                    + tPResusIM_bottom/Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))                & ! Resuspension
                                    + LimnoParam%uDErosIMW/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))*Max(HydroParam%Pcri,HydroParam%Dzi(iLayer,iElem))/V(HydroParam%eta(iElem)-HydroParam%etaplus(iElem),HydroParam%hb(iElem))            !Erosion
            
        EndDo !Loop Layer
    EndDo !Loop Cell
      
    !-------------------------------------------------------------------------------------!
    !                       Solver Transport Equation for IM               	    	      !
    !-------------------------------------------------------------------------------------!
    If (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC Scheme
        Call UpWindCWC(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDIMW,LimnoParam%sDIMWP,dt,dtday,HydroParam,MeshParam) 
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==0) Then ! UpWind CWC High Resolution Scheme
        Call UpWindCWC_HR(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDIMW,LimnoParam%sDIMWP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag) 
    ElseIf (LimnoParam%iTranspFlag == 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC with Global Time Stepping Scheme
        Call UpWindCWC_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDIMW,LimnoParam%sDIMWP,dt,dtday,HydroParam,MeshParam)
    ElseIf (LimnoParam%iTranspFlag > 0.and.LimnoParam%iAdaptativeTimeStep==1) Then ! UpWind CWC High Resolution with Global Time Stepping Scheme
        Call UpWindCWC_HR_GATS(index,HydroParam%uLoadVarEst,HydroParam%dVarEst(:,:,1),LimnoParam%sDIMW,LimnoParam%sDIMWP,dt,dtday,HydroParam,MeshParam,LimnoParam%iTranspFlag)
    EndiF        
    
Return
End
