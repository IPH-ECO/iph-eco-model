Subroutine Tension(SurfTensionFlag,BottomTensionFlag,iEdge,kSurf,kBottom,g,ub,vb,us,vs,Chezy,rhoair,rho0,CV,WindEdge,WindEdgeXY,GammaT,GammaB) 

    ! Calculate the Shear Stress on the surface and bottom layers
    ! Based on:

    ! Input:
    ! Hydrodynamic Features
    ! Output:
    ! Shear Stress in the Surface (GammaT)
    ! Shear Stress in the Bottom (GammaB)
    
    ! List of Modifications:
    !   15.06.2015: Routine Implementation       (Carlos Ruberto)
    ! Programmer: Carlos Ruberto    
    
    Implicit None
    Integer:: SurfTensionFlag,BottomTensionFlag,iEdge
    Integer:: kSurf,kBottom
    Real:: g,ub,vb,us,vs,Chezy,rhoair,rho0,GammaT,GammaB,WindEdge(2),CV,WindEdgeXY(2)

    ! 1. Shear stress in water surface
    If (SurfTensionFlag==1) Then ! Formulation 1 - Based on Wind Drag coefficient
        GammaT = CV*(WindEdgeXY(1)**2.+WindEdgeXY(2)**2.)**0.5 
        !GammaT = 0.013/(rho0*sqrt((WindEdgeXY(1)-us)**2. + (WindEdgeXY(2)-vs)**2.)) 
    ElseIf (SurfTensionFlag==0) Then ! Formulation 2 - Based on Air density
        GammaT = rhoair/rho0*(0.63+0.066*sqrt(WindEdgeXY(1)**2.+WindEdgeXY(2)**2.))*0.001*sqrt((WindEdgeXY(1)-us)**2. + (WindEdgeXY(2)-vs)**2.)
    EndIf
    
    ! 2. Shear stress in the bottom
    If (BottomTensionFlag==0) Then
        
        If (Chezy > 0) Then 
            GammaB = (g*(ub**2.+vb**2.)**0.5)/(Chezy**2.)
        Else
            GammaB = 0.d0
        EndIf
        
    EndIf

    
    

    
    
    Return
End Subroutine Tension