!> This subroutine reads the simulation parameters. 
Subroutine ReadGridData(GridDataConfig,MeshParam,HydroParam)
    
    Use domain_types
    Use Hydrodynamic
    Use MeshVars
    
    Implicit none
    Integer i, iElem
    type(GridDataConfiguration) :: GridDataConfig    
    type(MeshGridParam) :: MeshParam
    type(GridData), pointer :: Layers(:)
    Real, pointer :: weights(:)
    type(HydrodynamicParam) :: HydroParam
    
    call c_f_pointer(GridDataConfig%layers, Layers, [GridDataConfig%numberOfLayers])
    
    MeshParam%iWindRed = 0
    MeshParam%iWetland = 0
    MeshParam%id50 = 0
    MeshParam%iOMfraction = 0
    MeshParam%iBedrock = 0
    MeshParam%ieta0 = 0
    
    MeshParam%CREDV = 1.d0
    MeshParam%BANHADO = 0
    MeshParam%d50 = 0.d0
    MeshParam%OMfraction = 0.d0
    MeshParam%eta0 = 0.d0
    
    Do i = 1,GridDataConfig%numberOfLayers
        If (Layers(i)%typeid == 1) Then !< Sediment level (Bathymetry) 
            call c_f_pointer(Layers(i)%weights, weights, [Layers(i)%numberOfElements])
            Do iElem = 1, MeshParam%nElem
                HydroParam%hb(iElem) = weights(iElem)
            EndDo
        ElseIf (Layers(i)%typeid == 2) Then !< Roughness 
            If (HydroParam%iRoughForm == 0) Then ! roughnessChezyConstant
                HydroParam%Rug = HydroParam%RugChezyConst
            ElseIf (HydroParam%iRoughForm == 1) Then ! roughnessManningConstant
                Do iElem = 1, MeshParam%nElem
                    HydroParam%Rug(iElem) = HydroParam%RugManConst !
                EndDo
            ElseIf (HydroParam%iRoughForm == 2) Then ! roughnessWhiteColebrookConstant
                Do iElem = 1, MeshParam%nElem
                    HydroParam%Rug(iElem) = HydroParam%RugWCConst !18.*log10(12.*Max(HydroParam%Pcri,HydroParam%eta(iElem)-HydroParam%hb(iElem))/(HydroParam%RugWCConst/30.))
                EndDo
            ElseIf (HydroParam%iRoughForm == 3) Then ! roughnessChezyUseGridData
                call c_f_pointer(Layers(i)%weights, weights, [Layers(i)%numberOfElements])
                Do iElem = 1, MeshParam%nElem
                    HydroParam%Rug(iElem) = weights(iElem)
                EndDo
            ElseIf (HydroParam%iRoughForm == 4) Then ! roughnessManningUseGridData
                call c_f_pointer(Layers(i)%weights, weights, [Layers(i)%numberOfElements])
                Do iElem = 1, MeshParam%nElem
                    HydroParam%Rug(iElem) = weights(iElem) !Max(HydroParam%Pcri,HydroParam%eta(iElem)-HydroParam%hb(iElem))**(1./6.)/weights(iElem)
                EndDo
            ElseIf (HydroParam%iRoughForm == 5) Then ! roughnessWhiteColebrookUseGridData
                call c_f_pointer(Layers(i)%weights, weights, [Layers(i)%numberOfElements])
                Do iElem = 1, MeshParam%nElem
                    HydroParam%Rug(iElem) = weights(iElem) !18.*log10(12.*Max(HydroParam%Pcri,HydroParam%eta(iElem)-HydroParam%hb(iElem))/(weights(iElem)/30.)) 
                EndDo
            EndIf
        ElseIf (Layers(i)%typeid == 3) Then !< Wind reduction coefficient 
            MeshParam%iWindRed = 1
            call c_f_pointer(Layers(i)%weights, weights, [Layers(i)%numberOfElements])
            Do iElem = 1, MeshParam%nElem
                MeshParam%CREDV(iElem) = weights(iElem)
            EndDo
        ElseIf (Layers(i)%typeid == 4) Then !< Wetland 
            MeshParam%iWetland = 1
            call c_f_pointer(Layers(i)%weights, weights, [Layers(i)%numberOfElements])
            Do iElem = 1, MeshParam%nElem
                MeshParam%BANHADO(iElem) = weights(iElem)
            EndDo
        ElseIf (Layers(i)%typeid == 5) Then !< D50 GRAIN SIZE 
            MeshParam%id50 = 1
            call c_f_pointer(Layers(i)%weights, weights, [Layers(i)%numberOfElements])
            Do iElem = 1, MeshParam%nElem
                MeshParam%d50(iElem) = weights(iElem)
            EndDo
            
        ElseIf (Layers(i)%typeid == 6) Then !< FRACTION OF ORGANIC MATTER 
            MeshParam%iOMfraction = 1
            call c_f_pointer(Layers(i)%weights, weights, [Layers(i)%numberOfElements])
            Do iElem = 1, MeshParam%nElem
                MeshParam%OMfraction(iElem) = weights(iElem)
            EndDo
        ElseIf (Layers(i)%typeid == 7) Then !< Bedrock elevation 
            MeshParam%iBedrock = 1
            call c_f_pointer(Layers(i)%weights, weights, [Layers(i)%numberOfElements])
            Do iElem = 1, MeshParam%nElem
                HydroParam%sb(iElem) = weights(iElem)
            EndDo
            
        ElseIf (Layers(i)%typeid == 8) Then !< Free-surface elevation 
            MeshParam%ieta0 = 1
            call c_f_pointer(Layers(i)%weights, weights, [Layers(i)%numberOfElements])
            Do iElem = 1, MeshParam%nElem
                MeshParam%eta0(iElem) = weights(iElem)
            EndDo
        EndIf
    EndDo
    
    If (MeshParam%iBedrock == 0) Then
            HydroParam%sb =  HydroParam%hb
    EndIf

   
End Subroutine ReadGridData
    