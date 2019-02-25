  !> ompute Water Density based on Temperature and Salt
  !>\param a a vector to be multiplied
  !>\param b a vector to be multiplied
  !>\return c the resultant vector
  !>\author  Rafael Cavalcanti
  !>\attention  List of Modifications: \n
  !> - 04.03.2014: Routine Implementation (Rafael Cavalcanti)

Subroutine WaterDensity(HydroParam,MeshParam,LimnoParam)
    
    ! Compute Water Density based on Temperature and Salt
    
    ! Based on: 
    ! [1] Read, J.S.; Hamilton, D.P.; Jones, I.D.; Muraoka, K.; Winslow, L.A.; Kroiss, R.; Wu, C.H.; Gaiser, E. 2011. Derivation of lake mixing and stratification indices from high-resolution lake buoy data
    !   Environmental Modelling & Software, 26, pp. 1325-1336.
    ! [2] Ji, Z-G. Hydrodynamics and Water Quality: modeling rivers, lakes, and estuaries. 
    !   John Wiley & Sons, 2008. (p. 15)
    ! [3] Chen,C.-T.; Millero, F.J. The use and misuse of pure water PVT properties for lake waters.
    !   Letters to Nature - Nature, vol. 266, pp. 707-708, 1977.

    ! Input:
    ! T -> Water Temperature (°C)
    ! S -> Water Salt Content (parts per thousand)
    
    ! Output:
    ! Rho -> Water Density (Kg/m³)
    
    ! List of Modifications:
    !   13.12.2016: Routine Implementation       (J. Rafael Cavalcanti)
    !   13.09.2017: Chen and Millero Density     (J. Rafael Cavalcanti)

    ! Programmer: J. Rafael Cavalcanti
    Use MeshVars
    Use Hydrodynamic
    Use LimnologyVars
    

    Implicit None
    Integer, Parameter:: iDens=2
    Integer:: iElem, iLayer
    Real:: Spsu
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    
    Select Case (iDens)
    Case(1)         ! UNESCO [1,2]
        Do iElem = 1,MeshParam%nElem
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                HydroParam%sDRhoW(iLayer,iElem) = 0. 
                ! 1. Water Temperature Effects
                HydroParam%sDRhoW(iLayer,iElem) = 999.842594 + 6.793952*1.E-2*LimnoParam%sDTempW(iLayer,iElem) - 9.095290*1.E-3*LimnoParam%sDTempW(iLayer,iElem)**2. + 1.001685*1.E-4*LimnoParam%sDTempW(iLayer,iElem)**3. - 1.120083*1.E-6*LimnoParam%sDTempW(iLayer,iElem)**4. + 6.536332*1.E-9*LimnoParam%sDTempW(iLayer,iElem)**5.
                ! 2. Water Salinity Effects
                HydroParam%sDRhoW(iLayer,iElem) = HydroParam%sDRhoW(iLayer,iElem) + LimnoParam%sDSal(iLayer,iElem)*( 0.824493 - 4.0899*1.E-3*LimnoParam%sDTempW(iLayer,iElem) + 7.6438*1.E-5*LimnoParam%sDTempW(iLayer,iElem)**2. - 8.2467*1.E-7*LimnoParam%sDTempW(iLayer,iElem)**3. + 5.3875*1.E-9*LimnoParam%sDTempW(iLayer,iElem)**4. ) + ( -5.72466*1.E-3 + 1.0227*1.E-4*LimnoParam%sDTempW(iLayer,iElem) - 1.6546*1.E-6*LimnoParam%sDTempW(iLayer,iElem)**2. )*LimnoParam%sDSal(iLayer,iElem)**(3/2.) + ( 4.8314*1.E-4 )*LimnoParam%sDSal(iLayer,iElem)**2.
                ! 3. Turbidity Effects (Total Suspended Solid and Total Dissolved Solids)
                Continue
            EndDo       ! Layer Loop
        End Do       ! Element Loop
    Case(2)         ! Chen and Millero [3]
        Do iElem = 1,MeshParam%nElem
            Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
                HydroParam%sDRhoW(iLayer,iElem) = 0. 
                ! 1. Water Temperature Effects
                HydroParam%sDRhoW(iLayer,iElem) = 0.999835 + 6.7914E-5*LimnoParam%sDTempW(iLayer,iElem) - 9.0894E-6*LimnoParam%sDTempW(iLayer,iElem)**2. + 1.0171E-7*LimnoParam%sDTempW(iLayer,iElem)**3. - 1.2846E-9*LimnoParam%sDTempW(iLayer,iElem)**4. + 1.1592E-11*LimnoParam%sDTempW(iLayer,iElem)**5. - 5.0125E-14*LimnoParam%sDTempW(iLayer,iElem)**6.
                ! 2. Water Salinity Effects
                Spsu = LimnoParam%sDSal(iLayer,iElem)*0.995307
                HydroParam%sDRhoW(iLayer,iElem) = HydroParam%sDRhoW(iLayer,iElem) + Spsu*( 8.221E-4 - 3.87E-6*LimnoParam%sDTempW(iLayer,iElem) + 4.99E-8*LimnoParam%sDTempW(iLayer,iElem)**2 )
                ! 3. Turbidity Effects (Total Suspended Solid and Total Dissolved Solids)
                Continue
                
                HydroParam%sDRhoW(iLayer,iElem) = HydroParam%sDRhoW(iLayer,iElem)*1000.
            End Do       ! Layer Loop
        End Do       ! Element Loop 
    Case Default
        Print*, '**************************************************************'
        Print*, 'Something Wrong...'
        Print*, 'UnIPHECO could not Recognize the Density Equation Option...'
        Print*, 'It is not your fault. Talk to the development team to solve it.'
        Pause
        Stop
    End Select
     Do iElem = 1,MeshParam%nElem
        Do iLayer =1,  HydroParam%ElSmallm(iElem)-1
            HydroParam%sDRhoW(iLayer,iElem) = HydroParam%sDRhoW( HydroParam%ElSmallm(iElem),iElem)
        enddo
     enddo
     
    
    Return
End Subroutine WaterDensity