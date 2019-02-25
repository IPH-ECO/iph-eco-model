!> This subroutine reads the hydrodynamic parameters. 
Subroutine ReadWQIniCond(HydroParam,MeshParam,LimnoParam,simParam)
	
    Use Hydrodynamic
	Use LimnologyVars
    Use MeshVars
    Use SimulationModel
    
	Implicit None
	
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    type(MeshGridParam) :: MeshParam
    type(SimulationParam) :: simParam
    Integer:: bb,gg,mm,zz,ben,ff,pom,dom,iElem,iLayer
    Real:: aAlkW,aAlkS
	Real:: Hion
	Real:: k1,k2,kwater
	Real:: a0,a1,a2,kzeta, k_t,a_t
    Real:: dt
    Integer:: i,ano2,diajulcont1,dayoftheyear
    Real:: delta,cos_hss,sunSetHourAngle,dayPhotoFration,idate(6)
    character(len=200):: text
    
    HydroParam%Locdt = simParam%dt

    If (LimnoParam%iTempW == 1) Then
	    LimnoParam%sDTempW = LimnoParam%sDTempW0           ! Water Temperature (Current Time Step)
	    LimnoParam%sDTempWP = LimnoParam%sDTempW            ! Water Temperatute (Past Time Step)
        !LimnoParam%sDTempW(1:5,:)  = 12
        !LimnoParam%sDTempW(6,:)  = 20
        
    Else
        LimnoParam%sDTempW = LimnoParam%WtempRef
        LimnoParam%sDTempWP = LimnoParam%sDTempW
    EndIf

    If (LimnoParam%iSal == 1) Then
        
	    LimnoParam%sDSal  = LimnoParam%sDSal0 ! Salinity (Current Time Step)
	    LimnoParam%sDSalP = LimnoParam%sDSal ! Salinity (Past Time Step)
        
        !!Lock exchange case
        !Do iElem = 1,MeshParam%nElem
        !    k_t = HydroParam%pi
        !    a_t = 0.1d0/k_t
        !    kzeta=k_t*a_t*((1-(k_t*a_t)**2./64.)*COS(k_t*(MeshParam%xb(iElem)))-(k_t*a_t)**2./8.*COS(3.*k_t*(MeshParam%xb(iElem))))
        !    Do iLayer = HydroParam%ElSmallm(iElem), HydroParam%ElCapitalM(iElem)
        !        LimnoParam%sDSal(iLayer,iElem)=39.5082/2.-39.5082/2.*TANH(2.*(1./TANH(0.99))/(0.05*HydroParam%pi)*(k_t*HydroParam%Zb(iLayer,iElem)-kzeta+k_t*1./2.))
        !        HydroParam%sDRhoW(iLayer,iElem)=985.+((-30./2.)*TANH(2.*(1./TANH(0.99))/(0.05*k_t)*(k_t*HydroParam%Zb(iLayer,iElem)-kzeta+k_t*1./2.)))
        !    EndDo
            
            
        !    If (MeshParam%xb(iElem)<15) Then
	       !     LimnoParam%sDSal(:,iElem)  = 1 ! Salinity (Current Time Step)
	       !     LimnoParam%sDSalP(:,iElem) = 1 ! Salinity (Past Time Step)
	       !     !LimnoParam%sDTempW(:,iElem)  = 15.8 ! Salinity (Current Time Step)
	       !     !LimnoParam%sDTempWP(:,iElem) = 15.8 ! Salinity (Past Time Step)
        !    Else
	       !     LimnoParam%sDSal(:,iElem)  = 0 ! Salinity (Current Time Step)
	       !     LimnoParam%sDSalP(:,iElem) = 0 ! Salinity (Past Time Step)
	       !     !LimnoParam%sDTempW(:,iElem)  = 20 ! Salinity (Current Time Step)
	       !     !LimnoParam%sDTempWP(:,iElem) = 20 ! Salinity (Past Time Step)
        !    EndIf
        !EndDo
        
    Else
        LimnoParam%sDSal  = 0.d0
        LimnoParam%sDSalP = 0.d0
    EndIf
        
!------------------------------------------------------------------------------------------------------------------!	
!  M O D U L O    E C O L O G I C O  =>  B I O T I C O / E C O L O G I C O                                         !	
!------------------------------------------------------------------------------------------------------------------!	
!------------------------------------------------------------------------------------------------------------------!	
!  B I O T I C O                                                                                                   !	
!------------------------------------------------------------------------------------------------------------------!	
!------------------------------------------------------------------------------------------------------------------!	
!  B A C T E R I O P L A N C T O N                                                                                 !	
!------------------------------------------------------------------------------------------------------------------!
	If (LimnoParam%InclBac==1) Then
		Do bb=1,LimnoParam%nbac
	        !Water
			LimnoParam%sDBacW(:,:,bb) = LimnoParam%sDBacW0(bb)
			LimnoParam%sCBacW(:,:,bb) = LimnoParam%cCDBacRef(bb) * LimnoParam%sDBacW(:,:,bb)  
			LimnoParam%sNBacW(:,:,bb) = LimnoParam%cNDBacRef(bb) * LimnoParam%sDBacW(:,:,bb)  
			LimnoParam%sPBacW(:,:,bb) = LimnoParam%cPDBacRef(bb) * LimnoParam%sDBacW(:,:,bb)  
			LimnoParam%rCDBacW(:,:,bb) = LimnoParam%cCDBacRef(bb)
			LimnoParam%rNDBacW(:,:,bb) = LimnoParam%cNDBacRef(bb)
			LimnoParam%rPDBacW(:,:,bb) = LimnoParam%cPDBacRef(bb) 
			!Sediment
			LimnoParam%sDBacS(:,bb) = LimnoParam%sDBacS0(bb)
			LimnoParam%sCBacS(:,bb) = LimnoParam%cCDBacRef(bb) * LimnoParam%sDBacS(:,bb)  
			LimnoParam%sNBacS(:,bb) = LimnoParam%cNDBacRef(bb) * LimnoParam%sDBacS(:,bb)  
			LimnoParam%sPBacS(:,bb) = LimnoParam%cPDBacRef(bb) * LimnoParam%sDBacS(:,bb)  
			LimnoParam%rCDBacS(:,bb) = LimnoParam%cCDBacRef(bb) 
			LimnoParam%rNDBacS(:,bb) = LimnoParam%cNDBacRef(bb)
			LimnoParam%rPDBacS(:,bb) = LimnoParam%cPDBacRef(bb)
		EndDo
    Else
	    !Water
		LimnoParam%sDBacW = 0.d0
		LimnoParam%sCBacW = 0.d0
		LimnoParam%sNBacW = 0.d0  
		LimnoParam%sPBacW = 0.d0 
		LimnoParam%rCDBacW = 0.d0
		LimnoParam%rNDBacW = 0.d0
		LimnoParam%rPDBacW = 0.d0
		!Sediment
		LimnoParam%sDBacS = 0.d0
		LimnoParam%sCBacS = 0.d0 
		LimnoParam%sNBacS = 0.d0  
		LimnoParam%sPBacS = 0.d0  
		LimnoParam%rCDBacS = 0.d0
		LimnoParam%rNDBacS = 0.d0
		LimnoParam%rPDBacS = 0.d0
    EndIf
	!Processes
	LimnoParam%wDAssBac = 0.d0  
	LimnoParam%wCAssBac = 0.d0  
	LimnoParam%wNAssBac = 0.d0  
	LimnoParam%wPAssBac = 0.d0  
	LimnoParam%wSiAssBac = 0.d0
	LimnoParam%wDAssBacDom = 0.d0  
	LimnoParam%wCAssBacDom = 0.d0  
	LimnoParam%wNAssBacDom = 0.d0  
	LimnoParam%wPAssBacDom = 0.d0 
	LimnoParam%wSiAssBacDom = 0.d0  
	LimnoParam%wCRespBac2CO2 = 0.d0  
	LimnoParam%wCRespBac2CH4 = 0.d0  
	LimnoParam%wDRespBac = 0.d0  
    LimnoParam%wCRespBac = 0.d0
    LimnoParam%wNExcrBac = 0.d0
    LimnoParam%wPExcrBac = 0.d0
	LimnoParam%wPUptBac = 0.d0  
	LimnoParam%wNUptBac = 0.d0  
	LimnoParam%wDMortBac = 0.d0  
	LimnoParam%wCMortBac = 0.d0  
	LimnoParam%wNMortBac = 0.d0  
	LimnoParam%wPMortBac = 0.d0  
	LimnoParam%tDSetBac = 0.d0  
	LimnoParam%tCSetBac = 0.d0  
	LimnoParam%tNSetBac = 0.d0  
	LimnoParam%tPSetBac = 0.d0
    LimnoParam%tDMortBac = 0.d0
    LimnoParam%tCMortBac = 0.d0
    LimnoParam%tNMortBac = 0.d0
    LimnoParam%tPMortBac = 0.d0
    LimnoParam%tPUptBac = 0.d0
    LimnoParam%tNUptBac = 0.d0
    LimnoParam%tNExcrBac = 0.d0
    LimnoParam%tPExcrBac = 0.d0
    LimnoParam%tDRespBac = 0.d0
    LimnoParam%tCRespBac = 0.d0
    LimnoParam%tCRespBac2CO2 = 0.d0
    LimnoParam%tCRespBac2CH4 = 0.d0
    LimnoParam%tDAssBac = 0.d0
    LimnoParam%tCAssBac = 0.d0
    LimnoParam%tNAssBac = 0.d0
    LimnoParam%tPAssBac = 0.d0
    LimnoParam%tDAssBacDom = 0.d0
    LimnoParam%tCAssBacDom = 0.d0 
    LimnoParam%tNAssBacDom = 0.d0 
    LimnoParam%tPAssBacDom = 0.d0
    
    
!------------------------------------------------------------------------------------------------------------------!	
!  F I T O P L A N C T O N                                                                                         !	
!------------------------------------------------------------------------------------------------------------------!
	If (LimnoParam%InclPhyt==1) Then
		Do gg=1,LimnoParam%nphy
	        !Water
			LimnoParam%sDPhytW(:,:,gg)  = LimnoParam%sDPhytW0(gg) 
			LimnoParam%sCPhytW(:,:,gg)  = LimnoParam%cCDPhytOpt(gg) * LimnoParam%sDPhytW(:,:,gg)  
			LimnoParam%sNPhytW(:,:,gg)  = LimnoParam%cNDPhytOpt(gg) * LimnoParam%sDPhytW(:,:,gg)  
			LimnoParam%sPPhytW(:,:,gg)  = LimnoParam%cPDPhytOpt(gg) * LimnoParam%sDPhytW(:,:,gg)  
			LimnoParam%sSiPhytW(:,:,gg) = LimnoParam%cSiDPhytOpt(gg) * LimnoParam%sDPhytW(:,:,gg)
			LimnoParam%rCDPhytW(:,:,gg) = LimnoParam%cCDPhytOpt(gg)
			LimnoParam%rNDPhytW(:,:,gg) = LimnoParam%cNDPhytOpt(gg) 
			LimnoParam%rPDPhytW(:,:,gg) = LimnoParam%cPDPhytOpt(gg) 
			LimnoParam%rSiDPhytW(:,:,gg)= LimnoParam%cSiDPhytOpt(gg)
			!Sediment
			LimnoParam%sDPhytS(:,gg) = LimnoParam%sDPhytS0(gg)
			LimnoParam%sCPhytS(:,gg) = LimnoParam%cCDPhytOpt(gg) * LimnoParam%sDPhytS(:,gg)  
			LimnoParam%sNPhytS(:,gg) = LimnoParam%cNDPhytOpt(gg) * LimnoParam%sDPhytS(:,gg)  
			LimnoParam%sPPhytS(:,gg) = LimnoParam%cPDPhytOpt(gg) * LimnoParam%sDPhytS(:,gg)  
			LimnoParam%sSiPhytS(:,gg)= LimnoParam%cSiDPhytOpt(gg)* LimnoParam%sDPhytS(:,gg)
			LimnoParam%rCDPhytS(:,gg) = LimnoParam%cCDPhytOpt(gg) 
			LimnoParam%rNDPhytS(:,gg) = LimnoParam%cNDPhytOpt(gg) 
			LimnoParam%rPDPhytS(:,gg) = LimnoParam%cPDPhytOpt(gg) 
			LimnoParam%rSiDPhytS(:,gg)= LimnoParam%cSiDPhytOpt(gg)

		EndDo
    Else
	    !Zona Pelagica
		LimnoParam%sDPhytW = 0.d0 
		LimnoParam%sCPhytW = 0.d0
		LimnoParam%sNPhytW = 0.d0
		LimnoParam%sPPhytW = 0.d0 
		LimnoParam%sSiPhytW = 0.d0
		LimnoParam%rCDPhytW = 0.d0
		LimnoParam%rNDPhytW = 0.d0
		LimnoParam%rPDPhytW = 0.d0
		LimnoParam%rSiDPhytW= 0.
		!Sedimento
		LimnoParam%sDPhytS = 0.d0
		LimnoParam%sCPhytS = 0.d0
		LimnoParam%sNPhytS = 0.d0
		LimnoParam%sPPhytS = 0.d0
		LimnoParam%sSiPhytS = 0.d0
		LimnoParam%rCDPhytS = 0.d0
		LimnoParam%rNDPhytS = 0.d0
		LimnoParam%rPDPhytS = 0.d0
		LimnoParam%rSiDPhytS = 0.d0
    EndIf
	!Processes
	!Consumo de nutrientes
	LimnoParam%wCUptPhyt = 0.d0
	LimnoParam%wNUptPhyt = 0.d0
	LimnoParam%wPUptPhyt = 0.d0
	LimnoParam%wSiUptPhyt = 0.d0 
	LimnoParam%wNUptNH4Phyt = 0.d0
	LimnoParam%wNUptNO3Phyt = 0.d0                          
	!Produção Primária
	LimnoParam%wDAssPhyt = 0.d0
	!Respiração 
	LimnoParam%wDRespbPhyt = 0.d0
	LimnoParam%wDResppPhyt = 0.d0
	LimnoParam%wDRespPhyt = 0.d0
	LimnoParam%tDRespPhytS = 0.d0
    LimnoParam%wCExcrPhytW = 0.d0
    LimnoParam%wPExcrPhytW = 0.d0
    LimnoParam%wNExcrPhytW = 0.d0
    LimnoParam%tCExcrPhytS = 0.d0
    LimnoParam%tPExcrPhytS = 0.d0
    LimnoParam%tNExcrPhytS = 0.d0
	!Resuspensão
	LimnoParam%tDResusPhyt = 0.d0
	LimnoParam%tCResusPhyt = 0.d0
	LimnoParam%tPResusPhyt = 0.d0
	LimnoParam%tNResusPhyt = 0.d0
	LimnoParam%tSiResusPhyt = 0.d0
	!Sedimentação
	LimnoParam%tDSetPhyt = 0.d0
	!Lise
	LimnoParam%wDLisPhyt = 0.d0
	LimnoParam%wCLisPhyt = 0.d0
	LimnoParam%wNLisPhyt = 0.d0
	LimnoParam%wPLisPhyt = 0.d0
	LimnoParam%wSiLisPhyt = 0.d0
	LimnoParam%tDLisPhytS = 0.d0
	LimnoParam%tCLisPhytS = 0.d0
	LimnoParam%tNLisPhytS = 0.d0
	LimnoParam%tPLisPhytS = 0.d0
	LimnoParam%tSiLisPhytS = 0.d0
	LimnoParam%tDLisPomPhytS = 0.d0
	LimnoParam%tCLisPomPhytS = 0.d0
	LimnoParam%tNLisPomPhytS = 0.d0
	LimnoParam%tPLisPomPhytS = 0.d0
	LimnoParam%tSiLisPomPhytS = 0.d0
	LimnoParam%tDLisDomPhytS = 0.d0
	LimnoParam%tCLisDomPhytS = 0.d0
	LimnoParam%tNLisDomPhytS = 0.d0
	LimnoParam%tPLisDomPhytS = 0.d0
	LimnoParam%tSiLisDomPhytS = 0.d0
    
!------------------------------------------------------------------------------------------------------------------!	
!  M A C R O F I T A S                                                                                             !	
!------------------------------------------------------------------------------------------------------------------!	
    If (LimnoParam%InclMacr==1) Then
	    Do mm=1,LimnoParam%nmac
			LimnoParam%rCDMac(:,mm) = LimnoParam%cCDMac0(mm)
			LimnoParam%rNDMac(:,mm) = LimnoParam%cNDMac0(mm)
			LimnoParam%rPDMac(:,mm) = LimnoParam%cPDMac0(mm)
                    
		    LimnoParam%sDMac(:,mm) = LimnoParam%sDMac0(mm)			
		    LimnoParam%sCMac(:,mm) = LimnoParam%sDMac(:,mm)*LimnoParam%rCDMac(:,mm)
		    LimnoParam%sNMac(:,mm) = LimnoParam%sDMac(:,mm)*LimnoParam%rNDMac(:,mm)
		    LimnoParam%sPMac(:,mm) = LimnoParam%sDMac(:,mm)*LimnoParam%rPDMac(:,mm)
		EndDo
    Else
   		LimnoParam%sDMac = 0.d0	
		LimnoParam%sCMac = 0.d0
		LimnoParam%sNMac = 0.d0
		LimnoParam%sPMac = 0.d0
		LimnoParam%rCDMac = 0.d0
		LimnoParam%rNDMac = 0.d0
		LimnoParam%rPDMac = 0.d0
    EndIf
	!Processes
	!Consumo de Nutrientes
	LimnoParam%tCUptMacW = 0.d0
	LimnoParam%tNUptMacW = 0.d0
	LimnoParam%tPUptMacW = 0.d0
	LimnoParam%tNUptNH4MacW = 0.d0
	LimnoParam%tNUptNO3MacW = 0.d0
	LimnoParam%tCUptMacS = 0.d0
	LimnoParam%tNUptMacS = 0.d0
	LimnoParam%tPUptMacS = 0.d0
	LimnoParam%tNUptNH4MacS = 0.d0
	LimnoParam%tNUptNO3MacS = 0.d0
	!Produção Primária
	LimnoParam%tDProdMac = 0.d0
	!Exducao
	LimnoParam%tCExdMac = 0.d0
	LimnoParam%tNExdMac = 0.d0
	LimnoParam%tPExdMac = 0.d0
	LimnoParam%tCExdMacW = 0.d0
	LimnoParam%tNExdMacW = 0.d0
	LimnoParam%tPExdMacW = 0.d0
	LimnoParam%tCExdMacS = 0.d0
	LimnoParam%tNExdMacS = 0.d0
	LimnoParam%tPExdMacS = 0.d0			
	!Respiração
    LimnoParam%tDRespMac = 0.d0
	LimnoParam%tCRespMac = 0.d0
	LimnoParam%tDRespMacW = 0.d0
	LimnoParam%tDRespMacS = 0.d0
	!Mortalidade
    LimnoParam%tDMortMac = 0.d0
    LimnoParam%tCMortMac = 0.d0
    LimnoParam%tNMortMac = 0.d0
    LimnoParam%tPMortMac = 0.d0
	LimnoParam%tDMortMacW = 0.d0
	LimnoParam%tCMortMacW = 0.d0
	LimnoParam%tNMortMacW = 0.d0
	LimnoParam%tPMortMacW = 0.d0
	LimnoParam%tDMortMacS = 0.d0
	LimnoParam%tCMortMacS = 0.d0
	LimnoParam%tNMortMacS = 0.d0
	LimnoParam%tPMortMacS = 0.d0
	!Consumo por aves
	LimnoParam%tDGrazMacBird = 0.d0
	LimnoParam%tCGrazMacBird = 0.d0
	LimnoParam%tNGrazMacBird = 0.d0
	LimnoParam%tPGrazMacBird = 0.d0
	LimnoParam%tDEgesBird = 0.d0
	LimnoParam%tCEgesBird = 0.d0
	LimnoParam%tNEgesBird = 0.d0
	LimnoParam%tPEgesBird = 0.d0	
!------------------------------------------------------------------------------------------------------------------!	
!  Z O O P L A N C T O N                                                                                           !	
!------------------------------------------------------------------------------------------------------------------!	
	If (LimnoParam%InclZoo==1) Then
		Do zz=1,LimnoParam%nzoo
			LimnoParam%sDZoo(:,:,zz) = LimnoParam%sDZoo0(zz)				
			LimnoParam%sCZoo(:,:,zz) = LimnoParam%cCDZooRef(zz) * LimnoParam%sDZoo(:,:,zz)
			LimnoParam%sPZoo(:,:,zz) = LimnoParam%cPDZooRef(zz) * LimnoParam%sDZoo(:,:,zz)
			LimnoParam%sNZoo(:,:,zz) = LimnoParam%cNDZooRef(zz) * LimnoParam%sDZoo(:,:,zz)
			LimnoParam%rCDZoo(:,:,zz)= LimnoParam%cCDZooRef(zz)
			LimnoParam%rPDZoo(:,:,zz)= LimnoParam%cPDZooRef(zz)
			LimnoParam%rNDZoo(:,:,zz)= LimnoParam%cNDZooRef(zz)	
		EndDo
	Else
 	    LimnoParam%sDZoo = 0.d0				
		LimnoParam%sCZoo = 0.d0
		LimnoParam%sPZoo = 0.d0
		LimnoParam%sNZoo = 0.d0
		LimnoParam%rCDZoo= 0.d0
		LimnoParam%rPDZoo= 0.d0
		LimnoParam%rNDZoo= 0.d0	
    EndIf
	!Processes
	LimnoParam%wDGrzZoo = 0.d0
	LimnoParam%wCGrzZoo = 0.d0
	LimnoParam%wNGrzZoo = 0.d0
	LimnoParam%wPGrzZoo = 0.d0
	LimnoParam%wDGrzZooZoo = 0.d0
	LimnoParam%wCGrzZooZoo = 0.d0
	LimnoParam%wNGrzZooZoo = 0.d0
	LimnoParam%wPGrzZooZoo = 0.d0
	LimnoParam%wDGrzZooPhyt = 0.d0
	LimnoParam%wCGrzZooPhyt = 0.d0
	LimnoParam%wNGrzZooPhyt = 0.d0
	LimnoParam%wPGrzZooPhyt = 0.d0
	LimnoParam%wSiGrzZooPhyt = 0.d0
	LimnoParam%wDGrzZooBac = 0.d0
	LimnoParam%wCGrzZooBac = 0.d0
	LimnoParam%wNGrzZooBac = 0.d0
	LimnoParam%wPGrzZooBac = 0.d0
	LimnoParam%wDGrzZooPom = 0.d0
	LimnoParam%wCGrzZooPom = 0.d0
	LimnoParam%wNGrzZooPom = 0.d0
	LimnoParam%wPGrzZooPom = 0.d0
	LimnoParam%wDRespZoo = 0.d0
	LimnoParam%wCRespZoo = 0.d0
	LimnoParam%wDEgesZoo = 0.d0
	LimnoParam%wCEgesZoo = 0.d0
	LimnoParam%wNEgesZoo = 0.d0
	LimnoParam%wPEgesZoo = 0.d0
	LimnoParam%wDMortZoo = 0.d0
	LimnoParam%wCMortZoo = 0.d0
	LimnoParam%wNMortZoo = 0.d0
	LimnoParam%wPMortZoo = 0.d0
	LimnoParam%wNExcrZoo = 0.d0
	LimnoParam%wPExcrZoo = 0.d0
    LimnoParam%wDMessZooOM = 0.d0
    LimnoParam%wCMessZooOM = 0.d0
    LimnoParam%wNMessZooOM = 0.d0
    LimnoParam%wPMessZooOM = 0.d0
    LimnoParam%wDAssZoo = 0.d0
    LimnoParam%wCAssZoo = 0.d0
    LimnoParam%wNAssZoo = 0.d0
    LimnoParam%wPAssZoo = 0.d0
    
!------------------------------------------------------------------------------------------------------------------!	
!  Z O O B E N T O S                                                                                               !	
!------------------------------------------------------------------------------------------------------------------!	
	If (LimnoParam%InclBent==1) Then
		Do ben=1,LimnoParam%nben
			LimnoParam%sDBent(:,ben) = LimnoParam%sDBent0(ben)				
			LimnoParam%sCBent(:,ben) = LimnoParam%cCDBentRef(ben) * LimnoParam%sDBent(:,ben)
			LimnoParam%sPBent(:,ben) = LimnoParam%cPDBentRef(ben) * LimnoParam%sDBent(:,ben)
			LimnoParam%sNBent(:,ben) = LimnoParam%cNDBentRef(ben) * LimnoParam%sDBent(:,ben)
			LimnoParam%rCDBent(:,ben)= LimnoParam%cCDBentRef(ben)
			LimnoParam%rPDBent(:,ben)= LimnoParam%cPDBentRef(ben)
			LimnoParam%rNDBent(:,ben)= LimnoParam%cNDBentRef(ben)		
		EndDo
	Else
		LimnoParam%sDBent = 0.d0			
		LimnoParam%sCBent = 0.d0
		LimnoParam%sPBent = 0.d0
		LimnoParam%sNBent = 0.d0
		LimnoParam%rCDBent= 0.d0
		LimnoParam%rPDBent= 0.d0
		LimnoParam%rNDBent= 0.d0
    EndIf
	!Processes
    LimnoParam%tDGrzBent = 0.d0
    LimnoParam%tCGrzBent = 0.d0
    LimnoParam%tNGrzBent = 0.d0
    LimnoParam%tPGrzBent = 0.d0
    LimnoParam%tDGrzBentPhyt = 0.d0
    LimnoParam%tCGrzBentPhyt = 0.d0
    LimnoParam%tNGrzBentPhyt = 0.d0
    LimnoParam%tPGrzBentPhyt = 0.d0
    LimnoParam%tSiGrzBentPhyt = 0.d0
    LimnoParam%tDGrzBentBac = 0.d0
    LimnoParam%tCGrzBentBac = 0.d0
    LimnoParam%tNGrzBentBac = 0.d0
    LimnoParam%tPGrzBentBac = 0.d0
    LimnoParam%tDGrzBentPom = 0.d0
    LimnoParam%tCGrzBentPom = 0.d0
    LimnoParam%tNGrzBentPom = 0.d0
    LimnoParam%tPGrzBentPom = 0.d0
    LimnoParam%tDRespBent = 0.d0
    LimnoParam%tCRespBent = 0.d0
    LimnoParam%tDEgesBent = 0.d0
    LimnoParam%tCEgesBent = 0.d0
    LimnoParam%tNEgesBent = 0.d0
    LimnoParam%tPEgesBent = 0.d0
    LimnoParam%tDMortBent = 0.d0
    LimnoParam%tCMortBent = 0.d0
    LimnoParam%tNMortBent  = 0.d0
    LimnoParam%tPMortBent = 0.d0
    LimnoParam%tNExcBent = 0.d0
    LimnoParam%tPExcBent = 0.d0
    LimnoParam%tDMessBentOM = 0.d0
    LimnoParam%tCMessBentOM = 0.d0
    LimnoParam%tNMessBentOM = 0.d0
    LimnoParam%tPMessBentOM = 0.d0
    
!------------------------------------------------------------------------------------------------------------------!	
!  P E I X E S                                                                                                     !	
!------------------------------------------------------------------------------------------------------------------!	
	If (LimnoParam%InclFish==1) Then
		Do ff=1,LimnoParam%nfish
			LimnoParam%sDFiAd(:,:,ff) = LimnoParam%sDFiAd0(ff)
			LimnoParam%sCFiAd(:,:,ff) = LimnoParam%cCDFishRef(ff) * LimnoParam%sDFiAd(:,:,ff)
			LimnoParam%sPFiAd(:,:,ff) = LimnoParam%cPDFishRef(ff) * LimnoParam%sDFiAd(:,:,ff)
			LimnoParam%sNFiAd(:,:,ff) = LimnoParam%cNDFishRef(ff) * LimnoParam%sDFiAd(:,:,ff)
            LimnoParam%sDFiJv(:,:,ff) = LimnoParam%sDFiJv0(ff)
			LimnoParam%sCFiJv(:,:,ff) = LimnoParam%cCDFishRef(ff) * LimnoParam%sDFiJv(:,:,ff)
			LimnoParam%sPFiJv(:,:,ff) = LimnoParam%cPDFishRef(ff) * LimnoParam%sDFiJv(:,:,ff)
			LimnoParam%sNFiJv(:,:,ff) = LimnoParam%cNDFishRef(ff) * LimnoParam%sDFiJv(:,:,ff)
			LimnoParam%rCDFiAd(:,:,ff)= LimnoParam%cCDFishRef(ff)
			LimnoParam%rPDFiAd(:,:,ff)= LimnoParam%cPDFishRef(ff)
			LimnoParam%rNDFiAd(:,:,ff)= LimnoParam%cNDFishRef(ff)		
			LimnoParam%rCDFiJv(:,:,ff)= LimnoParam%cCDFishRef(ff)
			LimnoParam%rPDFiJv(:,:,ff)= LimnoParam%cPDFishRef(ff)
			LimnoParam%rNDFiJv(:,:,ff)= LimnoParam%cNDFishRef(ff)		
		EndDo
	Else
		LimnoParam%sDFiAd = 0.d0		
		LimnoParam%sCFiAd = 0.d0
		LimnoParam%sPFiAd = 0.d0
		LimnoParam%sNFiAd = 0.d0
		LimnoParam%sDFiJv = 0.d0		
		LimnoParam%sCFiJv = 0.d0
		LimnoParam%sPFiJv = 0.d0
		LimnoParam%sNFiJv = 0.d0
        LimnoParam%rCDFiAd= 0.d0
		LimnoParam%rPDFiAd= 0.d0
		LimnoParam%rNDFiAd= 0.d0	
        LimnoParam%rCDFiJv= 0.d0
		LimnoParam%rPDFiJv= 0.d0
		LimnoParam%rNDFiJv= 0.d0		
    ENDIF
    !Processes
    LimnoParam%wDMigrFiAd = 0.d0
    LimnoParam%wCMigrFiAd = 0.d0
    LimnoParam%wNMigrFiAd = 0.d0
    LimnoParam%wPMigrFiAd = 0.d0
    LimnoParam%wDMigrFiJv = 0.d0
    LimnoParam%wCMigrFiJv = 0.d0
    LimnoParam%wNMigrFiJv = 0.d0
    LimnoParam%wPMigrFiJv = 0.d0
            
    LimnoParam%wDReprFish = 0.d0
    LimnoParam%wCReprFish = 0.d0   
    LimnoParam%wNReprFish = 0.d0
    LimnoParam%wPReprFish = 0.d0
    LimnoParam%wDAgeFish = 0.d0      
    LimnoParam%wCAgeFish = 0.d0      
    LimnoParam%wNAgeFish = 0.d0
    LimnoParam%wPAgeFish = 0.d0 
    LimnoParam%wDGrzFiAdFiAd = 0.d0
    LimnoParam%wCGrzFiAdFiAd = 0.d0
    LimnoParam%wNGrzFiAdFiAd = 0.d0
    LimnoParam%wPGrzFiAdFiAd = 0.d0
    LimnoParam%wDGrzFiAdFiJv = 0.d0
    LimnoParam%wCGrzFiAdFiJv = 0.d0
    LimnoParam%wNGrzFiAdFiJv = 0.d0
    LimnoParam%wPGrzFiAdFiAd = 0.d0
    LimnoParam%wDGrzFiAdZoo = 0.d0 
    LimnoParam%wCGrzFiAdZoo = 0.d0
    LimnoParam%wNGrzFiAdZoo = 0.d0
    LimnoParam%wPGrzFiAdZoo = 0.d0
    LimnoParam%wDGrzFiJvZoo = 0.d0 
    LimnoParam%wCGrzFiJvZoo = 0.d0
    LimnoParam%wNGrzFiJvZoo = 0.d0
    LimnoParam%wPGrzFiJvZoo = 0.d0
    LimnoParam%wDGrzFiAdPhyt = 0.d0
    LimnoParam%wCGrzFiAdPhyt = 0.d0
    LimnoParam%wNGrzFiAdPhyt = 0.d0
    LimnoParam%wPGrzFiAdPhyt = 0.d0
    LimnoParam%wDGrzFiJvPhyt = 0.d0
    LimnoParam%wCGrzFiJvPhyt = 0.d0
    LimnoParam%wNGrzFiJvPhyt = 0.d0
    LimnoParam%wPGrzFiJvPhyt = 0.d0
    !tDGrzFishBac = 0.d0 
    !tCGrzFishBac = 0.d0
    !tNGrzFishBac = 0.d0
    !tPGrzFishBac = 0.d0
    !tDGrzFishBent = 0.d0
    !tCGrzFishBent = 0.d0
    !tNGrzFishBent = 0.d0
    !tPGrzFishBent = 0.d0
    LimnoParam%wDGrzFiAdPom = 0.d0 
    LimnoParam%wCGrzFiAdPom = 0.d0 
    LimnoParam%wNGrzFiAdPom = 0.d0
    LimnoParam%wPGrzFiAdPom = 0.d0
    LimnoParam%wDGrzFiJvPom = 0.d0 
    LimnoParam%wCGrzFiJvPom = 0.d0 
    LimnoParam%wNGrzFiJvPom = 0.d0
    LimnoParam%wPGrzFiJvPom = 0.d0
    LimnoParam%wDHarvFiAd = 0.d0     
    LimnoParam%wCHarvFiAd = 0.d0     
    LimnoParam%wNHarvFiAd = 0.d0
    LimnoParam%wPHarvFiAd = 0.d0
    LimnoParam%wDHarvFiJv = 0.d0     
    LimnoParam%wCHarvFiJv = 0.d0     
    LimnoParam%wNHarvFiJv = 0.d0
    LimnoParam%wPHarvFiJv = 0.d0
    LimnoParam%wDRespFiAd = 0.d0     
    LimnoParam%wDRespFiJv = 0.d0     
    LimnoParam%wCRespFiAd = 0.d0
    LimnoParam%wCRespFiJv = 0.d0
    LimnoParam%wDEgesFiAd = 0.d0     
    LimnoParam%wCEgesFiAd = 0.d0     
    LimnoParam%wNEgesFiAd = 0.d0
    LimnoParam%wPEgesFiAd = 0.d0
    LimnoParam%wDEgesFiJv = 0.d0     
    LimnoParam%wCEgesFiJv = 0.d0     
    LimnoParam%wNEgesFiJv = 0.d0
    LimnoParam%wPEgesFiJv = 0.d0
    LimnoParam%wDMortFiAd = 0.d0     
    LimnoParam%wCMortFiAd = 0.d0     
    LimnoParam%wNMortFiAd = 0.d0
    LimnoParam%wPMortFiAd = 0.d0
    LimnoParam%wDMortFiJv = 0.d0     
    LimnoParam%wCMortFiJv = 0.d0     
    LimnoParam%wNMortFiJv = 0.d0
    LimnoParam%wPMortFiJv = 0.d0
    LimnoParam%wNExcrFiAd = 0.d0
    LimnoParam%wPExcrFiAd = 0.d0
    LimnoParam%wNExcrFiJv = 0.d0
    LimnoParam%wPExcrFiJv = 0.d0
    LimnoParam%wDMessFiAdOM = 0.d0
    LimnoParam%wCMessFiAdOM = 0.d0
    LimnoParam%wNMessFiAdOM = 0.d0
    LimnoParam%wPMessFiAdOM = 0.d0
    LimnoParam%wDMessFiJvOM = 0.d0
    LimnoParam%wCMessFiJvOM = 0.d0
    LimnoParam%wNMessFiJvOM = 0.d0
    LimnoParam%wPMessFiJvOM = 0.d0    
!------------------------------------------------------------------------------------------------------------------!	
!  A B I O T I C O                                                                                                 !	
!------------------------------------------------------------------------------------------------------------------!	
!------------------------------------------------------------------------------------------------------------------!	
!  Z O N A    P E L A G I C A                                                                                      !	
!------------------------------------------------------------------------------------------------------------------!	
!1)Materia Organica
 	!PARTICULATE ORGANIC MATTER (POM)
    If (LimnoParam%InclMatOrg == 1) Then
 	    Do pom=1,LimnoParam%npom
	        LimnoParam%sDPomW(:,:,pom) = LimnoParam%sDPomW0(pom)
            LimnoParam%sCPomW(:,:,pom) = LimnoParam%cCDPomRef(pom)*LimnoParam%sDPomW(:,:,pom)			                        !inicial N no humus (gN/m²)
            LimnoParam%sNPomW(:,:,pom) = LimnoParam%cNDPomRef(pom)*LimnoParam%sDPomW(:,:,pom)
            LimnoParam%sPPomW(:,:,pom) = LimnoParam%cPDPomRef(pom)*LimnoParam%sDPomW(:,:,pom)
	
	        LimnoParam%rCDPomW(:,:,pom) = LimnoParam%cCDPomRef(pom)				                                    !fração inicial de P no humus (gP/gD)
	        LimnoParam%rNDPomW(:,:,pom) = LimnoParam%cNDPomRef(pom)				                                    !fração inicial de P no humus (gP/gD)
	        LimnoParam%rPDPomW(:,:,pom) = LimnoParam%cPDPomRef(pom)				                                    !fração inicial de N no humus (gN/gD)
        Enddo
        
 	    Do dom=1,LimnoParam%ndom
            If (LimnoParam%InclMatOrgSplit == 1) Then
                LimnoParam%sDDomW(:,:,dom) = LimnoParam%sDDomW0(dom)           						                    !total de humos inicial no topo de camada (gD/m²)
                LimnoParam%sCDomW(:,:,dom) = LimnoParam%cCDDomRef(dom)*LimnoParam%sDDomW(:,:,dom)
                LimnoParam%sNDomW(:,:,dom) = LimnoParam%cNDDomRef(dom)*LimnoParam%sDDomW(:,:,dom)
                LimnoParam%sPDomW(:,:,dom) = LimnoParam%cPDDomRef(dom)*LimnoParam%sDDomW(:,:,dom)
                !LimnoParam%sSiDomW(:,:,dom)= LimnoParam%cSiDDomRef(dom)*LimnoParam%sDDomW(:,:,dom)
	
	            LimnoParam%rCDDomW(:,:,dom) = LimnoParam%cCDDomRef(dom)				                                    !fração inicial de P no humus (gP/gD)
	            LimnoParam%rNDDomW(:,:,dom) = LimnoParam%cNDDomRef(dom)				                                    !fração inicial de P no humus (gP/gD)
                LimnoParam%rPDDomW(:,:,dom) = LimnoParam%cPDDomRef(dom)				                                    !fração inicial de N no humus (gN/gD)
                LimnoParam%rSiDDomW(:,:,dom)= LimnoParam%cSiDDomRef(dom)
            Else
                LimnoParam%sDDomW(:,:,dom) = LimnoParam%sDPomW0(dom)           						                    !total de humos inicial no topo de camada (gD/m²)
                LimnoParam%sCDomW(:,:,dom) = LimnoParam%cCDPomRef(dom)*LimnoParam%sDPomW(:,:,dom)
                LimnoParam%sNDomW(:,:,dom) = LimnoParam%cNDPomRef(dom)*LimnoParam%sDPomW(:,:,dom)
                LimnoParam%sPDomW(:,:,dom) = LimnoParam%cPDPomRef(dom)*LimnoParam%sDPomW(:,:,dom)
                !LimnoParam%sSiDomW(:,:,dom)= LimnoParam%cSiDPomRef(dom)*LimnoParam%sDPomW(:,:,dom)
	
	            LimnoParam%rCDDomW(:,:,dom) = LimnoParam%cCDPomRef(dom)				                                    !fração inicial de P no humus (gP/gD)
	            LimnoParam%rNDDomW(:,:,dom) = LimnoParam%cNDPomRef(dom)				                                    !fração inicial de P no humus (gP/gD)
                LimnoParam%rPDDomW(:,:,dom) = LimnoParam%cPDPomRef(dom)				                                    !fração inicial de N no humus (gN/gD)
                !LimnoParam%rSiDDomW(:,:,dom)= LimnoParam%cSiDPomRef(dom)
            EndIf
        EndDo        
    Else
        LimnoParam%sDPomW = 0.d0
        LimnoParam%sCPomW = 0.d0
        LimnoParam%sNPomW = 0.d0
        LimnoParam%sPPomW = 0.d0
        LimnoParam%rCDPomW = 0.d0
        LimnoParam%rNDPomW = 0.d0
        LimnoParam%rPDPomW = 0.d0
        LimnoParam%sDDomW = 0.d0
        LimnoParam%sCDomW = 0.d0
        LimnoParam%sNDomW = 0.d0
        LimnoParam%sPDomW = 0.d0
        LimnoParam%rCDDomW = 0.d0
        LimnoParam%rNDDomW = 0.d0
        LimnoParam%rPDDomW = 0.d0
    EndIf
    
 	!Processes
 	!Hidrolise
    LimnoParam%wDHidPomW = 0.d0
	LimnoParam%wCHidPomW = 0.d0
	LimnoParam%wNHidPomW = 0.d0
    LimnoParam%wPHidPomW = 0.d0
    LimnoParam%aTLimHidPom = 0.d0
    LimnoParam%aO2LimHidPom = 0.d0
    LimnoParam%aBacLimHidPom = 0.d0
    LimnoParam%kCorHidPom = 0.d0
    

	!Processes
    
	!Infiltração
	LimnoParam%tDInfDomW = 0.d0
	LimnoParam%tCInfDomW = 0.d0
	LimnoParam%tNInfDomW = 0.d0
	LimnoParam%tPInfDomW = 0.d0
	!Mineração - IF InclBac=0
	LimnoParam%wDMinAerW = 0.d0
	LimnoParam%wCMinAerW = 0.d0
	LimnoParam%wNMinAerW = 0.d0
	LimnoParam%wPMinAerW = 0.d0
	LimnoParam%wSiMinAerW = 0.d0
    LimnoParam%wCMinAer2CO2W = 0.d0
    LimnoParam%wCMinAer2CH4W = 0.d0
	LimnoParam%wDMinAnaerW = 0.d0
	LimnoParam%wCMinAnaerW = 0.d0
	LimnoParam%wNMinAnaerW = 0.d0
	LimnoParam%wPMinAnaerW = 0.d0
	LimnoParam%wSiMinAnaerW = 0.d0
    
    !  
!2)Inorganci components
    If (LimnoParam%InclMatInorg == 1) Then
	    LimnoParam%sDIMW=LimnoParam%sDIMW0		        !Materia inorganica na água (mgDW/l)
        LimnoParam%sPAIMW=LimnoParam%sPAIMW0		        !fosforo absorvido na matéria inorganica (mgP/l)
    Else
        LimnoParam%sDIMW=0.
        LimnoParam%sPAIMW = 0.d0
    EndIf
	LimnoParam%sPO4W=LimnoParam%sPO4W0		        !Ortofosfato da água (mgP/l)
	
	LimnoParam%sNO3W=LimnoParam%sNO3W0		        !NO3 na água (mgN/l)
	LimnoParam%sNH4W=LimnoParam%sNH4W0		        !NH4 na água (mgN/l)
	LimnoParam%sSiO2W=LimnoParam%sSiO2W0		        !silica dissolvida da água (mgSi/l)
	LimnoParam%sCFW=LimnoParam%sCFW0
	
	
	!Processes
	!Desnitrificação
    LimnoParam%wNDenitW = 0.d0
	!Nitrificação
    LimnoParam%wNNitrW = 0.d0
	!Adsorção do Fósforo
    LimnoParam%wPSorpIMW = 0.d0
    LimnoParam%tPSorpIMS = 0.d0
!3)Gases
    If (LimnoParam%iOxyBOD==1) Then
	    LimnoParam%sO2W=LimnoParam%sO2W0			        !oxigenio na água (mgO2/l)
        LimnoParam%sDBOW=LimnoParam%sDBOW0
    Else
        LimnoParam%sO2W=10.d0
        LimnoParam%sDBOW=0.d0
    EndIf
    
    If (LimnoParam%iCarbon==1) Then
        LimnoParam%sDicW=LimnoParam%sDicW0
        LimnoParam%sCH4W=LimnoParam%sCH4W0
    Else
        LimnoParam%sDicW =0.d0
        LimnoParam%sCH4W = 0.d0
    EndIf
    
	!ZERAR FLUXOS PARA TEMPODIA=TEMPODIA0
	!Difusao Agua-Ar
    LimnoParam%tO2Reaer(:) = 0.d0
    LimnoParam%tCO2Reaer(:) = 0.d0

!4)Variaveis Secundarias
    LimnoParam%spHW=LimnoParam%spHW0
    LimnoParam%sAlkW=LimnoParam%sAlkW0
	LimnoParam%sGPPW=0.d0
	LimnoParam%sRespW=0.d0
	LimnoParam%sNEPW=0.d0

!5)Carbonate System    
    Hion =  10.0 **(-LimnoParam%spHW0)
    aAlkW = LimnoParam%sAlkW0/50000. !mg/L -> eq/L
    
    k1=10**-6.3
    k2=10**-10.3
    kwater=10**-14.

    a0 = Hion*Hion/(Hion*Hion + k1*Hion + k1*k2)
    a1 = k1*Hion / (Hion*Hion + k1*Hion + k1*k2)
    a2 = k1*k2   / (Hion*Hion + k1*Hion + k1*k2)  
        
    LimnoParam%sDicW = ((aAlkW - kwater/Hion + Hion) / (a1 + 2*a2))*12000. !mg/L
    
    LimnoParam%sH2CO3W  = LimnoParam%sDicW  *a0  !mg/L
    LimnoParam%sHCO3W   = LimnoParam%sDicW * a1  !mg/L
    LimnoParam%sCO3W    = LimnoParam%sDicW * a2  !mg/L	 
	
    
    !  Day_length
    Call unix2c(simParam%time, idate)
    !print *,idate
    !Day0*19                                        ! Initial time of simulation (AAAA/MM/DD HH:MM:SS)
	CALL JULIANDAY(floor(idate(1)),1,1,diajulcont1) !January 1st if the current year
    
    CALL JULIANDAY(floor(idate(1)),floor(idate(2)),floor(idate(3)),dayoftheyear)
    simParam%Simday = dayoftheyear
	dayoftheyear = int(dayoftheyear-diajulcont1)+1
                
	IF (dayoftheyear<1.OR.dayoftheyear>366) THEN
        write(*,*) ' ERROR: Day of year is out of bounds.'
        write(*,*) ' dayoftheyear =',dayoftheyear
        STOP 
    ENDIF
    simParam%Julday = dayoftheyear
    
				
	delta = 23.45d0*(HydroParam%PI/180)*COS(2*HydroParam%PI/365*(172-dayoftheyear))
	cos_hss = -tan(HydroParam%LAT*HydroParam%PI/180.)*tan(delta)
                
	IF (abs(cos_hss) > 1.0) THEN
        write(*,*) ' ERROR: cos_hss (cosine of the sunset hour angle) is out of range.'
        write(*,*) ' cos_hss =',cos_hss
        STOP 
    ENDIF
					    
	sunSetHourAngle = acos(cos_hss)
	dayPhotoFration = sunSetHourAngle/HydroParam%PI
	LimnoParam%ufDay = dayPhotoFration  !cfDayAve-cfDayVar*(-cos(2*PI*(DiaJul+TempLag)/365.)) !Fotoperíodo (comprimento do dia) (h/24h)
    
!------------------------------------------------------------------------------------------------------------------!	
!  S E D I M E N T O                                                                                               !	
!------------------------------------------------------------------------------------------------------------------!	

   
    !1)Matriz do Sedimento Topo
	LimnoParam%bRhoSolidS0=LimnoParam%fDOrgS0*LimnoParam%cRhoOM+(1.0-LimnoParam%fDOrgS0)*LimnoParam%cRhoIM										  !Densidade inicial do material sólido (g/m³ de sólido)
	LimnoParam%bPorS=(1.-LimnoParam%fDTotS0)*(LimnoParam%bRhoSolidS0/1000000.)/(LimnoParam%fDTotS0+(1.-LimnoParam%fDTotS0)*(LimnoParam%bRhoSolidS0/1000000.)) !porosidade (m³ água/m³ de sedimento)
	LimnoParam%bPorCorS=LimnoParam%bPorS**(1.0+LimnoParam%bPorS)															  !Porosidade do sedimento corrigida pela tortuosidade
	LimnoParam%bRhoTotS0=LimnoParam%bRhoSolidS0*(1.0-LimnoParam%bPorS)													  !Densidade aparente do sedimento (g solido/m³ de sedimento)
	LimnoParam%bDTotS0=LimnoParam%bRhoTotS0*LimnoParam%cDepthS															  !total de peso seco inicial no topo de camada (gD/m²)
    
    LimnoParam%tDResusTot=0.

!2)Materia Organica
 	!PARTICULATE ORGANIC MATTER (POM)
    If (LimnoParam%InclMatOrg == 1) Then
 	    Do pom=1,LimnoParam%npom
	        If(pom==LimnoParam%LABIL)Then
	            LimnoParam%sDPomS0(pom) = LimnoParam%fDLabilS0*LimnoParam%fDOrgS0*LimnoParam%bDTotS0    						 !total de detritos inicial no topo de camada (gD/m²)
            Else
                LimnoParam%sDPomS0(pom) = (1.0-LimnoParam%fDLabilS0)*LimnoParam%fDOrgS0*LimnoParam%bDTotS0						  !total de humos inicial no topo de camada (gD/m²)
            EndIf
        
            LimnoParam%sDPomS(:,pom) = LimnoParam%sDPomS0(pom)
	        LimnoParam%sCPomS(:,pom) = LimnoParam%cCDPomRef(pom)*LimnoParam%sDPomS(:,pom)				         !inicial N no humus (gN/m²)
            LimnoParam%sNPomS(:,pom) = LimnoParam%cNDPomRef(pom)*LimnoParam%sDPomS(:,pom)
	        LimnoParam%sPPomS(:,pom) = LimnoParam%cPDPomRef(pom)*LimnoParam%sDPomS(:,pom)
	
	        LimnoParam%rCDPomS(:,pom) = LimnoParam%cCDPomRef(pom)				                     !fração inicial de P no humus (gP/gD)
	        LimnoParam%rNDPomS(:,pom) = LimnoParam%cNDPomRef(pom)				                     !fração inicial de P no humus (gP/gD)
	        LimnoParam%rPDPomS(:,pom) = LimnoParam%cPDPomRef(pom)				                     !fração inicial de N no humus (gN/gD)
        EndDo	
 	    Do dom=1,LimnoParam%ndom
            If (LimnoParam%InclMatOrgSplit == 1) Then
                LimnoParam%sDDomS(:,dom) = LimnoParam%sDDomS0(dom)           						     !total de humos inicial no topo de camada (gD/m²)
   	            LimnoParam%sCDomS(:,dom) = LimnoParam%cCDDomRef(dom)*LimnoParam%sDDomS(:,dom)			         !inicial N no humus (gN/m²)
                LimnoParam%sNDomS(:,dom) = LimnoParam%cNDDomRef(dom)*LimnoParam%sDDomS(:,dom)
                LimnoParam%sPDomS(:,dom) = LimnoParam%cPDDomRef(dom)*LimnoParam%sDDomS(:,dom)
                !LimnoParam%sSiDomS(:,dom)= LimnoParam%cSiDDomRef(dom)*LimnoParam%sDDomS(:,dom)

	            LimnoParam%rCDDomS(:,dom) = LimnoParam%cCDDomRef(dom)				                     !fração inicial de P no humus (gP/gD)
	            LimnoParam%rNDDomS(:,dom) = LimnoParam%cNDDomRef(dom)				                     !fração inicial de P no humus (gP/gD)
	            LimnoParam%rPDDomS(:,dom) = LimnoParam%cPDDomRef(dom)				                     !fração inicial de N no humus (gN/gD)
	            !LimnoParam%rSiDDomS(:,dom)= LimnoParam%cSiDDomRef(dom)				                     !fração inicial de N no humus (gN/gD)
            Else
                LimnoParam%sDDomS(:,dom) = LimnoParam%sDPomS0(dom)           						     !total de humos inicial no topo de camada (gD/m²)
   	            LimnoParam%sCDomS(:,dom) = LimnoParam%cCDPomRef(dom)*LimnoParam%sDPomS(:,dom)			         !inicial N no humus (gN/m²)
                LimnoParam%sNDomS(:,dom) = LimnoParam%cNDPomRef(dom)*LimnoParam%sDPomS(:,dom)
                LimnoParam%sPDomS(:,dom) = LimnoParam%cPDPomRef(dom)*LimnoParam%sDPomS(:,dom)
                !LimnoParam%sSiDomS(:,dom)= LimnoParam%cSiDPomRef(dom)*LimnoParam%sDPomS(:,dom)

	            LimnoParam%rCDDomS(:,dom) = LimnoParam%cCDPomRef(dom)				                     !fração inicial de P no humus (gP/gD)
	            LimnoParam%rNDDomS(:,dom) = LimnoParam%cNDPomRef(dom)				                     !fração inicial de P no humus (gP/gD)
	            LimnoParam%rPDDomS(:,dom) = LimnoParam%cPDPomRef(dom)				                     !fração inicial de N no humus (gN/gD)
	            !LimnoParam%rSiDDomS(:,dom)= LimnoParam%cSiDPomRef(dom)				                     !fração inicial de N no humus (gN/gD)
            EndIf
        ENDDO	
        
    Else
        LimnoParam%sDPomS = 0.d0
        LimnoParam%sCPomS = 0.d0
        LimnoParam%sNPomS = 0.d0
        LimnoParam%sPPomS = 0.d0
        LimnoParam%rCDPomS = 0.d0
        LimnoParam%rNDPomS = 0.d0
        LimnoParam%rPDPomS = 0.d0
        LimnoParam%sDDomS = 0.d0
        LimnoParam%sCDomS = 0.d0
        LimnoParam%sNDomS = 0.d0
        LimnoParam%sPDomS = 0.d0
        LimnoParam%rCDDomS = 0.d0
        LimnoParam%rNDDomS = 0.d0
        LimnoParam%rPDDomS = 0.d0
    EndIf
    
	!ZERAR FLUXOS PARA TEMPODIA=TEMPODIA0
	!Hidrolise
    LimnoParam%tDHidPomS = 0.d0
    LimnoParam%tCHidPomS = 0.d0
    LimnoParam%tNHidPomS = 0.d0
    LimnoParam%tPHidPomS = 0.d0	
	!Erosão
    LimnoParam%uDErosIMS = 0.d0
    LimnoParam%uDErosIMW = 0.d0
    LimnoParam%tDErosPom = 0.d0
    LimnoParam%tCErosPom = 0.d0
    LimnoParam%tNErosPom = 0.d0
    LimnoParam%tPErosPom = 0.d0
	!Enterro
    LimnoParam%tDBurPom = 0.d0
    LimnoParam%tCBurPom = 0.d0
    LimnoParam%tNBurPom = 0.d0
    LimnoParam%tPBurPom = 0.d0
	!Resuspensao
    LimnoParam%tDResusPom = 0.d0
    LimnoParam%tCResusPom = 0.d0
    LimnoParam%tPResusPom = 0.d0
    LimnoParam%tNResusPom = 0.d0
	!Sedimentação
    LimnoParam%tDSetPom = 0.d0
    LimnoParam%tCSetPom = 0.d0
    LimnoParam%tNSetPom = 0.d0
    LimnoParam%tPSetPom = 0.d0! 	
 	!DISSOLVED ORGANIC MATTER (DOM)
    
	!Processes
	!Infiltração
    LimnoParam%tDInfDomS = 0.d0
    LimnoParam%tCInfDomS = 0.d0
    LimnoParam%tNInfDomS = 0.d0
    LimnoParam%tPInfDomS = 0.d0
	!Mineração - IF InclBac=0
    LimnoParam%tDMinAerS = 0.d0
    LimnoParam%tCMinAerS = 0.d0
    LimnoParam%tNMinAerS = 0.d0
    LimnoParam%tPMinAerS = 0.d0    
    LimnoParam%tCMinAer2CO2S = 0.d0
    LimnoParam%tCMinAer2CH4S = 0.d0
    LimnoParam%tO2MinDetS = 0.d0
        
    LimnoParam%tDMinAnaerS = 0.d0
    LimnoParam%tCMinAnaerS = 0.d0
    LimnoParam%tNMinAnaerS = 0.d0
    LimnoParam%tPMinAnaerS = 0.d0    
    LimnoParam%tCMinAnaer2CO2S = 0.d0
    LimnoParam%tCMinAnaer2CH4S = 0.d0
        
	!Difusao Agua-Sed
    LimnoParam%tDDifDom = 0.d0
    LimnoParam%tCDifDom = 0.d0
    LimnoParam%tNDifDom = 0.d0
    LimnoParam%tPDifDom = 0.d0
	!Resuspensao
    LimnoParam%tDResusDom = 0.d0
    LimnoParam%tCResusDom = 0.d0
    LimnoParam%tNResusDom = 0.d0
    LimnoParam%tPResusDom = 0.d0
	!Enterro
    LimnoParam%tDBurDom = 0.d0
    LimnoParam%tCBurDom = 0.d0
    LimnoParam%tNBurDom = 0.d0
    LimnoParam%tPBurDom = 0.d0
        
!3)Materia Inorganica
    If (LimnoParam%InclMatInorg == 1) Then
	    LimnoParam%sDIMS0=LimnoParam%bDTotS0-SUM(LimnoParam%sDPomS0(:)) 												 !total de materia inorganica inicial no topo de camada (gD/m²)
	    LimnoParam%sDIMS=LimnoParam%sDIMS0
        LimnoParam%sPAIMS0=LimnoParam%fPAdsS0*LimnoParam%fPInorgS0*LimnoParam%bDTotS0											 !P absorvido na materia inorganica no sedimento (gP/m²)
        LimnoParam%sPAIMS=LimnoParam%sPAIMS0
        LimnoParam%rPDIMS=LimnoParam%sPAIMS0/LimnoParam%sDIMS0
    Else
        LimnoParam%sDIMS=0.
        LimnoParam%sPAIMS=0.
        LimnoParam%rPDIMS=0.
    EndIf
    
    If (LimnoParam%iCarbon==1) Then
        LimnoParam%sDicS=LimnoParam%sDicS0			                                                         !NH4 no sedimento (gN/m²)    
        LimnoParam%oDicS=LimnoParam%sDicS0/LimnoParam%cDepthS/LimnoParam%bPorS
        LimnoParam%sCH4S=LimnoParam%sCH4S0
    Else
        LimnoParam%sDicS = 0.d0
        LimnoParam%oDicS = 0.d0
        LimnoParam%sCH4S = 0.d0
    EndIf
    
        
    LimnoParam%sPO4S0=(1.0-LimnoParam%fPAdsS0)*LimnoParam%fPInorgS0*LimnoParam%bDTotS0										 !P dissolvido no sedimento (gP/m²)
	LimnoParam%sPO4S=LimnoParam%sPO4S0
	LimnoParam%sNO3S=LimnoParam%sNO3S0
	LimnoParam%sNH4S=LimnoParam%sNH4S0
    
	LimnoParam%oNO3S=LimnoParam%sNO3S0/LimnoParam%cDepthS/LimnoParam%bPorS	                                                 !conc. de N-NO3 dissolvido na água intersticial (gN/m³)
	LimnoParam%oNH4S=LimnoParam%sNH4S0/LimnoParam%cDepthS/LimnoParam%bPorS	                                                 !conc. de N-NH4 dissolvido na água intersticial (gN/m³)
	LimnoParam%oPO4S=LimnoParam%sPO4S0/LimnoParam%cDepthS/LimnoParam%bPorS	                                                 !conc. de P-PO4 dissolvido na água intersticial (gP/m³)
    
	LimnoParam%aDTotS=LimnoParam%sDIMS0+SUM(LimnoParam%sDPomS0(:))	
	!ZERAR FLUXOS PARA TEMPODIA=TEMPODIA0
	!Ressuspencao
	LimnoParam%tPResusPO4 = 0.d0
	LimnoParam%tNResusNH4 = 0.d0 
	LimnoParam%tNResusNO3 = 0.d0 
    LimnoParam%tCResusDic = 0.d0 
    LimnoParam%tDResusIM  = 0.d0
    LimnoParam%tPResusAIM = 0.d0
    !Sedimentação
    LimnoParam%tDSetIM = 0.d0
    LimnoParam%tPSetAIM = 0.d0
	!Desnitrificação
    LimnoParam%tNDenitS = 0.d0
	!Nitrificação
    LimnoParam%tNNitrS = 0.d0
    LimnoParam%tO2NitrS = 0.d0
	!Adsorção do Fósforo
    LimnoParam%tPSorpIMS = 0.d0
	
!4)Gases
    LimnoParam%sO2S = LimnoParam%sO2W0/2.
    

!5)Variaveis Secundarias
    LimnoParam%spHS=LimnoParam%spHS0
    LimnoParam%sAlkS=LimnoParam%sAlkS0
    
!4)Carbonate System 
    Hion =  10.0 **(-LimnoParam%spHS0)
    aAlkS = LimnoParam%sAlkS0/50000. !mg/L -> eq/L
    
    k1=10**-6.3
    k2=10**-10.3
    kwater=10**-14.

    a0 = Hion*Hion/(Hion*Hion + k1*Hion + k1*k2)
    a1 = k1*Hion / (Hion*Hion + k1*Hion + k1*k2)
    a2 = k1*k2   / (Hion*Hion + k1*Hion + k1*k2)  
        
    LimnoParam%sDicS = (aAlkS - kwater/Hion + Hion) / (a1 + 2*a2)*12000. !mol/L
    
    LimnoParam%sH2CO3S = LimnoParam%sDicS  *a0 
    LimnoParam%sHCO3S  = LimnoParam%sDicS * a1 
    LimnoParam%sCO3S   = LimnoParam%sDicS * a2    
    
    !3.Sediment
    !3.1 Diffusion in the sediment
	LimnoParam%tPDifPO4 = 0.d0 
    LimnoParam%tNDifNH4 = 0.d0
    LimnoParam%tNDifNO3 = 0.d0
    LimnoParam%tCDifDic = 0.d0
    LimnoParam%tO2Dif = 0.d0
	LimnoParam%tDDifDom = 0.d0
    LimnoParam%tCDifDom = 0.d0
    LimnoParam%tNDifDom = 0.d0
    LimnoParam%tPDifDom = 0.d0
    
    !Autorecovery
    If (simParam%it > 0) Then
        Do i = 1, simParam%nOutput
            text = trim(simParam%wqoOutputParameters(i)%name)
            If (trim(text)=='sDTempW') Then
                LimnoParam%sDTempW = simParam%wqosave(:,:,i)
            EndIf
            If (trim(text)=='sDSal') Then
                LimnoParam%sDSal = simParam%wqosave(:,:,i)
            EndIf
        EndDo  
    EndIf
    
	 !!Water Density as a Function of Water Temperature and Salinity
    
    Call WaterDensity(HydroParam,MeshParam,LimnoParam)
    
    Call UpdateLimnVars(HydroParam,MeshParam,LimnoParam,simParam)
    
    Call UpdateSedVars(HydroParam,MeshParam,LimnoParam)
    
End Subroutine ReadWQIniCond
    