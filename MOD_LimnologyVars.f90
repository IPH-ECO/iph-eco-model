Module LimnologyVars
    
    ! Declare the Water Quality Variables
    
    ! List of Modifications: 
    !   -> 30.12.2014: Routine Implementation      (Rafael Cavalcanti)
    ! Programmer: Rafael Cavalcanti
    use domain_types
    
    Implicit None
    
    type LimnologyParam
        !General 
        Integer:: iLimn !< Limnology module flag: iLimn = 0 (off); iLimn = 1 (on)
        Integer:: iSed !< Sediment module flag: iSed = 0 (off); iSed = 1 (on)
        !Integer,ALLOCATABLE:: IndexWQ(:,:) !< Cell Index of Water Quality boundary condition (IndexWQ(:,1) -> Water Temperature; IndexWQ(:,2) -> Salinity; IndexWQ(:,3) -> Inorganic Matter; IndexWQ(:,4) -> Dissolved Organic Matter; IndexWQ(:,4) -> Particulated Organic Matter)
        Real, Allocatable:: dVarEstS(:,:)

    
        !1. Water Temperature
        !1.1 Flags
        Integer:: iTempW !< Water temperatue module flag: iTempW = 0 (off); iTempW = 1 (on)
        !1.2 State variables     
        Real,ALLOCATABLE:: sDTempW(:,:) !<Water Temperature in the current time
        Real,ALLOCATABLE:: sDTempWP(:,:) !<Water Temperature in the past time
        !1.3 Initial condition    
        Real sDTempW0 !<Initial condition of Water Temperature
        !1.4 Boundary condition 
        Real,ALLOCATABLE:: uDLoadTemp(:,:) !< Water temperature loading value 
        Integer:: NuDLoadTemp !< Number of boundary conditions of Water temperature
        Integer,ALLOCATABLE:: TempnTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: TempTime(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: TempValue(:,:) !< value of each time-serie
        Integer,ALLOCATABLE:: TempSmallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: TempCapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NTempIndex(:,:)
        
        !1.5 Parameters    
        Real:: a_param                    !< Water Temperature Model Coefficient (0.5 to 0.7)
        Real:: tau_param                  !< Stefan-Boltzmann Constant = 11.7 E-8 cal/cm²/d/K
        Real:: e_param                    !< Emissivity of the Radiating Body
        Real:: c1_param                   !< Bowen's Coefficient (0.47 mmHg/ºC)
        Real:: cd_param                   !< Specific Heat (cal/g/ºC)
        Real:: WtempRef
    
        !2. Salinity
        !2.1 Flags
        Integer:: iSal !< Salinity module flag: iSal = 0 (off); iSal = 1 (on)
        !2.2 State variables     
	    Real,ALLOCATABLE:: sDSal(:,:) !<Salinity in the current time
        Real,ALLOCATABLE:: sDSalP(:,:) !<Salinity in the past time
        !2.3 Initial condition    
	    Real sDSal0 !<Initial condition of Salinity
        !2.4 Boundary condition 
        Real,ALLOCATABLE:: uDLoadSal(:,:) !< Salinity loading value 
        Integer:: NuDLoadSal !< Number of cells with boundary condition of Salinity
        Integer,ALLOCATABLE:: SalnTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: SalTime(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: SalValue(:,:) !< value of each time-serie
        Integer,ALLOCATABLE:: SalSmallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: SalCapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NSalIndex(:,:)
        

        Integer:: iLimno !< Water temperatue module flag: iLimno = 0 (off); iLimno = 1 (on
        !3. Inorganic Matter
        !3.1 Flags
        Integer:: InclMatInorg !< Inorganic Matter module flag: InclMatInorg = 0 (off); InclMatInorg = 1 (on)
        Integer:: InclMatInorgSplit !< Split Inorganic Matter flag: InclMatInorgSplit = 0 (no); InclMatInorgSplit = 1 (yes)
        !3.2 State variables  
	    Real,ALLOCATABLE:: sDIMW(:,:) !<Inorganic Matter in the water at current time
        Real,ALLOCATABLE:: sDIMWP(:,:) !<Inorganic Matter in the water at past time
	    Real,ALLOCATABLE:: sDIMS(:) !<Inorganic Matter in the sediment at current time
        Real,ALLOCATABLE:: sDIMSP(:) !<Inorganic Matter in the sediment at past time
        !3.3 Initial condition    
	    Real sDIMW0    !<Initial condition of Inorganic Matter in the water (mgDW/l)
        Real sDIMS0    !<Initial condition of Inorganic Matter in the sediment (gD/m²)
        !3.4 Boundary condition 
        Real,ALLOCATABLE:: uDLoadIM(:,:) !< Inorganic Matter loading value 
        Integer:: NuDLoadIM !< Number of cells with boundary condition of Inorganic Matter
        Integer,ALLOCATABLE:: IMnTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: IMTime(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: IMValue(:,:) !< value of each time-serie
        Integer,ALLOCATABLE:: IMSmallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: IMCapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NIMIndex(:,:)
   
        !4/5. Dissolved and Particulted Organic Matter
        !4/5.1 Flags
        Integer:: InclMatOrg !< Organic Matter flag: InclMatOrg = 0 (no); InclMatOrg = 1 (yes)
        Integer:: InclMatOrgSplit !< Split Organic Matter flag: InclMatOrgSplit = 0 (no); InclMatOrgSplit = 1 (yes)
        Integer:: InclPomDomSplit !< Include Particulated Organic Matter flag: InclPomSplit = 0 (no); InclPomSplit = 1 (yes)
        !4/5.2 State variables     
	    Real,ALLOCATABLE:: sDPomW(:,:,:), sCPomW(:,:,:), sNPomW(:,:,:), sPPomW(:,:,:), sSiPomW(:,:,:) !<Fractions of Particulated Organic Matter in the water at current time
	    Real,ALLOCATABLE:: sDDomW(:,:,:), sCDomW(:,:,:), sNDomW(:,:,:), sPDomW(:,:,:), sSiDomW(:,:,:) !<Fractions of Dissolved Organic Matter in the water at current time
	    Real,ALLOCATABLE:: sDPomS(:,:),   sCPomS(:,:),   sNPomS(:,:),   sPPomS(:,:),   sSiPomS(:,:) !<Fractions of Particulated Organic Matter in the sediment at current time
	    Real,ALLOCATABLE:: sDDomS(:,:),   sCDomS(:,:),   sNDomS(:,:),   sPDomS(:,:),   sSiDomS(:,:) !<Fractions of Dissolved Organic Matter in the sediment at current time
	    Real,ALLOCATABLE:: sDPomWP(:,:,:), sCPomWP(:,:,:), sNPomWP(:,:,:), sPPomWP(:,:,:), sSiPomWP(:,:,:) !<Fractions of Particulated Organic Matter in the water at past time
	    Real,ALLOCATABLE:: sDDomWP(:,:,:), sCDomWP(:,:,:), sNDomWP(:,:,:), sPDomWP(:,:,:), sSiDomWP(:,:,:) !<Fractions of Dissolved Organic Matter in the water at past time
	    Real,ALLOCATABLE:: sDPomSP(:,:),   sCPomSP(:,:),   sNPomSP(:,:),   sPPomSP(:,:),   sSiPomSP(:,:) !<Fractions of Particulated Organic Matter in the sediment at past time
	    Real,ALLOCATABLE:: sDDomSP(:,:),   sCDomSP(:,:),   sNDomSP(:,:),   sPDomSP(:,:),   sSiDomSP(:,:) !<Fractions of Dissolved Organic Matter in the sediment at past time
        !4/5.3 Initial condition    
	    Real,ALLOCATABLE:: sDDomW0(:),sDDomS0(:)        !<Initial condition of Dissolved Organic Matter in the water and in the sediment (mgDW/l)
	    Real,ALLOCATABLE:: sDPomW0(:),sDPomS0(:)        !<Initial condition of Particulated Organic Matter in the water and in the sediment (gD/m²)
        !4/5.4 Boundary condition 
        Real,ALLOCATABLE:: uDLoadDom(:,:) !< Dissolved Organic Matter loading value
        Real,ALLOCATABLE:: uDLoadPom(:,:) !< Particulated Organic Matter loading value
        Integer:: NuDLoadDom !< Number of cells with boundary condition of Dissolved Organic Matter
        Integer,ALLOCATABLE:: DomnTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: DomTime(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: DomValue(:,:) !< value of each time-serie
        Integer:: NuDLoadPom !< Number of cells with boundary condition of Particulated Organic Matter
        Integer,ALLOCATABLE:: PomnTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: PomTime(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: PomValue(:,:) !< value of each time-serie
        Integer,ALLOCATABLE:: DomSmallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: DomCapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NDomIndex(:,:)
        Integer,ALLOCATABLE:: PomSmallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: PomCapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NPomIndex(:,:)
    
        !6. PO4
        !6.1 State variables  
	    Real,ALLOCATABLE:: sPO4W(:,:) !<PO4 in the water at current time
        Real,ALLOCATABLE:: sPO4WP(:,:) !<PO4 in the water at past time
	    Real,ALLOCATABLE:: sPO4S(:) !<PO4 in the sediment at current time
        Real,ALLOCATABLE:: sPO4SP(:) !<PO4 in the sediment at past time
        !6.2 Initial condition    
	    Real sPO4W0     !<Initial condition of PO4 in the water (mgP/l)
        Real sPO4S0     !<Initial condition of PO4 in the sediment (mgP/l)
        !6.3 Boundary condition 
        Real,ALLOCATABLE:: uPLoadPO4(:,:) !< PO4 loading value 
        Integer:: NuPLoadPO4 !< Number of cells with boundary condition of PO4
        Integer,ALLOCATABLE:: PO4nTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: PO4Time(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: PO4Value(:,:) !< value of each time-serie
        Integer,ALLOCATABLE:: PO4Smallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: PO4CapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NPO4Index(:,:)
   
        !7. PAIM
        !7.1 State variables  
	    Real,ALLOCATABLE:: sPAIMW(:,:) !<PAIM in the water at current time
        Real,ALLOCATABLE:: sPAIMWP(:,:) !<PAIM in the water at past time
	    Real,ALLOCATABLE:: sPAIMS(:) !<PAIM in the sediment at current time
        Real,ALLOCATABLE:: sPAIMSP(:) !<PAIM in the sediment at past time
        !7.2 Initial condition    
	    Real sPAIMW0     !<Initial condition of PAIM in the water (mgP/l)
        Real sPAIMS0     !<Initial condition of PAIM in the sediment (mgP/l)
        !7.3 Boundary condition 
        Real,ALLOCATABLE:: uPLoadPAIM(:,:) !< PAIM loading value 
        Integer:: NuPLoadPAIM !< Number of cells with boundary condition of PAIM
        Integer,ALLOCATABLE:: PAIMnTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: PAIMTime(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: PAIMValue(:,:) !< value of each time-serie
        Integer,ALLOCATABLE:: PAIMSmallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: PAIMCapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NPAIMIndex(:,:)
    
        !8. NH4
        !8.1 State variables  
	    Real,ALLOCATABLE:: sNH4W(:,:) !<NH4 in the water at current time
        Real,ALLOCATABLE:: sNH4WP(:,:) !<NH4 in the water at past time
	    Real,ALLOCATABLE:: sNH4S(:) !<NH4 in the sediment at current time
        Real,ALLOCATABLE:: sNH4SP(:) !<NH4 in the sediment at past time
        !8.2 Initial condition    
	    Real sNH4W0     !<Initial condition of NH4 in the water (mgP/l)
        Real sNH4S0     !<Initial condition of NH4 in the sediment (mgP/l)
        !8.3 Boundary condition 
        Real,ALLOCATABLE:: uNLoadNH4(:,:) !< NH4 loading value 
        Integer:: NuNLoadNH4 !< Number of cells with boundary condition of NH4
        Integer,ALLOCATABLE:: NH4nTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: NH4Time(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: NH4Value(:,:) !< value of each time-serie
        Integer,ALLOCATABLE:: NH4Smallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NH4CapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NNH4Index(:,:)

        !9. NO3
        !9.1 State variables  
	    Real,ALLOCATABLE:: sNO3W(:,:) !<NO3 in the water at current time
        Real,ALLOCATABLE:: sNO3WP(:,:) !<NO3 in the water at past time
	    Real,ALLOCATABLE:: sNO3S(:) !<NO3 in the sediment at current time
        Real,ALLOCATABLE:: sNO3SP(:) !<NO3 in the sediment at past time
        !9.2 Initial condition    
	    Real sNO3W0     !<Initial condition of NO3 in the water (mgP/l)
        Real sNO3S0     !<Initial condition of NO3 in the sediment (mgP/l)
        !9.3 Boundary condition 
        Real,ALLOCATABLE:: uNLoadNO3(:,:) !< NO3 loading value 
        Integer:: NuNLoadNO3 !< Number of cells with boundary condition of NO3
        Integer,ALLOCATABLE:: NO3nTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: NO3Time(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: NO3Value(:,:) !< value of each time-serie
        Integer,ALLOCATABLE:: NO3Smallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NO3CapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NNO3Index(:,:)
        
        !10. SiO2
        !10.1 State variables  
	    Real,ALLOCATABLE:: sSiO2W(:,:) !<SiO2 in the water at current time
        Real,ALLOCATABLE:: sSiO2WP(:,:) !<SiO2 in the water at past time
        !10.2 Initial condition    
	    Real sSiO2W0     !<Initial condition of SiO2 in the water (mgP/l)
        !10.3 Boundary condition 
        Real,ALLOCATABLE:: uSiLoadSiO2(:,:) !< SiO2 loading value 
        Integer:: NuSiLoadSiO2 !< Number of cells with boundary condition of SiO2
        Integer,ALLOCATABLE:: SiO2nTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: SiO2Time(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: SiO2Value(:,:) !< value of each time-serie
        Integer,ALLOCATABLE:: SiO2Smallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: SiO2CapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NSiO2Index(:,:)
        
        !11. O2
        !11.1 State variables
        Integer:: iOxyBOD, iMetabO2   
        Real,ALLOCATABLE:: sO2W(:,:) !<O2 in the water at current time
        Real,ALLOCATABLE:: sO2WP(:,:) !<O2 in the water at past time
        !11.2 Initial condition    
        Real sO2W0           !<Initial condition of O2 in the water (mgP/l)
        !11.3 Boundary condition 
        Real,ALLOCATABLE:: uO2LoadO2(:,:) !< O2 loading value
        Integer:: NuO2LoadO2 !< Number of cells with boundary condition of O2
        Integer,ALLOCATABLE:: O2nTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: O2Time(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: O2Value(:,:) !< value of each time-serie
        Integer,ALLOCATABLE:: O2Smallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: O2CapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NO2Index(:,:)

        !12. BOD
        !12.1 State variables  
        Real,ALLOCATABLE:: sDBOW(:,:) !<BOD in the water at current time
        Real,ALLOCATABLE:: sDBOWP(:,:) !<BOD in the water at past time
        !12.2 Initial condition    
        Real sDBOW0           !<Initial condition of BOD in the water (mgP/l)
        !12.3 Boundary condition 
        Real,ALLOCATABLE:: uDLoadDBO(:,:) !< BOD loading value
        Integer:: NuDBOLoadDBO !< Number of cells with boundary condition of BOD
        Integer,ALLOCATABLE:: DBOnTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: DBOTime(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: DBOValue(:,:) !< value of each time-serie
        Integer,ALLOCATABLE:: DBOSmallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: DBOCapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NDBOIndex(:,:)
    
        !13. DIC
        !13.1 State variables  
        Integer:: iCarbon, iGHG, iMetabC 
        Real,ALLOCATABLE:: sDicW(:,:) !<DIC in the water at current time
        Real,ALLOCATABLE:: sDicWP(:,:) !<DIC in the water at past time
        !13.2 Initial condition    
        Real sDicW0           !<Initial condition of DIC in the water (mgP/l)
        !13.3 Boundary condition 
        Real,ALLOCATABLE:: uDLoadDic(:,:) !< DIC loading value
        Integer:: NuDicLoadDic !< Number of cells with boundary condition of DIC
        Integer,ALLOCATABLE:: DicnTime(:) !< Number of lines in each time-serie
        Real,ALLOCATABLE:: DicTime(:,:) !< unix time of each time-serie
        Real,ALLOCATABLE:: DicValue(:,:) !< value of each time-serie
        Integer,ALLOCATABLE:: DicSmallm(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: DicCapitalM(:) !< Number of lines in each time-serie
        Integer,ALLOCATABLE:: NDicIndex(:,:)
    
        !Mass Transport numerical scheme
        Integer:: iTranspFlag
        Integer:: iAdaptativeTimeStep
    
    
    !------------------------------------------------------------------------------------------------------------------!	
    !  M O D U L O    E C O L O G I C O  =>  B I O T I C O / A B I O T I C O                                           !	
    !------------------------------------------------------------------------------------------------------------------!	
    !------------------------------------------------------------------------------------------------------------------!	
    !  B I O T I C O                                                                                                   !	
    !------------------------------------------------------------------------------------------------------------------!	
    !------------------------------------------------------------------------------------------------------------------!	
    !  B A C T E R I O P L A N C T O N                                                                                 !	
    !------------------------------------------------------------------------------------------------------------------!
    !1)Configuracoes Iniciais
	    INTEGER	InclBac
	    INTEGER nbac      !numero de grupos funcionais de bacterioplancton simulados numa aplicacao
        INTEGER nbacaer
        INTEGER nbacanaer
	    INTEGER,ALLOCATABLE:: abac(:)
	    INTEGER:: AEROBICA   = 1
	    INTEGER:: ANAEROBICA = 2
	    INTEGER:: BACTERIA   = 3
    !2)Parametros	
	    Real,ALLOCATABLE:: cThetaBac(:)              !parametro de limitacao para temperatura
	    Real pHmax,pHmin                                    !parametros de limitacao para pH
	    Real hO2Bac,fbacAn,kbacAn                           !parametros de limitacao para OD
	    Real,ALLOCATABLE:: hAssBac(:)                       !parametros de limitacao para disponibilidade de DOM - conc meia-saturacao (mg/L)
	    Real,ALLOCATABLE:: cCDBacRef(:),cPDBacRef(:),cNDBacRef(:)                  !razoes de referencia para nutrientes internos bacteria
	    Real,ALLOCATABLE:: cPrefBacDom(:,:)                 !preferencia da bac por DOM
	    Real,ALLOCATABLE:: Ydom2bac(:,:)                !eficiencia de conversao de DOM em BAC
	    Real,ALLOCATABLE:: MuBac(:)                        !taxa de crescimento (1/d)
	    Real,ALLOCATABLE:: kResbBac(:),fBac2CO2(:)   !parametros de respiracao basal, atividade e fracao CO2
	    Real,ALLOCATABLE:: hPO4UptBac(:),hNH4UptBac(:)        !conc meia-saturacao p/ consumo de nutrientes
	    Real,ALLOCATABLE:: kMortBac(:)                   !mortalidade
	    Real,ALLOCATABLE:: hSedBac(:)                                        !meia-saturacao para sedimentacao
        Real,ALLOCATABLE:: hBacMin(:)
	    Real fSiBac
    !3)Condicoes Iniciais 
	    Real,ALLOCATABLE:: sDBacW0(:) !Bacterioplancton na água (mgDW/l)
	    Real,ALLOCATABLE:: sDBacS0(:) !Bacterioplancton na sedimento (mgDW/l)
    !4)Fluxos
	    Real,ALLOCATABLE:: wDAssBac(:,:,:),     wCAssBac(:,:,:),     wNAssBac(:,:,:),     wPAssBac(:,:,:),     wSiAssBac(:,:,:)
        Real,ALLOCATABLE:: tDAssBac(:,:),     tCAssBac(:,:),     tNAssBac(:,:),     tPAssBac(:,:)
	    Real,ALLOCATABLE:: wDAssBacDom(:,:,:,:),wCAssBacDom(:,:,:,:),wNAssBacDom(:,:,:,:),wPAssBacDom(:,:,:,:),wSiAssBacDom(:,:,:,:)
        Real,ALLOCATABLE:: tDAssBacDom(:,:,:),tCAssBacDom(:,:,:),tNAssBacDom(:,:,:),tPAssBacDom(:,:,:)
        Real:: wDRespbBac, wDRespaBac,tDRespbBac, tDRespaBac
	    Real,ALLOCATABLE:: wDRespBac(:,:,:),    wCRespBac(:,:,:) , wCRespBac2CO2(:,:,:), wCRespBac2CH4(:,:,:)
        Real,ALLOCATABLE:: tDRespBac(:,:),    tCRespBac(:,:) , tCRespBac2CO2(:,:), tCRespBac2CH4(:,:)
        Real,ALLOCATABLE:: wNExcrBac(:,:,:), wPExcrBac(:,:,:)
        Real,ALLOCATABLE:: tNExcrBac(:,:), tPExcrBac(:,:)
	    Real,ALLOCATABLE:: wPUptBac(:,:,:),     wNUptBac(:,:,:)
        Real,ALLOCATABLE:: tPUptBac(:,:),     tNUptBac(:,:)
	    Real,ALLOCATABLE:: wDMortBac(:,:,:),    wCMortBac(:,:,:),    wNMortBac(:,:,:),    wPMortBac(:,:,:)
        Real,ALLOCATABLE:: tDMortBac(:,:),    tCMortBac(:,:),    tNMortBac(:,:),    tPMortBac(:,:)
	    Real,ALLOCATABLE:: tDSetBac(:,:,:),     tCSetBac(:,:,:),     tNSetBac(:,:,:),     tPSetBac(:,:,:)
    !5)Funcoes Auxiliares
	    Real:: aTLimMinBac
	    Real:: apHLimMinBac
	    Real:: aO2LimMinBac
	    Real:: aNutLimMinBac
	    Real:: aCorMinBac
	    Real,ALLOCATABLE:: oDFoodBacDom(:)
        Real,ALLOCATABLE:: oDSetBacSum(:)
	    Real:: oDFoodBac
	    Real:: aDFoodSatBac
        Real:: oDSetBac
	    !Real,ALLOCATABLE:: aDAssBacSub(:,:,:)
	    Real,ALLOCATABLE:: fDAssBacDom(:)
	    Real:: kCorResbBac
	    Real:: kCorMortBac
    !7)Variaveis de Estado	
	    Real,ALLOCATABLE:: sDBacW(:,:,:), sCBacW(:,:,:), sNBacW(:,:,:), sPBacW(:,:,:)
	    Real,ALLOCATABLE:: sDBacS(:,:),   sCBacS(:,:),   sNBacS(:,:),   sPBacS(:,:)   
	    Real,ALLOCATABLE:: sDBacWP(:,:,:),sCBacWP(:,:,:),sNBacWP(:,:,:),sPBacWP(:,:,:)
	    Real,ALLOCATABLE:: sDBacSP(:,:),  sCBacSP(:,:),  sNBacSP(:,:),  sPBacSP(:,:)
	    Real,ALLOCATABLE:: rCDBacW(:,:,:),rNDBacW(:,:,:),rPDBacW(:,:,:)                         
	    Real,ALLOCATABLE:: rCDBacS(:,:),  rNDBacS(:,:),  rPDBacS(:,:)                          

    !------------------------------------------------------------------------------------------------------------------!	
    !  F I T O P L A N C T O N                                                                                         !	
    !------------------------------------------------------------------------------------------------------------------!
    !1)Configuracoes Iniciais
	    INTEGER	InclPhyt  !Incluir os três grupos de algas? Sim=1, Não=0
	    INTEGER nphy,nphydino,nphyfcyano,nphymcyano,nphychloro,nphycryto,nphymdiato,nphyfdiato     !numero de grupos funcionais de fitoplancton simulados numa aplicacao
	    INTEGER,ALLOCATABLE:: aphy(:)
	    INTEGER LightMethodPhyt !
        INTEGER:: DINOF = 1
	    INTEGER:: CYANO = 2
	    INTEGER:: NODUL = 3
	    INTEGER:: CHLOR = 4
	    INTEGER:: CRYPT = 5
	    INTEGER:: MDIAT = 6
	    INTEGER:: FDIAT = 7
        
    !2)Parametros(pg)	
        !Produção primária
	    Real cExtWat								 !coeficiente de atenuação da luz na água (1/m)
	    Real cExtSpIM
	    Real,ALLOCATABLE:: cExtSpPom(:) !gpom
	    Real,ALLOCATABLE:: cExtSpDom(:) !gdom             !coeficiente de atenuação da luz na água devido ao detritos (m²/gD)
	    Real,ALLOCATABLE:: cExtSpZoo(:) !gzoo             !coeficiente de atenuação da luz na água devido ao detritos (m²/gD)
	    Real fPAR									 !Fração de radiação fotosinteticamente ativa (PAR)
	    Real fRefl									 !Fração da luz refletida na superfície (-)
	    Real:: cPACoefMin								 !coeficiente minimo de Poole-Atkins (1/dia)
	    Real:: cPACoefMax								 !coeficiente máximo de Poole-Atkins (1/dia)
	    Real:: hPACoef					 !constante de meia-saturação para o coef. de Poole-Atkins com materia organica (g/m²)
	    Real cfDayAve								 !fotoperíodo médio anual (h/24h)
	    Real cfDayVar								 !variação anual do fotoperíodo (h/24h)
	    Real,ALLOCATABLE:: cExtSpPhyt(:)			 !coeficiente de atenuação da luz na água devido a algas total (m²/gD)
	    Real,ALLOCATABLE:: hLRefPhyt(:)			 !meia-saturação do PAR para algas total a 20oC (W/m²)
	    Real,ALLOCATABLE:: cLOptRefPhyt(:)		     !otimo PAR (W/m²)
	    Real,ALLOCATABLE:: cMuMaxPhyt(:)			 !taxa máxima de crescimento das algas total (1/dia)
	    Real,ALLOCATABLE:: cChDPhytMax(:)		     !razão maxima Clorofila/D na algas total (mgChl/mgDW)
	    Real,ALLOCATABLE:: cChDPhytMin(:)  		 !razão minima Clorofila/D na algas total (mgChl/mgDW)
	    Real,ALLOCATABLE:: hSiAssDiat(:)								 !meia-saturação de conc de Si para cresc. algal (mgSi/l)
        !Consumo de Nutrientes
	    Real cTmRef                                  !temperatura de referencia (oC)
	    Real,ALLOCATABLE:: cSigTmPhyt(:)            !temperatura constante para Algas total-sigma na curva Gaussian (oC)
	    Real,ALLOCATABLE:: cTmOptPhyt(:)            !temperatura otima para algas (oC)
	    Real,ALLOCATABLE:: cVCUptMaxPhyt(:)         !Capacidade máxima de consumo de N pela algas total (mgN/mgDW/dia)
	    Real,ALLOCATABLE:: cVNUptMaxPhyt(:)         !Capacidade máxima de consumo de N pela algas total (mgN/mgDW/dia)
	    Real,ALLOCATABLE:: cVPUptMaxPhyt(:)         !Capacidade máxima de consumo de P pela algas total (mgP/mgDW/dia)
	    Real,ALLOCATABLE:: cVSiUptMaxPhyt(:)        !Capacidade máxima de consumo de Si pela diatomáceas (mgN/mgDW/dia)
	    Real,ALLOCATABLE:: cCDPhytMin(:)
	    Real,ALLOCATABLE:: cCDPhytMax(:)
	    Real,ALLOCATABLE:: cNDPhytMax(:)            !razão máxima N/D pela algas (mgN/mgDW)
	    Real,ALLOCATABLE:: cNDPhytMin(:)            !razão mínima N/D pela algas (mgN/mgDW)
	    Real,ALLOCATABLE:: cPDPhytMax(:)            !razão máxima P/D pela algas (mgP/mgDW)
	    Real,ALLOCATABLE:: cPDPhytMin(:)            !razão mínima P/D pela algas (mgP/mgDW)
	    Real,ALLOCATABLE:: cSiDPhytMax(:)           !razão máxima Si/D pela diatomaceas (mgSi/mgDW)
	    Real,ALLOCATABLE:: cSiDPhytMin(:)           !razão mínima Si/D pela diatomaceas (mgSi/mgDW)
	    Real,ALLOCATABLE:: cAffCUptPhyt(:)          !Afinidade de consumo de C para algas total (1/mgDW/dia)
        Real,ALLOCATABLE:: cAffNUptPhyt(:)          !Afinidade de consumo de N para algas total (1/mgDW/dia)
	    Real,ALLOCATABLE:: cAffPUptPhyt(:)          !Afinidade de consumo de P para algas total (1/mgDW/dia)
	    Real,ALLOCATABLE:: cAffSiUptPhyt(:)         !Afinidade de consumo de Si para diatomaceas (1/mgDW/dia)
	    !Exudacao
	    !Real,ALLOCATABLE:: cExdPhytW(:)             !Fracao da prod primaria exudada como DOM(-)
        !Respiração
	    Real,ALLOCATABLE:: kDRespbPhyt(:)           !taxa de respiração basal para algas (1/dia)
        Real,ALLOCATABLE:: alfap(:)                 !parcela da prod primaria respirada como CO2 (-)
        !Sedimentação
	    Real,ALLOCATABLE:: cVSetPhyt(:)             !velocidade de sedimentação para as algas total (m/dia)
        !Lise
	    Real,ALLOCATABLE:: cCDPhytOpt(:)
	    Real,ALLOCATABLE:: cPDPhytOpt(:)
	    Real,ALLOCATABLE:: cNDPhytOpt(:)
	    Real,ALLOCATABLE:: cSiDPhytOpt(:)	
	    Real,ALLOCATABLE:: cLisPhytW(:)
	    !Real,ALLOCATABLE:: hLisNutPhyt(:)          !taxa de mortalidade das algas na água (1/dia)
	    Real,ALLOCATABLE:: cLisPhytS(:)            !taxa de mortalidade das algas no Sed (1/dia)
        
    !3)Condicoes Iniciais 
	    Real,ALLOCATABLE:: sDPhytW0(:)             !fitoplancton na água (mgDW/l)
	    Real,ALLOCATABLE:: sDPhytS0(:)             !fitoplancton no sedimento (mgDW/l)
    !4)Fluxos - Agua(K,gg),Sed(K,gg)
	    !Consumo de nutrientes
	    Real,ALLOCATABLE:: wCUptPhyt(:,:,:),wNUptPhyt(:,:,:),wPUptPhyt(:,:,:),wSiUptPhyt(:,:,:) !fluxo de consumo total (gP/m³/d)
	    Real,ALLOCATABLE:: wNUptNH4Phyt(:,:,:),wNUptNO3Phyt(:,:,:)                          !fluxo de consumo total de NO3 (gN/m³/d)
	    !Produção Primária
	    Real,ALLOCATABLE:: wDAssPhyt(:,:,:)
	    !Exudacao/Respiração 
	    Real:: wDRespbPhyt,wDResppPhyt !basal, fotossintetica
	    Real,ALLOCATABLE:: wDRespPhyt(:,:,:), wCExcrPhytW(:,:,:), wPExcrPhytW(:,:,:), wNExcrPhytW(:,:,:) 
	    Real,ALLOCATABLE:: tDRespPhytS(:,:),  tCExcrPhytS(:,:), tPExcrPhytS(:,:), tNExcrPhytS(:,:) 
	    !Resuspensão
        Real:: akResusPhytRef
	    Real,ALLOCATABLE:: tDResusPhyt(:,:),tCResusPhyt(:,:),tPResusPhyt(:,:),tNResusPhyt(:,:),tSiResusPhyt(:,:)
	    !Sedimentação
        Real,ALLOCATABLE:: uCorVSetPhyt(:),StokesVelPhyt(:,:),ppPhyt(:),diPhyt(:)
	    Real,ALLOCATABLE:: tDSetPhyt(:,:,:) !,tCSetPhyt(:,:),tPSetPhyt(:,:),tNSetPhyt(:,:),tSiSetPhyt(:,:)
	    !Lise
	    Real,ALLOCATABLE:: wDLisPhyt(:,:,:),   wCLisPhyt(:,:,:),   wNLisPhyt(:,:,:),   wPLisPhyt(:,:,:),   wSiLisPhyt(:,:,:)
        !Real,ALLOCATABLE:: wPLisPhytPO4W(:,:,:),wNLisPhytNH4W(:,:,:),wCLisPhytCO2W(:,:,:)
	    !Real,ALLOCATABLE:: wDLisPomPhyt(:,:,:),wCLisPomPhyt(:,:,:),wNLisPomPhyt(:,:,:),wPLisPomPhyt(:,:,:),wSiLisPomPhyt(:,:,:) 
	    !Real,ALLOCATABLE:: wDLisDomPhyt(:,:,:),wCLisDomPhyt(:,:,:),wNLisDomPhyt(:,:,:),wPLisDomPhyt(:,:,:),wSiLisDomPhyt(:,:,:)
	    Real,ALLOCATABLE:: tDLisPhytS(:,:),    tCLisPhytS(:,:),    tNLisPhytS(:,:),    tPLisPhytS(:,:),    tSiLisPhytS(:,:)
	    Real,ALLOCATABLE:: tDLisPomPhytS(:,:), tCLisPomPhytS(:,:), tNLisPomPhytS(:,:), tPLisPomPhytS(:,:), tSiLisPomPhytS(:,:)
	    Real,ALLOCATABLE:: tDLisDomPhytS(:,:), tCLisDomPhytS(:,:), tNLisDomPhytS(:,:), tPLisDomPhytS(:,:), tSiLisDomPhytS(:,:)
    !5)Funcoes Auxiliares(gg)
	    !Consumo de nutrientes
	    Real:: uFunTmPhyt                  !função de temperatura para algas (-)
	    Real:: aVCUptMaxCorPhyt            !taxa de consumo máxima (mgC/mgD/d)
	    Real:: aVNUptMaxCorPhyt            !taxa de consumo máxima (mgN/mgD/d)
	    Real:: aVPUptMaxCorPhyt            !taxa de consumo máxima (mgP/mgD/d)
	    Real:: aVSiUptMaxCorPhyt           !taxa de consumo máxima (mgSi/mgD/d)
	    Real:: aVCUptPhyt                  !taxa de consumo de C especifica (mgC/mgD/d)
	    Real:: aVNUptPhyt                  !taxa de consumo de N especifica (mgN/mgD/d)
	    Real:: aVPUptPhyt                  !taxa de consumo de P especifica (mgP/mgD/d)
	    Real:: aVSiUptPhyt                 !taxa de consumo de Si especifica (mgSi/mgD/d)
	    Real:: ahPUptPhyt                  !concentração SRP de meia-saturação (mgP/l)
	    Real:: ahNUptPhyt                  !concentração SRN de meia-saturação (mgN/l)
	    Real:: afNH4UptPhyt                !Fração de N absorvido como amonia
	    !Produção Primária
	    Real ufDay
	    Real aExtVeg,aPACoef
	    Real,ALLOCATABLE:: aExtIM(:,:),aExtPom(:,:,:),aExtDom(:,:,:),aExtZoo(:,:,:),aExtPhyt(:,:,:),aExtCoefOpen(:,:)
	    Real,ALLOCATABLE:: aExtCoef(:,:),aSecchi(:,:)
	    Real,ALLOCATABLE:: uLPAR0(:),aLPARBot(:)
	    Real:: aLLimPhyt,uhLPhyt
	    Real:: aMuTmLPhyt
	    Real aCLimPhyt,aPLimPhyt,aNLimPhyt,aSiLimPhyt
	    Real aNutLimPhyt,aMuPhyt
	    Real,ALLOCATABLE:: rChDPhyt(:)
	    !Respiração
	    Real ukDRespbTmPhyt
	    !Lise
	    !Real,ALLOCATABLE:: aCLisPhyt(:), aPLisPhyt(:), aNLisPhyt(:), aSiLisPhyt(:)
	    !Real,ALLOCATABLE:: aCLisPhytS(:),aPLisPhytS(:),aNLisPhytS(:),aSiLisPhytS(:)
	    !Real,ALLOCATABLE:: aNutLisPhyt(:), aLisPomPhyt(:)
	    !Real,ALLOCATABLE:: aNutLisPhytS(:),aLisPomPhytS(:)
    !6)Derivadas	
	    Real,ALLOCATABLE::  dDPhytW(:),      dCPhytW(:),       dNPhytW(:),       dPPhytW(:),       dSiPhytW(:) 
	    Real,ALLOCATABLE::  dDPhytS(:),      dCPhytS(:),       dNPhytS(:),       dPPhytS(:),       dSiPhytS(:)
    !7)Variaveis de Estado	
	    Real,ALLOCATABLE:: sDPhytW(:,:,:), sCPhytW(:,:,:), sNPhytW(:,:,:), sPPhytW(:,:,:), sSiPhytW(:,:,:)
        Real,ALLOCATABLE:: sDPhytS(:,:),   sCPhytS(:,:),   sNPhytS(:,:),   sPPhytS(:,:),   sSiPhytS(:,:)
	    Real,ALLOCATABLE:: sDPhytWP(:,:,:),sCPhytWP(:,:,:),sNPhytWP(:,:,:),sPPhytWP(:,:,:),sSiPhytWP(:,:,:)
        Real,ALLOCATABLE:: sDPhytSP(:,:),  sCPhytSP(:,:),  sNPhytSP(:,:),  sPPhytSP(:,:),  sSiPhytSP(:,:)
	    Real,ALLOCATABLE:: rCDPhytW(:,:,:),rNDPhytW(:,:,:),rPDPhytW(:,:,:),rSiDPhytW(:,:,:)     
	    Real,ALLOCATABLE:: rCDPhytS(:,:),  rNDPhytS(:,:),  rPDPhytS(:,:),  rSiDPhytS(:,:)       

    !------------------------------------------------------------------------------------------------------------------!	
    !  M A C R O F I T A S                                                                                             !	
    !------------------------------------------------------------------------------------------------------------------!	
    !1)Configuracoes Iniciais
	    INTEGER InclMacr
	    INTEGER nmac,nmacElod,nmacCharo,nmacCera,nmacLemna,nmacNymp,nmacHelo,nmacGeral      !numero de grupos funcionais de macrofitas simulados numa aplicacao
	    INTEGER,ALLOCATABLE:: amac(:)
	    INTEGER:: ELOD = 1   !Submersas
        INTEGER:: CHARO = 2   !Submersas
	    INTEGER:: CERA = 3   !Nao enraizadas
        INTEGER:: LEMNA = 4   !Nao enraizadas
        INTEGER:: NYMP = 5   !Flutuantes
	    INTEGER:: HELO = 6   !Emergentes
	    INTEGER:: MACR = 7   !Geral
    
    
    
    !2)Parametros(mg)
 	    !Produção primária
        Real,ALLOCATABLE:: cTmInitVeg(:) !temperature for initial growth
        Real,ALLOCATABLE:: fRootVegWin(:)
        Real,ALLOCATABLE:: fRootVegSum(:)
        Real,ALLOCATABLE:: fEmergVeg(:)
        Real,ALLOCATABLE:: cDLayerVeg(:)
        Real,ALLOCATABLE:: fFloatVeg(:)
        Real,ALLOCATABLE:: cCovSpVeg(:)
	    !Real,ALLOCATABLE:: cExtSpMac    !coeficiente de atenuação da luz na água devido a vegetação (m²/gD)
 	   ! Real,ALLOCATABLE:: cCovSpMac    !cobertura por gD/m² de elodea (%)
	    Real,ALLOCATABLE:: cMuMaxMac(:)    !taxa máxima de crescimento da Elodea (1/dia)
	    !Real,ALLOCATABLE:: fDepth1Mac   !limite maximo da planta, como fração da profundidade da água
	    !Real,ALLOCATABLE:: fDepth2Mac   !limite minimo da planta, como fração da profundidade da água
	    Real,ALLOCATABLE:: hLRefMac(:)     !meia-saturação da luz à 20oC (W/m² PAR)
	    Real,ALLOCATABLE:: cDCarrMac(:)    !Capacidade de crescimento maxima (gD/m²)
	    !Real,ALLOCATABLE:: cDLayerMac   !biomassa de uma simples camada de folhas flutuantes (gD/m²)
	    !Real,ALLOCATABLE:: fEmergMac
        Real,ALLOCATABLE:: fSedUptVegMax(:)
        Real,ALLOCATABLE:: fSedUptVegCoef(:)
        Real,ALLOCATABLE:: fSedUptVegExp(:)
        !Consumo de Nutrientes
	    Real,ALLOCATABLE:: cQ10ProdMac(:)    !expoente de produção para a temperatura (-)
	    Real,ALLOCATABLE:: cVCUptMaxMac(:)   !Capacidade máxima de consumo de C (mgN/mgD/dia)
	    Real,ALLOCATABLE:: cVNUptMaxMac(:)   !Capacidade máxima de consumo de N (mgN/mgD/dia)
	    Real,ALLOCATABLE:: cVPUptMaxMac(:)   !Capacidade máxima de consumo de P (mgP/mgD/dia)
	    Real,ALLOCATABLE:: cAffCUptMac(:)    !Afinidade máxima de consumo de C (1/mgD/dia)
	    Real,ALLOCATABLE:: cAffNUptMac(:)    !Afinidade máxima de consumo de N (1/mgD/dia)
	    Real,ALLOCATABLE:: cAffPUptMac(:)    !Afinidade máxima de consumo de P (1/mgD/dia)
	    Real,ALLOCATABLE:: hCRoot(:),hCShoot(:) !conc de meia-saturacao para a assimilacao de C roots e shoot (mgC/L,gC/m3)
	    Real,ALLOCATABLE:: hNShoot(:),hNRoot(:) !conc de meia-saturacao para a assimilacao de N roots e shoot (mgN/L,gN/m3)
	    Real,ALLOCATABLE:: hPShoot(:),hPRoot(:) !conc de meia-saturacao para a assimilacao de P roots e shoot (mgP/L,gP/m3)
	    Real,ALLOCATABLE:: hNSatMac(:)       !conc de satutacao de N (mgN/L,gN/m3)
	    Real,ALLOCATABLE:: cCDMacMax(:)      !razão C/D máxima (mgC/mgD)
	    Real,ALLOCATABLE:: cCDMacMin (:)     !razão C/D minima (mgC/mgD)
	    Real,ALLOCATABLE:: cNDMacMax(:)     !razão N/D máxima (mgN/mgD)
	    Real,ALLOCATABLE:: cNDMacMin(:)      !razão N/D minima (mgN/mgD)
	    Real,ALLOCATABLE:: cPDMacMax(:)      !razão P/D máxima (mgP/mgD)
	    Real,ALLOCATABLE:: cPDMacMin(:)      !razão P/D minima (mgP/mgD)
        
        !Exudacao
	    Real,ALLOCATABLE:: fExdMac(:)        !parcela da prod primaria exudada como DOM (-)
        !Respiração
        Real,ALLOCATABLE:: cQ10RespMac(:)    !expoente de produção para a temperatura (-)
	    Real,ALLOCATABLE:: kDRespMac(:)      !taxa de respiração no escuro (1/dia)
	    Real,ALLOCATABLE:: fPhotoRespMac(:)  !parcela da prod primaria respirada como CO2 (-)
        !Mortalidade
	    Real,ALLOCATABLE:: cSigTmMac(:),cTmOptMac(:) !Parametros da curva de gauss para representacao do aumento da mortalidade de MAC durante o inverno (celsius)
	    Real,ALLOCATABLE:: kMortMac(:)            !taxa de mortalidade de vegetação
        !Consumo por aves
	    Real,ALLOCATABLE:: hDMacBird(:)      !meia-saturação de biomassa
	    Real,ALLOCATABLE:: cPrefMacBird(:)   !comestível por aves
	    Real cBirdsPerha                      !Numéros de aves por ha do lago
	    Real,ALLOCATABLE:: cDGrazPerBird(:)                    !consumo diario de macrofitas pelas aves
	    Real,ALLOCATABLE:: fDAssBird(:)                    !Eficiencia de assimilação das aves
	    !Shoot<->Root
	    Real,ALLOCATABLE:: kroot2shoot(:)    !taxa de realocacao de biomassa Shoot<->Root (1/dia)
    !3)Condicoes Iniciais (mm)
	    Real,ALLOCATABLE:: sDMac0(:)         ! biomassa de MAC(gD/m²)
	    Real,ALLOCATABLE:: cCDMac0(:)        ! (gC/gD)
	    Real,ALLOCATABLE:: cNDMac0(:)        ! (gN/gD)        
	    Real,ALLOCATABLE:: cPDMac0(:)        ! (gP/gD)       
	    Real,ALLOCATABLE:: hmac(:)           ! prof máxima com presenca de MAC
    !4)Fluxos (mm)
	    !Consumo de Nutrientes
	    Real,ALLOCATABLE:: tCUptMacW(:,:),tNUptMacW(:,:),tPUptMacW(:,:),tNUptNH4MacW(:,:),tNUptNO3MacW(:,:)
	    Real,ALLOCATABLE:: tCUptMacS(:,:),tNUptMacS(:,:),tPUptMacS(:,:),tNUptNH4MacS(:,:),tNUptNO3MacS(:,:)
	    !Produção Primária
        Real:: tDEnvMac, tDEnvProdMac,tDEnvMortMac
        Real,ALLOCATABLE:: tDProdMac(:,:)
	    !Exudacao
        Real,ALLOCATABLE:: tCExdMac(:,:),tNExdMac(:,:),tPExdMac(:,:)
	    Real,ALLOCATABLE:: tCExdMacW(:,:),tNExdMacW(:,:),tPExdMacW(:,:)
	    Real,ALLOCATABLE:: tCExdMacS(:,:),tNExdMacS(:,:),tPExdMacS(:,:)		
	    !Respiração
        Real,ALLOCATABLE:: tDRespMac(:,:),tCRespMac(:,:)
	    Real,ALLOCATABLE:: tDRespMacW(:,:)
	    Real,ALLOCATABLE:: tDRespMacS(:,:)
        
	    !Mortalidade
        Real,ALLOCATABLE:: tDMortMac(:,:),tCMortMac(:,:),tNMortMac(:,:),tPMortMac(:,:)
	    Real,ALLOCATABLE:: tDMortMacW(:,:),tCMortMacW(:,:),tNMortMacW(:,:),tPMortMacW(:,:)
	    Real,ALLOCATABLE:: tDMortMacS(:,:),tCMortMacS(:,:),tNMortMacS(:,:),tPMortMacS(:,:)
	    !Consumo por aves
	    Real,ALLOCATABLE:: tDGrazMacBird(:,:),tCGrazMacBird(:,:),tNGrazMacBird(:,:),tPGrazMacBird(:,:)
	    Real,ALLOCATABLE:: tDEgesBird(:,:),tCEgesBird(:,:),tNEgesBird(:,:),tPEgesBird(:,:)
        !oxygen
        Real,ALLOCATABLE:: tO2UptNO3MacW(:,:)
        Real,ALLOCATABLE:: tO2RespMacW(:,:)
        Real,ALLOCATABLE:: tO2RespMacS(:,:)
        Real,ALLOCATABLE:: tO2ProdMac(:,:)
        Real,ALLOCATABLE:: tO2ProdMacW(:,:)
        Real,ALLOCATABLE:: tO2ProdMacS(:,:)
    !5)Funcoes Auxiliares (mm)
	    !Consumo de Nutrientes
        Real:: uFunTmProdMac, uFunTmRespMac
 	    Real:: afCovSurfVeg
        Real:: aDayInitVeg,bfShootVeg,bfRootVeg,aDRootVeg,aDShootVeg,bfSubVeg,aDSubVeg,aDEmergVeg,aDFloatVeg,afCovEmergVeg,aCovVeg
        Real:: afCUptVegS,afNUptVegS,afPUptVegS
	    Real:: aVCUptMaxCorMac,aVNUptMaxCorMac,aVPUptMaxCorMac
	    Real:: aVCUptMacW,aVNUptMacW,aVPUptMacW
	    Real:: aVCUptMacS,aVNUptMacS,aVPUptMacS
	    Real:: ahCUptMac,ahNUptMac,ahPUptMac
	    Real:: afNH4UptMacW
	    Real:: afNH4UptMacS
	    !Produção Primária
	    Real:: uMuMaxTmMac, aLPAR1Mac, aLPAR2Mac, uhLMac, aFunLSubMac
	    Real:: aMuTmLMac, aCLimMac,aPLimMac, aNLimMac, aNutLimMac
	    Real:: aMuMac, akDIncrMac
        
        
	    !Respiração
	    Real,ALLOCATABLE:: ukDRespTmMac
	    !Mortalidade
	    Real:: bkMortMac
    !7)Variaveis de Estado (I,J,mm)	
	    !Total
	    Real,ALLOCATABLE:: sDMac(:,:),   sCMac(:,:),   sNMac(:,:),   sPMac(:,:)
        Real,ALLOCATABLE:: sDMacP(:,:),   sCMacP(:,:),   sNMacP(:,:),   sPMacP(:,:)
        !Razao P/D de macrofitas na água (gP/gD)
	    Real,ALLOCATABLE:: rCDMac(:,:), rNDMac(:,:), rPDMac(:,:)    
    !------------------------------------------------------------------------------------------------------------------!	
    !  Z O O P L A N C T O N                                                                                           !	
    !------------------------------------------------------------------------------------------------------------------!	
    !1)Configuracoes Iniciais
	    INTEGER InclZoo  !Incluir zooplântons? Sim=1, Não=0
	    INTEGER InclZooDist
	    INTEGER zzaux,zgaux
	    INTEGER nzoo,nzoomicro,nzoomeso,nzoomacro     !numero de grupos funcionais de zooplancton simulados numa aplicacao
	    INTEGER,ALLOCATABLE:: azoo(:)
	    INTEGER:: MACRO = 1
	    INTEGER:: MESO  = 2
	    INTEGER:: MICRO = 3
	    INTEGER:: ZOOP  = 4
    !2)Parametros(zg)
	    Real,ALLOCATABLE:: cSigTmZoo(:),cTmOptZoo(:)
	    Real,ALLOCATABLE:: fDAssZoo(:)
	    Real,ALLOCATABLE:: cPrefZooZoo(:,:) !(gzoo,gzoo)
	    Real,ALLOCATABLE:: cPrefZooPhyt(:,:) !(gzoo,gphy)
	    Real,ALLOCATABLE:: cPrefZooBac(:,:) !(gzoo,gbaz)
	    Real,ALLOCATABLE:: cPrefZooPom(:,:) !(gzoo,gpom)
	    Real,ALLOCATABLE:: cDCarrZoo(:) !(gzoo)
	    Real,ALLOCATABLE:: hFilt(:) !(gzoo)
	    Real,ALLOCATABLE:: cFiltMax(:) !(gzoo)
	    Real,ALLOCATABLE:: kDRespZoo(:) !(gzoo)
	    !Real,ALLOCATABLE:: fDMessyZoo(:) !(gzoo)
	    Real,ALLOCATABLE:: kPelZoo(:) !(gzoo)
	    Real,ALLOCATABLE:: kMortZoo(:) !(gzoo)
	    Real,ALLOCATABLE:: cCDZooRef(:) !(gzoo)
	    Real,ALLOCATABLE:: cNDZooRef(:) !(gzoo)
	    Real,ALLOCATABLE:: cPDZooRef(:) !(gzoo)
	    Real zooMinBiomass
	    Real,ALLOCATABLE:: ODmin(:)	
    !3)Condicoes Iniciais(zz)
	    Real,ALLOCATABLE:: sDZoo0(:)     !zooplancton (mgDW/l)
    !4)Fluxos(K,zz) ou (K,zz,presa)
	    Real,ALLOCATABLE:: wDGrzZoo(:,:),      wCGrzZoo(:,:),      wNGrzZoo(:,:),      wPGrzZoo(:,:)
	    Real,ALLOCATABLE:: wDAssZoo(:,:,:),      wCAssZoo(:,:,:),      wNAssZoo(:,:,:),      wPAssZoo(:,:,:)
	    Real:: wDEnvZoo
        Real:: wDConsZoo,wCConsZoo,wNConsZoo,wPConsZoo
	    Real,ALLOCATABLE:: wDGrzZooZoo(:,:,:,:), wCGrzZooZoo(:,:,:,:), wNGrzZooZoo(:,:,:,:), wPGrzZooZoo(:,:,:,:)
	    Real,ALLOCATABLE:: wDGrzZooPhyt(:,:,:,:),wCGrzZooPhyt(:,:,:,:),wNGrzZooPhyt(:,:,:,:),wPGrzZooPhyt(:,:,:,:),wSiGrzZooPhyt(:,:,:,:)
	    Real,ALLOCATABLE:: wDGrzZooBac(:,:,:,:), wCGrzZooBac(:,:,:,:), wNGrzZooBac(:,:,:,:), wPGrzZooBac(:,:,:,:)
	    Real,ALLOCATABLE:: wDGrzZooPom(:,:,:,:), wCGrzZooPom(:,:,:,:), wNGrzZooPom(:,:,:,:), wPGrzZooPom(:,:,:,:)
	    Real,ALLOCATABLE:: wDRespZoo(:,:,:),wCRespZoo(:,:,:)
	    Real,ALLOCATABLE:: wDEgesZoo(:,:,:),     wCEgesZoo(:,:,:),     wNEgesZoo(:,:,:),     wPEgesZoo(:,:,:)
	    Real,ALLOCATABLE:: wDMortZoo(:,:,:),     wCMortZoo(:,:,:),     wNMortZoo(:,:,:),     wPMortZoo(:,:,:)
	    Real,ALLOCATABLE:: wNExcrZoo(:,:,:),      wPExcrZoo(:,:,:)
	    Real,ALLOCATABLE:: wDMessZooOM(:,:,:),   wCMessZooOM(:,:,:),   wNMessZooOM(:,:,:),   wPMessZooOM(:,:,:)
    !5)Funcoes Auxiliares(zz)	
	    Real:: aFilt
	    Real:: uFunTmZoo,aDSatZoo
	    Real:: oDFoodZooZoo,oDFoodZooPhyt,oDFoodZooBac,oDFoodZooPom
	    Real:: oDFoodZoo,oCFoodZoo,oNFoodZoo,oPFoodZoo
	    Real,ALLOCATABLE:: aDSatFoodZoo(:)
	    Real,ALLOCATABLE:: fDGrzZooZoo(:,:),fDGrzZooPhyt(:,:),fDGrzZooBac(:,:),fDGrzZooPom(:,:)
	    Real:: ukDRespTmZoo,aCorDRespZoo
	    Real:: ukDAssTmZoo
	    Real:: ukDIncrZoo
	    Real,ALLOCATABLE:: akCExcrZoo(:),akNExcrZoo(:),akPExcrZoo(:)
	    Real:: afCAssZoo ,afNAssZoo ,afPAssZoo
	    Real:: rCDFoodZoo,rNDFoodZoo,rPDFoodZoo
    !7)Variaveis de Estado(I,J,K,zz)	
	    Real,ALLOCATABLE:: sDZoo(:,:,:),  sCZoo(:,:,:),  sNZoo(:,:,:),  sPZoo(:,:,:)                        
	    Real,ALLOCATABLE:: sDZooP(:,:,:), sCZooP(:,:,:), sNZooP(:,:,:), sPZooP(:,:,:)                       
	    Real,ALLOCATABLE:: rCDZoo(:,:,:), rNDZoo(:,:,:), rPDZoo(:,:,:)

    !------------------------------------------------------------------------------------------------------------------!	
    !  Z O O B E N T O S                                                                                               !	
    !------------------------------------------------------------------------------------------------------------------!	
    !1)Configuracoes Iniciais
	    INTEGER InclBent !Incluir Macrobentos? Sim=1, Não=0
	    INTEGER nben,nbenmicro,nbenmeso,nbenmacro    !numero de grupos funcionais de zoobentos simulados numa aplicacao
	    INTEGER,ALLOCATABLE:: aben(:)
	    INTEGER:: BENTOS=1
    !2)Parametros(beng)	
	    Real,ALLOCATABLE:: cSigTmBent(:),cTmOptBent(:)
	    Real,ALLOCATABLE:: cPrefBentPhyt(:,:) !(gben,gphy)
	    Real,ALLOCATABLE:: cPrefBentBac(:,:) !(gben,gbac)
	    Real,ALLOCATABLE:: cPrefBentPom(:,:) !gben,gpom
        Real,ALLOCATABLE:: hDFoodBent(:)
        Real,ALLOCATABLE:: kDAssBent(:)
        Real,ALLOCATABLE:: cDCarrBent(:)
        Real,ALLOCATABLE:: fDAssBent(:)
	    Real,ALLOCATABLE:: kDRespBent(:)
	    Real,ALLOCATABLE:: kPelBent(:)
	    Real,ALLOCATABLE:: kMortBent(:)
	    Real,ALLOCATABLE:: cCDBentRef(:)
	    Real,ALLOCATABLE:: cNDBentRef(:)
	    Real,ALLOCATABLE:: cPDBentRef(:)
    !3)Condicoes Iniciais(ben) 
	    Real,ALLOCATABLE:: sDBent0(:)     !Bentos (mgDW/m²)
    !4)Fluxos(ben)
	    Real,ALLOCATABLE:: tDGrzBent(:),      tCGrzBent(:),      tNGrzBent(:),      tPGrzBent(:)
	    Real,ALLOCATABLE:: tDGrzBentPhyt(:,:,:),tCGrzBentPhyt(:,:,:),tNGrzBentPhyt(:,:,:),tPGrzBentPhyt(:,:,:),tSiGrzBentPhyt(:,:,:)
	    Real,ALLOCATABLE:: tDGrzBentBac(:,:,:), tCGrzBentBac(:,:,:), tNGrzBentBac(:,:,:), tPGrzBentBac(:,:,:)
	    Real,ALLOCATABLE:: tDGrzBentPom(:,:,:), tCGrzBentPom(:,:,:), tNGrzBentPom(:,:,:), tPGrzBentPom(:,:,:)
	    Real,ALLOCATABLE:: tDRespBent(:,:),     tCRespBent(:,:)
	    Real,ALLOCATABLE:: tDEgesBent(:,:),     tCEgesBent(:,:),     tNEgesBent(:,:),     tPEgesBent(:,:)
	    Real,ALLOCATABLE:: tDMortBent(:,:),     tCMortBent(:,:),     tNMortBent(:,:),     tPMortBent(:,:)
	    Real,ALLOCATABLE:: tNExcBent(:,:),      tPExcBent(:,:)
	    Real,ALLOCATABLE:: tDMessBentOM(:,:),   tCMessBentOM(:,:),   tNMessBentOM(:,:),   tPMessBentOM(:,:)
    !5)Funcoes Auxiliares(ben)
	    Real:: uFunTmBent
	    Real:: oDFoodBentPhyt,oDFoodBentBac,oDFoodBentPom
	    Real:: oDFoodBent,oCFoodBent,oNFoodBent,oPFoodBent
	    Real:: aDSatFoodBent
	    Real:: fDGrzBentPhyt,fDGrzBentBac,fDGrzBentPom
	    Real:: ukDRespTmBent,aCorDRespBent
	    Real:: aDSatBent
	    Real:: ukDIncrBent,tDEnvBent,tDConsBent,tDAssBent,tCAssBent,tNAssBent,tPAssBent
        Real:: rCDFoodBent,tCConsBent,afCAssBent,rNDFoodBent,tNConsBent,afNAssBent,rPDFoodBent,tPConsBent,afPAssBent
        
    !7)Variaveis de Estado(I,J,ben)	
	    Real,ALLOCATABLE:: sDBent(:,:),  sCBent(:,:),  sNBent(:,:),  sPBent(:,:)             !Zoobentos (g/m²)
	    Real,ALLOCATABLE:: sDBentP(:,:), sCBentP(:,:), sNBentP(:,:), sPBentP(:,:)            !Zoobentos (g/m²)
	    Real,ALLOCATABLE:: rCDBent(:,:), rNDBent(:,:), rPDBent(:,:)                            !Razao P/D e N/D de Zooplancton e Bentos na água (gP/gD)
	
    !------------------------------------------------------------------------------------------------------------------!	
    !  P E I X E S                                                                                                     !	
    !------------------------------------------------------------------------------------------------------------------!	
    !1)Configuracoes Iniciais
	    INTEGER InclFish  !Incluir Peixes? Sim=1, Não=0
        INTEGER iFishMov  !Include Fish movement?
        INTEGER iFishRefuge  !Include Fish refuge?
        INTEGER iFishStage  !Include Stage of Live (Juvenile and Adult)?
	    INTEGER ffaux,fgaux
	    INTEGER nfish,nfishomniv,nfishplank,nfishpisc     !numero de grupos funcionais de peixes simulados numa aplicacao
	    INTEGER,ALLOCATABLE:: afish(:)
    !2)Parametros
	    !Migration
	    Real,ALLOCATABLE:: kMigrFiAd(:)
        Real,ALLOCATABLE:: cDFiAdIn(:)
	    Real,ALLOCATABLE:: kMigrFiJv(:)
        Real,ALLOCATABLE:: cDFiJvIn(:)
        !Reproduction/Aging
        Real,ALLOCATABLE:: cDayReprFish(:)
        Real,ALLOCATABLE:: fReprFish(:)
	    Real,ALLOCATABLE:: fAgeFish(:)

        Real,ALLOCATABLE:: cRelVegFish(:)
        Real,ALLOCATABLE::cSigTmFish(:),cTmOptFish(:)
    !	Real cTmRef
	    Real,ALLOCATABLE:: cPrefFiAdFiAd(:,:) !(gfish,gfish)
        Real,ALLOCATABLE:: cPrefFiAdFiJv(:,:) !(gfish,gfish)
	    Real,ALLOCATABLE::  cPrefFiAdZoo(:,:) !(gfish,gzoo)
        Real,ALLOCATABLE::  cPrefFiJvZoo(:,:) !(gfish,gzoo)
	    Real,ALLOCATABLE::  cPrefFiAdPhyt(:,:) !(gfish,gphy)
        Real,ALLOCATABLE::  cPrefFiJvPhyt(:,:) !(gfish,gphy)
	    !Real,ALLOCATABLE:: cPrefFishBac(:,:) !(gfish,gbac)
	    Real,ALLOCATABLE::  cPrefFiAdBent(:,:) !(gfish,gben)
        Real,ALLOCATABLE::  cPrefFiJvBent(:,:) !(gfish,gben)
	    Real,ALLOCATABLE::  cPrefFiAdPom(:,:) !(gfish,gpom)
        Real,ALLOCATABLE::  cPrefFiJvPom(:,:) !(gfish,gpom)
	    Real,ALLOCATABLE:: hDFiAd(:)
        Real,ALLOCATABLE:: hDFiJv(:)
	    Real,ALLOCATABLE:: kDAssFiAd(:)
        Real,ALLOCATABLE:: kDAssFiJv(:)
        Real,ALLOCATABLE:: fDAssFiAd(:)
        Real,ALLOCATABLE:: fDAssFiJv(:)
	    Real,ALLOCATABLE:: kHarvFiAd(:)
        Real,ALLOCATABLE:: kHarvFiJv(:)
	    Real,ALLOCATABLE:: kDRespFiAd(:)
        Real,ALLOCATABLE:: kDRespFiJv(:)
	    !Real,ALLOCATABLE:: fDMessyFiAd(:)
        !Real,ALLOCATABLE:: fDMessyFiJv(:)
	    Real,ALLOCATABLE:: kPelFiAd(:)
        Real,ALLOCATABLE:: kPelFiJv(:)
	    Real,ALLOCATABLE:: kMortFiAd(:)
        Real,ALLOCATABLE:: kMortFiJv(:)
        Real,ALLOCATABLE:: cDCarrFish(:)
	    Real,ALLOCATABLE:: cCDFishRef(:)
	    Real,ALLOCATABLE:: cNDFishRef(:)
	    Real,ALLOCATABLE:: cPDFishRef(:)
	    Real :: fishMinBiomass
    !3)Condicoes Iniciais 
	    Real,ALLOCATABLE:: sDFiAd0(:)     !Piscivoro jovem (mgDW/m²)
        Real,ALLOCATABLE:: sDFiJv0(:)     !Piscivoro jovem (mgDW/m²)
    !4)Fluxos
        Real:: wDConsFiAd,wCConsFiAd,wNConsFiAd,wPConsFiAd
        Real:: wDConsFiJv,wCConsFiJv,wNConsFiJv,wPConsFiJv
	    Real,ALLOCATABLE:: wDMigrFiAd(:)   ,wCMigrFiAd(:)   ,wNMigrFiAd(:)   ,wPMigrFiAd(:)
        Real,ALLOCATABLE:: wDMigrFiJv(:)   ,wCMigrFiJv(:)   ,wNMigrFiJv(:)   ,wPMigrFiJv(:)
	    Real,ALLOCATABLE:: wDReprFish(:)     ,wCReprFish(:)     ,wNReprFish(:)     ,wPReprFish(:)
	    Real,ALLOCATABLE:: wDAgeFish(:)      ,wCAgeFish(:)      ,wNAgeFish(:)      ,wPAgeFish(:) 
        Real,ALLOCATABLE:: wDAssFiAd(:,:,:),wCAssFiAd(:,:,:),wNAssFiAd(:,:,:),wPAssFiAd(:,:,:)
        Real,ALLOCATABLE:: wDAssFiJv(:,:,:),wCAssFiJv(:,:,:),wNAssFiJv(:,:,:),wPAssFiJv(:,:,:)
	    Real,ALLOCATABLE:: wDGrzFiAdFiAd(:,:,:,:),wCGrzFiAdFiAd(:,:,:,:),wNGrzFiAdFiAd(:,:,:,:),wPGrzFiAdFiAd(:,:,:,:)
        Real,ALLOCATABLE:: wDGrzFiAdFiJv(:,:,:,:),wCGrzFiAdFiJv(:,:,:,:),wNGrzFiAdFiJv(:,:,:,:),wPGrzFiAdFiJv(:,:,:,:)
	    Real,ALLOCATABLE:: wDGrzFiAdZoo(:,:,:,:) ,wCGrzFiAdZoo(:,:,:,:) ,wNGrzFiAdZoo(:,:,:,:) ,wPGrzFiAdZoo(:,:,:,:)
        Real,ALLOCATABLE:: wDGrzFiJvZoo(:,:,:,:) ,wCGrzFiJvZoo(:,:,:,:) ,wNGrzFiJvZoo(:,:,:,:) ,wPGrzFiJvZoo(:,:,:,:)
	    Real,ALLOCATABLE:: wDGrzFiAdBent(:,:,:,:) ,wCGrzFiAdBent(:,:,:,:) ,wNGrzFiAdBent(:,:,:,:) ,wPGrzFiAdBent(:,:,:,:)
        Real,ALLOCATABLE:: wDGrzFiJvBent(:,:,:,:) ,wCGrzFiJvBent(:,:,:,:) ,wNGrzFiJvBent(:,:,:,:) ,wPGrzFiJvBent(:,:,:,:)
	    Real,ALLOCATABLE:: wDGrzFiAdPhyt(:,:,:,:),wCGrzFiAdPhyt(:,:,:,:),wNGrzFiAdPhyt(:,:,:,:),wPGrzFiAdPhyt(:,:,:,:)
        Real,ALLOCATABLE:: wDGrzFiJvPhyt(:,:,:,:),wCGrzFiJvPhyt(:,:,:,:),wNGrzFiJvPhyt(:,:,:,:),wPGrzFiJvPhyt(:,:,:,:)
	    !Real,ALLOCATABLE:: tDGrzFishBac(:,:) ,tCGrzFishBac(:,:) ,tNGrzFishBac(:,:) ,tPGrzFishBac(:,:)
	    !Real,ALLOCATABLE:: tDGrzFishBent(:,:),tCGrzFishBent(:,:),tNGrzFishBent(:,:),tPGrzFishBent(:,:)
	    Real,ALLOCATABLE:: wDGrzFiAdPom(:,:,:,:) ,wCGrzFiAdPom(:,:,:,:) ,wNGrzFiAdPom(:,:,:,:) ,wPGrzFiAdPom(:,:,:,:)
        Real,ALLOCATABLE:: wDGrzFiJvPom(:,:,:,:) ,wCGrzFiJvPom(:,:,:,:) ,wNGrzFiJvPom(:,:,:,:) ,wPGrzFiJvPom(:,:,:,:)
	    Real,ALLOCATABLE:: wDHarvFiAd(:)     ,wCHarvFiAd(:)     ,wNHarvFiAd(:)     ,wPHarvFiAd(:)
        Real,ALLOCATABLE:: wDHarvFiJv(:)     ,wCHarvFiJv(:)     ,wNHarvFiJv(:)     ,wPHarvFiJv(:)
	    Real,ALLOCATABLE:: wDRespFiAd(:,:,:)     ,wCRespFiAd(:,:,:)
        Real,ALLOCATABLE:: wDRespFiJv(:,:,:)     ,wCRespFiJv(:,:,:)
	    Real,ALLOCATABLE:: wDEgesFiAd(:,:,:)     ,wCEgesFiAd(:,:,:)     ,wNEgesFiAd(:,:,:)     ,wPEgesFiAd(:,:,:)
        Real,ALLOCATABLE:: wDEgesFiJv(:,:,:)     ,wCEgesFiJv(:,:,:)     ,wNEgesFiJv(:,:,:)     ,wPEgesFiJv(:,:,:)
	    Real,ALLOCATABLE:: wDMortFiAd(:,:,:)     ,wCMortFiAd(:,:,:)     ,wNMortFiAd(:,:,:)     ,wPMortFiAd(:,:,:)
        Real,ALLOCATABLE:: wDMortFiJv(:,:,:)     ,wCMortFiJv(:,:,:)     ,wNMortFiJv(:,:,:)     ,wPMortFiJv(:,:,:)
	    Real,ALLOCATABLE:: wNExcrFiAd(:,:,:)      ,wPExcrFiAd(:,:,:)
        Real,ALLOCATABLE:: wNExcrFiJv(:,:,:)      ,wPExcrFiJv(:,:,:)
        Real,ALLOCATABLE:: wDMessFiAdOM(:,:,:),   wCMessFiAdOM(:,:,:),   wNMessFiAdOM(:,:,:),   wPMessFiAdOM(:,:,:)
        Real,ALLOCATABLE:: wDMessFiJvOM(:,:,:),   wCMessFiJvOM(:,:,:),   wNMessFiJvOM(:,:,:),   wPMessFiJvOM(:,:,:)
    !5)Funcoes Auxiliares
	    Real aDFish
	    Real uFunTmFish
        Real:: oDFoodFiAd,oCFoodFiAd,oNFoodFiAd,oPFoodFiAd
        Real:: oDFoodFiJv,oCFoodFiJv,oNFoodFiJv,oPFoodFiJv
	    Real:: oDFoodFiAdFish,oDFoodFiAdZoo,oDFoodFiAdPhyt,oDFoodFiAdBent,oDFoodFiAdPom
        Real:: oDFoodFiJvZoo,oDFoodFiJvPhyt,oDFoodFiJvBent,oDFoodFiJvPom
	    Real:: oDFoodFish
	    Real:: ukDIncrFiAd,tDEnvFiAd,aDSatFiAd
        Real:: ukDIncrFiJv,tDEnvFiJv,aDSatFiJv
        Real:: rCDFoodFiAd,rNDFoodFiAd,rPDFoodFiAd
        Real:: rCDFoodFiJv,rNDFoodFiJv,rPDFoodFiJv
        Real:: afCAssFiAd ,afNAssFiAd ,afPAssFiAd
        Real:: afCAssFiJv ,afNAssFiJv ,afPAssFiJv
	    Real,ALLOCATABLE:: fDGrzFishFish(:,:) ,fDGrzFishZoo(:,:), fDGrzFishPhyt(:,:), fDGrzFishBac(:,:), fDGrzFishBent(:,:), fDGrzFishPom(:,:)
    !7)Variaveis de Estado	
	    Real,ALLOCATABLE:: sDFiAd(:,:,:),   sCFiAd(:,:,:),   sNFiAd(:,:,:),   sPFiAd(:,:,:)        !Peixes Piscivoros (g/m²)
        Real,ALLOCATABLE:: sDFiJv(:,:,:),   sCFiJv(:,:,:),   sNFiJv(:,:,:),   sPFiJv(:,:,:)        !Peixes Piscivoros (g/m²)
	    Real,ALLOCATABLE:: sDFiAdP(:,:,:),  sCFiAdP(:,:,:),  sNFiAdP(:,:,:),  sPFiAdP(:,:,:)       !Peixes Piscivoros (g/m²)
        Real,ALLOCATABLE:: sDFiJvP(:,:,:),  sCFiJvP(:,:,:),  sNFiJvP(:,:,:),  sPFiJvP(:,:,:)       !Peixes Piscivoros (g/m²)
	    Real,ALLOCATABLE:: rCDFiAd(:,:,:),  rNDFiAd(:,:,:),  rPDFiAd(:,:,:)
        Real,ALLOCATABLE:: rCDFiJv(:,:,:),  rNDFiJv(:,:,:),  rPDFiJv(:,:,:)

    !------------------------------------------------------------------------------------------------------------------!	
    !  A B I O T I C O                                                                                                 !	
    !------------------------------------------------------------------------------------------------------------------!	
    !------------------------------------------------------------------------------------------------------------------!	
    !  C O N F I G U R A C O E S     I N I C I A I S                                                                   !	
    !------------------------------------------------------------------------------------------------------------------!	
        INTEGER ndom,npom
	    INTEGER domg
	    INTEGER pomg
	    INTEGER:: LABIL=1
	    INTEGER:: REFRATARIA=2
	    Real :: pomMinConc = 0.001
	    Real :: domMinConc = 0.001	
    
    
    
	    INTEGER :: nsed
	    INTEGER ConstDepth  !A profundidade é constante?    Sim=0, Não=1
    

    !------------------------------------------------------------------------------------------------------------------!	
    !  P A R A M E T R O S                                                                                             !	
    !------------------------------------------------------------------------------------------------------------------!	
    !Link BIOTICO<->ABIOTICO
	    !Particulate Organic Matter (POM)
	    Real,ALLOCATABLE:: cCDPomRef(:)
	    Real,ALLOCATABLE:: cNDPomRef(:)
	    Real,ALLOCATABLE:: cPDPomRef(:)
	    !Dissolved Organic Matter (DOM)
	    Real,ALLOCATABLE:: cCDDomRef(:)
	    Real,ALLOCATABLE:: cNDDomRef(:)
	    Real,ALLOCATABLE:: cPDDomRef(:)
	    Real,ALLOCATABLE:: cSiDDomRef(:)
	    !Real:: fDExdDomPhyt
	    Real:: fDLisDomPhyt
	    Real:: fDMessDomZoo
	    Real:: fDExcDomZoo
	    Real:: fDMessDomBen
	    Real:: fDEgesDomZoo
	    Real:: fDMortDomZoo
	    Real:: fDExcDomBen
	    Real:: fDEgesDomBen
	    Real:: fDMortDomBen
	    Real:: fDMessDomFish
	    Real:: fDExcDomFish
	    Real:: fDEgesDomFish
	    Real:: fDMortDomFish
	    Real:: fDMortDomMac
	    Real:: fDAssDomBird
        !Sedimento
	    Real fLutum						!Fração de lodo na matéria inorganica (g/g)
	    Real fLutumRef					!Fração de lodo na matéria inorganica (g/g)
	    Real cRhoOM						!Densidade da matéria organica no solido (g/m³ de sólido)
	    Real cRhoIM						!Densidade da matéria inorganica no solido (g/m³ de sólido) (diferentes para argila (1.95x10^6) e areia(1.7x10^6))
	    Real cDepthS					!Profundidade do topo de camada de sedimento (m)
        !Infiltracao
	    Real cQinf                        !taxa de infiltração (mm/d)
	    Real cPBackLoad                   !Carga direta de P na superficie da água (mgP/l) - deposição atmosferica na holanda
	    Real cNBackLoad                   !Carga direta de N na superficie da água (mgN/l) - deposição atmosferica na holanda
	    Real cPO4Ground                   !Concentração de PO4 em camadas profundas (mgP/l) 
	    Real cNH4Ground                   !Concentração de NH4 em camadas profundas (mgN/l)
	    Real cNO3Ground                   !Concentração de NO3 em camadas profundas (mgN/l)
	    Real cDicGround                   !Concentração de Dic em camadas profundas (mgN/l)
	    Real,ALLOCATABLE:: cDomGround(:) !gdom !Concentração de Dom em camadas profundas (mgN/l)
	    !Erosao
	    Real cDErosTot                     !Fluxo total de erosão (g/m²/d)
	    Real fSedErosIM                    !fração de materia inorganica erodida instantaneamente (-)
	    Real fDOrgSoil                     !fração de matéria organica no solo (-)
	    !Real,ALLOCATABLE:: fDErosPom(:)  !gpom !fracao refrataria e particulada (-)
    
        !Resuspension
        Integer:: ResuspMethodFlag                      !< Resuspension Method flag: ResuspMethodFlag = 0 (No resuspension); ResuspMethodFlag = 1 (PCLake); ResuspMethodFlag = 2 (Garcia and Parker)
	    Real caRes						!coeficientes ALFA da equação de regressão 
	    Real cbRes						!coeficientes BETA da equação de regressão 
	    Real cSuspRef					!função de referencia da materia suspensa (-)
	    Real cSuspMin					!minimo valor da função logistica
	    Real cSuspMax					!maximo valor da função logistica
	    Real cSuspSlope					!inclinação da função logistica
	    Real hDepthSusp					!valor de meia saturaçao da profundidade para a função logista
	    Real cFetchRef					!fetch de referencia
	    Real cFetch						!fetch na lagoa Mangueira
	    Real kTurbFish					!resuspensao relativa aos peixes adultos (g de sedimento/g de peixe)/dia
	    Real kVegResus					!coeficiente de atenuação da resuspenção devido a vegetação (m²/gD)
	    Real kResusPhytMax				!resuspensao maxima do fitoplancton (/d)
	    Real cResusPhytExp				!parametro exponencial para resuspensao do fitoplancton (d.m²/gD)
	    Real cThetaResus
    
        !Sedimentation
        Integer:: SetMethodFlag                         !< Sedimentation Method flag: SetMethodFlag = 0 (No resuspension); SetMethodFlag = 1 (PCLake); SetMethodFlag = 2 (Stokes)
        Real aFunTauSet
        Real uCorVSetIM
	    Real cThetaSet					!parametro de temperatura para sedimentação (1/e^oC)
	    Real cVSetPom					!velocidade de sedimentação maxima da materia organica particulada (m/day)
	    Real cVSetIM					!velocidade de sedimentação maxima da materia inorganica inerte (m/day)
	    Real ppOM						!densidade da OM
	    Real ppIM						!densidade da IM 
	    Real diOM						!diametro da OM (m)
	    Real diIM						!diametro da IM (m)
        !Hidrolise - POM -> DOM
	    Real cThetaHidW
	    Real hBacHid
	    Real,ALLOCATABLE:: kHidPomW(:) !gpom
	    Real cThetaHidS
	    Real,ALLOCATABLE:: kHidPomS(:) !gpom
        !Mineração - IF InclBac=0 - DOM->CO2,NH4,PO4,SiO2
	    Real cThetaMinAerW,cThetaMinAnaerW                     !Constante de temperatura exponecial de mineração na água
	    Real,ALLOCATABLE:: kMinAerDomW(:),kMinAnaerDomW(:) !gdom	!taxa de mineralizacao de DOM na água
	    Real cThetaMinAerS,cThetaMinAnaerS                     !Constante de temperatura exponecial de mineração no Sed
	    Real,ALLOCATABLE:: kMinAerDomS(:),kMinAnaerDomS(:) !gdom	!taxa de mineralizacao de DOM no Sed
        !DBO e Coliformes
	    Real Kdbo1
	    Real Kdbo3
	    Real Kcol
        !Desnitrificação
	    Real hNO3Denit					!meia-saturação quadratica da conc. de NO3 para desnitrificação
	    Real NO3PerC					!moles de NO3 denitrificados por mol de C mineralizado
        !Nitrificação
	    Real cThetaNitr                 !Constante de temperatura exponecial de nitrificação
	    Real kNitrW						!taxa de nitrificação constante na água (1/d)
	    Real kNitrS                     !taxa de nitrificação constante no Sed (1/d)
	    Real hO2Nitr					!meia-saturação quadratica da conc. de O2 para Nitrificação mgO2/l
	    Real O2PerNH4					!moléculas de O2 usadas por mol de NH4 nitrificado (-)
        !Adsorção do P
	    Real cRelPAdsD					!P máximo absorvido per g matéria organica (gP/gD)
	    Real cRelPAdsFe					!P máximo absorvido per g Fe (gP/gFe)
	    Real cRelPAdsAl					!P máximo absorvido per g matéria organica (gP/gAl)
	    Real fFeDIM						!Conteúdo de Fe na matéria inorganica (gFe/gD)
	    Real fAlDIM						!Conteúdo de Al na matéria inorganica (gAl/gD)
	    Real fRedMax					!fator de redução máxima da afinidade de absorção P
	    Real cKPAdsOx					!afinidade de absorção P em condições oxidadas
	    Real kPSorp						!taxa de absorção do fósforo (1/d)
	    Real hAdsP
        !Imobilizacao do P
	    Real kPChemPO4					!taxa constante (1/d)
	    Real cPO4Max					!Conc. Maxima de PO4 no sedimento
        !Difusao Agua-Sed
        Real fDepthDifS					!Fraction of Sediment's depth with diffusion (-) 
	    Real cThetaDif					!coeficiente de temperatura devido a difusao (1/e^oC)
	    Real,ALLOCATABLE:: kPDifDom(:)					!constante de difusao das moleculas de Dom (m²/d)
	    Real kCDifDic					!constante de difusao das moleculas de Dic (m²/d)
	    Real kPDifPO4					!constante de difusao das moleculas de PO4 (m²/d)
	    Real kNDifNO3					!constante de difusao das moleculas de NO3 (m²/d)
	    Real kNDifNH4					!constante de difusao das moleculas de NH4 (m²/d)
 	    Real kO2Dif						!Constante de difusão das moleculas de O2 (m²/d)
 	    Real cTurbDifNut				!fator de bioturbação para difusão
 	    Real cturbDifO2					!Fator de bioturbação para difusao de O2 (-)
        !Difusao Agua-Ar
	    !O2
	    Real cThetaReaer
	    Real uO2Sat
	    Real kReaer
        !CO2
	    Real pCO2
	    !Geral
	    Real cCPerDW					!Conteúdo de C na matéria organica (gC/gDW)
        Real O2PerNO3                   !Conteúdo de O2 em NO3 (gC/gDW)
	    Real molO2molC					!razão de peso molecular (gO2/gC)
	    Real molO2molN					!razão de peso molecular (gO2/gN)
	    Real molNmolC					!razão de peso molecular (gN/gC)
    !------------------------------------------------------------------------------------------------------------------!	
    !  C O N D I C A O     I N I C I A L                                                                               !	
    !------------------------------------------------------------------------------------------------------------------!	

        !Sedimento
	    Real fDTotS0     !fração inicial de peso seco no sedimento (g solido/g sedimento)
	    Real fDOrgS0     !fração inicial de materia organica no sedimento (g residuos de MO/g solido)
	    Real fDLabilS0   !fração inicial de detritos na matéria orgânica (g/g)
	    Real bRhoSolidS0 !Densidade inicial do material sólido (g/m³ de sólido)
	    Real bPorS       !porosidade (m³ água/m³ de sedimento)
	    Real bPorCorS    !Porosidade do sedimento corrigida pela tortuosidade
	    Real bRhoTotS0   !Densidade aparente do sedimento (g solido/m³ de sedimento)
	    Real bDTotS0     !total de peso seco inicial no topo de camada (gD/m²)
	    !Materia Organica
	    Real sCFW0                                  !CF (mgDW/l)
	    !Materia Inorganica
	
	    Real fPInorgS0        !fração inicial de P inorganico no sedimento (gP/gD)
	    Real fPAdsS0          !fração inicial absorvida de P inorganico no sedimento (-)
        Real sDicS0    !NH4 na água (mgN/l)
        Real spHW0,spHS0      !NH4 na água (mgN/l)
        Real sAlkW0,sAlkS0    !NH4 na água (mgN/l)
        !Gases
	    Real sCH4W0,sCH4S0
    !------------------------------------------------------------------------------------------------------------------!	
    ! F L U X O S   E   V A R S   A U X I L I A R E S                                                                  !	
    !------------------------------------------------------------------------------------------------------------------!	
	    !Sedimento
	    Real akO2DifCor,tSOD,aDepthOxySed,afOxySed
	    Real aO2LimS
	    Real,ALLOCATABLE :: aDepthS(:)
	    Real,ALLOCATABLE :: fDepthS(:)
	    !Infiltração
	    Real tPInfPO4W,tNInfNH4W,tNInfNO3W,tCInfDicW
	    Real tPInfPO4S,tNInfNH4S,tNInfNO3S,tCInfDicS
	    Real,ALLOCATABLE:: tDInfDomW(:),tCInfDomW(:),tNInfDomW(:),tPInfDomW(:)
	    Real,ALLOCATABLE:: tDInfDomS(:),tCInfDomS(:),tNInfDomS(:),tPInfDomS(:)
	    !Erosão
	    Real uDErosIM,uDErosIMW
        Real,ALLOCATABLE :: uDErosIMS(:)
	    Real,ALLOCATABLE :: tDErosPom(:,:,:),tCErosPom(:,:),tNErosPom(:,:),tPErosPom(:,:)
	    !Enterro
	    Real vDeltaS
	    Real tDIMS
	    Real,ALLOCATABLE :: tDPomS(:)
	    Real,ALLOCATABLE :: tDBurPom(:),tCBurPom(:),tNBurPom(:),tPBurPom(:)
	    Real,ALLOCATABLE :: tDBurDom(:),tCBurDom(:),tNBurDom(:),tPBurDom(:)
	    Real tDBurTot,tDBurIM,tDBurOM,tCBurDic,tPBurPO4,tPBurAIM,tNBurNH4,tNBurNO3
	    Real vDeltaW,aRelDeltaW
	    !Resuspensao
	    Real aFunDimSusp,tDResusTauDead,tDTurbFish,aFunVegResus,tDResusDead
	    Real tDResusPhytTot
	    Real,ALLOCATABLE :: tDResusTot(:),tDResusIM(:),tPResusAIM(:)
	    Real,ALLOCATABLE :: tPResusPO4(:),tNResusNH4(:),tNResusNO3(:),tCResusDic(:)
	    Real,ALLOCATABLE :: tDResusPom(:,:),tCResusPom(:,:),tPResusPom(:,:),tNResusPom(:,:)
	    Real,ALLOCATABLE :: tDResusDom(:,:),tCResusDom(:,:),tNResusDom(:,:),tPResusDom(:,:)
	    !Sedimentação
	    Real mu
	    Real uFunTmSet,uCorVSetDet
	    Real,ALLOCATABLE :: StokesVelIM(:),StokesVelPom(:,:)
	    Real,ALLOCATABLE :: tDSetIM(:,:),tPSetAIM(:,:)
	    Real,ALLOCATABLE :: tDSetPom(:,:,:),tCSetPom(:,:,:),tNSetPom(:,:,:),tPSetPom(:,:,:)
	    !Hidrolise
	    Real,ALLOCATABLE :: aTLimHidPom(:),aO2LimHidPom(:),aBacLimHidPom(:),apHLimHidPom(:)
	    Real,ALLOCATABLE :: kCorHidPom(:)	
	    Real,ALLOCATABLE :: wDHidPomW(:,:,:),wCHidPomW(:,:),wNHidPomW(:,:),wPHidPomW(:,:)
	    Real,ALLOCATABLE :: aTLimHidPomS(:)
	    Real,ALLOCATABLE :: tDHidPomS(:,:),tCHidPomS(:,:),tNHidPomS(:,:),tPHidPomS(:,:)    
        Real:: aO2LimHidPomS,kCorHidPomS
	    !Mineração - IF InclBac=0
	    Real aCorO2BOD 
	    Real,ALLOCATABLE :: wDMinAerW(:,:,:),wCMinAerW(:,:),wNMinAerW(:,:),wPMinAerW(:,:),wSiMinAerW(:,:)
        Real,ALLOCATABLE :: wCMinAer2CO2W(:,:,:),wCMinAer2CH4W(:,:)
	    Real,ALLOCATABLE :: wDMinAnaerW(:,:,:),wCMinAnaerW(:,:),wNMinAnaerW(:,:),wPMinAnaerW(:,:),wSiMinAnaerW(:,:)
        Real,ALLOCATABLE :: wCMinAnaer2CO2W(:,:),wCMinAnaer2CH4W(:,:)
	    Real,ALLOCATABLE :: aTLimMinAer(:),aO2LimMinAer(:),apHLimMinAer(:),aBacLimMinAer(:)
	    Real,ALLOCATABLE :: aTLimMinAnaer(:),aO2LimMinAnaer(:),apHLimMinAnaer(:),aBacLimMinAnaer(:)
	    Real,ALLOCATABLE :: aTLimMinDomS(:)
	    Real,ALLOCATABLE :: tDMinAerS(:,:),tCMinAerS(:,:),tNMinAerS(:,:),tPMinAerS(:,:)
        Real,ALLOCATABLE :: tDMinAnaerS(:,:),tCMinAnaerS(:,:),tNMinAnaerS(:,:),tPMinAnaerS(:,:)
	    Real,ALLOCATABLE :: tCMinAer2CO2S(:,:),tCMinAer2CH4S(:,:)
        Real,ALLOCATABLE :: tCMinAnaer2CO2S(:,:),tCMinAnaer2CH4S(:,:)
        Real:: tO2MinDetS,tDMinOxyDetS,aTLimMinAerS,aO2LimMinAerS,uFunTmMinS,aTLimMinAnaerS,aO2LimMinAnaerS
    
	    !Desnitrificação
	    Real wDDenitW,tDDenitS
        Real,ALLOCATABLE :: wNDenitW(:,:)
	    Real,ALLOCATABLE :: tNDenitS(:)    
	    !Nitrificação
        Real uFunTmNitr
	    Real aCorO2NitrW
	    Real,ALLOCATABLE :: wNNitrW(:,:)
	    Real,ALLOCATABLE :: tNNitrS(:)
        Real,ALLOCATABLE :: tO2NitrS(:)
	    Real aCorO2NitrS
	    !Adsorção do Fósforo
	    Real aPAdsMaxW,aKPAdsW,aPIsoAdsW,aPEqIMW
	    Real aPAdsMaxS,aKPAdsS,aPIsoAdsS,aPEqIMS
	    Real,ALLOCATABLE :: wPSorpIMW(:,:)
	    Real,ALLOCATABLE :: tPSorpIMS(:)
	    !Imobilização do Fósforo
	    Real tPChemPO4
	    !Difusao Agua-Sed
        Real uFunTmDif
	    Real aDepthDif
        Real,ALLOCATABLE:: tPDifPO4(:),tNDifNH4(:),tNDifNO3(:),tCDifDic(:),tO2Dif(:)
	    Real,ALLOCATABLE:: tDDifDom(:,:),tCDifDom(:,:),tNDifDom(:,:),tPDifDom(:,:)
	    !Difusao Agua-Ar
	    !O2
	    Real,ALLOCATABLE:: tO2Reaer(:)
	    !CO2
	    Real kCO2Exch,pkHenry,kHenry,uCO2Sat
	    Real,ALLOCATABLE:: tCO2Reaer(:)
    !------------------------------------------------------------------------------------------------------------------!	
    !  D E R I V A D A S                                                                                               !	
    !------------------------------------------------------------------------------------------------------------------!	
	    !Real,ALLOCATABLE:: dDDomW(:),  dCDomW(:),  dNDomW(:),  dPDomW(:), dSiDomW(:)
	    !Real,ALLOCATABLE:: dDDomS(:),  dCDomS(:),  dNDomS(:),  dPDomS(:), dSiDomS(:)  
	    !Real,ALLOCATABLE:: dDPomW(:),  dCPomW(:),  dNPomW(:),  dPPomW(:), dSiPomW(:)
	    !Real,ALLOCATABLE:: dDPomS(:),  dCPomS(:),  dNPomS(:),  dPPomS(:), dSiPomS(:)
	    !Real dDBOW 
	    !Real dCFW
	    !Real dDIMW,dPO4W,dNH4W,dPAIMW,dNO3W,dSiO2W  
	    !Real dDIMS,dPO4S,dNH4S,dPAIMS,dNO3S                                                         
	    !Real dO2W, dDicW,dCH4W
	    !Real dDicS,dCH4S
    !------------------------------------------------------------------------------------------------------------------!	
    !  V A R I A V E I S    D E    E S T A D O                                                                         !	
    !------------------------------------------------------------------------------------------------------------------!
    !Variáveis de Estado (No tempo ATUAL)
	    !Bioticos
        Real,ALLOCATABLE:: sCFW(:,:)
	    !Gases
	    Real,ALLOCATABLE:: sCH4W(:,:)
	    Real,ALLOCATABLE::           sDicS(:),  sCH4S(:)

    !Variáveis de Estado (No tempo PASSADO)
	    !Materia Organica
        Real,ALLOCATABLE:: sCFWP(:,:)
	    !Gases
	    Real,ALLOCATABLE:: sCH4WP(:,:) 
	    Real,ALLOCATABLE:: sDicSP(:),  sCH4SP(:)

    !Razoes de nutrientes internos Abióticas
	    Real,ALLOCATABLE:: rCDPomW(:,:,:), rNDPomW(:,:,:), rPDPomW(:,:,:), rSiDPomW(:,:,:)
	    Real,ALLOCATABLE:: rCDDomW(:,:,:), rNDDomW(:,:,:), rPDDomW(:,:,:), rSiDDomW(:,:,:)       ! Razões de nutrientes na água
	    Real,ALLOCATABLE:: rCDPomS(:,:),   rNDPomS(:,:),   rPDPomS(:,:),   rSiDPomS(:,:)
	    Real,ALLOCATABLE:: rCDDomS(:,:),   rNDDomS(:,:),   rPDDomS(:,:),   rSiDDomS(:,:)         ! Razões de nutrientes na água
	    Real,ALLOCATABLE::                   rPDIMW(:,:)
	    Real,ALLOCATABLE::                   rPDIMS(:)

    !Variaveis Secundarias	
        Real,ALLOCATABLE:: oDPhytW(:,:),oDOMW(:,:),oDSestW(:,:)                                                    ! Variaveis de peso seco na água
	    Real,ALLOCATABLE:: oCPhytW(:,:),oCOMW(:,:),oCSestW(:,:)
	    Real,ALLOCATABLE:: oNPhytW(:,:),oNOMW(:,:),oNSestW(:,:),oNTotW(:,:),oNkjW(:,:),oNDissW(:,:)          ! Variaveis de N na água
	    Real,ALLOCATABLE:: oPPhytW(:,:),oPOMW(:,:),oPSestW(:,:),oPTotW(:,:),oPInorgW(:,:)                      ! Variaveis de P na água
		
	    Real,ALLOCATABLE:: aDTotS(:),aDPhytS(:),aRhoTotS(:),aRhoSolidS(:),afDTotS(:),afDOrgS(:),afPOMS(:,:),afPOMTotS(:)                                    ! Variaveis abiotica no sedimento (peso seco)
	    Real,ALLOCATABLE:: aCPhytS(:)
	    Real,ALLOCATABLE:: aNDissS(:),aNPhytS(:),aNkjAvailS(:),aNkjS(:),aNTotAvailS(:),aNTotS(:),afNInorS(:),afNTotS(:),oNO3S(:),oNH4S(:),oNDissS(:) ! Variaveis abiotica no sedimento  (N)
	    Real,ALLOCATABLE:: aPInorgS(:),aPPhytS(:),aPTotAvailS(:),aPTotS(:)    ,afPInorgS(:)  ,afPO4S(:) ,oPO4S(:)                                                      ! Variaveis abiotica no sedimento  (P)
	    Real,ALLOCATABLE:: oDDomS(:,:)
	    Real,ALLOCATABLE:: oDicS(:), oHCO3S(:) ,oH2CO3S(:),  oCO3S(:)	
	    Real,ALLOCATABLE:: sAlkW(:,:), spHW(:,:),sH2CO3W(:,:),sHCO3W(:,:),sCO3W(:,:),pCO2W(:,:)
        Real,ALLOCATABLE:: sAlkS(:),   spHS(:),  sH2CO3S(:),  sHCO3S(:),  sCO3S(:)
      
        Real,ALLOCATABLE:: sO2S(:) 
        Real,ALLOCATABLE:: sGPPW(:,:),sRespW(:,:),sNEPW(:,:)
        
    contains
        procedure :: initializeLimnoParam
    end type
    
    contains
    
    subroutine initializeLimnoParam(this, wqConfiguration)     
        
        Integer:: i,j,k,pindex
        type(WaterQualityConfiguration), pointer :: wqConfiguration
        type(WaterQualityParameter), dimension(:), pointer :: wqParam
        type(WaterQualityGroup), dimension(:), pointer :: wqgroups
        type(FoodMatrixValue), dimension(:), pointer :: wqFoodMatrix
        real(c_double), dimension(:), pointer :: wqvalues
        integer(c_int) :: groupQuantity
        character(len=255):: text
        class(LimnologyParam) :: this
        
        ! Reading water quality parameters 
        call c_f_pointer(wqConfiguration%parameters, wqParam, [wqConfiguration%numberOfParameters])
        
        Do i = 1, wqConfiguration%numberOfParameters
            
            text = trim(wqParam(i)%name)
            
            !Structure Flags
            If (trim(text) == 'heatBudgetModeling') Then 
                this%iTempW = wqParam(i)%value
            ElseIf (trim(text) == 'salinityModeling') Then 
                this%iSal = wqParam(i)%value
            ElseIf (trim(text) == 'limnologyModeling') Then 
                this%iLimno = wqParam(i)%value
            ElseIf (trim(text) == 'iOxyBOD') Then 
                this%iOxyBOD = wqParam(i)%value
            ElseIf (trim(text) == 'iMetabO2') Then 
                this%iMetabO2 = wqParam(i)%value
            ElseIf (trim(text) == 'iCarbon') Then 
                this%iCarbon = wqParam(i)%value
            ElseIf (trim(text) == 'iGHG') Then 
                this%iGHG = wqParam(i)%value
            ElseIf (trim(text) == 'iMetabC') Then 
                this%iMetabC = wqParam(i)%value
            ElseIf (trim(text) == 'bacterioplanktonModeling') Then 
                this%InclBac = wqParam(i)%value
            ElseIf (trim(text) == 'iMatInorg') Then 
                this%InclMatInorg = wqParam(i)%value
            ElseIf (trim(text) == 'iSuspIMDist') Then 
                this%InclMatInorgSplit = wqParam(i)%value
            ElseIf (trim(text) == 'organicMatterModeling') Then 
                this%InclMatOrg = wqParam(i)%value
            ElseIf (trim(text) == 'iMatOrgSplit') Then 
                this%InclMatOrgSplit = wqParam(i)%value
            ElseIf (trim(text) == 'iPomDomSplit') Then 
                this%InclPomDomSplit = wqParam(i)%value
            ElseIf (trim(text) == 'phytoplanktonModeling') Then 
                this%InclPhyt = wqParam(i)%value
            ElseIf (trim(text) == 'chalker') Then
                If (wqParam(i)%value==1) Then
                    this%LightMethodPhyt = 0
                EndIf
            ElseIf (trim(text) == 'klepper') Then
                If (wqParam(i)%value==1) Then
                    this%LightMethodPhyt = 1
                EndIf
            ElseIf (trim(text) == 'zooplanktonModeling') Then 
                this%InclZoo = wqParam(i)%value
            ElseIf (trim(text) == 'iBen') Then 
                this%InclBent = wqParam(i)%value
            ElseIf (trim(text) == 'iMac') Then 
                this%InclMacr = wqParam(i)%value
            ElseIf (trim(text) == 'fishModeling') Then 
                this%InclFish = wqParam(i)%value
            ElseIf (trim(text) == 'iFishMov') Then 
                this%iFishMov = wqParam(i)%value
            ElseIf (trim(text) == 'fishModelingRefuge') Then 
                this%iFishRefuge = wqParam(i)%value
            ElseIf (trim(text) == 'fishModelingStage') Then 
                this%iFishStage = wqParam(i)%value
            ElseIf (trim(text) == 'sedimentation') Then 
                this%iSed = wqParam(i)%value
            ElseIf (trim(text) == 'noResuspension') Then 
                If (wqParam(i)%value==1) Then
                    this%ResuspMethodFlag = 0
                EndIf
            ElseIf (trim(text) == 'pcLake1') Then 
                If (wqParam(i)%value==1) Then
                    this%ResuspMethodFlag = 1
                EndIf
            ElseIf (trim(text) == 'garciaParker') Then 
                If (wqParam(i)%value==1) Then
                    this%ResuspMethodFlag = 1 !Not implemented yet
                EndIf
            ElseIf (trim(text) == 'noSedimentation') Then 
                If (wqParam(i)%value==1) Then
                    this%SetMethodFlag = 0
                EndIf
            ElseIf (trim(text) == 'pcLake2') Then 
                If (wqParam(i)%value==1) Then
                    this%SetMethodFlag = 1
                EndIf
            ElseIf (trim(text) == 'stokes') Then 
                If (wqParam(i)%value==1) Then
                    this%SetMethodFlag = 0
                EndIf
            ElseIf (trim(text) == 'iAdaptativeTimeStep') Then 
                this%iAdaptativeTimeStep = wqParam(i)%value
            ElseIf (trim(text) == 'nsUpwind') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 0
                EndIf
            ElseIf (trim(text) == 'nsLaxWendroff') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 1
                EndIf
            ElseIf (trim(text) == 'nsSuperbee') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 2
                EndIf
            ElseIf (trim(text) == 'nsVanLeer') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 3
                EndIf
            ElseIf (trim(text) == 'nsMinMod') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 4
                EndIf
            ElseIf (trim(text) == 'nsMuscl') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 5
                EndIf
            ElseIf (trim(text) == 'nsSuperC') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 6
                EndIf
            ElseIf (trim(text) == 'nsUltimateQuickest') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 7
                EndIf
            ElseIf (trim(text) == 'nsHyperC') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 8
                EndIf
            ElseIf (trim(text) == 'nsOSPRE') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 9
                EndIf
            ElseIf (trim(text) == 'nsSPL13') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 10
                EndIf
            ElseIf (trim(text) == 'nsHCUI') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 11
                EndIf
            ElseIf (trim(text) == 'nsSmart') Then 
                If (wqParam(i)%value==1) Then
                    this%iTranspFlag = 12
                EndIf
                !Group Size
            ElseIf (trim(text) == 'aerobicGroup') Then 
                this%nbacaer = wqParam(i)%value
            ElseIf (trim(text) == 'anaerobicGroup') Then 
                this%nbacanaer = wqParam(i)%value
            ElseIf (trim(text) == 'dinoflagellatesGroup') Then 
                this%nphydino = wqParam(i)%value
            ElseIf (trim(text) == 'freshwaterCyanobacteriaGroup') Then 
                this%nphyfcyano = wqParam(i)%value
            ElseIf (trim(text) == 'marineCyanobacteriaGroup') Then 
                this%nphymcyano = wqParam(i)%value
            ElseIf (trim(text) == 'chlorophytesGroup') Then 
                this%nphychloro = wqParam(i)%value
            ElseIf (trim(text) == 'cryptophytesGroup') Then 
                this%nphycryto = wqParam(i)%value
            ElseIf (trim(text) == 'marineDiatomsGroup') Then 
                this%nphymdiato = wqParam(i)%value
            ElseIf (trim(text) == 'freshwaterDiatomsGroup') Then 
                this%nphyfdiato = wqParam(i)%value
            ElseIf (trim(text) == 'microZooGroup') Then 
                this%nzoomicro = wqParam(i)%value
            ElseIf (trim(text) == 'mesoZooGroup') Then 
                this%nzoomeso = wqParam(i)%value
            ElseIf (trim(text) == 'macroZooGroup') Then 
                this%nzoomacro = wqParam(i)%value
            ElseIf (trim(text) == 'microZoobenthosGroup') Then 
                this%nbenmicro = wqParam(i)%value
            ElseIf (trim(text) == 'mesoZoobenthosGroup') Then 
                this%nbenmeso = wqParam(i)%value
            ElseIf (trim(text) == 'macroZoobenthosGroup') Then 
                this%nbenmacro = wqParam(i)%value
            ElseIf (trim(text) == 'planktivorousGroup') Then 
                this%nfishplank = wqParam(i)%value
            ElseIf (trim(text) == 'omnivorousGroup') Then 
                this%nfishomniv = wqParam(i)%value
            ElseIf (trim(text) == 'piscivorousGroup') Then 
                this%nfishpisc = wqParam(i)%value
            ElseIf (trim(text) == 'elodeidGroup') Then 
                this%nmacElod = wqParam(i)%value
            ElseIf (trim(text) == 'charophytesGroup') Then 
                this%nmacCharo = wqParam(i)%value
            ElseIf (trim(text) == 'ceratophyllumGroup') Then 
                this%nmacCera = wqParam(i)%value
            ElseIf (trim(text) == 'lemnaceaeGroup') Then 
                this%nmacLemna = wqParam(i)%value
            ElseIf (trim(text) == 'nymphaeasGroup') Then 
                this%nmacNymp = wqParam(i)%value
            ElseIf (trim(text) == 'helophytesGroup') Then 
                this%nmacHelo = wqParam(i)%value
            ElseIf (trim(text) == 'generalGroup') Then 
                this%nmacGeral = wqParam(i)%value
            
            EndIf
            
        EndDo
        
        !Allocate arrays
        
        !Organic Matter   
        If (this%InclMatOrg == 1) Then
            If (this%InclMatOrgSplit == 1) Then
                If (this%InclPomDomSplit == 1) Then
                    this%npom = 2
                    this%ndom = 2
                    ALLOCATE (this%sDPomW0(this%ndom))
                    ALLOCATE (this%sDPomS0(this%ndom))
                    ALLOCATE (this%sDDomW0(this%ndom))
                    ALLOCATE (this%sDDomS0(this%ndom))
                Else
                    this%npom = 1
                    this%ndom = 1
                    ALLOCATE (this%sDPomW0(this%ndom))
                    ALLOCATE (this%sDPomS0(this%ndom))
                    ALLOCATE (this%sDDomW0(this%ndom))
                    ALLOCATE (this%sDDomS0(this%ndom))
                EndIf
            Else
                this%npom = 1
                this%ndom = 1
                ALLOCATE (this%sDPomW0(this%ndom))
                ALLOCATE (this%sDPomS0(this%ndom))
                ALLOCATE (this%sDDomW0(this%ndom))
                ALLOCATE (this%sDDomS0(this%ndom))
            EndIf   
            
            ALLOCATE (this%kPDifDom(this%ndom))
            ALLOCATE (this%cCDPomRef(this%npom))
            ALLOCATE (this%cNDPomRef(this%npom))
            ALLOCATE (this%cPDPomRef(this%npom))
            ALLOCATE (this%cExtSpPom(this%npom))
	        ALLOCATE (this%cDomGround(this%ndom))  
	        ALLOCATE (this%kHidPomW(this%npom))  
	        ALLOCATE (this%kHidPomS(this%npom))  
            ALLOCATE (this%kMinAerDomW(this%ndom))  
	        ALLOCATE (this%kMinAnaerDomW(this%ndom))  
	        ALLOCATE (this%kMinAerDomS(this%ndom))
            ALLOCATE (this%kMinAnaerDomS(this%ndom))
            ALLOCATE (this%cCDDomRef(this%ndom))
            ALLOCATE (this%cNDDomRef(this%ndom))
            ALLOCATE (this%cPDDomRef(this%ndom))
            ALLOCATE (this%cSiDDomRef(this%ndom))
            ALLOCATE (this%cExtSpDom(this%ndom))               
        Else
            this%npom = 1
            this%ndom = 1
            ALLOCATE (this%sDPomW0(this%ndom))
            ALLOCATE (this%sDPomS0(this%ndom))
            ALLOCATE (this%sDDomW0(this%ndom))
            ALLOCATE (this%sDDomS0(this%ndom))
            ALLOCATE (this%kPDifDom(this%ndom))
            ALLOCATE (this%cCDPomRef(this%npom))
            ALLOCATE (this%cNDPomRef(this%npom))
            ALLOCATE (this%cPDPomRef(this%npom))
            ALLOCATE (this%cExtSpPom(this%npom))
	        ALLOCATE (this%cDomGround(this%ndom))  
	        ALLOCATE (this%kHidPomW(this%npom))  
	        ALLOCATE (this%kHidPomS(this%npom))  
            ALLOCATE (this%kMinAerDomW(this%ndom))  
	        ALLOCATE (this%kMinAnaerDomW(this%ndom))  
	        ALLOCATE (this%kMinAerDomS(this%ndom))
            ALLOCATE (this%kMinAnaerDomS(this%ndom))
            ALLOCATE (this%cCDDomRef(this%ndom))
            ALLOCATE (this%cNDDomRef(this%ndom))
            ALLOCATE (this%cPDDomRef(this%ndom))
            ALLOCATE (this%cSiDDomRef(this%ndom))
            ALLOCATE (this%cExtSpDom(this%ndom))  
        EndIf
        
        !Bacterioplankton    
        If (this%InclBac==0) Then
            this%nbac = 0
        Else
            this%nbac = this%nbacaer + this%nbacanaer
            ALLOCATE (this%abac(this%nbac)) ! Bacterioplankton functional groups (= 1 -> Aerobic; = 2 -> Anaerobic)
            this%abac(1:this%nbacaer) = 1
            this%abac(this%nbacaer+1:this%nbac) = 2    
            ALLOCATE (this%sDBacW0(this%nbac))
            ALLOCATE (this%sDBacS0(this%nbac))
            ALLOCATE (this%cThetaBac(this%nbac))
            ALLOCATE (this%MuBac(this%nbac))
            ALLOCATE (this%hAssBac(this%nbac))
            ALLOCATE (this%kResbBac(this%nbac))
            ALLOCATE (this%fBac2CO2(this%nbac))
            ALLOCATE (this%hPO4UptBac(this%nbac))
            ALLOCATE (this%hNH4UptBac(this%nbac))
            ALLOCATE (this%kMortBac(this%nbac))
            ALLOCATE (this%hSedBac(this%nbac))
            ALLOCATE (this%cCDBacRef(this%nbac))
            ALLOCATE (this%cNDBacRef(this%nbac))
            ALLOCATE (this%cPDBacRef(this%nbac))
            ALLOCATE (this%hBacMin(this%nbac))
        EndIf
        
        !Phytoplankton    
        If (this%InclPhyt==0) Then
            this%nphy = 0
        Else
            this%nphy = this%nphydino + this%nphyfcyano + this%nphymcyano + this%nphychloro + this%nphycryto + this%nphymdiato + this%nphyfdiato
            ALLOCATE (this%aphy(this%nphy)) 
            pindex = 0
            Do i = 1,this%nphydino
                pindex = pindex + 1
                this%aphy(pindex) = 1
            EndDo
            Do i = 1,this%nphyfcyano
                pindex = pindex + 1
                this%aphy(pindex) = 2
            EndDo
            Do i = 1,this%nphymcyano
                pindex = pindex + 1
                this%aphy(pindex) = 3
            EndDo
            Do i = 1,this%nphychloro
                pindex = pindex + 1
                this%aphy(pindex) = 4
            EndDo
            Do i = 1,this%nphycryto
                pindex = pindex + 1
                this%aphy(pindex) = 5
            EndDo
            Do i = 1,this%nphymdiato
                pindex = pindex + 1
                this%aphy(pindex) = 6
            EndDo
            Do i = 1,this%nphyfdiato
                pindex = pindex + 1
                this%aphy(pindex) = 7
            EndDo
            ALLOCATE (this%sDPhytW0(this%nphy))
            ALLOCATE (this%sDPhytS0(this%nphy))
            ALLOCATE (this%cCDPhytOpt(this%nphy))
            ALLOCATE (this%cNDPhytOpt(this%nphy))
            ALLOCATE (this%cPDPhytOpt(this%nphy))
            ALLOCATE (this%cSiDPhytOpt(this%nphy))
            ALLOCATE (this%cExtSpPhyt(this%nphy))
            ALLOCATE (this%hLRefPhyt(this%nphy))
            ALLOCATE (this%cLOptRefPhyt(this%nphy))
            ALLOCATE (this%cMuMaxPhyt(this%nphy))
            ALLOCATE (this%cChDPhytMax(this%nphy))
            ALLOCATE (this%cChDPhytMin(this%nphy))
            ALLOCATE (this%hSiAssDiat(this%nphy))
            ALLOCATE (this%cSigTmPhyt(this%nphy))
            ALLOCATE (this%cTmOptPhyt(this%nphy))
            ALLOCATE (this%cVCUptMaxPhyt(this%nphy))
            ALLOCATE (this%cVNUptMaxPhyt(this%nphy))
            ALLOCATE (this%cVPUptMaxPhyt(this%nphy))
            ALLOCATE (this%cVSiUptMaxPhyt(this%nphy))
            ALLOCATE (this%cCDPhytMax(this%nphy))
            ALLOCATE (this%cCDPhytMin(this%nphy))
            ALLOCATE (this%cNDPhytMax(this%nphy))
            ALLOCATE (this%cNDPhytMin(this%nphy))
            ALLOCATE (this%cPDPhytMax(this%nphy))
            ALLOCATE (this%cPDPhytMin(this%nphy))
            ALLOCATE (this%cSiDPhytMax(this%nphy))
            ALLOCATE (this%cSiDPhytMin(this%nphy))
            ALLOCATE (this%cAffCUptPhyt(this%nphy))
            ALLOCATE (this%cAffNUptPhyt(this%nphy))
            ALLOCATE (this%cAffPUptPhyt(this%nphy))
            ALLOCATE (this%cAffSiUptPhyt(this%nphy))      
            ALLOCATE (this%kDRespbPhyt(this%nphy))
            ALLOCATE (this%alfap(this%nphy))
            ALLOCATE (this%cVSetPhyt(this%nphy))
            ALLOCATE (this%ppPhyt(this%nphy)) !SetMethodFlag = 1
            ALLOCATE (this%diPhyt(this%nphy)) !SetMethodFlag = 1
            ALLOCATE (this%cLisPhytW(this%nphy))
            ALLOCATE (this%cLisPhytS(this%nphy))        
        EndIf    
        
        !Zooplankton    
        If (this%InclZoo==0) Then
            this%nzoo = 0
        Else
            this%nzoo = this%nzoomicro + this%nzoomeso + this%nzoomacro
            ALLOCATE (this%azoo(this%nzoo))
            pindex = 0
            Do i = 1,this%nzoomicro
                pindex = pindex + 1
                this%azoo(pindex) = 1
            EndDo
            Do i = 1,this%nzoomeso
                pindex = pindex + 1
                this%azoo(pindex) = 2
            EndDo
            Do i = 1,this%nzoomeso
                pindex = pindex + 1
                this%azoo(pindex) = 3
            EndDo
            ALLOCATE (this%sDZoo0(this%nzoo))
            ALLOCATE (this%cCDZooRef(this%nzoo))
            ALLOCATE (this%cNDZooRef(this%nzoo))
            ALLOCATE (this%cPDZooRef(this%nzoo))
            ALLOCATE (this%fDAssZoo(this%nzoo))
            ALLOCATE (this%cSigTmZoo(this%nzoo))
            ALLOCATE (this%cTmOptZoo(this%nzoo))
            ALLOCATE (this%cDCarrZoo(this%nzoo))
            ALLOCATE (this%hFilt(this%nzoo))
            ALLOCATE (this%cFiltMax(this%nzoo))
            ALLOCATE (this%kDRespZoo(this%nzoo))
            ALLOCATE (this%kPelZoo(this%nzoo))
            ALLOCATE (this%kMortZoo(this%nzoo))            
        EndIf
        
        !Zoobenthos    
        If (this%InclBent==0) Then
            this%nben = 0
        Else
            this%nben = this%nbenmicro + this%nbenmeso + this%nbenmacro
            ALLOCATE (this%aben(this%nben))
            pindex = 0
            Do i = 1,this%nbenmicro
                pindex = pindex + 1
                this%aben(pindex) = 1
            EndDo
            Do i = 1,this%nbenmeso
                pindex = pindex + 1
                this%aben(pindex) = 2
            EndDo
            Do i = 1,this%nbenmacro
                pindex = pindex + 1
                this%aben(pindex) = 3
            EndDo
            ALLOCATE (this%sDBent0(this%nben))
            ALLOCATE (this%cSigTmBent(this%nben))
            ALLOCATE (this%cTmOptBent(this%nben))
            ALLOCATE (this%hDFoodBent(this%nben))
            ALLOCATE (this%kDAssBent(this%nben))
            ALLOCATE (this%cDCarrBent(this%nben))
            ALLOCATE (this%fDAssBent(this%nben))
            ALLOCATE (this%kDRespBent(this%nben))
            ALLOCATE (this%kPelBent(this%nben))
            ALLOCATE (this%kMortBent(this%nben))
            ALLOCATE (this%cCDBentRef(this%nben))
            ALLOCATE (this%cNDBentRef(this%nben))
            ALLOCATE (this%cPDBentRef(this%nben))

        EndIf
        
        !Fish    
        If (this%InclFish==0) Then
            this%nfish = 0
        Else
            this%nfish = this%nfishomniv + this%nfishplank + this%nfishpisc
            ALLOCATE (this%afish(this%nfish))
            pindex = 0
            Do i = 1,this%nfishpisc
                pindex = pindex + 1
                this%afish(pindex) = 1
            EndDo
            Do i = 1,this%nfishomniv
                pindex = pindex + 1
                this%afish(pindex) = 2
            EndDo
            Do i = 1,this%nfishplank
                pindex = pindex + 1
                this%afish(pindex) = 3
            EndDo
            ALLOCATE (this%sDFiAd0(this%nfish))
            ALLOCATE (this%sDFiJv0(this%nfish))
            ALLOCATE (this%kMigrFiAd(this%nfish))
            ALLOCATE (this%cDFiAdIn(this%nfish))
            ALLOCATE (this%kMigrFiJv(this%nfish))
            ALLOCATE (this%cDFiJvIn(this%nfish))
            ALLOCATE (this%cDayReprFish(this%nfish))
            ALLOCATE (this%fReprFish(this%nfish))
            ALLOCATE (this%fAgeFish(this%nfish))
            ALLOCATE (this%cRelVegFish(this%nfish))
            ALLOCATE (this%cSigTmFish(this%nfish))
            ALLOCATE (this%cTmOptFish(this%nfish))
            ALLOCATE (this%hDFiAd(this%nfish))
            ALLOCATE (this%hDFiJv(this%nfish))
            ALLOCATE (this%kDAssFiAd(this%nfish))
            ALLOCATE (this%kDAssFiJv(this%nfish))
            ALLOCATE (this%fDAssFiAd(this%nfish))
            ALLOCATE (this%fDAssFiJv(this%nfish))
            ALLOCATE (this%kHarvFiAd(this%nfish))
            ALLOCATE (this%kHarvFiJv(this%nfish))
            ALLOCATE (this%kDRespFiAd(this%nfish))
            ALLOCATE (this%kDRespFiJv(this%nfish))
            ALLOCATE (this%kPelFiAd(this%nfish))
            ALLOCATE (this%kPelFiJv(this%nfish))
            ALLOCATE (this%kMortFiAd(this%nfish))
            ALLOCATE (this%kMortFiJv(this%nfish))
            ALLOCATE (this%cDCarrFish(this%nfish))
            ALLOCATE (this%cCDFishRef(this%nfish))
            ALLOCATE (this%cNDFishRef(this%nfish))
            ALLOCATE (this%cPDFishRef(this%nfish))            
        EndIf
        
        !Macrophytes    
        If (this%InclMacr==0) Then
            this%nmac = 0
        Else
            this%nmac = this%nmacElod + this%nmacCharo + this%nmacCera + this%nmacLemna + this%nmacNymp + this%nmacHelo + this%nmacGeral
            ALLOCATE (this%amac(this%nmac))
            pindex = 0
            Do i = 1,this%nmacElod
                pindex = pindex + 1
                this%amac(pindex) = 1
            EndDo
            Do i = 1,this%nmacCharo
                pindex = pindex + 1
                this%amac(pindex) = 2
            EndDo
            Do i = 1,this%nmacCera
                pindex = pindex + 1
                this%amac(pindex) = 3
            EndDo
            Do i = 1,this%nmacLemna
                pindex = pindex + 1
                this%amac(pindex) = 4
            EndDo
            Do i = 1,this%nmacNymp
                pindex = pindex + 1
                this%amac(pindex) = 5
            EndDo
            Do i = 1,this%nmacHelo
                pindex = pindex + 1
                this%amac(pindex) = 6
            EndDo
            Do i = 1,this%nmacGeral
                pindex = pindex + 1
                this%amac(pindex) = 7
            EndDo
            ALLOCATE (this%sDMac0(this%nmac))
            ALLOCATE (this%cCDMac0(this%nmac))
            ALLOCATE (this%cNDMac0(this%nmac))
            ALLOCATE (this%cPDMac0(this%nmac))
            ALLOCATE (this%cQ10ProdMac(this%nmac))
            ALLOCATE (this%cMuMaxMac(this%nmac))
            ALLOCATE (this%hLRefMac(this%nmac))
            ALLOCATE (this%cTmInitVeg(this%nmac))
            ALLOCATE (this%fRootVegWin(this%nmac))
            ALLOCATE (this%fRootVegSum(this%nmac))
            ALLOCATE (this%fEmergVeg(this%nmac))
            ALLOCATE (this%fFloatVeg(this%nmac))
            ALLOCATE (this%cDLayerVeg(this%nmac))
            ALLOCATE (this%cDCarrMac(this%nmac))
            ALLOCATE (this%cCovSpVeg(this%nmac))
            ALLOCATE (this%cCDMacMax(this%nmac))
            ALLOCATE (this%cNDMacMax(this%nmac))
            ALLOCATE (this%cPDMacMax(this%nmac))
            ALLOCATE (this%cCDMacMin(this%nmac))
            ALLOCATE (this%cNDMacMin(this%nmac))
            ALLOCATE (this%cPDMacMin(this%nmac))
            ALLOCATE (this%kDRespMac(this%nmac))
            ALLOCATE (this%cQ10RespMac(this%nmac))
            ALLOCATE (this%kMortMac(this%nmac))
            ALLOCATE (this%cVCUptMaxMac(this%nmac))
            ALLOCATE (this%cVNUptMaxMac(this%nmac))
            ALLOCATE (this%cVPUptMaxMac(this%nmac))
            ALLOCATE (this%cAffCUptMac(this%nmac))
            ALLOCATE (this%cAffNUptMac(this%nmac))
            ALLOCATE (this%cAffPUptMac(this%nmac))
            ALLOCATE (this%fSedUptVegMax(this%nmac))
            ALLOCATE (this%fSedUptVegCoef(this%nmac))
            ALLOCATE (this%fSedUptVegExp(this%nmac))
            ALLOCATE (this%cPrefMacBird(this%nmac))
            ALLOCATE (this%hDMacBird(this%nmac))
            ALLOCATE (this%cDGrazPerBird(this%nmac))
            ALLOCATE (this%fDAssBird(this%nmac))            
        EndIf
        
        
        
        
        Do i = 1, wqConfiguration%numberOfParameters
        
            text = trim(wqParam(i)%name)
        
                !Heat Budget Parameters
            If (trim(text) == 'cExtWat') Then
                this%cExtWat = wqParam(i)%value
            ElseIf (trim(text) == 'tau_param') Then 
               this%tau_param = wqParam(i)%value
            ElseIf (trim(text) == 'a_param') Then 
               this%a_param = wqParam(i)%value
            ElseIf (trim(text) == 'e_param') Then 
               this%e_param = wqParam(i)%value
            ElseIf (trim(text) == 'c1_param') Then 
               this%c1_param = wqParam(i)%value
            ElseIf (trim(text) == 'Cd') Then 
               this%cd_param = wqParam(i)%value
            ElseIf (trim(text) == 'Cd') Then 
               this%cd_param = wqParam(i)%value
               
            !General Parameters
            ElseIf (trim(text) == 'cCPerDW') Then 
               this%cCPerDW = wqParam(i)%value
            ElseIf (trim(text) == 'O2PerNO3') Then 
               this%O2PerNO3 = wqParam(i)%value
            ElseIf (trim(text) == 'molO2molC') Then 
               this%molO2molC = wqParam(i)%value
            ElseIf (trim(text) == 'molO2molN') Then 
               this%molO2molN = wqParam(i)%value
            ElseIf (trim(text) == 'molNmolC') Then 
               this%molNmolC = wqParam(i)%value
            ElseIf (trim(text) == 'O2PerNH4') Then 
               this%O2PerNH4 = wqParam(i)%value
            ElseIf (trim(text) == 'NO3PerC') Then 
               this%NO3PerC = wqParam(i)%value
            ElseIf (trim(text) == 'pHmax') Then 
               this%pHmax = wqParam(i)%value
            ElseIf (trim(text) == 'pHmin') Then 
               this%pHmin = wqParam(i)%value
            ElseIf (trim(text) == 'fPAR') Then 
               this%fPAR = wqParam(i)%value
            !Inorganic Matter Parameters
            ElseIf (trim(text) == 'cExtSpIM') Then
                this%cExtSpIM = wqParam(i)%value
            ElseIf (trim(text) == 'cVSetIM') Then
                this%cVSetIM = wqParam(i)%value
            ElseIf (trim(text) == 'ppIM') Then
                this%ppIM = wqParam(i)%value
            ElseIf (trim(text) == 'diIM') Then
                this%diIM = wqParam(i)%value
          
            !this%InclMatOrg == 1
            !this%InclMatOrgSplit == 1
            !this%InclPomDomSplit == 1

            !Organic Matter Parameters
                !Case 1
            ElseIf (trim(text) == 'cExtSpOM') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%cExtSpPom(1) = wqParam(i)%value
                    this%cExtSpDom(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cCDOMRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%cCDDomRef(1) = wqParam(i)%value
                    this%cCDPomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cNDOMRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%cNDDomRef(1) = wqParam(i)%value
                    this%cNDPomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cPDOMRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%cPDDomRef(1) = wqParam(i)%value
                    this%cPDPomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cOMGround') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%cDomGround(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cVSetOM') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%cVSetPom = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'ppOM') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%ppOM = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'diOM') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%diOM = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaMinAerW1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%cThetaMinAerW = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaMinAerS1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%cThetaMinAerS = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaMinAnaerW1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%cThetaMinAnaerW = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaMinAnaerS1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%cThetaMinAnaerS = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kMinAerOMW1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%kMinAerDomW(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kMinAerOMS1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%kMinAerDomS(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kMinAnaerOMW1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%kMinAnaerDomW(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kMinAnaerOMS1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 0.and.this%InclPomDomSplit == 0) Then
                    this%kMinAnaerDomS(1) = wqParam(i)%value
                EndIf
                !Case 2
            ElseIf (trim(text) == 'cExtSpPom') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cExtSpPom(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cExtSpDom') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cExtSpDom(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cCDPomRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cCDPomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cNDPomRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cNDPomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cPDPomRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cPDPomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cCDDomRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cCDDomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cNDDomRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cNDDomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cPDDomRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cPDDomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cDomGround') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cDomGround(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cVSetPom1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cVSetPom = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'ppOM2') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%ppOM = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'diOM2') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%diOM = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaHidW1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cThetaHidW = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaHidS1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cThetaHidS = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'hBacHid1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%hBacHid = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kHidPomW1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%kHidPomW(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kHidPomS1') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%kHidPomS(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaMinAerW2') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cThetaMinAerW = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaMinAerS2') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cThetaMinAerS = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaMinAnaerW2') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cThetaMinAnaerW = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaMinAnaerS2') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%cThetaMinAnaerS = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kMinAerDomW') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%kMinAerDomW(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kMinAerDomS') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%kMinAerDomS(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kMinAnaerDomW') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%kMinAnaerDomW(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kMinAnaerDomS') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 0) Then
                    this%kMinAnaerDomS(1) = wqParam(i)%value
                EndIf
                !Case 3
            ElseIf (trim(text) == 'fDLisDomPhyt') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDLisDomPhyt = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDMortDomMac') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDMortDomMac = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDMessDomZoo') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDMessDomZoo = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDExcDomZoo') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDExcDomZoo = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDEgesDomZoo') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDEgesDomZoo = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDMortDomZoo') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDMortDomZoo = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDMessDomBent') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDMessDomBen = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDEgesDomBent') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDEgesDomBen = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDMortDomBent') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDMortDomBen = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDMessDomFish') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDMessDomFish = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDExcDomFish') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDExcDomFish = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDEgesDomFish') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDEgesDomFish = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDMortDomFish') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDMortDomFish = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fDEgesDomBird') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%fDAssDomBird = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cExtSpPomL') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cExtSpPom(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cExtSpPomR') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cExtSpPom(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cExtSpDomL') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cExtSpDom(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cExtSpDomR') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cExtSpDom(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cCDPomLRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cCDPomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cCDPomRRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cCDPomRef(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cNDPomLRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cNDPomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cNDPomRRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cNDPomRef(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cPDPomLRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cPDPomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cPDPomRRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cPDPomRef(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cCDDomLRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cCDDomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cCDDomRRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cCDDomRef(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cNDDomLRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cNDDomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cNDDomRRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cNDDomRef(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cPDDomLRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cPDDomRef(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cPDDomRRef') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cPDDomRef(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cDomLGround') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cDomGround(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cDomRGround') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cDomGround(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cVSetPom2') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cVSetPom = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'ppOM3') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%ppOM = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'diOM3') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%diOM = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaHidW2') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cThetaHidW = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaHidS2') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cThetaHidS = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'hBacHid2') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%hBacHid = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kHidPomLW') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%kHidPomW(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kHidPomRW') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%kHidPomW(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kHidPomLS') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%kHidPomS(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kHidPomRS') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%kHidPomS(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaMinAerW3') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cThetaMinAerW = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaMinAerS3') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cThetaMinAerS = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaMinAnaerW3') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cThetaMinAnaerW = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'cThetaMinAnaerS3') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%cThetaMinAnaerS = wqParam(i)%value                
                EndIf
            ElseIf (trim(text) == 'kMinAerDomLW') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%kMinAerDomW(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kMinAerDomRW') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%kMinAerDomW(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kMinAerDomLS') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%kMinAerDomS(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kMinAerDomRS') Then
                If (this%InclMatOrg == 1.and.this%InclMatOrgSplit == 1.and.this%InclPomDomSplit == 1) Then
                    this%kMinAerDomS(2) = wqParam(i)%value
                EndIf
                !Sediment Parameters
            ElseIf (trim(text) == 'nsed') Then 
               this%nsed = wqParam(i)%value
            ElseIf (trim(text) == 'cDepthS') Then 
               this%cDepthS = wqParam(i)%value
            ElseIf (trim(text) == 'cRhoIM') Then 
               this%cRhoIM = wqParam(i)%value
            ElseIf (trim(text) == 'cRhoOM') Then 
               this%cRhoOM = wqParam(i)%value
            ElseIf (trim(text) == 'fDepthDifS') Then 
               this%fDepthDifS = wqParam(i)%value
            ElseIf (trim(text) == 'cThetaDif') Then 
               this%cThetaDif = wqParam(i)%value
            ElseIf (trim(text) == 'kPDifDomR') Then 
                If (this%ndom>=1) Then
                    this%kPDifDom(1) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kPDifDomL') Then 
                If (this%ndom==2) Then
                    this%kPDifDom(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kCDifDic') Then 
               this%kCDifDic = wqParam(i)%value
            ElseIf (trim(text) == 'kPDifPO4') Then 
               this%kPDifPO4 = wqParam(i)%value
            ElseIf (trim(text) == 'kNDifNO3') Then 
               this%kNDifNO3 = wqParam(i)%value
            ElseIf (trim(text) == 'kNDifNH4') Then 
               this%kNDifNH4 = wqParam(i)%value
            ElseIf (trim(text) == 'kO2Dif') Then 
               this%kO2Dif = wqParam(i)%value
            ElseIf (trim(text) == 'cTurbDifNut') Then 
               this%cTurbDifNut = wqParam(i)%value
            ElseIf (trim(text) == 'cturbDifO2') Then 
               this%cturbDifO2 = wqParam(i)%value
            ElseIf (trim(text) == 'cDErosTot') Then 
               this%cDErosTot = wqParam(i)%value
            ElseIf (trim(text) == 'fDOrgSoil') Then 
               this%fDOrgSoil = wqParam(i)%value
            ElseIf (trim(text) == 'fSedErosIM') Then 
               this%fSedErosIM = wqParam(i)%value
            ElseIf (trim(text) == 'cSuspRef') Then 
               this%cSuspRef = wqParam(i)%value
            ElseIf (trim(text) == 'cSuspMin') Then 
               this%cSuspMin = wqParam(i)%value
            ElseIf (trim(text) == 'cSuspMax') Then 
               this%cSuspMax = wqParam(i)%value
            ElseIf (trim(text) == 'cSuspSlope') Then 
               this%cSuspSlope = wqParam(i)%value
            ElseIf (trim(text) == 'hDepthSusp') Then 
               this%hDepthSusp = wqParam(i)%value
            ElseIf (trim(text) == 'cFetchRef') Then 
               this%cFetchRef = wqParam(i)%value
            ElseIf (trim(text) == 'kTurbFish') Then 
               this%kTurbFish = wqParam(i)%value
            ElseIf (trim(text) == 'kVegResus') Then 
               this%kVegResus = wqParam(i)%value
            ElseIf (trim(text) == 'fDTotS0') Then 
               this%fDTotS0 = wqParam(i)%value
            ElseIf (trim(text) == 'fDOrgS0') Then 
               this%fDOrgS0 = wqParam(i)%value
            ElseIf (trim(text) == 'fDLabilS0') Then 
               this%fDLabilS0 = wqParam(i)%value
            ElseIf (trim(text) == 'fPAdsS0') Then 
               this%fPAdsS0 = wqParam(i)%value
            ElseIf (trim(text) == 'fPInorgS0') Then 
               this%fPInorgS0 = wqParam(i)%value
            ElseIf (trim(text) == 'fLutum') Then 
               this%fLutum = wqParam(i)%value
            ElseIf (trim(text) == 'fLutumRef') Then 
               this%fLutumRef = wqParam(i)%value
               
            !Nutrient&Gases Parameters
            ElseIf (trim(text) == 'cThetaNitr') Then 
               this%cThetaNitr = wqParam(i)%value
            ElseIf (trim(text) == 'kNitrS') Then 
               this%kNitrS = wqParam(i)%value
            ElseIf (trim(text) == 'kNitrW') Then 
               this%kNitrW = wqParam(i)%value
            ElseIf (trim(text) == 'hO2Nitr') Then 
               this%hO2Nitr = wqParam(i)%value
            ElseIf (trim(text) == 'hNO3Denit') Then 
               this%hNO3Denit = wqParam(i)%value
            ElseIf (trim(text) == 'cRelPAdsD') Then 
               this%cRelPAdsD = wqParam(i)%value
            ElseIf (trim(text) == 'cRelPAdsFe') Then 
               this%cRelPAdsFe = wqParam(i)%value
            ElseIf (trim(text) == 'cRelPAdsAl') Then 
               this%cRelPAdsAl = wqParam(i)%value
            ElseIf (trim(text) == 'fFeDIM') Then 
               this%fFeDIM = wqParam(i)%value
            ElseIf (trim(text) == 'fAlDIM') Then 
               this%fAlDIM = wqParam(i)%value
            ElseIf (trim(text) == 'fRedMax') Then 
               this%fRedMax = wqParam(i)%value
            ElseIf (trim(text) == 'cKPAdsOx') Then 
               this%cKPAdsOx = wqParam(i)%value
            ElseIf (trim(text) == 'kPSorp') Then 
               this%kPSorp = wqParam(i)%value
            ElseIf (trim(text) == 'kPChemPO4') Then 
               this%kPChemPO4 = wqParam(i)%value
            ElseIf (trim(text) == 'cPO4Max') Then 
               this%cPO4Max = wqParam(i)%value
            ElseIf (trim(text) == 'cThetaReaer') Then 
               this%cThetaReaer = wqParam(i)%value
        
            !Bacterioplankton Parameters
            ElseIf (trim(text) == 'cThetaBac') Then 
                If (this%InclBac == 1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cThetaBac(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'hO2Bac') Then
                If (this%InclBac == 1) Then
                    this%hO2Bac = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'fbacAn') Then 
                If (this%InclBac == 1) Then
                    this%fbacAn = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'kbacAn') Then 
                If (this%InclBac == 1) Then
                    this%kbacAn = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'MuBac') Then 
                If (this%InclBac == 1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%MuBac(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'hAssBac') Then 
                If (this%InclBac == 1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%hAssBac(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kResbBac') Then 
                If (this%InclBac == 1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kResbBac(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cCDBacRef') Then 
                If (this%InclBac == 1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cCDBacRef(pindex) = wqvalues(k)
                        EndDo
                    EndDo 
                EndIf
            ElseIf (trim(text) == 'cNDBacRef') Then 
                If (this%InclBac == 1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cNDBacRef(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cPDBacRef') Then 
                If (this%InclBac == 1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cPDBacRef(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kMortBac') Then 
                If (this%InclBac == 1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kMortBac(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'hSedBac') Then
                If (this%InclBac == 1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%hSedBac(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
                
            !Phytoplankton Parameters
            ElseIf (trim(text) == 'cExtSpPhyt') Then 
                If (this%InclPhyt==1) Then                
                call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                pindex = 0
                Do j = 1, wqParam(i)%numberOfGroups
                    call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                    Do k = 1, wqgroups(j)%numberOfValues
                        pindex = pindex + 1
                        this%cExtSpPhyt(pindex) = wqvalues(k)
                    EndDo
                EndDo
                EndIf
            ElseIf (trim(text) == 'hLRefPhyt') Then
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%hLRefPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo  
                EndIf
            ElseIf (trim(text) == 'cLOptRefPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cLOptRefPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo
                EndIf
            ElseIf (trim(text) == 'cMuMaxPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cMuMaxPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo
                EndIf
            ElseIf (trim(text) == 'cChDPhytMin') Then
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cChDPhytMin(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cChDPhytMax') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cChDPhytMax(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cSigTmPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cSigTmPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cTmOptPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cTmOptPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo  
                EndIf
            ElseIf (trim(text) == 'cVCUptMaxPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cVCUptMaxPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cVNUptMaxPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cVNUptMaxPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo  
                EndIf
            ElseIf (trim(text) == 'cVPUptMaxPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cVPUptMaxPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cVSiUptMaxPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cVSiUptMaxPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cCDPhytMax') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cCDPhytMax(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cCDPhytMin') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cCDPhytMin(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cNDPhytMax') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cNDPhytMax(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cNDPhytMin') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cNDPhytMin(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cPDPhytMax') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cPDPhytMax(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cPDPhytMin') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cPDPhytMin(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cSiDPhytMax') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cSiDPhytMax(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cSiDPhytMin') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cSiDPhytMin(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cAffCUptPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cAffCUptPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cAffNUptPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cAffNUptPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cAffPUptPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cAffPUptPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cAffSiUptPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cAffSiUptPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'hSiAssDiat') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%hSiAssDiat(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kDRespbPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kDRespbPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'alfap') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%alfap(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cLisPhytW') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cLisPhytW(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cVSetPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cVSetPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'ppPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%ppPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'diPhyt') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%diPhyt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cPACoefMin') Then 
                this%cPACoefMin = wqParam(i)%value
            ElseIf (trim(text) == 'cPACoefMax') Then 
                this%cPACoefMax = wqParam(i)%value
            ElseIf (trim(text) == 'hPACoef') Then 
                this%hPACoef = wqParam(i)%value
                
            !Zooplankton Parameters
            ElseIf (trim(text) == 'fDAssZoo') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fDAssZoo(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cSigTmZoo') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cSigTmZoo(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cTmOptZoo') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cTmOptZoo(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cDCarrZoo') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cDCarrZoo(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'hFilt') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%hFilt(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cFiltMax') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cFiltMax(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kPelZoo') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kPelZoo(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kDRespZoo') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kDRespZoo(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kMortZoo') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kMortZoo(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cCDZooRef') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cCDZooRef(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cNDZooRef') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cNDZooRef(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cPDZooRef') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cPDZooRef(pindex) = wqvalues(k)
                        EndDo
                    EndDo  
                EndIf
                
            !Zoobenthos Parameters
            ElseIf (trim(text) == 'fDAssBent') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fDAssBent(pindex) = wqvalues(k)
                        EndDo
                    EndDo 
                EndIf
            ElseIf (trim(text) == 'cSigTmBent') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cSigTmBent(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cTmOptBen') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cTmOptBent(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cDCarrBent') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cDCarrBent(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kDAssBent') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kDAssBent(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'hDFoodBent') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%hDFoodBent(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kPelBent') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kPelBent(pindex) = wqvalues(k)
                        EndDo
                    EndDo  
                EndIf
            ElseIf (trim(text) == 'kDRespBent') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kDRespBent(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kMortBent') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kMortBent(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cCDBentRef') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cCDBentRef(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cNDBentRef') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cNDBentRef(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cPDBentRef') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cPDBentRef(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
                
            !Fish Parameters
            ElseIf (trim(text) == 'kMigrFiAd') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kMigrFiAd(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cDFiAdIn') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cDFiAdIn(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kMigrFiJv') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kMigrFiJv(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cDFiJvIn') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cDFiJvIn(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cDayReprFish') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cDayReprFish(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'fReprFish') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fReprFish(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'fAgeFish') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fAgeFish(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cRelVegFish') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cRelVegFish(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cSigTmFish') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cSigTmFish(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cTmOptFish') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cTmOptFish(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'hDFiAd') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%hDFiAd(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'hDFiJv') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%hDFiJv(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kDAssFiAd') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kDAssFiAd(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kDAssFiJv') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kDAssFiJv(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'fDAssFiAd') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fDAssFiAd(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'fDAssFiJv') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fDAssFiJv(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cDCarrFish') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cDCarrFish(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kHarvFiAd') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kHarvFiAd(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kHarvFiJv') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kHarvFiJv(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kDRespFiAd') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kDRespFiAd(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kDRespFiJv') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kDRespFiJv(pindex) = wqvalues(k)
                        EndDo
                    EndDo  
                EndIf
            ElseIf (trim(text) == 'kPelFiAd') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kPelFiAd(pindex) = wqvalues(k)
                        EndDo
                    EndDo  
                EndIf
            ElseIf (trim(text) == 'kPelFiJv') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kPelFiJv(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kMortFiAd') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kMortFiAd(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'kMortFiJv') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kMortFiJv(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cCDFishRef') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cCDFishRef(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cNDFishRef') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cNDFishRef(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cPDFishRef') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cPDFishRef(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
                
            !Macrophytes Parameters
            ElseIf (trim(text) == 'cQ10ProdMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cQ10ProdMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'hLRefMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%hLRefMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo    
                EndIf
            ElseIf (trim(text) == 'cMuMaxMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cMuMaxMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo            
                EndIf
            ElseIf (trim(text) == 'cTmInitVeg') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cTmInitVeg(pindex) = wqvalues(k)
                        EndDo
                    EndDo    
                EndIf
            ElseIf (trim(text) == 'fRootVegWin') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fRootVegWin(pindex) = wqvalues(k)
                        EndDo
                    EndDo      
                EndIf
            ElseIf (trim(text) == 'fRootVegSum') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fRootVegSum(pindex) = wqvalues(k)
                        EndDo
                    EndDo       
                EndIf
            ElseIf (trim(text) == 'fEmergVeg') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fEmergVeg(pindex) = wqvalues(k)
                        EndDo
                    EndDo  
                EndIf
            ElseIf (trim(text) == 'fFloatVeg') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fFloatVeg(pindex) = wqvalues(k)
                        EndDo
                    EndDo         
                EndIf
            ElseIf (trim(text) == 'cDLayerVeg') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cDLayerVeg(pindex) = wqvalues(k)
                        EndDo
                    EndDo        
                EndIf
            ElseIf (trim(text) == 'cDCarrMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cDCarrMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo       
                EndIf
            ElseIf (trim(text) == 'cCovSpVeg') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cCovSpVeg(pindex) = wqvalues(k)
                        EndDo
                    EndDo        
                EndIf
            ElseIf (trim(text) == 'cVCUptMaxMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cVCUptMaxMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo       
                EndIf
            ElseIf (trim(text) == 'cVNUptMaxMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cVNUptMaxMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo    
                EndIf
            ElseIf (trim(text) == 'cVPUptMaxMac') Then
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cVPUptMaxMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo  
                EndIf
            ElseIf (trim(text) == 'cCDMacMax') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cCDMacMax(pindex) = wqvalues(k)
                        EndDo
                    EndDo    
                EndIf
            ElseIf (trim(text) == 'cCDMacMin') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cCDMacMin(pindex) = wqvalues(k)
                        EndDo
                    EndDo 
                EndIf
            ElseIf (trim(text) == 'cNDMacMin') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cNDMacMin(pindex) = wqvalues(k)
                        EndDo
                    EndDo       
                EndIf
            ElseIf (trim(text) == 'cNDMacMax') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cNDMacMax(pindex) = wqvalues(k)
                        EndDo
                    EndDo       
                EndIf
            ElseIf (trim(text) == 'cPDMacMin') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cPDMacMin(pindex) = wqvalues(k)
                        EndDo
                    EndDo       
                EndIf
            ElseIf (trim(text) == 'cPDMacMax') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cPDMacMax(pindex) = wqvalues(k)
                        EndDo
                    EndDo        
                EndIf
            ElseIf (trim(text) == 'cAffCUptMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cAffCUptMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cAffNUptMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cAffNUptMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'cAffPUptMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cAffPUptMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo    
                EndIf
            ElseIf (trim(text) == 'fSedUptVegMax') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fSedUptVegMax(pindex) = wqvalues(k)
                        EndDo
                    EndDo    
                EndIf
            ElseIf (trim(text) == 'fSedUptVegCoef') Then 
                    If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fSedUptVegCoef(pindex) = wqvalues(k)
                        EndDo
                    EndDo      
                EndIf
            ElseIf (trim(text) == 'fSedUptVegExp') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fSedUptVegExp(pindex) = wqvalues(k)
                        EndDo
                    EndDo  
                EndIf
            ElseIf (trim(text) == 'cQ10RespMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cQ10RespMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo      
                EndIf
            ElseIf (trim(text) == 'kDRespMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kDRespMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo     
                EndIf
            ElseIf (trim(text) == 'kMortMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%kMortMac(pindex) = wqvalues(k)
                        EndDo
                    EndDo 
                EndIf
            ElseIf (trim(text) == 'hDMacBird') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%hDMacBird(pindex) = wqvalues(k)
                        EndDo
                    EndDo 
                EndIf
            ElseIf (trim(text) == 'cDGrazPerBird') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%cDGrazPerBird(pindex) = wqvalues(k)
                        EndDo
                    EndDo    
                EndIf
            ElseIf (trim(text) == 'fDAssBird') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%fDAssBird(pindex) = wqvalues(k)
                        EndDo
                    EndDo  
                EndIf
                
            !Initial Conditions
            ElseIf (trim(text) == 'sDTempW') Then 
               this%sDTempW0 = wqParam(i)%value
            ElseIf (trim(text) == 'sDSal') Then 
               this%sDSal0 = wqParam(i)%value
            ElseIf (trim(text) == 'sO2W') Then 
               this%sO2W0 = wqParam(i)%value
            ElseIf (trim(text) == 'spHW') Then 
               this%spHW0 = wqParam(i)%value
            ElseIf (trim(text) == 'spHS') Then 
               this%spHS0 = wqParam(i)%value
            ElseIf (trim(text) == 'sAlkW') Then 
               this%sAlkW0 = wqParam(i)%value
            ElseIf (trim(text) == 'sAlkS') Then 
               this%sAlkS0 = wqParam(i)%value
            ElseIf (trim(text) == 'sDicW') Then 
               this%sDicW0 = wqParam(i)%value
            ElseIf (trim(text) == 'sDicS') Then 
               this%sDicS0 = wqParam(i)%value
            ElseIf (trim(text) == 'sCH4W') Then 
               this%sCH4W0 = wqParam(i)%value
            ElseIf (trim(text) == 'sCH4S') Then 
               this%sCH4S0 = wqParam(i)%value
            ElseIf (trim(text) == 'sPO4W') Then 
               this%sPO4W0 = wqParam(i)%value
            ElseIf (trim(text) == 'sPAIMW') Then 
               this%sPAIMW0 = wqParam(i)%value
            ElseIf (trim(text) == 'sNO3W') Then 
               this%sNO3W0 = wqParam(i)%value
            ElseIf (trim(text) == 'sNO3S') Then 
               this%sNO3S0 = wqParam(i)%value
            ElseIf (trim(text) == 'sNH4W') Then 
               this%sNH4W0 = wqParam(i)%value
            ElseIf (trim(text) == 'sNH4S') Then 
               this%sNH4S0 = wqParam(i)%value
            ElseIf (trim(text) == 'sSiO2W') Then 
               this%sSiO2W0 = wqParam(i)%value
            ElseIf (trim(text) == 'sDIMW') Then 
               this%sDIMW0 = wqParam(i)%value
            ElseIf (trim(text) == 'sOMW') Then 
               this%sDPomW0(1) = wqParam(i)%value
               this%sDDomW0(1) = wqParam(i)%value
            ElseIf (trim(text) == 'sOMS') Then 
               this%sDPomS0(1) = wqParam(i)%value
               this%sDDomS0(1) = wqParam(i)%value
            ElseIf (trim(text) == 'sDPomW') Then 
               this%sDPomW0(1) = wqParam(i)%value
            ElseIf (trim(text) == 'sDPomS') Then 
               this%sDPomS0(1) = wqParam(i)%value
            ElseIf (trim(text) == 'sDDomW') Then 
               this%sDDomW0(1) = wqParam(i)%value
            ElseIf (trim(text) == 'sDDomS') Then 
               this%sDDomS0(1) = wqParam(i)%value
            ElseIf (trim(text) == 'sDPomLW') Then 
               this%sDPomW0(1) = wqParam(i)%value
            ElseIf (trim(text) == 'sDPomRW') Then 
                If (this%InclPomDomSplit==1) Then
                    this%sDPomW0(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'sDPomLS') Then 
               this%sDPomS0(1) = wqParam(i)%value
            ElseIf (trim(text) == 'sDPomRS') Then 
                If (this%InclPomDomSplit==1) Then
                    this%sDPomS0(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'sDDomLW') Then 
               this%sDDomW0(1) = wqParam(i)%value
            ElseIf (trim(text) == 'sDDomRW') Then 
                If (this%InclPomDomSplit==1) Then
                    this%sDDomW0(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'sDDomLS') Then 
               this%sDDomS0(1) = wqParam(i)%value
            ElseIf (trim(text) == 'sDDomRS') Then 
                If (this%InclPomDomSplit==1) Then
                    this%sDDomS0(2) = wqParam(i)%value
                EndIf
            ElseIf (trim(text) == 'sDBacW') Then 
                If (this%InclBac==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%sDBacW0(pindex) = wqvalues(k)
                        EndDo
                    EndDo    
                EndIf
            ElseIf (trim(text) == 'sDBacS') Then 
                If (this%InclBac==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%sDBacS0(pindex) = wqvalues(k)
                        EndDo
                    EndDo     
                EndIf
            ElseIf (trim(text) == 'sDPhytW') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%sDPhytW0(pindex) = wqvalues(k)
                        EndDo
                    EndDo 
                EndIf
            ElseIf (trim(text) == 'sDPhytS') Then 
                If (this%InclPhyt==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%sDPhytS0(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'sDMac') Then 
                If (this%InclMacr==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%sDMac0(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'sDZoo') Then 
                If (this%InclZoo==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%sDZoo0(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'sDBent') Then 
                If (this%InclBent==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%sDBent0(pindex) = wqvalues(k)
                        EndDo
                    EndDo 
                EndIf
            ElseIf (trim(text) == 'sDFiAd') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%sDFiAd0(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
            ElseIf (trim(text) == 'sDFiJv') Then 
                If (this%InclFish==1) Then
                    call c_f_pointer(wqParam(i)%groups, wqgroups, [wqParam(i)%numberOfGroups])
                    pindex = 0
                    Do j = 1, wqParam(i)%numberOfGroups
                        call c_f_pointer(wqgroups(j)%values, wqvalues, [wqgroups(j)%numberOfValues])
                        Do k = 1, wqgroups(j)%numberOfValues
                            pindex = pindex + 1
                            this%sDFiJv0(pindex) = wqvalues(k)
                        EndDo
                    EndDo   
                EndIf
               
            EndIf
        EndDo
        
        !Reading foodweb matrix
        call c_f_pointer(wqConfiguration%foodMatrix, wqFoodMatrix, [wqConfiguration%foodMatrixSize])
        
        ALLOCATE (this%cPrefZooZoo(this%nzoo,this%nzoo))
        ALLOCATE (this%cPrefZooPhyt(this%nzoo,this%nphy))
        ALLOCATE (this%cPrefZooBac(this%nzoo,this%nbac))
        
        Do i = 1, wqConfiguration%foodMatrixSize
            ! MicroZooplankton Preferences
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='om') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='pom') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,1) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='pomLabile') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,1) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='pomRefractory') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,2) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='anaerobicBacteria') Then
                this%cPrefZooBac(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup+this%nbacaer) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='aerobicBacteria') Then
                this%cPrefZooBac(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='dinoflagellates') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='chlorophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='cryptophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='marineDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            ! MesoZooplankton Preferences
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='om') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='pom') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,1) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='pomLabile') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,1) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='pomRefractory') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,2) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='anaerobicBacteria') Then
                this%cPrefZooBac(wqFoodMatrix(i)%predatorGroup+this%nzoomicro,wqFoodMatrix(i)%preyGroup+this%nbacaer) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='aerobicBacteria') Then
                this%cPrefZooBac(wqFoodMatrix(i)%predatorGroup+this%nzoomicro,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='dinoflagellates') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='chlorophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='cryptophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='marineDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            ! MacroZooplankton Preferences
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='om') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='pom') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,1) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='pomLabile') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,1) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='pomRefractory') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,2) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='anaerobicBacteria') Then
                this%cPrefZooBac(wqFoodMatrix(i)%predatorGroup+this%nzoomicro+this%nzoomeso,wqFoodMatrix(i)%preyGroup+this%nbacaer) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='aerobicBacteria') Then
                this%cPrefZooBac(wqFoodMatrix(i)%predatorGroup+this%nzoomicro+this%nzoomeso,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='dinoflagellates') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='chlorophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='cryptophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZooGroup'.and.trim(wqFoodMatrix(i)%prey)=='marineDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf    
            ! MicroZoobenthos Preferences
            If (trim(wqFoodMatrix(i)%predator)=='microZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='anaerobicBacteria') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,2) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='aerobicBacteria') Then
                this%cPrefZooPom(wqFoodMatrix(i)%predatorGroup,2) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='dinoflagellates') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='chlorophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='cryptophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='microZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='marineDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            ! MesoZoobenthos Preferences
            If (trim(wqFoodMatrix(i)%predator)=='mesoZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='anaerobicBacteria') Then
                this%cPrefZooBac(wqFoodMatrix(i)%predatorGroup,2) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='aerobicBacteria') Then
                this%cPrefZooBac(wqFoodMatrix(i)%predatorGroup,2) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='dinoflagellates') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='chlorophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='cryptophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='mesoZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='marineDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            ! MacroZoobenthos Preferences
            If (trim(wqFoodMatrix(i)%predator)=='macroZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='anaerobicBacteria') Then
                this%cPrefZooBac(wqFoodMatrix(i)%predatorGroup,2) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='aerobicBacteria') Then
                this%cPrefZooBac(wqFoodMatrix(i)%predatorGroup,2) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='dinoflagellates') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='chlorophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='cryptophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='macroZoobenthos'.and.trim(wqFoodMatrix(i)%prey)=='marineDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            ! Adult planktivorous fish Preferences
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='dinoflagellates') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='chlorophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='cryptophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='marineDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='microZooGroup') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='mesoZooGroup') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='macroZooGroup') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='microZoobenthos') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='mesoZoobenthos') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='macroZoobenthos') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='planktivorousJuvenile') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            ! Juvenile planktivorous fish Preferences
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='dinoflagellates') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='chlorophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='cryptophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='marineDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='microZooGroup') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='mesoZooGroup') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='macroZooGroup') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='microZoobenthos') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='mesoZoobenthos') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='planktivorousJuvenile'.and.trim(wqFoodMatrix(i)%prey)=='macroZoobenthos') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            ! Adult omnivorous fish Preferences
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='dinoflagellates') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='marineCyanobacteria') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='chlorophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='cryptophytes') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='freshwaterDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='marineDiatoms') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='microZooGroup') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='mesoZooGroup') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='macroZooGroup') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='microZoobenthos') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='mesoZoobenthos') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='macroZoobenthos') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='planktivorousJuvenile') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='planktivorousAdult') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            If (trim(wqFoodMatrix(i)%predator)=='omnivorousAdult'.and.trim(wqFoodMatrix(i)%prey)=='omnivorousJuvenile') Then
                this%cPrefZooPhyt(wqFoodMatrix(i)%predatorGroup,wqFoodMatrix(i)%preyGroup) = wqFoodMatrix(i)%value
            EndIf
            
            
        End Do        
        
        
        !this%cPrefZooZoo(:,1) = (/ 0.0,  0.0,   0.0 /)
        !this%cPrefZooZoo(:,2) = (/ 0.0,  0.0,   0.0 /)
        !this%cPrefZooZoo(:,3) = (/ 0.9,  0.75,   0.0 /)
        !)
        !this%cPrefZooPhyt(:,1) = (/ 0.4,  0.4,   0.4 /)
        !!this%cPrefZooPhyt(:,2) = (/ 0.8,  0.8,   0.8 /)
        !!this%cPrefZooPhyt(:,3) = (/ 0.4,  0.4,   0.4 /)
        !
        !this%cPrefZooBac(:,1) = (/ 0.8,  0.8,   0.8 /)
        !this%cPrefZooBac(:,2) = (/ 0.4,  0.4,   0.4 /)
        !this%cPrefZooBac(:,3) = (/ 0.4,  0.4,   0.4 /)
        !ALLOCATE (this%cPrefZooPom(this%nzoo,this%npom))
        !this%cPrefZooPom(:,1) = (/ 0.2,  0.2,   0.2 /)
        !!this%cPrefZooPom(:,2) = (/ 0.1,  0.1,   0.1 /)        
        
     !   this%iSal = 1
     !   !0. Numerical parameters 
     !   this%TranspFlag = 3
     !   
     !   this%cExtWat = 2.2 !<background extinction ok
     !   !1. Water temperature
     !   this%iTempW = 1 ! ok
     !   this%WtempRef = 20. ! ok
     !   this%a_param = 0.69 ! ok
     !   this%tau_param = 0.000000117! ok
     !   this%e_param = 0.8  ! ok
     !   this%c1_param = 0.47 ! ok
     !   this%cd_param = 4182 ! ok
     !   this%sDTempW0 = 12. ! ok
     !   
     !   !2.X General parameters
     !   this%cCPerDW = 0.4 !<C content of organic matter ok
     !   this%O2PerNO3 = 1.5 !<mol O2 formed per mol NO3- ammonified
	    !this%molO2molC = 2.6667 !< 32/12 [gO2/gC], ratio of mol.weights
	    !this%molO2molN = 2.2857; !<= 32/14 [gO2/gN], ratio of mol.weights
	    !this%molNmolC = 1.1667; !<= 14/12 [gN/gC], ratio of mol.weights 
     !   this%O2PerNH4 = 2.0 !<=mol O2 used per mol NH4+ nitrified [-]
     !   this%NO3PerC = 0.8 !<mol NO3 denitrified per mol C mineralised
     !   
     !   !2.X pH Modelling
     !   this%pHmax = 9.0
     !   this%pHmin = 5.0
     !   
     !   ! 2. Limnological Parameters
     !   ! 2.1 Organic Matter
     !   ! 2.1.1 Initial configuration
     !   this%InclMatOrgSplit = 0
     !   this%InclPomDomSplit = 0
     !   
     !   If (this%InclMatOrgSplit == 1) Then
     !   
     !       If (this%InclPomDomSplit == 1) Then
     !           this%npom = 2
     !           ALLOCATE (this%sDPomW0(this%npom))
     !           ALLOCATE (this%sDPomS0(this%npom))
     !           this%sDPomW0(:) = (/ 0.10, 0.00 /)     ! Defining Initial condition of POM functional groups in water
     !           this%sDPomS0(:) = (/ 0.50, 0.00 /)     ! Defining Initial condition of POM functional groups in sediment
     !       Else
     !           this%npom = 1
     !           ALLOCATE (this%sDPomW0(this%npom))
     !           ALLOCATE (this%sDPomS0(this%npom))
     !           this%sDPomW0(:) = (/ 0.10 /)     ! Defining Initial condition of POM functional groups in water
     !           this%sDPomS0(:) = (/ 0.50 /)     ! Defining Initial condition of POM functional groups in sediment
     !       EndIf
     !   
     !       If (this%InclPomDomSplit == 1) Then
     !           this%ndom = 2
     !           ALLOCATE (this%sDDomW0(this%ndom))
     !           ALLOCATE (this%sDDomS0(this%ndom))
     !           this%sDDomW0(:) = (/ 0.10, 0.00 /)     ! Defining Initial condition of DOM functional groups in water
     !           this%sDDomS0(:) = (/ 0.50, 0.00 /)     ! Defining Initial condition of DOM functional groups in sediment
     !       Else
     !           this%ndom = 1
     !           ALLOCATE (this%sDDomW0(this%ndom))
     !           ALLOCATE (this%sDDomS0(this%ndom))
     !           this%sDDomW0(:) = (/ 0.10 /)     ! Defining Initial condition of DOM functional groups in water
     !           this%sDDomS0(:) = (/ 0.50 /)     ! Defining Initial condition of DOM functional groups in sediment
     !       EndIf
     !   Else
     !       this%npom = 1
     !       this%ndom = 1
     !       ALLOCATE (this%sDPomW0(this%ndom))
     !       ALLOCATE (this%sDPomS0(this%ndom))
     !       ALLOCATE (this%sDDomW0(this%ndom))
     !       ALLOCATE (this%sDDomS0(this%ndom))
     !       this%sDPomW0(:) = (/ 0.10 /)     ! Defining Initial condition of DOM functional groups in water
     !       this%sDPomS0(:) = (/ 0.50 /)     ! Defining Initial condition of DOM functional groups in sediment
     !       this%sDDomW0(:) = this%sDPomW0(:)     ! Defining Initial condition of DOM functional groups in water
     !       this%sDDomS0(:) = this%sDPomS0(:)     ! Defining Initial condition of DOM functional groups in sediment
     !   EndIf
     !
     !   ! 2.1.2 Particulate Organic Matter (POM)
     !   ! 2.1.2.1 Parameters of Particulate Organic Matter (POM)
     !   ALLOCATE (this%cCDPomRef(this%npom))
     !   ALLOCATE (this%cNDPomRef(this%npom))
     !   ALLOCATE (this%cPDPomRef(this%npom))
     !   ALLOCATE (this%cExtSpPom(this%npom))
     !   this%cCDPomRef(:)     = (/ 0.86,  0.86 /)
     !   this%cNDPomRef(:)     = (/ 0.13,  0.13 /) 
     !   this%cPDPomRef(:)     = (/ 0.008, 0.008 /)
     !   this%cExtSpPom(:) = (/ 0.0, 0.0 /) !(/ 0.12, 0.12 /)
     !   !fDLisPomPhyt(:)  = (/ 1.00, 0.00 /) 
     !   !fDMessPomZoo(:)  = (/ 1.00, 0.00 /)
     !   !fDEgesPomZoo(:)  = (/ 1.00, 0.00 /)
     !   !fDMortPomZoo(:)  = (/ 1.00, 0.00 /)
     !   !fDMessPomBen(:)  = (/ 1.00, 0.00 /)
     !   !fDEgesPomBen(:)  = (/ 1.00, 0.00 /)
     !   !fDMortPomBen(:)  = (/ 1.00, 0.00 /)
     !   !fDMessPomFish(:) = (/ 1.00, 0.00 /)
     !   !fDEgesPomFish(:) = (/ 1.00, 0.00 /)
     !   !fDMortPomFish(:) = (/ 0.20, 0.80 /)
     !   !fDMortPomMac(:)  = (/ 0.80, 0.00 /)
     !   !fDMortPomMacW2S(:) = (/ 0.80, 0.00 /)  
     !   !fDAssPomBird(:)  = (/ 0.80, 0.00 /)
     !   ! 2.1.2.2 Infiltration
	    !ALLOCATE (this%cDomGround(this%ndom))  
     !   this%cDomGround(:) = (/ 0.0, 0.0 /)
     !   ! 2.1.2.4 Sedimentation
     !   this%cVSetPom  = 0.50 !SetMethodFlag = 0               
     !   this%ppOM      = 1030 !SetMethodFlag = 1
     !   this%diOM      = 0.000005 !SetMethodFlag = 1		
     !   ! 2.1.2.5 Hydrolises
	    !ALLOCATE (this%kHidPomW(this%npom))  
	    !ALLOCATE (this%kHidPomS(this%npom))  
     !   this%cThetaHidW    = 1.07		
     !   this%hBacHid       = 0.02		
     !   this%kHidPomW(:) = (/ 0.08,   0.0002 /)   
     !   this%cThetaHidS    = 1.07		
     !   this%kHidPomS(:) = (/ 0.0002, 0.0002 /) 
     !   ! 2.1.2.6 Mineralization (IF InclBac=0)
     !   this%cThetaMinAerW    = 1.07
     !   this%cThetaMinAnaerW    = 1.07	
     !   this%cThetaMinAerS    = 1.07	
     !   this%cThetaMinAnaerS    = 1.07
     !   ALLOCATE (this%kMinAerDomW(this%ndom))  
	    !ALLOCATE (this%kMinAnaerDomW(this%ndom))  
	    !ALLOCATE (this%kMinAerDomS(this%ndom))
     !   ALLOCATE (this%kMinAnaerDomS(this%ndom))
     !   this%kMinAerDomW(:)   = (/ 0.10,   0.002 /)
     !   this%kMinAnaerDomW(:)   = (/ 0.10,   0.002 /)
     !   this%kMinAerDomS(:)   = (/ 0.01,   0.002 /)
     !   this%kMinAnaerDomS(:)   = (/ 0.01,   0.002 /)
     !   ! 2.1.3 Dissolved Organic Matter (DOM)
     !   ! 2.1.3.1 Parameters of Dissolved Organic Matter (DOM)
     !   ALLOCATE (this%cCDDomRef(this%ndom))
     !   ALLOCATE (this%cNDDomRef(this%ndom))
     !   ALLOCATE (this%cPDDomRef(this%ndom))
     !   ALLOCATE (this%cSiDDomRef(this%ndom))
     !   ALLOCATE (this%cExtSpDom(this%ndom))
     !   this%cCDDomRef(:)     = (/ 0.86,   0.86 /) 
     !   this%cNDDomRef(:)     = (/ 0.13,   0.13 /)  
     !   this%cPDDomRef(:)     = (/ 0.008,  0.0025 /) 
     !   this%cSiDDomRef(:)    = (/ 0.01, 0.01 /) 
     !   this%cExtSpDom(:) = (/ 0.0, 0.0 /) !(/ 0.05, 0.05 /)
     !   this%fDMessDomZoo = 0.01
     !   this%fDExcDomZoo = 1.0
     !   this%fDEgesDomZoo  = 0.10
     !   this%fDMortDomZoo  = 0.10
     !   this%fDMessDomBen  = 0.10
     !   this%fDEgesDomBen = 0.10
     !   this%fDMortDomBen  = 0.10
     !   this%fDMessDomFish = 0.10
     !   this%fDExcDomFish  = 0.10
     !   this%fDEgesDomFish = 0.10
     !   this%fDMortDomFish = 0.10
     !   this%fDMortDomMac  = 0.20
     !   this%fDAssDomBird  = 0.20 
     !   this%fDLisDomPhyt  = 0.20 !Labile fractions of died_Algae (Autoctone Organic Matter)
     !
     !   ! 2.2 Inorganic Matter
     !   ! 2.2.1 Initial configuration
     !   this%InclMatInorg = 1
     !   this%InclMatInorgSplit = 1
     !   this%sDIMW0    = 5.0
     !   ! 2.2.2 Parameters
     !   this%cExtSpIM  = 0.0 !0.05
     !   !Sedimentation
     !   this%cVSetIM = 0.5 !SetMethodFlag = 0
     !   this%ppIM      = 2650 !SetMethodFlag = 1
     !   this%diIM      = 0.000003 !SetMethodFlag = 1	
     !
     !   ! 2.3 Inorganic Nutrients
     !   this%sO2W0    = 10.
     !   this%sPO4W0    = 0.005
     !   this%sPAIMW0   = 0.000
     !   this%sNO3W0    = 1.000
     !   this%sNO3S0    = 0.002      
     !   this%sNH4W0    = 0.001
     !   this%sNH4S0    = 0.020      
     !   this%sSiO2W0   = 3.0
     !   this%sDicW0    = 0.0035
     !   this%sDicS0    = 0.0035
     !   this%spHW0     = 8.0
     !   this%spHS0     = 8.0
     !   this%sAlkW0    = 100
     !   this%sAlkS0    = 100  
     !   this%sCH4W0    = 0.00
     !   this%sCH4S0    = 0.001      
     !   
     !
     !   ! 2.2 Bacterioplakton
     !   ! 2.2.1 Initial configuration
     !   this%InclBac = 0                             ! Include Bacterioplakton: 0 = No; 1 = Yes
     !   this%nbac = 3  
     !   ALLOCATE (this%abac(this%nbac))
     !   this%abac(:) = (/ 1, 1, 2 /)                 ! Bacterioplankton functional groups (= 1 -> Aerobic; = 2 -> Anaerobic)
     !   ALLOCATE (this%sDBacW0(this%nbac))
     !   ALLOCATE (this%sDBacS0(this%nbac))
     !   this%sDBacW0(:) = (/ 0.03, 0.10, 0.10 /)     ! Initial condition of Bacterioplakton functional groups in water
     !   this%sDBacS0(:) = (/ 0.20, 0.10, 0.10 /)     ! Initial condition of Bacterioplakton functional groups in sediment
     !   ! 2.2.2 Parameters of Bacterioplakton
     !   ALLOCATE (this%cThetaBac(this%nbac))
     !   ALLOCATE (this%MuBac(this%nbac))
     !   ALLOCATE (this%hAssBac(this%nbac))
     !   ALLOCATE (this%kResbBac(this%nbac))
     !   ALLOCATE (this%fBac2CO2(this%nbac))
     !   ALLOCATE (this%hPO4UptBac(this%nbac))
     !   ALLOCATE (this%hNH4UptBac(this%nbac))
     !   ALLOCATE (this%kMortBac(this%nbac))
     !   ALLOCATE (this%cCDBacRef(this%nbac))
     !   ALLOCATE (this%cNDBacRef(this%nbac))
     !   ALLOCATE (this%cPDBacRef(this%nbac))
     !   ALLOCATE (this%hBacMin(this%nbac))
     !   this%cThetaBac(:) = (/ 1.08, 1.08, 1.08 /) 
     !   this%hO2Bac = 1.0 !<half-sat. oxygen conc. for BOD
     !   this%hBacMin(:) = (/ 0.05, 0.05, 0.01 /)  !<half-sat. of bacterioplankton for mineralization
     !   this%fbacAn = 0.1
     !   this%kbacAn = 1.0 
     !   this%cCDBacRef(:) = (/ 0.45, 0.45, 0.45 /) 
     !   this%cNDBacRef(:) = (/ 0.09, 0.09, 0.09 /)
     !   this%cPDBacRef(:) = (/ 0.01, 0.01, 0.01 /)
     !   this%MuBac(:)      =  (/ 4.0,  4.0,  4.0 /)
     !   this%hAssBac(:)    =  (/ 0.0,  0.0,  0.0 /)
     !   this%kResbBac(:)   = (/ 0.01, 0.01, 0.01 /)
     !   this%fBac2CO2(:) = (/ 0.9, 0.9, 0.9 /)
     !   this%hPO4UptBac(:) =  (/ 0.003, 0.003, 0.003 /) 
     !   this%hNH4UptBac(:) =  (/ 0.007, 0.007, 0.007 /)
     !   this%kMortBac(:)   =  (/ 0.05,  0.05,  0.05 /)
     !   this%hSedBac(:) = (/ 0.3, 0.3, 0.3 /)
     !   this%fSiBac  = 0.0 
     !   ! Bacteria Preferences
     !   ALLOCATE (this%cPrefBacDom(this%nbac,this%ndom))
     !   this%cPrefBacDom(:,1)   =  (/ 0.8,  0.8,  0.8 /) !Bacteria preference for DOM Labil
     !   !this%cPrefBacDom(:,2)   =  (/ 0.2,  0.2,  0.2 /) !Bacteria preference for DOM Refratory
     !
     !   ! 2.3 Phytoplakton
     !   ! 2.3.1 Initial configuration
     !   this%InclPhyt  = 1                               ! Include Phytoplakton: 0 = No; 1 = Yes
     !   this%nphy = 1
     !   ALLOCATE (this%aphy(this%nphy))
     !   this%aphy(:) = (/ 4, 7, 4 /)                 !(= 1 -> Dinoflagellates ; = 2 -> Freshwater cyanobacteria; = 3 -> Marine/estuarine cyanobacteria; = 4 -> Chlorophytes; = 5 -> Cryptophytes; = 6 -> Marine/estuarine diatoms; = 7 -> Freshwater diatoms
     !   ALLOCATE (this%LightMethodPhyt(this%nphy))
     !   this%LightMethodPhyt = (/ 1, 1, 0 /) 
     !   ALLOCATE (this%sDPhytW0(this%nphy))
     !   ALLOCATE (this%sDPhytS0(this%nphy))
     !   this%sDPhytW0(:) = (/ 0.10, 0.10, 0.10 /)     ! Initial condition of Phytoplakton functional groups in water
     !   this%sDPhytS0(:) = (/ 0.50, 0.50, 0.50 /)     ! Initial condition of Phytoplakton functional groups in sediment
     !   ! 2.3.2 Parameters of Phytoplakton
     !   !Ratios
     !   ALLOCATE (this%cCDPhytOpt(this%nphy))
     !   ALLOCATE (this%cNDPhytOpt(this%nphy))
     !   ALLOCATE (this%cPDPhytOpt(this%nphy))
     !   ALLOCATE (this%cSiDPhytOpt(this%nphy))
     !   this%cCDPhytOpt(:)  = (/ 0.86, 0.86, 0.86 /)
     !   this%cNDPhytOpt(:)  = (/ 0.13, 0.13, 0.13 /)
     !   this%cPDPhytOpt(:)  = (/ 0.008, 0.008, 0.008 /)
     !   this%cSiDPhytOpt(:) = (/ 0.38, 0.00, 0.00 /)
     !
     !   !Primary Production
     !   ALLOCATE (this%cExtSpPhyt(this%nphy))
     !   ALLOCATE (this%hLRefPhyt(this%nphy))
     !   ALLOCATE (this%cLOptRefPhyt(this%nphy))
     !   ALLOCATE (this%cMuMaxPhyt(this%nphy))
     !   ALLOCATE (this%cChDPhytMax(this%nphy))
     !   ALLOCATE (this%cChDPhytMin(this%nphy))
     !   this%fPAR  = 0.45									 
     !   this%cExtSpPhyt(:)   =  (/ 0.0,  0.0,  0.0 /) !(/ 0.20,  0.25,  0.25 /)
     !   this%hLRefPhyt(:)    =   (/ 0.0,  17.0,   0.0 /)
     !   this%cLOptRefPhyt(:) =  (/ 20.0,   0.0,  13.6 /)
     !   this%cMuMaxPhyt(:)   =  (/ 2.00,   1.5,  0.60 /)
     !   this%cChDPhytMax(:)  = (/ 0.006, 0.020, 0.015 /)
     !   this%cChDPhytMin(:)  = (/ 0.001, 0.010, 0.005 /)
     !   this%hSiAssDiat = 	    0.09
     !   !Schecci Disk
     !   this%cPACoefMin = 1.5								 
     !   this%cPACoefMax = 2.8							
     !   this%hPACoef = 3	
     !   ! Nutrient Uptake
     !   ALLOCATE (this%cSigTmPhyt(this%nphy))
     !   ALLOCATE (this%cTmOptPhyt(this%nphy))
     !   ALLOCATE (this%cVCUptMaxPhyt(this%nphy))
     !   ALLOCATE (this%cVNUptMaxPhyt(this%nphy))
     !   ALLOCATE (this%cVPUptMaxPhyt(this%nphy))
     !   ALLOCATE (this%cVSiUptMaxPhyt(this%nphy))
     !   ALLOCATE (this%cCDPhytMax(this%nphy))
     !   ALLOCATE (this%cCDPhytMin(this%nphy))
     !   ALLOCATE (this%cNDPhytMax(this%nphy))
     !   ALLOCATE (this%cNDPhytMin(this%nphy))
     !   ALLOCATE (this%cPDPhytMax(this%nphy))
     !   ALLOCATE (this%cPDPhytMin(this%nphy))
     !   ALLOCATE (this%cSiDPhytMax(this%nphy))
     !   ALLOCATE (this%cSiDPhytMin(this%nphy))
     !   ALLOCATE (this%cAffCUptPhyt(this%nphy))
     !   ALLOCATE (this%cAffNUptPhyt(this%nphy))
     !   ALLOCATE (this%cAffPUptPhyt(this%nphy))
     !   ALLOCATE (this%cAffSiUptPhyt(this%nphy))
     !   this%cSigTmPhyt(:) =         (/ 20,     15,     12 /)
     !   this%cTmOptPhyt(:) =         (/ 18,     25,     25 /)
     !   this%cVCUptMaxPhyt(:) =    (/ 0.86,   0.86,   0.86 /) 
     !   this%cVNUptMaxPhyt(:) =    (/ 0.15,   0.15,   0.15 /)
     !   this%cVPUptMaxPhyt(:) =    (/ 0.01,   0.01,   0.04 /)
     !   this%cVSiUptMaxPhyt(:) =   (/ 0.20,   0.00,   0.00 /)
     !   this%cCDPhytMax(:) =       (/ 0.60,   0.60,   0.60 /)
     !   this%cCDPhytMin(:) =       (/ 0.40,   0.40,   0.4 /)
     !   this%cNDPhytMax(:) =       (/ 0.07,   0.50,   0.20 /)
     !   this%cNDPhytMin(:) =       (/ 0.01,   0.02,   0.03 /)
     !   this%cPDPhytMax(:) =      (/ 0.005,  0.015,  0.025 /)
     !   this%cPDPhytMin(:) =     (/ 0.0005, 0.0015, 0.0025 /)
     !   this%cSiDPhytMax(:) =      (/ 0.20,   0.00,   0.00 /)
     !   this%cSiDPhytMin(:) =      (/ 0.04,   0.00,   0.00 /)
     !   this%cAffCUptPhyt(:) =     (/ 0.00,   0.00,   0.00 /)
     !   this%cAffNUptPhyt(:) =     (/ 0.08,   0.20,   0.20 /)
     !   this%cAffPUptPhyt(:) =     (/ 0.20,   0.20,   0.80 /)
     !   this%cAffSiUptPhyt(:) =    (/ 0.50,   0.00,   0.00 /)
     !   !Respiration
     !   ALLOCATE (this%kDRespbPhyt(this%nphy))
     !   ALLOCATE (this%alfap(this%nphy))
     !   this%kDRespbPhyt(:) =      (/ 0.01,   0.07,   0.03 /)
     !   this%alfap(:) =            (/ 0.00,   0.00,   0.00 /)
     !   !Sedimentation
     !   ALLOCATE (this%cVSetPhyt(this%nphy))
     !   this%cVSetPhyt(:) =        (/ 0.15,   0.15,   0.10 /) !SetMethodFlag = 0
     !   ALLOCATE (this%ppPhyt(this%nphy)) !SetMethodFlag = 1
     !   ALLOCATE (this%diPhyt(this%nphy)) !SetMethodFlag = 1
     !   this%ppPhyt(:) =        (/ 1300,   1300,   1300 /) !SetMethodFlag = 1
     !   this%diPhyt(:) =        (/ 0.000005,   0.000005,   0.000005 /) !SetMethodFlag = 1
     !   !Lise
     !   ALLOCATE (this%cLisPhytW(this%nphy))
     !   !IF (.NOT. ALLOCATED(hLisNutPhyt))    ALLOCATE (hLisNutPhyt(nphy))
     !   ALLOCATE (this%cLisPhytS(this%nphy))
     !   this%cLisPhytW(:)   = (/ 0.01, 0.1, 0.1 /) ! mortality considering predators off
     !   this%cLisPhytS(:)   = (/ 0.05, 0.1, 0.1 /) ! mortality considering predators off
     !   !Exudation
     !   !IF (.NOT. ALLOCATED(cExdPhytW))    ALLOCATE (cExdPhytW(nphy))
     !   !cExdPhytW(:)   = (/ 0.10, 0.10, 0.10 /)
     !   !If (InclMatOrgSplit == 1) Then
     !   !    fDExdDomPhyt  = 1.00 !Labile fractions of died_Algae (Autoctone Organic Matter)
     !   !Else
     !   !    fDExdDomPhyt  = 1.00 !Labile fractions of died_Algae (Autoctone Organic Matter)
     !   !EndIf
     !
     !   ! 2.4 Zooplankton
     !   ! 2.4.1 Initial configuration
     !   this%InclZoo  = 0                             ! Include Zooplankton: 0 = No; 1 = Yes
     !   this%InclZooDist = 1                              ! Include specimen distinction among Zooplankton: 0 = No; 1 = Yes
     !   this%nzoo = 3
     !   ALLOCATE (this%azoo(this%nzoo))
     !   this%azoo(:) = (/ 1, 1, 2 /)                 !(= 1 -> Micro ; = 2 -> Meso; = 3 -> Macro; = 4 -> General)
     !   ALLOCATE (this%sDZoo0(this%nzoo))
     !   this%sDZoo0(:) = (/ 0.01, 0.01, 0.01 /)      ! Initial condition of Zooplankton functional groups
     !   ! 2.3.2 Parameters of Zooplankton
     !   !Ratios
     !   ALLOCATE (this%cCDZooRef(this%nzoo))
     !   ALLOCATE (this%cNDZooRef(this%nzoo))
     !   ALLOCATE (this%cPDZooRef(this%nzoo))
     !   this%cCDZooRef(:) =         (/  0.45,  0.45,   0.45 /)
     !   this%cNDZooRef(:) =         (/  0.08,  0.08,   0.08 /)
     !   this%cPDZooRef(:) =         (/ 0.005, 0.005,  0.005 /)  
     !   !Parameters
     !   ALLOCATE (this%fDAssZoo(this%nzoo))
     !   ALLOCATE (this%cSigTmZoo(this%nzoo))
     !   ALLOCATE (this%cTmOptZoo(this%nzoo))
     !   ALLOCATE (this%cDCarrZoo(this%nzoo))
     !   ALLOCATE (this%hFilt(this%nzoo))
     !   ALLOCATE (this%cFiltMax(this%nzoo))
     !   ALLOCATE (this%kDRespZoo(this%nzoo))
     !   !ALLOCATE (this%fDMessyZoo(this%nzoo))
     !   ALLOCATE (this%kPelZoo(this%nzoo))
     !   ALLOCATE (this%kMortZoo(this%nzoo))
     !   this%fDAssZoo(:) =          (/  0.35,  0.35,  0.35 /)
     !   this%cSigTmZoo(:) =          (/  13.0,  13.0,  13.0 /)
     !   this%cTmOptZoo(:) =          (/  25.0,  25.0,  25.0 /)
     !   this%cDCarrZoo(:) =          (/  15.0,  15.0,   15.0 /)
     !   this%hFilt(:) =              (/  0.50,  0.30,   0.30 /)
     !   this%cFiltMax(:) =           (/  0.50,  2.00,  10.00 /)
     !   this%kDRespZoo(:) =          (/  0.01,  0.02,   0.02 /)
     !   !this%fDMessyZoo(:) =         (/  0.20,  0.20,   0.20 /)
     !   this%kPelZoo(:) =           (/  0.04,  0.01,   0.01 /)
     !   this%kMortZoo(:) =           (/  0.04,  0.04,   0.04 /)
     !   
     !   
     !   ! 2.5 Zoobenthos
     !   this%InclBent  = 0
     !   this%nben = 3
     !   ALLOCATE (this%aben(this%nben))
     !   this%aben(:) = (/ 4, 2, 2 /)                 !(= 1 -> Micro ; = 2 -> Meso; = 3 -> Macro; = 4 -> General)
     !   ALLOCATE (this%sDBent0(this%nben))
     !   this%sDBent0(:) = (/ 0.1, 0.1, 0.1 /)
     !   
     !   !Parametros
     !   ALLOCATE (this%cSigTmBent(this%nben))
     !   ALLOCATE (this%cTmOptBent(this%nben))
     !   ALLOCATE (this%hDFoodBent(this%nben))
     !   ALLOCATE (this%kDAssBent(this%nben))
     !   ALLOCATE (this%cDCarrBent(this%nben))
     !   ALLOCATE (this%fDAssBent(this%nben))
     !   ALLOCATE (this%kDRespBent(this%nben))
     !   ALLOCATE (this%kPelBent(this%nben))
     !   ALLOCATE (this%kMortBent(this%nben))
     !   ALLOCATE (this%cCDBentRef(this%nben))
     !   ALLOCATE (this%cNDBentRef(this%nben))
     !   ALLOCATE (this%cPDBentRef(this%nben))
     !   
     !   
     !   this%cSigTmBent(:) = (/ 16.,  16.,  16. /) !< temperature constant of zoobenthos (sigma in Gaussian curve)
     !   this%cTmOptBent(:) = (/ 25.,  25.,  25. /) !< optimum temp. of zoobenthos
     !   this%hDFoodBent(:) = (/ 200.,  200.,  200. /) !< half-saturating food for zoobenthos (g/m2)
     !   this%kDAssBent(:) = (/ 0.1,  0.1,  0.1 /) !< maximum assimilation rate (1/day)
     !   this%cDCarrBent(:) = (/ 10.,  10.,  10. /) !< Carring capacity (mgDW/L)
     !   this%fDAssBent(:) = (/ 0.3,  0.3,  0.3 /) !< C ass. efficiency of zoobenthos (-)
     !   
     !   !this%hFiltBent(:) = (/ 0.5,  0.5,  0.5 /) !< half-sat. food conc. for filtering (mgDW/L)
     !   !this%cFiltMaxBent(:) = (/ 2.5,  2.5,  2.5 /) !< Max half-sat. food conc. for filtering (mgDW/L)
     !   this%kDRespBent(:) = (/ 0.005,  0.005,  0.005 /) !< maint. respiration constant of zoobenthos
     !   !this%fDMessyBent(:) = (/ 0.00,  0.00,  0.00 /) !< 
     !   this%kPelBent(:) = (/ 0.01,  0.01,  0.01 /) !< 
     !   this%kMortBent(:) = (/ 0.005,  0.005,  0.005 /) !< mortality constant of zoobenthos
     !   this%cCDBentRef(:) = (/ 0.45,  0.45,  0.45 /) !<
     !   this%cNDBentRef(:) = (/ 0.08,  0.08,  0.08 /) !<
     !   this%cPDBentRef(:) = (/ 0.01,  0.01,  0.01 /) !<
     !
     !   ALLOCATE (this%cPrefBentPhyt(this%nben,this%nphy))
     !   this%cPrefBentPhyt(:,1) = (/ 0.9,  0.9,   0.9 /)
     !   !this%cPrefBentPhyt(:,2) = (/ 0.9,  0.9,   0.9 /)
     !   !this%cPrefBentPhyt(:,3) = (/ 0.9,  0.75,   0.9 /)
     !   ALLOCATE (this%cPrefBentBac(this%nben,this%nbac))
     !   this%cPrefBentBac(:,1) = (/ 0.1,  0.1,   0.1 /)
     !   this%cPrefBentBac(:,2) = (/ 0.1,  0.1,   0.1 /)
     !   this%cPrefBentBac(:,3) = (/ 0.1,  0.1,   0.1 /)
     !   ALLOCATE (this%cPrefBentPom(this%nben,this%npom))
     !   this%cPrefBentPom(:,1) = (/ 0.3,  0.4,   0.3 /)
     !   !this%cPrefBentPom(:,2) = (/ 0.3,  0.4,   0.3 /)        
     !   
     !   
     !   
     !   ! 2.X Fish
     !   ! 2.X.1 Initial configuration
     !   this%InclFish  = 0
     !   this%nfish = 3
     !   ALLOCATE (this%afish(this%nfish))
     !   this%afish(:) = (/ 2, 3, 3 /)                 ! Fish functional groups (= 1 -> Piscivorous ; = 2 -> Omnivorous; = 3 -> Planktivorous)
     !   ALLOCATE (this%sDFiAd0(this%nfish))
     !   this%sDFiAd0(:) = (/ 2.0, 0.1, 3.0 /)
     !   ALLOCATE (this%sDFiJv0(this%nfish))
     !   this%sDFiJv0(:) = (/ 0.1, 0., 0.1 /)
     !
     !   ! 2.X.2 Parameters of Fish
     !   ALLOCATE (this%kMigrFiAd(this%nfish))
     !   ALLOCATE (this%cDFiAdIn(this%nfish))
     !   ALLOCATE (this%kMigrFiJv(this%nfish))
     !   ALLOCATE (this%cDFiJvIn(this%nfish))
     !
     !   !IF (.NOT. ALLOCATED(MigrOutTime))    ALLOCATE (MigrOutTime(nfish))
     !   !IF (.NOT. ALLOCATED(cDayMigrInFish))    ALLOCATE (cDayMigrInFish(nfish))
     !   !IF (.NOT. ALLOCATED(cDayMigrOutFish))    ALLOCATE (cDayMigrOutFish(nfish))
     !   ALLOCATE (this%cDayReprFish(this%nfish))
     !   ALLOCATE (this%fReprFish(this%nfish))
     !   ALLOCATE (this%fAgeFish(this%nfish))
     !   ALLOCATE (this%cRelVegFish(this%nfish))
     !   ALLOCATE (this%cSigTmFish(this%nfish))
     !   ALLOCATE (this%cTmOptFish(this%nfish))
     !   ALLOCATE (this%hDFiAd(this%nfish))
     !   ALLOCATE (this%hDFiJv(this%nfish))
     !   ALLOCATE (this%kDAssFiAd(this%nfish))
     !   ALLOCATE (this%kDAssFiJv(this%nfish))
     !   ALLOCATE (this%fDAssFiAd(this%nfish))
     !   ALLOCATE (this%fDAssFiJv(this%nfish))
     !   ALLOCATE (this%kHarvFiAd(this%nfish))
     !   ALLOCATE (this%kHarvFiJv(this%nfish))
     !   ALLOCATE (this%kDRespFiAd(this%nfish))
     !   ALLOCATE (this%kDRespFiJv(this%nfish))
     !   !ALLOCATE (this%fDMessyFiAd(this%nfish))
     !   !ALLOCATE (this%fDMessyFiJv(this%nfish))
     !   ALLOCATE (this%kPelFiAd(this%nfish))
     !   ALLOCATE (this%kPelFiJv(this%nfish))
     !   ALLOCATE (this%kMortFiAd(this%nfish))
     !   ALLOCATE (this%kMortFiJv(this%nfish))
     !   ALLOCATE (this%cDCarrFish(this%nfish))
     !   ALLOCATE (this%cCDFishRef(this%nfish))
     !   ALLOCATE (this%cNDFishRef(this%nfish))
     !   ALLOCATE (this%cPDFishRef(this%nfish))
     !   this%kMigrFiAd(:)            = (/ 0.0,  0.001,  0.001 /) !< Adult fish migration rate (day-1)
     !   this%cDFiAdIn(:)            = (/ 0.005,  0.005,  0.005 /) !< External adult fish density (gD/m2)
     !   this%kMigrFiJv(:)            = (/ 0.0,  0.001,  0.001 /) !< Juvenile fish migration rate (day-1)
     !   this%cDFiJvIn(:)            = (/ 0.005,  0.005,  0.005 /) !< External juvenile fish density (gD/m2)
     !   this%cDayReprFish(:)           = (/ 0.00,  0.00,   120. /) 
     !   this%fReprFish(:)              = (/ 0.00,  0.00,  0.02 /)
     !   this%fAgeFish(:)               = (/ 0.00,   0.5,  0.00 /)
     !   this%cRelVegFish(:)            = (/ 0.03,   0.0, 0.009 /)
     !   this%cSigTmFish(:)             = (/   10.,    10.,    10. /)
     !   this%cTmOptFish(:)             = (/   25.,    25.,    25. /)
     !   this%hDFiAd(:)                 = (/ 5.00,  1.25,  5.00 /)
     !   this%hDFiJv(:)                 = (/ 5.00,  1.25,  5.00 /)
     !   this%kDAssFiAd(:)              = (/ 0.15,  0.12,  0.06 /)
     !   this%kDAssFiJv(:)              = (/ 0.15,  0.12,  0.06 /)
     !   this%fDAssFiAd(:)              = (/ 0.5,  0.5,  0.5 /)
     !   this%fDAssFiJv(:)              = (/ 0.35,  0.35,  0.35 /)
     !   this%kHarvFiAd(:)              = (/ 0.00,  0.00,  0.00 /)
     !   this%kHarvFiJv(:)              = (/ 0.00,  0.00,  0.00 /)
     !   this%kDRespFiAd(:)             = (/ 0.004, 0.010, 0.004 /) 
     !   this%kDRespFiJv(:)             = (/ 0.004, 0.010, 0.004 /) 
     !   this%kPelFiAd(:)              = (/ 0.02,  0.02,  0.02 /)
     !   this%kPelFiJv(:)              = (/ 0.02,  0.02,  0.02 /)
     !   this%kMortFiAd(:)              = (/ 0.001, 0.001, 0.001 /)
     !   this%kMortFiJv(:)              = (/ 0.001, 0.001, 0.001 /)
     !   this%cDCarrFish(:)              = (/ 15, 15, 15 /)
     !   this%cCDFishRef(:)             = (/ 0.45,  0.45,  0.45 /) 
     !   this%cNDFishRef(:)             = (/ 0.08,  0.08,  0.08 /)
     !   this%cPDFishRef(:)             = (/ 0.005, 0.005, 0.005 /)    
     !   ! Fish Preferences
     !   ALLOCATE (this%cPrefFiAdFiAd(this%nfish,this%nfish))
     !   this%cPrefFiAdFiAd(:,1) = (/ 0.0,  0.0,   0.0 /)
     !   this%cPrefFiAdFiAd(:,2) = (/ 0.0,  0.0,   0.0 /)
     !   this%cPrefFiAdFiAd(:,3) = (/ 0.9,  0.75,   0.0 /)
     !   ALLOCATE (this%cPrefFiAdFiJv(this%nfish,this%nfish))
     !   this%cPrefFiAdFiJv(:,1) = (/ 0.0,  0.0,   0.0 /)
     !   this%cPrefFiAdFiJv(:,2) = (/ 0.0,  0.0,   0.0 /)
     !   this%cPrefFiAdFiJv(:,3) = (/ 0.9,  0.75,   0.0 /)
     !   ALLOCATE (this%cPrefFiAdZoo(this%nfish,this%nzoo))
     !   this%cPrefFiAdZoo(:,1) = (/ 0.0,  0.0,   0.0 /)
     !   this%cPrefFiAdZoo(:,2) = (/ 0.9,  0.9,   0.9 /)
     !   this%cPrefFiAdZoo(:,3) = (/ 0.9,  0.75,   0.9 /)
     !   ALLOCATE (this%cPrefFiJvZoo(this%nfish,this%nzoo))
     !   this%cPrefFiJvZoo(:,1) = (/ 0.0,  0.0,   0.0 /)
     !   this%cPrefFiJvZoo(:,2) = (/ 0.9,  0.9,   0.9 /)
     !   this%cPrefFiJvZoo(:,3) = (/ 0.9,  0.75,   0.9 /)
     !   ALLOCATE (this%cPrefFiAdBent(this%nfish,this%nben))
     !   this%cPrefFiAdBent(:,1) = (/ 0.0,  0.0,   0.0 /)
     !   this%cPrefFiAdBent(:,2) = (/ 0.9,  0.9,   0.9 /)
     !   this%cPrefFiAdBent(:,3) = (/ 0.9,  0.75,   0.9 /)
     !   ALLOCATE (this%cPrefFiJvBent(this%nfish,this%nben))
     !   this%cPrefFiJvBent(:,1) = (/ 0.0,  0.0,   0.0 /)
     !   this%cPrefFiJvBent(:,2) = (/ 0.9,  0.9,   0.9 /)
     !   this%cPrefFiJvBent(:,3) = (/ 0.9,  0.75,   0.9 /)
     !   ALLOCATE (this%cPrefFiAdPhyt(this%nfish,this%nphy))
     !   this%cPrefFiAdPhyt(:,1) = (/ 0.1,  0.1,   0.1 /)
     !   !this%cPrefFiAdPhyt(:,2) = (/ 0.1,  0.1,   0.1 /)
     !   !this%cPrefFiAdPhyt(:,3) = (/ 0.1,  0.1,   0.1 /)
     !   ALLOCATE (this%cPrefFiJvPhyt(this%nfish,this%nphy))
     !   this%cPrefFiJvPhyt(:,1) = (/ 0.8,  0.8,   0.8 /)
     !   !this%cPrefFiJvPhyt(:,2) = (/ 0.8,  0.8,   0.8 /)
     !   !this%cPrefFiJvPhyt(:,3) = (/ 0.8,  0.8,   0.8 /)
     !   ALLOCATE (this%cPrefFiAdPom(this%nfish,this%npom))
     !   this%cPrefFiAdPom(:,1) = (/ 0.1,  0.1,   0.1 /)
     !   ALLOCATE (this%cPrefFiJvPom(this%nfish,this%npom))
     !   this%cPrefFiJvPom(:,1) = (/ 0.8,  0.8,   0.8 /)
     !   
     !   !this%cPrefFiAdPom(:,2) = (/ 0.1,  0.1,   0.1 /)
     !   !cPrefFiAdPom(:,3) = (/ 0.1,  0.1,   0.1 /)
     !
     !   ! 2.X Resuspension
     !   this%cSuspRef      = 0.50 !<reference suspended matter function		
     !   this%cSuspMin      = 6.10 !<minimum value of logistic function		
     !   this%cSuspMax      = 25.2 !<maximum value of logistic function		
     !   this%cSuspSlope    = 2.10 !<slope of logistic function	
     !   this%hDepthSusp    = 2.00 !<‘half-sat. value’ of depth in logistic function		
     !   this%cFetchRef     = 1000.0 !<‘reference’ fetch		
     !   this%kTurbFish     = 1.00 !<relative resuspension by adult fish browsing
     !   this%kVegResus     = 0.01 !<relative resuspension reduction per g vegetation
     !
     !   
     !   ! Nutrients
     !   !3.1 Nitrification
	    !this%cThetaNitr = 1.08 !Constante de temperatura exponecial de nitrificação
	    !this%kNitrW = 0.1 !<nitrification rate constant in water [d-1]
     !   this%kNitrS = 1.0 !<nitrification rate constant in sediment [d-1]
	    !this%hO2Nitr = 2.0  !meia-saturação quadratica da conc. de O2 para Nitrificação mgO2/l
     !   !3.1 Desnitrification
	    !this%hNO3Denit = 2. !<quadratic half-sat. NO3 conc. for denitrification (mgN/l)    
     !   !3.2 Absorbed P in sediment,oxygen dependent
     !   this%cRelPAdsD=0.00003 !<max. P adsorption per g DW (gP/gD) 
     !   this%cRelPAdsFe=0.065 !<max. P adsorption per g Fe (gP/gFe) 
     !   this%cRelPAdsAl=0.134 !<max. P adsorption per g Al (gP/gAl) 
     !   this%fFeDIM=0.01 !<Fe content of inorg. matter (gFe/gD) 
     !   this%fAlDIM=0.01 !<Al content of inorg. matter (gAl/gD) 
     !   this%fRedMax=0.9 !< max. reduction factor of P adsorption affinity
     !   this%cKPAdsOx=0.6 !<P adsorption affinity at oxidized conditions 
     !   this%kPSorp=0.05 !taxa de absorção do fósforo (1/d)     
     !   !3.3 P immobilization    
     !   this%kPChemPO4=0.03 !<chem. PO4 loss rate (1/d) 
     !   this%cPO4Max=2.0 !<max. SRP conc. in pore water
     !   !3.4 Reaeration        
     !   this%cThetaReaer = 1.1
     !   
     !   
     !   
     !   !3. Sediment
     !   this%nsed = 2 !<number of layers in the sediment
     !   this%cDepthS = 0.1              !<sediment depth (m)
     !   this%cRhoIM = 2500000.0 !<g/m3 solid ; density of sediment IM
     !   this%cRhoOM = 1400000.0 !<g/m3 ; density of sediment detritus
     !   !3.1 Difussion in the sediment
     !   this%fDepthDifS = 0.5					!< Fraction of Sediment's depth with diffusion (-) 
	    !this%cThetaDif = 1.02					!< Temperature coefficient for diffusion
     !   Allocate (this%kPDifDom(this%ndom))
	    !this%kPDifDom(:) = (/ 0.000072,  0.000072 /) !<DOM diffusion constant (m²/d)
	    !this%kCDifDic = 0.000072 !<DIC diffusion constant (m²/d)
	    !this%kPDifPO4 = 0.000072 !<PO4 diffusion constant (m²/d)
	    !this%kNDifNO3 = 0.000086 !<NO3 diffusion constant (m²/d)
	    !this%kNDifNH4 = 0.000072 !<NH4 diffusion constant (m²/d)
 	   ! this%kO2Dif = 0.000026	!<NH4 diffusion constant (m²/d)
 	   ! this%cTurbDifNut = 1.0  !<bioturbation factor for diffusion for nutrients
 	   ! this%cturbDifO2 = 1.0   !<bioturbation factor for diffusion for O2 (-)
     !
     !
     !   ! 3.4 Erosion
     !   this%cDErosTot = 0.1   !< Erosion input (tentative) (g/m2/day)
     !   this%fDOrgSoil = 0.1 ! - ; fraction soil organic matter
     !   this%fSedErosIM = 0.95 !< instantly sedimentating fraction of IM
     !   
     !   !3.5 Fractions
	    !this%fDTotS0 = 0.3              !<initial dry-weight fraction in sediment (g solid g-1 sediment)
	    !this%fDOrgS0 = 0.001              !<initial organic fraction of sediment DW (g residuos de MO/g solido)
	    !this%fDLabilS0 =0.3             !<initial labile fraction in the sediment OM
     !   this%fPAdsS0   = 0.99 !< initial adsorbed fraction of inorg. P in sed.   
     !   this%fPInorgS0 = 0.000001 !< initial inorg. P fraction in sed
	    !this%fLutum	= 0.2					!<fraction lutum in the inorganic matter (g/g)
	    !this%fLutumRef = 0.2				!<reference lutum fraction (g/g)
     !   
     !   !3.6 Macrophytes
     !   this%InclMacr = 0
     !   this%nmac = 3
     !   ALLOCATE (this%amac(this%nmac ))
     !   this%amac(:) = (/ 1, 4, 5 /)                 !(= 1 -> Elodea ; = 2 -> Ceraphytes; = 3 -> Helophytes; = 4 -> Nymphytes; = 5 -> General)
     !   ALLOCATE (this%sDMac0(this%nmac))
     !   this%sDMac0(:) = (/ 10., 10., 10. /)
     !   ALLOCATE (this%cCDMac0(this%nmac))
     !   ALLOCATE (this%cNDMac0(this%nmac))
     !   ALLOCATE (this%cPDMac0(this%nmac))
     !   this%cCDMac0(:) = (/ 0.45, 0.45, 0.45 /)
     !   this%cNDMac0(:) = (/ 0.2, 0.2, 0.2 /)
     !   this%cPDMac0(:) = (/ 0.02, 0.02, 0.02 /)
     !   
     !   !2)Parametros(mg)
 	   ! !Produção primária
     !   ALLOCATE (this%cQ10ProdMac(this%nmac))
     !   ALLOCATE (this%cMuMaxMac(this%nmac))
     !   ALLOCATE (this%hLRefMac(this%nmac))
     !   ALLOCATE (this%cTmInitVeg(this%nmac))
     !   ALLOCATE (this%fRootVegWin(this%nmac))
     !   ALLOCATE (this%fRootVegSum(this%nmac))
     !   ALLOCATE (this%fEmergVeg(this%nmac))
     !   ALLOCATE (this%fFloatVeg(this%nmac))
     !   ALLOCATE (this%cDLayerVeg(this%nmac))
     !   ALLOCATE (this%cDCarrMac(this%nmac))
     !   ALLOCATE (this%cCovSpVeg(this%nmac))
     !   ALLOCATE (this%cCDMacMax(this%nmac))
     !   ALLOCATE (this%cNDMacMax(this%nmac))
     !   ALLOCATE (this%cPDMacMax(this%nmac))
     !   ALLOCATE (this%cCDMacMin(this%nmac))
     !   ALLOCATE (this%cNDMacMin(this%nmac))
     !   ALLOCATE (this%cPDMacMin(this%nmac))
     !   ALLOCATE (this%kDRespMac(this%nmac))
     !   ALLOCATE (this%cQ10RespMac(this%nmac))
     !   ALLOCATE (this%kMortMac(this%nmac))
     !   ALLOCATE (this%cVCUptMaxMac(this%nmac))
     !   ALLOCATE (this%cVNUptMaxMac(this%nmac))
     !   ALLOCATE (this%cVPUptMaxMac(this%nmac))
     !   ALLOCATE (this%cAffCUptMac(this%nmac))
     !   ALLOCATE (this%cAffNUptMac(this%nmac))
     !   ALLOCATE (this%cAffPUptMac(this%nmac))
     !   ALLOCATE (this%fSedUptVegMax(this%nmac))
     !   ALLOCATE (this%fSedUptVegCoef(this%nmac))
     !   ALLOCATE (this%fSedUptVegExp(this%nmac))
     !   ALLOCATE (this%cPrefMacBird(this%nmac))
     !   ALLOCATE (this%hDMacBird(this%nmac))
     !   ALLOCATE (this%cDGrazPerBird(this%nmac))
     !   ALLOCATE (this%fDAssBird(this%nmac))
     !   
     !   this%cQ10ProdMac(:) = (/ 1.2,  1.2,  1.2 /) !< temperature quotient of production
     !   this%cMuMaxMac(:) = (/ 0.2,  0.2,  0.2 /) !< maximum growth rate of vegetation at 20oC (g/g shoot/day)
     !   this%hLRefMac(:) = (/ 17.,  17.,  17. /) !< half-sat. light at 20 oC (W/m2 PAR)
     !   this%cTmInitVeg(:) = (/ 9.0,  9.0,  9.0 /) !< temperature for initial growth
     !   this%fRootVegWin(:) = (/ 0.6,  0.6,  0.6 /) !< root fraction outside growing season
     !   this%fRootVegSum(:) = (/ 0.1,  0.1,  0.1 /) !<root fraction inside growing season
     !   this%fEmergVeg(:) = (/ 0.,  0.,  0. /) !<g emergent / g shoot ; emergent fraction of shoot
     !   this%fFloatVeg(:) = (/ 0.,  0.,  0. /) !<g floating / g shoot ; floating fraction of shoot
     !   this%cDLayerVeg(:) = (/ 0.,  0.,  0. /) !<biomass of a single layer floating leaves (gD/m2)
     !   this%cDCarrMac(:) = (/ 400., 400., 400. /) !<max. vegetation standing crop (gDW/m2)
     !   this%cCovSpVeg(:) = (/ 0.5, 0.5, 0.5 /) !<% cover per gD/m2 ; specific cover
     !   this%cCDMacMax(:) = (/ 0.6, 0.6, 0.6 /) !<Maximum C/D ratio vegetation
     !   this%cCDMacMin(:) = (/ 0.4, 0.4, 0.4 /) !<Minimum C/D ratio vegetation
     !   this%cNDMacMax(:) = (/0.035, 0.035, 0.035 /) !<Maximum N/D ratio vegetation
     !   this%cNDMacMin(:) = (/0.01, 0.01, 0.01 /) !<Minimum N/D ratio vegetation
     !   this%cPDMacMax(:) = (/ 0.0035, 0.0035, 0.0035 /) !<Maximum P/D ratio vegetation
     !   this%cPDMacMin(:) = (/ 0.0008, 0.0008, 0.0008 /) !<Minimum P/D ratio vegetation
     !   this%kDRespMac(:) = (/ 0.02, 0.02, 0.02 /) !<dark respiration rate of vegetation (1/day)
     !   this%cQ10RespMac(:) = (/ 2.0, 2.0, 2.0 /) !<temperature quotient of respiration
     !   this%kMortMac(:) = (/ 0.010, 0.010, 0.010 /) !<vegetation mortality rate (1/day)
     !   this%cVCUptMaxMac(:) = (/ 0.4, 0.4, 0.4 /) !<maximum C uptake capacity of vegetation (mgC/mgD/day)
     !   this%cVNUptMaxMac(:) = (/ 0.10, 0.10, 0.10 /) !<maximum N uptake capacity of vegetation (mgN/mgD/day)
     !   this%cVPUptMaxMac(:) = (/ 0.010, 0.010, 0.010 /) !<maximum P uptake capacity of vegetation (mgP/mgD/day)
     !   this%cAffCUptMac(:) = (/ 0.10, 0.10, 0.10 /) !<initial C uptake affinity vegetation (L/mgD/day) 
     !   this%cAffNUptMac(:) = (/ 0.20, 0.20, 0.20 /) !<initial N uptake affinity vegetation (L/mgD/day) 
     !   this%cAffPUptMac(:) = (/ 0.20, 0.20, 0.20 /) !<initial P uptake affinity vegetation (L/mgD/day) 
     !   this%fSedUptVegMax(:) = (/ 0.998, 0.998, 0.998 /) !< maximum_sediment_fraction_of_nutrient_uptake 
     !   this%fSedUptVegCoef(:) = (/ 2.66 , 2.66 , 2.66  /) !<sigm._regr._coeff._for_sediment_fraction_of_nutrient_uptake
     !   this%fSedUptVegExp(:) = (/ -0.83, -0.83, -0.83 /) !<exponent_in_sigm._regr._for_sediment_fraction_of_nutrient_uptake
     !   this%cPrefMacBird(:) = (/ 1., 1., 1. /) !<edibility for birds
     !   this%hDMacBird(:) = (/ 5., 5., 5. /) !<half-sat. vegetation biomass (gDW/m2)
     !   this%cDGrazPerBird(:) = (/ 45., 45., 45. /) !<daily grazing of birds (gD/coot/day)
     !   this%fDAssBird(:) = (/ 0.5, 0.5, 0.5 /) !<birds assim. efficiency
       
        
    
    end subroutine
End Module LimnologyVars