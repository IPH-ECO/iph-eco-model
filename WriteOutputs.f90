    !> This subroutine reads the hydrodynamic parameters.
    Subroutine WriteOutputs(simParam,HydroParam,MeshParam,LimnoParam,MeteoParam,iNewton,innerNewton)

    Use SimulationModel
    Use Hydrodynamic
    Use MeshVars
    Use LimnologyVars
    Use ParticleTracking
    Use Meteorological
    Implicit none

    type(SimulationParam) :: simParam
    type(MeshGridParam) :: MeshParam
    type(HydrodynamicParam) :: HydroParam
    type(LimnologyParam) :: LimnoParam
    type(ParticleParam) :: PartParam
    type(MeteorologicalParam) :: MeteoParam

    Character(200):: FileName,Basename
    Integer:: i,c,icell,iElem,iLayer,face,iNewton, innerNewton
    Real:: V,eair,upup,up


    If (simParam%outputTimeStep /= 0) Then
        !If (mod(simParam%time-simParam%IniTime,int(simParam%outputTimeStep*simParam%dt))==0.) Then
        If (mod(simParam%it,simParam%outputTimeStep)==0.or.simParam%it==1) Then
            Call VTKOutput(simParam,HydroParam,MeshParam,LimnoParam)
            
             
            !If (simParam%it == 1) Then  
            !    Basename = 'EtaIphEco'
            !    !Write(FileName,'(i10)') EtaIPHECO 
            !    Open(94,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !    Write(94,'(1I15,100F30.20)') simParam%it, HydroParam%eta(348:360)!hydroParam%eta(540),HydroParam%eta(580),HydroParam%eta(628),HydroParam%eta(692),HydroParam%eta(760),HydroParam%eta(840) !hydroParam%eta(348:360)!
            !    
            !    Basename = 'Continuidade'
            !    Open(104,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !    Write(104,'(1I15,100F30.20)') simParam%it,HydroParam%SumVer,HydroParam%SumVerAcum
            !     
            !    !Basename = 'Time'
            !    !Open(666,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !    !Write(666,'(1I15,100F30.20)') simParam%it,simParam%start(simParam%it),simParam%finish(simParam%it), (simParam%finish(simParam%it)-simParam%start(simParam%it))
            !    
            !    !  
            !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !    
            !    !Gravando densidade
            !    icell = 1
            !    Basename = 'density'
            !    Write(FileName,'(i10)') icell
            !    Open(665,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !    Do iLayer=1,12!iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !        If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !            Write(665,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%sDRhoW(iLayer,348:360)!,HydroParam%sDRhoW(iLayer,2),HydroParam%sDRhoW(iLayer,3),HydroParam%sDRhoW(iLayer,4),HydroParam%sDRhoW(iLayer,5),HydroParam%sDRhoW(iLayer,6),HydroParam%sDRhoW(iLayer,7),HydroParam%sDRhoW(iLayer,8),HydroParam%sDRhoW(iLayer,9),HydroParam%sDRhoW(iLayer,10),HydroParam%sDRhoW(iLayer,11),HydroParam%sDRhoW(iLayer,12),HydroParam%sDRhoW(iLayer,13),HydroParam%sDRhoW(iLayer,14),HydroParam%sDRhoW(iLayer,15),HydroParam%sDRhoW(iLayer,16),HydroParam%sDRhoW(iLayer,17),HydroParam%sDRhoW(iLayer,18),HydroParam%sDRhoW(iLayer,19),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,31),HydroParam%sDRhoW(iLayer,32),HydroParam%sDRhoW(iLayer,33),HydroParam%sDRhoW(iLayer,34),HydroParam%sDRhoW(iLayer,35),HydroParam%sDRhoW(iLayer,36),HydroParam%sDRhoW(iLayer,37),HydroParam%sDRhoW(iLayer,38),HydroParam%sDRhoW(iLayer,39),HydroParam%sDRhoW(iLayer,40),HydroParam%sDRhoW(iLayer,41),HydroParam%sDRhoW(iLayer,42),HydroParam%sDRhoW(iLayer,43),HydroParam%sDRhoW(iLayer,44),HydroParam%sDRhoW(iLayer,45),HydroParam%sDRhoW(iLayer,46),HydroParam%sDRhoW(iLayer,47),HydroParam%sDRhoW(iLayer,48),HydroParam%sDRhoW(iLayer,49),HydroParam%sDRhoW(iLayer,50),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,60),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,61),HydroParam%sDRhoW(iLayer,62),HydroParam%sDRhoW(iLayer,63),HydroParam%sDRhoW(iLayer,64),HydroParam%sDRhoW(iLayer,65),HydroParam%sDRhoW(iLayer,66),HydroParam%sDRhoW(iLayer,67),HydroParam%sDRhoW(iLayer,68),HydroParam%sDRhoW(iLayer,69),HydroParam%sDRhoW(iLayer,70),HydroParam%sDRhoW(iLayer,71),HydroParam%sDRhoW(iLayer,72),HydroParam%sDRhoW(iLayer,73),HydroParam%sDRhoW(iLayer,74),HydroParam%sDRhoW(iLayer,75),HydroParam%sDRhoW(iLayer,76),HydroParam%sDRhoW(iLayer,77),HydroParam%sDRhoW(iLayer,78),HydroParam%sDRhoW(iLayer,79),HydroParam%sDRhoW(iLayer,80)
            !        Else
            !            Write(665,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%sDRhoW(iLayer,348:360)!,HydroParam%sDRhoW(iLayer,2),HydroParam%sDRhoW(iLayer,3),HydroParam%sDRhoW(iLayer,4),HydroParam%sDRhoW(iLayer,5),HydroParam%sDRhoW(iLayer,6),HydroParam%sDRhoW(iLayer,7),HydroParam%sDRhoW(iLayer,8),HydroParam%sDRhoW(iLayer,9),HydroParam%sDRhoW(iLayer,10),HydroParam%sDRhoW(iLayer,11),HydroParam%sDRhoW(iLayer,12),HydroParam%sDRhoW(iLayer,13),HydroParam%sDRhoW(iLayer,14),HydroParam%sDRhoW(iLayer,15),HydroParam%sDRhoW(iLayer,16),HydroParam%sDRhoW(iLayer,17),HydroParam%sDRhoW(iLayer,18),HydroParam%sDRhoW(iLayer,19),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,31),HydroParam%sDRhoW(iLayer,32),HydroParam%sDRhoW(iLayer,33),HydroParam%sDRhoW(iLayer,34),HydroParam%sDRhoW(iLayer,35),HydroParam%sDRhoW(iLayer,36),HydroParam%sDRhoW(iLayer,37),HydroParam%sDRhoW(iLayer,38),HydroParam%sDRhoW(iLayer,39),HydroParam%sDRhoW(iLayer,40),HydroParam%sDRhoW(iLayer,41),HydroParam%sDRhoW(iLayer,42),HydroParam%sDRhoW(iLayer,43),HydroParam%sDRhoW(iLayer,44),HydroParam%sDRhoW(iLayer,45),HydroParam%sDRhoW(iLayer,46),HydroParam%sDRhoW(iLayer,47),HydroParam%sDRhoW(iLayer,48),HydroParam%sDRhoW(iLayer,49),HydroParam%sDRhoW(iLayer,50),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,60),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,61),HydroParam%sDRhoW(iLayer,62),HydroParam%sDRhoW(iLayer,63),HydroParam%sDRhoW(iLayer,64),HydroParam%sDRhoW(iLayer,65),HydroParam%sDRhoW(iLayer,66),HydroParam%sDRhoW(iLayer,67),HydroParam%sDRhoW(iLayer,68),HydroParam%sDRhoW(iLayer,69),HydroParam%sDRhoW(iLayer,70),HydroParam%sDRhoW(iLayer,71),HydroParam%sDRhoW(iLayer,72),HydroParam%sDRhoW(iLayer,73),HydroParam%sDRhoW(iLayer,74),HydroParam%sDRhoW(iLayer,75),HydroParam%sDRhoW(iLayer,76),HydroParam%sDRhoW(iLayer,77),HydroParam%sDRhoW(iLayer,78),HydroParam%sDRhoW(iLayer,79),HydroParam%sDRhoW(iLayer,80)
            !        EndIf
            !    EndDo
            !    
            !    !Gravando temperatura
            !    icell = 1
            !    Basename = 'Temperature'
            !    Write(FileName,'(i10)') icell
            !    Open(664,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !    Do iLayer=1,12!iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !        If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !            Write(664,'(I10,A10,100F30.20)') simParam%it,'Sup',LimnoParam%sDTempW(iLayer,348:360)!,HydroParam%sDRhoW(iLayer,2),HydroParam%sDRhoW(iLayer,3),HydroParam%sDRhoW(iLayer,4),HydroParam%sDRhoW(iLayer,5),HydroParam%sDRhoW(iLayer,6),HydroParam%sDRhoW(iLayer,7),HydroParam%sDRhoW(iLayer,8),HydroParam%sDRhoW(iLayer,9),HydroParam%sDRhoW(iLayer,10),HydroParam%sDRhoW(iLayer,11),HydroParam%sDRhoW(iLayer,12),HydroParam%sDRhoW(iLayer,13),HydroParam%sDRhoW(iLayer,14),HydroParam%sDRhoW(iLayer,15),HydroParam%sDRhoW(iLayer,16),HydroParam%sDRhoW(iLayer,17),HydroParam%sDRhoW(iLayer,18),HydroParam%sDRhoW(iLayer,19),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,31),HydroParam%sDRhoW(iLayer,32),HydroParam%sDRhoW(iLayer,33),HydroParam%sDRhoW(iLayer,34),HydroParam%sDRhoW(iLayer,35),HydroParam%sDRhoW(iLayer,36),HydroParam%sDRhoW(iLayer,37),HydroParam%sDRhoW(iLayer,38),HydroParam%sDRhoW(iLayer,39),HydroParam%sDRhoW(iLayer,40),HydroParam%sDRhoW(iLayer,41),HydroParam%sDRhoW(iLayer,42),HydroParam%sDRhoW(iLayer,43),HydroParam%sDRhoW(iLayer,44),HydroParam%sDRhoW(iLayer,45),HydroParam%sDRhoW(iLayer,46),HydroParam%sDRhoW(iLayer,47),HydroParam%sDRhoW(iLayer,48),HydroParam%sDRhoW(iLayer,49),HydroParam%sDRhoW(iLayer,50),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,60),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,61),HydroParam%sDRhoW(iLayer,62),HydroParam%sDRhoW(iLayer,63),HydroParam%sDRhoW(iLayer,64),HydroParam%sDRhoW(iLayer,65),HydroParam%sDRhoW(iLayer,66),HydroParam%sDRhoW(iLayer,67),HydroParam%sDRhoW(iLayer,68),HydroParam%sDRhoW(iLayer,69),HydroParam%sDRhoW(iLayer,70),HydroParam%sDRhoW(iLayer,71),HydroParam%sDRhoW(iLayer,72),HydroParam%sDRhoW(iLayer,73),HydroParam%sDRhoW(iLayer,74),HydroParam%sDRhoW(iLayer,75),HydroParam%sDRhoW(iLayer,76),HydroParam%sDRhoW(iLayer,77),HydroParam%sDRhoW(iLayer,78),HydroParam%sDRhoW(iLayer,79),HydroParam%sDRhoW(iLayer,80)
            !        Else
            !            Write(664,'(2I10,100F30.20)') simParam%it,iLayer,LimnoParam%sDTempW(iLayer,348:360)!,HydroParam%sDRhoW(iLayer,2),HydroParam%sDRhoW(iLayer,3),HydroParam%sDRhoW(iLayer,4),HydroParam%sDRhoW(iLayer,5),HydroParam%sDRhoW(iLayer,6),HydroParam%sDRhoW(iLayer,7),HydroParam%sDRhoW(iLayer,8),HydroParam%sDRhoW(iLayer,9),HydroParam%sDRhoW(iLayer,10),HydroParam%sDRhoW(iLayer,11),HydroParam%sDRhoW(iLayer,12),HydroParam%sDRhoW(iLayer,13),HydroParam%sDRhoW(iLayer,14),HydroParam%sDRhoW(iLayer,15),HydroParam%sDRhoW(iLayer,16),HydroParam%sDRhoW(iLayer,17),HydroParam%sDRhoW(iLayer,18),HydroParam%sDRhoW(iLayer,19),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,31),HydroParam%sDRhoW(iLayer,32),HydroParam%sDRhoW(iLayer,33),HydroParam%sDRhoW(iLayer,34),HydroParam%sDRhoW(iLayer,35),HydroParam%sDRhoW(iLayer,36),HydroParam%sDRhoW(iLayer,37),HydroParam%sDRhoW(iLayer,38),HydroParam%sDRhoW(iLayer,39),HydroParam%sDRhoW(iLayer,40),HydroParam%sDRhoW(iLayer,41),HydroParam%sDRhoW(iLayer,42),HydroParam%sDRhoW(iLayer,43),HydroParam%sDRhoW(iLayer,44),HydroParam%sDRhoW(iLayer,45),HydroParam%sDRhoW(iLayer,46),HydroParam%sDRhoW(iLayer,47),HydroParam%sDRhoW(iLayer,48),HydroParam%sDRhoW(iLayer,49),HydroParam%sDRhoW(iLayer,50),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,60),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,61),HydroParam%sDRhoW(iLayer,62),HydroParam%sDRhoW(iLayer,63),HydroParam%sDRhoW(iLayer,64),HydroParam%sDRhoW(iLayer,65),HydroParam%sDRhoW(iLayer,66),HydroParam%sDRhoW(iLayer,67),HydroParam%sDRhoW(iLayer,68),HydroParam%sDRhoW(iLayer,69),HydroParam%sDRhoW(iLayer,70),HydroParam%sDRhoW(iLayer,71),HydroParam%sDRhoW(iLayer,72),HydroParam%sDRhoW(iLayer,73),HydroParam%sDRhoW(iLayer,74),HydroParam%sDRhoW(iLayer,75),HydroParam%sDRhoW(iLayer,76),HydroParam%sDRhoW(iLayer,77),HydroParam%sDRhoW(iLayer,78),HydroParam%sDRhoW(iLayer,79),HydroParam%sDRhoW(iLayer,80)
            !        EndIf
            !    EndDo
            !    
                !Gravando u!CAYO
                !!BENCH 01:
                icell = 3!41
                face = 9 !123
                Basename = 'U'
                Write(FileName,'(i10)') icell
                Open(96,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
                iLayer = 1    
                Write(96,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%u(iLayer,337), HydroParam%H(337), HydroParam%u(iLayer,457),HydroParam%H(457),HydroParam%u(iLayer,757),HydroParam%H(757)
               ! Write(96,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%Vol(icell-1),HydroParam%Vol(icell),HydroParam%Vol(icell+1), HydroParam%u(iLayer,face), HydroParam%H(face), HydroParam%u(iLayer,face),HydroParam%H(face),HydroParam%u(iLayer,face),HydroParam%H(face)
                !iLayer = 2
                !Write(96,'(I10,A10,100F30.20)') simParam%it,'Sup', HydroParam%u(iLayer,face), HydroParam%H(face),simParam%dt, HydroParam%u(iLayer,122), HydroParam%H(122),  HydroParam%ub(iLayer,2,icell), HydroParam%eta(icell)!HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !
            
                !
                !!!BENCH02:
                !icell = 41!2!41
                !face = 123!6!123!122
                !Basename = 'U'
                !Write(FileName,'(i10)') icell
                !Open(96,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
                !
                !iLayer = 1
                !
                ! up  = MeshParam%Right(MeshParam%Edge(1,icell))
                ! upup = MeshParam%Right(MeshParam%Edge(1,up))
                !
                !!!Write(96,'(I10,A10,100F30.20)') simParam%it,'Sup', HydroParam%u(iLayer,face), HydroParam%H(face),simParam%dt, HydroParam%u(iLayer,122), HydroParam%H(122),  HydroParam%ub(iLayer,2,icell), HydroParam%eta(icell)!HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !!!Write(96,'(I10,A10,100F30.20)') simParam%it,'Sup', HydroParam%Vol(icell-1),HydroParam%Vol(icell),HydroParam%Vol(icell+1), HydroParam%u(iLayer,face), HydroParam%H(face),simParam%dt, HydroParam%u(iLayer,122), HydroParam%H(122),  HydroParam%ub(iLayer,2,icell), HydroParam%eta(icell)!HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !Write(96,'(I10,A10,100F30.20)') simParam%it,'Sup', HydroParam%Vol(icell-1),HydroParam%Vol(icell),HydroParam%Vol(icell+1),HydroParam%Vol(up), HydroParam%Vol(upup), HydroParam%u(iLayer,face), HydroParam%H(face),simParam%dt, HydroParam%ub(iLayer,2,icell), HydroParam%eta(icell),HydroParam%eta(up),HydroParam%eta(upup), iNewton, innerNewton !HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !!Write(96,'(I10,A10,100F30.20)') simParam%it,'Sup', HydroParam%Vol(icell-1),HydroParam%Vol(icell),HydroParam%Vol(icell+1),HydroParam%Vol(up), HydroParam%u(iLayer,face), HydroParam%H(face),simParam%dt, HydroParam%ub(iLayer,2,icell), HydroParam%eta(icell),HydroParam%eta(up)!HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !!Write(96,'(I10,A10,100F30.20)') simParam%it,'Sup', HydroParam%Vol(icell), HydroParam%u(iLayer,face), HydroParam%H(face),simParam%dt, HydroParam%ub(iLayer,2,icell), HydroParam%eta(icell)!HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !!
                !!!
                
                
                !
                !!!BENCH02 TESTE:
                !icell = 2!41
                !face = 6!123!122
                !Basename = 'U'
                !Write(FileName,'(i10)') icell
                !Open(96,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
                !
                !iLayer = 1
                !
                ! up  = MeshParam%Right(MeshParam%Edge(1,icell))
                !Write(96,'(I10,A10,100F30.20)') simParam%it,'Sup', HydroParam%Vol(icell-1),HydroParam%Vol(icell),HydroParam%Vol(icell+1),HydroParam%Vol(up), HydroParam%u(iLayer,face), HydroParam%H(face),simParam%dt, HydroParam%ub(iLayer,2,icell), HydroParam%eta(icell),HydroParam%eta(up)!HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !
                !
                !
                
                !Do iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell) !HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
                !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
                !        Write(96,'(I10,A10,100F30.20)') simParam%it,'Sup', HydroParam%u(iLayer,face), HydroParam%H(face),simParam%dt, HydroParam%u(iLayer,122), HydroParam%H(122),  HydroParam%ub(iLayer,2,icell), HydroParam%eta(icell)!HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !        !Write(96,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%u(iLayer,face), HydroParam%H(face), HydroParam%u(iLayer,122), HydroParam%H(122),  HydroParam%ub(iLayer,2,icell), HydroParam%eta(icell)!HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !    Else
                !        Write(96,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%u(iLayer,face), HydroParam%H(face),simParam%dt, HydroParam%u(iLayer,122), HydroParam%H(122), HydroParam%ub(iLayer,2,icell), HydroParam%eta(icell)!HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !        
                !        !Write(96,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%u(iLayer,face), HydroParam%H(face), HydroParam%u(iLayer,122), HydroParam%H(122), HydroParam%ub(iLayer,2,icell), HydroParam%eta(icell)!HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !    EndIf
                !EndDo
            !    
            !    !Gravando V
            !    icell = 1
            !    Basename = 'V'
            !    Write(FileName,'(i10)') icell
            !    Open(100,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !    Do iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !        If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !            Write(100,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%ub(iLayer,2,348:360)!HydroParam%ub(iLayer,2,1),HydroParam%ub(iLayer,2,2),HydroParam%ub(iLayer,2,3),HydroParam%ub(iLayer,2,4),HydroParam%ub(iLayer,2,5),HydroParam%ub(iLayer,2,6),HydroParam%ub(iLayer,2,7),HydroParam%ub(iLayer,2,8),HydroParam%ub(iLayer,2,9),HydroParam%ub(iLayer,2,10),HydroParam%ub(iLayer,2,11),HydroParam%ub(iLayer,2,12),HydroParam%ub(iLayer,2,13),HydroParam%ub(iLayer,2,14),HydroParam%ub(iLayer,2,15),HydroParam%ub(iLayer,2,16),HydroParam%ub(iLayer,2,17),HydroParam%ub(iLayer,2,18),HydroParam%ub(iLayer,2,19),HydroParam%ub(iLayer,2,20)
            !        Else
            !            Write(100,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%ub(iLayer,2,348:360)!HydroParam%ub(iLayer,2,1),HydroParam%ub(iLayer,2,2),HydroParam%ub(iLayer,2,3),HydroParam%ub(iLayer,2,4),HydroParam%ub(iLayer,2,5),HydroParam%ub(iLayer,2,6),HydroParam%ub(iLayer,2,7),HydroParam%ub(iLayer,2,8),HydroParam%ub(iLayer,2,9),HydroParam%ub(iLayer,2,10),HydroParam%ub(iLayer,2,11),HydroParam%ub(iLayer,2,12),HydroParam%ub(iLayer,2,13),HydroParam%ub(iLayer,2,14),HydroParam%ub(iLayer,2,15),HydroParam%ub(iLayer,2,16),HydroParam%ub(iLayer,2,17),HydroParam%ub(iLayer,2,18),HydroParam%ub(iLayer,2,19),HydroParam%ub(iLayer,2,20)
            !        EndIf
            !    EndDo
            !    
            !    !Gravando W
            !    icell = 5   !Gravando velocidade vertical
            !    Basename = 'W'
            !    Write(FileName,'(i10)') icell
            !    Open(95,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !    icell = 5
            !    Do iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !        If (iLayer==HydroParam%ElCapitalM(icell)+1) Then
            !            Write(95,'(I10,A10,100F30.20)') simParam%it,'Sup', HydroParam%ub(iLayer,3,412), HydroParam%ub(iLayer,3,348:360)!HydroParam%ub(iLayer,3,1),HydroParam%ub(iLayer,3,2),HydroParam%ub(iLayer,3,3),HydroParam%ub(iLayer,3,4),HydroParam%ub(iLayer,3,5),HydroParam%ub(iLayer,3,6),HydroParam%ub(iLayer,3,7),HydroParam%ub(iLayer,3,8),HydroParam%ub(iLayer,3,9),HydroParam%ub(iLayer,3,10),HydroParam%ub(iLayer,3,11),HydroParam%ub(iLayer,3,12),HydroParam%ub(iLayer,3,13),HydroParam%ub(iLayer,3,14),HydroParam%ub(iLayer,3,15),HydroParam%ub(iLayer,3,16),HydroParam%ub(iLayer,3,17),HydroParam%ub(iLayer,3,18),HydroParam%ub(iLayer,3,19),HydroParam%ub(iLayer,3,20)
            !        Else
            !            Write(95,'(2I10,100F30.20)') simParam%it,iLayer, HydroParam%ub(iLayer,3,412), HydroParam%ub(iLayer,3,348:360)!HydroParam%ub(iLayer,3,1),HydroParam%ub(iLayer,3,2),HydroParam%ub(iLayer,3,3),HydroParam%ub(iLayer,3,4),HydroParam%ub(iLayer,3,5),HydroParam%ub(iLayer,3,6),HydroParam%ub(iLayer,3,7),HydroParam%ub(iLayer,3,8),HydroParam%ub(iLayer,3,9),HydroParam%ub(iLayer,3,10),HydroParam%ub(iLayer,3,11),HydroParam%ub(iLayer,3,12),HydroParam%ub(iLayer,3,13),HydroParam%ub(iLayer,3,14),HydroParam%ub(iLayer,3,15),HydroParam%ub(iLayer,3,16),HydroParam%ub(iLayer,3,17),HydroParam%ub(iLayer,3,18),HydroParam%ub(iLayer,3,19),HydroParam%ub(iLayer,3,20)
            !        EndIf
            !    EndDo 
            !!    
            !!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!    !Gravando Fu
            !!    icell = 16
            !!    Basename = 'FU'
            !!    Write(FileName,'(i10)') icell
            !!    Open(201,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !!    Do iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !!        If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !!            Write(201,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%Fub(iLayer,1,1),HydroParam%Fub(iLayer,1,2),HydroParam%Fub(iLayer,1,3),HydroParam%Fub(iLayer,1,4),HydroParam%Fub(iLayer,1,5),HydroParam%Fub(iLayer,1,6),HydroParam%Fub(iLayer,1,7),HydroParam%Fub(iLayer,1,8),HydroParam%Fub(iLayer,1,9),HydroParam%Fub(iLayer,1,10),HydroParam%Fub(iLayer,1,11),HydroParam%Fub(iLayer,1,12),HydroParam%Fub(iLayer,1,13),HydroParam%Fub(iLayer,1,14),HydroParam%Fub(iLayer,1,15),HydroParam%Fub(iLayer,1,16),HydroParam%Fub(iLayer,1,17),HydroParam%Fub(iLayer,1,18),HydroParam%Fub(iLayer,1,19),HydroParam%Fub(iLayer,1,20)
            !!        Else
            !!            Write(201,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%Fub(iLayer,1,1),HydroParam%Fub(iLayer,1,2),HydroParam%Fub(iLayer,1,3),HydroParam%Fub(iLayer,1,4),HydroParam%Fub(iLayer,1,5),HydroParam%Fub(iLayer,1,6),HydroParam%Fub(iLayer,1,7),HydroParam%Fub(iLayer,1,8),HydroParam%Fub(iLayer,1,9),HydroParam%Fub(iLayer,1,10),HydroParam%Fub(iLayer,1,11),HydroParam%Fub(iLayer,1,12),HydroParam%Fub(iLayer,1,13),HydroParam%Fub(iLayer,1,14),HydroParam%Fub(iLayer,1,15),HydroParam%Fub(iLayer,1,16),HydroParam%Fub(iLayer,1,17),HydroParam%Fub(iLayer,1,18),HydroParam%Fub(iLayer,1,19),HydroParam%Fub(iLayer,1,20)
            !!        EndIf
            !!    EndDo
            !!    
            !!    !Gravando FV
            !!    icell = 16
            !!    Basename = 'FV'
            !!    Write(FileName,'(i10)') icell
            !!    Open(202,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !!    Do iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !!        If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !!            Write(202,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%Fub(iLayer,2,1),HydroParam%Fub(iLayer,2,2),HydroParam%Fub(iLayer,2,3),HydroParam%Fub(iLayer,2,4),HydroParam%Fub(iLayer,2,5),HydroParam%Fub(iLayer,2,6),HydroParam%Fub(iLayer,2,7),HydroParam%Fub(iLayer,2,8),HydroParam%Fub(iLayer,2,9),HydroParam%Fub(iLayer,2,10),HydroParam%Fub(iLayer,2,11),HydroParam%Fub(iLayer,2,12),HydroParam%Fub(iLayer,2,13),HydroParam%Fub(iLayer,2,14),HydroParam%Fub(iLayer,2,15),HydroParam%Fub(iLayer,2,16),HydroParam%Fub(iLayer,2,17),HydroParam%Fub(iLayer,2,18),HydroParam%Fub(iLayer,2,19),HydroParam%Fub(iLayer,2,20)
            !!        Else
            !!            Write(202,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%Fub(iLayer,2,1),HydroParam%Fub(iLayer,2,2),HydroParam%Fub(iLayer,2,3),HydroParam%Fub(iLayer,2,4),HydroParam%Fub(iLayer,2,5),HydroParam%Fub(iLayer,2,6),HydroParam%Fub(iLayer,2,7),HydroParam%Fub(iLayer,2,8),HydroParam%Fub(iLayer,2,9),HydroParam%Fub(iLayer,2,10),HydroParam%Fub(iLayer,2,11),HydroParam%Fub(iLayer,2,12),HydroParam%Fub(iLayer,2,13),HydroParam%Fub(iLayer,2,14),HydroParam%Fub(iLayer,2,15),HydroParam%Fub(iLayer,2,16),HydroParam%Fub(iLayer,2,17),HydroParam%Fub(iLayer,2,18),HydroParam%Fub(iLayer,2,19),HydroParam%Fub(iLayer,2,20)
            !!        EndIf
            !!    EndDo
            !!    
            !!    !Gravando FW
            !!    icell = 5   !Gravando velocidade vertical
            !!    Basename = 'FW'
            !!    Write(FileName,'(i10)') icell
            !!    Open(203,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !!    icell = 5
            !!    Do iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !!        If (iLayer==HydroParam%ElCapitalM(icell)+1) Then
            !!            Write(203,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%Fub(iLayer,3,1),HydroParam%Fub(iLayer,3,2),HydroParam%Fub(iLayer,3,3),HydroParam%Fub(iLayer,3,4),HydroParam%Fub(iLayer,3,5),HydroParam%Fub(iLayer,3,6),HydroParam%Fub(iLayer,3,7),HydroParam%Fub(iLayer,3,8),HydroParam%Fub(iLayer,3,9),HydroParam%Fub(iLayer,3,10),HydroParam%Fub(iLayer,3,11),HydroParam%Fub(iLayer,3,12),HydroParam%Fub(iLayer,3,13),HydroParam%Fub(iLayer,3,14),HydroParam%Fub(iLayer,3,15),HydroParam%Fub(iLayer,3,16),HydroParam%Fub(iLayer,3,17),HydroParam%Fub(iLayer,3,18),HydroParam%Fub(iLayer,3,19),HydroParam%Fub(iLayer,3,20)
            !!        Else
            !!            Write(203,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%Fub(iLayer,3,1),HydroParam%Fub(iLayer,3,2),HydroParam%Fub(iLayer,3,3),HydroParam%Fub(iLayer,3,4),HydroParam%Fub(iLayer,3,5),HydroParam%Fub(iLayer,3,6),HydroParam%Fub(iLayer,3,7),HydroParam%Fub(iLayer,3,8),HydroParam%Fub(iLayer,3,9),HydroParam%Fub(iLayer,3,10),HydroParam%Fub(iLayer,3,11),HydroParam%Fub(iLayer,3,12),HydroParam%Fub(iLayer,3,13),HydroParam%Fub(iLayer,3,14),HydroParam%Fub(iLayer,3,15),HydroParam%Fub(iLayer,3,16),HydroParam%Fub(iLayer,3,17),HydroParam%Fub(iLayer,3,18),HydroParam%Fub(iLayer,3,19),HydroParam%Fub(iLayer,3,20)
            !!        EndIf
            !!    EndDo 
            !!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
            !else
            !    Write(94,'(1I15,100F30.20)') simParam%it, HydroParam%eta(348:360)!HydroParam%eta(540),HydroParam%eta(580),HydroParam%eta(628),HydroParam%eta(692),HydroParam%eta(760),HydroParam%eta(840)!HydroParam%eta(348:360)
            !    
            !    Write(104,'(1I15,100F30.20)') simParam%it,HydroParam%SumVer,HydroParam%SumVerAcum
            !!    
            !   ! Write(666,'(1I15,100F30.20)') simParam%it, simParam%start(simParam%it),simParam%finish(simParam%it), (simParam%finish(simParam%it)-simParam%start(simParam%it))
            !!    
            !!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!    !
            !    icell = 1
            !    Do iLayer=1,12!iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !        If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !            Write(665,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%sDRhoW(iLayer,348:360)!,HydroParam%sDRhoW(iLayer,2),HydroParam%sDRhoW(iLayer,3),HydroParam%sDRhoW(iLayer,4),HydroParam%sDRhoW(iLayer,5),HydroParam%sDRhoW(iLayer,6),HydroParam%sDRhoW(iLayer,7),HydroParam%sDRhoW(iLayer,8),HydroParam%sDRhoW(iLayer,9),HydroParam%sDRhoW(iLayer,10),HydroParam%sDRhoW(iLayer,11),HydroParam%sDRhoW(iLayer,12),HydroParam%sDRhoW(iLayer,13),HydroParam%sDRhoW(iLayer,14),HydroParam%sDRhoW(iLayer,15),HydroParam%sDRhoW(iLayer,16),HydroParam%sDRhoW(iLayer,17),HydroParam%sDRhoW(iLayer,18),HydroParam%sDRhoW(iLayer,19),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,31),HydroParam%sDRhoW(iLayer,32),HydroParam%sDRhoW(iLayer,33),HydroParam%sDRhoW(iLayer,34),HydroParam%sDRhoW(iLayer,35),HydroParam%sDRhoW(iLayer,36),HydroParam%sDRhoW(iLayer,37),HydroParam%sDRhoW(iLayer,38),HydroParam%sDRhoW(iLayer,39),HydroParam%sDRhoW(iLayer,40),HydroParam%sDRhoW(iLayer,41),HydroParam%sDRhoW(iLayer,42),HydroParam%sDRhoW(iLayer,43),HydroParam%sDRhoW(iLayer,44),HydroParam%sDRhoW(iLayer,45),HydroParam%sDRhoW(iLayer,46),HydroParam%sDRhoW(iLayer,47),HydroParam%sDRhoW(iLayer,48),HydroParam%sDRhoW(iLayer,49),HydroParam%sDRhoW(iLayer,50),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,60),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,61),HydroParam%sDRhoW(iLayer,62),HydroParam%sDRhoW(iLayer,63),HydroParam%sDRhoW(iLayer,64),HydroParam%sDRhoW(iLayer,65),HydroParam%sDRhoW(iLayer,66),HydroParam%sDRhoW(iLayer,67),HydroParam%sDRhoW(iLayer,68),HydroParam%sDRhoW(iLayer,69),HydroParam%sDRhoW(iLayer,70),HydroParam%sDRhoW(iLayer,71),HydroParam%sDRhoW(iLayer,72),HydroParam%sDRhoW(iLayer,73),HydroParam%sDRhoW(iLayer,74),HydroParam%sDRhoW(iLayer,75),HydroParam%sDRhoW(iLayer,76),HydroParam%sDRhoW(iLayer,77),HydroParam%sDRhoW(iLayer,78),HydroParam%sDRhoW(iLayer,79),HydroParam%sDRhoW(iLayer,80)
            !        Else
            !            Write(665,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%sDRhoW(iLayer,348:360)!,HydroParam%sDRhoW(iLayer,2),HydroParam%sDRhoW(iLayer,3),HydroParam%sDRhoW(iLayer,4),HydroParam%sDRhoW(iLayer,5),HydroParam%sDRhoW(iLayer,6),HydroParam%sDRhoW(iLayer,7),HydroParam%sDRhoW(iLayer,8),HydroParam%sDRhoW(iLayer,9),HydroParam%sDRhoW(iLayer,10),HydroParam%sDRhoW(iLayer,11),HydroParam%sDRhoW(iLayer,12),HydroParam%sDRhoW(iLayer,13),HydroParam%sDRhoW(iLayer,14),HydroParam%sDRhoW(iLayer,15),HydroParam%sDRhoW(iLayer,16),HydroParam%sDRhoW(iLayer,17),HydroParam%sDRhoW(iLayer,18),HydroParam%sDRhoW(iLayer,19),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,31),HydroParam%sDRhoW(iLayer,32),HydroParam%sDRhoW(iLayer,33),HydroParam%sDRhoW(iLayer,34),HydroParam%sDRhoW(iLayer,35),HydroParam%sDRhoW(iLayer,36),HydroParam%sDRhoW(iLayer,37),HydroParam%sDRhoW(iLayer,38),HydroParam%sDRhoW(iLayer,39),HydroParam%sDRhoW(iLayer,40),HydroParam%sDRhoW(iLayer,41),HydroParam%sDRhoW(iLayer,42),HydroParam%sDRhoW(iLayer,43),HydroParam%sDRhoW(iLayer,44),HydroParam%sDRhoW(iLayer,45),HydroParam%sDRhoW(iLayer,46),HydroParam%sDRhoW(iLayer,47),HydroParam%sDRhoW(iLayer,48),HydroParam%sDRhoW(iLayer,49),HydroParam%sDRhoW(iLayer,50),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,60),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,61),HydroParam%sDRhoW(iLayer,62),HydroParam%sDRhoW(iLayer,63),HydroParam%sDRhoW(iLayer,64),HydroParam%sDRhoW(iLayer,65),HydroParam%sDRhoW(iLayer,66),HydroParam%sDRhoW(iLayer,67),HydroParam%sDRhoW(iLayer,68),HydroParam%sDRhoW(iLayer,69),HydroParam%sDRhoW(iLayer,70),HydroParam%sDRhoW(iLayer,71),HydroParam%sDRhoW(iLayer,72),HydroParam%sDRhoW(iLayer,73),HydroParam%sDRhoW(iLayer,74),HydroParam%sDRhoW(iLayer,75),HydroParam%sDRhoW(iLayer,76),HydroParam%sDRhoW(iLayer,77),HydroParam%sDRhoW(iLayer,78),HydroParam%sDRhoW(iLayer,79),HydroParam%sDRhoW(iLayer,80)
            !        EndIf
            !    EndDo
            !    
            !    icell = 1
            !    Do iLayer=1,12
            !        If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !            Write(664,'(I10,A10,100F30.20)') simParam%it,'Sup',LimnoParam%sDTempW(iLayer,348:360)!,HydroParam%sDRhoW(iLayer,2),HydroParam%sDRhoW(iLayer,3),HydroParam%sDRhoW(iLayer,4),HydroParam%sDRhoW(iLayer,5),HydroParam%sDRhoW(iLayer,6),HydroParam%sDRhoW(iLayer,7),HydroParam%sDRhoW(iLayer,8),HydroParam%sDRhoW(iLayer,9),HydroParam%sDRhoW(iLayer,10),HydroParam%sDRhoW(iLayer,11),HydroParam%sDRhoW(iLayer,12),HydroParam%sDRhoW(iLayer,13),HydroParam%sDRhoW(iLayer,14),HydroParam%sDRhoW(iLayer,15),HydroParam%sDRhoW(iLayer,16),HydroParam%sDRhoW(iLayer,17),HydroParam%sDRhoW(iLayer,18),HydroParam%sDRhoW(iLayer,19),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,31),HydroParam%sDRhoW(iLayer,32),HydroParam%sDRhoW(iLayer,33),HydroParam%sDRhoW(iLayer,34),HydroParam%sDRhoW(iLayer,35),HydroParam%sDRhoW(iLayer,36),HydroParam%sDRhoW(iLayer,37),HydroParam%sDRhoW(iLayer,38),HydroParam%sDRhoW(iLayer,39),HydroParam%sDRhoW(iLayer,40),HydroParam%sDRhoW(iLayer,41),HydroParam%sDRhoW(iLayer,42),HydroParam%sDRhoW(iLayer,43),HydroParam%sDRhoW(iLayer,44),HydroParam%sDRhoW(iLayer,45),HydroParam%sDRhoW(iLayer,46),HydroParam%sDRhoW(iLayer,47),HydroParam%sDRhoW(iLayer,48),HydroParam%sDRhoW(iLayer,49),HydroParam%sDRhoW(iLayer,50),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,60),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,61),HydroParam%sDRhoW(iLayer,62),HydroParam%sDRhoW(iLayer,63),HydroParam%sDRhoW(iLayer,64),HydroParam%sDRhoW(iLayer,65),HydroParam%sDRhoW(iLayer,66),HydroParam%sDRhoW(iLayer,67),HydroParam%sDRhoW(iLayer,68),HydroParam%sDRhoW(iLayer,69),HydroParam%sDRhoW(iLayer,70),HydroParam%sDRhoW(iLayer,71),HydroParam%sDRhoW(iLayer,72),HydroParam%sDRhoW(iLayer,73),HydroParam%sDRhoW(iLayer,74),HydroParam%sDRhoW(iLayer,75),HydroParam%sDRhoW(iLayer,76),HydroParam%sDRhoW(iLayer,77),HydroParam%sDRhoW(iLayer,78),HydroParam%sDRhoW(iLayer,79),HydroParam%sDRhoW(iLayer,80)
            !        Else
            !            Write(664,'(2I10,100F30.20)') simParam%it,iLayer,LimnoParam%sDTempW(iLayer,348:360)!,HydroParam%sDRhoW(iLayer,2),HydroParam%sDRhoW(iLayer,3),HydroParam%sDRhoW(iLayer,4),HydroParam%sDRhoW(iLayer,5),HydroParam%sDRhoW(iLayer,6),HydroParam%sDRhoW(iLayer,7),HydroParam%sDRhoW(iLayer,8),HydroParam%sDRhoW(iLayer,9),HydroParam%sDRhoW(iLayer,10),HydroParam%sDRhoW(iLayer,11),HydroParam%sDRhoW(iLayer,12),HydroParam%sDRhoW(iLayer,13),HydroParam%sDRhoW(iLayer,14),HydroParam%sDRhoW(iLayer,15),HydroParam%sDRhoW(iLayer,16),HydroParam%sDRhoW(iLayer,17),HydroParam%sDRhoW(iLayer,18),HydroParam%sDRhoW(iLayer,19),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,31),HydroParam%sDRhoW(iLayer,32),HydroParam%sDRhoW(iLayer,33),HydroParam%sDRhoW(iLayer,34),HydroParam%sDRhoW(iLayer,35),HydroParam%sDRhoW(iLayer,36),HydroParam%sDRhoW(iLayer,37),HydroParam%sDRhoW(iLayer,38),HydroParam%sDRhoW(iLayer,39),HydroParam%sDRhoW(iLayer,40),HydroParam%sDRhoW(iLayer,41),HydroParam%sDRhoW(iLayer,42),HydroParam%sDRhoW(iLayer,43),HydroParam%sDRhoW(iLayer,44),HydroParam%sDRhoW(iLayer,45),HydroParam%sDRhoW(iLayer,46),HydroParam%sDRhoW(iLayer,47),HydroParam%sDRhoW(iLayer,48),HydroParam%sDRhoW(iLayer,49),HydroParam%sDRhoW(iLayer,50),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,60),HydroParam%sDRhoW(iLayer,51),HydroParam%sDRhoW(iLayer,52),HydroParam%sDRhoW(iLayer,53),HydroParam%sDRhoW(iLayer,54),HydroParam%sDRhoW(iLayer,55),HydroParam%sDRhoW(iLayer,56),HydroParam%sDRhoW(iLayer,57),HydroParam%sDRhoW(iLayer,58),HydroParam%sDRhoW(iLayer,59),HydroParam%sDRhoW(iLayer,20),HydroParam%sDRhoW(iLayer,61),HydroParam%sDRhoW(iLayer,62),HydroParam%sDRhoW(iLayer,63),HydroParam%sDRhoW(iLayer,64),HydroParam%sDRhoW(iLayer,65),HydroParam%sDRhoW(iLayer,66),HydroParam%sDRhoW(iLayer,67),HydroParam%sDRhoW(iLayer,68),HydroParam%sDRhoW(iLayer,69),HydroParam%sDRhoW(iLayer,70),HydroParam%sDRhoW(iLayer,71),HydroParam%sDRhoW(iLayer,72),HydroParam%sDRhoW(iLayer,73),HydroParam%sDRhoW(iLayer,74),HydroParam%sDRhoW(iLayer,75),HydroParam%sDRhoW(iLayer,76),HydroParam%sDRhoW(iLayer,77),HydroParam%sDRhoW(iLayer,78),HydroParam%sDRhoW(iLayer,79),HydroParam%sDRhoW(iLayer,80)
            !        EndIf
            !    EndDo
            !    
                !Gravando u !CAYO
                !icell = 1
                !Do iLayer=1,2 !HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
                !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
                !        Write(96,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%u(iLayer,337), HydroParam%H(337), HydroParam%u(iLayer,457),HydroParam%H(457),HydroParam%u(iLayer,757),HydroParam%H(757) !HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !    Else
                !        Write(96,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%u(iLayer,875), HydroParam%u(iLayer,6),  HydroParam%ub(iLayer,1,348:360)!HydroParam%u(iLayer,874),HydroParam%u(iLayer,875),HydroParam%u(iLayer,856),HydroParam%u(iLayer,876),HydroParam%Fu(iLayer,874),HydroParam%Fu(iLayer,875),HydroParam%Fu(iLayer,856),HydroParam%Fu(iLayer,876) !HydroParam%ub(iLayer,1,1),HydroParam%ub(iLayer,1,2),HydroParam%ub(iLayer,1,3),HydroParam%ub(iLayer,1,4),HydroParam%ub(iLayer,1,5),HydroParam%ub(iLayer,1,6),HydroParam%ub(iLayer,1,7),HydroParam%ub(iLayer,1,8),HydroParam%ub(iLayer,1,9),HydroParam%ub(iLayer,1,10),HydroParam%ub(iLayer,1,11),HydroParam%ub(iLayer,1,12),HydroParam%ub(iLayer,1,13),HydroParam%ub(iLayer,1,14),HydroParam%ub(iLayer,1,15),HydroParam%ub(iLayer,1,16),HydroParam%ub(iLayer,1,17),HydroParam%ub(iLayer,1,18),HydroParam%ub(iLayer,1,19),HydroParam%ub(iLayer,1,20)
                !    EndIf
                !EndDo
            !    
            !    !Gravando v
            !    icell = 1 
            !    Do iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !        If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !            Write(100,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%ub(iLayer,2,348:360)!HydroParam%ub(iLayer,2,1),HydroParam%ub(iLayer,2,2),HydroParam%ub(iLayer,2,3),HydroParam%ub(iLayer,2,4),HydroParam%ub(iLayer,2,5),HydroParam%ub(iLayer,2,6),HydroParam%ub(iLayer,2,7),HydroParam%ub(iLayer,2,8),HydroParam%ub(iLayer,2,9),HydroParam%ub(iLayer,2,10),HydroParam%ub(iLayer,2,11),HydroParam%ub(iLayer,2,12),HydroParam%ub(iLayer,2,13),HydroParam%ub(iLayer,2,14),HydroParam%ub(iLayer,2,15),HydroParam%ub(iLayer,2,16),HydroParam%ub(iLayer,2,17),HydroParam%ub(iLayer,2,18),HydroParam%ub(iLayer,2,19),HydroParam%ub(iLayer,2,20)
            !        Else
            !            Write(100,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%ub(iLayer,2,348:360)!HydroParam%ub(iLayer,2,1),HydroParam%ub(iLayer,2,2),HydroParam%ub(iLayer,2,3),HydroParam%ub(iLayer,2,4),HydroParam%ub(iLayer,2,5),HydroParam%ub(iLayer,2,6),HydroParam%ub(iLayer,2,7),HydroParam%ub(iLayer,2,8),HydroParam%ub(iLayer,2,9),HydroParam%ub(iLayer,2,10),HydroParam%ub(iLayer,2,11),HydroParam%ub(iLayer,2,12),HydroParam%ub(iLayer,2,13),HydroParam%ub(iLayer,2,14),HydroParam%ub(iLayer,2,15),HydroParam%ub(iLayer,2,16),HydroParam%ub(iLayer,2,17),HydroParam%ub(iLayer,2,18),HydroParam%ub(iLayer,2,19),HydroParam%ub(iLayer,2,20)
            !        EndIf
            !    EndDo
            !    
            !    !Gravando W
            !    icell = 5
            !    Do iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !        If (iLayer==HydroParam%ElCapitalM(icell)+1) Then
            !            Write(95,'(I10,A10,100F30.20)') simParam%it,'Sup', HydroParam%ub(iLayer,3,412), HydroParam%ub(iLayer,3,348:360)!HydroParam%ub(iLayer,3,1),HydroParam%ub(iLayer,3,2),HydroParam%ub(iLayer,3,3),HydroParam%ub(iLayer,3,4),HydroParam%ub(iLayer,3,5),HydroParam%ub(iLayer,3,6),HydroParam%ub(iLayer,3,7),HydroParam%ub(iLayer,3,8),HydroParam%ub(iLayer,3,9),HydroParam%ub(iLayer,3,10),HydroParam%ub(iLayer,3,11),HydroParam%ub(iLayer,3,12),HydroParam%ub(iLayer,3,13),HydroParam%ub(iLayer,3,14),HydroParam%ub(iLayer,3,15),HydroParam%ub(iLayer,3,16),HydroParam%ub(iLayer,3,17),HydroParam%ub(iLayer,3,18),HydroParam%ub(iLayer,3,19),HydroParam%ub(iLayer,3,20)
            !        Else
            !            Write(95,'(2I10,100F30.20)') simParam%it,iLayer, HydroParam%ub(iLayer,3,412), HydroParam%ub(iLayer,3,348:360)!HydroParam%ub(iLayer,3,1),HydroParam%ub(iLayer,3,2),HydroParam%ub(iLayer,3,3),HydroParam%ub(iLayer,3,4),HydroParam%ub(iLayer,3,5),HydroParam%ub(iLayer,3,6),HydroParam%ub(iLayer,3,7),HydroParam%ub(iLayer,3,8),HydroParam%ub(iLayer,3,9),HydroParam%ub(iLayer,3,10),HydroParam%ub(iLayer,3,11),HydroParam%ub(iLayer,3,12),HydroParam%ub(iLayer,3,13),HydroParam%ub(iLayer,3,14),HydroParam%ub(iLayer,3,15),HydroParam%ub(iLayer,3,16),HydroParam%ub(iLayer,3,17),HydroParam%ub(iLayer,3,18),HydroParam%ub(iLayer,3,19),HydroParam%ub(iLayer,3,20)
            !        EndIf
            !    EndDo 
            !!    
            !!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!    !Gravando Fu
            !!    icell = 16
            !!    Do iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !!        If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !!            Write(201,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%Fub(iLayer,1,1),HydroParam%Fub(iLayer,1,2),HydroParam%Fub(iLayer,1,3),HydroParam%Fub(iLayer,1,4),HydroParam%Fub(iLayer,1,5),HydroParam%Fub(iLayer,1,6),HydroParam%Fub(iLayer,1,7),HydroParam%Fub(iLayer,1,8),HydroParam%Fub(iLayer,1,9),HydroParam%Fub(iLayer,1,10),HydroParam%Fub(iLayer,1,11),HydroParam%Fub(iLayer,1,12),HydroParam%Fub(iLayer,1,13),HydroParam%Fub(iLayer,1,14),HydroParam%Fub(iLayer,1,15),HydroParam%Fub(iLayer,1,16),HydroParam%Fub(iLayer,1,17),HydroParam%Fub(iLayer,1,18),HydroParam%Fub(iLayer,1,19),HydroParam%Fub(iLayer,1,20)
            !!        Else
            !!            Write(201,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%Fub(iLayer,1,1),HydroParam%Fub(iLayer,1,2),HydroParam%Fub(iLayer,1,3),HydroParam%Fub(iLayer,1,4),HydroParam%Fub(iLayer,1,5),HydroParam%Fub(iLayer,1,6),HydroParam%Fub(iLayer,1,7),HydroParam%Fub(iLayer,1,8),HydroParam%Fub(iLayer,1,9),HydroParam%Fub(iLayer,1,10),HydroParam%Fub(iLayer,1,11),HydroParam%Fub(iLayer,1,12),HydroParam%Fub(iLayer,1,13),HydroParam%Fub(iLayer,1,14),HydroParam%Fub(iLayer,1,15),HydroParam%Fub(iLayer,1,16),HydroParam%Fub(iLayer,1,17),HydroParam%Fub(iLayer,1,18),HydroParam%Fub(iLayer,1,19),HydroParam%Fub(iLayer,1,20)
            !!        EndIf
            !!    EndDo
            !!    
            !!    !Gravando Fv
            !!    icell = 16
            !!    Do iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !!        If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !!            Write(202,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%Fub(iLayer,2,1),HydroParam%Fub(iLayer,2,2),HydroParam%Fub(iLayer,2,3),HydroParam%Fub(iLayer,2,4),HydroParam%Fub(iLayer,2,5),HydroParam%Fub(iLayer,2,6),HydroParam%Fub(iLayer,2,7),HydroParam%Fub(iLayer,2,8),HydroParam%Fub(iLayer,2,9),HydroParam%Fub(iLayer,2,10),HydroParam%Fub(iLayer,2,11),HydroParam%Fub(iLayer,2,12),HydroParam%Fub(iLayer,2,13),HydroParam%Fub(iLayer,2,14),HydroParam%Fub(iLayer,2,15),HydroParam%Fub(iLayer,2,16),HydroParam%Fub(iLayer,2,17),HydroParam%Fub(iLayer,2,18),HydroParam%Fub(iLayer,2,19),HydroParam%Fub(iLayer,2,20)
            !!        Else
            !!            Write(202,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%Fub(iLayer,2,1),HydroParam%Fub(iLayer,2,2),HydroParam%Fub(iLayer,2,3),HydroParam%Fub(iLayer,2,4),HydroParam%Fub(iLayer,2,5),HydroParam%Fub(iLayer,2,6),HydroParam%Fub(iLayer,2,7),HydroParam%Fub(iLayer,2,8),HydroParam%Fub(iLayer,2,9),HydroParam%Fub(iLayer,2,10),HydroParam%Fub(iLayer,2,11),HydroParam%Fub(iLayer,2,12),HydroParam%Fub(iLayer,2,13),HydroParam%Fub(iLayer,2,14),HydroParam%Fub(iLayer,2,15),HydroParam%Fub(iLayer,2,16),HydroParam%Fub(iLayer,2,17),HydroParam%Fub(iLayer,2,18),HydroParam%Fub(iLayer,2,19),HydroParam%Fub(iLayer,2,20)
            !!        EndIf
            !!    EndDo
            !!    
            !!    !Gravando FW
            !!    icell = 5
            !!    Do iLayer=HydroParam%ElSmallm(icell),HydroParam%ElCapitalM(icell)
            !!        If (iLayer==HydroParam%ElCapitalM(icell)+1) Then
            !!            Write(203,'(I10,A10,100F30.20)') simParam%it,'Sup',HydroParam%Fub(iLayer,3,1),HydroParam%Fub(iLayer,3,2),HydroParam%Fub(iLayer,3,3),HydroParam%Fub(iLayer,3,4),HydroParam%Fub(iLayer,3,5),HydroParam%Fub(iLayer,3,6),HydroParam%Fub(iLayer,3,7),HydroParam%Fub(iLayer,3,8),HydroParam%Fub(iLayer,3,9),HydroParam%Fub(iLayer,3,10),HydroParam%Fub(iLayer,3,11),HydroParam%Fub(iLayer,3,12),HydroParam%Fub(iLayer,3,13),HydroParam%Fub(iLayer,3,14),HydroParam%Fub(iLayer,3,15),HydroParam%Fub(iLayer,3,16),HydroParam%Fub(iLayer,3,17),HydroParam%Fub(iLayer,3,18),HydroParam%Fub(iLayer,3,19),HydroParam%Fub(iLayer,3,20)
            !!        Else
            !!            Write(203,'(2I10,100F30.20)') simParam%it,iLayer,HydroParam%Fub(iLayer,3,1),HydroParam%Fub(iLayer,3,2),HydroParam%Fub(iLayer,3,3),HydroParam%Fub(iLayer,3,4),HydroParam%Fub(iLayer,3,5),HydroParam%Fub(iLayer,3,6),HydroParam%Fub(iLayer,3,7),HydroParam%Fub(iLayer,3,8),HydroParam%Fub(iLayer,3,9),HydroParam%Fub(iLayer,3,10),HydroParam%Fub(iLayer,3,11),HydroParam%Fub(iLayer,3,12),HydroParam%Fub(iLayer,3,13),HydroParam%Fub(iLayer,3,14),HydroParam%Fub(iLayer,3,15),HydroParam%Fub(iLayer,3,16),HydroParam%Fub(iLayer,3,17),HydroParam%Fub(iLayer,3,18),HydroParam%Fub(iLayer,3,19),HydroParam%Fub(iLayer,3,20)
            !!        EndIf
            !!    EndDo 
            !!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            !    
            !endif
            
            !If (simParam%it == 8000) Then
            !    pause
            !endif            
            
            !!< Define cell to print
            !icell = 1781     !CONDI플O DE CONTORNO 1 
            !If (simParam%TIME == simParam%IniTime) Then 
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(99,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        Write(99,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    Else
            !        Write(99,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    EndIf
            !     If ( isNaN(LimnoParam%sDTempW(iLayer,icell)) ) Then
            !          Continue
            !     End If
            !EndDo          
            ! 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !icell = 1353     !CONDI플O DE CONTORNO 2
            !If (simParam%TIME == simParam%IniTime) Then 
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(98,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        Write(98,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    Else
            !        Write(98,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    EndIf
            !EndDo   
            !
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
            !icell = 1185   !CONDI플O DE CONTORNO 3
            !If (simParam%TIME == simParam%IniTime) Then 
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(97,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        Write(97,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    Else
            !        Write(97,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    EndIf
            !EndDo   
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !icell = 160   !CONDI플O DE CONTORNO 4
            !If (simParam%TIME == simParam%IniTime) Then 
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(96,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        Write(96,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    Else
            !        Write(96,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    EndIf
            !EndDo   
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !icell = 4    !CONDI플O DE CONTORNO 5
            !If (simParam%TIME == simParam%IniTime) Then 
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(95,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        Write(95,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    Else
            !        Write(95,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    EndIf
            !EndDo   
            !
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !icell = 28   !CONDI플O DE CONTORNO 6
            !If (simParam%TIME == simParam%IniTime) Then 
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(94,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        Write(94,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    Else
            !        Write(94,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    EndIf
            !EndDo   
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !icell = 460   !CONDI플O DE CONTORNO 7
            !If (simParam%TIME == simParam%IniTime) Then 
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(93,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        Write(93,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    Else
            !        Write(93,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    EndIf
            !EndDo  
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !icell = 3683   !CONDI플O DE CONTORNO VAZ? DE SA?A
            !If (simParam%TIME == simParam%IniTime) Then 
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(92,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        Write(92,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    Else
            !        Write(92,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    EndIf
            !EndDo   
            !
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !icell = 3777    !CONDI플O DE CONTORNO VERTEDOR
            !If (simParam%TIME == simParam%IniTime) Then 
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(91,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        Write(91,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    Else
            !        Write(91,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell)!,LimnoParam%sO2W(iLayer,icell),LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPomW(iLayer,icell,1),LimnoParam%sDPhytW(iLayer,icell,1) !,MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param) !,LimnoParam%sPO4W(iLayer,icell),LimnoParam%sPAIMW(iLayer,icell),LimnoParam%sNH4W(iLayer,icell),LimnoParam%sNO3W(iLayer,icell),LimnoParam%sDPhytW(iLayer,icell,1) 
            !    EndIf
            !EndDo   
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !icell = 3670    ! Problema PERFIL
            !If (simParam%TIME == simParam%IniTime) Then 
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(90,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!--------------------------------------------------------------------------
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        eair = 4.596*exp(17.27*MeteoParam%AirTemp(icell)/(237.3+MeteoParam%AirTemp(icell)))
            !        eair= MeteoParam%RelHum(icell)*eair/100.
            !        !Write(97,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell), LimnoParam%sDTempW(iLayer,icell) !,Sum( HydroParam%VerEddyVisc(iLayer,MeshParam%Edge(1:4,icell)) )/4., Sum( HydroParam%VerEddyDiff(iLayer,MeshParam%Edge(1:4,icell)) )/4., Sum( HydroParam%Rich(iLayer,MeshParam%Edge(1:4,icell)) )/4., LimnoParam%sDTempW(iLayer,icell), MeshParam%CREDV(iElem)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(iElem)*(exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer+1,iElem)))), 0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(iElem)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB), 0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,iElem)+273.)**4., 0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.))*(LimnoParam%sDTempW(iLayer,iElem)-MeteoParam%AirTemp(iElem)), 0.48426*(19.0+0.95*(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,iElem)/(237.3+LimnoParam%sDTempW(iLayer,iElem)))-eair)
            !        Write(90,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell), LimnoParam%sDTempW(iLayer,icell),HydroParam%u(iLayer,MeshParam%Edge(1,icell)),HydroParam%u(iLayer,MeshParam%Edge(2,icell)),HydroParam%u(iLayer,MeshParam%Edge(3,icell)),HydroParam%u(iLayer,MeshParam%Edge(4,icell)),HydroParam%w(iLayer,icell),MeshParam%EdgeBary(1,MeshParam%Edge(1,icell)),MeshParam%EdgeBary(2,MeshParam%Edge(1,icell)),MeshParam%EdgeBary(1,MeshParam%Edge(2,icell)),MeshParam%EdgeBary(2,MeshParam%Edge(2,icell)),MeshParam%EdgeBary(1,MeshParam%Edge(3,icell)),MeshParam%EdgeBary(2,MeshParam%Edge(3,icell)),MeshParam%EdgeBary(1,MeshParam%Edge(4,icell)),MeshParam%EdgeBary(2,MeshParam%Edge(4,icell))  !, MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell)))), 0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB), 0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4., 0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell)), 0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)
            !    Else
            !        Write(90,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell), LimnoParam%sDTempW(iLayer,icell),HydroParam%u(iLayer,MeshParam%Edge(1,icell)),HydroParam%u(iLayer,MeshParam%Edge(2,icell)),HydroParam%u(iLayer,MeshParam%Edge(3,icell)),HydroParam%u(iLayer,MeshParam%Edge(4,icell)),HydroParam%w(iLayer,icell) !,(MeshParam%CREDV(icell)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(icell)*(exp(-LimnoParam%aExtCoef(iLayer,icell)*(HydroParam%Ze(HydroParam%ElCapitalM(icell)+1,icell) - HydroParam%Ze(iLayer+1,icell))))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param)+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(icell)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param)-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,icell)+273.)**4./(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param)-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(LimnoParam%sDTempW(iLayer,icell)-MeteoParam%AirTemp(icell))/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param)-0.48426*(19.0+0.95*(HydroParam%Windix(icell)**2.+HydroParam%Windiy(icell)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,icell)/(237.3+LimnoParam%sDTempW(iLayer,icell)))-eair)/(HydroParam%sDRhoW(iLayer,icell)*LimnoParam%cd_param)+0.*0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer-1,icell)+273.)**4./(HydroParam%sDRhoW(iLayer-1,icell)*LimnoParam%cd_param*(Max(1.,HydroParam%DZi(iLayer-1,icell)))))
            !    EndIf
            !EndDo
            
            

            !< Print Particles Position
            !Basename = 'Cells'
            !Write(FileName,'(i10)') simParam%TIME
            !Open(99,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !
            !!< Saving position on text file
            !Do i = 1,MeshParam%nElem
            !    Write(99,'(F20.2,F20.2,I20)') MeshParam%xb(i), MeshParam%yb(i), i
            !EndDo
            !Close(99)
            !
            !
            !!< Print Particles Position
            !Basename = 'ParticleTracking'
            !Write(FileName,'(i10)') simParam%TIME
            !Open(99,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !
            !!< Saving position on text file
            !Do i=1,PartParam%nPart
            !    Write(99,'(I10,I10,F11.2,F11.2,F11.4)') PartParam%npartElem0(i),PartParam%npartElem(i),PartParam%xpart(i),PartParam%ypart(i),PartParam%RTime(i)
            !EndDo
            !Close(99)
            
            !!< Print Particles Position
            !icell = 1
            !If (simParam%TIME == simParam%IniTime) Then
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(99,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !
            !Write(99,'(I10,100F30.20)') simParam%TIME,HydroParam%eta(icell),LimnoParam%sDTempW(HydroParam%ElCapitalM(icell),icell) !HydroParam%ElSmallM(icell):HydroParam%ElSmallM(icell):
            !!EndDo
            
            !< Define cell to print
            !icell = 15257 !PAC ! 15374 !Munda?
            !iElem = 15257
            !If (simParam%TIME == simParam%IniTime) Then 
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(99,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            
            !icell =  19281 !3761 !PAC (150 m)
            !iElem = 19281 !3761
            !If (simParam%TIME == simParam%IniTime) Then
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(99,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        Write(99,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDSal(iLayer,icell) !LimnoParam%sDTempW(iLayer,icell),LimnoParam%sPO4W(iLayer,icell), LimnoParam%sPAIMW(iLayer,icell), LimnoParam%sNH4W(iLayer,icell), LimnoParam%sNO3W(iLayer,icell),LimnoParam%sO2W(iLayer,icell), LimnoParam%sDPhytW(iLayer,icell,1),LimnoParam%sDDomW(iLayer,icell,1),LimnoParam%cPBackLoad,LimnoParam%wPSorpIMW(iLayer,icell),LimnoParam%tPResusPO4(icell),LimnoParam%tPDifPO4(icell),LimnoParam%rPDDomW(iLayer,icell,1) * LimnoParam%wDMinAerW(iLayer,icell,1),LimnoParam%rPDDomW(iLayer,icell,1) * LimnoParam%wDMinAnaerW(iLayer,icell,1),SUM(LimnoParam%wPUptPhyt(iLayer,icell,:)),SUM(LimnoParam%wPExcrPhytW(iLayer,icell,:)),SUM(LimnoParam%wPUptBac(iLayer,icell,:)),SUM(LimnoParam%wPExcrZoo(iLayer,icell,:)),SUM( LimnoParam%wPExcrFiAd(iLayer,icell,:)),SUM( LimnoParam%wPExcrFiJv(iLayer,icell,:)),SUM(LimnoParam%tPUptMacW(icell,:)),LimnoParam%sPDomW(iLayer,icell,1),LimnoParam%rPDDomW(iLayer,icell,1),LimnoParam%rNDDomW(iLayer,icell,1),LimnoParam%wDMinAerW(iLayer,icell,1)
            !    Else
            !        Write(99,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDSal(iLayer,icell) !LimnoParam%sDTempW(iLayer,icell),LimnoParam%sPO4W(iLayer,icell), LimnoParam%sPAIMW(iLayer,icell), LimnoParam%sNH4W(iLayer,icell), LimnoParam%sNO3W(iLayer,icell),LimnoParam%sO2W(iLayer,icell), LimnoParam%sDPhytW(iLayer,icell,1),LimnoParam%sDDomW(iLayer,icell,1),LimnoParam%cPBackLoad,LimnoParam%wPSorpIMW(iLayer,icell),LimnoParam%tPResusPO4(icell),LimnoParam%tPDifPO4(icell),LimnoParam%rPDDomW(iLayer,icell,1) * LimnoParam%wDMinAerW(iLayer,icell,1),LimnoParam%rPDDomW(iLayer,icell,1) * LimnoParam%wDMinAnaerW(iLayer,icell,1),SUM(LimnoParam%wPUptPhyt(iLayer,icell,:)),SUM(LimnoParam%wPExcrPhytW(iLayer,icell,:)),SUM(LimnoParam%wPUptBac(iLayer,icell,:)),SUM(LimnoParam%wPExcrZoo(iLayer,icell,:)),SUM( LimnoParam%wPExcrFiAd(iLayer,icell,:)),SUM( LimnoParam%wPExcrFiJv(iLayer,icell,:)),SUM(LimnoParam%tPUptMacW(icell,:)),LimnoParam%sPDomW(iLayer,icell,1),LimnoParam%rPDDomW(iLayer,icell,1),LimnoParam%rNDDomW(iLayer,icell,1),LimnoParam%wDMinAerW(iLayer,icell,1)
            !    EndIf
            !EndDo    
            !!
            !!< Define cell to print
            !icell = 1379 !PCE (150 m)
            !If (simParam%TIME == simParam%IniTime) Then
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(98,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        Write(98,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell),LimnoParam%sPO4W(iLayer,icell), LimnoParam%sPAIMW(iLayer,icell), LimnoParam%sNH4W(iLayer,icell), LimnoParam%sNO3W(iLayer,icell),LimnoParam%sO2W(iLayer,icell), LimnoParam%sDPhytW(iLayer,icell,1),LimnoParam%sDDomW(iLayer,icell,1)
            !    Else
            !        Write(98,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell),LimnoParam%sPO4W(iLayer,icell), LimnoParam%sPAIMW(iLayer,icell), LimnoParam%sNH4W(iLayer,icell), LimnoParam%sNO3W(iLayer,icell),LimnoParam%sO2W(iLayer,icell), LimnoParam%sDPhytW(iLayer,icell,1),LimnoParam%sDDomW(iLayer,icell,1)
            !    EndIf
            !EndDo            
            !!
            !!< Define cell to print
            !icell = 2 !PCA (150 m)
            !iElem = 2
            !If (simParam%TIME == simParam%IniTime) Then
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(97,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        eair = 4.596*exp(17.27*MeteoParam%AirTemp(iElem)/(237.3+MeteoParam%AirTemp(iElem)))
            !        eair= MeteoParam%RelHum(iElem)*eair/100.
            !        Write(97,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell),LimnoParam%sPO4W(iLayer,icell), LimnoParam%sPAIMW(iLayer,icell), LimnoParam%sNH4W(iLayer,icell), LimnoParam%sNO3W(iLayer,icell),LimnoParam%sO2W(iLayer,icell), LimnoParam%sDPhytW(iLayer,icell,1),LimnoParam%sDDomW(iLayer,icell,1)
            !        !Write(97,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell),MeshParam%CREDV(iElem)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(iElem)*(exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer+1,iElem))))/(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(iElem)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,iElem)+273.)**4./(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.))*(LimnoParam%sDTempW(iLayer,iElem)-MeteoParam%AirTemp(iElem))/(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,iElem)/(237.3+LimnoParam%sDTempW(iLayer,iElem)))-eair)/(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param)
            !    Else
            !        Write(97,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell),LimnoParam%sPO4W(iLayer,icell), LimnoParam%sPAIMW(iLayer,icell), LimnoParam%sNH4W(iLayer,icell), LimnoParam%sNO3W(iLayer,icell),LimnoParam%sO2W(iLayer,icell), LimnoParam%sDPhytW(iLayer,icell,1),LimnoParam%sDDomW(iLayer,icell,1)
            !    EndIf
            !EndDo    
            
            
            
            !icell = 161 !PCA
            !iElem = 161
            !If (simParam%TIME == simParam%IniTime) Then
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(98,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        eair = 4.596*exp(17.27*MeteoParam%AirTemp(iElem)/(237.3+MeteoParam%AirTemp(iElem)))
            !        eair= MeteoParam%RelHum(iElem)*eair/100.
            !        Write(98,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell), HydroParam%u(iLayer,MeshParam%Edge(1:4,icell)),HydroParam%Fu(iLayer,MeshParam%Edge(1:4,icell)),HydroParam%w(iLayer,icell),HydroParam%sDRhoW(iLayer,icell),HydroParam%sDRhoW(iLayer,icell+1) !,LimnoParam%sDTempW(iLayer,icell),LimnoParam%sPO4W(iLayer,icell), LimnoParam%sPAIMW(iLayer,icell), LimnoParam%sNH4W(iLayer,icell), LimnoParam%sNO3W(iLayer,icell),LimnoParam%sO2W(iLayer,icell), LimnoParam%sDPhytW(iLayer,icell,1)
            !        !Write(97,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell),MeshParam%CREDV(iElem)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(iElem)*(exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer+1,iElem))))/(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(iElem)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,iElem)+273.)**4./(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.))*(LimnoParam%sDTempW(iLayer,iElem)-MeteoParam%AirTemp(iElem))/(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,iElem)/(237.3+LimnoParam%sDTempW(iLayer,iElem)))-eair)/(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param)
            !    Else
            !        Write(98,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell), HydroParam%u(iLayer,MeshParam%Edge(1:4,icell)),HydroParam%Fu(iLayer,MeshParam%Edge(1:4,icell)),HydroParam%w(iLayer,icell),HydroParam%sDRhoW(iLayer,icell),HydroParam%sDRhoW(iLayer,icell+1) !,LimnoParam%sDTempW(iLayer,icell),LimnoParam%sPO4W(iLayer,icell), LimnoParam%sPAIMW(iLayer,icell), LimnoParam%sNH4W(iLayer,icell), LimnoParam%sNO3W(iLayer,icell),LimnoParam%sO2W(iLayer,icell), LimnoParam%sDPhytW(iLayer,icell,1)
            !    EndIf
            !EndDo                
            !icell = 162 !PCA
            !iElem = 162
            !If (simParam%TIME == simParam%IniTime) Then
            !    Basename = 'Time-Series'
            !    Write(FileName,'(i10)') icell
            !    Open(99,FILE=trim(simParam%OutputPath)//'/'//trim(Basename)//trim(FileName)//'.txt',STATUS='UNKNOWN',ACTION='WRITE') !< File to save the positions
            !EndIf
            !
            !!< Saving position on text file
            !!Do i=1,PartParam%nPart
            !Do iLayer=HydroParam%ElSmallm(icell), HydroParam%ElCapitalM(icell)
            !    If (iLayer==HydroParam%ElCapitalM(icell)) Then
            !        eair = 4.596*exp(17.27*MeteoParam%AirTemp(iElem)/(237.3+MeteoParam%AirTemp(iElem)))
            !        eair= MeteoParam%RelHum(iElem)*eair/100.
            !        Write(99,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell), HydroParam%u(iLayer,MeshParam%Edge(1:4,icell)),HydroParam%Fu(iLayer,MeshParam%Edge(1:4,icell)),HydroParam%w(iLayer,icell),HydroParam%sDRhoW(iLayer,icell),HydroParam%sDRhoW(iLayer,icell+1) !,LimnoParam%sDTempW(iLayer,icell),LimnoParam%sPO4W(iLayer,icell), LimnoParam%sPAIMW(iLayer,icell), LimnoParam%sNH4W(iLayer,icell), LimnoParam%sNO3W(iLayer,icell),LimnoParam%sO2W(iLayer,icell), LimnoParam%sDPhytW(iLayer,icell,1)
            !        !Write(97,'(I10,A10,100F30.20)') simParam%TIME,'Sup',HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell),LimnoParam%sDTempW(iLayer,icell),MeshParam%CREDV(iElem)*(1.-HydroParam%ALB)*MeteoParam%SolarRad(iElem)*(exp(-LimnoParam%aExtCoef(iLayer,iElem)*(HydroParam%Ze(HydroParam%ElCapitalM(iElem)+1,iElem) - HydroParam%Ze(iLayer+1,iElem))))/(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param),+0.48426*LimnoParam%tau_param*(MeteoParam%AirTemp(iElem)+273.)**4.*(LimnoParam%a_param+0.031*sqrt(eair))*(1.-HydroParam%ALB)/(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param),-0.48426*LimnoParam%e_param*LimnoParam%tau_param*(LimnoParam%sDTempWP(iLayer,iElem)+273.)**4./(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param),-0.48426*LimnoParam%c1_param*(19.0+0.95*(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.))*(LimnoParam%sDTempW(iLayer,iElem)-MeteoParam%AirTemp(iElem))/(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param),-0.48426*(19.0+0.95*(HydroParam%Windix(iElem)**2.+HydroParam%Windiy(iElem)**2.))*(4.596*exp(17.27*LimnoParam%sDTempW(iLayer,iElem)/(237.3+LimnoParam%sDTempW(iLayer,iElem)))-eair)/(HydroParam%sDRhoW(iLayer,iElem)*LimnoParam%cd_param)
            !    Else
            !        Write(99,'(2I10,100F30.20)') simParam%TIME,iLayer,HydroParam%ub(iLayer,1,icell),HydroParam%ub(iLayer,2,icell),HydroParam%ub(iLayer,3,icell),V(HydroParam%eta(icell),HydroParam%hb(icell)),HydroParam%Ze(iLayer+1,icell), HydroParam%u(iLayer,MeshParam%Edge(1:4,icell)),HydroParam%Fu(iLayer,MeshParam%Edge(1:4,icell)),HydroParam%w(iLayer,icell),HydroParam%sDRhoW(iLayer,icell),HydroParam%sDRhoW(iLayer,icell+1) !,LimnoParam%sDTempW(iLayer,icell),LimnoParam%sPO4W(iLayer,icell), LimnoParam%sPAIMW(iLayer,icell), LimnoParam%sNH4W(iLayer,icell), LimnoParam%sNO3W(iLayer,icell),LimnoParam%sO2W(iLayer,icell), LimnoParam%sDPhytW(iLayer,icell,1)
            !    EndIf
            !EndDo    
             
            
        EndIf
    EndIf

    Return

    End Subroutine WriteOutputs
