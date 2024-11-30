Subroutine DestroyHydro(HydroParam)

    Use Hydrodynamic
    
	Implicit none
   
    type(HydrodynamicParam) :: HydroParam
    
    Deallocate(HydroParam%Z)
    Deallocate(HydroParam%Ze)
    Deallocate(HydroParam%Zb)
    Deallocate(HydroParam%DZj)
    Deallocate(HydroParam%DZjt)
    Deallocate(HydroParam%DZi)
    Deallocate(HydroParam%DZit)
    Deallocate(HydroParam%iADZ)
    Deallocate(HydroParam%iAG)
    Deallocate(HydroParam%DZiADZ)
    Deallocate(HydroParam%DZiAG)
    Deallocate(HydroParam%Smallm)
    Deallocate(HydroParam%CapitalM)
    Deallocate(HydroParam%ElSmallm)
    Deallocate(HydroParam%ElCapitalM)
    
    Deallocate(HydroParam%psij)
    Deallocate(HydroParam%rj)
    
    Deallocate(HydroParam%ubt)
    Deallocate(HydroParam%uxyt)
    Deallocate(HydroParam%uNodet)  
    Deallocate(HydroParam%ugt)
    Deallocate(HydroParam%vgt)
    Deallocate(HydroParam%wgt)
    Deallocate(HydroParam%ubVt)
    Deallocate(HydroParam%uxyLt)
    Deallocate(HydroParam%wfct)
    
    
    Deallocate(HydroParam%etaInf)
    Deallocate(HydroParam%etaInfn)
    Deallocate(HydroParam%etaplus)
    Deallocate(HydroParam%peta)
    Deallocate(HydroParam%petan)
    Deallocate(HydroParam%eta)
    Deallocate(HydroParam%etan)
    Deallocate(HydroParam%hb)
    Deallocate(HydroParam%H)
    Deallocate(HydroParam%hj)
    Deallocate(HydroParam%u)
    Deallocate(HydroParam%ut)
    Deallocate(HydroParam%uxyback)
    Deallocate(HydroParam%uArrow)
    Deallocate(HydroParam%uNode)
    Deallocate(HydroParam%uxy)
    Deallocate(HydroParam%Wu)
    Deallocate(HydroParam%ub)
    Deallocate(HydroParam%w)
    Deallocate(HydroParam%Hu)
    Deallocate(HydroParam%P)
    Deallocate(HydroParam%Aeta)
    Deallocate(HydroParam%f)
    Deallocate(HydroParam%Deta)
    Deallocate(HydroParam%rhs)
    Deallocate(HydroParam%Gu)
    Deallocate(HydroParam%Fu)
    Deallocate(HydroParam%Fv)
    Deallocate(HydroParam%Fuv)
    Deallocate(HydroParam%Fvu)
    Deallocate(HydroParam%FuxyNode)
    Deallocate(HydroParam%Fub)
    Deallocate(HydroParam%Rug)
    Deallocate(HydroParam%uIniVec)
    Deallocate(HydroParam%sDRhoW)
    
    Deallocate(HydroParam%HorViscosity)
    Deallocate(HydroParam%HorDiffusivity)
    Deallocate(HydroParam%VerEddyVisc)
    Deallocate(HydroParam%VerEddyDiff)
    
    Deallocate(HydroParam%WindVel)
    Deallocate(HydroParam%WindXY)
    Deallocate(HydroParam%Windix)
    Deallocate(HydroParam%Windiy)
    Deallocate(HydroParam%InFlowValue)
    Deallocate(HydroParam%WaterLevelValue)
    Deallocate(HydroParam%InFlowTime)
    Deallocate(HydroParam%WaterLevelTime)
    Deallocate(HydroParam%InFlownTime)
    Deallocate(HydroParam%WaterLevelnTime)
    Deallocate(HydroParam%InFlowSmallm)
    Deallocate(HydroParam%InFlowCapitalM)
    Deallocate(HydroParam%IndexWaterLevel)
    Deallocate(HydroParam%IndexInflow)
    Deallocate(HydroParam%WaterLevel)
    Deallocate(HydroParam%fetch_m)
    
    Deallocate(HydroParam%SScalar)
    Deallocate(HydroParam%SScalar2D)
    Deallocate(HydroParam%SVector)
    
    Deallocate(HydroParam%uLoadVarEst)
    Deallocate(HydroParam%dVarEst)
    Deallocate(HydroParam%Locdt)
    Deallocate(HydroParam%DZitau)
    Deallocate(HydroParam%DZistau)
    
    Deallocate(HydroParam%Smallms)
    Deallocate(HydroParam%CapitalMs)
    Deallocate(HydroParam%ElSmallms)
    Deallocate(HydroParam%ElCapitalMs)    
    
    Deallocate(HydroParam%etak)
    Deallocate(HydroParam%etam)
    
    Deallocate(HydroParam%DZsj) !CAYO
    Deallocate(HydroParam%DZsjt)!CAYO
    Deallocate(HydroParam%DZhj)!CAYO
    Deallocate(HydroParam%DZhjt)!CAYO
    Deallocate(HydroParam%DZsi) !CAYO
    Deallocate(HydroParam%DZsit)!CAYO
    Deallocate(HydroParam%DZhi)!CAYO
    Deallocate(HydroParam%DZhit)!CAYO
    Deallocate(HydroParam%Vol)!CAYO
    Deallocate(HydroParam%DZK) !Sediment Layer
    Deallocate(HydroParam%Gusub) !Sediment Layer
    Deallocate(HydroParam%PsiCrit) !Sediment Layer
    
    Deallocate(HydroParam%etak)
    Deallocate(HydroParam%Qk)
    Deallocate(HydroParam%Ci)
    Deallocate(HydroParam%d)    
    Deallocate(HydroParam%Vol2)
    Deallocate(HydroParam%Vol1)
    
    Deallocate(HydroParam%us)!CAYO
    Deallocate(HydroParam%ust)!CAYO
    Deallocate(HydroParam%ustang)!CAYO
    Deallocate(HydroParam%um)!CAYO
    Deallocate(HydroParam%umt)!CAYO
    Deallocate(HydroParam%umtang)!CAYO
    Deallocate(HydroParam%wm)
    Deallocate(HydroParam%wmt)
    Deallocate(HydroParam%uxysub)
    Deallocate(HydroParam%ubsub)
    
    Deallocate(HydroParam%utangNodes)
    
    Return
End Subroutine DestroyHydro
	
