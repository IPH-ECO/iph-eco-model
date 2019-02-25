!*******************************************************************************
!-------------------------------------------------------------------------------
!  uTempFunction:

!  Original implemented by Fenjuan Hu(feh@bios.au.dk)
!-------------------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_pclake_foodweb_water
! !INTERFACE:
    Module uTempFunction
! !DESCRIPTION:

! !USES:
    Use LimnologyVars
   
   
   implicit none

!  default: all is private.
      private
      
      PUBLIC uFunTmAbio
      PUBLIC uFunTmBio
      PUBLIC uFunTmVeg

   contains
   
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
   pure Real function uFunTmAbio(uTm,cTheta)
!
!    !DESCRIPTION:
!    The temperature function of abiotic processes is defined here
!    This function is applied for mineralization both in water and sediment
!    nitrification both in water and sediment, diffusion and sedimentation
!   
!    !INPUT PARAMETERS:

      Real, intent(in) :: uTm,cTheta
      !real(rk)              :: cTmRef=20.0_rk
!    !REVISION HISTORY:
!     Original author(s): Fenjuan Hu
!
!EOP
!-----------------------------------------------------------------------
!BOC
   uFunTmAbio = cTheta**(uTm-20.0)

   end function uFunTmAbio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
   pure Real function uFunTmBio(uTm,cSigTm,cTmOpt)
!
!    !DESCRIPTION:
!    The temperature function of biotic processes is defined here
!    This funcition is applied for all groups of phytoplankton both in water and sediment
!    zooplankton in water column, fish and zoobenthos
!   
!    !INPUT PARAMETERS:
      Real, intent(in) :: uTm,cSigTm,cTmOpt
      !real(rk)              :: cTmRef=20.0_rk
!   
!    !REVISION HISTORY:
!     Original author(s): Fenjuan Hu
!
!EOP
!-----------------------------------------------------------------------
!BOC
   uFunTmBio = exp(-0.5/cSigTm**2 *((uTm-cTmOpt)**2- (20.0-cTmOpt)**2))

   end function uFunTmBio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
   pure Real function uFunTmVeg(uTm,cQ10Veg)
!
!    !DESCRIPTION:
!    The temperature function of vegetation processes is defined here
!    This function is applied for macrophytes 
!   
!    !INPUT PARAMETERS:
      Real, intent(in) :: uTm,cQ10Veg
      !real(rk)              :: cTmRef=20.0_rk
!   
!    !REVISION HISTORY:
!     Original author(s): Fenjuan Hu
!
!EOP
!-----------------------------------------------------------------------
!BOC
   uFunTmVeg = ((cQ10Veg )** (0.1 * (uTm - 20.0)))

   end function uFunTmVeg

   end module uTempFunction

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
