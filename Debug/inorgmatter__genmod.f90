        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:50:19 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INORGMATTER__genmod
          INTERFACE 
            SUBROUTINE INORGMATTER(HYDROPARAM,MESHPARAM,LIMNOPARAM,DT,  &
     &DTDAY)
              USE LIMNOLOGYVARS
              USE MESHVARS
              USE HYDRODYNAMIC
              TYPE (HYDRODYNAMICPARAM) :: HYDROPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              REAL(KIND=4) :: DT
              REAL(KIND=4) :: DTDAY
            END SUBROUTINE INORGMATTER
          END INTERFACE 
        END MODULE INORGMATTER__genmod
