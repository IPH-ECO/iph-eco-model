        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:51:11 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ALLOCATEWQVARS__genmod
          INTERFACE 
            SUBROUTINE ALLOCATEWQVARS(LIMNOPARAM,MESHPARAM)
              USE MESHVARS
              USE LIMNOLOGYVARS
              TYPE (LIMNOLOGYPARAM) :: LIMNOPARAM
              TYPE (MESHGRIDPARAM) :: MESHPARAM
            END SUBROUTINE ALLOCATEWQVARS
          END INTERFACE 
        END MODULE ALLOCATEWQVARS__genmod
