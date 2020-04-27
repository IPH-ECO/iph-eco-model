        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 27 10:43:08 2020
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
