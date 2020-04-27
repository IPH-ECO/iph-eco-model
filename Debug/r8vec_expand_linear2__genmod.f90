        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:49:47 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE R8VEC_EXPAND_LINEAR2__genmod
          INTERFACE 
            SUBROUTINE R8VEC_EXPAND_LINEAR2(N,X,BEFORE,FAT,AFTER,XFAT)
              INTEGER(KIND=4) :: AFTER
              INTEGER(KIND=4) :: FAT
              INTEGER(KIND=4) :: BEFORE
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: XFAT(BEFORE+(N-1)*(FAT+1)+1+AFTER)
            END SUBROUTINE R8VEC_EXPAND_LINEAR2
          END INTERFACE 
        END MODULE R8VEC_EXPAND_LINEAR2__genmod
