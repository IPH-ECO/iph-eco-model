        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 24 15:49:47 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE R8VEC_EXPAND_LINEAR__genmod
          INTERFACE 
            SUBROUTINE R8VEC_EXPAND_LINEAR(N,X,FAT,XFAT)
              INTEGER(KIND=4) :: FAT
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: XFAT((N-1)*(FAT+1)+1)
            END SUBROUTINE R8VEC_EXPAND_LINEAR
          END INTERFACE 
        END MODULE R8VEC_EXPAND_LINEAR__genmod
