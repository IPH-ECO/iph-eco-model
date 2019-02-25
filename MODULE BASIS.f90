!***********************************************************
!* Fortran 90 version of basis_r.cpp, reduced and modified *
!* version of basis.c by J-P Moreau, Paris                 *
!* ------------------------------------------------------- *
!* Reference for original basis.c:                         *
!*                                                         *
!* "Numerical Algorithms with C, by Gisela Engeln-Muellges *
!*  and Frank Uhlig, Springer-Verlag, 1996".               *
!***********************************************************
  MODULE BASIS 
 
  PARAMETER(HALF=0.5d0,ONE=1.d0,TWO=2.d0,THREE=3.d0)

  CONTAINS

  Subroutine NormMax(vektor, n, normaxi)
  DOUBLE PRECISION normaxi, vektor(0:n-1) 
  !************************************************************************
  !* Return the maximum norm of a (0..n-1) vector  v.                     *
  !*                                                                      *
  !* global names used:                                                   *
  !* ==================                                                   *
  !*  None.                                                               *
  !************************************************************************
  DOUBLE PRECISION norm                                              ! local max
  DOUBLE PRECISION betrag                             ! magnitude of a component 
  norm=0.d0
  do i=0, n-1 
    betrag=abs(vektor(i))
    if (betrag > norm) norm = betrag
  end do   
  normaxi = norm
  return
  end subroutine NormMax 
  
  Subroutine CopyMat(ziel, quelle, n, m)
  DOUBLE PRECISION ziel(0:n-1,0:m-1),quelle(0:n-1,0:m-1)
  !************************************************************************
  !* copy the n*m elements of the matrix quelle into the matrix ziel.     *
  !*                                                                      *
  !* global names used:                                                   *
  !* ==================                                                   *
  !* none                                                                 *
  !************************************************************************
  do i=0, n-1
    do j=0, m-1 
      ziel(i,j) = quelle(i,j)
    end do
  end do
  return
  end subroutine CopyMat 

  !*****************************************
  ! Read a vector from unit u from index = 0
  !*****************************************
  Subroutine ReadVec (u, n, x) 
  INTEGER u
  DOUBLE PRECISION x(0:n-1)
  read(u,*) (x(i),i=0,n-1)
  return
  end subroutine ReadVec
  
  !*****************************************
  ! Read a vector from unit u from index = 1
  !*****************************************
  Subroutine ReadVec1 (u, n, x) 
  INTEGER u
  DOUBLE PRECISION x(1:n)
  read(u,*) (x(i),i=1,n)
  return
  end subroutine ReadVec1

  !*****************************
  ! Put all the components of a
  ! vector to a given value val
  !*****************************
  Subroutine SetVec(n,x,val)
  DOUBLE PRECISION x(0:n-1)
  DOUBLE PRECISION val
  do i = 0, n-1
    x(i) = val
  end do
  return
  end subroutine SetVec       

  !********************************
  ! Put all the components of a
  ! nxm matrix to a given value val
  !********************************     
  Subroutine SetMat(n,m,a,val)
  DOUBLE PRECISION  a(0:n-1,0:m-1)
  DOUBLE PRECISION val
  do i = 0, n-1
    do j = 0, m-1
      a(i,j) = val
    end do
  end do
  return
  end subroutine SetMat
     
  Subroutine WriteVec(u, n, x)
  INTEGER u, n
  DOUBLE PRECISION x(0:n-1)
!*====================================================================*
!*                                                                    *
!*  Put out vector x of length n to  output file (10 items per line)  *
!*  Index starts at ZERO.                                             *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*  Input parameters:                                                 *
!*  ================                                                  *
!*                                                                    *
!*      n        integer n;                                           *
!*               lenght of vector                                     *
!*      x        DOUBLE PRECISION*8 x(n)                                          *
!*               vector                                               *
!*                                                                    *
!*====================================================================*
  write(u,10) (x(i),i=0,n-1)
  return
10 format(10E13.6)
  end subroutine WriteVec
   
  Subroutine WriteVec1 (u, n, x)
  INTEGER u, n
  DOUBLE PRECISION x(1:n)
!*====================================================================*
!*                                                                    *
!*  Put out vector x of length n to  output file (10 items per line)  *
!*  Index starts at ONE.                                              *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*  Input parameters:                                                 *
!*  ================                                                  *
!*                                                                    *
!*      n        integer n;                                           *
!*               lenght of vector                                     *
!*      x        DOUBLE PRECISION*8 x(n)                                          *
!*               vector                                               *
!*                                                                    *
!*====================================================================*
  write(u,10) (x(i),i=1,n)
  return
10 format(10E13.6)
  end subroutine WriteVec1

  Subroutine ReadMat (u,n, m, a)
!*====================================================================*
!*                                                                    *
!*  Read an n x m matrix from input file.                             *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*  Input parameters :                                                *
!*  ==================                                                *
!*      u        integer u; # of unit to read                         *
!*      n        integer m; ( m > 0 )                                 *
!*               number of rows of matrix                             *
!*      m        integer n; ( n > 0 )                                 *
!*               column number of  matrix                             *
!*                                                                    *
!*   Output parameter:                                                *
!*   ================                                                 *
!*      a        DOUBLE PRECISION*8 a(n,m)                                        *
!*               matrix to read                                       *
!*                                                                    *
!*   ATTENTION : WE do not allocate storage for a here.               *
!*                                                                    *
!*====================================================================*
  INTEGER u, n, m
  DOUBLE PRECISION a(0:n-1,0:m-1)
  Read(u,*) ((a(i,j),j=0,m-1),i=0,n-1)
  return
  end subroutine ReadMat

  Subroutine ReadMat1(u,n, m, a)
!*====================================================================*
!*                                                                    *
!*  Read an n x m matrix from input file (index begins at 1)          *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*  Input parameters :                                                *
!*  ==================                                                *
!*      u        integer u; # of unit to read                         *
!*      n        integer m; ( m > 0 )                                 *
!*               number of rows of matrix                             *
!*      m        integer n; ( n > 0 )                                 *
!*               column number of  matrix                             *
!*                                                                    *
!*   Output parameter:                                                *
!*   ================                                                 *
!*      a        DOUBLE PRECISION*8 a(n,m)                                        *
!*               matrix to read                                       *
!*                                                                    *
!*   ATTENTION : WE do not allocate storage for a here.               *
!*                                                                    *
!*====================================================================*
  INTEGER u, n, m
  DOUBLE PRECISION a(n,m)
  Read(u,*) ((a(i,j),j=1,m),i=1,n)
  return
  end subroutine ReadMat1

  Subroutine WriteMat(u,n, m, a)
!*====================================================================*
!*                                                                    *
!*  Put out an m x n matrix in output file ( 10 items per line )      *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*  Input parameters:                                                 *
!*  ================                                                  *
!*      u        integer u; # of unit to read                         *
!*      n        int m; ( m > 0 )                                     *
!*               row number of matrix                                 *
!*      m        int n; ( n > 0 )                                     *
!*               column number of matrix                              *
!*      a        DOUBLE PRECISION*8 a(n,n)                                        *
!*               matrix to write                                               *
!*====================================================================*
  INTEGER u
  DOUBLE PRECISION a(0:n-1,0:m-1)
  do i=0, n-1
    write(u,10) (a(i,j),j=0,m-1)
  enddo
  return
  10 format(10E14.6)
  end subroutine WriteMat

  Subroutine WriteMat1(u,n, m, a)
!*====================================================================*
!*                                                                    *
!*  Put out an m x n matrix in output file ( 10 items per line )      *
!*  (index begins at one)                                             *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*  Input parameters:                                                 *
!*  ================                                                  *
!*      u        integer u; # of unit to read                         *
!*      n        int m; ( m > 0 )                                     *
!*               row number of matrix                                 *
!*      m        int n; ( n > 0 )                                     *
!*               column number of matrix                              *
!*      a        DOUBLE PRECISION*8 a(n,n)                                        *
!*               matrix to write                                      *
!*====================================================================*
  INTEGER u
  DOUBLE PRECISION a(n,m)
  do i=1, n
    write(u,10) (a(i,j),j=1,m)
  enddo
  return
  10 format(10E14.6)
  end subroutine WriteMat1

  Subroutine WriteHead (u, nom)
!*====================================================================*
!*                                                                    *
!*  Put out header with text in unit u           .                    *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*  Input parameters:                                                 *
!*  ================                                                  *
!*      u        integer u; # of unit to read                         *
!*      nom      character*(*) nom                                    *
!*               text of headertext (last byte is a 0)                *
!*                                                                    *
!*   Return value :                                                   *
!*   =============                                                    *
!*      None                                                          *
!*                                                                    *
!*====================================================================*
  INTEGER u
  CHARACTER*(*) nom
  Write(u,*) '----------------------------------------------------------'
  Write(u,*) nom
  Write(u,*) '----------------------------------------------------------'
  return
  end subroutine WriteHead

  Subroutine WriteEnd(u)
!*====================================================================*
!*                                                                    *
!*  Put out end of writing onto unit u.                               *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Return value :                                                   *
!*   =============                                                    *
!*      None                                                          *
!*                                                                    *
!*====================================================================*
  INTEGER u
  Write(u,*) ' '
  Write(u,*) '----------------------------------------------------------'
  return
  end subroutine WriteEnd

!***********************
! Swap two DOUBLE PRECISION values
!***********************
  Subroutine Swap(a,b)
  DOUBLE PRECISION a,b,temp
  temp=a
  a=b
  b=temp
  return
  end subroutine Swap

END MODULE

! end of file basis.f90
 