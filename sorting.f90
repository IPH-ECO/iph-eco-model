! --------------------------------------------------------------------
! PROGRAM  Sorting:
!    This program can sort a set of numbers.  The method used is 
! usually referred to as "selection" method.
! --------------------------------------------------------------------

Module  Sorting


CONTAINS

! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

   INTEGER FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      Real, DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, End
      INTEGER                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)		! assume the first is the min
      Location = Start			! record its position
      DO i = Start+1, End		! start with next elements
         IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
            Minimum  = x(i)		!      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location        	! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMaximum():
!    This function returns the location of the maximum in the section
! between Start and End.
! --------------------------------------------------------------------

   INTEGER FUNCTION  FindMaximum(x, Start, End)
      IMPLICIT  NONE
      Real, DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, End
      Real                               :: Maximum
      INTEGER                            :: Location
      INTEGER                            :: i

      Maximum  = x(Start)		! assume the first is the max
      Location = Start			! record its position
      DO i = Start+1, End		! start with next elements
         IF (x(i) > Maximum) THEN	!   if x(i) greater than the max?
            Maximum  = x(i)		!      Yes, a new maximum found
            Location = i                !      record its position
         END IF
      END DO
      FindMaximum = Location        	! return the position
   END FUNCTION  FindMaximum
   
! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      Real, INTENT(INOUT) :: a, b
      Real                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Size)
      IMPLICIT  NONE
      Real, DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1			! except for the last
         Location = FindMinimum(x, i, Size)	! find min from this to last
         CALL  Swap(x(i), x(Location))	! swap this and the minimum
      END DO
   END SUBROUTINE  Sort
    
   SUBROUTINE  SortDecreasing(x, Size)
      IMPLICIT  NONE
      Real, DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1			! except for the last
         Location = FindMaximum(x, i, Size)	! find min from this to last
         CALL  Swap(x(i), x(Location))	! swap this and the minimum
      END DO
   END SUBROUTINE  SortDecreasing
   
END Module  Sorting
