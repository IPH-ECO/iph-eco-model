  !> Compute the Cross Product Between Two Vectors
  !>\param a a vector to be multiplied
  !>\param b a vector to be multiplied
  !>\return c the resultant vector
  !>\author  Rafael Cavalcanti
  !>\attention  List of Modifications: \n
  !> - 04.03.2014: Routine Implementation (Rafael Cavalcanti)

Subroutine Cross(a,b,c)
    
    Implicit None
    Real, intent(in) ::  a(1:3) !<input vector
    Real, intent(in) ::  b(1:3) !<input vector
    Real, intent(inout) :: c(1:3) !<output vector
    
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
    
    Return
End Subroutine Cross