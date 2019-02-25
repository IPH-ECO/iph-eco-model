  !>Compute the norm of a Vector
  !>\author Rafael Cavalcanti
  !>\param[in] a the Vector from which the norm is computed
  !>\param[in] Dim  dimension of a vector
  !>\return res the norm of the vector
  !>\attention List of Modifications:\n
  !! - 04.03.2014: Routine Implementation (Rafael Cavalcanti)

    
Subroutine Norm(a,Dim,res)
    
    Implicit None
    Integer, intent(in) :: Dim !<dimension of the a Vector
    Real, intent(in)   :: a(Dim) !< a a Vector
    Real, intent(inout) :: res  !< norm of the Vector
    !local variables
    Real :: Sum
    Integer :: iElem
        
    Sum = 0.
    Do iElem = 1,Dim
        Sum = Sum + a(iElem)**2
    EndDo
    
    res = sqrt(Sum)
    
    Return
End Subroutine Norm