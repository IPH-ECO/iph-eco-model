subroutine solve_tridiag(a,b,c,v,x,n)
      implicit none
!      a - sub-diagonal (means it is the diagonal below the main diagonal)
!      b - the main diagonal
!      c - sup-diagonal (means it is the diagonal above the main diagonal)
!      v - right part
!      x - the answer
!      n - number of equations
!      Taken from the following Wikipedia page on 25 January 2012:
!        <http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm>.
 
        integer,intent(in) :: n
        Real,dimension(n),intent(in) :: a,b,c,v
        Real,dimension(n),intent(out) :: x
        Real,dimension(n) :: bp,vp,a_aux,b_aux,c_aux,v_aux,x_aux,d_aux
        Real :: w
        integer i
        
        !Flipping matrix
        do i = 1,n
         a_aux(i) = a(n+1-i)
         b_aux(i) = b(n+1-i)
         c_aux(i) = c(n+1-i)
         v_aux(i) = v(n+1-i)
        enddo
 
        d_aux = 0   
        x = d_aux
        w = b_aux(1)
        x(1) = v_aux(1)/w
        do i=2,n
            d_aux(i-1) = c_aux(i-1)/w;
            w = b_aux(i) - a_aux(i)*d_aux(i-1);
            x(i) = ( v_aux(i) - a_aux(i)*x(i-1) )/w;
        enddo
        do i=n-1, 1, -1
           x(i) = x(i) - d_aux(i)*x(i+1);
        enddo

        !Flipping matrix
        do i = 1,n
         x_aux(i) = x(n+1-i)
        enddo
        x = x_aux
        
    end subroutine solve_tridiag
    
Subroutine TridT(DIM,Mat,D,Sol) 

    ! Tridiagonal Matrix solver for scalar equation using the double-sweep method
    
    ! Input:
    ! DIM   - Number of Equations (and Variables to Solve)
    ! Mat   - Coefficient Matrix: Upper (3), Lower (1), and Central (2) Diagonal
    ! D     - Right-Hand Side Matrix
    
    ! Output:
    ! Sol - System Solution
    
    ! List of Modifications: 
    !   -> 05.09.2017: Routine Implementation     (J. Rafael Cavalcanti)

    ! Programmer: J. Rafael Cavalcanti

    Implicit None
    Integer, Intent(In):: DIM
    Integer:: k
    Real, Dimension(DIM), Intent(In):: D
    Real, Dimension(DIM), Intent(Out):: Sol
    Real, Dimension(3,DIM), Intent(In):: Mat
    Real, Dimension(DIM):: a, b, c, e, f

    ! 1. Load diagonals of coefficient matrix into 1-D arrays and rename RHS vector
    a(1:DIM) = Mat(1,1:DIM); b(1:DIM) = Mat(2,1:DIM); c(1:DIM) = Mat(3,1:DIM)
    
    ! 2. Initialize Variables for Forward Sweep
    e(1) = -c(1)/b(1); f(1) = D(1)/b(1)
    
    ! 3. Forward Sweep
    If ( DIM == 2 ) Then
        GoTo 1
    Else
        Continue
    End If
    
    Do k = 2, DIM-1
        e(k) = -c(k)/( b(k) + a(k)*e(k-1) )
        f(k) = ( D(k) - a(k)*f(k-1) )/( b(k) + a(k)*e(k-1) )
    End Do
    
    ! 4. Compute Scalar in Bottom Layer
1   Sol(DIM) = ( D(DIM) - a(DIM)*f(DIM-1) )/( b(DIM) + a(DIM)*e(DIM-1) )
    
    ! 5. Backward Sweep
    Do k = DIM-1, 1, -1
        Sol(k) = e(k)*Sol(k+1) + f(k)
    End Do
    
    Return
End Subroutine TridT
