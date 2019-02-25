Subroutine Thomas(DIM,a,b,c,d,e) 

    Implicit None
    Integer:: i, N, DIM
    Double Precision:: gamma
    Double Precision, Dimension(DIM):: a,b,c,d,e,c_aux,d_aux
    
    N = size(d,1)      ! number of unknowns 
    If (N == 0) then
        Continue
    EndIf
    ! Part I 
    c_aux(1) = c(1)/b(1)
    d_aux(1) = d(1)/b(1)
    Do i = 2,N
        gamma = 1/( b(i) - c_aux(i-1)*a(i))
        c_aux(i)  = c(i)*gamma
        d_aux(i)  = (d(i)-a(i)*d(i-1))*gamma
    EndDo
    
    ! Part II 
    e = 0.
    e(N) = d_aux(N) 
    Do i = N-1,1,-1
        e(i) = d_aux(i) - c_aux(i)*e(i+1)
    EndDo
    
    Return
End Subroutine Thomas
