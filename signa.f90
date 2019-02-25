  !> Compute signed area formed by vertices 1,2,3    
    Function signa(x1,x2,x3,y1,y2,y3)

    Implicit none
    Real:: x1,x2,x3,y1,y2,y3,signa

    signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
      
    
    End Function signa     
