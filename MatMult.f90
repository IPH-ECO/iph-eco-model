SUBROUTINE MatMult(A,B,C,N,M)
    ! This routine multiplies two Real matrixes
    ! Called in routines: ATUALIZA and COEFICIENTES
    Integer, intent(in) :: N,M
    Real, intent(in)  :: A(N,N),B(N,M)
    Real, intent(out) ::C(N,M)
    Real              :: SUM
    DO I = 1,N                                                                  
        DO J = 1,M                                                                
            SUM = 0.                                                                
            DO K = 1,N                                                              
                SUM = SUM + A(I,K)*B(K,J)                                               
            ENDDO                                                                 
            C(I,J) = SUM                                                            
        ENDDO                                                                   
    ENDDO                                                                     
    RETURN
END SUBROUTINE MATMULT