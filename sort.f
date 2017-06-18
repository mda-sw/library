	DIMENSION X(4)
	N=4
	X(1) = 4
	X(2) = 7
	X(3) = 2
	X(4) = 5
        write(6,*) x
	CALL TRIAGE(N,X)
        write(6,*) x
	END
      SUBROUTINE TRIAGE (N,X)
      DIMENSION X(N)
      DO 1 J=2,N                 
      JJ=J-1  
      XJ=X(J) 
      DO 2 K=1,JJ     
      IF(XJ.LT.X(K)) GOTO 3       
    2 CONTINUE     
      GOTO 1       
    3 CONTINUE
      LL=J-K      
      DO  4 L=1,LL 
      JML=J-L      
      X(JML+1)=X(JML)    
    4 CONTINUE
      X(K)=XJ 
    1 CONTINUE            
      RETURN         
      END                                                                       
