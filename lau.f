C - NOTE - SEE ACCOMPANYING DATA SET IN LAU.DAT
C Program for graph partitioning due to Lau.
C Ref: Combinatorial Heuristic Algorithms with Fortran, H.T. Lau,
C Springer-Verlag, 1986 (pp. 98-105).
C Alg: Partition set of 2n nodes, V, with an associated 2n*2n symmetric
C      cost matrix, c(i,j), into two sets P and Q = V-P such that:
C      sum(p member of P, q member of Q) c(p,q)
C      is minimized.
C Input and tested on one case by F. Murtagh (Feb. 1988).
	INTEGER IP(5),IQ(5),KP(5),KQ(5)
	REAL COST(10,10),WK1(5),WK2(5),WK3(5)
	LOGICAL IWK4(5),IWK5(5), INIT
	OPEN(UNIT=21,STATUS='OLD',FILE='LAU.DAT')  ! STATEMENT OF MINE - FM.
C
	READ(21,10) N2
  10	FORMAT(I3)
	DO 20 I = 1, N2
  20	   READ(21,30) (COST(I,J),J=1,N2)
  30	   FORMAT(10F3.0)
	N = N2/2
	ICDIM = 10
	INIT = .TRUE.
	CALL PARTIT(N2,N,COST,ICDIM,INIT,IP,IQ,KP,KQ,TCOST,
     X              WK1,WK2,WK3,IWK4,IWK5)
	WRITE(6,40) (KP(I),I=1,N)
  40	FORMAT(/' FIRST SET : ',5I3)
	WRITE(6,50) (KQ(I),I=1, N)
  50    FORMAT(/' SECOND SET :',5I3)
	WRITE(6,60) TCOST
  60    FORMAT(/' TOTAL COST =',F10.1)
	STOP
	END
C--------------------------------------------------------------------
	SUBROUTINE PARTIT (N2,N,COST,ICDIM,INIT,IP,IQ,KP,KQ,TCOST,
     X			   WK1,WK2,WK3,IWK4,IWK5)
C
	INTEGER IP(N),IQ(N),KP(N),KQ(N)
	REAL COST(ICDIM,1),WK1(N),WK2(N),WK3(N)
	LOGICAL IWK4(N),IWK5(N),INIT
C
	IF (INIT) THEN
	   DO 10 I = 1, N
	      IP(I) = I
	      IQ(I) = I+N
  10	   CONTINUE
	ENDIF
C
  20 	DO 30 I = 1, N
	   IWK4(I) = .TRUE.
	   IWK5(I) = .TRUE.
  30	CONTINUE
	TCOST = 0.
	DO 40 I = 1, N
	   DO 40 J = 1, N
  40	   TCOST = TCOST + COST(IP(I),IQ(J))
	SMALL = -2.0*TCOST
C
	DO 70 I = 1, N
	   TOT1 = 0.
	   DO 50 J = 1, N
  50	      TOT1 = TOT1 + COST(IP(I),IQ(J))
C
	TOT2 = 0.
	DO 60 K = 1, N
  60	   TOT2 = TOT2 + COST(IP(I),IP(K))
C
	WK1(I) = TOT1 - TOT2
  70  	CONTINUE
	DO 100 I = 1, N
C
	TOT1 = 0.
	DO 80 J = 1, N
  80	   TOT1 = TOT1 + COST(IQ(I),IP(J))
C
	TOT2 = 0.
	DO 90 K = 1, N
  90	   TOT2 = TOT2 + COST(IQ(I),IQ(K))
C
	WK2(I) = TOT1 - TOT2
 100	CONTINUE
	DO 140 I = 1, N
C
	   TMAX = SMALL
	   DO 120 J = 1, N
	      IF (IWK4(J)) THEN
		 DO 110 K = 1, N
		    IF (IWK5(K)) THEN
		    GAIN = WK1(J)+WK2(K)-2.*COST(IP(J),IQ(K))
		    IF (GAIN.GT.TMAX) THEN
		       TMAX = GAIN
		       IA = IP(J)
		       IB = IQ(K)
		       IND1 = J
		       IND2 = K
		    ENDIF
	       ENDIF
 110	    CONTINUE
	 ENDIF
 120     CONTINUE
C
	WK3(I) = TMAX	
	KP(I) = IA
	KQ(I) = IB
	IWK4(IND1) = .FALSE.
	IWK5(IND2) = .FALSE.
C
	DO 130 J = 1, N
	   IF (IWK4(J))
     X        WK1(J) = WK1(J)+2.0*COST(IP(J),IA)-2.0*COST(IP(J),IB)
           IF (IWK5(J))
     X        WK2(J) = WK2(J)+2.0*COST(IQ(J),IB)-2.0*COST(IQ(J),IA)
 130	CONTINUE
 140	CONTINUE
C
	TMAX = SMALL
	DO 160 I = 1,N
	   TOT1 = 0.
	   DO 150 J = 1, I
 150	      TOT1 = TOT1 + WK3(J)
	   IF (TOT1.GT.TMAX) THEN
	      TMAX = TOT1
	      K = I
	   ENDIF
 160	CONTINUE
C
 	IF (TMAX.GT.0) THEN
	   DO 170 I = 1, K
	      IP(I) = KQ(I)
	      IQ(I) = KP(I)
 170	   CONTINUE
	   K1 = K+1
	   DO 180 I = K1,N
	      IP(I) = KP(I)
	      IQ(I) = KQ(I)
 180	   CONTINUE
	   GOTO 20
        ENDIF
C
	RETURN
	END
