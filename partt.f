	DIMENSION A(18,16), MEMGP(18), NUMGP(18), COMP(2)
	DIMENSION GPCEN(2,16)
C
	OPEN(UNIT=21,STATUS='OLD',FILE='spectr.dat')
C
	N = 18
	M = 16
	NG = 2
	DO I = 1, N
	   READ(21,1000) (A(I,J),J=1,M)
 1000	   FORMAT(8F7.1)
	ENDDO
C
	IERR = 0
	NGP0 = 1
	NTRIES = 15
	ISEED = 37519
	DO I = 1, NTRIES
	   CALL RANDP(N,NG,MEMGP,ISEED)
C  In the following, we optimise either by assigning (using
C  minimum distances to group centres) or by exchanging.
C	   CALL MINDST(A,N,M,MEMGP,NGP0,NUMGP,GPCEN,NG,COMP,CTOT,ITER,
C     X	                             IERR)
	   CALL EXCH(A,N,M,MEMGP,NGP0,NUMGP,GPCEN,NG,COMP,CTOT,ITER,
     X	                             IERR)
	   IF (IERR.EQ.2) THEN
	      WRITE(6,*) ' A GROUP HAS LESS THAN',NGP0,' MEMBERS.',
     X		         ' REDUCE THE NUMBER OF GROUPS AND TRY AGAIN.'
	      STOP
	   ENDIF
	   IF (IERR.EQ.1) THEN
	      WRITE(6,*) ' AN INVALID GROUP ASSIGNMENT HAS BEEN DETECTED',
     X			 ' (LESS THAN 1 OR GREATER THAN THE NO. OF GPS.).',
     X		         ' IS THE NUMBER OF GROUPS CORRECTLY SPECIFIED?'
	      STOP
	   ENDIF
	   IF (IERR.NE.0) THEN
	      WRITE(6,*) ' IERR =',IERR
	      STOP
	   ENDIF
	   WRITE(6,*) ' SUM OF VARIANCES:',CTOT
	   WRITE(6,*) ' NO. OF ITERATIONS:',ITER
	   WRITE(6,*) ' ASSIGNMENTS:', (MEMGP(K),K=1,N)
	ENDDO
C
	END
