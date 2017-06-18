c        REAL A(12,12), W1(12), W2(12), A2(12,12)
C
c	OPEN(UNIT=21,STATUS='OLD',FILE='CITIES.DAT')
C
C	   Get input data.
C
c	N = 12
c	DO I = 1, N
c	   READ(21,100)(A(I,J),J=I+1,N)
c  100	   FORMAT(4X,12F4.0)
c	   DO J = I+1, N
c	      A(J,I) = A(I,J)
c	   ENDDO
c	ENDDO
C
c	DO I = 1,3
c	WRITE (6,150) (A(I,J),J=1,N)
c  150	FORMAT(' 1ST 3 LINES OF INPUT DATA:',12F5.0)
c	ENDDO
C
        REAL*4 DIST(100,100)
        REAL*4 A(100,2), w1(100), w2(100), a2(100,100)
C
        N = 100
        M = 2
        NEUR = 2*N
        IX = 137897
        DO 111 I = 1, 80
           A(I,1) = RAN(IX)+2.0
           A(I,2) = RAN(IX)+2.0
  111   CONTINUE
        DO 112 I = 81,100
           A(I,1) = RAN(IX)+3.0
           A(I,2) = RAN(IX)+3.0
  112   CONTINUE
C
C-----  Construct distances  ----------------------------------------
C
        VALM = 0.0
        DO 3 I = 1,N
           DO 2 I2 = 1,N
              DIST(I,I2) = 0.0
                 DO 1 J = 1,M
                    DIST(I,I2)=DIST(I,I2)+(A(I,J)-A(I2,J))**2
   1             CONTINUE
              DIST(I,I2)=SQRT(DIST(I,I2))
           IF (DIST(I,I2).GT.VALM) VALM = DIST(I,I2)
   2       CONTINUE
   3    CONTINUE
C
C-----  Normalize to maximum distance = unity  ----------------------
C
        DO 5 I = 1,N
           DO 4 J = 1,N
              DIST(I,J) = DIST(I,J)/VALM
   4       CONTINUE
   5    CONTINUE
c
c
	IERR = 0
	IPRINT = 2
	CALL CMDS(N,dist,IPRINT,
     X			W1,W2,A2,IERR)
	IF (IERR.NE.0) GOTO 9000
C
	GOTO 9900
 9000	WRITE (6,*) ' ABNORMAL END: IERR =', IERR
 9900	CONTINUE
	END
