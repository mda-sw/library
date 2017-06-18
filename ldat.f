	REAL DATA(18,16) 
        REAL TOTAL(16,16), W1(16), W2(16), W3(16)
	REAL MEAN(16), MGP(2,16)
	INTEGER GP(18),nog(2), IW1(16), IW2(16)
C
	OPEN(UNIT=21,STATUS='OLD',FILE='spectr2.dat')
C
C	   Get input data.
C
	N = 18
	M = 16
	DO I = 1, N
	   READ(21,100)(DATA(I,J),J=1,M),GP(I)
  100	   FORMAT(8F7.1,/,8F7.1,/,I1)
	ENDDO
        do i = 1, 13
          gp(i)=1
        enddo
        do i = 14, 18
          gp(i)=2
        enddo
C
	IERR = 0
	NG = 2
	IPRINT = 3
	CALL LDA(N,M,DATA,GP,IPRINT,MEAN,MGP,TOTAL,W1,W2,W3,nog,
     X         IW1, IW2,IERR)
	IF (IERR.NE.0) GOTO 9000
C
	GOTO 9900
 9000	WRITE (6,*) ' ABNORMAL END: IERR =', IERR
 9900	CONTINUE
	END
