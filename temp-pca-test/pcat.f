	REAL DATA(150,4)
        REAL ARRAY1(4,4), ARRAY2(4,4), VECT1(4), VECT2(4)
        real vect3(150), vect4(150)
C
	OPEN(UNIT=21,STATUS='OLD',FILE='iris.dat')
C
	N = 150
	M = 4
	DO I = 1, N
	   READ(21,*)(DATA(I,J),J=1,M)
c  100	   FORMAT(8F7.1)
 	ENDDO

C	CALL OUTMAT(N,M,DATA)
C
	METHOD = 3
	IPRINT = 3
	CALL PCA(N,M,DATA,METHOD,IPRINT,ARRAY1,VECT1,VECT2,
     X              vect3,vect4,ARRAY2,IERR)
             
	IF (IERR.NE.0) GOTO 9000
C
	GOTO 9900
 9000	WRITE (6,*) ' ABNORMAL END: IERR =', IERR
 9900	CONTINUE
	END