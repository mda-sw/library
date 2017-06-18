C       Driver routine for linear regression with errors in both variables.
C       Input data given in DATA statements.
C       Reference: Ripley and Thompson (see routine RIPLEY).

	REAL*4		X(5000),Y(5000),WTX(5000),WTY(5000),W(5000)
        REAL*4          WW(5000)
        COMMON          N,X,Y,WTX,WTY,W,WW
        DATA            X /0.0,0.9,1.8,2.6,3.3,4.4,5.2,6.1,6.5,7.4,
     +                     4990*0.0/
        DATA            Y /5.9,5.4,4.4,4.6,3.5,3.7,2.8,2.8,2.4,1.5, 
     +                     4990*0.0/
        DATA            WTX /1000.,1000.,500.,800.,200.,80.,60.,20.,
     +                     1.8,1.0,4990*0.0/
        DATA            WTY /1.0,1.8,4.0,8.0,20.,20.,70.,70.,100.,500.,
     +                      4990*0.0/
C
	N = 10
C
C       Weights are given; determine sigmas:
	DO I=1,N			              
	   WTX(I) = 1./WTX(I)
	   WTY(I) = 1./WTY(I)
	ENDDO	                              
C							
C
	ISTAT = 0
	CALL Ripley(ALPH,BETA,SEALPHA,SEBETA)
  90	WRITE(6,1001) N,ITER,BETA,ALPH,SEBETA,SEALPHA	
1001	FORMAT(1X,'N =',I4,' Niter =',I5,' Slope =',G16.8,		
     x  ' Intercept =',G16.8,/,' Errors: se(slope) =',G11.4,	
     X  ' se(intercept) =',G11.4)
C              
	END
C-----------------------------------------------------------------------------
C     Program to fit functional relationships Y=X=U+DELTA,
C     V=ALPHA+BETA*U with known variances for DELTA and ETA.
C     Reference: B.D. Ripley and M. Thompson, Regression techniques for
C                the detection of analytical bias, Analyst, 112, 377-383,
C                1987.
C
C     B.D. RIPLEY      1981,1986
C-----------------------------------------------------------------------------
      SUBROUTINE RIPLEY(ALPH,BETA,SEALPH,SEBETA)
      REAL X(5000), Y(5000), VARX(5000), VARY(5000), WT(5000), U(5000)
      CHARACTER MMM*3
      LOGICAL INTER
      COMMON N,X,Y,VARX,VARY,U,WT
      COMMON /CNTL/ INTER,ALPHA
C     The following is commented out because data is given in DATA statement
C     in calling routine.  This subroutine can be changed to a main routine,
C     and, if so, activating the following allows a data file name to be
C     given to the program.
C      CALL DATIN (X,Y,N,VARX,VARY)
      IF (N.LE.2 .OR. N.GT.5000) THEN
        PRINT 1010,N
        STOP
      ENDIF
C      The $ below is DEC-specific and allows answers on the prompt line
C      delete the ,$ or replace by a similar feature.
C      See also line in DATIN.
      PRINT '(1X,A,$)','Line through (0,0) '
      READ '(A3)',MMM
      INTER = MMM.NE.'YES'
      BETA = GUESS (X,Y,N)
      CALL MLE (ITER,BETA)
      CALL SE (U,N,WT,SEALPH,SEBETA)
      PRINT 1000,ALPHA,BETA,SEALPH,SEBETA
      CALL RESCHK (X,Y,N,WT,U,ALPHA,BETA)
 1000 FORMAT(/8X,'Intercept',3X,'Slope'/' value',2X,G10.4,1X,G10.4/
     + ' S.E. ',2X,G10.4,1X,G10.4/)
 1010 FORMAT (1X,I4,' observations out of range 3 to 100')
      ALPH = ALPHA
      END
C
      SUBROUTINE DATIN (X,Y,N,VX,VY)
C-----------------------------------------------------------------------------
C     read in data from file data in colums for X,Y,VARX,VARY
C     X and Y are the measurements by the two methods
C     VARX and VARY are the variances of the X and Y measurements
C      respectively
C-----------------------------------------------------------------------------
      REAL X(5000), Y(5000), VX(5000), VY(5000)
      CHARACTER NAME*20
C      non-standard feature in line below
      PRINT '(1X,A,$)','Data file name'
      READ '(A)',NAME
      OPEN (UNIT=1,FILE=NAME,STATUS='OLD',IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
        PRINT *,'File name error'
        STOP
      ENDIF
      M = 1
   10 READ(1,*,END=20) X(M),Y(M),VX(M),VY(M)
      M = M+1
      IF (M.LE.101) GO TO 10
   20 N = M-1
      RETURN
      END
C
      FUNCTION WTSUM (BETA)
C-----------------------------------------------------------------------------
C      Calculate weighted sum of squares to be minimized, given beta
C-----------------------------------------------------------------------------
      COMMON N, X(5000), Y(5000), VARX(5000), VARY(5000), U(5000), 
     .                                                   WT(5000)
      COMMON /CNTL/ INTER, ALPHA
      LOGICAL INTER
      XSUM = 0.0
      YSUM = 0.0
      WS = 0.0
      DO 10 I = 1,N
        A = 1.0/(VARY(I)+BETA*BETA*VARX(I))
        XSUM = XSUM+A*X(I)
        YSUM = YSUM+A*Y(I)
        WS = WS+A
   10   WT(I) = A
      ALPHA = 0.0
      IF (INTER) ALPHA = (YSUM-BETA*XSUM)/WS
      DO 20 I = 1,N
        B = VARY(I)*X(I)+BETA*VARX(I)*(Y(I)-ALPHA)
   20   U(I) = B*WT(I)
      A = 0.0
      DO 30 I = 1,N
        B = X(I)-U(I)
        C = Y(I)-ALPHA-BETA*U(I)
   30   A = A+B*B/VARX(I)+C*C/VARY(I)
      WTSUM = A
      RETURN
      END
C
      SUBROUTINE MLE (ITER,BETA)
C-----------------------------------------------------------------------------
C      Minimizes sum from WTSUM
C-----------------------------------------------------------------------------
      EXTERNAL WTSUM
      STEP = 0.2*ABS(BETA)
      CALL OPT (WTSUM,BETA,STEP)
      RETURN
      END
C
      FUNCTION GUESS (X,Y,N)
C-----------------------------------------------------------------------------
C      Gives guess of slope beta via least squares
C-----------------------------------------------------------------------------
      REAL X(N), Y(N)
      SX = 0.0
      SY = 0.0
      DO 10 I = 1,N
        SX = SX+X(I)
   10   SY = SY+Y(I)
      XBAR = SX/N
      YBAR = SY/N
      SXX = 0.0
      SYY = 0.0
      DO 20 I = 1,N
        SXX = SXX+(X(I)-XBAR)**2
   20   SYY = SYY+(X(I)-XBAR)*(Y(I)-YBAR)
      GUESS = SYY/SXX
      RETURN
      END
C
      SUBROUTINE SE (X,N,WT,SEALPH,SEBETA)
C-----------------------------------------------------------------------------
C      Calculate s.e.'s for alpha and beta
C-----------------------------------------------------------------------------
      REAL X(N), WT(N)
      LOGICAL INTER
      COMMON /CNTL/ INTER, ALPHA
      XSUM = 0.0
      WS = 0.0
      DO 10 I = 1,N
        A = WT(I)
        WS = WS+A
   10   XSUM = XSUM+A*X(I)
      XBAR = XSUM/WS
      A1 = 0.0
      A2 = 0.0
      DO 20 I = 1,N
        B = X(I)
        A = WT(I)
        A1 = A1+A*B*B
        B = B-XBAR
   20   A2 = A2+A*B*B
      IF (INTER) THEN
        SEALPH = SQRT(A1/(A2*WS))
        SEBETA = SQRT(1.0/A2)
      ELSE
        SEALPH = 0.0
        SEBETA = SQRT(1.0/A1)
      ENDIF
      RETURN
      END
C
      SUBROUTINE RESCHK (X,Y,N,WT,U,ALPHA,BETA)
C-----------------------------------------------------------------------------
C      Calculate scaled residuals
C-----------------------------------------------------------------------------
      REAL X(N), Y(N), WT(N), U(N), RES(5000)
      CHARACTER ANS*3
      SUM = 0.0
      DO 10 I = 1,N
        A = Y(I)-ALPHA-BETA*X(I)
        RES(I) = A*SQRT(WT(I))
   10   SUM = SUM + RES(I)**2
      PRINT 1000
      PRINT 1010,(RES(I),I=1,N)
      PRINT 1020,SUM/N
      PRINT 1030
      READ '(A)',ANS
      CALL RESPLT (N,RES,U)
      RETURN
 1000 FORMAT('   Scaled residuals'/)
 1010 FORMAT(10F6.2)
 1020 FORMAT(/1X,'Aver. squared redidual ',F6.2)
 1030 FORMAT(/1X,'Press RETURN to plot',
     + ' scaled residuals vs fitted values')
      END
C
      SUBROUTINE RESPLT (N,RES,FITTED)
C-----------------------------------------------------------------------------
C      plot scaled residuals vs fitted values
C-----------------------------------------------------------------------------
      REAL RES(N), FITTED(N)
      XU = FITTED(1)
      DO 10 I = 2,N
   10   XU = MAX(XU,FITTED(I))
      MX = INT(ALOG10(XU)+10.01)-10
      XU1 = 10.0**MX
      XU = XU1*(INT(XU/XU1+0.99))
      YM = 0.0
      DO 20 I = 1,N
   20   YM = MAX(YM,ABS(RES(I)))
      YM = INT(YM+0.99)
      CALL STPLT (0.0,XU,-YM,YM,50,20)
      DO 30 I = 1,N
   30   CALL ADDPT (FITTED(I),RES(I),'*')
      CALL ENPLT (6)
      RETURN
      END
C
      SUBROUTINE OPT (FUNC,START,STEP)
C
C       Optimization routine from
C       J.C. Nash - Computer Numerical Methods
C       for Computers, Adam Hilger,1979.
C
C       START is initial guess
C       STEP is guess at step size
C
      PARAMETER (BIG = 1.0E20)
      DATA A1,A2/1.5,-0.25/
      IFN = 0
      B = START
    1 IFN = IFN+1
      P = FUNC(B)
    2 S1 = P
      SO = -BIG
      X1 = 0.0
      BMIN = B
    3 X2 = X1 +STEP
      B = BMIN+X2
      IF (B .EQ. BMIN+X1) GOTO 100
      IFN = IFN+1
      P = FUNC(B)
      IF (P .LT. S1) GOTO 10
      IF (SO .GE. S1) GOTO 11
      S0 = P
      X0 = X2
      STEP = A2*STEP
      GO TO 3
   10 X0 = X1
      S0 = S1
      X1 = X2
      S1 = P
      STEP = STEP * A1
      GO TO 3
   11 X0 = X0-X1
      S0 = (S0-S1)*STEP
      P = (P-S1)*X0
      IF (P .EQ. S0) GOTO 18
      STEP = 0.5*(P*X0-S0*STEP)/(P-S0)
      X2 = X1+STEP
      B = BMIN+X2
      IF (B .EQ. BMIN+X1) GOTO 20
      IFN = IFN+1
      P = FUNC(B)
      IF (P .LT. S1) GO TO 19
   18 B = BMIN+X1
      P = S1
      GOTO 20
   19 X1 = X2
   20 STEP = STEP * A2
      GO TO 2
  100 CONTINUE
      START = B
      RETURN
      END
C
      SUBROUTINE STPLT (XL,XU,YL,YU,NX,NY)
      COMMON /SPL/ NR,NC,XL1,YL1,XI,YI,IF,N1
      COMMON /SPLC/ IA(1000), VF
      CHARACTER IA*10, VF*15, FV*9, IT*10
      IF = 0
      IF (NX.LE.0) GO TO 60
      IF (NY.LE.0) GO TO 60
      IF (XL.EQ.XU) GO TO 60
      IF (YL.EQ.YU) GO TO 60
      XL1 = XL
      YL = YL
      NR = NY+1
      NC = NX/10+1
      N1 = NX+11-10*NC
      XI = (XU-XL)/REAL(NX)
      YI = (YU-YL)/REAL(NY)
      DO 10 J = 1,1000
   10   IA(J) = '          '
      XM = MAX(ABS(XL),ABS(XU))
      YM = MAX(ABS(YL),ABS(YU))
      MX = INT(ALOG10(XM)+10.01)-10
      MY = INT(ALOG10(YM)+10.01)-10
      NDX = 2-MX
      NDY = 2-MY
      NDX = MAX(0,MIN(NDX,7))
      NDY = MAX(0,MIN(NDY,5))
      WRITE (FV,1000) NDY
      NR2 = NR+2
      NC2 = NC+2
      IF (NR2*NC2.GT.1000) GO TO 60
      YI1 = YI/2.0
      RF = YL1+REAL(NR)*YI+YI1
      DO 20 I = 1,NR2
        RF = RF-YI
        IF (I.EQ.2.OR.I.EQ.NR2) RF = RF+YI1
        WRITE (IT,FV) RF
        JL = 1+(I-1)*NC2
        JU = I*NC2
        IA(JU) = IT
        IF (I.NE.1.AND.I.NE.NR2) GO TO 20
        IF (I.EQ.1) IT(1:1) = '>'
        IF (I.EQ.NR2) IT(1:1) = '<'
   20   IA(JL) = IT
      WRITE (VF,1010) NDX+1,NDX
      DO 30 I = 1,NR
        J = I*(NC+2)+1
        K = J+NC+1
        IA(J)(10:10) = ':'
   30   IA(K)(1:1) = ':'
      DO 40 I = 1,NC
        JJ = 1+(NC+2)*(NR+1)
        IA(I+1) = '+---------'
   40   IA(I+JJ) = '+---------'
      RETURN
   60 WRITE (6,1020)
      WRITE (6,1030) XL,XU,YL,YU,NX,NY
      IF = 1
      RETURN
 1000 FORMAT ('(1X,F7.',I1,')')
 1010 FORMAT ('(1X,',I1,'X,12F10.',I1,')')
 1020 FORMAT (1X,'FAULT IN PLOT')
 1030 FORMAT (1X,'ARGUMENTS WERE'/1X,4F8.3,2I4)
      END
      SUBROUTINE ADDPT (X,Y,CHR)
      COMMON /SPL/ NR,NC,XL1,YL1,XI,YI,IF,N1
      COMMON /SPLC/ IA,VF
      CHARACTER IA(1000)*10, VF*15,CHR*1
      IF (IF.GT.0) RETURN
      K2 = (Y-YL1)/YI+0.5
      K2 = NR-K2
      IF (K2.LT.1) K2 = 0
      IF (K2.GT.NR) K2 = NR+1
      X2 = (X-XL1)/XI
      K1 = X2+0.5
      IF (K1.LT.0) K1 = -1
      IF (K1.GE.10*NC) K1 = 10*NC
      K11 = (K1+10)/10
      K12 = K1-10*(K11-1)
      K12 = K12+1
      J = (NC+2)*K2+K11+1
      IA(J)(K12:K12) = CHR
      RETURN
      END
      SUBROUTINE ENPLT (M)
      COMMON /SPL/ NR,NC,XL1,YL1,XI,YI,IF,N1
      COMMON /SPLC/ IA,VF
      CHARACTER IA(1000)*10, VF*15, OF*18
      DIMENSION XV(20)
      IF (IF.GT.0) RETURN
      DO 10 I = 1,NC
   10   XV(I) = XL1+REAL(I-1)*XI*10.0
      NR2 = NR+2
      NC2 = NC+2
      OF = '(1X,  A10,A  ,A10)'
      WRITE (OF(5:6),'(I2)') NC
      WRITE (OF(12:13),'(I2)') N1
      WRITE (M,VF) (XV(I),I=1,NC)
      DO 20 I = 1,NR2
        JL = 1+(I-1)*NC2
        JU = I*NC2
   20  WRITE (M,OF) (IA(J),J=JL,JU)
      WRITE (M,VF) (XV(I),I=1,NC)
      RETURN
      END
