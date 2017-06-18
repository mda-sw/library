C       Driver routine for linear regression with errors in both variables.
C       Input data given in DATA statements.
C       Reference: York (see routine REGRWT).

	REAL*4		X(5000),Y(5000),WTX(5000),WTY(5000),W(5000)
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
	ISTAT = 0
	CALL REGRWT(N,X,Y,WTX,WTY,W,ITER,SL,XIN,SIGMA,SIGMAI,
     X                                           XBAR,YBAR,ISTAT)
        IF (ISTAT.NE.0) THEN
	   WRITE(6,*) ' ON RETURN FROM REGRWT, ISTAT =',ISTAT
	   GOTO 9000
	ENDIF
  90	WRITE(6,1001) N,ITER,SL,XIN,SIGMA,SIGMAI,XBAR,YBAR	
1001	FORMAT(1X,'N =',I4,' Niter =',I5,' Slope =',G16.8,		
     x  ' Intercept =',G16.8,/,' Errors: sigma(slope) =',G11.4,	
     X  ' sigma(intercept) =',G11.4,/,' Centroid =',G11.4,',',G11.4)
C
C
 9000	CONTINUE
	END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C.IDENTIFICATION
C
C  Program REGRWT.FOR
C  
C.AUTHOR
C
C  F. Murtagh, ST-ECF, Garching.                   Version 1.0 9-July-1986
C  (Following discussions with J. Melnick and L. Lucy.)
C
C.KEYWORDS
C
C  Least squares fitting, linear regression.
C
C.INPUT/OUTPUT
C
C.ALGORITHM
C
C  Error weights in both dep. and indep. vbes.
C  Following initial guess at slope (not using weights), iterative
C  improvement of slope is carried out until convergence.
C  
C  Reference: D. York, "Least squares fitting of a straight line",
C                      Canadian Journal of Physics, 44, 1079-1086,
C                      1966.
C  See also:  M. Lybanon, "A better least squares method when both
C                      variables have uncertainties", Americal 
C                      Journal of Physics, 52, 22-26, 1984.
C  and:       F.S. Acton, "Analysis of Straight Line Data", Dover,
C                      New York, 1959, Chapter 5.
C
C.INPUTS
C
C  N          ... dimension of X, Y, WTX, WTY, W.
C  X, Y       ... vectors of independent and dependent variables.
C  WTX, WTY   ... vectors of weights associated with X and Y.
C  W          ... work vector.
C 
C.OUTPUTS
C
C  ITER       ... number of iterations taked by the optimisation algorithm.
C  SL1        ... slope of least squares solution.
C  XIN        ... intercept of least squares solution.
C  SIGMA      ... uncertainty in slope (sigma).
C  SIGMAI     ... uncertainty in intercept (sigma).
C  XBAR, YBAR ... centroid of points.
C  ISTAT      ... error indicator (= 1 if no convergence; normally 0).
C
C---------------------------------------------------------------------
	SUBROUTINE      REGRWT(N,X,Y,WTX,WTY,W,ITER,SL1,XIN,
     X                         SIGMA,SIGMAI,XBAR,YBAR,ISTAT)
	REAL*4		X(N),Y(N),WTX(N),WTY(N),W(N)
	REAL*4          SL0,SL1,XIN,SIGMA,SIGMAI,XBAR,YBAR,A1,A2,B1,B2,
     X                  C1,SUM,PHI,RPHI,U,V,SUMI,EPS
	DATA            EPS /1.0E-10/
C
C
  40	CALL SOLVE(X,Y,WTX,WTY,N,SL0)			! Initial guess
C							! for the slope, SL0.
	ITER= 0
C
C Reference for the following: York, Can. J. Physics, 1966.
C In the following, the "least squares cubic" is solved iteratively.
C
  55	DO I = 1, N
	   W(I) = WTX(I)*WTY(I)/(WTX(I)+SL0*SL0*WTY(I))   
	ENDDO		
C
	XBAR = 0.0
	YBAR = 0.0
	SUM  = 0.0
	DO I = 1, N
	   XBAR = XBAR + W(I)*X(I)
	   YBAR = YBAR + W(I)*Y(I)
	   SUM  = SUM  + W(I)
	ENDDO
	XBAR = XBAR/SUM
	YBAR = YBAR/SUM
C
C In York's notation, A1, B1 and C1 here become, respectively,
C alpha, beta and gamma.
C
	A1 = 0.0
	A2 = 0.0
	B1 = 0.0
	B2 = 0.0
	C1 = 0.0
	DO I = 1, N
	   U = X(I) - XBAR
	   V = Y(I) - YBAR
	   A1 = A1 + W(I)*W(I)*U*V/WTX(I)
	   A2 = A2 + W(I)*W(I)*U*U/WTX(I)
	   B1 = B1 + W(I)*W(I)*V*V/WTX(I)
	   B2 = B2 + W(I)*U*U
	   C1 = C1 + W(I)*U*V
	ENDDO
	A1 = 2.0*A1/(3.0*A2)
	B1 = (B1 - B2)/(3.0*A2)
	C1 = - C1/A2
C
	IF ((A1*A1 - B).LE.EPS) THEN
	   ISTAT = 2
	   RETURN
	ENDIF
	PHI = (A1*A1*A1 - 1.5*A1*B1 + 0.5*C1)/
     X                          ((A1*A1 - B1)**1.5)
C       RPHI is angle given by arccos of the above expr. (in radians).
	RPHI = ACOS(PHI)
	RPHI =  RPHI/3.0 + 4.0*3.1415926/3.0 
C
	SL1 = A1 + 2.0*SQRT(A1*A1 - B1)*COS(RPHI)	
C
	ITER = ITER+1
	IF(ABS(SL1-SL0).GT.EPS) THEN	
	   SL0 = SL1						
	   IF (ITER.GT.100) THEN
	      ISTAT = 1
	      RETURN
	   ENDIF
	   GOTO 55						
	ENDIF
C
C  Intercept.
C
	XIN = YBAR - SL1*XBAR
C
C  Errors.
C
	SIGMA = 0.0
	SIGMAI = 0.0
	SUM = 0.0
	SUMI = 0.0
	DO I = 1, N
	   U = X(I) - XBAR
	   V = Y(I) - YBAR
	   SIGMA = SIGMA + W(I)*(SL1*U-V)*(SL1*U-V)	   
	   SUM = SUM + W(I)*U*U
	   SIGMAI = SIGMAI + W(I)*X(I)*X(I)
	   SUMI = SUMI + W(I)
	ENDDO
	SIGMA = SIGMA/(FLOAT(N-2)*SUM)
	SIGMAI = SIGMA*SIGMAI/SUMI
	SIGMA = SQRT(SIGMA)
	SIGMAI = SQRT(SIGMAI)
C
C
	RETURN
	END
C----------------------------------------------------------------
	SUBROUTINE SOLVE(X,Y,ERRX,ERRY,N,SL1)
C
C Get estimate of slope, by ignoring weights;
C Equation (10), p. 135, of Acton, F.S., "Analysis of Straight-
C Line Data", Dover, 1959, is used.
C
	REAL*4		X(1),Y(1),ERRX(1),ERRY(1)
	REAL*4		SIGMA(5),AVE(3)
	REAL*8		W(30),RMS
C
	EQUIVALENCE	(SIGMA(1),SXX),(SIGMA(2),SYY),(SIGMA(3),SXY),
     X                  (SIGMA(4),SU) ,(SIGMA(5),SV)
C
C	
  	CALL W_MEAN(X,Y,N,ERRX,ERRY,AVE,SIGMA)		
	RHO = 0.0
	R  = RHO*SQRT(SV/SU)					
	R2 = SV/SU						
	A  = SXY - R*SXX					
	B  = SYY - R2*SXX				
	C  = -R2*SXY+R*SYY			
	SL1 = 0.5*(B+SQRT(B*B-4.*A*C))/A		! Slope
C	XIN = AVE(4) - SL1*AVE(3)			! Intercept
C
	RETURN
	END
C-----------------------------------------------------------------
	SUBROUTINE W_MEAN(X,Y,N,ERRX,ERRY,AVE,SIGMA)
	REAL*4		X(1),Y(1),ERRX(1),ERRY(1)
	REAL*4		AVE(1),SIGMA(1)
	IF(N.LT.2.) RETURN
C
C
	AVEX = 0.
	AVEY = 0.
	WAVX = 0.
	WAVY = 0.
	WX = 0.
	WY = 0.
	SU = 0.
	SV = 0.
	DO I=1,N
	   WTX = 1./ERRX(I)**2
	   WTY = 1./ERRY(I)**2
	   WX = WX + WTX 
	   WY = WY + WTY
	   WAVX = WAVX + X(I)*WTX
	   WAVY = WAVY + Y(I)*WTY
	   AVEX = AVEX + X(I)
	   AVEY = AVEY + Y(I)
	   SU = SU + ERRX(I)*ERRX(I)
	   SV = SV + ERRY(I)*ERRY(I)
	ENDDO
	AVEX = AVEX/N
	AVEY = AVEY/N
	WAVX = WAVX/WX
	WAVY = WAVY/WY
C
	SXX = 0.
	SYY = 0.
	SXY = 0.
	DO I=1,N
	   SXX = SXX + (X(I)-AVEX)*(X(I)-AVEX)
	   SYY = SYY + (Y(I)-AVEY)*(Y(I)-AVEY)
	   SXY = SXY + (X(I)-AVEX)*(Y(I)-AVEY)
	ENDDO
C
C
	AVE(1) = AVEX
	AVE(2) = AVEY
	AVE(3) = WAVX
	AVE(4) = WAVY
	SIGMA(1) = SXX
	SIGMA(2) = SYY
	SIGMA(3) = SXY
	SIGMA(4) = SU/(N-1.)
	SIGMA(5) = SV/(N-1.)
	RETURN
	END
