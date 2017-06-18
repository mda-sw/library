C       Driver routine for linear regression with errors in both variables.
C       Input data given in DATA statements.
C       Reference: Fasano and Vio (see routine FV).

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
C       Det. std. devs. from wts. input.
	DO I=1,N			                
	   WTX(I) = 1./sqrt(WTX(I))
	   WTY(I) = 1./sqrt(WTY(I))
	ENDDO
C							
C
	ISTAT = 0
	CALL fv(N,X,Y,WTX,WTY,0.01,0,XIN,SL,VAR,SIGMAI,SIGMAS,
     X                                      SIGMAV,RES,NITER)
  90	WRITE(6,1001) N,NITER,SL,XIN,SIGMAS,SIGMAI	
1001	FORMAT(1X,'N =',I4,' Niter =',I5,' Slope =',G16.8,		
     x  ' Intercept =',G16.8,/,' Errors: sigma1(slope) =',G11.4,	
     X  ' sigma1(intercept) =',G11.4)
C
C
	END
C
	subroutine fv(np,x,y,sx,sy,
     !	tol,ich,a,b,v,ea,eb,ev,s,niter)

c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c Fits a straight line to a distribution of data points with
c errors on both coordinates.
c
c np --> Number of data points (input)
c x,y --> Vectors of length np containing the coordinates 
c         of the data points  (input)
c sx,sy --> Vectors of length np containing the standard
c           deviations of the data points (input)
c tol --> Allowed tolerance in estimating the slope (input).
c          The iteration stops when the estimate of the slope
c	   differs from the previous one less than 'tol'. 
c ich --> control index which allows (if ich=0) to avoid the natural
c         variance estimation (input)
c
c a --> intercept (output)
c b --> slope (output)
c v --> natural variance along the Y axis (output)
c ea --> expected error (1 sigma) of the intercept (output)
c eb --> expected error (1 sigma) of the slope (output)
c ev --> expected error (1 sigma) of the natural variance (output)
c s --> weighted r.m.s. of the residuals (output)
c niter --> number of iterations (output)
c
c Authors: G. Fasano and R. Vio, Astronomical Observatory of Padua,
c          Padua, Italy.  
c Reference: Newsletter of Working Group for Modern Astronomical
c          Methodology, Issue 7, Sept. 1988, pp. 2-7.
c          (Newsletter edited and distributed by F. Murtagh and A. Heck.)
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	dimension x(1),y(1),sx(1),sy(1),w(10000)

	niter=1                      ! number of iterations
	anp=float(np)
	vn=0.

c -----------------------------------------------------------------
c    First estimation of the slope and of the intercept by means of 
c    the unweighted regression 

	xsum=0.
	x2sum=0.
	ysum=0.
	xysum=0.

	do i=1,np
	 xsum=xsum+x(i)
	 x2sum=x2sum+x(i)*x(i)
	 ysum=ysum+y(i)
	 xysum=xysum+x(i)*y(i)
	enddo

	delta=anp*x2sum-xsum*xsum
	
c ------------------------------------------------------------------
c   initial guesses of the slope and of the intercept

	b=(anp*xysum-xsum*ysum)/delta	
	a=(x2sum*ysum-xsum*xysum)/delta

c --------------------------------------------------------------------
c ********************************************************************
c    begin the iterative loop

 2	xb=0.
	yb=0.
	b1=0.
	b2=0.
	b3=0.
	s1=0.
	s2=0.
	wsum=0.
	v1=0.
	v2=0.
	v3=0.

c ------------------------------------------------------------------
c   computing the weighting factors

	do i=1,np
	 w(i)=1./(sx(i)*sx(i)*b*b+sy(i)*sy(i)+vn)
 	enddo

c ------------------------------------------------------------------
c   computing the barycenter coordinates

	do i=1,np
	 xb=xb+w(i)*x(i)
	 yb=yb+w(i)*y(i)
	 wsum=wsum+w(i)
	enddo

	xb=xb/wsum
	yb=yb/wsum

c ----------------------------------------------------------------------
c    estimating the slope
c ----------------------------------------------------------------------
c   computing the coefficients of the quadric

	do i=1,np
	 b1=b1+w(i)*w(i)*(sx(i)*sx(i)+vn/(1.+b*b))*(x(i)-xb)*(y(i)-yb)
	 b2=b2+w(i)*w(i)*((x(i)-xb)**2.*(sy(i)*sy(i)+vn/(1.+b*b))
     $   -(y(i)-yb)**2.*(sx(i)*sx(i)+vn/(1.+b*b)))
	 b3=b3-w(i)*w(i)*(sy(i)*sy(i)+vn/(1.+b*b))*(x(i)-xb)*(y(i)-yb)
	enddo

c ----------------------------------------------------------------------
c   if the quadric reduces to a linear form then....

	if(b1.eq.0.) then
	 b=-b3/b2
	 goto 20
	endif

c   ----------------------------------------------------------------------
c   discriminant of the quadric 

	discr=b2*b2-4.*b1*b3      

	if(discr.le.0.) then
	 sqdis=0.
	 goto 30
	endif

	sqdis=sqrt(discr)

c ------------------------------------------------------------------------
c   computing the two solutions of the quadric

 30	b_1=(-b2+sqdis)/(2.*b1)
	b_2=(-b2-sqdis)/(2.*b1)

c -----------------------------------------------------------------------
c   choosing the solution with the lowest weighted r.m.s.

	do i=1,np
	 s1=s1+(y(i)-b_1*x(i)-a)**2./(b_1*b_1*sx(i)*sx(i)+
     !	 sy(i)*sy(i)+vn)	
	 s2=s2+(y(i)-b_2*x(i)-a)**2./(b_2*b_2*sx(i)*sx(i)+
     !	 sy(i)*sy(i)+vn)	
	enddo

	if(s1.le.s2) then
	 b=b_1
	 icon=1
	 s=sqrt(s1/wsum)
	 else
	 b=b_2
	 icon=2
	 s=sqrt(s2/wsum)
	endif

 20     continue

c ---------------------------------------------------------------------
c  estimating the natural variance
c ---------------------------------------------------------------------

	if(ich.eq.0) goto 100
 
 	do i=1,np
	 v1=v1+w(i)*(y(i)-yb-b*(x(i)-xb))**2.
	 v2=v2+w(i)*(sy(i)*sy(i)+b*b*sx(i)*sx(i))
	 v3=v3+w(i)*w(i)
	enddo

	vn=(v1*anp/wsum/(anp-2.)-v2/wsum)
	v=vn

	if(vn.lt.0.) vn=0.

 100	continue
c -------------------------------------------------------------------
c   estimating the intercept

	a=yb-b*xb           

c ----------------------------------------------------------------------
c   testing if the slope variation is within the tolerance

  	if(abs(b-bcont).gt.tol) then
	 niter=niter+1
	 bcont=b
	 goto 2
	endif

c ----------------------------------------------------------------------
c   close the iterative loop
c ***********************************************************************
c ----------------------------------------------------------------------
c   estimating the errors of the parameters
c -----------------------------------------------------------------------
c   inizializations

	e1=0.
	e2=0.
	e3=0.
	e4=0.
	e5=0.
	e6=0.
	an2=anp/(anp-2.)	

	do i=1,np
	 delyi=y(i)-b*x(i)-a
	 sigi=sx(i)*sx(i)*b*b+sy(i)*sy(i)
	 e1=e1+w(i)*(x(i)-xb)**2.
	 e2=e2+w(i)*x(i)*x(i)
	 e3=e3+delyi*delyi*delyi*delyi
	 e4=e4+delyi**2.
	 e5=e5+sigi**2.
	 e6=e6+sigi
	enddo

	e3=e3/anp
	e4=(e4/anp)**2.
	e5=e5/anp
	e6=(e6/anp)**2.

	eb=sqrt(1./e1)            ! error of the slope
	ea=sqrt(e2/(wsum*e1))     ! error of the intercept

	if(ich.eq.0) goto 40
	ev=sqrt((an2**2.*(e3-e4)+e5-e6)/float(np-1))
                                  ! error of the natural variance          
 40	return
	end
