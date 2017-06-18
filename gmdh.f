c This main program reads in the data arrays necessary to find the
C Ivakhnenko polynomial (ip).  Subroutine GMDH finds the IP and stores
c it in the arrays TREE, ITREE, and ITER.  These arrays are then used
c by subsroutine COMP to evaluate the IP at the original data points.
c The user can easily modify this program by either adding a 
c preprocessor subroutine to transform the data before GMDH or by
c adding another read statement to evaluate the IP at other points.
c
c Input data is m = number of variables
c               n = number of data points
c              nt = number of data points in the training set (the first
c                   nt data points are assumed to be in the training set)
c           niter = number of levels GMDH performs before stopping
c                   (if niter = 0 GMDH decides itself)
c           nprnt = 0 if input is not printed
c                 = 1 if input is printed
c              pi = fractional increase in the number of variables
c                   at each iteration (should be between 0 and 1)
c               x = n by m array of independent variables
c               y = n array of dep variables
c
c Order of input: cards 1-10 hollerith info to be printed
c card 11  m, n, nt, pi    format (3i5,f10.2)
c cards 12-  the i'th card contains  y(i), (x(i,j),j=1,m)  i = 1, n
c
	dimension x(100,30), y(100), itree(10,75), tree(10,75,6)
	dimension zz(30), ev(100,30), ysave(100), itr(256), data(80)
c
c       read no. of variables, no. of observations, no. of observations
c       in the training set, rms value
	write (3,4)
	do 60 i = 1, 10
	   read (1,42) data
  42       format (80a1)
	   write (3,45) data
  45       format (10x,80a1)
  60    continue
        read (1,1) m,n,nt,niter,nprnt,pi
   1    format (5i5,5f10.2)
 	do 2 i = 1, n
c          read dep variable, m indep variables
	   read (1,3) y(i), (x(i,j),j=1,m)
   2    ysave(i) = y(i)
   3    format (8f10.2)
   4    format (1h1)
	write (3,43)
  43    format (' --------------------------------------------------------------
     x--------------------------')
	write (3,5) m,n
   5    format (' number of indep var = ',i4/' number of obs = ',i4)
	write (3,6) nt
   6    format (' number of var in training set = ', i4)
	write (3,7) pi
   7    format (' fractional increase in variables = ', f8.3//)
	write (3, 43)
	if (nprnt.eq.0.) goto 32
	write (3, 30)
  30    format (' print input data'//)
	write (3,8)
   8    format (3x,' obs   y     (x(i),i=1,m)'/)
	do 9 i = 1, n
   9       write (3,10) i,y(i),(x(i,j),j=1,m)
  10    format (i6,10f7.2)
        write (3,4)
c       save the data array for later
  32    do 11 i = 1, n
	   do 11 j = 1, m
  11          ev(i,j) = x(i,j)
c       call gmdh subroutine
	call gmdh(n,m,nt,pi,x,y,itree,tree,iter,niter)
	write (3,40)
  40    format (' case no.   observed value      estimate          error
     x            percent error'/)
	do 12 i = 1, n
	   do 13 j = 1, m
  13          zz(j) = ev(i,j)
c          call subroutine to evaluate the Ivakhnenko polynomial
	   call comp(zz,yy,itree,tree,iter,itr)
	   er = ysave(i) - yy
	   perer = 100.*er/ysave(i)
  12    write (3,14) i,ysave(i),yy,er,perer
  14    format (i5,3x,4e18.8)
        write (3,43)
        write (3,4)
c       The Ivakhnenko polynomial is printed only if it is a simple
c       quadratic
        if (iter.gt.1) stop
        write (3,70) (tree(1,1,j), j=1,6)
  70    format(//' Ivakhnenko polynomial '//
     x   ' y = a + b*u + c*v + d*u*u + e*v*v + f*u*v   '/
     x   ' a = ',e12.4/' b = ',e12.4/' c = ',e12.4/
     x   ' d = ',e12.4/' e = ',e12.4/' f = ',e12.4//)
        write (3,71) itr(2),itr(3)
  71    format (' u = x(',i2,' )      v = x(',i2,' )')
        write (3,4)
        stop
        end
c
c-----------------------------------------------------------------------------
	subroutine gmdh(n,m,nt,pi,x,y,itree,tree,iter,niter)
	dimension x(100,30), y(100), tree(10,75,6), itree(10,75)
	dimension poly(6,75), xtx(6,7), xty(6), index(435), work(100,75)
	dimension d(435), ind(435), zzz(6), ma(20), wrk(100)
	dimension xwork(100), ywork(100)
	rms = pi
	ntp1 = nt+ 1
	nc = n - nt
	niter = niter + 1
	mm = m
	iter = 1
	dmin = 1.0e20
  28    l = 1
	mm1 = m-1
c
c       begin loop to compute new variables
c
	do 4 i = 1, mm1
	ip1 = i + 1
	   do 4 j = ip1, m
	      do 5 ii = 1, 6
	   	 xty(ii) = 0.0
	         do 5 jj = 1, 6
   5             xtx(ii,jj) = 0.
              xtx(1,1) = nt
              do 7 k = 1, nt
                 xtx(1,2) = xtx(1,2) + x(k,i)
                 xtx(1,3) = xtx(1,3) + x(k,j)
                 xtx(1,4) = xtx(1,4) + x(k,i)**2
                 xtx(1,5) = xtx(1,5) + x(k,j)**2
                 xtx(1,6) = xtx(1,6) + x(k,i)*x(k,j)
                 xtx(2,2) = xtx(2,2) + x(k,i)**2
                 xtx(2,3) = xtx(2,3) + x(k,i)*x(k,j)
                 xtx(2,4) = xtx(2,4) + x(k,i)**3
                 xtx(2,5) = xtx(2,5) + x(k,i)*x(k,j)**2
                 xtx(2,6) = xtx(2,6) + x(k,i)*x(k,i)*x(k,j)
                 xtx(3,3) = xtx(3,3) + x(k,j)**2
                 xtx(3,4) = xtx(3,4) + x(k,i)*x(k,i)*x(k,j)
                 xtx(3,5) = xtx(3,5) + x(k,j)**3
                 xtx(3,6) = xtx(3,6) + x(k,i)*x(k,j)**2
                 xtx(4,4) = xtx(4,4) + x(k,i)**4
                 xtx(4,5) = xtx(4,5) + (x(k,i)*x(k,j))**2
                 xtx(4,6) = xtx(4,6) + x(k,j)*x(k,i)**3
                 xtx(5,5) = xtx(5,5) + x(k,j)**4
                 xtx(5,6) = xtx(5,6) + x(k,i)*x(k,j)**3
                 xtx(6,6) = xtx(6,6) + (x(k,i)*x(k,j))**2
                 xty(1) = xty(1) + y(k)
                 xty(2) = xty(2) + x(k,i)*y(k)
                 xty(3) = xty(3) + x(k,j)*y(k)
                 xty(4) = xty(4) + x(k,i)*x(k,i)*y(k)
                 xty(5) = xty(5) + x(k,j)*x(k,j)*y(k)
  7              xty(6) = xty(6) + x(k,i)*x(k,j)*y(k)
              do 45 ii = 2, 6
                 im1 = ii-1
                 do 45 jj = 1, im1
  45          xtx(ii,jj) = xtx(jj,ii)
c
c	      solve regression equation for variables i and j
c
	      call sys(xtx,zzz,xty,6,iflag)
              if (iflag.eq.1) goto  41
              do 40 iii = 1, 6
  40          poly(iii,l) = zzz(iii)
              do 15 k = 1, n
                 ww = poly(1,l) + poly(2,l)*x(k,i) + poly(3,l)*x(k,j)
		 ww = ww + poly(4,l)*x(k,i)**2 + poly(5,l)*x(k,j)**2
		 ww = ww + poly(6,l)*x(k,i)*x(k,j)
  15          work(k,l) = ww
              ind(l) = 100*(i+10) + (j+10)
              if (l.eq.75) goto 385
              l = l+1
  41       continue
   4  	continue
c
c       completed construction of m*(m-1)/2 new variables
c
	l = l-1
 385    do 120 i = 1, nc
 120    ywork(i) = y(nt+i)
        do 122 j = 1,l
           do 124 i = 1, nc
 124       xwork(i) = work(nt+i,j)
c
c       compute the goodness of fit statistics
c
	call stat(nc,ywork,xwork,st)
	d(j) = st
 122	index(j) = j
c
c	sort the values of the statistics from low to high
c
	call sort(d,l,index)
	m = m + rms*m
	if (m.gt.l) m = l
  53 	continue
c	the largest number of var is set to 75
	if (m.gt.75) m = 75
	if (m.lt.mm) m = mm
	do 22 j = 1, m
	   itree(iter,j) = ind(index(j))
	   do 23 k = 1, 6
  23       tree(iter,j,k) = poly(k,index(j))
  22    continue
c
c       test for convergence of gmdh algorithm
c
	if (niter.eq.1) goto 55
	if (iter.eq.niter) goto 60
	goto 56
  55    test = d(1) - dmin + .0005
	if (test.ge.0.) goto 60
  56    write (3,300)
 300    format(' ----------------------------------------------------------
     x-------------------------')
	write  (3,200) iter
 200 	format (' Level number = ',i3)
	write  (3,203) m,d(1)
 203 	format (' No. variables saved = ',i3/' rmin value (summed ',
     x	'over checking set) = ',e12.5)
	dmin = d(1)
	ma(iter) = m
	iter = iter + 1
c  	max number of iteration is I1
	do 26 i = 1, n
	   do 26 j = 1, m
  26	x(i,j) = work(i,index(j))
	sum = 0.
	do 190 i = 1, nt
 190	sum = sum + y(i)
	sum = sum/nt
	sum1 = 0.
	sum2 = 0.
	do 191 i = 1, nt
	   sum1 = sum1 + (sum - x(i,1))**2
 191	sum2 = sum2 + (y(i) - sum)**2
	sum = sum1/sum2
	write (3,92) sum
	write (3,300)
	goto 28
  60	if (iter.eq.1) goto 65
	iter = iter - 1
  65	m = mm
	write  (3,130) iter
 130	format ( /' GMDH converged after ',i3,' generation(s)')
	sum = 0.
	do 90 i = 1, nt
  90	sum = sum + y(i)
	sum = sum/nt
	sum1 = 0.
	sum2 = 0.
	do 91 i = 1, nt
	   sum1 = sum1 + (sum-x(i,1))**2
  91	sum2 = sum2 + (y(i) - sum)**2
	sum = sum1/sum2
	write  (3, 92) sum
  92	format (' Multiple correlation (summed over training set) = '
     x	,f6.2)
	write (3,300)
	return
	end
c---------------------------------------------------------------------------
	subroutine sort(a,la,ir)
	integer la,ir(la)
	real a(la)
	integer iu(21),il(21),i,m,j,k,ij,it,l,itt
	real t,tt,r
	if (la.le.0) return
	do 5 i = 1, la
	   if (a(i).lt.0) a(i) = -a(i)
  5	continue
	m = 1
	i = 1
	j = la
	r = .375
  10	if (i.eq.j) goto 55
  15	if (r.gt..5898437) goto 20
	r = r+3.90625e-2
	goto 25
  20	r = r - .21875
  25 	k = i
c
	ij = i + (j-i)*r
	t = a(ij)
	it = ir(ij)
c
	if (a(i).le.t) goto 30
	a(ij) = a(i)
	a(i) = t
	t = a(ij)
	ir(ij) = ir(i)
	ir(i) = it
	it = ir(ij)
  30	l = j
c
	if (a(j).ge.t) goto 40
	a(ij) = a(j)
	a(j) = t
	t = a(ij)
	ir(ij) = ir(j)
	ir(j) = it
	it = ir(ij)
c
	if (a(i).le.t) goto 40
	a(ij) = a(i)
	a(i) = t
	t = a(ij)
	ir(ij) = ir(i)
	ir(i) = it
	it = ir(ij)
	goto 40
 35	if (a(l).eq.a(k)) goto 40
	tt = a(l)
	a(l) = a(k)
	a(k) = tt
	itt = ir(l)
	ir(l) = ir(k)
	ir(k) = itt
c
  40	l = l-1
	if (a(l).gt.t) goto 40
c
  45	k = k+1
	if (a(k).lt.t) goto 45
c
	if (k.le.l) goto 35
c
	if (l-i.le.j-k) goto 50
	 il(m) = i
	iu(m) = l
	i = k
	m = m+1
	goto 60
  50	il(m) = k
	iu(m) = j
	j = l
	m = m + 1
	goto 60
c
  55 	m = m -1
	 if(m.eq.0) return
	i = il(m)
	j = iu(m)
  60	if (j-i.ge.11) goto 25
	if (i.eq.1) goto 10
	i = i-1
  65	i = i+1
	if (i.eq.j) goto 55
	t = a(i+1)
	it = ir(i+1)
	if (a(i).le.t) goto 65
	k = i
  70	a(k+1) = a(k)
	ir(k+1) = ir(k)
	k=k-1
	if (t.lt.a(k)) goto 70
	a(k+1) = t
	ir(k+1) = it
	goto 65
	end
c--------------------------------------------------------------------------
	subroutine sys(a,x,b,n,iflag)
	dimension a(6,7),x(6),index(6),b(6)
	np1 = n + 1
	do 1 i = 1, n
	   index(i) = i
  1	a(i,np1) = b(i)
	do 5 k = 1, n
	   if (k.eq.n) goto 2
	   call inter(a,n,k,index)
  2	   kp1 = k + 1
	   do 4 j = kp1, np1
	      if (abs(a(k,k)).lt..000001) goto 8
	      a(k,j) = a(k,j)/a(k,k)
	      do 4 i = 1, n
	         if (k-i) 3,4,3
  3	         a(i,j) = a(i,j) -a(i,k)*a(k,j)
  4	      continue
  5	continue
	do 7 i = 1,n
	   do 6 j = 1, n
	      if (index(j).ne.i) goto 6
	      x(i) = a(j,np1)
	      goto 7
  6	   continue
  7	continue
	iflag = 0
	return
  8	iflag = 1
	return 
	end
c---------------------------------------------------------------------------
	subroutine inter(a,n,k,index)
	dimension a(6,7),index(6)
	np1 = n + 1
	nr = k
	nc = k
	ab = abs(a(k,k))
	do 2 i = k, n
	   do 2 j = k, n
	      if (abs(a(i,j))-ab) 2,2,1
  1	      ab = abs(a(i,j))
  	      nr = i
	      nc = j
  2	continue
	do 3 j = k, np1
	   de = a(nr,j)
	   a(nr,j) = a(k,j)
  3	a(k,j) = de
	do 4 i = 1, n
	   de = a(i,nc)
	   a(i,nc) = a(i,k)
  4	a(i,k) = de
	is = index(nc)
	index(nc) = index(k)
	index(k) = is
	return
	end
c-------------------------------------------------------------------------
	subroutine comp(x,y,itree,tree,iter,itr)
	dimension x(30),itree(10,75),tree(10,75,6),work(256)
	dimension itr(256)
	it = iter
	itr(1) = 1
	i = 1
  6 	l = 0
	nn = 2**(i-1)
 	n1 = 2**i
	nz = 2**(i+1)-1
	j = n1
  4	jj = itr(nn+l)
	xx = itree(iter,jj)
	itr(j) = itree(iter,jj)/100-10
	iz = itree(iter,jj)/100
	itr(j+1) = xx -100*iz - 10
	j = j + 2
	if (j.gt.nz) goto 3
	l = l + 1
	goto 4
  3	if (iter.eq.1) goto 5
	iter = iter -1
	i = i + 1
	goto 6
  5	iter = it
	nz = 2**(iter-1)
	nzz = nz
	n1 = 2**iter
	do 8 j = 1, nzz
	   jj1 = itr(nz)
	   jj2 = itr(n1)
	   jj3 = itr(n1+1)
	   wk = tree(1,jj1,1) + tree(1,jj1,2)*x(jj2) + 
     x                                           tree(1,jj1,3)*x(jj3)
	   wk = wk + tree(1,jj1,4)*x(jj2)**2 + tree(1,jj1,5)*x(jj3)**2
	   wk = wk + tree(1,jj1,6)*x(jj2)*x(jj3)
	   work(j) = wk
	   nz = nz + 1
  8	n1 = n1 + 2
	iter = iter -1
	if (iter.ne.0) goto 12
	y = work(1)
	iter = it
 	return
  12    i = 2
  18    nz = 2**(iter-1)
  	n1 = 2**iter
	nzz = nz
	n11 = n1
	do 13 j = 1, nzz
	   jj = 2*j-1
	   jj1 = itr(nz)
	   jj2 = itr(n1)
	   jj3 = itr(n1+1)
	   wk = tree(i,jj1,1) + tree(i,jj1,2)*work(jj)
	   wk = wk + tree(i,jj1,3)*work(jj+1) + tree(i,jj1,4)*work(jj)**2
	   wk = wk + tree(i,jj1,5)*work(jj+1)**2
	   wk = wk + tree(i,jj1,6)*work(jj)*work(jj+1)
	   work(n11+j) = wk
	   nz = nz + 1
  13    n1 = n1 + 2
	iter = iter - 1
        if (iter.eq.0) goto 15
	do 14 l = 1, nzz
  14	work(l) = work(n11+l)
	i = i + 1
	goto 18
  15	y = work(3)
	iter = it
	return
        end
c---------------------------------------------------------------------------
	subroutine stat(n,y,x,rms)
	dimension y(50),x(50)
	ss = 0.
	sss = 0.
	do 5 i = 1, n
	   ss = ss + (y(i)-x(i))**2
  5	sss = sss + y(i)**2
	rms = ss/sss
	return
	end
