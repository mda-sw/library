      dimension x(150,4), y(150,2), dstar(150,150), d(150,150)
      open (unit=21, status='old', file='iris.dat')
c
c--   constants for the problem
      n = 150
      m = 4
      p = 2
      maxit = 100
c
c--   get input data
      do i = 1, n
         read (21,100) (x(i,k),k=1,m)
      enddo
 100  format(4f4.1)
c
c--   get distances
      do i = 1, n
         do j = 1, n
            dstar(i,j) = 0.0
            do k = 1, m
               dstar(i,j) = dstar(i,j) + (x(i,k)-x(j,k))**2
            enddo
            dstar(i,j) = sqrt(dstar(i,j))
         enddo
      enddo
c
c--   choose initial y-configuration 
      do k = 1, p
         do i = 1, n
            y(i,k) = float(i)+float(k)/float(p+1)
         enddo
      enddo
c
c--   determine distances in new space
      do i = 1, n
         do j = 1, n
            d(i,j) = 0.0
            do k = 1, p
               d(i,j) = d(i,j) + 
     .                   (y(i,k)-y(j,k))*(y(i,k)-y(j,k))
            enddo
            d(i,j) = sqrt(d(i,j))
         enddo
      enddo
c
c--   now off to subr.
      tol = 0.002
      call smap(n,m,p,x,y,dstar,d,maxit,iter,tol,err)
c
c--
      write(6,*) ' err, iter:', err, iter
      end

c-----------------------------------------------------------------------------
      subroutine smap(n,m,p,x,y,dstar,d,maxit,iter,tol,err)
      dimension x(n,m), y(n,p), dstar(n,n), d(n,n)
c
      iter = 0
      alpha = 0.5
c
      c = 0.0
      do i = 1, n-1
         do j = i+1, n
            c = c + dstar(i,j)
         enddo
      enddo
c
 100  continue
      iter = iter + 1

      do i = 1, n
         do k = 1, p
c
            dedy   = 0.0
            d2edy2 = 0.0
            do j   = 1, n
c
               if (j.eq.i) goto 200
               if (dstar(i,j).le.0.001) goto 200
               if (d(i,j).le.0.001) goto 200
c
               xnumer1 = y(i,k)-y(j,k)
               xnumsq  = xnumer1*xnumer1
               xnumer2 = dstar(i,j)-d(i,j)
               denom1 = dstar(i,j)*d(i,j)
c
               dedy = dedy + xnumer1*xnumer2/denom1
c
               d2edy2 = d2edy2 + (1.0/denom1)*
     .            (xnumer2-(xnumsq/d(i,j))*(1.0+xnumer2/d(i,j)))
c
 200           continue
            enddo
            dedy = -dedy
            d2edy2 = abs(-d2edy2)
            deltay = dedy/d2edy2
c
            y(i,k) = y(i,k) - alpha*deltay
         enddo
       enddo
c
       do i = 1, n-1
          do j = i+1, n
             d(i,j) = 0.0
             do k = 1, p
                d(i,j) = d(i,j) + (y(i,k)-y(j,k))*(y(i,k)-y(j,k))
             enddo
             d(i,j) = sqrt(d(i,j))
             d(j,i) = d(i,j)
          enddo
       enddo
c
       err = 0.0
       do i = 1, n-1
          do j = i+1, n
             if (dstar(i,j).gt.0.001)
     .         err = err + 
     .         ((dstar(i,j)-d(i,j))*(dstar(i,j)-d(i,j)))/dstar(i,j)
          enddo
       enddo
       err = err/c
c
       write(6,*) ' Error: ',err
       if (err.le.tol) return
       if (iter.ge.maxit) return
       goto 100
c
       end
