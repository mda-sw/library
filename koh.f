      dimension a(150,4)
c      dimension grid(3,3,4), nfreq(3,3)
      dimension grid(4,4,4), nfreq(4,4)
      dimension mtchx(150), mtchy(150)
c      dimension ioutrep(3,3), xmax(4), xmin(4)
      dimension ioutrep(4,4), xmax(4), xmin(4)

      open (unit=21,status='old',file='iris_norm.dat')

      n = 150
      iseq = 0
      do j = 1, 4
         xmax(j) = -10000000.0
         xmin(j) =  10000000.0
      enddo


      do i = 1, n
         read (21, 100) (a(i,j), j=1,4)

         do j = 1, 4
            if (a(i,j).gt.xmax(j)) xmax(j) = a(i,j)
            if (a(i,j).lt.xmin(j)) xmin(j) = a(i,j)
         enddo

      enddo
 100  format(4f5.2)


      write (6,*) ' Maxima:'
      write (6,*) (xmax(j),j=1,4)
      write (6,*) ' Minima:'
      write (6,*) (xmin(j),j=1,4)



      write (6,*) ' No. of cases:',n

      write (6,*) ' Test -- sample of data array to be processed:'
      do i = 1, 10
         write (6,*) (a(i,j),j=1,4)
      enddo


c     set parameter values
      m  = 4
      nx = 4
      ny = 4
c      ix = 33731
      ix = 17
      nsteps = 50 * nx * ny 
c      nepoks = nsteps/n
      nepoks = 60
      rmax = max(nx,ny)
      rmin = 1.0
      do i = 1, nx
         do j = 1, ny
            ioutrep(i,j) = 0
         enddo
      enddo

c     initialize grid

      do i = 1, nx
         do j = 1, ny
            do k = 1, m
               grid(i,j,k) = ran(ix)
            enddo
         enddo
      enddo

c     for each object in each epoch...
 
      iter = 0
      discrep = 0.0
      do ne = 1, nepoks
         write(6,*) ' Starting epoch -- discrep',ne,' -- ',discrep
         discrep = 0.0
      do i = 1, n
         iter = iter + 1
         if (iter.eq.1000) write(6,*) ' Iter 1000'
         if (iter.eq.5000) write(6,*) ' Iter 5000'
         if (iter.eq.10000) write(6,*) ' Iter 10000'
         if (iter.eq.15000) write(6,*) ' Iter 15000'

c        det best match grid unit

         ihit = 0
         jhit = 0
         dhit = 100000.0
         do ig = 1, nx  
            do jg = 1, ny
               d = 0.0
               neff = 0
               do k = 1, m
                  d = d + (a(i,k)-grid(ig,jg,k))**2
               enddo
               d = d/float(m)
               if (d.lt.dhit) then
                  dhit = d
                  ihit = ig
                  jhit = jg
               endif
            enddo
         enddo
         discrep = discrep + dhit

c        we now have hit - i.e. best match grid unit

c        define neighborhood: center = (ihit, jhit), radius = ??, 
c        max limits: 1 <= ihit <= nx, 1 <= jhit <= ny

c        max radius = max(nx,ny)
c        min radius = 1
c        shrinks linearly over 1000 iterations

c        note following enlarges rmax a bit
         ir = max(rmax*float(1001-iter)/1000.0 + 0.9999999999,1)

         alpha = max(0.9*(1.0-float(iter)/1000.0),0.01)

         do in = ihit-ir, ihit+ir
            do jn = jhit-ir, jhit+ir
               if (in.lt.1.or.in.gt.nx) goto 700
               if (jn.lt.1.or.jn.gt.ny) goto 700

               do k = 1, m
                  grid(in,jn,k) = grid(in,jn,k) + alpha*
     x            (a(i,k)-grid(in,jn,k))
               enddo

 700           continue

           enddo
        enddo

c       that's it

      enddo

      enddo

c     final best match        
               
      do i = 1, n
         ihit = 0
         jhit = 0
         dhit = 100000.0
         do ig = 1, nx  
            do jg = 1, ny
               d = 0.0
               neff = 0
               do k = 1, m
                  d = d + (a(i,k)-grid(ig,jg,k))**2
               enddo
               d = d/float(m)
               if (d.lt.dhit) then
                  dhit = d
                  ihit = ig
                  jhit = jg
               endif
            enddo
         enddo

         mtchx(i) = ihit
         mtchy(i) = jhit
         ioutrep(ihit,jhit) = i

      enddo


c     count # spectra associated with nodes
      do i = 1, nx
         do j = 1, ny
            nfreq(i,j) = 0
            do k = 1, n
               if (mtchx(k).eq.i.and.mtchy(k).eq.j) 
     .              nfreq(i,j) = nfreq(i,j) + 1
            enddo
         enddo
      enddo

      do i = 1, n
         write(22,*) i, mtchx(i), mtchy(i)
      enddo


c     finito


      write(6,*) ' Cells assigned to all cases to unit 22.'
      write(6,*) ' Grid cell freqs. and unit/vector values to unit 24.'
      do i = 1, nx
         do j = 1, ny
            ij = j + (i-1)*nx
            write(24,780) i, j, ij, nfreq(i,j), (grid(i,j,k),k=1,4)
         enddo
      enddo
 780  format(2i3,i5,i8,4f5.2)


      end

c----------------------------------------------------------------------------


